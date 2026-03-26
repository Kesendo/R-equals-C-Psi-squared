using System.Numerics;

namespace RCPsiSquared.Propagate;

/// <summary>
/// Matrix-free RK4 propagation of the Lindblad master equation on raw Complex[] arrays.
/// Uses accumulating RK4 to minimize work buffers: 3 internal arrays + rho = 4 total.
///
/// Memory budget at N=15 (d=32768):
///   rho:      16 GB (passed in)
///   _accum:   16 GB
///   _kBuf:    16 GB
///   _tmp:     16 GB
///   _dephMask: 8 GB (double[])
///   Total:    ~72 GB (fits in 128 GB)
/// </summary>
public sealed class MatrixFreePropagator
{
    private readonly MatrixFreeHamiltonian _H;
    private readonly int _d;
    private readonly int _d2;
    private readonly double[] _dephMask;

    private readonly Complex[] _accum;
    private readonly Complex[] _kBuf;
    private readonly Complex[] _tmp;

    public MatrixFreePropagator(MatrixFreeHamiltonian H, double[] gammas, int nQubits)
    {
        _H = H;
        _d = 1 << nQubits;
        _d2 = _d * _d;

        Console.Error.WriteLine($"MatrixFreePropagator: d={_d}, arrays={_d2:N0} elements ({(long)_d2 * 16 / 1e9:F1} GB each)");

        // Precompute dephasing mask (column-major)
        _dephMask = new double[_d2];
        Parallel.For(0, _d, j =>
        {
            int jd = j * _d;
            for (int i = 0; i < _d; i++)
            {
                double rate = 0;
                int xor = i ^ j;
                for (int k = 0; k < nQubits; k++)
                {
                    if (((xor >> (nQubits - 1 - k)) & 1) != 0)
                        rate -= 2 * gammas[k];
                }
                _dephMask[jd + i] = rate;
            }
        });

        Console.Error.WriteLine("Allocating work buffers...");
        _accum = new Complex[_d2];
        _kBuf = new Complex[_d2];
        _tmp = new Complex[_d2];
        Console.Error.WriteLine("MatrixFreePropagator ready.");
    }

    /// <summary>
    /// Evaluate drho/dt = -i[H, rho] + dephasing.
    /// Reads from input, writes to output. No extra allocations.
    /// </summary>
    public void EvalRHS(Complex[] input, Complex[] output)
    {
        // Step 1: output = [H, input]
        _H.ApplyCommutator(input, output);

        // Step 2: output = -i * [H, input] + deph * input
        int d2 = _d2;
        var mask = _dephMask;
        Parallel.For(0, d2, idx =>
        {
            var comm = output[idx];
            // -i * comm = (comm.Imag, -comm.Real)
            output[idx] = new Complex(comm.Imaginary, -comm.Real)
                        + mask[idx] * input[idx];
        });
    }

    /// <summary>
    /// dest = a + scale * b, element-wise.
    /// </summary>
    private void AddScaled(Complex[] a, double scale, Complex[] b, Complex[] dest)
    {
        int d2 = _d2;
        Parallel.For(0, d2, i =>
        {
            dest[i] = a[i] + scale * b[i];
        });
    }

    /// <summary>
    /// a += scale * b, in-place.
    /// </summary>
    private void AddScaledInPlace(Complex[] a, double scale, Complex[] b)
    {
        int d2 = _d2;
        Parallel.For(0, d2, i =>
        {
            a[i] += scale * b[i];
        });
    }

    /// <summary>
    /// dest = scale * src.
    /// </summary>
    private void CopyScaled(Complex[] dest, double scale, Complex[] src)
    {
        int d2 = _d2;
        Parallel.For(0, d2, i =>
        {
            dest[i] = scale * src[i];
        });
    }

    /// <summary>
    /// Enforce Hermiticity: rho = (rho + rho^dag) / 2.
    /// </summary>
    private void EnforceHermiticity(Complex[] rho)
    {
        int d = _d;
        Parallel.For(0, d, i =>
        {
            // Diagonal: keep real part only
            rho[i * d + i] = new Complex(rho[i * d + i].Real, 0);

            for (int j = i + 1; j < d; j++)
            {
                // rho[i,j] = buf[j*d + i], rho[j,i] = buf[i*d + j]
                var ij = rho[j * d + i];
                var ji = rho[i * d + j];
                var avg = (ij + Complex.Conjugate(ji)) * 0.5;
                rho[j * d + i] = avg;
                rho[i * d + j] = Complex.Conjugate(avg);
            }
        });
    }

    /// <summary>
    /// Accumulating RK4 integration.
    /// Instead of storing k1-k4 separately (4 buffers),
    /// we accumulate the weighted sum incrementally (1 buffer).
    /// </summary>
    public void Propagate(Complex[] rho, double tMax, double dt,
        double[] measureTimes, Action<double, Complex[]> callback)
    {
        int measureIdx = 0;
        double t = 0;

        if (measureIdx < measureTimes.Length && measureTimes[measureIdx] <= dt / 2)
        {
            callback(0, rho);
            measureIdx++;
        }

        double dt2 = dt / 2;
        int nSteps = (int)(tMax / dt) + 1;

        for (int step = 0; step < nSteps; step++)
        {
            // Stage 1: k = f(rho), accum = (1/6) * k
            EvalRHS(rho, _kBuf);
            CopyScaled(_accum, 1.0 / 6.0, _kBuf);

            // Stage 2: tmp = rho + dt/2 * k, k = f(tmp), accum += (2/6) * k
            AddScaled(rho, dt2, _kBuf, _tmp);
            EvalRHS(_tmp, _kBuf);
            AddScaledInPlace(_accum, 2.0 / 6.0, _kBuf);

            // Stage 3: tmp = rho + dt/2 * k, k = f(tmp), accum += (2/6) * k
            AddScaled(rho, dt2, _kBuf, _tmp);
            EvalRHS(_tmp, _kBuf);
            AddScaledInPlace(_accum, 2.0 / 6.0, _kBuf);

            // Stage 4: tmp = rho + dt * k, k = f(tmp), accum += (1/6) * k
            AddScaled(rho, dt, _kBuf, _tmp);
            EvalRHS(_tmp, _kBuf);
            AddScaledInPlace(_accum, 1.0 / 6.0, _kBuf);

            // Final: rho += dt * accum
            AddScaledInPlace(rho, dt, _accum);
            EnforceHermiticity(rho);

            t = (step + 1) * dt;

            while (measureIdx < measureTimes.Length &&
                   Math.Abs(t - measureTimes[measureIdx]) < dt / 2)
            {
                callback(t, rho);
                measureIdx++;
            }
        }
    }
}
