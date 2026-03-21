using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Propagate;

/// <summary>
/// RK4 propagation of the Lindblad master equation with zero-allocation hot path.
/// drho/dt = -i[H, rho] + sum_k gamma_k * (Z_k rho Z_k - rho)
///
/// Z-dephasing dissipator is element-wise (no matrix multiply).
/// Only the commutator [H, rho] requires BLAS zgemm (2 calls per evaluation).
/// All intermediate matrices are pre-allocated in the constructor.
/// </summary>
public class LindbladPropagator
{
    private readonly Matrix<Complex> _H;
    private readonly int _d;
    private readonly int _d2; // d*d

    // Precomputed dephasing mask as flat array for cache-friendly access
    // mask[i*d+j] = -2 * sum_{k where bit_k(i) != bit_k(j)} gamma_k
    private readonly double[] _dephMask;

    // Pre-allocated work buffers (created once, reused every step)
    private readonly Matrix<Complex> _Hrho;   // H * rho
    private readonly Matrix<Complex> _rhoH;   // rho * H
    private readonly Matrix<Complex> _drho;   // EvalRHS output
    private readonly Matrix<Complex> _k1, _k2, _k3, _k4;
    private readonly Matrix<Complex> _tmp;    // scratch for rho + dt/2 * k etc.

    public LindbladPropagator(Matrix<Complex> H, double[] gammas, int nQubits)
    {
        _H = H;
        _d = 1 << nQubits;
        _d2 = _d * _d;

        // Precompute dephasing mask (flat array, column-major to match MathNet storage)
        _dephMask = new double[_d2];
        for (int j = 0; j < _d; j++)
            for (int i = 0; i < _d; i++)
            {
                double rate = 0;
                int xor = i ^ j;
                for (int k = 0; k < nQubits; k++)
                {
                    if (((xor >> (nQubits - 1 - k)) & 1) != 0)
                        rate -= 2 * gammas[k];
                }
                _dephMask[j * _d + i] = rate; // column-major
            }

        // Pre-allocate all work matrices
        _Hrho = DenseMatrix.Create(_d, _d, Complex.Zero);
        _rhoH = DenseMatrix.Create(_d, _d, Complex.Zero);
        _drho = DenseMatrix.Create(_d, _d, Complex.Zero);
        _k1 = DenseMatrix.Create(_d, _d, Complex.Zero);
        _k2 = DenseMatrix.Create(_d, _d, Complex.Zero);
        _k3 = DenseMatrix.Create(_d, _d, Complex.Zero);
        _k4 = DenseMatrix.Create(_d, _d, Complex.Zero);
        _tmp = DenseMatrix.Create(_d, _d, Complex.Zero);
    }

    /// <summary>
    /// Evaluate drho/dt into pre-allocated output buffer. Zero allocations.
    /// </summary>
    public void EvalRHS(Matrix<Complex> rho, Matrix<Complex> output)
    {
        _H.Multiply(rho, _Hrho);
        rho.Multiply(_H, _rhoH);

        var outVals = ((DenseMatrix)output).Values;
        var hrVals = ((DenseMatrix)_Hrho).Values;
        var rhVals = ((DenseMatrix)_rhoH).Values;
        var rhoVals = ((DenseMatrix)rho).Values;
        var mask = _dephMask;
        int d2 = _d2;

        Parallel.For(0, d2, idx =>
        {
            var diff = hrVals[idx] - rhVals[idx];
            outVals[idx] = new Complex(diff.Imaginary, -diff.Real)
                         + mask[idx] * rhoVals[idx];
        });
    }

    /// <summary>
    /// rhs = EvalRHS(rho) — convenience wrapper using the public API.
    /// </summary>
    public Matrix<Complex> EvalRHS(Matrix<Complex> rho)
    {
        var result = DenseMatrix.Create(_d, _d, Complex.Zero);
        EvalRHS(rho, result);
        return result;
    }

    /// <summary>
    /// dest = a + scale * b, element-wise, zero allocations.
    /// </summary>
    private void AddScaled(Matrix<Complex> a, double scale, Matrix<Complex> b,
                           Matrix<Complex> dest)
    {
        var av = ((DenseMatrix)a).Values;
        var bv = ((DenseMatrix)b).Values;
        var dv = ((DenseMatrix)dest).Values;
        int d2 = _d2;

        Parallel.For(0, d2, i =>
        {
            dv[i] = av[i] + scale * bv[i];
        });
    }

    /// <summary>
    /// Enforce Hermiticity in-place: rho = (rho + rho†) / 2.
    /// </summary>
    private void EnforceHermiticity(Matrix<Complex> rho)
    {
        for (int i = 0; i < _d; i++)
        {
            // Diagonal: keep real part
            rho[i, i] = new Complex(rho[i, i].Real, 0);
            for (int j = i + 1; j < _d; j++)
            {
                var avg = (rho[i, j] + Complex.Conjugate(rho[j, i])) / 2.0;
                rho[i, j] = avg;
                rho[j, i] = Complex.Conjugate(avg);
            }
        }
    }

    /// <summary>
    /// RK4 integration from t=0 to tMax. Zero allocations in the hot loop.
    /// </summary>
    public void Propagate(Matrix<Complex> rho0, double tMax, double dt,
        double[] measureTimes, Action<double, Matrix<Complex>> callback)
    {
        // Copy rho0 into working matrix
        var rho = rho0.Clone();
        int measureIdx = 0;
        double t = 0;

        if (measureIdx < measureTimes.Length && measureTimes[measureIdx] <= dt / 2)
        {
            callback(0, rho);
            measureIdx++;
        }

        double dt2 = dt / 2, dt6 = dt / 6;
        int nSteps = (int)(tMax / dt) + 1;

        for (int step = 0; step < nSteps; step++)
        {
            // k1 = f(rho)
            EvalRHS(rho, _k1);

            // tmp = rho + dt/2 * k1;  k2 = f(tmp)
            AddScaled(rho, dt2, _k1, _tmp);
            EvalRHS(_tmp, _k2);

            // tmp = rho + dt/2 * k2;  k3 = f(tmp)
            AddScaled(rho, dt2, _k2, _tmp);
            EvalRHS(_tmp, _k3);

            // tmp = rho + dt * k3;  k4 = f(tmp)
            AddScaled(rho, dt, _k3, _tmp);
            EvalRHS(_tmp, _k4);

            var rv = ((DenseMatrix)rho).Values;
            var k1v = ((DenseMatrix)_k1).Values;
            var k2v = ((DenseMatrix)_k2).Values;
            var k3v = ((DenseMatrix)_k3).Values;
            var k4v = ((DenseMatrix)_k4).Values;
            var dt6local = dt6;
            Parallel.For(0, _d2, i =>
            {
                rv[i] += dt6local * (k1v[i] + 2 * k2v[i] + 2 * k3v[i] + k4v[i]);
            });

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
