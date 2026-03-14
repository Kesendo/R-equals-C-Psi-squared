using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Compute;

/// <summary>
/// Constructs the Lindblad Liouvillian superoperator.
/// L(rho) = -i[H,rho] + sum_k (L_k rho L_k^dag - 1/2 {L_k^dag L_k, rho})
/// In vectorized form: L_mat |rho>> where |rho>> = vec(rho).
/// </summary>
public static class Liouvillian
{
    /// <summary>
    /// Build Liouvillian using MathNet KroneckerProduct. Fast for N &lt;= 6.
    /// </summary>
    public static Matrix<Complex> Build(int nQubits, Bond[] bonds, double[] gammaPerQubit)
    {
        int d = 1 << nQubits;
        var H = Topology.BuildHamiltonian(nQubits, bonds);
        var Id = DenseMatrix.CreateIdentity(d);

        var L = new Complex(0, -1) * (H.KroneckerProduct(Id) - Id.KroneckerProduct(H.Transpose()));

        for (int k = 0; k < nQubits; k++)
        {
            var Lk = Math.Sqrt(gammaPerQubit[k]) * PauliOps.At(PauliOps.Z, k, nQubits);
            var LkDag = Lk.ConjugateTranspose();
            var LdL = LkDag * Lk;

            L += Lk.KroneckerProduct(Lk.Conjugate())
               - 0.5 * (LdL.KroneckerProduct(Id) + Id.KroneckerProduct(LdL.Transpose()));
        }
        return L;
    }

    /// <summary>
    /// Build Liouvillian directly into a raw array for MKL.
    /// Returns column-major Complex[] ready for z_eigen.
    /// Avoids MathNet Matrix overhead entirely.
    /// </summary>
    public static Complex[] BuildDirectRaw(int nQubits, Bond[] bonds, double[] gammaPerQubit, Action<string>? log = null)
    {
        int d = 1 << nQubits;
        int d2 = d * d;

        var H = Topology.BuildHamiltonian(nQubits, bonds);

        // Column-major: element at row r, col c = data[c * d2 + r]
        var data = new Complex[(long)d2 * d2];

        log?.Invoke($"Raw direct build: {d2}x{d2} column-major ({(long)d2 * d2 * 16 / 1e9:F2} GB)");

        var minusI = new Complex(0, -1);

        // Hamiltonian: -i(H kron I - I kron H^T)
        // Element at superop row (i*d+j), col (k*d+l):
        //   = -i * H[i,k] * delta(j,l) + i * delta(i,k) * H[l,j]

        log?.Invoke("Filling Hamiltonian...");
        for (int i = 0; i < d; i++)
        {
            for (int k = 0; k < d; k++)
            {
                var hik = H[i, k];
                if (hik == Complex.Zero) continue;
                var val = minusI * hik;
                // For all j: L[i*d+j, k*d+j] += val
                for (int j = 0; j < d; j++)
                {
                    int row = i * d + j;
                    int col = k * d + j;
                    data[(long)col * d2 + row] += val;
                }
            }
        }

        for (int j = 0; j < d; j++)
        {
            for (int l = 0; l < d; l++)
            {
                var hlj = H[l, j];
                if (hlj == Complex.Zero) continue;
                var val = minusI * hlj; // +i * H[l,j] (note the sign)
                // For all i: L[i*d+j, i*d+l] -= val
                for (int i = 0; i < d; i++)
                {
                    int row = i * d + j;
                    int col = i * d + l;
                    data[(long)col * d2 + row] -= val;
                }
            }
        }

        log?.Invoke("Filling dephasing diagonal...");
        for (int i = 0; i < d; i++)
        {
            for (int j = 0; j < d; j++)
            {
                int xor = i ^ j;
                double rate = 0;
                for (int m = 0; m < nQubits; m++)
                    if (((xor >> m) & 1) == 1)
                        rate += gammaPerQubit[m];
                if (rate > 0)
                {
                    int idx_row = i * d + j;
                    data[(long)idx_row * d2 + idx_row] -= 2.0 * rate;
                }
            }
        }

        log?.Invoke("Raw build complete.");
        return data;
    }

    public static List<double> GetOscillatoryRates(Matrix<Complex> L, double threshold = 0.05)
    {
        var evals = L.Evd().EigenValues;
        return ExtractRates(evals, threshold);
    }

    /// <summary>
    /// Eigenvalues directly from raw column-major array via MKL.
    /// No MathNet Matrix involved at all. Minimum memory path.
    /// </summary>
    public static List<double> GetOscillatoryRatesMklRaw(Complex[] columnMajorData, int n, double threshold = 0.05)
    {
        var evals = MklDirect.EigenvaluesRaw(columnMajorData, n);
        return ExtractRates(evals, threshold);
    }

    private static List<double> ExtractRates(IEnumerable<Complex> evals, double threshold)
    {
        var rates = new List<double>();
        foreach (var ev in evals)
        {
            if (Math.Abs(ev.Imaginary) > threshold)
            {
                double rate = -ev.Real;
                if (rate > 0.0001)
                    rates.Add(Math.Round(rate, 6));
            }
        }
        rates.Sort();
        return rates;
    }
}
