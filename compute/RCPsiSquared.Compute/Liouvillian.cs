using System.Numerics;
using System.Runtime.InteropServices;
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

    /// <summary>
    /// Build Liouvillian into unmanaged (native) memory for N >= 8.
    /// Returns an IntPtr to column-major Complex data.
    /// Caller must free with NativeMemory.Free.
    /// For N=8: 65536² × 16 bytes ≈ 64 GB.
    ///
    /// Optimizations over BuildDirectRaw:
    /// - Parallel.For on Hamiltonian fill (concurrent page fault servicing)
    /// - Parallel page warm-up after build (pre-faults all pages for LAPACK)
    /// </summary>
    public static unsafe IntPtr BuildDirectNative(int nQubits, Bond[] bonds, double[] gammaPerQubit, Action<string>? log = null)
    {
        int d = 1 << nQubits;
        long d2 = d * d;
        long totalElements = d2 * d2;
        long totalBytes = totalElements * sizeof(Complex);

        log?.Invoke($"Native build: {d2}x{d2} ({totalElements:N0} elements, {totalBytes / 1e9:F1} GB)");

        var ptr = NativeMemory.AllocZeroed((nuint)totalBytes);
        var data = (Complex*)ptr;

        var H = Topology.BuildHamiltonian(nQubits, bonds);
        var minusI = new Complex(0, -1);

        // --- Hamiltonian: -i(H⊗I - I⊗H^T) ---
        // Loop 1: -i * H[i,k] * δ(j,l) → parallelize over i (each i writes distinct rows)
        log?.Invoke($"Filling Hamiltonian (native, parallel over {Environment.ProcessorCount} cores)...");

        Parallel.For(0, d, i =>
        {
            for (int k = 0; k < d; k++)
            {
                var hik = H[i, k];
                if (hik == Complex.Zero) continue;
                var val = minusI * hik;
                for (int j = 0; j < d; j++)
                {
                    long row = (long)i * d + j;
                    long col = (long)k * d + j;
                    data[col * d2 + row] += val;
                }
            }
        });

        // Loop 2: +i * δ(i,k) * H[l,j] → parallelize over j (each j writes distinct rows)
        Parallel.For(0, d, j =>
        {
            for (int l = 0; l < d; l++)
            {
                var hlj = H[l, j];
                if (hlj == Complex.Zero) continue;
                var val = minusI * hlj;
                for (int i = 0; i < d; i++)
                {
                    long row = (long)i * d + j;
                    long col = (long)i * d + l;
                    data[col * d2 + row] -= val;
                }
            }
        });

        // --- Dephasing: diagonal in |i⟩⟨j| basis ---
        log?.Invoke("Filling dephasing (native, parallel)...");
        Parallel.For(0, d, i =>
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
                    long idx = (long)i * d + j;
                    data[idx * d2 + idx] -= 2.0 * rate;
                }
            }
        });

        // --- Parallel page warm-up ---
        // The build only touches ~25% of pages (sparse Hamiltonian + diagonal).
        // Pre-fault the remaining pages so LAPACK doesn't stall on cold page faults.
        // One read per 4KB page across all cores.
        if (totalBytes > 1_000_000_000)  // only for large matrices
        {
            long pageSize = 4096;
            long totalPages = totalBytes / pageSize;
            log?.Invoke($"Warming {totalPages:N0} pages ({totalBytes / 1e9:F1} GB) across {Environment.ProcessorCount} cores...");

            var bytePtr = (byte*)ptr;
            // Volatile read forces the page to be faulted in
            Parallel.For(0L, totalPages, page =>
            {
                _ = System.Threading.Volatile.Read(ref bytePtr[page * pageSize]);
            });
            log?.Invoke("Page warm-up complete.");
        }

        log?.Invoke("Native build complete.");
        return (IntPtr)data;
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

    /// <summary>
    /// Eigenvalues ONLY from raw column-major array via direct LAPACK.
    /// No eigenvectors. Memory: input + O(n) workspace.
    /// For N=8: ~65 GB instead of ~192 GB.
    /// </summary>
    public static List<double> GetOscillatoryRatesEigenvaluesOnly(Complex[] columnMajorData, int n, Action<string>? log = null, double threshold = 0.05)
    {
        var evals = MklDirect.EigenvaluesOnlyRaw(columnMajorData, n, log);
        return ExtractRates(evals, threshold);
    }

    /// <summary>
    /// Eigenvalues ONLY from native pointer.
    /// For N >= 8 where managed arrays exceed .NET limits.
    /// Auto-selects LP64 or ILP64 based on n.
    /// </summary>
    public static List<double> GetOscillatoryRatesNative(IntPtr matrixPtr, int n, Action<string>? log = null, double threshold = 0.05)
    {
        var evals = MklDirect.EigenvaluesOnlyNative(matrixPtr, n, log);
        return ExtractRates(evals, threshold);
    }

    /// <summary>
    /// Force ILP64 path for smoke-testing at small n.
    /// </summary>
    public static List<double> GetOscillatoryRatesNativeIlp64(IntPtr matrixPtr, int n, Action<string>? log = null, double threshold = 0.05)
    {
        var evals = MklDirect.EigenvaluesOnlyNativeIlp64(matrixPtr, n, log);
        return ExtractRates(evals, threshold);
    }

    // ---- Cavity mode analysis (all eigenvalues, gamma=0) ----

    public record CavityModes(
        int Total,
        int Stationary,
        int Oscillating,
        int Anomalous,          // |Re| > eps (should be 0 for gamma=0)
        double[] Frequencies,   // sorted unique |Im| values
        double MaxAbsReal       // worst-case |Re|, for diagnostics
    );

    public static CavityModes GetCavityModes(Matrix<Complex> L, double eps = 1e-10, double freqTol = 1e-8)
    {
        var evals = L.Evd().EigenValues;
        return ClassifyCavityModes(evals, eps, freqTol);
    }

    public static CavityModes GetCavityModesMklRaw(Complex[] columnMajorData, int n, double eps = 1e-10, double freqTol = 1e-8)
    {
        var evals = MklDirect.EigenvaluesRaw(columnMajorData, n);
        return ClassifyCavityModes(evals, eps, freqTol);
    }

    private static CavityModes ClassifyCavityModes(IEnumerable<Complex> evals, double eps, double freqTol)
    {
        int total = 0, stationary = 0, oscillating = 0, anomalous = 0;
        double maxAbsReal = 0;
        var freqs = new List<double>();

        foreach (var ev in evals)
        {
            total++;
            double absRe = Math.Abs(ev.Real);
            double absIm = Math.Abs(ev.Imaginary);
            if (absRe > maxAbsReal) maxAbsReal = absRe;

            if (absRe > eps)
            {
                anomalous++;
            }
            else if (absIm < eps)
            {
                stationary++;
            }
            else
            {
                oscillating++;
                freqs.Add(absIm);
            }
        }

        // Deduplicate frequencies with tolerance
        freqs.Sort();
        var unique = new List<double>();
        foreach (var f in freqs)
        {
            if (unique.Count == 0 || Math.Abs(f - unique[^1]) > freqTol)
                unique.Add(f);
        }

        return new CavityModes(total, stationary, oscillating, anomalous, unique.ToArray(), maxAbsReal);
    }

    // ---- All eigenvalues (for RMT analysis) ----

    public static Complex[] GetAllEigenvalues(Matrix<Complex> L)
    {
        var evd = L.Evd();
        var vec = evd.EigenValues;
        var result = new Complex[vec.Count];
        for (int i = 0; i < vec.Count; i++)
            result[i] = vec[i];
        return result;
    }

    public static Complex[] GetAllEigenvaluesMklRaw(Complex[] columnMajorData, int n)
    {
        var evals = MklDirect.EigenvaluesRaw(columnMajorData, n);
        return evals is Complex[] arr ? arr : evals.ToArray();
    }

    // ---- Eigenvalues + eigenvectors (for Pauli projection) ----

    /// <summary>
    /// All eigenvalues and right eigenvectors via MathNet Evd. For N ≤ 6.
    /// vectors[j * d2 + i] = i-th component of the j-th eigenvector (column-major).
    /// </summary>
    public static (Complex[] values, Complex[] vectors) GetAllEigenvaluesAndVectors(Matrix<Complex> L)
    {
        var evd = L.Evd();
        var eigenValues = evd.EigenValues;
        var eigenVectors = evd.EigenVectors;

        int n = eigenValues.Count;
        var values = new Complex[n];
        for (int i = 0; i < n; i++)
            values[i] = eigenValues[i];

        var vectors = new Complex[(long)n * n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                vectors[(long)j * n + i] = eigenVectors[i, j];

        return (values, vectors);
    }

    /// <summary>
    /// All eigenvalues and right eigenvectors via direct MKL z_eigen. For N ≤ 7.
    /// vectors[j * n + i] = i-th component of the j-th eigenvector (column-major).
    /// </summary>
    public static (Complex[] values, Complex[] vectors) GetAllEigenvaluesAndVectorsMklRaw(Complex[] columnMajorData, int n)
    {
        return MklDirect.EigenvaluesAndVectorsRaw(columnMajorData, n);
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
