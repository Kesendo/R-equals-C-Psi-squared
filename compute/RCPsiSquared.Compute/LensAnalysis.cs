using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Compute;

/// <summary>
/// Lens pipeline: extract optimal single-excitation state from the slow
/// Liouvillian eigenmode. Bit convention: qubit 0 = LSB (2^0).
/// Reshape convention: vec[i*d + j] = V[i,j] (row-major, matches Python).
/// </summary>
public static class LensAnalysis
{
    public record SlowModeInfo(
        Complex Eigenvalue,
        int OriginalIndex,
        Complex[] RightEigenvector
    );

    public record LensResult(
        double[] Amplitudes,
        double MaxEigenvalue,
        double SlowModeProjection
    );

    public record LensSurveyResult(
        string Label,
        int N,
        string Topology,
        double[] Gammas,
        SlowModeInfo[] SlowModes,
        double[] SEFrobRatios,
        LensResult? LensState,
        bool SecondModeAccessible,
        double SecondModeRate
    );

    // ================================================================
    // Method 1: FindSlowModes
    // ================================================================
    public static List<SlowModeInfo> FindSlowModes(
        Complex[] eigenvalues, Complex[] vectors, int d2, int nSlowModes = 3)
    {
        var candidates = new List<(int idx, double rate)>();
        for (int j = 0; j < d2; j++)
        {
            double rate = -eigenvalues[j].Real;
            if (Math.Abs(eigenvalues[j].Real) > 1e-10)
                candidates.Add((j, rate));
        }
        candidates.Sort((a, b) => a.rate.CompareTo(b.rate));

        var result = new List<SlowModeInfo>();
        for (int i = 0; i < Math.Min(nSlowModes, candidates.Count); i++)
        {
            int j = candidates[i].idx;
            var vec = new Complex[d2];
            for (int k = 0; k < d2; k++)
                vec[k] = vectors[j * d2 + k]; // column-major: j-th eigvec
            result.Add(new SlowModeInfo(eigenvalues[j], j, vec));
        }
        return result;
    }

    // ================================================================
    // Method 2: ComputeLeftEigenvectors via full R^{-1}
    // ================================================================
    // The projection c_k = (row k of R^{-1}) . vec(rho)
    // The "left eigenvector" for the SE block is row k of R^{-1},
    // reshaped to d x d.  No conjugation needed because:
    //   c_slow = R_inv[k, :] @ rho_vec  (direct dot product)
    //
    // Returns R_inv as a dense matrix for reuse across multiple modes.
    public static Matrix<Complex> ComputeRInverse(Complex[] vectors, int d2)
    {
        var R = DenseMatrix.Create(d2, d2, (i, j) => vectors[j * d2 + i]);
        return R.Inverse();
    }

    /// <summary>
    /// Extract row k of R_inv as a flat array (the "left co-vector").
    /// This is the vector whose dot product with rho_vec gives the projection.
    /// </summary>
    public static Complex[] GetLeftCovector(Matrix<Complex> rInv, int k, int d2)
    {
        var row = new Complex[d2];
        for (int i = 0; i < d2; i++)
            row[i] = rInv[k, i];
        return row;
    }

    // ================================================================
    // Method 3: ExtractSESectorBlock
    // ================================================================
    public static (Complex[,] block, double frobRatio) ExtractSESectorBlock(
        Complex[] eigenvector, int nQubits)
    {
        int d = 1 << nQubits;

        // SE basis indices: qubit k excited = 2^(N-1-k) (qubit 0 = MSB).
        // This matches PauliOps.At() convention where qubit 0 is placed
        // leftmost in the Kronecker product (highest bit position).
        // NOTE: Build() uses MSB; BuildDirectRaw() uses LSB. They differ.
        // We use Build() here (N <= 6), so MSB convention applies.
        int[] seIdx = new int[nQubits];
        for (int k = 0; k < nQubits; k++)
            seIdx[k] = 1 << (nQubits - 1 - k);

        // Reshape: V[i,j] = eigenvector[i*d + j]
        // Extract N x N block
        var block = new Complex[nQubits, nQubits];
        for (int a = 0; a < nQubits; a++)
            for (int b = 0; b < nQubits; b++)
                block[a, b] = eigenvector[seIdx[a] * d + seIdx[b]];

        // Frobenius norms
        double blockNorm = 0, fullNorm = 0;
        for (int a = 0; a < nQubits; a++)
            for (int b = 0; b < nQubits; b++)
            {
                double m = Complex.Abs(block[a, b]);
                blockNorm += m * m;
            }
        for (int i = 0; i < d * d; i++)
        {
            double m = Complex.Abs(eigenvector[i]);
            fullNorm += m * m;
        }
        blockNorm = Math.Sqrt(blockNorm);
        fullNorm = Math.Sqrt(fullNorm);

        double ratio = fullNorm > 1e-30 ? blockNorm / fullNorm : 0.0;
        return (block, ratio);
    }

    // ================================================================
    // Method 4: ExtractLensState
    // ================================================================
    public static LensResult ExtractLensState(Complex[,] seBlock, int nQubits)
    {
        // Build Q = seBlock^T (because c_slow = a^H Q a where Q = B^T)
        var Q = DenseMatrix.Create(nQubits, nQubits,
            (i, j) => seBlock[j, i]);

        // Hermitize: M = (Q + Q^H) / 2
        var M = DenseMatrix.Create(nQubits, nQubits,
            (i, j) => (Q[i, j] + Complex.Conjugate(Q[j, i])) / 2.0);

        var evd = M.Evd();
        var evals = evd.EigenValues;
        var evecs = evd.EigenVectors;

        // Find eigenvector with largest |eigenvalue|
        int bestIdx = 0;
        double bestAbs = 0;
        for (int k = 0; k < nQubits; k++)
        {
            double a = Complex.Abs(evals[k]);
            if (a > bestAbs) { bestAbs = a; bestIdx = k; }
        }

        // Extract amplitudes (magnitudes, normalized)
        var amps = new double[nQubits];
        double norm = 0;
        for (int k = 0; k < nQubits; k++)
        {
            amps[k] = Complex.Abs(evecs[k, bestIdx]);
            norm += amps[k] * amps[k];
        }
        norm = Math.Sqrt(norm);
        for (int k = 0; k < nQubits; k++)
            amps[k] /= norm;

        // Compute |c_slow| = |a^H Q a|
        var a_vec = DenseVector.Create(nQubits,
            i => new Complex(amps[i], 0));
        var proj = a_vec.ConjugateDotProduct(Q * a_vec);

        return new LensResult(amps, evals[bestIdx].Real, Complex.Abs(proj));
    }

    // ================================================================
    // Method 5: RunFullLensPipeline
    // ================================================================
    public static LensSurveyResult RunFullLensPipeline(
        int nQubits, Bond[] bonds, double[] gammaPerQubit,
        string label, string topology, Action<string>? log = null)
    {
        log ??= _ => { };
        int d = 1 << nQubits;
        int d2 = d * d;

        // Build Liouvillian
        log($"  Building L for N={nQubits} ({label})...");
        var L = Liouvillian.Build(nQubits, bonds, gammaPerQubit);

        // Eigendecomposition
        log($"  Eigendecomposition (d2={d2})...");
        var (eigenvalues, vectors) = Liouvillian.GetAllEigenvaluesAndVectors(L);

        // Find many slow modes to reach beyond inaccessible clusters
        // (at N=5, there are 10+ modes slower than the first SE-accessible one)
        int nScan = Math.Min(d2 / 4, 50);
        var slowModes = FindSlowModes(eigenvalues, vectors, d2, nScan);
        log($"  Found {slowModes.Count} slow modes: " +
            string.Join(", ", slowModes.Select(m =>
                $"Re={m.Eigenvalue.Real:F4}")));

        // Compute R^{-1} once (expensive but correct)
        log($"  Computing R^{{-1}} ({d2}x{d2})...");
        var rInv = ComputeRInverse(vectors, d2);

        // Scan slow modes for SE content using left co-vectors (rows of R_inv)
        var processedModes = new List<SlowModeInfo>();
        var processedSE = new List<double>();
        LensResult? lensState = null;
        int firstNonSE = -1;

        for (int mi = 0; mi < slowModes.Count && (lensState == null || firstNonSE < 0); mi++)
        {
            var mode = slowModes[mi];
            var leftCovec = GetLeftCovector(rInv, mode.OriginalIndex, d2);

            // Bi-orthogonality check: leftCovec . rightVec should be ~1
            Complex biorth = Complex.Zero;
            for (int i = 0; i < d2; i++)
                biorth += leftCovec[i] * mode.RightEigenvector[i];

            var (block, frobRatio) = ExtractSESectorBlock(leftCovec, nQubits);

            if (frobRatio > 0.01 && lensState == null)
            {
                log($"  Mode {mi} (Re={mode.Eigenvalue.Real:F4}): SE={frobRatio:F4} biorth={Complex.Abs(biorth):F4} -> LENS");
                processedModes.Add(mode);
                processedSE.Add(frobRatio);
                lensState = ExtractLensState(block, nQubits);
                log($"    |c_slow|={lensState.SlowModeProjection:F4}");
                log($"    Amplitudes: [{string.Join(", ", lensState.Amplitudes.Select(a => a.ToString("F4")))}]");
            }
            else if (frobRatio <= 0.01 && firstNonSE < 0)
            {
                firstNonSE = mi;
                processedModes.Add(mode);
                processedSE.Add(frobRatio);
                log($"  Mode {mi} (Re={mode.Eigenvalue.Real:F4}): SE={frobRatio:E2} -> inaccessible");
            }
        }

        // Build final arrays: lens mode at index 0, inaccessible at index 1
        var finalModes = new List<SlowModeInfo>();
        var seFrobRatios = new List<double>();

        // Find the lens mode and inaccessible mode from processedModes
        for (int i = 0; i < processedModes.Count; i++)
            if (processedSE[i] > 0.01) { finalModes.Insert(0, processedModes[i]); seFrobRatios.Insert(0, processedSE[i]); }
            else { finalModes.Add(processedModes[i]); seFrobRatios.Add(processedSE[i]); }

        bool secondAccessible = seFrobRatios.Count > 1 && seFrobRatios[1] > 0.01;
        double secondRate = firstNonSE >= 0 ? slowModes[firstNonSE].Eigenvalue.Real : double.NaN;

        return new LensSurveyResult(
            label, nQubits, topology, gammaPerQubit,
            finalModes.ToArray(), seFrobRatios.ToArray(),
            lensState, secondAccessible, secondRate);
    }
}
