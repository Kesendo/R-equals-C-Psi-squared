using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Compute;

/// <summary>
/// Lens pipeline: extract optimal single-excitation state from the slow
/// Liouvillian eigenmode.
///
/// Bit convention: qubit 0 = MSB (all paths).
///   Build() uses PauliOps.At() Kronecker order (MSB).
///   BuildDirectRaw()/BuildDirectNative() use MSB-consistent dephasing.
///   SE index for qubit k: 2^(N-1-k).
///
/// Reshape convention (all paths): vec[i*d + j] = V[i,j] (row-major).
/// </summary>
public static class LensAnalysis
{
    // SlowModeInfo stores only the column index into the vectors array,
    // not a copy of the 268M-element eigenvector. Read elements on demand.
    public record SlowModeInfo(Complex Eigenvalue, int OriginalIndex);
    public record LensResult(double[] Amplitudes, double MaxEigenvalue, double SlowModeProjection);
    public record LensSurveyResult(
        string Label, int N, string Topology, double[] Gammas,
        SlowModeInfo[] SlowModes, double[] SEFrobRatios,
        LensResult? LensState, bool SecondModeAccessible, double SecondModeRate);

    public static int[] SEIndices(int nQubits, bool msb)
    {
        var idx = new int[nQubits];
        for (int k = 0; k < nQubits; k++)
            idx[k] = msb ? 1 << (nQubits - 1 - k) : 1 << k;
        return idx;
    }

    // ================================================================
    // FindSlowModes - stores only indices, no vector copies
    // ================================================================
    public static List<SlowModeInfo> FindSlowModes(
        Complex[] eigenvalues, int d2, int nSlowModes = 3)
    {
        var candidates = new List<(int idx, double rate)>();
        for (int j = 0; j < d2; j++)
        {
            if (Math.Abs(eigenvalues[j].Real) > 1e-10)
                candidates.Add((j, -eigenvalues[j].Real));
        }
        candidates.Sort((a, b) => a.rate.CompareTo(b.rate));

        var result = new List<SlowModeInfo>();
        for (int i = 0; i < Math.Min(nSlowModes, candidates.Count); i++)
            result.Add(new SlowModeInfo(eigenvalues[candidates[i].idx], candidates[i].idx));
        return result;
    }

    // ================================================================
    // Read N elements from a column of the vectors array (no full copy)
    // ================================================================
    static Complex[] ReadSEBlock(Complex[] vectors, int d2, int colIndex, int[] seIdx, int d)
    {
        // Extract only the SE-indexed elements from column colIndex
        // vectors[colIndex * d2 + row] = element at row for the colIndex-th eigenvector
        int n = seIdx.Length;
        var block = new Complex[n * n];
        long colOff = (long)colIndex * d2;
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
                block[a * n + b] = vectors[colOff + seIdx[a] * d + seIdx[b]];
        return block;
    }

    // ================================================================
    // Left eigenvectors
    // ================================================================

    /// <summary>Full R^{-1} for N &lt;= 6.</summary>
    public static Matrix<Complex> ComputeRInverse(Complex[] vectors, int d2)
    {
        var storage = new Complex[vectors.Length];
        Array.Copy(vectors, storage, vectors.Length);
        return new DenseMatrix(d2, d2, storage).Inverse();
    }

    public static Complex[] GetLeftCovector(Matrix<Complex> rInv, int k, int d2)
    {
        var row = new Complex[d2];
        for (int i = 0; i < d2; i++)
            row[i] = rInv[k, i];
        return row;
    }

    /// <summary>
    /// Targeted solve for N >= 7 via direct LAPACK zgesv (bypasses MathNet 2GB limit).
    /// Solves R^T y = e_k to get row k of R^{-1}.
    /// </summary>
    public static Complex[] ComputeLeftCovectorTargeted(
        Complex[] vectors, int d2, int targetIndex, Action<string>? log = null)
    {
        return MklDirect.SolveRTranspose(vectors, d2, targetIndex, log);
    }

    // ================================================================
    // ExtractSESectorBlock from a flat eigenvector, with precomputed fullNormSq
    // ================================================================
    public static (Complex[,] block, double frobRatio) ExtractSESectorBlock(
        Complex[] eigenvector, int nQubits, int[] seIdx, double fullNormSq)
    {
        int d = 1 << nQubits;
        var block = new Complex[nQubits, nQubits];
        for (int a = 0; a < nQubits; a++)
            for (int b = 0; b < nQubits; b++)
                block[a, b] = eigenvector[seIdx[a] * d + seIdx[b]];

        double blockNormSq = 0;
        for (int a = 0; a < nQubits; a++)
            for (int b = 0; b < nQubits; b++)
                blockNormSq += Complex.Abs(block[a, b]) * Complex.Abs(block[a, b]);

        double ratio = fullNormSq > 1e-60
            ? Math.Sqrt(blockNormSq / fullNormSq) : 0.0;
        return (block, ratio);
    }

    // Compute Frobenius norm squared of a flat vector
    static double FrobNormSq(Complex[] v)
    {
        double s = 0;
        for (int i = 0; i < v.Length; i++)
            s += v[i].Real * v[i].Real + v[i].Imaginary * v[i].Imaginary;
        return s;
    }

    // SE block Frobenius ratio from the vectors array directly (no full copy)
    static double SEFrobRatioFromVectors(Complex[] vectors, int d2, int colIndex,
        int[] seIdx, int d, int nQubits)
    {
        long colOff = (long)colIndex * d2;
        double blockSq = 0;
        for (int a = 0; a < nQubits; a++)
            for (int b = 0; b < nQubits; b++)
            {
                var c = vectors[colOff + seIdx[a] * d + seIdx[b]];
                blockSq += c.Real * c.Real + c.Imaginary * c.Imaginary;
            }
        double fullSq = 0;
        for (int i = 0; i < d2; i++)
        {
            var c = vectors[colOff + i];
            fullSq += c.Real * c.Real + c.Imaginary * c.Imaginary;
        }
        return fullSq > 1e-60 ? Math.Sqrt(blockSq / fullSq) : 0.0;
    }

    // ================================================================
    // ExtractLensState
    // ================================================================
    public static LensResult ExtractLensState(Complex[,] seBlock, int nQubits)
    {
        var Q = DenseMatrix.Create(nQubits, nQubits, (i, j) => seBlock[j, i]);
        var M = DenseMatrix.Create(nQubits, nQubits,
            (i, j) => (Q[i, j] + Complex.Conjugate(Q[j, i])) / 2.0);

        var evd = M.Evd();
        var evals = evd.EigenValues;
        var evecs = evd.EigenVectors;

        int bestIdx = 0;
        double bestAbs = 0;
        for (int k = 0; k < nQubits; k++)
        {
            double a = Complex.Abs(evals[k]);
            if (a > bestAbs) { bestAbs = a; bestIdx = k; }
        }

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

        var a_vec = DenseVector.Create(nQubits, i => new Complex(amps[i], 0));
        var proj = a_vec.ConjugateDotProduct(Q * a_vec);
        return new LensResult(amps, evals[bestIdx].Real, Complex.Abs(proj));
    }

    // ================================================================
    // RunFullLensPipeline
    // ================================================================
    public static LensSurveyResult RunFullLensPipeline(
        int nQubits, Bond[] bonds, double[] gammaPerQubit,
        string label, string topology, Action<string>? log = null)
    {
        log ??= _ => { };
        int d = 1 << nQubits;
        int d2 = d * d;
        // All build paths use MSB convention after the dephasing fix
        var seIdx = SEIndices(nQubits, msb: true);

        // Build and eigendecompose
        Complex[] eigenvalues, vectors;
        if (nQubits <= 6)
        {
            log($"  Building L (N={nQubits}, MSB, MathNet)...");
            var L = Liouvillian.Build(nQubits, bonds, gammaPerQubit);
            log($"  Eigendecomposition (d2={d2})...");
            (eigenvalues, vectors) = Liouvillian.GetAllEigenvaluesAndVectors(L);
        }
        else
        {
            int cores = Environment.ProcessorCount;
            log($"  Building L (N={nQubits}, LSB, {cores} cores)...");
            var rawL = Liouvillian.BuildDirectRaw(nQubits, bonds, gammaPerQubit,
                msg => log($"    {msg}"));
            MklDirect.ConfigureThreads(cores, msg => log($"    {msg}"));
            log($"  Eigendecomposition LAPACK zgeev (d2={d2})...");
            (eigenvalues, vectors) = MklDirect.EigenvaluesAndVectorsDirectRaw(
                rawL, d2, msg => log($"    {msg}"));
        }

        // Find slow modes (indices only, no vector copies)
        int nScan = Math.Min(d2 / 4, 50);
        var slowModes = FindSlowModes(eigenvalues, d2, nScan);
        log($"  {slowModes.Count} slow modes found");

        var processedModes = new List<SlowModeInfo>();
        var processedSE = new List<double>();
        LensResult? lensState = null;
        int firstNonSE = -1;

        if (nQubits <= 6)
        {
            // N<=6: full R_inv, scan all modes
            log($"  Computing R^{{-1}} ({d2}x{d2})...");
            var rInv = ComputeRInverse(vectors, d2);

            for (int mi = 0; mi < slowModes.Count && (lensState == null || firstNonSE < 0); mi++)
            {
                var mode = slowModes[mi];
                var leftCovec = GetLeftCovector(rInv, mode.OriginalIndex, d2);
                double leftNormSq = FrobNormSq(leftCovec);
                var (block, frobRatio) = ExtractSESectorBlock(leftCovec, nQubits, seIdx, leftNormSq);

                if (frobRatio > 0.01 && lensState == null)
                {
                    log($"  Mode {mi} (Re={mode.Eigenvalue.Real:F4}): SE={frobRatio:F4} -> LENS");
                    processedModes.Add(mode);
                    processedSE.Add(frobRatio);
                    lensState = ExtractLensState(block, nQubits);
                    log($"    Amplitudes: [{string.Join(", ", lensState.Amplitudes.Select(a => a.ToString("F4")))}]");
                }
                else if (frobRatio <= 0.01 && firstNonSE < 0)
                    firstNonSE = mi;
            }
        }
        else
        {
            // N>=7: pre-screen on right eigvec (reads 49 elements, no copy),
            // then targeted LU solve only for candidates
            log($"  Pre-screening {slowModes.Count} modes...");
            var seCandidates = new List<int>();
            for (int mi = 0; mi < slowModes.Count; mi++)
            {
                double rf = SEFrobRatioFromVectors(vectors, d2, slowModes[mi].OriginalIndex,
                    seIdx, d, nQubits);
                if (rf > 0.001) seCandidates.Add(mi);
                if (rf <= 0.001 && firstNonSE < 0) firstNonSE = mi;
            }
            log($"  {seCandidates.Count} SE candidates, first non-SE={firstNonSE}");

            // Build DenseMatrix once for all targeted solves (single 4GB copy)
            DenseMatrix? R = null;
            foreach (int mi in seCandidates)
            {
                if (lensState != null) break;
                var mode = slowModes[mi];
                log($"  Targeted solve for mode {mi} (Re={mode.Eigenvalue.Real:F4})...");
                var leftCovec = ComputeLeftCovectorTargeted(
                    vectors, d2, mode.OriginalIndex, msg => log($"    {msg}"));
                double leftNormSq = FrobNormSq(leftCovec);
                var (block, frobRatio) = ExtractSESectorBlock(leftCovec, nQubits, seIdx, leftNormSq);

                if (frobRatio > 0.01)
                {
                    log($"    SE={frobRatio:F4} -> LENS");
                    processedModes.Add(mode);
                    processedSE.Add(frobRatio);
                    lensState = ExtractLensState(block, nQubits);
                    log($"    Amplitudes: [{string.Join(", ", lensState.Amplitudes.Select(a => a.ToString("F4")))}]");
                }
                else
                    log($"    SE={frobRatio:E2} (skipping)");
            }
        }

        // Result: lens mode first, inaccessible second
        var finalModes = new List<SlowModeInfo>(processedModes);
        var seFrobRatios = new List<double>(processedSE);
        if (firstNonSE >= 0)
        {
            finalModes.Add(slowModes[firstNonSE]);
            seFrobRatios.Add(0.0);
        }

        return new LensSurveyResult(
            label, nQubits, topology, gammaPerQubit,
            finalModes.ToArray(), seFrobRatios.ToArray(),
            lensState,
            seFrobRatios.Count > 1 && seFrobRatios[1] > 0.01,
            firstNonSE >= 0 ? slowModes[firstNonSE].Eigenvalue.Real : double.NaN);
    }
}
