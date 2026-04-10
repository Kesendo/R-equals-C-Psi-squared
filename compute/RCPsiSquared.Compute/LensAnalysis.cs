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
    /// Compute left co-vector via targeted LU solve (no full inverse).
    /// Solves R^T y = e_k to get row k of R^{-1}.
    /// Suitable for large d2 where full inv(R) is infeasible.
    /// </summary>
    public static Complex[] ComputeLeftCovectorTargeted(
        Complex[] vectors, int d2, int targetIndex)
    {
        var R = DenseMatrix.Create(d2, d2, (i, j) => vectors[j * d2 + i]);
        var ek = DenseVector.Create(d2, 0.0);
        ek[targetIndex] = Complex.One;
        // row k of R^{-1} = solution of R^T y = e_k
        var y = R.Transpose().Solve(ek);
        var left = new Complex[d2];
        for (int i = 0; i < d2; i++)
            left[i] = y[i];
        return left;
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

        // Build Liouvillian and eigendecompose
        Complex[] eigenvalues, vectors;
        if (nQubits <= 6)
        {
            log($"  Building L for N={nQubits} ({label}) via MathNet...");
            var L = Liouvillian.Build(nQubits, bonds, gammaPerQubit);
            log($"  Eigendecomposition (d2={d2})...");
            (eigenvalues, vectors) = Liouvillian.GetAllEigenvaluesAndVectors(L);
        }
        else
        {
            log($"  Building L for N={nQubits} ({label}) via BuildDirectRaw...");
            var rawL = Liouvillian.BuildDirectRaw(nQubits,
                Topology.Chain(nQubits, Enumerable.Repeat(1.0, nQubits - 1).ToArray()),
                gammaPerQubit, msg => log($"    {msg}"));
            // Re-build with actual bonds if not chain
            if (topology != "chain")
            {
                rawL = Liouvillian.BuildDirectRaw(nQubits, bonds, gammaPerQubit);
            }
            log($"  Eigendecomposition MKL (d2={d2})...");
            (eigenvalues, vectors) = Liouvillian.GetAllEigenvaluesAndVectorsMklRaw(rawL, d2);
        }

        // Find many slow modes to reach beyond inaccessible clusters
        int nScan = Math.Min(d2 / 4, 50);
        var slowModes = FindSlowModes(eigenvalues, vectors, d2, nScan);
        log($"  Found {slowModes.Count} slow modes: " +
            string.Join(", ", slowModes.Take(5).Select(m =>
                $"Re={m.Eigenvalue.Real:F4}")) +
            (slowModes.Count > 5 ? $" ... ({slowModes.Count} total)" : ""));

        // Compute left co-vectors
        // For N<=6: full R^{-1} (fast, cached for multiple modes)
        // For N>=7: targeted LU solve per mode (avoids 16k x 16k inverse)
        bool useTargeted = d2 > 4096;
        Matrix<Complex>? rInv = null;
        if (!useTargeted)
        {
            log($"  Computing R^{{-1}} ({d2}x{d2})...");
            rInv = ComputeRInverse(vectors, d2);
        }
        else
        {
            log($"  Using targeted LU solve for left co-vectors (d2={d2}, too large for full inverse)...");
        }

        // Scan slow modes for SE content
        var processedModes = new List<SlowModeInfo>();
        var processedSE = new List<double>();
        LensResult? lensState = null;
        int firstNonSE = -1;

        if (!useTargeted)
        {
            // N<=6: scan all modes using cached R_inv (cheap, each is just a row lookup)
            log($"  Scanning {slowModes.Count} modes via R_inv rows...");
            for (int mi = 0; mi < slowModes.Count && (lensState == null || firstNonSE < 0); mi++)
            {
                var mode = slowModes[mi];
                var leftCovec = GetLeftCovector(rInv!, mode.OriginalIndex, d2);
                var (block, frobRatio) = ExtractSESectorBlock(leftCovec, nQubits);

                if (frobRatio > 0.01 && lensState == null)
                {
                    Complex biorth = Complex.Zero;
                    for (int i = 0; i < d2; i++)
                        biorth += leftCovec[i] * mode.RightEigenvector[i];
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
                }
            }
        }
        else
        {
            // N>=7: pre-screen on right eigvec, then targeted LU solve for ALL candidates
            log($"  Pre-screening {slowModes.Count} modes via right-eigvec SE block...");
            var seCandidates = new List<int>();
            for (int mi = 0; mi < slowModes.Count; mi++)
            {
                var (_, rightFrob) = ExtractSESectorBlock(slowModes[mi].RightEigenvector, nQubits);
                if (rightFrob > 0.001) seCandidates.Add(mi); // low threshold for right eigvec
                if (rightFrob <= 0.001 && firstNonSE < 0) firstNonSE = mi;
            }
            log($"  Pre-screen: {seCandidates.Count} SE candidates, first non-SE={firstNonSE}");

            foreach (int mi in seCandidates)
            {
                if (lensState != null) break;
                var mode = slowModes[mi];
                log($"  Targeted LU solve for mode {mi} (Re={mode.Eigenvalue.Real:F4})...");
                var leftCovec = ComputeLeftCovectorTargeted(vectors, d2, mode.OriginalIndex);
                var (block, frobRatio) = ExtractSESectorBlock(leftCovec, nQubits);

                if (frobRatio > 0.01)
                {
                    Complex biorth = Complex.Zero;
                    for (int i = 0; i < d2; i++)
                        biorth += leftCovec[i] * mode.RightEigenvector[i];
                    log($"    SE={frobRatio:F4} biorth={Complex.Abs(biorth):F4} -> LENS");
                    processedModes.Add(mode);
                    processedSE.Add(frobRatio);
                    lensState = ExtractLensState(block, nQubits);
                    log($"    |c_slow|={lensState.SlowModeProjection:F4}");
                    log($"    Amplitudes: [{string.Join(", ", lensState.Amplitudes.Select(a => a.ToString("F4")))}]");
                }
                else
                {
                    log($"    SE={frobRatio:E2} (right eigvec was misleading, skipping)");
                }
            }
        }

        // Build final arrays: lens mode at index 0, inaccessible at index 1
        var finalModes = new List<SlowModeInfo>(processedModes);
        var seFrobRatios = new List<double>(processedSE);

        // Add inaccessible mode info
        if (firstNonSE >= 0)
        {
            finalModes.Add(slowModes[firstNonSE]);
            seFrobRatios.Add(0.0);
        }

        bool secondAccessible = seFrobRatios.Count > 1 && seFrobRatios[1] > 0.01;
        double secondRate = firstNonSE >= 0 ? slowModes[firstNonSE].Eigenvalue.Real : double.NaN;

        return new LensSurveyResult(
            label, nQubits, topology, gammaPerQubit,
            finalModes.ToArray(), seFrobRatios.ToArray(),
            lensState, secondAccessible, secondRate);
    }
}
