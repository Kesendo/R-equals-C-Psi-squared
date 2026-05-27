using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Threading.Tasks;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>F112 non-Hermitian extension verifier: enumerates the Pauli-string basis
/// at length N and checks the open identity Im⟨L_{σ_α,-i}, L_{σ_β,-i}⟩ = 0 for every
/// pair. If all pairs vanish bit-exact (within tolerance), F112 non-Hermitian extension
/// is Tier1Derived at this N by bilinearity + basis spanning.
///
/// <para>C# port of <c>simulations/_f112_open_identity_basis_enum.py</c>. Existing Python
/// anchor: N=2, 3, 4 bit-exact (35,112 distinct ordered pairs total). This primitive
/// extends to N=5 (524,800 pairs, ~16 GB working memory).</para>
///
/// <para>Per pair cost: one Frobenius inner product over a 4^N × 4^N complex matrix.
/// Per-α pre-cache cost: 4^N L matrices each of size (4^N)² complex = 16·4^(2N) bytes.
/// At N=5: 1024 matrices × 16 MB each = 16 GB total cache.</para></summary>
public static class F112NonHermitianBasisEnumeration
{
    /// <summary>Build L_H = -i[H, ·] in the Pauli basis. Mirrors the Python
    /// <c>build_L_H_pauli</c>: L_vec = -i (H ⊗ I − I ⊗ H^T) is the d² × d² vec-basis
    /// commutator superoperator; rotated into the 4^N-dim Pauli basis via
    /// T^† L_vec T / 2^N where T = <see cref="PauliBasis.VecToPauliBasisTransform"/>.
    /// Returns a 4^N × 4^N complex matrix.</summary>
    public static ComplexMatrix BuildLHInPauliBasis(ComplexMatrix H, int N)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        int d = 1 << N;
        if (H.RowCount != d || H.ColumnCount != d)
            throw new ArgumentException($"H must be {d}x{d} (2^N x 2^N); got {H.RowCount}x{H.ColumnCount}", nameof(H));

        var Id = Matrix<Complex>.Build.DenseIdentity(d);
        // L_vec = -i (H ⊗ I − I ⊗ H^T)
        var lVec = -Complex.ImaginaryOne * (H.KroneckerProduct(Id) - Id.KroneckerProduct(H.Transpose()));
        var T = PauliBasis.VecToPauliBasisTransform(N);
        return T.ConjugateTranspose() * lVec * T / (1 << N); // 2^N
    }

    /// <summary>Project a 4^N × 4^N operator-space matrix M onto the eigenspace of the
    /// Π-conjugation map X ↦ Π X Π⁻¹ at the requested eigenvalue λ ∈ {+1, -1, +i, -i}.
    /// Π is order-4 on operator space (Π⁴ = I) so the standard idempotent projector is
    /// P_λ(M) = (1/4) Σ_{k=0..3} λ^{-k} Π^k M Π^{-k}.</summary>
    public static ComplexMatrix ProjectOntoPiEigenspace(ComplexMatrix M, ComplexMatrix pi, Complex targetEigenvalue)
    {
        if (M is null) throw new ArgumentNullException(nameof(M));
        if (pi is null) throw new ArgumentNullException(nameof(pi));

        var piInv = pi.ConjugateTranspose();
        var result = Matrix<Complex>.Build.Dense(M.RowCount, M.ColumnCount);
        var curPi = Matrix<Complex>.Build.DenseIdentity(pi.RowCount);
        var curPiInv = Matrix<Complex>.Build.DenseIdentity(pi.RowCount);
        Complex lambdaToK = Complex.One;
        for (int k = 0; k < 4; k++)
        {
            var coef = (1.0 / lambdaToK) / 4.0;
            result = result + coef * (curPi * M * curPiInv);
            curPi = curPi * pi;
            curPiInv = curPiInv * piInv;
            lambdaToK = lambdaToK * targetEigenvalue;
        }
        return result;
    }

    /// <summary>Frobenius inner product ⟨A, B⟩ = Σ A[i,j]* · B[i,j] = Tr(A† B).
    ///
    /// <para>Uses the underlying dense column-major <c>Complex[]</c> storage via
    /// <see cref="DenseMatrix.Values"/> for a tight inner loop. Frobenius inner product
    /// is invariant under any consistent linearization, so the flat array can be walked
    /// in storage order without re-permuting (i, j) indices. This skips MathNet's
    /// per-element method call + bounds check; at N=5 (1024×1024 matrices × 524k pairs)
    /// the difference is hours vs minutes single-threaded.</para>
    ///
    /// <para>Throws <see cref="InvalidCastException"/> if the inputs are not dense
    /// column-major matrices. Deliberately no fallback: a silent slow-path would
    /// regress the N=5 forecast.</para></summary>
    public static Complex FrobeniusInner(ComplexMatrix A, ComplexMatrix B)
    {
        if (A is null) throw new ArgumentNullException(nameof(A));
        if (B is null) throw new ArgumentNullException(nameof(B));
        if (A.RowCount != B.RowCount || A.ColumnCount != B.ColumnCount)
            throw new ArgumentException($"shape mismatch: A is {A.RowCount}x{A.ColumnCount}, B is {B.RowCount}x{B.ColumnCount}");

        var aData = ((DenseMatrix)A).Values;
        var bData = ((DenseMatrix)B).Values;

        Complex sum = Complex.Zero;
        for (int i = 0; i < aData.Length; i++)
            sum += Complex.Conjugate(aData[i]) * bData[i];
        return sum;
    }

    /// <summary>Result record from a single F112 non-Hermitian enumeration run.
    /// Mirrors the Python script's stdout summary.</summary>
    public sealed record EnumerationResult(
        int N,
        int TotalPairs,
        int NonzeroCount,
        double MaxImaginary,
        double MeanAbsImaginary,
        TimeSpan Elapsed,
        IReadOnlyList<(string Alpha, string Beta, double Imag)> NonzeroExamples);

    /// <summary>Enumerate all upper-triangular unordered pairs (σ_α, σ_β) of Pauli
    /// strings at length N, build L_α,-i and L_β,-i for each, compute the Frobenius
    /// inner product ⟨L_α,-i, L_β,-i⟩, and report the maximum |Im|, count of pairs
    /// above tolerance, and up to 10 non-zero examples.
    ///
    /// <para>If all pairs give Im &lt; tolerance, the F112 non-Hermitian extension is
    /// proven at this N via bilinearity + Pauli-basis spanning.</para>
    ///
    /// <para>Memory: pre-caches 4^N L_α,-i matrices, each 4^N × 4^N complex = 16·4^(2N)
    /// bytes. At N=5: 16 GB total cache. At N=6: 1 TB (infeasible without sparse rep).</para>
    /// </summary>
    public static EnumerationResult Enumerate(int N, double tolerance = 1e-10)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1");
        long count = 1L << (2 * N);  // 4^N
        if (count > int.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(N), N, "N too large: 4^N overflows int32");
        int stringCount = (int)count;

        long totalPairs = (long)stringCount * (stringCount + 1) / 2;
        if (totalPairs > int.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(N), N,
                $"N too large: upper-triangular pair count {totalPairs} overflows int32 " +
                "(consider widening EnumerationResult.TotalPairs to long for N ≥ 8).");

        var stopwatch = System.Diagnostics.Stopwatch.StartNew();
        var pi = PiOperator.BuildFull(N, PauliLetter.Z);

        // Pre-compute L_α,-i for every Pauli string α. Parallelized at FULL DOP
        // (Environment.ProcessorCount). The standard codebase convention
        // (LiouvillianBlockSpectrum.cs:205) caps outer to ProcessorCount/4 to
        // leave cores for MKL's internal threading on inner BLAS calls. Welle
        // 10b benchmarks at N=5 showed outerDop=6 saturated at 6.1 cores
        // (steady 611%) for 161.8 min of cache build, indicating MKL was NOT
        // adding meaningful inner parallelism on these 1024x1024 ZGEMM calls.
        // Lifting to outerDop=ProcessorCount=24 trades the (apparently unused)
        // inner MKL threading for direct outer parallelism, with the bulk of
        // the per-iteration work being managed KroneckerProduct (not BLAS) anyway.
        var cache = new ComplexMatrix[stringCount];
        var cachePo = new ParallelOptions { MaxDegreeOfParallelism = Environment.ProcessorCount };
        Parallel.For(0, stringCount, cachePo, k =>
        {
            var letters = PauliIndex.FromFlat(k, N);
            var sigma = PauliString.Build(letters);
            var L = BuildLHInPauliBasis(sigma, N);
            cache[k] = ProjectOntoPiEigenspace(L, pi, -Complex.ImaginaryOne);
        });

        // Iterate upper-triangular pairs in parallel. Each Parallel.ForEach iteration
        // consumes a Range partition of outer indices and accumulates into thread-local
        // state; once the iteration finishes, localFinally merges per-thread totals
        // into the globals under one lock (called ~ProcessorCount times, no contention).
        //
        // Why Parallel.ForEach + thread-locals instead of BLAS ZDOTC: MathNet's Complex
        // DenseVector.DoConjugateDotProduct is a managed loop with a "TODO: provide
        // native zdotc routine" comment in the source — there is no BLAS path for
        // Complex inner products. The flat-array hot loop in FrobeniusInner (Task 3 fix,
        // commit 2ea7e81) is already the best per-pair primitive available. Parallelism
        // over the outer index is the only remaining lever.
        //
        // Why full DOP (not ProcessorCount/4): the inner work is a managed memory-bound
        // loop, NOT MKL-multithreaded BLAS. The codebase's ProcessorCount/4 throttle
        // (LiouvillianBlockSpectrum.cs:205) applies only when inner work dispatches to
        // multi-threaded MKL; oversubscription is not a concern for managed loops.
        //
        // Why chunked Partitioner: the per-outer-index workload is decreasing-triangular
        // (a=0 → stringCount pair-evals, a=stringCount-1 → 1). Default range partitioning
        // packs the heaviest iterations into one chunk and one core finishes late. A small
        // chunk size + work-stealing balances naturally.
        long globalNonzeroCount = 0;
        double globalMaxIm = 0.0;
        double globalSumAbsIm = 0.0;
        var globalExamples = new ConcurrentBag<(string Alpha, string Beta, double Imag)>();
        object lockObj = new();

        var partitioner = Partitioner.Create(0, stringCount, rangeSize: 8);

        Parallel.ForEach(partitioner,
            localInit: () => (LocalMaxIm: 0.0, LocalSumAbsIm: 0.0, LocalNonzero: 0L),
            body: (range, _, local) =>
            {
                for (int a = range.Item1; a < range.Item2; a++)
                {
                    var cacheA = cache[a];
                    for (int b = a; b < stringCount; b++)
                    {
                        var inner = FrobeniusInner(cacheA, cache[b]);
                        double absIm = Math.Abs(inner.Imaginary);
                        if (absIm > local.LocalMaxIm) local.LocalMaxIm = absIm;
                        local.LocalSumAbsIm += absIm;
                        if (absIm > tolerance)
                        {
                            local.LocalNonzero++;
                            // Race-tolerant cap: multiple threads may all see Count < 10
                            // and add concurrently; worst case ~ProcessorCount extra
                            // entries past the cap. Final .Take(10) below truncates.
                            // Acceptable: contract is "up to 10 examples", not exactly 10.
                            if (globalExamples.Count < 10)
                            {
                                string alphaName = PauliLabel.Format(PauliIndex.FromFlat(a, N));
                                string betaName = PauliLabel.Format(PauliIndex.FromFlat(b, N));
                                globalExamples.Add((alphaName, betaName, inner.Imaginary));
                            }
                        }
                    }
                }
                return local;
            },
            localFinally: local =>
            {
                lock (lockObj)
                {
                    if (local.LocalMaxIm > globalMaxIm) globalMaxIm = local.LocalMaxIm;
                    globalSumAbsIm += local.LocalSumAbsIm;
                    globalNonzeroCount += local.LocalNonzero;
                }
            });

        stopwatch.Stop();
        double meanAbsIm = globalSumAbsIm / totalPairs;
        var examples = globalExamples.Take(10).ToList();
        return new EnumerationResult(N, (int)totalPairs, (int)globalNonzeroCount,
            globalMaxIm, meanAbsIm, stopwatch.Elapsed, examples);
    }

    /// <summary>Sparse-path counterpart to <see cref="Enumerate"/>. Same contract
    /// (same EnumerationResult shape, same per-pair Frobenius inner products with
    /// Im &lt; tolerance check), but L_σ matrices are cached as SparseLSigma
    /// (~24 KB per matrix at N=5, ~50 KB at N=6, ~200 KB at N=7) instead of dense
    /// 4^N × 4^N (16 MB at N=5, 256 MB at N=6, 4 GB at N=7).
    ///
    /// <para>Unlocks N=6 as empirical anchor for the universal-N F112 closure (Welle 11):
    /// sparse cache at N=6 fits in ~200 MB instead of 1 TB. Per-pair Frobenius inner via
    /// two-pointer merge runs in O(nnz_α + nnz_β) ≈ O(2·4^N) per pair, vs O(4^(2N)) for
    /// dense. Parallel.ForEach over the partitioned outer loop follows the same pattern
    /// as Enumerate (full DOP).</para>
    ///
    /// <para>Bit-exact-matches Enumerate at N=2, 3, 4 (verified by
    /// EnumerateSparse_MatchesDenseAtSmallN parity tests). At N=5+ the sparse path is
    /// the only practical path.</para></summary>
    public static EnumerationResult EnumerateSparse(int N, double tolerance = 1e-10)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1");
        long count = 1L << (2 * N);  // 4^N
        if (count > int.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(N), N, "N too large: 4^N overflows int32");
        int stringCount = (int)count;

        long totalPairs = (long)stringCount * (stringCount + 1) / 2;
        if (totalPairs > int.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(N), N,
                $"N too large: upper-triangular pair count {totalPairs} overflows int32 " +
                "(consider widening EnumerationResult.TotalPairs to long for N ≥ 8).");

        var stopwatch = System.Diagnostics.Stopwatch.StartNew();

        // Pre-compute SparseLSigma + ProjectSparseOntoPiMinusI for every Pauli string.
        // Parallelized at full DOP; per-iteration work is managed (no MKL inner) so no
        // oversubscription concern. See Enumerate's cache-build comment block.
        var cache = new SparseLSigma[stringCount];
        var cachePo = new ParallelOptions { MaxDegreeOfParallelism = Environment.ProcessorCount };
        Parallel.For(0, stringCount, cachePo, k =>
        {
            var letters = PauliIndex.FromFlat(k, N);
            var sparseL = BuildSparseLSigma(letters, N);
            cache[k] = ProjectSparseOntoPiMinusI(sparseL);
        });

        // Iterate upper-triangular pairs in parallel via the same Partitioner + thread-local
        // accumulator pattern as Enumerate.
        long globalNonzeroCount = 0;
        double globalMaxIm = 0.0;
        double globalSumAbsIm = 0.0;
        var globalExamples = new ConcurrentBag<(string Alpha, string Beta, double Imag)>();
        object lockObj = new();

        var partitioner = Partitioner.Create(0, stringCount, rangeSize: 8);

        Parallel.ForEach(partitioner,
            localInit: () => (LocalMaxIm: 0.0, LocalSumAbsIm: 0.0, LocalNonzero: 0L),
            body: (range, _, local) =>
            {
                for (int a = range.Item1; a < range.Item2; a++)
                {
                    var cacheA = cache[a];
                    for (int b = a; b < stringCount; b++)
                    {
                        var inner = FrobeniusInnerSparse(cacheA, cache[b]);
                        double absIm = Math.Abs(inner.Imaginary);
                        if (absIm > local.LocalMaxIm) local.LocalMaxIm = absIm;
                        local.LocalSumAbsIm += absIm;
                        if (absIm > tolerance)
                        {
                            local.LocalNonzero++;
                            if (globalExamples.Count < 10)
                            {
                                string alphaName = PauliLabel.Format(PauliIndex.FromFlat(a, N));
                                string betaName = PauliLabel.Format(PauliIndex.FromFlat(b, N));
                                globalExamples.Add((alphaName, betaName, inner.Imaginary));
                            }
                        }
                    }
                }
                return local;
            },
            localFinally: local =>
            {
                lock (lockObj)
                {
                    if (local.LocalMaxIm > globalMaxIm) globalMaxIm = local.LocalMaxIm;
                    globalSumAbsIm += local.LocalSumAbsIm;
                    globalNonzeroCount += local.LocalNonzero;
                }
            });

        stopwatch.Stop();
        double meanAbsIm = globalSumAbsIm / totalPairs;
        var examples = globalExamples.Take(10).ToList();
        return new EnumerationResult(N, (int)totalPairs, (int)globalNonzeroCount,
            globalMaxIm, meanAbsIm, stopwatch.Elapsed, examples);
    }

    /// <summary>Build the sparse Pauli-basis representation of L_σ = -i[σ, ·]
    /// directly from σ's per-letter encoding, without going through dense matrices.
    /// Iterates τ ∈ Pauli basis (4^N strings) in flat-index order, computes per-letter
    /// commute/anticommute pattern, accumulates the result phase, and emits a nonzero
    /// entry only when σ and τ anticommute as full strings (odd parity over sites).
    /// Result is sorted by column index ascending.
    ///
    /// <para>Basis convention: matches <see cref="BuildLHInPauliBasis"/>. The dense
    /// path combines a vec_F (column-stacked) transform T with a vec_R-style
    /// commutator superoperator (-i·(H⊗I − I⊗Hᵀ)), which implicitly maps
    /// σ_k ↦ σ_kᵀ in the chosen Pauli basis. Single-site σ_Iᵀ = σ_I, σ_Xᵀ = σ_X,
    /// σ_Zᵀ = σ_Z but σ_Yᵀ = −σ_Y, so every nonzero element picks up the sign
    /// factor (−1)^(n_Y(τ) + n_Y(σ·τ)) where n_Y counts Y letters in the string.
    /// We apply that sign here so the sparse rep is bit-exact identical to the
    /// dense rep, enabling downstream Frobenius inner products to agree
    /// element-by-element with the existing F112 enumeration cache.</para></summary>
    public static SparseLSigma BuildSparseLSigma(IReadOnlyList<PauliLetter> sigma, int N)
    {
        if (sigma is null) throw new ArgumentNullException(nameof(sigma));
        if (sigma.Count != N) throw new ArgumentException($"sigma has {sigma.Count} letters; expected N={N}", nameof(sigma));
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be >= 1");

        long dim = 1L << (2 * N);  // 4^N
        if (dim > int.MaxValue) throw new ArgumentOutOfRangeException(nameof(N), N, "4^N overflows int32 at N >= 16");

        // Pre-allocate worst case (2 * 4^(N-1) nonzeros), trim at end.
        int maxNnz = (int)(2 * (dim / 4));
        var cols = new List<int>(maxNnz);
        var rows = new List<int>(maxNnz);
        var values = new List<Complex>(maxNnz);

        for (int col = 0; col < (int)dim; col++)
        {
            var tau = PauliIndex.FromFlat(col, N);
            var (rowLetters, phase, anticommute) = PauliProductWithPhase(sigma, tau);
            if (!anticommute) continue;
            int row = (int)PauliIndex.ToFlat(rowLetters);
            int yParity = (CountY(tau) + CountY(rowLetters)) & 1;
            double basisSign = yParity == 0 ? 1.0 : -1.0;
            cols.Add(col);
            rows.Add(row);
            // L_σ[τ] = -2i · (σ·τ), with σ·τ = phase · rowLetters; basisSign matches
            // BuildLHInPauliBasis's σ_k ↦ σ_kᵀ convention (Y picks up −1).
            values.Add(new Complex(0, -2) * phase * basisSign);
        }

        return new SparseLSigma(N, (int)dim, cols.ToArray(), rows.ToArray(), values.ToArray());
    }

    /// <summary>Project a SparseLSigma onto the Π-conjugation eigenspace at λ = −i:
    /// P_{−i}(M) = (1/4) Σ_{k=0..3} i^k · Π^k M Π^{−k}. Returns a new SparseLSigma
    /// with the same N and Dim; nnz can grow up to 4 × input nnz before sparsification
    /// collapses overlapping (row, col) entries.
    ///
    /// <para>Since Π is a signed permutation in Pauli basis (Π · σ_τ = phase(τ) ·
    /// σ_{π(τ)}), the conjugation Π · M · Π⁻¹ permutes M's (row, col) indices to
    /// (π(row), π(col)) and multiplies each entry by phase_row · conj(phase_col).
    /// The k-sum accumulates four phase-weighted copies and merges duplicates.</para></summary>
    public static SparseLSigma ProjectSparseOntoPiMinusI(SparseLSigma sparseL)
    {
        if (sparseL is null) throw new ArgumentNullException(nameof(sparseL));
        int N = sparseL.N;
        int dim = sparseL.Dim;
        int nnz = sparseL.Nnz;

        // Pre-compute the Π-action on Pauli-string indices: permutation[k][i] = π^k(i)
        // and phase[k][i] = the phase that Π^k applies when permuting index i.
        // We accumulate by composing per-letter π action over sites.
        var (perm, phases) = GetPiOrbitTables(N);  // arrays of length 4 × dim

        // Accumulator: (row, col) -> Complex sum.
        // Using a Dictionary keyed by (long)row * dim + col for fast merging.
        var accumulator = new Dictionary<long, Complex>(nnz * 4);

        var iPowers = new Complex[] { Complex.One, Complex.ImaginaryOne, -Complex.One, -Complex.ImaginaryOne };
        // i^k for k = 0, 1, 2, 3

        for (int k = 0; k < 4; k++)
        {
            var permK = perm[k];
            var phaseK = phases[k];
            var coef = iPowers[k] * 0.25;
            for (int idx = 0; idx < nnz; idx++)
            {
                int origRow = sparseL.RowIndices[idx];
                int origCol = sparseL.ColIndices[idx];
                Complex origVal = sparseL.Values[idx];
                int newRow = permK[origRow];
                int newCol = permK[origCol];
                Complex newVal = coef * phaseK[origRow] * Complex.Conjugate(phaseK[origCol]) * origVal;
                long key = (long)newRow * dim + newCol;
                if (accumulator.TryGetValue(key, out var existing))
                    accumulator[key] = existing + newVal;
                else
                    accumulator[key] = newVal;
            }
        }

        // Sparsify: drop entries below tolerance, sort by col then row for FrobeniusInnerSparse.
        const double dropTol = 1e-15;
        var keptEntries = new List<(int Col, int Row, Complex Val)>(accumulator.Count);
        foreach (var (key, val) in accumulator)
        {
            if (val.Magnitude < dropTol) continue;
            int row = (int)(key / dim);
            int col = (int)(key % dim);
            keptEntries.Add((col, row, val));
        }
        keptEntries.Sort((a, b) => a.Col != b.Col ? a.Col.CompareTo(b.Col) : a.Row.CompareTo(b.Row));

        var newCols = new int[keptEntries.Count];
        var newRows = new int[keptEntries.Count];
        var newValues = new Complex[keptEntries.Count];
        for (int i = 0; i < keptEntries.Count; i++)
        {
            newCols[i] = keptEntries[i].Col;
            newRows[i] = keptEntries[i].Row;
            newValues[i] = keptEntries[i].Val;
        }
        return new SparseLSigma(N, dim, newCols, newRows, newValues);
    }

    /// <summary>Frobenius inner product ⟨spA, spB⟩ = Σ_{(r,c)} conj(spA[r,c]) · spB[r,c]
    /// over two SparseLSigma representations. Walks both COO-style arrays in parallel
    /// via two-pointer merge on the (col, row) lexicographic order. Cost is
    /// O(nnz_A + nnz_B); both inputs must be sorted (ColIndices ascending, then
    /// RowIndices ascending within a column).</summary>
    public static Complex FrobeniusInnerSparse(SparseLSigma spA, SparseLSigma spB)
    {
        if (spA is null) throw new ArgumentNullException(nameof(spA));
        if (spB is null) throw new ArgumentNullException(nameof(spB));
        if (spA.N != spB.N) throw new ArgumentException($"N mismatch: spA.N={spA.N}, spB.N={spB.N}");

        Complex sum = Complex.Zero;
        int i = 0, j = 0;
        int nnzA = spA.Nnz, nnzB = spB.Nnz;
        while (i < nnzA && j < nnzB)
        {
            int colA = spA.ColIndices[i], colB = spB.ColIndices[j];
            if (colA < colB) { i++; continue; }
            if (colA > colB) { j++; continue; }
            // Same column: compare row.
            int rowA = spA.RowIndices[i], rowB = spB.RowIndices[j];
            if (rowA < rowB) { i++; continue; }
            if (rowA > rowB) { j++; continue; }
            // Match. Accumulate.
            sum += Complex.Conjugate(spA.Values[i]) * spB.Values[j];
            i++;
            j++;
        }
        return sum;
    }

    /// <summary>Cache of (Π^k permutation, Π^k phase) tables per N. Computes the
    /// Π-action on Pauli-string indices for k = 0, 1, 2, 3 by composing the per-letter
    /// PiOperator.ActOnLetter step-by-step.</summary>
    private static readonly ConcurrentDictionary<int, (int[][] Perm, Complex[][] Phase)> _piOrbitCache = new();

    private static (int[][] Perm, Complex[][] Phase) GetPiOrbitTables(int N) =>
        _piOrbitCache.GetOrAdd(N, BuildPiOrbitTables);

    private static (int[][] Perm, Complex[][] Phase) BuildPiOrbitTables(int N)
    {
        int dim = 1 << (2 * N);
        var perm = new int[4][];
        var phases = new Complex[4][];

        // k = 0: identity permutation, phase 1 for all.
        perm[0] = Enumerable.Range(0, dim).ToArray();
        phases[0] = Enumerable.Repeat(Complex.One, dim).ToArray();

        // For k = 1: apply Π once via ActOnLetter per site.
        perm[1] = new int[dim];
        phases[1] = new Complex[dim];
        for (int i = 0; i < dim; i++)
        {
            var letters = PauliIndex.FromFlat(i, N);
            var newLetters = new PauliLetter[N];
            Complex totalPhase = Complex.One;
            for (int l = 0; l < N; l++)
            {
                var (newL, phase) = PiOperator.ActOnLetter(letters[l], PauliLetter.Z);
                newLetters[l] = newL;
                totalPhase *= phase;
            }
            perm[1][i] = (int)PauliIndex.ToFlat(newLetters);
            phases[1][i] = totalPhase;
        }

        // For k = 2, 3: compose perm[k-1] with perm[1].
        // Π^k(i) = perm[1][perm[k-1][i]], phase = phases[k-1][i] * phases[1][perm[k-1][i]]
        for (int k = 2; k < 4; k++)
        {
            perm[k] = new int[dim];
            phases[k] = new Complex[dim];
            for (int i = 0; i < dim; i++)
            {
                int prev = perm[k - 1][i];
                perm[k][i] = perm[1][prev];
                phases[k][i] = phases[k - 1][i] * phases[1][prev];
            }
        }

        return (perm, phases);
    }

    /// <summary>Count the number of Y letters in a Pauli string. Used by
    /// <see cref="BuildSparseLSigma"/> to apply the σ_k ↦ σ_kᵀ basis-sign
    /// correction that aligns the sparse rep with the dense
    /// <see cref="BuildLHInPauliBasis"/> output.</summary>
    private static int CountY(IReadOnlyList<PauliLetter> letters)
    {
        int count = 0;
        for (int i = 0; i < letters.Count; i++)
            if (letters[i] == PauliLetter.Y) count++;
        return count;
    }

    /// <summary>Compute the per-site Pauli product σ·τ as a Pauli string plus
    /// the accumulated phase ∈ {±1, ±i} plus the parity of (σ, τ) anticommutation
    /// over sites (true = odd = strings anticommute = nonzero commutator).</summary>
    private static (PauliLetter[] ResultLetters, Complex Phase, bool Anticommute) PauliProductWithPhase(
        IReadOnlyList<PauliLetter> sigma, IReadOnlyList<PauliLetter> tau)
    {
        int N = sigma.Count;
        var resultLetters = new PauliLetter[N];
        Complex totalPhase = Complex.One;
        int anticommuteCount = 0;

        for (int l = 0; l < N; l++)
        {
            var (letter, phase, antiBit) = SingleSitePauliProduct(sigma[l], tau[l]);
            resultLetters[l] = letter;
            totalPhase *= phase;
            anticommuteCount += antiBit;
        }

        return (resultLetters, totalPhase, (anticommuteCount & 1) == 1);
    }

    /// <summary>Single-site Pauli product σ_l · τ_l. Returns the result letter,
    /// the multiplicative phase, and 1 if (σ_l, τ_l) anticommute else 0.
    /// Hand-coded 4×4 lookup; see Pauli matrix products (I·X=X, X·X=I, X·Y=iZ,
    /// Y·X=-iZ, etc.).</summary>
    private static (PauliLetter Result, Complex Phase, int AntiBit) SingleSitePauliProduct(
        PauliLetter a, PauliLetter b)
    {
        // 16-entry lookup. PauliLetter values: I=0, X=1, Z=2, Y=3.
        // Anticommute: both non-I and distinct.
        // Standard Pauli matrix products:
        //   X²=Y²=Z²=I (phase 1)
        //   XY=iZ, YX=-iZ, YZ=iX, ZY=-iX, ZX=iY, XZ=-iY
        //   I·anything = anything (phase 1)
        return (a, b) switch
        {
            (PauliLetter.I, _) => (b, Complex.One, 0),
            (_, PauliLetter.I) => (a, Complex.One, 0),
            (PauliLetter.X, PauliLetter.X) => (PauliLetter.I, Complex.One, 0),
            (PauliLetter.Y, PauliLetter.Y) => (PauliLetter.I, Complex.One, 0),
            (PauliLetter.Z, PauliLetter.Z) => (PauliLetter.I, Complex.One, 0),
            (PauliLetter.X, PauliLetter.Y) => (PauliLetter.Z, Complex.ImaginaryOne, 1),
            (PauliLetter.Y, PauliLetter.X) => (PauliLetter.Z, -Complex.ImaginaryOne, 1),
            (PauliLetter.Y, PauliLetter.Z) => (PauliLetter.X, Complex.ImaginaryOne, 1),
            (PauliLetter.Z, PauliLetter.Y) => (PauliLetter.X, -Complex.ImaginaryOne, 1),
            (PauliLetter.Z, PauliLetter.X) => (PauliLetter.Y, Complex.ImaginaryOne, 1),
            (PauliLetter.X, PauliLetter.Z) => (PauliLetter.Y, -Complex.ImaginaryOne, 1),
            _ => throw new ArgumentException($"unexpected Pauli letter combo ({a}, {b})"),
        };
    }
}
