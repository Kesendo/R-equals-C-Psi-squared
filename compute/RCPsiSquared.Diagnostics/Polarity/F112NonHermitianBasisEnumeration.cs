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

        // Pre-compute L_α,-i for every Pauli string α. Parallelized with
        // outerDop = ProcessorCount/4: the inner work (KroneckerProduct +
        // ConjugateTranspose + Multiplication on 4^N × 4^N matrices + 4 matmuls in
        // ProjectOntoPiEigenspace) dispatches to MKL-multithreaded BLAS for the
        // matrix multiplications. Capping outer parallelism leaves cores free for
        // MKL's internal threading on the larger N=5 matmuls (1024^3 ZGEMM each).
        // Same pattern as LiouvillianBlockSpectrum.cs:205. Empirical observation
        // (Welle 10b SLOW N=5 run): cache build at outerDop=6 took 161.8 min;
        // outerDop=ProcessorCount=24 trials in Welle 10c may show better scaling
        // since MKL ZGEMM on 1024x1024 doesn't saturate 24 cores anyway.
        var cache = new ComplexMatrix[stringCount];
        int cacheOuterDop = Math.Max(1, Environment.ProcessorCount / 4);
        var cachePo = new ParallelOptions { MaxDegreeOfParallelism = cacheOuterDop };
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
}
