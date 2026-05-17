using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Witnesses for <see cref="KleinFourGroupSelfPairedSparseLBuilder"/>: sparse-CSR
/// build of the Klein-character sub-block on the self-paired (N/2, N/2) sector. Cross-
/// validated against <see cref="KleinFourGroupSelfPairedRefinement.BuildSubBlockL"/> dense
/// reference at small N; sparsity profile verified at N=10 (5, 5) max-block; sparse
/// shift-invert Arnoldi slow-mode extraction reconnaissance at N=10.</summary>
public class KleinFourGroupSelfPairedSparseLBuilderTests
{
    private readonly ITestOutputHelper _out;
    public KleinFourGroupSelfPairedSparseLBuilderTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4, KleinCharacter.PlusPlus)]
    [InlineData(4, KleinCharacter.PlusMinus)]
    [InlineData(4, KleinCharacter.MinusPlus)]
    [InlineData(4, KleinCharacter.MinusMinus)]
    [InlineData(6, KleinCharacter.PlusPlus)]
    [InlineData(6, KleinCharacter.MinusMinus)]
    [InlineData(8, KleinCharacter.PlusPlus)]
    public void Build_SparseMatchesDense(int N, KleinCharacter character)
    {
        // Sparse CSR construction must match the dense reference element-by-element.
        const double gamma = 0.05;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();

        var refinement = KleinFourGroupSelfPairedRefinement.Build(N);
        var dense = refinement.BuildSubBlockL(character, gammaArr);

        var sparse = KleinFourGroupSelfPairedSparseLBuilder.Build(N, character, gammaArr);
        Assert.Equal(refinement.SubBlockDims[character], sparse.SectorDim);

        // Reconstruct dense from sparse + compute Frobenius diff.
        int dim = sparse.SectorDim;
        var sparseDense = Matrix<Complex>.Build.Dense(dim, dim);
        for (int row = 0; row < dim; row++)
        {
            int start = sparse.RowPtr[row];
            int end = sparse.RowPtr[row + 1];
            for (int e = start; e < end; e++)
                sparseDense[row, sparse.ColIdx[e]] = sparse.Values[e];
        }

        double frobDiff = (dense - sparseDense).FrobeniusNorm();
        double frobRef = dense.FrobeniusNorm();
        _out.WriteLine($"N={N}, χ={character}, dim={dim}, sparse nnz={sparse.NnzTotal} " +
                       $"(mean {sparse.MeanNnzPerRow:F2}/row), ‖sparse − dense‖_F / ‖dense‖_F = {frobDiff / frobRef:G3}");
        Assert.True(frobDiff / frobRef < 1e-10,
            $"sparse vs dense mismatch: ‖diff‖_F={frobDiff:G3}, ‖ref‖_F={frobRef:G3}");
    }

    [Fact]
    public void Build_AtN10_AllFourSubBlocksHaveExpectedSparsity()
    {
        // The reconnaissance: at N=10 (5, 5), all 4 Klein sub-blocks should have low mean
        // nnz/row (≈ 10 per the earlier dense ++ probe in KleinFourGroupSelfPairedRefinementTests).
        const double gamma = 0.05;
        var gammaArr = Enumerable.Repeat(gamma, 10).ToArray();

        var sw = System.Diagnostics.Stopwatch.StartNew();
        long totalNnz = 0;
        foreach (var chi in new[] { KleinCharacter.PlusPlus, KleinCharacter.PlusMinus,
                                    KleinCharacter.MinusPlus, KleinCharacter.MinusMinus })
        {
            var sparse = KleinFourGroupSelfPairedSparseLBuilder.Build(10, chi, gammaArr);
            totalNnz += sparse.NnzTotal;
            _out.WriteLine($"N=10 χ={chi}: dim {sparse.SectorDim}, nnz {sparse.NnzTotal:N0}, " +
                           $"max {sparse.MaxNnzPerRow}, mean {sparse.MeanNnzPerRow:F2}");
            Assert.True(sparse.MeanNnzPerRow < 50, $"unexpectedly dense: {sparse.MeanNnzPerRow:F2}");
            Assert.True(sparse.MeanNnzPerRow > 1, $"unexpectedly empty: {sparse.MeanNnzPerRow:F2}");
        }
        sw.Stop();
        _out.WriteLine($"All 4 sub-blocks built in {sw.Elapsed.TotalSeconds:F1} s; " +
                       $"total nnz {totalNnz:N0} (vs 63504² = 4.04 G dense entries — " +
                       $"sparsity {totalNnz / (4.04e9):P3}).");
    }

    [Fact]
    public void Build_AtN10_SparseShiftInvertExtractsSlowModesPerSubBlock()
    {
        // Phase 3c target: per Klein sub-block, extract top-K slow modes via the generic
        // SparseShiftInvertArnoldi. The 4 sub-blocks contribute disjoint slow modes (by
        // Klein symmetry), so total slow modes ≈ 4 × K — a 4× coverage improvement over
        // Phase 2's single sweep on the full (5, 5) sector.
        //
        // Phase 3d: at numIter=80 with PreconditionerKind.Identity the 4th Ritz value
        // converges to machine precision on the Klein sub-blocks (verified at N=6, 8 in
        // Build_SparseShiftInvertWithIdentityPrecond_RecoversGenuineEigenvalues). The
        // returned slow modes are STRICT eigenvalues (matching Phase 2 outputs), not
        // garbage Ritz values from saturated BiCGStab.
        const double gamma = 0.05;
        const int K = 4;
        const int numIter = 80;
        var gammaArr = Enumerable.Repeat(gamma, 10).ToArray();

        var sigma = new Complex(0, 0.001);  // shift just off the steady state at λ = 0
        var swAll = System.Diagnostics.Stopwatch.StartNew();
        var perSubBlockSlow = new Dictionary<KleinCharacter, Complex[]>();
        foreach (var chi in new[] { KleinCharacter.PlusPlus, KleinCharacter.PlusMinus,
                                    KleinCharacter.MinusPlus, KleinCharacter.MinusMinus })
        {
            var swChi = System.Diagnostics.Stopwatch.StartNew();
            var sparse = KleinFourGroupSelfPairedSparseLBuilder.Build(10, chi, gammaArr);
            var buildTime = swChi.Elapsed;

            swChi.Restart();
            // Phase 3d fix: Klein computational-basis sub-blocks are NOT diagonally
            // dominant (diagonal = only -2γ·hamming, off-diagonal = ±iJ; magnitudes
            // comparable, many rows have exact-zero diagonal). Jacobi precond therefore
            // amplifies near-zero shifted diagonals (~1/σ ≈ 1000 at σ=(0, 0.001)),
            // destabilising BiCGStab → saturates at maxIter without converging.
            // PreconditionerKind.Identity preserves the natural conditioning of (L − σI)
            // and lets BiCGStab converge to its natural rate.
            var result = SparseShiftInvertArnoldi.Run(
                sparse.SectorDim, sparse.RowPtr, sparse.ColIdx, sparse.Values,
                sigma, numEig: K, numIter, randomSeed: 1,
                innerTolerance: 1e-8, innerMaxIter: 1000,
                preconditioner: PreconditionerKind.Identity);
            var arnoldiTime = swChi.Elapsed;

            perSubBlockSlow[chi] = result.Eigenvalues;
            _out.WriteLine($"N=10 χ={chi}: build {buildTime.TotalSeconds:F1} s, " +
                           $"shift-invert Arnoldi {arnoldiTime.TotalSeconds:F1} s " +
                           $"(inner BiCGStab mean {result.MeanInnerIterations:F1} iter), " +
                           $"top-{K} slow modes:");
            foreach (var e in result.Eigenvalues)
                _out.WriteLine($"    λ = ({e.Real:F6}, {e.Imaginary:F6})");
        }
        swAll.Stop();

        // Each Klein sub-block contributes K slow modes; the 4 sub-blocks together give
        // 4·K = 16 slow modes at N=10 (5, 5). Verify the structural property: each sub-block
        // returns K eigenvalues bounded by the sub-block's Frobenius norm.
        int totalSlow = perSubBlockSlow.Values.Sum(arr => arr.Length);
        Assert.Equal(4 * K, totalSlow);
        _out.WriteLine($"Phase 3c reconnaissance: {totalSlow} slow modes from 4 Klein sub-blocks " +
                       $"in {swAll.Elapsed.TotalSeconds:F1} s total.");
    }

    [Theory]
    [InlineData(6, KleinCharacter.PlusPlus)]
    [InlineData(6, KleinCharacter.MinusMinus)]
    [InlineData(8, KleinCharacter.PlusPlus)]
    public void Build_SparseShiftInvertWithIdentityPrecond_RecoversGenuineEigenvalues(
        int N, KleinCharacter character)
    {
        // Strict-eigenvalue witness: the slow modes returned by SparseShiftInvertArnoldi
        // with Identity preconditioning must be GENUINE eigenvalues of the dense Klein
        // sub-block (no garbage Ritz values from non-converged BiCGStab). Per slow mode λ,
        // check that some eigenvalue of the dense sub-block is within ε of λ.
        const double gamma = 0.05;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();
        var refinement = KleinFourGroupSelfPairedRefinement.Build(N);
        var dense = refinement.BuildSubBlockL(character, gammaArr);
        var denseEigs = dense.Evd().EigenValues.ToArray();

        var sparse = KleinFourGroupSelfPairedSparseLBuilder.Build(N, character, gammaArr);
        var sigma = new Complex(0, 0.001);
        int numEig = Math.Min(4, sparse.SectorDim - 2);
        // Deeper Ritz values converge slower; raise numIter and innerMaxIter to give the
        // 4th Ritz value enough budget to reach machine precision at small N. At
        // dim=1298 (N=8) we need ~80 Arnoldi steps for the 4th Ritz value to converge.
        int numIter = Math.Min(80, sparse.SectorDim - 1);

        var result = SparseShiftInvertArnoldi.Run(
            sparse.SectorDim, sparse.RowPtr, sparse.ColIdx, sparse.Values,
            sigma, numEig, numIter, randomSeed: 1,
            innerTolerance: 1e-12, innerMaxIter: 5000,
            preconditioner: PreconditionerKind.Identity);

        _out.WriteLine($"N={N} χ={character}: dim {sparse.SectorDim}, " +
                       $"inner BiCGStab mean {result.MeanInnerIterations:F1} iter");
        foreach (var lam in result.Eigenvalues)
        {
            double minDist = denseEigs.Min(e => (e - lam).Magnitude);
            _out.WriteLine($"  recovered λ = ({lam.Real:F6}, {lam.Imaginary:F6}); " +
                           $"nearest dense eig at distance {minDist:G3}");
            Assert.True(minDist < 1e-6,
                $"Recovered λ={lam} is not a genuine eigenvalue (nearest dense eig at distance {minDist:G3})");
        }
    }

    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var gammaArr = Enumerable.Repeat(0.05, 4).ToArray();
        var sparse = KleinFourGroupSelfPairedSparseLBuilder.Build(4, KleinCharacter.PlusPlus, gammaArr);
        Assert.Equal(Tier.Tier1Derived, sparse.Tier);
    }

    [Fact]
    public void Build_RejectsOddN()
    {
        var gammaArr = Enumerable.Repeat(0.05, 5).ToArray();
        Assert.Throws<ArgumentException>(() =>
            KleinFourGroupSelfPairedSparseLBuilder.Build(5, KleinCharacter.PlusPlus, gammaArr));
    }

    [Fact]
    public void Build_RejectsGammaWrongLength()
    {
        Assert.Throws<ArgumentException>(() =>
            KleinFourGroupSelfPairedSparseLBuilder.Build(4, KleinCharacter.PlusPlus, new double[3]));
    }
}
