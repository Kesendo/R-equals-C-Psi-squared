using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2InterChannelProjectorTests
{
    private static CoherenceBlock C2Block(int N) => new(N: N, n: 1, gammaZero: 0.05);

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Build_AcrossN5To8_LandsTier1Derived(int N)
    {
        var p = C2InterChannelProjector.Build(C2Block(N));
        Assert.Equal(Tier.Tier1Derived, p.Tier);
        Assert.True(p.IsAnalyticallyDerived);
    }

    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 2)]
    [InlineData(7, 1)]
    [InlineData(8, 2)]
    public void RankTop_MatchesAnalyticalDegeneracy(int N, int expectedRank)
    {
        var p = C2InterChannelProjector.Build(C2Block(N));
        Assert.Equal(expectedRank, p.RankTop);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Sigma0_AgreesWithInterChannelSvdReference(int N)
    {
        var block = C2Block(N);
        var p = C2InterChannelProjector.Build(block);
        var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
        Assert.Equal(svd.Sigma0, p.Sigma0, precision: 12);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void PTopL_IsHermitian(int N)
    {
        var p = C2InterChannelProjector.Build(C2Block(N));
        var diff = p.PTopL - p.PTopL.ConjugateTranspose();
        Assert.True(diff.FrobeniusNorm() < 1e-12,
            $"PTopL must be Hermitian; got Frobenius residual {diff.FrobeniusNorm():G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void PTopL_IsIdempotent(int N)
    {
        var p = C2InterChannelProjector.Build(C2Block(N));
        var diff = p.PTopL * p.PTopL - p.PTopL;
        Assert.True(diff.FrobeniusNorm() < 1e-10,
            $"PTopL must be idempotent (P² = P); got Frobenius residual {diff.FrobeniusNorm():G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void PTopL_Trace_EqualsRankTop(int N)
    {
        var p = C2InterChannelProjector.Build(C2Block(N));
        Assert.Equal(p.RankTop, p.PTopL.Trace().Real, precision: 10);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void PTopR_IsHermitianAndIdempotent_WithTraceEqualToRank(int N)
    {
        var p = C2InterChannelProjector.Build(C2Block(N));
        var hermitianRes = (p.PTopR - p.PTopR.ConjugateTranspose()).FrobeniusNorm();
        var idempotentRes = (p.PTopR * p.PTopR - p.PTopR).FrobeniusNorm();
        Assert.True(hermitianRes < 1e-12, $"PTopR Hermitian residual {hermitianRes:G3}");
        Assert.True(idempotentRes < 1e-10, $"PTopR idempotent residual {idempotentRes:G3}");
        Assert.Equal(p.RankTop, p.PTopR.Trace().Real, precision: 10);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void LibraryStabilityResidual_BelowTolerance_AcrossN(int N)
    {
        var p = C2InterChannelProjector.Build(C2Block(N));
        Assert.True(p.LibraryStabilityResidual < C2InterChannelProjector.LibraryStabilityTolerance,
            $"library stability witness failed at N={N}: residual {p.LibraryStabilityResidual:G3} ≥ {C2InterChannelProjector.LibraryStabilityTolerance:G3}");
    }

    [Fact]
    public void OddN_NonDegenerate_PTopL_RecoversInterChannelSvdSingleVectorOuterProduct()
    {
        // At N=5 the top singular value is non-degenerate (RankTop=1), so PTopL must equal
        // |u_0⟩⟨u_0| (up to a global phase that cancels in the outer product) computed by
        // the existing single-vector InterChannelSvd primitive.
        var block = C2Block(5);
        var p = C2InterChannelProjector.Build(block);
        Assert.Equal(1, p.RankTop);

        var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
        var u0 = svd.U0InFullBlock;
        var pSingleVector = u0.OuterProduct(u0.Conjugate());

        var diff = p.PTopL - pSingleVector;
        Assert.True(diff.FrobeniusNorm() < 1e-10,
            $"odd-N PTopL must agree with single-vector outer product; residual {diff.FrobeniusNorm():G3}");
    }

    [Fact]
    public void EvenN_DegenerateTop_PTopL_DiffersFromSingleVectorOuterProduct_ButIsLibraryStable()
    {
        // At N=6 the top singular value is degenerate (RankTop=2). The single-vector
        // |u_0⟩⟨u_0| from InterChannelSvd is library-dependent (it picks one direction in the
        // 2D top eigenspace via MathNet's tiebreaker), while the full projector PTopL covers
        // both directions. The two should differ; PTopL must still be library-stable.
        var block = C2Block(6);
        var p = C2InterChannelProjector.Build(block);
        Assert.Equal(2, p.RankTop);
        Assert.True(p.IsAnalyticallyDerived);

        var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
        var u0 = svd.U0InFullBlock;
        var pSingleVector = u0.OuterProduct(u0.Conjugate());

        var diff = p.PTopL - pSingleVector;
        Assert.True(diff.FrobeniusNorm() > 0.1,
            $"at degenerate even N=6 PTopL must differ from single-vector projector " +
            $"(rank 2 ≠ rank 1); got residual {diff.FrobeniusNorm():G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void SpectralGapBelowTop_IsAtLeast1000xDegeneracyTolerance(int N)
    {
        // Empirically the gap below the σ_0 cluster is small for c=2 (0.008-0.03 across
        // N=5..8); the σ_0² spectrum has a near-degenerate σ_1² sitting a few percent
        // below the strictly-equal top eigenspace. What matters for Riesz reliability is
        // that the gap is unambiguous against the DegeneracyTolerance: σ_1 is excluded
        // from the top cluster cleanly. Require ≥ 1000× DegeneracyTolerance (= 1e-5) so
        // the rank classification is robust to floating-point noise.
        var p = C2InterChannelProjector.Build(C2Block(N));
        double minimumGap = 1000.0 * C2InterChannelProjector.DegeneracyTolerance;
        Assert.True(p.SpectralGapBelowTop > minimumGap,
            $"N={N}: spectral gap below σ_0 cluster should be > {minimumGap:G3} for unambiguous rank classification; got {p.SpectralGapBelowTop:G3}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Sigma0Squared_EqualsTraceOfPTopL_TimesV_VHerm(int N)
    {
        // σ_0² · RankTop = Tr(PTopL · V V†) for the rank-RankTop top eigenspace projection.
        var block = C2Block(N);
        var p = C2InterChannelProjector.Build(block);
        var pHd1 = HdSubspaceProjector.Build(block, 1);
        var pHd3 = HdSubspaceProjector.Build(block, 3);
        var vInter = pHd1.ConjugateTranspose() * block.Decomposition.MhTotal * pHd3;
        var vvh = pHd1 * (vInter * vInter.ConjugateTranspose()) * pHd1.ConjugateTranspose();
        Complex traced = (p.PTopL * vvh).Trace();
        double expected = p.RankTop * p.Sigma0 * p.Sigma0;
        Assert.Equal(expected, traced.Real, precision: 8);
    }

    [Fact]
    public void Build_RejectsNonC2Block_WithDefaultHdPair()
    {
        // c=3 N=5 → the default (hd1=1, hd2=3) constructor should refuse.
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Equal(3, block.C);
        Assert.Throws<ArgumentException>(() => C2InterChannelProjector.Build(block));
    }

    [Fact]
    public void Anchor_References_F86Item1()
    {
        var p = C2InterChannelProjector.Build(C2Block(5));
        Assert.Contains("PROOF_F86_QPEAK", p.Anchor);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(8)]
    public void Children_IncludeRankAndStabilityWitnessNodes(int N)
    {
        IInspectable claim = C2InterChannelProjector.Build(C2Block(N));
        var labels = claim.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains("RankTop", labels);
        Assert.Contains("LibraryStabilityResidual", labels);
        Assert.Contains("SpectralGapBelowTop", labels);
        Assert.Contains("IsAnalyticallyDerived", labels);
    }
}
