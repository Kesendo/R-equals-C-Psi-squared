using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86e Tier-1-derived identity: σ_0 = ‖[Π_HD1, M_H]‖ on the c=2 coherence block.
/// The F86e inter-channel SVD-top singular value equals the operator norm of the commutator
/// of the Hamming-distance-1 projector with the block Hamiltonian super-operator. Derived
/// from Π_HD1 + Π_HD3 = I on the c=2 block (HD ∈ {1, 3}) plus the lemma
/// ‖P·M·(1−P)‖ = ‖[P, M]‖. Verified bit-exact (Python residual ~1e-15) for N=4..12; this
/// class reproduces it through N=8. See <see cref="SigmaZeroCommutatorNormClaim"/>.</summary>
public class SigmaZeroCommutatorNormClaimTests
{
    private readonly ITestOutputHelper _out;
    public SigmaZeroCommutatorNormClaimTests(ITestOutputHelper output) => _out = output;

    /// <summary>The core identity: across N=5..8 the residual |σ_0 − ‖[Π_HD1, M_H]‖|
    /// vanishes to machine precision. The c=2 block is built with n=1 (popcount-1 ⊗
    /// popcount-2), the unique c=2 stratum.</summary>
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Sigma0_EqualsCommutatorNorm_OnC2Block(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        Assert.Equal(2, block.C); // sanity: n=1 is the c=2 stratum

        var claim = SigmaZeroCommutatorNormClaim.Build(block);

        _out.WriteLine($"N={N}: σ_0 = {claim.Sigma0:F12}, " +
                       $"‖[Π_HD1, M_H]‖ = {claim.CommutatorNorm:F12}, " +
                       $"residual = {claim.Residual:E3}, " +
                       $"projector-idempotency residual = {claim.ProjectorIdempotencyResidual:E3}");

        Assert.True(claim.Residual < 1e-10,
            $"σ_0 = ‖[Π_HD1, M_H]‖ identity broken at N={N}: " +
            $"σ_0={claim.Sigma0:F12}, ‖[Π_HD1, M_H]‖={claim.CommutatorNorm:F12}, residual={claim.Residual:E3}");

        // The square projector built as P·P† must be genuinely idempotent.
        Assert.True(claim.ProjectorIdempotencyResidual < 1e-10,
            $"Square HD=1 projector not idempotent at N={N}: ‖Π·Π − Π‖={claim.ProjectorIdempotencyResidual:E3}");
    }

    /// <summary>At N=7 the c=2 σ_0 hits the sweet-spot crossing 2√2 bit-exactly (the
    /// <see cref="SigmaZeroChromaticityScaling.SigmaZeroSweetSpotIdentity_C2N7"/> value).</summary>
    [Fact]
    public void Sigma0_AtN7_EqualsTwoSqrtTwo()
    {
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var claim = SigmaZeroCommutatorNormClaim.Build(block);

        double twoSqrt2 = 2.0 * Math.Sqrt(2.0);
        _out.WriteLine($"N=7: σ_0 = {claim.Sigma0:F12}, 2√2 = {twoSqrt2:F12}, " +
                       $"|σ_0 − 2√2| = {Math.Abs(claim.Sigma0 - twoSqrt2):E3}");

        Assert.Equal(twoSqrt2, claim.Sigma0, precision: 9);
    }

    /// <summary>The claim is registered as Tier1Derived (analytic proof, bit-exact verified).</summary>
    [Fact]
    public void Claim_IsTier1Derived()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var claim = SigmaZeroCommutatorNormClaim.Build(block);
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    /// <summary>The identity is c=2-specific: building on a c=1 block (n=0, HD ∈ {1} only)
    /// must throw, because Π_HD1 + Π_HD3 = I fails outside c=2.</summary>
    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var c1Block = new CoherenceBlock(N: 3, n: 0, gammaZero: 0.05);
        Assert.Equal(1, c1Block.C);
        Assert.Throws<ArgumentException>(() => SigmaZeroCommutatorNormClaim.Build(c1Block));

        var c3Block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Equal(3, c3Block.C);
        Assert.Throws<ArgumentException>(() => SigmaZeroCommutatorNormClaim.Build(c3Block));
    }

    /// <summary>The F86 knowledge base surfaces the claim under Tier 1 (derived) for a c=2
    /// block and omits it (null) for a c=1 block.</summary>
    [Fact]
    public void F86KnowledgeBase_ExposesClaim_ForC2_OmitsForC1()
    {
        var c2Kb = new F86KnowledgeBase(new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05));
        Assert.NotNull(c2Kb.SigmaZeroCommutatorNorm);
        Assert.Equal(Tier.Tier1Derived, c2Kb.SigmaZeroCommutatorNorm!.Tier);
        Assert.True(c2Kb.SigmaZeroCommutatorNorm.Residual < 1e-10);

        var c1Kb = new F86KnowledgeBase(new CoherenceBlock(N: 3, n: 0, gammaZero: 0.05));
        Assert.Null(c1Kb.SigmaZeroCommutatorNorm);
    }
}
