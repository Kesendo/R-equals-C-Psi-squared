using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class ArgmaxMaxvalPairClaimTests
{
    [Fact]
    public void Claim_IsTier1Derived()
    {
        var claim = new ArgmaxMaxvalPairClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Claim_AnchorReferences_OnTheHalfCoda()
    {
        var claim = new ArgmaxMaxvalPairClaim();
        Assert.Contains("ON_THE_HALF.md", claim.Anchor);
        Assert.Contains("Coda", claim.Anchor);
    }

    [Fact]
    public void Claim_AnchorReferences_ProofBlockCpsiQuarter_HalfOfHalfReading()
    {
        var claim = new ArgmaxMaxvalPairClaim();
        Assert.Contains("PROOF_BLOCK_CPSI_QUARTER.md", claim.Anchor);
        Assert.Contains("half of half", claim.Anchor);
    }

    [Fact]
    public void Claim_HasNamedChildren_BothSidesAndPairInvariance()
    {
        var claim = new ArgmaxMaxvalPairClaim();
        IInspectable c = claim;
        var labels = c.Children.Select(ch => ch.DisplayName).ToList();
        Assert.Contains("argmax side: BilinearApexClaim", labels);
        Assert.Contains("maxval side: QuarterAsBilinearMaxvalClaim", labels);
        Assert.Contains("pair invariance", labels);
        Assert.Contains("Tom's coda (2026-05-07)", labels);
        Assert.Contains("inheritance reading", labels);
    }

    [Fact]
    public void Pair_OneHalfAndOneQuarter_AreArgmaxMaxvalOfSameParabola()
    {
        // precision: 15 marks bit-exact algebraic identity, not numerical convergence.
        const double argmax = 0.5;
        const double expectedMaxval = 0.25;
        double observedMaxval = argmax * (1 - argmax);
        Assert.Equal(expectedMaxval, observedMaxval, precision: 15);
    }

    [Fact]
    public void Pair_QuarterIsHalfSquared_AlgebraicIdentity()
    {
        // Tom's coda 2026-05-07: "1/4 ist die Hälfte von 0.5" — the quadratic shadow.
        // precision: 15 marks bit-exact algebraic identity, not numerical convergence.
        const double half = 0.5;
        Assert.Equal(0.25, half * half, precision: 15);
    }

    [Fact]
    public void Claim_PairsBothBilinearApexAndQuarterMaxval_InSummary()
    {
        var claim = new ArgmaxMaxvalPairClaim();
        Assert.Contains("BilinearApexClaim", claim.Summary);
        Assert.Contains("QuarterAsBilinearMaxvalClaim", claim.Summary);
    }
}
