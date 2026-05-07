using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class QuarterAsBilinearMaxvalClaimTests
{
    [Fact]
    public void Claim_IsTier1Derived()
    {
        var claim = new QuarterAsBilinearMaxvalClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Claim_AnchorReferences_ProofBlockCpsiQuarter()
    {
        var claim = new QuarterAsBilinearMaxvalClaim();
        Assert.Contains("PROOF_BLOCK_CPSI_QUARTER.md", claim.Anchor);
    }

    [Fact]
    public void Claim_AnchorReferences_RoadmapQuarterBoundary()
    {
        var claim = new QuarterAsBilinearMaxvalClaim();
        Assert.Contains("PROOF_ROADMAP_QUARTER_BOUNDARY.md", claim.Anchor);
    }

    [Fact]
    public void Claim_HasNamedChildren_CalculusMandelbrotTheorem2()
    {
        var claim = new QuarterAsBilinearMaxvalClaim();
        IInspectable c = claim;
        var labels = c.Children.Select(ch => ch.DisplayName).ToList();
        Assert.Contains("calculus identity", labels);
        Assert.Contains("Mandelbrot cardioid cusp", labels);
        Assert.Contains("Theorem 2 c-block ceiling", labels);
        Assert.Contains("d = 2 Pauli normalization", labels);
        Assert.Contains("companion to BilinearApexClaim", labels);
    }

    [Fact]
    public void BilinearForm_AtP_OneHalf_EqualsExactlyOneQuarter()
    {
        // Algebraic identity: max p(1-p) on [0,1] is 1/4 at p = 1/2.
        const double p = 0.5;
        double bilinear = p * (1 - p);
        Assert.Equal(0.25, bilinear, precision: 15);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(0.1)]
    [InlineData(0.25)]
    [InlineData(0.4)]
    [InlineData(0.6)]
    [InlineData(0.75)]
    [InlineData(0.9)]
    [InlineData(1.0)]
    public void BilinearForm_OffPeak_IsStrictlyLessThanOneQuarter(double p)
    {
        // Off-peak: p ∈ [0,1] \ {1/2} satisfies p(1-p) < 1/4 strictly.
        double bilinear = p * (1 - p);
        Assert.True(bilinear < 0.25,
            $"p(1-p) at p={p} should be < 1/4; got {bilinear}");
    }

    [Fact]
    public void BilinearForm_NumericalArgmax_IsOneHalf_WithMaxvalOneQuarter()
    {
        // Operational verification of the argmax/maxval pair.
        const int samples = 10001;
        double bestP = double.NaN;
        double bestValue = double.NegativeInfinity;
        for (int i = 0; i < samples; i++)
        {
            double p = (double)i / (samples - 1);
            double v = p * (1 - p);
            if (v > bestValue)
            {
                bestValue = v;
                bestP = p;
            }
        }
        Assert.Equal(0.5, bestP, precision: 4);
        Assert.Equal(0.25, bestValue, precision: 8);
    }

    [Fact]
    public void MandelbrotDiscriminant_AtCPsi_OneQuarter_VanishesExactly()
    {
        // R = CΨ² discriminant 1 − 4·CΨ vanishes at CΨ = 1/4 = (1/2)².
        const double cpsi = 0.25;
        double discriminant = 1.0 - 4.0 * cpsi;
        Assert.Equal(0.0, discriminant, precision: 15);
    }

    [Fact]
    public void MandelbrotDiscriminant_BelowCPsi_OneQuarter_IsPositive()
    {
        // CΨ < 1/4: discriminant is positive (two real fixed points pre-cusp).
        const double cpsi = 0.20;
        double discriminant = 1.0 - 4.0 * cpsi;
        Assert.True(discriminant > 0,
            $"discriminant at CΨ=0.20 should be positive; got {discriminant}");
    }
}

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
        // Operational closure: the BilinearApexClaim argmax (1/2) and the
        // QuarterAsBilinearMaxvalClaim maxval (1/4) come from the same parabola p(1-p).
        const double argmax = 0.5;
        const double expectedMaxval = 0.25;
        double observedMaxval = argmax * (1 - argmax);
        Assert.Equal(expectedMaxval, observedMaxval, precision: 15);
    }

    [Fact]
    public void Pair_QuarterIsHalfSquared_AlgebraicIdentity()
    {
        // 1/4 = (1/2)² — Tom's coda 2026-05-07 "1/4 ist die Hälfte von 0.5"
        // (the quarter is the half's quadratic shadow).
        const double half = 0.5;
        Assert.Equal(0.25, half * half, precision: 15);
    }

    [Fact]
    public void Claim_PairsBothBilinearApexAndQuarterMaxval_InSummary()
    {
        // The synthesis claim explicitly names both paired primitives in its summary.
        var claim = new ArgmaxMaxvalPairClaim();
        Assert.Contains("BilinearApexClaim", claim.Summary);
        Assert.Contains("QuarterAsBilinearMaxvalClaim", claim.Summary);
    }
}
