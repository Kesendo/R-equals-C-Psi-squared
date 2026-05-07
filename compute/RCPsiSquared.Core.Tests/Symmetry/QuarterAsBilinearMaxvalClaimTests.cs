using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
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
    public void Claim_HasNamedChildren_CalculusMandelbrotTheorem2MaximallyMixedCompanion()
    {
        var claim = new QuarterAsBilinearMaxvalClaim();
        IInspectable c = claim;
        var labels = c.Children.Select(ch => ch.DisplayName).ToList();
        Assert.Contains("calculus identity", labels);
        Assert.Contains("Mandelbrot cardioid cusp", labels);
        Assert.Contains("Theorem 2 c-block ceiling", labels);
        Assert.Contains("maximally-mixed purity", labels);
        Assert.Contains("companion to BilinearApexClaim", labels);
    }

    [Fact]
    public void BilinearForm_AtP_OneHalf_EqualsExactlyOneQuarter()
    {
        // precision: 15 marks bit-exact algebraic identity, not numerical convergence.
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
        double bilinear = p * (1 - p);
        Assert.True(bilinear < 0.25,
            $"p(1-p) at p={p} should be < 1/4; got {bilinear}");
    }

    [Fact]
    public void BilinearForm_NumericalArgmax_IsOneHalf_WithMaxvalOneQuarter()
    {
        // Operational verification: ParabolicPeakFinder gives sub-grid precision via
        // 3-point parabolic refinement, so a coarse 21-point grid suffices.
        const int samples = 21;
        var qGrid = new double[samples];
        var kCurve = new double[samples];
        for (int i = 0; i < samples; i++)
        {
            qGrid[i] = (double)i / (samples - 1);
            kCurve[i] = qGrid[i] * (1 - qGrid[i]);
        }
        var peak = ParabolicPeakFinder.Find(qGrid, kCurve);
        Assert.Equal(0.5, peak.QPeak, precision: 12);
        Assert.Equal(0.25, peak.KMax, precision: 12);
    }

    [Fact]
    public void MandelbrotDiscriminant_AtCPsi_OneQuarter_VanishesExactly()
    {
        // precision: 15 marks bit-exact algebraic identity, not numerical convergence.
        const double cpsi = 0.25;
        double discriminant = 1.0 - 4.0 * cpsi;
        Assert.Equal(0.0, discriminant, precision: 15);
    }

    [Fact]
    public void MandelbrotDiscriminant_BelowCPsi_OneQuarter_IsPositive()
    {
        const double cpsi = 0.20;
        double discriminant = 1.0 - 4.0 * cpsi;
        Assert.True(discriminant > 0,
            $"discriminant at CΨ=0.20 should be positive; got {discriminant}");
    }
}
