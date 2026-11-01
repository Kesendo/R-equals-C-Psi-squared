using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class UniversalCarrierClaimTests
{
    // UniversalCarrierClaim(AbsorptionTheoremClaim, Pi2DyadicLadderClaim), as in
    // ClockHandLadderClaim.Build().
    private static UniversalCarrierClaim BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var absorption = new AbsorptionTheoremClaim(ladder);
        return new UniversalCarrierClaim(absorption, ladder);
    }

    [Fact]
    public void DefaultTPeak_IsTheInverseRungTwoRate()
    {
        // t_peak = 1/(2γ₀·2): the Absorption cost of a basis pair disagreeing on two sites.
        var claim = BuildClaim();
        Assert.True(claim.DefaultTPeak == 1.0 / (4.0 * UniversalCarrierClaim.DefaultGammaZero),
            $"DefaultTPeak {claim.DefaultTPeak:R}");
        Assert.True(claim.DefaultTPeak == 5.0, $"DefaultTPeak {claim.DefaultTPeak:R}");
    }

    [Fact]
    public void Claim_BreadcrumbsToTheLiveTempoCertification()
    {
        var claim = BuildClaim();
        var kids = ((IInspectable)claim).Children.ToList();
        Assert.Contains(kids, k => k.DisplayName.Contains("live certification")
            && k.Summary.Contains("--tempo-ratio"));
    }

    [Fact]
    public void Claim_BreadcrumbsToTheSpeedOfLightQuestion()
    {
        var claim = BuildClaim();
        var kids = ((IInspectable)claim).Children.ToList();
        Assert.Contains(kids, k => k.DisplayName.Contains("speed-of-light question")
            && k.Summary.Contains("R*=1/C")
            && k.Summary.Contains("γ₀=1/(2·T₂)")
            && k.Summary.Contains("GAMMA_IS_LIGHT"));
    }
}
