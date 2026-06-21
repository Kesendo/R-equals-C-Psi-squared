using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class UniversalCarrierClaimTests
{
    // Mirrors the canonical construction recipe in ApproachFamilyCarrierClaim.Build():
    // UniversalCarrierClaim(AbsorptionTheoremClaim, Pi2DyadicLadderClaim, PolynomialDiscriminantAnchorClaim),
    // with the shared Pi2-Foundation roots built via their parameterless / typed ctors.
    private static UniversalCarrierClaim BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var polynomial = new PolynomialFoundationClaim();
        var qubit = new QubitDimensionalAnchorClaim();
        var absorption = new AbsorptionTheoremClaim(ladder);
        var discriminant = new PolynomialDiscriminantAnchorClaim(polynomial, qubit, ladder);
        return new UniversalCarrierClaim(absorption, ladder, discriminant);
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
