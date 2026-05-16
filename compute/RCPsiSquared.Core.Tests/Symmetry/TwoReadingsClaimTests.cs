using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Witnesses for <see cref="TwoReadingsClaim"/>: the typed Claim that names the
/// recurring "two readings of one underlying object" meta-pattern in the R=CΨ² framework.
/// Parent is <see cref="PolynomialFoundationClaim"/> (d=2 trunk pair-making).</summary>
public class TwoReadingsClaimTests
{
    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var claim = new TwoReadingsClaim(new PolynomialFoundationClaim());
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Build_ParentIsPolynomialFoundationClaim()
    {
        var poly = new PolynomialFoundationClaim();
        var claim = new TwoReadingsClaim(poly);
        Assert.Same(poly, claim.Polynomial);
    }

    [Fact]
    public void Build_RejectsNullParent()
    {
        Assert.Throws<ArgumentNullException>(() => new TwoReadingsClaim(null!));
    }

    [Fact]
    public void Build_SummaryNamesAllSevenLayers()
    {
        var claim = new TwoReadingsClaim(new PolynomialFoundationClaim());
        // The summary should mention at least the seven layer-names that the claim
        // collapses into the meta-pattern: number/angle, argmax/maxval, M/Π·M·Π⁻¹,
        // bra/ket, inside/outside, classical/quantum, inter-sectoral.
        Assert.Contains("number/angle", claim.Summary);
        Assert.Contains("argmax/maxval", claim.Summary);
        Assert.Contains("M/Π·M·Π⁻¹", claim.Summary);
        Assert.Contains("bra/ket", claim.Summary);
        Assert.Contains("inside/outside", claim.Summary);
        Assert.Contains("classical/quantum", claim.Summary);
        Assert.Contains("inter-sectoral", claim.Summary);
    }

    [Fact]
    public void Build_AnchorCitesCanonicalSyntheses()
    {
        var claim = new TwoReadingsClaim(new PolynomialFoundationClaim());
        Assert.Contains("ON_BOTH_SIDES_OF_THE_MIRROR", claim.Anchor);
        Assert.Contains("ON_THE_HALF", claim.Anchor);
        Assert.Contains("PRIMORDIAL_QUBIT", claim.Anchor);
        Assert.Contains("THE_OTHER_SIDE", claim.Anchor);
    }
}
