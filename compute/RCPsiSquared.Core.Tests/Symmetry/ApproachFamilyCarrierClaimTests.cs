using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class ApproachFamilyCarrierClaimTests
{
    private static ApproachFamilyCarrierClaim BuildClaim() => ApproachFamilyCarrierClaim.Build();

    [Fact]
    public void Tier_IsTier1Derived() => Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);

    [Fact]
    public void Constants_AreCanonical()
    {
        Assert.Equal(4.0, ApproachFamilyCarrierClaim.CarrierRateCoefficient);
        Assert.Equal(12.0, ApproachFamilyCarrierClaim.HarmonicRateCoefficient);
        Assert.Equal(0.75, ApproachFamilyCarrierClaim.CrossingThresholdS);
        // the 3:1 odd-harmonic ratio
        Assert.Equal(3.0, ApproachFamilyCarrierClaim.HarmonicRateCoefficient / ApproachFamilyCarrierClaim.CarrierRateCoefficient, 12);
    }

    [Fact]
    public void TwoMathematicalParents_AreExposed()
    {
        var c = BuildClaim();
        Assert.NotNull(c.Absorption);
        Assert.NotNull(c.F25);
        Assert.Null(typeof(ApproachFamilyCarrierClaim).GetProperty("Carrier"));
        Assert.Null(typeof(ApproachFamilyCarrierClaim).GetProperty("C2Ptf"));
        Assert.Null(typeof(ApproachFamilyCarrierClaim).GetProperty("TwoReadings"));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var c = BuildClaim();
        var constructor = Assert.Single(typeof(ApproachFamilyCarrierClaim).GetConstructors());
        Assert.Equal(2, constructor.GetParameters().Length);
        Assert.Throws<ArgumentNullException>(() => new ApproachFamilyCarrierClaim(null!, c.F25));
        Assert.Throws<ArgumentNullException>(() => new ApproachFamilyCarrierClaim(c.Absorption, null!));
    }

    [Fact]
    public void Summary_NamesNonAncestralComparisonsAndExactBellSpecialization()
    {
        var s = BuildClaim().Summary.ToLowerInvariant();
        Assert.Contains("prose comparison", s);
        Assert.Contains("bell+", s);
        Assert.Contains("f25", s);
        Assert.Contains("absorption theorem", s);
        Assert.Contains("n_diff=2", s);
        Assert.DoesNotContain("mode", s);
        Assert.DoesNotContain("excites", s);
    }

    [Fact]
    public void Claim_RestrictsTheLateTimeCarrierToNonzeroMembers()
    {
        var name = BuildClaim().Name.ToLowerInvariant();
        Assert.Contains("every nonzero member", name);
        Assert.Contains("0<s≤1", name);
        Assert.DoesNotContain("all members share", name);
    }

    [Fact]
    public void Claim_LabelsTheInitialReadoutAndTemporalCrossingExactly()
    {
        var claim = BuildClaim();
        var text = $"{claim.Name} {claim.Summary}";

        Assert.Contains("s is the pure-state concurrence", text, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("one third", text, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("iff γ>0 and s>3/4", text, StringComparison.Ordinal);
        Assert.Contains("γ=0", text, StringComparison.Ordinal);
        Assert.Contains("constant", text, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("entanglement itself", text, StringComparison.OrdinalIgnoreCase);
    }
}
