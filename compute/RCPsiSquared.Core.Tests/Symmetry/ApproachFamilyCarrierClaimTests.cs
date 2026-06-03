using RCPsiSquared.Core.F86.Item1Derivation;
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
    public void FourTypedParents_AreExposed()
    {
        var c = BuildClaim();
        Assert.NotNull(c.Carrier);
        Assert.NotNull(c.C2Ptf);
        Assert.NotNull(c.TwoReadings);
        Assert.NotNull(c.F25);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var c = BuildClaim();
        Assert.Throws<ArgumentNullException>(() => new ApproachFamilyCarrierClaim(null!, c.C2Ptf, c.TwoReadings, c.F25));
        Assert.Throws<ArgumentNullException>(() => new ApproachFamilyCarrierClaim(c.Carrier, null!, c.TwoReadings, c.F25));
        Assert.Throws<ArgumentNullException>(() => new ApproachFamilyCarrierClaim(c.Carrier, c.C2Ptf, null!, c.F25));
        Assert.Throws<ArgumentNullException>(() => new ApproachFamilyCarrierClaim(c.Carrier, c.C2Ptf, c.TwoReadings, null!));
    }

    [Fact]
    public void Summary_NamesTheKinshipNotIdentity()
    {
        // The honesty seam: the C2 connection is a kinship/sibling, explicitly NOT a same-object identity.
        var s = BuildClaim().Summary.ToLowerInvariant();
        Assert.Contains("kinship not identity", s);
    }
}
