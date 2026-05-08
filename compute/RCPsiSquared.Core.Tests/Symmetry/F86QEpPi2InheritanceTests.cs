using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F86QEpPi2InheritanceTests
{
    private const double PinnedGEff = 1.74;   // Endpoint g_eff from PolarityInheritanceLink

    private static F86QEpPi2Inheritance Build() =>
        new F86QEpPi2Inheritance(new Pi2DyadicLadderClaim(), new QEpLaw(PinnedGEff));

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void TwoFactor_IsExactlyTwo_FromLadderTermZero()
    {
        // F86 Q_EP's "2" numerator is exactly a_0 = 2 on the Pi2 dyadic ladder.
        var f = Build();
        Assert.Equal(2.0, f.TwoFactor, precision: 14);
    }

    [Fact]
    public void GEff_AgreesWithParentQEpLaw()
    {
        var f = Build();
        Assert.Equal(PinnedGEff, f.GEff, precision: 14);
    }

    [Fact]
    public void LiveQEp_AgreesWithParentQEpLaw()
    {
        // Composition drift check: TwoFactor / GEff bit-exact equals QEpLaw.Value.
        var f = Build();
        var parent = new QEpLaw(PinnedGEff);
        Assert.Equal(parent.Value, f.LiveQEp, precision: 14);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F86QEpPi2Inheritance(null!, new QEpLaw(PinnedGEff)));
        Assert.Throws<ArgumentNullException>(() =>
            new F86QEpPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }

    [Fact]
    public void Anchor_References_F86Proof_AndPi2Ladder()
    {
        var f = Build();
        Assert.Contains("PROOF_F86_QPEAK.md", f.Anchor);
        Assert.Contains("QEpLaw.cs", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
    }
}
