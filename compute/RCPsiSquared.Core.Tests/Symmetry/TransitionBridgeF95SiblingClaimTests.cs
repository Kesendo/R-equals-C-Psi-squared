using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class TransitionBridgeF95SiblingClaimTests
{
    private static TransitionBridgeF95SiblingClaim BuildClaim() =>
        TransitionBridgeF95SiblingClaim.Build();

    [Fact]
    public void Tier_IsTier1Derived()
    {
        // The algebra is Tier 1: both bridges are the F95 angle at a quadratic's discriminant zero,
        // bit-exact. ("TransitionBridge" the label and the quantum-classical reading are Tier-4, not
        // asserted by this claim.)
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void EpClockAngle_EqualsEpF95Angle_BitExact_AboveTheEp()
    {
        // The heart of the siblinghood: above the EP the F86 2-level eigenvalue is −4γ₀ ± iω, so its
        // clock Rotation arctan(ω/gap) is EXACTLY the F95 angle arctan(√(c/b²−1)) of its quadratic
        // (b=4γ₀, c=12γ₀²+J²g_eff²). Checked at two (g_eff, Q) above the EP, both γ₀-units.
        var c = BuildClaim();
        foreach (var (gEff, q) in new[] { (4.0 / 3.0, 2.5), (0.8, 3.0), (4.0 / 3.0, 2.0) })
        {
            Assert.True(q > TransitionBridgeF95SiblingClaim.QEp(gEff), "test point must be above the EP");
            double clock = c.EpClockAngle(1.0, q, gEff);
            double f95 = c.EpF95Angle(1.0, q, gEff);
            Assert.False(double.IsNaN(clock));
            Assert.False(double.IsNaN(f95));
            Assert.Equal(f95, clock, 12);
            Assert.True(c.EpClockAngleEqualsF95Angle(1.0, q, gEff));
        }
    }

    [Fact]
    public void EpAngle_IsUndefinedAtOrBelowTheEp()
    {
        // At or below the EP the eigenvalue is real (pure decay, no rotation): no angle.
        var c = BuildClaim();
        double gEff = 4.0 / 3.0; // Q_EP = 1.5
        Assert.True(double.IsNaN(c.EpClockAngle(1.0, 1.0, gEff)), "below the EP: no rotation");
        Assert.True(double.IsNaN(c.EpF95Angle(1.0, 1.0, gEff)), "below the EP: c < b², no F95 angle");
        Assert.False(c.EpClockAngleEqualsF95Angle(1.0, 1.0, gEff));
    }

    [Fact]
    public void CuspAngle_HasTheF95AnchorsAtHalf()
    {
        // The cusp side is F95 at b = ½: θ(CΨ) = arctan(√(4CΨ − 1)). 30° at the Bell+ start CΨ=1/3,
        // 45° at the anchor CΨ=½; undefined (NaN) at or below the cusp ¼ (the classical, real-root side).
        var c = BuildClaim();
        Assert.Equal(Math.PI / 6.0, c.CuspAngle(1.0 / 3.0), 12); // 30°
        Assert.Equal(Math.PI / 4.0, c.CuspAngle(0.5), 12);       // 45°
        Assert.True(double.IsNaN(c.CuspAngle(0.20)));            // below ¼: no interior angle
    }

    [Fact]
    public void QEp_IsTwoOverGEff()
    {
        Assert.Equal(1.5, TransitionBridgeF95SiblingClaim.QEp(4.0 / 3.0), 12);
        Assert.Equal(2.5, TransitionBridgeF95SiblingClaim.QEp(0.8), 12);
    }

    [Fact]
    public void Constructor_RejectsNullParent()
    {
        Assert.Throws<ArgumentNullException>(() => new TransitionBridgeF95SiblingClaim(null!));
    }

    [Fact]
    public void TypedParent_F95_IsExposed()
    {
        Assert.NotNull(BuildClaim().F95);
    }

    [Fact]
    public void Anchor_References_F95_Ep_Cusp_And_FragileBridge()
    {
        var f = BuildClaim();
        Assert.Contains("PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md", f.Anchor);
        Assert.Contains("F86_EP_THROUGH_THE_CLOCK.md", f.Anchor);
        Assert.Contains("CRITICAL_SLOWING_AT_THE_CUSP.md", f.Anchor);
        Assert.Contains("FRAGILE_BRIDGE.md", f.Anchor);
    }
}
