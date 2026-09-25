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
        // The algebra is Tier 1. "TransitionBridge" and any physical interpretation are not.
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void EpClockAngle_EqualsEpF95AngleWithinTolerance_AboveTheEp()
    {
        // The heart of the siblinghood: above the EP the F86 2-level eigenvalue is −4γ₀ ± iω, so its
        // clock Rotation and F95 are the same algebraic expression. Floating evaluations are
        // compared with a tolerance.
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
    public void EpAngle_IsZeroAtTheEp_AndUndefinedBelowIt()
    {
        // At the discriminant zero the angle is exactly zero; below it the F95
        // complex-root angle is undefined.
        var c = BuildClaim();
        double gEff = 4.0 / 3.0; // Q_EP = 1.5
        Assert.True(double.IsNaN(c.EpClockAngle(1.0, 1.0, gEff)), "below the EP: no rotation");
        Assert.True(double.IsNaN(c.EpF95Angle(1.0, 1.0, gEff)), "below the EP: c < b², no F95 angle");
        Assert.False(c.EpClockAngleEqualsF95Angle(1.0, 1.0, gEff));
        Assert.Equal(0.0, c.EpClockAngle(1.0, 1.5, gEff));
        Assert.Equal(0.0, c.EpF95Angle(1.0, 1.5, gEff));
        Assert.True(c.EpClockAngleEqualsF95Angle(1.0, 1.5, gEff));
    }

    [Theory]
    [InlineData(4.0 / 3.0)]
    [InlineData(0.8)]
    public void EpAngle_AtQEpApiBoundary_IsExactlyZero(double gEff)
    {
        var c = BuildClaim();
        double qEp = TransitionBridgeF95SiblingClaim.QEp(gEff);
        Assert.Equal(0.0, c.EpClockAngle(1.0, qEp, gEff));
        Assert.Equal(0.0, c.EpF95Angle(1.0, qEp, gEff));
        Assert.True(c.EpClockAngleEqualsF95Angle(1.0, qEp, gEff));
    }

    [Fact]
    public void CuspAngle_HasTheF95AnchorsAtHalf()
    {
        // The cusp side is F95 at b = ½: θ(CΨ) = arctan(√(4CΨ − 1)). 30° at the Bell+ start CΨ=1/3,
        // 45° at the anchor CΨ=½; zero at the cusp and undefined below it.
        var c = BuildClaim();
        Assert.Equal(Math.PI / 6.0, c.CuspAngle(1.0 / 3.0), 12); // 30°
        Assert.Equal(Math.PI / 4.0, c.CuspAngle(0.5), 12);       // 45°
        Assert.Equal(0.0, c.CuspAngle(0.25));                    // discriminant zero
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
    public void Surface_SeparatesRecurrenceVariableFromLiouvillianDecayVariable()
    {
        var claim = BuildClaim();
        string surface = string.Join("\n", new[] { claim.Name, claim.DisplayName, claim.Summary }
            .Concat(claim.Children.Select(child => $"{child.DisplayName}\n{child.Summary}")));

        Assert.Contains("z_rec", surface, StringComparison.Ordinal);
        Assert.Contains("z_decay=−λ", surface, StringComparison.Ordinal);
        Assert.DoesNotContain("λ=−γ₀ alone", surface, StringComparison.Ordinal);
        Assert.DoesNotContain("bit-exact", surface, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Surface_NamesTwoDistinctQuadraticApplications_NotOneSiblingObject()
    {
        var claim = BuildClaim();
        string surface = string.Join("\n", new[] { claim.Name, claim.DisplayName, claim.Summary }
            .Concat(claim.Children.Select(child => $"{child.DisplayName}\n{child.Summary}")));

        Assert.Contains("two distinct positive-b quadratic applications", surface,
            StringComparison.OrdinalIgnoreCase);
    }

    [Theory]
    [InlineData(1.0)]
    [InlineData(0.05)]
    public void EpAnchorIsTheAbsorptionRungBetweenTheUncoupledRates(double gamma0)
    {
        Assert.True(BuildClaim().UncoupledRatesAreRungsOneAndThreeAroundTheAnchor(gamma0));
    }

    [Fact]
    public void EpRoots_PinPositiveDecaySignSumAndAnchorB()
    {
        const double gamma0 = 0.7;
        const double q = 3.0;
        const double gEff = 0.8;
        var claim = BuildClaim();

        var lambda = claim.EpLiouvillianRoots(gamma0, q, gEff);
        var decay = claim.EpDecayRoots(gamma0, q, gEff);

        Assert.Equal(-4.0 * gamma0, lambda.Plus.Real, precision: 13);
        Assert.Equal(-4.0 * gamma0, lambda.Minus.Real, precision: 13);
        Assert.Equal(-lambda.Plus.Real, decay.FromLambdaPlus.Real, precision: 13);
        Assert.Equal(-lambda.Plus.Imaginary, decay.FromLambdaPlus.Imaginary, precision: 13);
        Assert.Equal(-lambda.Minus.Real, decay.FromLambdaMinus.Real, precision: 13);
        Assert.Equal(-lambda.Minus.Imaginary, decay.FromLambdaMinus.Imaginary, precision: 13);
        Assert.Equal(8.0 * gamma0, (decay.FromLambdaPlus + decay.FromLambdaMinus).Real, precision: 13);
        Assert.Equal(0.0, (decay.FromLambdaPlus + decay.FromLambdaMinus).Imaginary, precision: 13);
        Assert.Equal(4.0 * gamma0, claim.EpAnchorB(gamma0), precision: 13);
        Assert.True(claim.EpAnchorB(gamma0) > 0.0);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-0.1)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void EpF95Surface_RejectsNonPositiveGamma0(double gamma0)
    {
        var claim = BuildClaim();
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EpAnchorB(gamma0));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EpLiouvillianRoots(gamma0, 3.0, 0.8));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EpDecayRoots(gamma0, 3.0, 0.8));
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-0.8)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void EpF95Surface_RejectsInvalidGEff(double gEff)
    {
        var claim = BuildClaim();
        Assert.Throws<ArgumentOutOfRangeException>(() => TransitionBridgeF95SiblingClaim.QEp(gEff));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EpLiouvillianRoots(1.0, 3.0, gEff));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EpClockAngle(1.0, 3.0, gEff));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EpF95Angle(1.0, 3.0, gEff));
    }

    [Theory]
    [InlineData(-0.1)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void EpF95Surface_RejectsInvalidQ(double q)
    {
        var claim = BuildClaim();
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EpLiouvillianRoots(1.0, q, 0.8));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EpClockAngle(1.0, q, 0.8));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EpF95Angle(1.0, q, 0.8));
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

    [Fact]
    public void FragileBridge_IsNotIdentifiedWithTheToyTwoLevelEp()
    {
        var claim = BuildClaim();
        string surface = string.Join("\n", new[] { claim.Name, claim.DisplayName, claim.Summary }
            .Concat(claim.Children.Select(child => $"{child.DisplayName}\n{child.Summary}")));

        Assert.Contains("F86 toy 2×2 EP", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("FRAGILE_BRIDGE spectral-abscissa axis departure", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("EP/Hopf/Jordan character OPEN", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("EP (FRAGILE_BRIDGE", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("FRAGILE_BRIDGE (the EP", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("Sigma-gamma Hopf", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("tracked axis departure", surface, StringComparison.OrdinalIgnoreCase);
    }
}
