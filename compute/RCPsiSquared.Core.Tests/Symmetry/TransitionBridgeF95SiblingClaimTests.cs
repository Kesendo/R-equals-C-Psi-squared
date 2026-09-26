using System.Numerics;
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
    public void EpAngle_IsUndefinedBelowTheEp()
    {
        var c = BuildClaim();
        double gEff = 0.5; // Q_EP = 4
        Assert.True(double.IsNaN(c.EpClockAngle(1.0, 3.0, gEff)), "below the EP: no rotation");
        Assert.True(double.IsNaN(c.EpF95Angle(1.0, 3.0, gEff)), "below the EP: c < b², no F95 angle");
        Assert.False(c.EpClockAngleEqualsF95Angle(1.0, 3.0, gEff));
    }

    [Theory]
    [InlineData(0.5, 1.0)]
    [InlineData(1.0, 1.0)]
    [InlineData(2.0, 0.05)]
    [InlineData(4.0, 0.05)]
    public void EpAngle_WhereTheEpIsADouble_IsComputedExactlyZero(double gEff, double gamma0)
    {
        // Q_EP = 2/g_eff is exact in binary for these g_eff, so q·g_eff − 2 is exactly 0.0 and
        // both routes return 0.0 from their own arithmetic: the clock from a zero branch of the
        // roots, the F95 route from c = b² exactly. Nothing is special-cased at Q_EP.
        var c = BuildClaim();
        double qEp = TransitionBridgeF95SiblingClaim.QEp(gEff);
        Assert.True(qEp * gEff == 2.0);
        Assert.Equal(0.0, c.EpClockAngle(gamma0, qEp, gEff));
        Assert.Equal(0.0, c.EpF95Angle(gamma0, qEp, gEff));
        Assert.True(c.EpClockAngleEqualsF95Angle(gamma0, qEp, gEff));
    }

    [Fact]
    public void EpAngle_AtTheDoubleNearestToQEp_ReadsTheRoundingOfTheQuotient()
    {
        // Where 2/g_eff is not a double, QEp returns the nearest one and the angle is the angle
        // of those inputs. For a given double g_eff the EP is q* = 2/g_eff exactly, so the offset
        // is the rounding of that quotient (exact g_eff = 3, 5, 6 miss the EP the same way):
        // 2/fl(4/3) rounds to 1.5 with 1.5·fl(4/3) < 2, just below the EP (both routes NaN);
        // 2/fl(0.8) rounds to 2.5 with 2.5·fl(0.8) > 2, just above, where the clock reads
        // √(s·(q·g + 2))/4 = 2^(−27.5) ≈ 5.3·10⁻⁹ and the √-sensitive F95 route rounds c back to b², 0.0.
        var c = BuildClaim();
        double below = TransitionBridgeF95SiblingClaim.QEp(4.0 / 3.0);
        Assert.True(double.IsNaN(c.EpClockAngle(1.0, below, 4.0 / 3.0)));
        Assert.True(double.IsNaN(c.EpF95Angle(1.0, below, 4.0 / 3.0)));

        double above = TransitionBridgeF95SiblingClaim.QEp(0.8);
        double clock = c.EpClockAngle(1.0, above, 0.8);
        double s = Math.FusedMultiplyAdd(above, 0.8, -2.0);
        Assert.True(s == Math.ScaleB(1.0, -53), $"s = {s:R}");
        // s = 2⁻⁵³ and q·g + 2 rounds to 4, so the angle is arctan(√(2⁻⁵¹)/4) = 2^(−27.5) to
        // within the few roundings of Sqrt, Complex.Sqrt and Atan (each within an ulp).
        Assert.True(Math.Abs(clock / Math.Pow(2.0, -27.5) - 1.0) <= 4.0 * Math.ScaleB(1.0, -52), $"clock = {clock:R}");
        Assert.Equal(0.0, c.EpF95Angle(1.0, above, 0.8));
        Assert.True(clock <= TransitionBridgeF95SiblingClaim.AngleLawMargin * TransitionBridgeF95SiblingClaim.AngleRoundingLaw(0.0));
        Assert.True(c.EpClockAngleEqualsF95Angle(1.0, above, 0.8));

        // Exact g_eff, inexact quotient: 3 and 6 land below, 5 above.
        foreach (double g in new[] { 3.0, 6.0 })
            Assert.True(Math.FusedMultiplyAdd(TransitionBridgeF95SiblingClaim.QEp(g), g, -2.0) < 0.0
                        && double.IsNaN(c.EpClockAngle(1.0, TransitionBridgeF95SiblingClaim.QEp(g), g)));
        Assert.True(c.EpClockAngle(1.0, TransitionBridgeF95SiblingClaim.QEp(5.0), 5.0) > 0.0);
    }

    /// <summary>The exact t = ((q·g)² − 4)/16 of the given doubles, rounded once to a double:
    /// q·g is an exact dyadic rational, so t is computed in BigInteger arithmetic.</summary>
    private static double ExactT(double q, double g)
    {
        static (BigInteger M, int E) Split(double x)
        {
            long bits = BitConverter.DoubleToInt64Bits(x);
            int exponent = (int)((bits >> 52) & 0x7FF);
            long mantissa = bits & ((1L << 52) - 1);
            Assert.True(exponent != 0 && x > 0.0, "normal positive doubles only");
            return (new BigInteger(mantissa | (1L << 52)), exponent - 1075);
        }
        var (mq, eq) = Split(q);
        var (mg, eg) = Split(g);
        BigInteger m = mq * mg;
        int k = -2 * (eq + eg);
        Assert.True(k > 0);
        BigInteger numerator = m * m - (new BigInteger(4) << k);
        return Math.ScaleB((double)numerator, -k - 4);
    }

    [Fact]
    public void F95Route_FollowsItsRoundingLaw_AcrossSeventeenDecades()
    {
        // No exact route exists for the F95 angle near the EP (c/b² − 1 is formed by a rounded
        // division of two nearly equal numbers), so its deviation is held to a LAW: against the
        // exact angle of the given doubles, arctan(√t), the worst deviation per decade of t is a
        // steady fraction of E(t) = ε(1+t)/(2√t) from 10⁻¹⁶ to 10¹ (the branch of AngleRoundingLaw
        // that holds for t above ε/4; numpy replicas give about 0.8 to 1.06 per decade, so the band
        // [0.6, 1.25] rejects a law off by a factor of two either way). The clock route, which
        // reads the factored discriminant, stays within the same bound. Below ε/4 see
        // NextToTheEp_TheSideIsExact_AndTheErrorIsUnderItsCeiling.
        var c = BuildClaim();
        var rng = new Random(7);
        for (int decade = -16; decade <= 0; decade++)
        {
            double worstF95 = 0.0, worstClock = 0.0;
            for (int sample = 0; sample < 400; sample++)
            {
                double gamma0 = Math.Pow(10.0, -3.0 + 4.0 * rng.NextDouble());
                double g = Math.Pow(10.0, -1.0 + 2.0 * rng.NextDouble());
                double tTarget = Math.Pow(10.0, decade + rng.NextDouble());
                double q = Math.Sqrt(16.0 * tTarget + 4.0) / g;
                double t = ExactT(q, g);
                Assert.True(t > 0.0);
                double exact = Math.Atan(Math.Sqrt(t));
                double law = TransitionBridgeF95SiblingClaim.AngleRoundingLaw(t);
                worstF95 = Math.Max(worstF95, Math.Abs(c.EpF95Angle(gamma0, q, g) - exact) / law);
                worstClock = Math.Max(worstClock, Math.Abs(c.EpClockAngle(gamma0, q, g) - exact) / law);
                Assert.True(c.EpClockAngleEqualsF95Angle(gamma0, q, g), $"routes disagree beyond the law at q={q:R}, g={g:R}");
            }
            Assert.InRange(worstF95, 0.6, 1.25);
            Assert.InRange(worstClock, 0.0, 1.25);
        }
    }

    [Fact]
    public void NextToTheEp_TheSideIsExact_AndTheErrorIsUnderItsCeiling()
    {
        // The last few doubles around 2/g_eff: t = ((q·g)² − 4)/16 is of order ε or below, where
        // c/b² − 1 is quantized and √(ε(1+t)) is a ceiling, not a law. Two things are gated:
        // the side of the EP read by FMA(q, g, −2) agrees with the sign of the EXACT t (NaN below,
        // 0.0 on it), and above it both routes stay within 1.25·AngleRoundingLaw(t), the sweep's
        // margin: that is the √(ε(1+t)) ceiling for t below ε/4 and the law branch ε(1+t)/(2√t)
        // for the steps that land above it (most of them).
        var c = BuildClaim();
        var rng = new Random(3);
        int below = 0, above = 0;
        for (int sample = 0; sample < 500; sample++)
        {
            double g = Math.Pow(10.0, -1.0 + 2.0 * rng.NextDouble());
            double gamma0 = Math.Pow(10.0, -3.0 + 4.0 * rng.NextDouble());
            double q0 = TransitionBridgeF95SiblingClaim.QEp(g);
            for (int k = -3; k <= 3; k++)
            {
                double q = q0;
                for (int step = 0; step < Math.Abs(k); step++)
                    q = k > 0 ? Math.BitIncrement(q) : Math.BitDecrement(q);
                double t = ExactT(q, g);
                double clock = c.EpClockAngle(gamma0, q, g), f95 = c.EpF95Angle(gamma0, q, g);
                if (t < 0.0)
                {
                    below++;
                    Assert.True(double.IsNaN(clock) && double.IsNaN(f95), $"below the EP at q={q:R}, g={g:R}");
                }
                else if (t == 0.0)
                {
                    Assert.True(clock == 0.0 && f95 == 0.0);
                }
                else
                {
                    above++;
                    double exact = Math.Atan(Math.Sqrt(t));
                    double ceiling = 1.25 * TransitionBridgeF95SiblingClaim.AngleRoundingLaw(t);
                    Assert.True(Math.Abs(f95 - exact) <= ceiling, $"F95 route over its ceiling at q={q:R}, g={g:R}");
                    Assert.True(Math.Abs(clock - exact) <= ceiling, $"clock route over the ceiling at q={q:R}, g={g:R}");
                }
            }
        }
        Assert.True(below > 0 && above > 0, "both sides of the EP must be visited");
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
    public void ToyQuadratic_DoubleRootSitsAtPositiveDecay_AndTheFragileBoundaryCarriesItsScope()
    {
        // Two different double roots. The toy quadratic's sits at the midpoint of its decay roots,
        // z_decay = b = 4γ₀ > 0 (checked here from the roots at J = 0); the FRAGILE_BRIDGE collision
        // sits at zero decay, Re λ* = 0, since Σγ = 0 centres its palindrome at 0 (recomputed from the
        // C# builders in Core.Tests.F86.FragileBridgeThresholdTests). The boundary node is held to its
        // scope, the chain length.
        var claim = BuildClaim();
        var decay = claim.EpDecayRoots(0.05, 0.0, 1.0);
        Assert.Equal(claim.EpAnchorB(0.05), (decay.FromLambdaPlus.Real + decay.FromLambdaMinus.Real) / 2.0, 14);
        Assert.True(claim.EpAnchorB(0.05) > 0.0);

        string boundary = claim.Children.Single(child => child.DisplayName == "FRAGILE_BRIDGE boundary").Summary;
        Assert.Contains("two qubits per chain", boundary, StringComparison.Ordinal);
    }
}
