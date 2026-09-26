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

    // The two floating routes below are deterministic functions of the same bit-identical inputs
    // (one f = e^(−4γt), or one s), so their difference is not noise to be bounded but a fixed
    // reading: the signed distance in units in the last place, pinned exactly per point. Both
    // routes evaluate positive terms only, with 4 + 5 (F25 against the family) and 5 + 1 (the
    // weights against s/3) correctly rounded operations, so a few ulp is the ceiling; the pins
    // record what the points actually give. The F25 pins read f through this platform's Math.Exp;
    // a different libm can move f and with it the offsets.
    private static long UlpOffset(double reference, double value) =>
        BitConverter.DoubleToInt64Bits(value) - BitConverter.DoubleToInt64Bits(reference);

    [Theory]
    [InlineData(0.05, 0.0)]
    [InlineData(0.05, 0.3)]
    [InlineData(0.05, 3.0)]
    [InlineData(0.05, 30.0)]
    [InlineData(1.0, 0.7)]
    public void BellPlusMember_HasTheWeightsOfF25_Exactly(double gamma, double t)
    {
        // s = 1 gives w₀ = 0.5/3 and w₁ = 1/6, the same binary64 as F25's 1/6; the family then
        // evaluates (1/6)·f + (1/6)·f·f·f in its own order, bit for bit.
        double f = Math.Exp(-4.0 * gamma * t);
        double expected = (1.0 / 6.0) * f + (1.0 / 6.0) * f * f * f;
        Assert.True(ApproachFamilyCarrierClaim.CPsi(1.0, gamma, t) == expected);
        Assert.True(0.5 / 3.0 == 1.0 / 6.0);
    }

    [Theory]
    [InlineData(0.05, 0.0, 0L)]
    [InlineData(0.05, 0.3, 1L)]
    [InlineData(0.05, 3.0, -1L)]
    [InlineData(0.05, 30.0, 0L)]
    [InlineData(1.0, 0.7, 0L)]
    public void BellPlusMember_AgainstF25ClosedForm_ReadsAFixedUlpOffset(double gamma, double t, long ulps)
    {
        // F25 evaluates f(1 + f²)/6, the family (1/6)f + (1/6)f³: the same polynomial, two orders.
        var c = BuildClaim();
        long offset = UlpOffset(c.F25.CPsiAtTime(gamma, t), ApproachFamilyCarrierClaim.CPsi(1.0, gamma, t));
        Assert.True(offset == ulps, $"ulp offset {offset}, pinned {ulps}");
    }

    [Theory]
    [InlineData(0.1, 1L)]
    [InlineData(0.5, 1L)]
    [InlineData(0.75, 0L)]
    [InlineData(0.9, 0L)]
    [InlineData(1.0, 0L)]
    public void InitialValue_IsOneThirdOfTheConcurrence_ToAFixedUlpOffset(double s, long ulps)
    {
        var (w0, w1) = ApproachFamilyCarrierClaim.Weights(s);
        long offset = UlpOffset(s / 3.0, w0 + w1);
        Assert.True(offset == ulps, $"ulp offset {offset}, pinned {ulps}");
        Assert.True(w0 > 0.0);
    }

    [Fact]
    public void ZeroMember_HasNoTermAtAll()
    {
        Assert.Equal((0.0, 0.0), ApproachFamilyCarrierClaim.Weights(0.0));
        Assert.True(ApproachFamilyCarrierClaim.CPsi(0.0, 0.05, 1.0) == 0.0);
    }

    [Theory]
    [InlineData(0.5, false)]
    [InlineData(0.74, false)]
    [InlineData(0.76, true)]
    [InlineData(1.0, true)]
    public void DownwardQuarterCrossing_OccursIffSAboveThreeQuarters(double s, bool crosses)
    {
        // CΨ decreases monotonically from s/3 to 0, so it crosses ¼ on the way down iff it
        // starts above ¼; at γ = 0 it never moves.
        const double gamma = 0.05;
        double start = ApproachFamilyCarrierClaim.CPsi(s, gamma, 0.0);
        double late = ApproachFamilyCarrierClaim.CPsi(s, gamma, 20.0 / gamma);
        Assert.True(late < 0.25);
        Assert.Equal(crosses, start > 0.25);
        double previous = start;
        for (int k = 1; k <= 200; k++)
        {
            double next = ApproachFamilyCarrierClaim.CPsi(s, gamma, 0.1 * k);
            Assert.True(next < previous);
            previous = next;
        }
        Assert.True(ApproachFamilyCarrierClaim.CPsi(s, 0.0, 50.0) == start);
    }

    [Fact]
    public void InvalidInputs_AreRejected()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => ApproachFamilyCarrierClaim.Weights(-0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => ApproachFamilyCarrierClaim.Weights(1.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => ApproachFamilyCarrierClaim.Weights(double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => ApproachFamilyCarrierClaim.CPsi(0.5, -1.0, 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => ApproachFamilyCarrierClaim.CPsi(0.5, 0.05, -1.0));
    }
}
