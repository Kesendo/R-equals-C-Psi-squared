using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F77MmSaturationPi2InheritanceTests
{
    private static F77MmSaturationPi2Inheritance BuildClaim() =>
        new F77MmSaturationPi2Inheritance(new Pi2DyadicLadderClaim(), new HalfAsStructuralFixedPointClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void SaturationBit_IsExactlyOne_FromSelfMirrorPivot()
    {
        // The 1-bit asymptotic saturation = a_1 = 1 = self-mirror pivot.
        Assert.Equal(1.0, BuildClaim().SaturationBit, precision: 14);
    }

    [Fact]
    public void LadderIndexForSaturation_IsSelfMirrorIndex()
    {
        var f = BuildClaim();
        Assert.Equal(1, f.LadderIndexForSaturation);
        Assert.True(f.LandsOnSelfMirrorPivot);
    }

    [Fact]
    public void TwoFactor_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().TwoFactor, precision: 14);
    }

    [Fact]
    public void HalfFactor_IsExactlyOneHalf()
    {
        Assert.Equal(0.5, BuildClaim().HalfFactor, precision: 14);
    }

    [Fact]
    public void InversionIdentityProduct_IsExactlyOne()
    {
        // The "2·(1/2) = 1" mechanism IS the Pi2 ladder inversion identity
        // a_0 · a_2 = 1 = a_1 = SaturationBit.
        var f = BuildClaim();
        Assert.Equal(1.0, f.InversionIdentityProduct, precision: 14);
    }

    [Fact]
    public void InversionIdentityHolds_BitExactly()
    {
        Assert.True(BuildClaim().InversionIdentityHolds);
    }

    [Fact]
    public void AsymptoticCorrectionCoefficient_MatchesClosedForm()
    {
        // 3 / (4 ln 2) = 1.0820...
        // Verified to 10⁻⁴ by N = 10⁴ in F77 closed form.
        var expected = 3.0 / (4.0 * Math.Log(2.0));
        Assert.Equal(expected, BuildClaim().AsymptoticCorrectionCoefficient, precision: 14);
        Assert.Equal(1.0820, BuildClaim().AsymptoticCorrectionCoefficient, precision: 4);
    }

    [Fact]
    public void FourFactorInCorrectionDenominator_IsExactlyFour()
    {
        // The "4" in 4(N+1) = a_{-1} on the dyadic ladder = d² for 1 qubit.
        Assert.Equal(4.0, BuildClaim().FourFactorInCorrectionDenominator, precision: 14);
    }

    [Theory]
    // Verified table from ANALYTICAL_FORMULAS F77:
    // | N   | k* | MM(0)    | (MM−1)·(N+1) |
    // | 101 | 82 | 1.01078  | 1.100        |
    // | 201 | 102| 1.00540  | 1.091        |
    // | 1001| 878| 1.00108  | 1.084        |
    // (Exact MM at finite N includes higher-order corrections beyond
    //  3/(4 ln 2) → use truncated O(1/(N+1)) form for direct comparison)
    [InlineData(101)]
    [InlineData(201)]
    [InlineData(501)]
    [InlineData(1001)]
    [InlineData(10001)]
    public void MmAtN_IsAboveOne_AndConvergesToOne(int N)
    {
        var f = BuildClaim();
        var mm = f.MmAtN(N);
        Assert.True(mm > 1.0);  // saturates from above
        Assert.True(mm < 1.02); // within 2% for N ≥ 101
    }

    [Theory]
    [InlineData(101)]
    [InlineData(1001)]
    [InlineData(10001)]
    public void RescaledDeviation_ConvergesToCorrectionCoefficient(int N)
    {
        // (MM(0)(N) − 1) · (N + 1) → 3/(4 ln 2) = 1.0820 as N → ∞
        // (Using truncated O(1/(N+1)) form, the rescaled deviation equals the
        // coefficient exactly modulo floating-point cancellation precision.)
        var f = BuildClaim();
        Assert.Equal(f.AsymptoticCorrectionCoefficient, f.RescaledDeviation(N), precision: 10);
    }

    [Fact]
    public void MmAtN_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MmAtN(0));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        Assert.Throws<ArgumentNullException>(() => new F77MmSaturationPi2Inheritance(null!, half));
        Assert.Throws<ArgumentNullException>(() => new F77MmSaturationPi2Inheritance(ladder, null!));
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var f = BuildClaim();
        Assert.NotNull(f.Half);
    }
}
