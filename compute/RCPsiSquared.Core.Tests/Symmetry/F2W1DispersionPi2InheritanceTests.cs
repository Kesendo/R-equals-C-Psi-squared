using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Spectrum;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F2W1DispersionPi2InheritanceTests
{
    private static F2W1DispersionPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var w1 = new W1Dispersion(N: 5, J: 1.0, gammaZero: 0.05);
        return new F2W1DispersionPi2Inheritance(ladder, w1);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void BandwidthPrefactor_IsExactlyFour()
    {
        // 4 = a_{-1} on dyadic ladder.
        Assert.Equal(4.0, BuildClaim().BandwidthPrefactor, precision: 14);
    }

    [Fact]
    public void HoppingFactor_IsExactlyTwo()
    {
        // 2 = a_0 (hopping coefficient).
        Assert.Equal(2.0, BuildClaim().HoppingFactor, precision: 14);
    }

    [Theory]
    // F2: ω_k = 4J·(1−cos(πk/N)) for k=1..N−1
    // N=2, J=1, k=1: 4·(1−cos(π/2)) = 4·1 = 4
    // N=3, J=1, k=1: 4·(1−cos(π/3)) = 4·(1−1/2) = 2
    // N=3, J=1, k=2: 4·(1−cos(2π/3)) = 4·(1−(−1/2)) = 6
    // N=4, J=1, k=1: 4·(1−cos(π/4)) = 4·(1−√2/2) ≈ 1.172
    [InlineData(2, 1.0, 1, 4.0)]
    [InlineData(3, 1.0, 1, 2.0)]
    [InlineData(3, 1.0, 2, 6.0)]
    [InlineData(4, 1.0, 1, 4.0 * (1.0 - 0.7071067811865476))]
    [InlineData(4, 0.5, 1, 2.0 * (1.0 - 0.7071067811865476))]   // J scaling
    public void Frequency_MatchesClosedForm(int N, double J, int k, double expected)
    {
        Assert.Equal(expected, BuildClaim().Frequency(N, J, k), precision: 12);
    }

    [Theory]
    [InlineData(2, 1)]
    [InlineData(3, 2)]
    [InlineData(5, 4)]
    [InlineData(10, 9)]
    public void ModeCount_EqualsNMinusOne(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().ModeCount(N));
    }

    [Fact]
    public void MatchesW1DispersionParent_HoldsExactly()
    {
        // F2's closed form must produce identical numerical values to W1Dispersion's
        // stored Frequencies for the parent's (N, J).
        Assert.True(BuildClaim().MatchesW1DispersionParent());
    }

    [Fact]
    public void Frequency_ScalesLinearlyInJ()
    {
        // ω_k ∝ J: doubling J should double the frequency.
        var f = BuildClaim();
        double w1J1 = f.Frequency(N: 5, J: 1.0, k: 2);
        double w1J2 = f.Frequency(N: 5, J: 2.0, k: 2);
        Assert.Equal(w1J1 * 2.0, w1J2, precision: 12);
    }

    [Fact]
    public void Frequency_AtKEqualsN_Throws()
    {
        // k=N is out of range (only k=1..N−1 are valid).
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Frequency(N: 5, J: 1.0, k: 5));
    }

    [Fact]
    public void Frequency_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Frequency(N: 1, J: 1.0, k: 1));
    }

    [Fact]
    public void Frequency_NegativeJ_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Frequency(N: 3, J: -1.0, k: 1));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var w1 = new W1Dispersion(N: 5, J: 1.0, gammaZero: 0.05);
        Assert.Throws<ArgumentNullException>(() =>
            new F2W1DispersionPi2Inheritance(null!, w1));
    }

    [Fact]
    public void Constructor_NullW1Dispersion_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F2W1DispersionPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
