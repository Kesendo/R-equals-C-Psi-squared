using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F41PalindromicTimePi2InheritanceTests
{
    private static F41PalindromicTimePi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(ladder, memoryLoop);
        return new F41PalindromicTimePi2Inheritance(ladder, f1);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void HoppingCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().HoppingCoefficient, precision: 14);
    }

    [Theory]
    // F41: ω_min = 4·J·sin²(π/(2N))
    // N=2, J=1: sin²(π/4) = 1/2 → ω_min = 4·1·(1/2) = 2
    // N=3, J=1: sin²(π/6) = 1/4 → ω_min = 4·1·(1/4) = 1
    // N=4, J=1: sin²(π/8) ≈ 0.1464 → ω_min ≈ 0.586
    [InlineData(2, 1.0, 2.0)]
    [InlineData(3, 1.0, 1.0)]
    [InlineData(3, 0.5, 0.5)]   // J scaling: ω_min ∝ J
    [InlineData(4, 1.0, 4.0 * 0.14644660940672627)]
    public void MinFrequency_MatchesClosedForm(int N, double J, double expected)
    {
        Assert.Equal(expected, BuildClaim().MinFrequency(N, J), precision: 12);
    }

    [Theory]
    // F41: t_Pi = π/(2·J·sin²(π/(2N)))
    // N=2, J=1: 2·sin²(π/4) = 1 → t_Pi = π
    // N=3, J=1: 2·sin²(π/6) = 1/2 → t_Pi = 2π ≈ 6.283
    // N=3, J=2: t_Pi = π / (4·(1/4)) = π
    [InlineData(2, 1.0, Math.PI)]
    [InlineData(3, 1.0, 2.0 * Math.PI)]
    [InlineData(3, 2.0, Math.PI)]
    public void PalindromicTime_MatchesClosedForm(int N, double J, double expected)
    {
        Assert.Equal(expected, BuildClaim().PalindromicTime(N, J), precision: 12);
    }

    [Fact]
    public void PalindromicTime_GrowsWithN()
    {
        // t_Pi ~ N²/J for large N; should increase with N at fixed J.
        var f = BuildClaim();
        double prev = f.PalindromicTime(2, 1.0);
        for (int N = 3; N <= 20; N++)
        {
            double curr = f.PalindromicTime(N, 1.0);
            Assert.True(curr > prev, $"t_Pi at N={N} ({curr}) should be > t_Pi at N-1 ({prev})");
            prev = curr;
        }
    }

    [Fact]
    public void PalindromicTime_InverselyScalesWithJ()
    {
        // t_Pi ∝ 1/J: doubling J should halve t_Pi.
        var f = BuildClaim();
        double tJ1 = f.PalindromicTime(5, 1.0);
        double tJ2 = f.PalindromicTime(5, 2.0);
        Assert.Equal(tJ1 / 2.0, tJ2, precision: 12);
    }

    [Theory]
    // Asymptotic: t_Pi → 2N²/(π·J) for large N.
    // N=20, J=1: asymptotic = 2·400/π ≈ 254.65
    // N=100, J=1: asymptotic = 2·10000/π ≈ 6366.2
    [InlineData(20, 1.0, 800.0 / Math.PI)]
    [InlineData(100, 1.0, 20000.0 / Math.PI)]
    public void AsymptoticPalindromicTime_MatchesClosedForm(int N, double J, double expected)
    {
        Assert.Equal(expected, BuildClaim().AsymptoticPalindromicTime(N, J), precision: 12);
    }

    [Fact]
    public void AsymptoticPalindromicTime_ApproachesActualAtLargeN()
    {
        // At N=100, the leading-order asymptotic should agree with the exact closed form
        // to ~O(1/N²) ≈ 0.01% relative error.
        var f = BuildClaim();
        double exact = f.PalindromicTime(100, 1.0);
        double asympt = f.AsymptoticPalindromicTime(100, 1.0);
        Assert.True(Math.Abs(exact - asympt) / asympt < 0.001,
            $"asymptotic-exact mismatch at N=100: exact={exact}, asympt={asympt}");
    }

    [Theory]
    [InlineData(2, 1.0)]
    [InlineData(3, 1.0)]
    [InlineData(5, 0.5)]
    [InlineData(7, 2.0)]
    public void PeriodFrequencyProductIsTwoPi_HoldsAcrossNJ(int N, double J)
    {
        Assert.True(BuildClaim().PeriodFrequencyProductIsTwoPi(N, J));
    }

    [Fact]
    public void PalindromicTime_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PalindromicTime(N: 1, J: 1.0));
    }

    [Fact]
    public void PalindromicTime_NonPositiveJ_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PalindromicTime(N: 3, J: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PalindromicTime(N: 3, J: -1.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(ladder, memoryLoop);
        Assert.Throws<ArgumentNullException>(() =>
            new F41PalindromicTimePi2Inheritance(null!, f1));
    }

    [Fact]
    public void Constructor_NullF1_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F41PalindromicTimePi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
