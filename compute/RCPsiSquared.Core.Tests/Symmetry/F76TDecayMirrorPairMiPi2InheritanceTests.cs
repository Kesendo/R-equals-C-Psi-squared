using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F76TDecayMirrorPairMiPi2InheritanceTests
{
    private static F76TDecayMirrorPairMiPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f71 = new F71MirrorSymmetryPi2Inheritance();
        var f75 = new F75MirrorPairMiPi2Inheritance(ladder, f71);
        return new F76TDecayMirrorPairMiPi2Inheritance(ladder, f75);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void DecayRateCoefficient_IsExactlyFour()
    {
        // 4 in e^{-4γ₀t} = a_{-1} on dyadic ladder.
        Assert.Equal(4.0, BuildClaim().DecayRateCoefficient, precision: 14);
    }

    [Fact]
    public void LadderIndexForDecayRate_IsNegativeOne()
    {
        Assert.Equal(-1, BuildClaim().LadderIndexForDecayRate);
    }

    [Theory]
    [InlineData(0.0, 0.0, 1.0)]              // λ(0) = 1
    [InlineData(0.05, 0.0, 1.0)]
    [InlineData(0.05, 0.1, 0.9802)]          // F76 0.93 envelope check
    [InlineData(0.025, 0.1, 0.99005)]        // halved γ₀
    [InlineData(0.1, 0.1, 0.9608)]           // doubled γ₀
    public void Lambda_MatchesClosedForm(double gammaZero, double t, double expected)
    {
        Assert.Equal(expected, BuildClaim().Lambda(gammaZero, t), precision: 4);
    }

    [Fact]
    public void Lambda_LargeT_GoesToZero()
    {
        Assert.Equal(0.0, BuildClaim().Lambda(0.05, 1000.0), precision: 6);
    }

    [Fact]
    public void RecoversF75AtZero_HoldsAcrossPRange()
    {
        var f = BuildClaim();
        Assert.True(f.RecoversF75AtZero(0.0));
        Assert.True(f.RecoversF75AtZero(0.1));
        Assert.True(f.RecoversF75AtZero(0.25));
        Assert.True(f.RecoversF75AtZero(0.4));
        Assert.True(f.RecoversF75AtZero(0.5));
    }

    [Fact]
    public void MIPairAtTime_AtTZero_EqualsF75()
    {
        // MI_pair(p, t=0) = MI_F75(p) for any p ∈ [0, 1/2]
        var ladder = new Pi2DyadicLadderClaim();
        var f71 = new F71MirrorSymmetryPi2Inheritance();
        var f75 = new F75MirrorPairMiPi2Inheritance(ladder, f71);
        var f76 = new F76TDecayMirrorPairMiPi2Inheritance(ladder, f75);

        for (double p = 0.05; p <= 0.5; p += 0.05)
        {
            Assert.Equal(f75.MIPerPair(p), f76.MIPairAtTime(p, 0.05, 0.0), precision: 12);
        }
    }

    [Fact]
    public void MIPairAtTime_AtBilinearApex_AndZeroT_IsTwoBits()
    {
        // p=1/2, t=0: F76 recovers F75 saturation = 2 bits.
        Assert.Equal(2.0, BuildClaim().MIPairAtTime(0.5, 0.05, 0.0), precision: 10);
    }

    [Fact]
    public void MIPairAtTime_DecreasesWithT()
    {
        // For fixed p, MI should decrease as t grows (decoherence).
        var f = BuildClaim();
        double mi0 = f.MIPairAtTime(0.25, 0.05, 0.0);
        double mi1 = f.MIPairAtTime(0.25, 0.05, 0.1);
        double mi2 = f.MIPairAtTime(0.25, 0.05, 1.0);

        Assert.True(mi1 < mi0);
        Assert.True(mi2 < mi1);
    }

    [Theory]
    // F76 verified table (γ₀=0.05, t=0.1):
    // | N | k | pure-dephasing MM/MM(0) |
    // | 5 | 2 | 0.936 |
    // | 7 | 2 | 0.932 |
    // | 9 | 4 | 0.934 |
    // | 11| 4 | 0.933 |
    // | 13| 4 | 0.928 |
    [InlineData(5, 2, 0.936)]
    [InlineData(7, 2, 0.932)]
    [InlineData(9, 4, 0.934)]
    public void EnvelopeRatioForBondingMode_MatchesF76VerifiedTable(int N, int k, double expected)
    {
        Assert.Equal(expected, BuildClaim().EnvelopeRatioForBondingMode(N, k, 0.05, 0.1), precision: 2);
    }

    [Fact]
    public void Lambda_NegativeT_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Lambda(0.05, -0.1));
    }

    [Fact]
    public void Lambda_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().Lambda(-0.05, 0.1));
    }

    [Fact]
    public void MIPairAtTime_OutOfPRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MIPairAtTime(-0.1, 0.05, 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MIPairAtTime(0.6, 0.05, 0.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f71 = new F71MirrorSymmetryPi2Inheritance();
        var f75 = new F75MirrorPairMiPi2Inheritance(ladder, f71);
        Assert.Throws<ArgumentNullException>(() =>
            new F76TDecayMirrorPairMiPi2Inheritance(null!, f75));
    }

    [Fact]
    public void Constructor_NullF75_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F76TDecayMirrorPairMiPi2Inheritance(ladder, null!));
    }
}
