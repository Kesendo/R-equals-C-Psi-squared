using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F43SectorSffPairingPi2InheritanceTests
{
    private static F43SectorSffPairingPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        return new F43SectorSffPairingPi2Inheritance(ladder, f1);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void XorRateCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().XorRateCoefficient, precision: 14);
    }

    [Theory]
    // PartnerSector: w → N − w
    [InlineData(0, 3, 3)]
    [InlineData(1, 3, 2)]
    [InlineData(2, 3, 1)]
    [InlineData(3, 3, 0)]
    [InlineData(0, 4, 4)]
    [InlineData(2, 4, 2)]   // self-paired (N=4, w=N/2=2)
    [InlineData(5, 10, 5)]  // self-paired (N=10, w=N/2=5)
    public void PartnerSector_EqualsNMinusW(int w, int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().PartnerSector(w, N));
    }

    [Theory]
    [InlineData(0, 3, false)]
    [InlineData(1, 3, false)]
    [InlineData(2, 4, true)]    // N/2 = 2 = w (even N)
    [InlineData(5, 10, true)]   // N/2 = 5 = w (even N)
    [InlineData(3, 6, true)]    // N/2 = 3 = w
    public void IsSelfPaired_TrueOnlyAtMirrorAxisForEvenN(int w, int N, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsSelfPaired(w, N));
    }

    [Fact]
    public void IsSelfPaired_NoneForOddN()
    {
        // Odd N has no integer w that equals N/2 (half-integer); no self-pairs.
        var f = BuildClaim();
        for (int w = 0; w <= 5; w++)
        {
            Assert.False(f.IsSelfPaired(w, N: 5));
        }
    }

    [Theory]
    [InlineData(2, 1.0)]
    [InlineData(3, 1.5)]
    [InlineData(4, 2.0)]
    [InlineData(7, 3.5)]
    public void MirrorAxis_EqualsNOver2(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().MirrorAxis(N), precision: 14);
    }

    [Theory]
    [InlineData(2, false)]
    [InlineData(3, true)]
    [InlineData(4, false)]
    [InlineData(5, true)]
    [InlineData(7, true)]
    [InlineData(10, false)]
    public void IsHalfIntegerMirrorRegime_TrueForOddN(int N, bool expected)
    {
        Assert.Equal(expected, BuildClaim().IsHalfIntegerMirrorRegime(N));
    }

    [Theory]
    // XOR sector decay rate: 2·N·γ (per F50 the per-site rate is 2γ for w=1; XOR
    // sector is w=N where ALL sites carry X or Y, so the rate scales with N).
    [InlineData(3, 0.05, 0.3)]
    [InlineData(5, 0.1, 1.0)]
    [InlineData(7, 0.05, 0.7)]
    public void XorSectorRate_Equals2NGamma(int N, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().XorSectorRate(N, gammaZero), precision: 14);
    }

    [Fact]
    public void XorSectorSffValue_IsExactlyOne()
    {
        // K_freq(w = N, t) = 1 exactly (delta-spike from full eigenvalue degeneracy).
        Assert.Equal(1.0, F43SectorSffPairingPi2Inheritance.XorSectorSffValue, precision: 14);
    }

    [Theory]
    [InlineData(0, 3)]
    [InlineData(1, 5)]
    [InlineData(2, 4)]
    [InlineData(3, 7)]
    public void PartnerSumEqualsN_HoldsForAnyW(int w, int N)
    {
        Assert.True(BuildClaim().PartnerSumEqualsN(w, N));
    }

    [Fact]
    public void PartnerSector_OutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PartnerSector(w: -1, N: 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PartnerSector(w: 4, N: 3));
    }

    [Fact]
    public void PartnerSector_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().PartnerSector(w: 0, N: 0));
    }

    [Fact]
    public void XorSectorRate_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().XorSectorRate(N: 3, gammaZero: -0.05));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        Assert.Throws<ArgumentNullException>(() =>
            new F43SectorSffPairingPi2Inheritance(null!, f1));
    }

    [Fact]
    public void Constructor_NullF1_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F43SectorSffPairingPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
