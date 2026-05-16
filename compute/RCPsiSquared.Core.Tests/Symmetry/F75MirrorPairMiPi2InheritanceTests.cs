using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F75MirrorPairMiPi2InheritanceTests
{
    private static F75MirrorPairMiPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var qubitAnchor = new QubitDimensionalAnchorClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubitAnchor);
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        return new F75MirrorPairMiPi2Inheritance(
            ladder,
            new F71MirrorSymmetryPi2Inheritance(new Pi2DyadicLadderClaim()),
            f65);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void MaxMIPerPair_IsExactlyTwoBits()
    {
        // MI saturates at 2 bits when p_ℓ = 1/2 (Bell-state mirror-pair).
        // = a_0 on the dyadic ladder = polynomial root d.
        Assert.Equal(2.0, BuildClaim().MaxMIPerPair, precision: 14);
    }

    [Fact]
    public void DomainUpperBound_IsExactlyOneHalf_FromBilinearApex()
    {
        Assert.Equal(0.5, BuildClaim().DomainUpperBound, precision: 14);
    }

    [Fact]
    public void TwoCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().TwoCoefficient, precision: 14);
    }

    [Fact]
    public void MirrorSigns_IsZ2Pair()
    {
        Assert.Equal(new[] { +1, -1 }, BuildClaim().MirrorSigns);
    }

    [Theory]
    [InlineData(2, 1)]    // ⌊2/2⌋ = 1
    [InlineData(3, 1)]    // ⌊3/2⌋ = 1
    [InlineData(4, 2)]
    [InlineData(5, 2)]
    [InlineData(7, 3)]
    [InlineData(13, 6)]
    public void MirrorPairCount_IsFloorNOver2(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().MirrorPairCount(N));
    }

    [Theory]
    [InlineData(2, 2.0)]    // 1 pair × 2 bits
    [InlineData(4, 4.0)]
    [InlineData(5, 4.0)]
    [InlineData(8, 8.0)]
    public void MaxTotalMM_IsPairCountTimesTwo(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().MaxTotalMM(N), precision: 12);
    }

    [Fact]
    public void MIPerPair_AtZero_IsZero()
    {
        Assert.Equal(0.0, BuildClaim().MIPerPair(0.0), precision: 14);
        Assert.True(BuildClaim().MIAtZeroIsZero());
    }

    [Fact]
    public void MIPerPair_AtBilinearApex_IsExactlyTwoBits()
    {
        // MI(1/2) = 2·h(1/2) − h(1) = 2·1 − 0 = 2 bits exactly
        Assert.Equal(2.0, BuildClaim().MIPerPair(0.5), precision: 12);
        Assert.True(BuildClaim().MIAtBilinearApexEqualsMaxMIPerPair());
    }

    [Theory]
    // f(p) = 2 h(p) − h(2p), at sample p values
    // h(0.1) ≈ 0.4690, h(0.2) ≈ 0.7219; f(0.1) = 2·0.4690 − 0.7219 = 0.2161
    // h(0.25) ≈ 0.8113, h(0.5) = 1; f(0.25) = 2·0.8113 − 1 = 0.6226
    // h(0.4) ≈ 0.9710, h(0.8) ≈ 0.7219; f(0.4) = 2·0.9710 − 0.7219 = 1.2200
    // h(1/2) = 1, h(1) = 0; f(1/2) = 2·1 − 0 = 2.0
    [InlineData(0.1, 0.2161)]
    [InlineData(0.25, 0.6226)]
    [InlineData(0.4, 1.2200)]
    [InlineData(0.5, 2.0)]
    public void MIPerPair_MatchesClosedForm(double p, double expected)
    {
        Assert.Equal(expected, BuildClaim().MIPerPair(p), precision: 3);
    }

    [Fact]
    public void MIPerPair_OutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MIPerPair(-0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MIPerPair(0.6));
    }

    [Theory]
    // F75 verified table (analytic MM(0) values):
    // N=5  k=1 → 0.800
    // N=5  k=2 → 1.245 (maximum at small N)
    // N=5  k=3 → 0.918
    // N=7  k=4 → 1.245 (resonance enhancement)
    [InlineData(5, 1, 0.800)]
    [InlineData(5, 2, 1.245)]
    [InlineData(5, 3, 0.918)]
    [InlineData(7, 4, 1.245)]
    public void BondingModeMMAtZero_MatchesF75VerifiedTable(int N, int k, double expected)
    {
        Assert.Equal(expected, BuildClaim().BondingModeMMAtZero(N, k), precision: 2);
    }

    [Theory]
    [InlineData(5, 2, 0)]   // sin(π·2·1/6)² · 2/6 = sin(π/3)² · 1/3 = 0.75·1/3 = 0.25
    [InlineData(5, 2, 2)]   // sin(π·2·3/6)² · 1/3 = sin(π)² · 1/3 = 0  (center site)
    public void BondingModePopulation_MatchesClosedForm(int N, int k, int site)
    {
        var f = BuildClaim();
        double expected = (2.0 / (N + 1)) * Math.Pow(Math.Sin(Math.PI * k * (site + 1) / (N + 1)), 2);
        Assert.Equal(expected, f.BondingModePopulation(N, k, site), precision: 12);
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, new QubitDimensionalAnchorClaim());
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        Assert.Throws<ArgumentNullException>(() =>
            new F75MirrorPairMiPi2Inheritance(null!, new F71MirrorSymmetryPi2Inheritance(new Pi2DyadicLadderClaim()), f65));
    }

    [Fact]
    public void Constructor_NullF71_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, new QubitDimensionalAnchorClaim());
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        Assert.Throws<ArgumentNullException>(() =>
            new F75MirrorPairMiPi2Inheritance(ladder, null!, f65));
    }

    [Fact]
    public void Constructor_NullF65_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F75MirrorPairMiPi2Inheritance(new Pi2DyadicLadderClaim(), new F71MirrorSymmetryPi2Inheritance(new Pi2DyadicLadderClaim()), null!));
    }

    [Fact]
    public void MirrorPairCount_NLessThan2_Throws()
    {
        // Delegated to F71.IndependentComponentCount which requires N ≥ 2.
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().MirrorPairCount(1));
    }
}
