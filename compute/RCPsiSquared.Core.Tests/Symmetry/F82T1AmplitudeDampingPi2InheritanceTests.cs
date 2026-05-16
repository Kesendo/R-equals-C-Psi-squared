using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F82T1AmplitudeDampingPi2InheritanceTests
{
    private static F82T1AmplitudeDampingPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f81 = new F81Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, mirror, memoryLoop);
        return new F82T1AmplitudeDampingPi2Inheritance(ladder, f81);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void Coefficient2_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().Coefficient2, precision: 14);
    }

    [Theory]
    [InlineData(2, 2.0)]   // 2^1 = a_0
    [InlineData(3, 4.0)]   // 2^2 = a_{-1}
    [InlineData(4, 8.0)]   // 2^3 = a_{-2}
    [InlineData(5, 16.0)]  // 2^4 = a_{-3}
    public void ScalingFactor_Is2PowerNMinus1(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().ScalingFactor(N), precision: 12);
    }

    [Theory]
    [InlineData(2, 0)]
    [InlineData(3, -1)]
    [InlineData(4, -2)]
    [InlineData(5, -3)]
    public void LadderIndexForScaling_Is2MinusN(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().LadderIndexForScaling(N));
    }

    [Theory]
    // F82 verified table from ANALYTICAL_FORMULAS:
    // N=3, γ_T1=0.05: violation = 0.05·√3·4 = 0.3464
    // N=3, γ_T1=0.10: violation = 0.10·√3·4 = 0.6928
    // N=3, γ_T1=1.00: violation = 1.00·√3·4 = 6.9282
    [InlineData(0.05, 3, 0.3464)]
    [InlineData(0.10, 3, 0.6928)]
    [InlineData(1.00, 3, 6.9282)]
    [InlineData(0.10, 2, 0.2828)]   // 0.1·√2·2 = 0.2828
    [InlineData(0.10, 4, 1.6000)]   // 0.1·2·8 = 1.6
    [InlineData(0.10, 5, 3.5777)]   // 0.1·√5·16 = 0.1·2.23607·16 = 3.57771
    public void T1DissipatorNormUniform_MatchesVerifiedTable(double gammaT1, int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().T1DissipatorNormUniform(gammaT1, N), precision: 4);
    }

    [Theory]
    // F82 non-uniform table:
    // N=3, single-site (0.10, 0, 0): √(0.01)·4 = 0.4000
    // N=3, two-site (0.10, 0.10, 0): √(0.02)·4 = 0.5657
    // N=3, non-uniform (0.05, 0.10, 0.15): √(0.035)·4 = 0.7483
    [InlineData(new[] { 0.10, 0.0, 0.0 }, 0.4000)]
    [InlineData(new[] { 0.10, 0.10, 0.0 }, 0.5657)]
    [InlineData(new[] { 0.05, 0.10, 0.15 }, 0.7483)]
    public void T1DissipatorNormNonUniform_MatchesVerifiedTable(double[] gammaT1Sites, double expected)
    {
        Assert.Equal(expected, BuildClaim().T1DissipatorNormNonUniform(gammaT1Sites), precision: 4);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void RecoversF81AtZeroT1_AcrossN(int N)
    {
        Assert.True(BuildClaim().RecoversF81AtZeroT1(N));
    }

    [Theory]
    [InlineData(0.6928, 3, 0.10)]
    [InlineData(0.3464, 3, 0.05)]
    [InlineData(1.6000, 4, 0.10)]
    public void EstimateT1FromViolation_RecoversInputGamma(double violation, int N, double expectedGammaT1)
    {
        Assert.Equal(expectedGammaT1, BuildClaim().EstimateT1FromViolation(violation, N), precision: 4);
    }

    [Theory]
    [InlineData(0.05, 3)]
    [InlineData(0.10, 5)]
    [InlineData(0.001, 4)]
    public void ForwardInverseRoundTrip_HoldsAcrossNAndGamma(double gammaT1, int N)
    {
        Assert.True(BuildClaim().ForwardInverseRoundTrip(gammaT1, N));
    }

    [Fact]
    public void T1DissipatorNormUniform_NLessThan2_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().T1DissipatorNormUniform(0.05, 1));
    }

    [Fact]
    public void T1DissipatorNormUniform_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().T1DissipatorNormUniform(-0.05, 3));
    }

    [Fact]
    public void T1DissipatorNormNonUniform_NullArray_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => BuildClaim().T1DissipatorNormNonUniform(null!));
    }

    [Fact]
    public void T1DissipatorNormNonUniform_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            BuildClaim().T1DissipatorNormNonUniform(new[] { 0.05, -0.10, 0.15 }));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f81 = new F81Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, mirror, memoryLoop);
        Assert.Throws<ArgumentNullException>(() =>
            new F82T1AmplitudeDampingPi2Inheritance(null!, f81));
    }

    [Fact]
    public void Constructor_NullF81_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F82T1AmplitudeDampingPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
