using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F84ThermalAmplitudeDampingPi2InheritanceTests
{
    private static F84ThermalAmplitudeDampingPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f81 = new F81Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, mirror, memoryLoop);
        var f82 = new F82T1AmplitudeDampingPi2Inheritance(ladder, f81);
        return new F84ThermalAmplitudeDampingPi2Inheritance(ladder, f82);
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
    [InlineData(2, 2.0)]
    [InlineData(3, 4.0)]
    [InlineData(5, 16.0)]
    public void ScalingFactor_DelegatesToF82(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().ScalingFactor(N), precision: 12);
    }

    [Theory]
    // F84 verified table from ANALYTICAL_FORMULAS (N=3):
    [InlineData(0.10, 0.00, 3, 0.6928)]    // cooling only (= F82)
    [InlineData(0.00, 0.10, 3, 0.6928)]    // heating only (same magnitude)
    [InlineData(0.10, 0.10, 3, 0.0000)]    // detailed balance
    [InlineData(0.10, 0.05, 3, 0.3464)]    // net cooling
    [InlineData(0.05, 0.10, 3, 0.3464)]    // net heating
    [InlineData(0.20, 0.05, 3, 1.0392)]    // strong cooling
    public void AmplitudeDampingNormUniform_MatchesF84VerifiedTable(double gCool, double gHeat, int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().AmplitudeDampingNormUniform(gCool, gHeat, N), precision: 4);
    }

    [Fact]
    public void AmplitudeDampingNormUniform_DetailedBalance_IsZero()
    {
        // γ_↓ = γ_↑ → Δγ = 0 → violation = 0 (T → ∞ limit)
        Assert.Equal(0.0, BuildClaim().AmplitudeDampingNormUniform(0.10, 0.10, 5), precision: 14);
    }

    [Theory]
    [InlineData(0.10, 3)]
    [InlineData(0.05, 5)]
    public void RecoversF82AtZeroHeating_AcrossNAndGamma(double gammaCool, int N)
    {
        Assert.True(BuildClaim().RecoversF82AtZeroHeating(gammaCool, N));
    }

    [Theory]
    [InlineData(0.10, 3)]
    [InlineData(0.05, 5)]
    public void DetailedBalanceGivesZeroViolation_AcrossNAndGamma(double gamma, int N)
    {
        Assert.True(BuildClaim().DetailedBalanceGivesZeroViolation(gamma, N));
    }

    [Fact]
    public void PauliChannelViolation_IsExactlyZero()
    {
        // Pauli-Channel Cancellation Lemma: D[Z], D[X], D[Y] all Π²-symmetric
        Assert.Equal(0.0, F84ThermalAmplitudeDampingPi2Inheritance.PauliChannelViolation, precision: 15);
    }

    [Fact]
    public void EstimateNetCoolingFromViolation_RecoversInputDeltaGamma()
    {
        // For Δγ = 0.05 net cooling at N=3: violation = 0.3464; inverse → 0.05
        var f = BuildClaim();
        double violation = f.AmplitudeDampingNormUniform(0.10, 0.05, 3);
        double recoveredDelta = f.EstimateNetCoolingFromViolation(violation, 3);
        Assert.Equal(0.05, recoveredDelta, precision: 4);
    }

    [Fact]
    public void AmplitudeDampingNormNonUniform_MatchesUniformAtUniformInput()
    {
        var f = BuildClaim();
        var cool = new[] { 0.10, 0.10, 0.10 };
        var heat = new[] { 0.05, 0.05, 0.05 };
        // Each Δγ = 0.05; sum² = 3·0.0025 = 0.0075; sqrt = 0.0866; ·4 = 0.3464
        Assert.Equal(f.AmplitudeDampingNormUniform(0.10, 0.05, 3),
                     f.AmplitudeDampingNormNonUniform(cool, heat),
                     precision: 4);
    }

    [Fact]
    public void AmplitudeDampingNormUniform_NegativeRates_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().AmplitudeDampingNormUniform(-0.05, 0.10, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().AmplitudeDampingNormUniform(0.10, -0.05, 3));
    }

    [Fact]
    public void AmplitudeDampingNormNonUniform_MismatchedArrays_Throws()
    {
        Assert.Throws<ArgumentException>(() =>
            BuildClaim().AmplitudeDampingNormNonUniform(new[] { 0.1, 0.2 }, new[] { 0.05 }));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f81 = new F81Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, mirror, memoryLoop);
        var f82 = new F82T1AmplitudeDampingPi2Inheritance(ladder, f81);
        Assert.Throws<ArgumentNullException>(() =>
            new F84ThermalAmplitudeDampingPi2Inheritance(null!, f82));
    }

    [Fact]
    public void Constructor_NullF82_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F84ThermalAmplitudeDampingPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
