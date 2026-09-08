using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F43BandSffPairingPi2InheritanceTests
{
    private static F43BandSffPairingPi2Inheritance BuildClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        return new F43BandSffPairingPi2Inheritance(ladder, f1);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void TypedClaim_DescribesReflectedBands_NotInvariantIntegerWeightSectors()
    {
        var claim = BuildClaim();
        string surface = string.Join("\n", new[] { claim.Name, claim.DisplayName, claim.Summary }
            .Concat(claim.Children.Select(child => child.Summary)));

        Assert.Contains("decay-rate band", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("sector SFF", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("sector-weight", surface, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void TypedClaim_DescribesConstantEndpointSff_NotADeltaSpike()
    {
        var claim = BuildClaim();
        string surface = string.Join("\n", new[] { claim.Name, claim.DisplayName, claim.Summary }
            .Concat(claim.Children.Select(child => child.Summary)));

        Assert.Contains("constant", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("positive uniform", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("γ=0 is excluded", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("delta-spike", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("XOR sector", surface, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void ReflectionCoefficient_IsExactlyTwo()
    {
        Assert.Equal(2.0, BuildClaim().ReflectionCoefficient, precision: 14);
    }

    [Theory]
    [InlineData(0.0, 5, 5.0)]
    [InlineData(1.25, 5, 3.75)]
    [InlineData(2.5, 5, 2.5)]
    [InlineData(4.75, 5, 0.25)]
    public void ReflectedAverageLight_AllowsFractionalBands(double averageLight, int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().ReflectedAverageLight(averageLight, N), precision: 14);
    }

    [Theory]
    [InlineData(0.0, 5, 0.05, 0.5)]
    [InlineData(0.13, 5, 0.05, 0.37)]
    [InlineData(0.25, 5, 0.05, 0.25)]
    public void ReflectedDecayRate_UsesFullPalindromeSpan(
        double decayRate,
        int N,
        double gamma,
        double expected)
    {
        Assert.Equal(expected, BuildClaim().ReflectedDecayRate(decayRate, N, gamma), precision: 14);
    }

    [Theory]
    [InlineData(4, 2.0)]
    [InlineData(5, 2.5)]
    public void MirrorAxis_IsAverageLightMidpoint(int N, double expected)
    {
        Assert.Equal(expected, BuildClaim().MirrorAxis(N), precision: 14);
    }

    [Fact]
    public void SelfReflection_AcceptsFractionalOddNMidpoint()
    {
        var claim = BuildClaim();
        Assert.True(claim.IsSelfReflected(2.5, N: 5));
        Assert.False(claim.IsSelfReflected(2.0, N: 5));
    }

    [Theory]
    [InlineData(3, 0.05, 0.3)]
    [InlineData(5, 0.1, 1.0)]
    [InlineData(7, 0.05, 0.7)]
    public void EndpointPartnerRate_Equals2NGamma(int N, double gamma, double expected)
    {
        Assert.Equal(expected, BuildClaim().EndpointPartnerRate(N, gamma), precision: 14);
    }

    [Theory]
    [InlineData(1, 2)]
    [InlineData(3, 4)]
    [InlineData(5, 6)]
    [InlineData(12, 13)]
    public void EndpointMultiplicity_IsNPlusOne(int N, int expected)
    {
        Assert.Equal(expected, BuildClaim().EndpointMultiplicity(N, gamma: 0.05));
    }

    [Theory]
    [InlineData(-11.0)]
    [InlineData(0.0)]
    [InlineData(0.125)]
    [InlineData(1000.0)]
    public void EndpointSff_IsConstantAtEveryFiniteTime(double time)
    {
        var claim = BuildClaim();
        Assert.Equal(1.0, claim.EndpointNormalizedFrequencySff(gamma: 0.05, time), precision: 14);
        Assert.Equal(36.0, claim.EndpointUnnormalizedFrequencySff(N: 5, gamma: 0.05, time), precision: 14);
    }

    [Fact]
    public void ExactEndpointApis_RejectZeroDephasingWhereTheEndpointCollapsesIntoTheKernel()
    {
        var claim = BuildClaim();
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EndpointPartnerRate(N: 2, gamma: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EndpointMultiplicity(N: 2, gamma: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EndpointNormalizedFrequencySff(gamma: 0.0, time: 1.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EndpointUnnormalizedFrequencySff(N: 2, gamma: 0.0, time: 1.0));
    }

    [Fact]
    public void InvalidBandInputs_Throw()
    {
        var claim = BuildClaim();
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.ReflectedAverageLight(-0.1, 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.ReflectedAverageLight(5.1, 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.ReflectedAverageLight(double.NaN, 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.ReflectedDecayRate(-0.1, 5, 0.05));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.ReflectedDecayRate(0.6, 5, 0.05));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.ReflectedDecayRate(0.1, 5, double.PositiveInfinity));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.EndpointNormalizedFrequencySff(gamma: 0.05, time: double.NaN));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        Assert.Throws<ArgumentNullException>(() =>
            new F43BandSffPairingPi2Inheritance(null!, f1));
    }

    [Fact]
    public void Constructor_NullF1_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F43BandSffPairingPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }
}
