using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F82T1AmplitudeDampingPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterPi2I4MemoryLoop()
            .RegisterF81Pi2Inheritance();

    [Fact]
    public void RegisterF82_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F82T1AmplitudeDampingPi2Inheritance>());
    }

    [Fact]
    public void RegisterF82_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F82T1AmplitudeDampingPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF82_AncestorsContainBothLadderAndF81()
    {
        var registry = BuildBaseRegistry()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F82T1AmplitudeDampingPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(F81Pi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterF82_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F82T1AmplitudeDampingPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Theory]
    [InlineData(2, 2.0)]
    [InlineData(3, 4.0)]
    [InlineData(5, 16.0)]
    public void RegisterF82_ScalingFactorAcrossRegistry(int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F82T1AmplitudeDampingPi2Inheritance>().ScalingFactor(N), precision: 12);
    }

    [Theory]
    [InlineData(0.05, 3, 0.3464)]
    [InlineData(0.10, 4, 1.6000)]
    public void RegisterF82_T1DissipatorNormUniformAcrossRegistry(double gammaT1, int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F82T1AmplitudeDampingPi2Inheritance>().T1DissipatorNormUniform(gammaT1, N), precision: 4);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(5)]
    public void RegisterF82_RecoversF81AtZeroT1AcrossN(int N)
    {
        var registry = BuildBaseRegistry()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F82T1AmplitudeDampingPi2Inheritance>().RecoversF81AtZeroT1(N));
    }

    [Fact]
    public void RegisterF82_ForwardInverseRoundTripAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF82T1AmplitudeDampingPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F82T1AmplitudeDampingPi2Inheritance>().ForwardInverseRoundTrip(0.075, 4));
    }
}
