using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F84ThermalAmplitudeDampingPi2InheritanceRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterPi2I4MemoryLoop()
            .RegisterF81Pi2Inheritance()
            .RegisterF82T1AmplitudeDampingPi2Inheritance();

    [Fact]
    public void RegisterF84_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF84ThermalAmplitudeDampingPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F84ThermalAmplitudeDampingPi2Inheritance>());
    }

    [Fact]
    public void RegisterF84_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF84ThermalAmplitudeDampingPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F84ThermalAmplitudeDampingPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF84_AncestorsContainF82()
    {
        var registry = BuildBaseRegistry()
            .RegisterF84ThermalAmplitudeDampingPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F84ThermalAmplitudeDampingPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F82T1AmplitudeDampingPi2Inheritance), ancestors);
        Assert.Contains(typeof(F81Pi2Inheritance), ancestors);
    }

    [Fact]
    public void RegisterF84_AncestorsTransitivelyReachPolynomialFoundation()
    {
        var registry = BuildBaseRegistry()
            .RegisterF84ThermalAmplitudeDampingPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F84ThermalAmplitudeDampingPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(PolynomialFoundationClaim), ancestors);
    }

    [Theory]
    [InlineData(0.10, 0.0, 3, 0.6928)]
    [InlineData(0.10, 0.10, 3, 0.0)]
    [InlineData(0.10, 0.05, 3, 0.3464)]
    public void RegisterF84_AmplitudeDampingTableAcrossRegistry(double gCool, double gHeat, int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF84ThermalAmplitudeDampingPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F84ThermalAmplitudeDampingPi2Inheritance>().AmplitudeDampingNormUniform(gCool, gHeat, N), precision: 4);
    }

    [Fact]
    public void RegisterF84_RecoversF82AtZeroHeatingAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF84ThermalAmplitudeDampingPi2Inheritance()
            .Build();

        Assert.True(registry.Get<F84ThermalAmplitudeDampingPi2Inheritance>().RecoversF82AtZeroHeating(0.10, 4));
    }
}
