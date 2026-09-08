using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F43BandSffPairingPi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            .RegisterF88bPopcountCoherence()
            .RegisterF88bStaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterF38BitAInvolutionInheritance()
            .RegisterF63BitAReference()
            .RegisterF38Pi2InvolutionPi2Inheritance()
            .RegisterF63LCommutesPi2Pi2Inheritance()
            .RegisterF61BitAParityPi2Inheritance()
            .RegisterF1Pi2Inheritance();

    [Fact]
    public void RegisterF43_AddsBandClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF43BandSffPairingPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F43BandSffPairingPi2Inheritance>());
    }

    [Fact]
    public void RegisterF43_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF43BandSffPairingPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F43BandSffPairingPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF43_AncestorsContainF1AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF43BandSffPairingPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F43BandSffPairingPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1Pi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
    }

    [Fact]
    public void RegisterF43_ReflectsFractionalAverageLightAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF43BandSffPairingPi2Inheritance()
            .Build();

        Assert.Equal(3.75,
            registry.Get<F43BandSffPairingPi2Inheritance>()
                .ReflectedAverageLight(averageLight: 1.25, N: 5),
            precision: 14);
    }

    [Fact]
    public void RegisterF43_EndpointIsNPlusOneZeroFrequencyModes()
    {
        var registry = BuildBaseRegistry()
            .RegisterF43BandSffPairingPi2Inheritance()
            .Build();

        var f43 = registry.Get<F43BandSffPairingPi2Inheritance>();
        Assert.Equal(6, f43.EndpointMultiplicity(N: 5, gamma: 0.05));
        Assert.Equal(1.0, f43.EndpointNormalizedFrequencySff(gamma: 0.05, time: 17.0), precision: 14);
        Assert.Equal(36.0, f43.EndpointUnnormalizedFrequencySff(N: 5, gamma: 0.05, time: 17.0), precision: 14);
        Assert.Throws<ArgumentOutOfRangeException>(() => f43.EndpointMultiplicity(N: 5, gamma: 0.0));
    }

    [Fact]
    public void RegisterF43_WithoutF1_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterPi2I4MemoryLoop()
                .RegisterF43BandSffPairingPi2Inheritance()
                .Build());
    }
}
