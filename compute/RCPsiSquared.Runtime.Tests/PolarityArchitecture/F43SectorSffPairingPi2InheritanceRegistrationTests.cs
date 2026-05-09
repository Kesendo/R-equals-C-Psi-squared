using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F43SectorSffPairingPi2InheritanceRegistrationTests
{
    private static ChainSystem DefaultChain() =>
        new ChainSystem(N: 5, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterPi2I4MemoryLoop()
            .RegisterF1Pi2Inheritance();

    [Fact]
    public void RegisterF43_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF43SectorSffPairingPi2Inheritance()
            .Build();

        Assert.True(registry.Contains<F43SectorSffPairingPi2Inheritance>());
    }

    [Fact]
    public void RegisterF43_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF43SectorSffPairingPi2Inheritance()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<F43SectorSffPairingPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF43_AncestorsContainF1AndPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterF43SectorSffPairingPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F43SectorSffPairingPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1Pi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
    }

    [Theory]
    [InlineData(1, 3, 2)]
    [InlineData(2, 4, 2)]   // self-paired
    [InlineData(0, 5, 5)]
    public void RegisterF43_PartnerSectorAcrossRegistry(int w, int N, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF43SectorSffPairingPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F43SectorSffPairingPi2Inheritance>().PartnerSector(w, N));
    }

    [Theory]
    [InlineData(3, 0.05, 0.3)]
    [InlineData(5, 0.1, 1.0)]
    public void RegisterF43_XorSectorRateAcrossRegistry(int N, double gammaZero, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF43SectorSffPairingPi2Inheritance()
            .Build();

        Assert.Equal(expected, registry.Get<F43SectorSffPairingPi2Inheritance>().XorSectorRate(N, gammaZero), precision: 14);
    }

    [Fact]
    public void RegisterF43_PartnerSumEqualsN_HoldsAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF43SectorSffPairingPi2Inheritance()
            .Build();

        var f43 = registry.Get<F43SectorSffPairingPi2Inheritance>();
        for (int N = 2; N <= 10; N++)
            for (int w = 0; w <= N; w++)
                Assert.True(f43.PartnerSumEqualsN(w, N));
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
                // Missing: RegisterF1Pi2Inheritance
                .RegisterF43SectorSffPairingPi2Inheritance()
                .Build());
    }
}
