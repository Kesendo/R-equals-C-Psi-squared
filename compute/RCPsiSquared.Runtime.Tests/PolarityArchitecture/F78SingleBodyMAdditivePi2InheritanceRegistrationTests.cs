using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F78SingleBodyMAdditivePi2InheritanceRegistrationTests
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
    public void RegisterF78_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF78SingleBodyMAdditivePi2Inheritance()
            .Build();
        Assert.True(registry.Contains<F78SingleBodyMAdditivePi2Inheritance>());
    }

    [Fact]
    public void RegisterF78_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF78SingleBodyMAdditivePi2Inheritance()
            .Build();
        Assert.Equal(Tier.Tier1Derived, registry.Get<F78SingleBodyMAdditivePi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF78_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF78SingleBodyMAdditivePi2Inheritance()
            .Build();
        var ancestors = registry.AncestorsOf<F78SingleBodyMAdditivePi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F1Pi2Inheritance), ancestors);
        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
        Assert.Contains(typeof(Pi2I4MemoryLoopClaim), ancestors);
    }

    [Fact]
    public void RegisterF78_TwoIFactorAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF78SingleBodyMAdditivePi2Inheritance()
            .Build();
        var f78 = registry.Get<F78SingleBodyMAdditivePi2Inheritance>();
        Assert.Equal(0.0, f78.TwoIFactor.Real, precision: 14);
        Assert.Equal(2.0, f78.TwoIFactor.Imaginary, precision: 14);
    }

    [Fact]
    public void RegisterF78_YEqualsZSvdSpectrumAcrossRegistry()
    {
        var registry = BuildBaseRegistry()
            .RegisterF78SingleBodyMAdditivePi2Inheritance()
            .Build();
        Assert.True(registry.Get<F78SingleBodyMAdditivePi2Inheritance>()
            .YEqualsZSvdSpectrum(cl: 1.0, gammaZero: 0.05));
    }

    [Fact]
    public void RegisterF78_WithoutF1_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterPi2I4MemoryLoop()
                .RegisterF78SingleBodyMAdditivePi2Inheritance()
                .Build());
    }
}
