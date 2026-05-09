using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.F71Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F86F71MirrorPi2InheritanceRegistrationTests
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
            .RegisterF71Family(N: 5)
            .RegisterF1Pi2Inheritance()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror()
            .RegisterHalfIntegerMirror(5)
            .RegisterF71MirrorSymmetryPi2Inheritance();

    [Fact]
    public void RegisterF86F71Mirror_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86F71MirrorPi2Inheritance()
            .Build();
        Assert.True(registry.Contains<F86F71MirrorPi2Inheritance>());
    }

    [Fact]
    public void RegisterF86F71Mirror_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86F71MirrorPi2Inheritance()
            .Build();
        Assert.Equal(Tier.Tier1Derived, registry.Get<F86F71MirrorPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF86F71Mirror_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF86F71MirrorPi2Inheritance()
            .Build();
        var ancestors = registry.AncestorsOf<F86F71MirrorPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F71MirrorSymmetryPi2Inheritance), ancestors);
        Assert.Contains(typeof(F86MirrorGeneralisationLink), ancestors);
    }

    [Theory]
    [InlineData(5, 0, 3)]
    [InlineData(6, 2, 2)]
    public void RegisterF86F71Mirror_MirrorPartnerBondAcrossRegistry(int N, int b, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF86F71MirrorPi2Inheritance()
            .Build();
        Assert.Equal(expected, registry.Get<F86F71MirrorPi2Inheritance>().MirrorPartnerBond(N, b));
    }

    [Fact]
    public void RegisterF86F71Mirror_WithoutF71Family_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterPi2I4MemoryLoop()
                .RegisterF1Pi2Inheritance()
                .RegisterF88PopcountCoherence()
                .RegisterF88StaticDyadicAnchor()
                .RegisterPi2OperatorSpaceMirror()
                .RegisterHalfIntegerMirror(5)
                .RegisterF71MirrorSymmetryPi2Inheritance()
                // Missing: RegisterF71Family
                .RegisterF86F71MirrorPi2Inheritance()
                .Build());
    }
}
