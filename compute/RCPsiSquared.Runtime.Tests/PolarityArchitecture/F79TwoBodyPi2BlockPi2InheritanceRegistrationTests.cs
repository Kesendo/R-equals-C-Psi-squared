using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F79TwoBodyPi2BlockPi2InheritanceRegistrationTests
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
            .RegisterF1Pi2Inheritance()
            .RegisterF88PopcountCoherence()
            .RegisterF88StaticDyadicAnchor()
            .RegisterPi2OperatorSpaceMirror();

    [Fact]
    public void RegisterF79_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterF79TwoBodyPi2BlockPi2Inheritance()
            .Build();
        Assert.True(registry.Contains<F79TwoBodyPi2BlockPi2Inheritance>());
    }

    [Fact]
    public void RegisterF79_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterF79TwoBodyPi2BlockPi2Inheritance()
            .Build();
        Assert.Equal(Tier.Tier1Derived, registry.Get<F79TwoBodyPi2BlockPi2Inheritance>().Tier);
    }

    [Fact]
    public void RegisterF79_AncestorsContainAllThreeParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterF79TwoBodyPi2BlockPi2Inheritance()
            .Build();
        var ancestors = registry.AncestorsOf<F79TwoBodyPi2BlockPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(KleinFourCellClaim), ancestors);
        Assert.Contains(typeof(Pi2OperatorSpaceMirrorClaim), ancestors);
        Assert.Contains(typeof(F1Pi2Inheritance), ancestors);
    }

    [Theory]
    [InlineData('X', 'Y', 1)]
    [InlineData('Y', 'Z', 0)]
    public void RegisterF79_Pi2ParityAcrossRegistry(char p, char q, int expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF79TwoBodyPi2BlockPi2Inheritance()
            .Build();
        Assert.Equal(expected, registry.Get<F79TwoBodyPi2BlockPi2Inheritance>().Pi2Parity(p, q));
    }

    [Theory]
    [InlineData(3, 32.0)]
    [InlineData(5, 512.0)]
    public void RegisterF79_Pi2BlockDimensionAcrossRegistry(int N, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterF79TwoBodyPi2BlockPi2Inheritance()
            .Build();
        Assert.Equal(expected, registry.Get<F79TwoBodyPi2BlockPi2Inheritance>().Pi2BlockDimension(N), precision: 12);
    }

    [Fact]
    public void RegisterF79_WithoutF1_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                .RegisterPi2DyadicLadder()
                .RegisterPi2I4MemoryLoop()
                .RegisterF88PopcountCoherence()
                .RegisterF88StaticDyadicAnchor()
                .RegisterPi2OperatorSpaceMirror()
                .RegisterF79TwoBodyPi2BlockPi2Inheritance()
                .Build());
    }
}
