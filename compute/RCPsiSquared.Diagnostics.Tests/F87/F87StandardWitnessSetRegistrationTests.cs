using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class F87StandardWitnessSetRegistrationTests
{
    private static ChainSystem DefaultChain(int N = 3) =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.Heisenberg, Topology: TopologyKind.Chain);

    [Fact]
    public void RegisterF87StandardWitnessSet_AddsClaim()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterF87Family()
            .RegisterF87StandardWitnessSet()
            .Build();

        Assert.True(registry.Contains<F87StandardWitnessSet>());
    }

    [Fact]
    public void RegisterF87StandardWitnessSet_HasFiveCanonicalWitnesses()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterF87Family()
            .RegisterF87StandardWitnessSet()
            .Build();

        var set = registry.Get<F87StandardWitnessSet>();
        Assert.Equal(5, set.Witnesses.Count);
    }

    [Fact]
    public void RegisterF87StandardWitnessSet_TierIsTier2Empirical()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterF87Family()
            .RegisterF87StandardWitnessSet()
            .Build();

        Assert.Equal(Tier.Tier2Empirical, registry.Get<F87StandardWitnessSet>().Tier);
    }

    [Fact]
    public void RegisterF87StandardWitnessSet_AncestorsContainsTrichotomyClassification_AndChain()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterF87Family()
            .RegisterF87StandardWitnessSet()
            .Build();

        var ancestors = registry.AncestorsOf<F87StandardWitnessSet>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
        Assert.Contains(typeof(ChainSystemPrimitive), ancestors);
    }

    [Theory]
    [InlineData(0, "XX+YY", TrichotomyClass.Truly)]
    [InlineData(1, "XX+YY+ZZ (Heisenberg)", TrichotomyClass.Truly)]
    [InlineData(2, "YZ+ZY (EQ-030 soft)", TrichotomyClass.Soft)]
    [InlineData(3, "XX+XY (mixed hard)", TrichotomyClass.Hard)]
    [InlineData(4, "XY+YX (bond-flip soft)", TrichotomyClass.Soft)]
    public void RegisterF87StandardWitnessSet_LazyClassification_MatchesExpected(
        int index, string expectedNameSuffix, TrichotomyClass expectedClass)
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .RegisterPi2Family()
            .RegisterF87Family()
            .RegisterF87StandardWitnessSet()
            .Build();

        var w = registry.Get<F87StandardWitnessSet>().Witnesses[index];
        Assert.Contains(expectedNameSuffix, w.Name);
        Assert.Equal(expectedClass, w.ExpectedClass);
        Assert.Equal(expectedClass, w.ActualClass);
        Assert.True(w.Matches, $"witness {w.Name} expected {w.ExpectedClass} got {w.ActualClass}");
    }

    [Fact]
    public void RegisterF87StandardWitnessSet_WithoutTrichotomyParent_Throws()
    {
        Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .RegisterF1Family(DefaultChain())
                .RegisterPi2Family()
                // Missing: RegisterF87Family
                .RegisterF87StandardWitnessSet()
                .Build());
    }
}
