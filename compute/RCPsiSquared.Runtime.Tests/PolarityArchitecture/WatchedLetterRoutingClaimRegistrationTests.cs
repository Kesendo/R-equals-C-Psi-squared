using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Registration gates for <see cref="WatchedLetterRoutingClaim"/> (the label layer,
/// typed): registered, Tier 1 derived, and carrying its two typed parent edges to
/// <see cref="AbsorptionTheoremClaim"/> (the −2γ·⟨n_XY⟩ price list) and
/// <see cref="Pi2KleinV4DephaseSwapGroup"/> (the letter swap).</summary>
public class WatchedLetterRoutingClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterPi2KleinV4DephaseSwapGroup();

    [Fact]
    public void RegisterWatchedLetterRouting_AddsClaim_Tier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterWatchedLetterRoutingClaim()
            .Build();
        Assert.True(registry.Contains<WatchedLetterRoutingClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<WatchedLetterRoutingClaim>().Tier);
    }

    [Fact]
    public void RegisterWatchedLetterRouting_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry()
            .RegisterWatchedLetterRoutingClaim()
            .Build();
        var ancestors = registry.AncestorsOf<WatchedLetterRoutingClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(Pi2KleinV4DephaseSwapGroup), ancestors);
    }

    [Fact]
    public void BuildDefault_ContainsWatchedLetterRoutingClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<WatchedLetterRoutingClaim>());
    }
}
