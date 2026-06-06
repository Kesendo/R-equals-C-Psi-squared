using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;
using Strategy = RCPsiSquared.Diagnostics.F87.PalindromeSoftCertifier.SoftStrategy;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="PalindromeSoftCertifierClaim"/>: it surfaces the §7.12
/// Liouvillian-free soft-certifier (sound, one-sided, structurally incomplete) as a registered,
/// registry-queryable Claim. The Claim asserts ONLY the settled facts: the certifier's soundness
/// against the spectral authority (<see cref="PauliPairTrichotomy"/>, now including the XX+XZ routing case)
/// and the k-body routed-soft ceiling (XZX+XZY+YZX is soft yet NotCertified, beyond the 2-body routing
/// table). Built from the full <see cref="KnowledgeRegistryFactory.BuildDefault"/> registry, so both typed
/// parents are present and the strength-inheritance check (4 ≥ 4 against F87DiagonalCellBipartiteWitnessSet,
/// 5 ≥ 4 against F87TrichotomyClassification) is exercised in production.</summary>
public class PalindromeSoftCertifierClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsPalindromeSoftCertifierClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<PalindromeSoftCertifierClaim>());
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_TypedParents_ArePresentInRegistry()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F87DiagonalCellBipartiteWitnessSet>(),
            "typed parent F87DiagonalCellBipartiteWitnessSet (the §7 criterion the linear strategy scales) must be registered");
        Assert.True(registry.Contains<F87TrichotomyClassification>(),
            "typed parent F87TrichotomyClassification (the spectral authority soundness is checked against) must be registered");
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_Ancestors_ContainBothTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<PalindromeSoftCertifierClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F87DiagonalCellBipartiteWitnessSet), ancestors);
        Assert.Contains(typeof(F87TrichotomyClassification), ancestors);
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_TierIsTier1Candidate()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Candidate, registry.Get<PalindromeSoftCertifierClaim>().Tier);
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_SoundnessBattery_AllCertifiedAndNotHard()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();

        Assert.True(claim.SoundnessBattery.Count >= 2,
            "soundness battery should hold at least XY+YX (ExcitationPairing) and XZ+ZX (ExcitationParity)");
        foreach (var c in claim.SoundnessBattery)
        {
            Assert.True(c.Certified, $"soundness case '{c.Name}': certifier must certify it");
            Assert.True(c.NotHard, $"soundness case '{c.Name}': spectral authority must agree (not hard)");
            Assert.True(c.Passes, $"soundness case '{c.Name}' did not pass");
        }
        Assert.Equal(claim.SoundnessBattery.Count, claim.SoundnessPassCount);
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_SoundnessBattery_IncludesTheRoutingCase()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();

        // The receded ceiling adds XX+XZ as a routing soundness case: certified by the non-diagonal
        // hidden-Q routing (the colourings cannot reach its non-bipartite basis-state graph) and not hard.
        var routing = claim.SoundnessBattery.Single(c => c.Name == "XX+XZ (Routing)");
        Assert.Equal(Strategy.Routing, routing.Strategy);
        Assert.True(routing.Certified, "XX+XZ must be certified (by Routing)");
        Assert.True(routing.NotHard, "XX+XZ must be not hard by the spectral authority");
        Assert.True(routing.Passes);
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_Ceiling_KBodyRoutedSoft_Soft_NotCertified()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();

        var ceiling = claim.Ceiling;
        Assert.Equal("XZX+XZY+YZX", ceiling.Name);
        Assert.True(ceiling.IsSoft, "ceiling witness XZX+XZY+YZX must be soft by the spectral authority");
        Assert.False(ceiling.Certified,
            "ceiling witness XZX+XZY+YZX must be NotCertified (the 2-body-gated certifier cannot reach a 3-body routed case)");
        Assert.True(ceiling.Holds, "the ceiling pair (soft, NotCertified) must hold");
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_SelfCheck_Passes()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();
        Assert.True(claim.SelfCheckPasses,
            "the full self-check (soundness battery all certified-and-not-hard; the XX+XZ ceiling triple) must pass");
    }
}
