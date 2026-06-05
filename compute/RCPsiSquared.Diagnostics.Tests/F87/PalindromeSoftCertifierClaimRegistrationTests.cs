using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="PalindromeSoftCertifierClaim"/>: it surfaces the §7.12
/// Liouvillian-free soft-certifier (sound, one-sided, structurally incomplete) as a registered,
/// registry-queryable Claim. The Claim asserts ONLY the settled facts: the certifier's soundness
/// against the spectral authority (<see cref="PauliPairTrichotomy"/>) and the proven non-bipartite-soft
/// ceiling (XX+XZ is soft, non-bipartite, NotCertified). Built from the full
/// <see cref="KnowledgeRegistryFactory.BuildDefault"/> registry, so both typed parents are present and
/// the strength-inheritance check (4 ≥ 4 against F87DiagonalCellBipartiteWitnessSet, 5 ≥ 4 against
/// F87TrichotomyClassification) is exercised in production.</summary>
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
    public void PalindromeSoftCertifierClaim_Ceiling_XXxz_Soft_NonBipartite_NotCertified()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();

        var ceiling = claim.Ceiling;
        Assert.True(ceiling.IsSoft, "ceiling witness XX+XZ must be soft by the spectral authority");
        Assert.False(ceiling.IsBipartite, "ceiling witness XX+XZ must have a non-bipartite basis-state graph");
        Assert.False(ceiling.Certified, "ceiling witness XX+XZ must be NotCertified (the proven incompleteness)");
        Assert.True(ceiling.Holds, "the ceiling triple (soft, non-bipartite, NotCertified) must hold");
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
