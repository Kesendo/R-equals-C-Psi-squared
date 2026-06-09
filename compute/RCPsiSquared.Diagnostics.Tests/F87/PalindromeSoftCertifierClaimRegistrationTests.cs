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
/// against the spectral authority (<see cref="PauliPairTrichotomy"/>, now including the XX+XZ (Routing) and
/// XIX+XXY+YXX (RoutingKBody) cases) and the 2 non-local Z-middle k-body routed-soft ceiling cases
/// (XZX+XZY+YZX, YZY+XZY+YZX are each soft yet NotCertified, admitting no per-site product Q, so even the
/// derived k-body routing declines them; the two I-heavy IXI+IIY+YII, IYI+IIX+XII are now LOCAL, certified by
/// the SingleSiteField strategy, the 4 to 2 step). Built
/// from the full <see cref="KnowledgeRegistryFactory.BuildDefault"/> registry, so both typed
/// parents are present and the strength-inheritance check (5 ≥ 4 against F87DiagonalCellBipartiteWitnessSet,
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
    public void PalindromeSoftCertifierClaim_SoundnessBattery_IncludesTheRoutingKBodyCase()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();

        // Stufe B adds XIX+XXY+YXX as a k-body routing soundness case: a routable 3-body set the 2-body
        // family table misses, certified by the derived per-term k-site routing (the P4 pattern) and not hard.
        var routingKBody = claim.SoundnessBattery.Single(c => c.Name == "XIX+XXY+YXX (RoutingKBody)");
        Assert.Equal(Strategy.RoutingKBody, routingKBody.Strategy);
        Assert.True(routingKBody.Certified, "XIX+XXY+YXX must be certified (by RoutingKBody)");
        Assert.True(routingKBody.NotHard, "XIX+XXY+YXX must be not hard by the spectral authority");
        Assert.True(routingKBody.Passes);
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_SoundnessBattery_IncludesTheSingleSiteFieldCase()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();

        // The 4 to 2 step adds IXI+IIY+YII as a single-site-field soundness case: a sum of single-site
        // transverse fields, certified by the SingleSiteField strategy (the per-site crossover product) and
        // not hard. An I-heavy case once on the ceiling, now local.
        var singleSite = claim.SoundnessBattery.Single(c => c.Name == "IXI+IIY+YII (SingleSiteField)");
        Assert.Equal(Strategy.SingleSiteField, singleSite.Strategy);
        Assert.True(singleSite.Certified, "IXI+IIY+YII must be certified (by SingleSiteField)");
        Assert.True(singleSite.NotHard, "IXI+IIY+YII must be not hard by the spectral authority");
        Assert.True(singleSite.Passes);
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_Ceiling_KBodyRoutedSoft_Soft_NotCertified()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();

        var ceiling = claim.Ceiling;
        Assert.Equal(2, ceiling.Count);
        var names = ceiling.Select(c => c.Name).ToList();
        Assert.Contains("XZX+XZY+YZX", names);
        Assert.Contains("YZY+XZY+YZX", names);
        Assert.DoesNotContain("IXI+IIY+YII", names);   // I-heavy now local (SingleSiteField), left the ceiling (4 to 2)
        Assert.DoesNotContain("IYI+IIX+XII", names);
        Assert.All(ceiling, c => Assert.True(c.IsSoft, $"ceiling witness {c.Name} must be soft"));
        Assert.All(ceiling, c => Assert.False(c.Certified, $"ceiling witness {c.Name} must be NotCertified (non-local)"));
        Assert.All(ceiling, c => Assert.True(c.Holds, $"ceiling pair {c.Name} (soft, NotCertified) must hold"));
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_SelfCheck_Passes()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();
        Assert.True(claim.SelfCheckPasses,
            "the full self-check (soundness battery all certified-and-not-hard; the 2 Z-middle k-body ceiling pairs) must pass");
    }
}
