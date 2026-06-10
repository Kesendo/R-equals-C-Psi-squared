using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;
using Strategy = RCPsiSquared.Diagnostics.F87.PalindromeSoftCertifier.SoftStrategy;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="PalindromeSoftCertifierClaim"/>: it surfaces the §7.12
/// Liouvillian-free soft-certifier (sound, one-sided, the structural ceiling closed at zero) as a
/// registered, registry-queryable Claim. The Claim asserts ONLY the settled facts: the certifier's
/// soundness against the spectral authority (<see cref="PauliPairTrichotomy"/>, including the XX+XZ
/// (Routing), XIX+XXY+YXX (RoutingKBody), IXI+IIY+YII (SingleSiteField), and the two Z-middle
/// XZX+XZY+YZX, YZY+XZY+YZX (RoutingWindowSummed, the golden period-4 router, Stufe B′, F116) cases) and
/// the EMPTY recomputed ceiling (no soft-yet-NotCertified member remains in the k=3 windowed family; the
/// arc closed 6 → 4 → 2 → 0, docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md). Built
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
    public void PalindromeSoftCertifierClaim_SoundnessBattery_IncludesTheRoutingWindowSummedCases()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();

        // The 2 to 0 step (F116): the former Z-middle ceiling pair joins the battery as certified-positive
        // witnesses, certified by the window-summed golden period-4 router (Stufe B′) and not hard.
        foreach (var name in new[] { "XZX+XZY+YZX (RoutingWindowSummed)", "YZY+XZY+YZX (RoutingWindowSummed)" })
        {
            var c = claim.SoundnessBattery.Single(x => x.Name == name);
            Assert.Equal(Strategy.RoutingWindowSummed, c.Strategy);
            Assert.True(c.Certified, $"{name} must be certified (by RoutingWindowSummed)");
            Assert.True(c.NotHard, $"{name} must be not hard by the spectral authority");
            Assert.True(c.Passes);
        }
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_Ceiling_IsClosedAtZero()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();

        // F116 closes the arc 6 → 4 → 2 → 0: the recomputed ceiling (soft-yet-NotCertified members of the
        // k=3 windowed family) is EMPTY. The former witnesses XZX+XZY+YZX, YZY+XZY+YZX did not vanish from
        // the record: they moved to the soundness battery (RoutingWindowSummed), asserted above.
        Assert.Empty(claim.Ceiling);
    }

    [Fact]
    public void PalindromeSoftCertifierClaim_SelfCheck_Passes()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeSoftCertifierClaim>();
        Assert.True(claim.SelfCheckPasses,
            "the full self-check (soundness battery all certified-and-not-hard, incl. the two Z-middle golden cases; the recomputed ceiling empty) must pass");
    }
}
