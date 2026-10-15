using System.Linq;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="GoldenRouterClaim"/>: it types the F116 golden/metallic
/// router (docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md, reviewed bit-exact 2026-06-22) as a standalone
/// Tier1Derived Claim, the "C# witness first" inverse gap (the live witness GoldenRouterWitness and the
/// reused helper <see cref="KBodyPalindromeRouting"/> already existed; only the tiered Claim wrapper was
/// missing). The Claim asserts ONLY the settled, bit-exact facts: the router existence + closed form
/// W L W⁻¹ = −L − 2σ at every N ≥ 3 for arbitrary site-dependent γ, certified live by the window-summed
/// router; and the metallic-family check r(c) live = closed form. Its two typed parents are
/// <see cref="F1PalindromeIdentity"/> (the global palindrome the router realizes locally for the
/// Z-middle ceiling class) and <see cref="WindowedConverseThresholdClaim"/> (the F87 two-reflection
/// chiral spine its two-sided form instantiates). Both Tier1Derived, so the strength-inheritance check
/// (5 ≥ 5 against each) is exercised in production.
///
/// <para>NOTE on the parent choice: the open-arc f116_golden_router_typed_claim originally SUGGESTED
/// parenting on PalindromeSoftCertifierClaim, but that is Tier1Candidate (4 &lt; 5) so it would violate
/// the parent ≥ child rule, AND it is backwards (the certifier USES this router as a helper, so this is
/// logically upstream of it). The certifier relationship is expressed as a see-cref, not a parent edge.</para></summary>
public class GoldenRouterClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsGoldenRouterClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<GoldenRouterClaim>());
    }

    [Fact]
    public void GoldenRouterClaim_TypedParents_ArePresentInRegistry()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F1PalindromeIdentity>(),
            "typed parent F1PalindromeIdentity (the palindrome form the router realizes) must be registered");
        Assert.True(registry.Contains<WindowedConverseThresholdClaim>(),
            "typed parent WindowedConverseThresholdClaim (the F87 two-reflection chiral spine) must be registered");
    }

    [Fact]
    public void GoldenRouterClaim_Ancestors_ContainBothTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<GoldenRouterClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(WindowedConverseThresholdClaim), ancestors);
    }

    [Fact]
    public void GoldenRouterClaim_TierIsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<GoldenRouterClaim>().Tier);
    }

    [Fact]
    public void GoldenRouterClaim_StrengthInheritance_BothParentsAtLeastAsStrong()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var child = registry.Get<GoldenRouterClaim>().Tier;
        Assert.True(TierStrength.IsAtLeastAsStrong(Tier.Tier1Derived, child),
            "F1PalindromeIdentity (Tier1Derived) must be at least as strong as the router claim");
        Assert.True(TierStrength.IsAtLeastAsStrong(
            registry.Get<WindowedConverseThresholdClaim>().Tier, child),
            "WindowedConverseThresholdClaim must be at least as strong as the router claim");
    }

    [Fact]
    public void GoldenRouterClaim_GoldenCase_CertifiedByWindowSummedRouter_PerTermDeclines()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<GoldenRouterClaim>();

        // The window-summed golden router certifies XZX+XZY+YZX (the constructive existence side), while
        // the per-term lens correctly declines it (cross-template cancellation is invisible to it). That
        // split IS the F116 result: a local period-4 palindromizer exists for the Z-middle ceiling case.
        Assert.Equal(KBodyPalindromeRouting.GoldenDescription, claim.GoldenRoutes);
        Assert.True(claim.GoldenPerTermDeclines,
            "the per-term lens must decline the golden case (the documented coverage gap, not non-locality)");
    }

    [Fact]
    public void GoldenRouterClaim_SiblingCase_CertifiedByMirrorRouter()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<GoldenRouterClaim>();

        // The X↔Y sibling YZY+XZY+YZX is certified by the mirror period-4 router (the X↔Y conjugate maps).
        Assert.Equal(KBodyPalindromeRouting.GoldenMirrorDescription, claim.SiblingRoutes);
    }

    [Fact]
    public void GoldenRouterClaim_MetallicMean_LiveMatchesClosedForm_AtGoldenAndSilver()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<GoldenRouterClaim>();

        // The genuine check: the live router-derived frame ratio (root of ‖{W,S}‖_F = 0) equals the
        // closed-form metallic mean r(c) = (c+√(c²+4))/2. Golden c=1 → φ, silver c=2 → 1+√2.
        var golden = claim.MetallicReadings.Single(r => r.C == 1.0);
        Assert.True(golden.Matches, $"golden c=1: live r={golden.LiveRatio} vs closed φ={golden.ClosedForm}");

        var silver = claim.MetallicReadings.Single(r => r.C == 2.0);
        Assert.True(silver.Matches, $"silver c=2: live r={silver.LiveRatio} vs closed 1+√2={silver.ClosedForm}");
    }

    [Fact]
    public void GoldenRouterClaim_LetterTable_EigenEquationIsLocus_UniformChiralExactlyOnIt()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<GoldenRouterClaim>();

        // The golden direction is an eigenvector of the window's letter table M_c = [[c,1],[1,0]], and the
        // uniform product â^⊗N is chiral exactly there (integer polynomial arithmetic, no rounding).
        Assert.Equal(new long[] { 0, 1, 2, 3 }, claim.LetterTableReadings.Select(t => t.C));
        foreach (var t in claim.LetterTableReadings)
        {
            Assert.True(t.EigenEquationIsLocus, $"c={t.C}: det[v, M_c v] must equal r² − cr − 1");
            Assert.True(t.ChiralResidualIsLocusTimesK, $"c={t.C}: the chiral residual must be (r²−cr−1)·K(r)");
        }
        Assert.True(claim.ReflectedTableLeavesRemainder,
            "the hard reflected table [[1,1],[−1,0]] must leave a remainder: the divisibility test fires off the line");
    }

    [Fact]
    public void GoldenRouterClaim_Rhythm_IsTheOneZPerWindowPattern_FourShifts()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<GoldenRouterClaim>();

        // The [a,a,b,b] rhythm is a^⊗N times the Z pattern [I,I,Z,Z]; one Z per window end pair admits
        // exactly its four cyclic shifts, and no uniform or period-2 pattern.
        Assert.Equal(new[] { "IIZZ", "IZZI", "ZIIZ", "ZZII" }, claim.RhythmSolutions);
    }

    [Fact]
    public void GoldenRouterClaim_TwoSided_GoldenChiralOnOneSide_ZPatternOnTheOther()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<GoldenRouterClaim>();

        // The dephaser needs opposite Z-parity on the two sides at every site, the Hamiltonian needs both
        // sides chiral: ρ ↦ a^⊗N ρ S routes L exactly (either side, open chain, ring at 4 | N), and equal
        // parity, a one-sided identity column, and the ring at N = 5 do not. The two groups are one gate: a
        // dropped commutator keeps the four positives true and is caught only by the controls.
        Assert.Equal(
            new[] { true, true, true, true, false, false, false },
            claim.SplitRouterReadings.Select(r => r.Routes));

        // On a ring the alternation a ↔ b closes along the parity cycles exactly when 4 | N.
        Assert.Equal(new[] { 0, 4, 0, 0, 0, 4, 0 }, claim.RingRhythmCounts);
    }

    [Fact]
    public void GoldenRouterClaim_SelfCheck_Passes()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<GoldenRouterClaim>();
        Assert.True(claim.SelfCheckPasses,
            "the full live self-check (window-summed golden + sibling routing certified, per-term declines, " +
            "metallic mean live = closed form at golden + silver, residual vanishes at the mean) must pass");
        Assert.Equal(claim.Battery.Count, claim.PassCount);
    }
}
