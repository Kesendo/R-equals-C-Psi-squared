using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.F1Family;

/// <summary>Wiring tests for <see cref="CollisionGapOddOrdersClaim"/> (F161), whose single typed parent is
/// <see cref="CrackedRingExactCurveClaim"/> (F160, the polynomial that is the whole input to the series and, in its
/// Theorem G, that series at first order). Built on the default registry, which is also where the statement's earned
/// fences are pinned: F160's own parents are four claims deep, so hand-building the chain in Core.Tests would test the
/// chain rather than this claim.</summary>
public class CollisionGapOddOrdersClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_HoldsTheClaim_Tier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<CollisionGapOddOrdersClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<CollisionGapOddOrdersClaim>().Tier);
    }

    [Fact]
    public void F160_IsTheTypedParent_AndF129sTypedHomeIsDeliberatelyNotOne()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<CollisionGapOddOrdersClaim>();

        var ancestors = registry.AncestorsOf<CollisionGapOddOrdersClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(CrackedRingExactCurveClaim), ancestors);
        Assert.Same(registry.Get<CrackedRingExactCurveClaim>(), claim.CrackedRing);

        // F2b and F65 are ancestors THROUGH F160 and are not re-added as edges here.
        Assert.Contains(typeof(F2bXyChainSpectrumPi2Inheritance), ancestors);

        // F129's typed home is Tier1Candidate, so a typed edge to it would break the tier rule for a Tier1Derived
        // child. It is carried by anchor instead, and executably by LevelCollisionCensus / CollisionFamilyInventory,
        // which the claim calls. If that claim is ever promoted, this assertion is the place the choice is revisited.
        Assert.Equal(Tier.Tier1Candidate, registry.Get<CrossTripleOrthogonalityClaim>().Tier);
        Assert.DoesNotContain(typeof(CrossTripleOrthogonalityClaim), ancestors);
    }

    [Fact]
    public void TheStatement_KeepsTheFencesThatWereEarned()
    {
        // Each of these was paid for. The letter n is the COMB modulus, not the site count (the crack code spends the
        // same letter the other way). Theorem E has ONE forced direction: n = 24 breaks the converse, and the first
        // draft read it as an equivalence off n = 30, a modulus that could not break it. The second-order count is a
        // LOWER BOUND past the census, the converse being checked at n <= 30 only. The ROT3 shape of the twelve is a
        // census observation, not a derivation. And only the PARITY of the multipliers is a theorem, not their range.
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<CollisionGapOddOrdersClaim>();
        Assert.Contains("COMB modulus", claim.Name);
        Assert.Contains("ONE forced direction", claim.Name);
        Assert.Contains("n = 24 breaks the converse", claim.Name);
        Assert.Contains("LOWER BOUND", claim.Name);
        Assert.Contains("READ and not derived", claim.Name);
        Assert.Contains("RANGE is not proved", claim.Name);
        Assert.Contains("docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md", claim.Anchor);
        Assert.Contains("simulations/collision_gap_odd_orders.py", claim.Anchor);
    }
}
