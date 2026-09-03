using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.F1Family;

/// <summary>Wiring tests for <see cref="BlindSeatSectorFactorisationClaim"/> (F162), whose two typed parents are
/// <see cref="SeatCutBlindnessClaim"/> (F157, which owns the locus polynomial this one factorises and states no
/// sign and no leading coefficient) and <see cref="CrackedRingExactCurveClaim"/> (F160, which owns the Cassini
/// identity the congruence turns on). Built on the default registry: F160's own parents run four claims deep, so
/// hand-building the chain here would test the chain rather than this claim.</summary>
public class BlindSeatSectorFactorisationClaimRegistrationTests
{
    [Fact]
    public void DefaultRegistry_HoldsTheClaim_Tier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<BlindSeatSectorFactorisationClaim>());
        Assert.Equal(Tier.Tier1Derived, registry.Get<BlindSeatSectorFactorisationClaim>().Tier);
    }

    [Fact]
    public void BothParents_AreTypedEdges_AndTheirTiersAreWhatMakeTheChildDerived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<BlindSeatSectorFactorisationClaim>();

        var ancestors = registry.AncestorsOf<BlindSeatSectorFactorisationClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(SeatCutBlindnessClaim), ancestors);
        Assert.Contains(typeof(CrackedRingExactCurveClaim), ancestors);
        Assert.Same(registry.Get<SeatCutBlindnessClaim>(), claim.BlindSeat);
        Assert.Same(registry.Get<CrackedRingExactCurveClaim>(), claim.CrackedRing);

        // The tier rule caps a child by its weakest typed parent. Both of these are Tier1Derived, which is what
        // lets this claim be one; if either is ever downgraded, the registry throws at wire time and this
        // assertion is where the choice is revisited.
        Assert.Equal(Tier.Tier1Derived, registry.Get<SeatCutBlindnessClaim>().Tier);
        Assert.Equal(Tier.Tier1Derived, registry.Get<CrackedRingExactCurveClaim>().Tier);

        // F2b is an ancestor THROUGH F160 and is not re-added as an edge here.
        Assert.Contains(typeof(F2bXyChainSpectrumPi2Inheritance), ancestors);
    }

    [Fact]
    public void TheStatement_KeepsTheFencesThatWereEarned()
    {
        // Each of these was paid for over six review rounds. A sign law is a statement about a NAMED argument
        // order and about nothing else, because a resultant is antisymmetric up to (-1)^(deg f * deg g) and the
        // routine the gate first read it through does not keep the order it is given. The exponent reads the fold
        // coordinate and NOT N, which is what the seat-indexed reading hid. It is one factor per NON-POLE root and
        // never "every factor simple": no reading below N = 10 repeats, so the narrower wording would have been
        // unfalsifiable in the range first swept. And squarefreeness is a hypothesis, not a convenience.
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<BlindSeatSectorFactorisationClaim>();
        Assert.Contains("FOLD HALF FIRST", claim.Name);
        Assert.Contains("NOT N", claim.Name);
        Assert.Contains("NON-POLE", claim.Name);
        Assert.Contains("NAMED argument order", claim.Name);
        Assert.Contains("Squarefreeness", claim.Name);
        Assert.Contains("docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md", claim.Anchor);
        Assert.Contains("docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md", claim.Anchor);
        Assert.Contains("simulations/blind_seat_two_axes_proof.py", claim.Anchor);
    }
}
