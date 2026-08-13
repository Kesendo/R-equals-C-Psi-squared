using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="BipartiteGaugeCriterionClaim"/> (F151): the bipartite gauge
/// 𝒟 = diag((−1)^(σ(a)+σ(b))) conjugates every coherence block to L(−q), and commutes with the site reflection
/// exactly when (p+q)(N−1) is even, which is what makes a block's R-sectors real matrix families on that parity.
/// Two typed parent edges, and the ORDER is the point:
///
/// <list type="bullet">
///   <item><see cref="ChiralKClaim"/>: the sublattice gauge one level down, K·H·K = −H on the single-excitation
///         site basis, proved in <c>PROOF_K_PARTNERSHIP</c>. 𝒟 is its Liouville-space lift, and every committed
///         claim that uses the sublattice grading already parents here (<c>DeadSetLawClaim</c>,
///         <c>MirrorOrderSortingClaim</c>, <c>SeedExistenceCountingClaim</c>). Parenting this claim on the FOLD
///         alone would have contradicted its own thesis, which is that gauge and fold are independent levers.</item>
///   <item><see cref="F89CrossFoldSimilarityClaim"/>: the fold lattice, where 𝒟 appears as ingredient (iv) and
///         whose §7 grading clause records that P, Q, 𝒟, T commute with R "up to a block-global sign" without
///         saying when the sign is −1. That is where the criterion is USED and the clause it completes.</item>
/// </list>
///
/// <para>Sharpening its first parent while descending from it: <see cref="ChiralKClaim"/> records that
/// "Z-dephasing trivially commutes with K, so K is a super-operator symmetry of the full Liouvillian". At the
/// block level the conjugation sends L(q) to L(−q), so it is an ANTI-symmetry of the pencil and a symmetry only
/// at q = 0; that is exactly the distinction this claim's criterion is built on.</para>
///
/// <para>Tier consistency: Tier 1 derived. The gauge identity is an exact entry-wise sign rearrangement (residual
/// 0.0 bit for bit) and the criterion follows in one line from σ(rev(c)) = w·(N−1) − σ(c). The parent is Tier 1
/// derived. What is NOT derived, and is scoped as such in the claim, is the converse: the criterion is sufficient
/// for conjugation-closure, and the non-closure on the odd parity is measured rather than proved.</para>
///
/// <para>Requires <see cref="F89CrossFoldSimilarityClaim"/> to be registered (the builder topo-resolves, so the
/// order of registration is free).</para></summary>
public static class BipartiteGaugeCriterionClaimRegistration
{
    public static ClaimRegistryBuilder RegisterBipartiteGaugeCriterionClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<BipartiteGaugeCriterionClaim>(b =>
            new BipartiteGaugeCriterionClaim(
                b.Get<ChiralKClaim>(),
                b.Get<F89CrossFoldSimilarityClaim>()));
}
