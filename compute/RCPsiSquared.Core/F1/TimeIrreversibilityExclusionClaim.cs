using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F1;

/// <summary>The typed ARGUMENT of TIME_IRREVERSIBILITY_EXCLUSION (2026-07-02). The VALUE
/// ‖{L_H, L_Dc}‖² = 4γ²(N−2)‖L_H‖² is already carried by the parent
/// <see cref="F49NonUniformCrossTermClaim"/>; this claim types the N=2-only orthogonality READING
/// that value implies, together with its load-bearing caveat.
///
/// <para>The anticommutator {L_H, L_Dc} of the Hamiltonian Liouvillian and the F1-centered
/// Z-dephasing Liouvillian vanishes EXACTLY at N=2 (every nonzero L_H entry connects Pauli strings
/// whose w_XY values sum to N, so at N=2 the centering cancels the anticommutator), giving the
/// Pythagorean split L_c² = L_H² + L_Dc²; it grows for N>2 with the γ- and topology-independent
/// relative R(N) = √((N−2)/(N·4^(N−1))).</para>
///
/// <para>The load-bearing caveat: the N=2 vanishing is a Frobenius-ORTHOGONALITY fact, NOT a
/// separability or reversibility criterion. Whether the flow factors is governed by the COMMUTATOR
/// [L_H, L_Dc], and ‖[L_H, L_Dc]‖ ≈ 22.6 ≠ 0 already at N=2. So the naive arrow-of-time reading a
/// vanishing anticommutator would suggest is EXCLUDED. Live witness (recomputes both norms across N
/// and cross-checks the anticommutator against the F49 closed form):
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/TimeIrreversibilityExclusionWitness.cs</c>,
/// <c>inspect --root time-exclusion</c>.</para></summary>
public sealed class TimeIrreversibilityExclusionClaim : Claim
{
    /// <summary>The typed parent: the F49 cross-term closed form whose (N−2) factor makes the
    /// anticommutator vanish at N=2 and grow afterwards.</summary>
    public F49NonUniformCrossTermClaim CrossTerm { get; }

    public TimeIrreversibilityExclusionClaim(F49NonUniformCrossTermClaim crossTerm)
        : base("time-irreversibility exclusion: {L_H, L_Dc} = 0 only at N=2 (orthogonality, not reversibility); [L_H, L_Dc] ≠ 0 there",
               Tier.Tier1Derived,
               "docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md")
    {
        CrossTerm = crossTerm ?? throw new ArgumentNullException(nameof(crossTerm));
    }

    /// <summary>The proof's relative cross-term R(N) = √((N−2)/(N·4^(N−1))): exactly 0 at N=2,
    /// then ≈ 1.83% (N=3), 2.07% (N=4); γ- and topology-independent.</summary>
    public static double RelativeCrossTerm(int n) => Math.Sqrt((n - 2.0) / (n * Math.Pow(4.0, n - 1)));

    public override string DisplayName =>
        "Time-irreversibility exclusion: {L_H, L_Dc} vanishes only at N=2 (orthogonality); the commutator does not";

    public override string Summary =>
        "the anticommutator {L_H, L_Dc} = 0 exactly at N=2 (Frobenius-orthogonality, the Pythagorean split " +
        "L_c² = L_H² + L_Dc²), growing for N>2 as R(N) = √((N−2)/(N·4^(N−1))) (γ- and topology-independent); but " +
        "the commutator [L_H, L_Dc] ≠ 0 already at N=2, so the vanishing is orthogonality, NOT a " +
        "separability/reversibility criterion, and the naive arrow-of-time reading is EXCLUDED. The value is the " +
        "parent F49 4γ²(N−2)‖L_H‖²; this node types the reading + caveat (live: TimeIrreversibilityExclusionWitness).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return CrossTerm;
            yield return new InspectableNode("the orthogonality (N=2 only)",
                summary: "every nonzero L_H entry connects Pauli strings whose w_XY values sum to N; at N=2 that sum " +
                         "is 2 = N for all 24 entries, so {L_H, L_D}_{ab} = (d_a + d_b)(L_H)_{ab} = −2Σγ(L_H)_{ab}, and " +
                         "adding 2Σγ·L_H (the F1 centering) gives {L_H, L_Dc} = 0. The Pythagorean split " +
                         "L_c² = L_H² + L_Dc² then holds exactly (‖·‖ = 0 verified at N=2).");
            yield return new InspectableNode("the break (N>2)",
                summary: "for N>2 the w_XY-sum rule fails for some entries; ‖{L_H, L_Dc}‖² = 4γ²(N−2)‖L_H‖² (F49), so " +
                         "the relative R(N) = √((N−2)/(N·4^(N−1))) is 0 at N=2, ≈ 1.83% (N=3), 2.07% (N=4), " +
                         "γ- and topology-independent (chain = ring = star = complete).");
            yield return new InspectableNode("the caveat: not reversibility",
                summary: "the N=2 vanishing is a Frobenius-ORTHOGONALITY fact, NOT a separability/reversibility " +
                         "criterion. Separability of the flow is governed by the COMMUTATOR [L_H, L_Dc], and " +
                         "‖[L_H, L_Dc]‖ ≈ 22.6 ≠ 0 already at N=2; so no arrow of time is read off the anticommutator, " +
                         "and the naive irreversibility claim is excluded (Tier-3 caveat of the proof).");
            yield return new InspectableNode("live witness",
                summary: "TimeIrreversibilityExclusionWitness (inspect --root time-exclusion) builds L_H and L_Dc " +
                         "across N via LindbladianBuilder, reports ‖{L_H, L_Dc}‖ (→ 0 at N=2) AND ‖[L_H, L_Dc]‖ " +
                         "(≠ 0 at N=2), and cross-checks the anticommutator against the F49 closed form fed a live " +
                         "per-bond norm (two independent computations meeting).");
        }
    }
}
