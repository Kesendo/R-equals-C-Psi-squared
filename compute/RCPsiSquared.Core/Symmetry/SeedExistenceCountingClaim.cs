using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The seed-existence counting theorem (Tier 1 derived;
/// <c>experiments/F89_SEED_EXISTENCE_REDUCTION.md</c>, all three pieces proved 2026-07-04): on
/// the (1,2) block pencil L(q) = A + q·C of the XY chain under uniform Z-dephasing, the
/// real-eigenvalue-count identity
///
/// <code>
///     r(0⁺) − r(∞) = N − 1      for every odd N
/// </code>
///
/// where r(0⁺) = nullity(P₋₂CP₋₂) + nullity(P₋₆CP₋₆) counts the modes that stay real past the
/// q = 0 lift-off and r(∞) = nullity(C) counts the asymptotically-real modes. Since N − 1 > 0,
/// the real count strictly drops at finite q > 0: a real↔complex transition, which at a simple
/// discriminant zero is a defective 2×2 Jordan exceptional point (the Kato lemma), the SEED the
/// codim-1 containment corollary of <c>PROOF_CODIM1_BY_ADDITIVITY</c> transports across the
/// diamond. Three counting lemmas carry the identity:
///
/// <para><b>(N2), the path count:</b> the −2 rung decomposes into N − 1 disjoint simple paths of
/// N vertices each (the q-hop blocked at the shared site, the swap glue); a signed tree gauges
/// unsigned, and a path P_N has a zero mode iff N is odd, so n₂ = (N−1)·[N odd]. The odd-N seat
/// that cannot be mirrored is the combinatorial engine of the whole seed.</para>
///
/// <para><b>(FF), the fusion-resonance count:</b> ker C = the intertwiners {ρ: H₂ρ = ρH₁}; in
/// the JW eigenmode basis C is diagonal on cubic monomials c†_a c†_b c_c with eigenvalue
/// −i(λ_a+λ_b−λ_c), λ_k = 2cos(kπ/(N+1)), so nullity(C) = ρ = #{({a,b},c): λ_a+λ_b = λ_c}.</para>
///
/// <para><b>(N1′), the ordering-sector theorem:</b> K₆ = −i·P₋₆CP₋₆ is three hard-core
/// particles (2 ket + 1 bra) on the OPEN chain with mutual exclusion (ket-onto-bra hops ARE the
/// coupling block, not K₆), so no passing: the bra's rank among the three is conserved and
/// K₆ = K_L ⊕ K_M ⊕ K_R, each of dimension C(N,3); the diagonal gauge (−1)^{z_bra} turns each
/// sector into exactly −H₃, giving the closed form spec(K₆) = 3 × {−(λ_a+λ_b+λ_c)} and
/// n₆ = 3·Z₃. The chiral pairing λ_{N+1−k} = −λ_k maps zero-sum triples 3-to-1 onto fusion
/// resonances; the degenerate case D = #{2λ_x+λ_y = 0, x ≠ y} vanishes for every N by a
/// cyclotomic-integrality norm bound (all Galois conjugates |σ(λ_x)| ≤ 1 with |Nm| ≥ 1 force
/// |λ_y| = 2, impossible). Hence n₆ = 3·Z₃ = ρ, both parities. The fusion resonances appear in
/// BOTH endpoints and cancel; the surplus is the path kernel alone.</para>
///
/// <para><b>Scope (kept honest):</b> this claim types the COUNTING identity. The seed-existence
/// CONCLUSION additionally needs the codim-2 β-exotic genericity item (a count-dropping
/// transition is defective unless the non-generic order-3 nilpotent-linear-term point), which
/// stays OPEN in the doc's Status. Two adversarial reviews (exact arithmetic in ℤ[t]/Φ₂ₘ with a
/// counterexample hunt to N = 200; a full-2^N spin rebuild with explicit JW strings) held every
/// step.</para>
///
/// <para><b>Typed parents.</b> <see cref="AbsorptionTheoremClaim"/> (the rate law behind the
/// dephasing diagonal A = −2·diag(n_diff), whose two values −2/−6 define the rungs being
/// counted) and <see cref="ChiralKClaim"/> (the chiral single-particle pairing
/// λ_{N+1−k} = −λ_k, the mirror that powers the triple↔resonance bijection and the spectral
/// inheritance). Verifier: <c>simulations/seed_existence_nullity_check.py</c> (gates F1/N2/FF/
/// N1P/SI, N = 3..13 both parities); live: <c>inspect --root seedcount</c>
/// (<c>SeedExistenceCountingWitness</c>).</para></summary>
public sealed class SeedExistenceCountingClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: the AT rate law behind A = −2·diag(n_diff).
    public AbsorptionTheoremClaim Absorption { get; }
    // Parent-edge marker for Schicht-1 wiring: the chiral pairing λ_{N+1−k} = −λ_k.
    public ChiralKClaim ChiralK { get; }

    public SeedExistenceCountingClaim(AbsorptionTheoremClaim absorption, ChiralKClaim chiralK)
        : base("The seed-existence counting theorem: on the (1,2) block pencil L(q) = A + qC of the XY chain " +
               "under uniform Z-dephasing, r(0+) - r(inf) = N - 1 for every odd N, via (N2) the -2 rung = " +
               "N-1 disjoint paths of N vertices (zero mode iff N odd), (FF) nullity(C) = the free-fermion " +
               "fusion-resonance count rho = #{lambda_a+lambda_b = lambda_c}, and (N1') the ordering-sector " +
               "theorem: K6 = three no-passing bra-rank sectors, each gauge-equivalent via (-1)^{z_bra} to " +
               "MINUS the 3-magnon block H3, so spec(K6) = 3 x {-(la+lb+lc)} and n6 = 3*Z3 = rho (D = 0 by " +
               "cyclotomic integrality), both parities; the resonances cancel between the endpoints and the " +
               "surplus is the odd-N path kernel alone; the drop forces a real-to-complex transition, i.e. a " +
               "defective seed at a simple discriminant zero (Kato); the seed-existence conclusion stays open " +
               "only at the codim-2 beta-exotic genericity",
               Tier.Tier1Derived,
               "experiments/F89_SEED_EXISTENCE_REDUCTION.md + " +
               "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md + " +
               "simulations/seed_existence_nullity_check.py")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        ChiralK = chiralK ?? throw new ArgumentNullException(nameof(chiralK));
    }

    public override string DisplayName =>
        "Seed-existence counting: r(0⁺) − r(∞) = N − 1 (odd N) on the (1,2) block, the census input as a theorem";

    public override string Summary =>
        "r(0⁺) − r(∞) = N − 1 for every odd N on the (1,2) block: (N2) the −2 rung = N−1 paths of N vertices " +
        "(zero mode iff N odd), (FF) nullity(C) = the fusion-resonance count ρ, (N1′) n₆ = 3·Z₃ = ρ (the " +
        "ordering-sector theorem, spec(K₆) = 3×{−(λ_a+λ_b+λ_c)}, D = 0 by cyclotomic integrality); the " +
        "resonances cancel, the surplus forces a defective seed at a simple discriminant zero (Kato); open " +
        $"ink: the codim-2 β-exotic genericity ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("(N2) the path count: n₂ = (N−1)·[N odd]",
                summary: "the −2 rung (dim N(N−1)) decomposes into N−1 disjoint simple paths of N vertices " +
                         "(q-hop blocked at the shared site + swap glue); every path is a tree, the ±1 signs " +
                         "gauge away, and P_N has a zero mode iff N is odd: the unmirrorable middle seat as " +
                         "combinatorics.");
            yield return new InspectableNode("(FF) nullity(C) = ρ, the fusion resonances",
                summary: "ker C = {ρ: H₂ρ = ρH₁}; in the eigenmode basis C is diagonal on c†_a c†_b c_c with " +
                         "eigenvalue −i(λ_a+λ_b−λ_c), so nullity(C) = #{λ_a+λ_b = λ_c}: number-theoretic in " +
                         "the cosines (3, 6, 9, 12, 21, 18 at N = 3..13, the N = 11 jump from the " +
                         "cos75°+cos45°−cos15° = 0 identity).");
            yield return new InspectableNode("(N1′) the ordering-sector theorem: n₆ = 3·Z₃ = ρ",
                summary: "no passing on the open chain conserves the bra's rank among the three hard-core " +
                         "particles: K₆ = K_L ⊕ K_M ⊕ K_R, each gauge-equivalent via diag((−1)^{z_bra}) to " +
                         "−H₃, so spec(K₆) = 3 × {−(λ_a+λ_b+λ_c)} (a closed form, not just a nullity); the " +
                         "chiral 3-to-1 bijection and the D = 0 cyclotomic norm bound turn 3·Z₃ into ρ. The " +
                         "'SUSY-like spectral inheritance' of the first landing is thereby explained: a " +
                         "per-block multiset theorem (chiral mirror + spectator pairs), NOT a partition of " +
                         "spec(C) (the images overlap).");
            yield return new InspectableNode("the scope boundary (what stays open)",
                summary: "the counting identity is the theorem; the seed-existence conclusion still needs the " +
                         "codim-2 β-exotic genericity item (order-3 nilpotent-linear-term points), OPEN in the " +
                         "doc's Status. The near-miss law min|2λ_x+λ_y| ≈ 2π³/(N+1)³ caps numerical D-scans " +
                         "(resonance tol 1e-7 safe to N ≈ 850, SVD floor 1e-9 to N ≈ 3959); only the exact " +
                         "argument makes D = 0 a theorem.");
            yield return new InspectableNode("live witness + verifier",
                summary: "inspect --root seedcount recomputes n₂/n₆/nullity(C) (SVD), ρ and Z₃ " +
                         "(combinatorial), the exact-zero cross-sector and gauge gates, and the two-sided " +
                         "nonzero controls at inspect time (odd N ≤ 9); " +
                         "simulations/seed_existence_nullity_check.py asserts the full gate set N = 3..13, " +
                         "both parities.");
        }
    }
}
