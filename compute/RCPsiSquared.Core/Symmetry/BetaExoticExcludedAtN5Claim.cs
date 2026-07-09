using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The β-exotic exclusion at N = 5 (Tier 1 derived;
/// <c>experiments/F89_SEED_EXISTENCE_REDUCTION.md</c>, section "The β-exotic is excluded at
/// N = 5 only", landed 2026-07-09 after three empty reviews): on the (1,2) coherence pencil
/// L(q) = A + q·C of the 5-site XY chain under uniform Z-dephasing, no branch locus q\* ≠ 0 of the
/// residual charpoly factor F_res carries the Puiseux exponent 3/2 of the β-exotic. Both R-parity
/// sectors, exactly, over ℚ(i).
///
/// <para><b>Why a multiplicity settles it.</b> At a branch locus, the order of vanishing of
/// disc_Λ(F_res) reads the Puiseux exponent of the branches meeting there: a defective EP2
/// (exponent ½) gives 1, a diabolic crossing gives 2, a cubic branch point (3×3 Jordan, exponent ⅓)
/// gives 2, and the β-exotic (normal form β(s) = [[0, s], [s², 0]], eigenvalues ±s^{3/2}) gives 3.
/// Every OTHER colliding pair at the same locus contributes a NON-NEGATIVE order, so a coincidental
/// extra collision can only raise the total (3 + 1, 3 + 2), never mask it. Hence "disc_Λ(F_res) has
/// no root of multiplicity ≥ 3 off q = 0" excludes the β-exotic outright. This is deliberately
/// weaker than squarefreeness: the diabolic loci are genuine DOUBLE roots and must survive, which is
/// exactly why the dead discriminant-Galois route is unnecessary here. A Galois group sees which
/// sheets swap; it never sees the ½-versus-3/2 exponent. A multiplicity does.</para>
///
/// <para><b>The certificate, and why one prime suffices.</b>
/// <see cref="F89PathK.FoldResultantCertificate"/> interpolates D(q) = disc_Λ(F_res)(q) mod p and
/// squarefree-splits it; at N = 5 the certified reading is [56, 26] (R-odd) and [56, 32] (R-even):
/// two layers, maximum root multiplicity 2. The lift is ONE-WAY: reduction modulo a prime ideal
/// π_p ∤ lc_q(D) is a degree-preserving ring homomorphism, so a factor (q − α)^m survives and
/// distinct roots can only MERGE; therefore max-mult(D mod p) ≥ max-mult(D), and reading
/// multiplicity 2 at a single CERTIFIED prime proves max-mult(D) ≤ 2 exactly. The prime must be good
/// at both ends of the q-axis: it attains the true deg_q D (no root escaped to q = ∞) and the true
/// q-valuation (no nonzero root collapsed onto q = 0, where the q-power strip would discard it). A
/// Hadamard bound on ‖D‖ certifies EACH of those two true values separately (one bound per statement,
/// not one for their union); a prime attaining BOTH is then searched for among the primes the run
/// samples anyway, and the certificate fails closed if none does.</para>
///
/// <para><b>H1 follows as a corollary at N = 5</b> (each forced seed has algebraic multiplicity exactly
/// 2, no 3×3 Jordan), by ELIMINATION, not by counting the order of the zero. The tempting shortcut
/// ("a count-drop changes the discriminant's sign, so the zero has odd order, so with max multiplicity
/// 2 the order is 1") silently assumes exactly one conjugate pair is born at q\*; two forced seeds
/// coinciding at one q\* would drop the real count by 4 at an order-2, sign-preserving zero. Do not use
/// it. What holds without any order parity: maximum multiplicity 2 already excludes the β-exotic
/// (ord 3), every Jordan block of size ≥ 4 (ord ≥ 3), and every semisimple degeneracy of three or more
/// branches (ord ≥ 3). What remains at a branch locus is a defective EP2 (ord 1), a diabolic crossing
/// (ord 2), an analytic-defective 2×2 (ord 2, e.g. [[0,1],[0,s]]), a cubic branch point (ord 2), or a
/// coincidence of these; and of that list ONLY the EP2 changes the real-eigenvalue count (the analytic
/// cases keep two real branches, a cubic keeps one real branch and one conjugate pair on both sides).
/// The bound also forbids an EP2 from SHARING its locus with any ord-2 structure, since 1 + 2 = 3. So a
/// locus carrying a count-drop carries EP2s and nothing else (one, or two coincident), each a defective
/// 2×2 Jordan block.</para>
///
/// <para>That step needs F_res to have REAL coefficients. At odd N the bipartite sign
/// T = diag((−1)^{a₀+a₁+b₀}) satisfies T K T = −K exactly (every entry of K is a single hop, which
/// flips the site-sum parity) and T A T = A, so T L(q) T = L(q)† is a similarity; T also commutes with
/// the reflection R (3(N+1) is even), so it restricts to each R-parity sector and each sector's
/// spectrum is self-conjugate. The AT spectrum inside the sector is separately self-conjugate, its
/// strand slopes κ coming in ± pairs (the chiral pairing of <see cref="ChiralKClaim"/>). Hence F_res's
/// roots are self-conjugate and the monic F_res is real. Checked from below at N = 5 and N = 7; the AT
/// step rests on the chiral pairing fixing the slopes, which is not derived in general.
/// <b>The β-exclusion itself needs none of this</b>: it is the multiplicity bound alone.</para>
///
/// <para><b>Scope, stated flatly.</b> This is a PER-N certificate, not a law. It retires N = 5.
/// N = 7 is NOT the same call: CertifyComplete also proves the R1 gcd, whose resultant runs against a
/// corner block of dimension 441 at N = 7 (25 at N = 5). Reaching N = 7 needs a D-only entry point.
/// It has not been built. The all-N item, the
/// codim-2 β-exotic genericity, reduced to the single scalar s₆ ≠ 0 at every forced seed, is
/// UNTOUCHED by this claim and remains open; see the parent claim's scope note and the doc's Status
/// item 2. The reading covers the branch loci of F_res only, which is sufficient because the AT
/// factor divided out of χ is REDUCING, not merely invariant: it is q-independent and invariant
/// for every q, so A and K each preserve it, and both being Hermitian each preserves its orthogonal
/// complement too. AT ⊕ F_res is genuinely block-diagonal, so no Jordan chain crosses the seam where
/// disc_Λ(F_res) would not see it, and (verified from below) all four N = 5 forced seeds lie wholly
/// in F_res.</para>
///
/// <para><b>Typed parent.</b> <see cref="SeedExistenceCountingClaim"/>: it forces the N − 1
/// real-to-complex count-drops whose CHARACTER this claim pins. Gate:
/// <c>DiscHasNoMultiplicityThreeRoot_ExcludesTheBetaExotic</c> (Category FOLDRESULTANT, ~3 s per
/// parity); live: <c>inspect --root betaexotic</c>
/// (<c>BetaExoticExclusionWitness</c>).</para></summary>
public sealed class BetaExoticExcludedAtN5Claim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: the theorem that forces the count-drops this qualifies.
    public SeedExistenceCountingClaim SeedExistence { get; }

    public BetaExoticExcludedAtN5Claim(SeedExistenceCountingClaim seedExistence)
        : base("The beta-exotic exclusion at N=5: on the (1,2) block pencil L(q) = A + qC of the 5-site XY chain " +
               "under uniform Z-dephasing, no branch locus q* != 0 of the residual charpoly factor F_res carries " +
               "the Puiseux exponent 3/2 of the beta-exotic (normal form [[0,s],[s^2,0]], eigenvalues +-s^{3/2}); " +
               "both R-parity sectors, exactly, over Q(i). The order of a disc_Lambda(F_res) zero reads the " +
               "exponent (defective EP2 -> 1, diabolic -> 2, cubic branch point -> 2, beta-exotic -> 3), and every " +
               "other colliding pair contributes non-negative order, so 'no root of multiplicity >= 3 off q = 0' " +
               "excludes the beta outright -- weaker than squarefreeness on purpose, since the diabolic loci are " +
               "genuine double roots. The certified squarefree-layer reading is [56,26] (R-odd) and [56,32] " +
               "(R-even): max multiplicity 2. The lift is one-way: reduction mod a prime ideal not dividing " +
               "lc_q(D) is a degree-preserving homomorphism, so roots can only merge and max-mult(D mod p) >= " +
               "max-mult(D); the layer prime is certified at both ends of the q-axis (true deg_q D, true " +
               "q-valuation). H1 (algebraic multiplicity exactly 2, no 3x3 Jordan) follows at N=5 by ELIMINATION, " +
               "not by the order of the zero: max-mult 2 already kills the beta (ord 3), every Jordan block of " +
               "size >= 4, and every semisimple degeneracy of >= 3 branches, and of what remains (EP2 ord 1, " +
               "diabolic ord 2, cubic branch point ord 2) ONLY the EP2 changes the real count. That needs F_res " +
               "real, which holds because T commutes with the reflection R at odd N (so each R-sector is " +
               "self-conjugate) and the AT slopes are chirally paired; checked at N=5, the AT step not derived " +
               "in general. The beta-exclusion itself needs only the multiplicity bound. SCOPE: a per-N " +
               "certificate, not a law; the all-N item (s6 != 0 at every forced seed) is untouched and open",
               Tier.Tier1Derived,
               "experiments/F89_SEED_EXISTENCE_REDUCTION.md + " +
               "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md")
    {
        SeedExistence = seedExistence ?? throw new ArgumentNullException(nameof(seedExistence));
    }

    public override string DisplayName =>
        "The β-exotic excluded exactly at N = 5: disc_Λ(F_res) has no root of multiplicity ≥ 3 off q = 0";

    public override string Summary =>
        "No branch locus of the N = 5 (1,2) residual factor carries the β-exotic's Puiseux exponent 3/2, both " +
        "R-parities, exactly over ℚ(i): a 3/2 point needs ord disc ≥ 3, and the certified squarefree layers are " +
        "[56, 26] / [56, 32]: maximum multiplicity 2. One-way lift from a prime certified at both ends of the " +
        "q-axis. H1 (algebraic multiplicity exactly 2) follows at N = 5 by elimination, on a real F_res. Per-N, " +
        "not all-N: the " +
        $"scalar s₆ ≠ 0 stays open ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the multiplicity reads the Puiseux exponent",
                summary: "ord_{q*} disc_Λ(F_res) = 1 for a defective EP2 (exponent ½), 2 for a diabolic crossing, " +
                         "2 for a cubic branch point (3×3 Jordan, exponent ⅓), and 3 for the β-exotic (exponent " +
                         "3/2): a 2-cycle with exponent e gives (λ₁−λ₂)² ~ (q−q*)^{2e}. Every other colliding pair " +
                         "contributes non-negative order, so a coincident collision raises the total (3+1, 3+2) " +
                         "and can never mask it. Max multiplicity 2 therefore excludes the β outright.");
            yield return new InspectableNode("why it is NOT the Galois route",
                summary: "the demand is strictly weaker than squarefreeness: the diabolic loci are genuine double " +
                         "roots and must survive. The Galois lever died because a full Galois group sees which " +
                         "sheets swap, never the ½-versus-3/2 exponent (counterexample: x³ − 3x + (2 + t³), full " +
                         "S₃, disc = −27t³(t³+4), an order-3 β-shaped count-drop). A multiplicity sees the exponent.");
            yield return new InspectableNode("the one-way lift, and the prime certified at both ends",
                summary: "reduction mod π_p ∤ lc_q(D) is a degree-preserving homomorphism: a factor (q−α)^m " +
                         "survives and distinct roots can only merge, so max-mult(D mod p) ≥ max-mult(D). One " +
                         "certified prime proves the bound exactly over ℚ(i). The layer prime must attain the true " +
                         "deg_q D (no root escaped to q = ∞) AND the true q-valuation (no nonzero root collapsed " +
                         "onto q = 0, where the q-power strip would discard it). A Hadamard bound certifies each " +
                         "of those two true values separately, one bound per statement and not one for their " +
                         "union; a prime attaining BOTH is searched for, and the certificate fails closed if none " +
                         "of the sampled primes does.");
            yield return new InspectableNode("H1 at N = 5 (algebraic multiplicity exactly 2), by elimination",
                summary: "DO NOT argue this by the order of the zero. The shortcut 'a count-drop flips the sign of " +
                         "the real discriminant, so the zero has odd order, so with max multiplicity 2 the order " +
                         "is 1' assumes exactly ONE conjugate pair is born at q*: two forced seeds coinciding " +
                         "there would drop the real count by 4 at an order-2, sign-PRESERVING zero. What holds " +
                         "instead: max multiplicity 2 excludes the β-exotic (ord 3), every Jordan block of size " +
                         "≥ 4 (ord ≥ 3), and every semisimple degeneracy of three or more branches (ord ≥ 3). It " +
                         "also forbids an EP2 from sharing its locus with a diabolic crossing or a cubic branch " +
                         "point, since 1 + 2 = 3. So a locus carrying a count-drop carries EP2s and nothing else " +
                         "(one, or two coincident), each a defective 2×2 Jordan block: algebraic multiplicity " +
                         "exactly 2. This needs F_res real: T = diag((−1)^{a₀+a₁+b₀}), the bipartite sign of the " +
                         "hop graph, gives T K T = −K and T A T = A, so T L(q) T = L(q)† is a similarity; at odd N " +
                         "T commutes with the reflection R, so each R-sector's spectrum is self-conjugate, and the " +
                         "AT spectrum is separately self-conjugate (chirally paired slopes). The AT step is " +
                         "checked at N = 5 and N = 7, not derived.");
            yield return new InspectableNode("the scope boundary (what stays open)",
                summary: "a PER-N certificate, not a law: it retires N = 5. N = 7 is NOT the same call: " +
                         "CertifyComplete also proves the R1 gcd, whose resultant runs against a corner " +
                         "block of dimension 441 at N = 7 (25 at N = 5), so deg_q R would be ≈ 23000 " +
                         "nodes per prime. The β-exclusion needs only D (bound ≈ 2756), so N = 7 needs a " +
                         "D-only entry point, not yet built. The all-N core, the " +
                         "codim-2 β-exotic genericity, reduced to the single scalar s₆ ≠ 0 at every forced seed, " +
                         "is untouched. Two traps recorded with it: the SIGN of s₆ is not gauge-invariant (the " +
                         "antilinear gauge has two branches, exchanged by r ↦ i·r), so never try to prove s₆ > 0; " +
                         "and |vᵀv| ≈ 0 certifies nothing, since an isotropic vector exists in every 2-dim complex " +
                         "span.");
        }
    }
}
