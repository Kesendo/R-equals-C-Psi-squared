using System;
using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F150 (Tier 1 derived on the full charpoly, Tier 1 candidate on the residual): the FOLD CHECKERBOARD,
/// a Λ-side coefficient-support parity on fold-fixed blocks, and its corollary that the discriminant is even in q
/// and therefore real.
///
/// <para>Take a block fixed by EITHER fold leg, that is ket weight = N/2 OR bra weight = N/2 at even N
/// (PROOF_CODIM1_BY_ADDITIVITY §"the fold lattice", consequence (b), which names both legs). The bare leg is not an
/// antiunitary but the LINEAR, q-REVERSING rearrangement P·S(q)·P = 2Λ₀ − S(−q) with Λ₀ = −2N, where S is the
/// ×2-CLEARED R-sector pencil and NOT the block itself (in the uncleared scale the same identity reads
/// P·L(q)·P = −2N·I − L(−q), so there Λ₀ = −N, the palindrome centre the rest of the repo uses). P commutes with
/// the site reflection, so it acts inside each R-sector. Taking determinants gives
/// χ_q(Λ₀ + u) = (−1)^d·χ_{−q}(Λ₀ − u) with d = deg_Λ, at every q, real or complex, over ℤ[i]. Writing
/// χ = Σ c_{j,k}·u^j·q^k that forces c_{j,k} = 0 unless j + k ≡ d (mod 2): a checkerboard SUPPORT. The reality half
/// is the separate, unconditional q-side antiparity, so the surviving cells are real on even j and purely imaginary
/// on odd j. The branch sign ε = (−1)^d is FORCED, not measured: the residual is monic, so the cell (d, 0) settles
/// it before any other coefficient is read. Two things about that sign are worth stating rather than leaving
/// for a reader to hit. The realness pattern quoted above, real on even j and imaginary on odd j, is the ε = +1
/// branch; at odd d the pattern flips, and <c>CheckerboardRead</c> carries both as Sign ±1. And every committed
/// row has EVEN d (resDeg 8, 42, 48, 224; on a fold-fixed block the dimension C(N,p)·C(N,N/2) is even because
/// C(2m,m) is), so ε = +1 throughout and the ε = −1 branch has never been exercised.</para>
///
/// <para>THE COROLLARY IS WHAT IT WAS DERIVED FOR. The law gives F_res(−u, −q) = ε·F_res(u, q), so the Λ-roots
/// negate under q ↦ −q while the discriminant is invariant under that negation: D = disc_Λ(F_res) is EVEN in q.
/// Composed with Y1 (conj-coeff(D)(q) = D(−q)) that forces conj-coeff(D) = D, so D is REAL on every fold-fixed
/// block. This derives the disc-reality prior the multi-sector monodromy work had been asserting, replacing a
/// priced multi-hour mod-p read with an algebraic fact. The device that consumes it halves its sample points by
/// interpolating D(q) = E(q²), realised at 1.76× rather than the naive 2×.</para>
///
/// <para>MEASURED. (1,2)@N=4 sign +1 with 6 odd and 15 even cells, at odd character, where all three committed
/// boolean reads are False, so the checkerboard is the whole of that block's structure and cannot pass vacuously;
/// (2,1)@N=4 sign +1, the KET-leg control, which exists because a first draft scoped the law to the bra leg and
/// this row refutes that; (1,3)@N=6 sign +1 with 0 odd cells and 247 (R-even) / 325 (R-odd) even ones, the
/// degenerate side where an even character empties the odd cells; (1,2)@N=6 the NEGATIVE CONTROL, fold-fixed on
/// neither leg, sign 0 on F_res and on the full charpoly. Confirmed at (1,4)@N=8 at resDeg 224, sign +1 on both,
/// odd cells occupied; that run gated the cell counts but did not record them, so the cell counts at (1,4)@N=8 do
/// not exist and must not be invented. The N=4 closed form had the law visible and unread all along:
/// disc(F_8) = const·q²⁴·(3q⁴+q²−1)²·P₁₀(q²) has every factor even in q.</para>
///
/// <para>SCOPE, and the third fence is the one that matters. Even N and fold-fixed blocks only; off them the Λ side
/// has no source and sign 0 is the correct answer, which is what the negative control measures. The identity is the
/// BARE leg, giving a Λ-reflection tied to q ↦ −q and NOT Λ-evenness at fixed q; reading it as the latter was a
/// first-draft error. On the full charpoly χ this is a theorem with no premise, but the descent to F_res = χ/AT
/// additionally needs the AT factor to carry the same parity, and that premise is MEASURED, not proved, which is
/// why a failure there is a finding rather than a bug. It says nothing about coefficient-realness of F_res itself,
/// a different object from realness of D: (1,4)@N=8 has D real and F_res coefficient-complex at once. Δ = 0 only,
/// inheriting the gauge's scope. And it is not a census: a support parity and a disc-realness are neither an
/// absence verdict nor a defective-locus count.</para>
///
/// <para>Do not conflate the two reads in the implementation: the D-device consumes the strictly weaker
/// <c>QNegationSupportParityExact</c>, which asks only for one j+k parity and not for the realness pattern, while
/// <c>CheckerboardRead</c> is the full ±1 law and would decline blocks whose support parity is perfectly good.</para>
///
/// <para>Gates <c>FoldCheckerboardParityTests</c> (4 theories × 2 parities = 8, Category EVEN_DISC13_SCOUT except
/// the (1,3)@N=6 row which is SLOW_EVEN_DISC13) and, for the (1,4)@N=8 confirmation, <c>Disc14N8ScoutTests</c>
/// (Category SLOW_EVEN_DISC14, over two hours). Implementation
/// <c>FoldResultantCertificate.CheckerboardRead</c>. Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> (F150) +
/// <c>experiments/F89_PATH_K_DIABOLIC.md</c> + the arc <c>sideways_spin_ladder</c>.</para></summary>
public sealed class FoldCheckerboardParityClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: the fold lattice supplies the leg this law takes determinants of.
    public F89CrossFoldSimilarityClaim FoldLattice { get; }
    // Parent-edge marker for Schicht-1 wiring: the Δ=0 gauge. NOT the source of the q-side antiparity (that is
    // Y1, the pencil reality, which is separate and unconditional): the edge is the Δ=0 scope this law inherits,
    // plus the fact that the gauge's character decides whether the odd-j cells are occupied at all.
    public BipartiteGaugeCriterionClaim Gauge { get; }

    /// <summary>Whether a block is fixed by EITHER fold leg, which is the law's hypothesis. Both legs, deliberately:
    /// scoping this to the bra leg was the first-draft error the (2,1)@N=4 control refutes.</summary>
    public static bool IsFoldFixed(int n, int wKet, int wBra) =>
        n % 2 == 0 && (wKet == n / 2 || wBra == n / 2);

    /// <summary>The forced branch sign ε = (−1)^d, d the Λ-degree. Monicity settles it at the cell (d, 0); the
    /// gates cross-check against the measured residual degree rather than pretending to measure the sign.</summary>
    public static int ForcedBranchSign(int lambdaDegree) => (lambdaDegree & 1) == 0 ? 1 : -1;

    public FoldCheckerboardParityClaim(F89CrossFoldSimilarityClaim foldLattice, BipartiteGaugeCriterionClaim gauge)
        : base("F150 fold checkerboard: on a block fixed by EITHER fold leg (ket or bra weight = N/2 at even N) the " +
               "bare leg is the LINEAR q-reversing rearrangement P S(q) P = 2*Lambda0 - S(-q), so taking determinants " +
               "chi_q(Lambda0+u) = (-1)^d chi_{-q}(Lambda0-u) at every q over Z[i], forcing the coefficient of " +
               "u^j q^k to vanish unless j+k = d (mod 2), a checkerboard SUPPORT with the survivors real on even j " +
               "and imaginary on odd j (branch sign eps = (-1)^d forced by monicity, not measured); corollary " +
               "F_res(-u,-q) = eps F_res(u,q), so D = disc_Lambda(F_res) is EVEN in q and, with Y1, REAL on every " +
               "fold-fixed block, which derives the B3 disc-reality prior algebraically. A theorem on the full " +
               "charpoly; on F_res the descent needs the AT factor to carry the same parity, which is MEASURED",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md + " +
               "experiments/F89_PATH_K_DIABOLIC.md + " +
               "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md")
    {
        FoldLattice = foldLattice ?? throw new ArgumentNullException(nameof(foldLattice));
        Gauge = gauge ?? throw new ArgumentNullException(nameof(gauge));
    }

    public override string DisplayName =>
        "F150: the fold checkerboard (support parity j+k = d mod 2, hence disc even in q and real)";

    public override string Summary =>
        "on a fold-fixed block (EITHER leg, even N) the bare leg is a linear q-reversing rearrangement, and its " +
        "determinant identity chi_q(Lambda0+u) = (-1)^d chi_{-q}(Lambda0-u) forces a checkerboard support " +
        "j+k = d (mod 2), real on even j and imaginary on odd j. Corollary: disc_Lambda(F_res) is even in q and, " +
        "with Y1, REAL on every fold-fixed block, closing the B3 prior algebraically and replacing a priced " +
        "multi-hour mod-p read. Measured sign +1 at (1,2)@N=4 (6 odd / 15 even cells, non-vacuous), (2,1)@N=4 (the " +
        "ket-leg control), (1,3)@N=6 (0 odd, 247/325 even) and (1,4)@N=8 at resDeg 224; sign 0 at the (1,2)@N=6 " +
        "negative control. A THEOREM on the full charpoly with no premise; the descent to F_res is a MEASUREMENT, " +
        $"and it is not a defective-locus count ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("both legs, and why the control exists",
                summary: "the hypothesis is ket weight = N/2 OR bra weight = N/2, which is what " +
                         "PROOF_CODIM1_BY_ADDITIVITY §7 consequence (b) has said all along. A first draft scoped it to the " +
                         "bra leg and would have passed every other row; the (2,1)@N=4 gate is the ket-leg control that " +
                         $"refutes that scoping. IsFoldFixed(4,2,1) = {IsFoldFixed(4, 2, 1)}, IsFoldFixed(6,1,2) = " +
                         $"{IsFoldFixed(6, 1, 2)} (the negative control, fixed on neither leg, sign 0).");
            yield return new InspectableNode("the branch sign is forced, not measured",
                summary: "eps = (-1)^d, and since the residual is monic the cell (d, 0) settles the branch before any other " +
                         "coefficient is read. Both branches are still evaluated in the implementation, as a self-check that " +
                         "they are mutually exclusive, and the gates cross-check the verdict against the measured residual " +
                         $"degree. ForcedBranchSign(224) = {ForcedBranchSign(224)} at (1,4)@N=8.");
            yield return new InspectableNode("theorem on chi, measurement on F_res",
                summary: "on the full characteristic polynomial the determinant identity has NO premise. The descent to " +
                         "F_res = chi/AT additionally needs the AT factor to carry the same parity, and that is the arc's " +
                         "ungated AT-swap premise. Reading both is what turns the premise into a measurement, and a failure " +
                         "there is a finding about the construction rather than a bug.");
            yield return new InspectableNode("what it does NOT say",
                summary: "nothing about coefficient-realness of F_res, which is a different object from realness of " +
                         "D = disc_Lambda(F_res): (1,4)@N=8 has D real and F_res coefficient-COMPLEX at the same time. It is " +
                         "also not a Lambda-evenness at fixed q, the first-draft misreading, but a Lambda-reflection tied to " +
                         "q -> -q. And it is not a census: a support parity and a disc-realness are neither an absence verdict " +
                         "nor a defective-locus count. The (1,4)@N=8 cell counts were gated but not recorded, so they do not " +
                         "exist and must not be invented.");
        }
    }
}
