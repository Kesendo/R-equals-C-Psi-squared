using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89d (Tier 1 derived): the (SE,DE) = (w1,w2) coherence block and its branch-locus-palindrome partner
/// (SE, w_{N−2}) = (w1, N−2) are related by an EXACT antiunitary similarity at the MATRIX level,
///
/// <code>
///     L(1, N−2)(q̄)  =  −P · conj(L(1,2)(q)) · Pᵀ  −  2N·I
/// </code>
///
/// where P is the bra-complement permutation (the F89c bra bit-flip ρ[a,b] → ρ[a,b̄], with n_diff(a,b̄) = N −
/// n_diff(a,b)). Verified to MACHINE ZERO for N = 4..9 at every q, real and complex (the partner block evaluated
/// at the conjugate coupling q̄, the F1 antiunitary form T L(q) T⁻¹ = −L(q̄) − 2σ).
///
/// <para>INTEGRABILITY-INDEPENDENT (the (q,Δ) extension, 2026-06-30): the identity holds for the FULL interacting
/// XXZ block L(q,Δ) at every anisotropy Δ, not just the integrable XY one (machine zero, N=4..9, all q, all Δ).
/// The reason is that the Δ·ZZ term is EVEN under the global bit-flip (Z_bZ_{b+1} ↦ (−Z_b)(−Z_{b+1}) = Z_bZ_{b+1},
/// so zz(b̄) = zz(b)), so the bra-complement carries it cleanly. The diabolics themselves DIE under Δ
/// (integrability-protected, the arc's Move 2), but the pairing STRUCTURE does not: a diabolic and its partner
/// turn defective in lockstep. The discriminant is bit-flip PARITY: a bit-flip-ODD perturbation breaks the fold (a
/// longitudinal Z-field Σ_k w_k Z_k has fe(b̄) = −fe(b), residual O(1)). So the cross-fold is a structural/
/// algebraic property of the Liouvillian, NOT a free-fermion artifact.</para>
///
/// <para>This is the matrix-level UPGRADE of <see cref="F89BranchLocusPalindromeClaim"/>'s spectrum-level cross-
/// fold (which states only spec(SE,DE) ↔ spec(SE,w_{N−2}) about σ = N). An antiunitary similarity preserves the
/// whole Jordan structure: complex conjugation, the permutation similarity P·Pᵀ, and the affine −(·) − 2N each
/// preserve the dimensions of generalized eigenspaces, so a SEMISIMPLE coalescence maps to a semisimple
/// coalescence. Hence every (SE,DE) diabolic at (q, λ) has a partner diabolic at (q̄, −λ̄ − 2N) in the
/// (SE,w_{N−2}) block with IDENTICAL character and coalescence gap, for all N and all q at once. The N=4 within-
/// block self-fold (the palindrome that placed one diabolic on the real axis, <see cref="F89BranchLocusPalindromeClaim"/>)
/// is the degenerate w_{N−2} = w2 partner = self case; for N ≥ 5 the N=4 on-line "zeros" become these cross-block
/// mirror partners.</para>
///
/// <para>Verified bit-exact: the N=7 real-q diabolic (q = 1.1264, λ = −4.942, full-block coalescence gap
/// 4.19·10⁻⁵) pairs with the (1,5) block at the fold image −λ − 2N = −9.058, same gap. Live witness
/// <c>inspect --root crossfold</c> (<c>CrossFoldSimilarityWitness</c>, on the shared
/// <c>WeightCoherenceBlock</c> builder + the bra-complement permutation). Anchors:
/// <c>experiments/F89_PATH_K_DIABOLIC.md</c> (the cross-fold section) + <c>docs/ANALYTICAL_FORMULAS.md</c> (F89d)
/// + <c>experiments/F89_BRANCH_LOCUS_PALINDROME.md</c>. The Move-4 grounding of the diabolic_over_higher_n
/// arc.</para></summary>
public sealed class F89CrossFoldSimilarityClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: the F1 mirror λ ↦ −λ̄ − 2σ the cross-fold realises across blocks.
    public F1PalindromeIdentity Palindrome { get; }
    // Parent-edge marker for Schicht-1 wiring: the spectrum-level cross-fold this claim upgrades to a similarity.
    public F89BranchLocusPalindromeClaim SpectrumFold { get; }

    /// <summary>The cross-fold centre −σ = −N·γ (the AT-rate Hamming-complement pair-sum 2γ·N halved): the Re λ
    /// of the (SE,w_{N−2}) partner reflects about this.</summary>
    public static double FoldCentre(double gamma, int n) => -n * gamma;

    /// <summary>The fold image of a REAL (SE,DE) eigenvalue λ at site count N: the partner block's eigenvalue
    /// −λ − 2N (the real-λ case of −λ̄ − 2N, γ = 1).</summary>
    public static double FoldImageReal(double lambda, int n) => -lambda - 2.0 * n;

    public F89CrossFoldSimilarityClaim(F1PalindromeIdentity palindrome, F89BranchLocusPalindromeClaim spectrumFold)
        : base("F89d cross-fold similarity: the (SE,DE)=(w1,w2) <-> (SE,w_{N-2})=(w1,N-2) fold is an EXACT antiunitary similarity L(1,N-2)(qbar) = -P conj(L(1,2)(q)) P^T - 2N I (machine zero, N=4..9, all q), upgrading the branch-locus palindrome's spectrum match to a Jordan-structure-preserving similarity, so every (SE,DE) diabolic at (q,lambda) pairs with a (SE,w_{N-2}) diabolic at (qbar,-lambdabar-2N) with identical character and gap; the N=4 self-fold is the degenerate partner=self case",
               Tier.Tier1Derived,
               "experiments/F89_PATH_K_DIABOLIC.md + " +
               "docs/ANALYTICAL_FORMULAS.md + " +
               "experiments/F89_BRANCH_LOCUS_PALINDROME.md")
    {
        Palindrome = palindrome ?? throw new ArgumentNullException(nameof(palindrome));
        SpectrumFold = spectrumFold ?? throw new ArgumentNullException(nameof(spectrumFold));
    }

    public override string DisplayName =>
        "F89d: (SE,DE) <-> (SE,w_{N-2}) cross-fold is an exact antiunitary similarity (the diabolics pair)";

    public override string Summary =>
        "L(1,N-2)(qbar,Delta) = -P conj(L(1,2)(q,Delta)) P^T - 2N I exactly (N=4..9, all q, all Delta): a Jordan-" +
        "preserving antiunitary similarity, so every (SE,DE) diabolic at (q,lambda) pairs with a (SE,w_{N-2}) " +
        "diabolic at (qbar, -lambdabar-2N), character and gap preserved; integrability-independent (holds for the " +
        $"full interacting XXZ block, the discriminant is bit-flip parity) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("cross-fold centre -σ = -N (γ=1, N=7)", FoldCentre(1.0, 7));
            yield return InspectableNode.RealScalar("fold image of the N=7 real-q diabolic λ=-4.942 (-λ-2N)",
                FoldImageReal(-4.942, 7));
            yield return new InspectableNode("the matrix-level upgrade (Jordan-preserving)",
                summary: "the branch-locus palindrome gives spec(SE,DE) <-> spec(SE,w_{N-2}) about σ=N; this claim shows the " +
                         "carrier is an EXACT antiunitary similarity L(1,N-2)(qbar) = -P conj(L(1,2)(q)) P^T - 2N I (P the bra-" +
                         "complement permutation), which preserves the whole Jordan structure, so the diabolics pair with identical " +
                         "character and gap. N=7: (q=1.1264, λ=-4.942) <-> (1,5) at -9.058, gap 4.19e-5 both. Live: inspect --root crossfold.");
            yield return new InspectableNode("integrability-independent: survives XXZ anisotropy (the (q,Δ) extension)",
                summary: "the identity holds for the FULL interacting XXZ block L(q,Δ) at every Δ (machine zero, N=4..9, all q): " +
                         "the Δ·ZZ term is EVEN under the global bit-flip (zz(b̄)=zz(b)), so the bra-complement carries it cleanly. " +
                         "The diabolics die under Δ (Move 2, integrability-protected), but the pairing structure does not: a " +
                         "diabolic and its partner turn defective in lockstep. Discriminant = bit-flip parity: a longitudinal " +
                         "Z-field (odd, fe(b̄)=−fe(b)) breaks the fold (residual O(1)). So the fold is structural, not free-fermion.");
        }
    }
}
