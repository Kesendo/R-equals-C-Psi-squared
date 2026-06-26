using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The F89 path-3 octic branch locus is a palindrome (Tier 1 derived; the mirror is forced, not
/// observed). As q = J/γ varies in the complex plane the eight octic rates braid and collide at exceptional
/// points (EPs) and one diabolic point; those collisions are mirror-symmetric about Re λ = −σ = −4γ (the
/// midpoint of the absorption rungs −2γ, −6γ, the F1 weight-complement pair n_diff = 1 ↔ 3 of the
/// <see cref="AbsorptionTheoremClaim"/>).
///
/// <para>It is FORCED by the F1 palindrome carried on the (SE,DE) block as an ANTIUNITARY symmetry
/// T = P·K (P the weight-complement rung swap, commuting with the J-hopping; K complex conjugation):
/// T L(q) T⁻¹ = −L(q̄) − 2σ (note the conjugate q̄ on the right: a same-q identity on the real axis, the
/// vertical fold, and a q → q̄ relation off it). The genuinely all-q holomorphic shadow is
/// F₈(λ, q) = F₈(−λ − 8, −q) (a q → −q identity). The spectral
/// action is antilinear, λ ↦ −λ̄ − 2σ (reflect Re about −σ, preserve Im), so the merged-eigenvalue locus
/// is invariant under it: every EP lies on the line Re λ = −σ or in a mirror pair across it, no orphan.
/// Verified on the committed octic literal to 4·10⁻¹³ (closes under λ↦−λ̄−2σ, fails the bare linear
/// λ↦−λ−2σ by 4.3) and over the whole complex-q plane to machine precision; the discriminant factorises
/// exactly as const·q²⁴·(3q⁴+q²−1)²·P₂₀(q).</para>
///
/// <para>Two separable gifts at the diabolic: it sits ON the line because its coalescing pair is
/// overlap-balanced (p = ½, dephasing scalar −4γ, the AT-midpoint, a free-fermion integrability fact,
/// <see cref="F89Path3OcticEpClaim"/>), and it is SILENT (semisimple) for the same integrability reason,
/// NOT because it is self-mirror: on-the-line does not imply silent (XXZ Δ ≠ 0 stays on the line yet
/// becomes defective; hypotheses/DIABOLIC_BY_INTEGRABILITY.md). The line is the palindrome's gift, the
/// silence is integrability's. Live witness <c>inspect --root branchpalindrome</c>
/// (<c>BranchLocusPalindromeWitness</c>); reading reflections/ON_WHO_WATCHES_WHOM.md. Scope, now checked
/// (foldlift probe, 2026-06-26): the block-internal self-fold is N_block=4 ONLY — the rung-swap P needs the
/// overlap/no-overlap multiplicity balance 2 = N−2 (DE = bar(DE) half-filling), true only at N=4 (residual
/// 3e-14 at N=4, ~1 with 0 on-line zeros at N=5,6,7). The GLOBAL palindrome Π L Π⁻¹ = −L − 2σ still lifts to
/// all N: its column bit-flip ρ[a,b]→ρ[a,bar b] pairs (SE,DE)=(w1,w2) with (SE,w_{N-2})=(w1,w_{N-2}) (the rungs
/// complement, n_diff(a,b)+n_diff(a,bar b)=N), and only the case where that partner IS (SE,DE) itself (w_{N-2}=w2,
/// i.e. N=4) gives the within-block self-fold. The `foldcross` probe confirms the cross-fold spec(SE,DE) ↔
/// spec(SE,w_{N-2}) holds to ~1e-13 about σ=N for N=4,5,6 (self-fold breaks at N≥5): the global mirror lifts as a
/// cross-block fold, the N=4 self-fold is the degenerate partner=self case. So the N=4 on-line zeros become
/// cross-block mirror partners for N≥5. The count 2=N−2 is distinct from the diabolic's overlap-fraction p=½.</para>
///
/// <para>Anchors: <c>experiments/F89_BRANCH_LOCUS_PALINDROME.md</c> +
/// <c>docs/proofs/MIRROR_SYMMETRY_PROOF.md</c> + <c>reflections/ON_WHO_WATCHES_WHOM.md</c>.</para></summary>
public sealed class F89BranchLocusPalindromeClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring (the F1 palindrome the branch locus inherits).
    public F1PalindromeIdentity Palindrome { get; }
    // Parent-edge marker for Schicht-1 wiring (the octic, the diabolic on the line, the AT-midpoint centre).
    public F89Path3OcticEpClaim Diabolic { get; }

    /// <summary>The palindrome centre Re λ = −σ = −N_block·γ for the path-3 (SE,DE) block (N_block = 4).</summary>
    public static double PalindromeCentre(double gamma, int nBlock) => -nBlock * gamma;

    public F89BranchLocusPalindromeClaim(F1PalindromeIdentity palindrome, F89Path3OcticEpClaim diabolic)
        : base("F89 path-3 octic branch locus is a palindrome: the EP/diabolic collisions are mirror-symmetric about Re λ = −σ = −4γ, FORCED by the F1 palindrome carried antiunitarily on the (SE,DE) block (T L(q) T⁻¹ = −L(q̄) − 2σ, a same-q fold at real q; all-q holomorphic identity F8(λ,q) = F8(−λ−8, −q); spectral action λ ↦ −λ̄ − 2σ; EP partner at q̄*); every EP on the line or in a mirror pair, no orphan",
               Tier.Tier1Derived,
               "experiments/F89_BRANCH_LOCUS_PALINDROME.md + " +
               "docs/proofs/MIRROR_SYMMETRY_PROOF.md + " +
               "reflections/ON_WHO_WATCHES_WHOM.md")
    {
        Palindrome = palindrome ?? throw new ArgumentNullException(nameof(palindrome));
        Diabolic = diabolic ?? throw new ArgumentNullException(nameof(diabolic));
    }

    public override string DisplayName =>
        "F89 octic branch locus is a palindrome: EPs mirror about Re λ = −4, forced by Π (antiunitary)";

    public override string Summary =>
        $"the EP/diabolic collisions mirror about Re λ = −σ = −4γ (octic closes under the antiunitary λ ↦ −λ̄ − 2σ, " +
        $"not the linear λ ↦ −λ − 2σ); every EP on the line or in a mirror pair, no orphan ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("palindrome centre Re λ = −σ (γ=1, N_block=4)", PalindromeCentre(1.0, 4));
            yield return new InspectableNode("the inherited mirror (antiunitary)",
                summary: "the (SE,DE) F1 palindrome is the antiunitary T = P·K: T L(q) T⁻¹ = −L(q̄) − 2σ (same-q fold at real q), so " +
                         "the branch locus is invariant under λ ↦ −λ̄ − 2σ (Re reflects about −σ, Im preserved); every EP " +
                         "on the line or in a mirror pair. Line = the palindrome's gift; the diabolic's silence = " +
                         "free-fermion integrability's, separate. Live: inspect --root branchpalindrome.");
        }
    }
}
