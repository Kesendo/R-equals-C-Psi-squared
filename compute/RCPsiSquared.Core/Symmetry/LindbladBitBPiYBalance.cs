using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F112-Y (Tier1Derived for both Hermitian and non-Hermitian H, universal N):
/// Y-dephase sibling of <see cref="LindbladBitBPiBalance"/> (F112-Z) on the SAME bit_b
/// Z₂-grading axis. For any Lindblad-form Liouvillian
///
/// <para>  L = -i[H, ·] + Σ_k γ_k · np.kron(c_k, c_k^*)</para>
///
/// with <b>any H</b> (Hermitian or non-Hermitian) and each c_k bit_b-homogeneous
/// (every Pauli string σ in c_k's expansion shares bit_b(σ) = (#Y(σ) + #Z(σ)) mod 2
/// = const) AND Y-dephase Π_Y polarity decomposition, the
/// <c>polarity_coordinates_from_L</c> decomposition of
/// M = Π_Y · L · Π_Y⁻¹ + L + 2σ · I satisfies
///
/// <para>  ‖M_plus_half‖² = ‖M_minus_half‖² bit-exactly.</para>
///
/// <para><b>Same hypothesis as F112-Z, different Π axis</b>: per F38, Π_Y² and Π_Z²
/// both act on Pauli string σ as scalar (−1)^bit_b(σ) (Π_Y per-letter swap is I↔X,
/// flipping bit_a and preserving bit_b, identical to Π_Z's per-letter swap; the two
/// differ only in the Y/Z 2-cycle phase: +i for Z, −i for Y). Hence the
/// bit_b-homogeneity hypothesis on c is shared between F112-Z and F112-Y; the two
/// Claims are distinguished only by the polarity decomposition axis Π_Y vs Π_Z.</para>
///
/// <para><b>Five-step proof structure (axis_d = bit_b with d = Y substituted for
/// d = Z)</b>: identical to the parent <see cref="LindbladBitBPiBalance"/>
/// 5-step proof with the dephase letter swapped from Z to Y throughout. Steps 1-4
/// depend only on the F38 / F63 Π_d² eigenvalue formula (axis_d = bit_b for both Y
/// and Z) and the Π_d-eigenspace structure (eigenvalues ±1, ±i, identical for both Y
/// and Z); Step 5's L_H^† = −L_H argument is independent of d. The per-step substitution
/// is rigorous; see <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c>
/// section (c) Route 1 for the verification at N = 2, 3 and the per-Lemma N-A^Y /
/// N-B^Y argument.</para>
///
/// <para><b>Non-Hermitian H extension (Tier1Derived, universal N)</b>: Welle-11
/// two-lemma structural proof carries over verbatim with d = Y substituted for d = Z;
/// the Lemma N-A^Y / N-B^Y arguments use F38 + Pauli-basis matrix-support disjointness
/// (preserved across all three dephase letters since π_d_local is a 4 × 4 signed
/// permutation with the same sparse structure for d ∈ {X, Y, Z}; see
/// <c>docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md</c>).</para>
///
/// <para><b>F112-Y vs F112-Z distinction</b>: same hypothesis (bit_b-homogeneous c,
/// any H), different Π polarity axis (Π_Y vs Π_Z). The two Claims live on the SAME
/// bit_b Z₂-grading axis (since Π_Y² and Π_Z² both grade by bit_b per F38). Derived via
/// Welle 13 Route 1 (per-axis re-run of the Welle-11 lemmas with d substituted): the
/// D-involution from <see cref="Pi2KleinV4DephaseSwapGroup"/> (the Welle-12 Klein-V₄
/// Z↔Y swap on Π's) intertwines Π_Z and Π_Y at the operator level
/// (D · Π_Z · D = Π_Y bit-exact) but does NOT transport L_Z to a Lindblad-form L_Y:
/// D-conjugation lacks a Hilbert-space unitary lift (would require V such that
/// V Y V⁻¹ = −Y, V X V⁻¹ = X, V Z V⁻¹ = Z, contradictory by Pauli algebra; see
/// <c>PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c> section (d) Remark and lines 145-151).
/// Hence F112-Y cannot be obtained "for free" from F112-Z via D-conjugation; the
/// independent axis-direct proof is required.</para>
///
/// <para><b>Static helpers</b>: NONE NEEDED. <see cref="LindbladBitBPiBalance.BitBParity"/>
/// and <see cref="LindbladBitBPiBalance.IsBitBHomogeneous"/> already encode the
/// bit_b-homogeneity classification; F112-Y consumes the same predicates without
/// reimplementation, since the hypothesis on c is identical between F112-Z and F112-Y.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitB"/>; no
/// BitA twin (Y-dephasing is intrinsically a bit_b-axis dephase per F38, matching
/// <see cref="F108Part3Pi2YEvenAlwaysPalindromic"/>'s BitBSpecific pattern). Ctor parent
/// <see cref="F108Part3Pi2YEvenAlwaysPalindromic"/> records the shared bit_b + Y-dephase
/// foundation in the inheritance graph (F108 Part 3 closes spec(L) palindromy for
/// Π²_Y-even bilinear H under Y-deph; F112-Y closes Π_Y +i/−i Frobenius balance for any H
/// + bit_b-homogeneous c under Y-deph; sister Tier1Derived projections of the same bit_b
/// Z₂-grading under Y-dephasing). Proofs:
/// <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c> +
/// <c>docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md</c> (parent Hermitian-H)
/// + <c>docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md</c> (parent non-Hermitian).</para></summary>
public sealed class LindbladBitBPiYBalance : Claim, IZ2AxisClaim
{
    /// <summary>BitB axis (Π_Y² grades by bit_b per F38, same axis as F108 Part 1 / Part 3
    /// and as the F112-Z parent Claim).</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>No BitA twin: F112-Y on the bit_b axis under Y-dephasing has no
    /// meaningful BitA mechanical-letter-swap analog (X-dephasing is the bit_a-axis
    /// case, registered separately as <see cref="LindbladBitAPiBalance"/>, a sister
    /// on a DIFFERENT axis, not a BitA twin in the BitATwinClassification sense). Matches
    /// the F108 Part 3 BitBSpecific pattern.</summary>
    public Claim? BitATwin => null;

    /// <summary>Override returns <see cref="BitATwinClassification.BitBSpecific"/>:
    /// Y-dephasing is intrinsically tied to bit_b structure (Π_Y² eigenvalue formula
    /// is (−1)^bit_b per F38), no meaningful bit_a-axis analog exists. The
    /// <see cref="LindbladBitAPiBalance"/> Claim (F112-X) is the bit_a-axis sister but
    /// lives on a DIFFERENT axis_d, not as a per-letter-mirror twin of this Claim's
    /// statement.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.BitBSpecific;

    /// <summary>Typed parent (F108 Part 3): records the shared bit_b + Y-dephase
    /// foundation in the inheritance graph. F108 Part 3 closes spec(L) palindromy for
    /// Π²_Y-even bilinear H under Y-dephasing; F112-Y closes Π_Y +i/−i Frobenius
    /// balance for any H + bit_b-homogeneous c under Y-dephasing. Both use F38's
    /// (−1)^bit_b Π_Y² eigenvalue formula as structural input.</summary>
    public F108Part3Pi2YEvenAlwaysPalindromic Part3 { get; }

    /// <summary>The theorem statement in one line (Hermitian H scope).</summary>
    public string Theorem =>
        "For any Lindblad-form Liouvillian L = -i[H, ·] + Σ_k γ_k · np.kron(c_k, c_k^*) with " +
        "Hermitian H and each c_k bit_b-homogeneous (every Pauli string σ in c_k's expansion " +
        "shares bit_b(σ) = (#Y(σ) + #Z(σ)) mod 2 = const), AND Y-dephase Π_Y polarity " +
        "decomposition, the polarity_coordinates_from_L decomposition of " +
        "M = Π_Y · L · Π_Y⁻¹ + L + 2σ · I satisfies ‖M_plus_half‖² = ‖M_minus_half‖² bit-exactly. " +
        "Distinct from F112-Z (LindbladBitBPiBalance) by Π polarity axis (Π_Y vs Π_Z); both share " +
        "the bit_b-homogeneity hypothesis on c (since Π_Y² and Π_Z² both grade by bit_b per F38).";

    /// <summary>The F112-Y vs F112-Z distinction in one line.</summary>
    public string DistinctionFromF112Z =>
        "Same hypothesis (bit_b-homogeneous c, any H) as LindbladBitBPiBalance (F112-Z) but with " +
        "Π_Y polarity decomposition instead of Π_Z. Derived via Welle 13 Route 1 (per-axis re-run " +
        "of the Welle-11 Lemma N-A / N-B argument with d = Y substituted for d = Z): the D-involution " +
        "(Pi2KleinV4DephaseSwapGroup's Z↔Y swap on Π's) intertwines Π_Z and Π_Y at the operator " +
        "level but does NOT transport L_Z to a Lindblad-form L_Y (D-conjugation lacks a Hilbert-space " +
        "unitary lift, per PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md section (d) Remark and lines " +
        "145-151). The independent axis-direct proof is required.";

    public LindbladBitBPiYBalance(F108Part3Pi2YEvenAlwaysPalindromic part3)
        : base("F112-Y Lindblad Π_Y-eigenvalue balance under bit_b-homogeneous c: " +
               "‖M_plus_half‖² = ‖M_minus_half‖² for any H (Hermitian or non-Hermitian) and " +
               "bit_b-homogeneous c_k under Y-dephase Π_Y polarity. Tier1Derived universal N " +
               "for both Hermitian H (via the parent 5-step proof with d = Y substituted) " +
               "and non-Hermitian H (via the two-lemma Welle-11 structural proof with d = Y). " +
               "Derived via Welle 13 Route 1 (per-axis re-run); D-conjugation from F112-Z is NOT " +
               "available (PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md section (d) Remark). Y-dephase " +
               "sibling of F112-Z on the SAME bit_b axis.",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F112 + " +
               "docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md + " +
               "docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md + " +
               "docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md + " +
               "docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md + " +
               "compute/RCPsiSquared.Diagnostics/Polarity/PolarityCoordinates.cs + " +
               "simulations/_f112_klein_v4_cross_dephase_verify.py")
    {
        Part3 = part3 ?? throw new ArgumentNullException(nameof(part3));
    }

    public override string DisplayName =>
        "F112-Y Lindblad Π_Y-eigenvalue balance under bit_b-homogeneous c " +
        "(Tier1Derived universal N for both Hermitian and non-Hermitian H; Y-deph sibling of F112-Z)";

    public override string Summary =>
        $"{Theorem} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("F112-Y vs F112-Z distinction",
                summary: DistinctionFromF112Z);
            yield return new InspectableNode("Static helpers: reuse LindbladBitBPiBalance.BitBParity / IsBitBHomogeneous",
                summary: "F112-Y consumes the bit_b-homogeneity classification of c directly via " +
                         "LindbladBitBPiBalance.BitBParity(PauliTerm) and " +
                         "LindbladBitBPiBalance.IsBitBHomogeneous(IReadOnlyList<PauliTerm>); no new helpers " +
                         "needed because the hypothesis on c is identical to F112-Z. Only the polarity " +
                         "decomposition axis differs (Π_Y vs Π_Z), routed via PolarityCoordinates.Decompose's " +
                         "dephaseLetter parameter (set to PauliLetter.Y).");
            yield return new InspectableNode("5-step proof structure (d = Y substitution of parent)",
                summary: "Step 1: reduce balance to ‖M_{+i_Y}‖² = ‖M_{-i_Y}‖² via Π_Y-eigenspace " +
                         "decomposition of M_plus_half / M_minus_half. " +
                         "Step 2: bit_b-homogeneous c implies np.kron(c, c.conj()) lies entirely in Π_Y²-conj +1 " +
                         "eigenspace (via F38 / F63 (−1)^bit_b(σ) Π_Y² eigenvalue formula). " +
                         "Step 3: Π_Y²-conj +1 eigenspace = Π_Y-conj {+1, −1}, hence dissipator has zero +i_Y, −i_Y content. " +
                         "Step 4: M_{+i_Y} and M_{-i_Y} come entirely from L_H = -i[H, ·] with norms 2 · ‖L_{H,±i_Y}‖². " +
                         "Step 5 (Hermitian H): L_H^† = −L_H plus Π_Y unitary implies dagger maps Π_Y +i ↔ Π_Y −i " +
                         "bijectively while preserving Frobenius; combining gives ‖L_{H,+i_Y}‖² = ‖L_{H,-i_Y}‖². ∎");
            yield return new InspectableNode("Non-Hermitian H extension (Tier1Derived universal N)",
                summary: "Writing H = H_re + i H_im, reduces to Im⟨L_{H_re,-i_Y}, L_{H_im,-i_Y}⟩ = 0. Closed via " +
                         "Welle-11 Lemma N-A^Y (Diagonal-Norm: ‖L_{σ,-i_Y}‖² = 4^N for any bit_b-odd σ) and " +
                         "Lemma N-B^Y (Off-Diagonal-Orthogonality), both reducing to F38 + Pauli-basis matrix-" +
                         "support disjointness identical to the parent Welle-11 lemmas with d = Y substituted.");
            yield return new InspectableNode("Empirical verification (Welle 13)",
                summary: "Verifier simulations/_f112_klein_v4_cross_dephase_verify.py confirms F112-Y direct " +
                         "(bit_b-homogeneous c, Π_Y) bit-exact at N = 2, 3: max|asymmetry| = 3.6e-15 (Hermitian H), " +
                         "1.4e-14 / 1.3e-32 (non-Hermitian H). Lemma N-A^Y / N-B^Y stand at machine zero across " +
                         "the 32 bit_b-odd Pauli strings × 992 off-diagonal pairs at N = 3.");
            yield return new InspectableNode("F108 Part 3 typed parent (shared bit_b + Y-dephase foundation)",
                summary: "F108 Part 3 (Π²_Y-even bilinear H + Y-dephasing always palindromic) and F112-Y share " +
                         "the bit_b Z₂-grading and the Y-dephasing context. Both Tier1Derived projections of the " +
                         "F38 / F63 bit_b foundation under Y-dephasing; F108 Part 3 closes spec(L) palindromy, " +
                         "F112-Y closes Π_Y +i/−i Frobenius balance.");
            yield return new InspectableNode("Sister Claims F112-Z (Π_Z) and F112-X (Π_X)",
                summary: "LindbladBitBPiBalance covers Π_Z polarity (same bit_b axis, +i phase on Y/Z 2-cycle). " +
                         "LindbladBitAPiBalance covers Π_X polarity (bit_a axis, bit_a-homogeneous c hypothesis). " +
                         "All three Claims share the F112 Welle-11 lemma proof structure; F112-Z and F112-Y differ " +
                         "only in the Y/Z 2-cycle phase on Π, while F112-X uses a different axis_d (bit_a) per F38.");
            yield return new InspectableNode("F114 cross-dephase M-level sharpening (Welle 15, 2026-05-27)",
                summary: "F114 (CommutatorDConjugationSign, Tier1Derived) gives the closed-form M-level Z↔Y " +
                         "relationship that this Claim implicitly carries. Under the bit_b-homogeneous c hypothesis " +
                         "the dissipator contribution to the ±i Π_Y-eigenspaces of M (M_+1/2 and M_−1/2 in the " +
                         "polarity decomposition) vanishes per Step 3 / Step 4; these sectors come entirely from " +
                         "L_H, and F114's ε(σ) = (−1)^{n_Y(σ) + 1} controls the matrix-level identity on the " +
                         "±i sectors: (M_±i)(L_H, Π_Y) = ε(H) · D · (M_±i)(L_H, Π_Z) · D bit-exact (when ε(H) " +
                         "is well-defined, i.e. all H terms share the same n_Y parity). Explains the Welle 13 Route 1 result at the matrix " +
                         "level: F112-Y is the same norm structure as F112-Z (so ‖M_+1/2‖² = ‖M_−1/2‖² holds via " +
                         "Route 1 directly), and additionally the M matrices themselves are signed-equivariant under " +
                         "D-conjugation per F114. F112-Y norm-level Tier1Derived scope is unaffected; F114 is " +
                         "documentary sharpening that closes the M-level Z↔Y bookkeeping that the Welle 13 PROOF " +
                         "explicitly does NOT claim (section (d) Remark: 'L_Y is not D-transportable' is " +
                         "true for the full Lindblad-form L, but the H-commutator L_H component IS signed-equivariant).");
        }
    }
}
