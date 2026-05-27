using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F112-X (Tier1Derived for both Hermitian and non-Hermitian H, universal N):
/// BitA-axis sister of <see cref="LindbladBitBPiBalance"/>. For any Lindblad-form
/// Liouvillian
///
/// <para>  L = -i[H, ·] + Σ_k γ_k · np.kron(c_k, c_k^*)</para>
///
/// with <b>any H</b> (Hermitian or non-Hermitian) and each c_k bit_a-homogeneous
/// (every Pauli string σ in c_k's expansion shares bit_a(σ) = (#X(σ) + #Y(σ)) mod 2
/// = const) AND X-dephase Π_X polarity decomposition, the
/// <c>polarity_coordinates_from_L</c> decomposition of
/// M = Π_X · L · Π_X⁻¹ + L + 2σ · I satisfies
///
/// <para>  ‖M_plus_half‖² = ‖M_minus_half‖² bit-exactly.</para>
///
/// <para><b>Why bit_a, not bit_b</b>: per F38 (PiOperator.SquaredEigenvalue),
/// Π_X² acts on Pauli string σ as scalar (−1)^bit_a(σ) (Π_X per-letter swap is I↔Z
/// flipping bit_b and preserving bit_a). The axis-relevant homogeneity for F112-X is
/// therefore bit_a, exactly as for F112-Z / F112-Y on the bit_b axis. The proof
/// structure is the Welle-11 Lemma N-A / N-B argument re-run with axis_d := bit_a
/// substituted for bit_b; see <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c>
/// section (c) Route 1 for the per-step substitution and section (d) Route 2 for the
/// independent Hadamard-transport derivation.</para>
///
/// <para><b>Five-step proof structure (axis_d = bit_a substitution of the parent
/// Hermitian-H proof)</b>:</para>
/// <list type="number">
///   <item>Π_X-eigenspace decomposition of M reduces the balance condition to
///         ‖M_{+i_X}‖² = ‖M_{-i_X}‖² (asymmetry = (1/2)(‖M_{+i_X}‖² − ‖M_{-i_X}‖²)).</item>
///   <item>bit_a-homogeneous c gives Π_X² · np.kron(c, c^*) · Π_X⁻² = ε² · np.kron(c, c^*)
///         with ε ∈ {+1, −1}; hence np.kron(c, c^*) lies entirely in the Π_X²-conj
///         +1 eigenspace, via the F38 / F63 Π_X² eigenvalue formula on Pauli strings.</item>
///   <item>Π_X²-conj +1 eigenspace = Π_X-conj {+1, −1}; the dissipator part of M has
///         zero +i_X and zero −i_X Π_X-conjugation content.</item>
///   <item>M_{+i_X} and M_{-i_X} come entirely from L_H = -i[H, ·], with norms
///         ‖M_{±i_X}‖² = 2 · ‖L_{H,±i_X}‖².</item>
///   <item>For Hermitian H: L_H^† = −L_H plus Π_X unitary implies dagger maps Π_X +i ↔ Π_X −i
///         bijectively while preserving Frobenius; combining gives
///         ‖L_{H,+i_X}‖² = ‖L_{H,-i_X}‖². ∎</item>
/// </list>
///
/// <para><b>Non-Hermitian H extension (Tier1Derived, universal N)</b>: writing
/// H = H_re + i H_im with both summands Hermitian, the equality reduces algebraically
/// to the identity Im⟨L_{H_re,-i_X}, L_{H_im,-i_X}⟩ = 0. The Welle-11 two-lemma
/// structural proof carries over verbatim with axis_d := bit_a substituted: Lemma N-A^X
/// (Diagonal-Norm: ‖L_{σ,-i_X}‖² = 4^N for any bit_a-odd Pauli string σ) and Lemma N-B^X
/// (Off-Diagonal-Orthogonality: ⟨L_{σ_α,-i_X}, L_{σ_β,-i_X}⟩ = 0 for σ_α ≠ σ_β both
/// bit_a-odd). Both lemmas reduce to F38 + Pauli-basis matrix-support disjointness, which
/// transfer cleanly between dephase letters (the per-site π_d_local matrix is a 4 × 4
/// signed permutation for every d ∈ {X, Y, Z} with the same sparse structure;
/// see <c>docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md</c>).</para>
///
/// <para><b>Two independent derivation routes (Welle 13)</b>: (Route 1) direct re-run of
/// the Welle-11 lemmas with axis_d := bit_a; (Route 2) Hadamard transport from F112-Z to
/// F112-X via Q_zx = H · D = Hadamard^⊗N on operator space. Q_zx is the operator-space
/// lift of the per-site Hadamard rotation U_H, which acts on Pauli letters by
/// X ↔ Z, Y → −Y, I → I, swapping bit_a ↔ bit_b. Hence a bit_b-homogeneous c for F112-Z
/// transforms to a bit_a-homogeneous c' for F112-X under simultaneous Q_zx-conjugation of
/// L and Π, preserving Frobenius norms. Routes 1 and 2 are complementary: Route 1 covers
/// all three d ∈ {X, Y, Z} structurally, Route 2 provides the explicit Π_Z ↔ Π_X
/// intertwiner via the Hilbert-space Hadamard^⊗N lift.</para>
///
/// <para><b>BitB sister Claim</b>: <see cref="LindbladBitBPiBalance"/> is the bit_b-axis
/// case (Π_Z or Π_Y polarity with bit_b-homogeneous c). The two Claims are sisters on
/// opposite Z₂-grading axes (bit_a vs bit_b), not in the mechanical-letter-swap sense of
/// BitATwinClassification — the existing BitB Claim's BitATwinStatus remains
/// <see cref="BitATwinClassification.BitBSpecific"/> because F112-X is a sibling on a
/// DIFFERENT axis (axis_d differs), not a per-letter-mirror twin of the same statement.
/// The mutual cross-reference is captured in this Claim's docstring and inspectables.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitA"/>; ctor
/// parent <see cref="F108Part2Pi2XEvenAlwaysPalindromic"/> records the shared bit_a +
/// X-dephase foundation in the inheritance graph (F108 Part 2's Π²_X-even bilinear set
/// {ZZ, XX, XY, YX, YY} sits on the same bit_a Z₂-grading that F112-X's c-homogeneity
/// hypothesis fixes). Proofs:
/// <c>docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c> +
/// <c>docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md</c> (parent Hermitian-H)
/// + <c>docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md</c> (parent non-Hermitian).</para></summary>
public sealed class LindbladBitAPiBalance : Claim, IZ2AxisClaim
{
    /// <summary>BitA axis (bit_a parity Z₂-grading of the Pauli group; same axis as
    /// F108 Part 2, anchored in F38 Π_X² eigenvalue formula).</summary>
    public Z2Axis Z2Axis => Z2Axis.BitA;

    /// <summary>BitA-axis Claims have no BitATwin slot (the twin concept lives on BitB-axis
    /// Claims pointing at BitA siblings; see <see cref="BitATwinClassification"/>).</summary>
    public Claim? BitATwin => null;

    /// <summary>Override returns <see cref="BitATwinClassification.NotApplicableForThisAxis"/>:
    /// matches the F108 Part 2 pattern for BitA-axis Claims. The structural BitB sister
    /// <see cref="LindbladBitBPiBalance"/> (F112-Z / F112-Y on the opposite axis) is referenced
    /// in this Claim's docstring and inspectables but is not the mechanical-letter-swap twin
    /// captured by <see cref="BitATwinClassification.Filled"/>: F112-X and F112-Z/F112-Y live
    /// on different axis_d values (bit_a vs bit_b), not on different per-letter mirrors of
    /// the same statement.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.NotApplicableForThisAxis;

    /// <summary>Typed parent (F108 Part 2): records the shared bit_a + X-dephase
    /// foundation in the inheritance graph. F112-X uses the F38 / F63 Π_X² eigenvalue
    /// formula on Pauli strings exactly as F108 Part 2 does; F108 Part 2's Π²_X-even
    /// bilinear set {ZZ, XX, XY, YX, YY} is the bit_a = 0 family on the same Z₂-grading
    /// F112-X's c-homogeneity hypothesis lives on.</summary>
    public F108Part2Pi2XEvenAlwaysPalindromic Part2 { get; }

    /// <summary>The theorem statement in one line (Hermitian H scope).</summary>
    public string Theorem =>
        "For any Lindblad-form Liouvillian L = -i[H, ·] + Σ_k γ_k · np.kron(c_k, c_k^*) with " +
        "Hermitian H and each c_k bit_a-homogeneous (every Pauli string σ in c_k's expansion " +
        "shares bit_a(σ) = (#X(σ) + #Y(σ)) mod 2 = const), AND X-dephase Π_X polarity " +
        "decomposition, the polarity_coordinates_from_L decomposition of " +
        "M = Π_X · L · Π_X⁻¹ + L + 2σ · I satisfies ‖M_plus_half‖² = ‖M_minus_half‖² bit-exactly.";

    /// <summary>The two independent derivation routes from
    /// <c>PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md</c>.</summary>
    public string DerivationRoutes =>
        "Route 1 (direct axis re-run): re-run of the Welle-11 Lemma N-A / N-B argument with " +
        "axis_d := bit_a substituted for bit_b; uses only F38 (Π_X² eigenvalue (−1)^bit_a(σ) " +
        "on Pauli strings) and Pauli-basis matrix-support disjointness (preserved across all " +
        "three dephase letters since π_d_local is a 4 × 4 signed permutation with the same " +
        "sparse structure for d ∈ {X, Y, Z}). Route 2 (Hadamard transport): Q_zx-conjugation " +
        "of L and Π via Q_zx = H · D = U_H^⊗N ⊗ (U_H^⊗N)^* (operator-space lift of the per-site " +
        "Hadamard rotation U_H). U_H acts on Pauli letters by X ↔ Z, Y → −Y, I → I, swapping " +
        "bit_a ↔ bit_b. A bit_b-homogeneous c for F112-Z thus transforms to a bit_a-homogeneous " +
        "c' for F112-X under simultaneous Q_zx-conjugation, with Frobenius norms preserved.";

    // ============================================================
    // Static helpers (delegating to PauliLetter.BitA())
    // ============================================================

    /// <summary>Returns the bit_a parity of <paramref name="term"/>: sum bit_a across
    /// the term's letters mod 2. By F38 / F63 Π_X² eigenvalue formula on Pauli strings,
    /// this is the Π_X²-conjugation eigenvalue exponent:
    /// Π_X² · σ · Π_X⁻² = (−1)^{BitAParity(σ)} · σ.</summary>
    public static int BitAParity(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        int sum = 0;
        foreach (var letter in term.Letters) sum += letter.BitA();
        return sum & 1;
    }

    /// <summary>Returns true iff every term in <paramref name="cTerms"/> shares the same
    /// <see cref="BitAParity"/> value. This is F112-X's structural hypothesis on the
    /// dissipator operators c_k: bit_a-homogeneous c implies, via Step 2 of the proof,
    /// that np.kron(c, c^*) lies entirely in the Π_X²-conjugation +1 eigenspace.
    ///
    /// <para>Edge cases: an empty list is vacuously homogeneous (returns true); a
    /// single-term list is trivially homogeneous (returns true); a null list throws.
    /// A term whose letters list is null or empty returns false (rejected as ill-formed).</para></summary>
    public static bool IsBitAHomogeneous(IReadOnlyList<PauliTerm> cTerms)
    {
        if (cTerms is null) throw new ArgumentNullException(nameof(cTerms));
        if (cTerms.Count == 0) return true;
        int? reference = null;
        foreach (var term in cTerms)
        {
            if (term is null) return false;
            if (term.Letters is null || term.Letters.Count == 0) return false;
            int a = BitAParity(term);
            if (reference is null) reference = a;
            else if (reference.Value != a) return false;
        }
        return true;
    }

    public LindbladBitAPiBalance(F108Part2Pi2XEvenAlwaysPalindromic part2)
        : base("F112-X Lindblad Π_X-eigenvalue balance under bit_a-homogeneous c: " +
               "‖M_plus_half‖² = ‖M_minus_half‖² for any H (Hermitian or non-Hermitian) and " +
               "bit_a-homogeneous c_k under X-dephase Π_X polarity. Tier1Derived universal N " +
               "for both Hermitian H (via the parent 5-step proof with axis_d = bit_a substituted) " +
               "and non-Hermitian H (via the two-lemma Welle-11 structural proof with axis_d = bit_a). " +
               "Two independent derivation routes per Welle 13: direct axis re-run (Route 1) " +
               "and Hadamard transport via Q_zx (Route 2). The bit_a-axis sister of F112-Z / F112-Y.",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F112 + " +
               "docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md + " +
               "docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md + " +
               "docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md + " +
               "docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md + " +
               "compute/RCPsiSquared.Diagnostics/Polarity/PolarityCoordinates.cs + " +
               "simulations/_f112_klein_v4_cross_dephase_verify.py")
    {
        Part2 = part2 ?? throw new ArgumentNullException(nameof(part2));
    }

    public override string DisplayName =>
        "F112-X Lindblad Π_X-eigenvalue balance under bit_a-homogeneous c " +
        "(Tier1Derived universal N for both Hermitian and non-Hermitian H; BitA sister of F112-Z / F112-Y)";

    public override string Summary =>
        $"{Theorem} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("5-step proof structure (axis_d = bit_a substitution)",
                summary: "Step 1: reduce balance to ‖M_{+i_X}‖² = ‖M_{-i_X}‖² via Π_X-eigenspace " +
                         "decomposition of M_plus_half / M_minus_half. " +
                         "Step 2: bit_a-homogeneous c implies np.kron(c, c.conj()) lies entirely in Π_X²-conj +1 " +
                         "eigenspace (via F38 / F63 Π_X² eigenvalue formula on Pauli strings; (−1)^bit_a(σ)). " +
                         "Step 3: Π_X²-conj +1 eigenspace = Π_X-conj {+1, −1}, hence dissipator has zero +i_X, −i_X content. " +
                         "Step 4: M_{+i_X} and M_{-i_X} come entirely from L_H = -i[H, ·] with norms 2 · ‖L_{H,±i_X}‖². " +
                         "Step 5 (Hermitian H): L_H^† = −L_H plus Π_X unitary implies dagger maps Π_X +i ↔ Π_X −i " +
                         "bijectively while preserving Frobenius; combining gives ‖L_{H,+i_X}‖² = ‖L_{H,-i_X}‖². ∎");
            yield return new InspectableNode("Derivation routes (Welle 13)",
                summary: DerivationRoutes);
            yield return new InspectableNode("Non-Hermitian H extension (Tier1Derived universal N)",
                summary: "Writing H = H_re + i H_im with both summands Hermitian, the equality reduces " +
                         "algebraically to Im⟨L_{H_re,-i_X}, L_{H_im,-i_X}⟩ = 0 for any Hermitian H_re, H_im. " +
                         "Closed via Welle-11 Lemma N-A^X (Diagonal-Norm: ‖L_{σ,-i_X}‖² = 4^N for any bit_a-odd σ) " +
                         "and Lemma N-B^X (Off-Diagonal-Orthogonality: ⟨L_{σ_α,-i_X}, L_{σ_β,-i_X}⟩ = 0 for " +
                         "σ_α ≠ σ_β both bit_a-odd), both reducing to F38 + Pauli-basis matrix-support disjointness.");
            yield return new InspectableNode("BitB sister F112-Z / F112-Y on opposite axis",
                summary: "LindbladBitBPiBalance covers the bit_b-axis case (Π_Z polarity with bit_b-homogeneous c). " +
                         "LindbladBitBPiYBalance covers the Π_Y polarity case (same bit_b axis, since Π_Y² also " +
                         "grades by bit_b per F38). F112-X (this Claim) covers the bit_a-axis case (Π_X polarity " +
                         "with bit_a-homogeneous c). The three Claims together close F112 across all dephase " +
                         "letters {Z, Y, X} and both Z₂-grading axes {bit_b, bit_a}.");
            yield return new InspectableNode("Empirical verification",
                summary: "Welle 13 verifier simulations/_f112_klein_v4_cross_dephase_verify.py confirms F112-X " +
                         "bit-exact at N = 2, 3: max|asymmetry| = 0.0 for Hermitian H + bit_a-homogeneous c " +
                         "and for non-Hermitian H + bit_a-homogeneous c (direct axis re-run); max|asymmetry| = " +
                         "3.6e-15 at N=2, 1.9e-30 at N=3 for the Hadamard-transport check (Z-config rotated → " +
                         "Π_X measurement). Lemma N-A^X / N-B^X stand at zero across the 32 bit_a-odd Pauli " +
                         "strings × 992 off-diagonal pairs at N=3.");
            yield return new InspectableNode("F108 Part 2 typed parent (shared bit_a + X-dephase foundation)",
                summary: "F108 Part 2 (Π²_X-even H + X-dephasing always palindromic) and F112-X share the bit_a " +
                         "Z₂-grading and the X-dephasing context. F108 Part 2's bilinear set {ZZ, XX, XY, YX, YY} " +
                         "is the bit_a = 0 family; F112-X's c-homogeneity hypothesis fixes c to one bit_a value " +
                         "(0 or 1). Both use the F38 Π_X² eigenvalue formula on Pauli strings as structural input.");
        }
    }
}
