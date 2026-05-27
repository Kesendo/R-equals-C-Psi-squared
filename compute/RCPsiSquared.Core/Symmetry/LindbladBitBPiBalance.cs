using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F112 (Tier1Derived for both Hermitian and non-Hermitian H, universal N): the structural identity behind the
/// `polarity_coordinates_from_L` diagnostic. For any Lindblad-form Liouvillian
///
/// <para>  L = -i[H, ·] + Σ_k γ_k · np.kron(c_k, c_k^*)</para>
///
/// with <b>Hermitian H</b> and each c_k bit_b-homogeneous (every Pauli string
/// σ in c_k's expansion shares bit_b(σ) = (#Y(σ) + #Z(σ)) mod 2 = const), the
/// `polarity_coordinates_from_L` decomposition of M = Π L Π⁻¹ + L + 2σ·I
/// satisfies
///
/// <para>  ‖M_plus_half‖² = ‖M_minus_half‖² bit-exactly.</para>
///
/// <para><b>Five-step proof structure</b> (Hermitian H, rigorous; see PROOF_F112):</para>
/// <list type="number">
///   <item>Π-eigenspace decomposition of M reduces the balance condition to
///         ‖M_{+i}‖² = ‖M_{-i}‖² (asymmetry = (1/2)(‖M_{+i}‖² − ‖M_{-i}‖²)).</item>
///   <item>bit_b-homogeneous c gives Π² · np.kron(c, c^*) · Π⁻² = ε² · np.kron(c, c^*)
///         with ε ∈ {+1, −1}; hence np.kron(c, c^*) lies entirely in the Π²-conj
///         +1 eigenspace, via the F38 / F63 Π² eigenvalue formula on Pauli strings.</item>
///   <item>Π²-conj +1 eigenspace = Π-conj {+1, −1}; the dissipator part of M has
///         zero +i and zero −i Π-conjugation content.</item>
///   <item>M_{+i} and M_{-i} come entirely from L_H = -i[H, ·], with norms
///         ‖M_{±i}‖² = 2 · ‖L_{H,±i}‖².</item>
///   <item>For Hermitian H: L_H^† = −L_H (anti-Hermitian as superoperator, Lemma B)
///         plus Π unitary implies dagger maps Π +i ↔ Π −i bijectively while preserving
///         Frobenius (Lemma A); combining gives ‖L_{H,+i}‖² = ‖L_{H,-i}‖². ∎</item>
/// </list>
///
/// <para><b>F87 orthogonality</b>: F87 (dissipator-resonance trichotomy) and F112 are
/// orthogonal axes on the shared bit_b Z₂-grading of the Pauli group. F87 lives in
/// spec(L) palindromy (Π conjugation acts on the spectrum); F112 lives in M_anti's
/// Π +i / −i split (Π conjugation acts on the operator). All three F87 classes
/// (truly, soft, hard) at N=3 under standard Z-dephasing give F112 asymmetry = 0
/// bit-exactly (`simulations/_polarity_probe_f87_connection.py`).</para>
///
/// <para><b>Non-Hermitian H extension (Tier1Derived, universal N)</b>: writing
/// H = H_re + i H_im with both summands Hermitian, the equality reduces algebraically
/// to the identity Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 for any Hermitian H_re, H_im.
/// This identity is Tier1Derived for all N via the two-lemma structural proof in
/// <c>docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md</c> (Welle 11, 2026-05-27):
/// Lemma N-A (Diagonal-Norm; the Welle-11 non-Hermitian-extension lemma, distinct from
/// the parent Hermitian-H Lemma A used in Step 5 above) gives ‖L_{σ,-i}‖² = 4^N for
/// any BitB-odd Pauli string σ via Π-eigenspace decomposition + matrix-support
/// disjointness; Lemma N-B (Off-Diagonal-Orthogonality; distinct from parent Lemma B)
/// gives ⟨L_{σ_α,-i}, L_{σ_β,-i}⟩ = 0 for σ_α ≠ σ_β
/// both BitB-odd via support disjointness across all four Π-orbit shifts. By
/// bilinearity + Pauli-basis spanning, F112 holds for any non-Hermitian H at any N.
/// The basis-enumeration anchor at N ≤ 5 (Welle 10a Python 559,912 pairs all bit-exact
/// 0; Welle 10b C# SLOW_F112 MaxImaginary &lt; 1e-10) is preserved as the empirical
/// motivation for the structural proof.</para>
///
/// <para><b>Diagnostic significance</b>: F112 makes the
/// `polarity_coordinates_from_L` asymmetry an exact witness; asymmetry ≠ 0 detects
/// c with cross-bit_b Pauli support, which sits OUTSIDE the F108 closure regime.
/// For any standard Lindblad system the diagnostic asymmetry is exactly 0 by
/// Steps 1-5.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitB"/>;
/// no BitA twin (the theorem is intrinsically about bit_b homogeneity, no meaningful
/// bit_a-axis analog exists), so
/// <see cref="BitATwinClassification.BitBSpecific"/>. Ctor parent
/// <see cref="F108Part1Pi2EvenAlwaysPalindromic"/> records the shared bit_b
/// foundation: F112 uses the F38 / F63 Π² eigenvalue formula on Pauli strings in
/// exactly the same way F108 does, and the bilinear set {XX, YY, YZ, ZY, ZZ}
/// that F108 palindromizes is exactly the bit_b = 0 (Π²-Z-even) family that lives
/// on the same Z₂-grading as F112's c constraint. Proof:
/// <c>docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md</c>.</para></summary>
public sealed class LindbladBitBPiBalance : Claim, IZ2AxisClaim
{
    /// <summary>BitB axis (bit_b parity Z₂-grading of the Pauli group; same axis
    /// as F108 Part 1 and F108 Part 3, anchored in F38 / F63).</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>No BitA twin: F112 is intrinsically a bit_b-axis theorem. The
    /// bit_b homogeneity of the c operators is the structural hypothesis; the
    /// bit_a parity (X+Y count) has no role in the proof or the conclusion, so
    /// no meaningful bit_a-axis analog exists.</summary>
    public Claim? BitATwin => null;

    /// <summary>Override returning <see cref="BitATwinClassification.BitBSpecific"/>:
    /// the algebraic content (bit_b homogeneity of c, Π²-conjugation +1 eigenspace
    /// containment, dagger anti-Hermiticity of L_H for Hermitian H) is intrinsically
    /// tied to bit_b structure, no meaningful bit_a-axis twin exists.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.BitBSpecific;

    /// <summary>Typed parent (F108 Part 1): records the shared bit_b foundation in
    /// the inheritance graph and surfaces in the inspector. F112 uses the F38 / F63
    /// Π² eigenvalue formula on Pauli strings exactly as F108 does, and F108's
    /// Π²-even bilinear set {XX, YY, YZ, ZY, ZZ} is the bit_b = 0 family on the same
    /// Z₂-grading.</summary>
    public F108Part1Pi2EvenAlwaysPalindromic Part1 { get; }

    /// <summary>The theorem statement in one line (Hermitian H scope).</summary>
    public string Theorem =>
        "For any Lindblad-form Liouvillian L = -i[H, ·] + Σ_k γ_k · np.kron(c_k, c_k^*) with " +
        "Hermitian H and each c_k bit_b-homogeneous (every Pauli string σ in c_k's expansion " +
        "shares bit_b(σ) = (#Y(σ) + #Z(σ)) mod 2 = const), the polarity_coordinates_from_L " +
        "decomposition of M = Π L Π⁻¹ + L + 2σ·I satisfies ‖M_plus_half‖² = ‖M_minus_half‖² " +
        "bit-exactly.";

    /// <summary>F87 ↔ F112 orthogonality on the shared bit_b Z₂-grading: F87 acts on
    /// spec(L), F112 acts on M_anti's Π-eigenspace decomposition.</summary>
    public string F87Orthogonality =>
        "F87 (dissipator-resonance trichotomy) and F112 are orthogonal axes on the shared " +
        "bit_b Z₂-grading of the Pauli group. F87 lives in spec(L) palindromy (Π conjugation " +
        "acts on the spectrum); F112 lives in M_anti's Π +i / −i split (Π conjugation acts on " +
        "the operator). All three F87 classes (truly, soft, hard) at N=3 under standard " +
        "Z-dephasing give F112 asymmetry = 0 bit-exactly.";

    /// <summary>Non-Hermitian H extension status: Tier1Derived universal N (Welle 11,
    /// 2026-05-27) via two-lemma structural proof in
    /// docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md.</summary>
    public string NonHermitianExtension =>
        "Writing H = H_re + i H_im with both summands Hermitian, the equality reduces " +
        "algebraically to the identity F(H_re, H_im) := Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 " +
        "for any Hermitian H_re, H_im. F is real-bilinear and antisymmetric under H_re ↔ H_im " +
        "exchange, so it is determined by its values on Pauli-string basis pairs. The " +
        "per-pair identity F(σ_α, σ_β) = 0 holds structurally at any N via two lemmas: " +
        "Lemma N-A (Diagonal-Norm; the Welle-11 non-Hermitian-extension lemma, distinct " +
        "from the parent Hermitian-H Lemma A used in Step 5 of the parent proof) gives " +
        "‖L_{σ,-i}‖² = 4^N for any BitB-odd σ; Lemma N-B " +
        "(Off-Diagonal-Orthogonality; distinct from parent Lemma B) gives " +
        "⟨L_{σ_α,-i}, L_{σ_β,-i}⟩ = 0 for σ_α ≠ σ_β " +
        "both BitB-odd. Both Welle-11 lemmas reduce to Pauli-basis matrix-support disjointness of " +
        "L_σ and Π^m L_σ Π^{-m} under the F38 / F63 Π² eigenvalue formula. F112 " +
        "non-Hermitian extension is therefore Tier1Derived for ALL N (proof in " +
        "docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md, Welle 11 2026-05-27). The " +
        "basis-enumeration anchor at N ≤ 5 (559,912 pairs all bit-exact 0; Welle 10a " +
        "Python simulations/_f112_open_identity_basis_enum.py and Welle 10b C# SLOW_F112 " +
        "compute/RCPsiSquared.Diagnostics/Polarity/F112NonHermitianBasisEnumeration.cs) " +
        "is preserved as the empirical motivation for the structural proof. The two-lemma " +
        "proof verifier simulations/_f112_universal_n_proof_verify.py confirms each step " +
        "within 1e-12 numpy double-precision tolerance at N = 1, 2, 3 (42 BitB-odd strings, " +
        "1050 off-diagonal pairs, 4368 all-pair F-values, all max deviations < 1e-12, " +
        "i.e. machine zero to numpy double precision). The Welle 10a Python enumeration " +
        "above is genuinely bit-exact at N ≤ 4 (rational matrix entries); the Welle 11 " +
        "verifier is numerical. See experiments/F112_NONHERMITIAN_BASIS_ENUMERATION.md.";

    // ============================================================
    // Static helpers (delegating to PauliLetter.BitB())
    // ============================================================

    /// <summary>Returns the bit_b parity of <paramref name="term"/>: sum bit_b across
    /// the term's letters mod 2. By the F38 / F63 Π² eigenvalue formula on Pauli
    /// strings, this is the Π²-Z conjugation eigenvalue exponent: Π² · σ · Π⁻² =
    /// (−1)^{BitBParity(σ)} · σ.</summary>
    public static int BitBParity(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        int sum = 0;
        foreach (var letter in term.Letters) sum += letter.BitB();
        return sum & 1;
    }

    /// <summary>Returns true iff every term in <paramref name="cTerms"/> shares the
    /// same <see cref="BitBParity"/> value. This is F112's structural hypothesis on
    /// the dissipator operators c_k: bit_b-homogeneous c implies, via Step 2 of the
    /// proof, that np.kron(c, c^*) lies entirely in the Π²-conjugation +1
    /// eigenspace.
    ///
    /// <para>Edge cases: an empty list is vacuously homogeneous (returns true);
    /// a single-term list is trivially homogeneous (returns true); a null list
    /// throws. A term whose letters list is null or empty returns false (rejected
    /// as ill-formed).</para></summary>
    public static bool IsBitBHomogeneous(IReadOnlyList<PauliTerm> cTerms)
    {
        if (cTerms is null) throw new ArgumentNullException(nameof(cTerms));
        if (cTerms.Count == 0) return true;
        int? reference = null;
        foreach (var term in cTerms)
        {
            if (term is null) return false;
            if (term.Letters is null || term.Letters.Count == 0) return false;
            int b = BitBParity(term);
            if (reference is null) reference = b;
            else if (reference.Value != b) return false;
        }
        return true;
    }

    /// <summary>Documents the F112 scope assertion: the Hamiltonian must be
    /// Hermitian for the Tier1Derived proof (Step 5 uses L_H^† = −L_H, which holds
    /// only when H = H^†). This method always returns true because any real-coefficient
    /// linear combination of Pauli strings IS Hermitian (Pauli strings are Hermitian
    /// matrices, and a real-coefficient Hermitian-summand sum stays Hermitian).
    /// The method exists to make the scope assertion explicit in code: when a
    /// caller routes a real-coefficient Hamiltonian into F112's pipeline, this is
    /// the documentation hook saying "yes, the Hermitian-H scope applies".
    ///
    /// <para>For non-real coefficients (e.g. complex weights in a Pauli decomposition),
    /// the Hermitian-H proof's Step 5 does not directly apply; that case is covered
    /// by the non-Hermitian extension, now Tier1Derived universal N via the two-lemma
    /// structural proof (Welle 11, 2026-05-27). See
    /// <see cref="NonHermitianExtension"/>.</para></summary>
    public static bool IsHermitianHamiltonian(IReadOnlyList<double> coefficients)
    {
        if (coefficients is null) throw new ArgumentNullException(nameof(coefficients));
        return true;
    }

    public LindbladBitBPiBalance(F108Part1Pi2EvenAlwaysPalindromic part1)
        : base("F112 Lindblad Π-eigenvalue balance under bit_b-homogeneous c: " +
               "‖M_plus_half‖² = ‖M_minus_half‖² for any H (Hermitian or non-Hermitian) and " +
               "bit_b-homogeneous c_k. Tier1Derived universal N for both Hermitian H (via the " +
               "parent 5-step proof) and non-Hermitian H (via the two-lemma structural proof " +
               "in docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md, Welle 11). The " +
               "structural identity behind polarity_coordinates_from_L; the diagnostic " +
               "asymmetry is an exact witness for c outside the bit_b-homogeneous regime.",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F112 + " +
               "docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md + " +
               "docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md + " +
               "simulations/framework/diagnostics/polarity_coordinates.py + " +
               "simulations/_polarity_proof_verify.py + simulations/_polarity_step5_stress.py + " +
               "simulations/_polarity_probe_f87_connection.py + " +
               "simulations/_f112_universal_n_proof_verify.py")
    {
        Part1 = part1 ?? throw new ArgumentNullException(nameof(part1));
    }

    public override string DisplayName =>
        "F112 Lindblad Π-eigenvalue balance under bit_b-homogeneous c " +
        "(Tier1Derived universal N for both Hermitian and non-Hermitian H)";

    public override string Summary =>
        $"{Theorem} {F87Orthogonality} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("5-step proof structure",
                summary: "Step 1: reduce balance to ‖M_{+i}‖² = ‖M_{-i}‖² via Π-eigenspace " +
                         "decomposition of M_plus_half / M_minus_half (asymmetry = (1/2)(‖M_{+i}‖² − ‖M_{-i}‖²)). " +
                         "Step 2: bit_b-homogeneous c implies np.kron(c, c.conj()) lies entirely in Π²-conj +1 " +
                         "eigenspace (via F38 / F63 Π² eigenvalue formula on Pauli strings; ε² = 1 for ε ∈ {+1, −1}). " +
                         "Step 3: Π²-conj +1 eigenspace = Π-conj {+1, −1}, hence dissipator has zero +i, −i content. " +
                         "Step 4: M_{+i} and M_{-i} come entirely from L_H = -i[H, ·] with norms 2 · ‖L_{H,±i}‖². " +
                         "Step 5 (Hermitian H): L_H^† = −L_H (Lemma B, anti-Hermitian superoperator), and dagger " +
                         "maps Π +i ↔ Π −i bijectively while preserving Frobenius (Lemma A); combining gives " +
                         "‖L_{H,+i}‖² = ‖L_{H,-i}‖². ∎");
            yield return new InspectableNode("Empirical anchor",
                summary: "14 probes in simulations/_polarity_probe_*.py, _polarity_proof_verify.py, " +
                         "_polarity_step5_stress.py: candidate-breakers (1-5), hand-engineered non-Lindblad " +
                         "L (6), random c with full Pauli rank (7-8), k_max boundary (9), exhaustive 136-pair " +
                         "N=2 enumeration (10), coefficient sweep (11), Z₂³-cell N=3, 4 scaling (12, 171 / 171 " +
                         "balanced within-cell), Π²-content verification (13, 100.00% Π²=+1 for bit_b-homogeneous " +
                         "c), and direct Π-eigenspace L_H projection across 30 random H (10 Hermitian + 10 " +
                         "non-Hermitian Pauli + 10 random complex matrix) at N=2, 3 (14, all bit-exact).");
            yield return new InspectableNode("F87 orthogonality", summary: F87Orthogonality);
            yield return new InspectableNode("Non-Hermitian H extension (Tier1Derived universal N)",
                summary: NonHermitianExtension +
                         " Reduced identity: Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 for any Hermitian H_re, H_im.");
            yield return new InspectableNode("BitA twin (BitBSpecific)",
                summary: "No BitA twin: F112 is intrinsically a bit_b-axis theorem. The bit_b homogeneity of " +
                         "the c operators is the structural hypothesis; the bit_a parity (X+Y count) has no role " +
                         "in the proof or the conclusion. BitATwinStatus = BitBSpecific.");
            yield return new InspectableNode("Shared bit_b axis with F108 + F87 + F38",
                summary: "F112 and F108 Part 1 / Part 3 share the bit_b Z₂-grading of the Pauli group. F108's " +
                         "Π²-even bilinear set {XX, YY, YZ, ZY, ZZ} is the bit_b = 0 family; F112's c-homogeneity " +
                         "hypothesis fixes c to one bit_b value. Both use the F38 Π² eigenvalue formula (Π² · σ · Π⁻² = " +
                         "(−1)^{BitBParity(σ)} · σ) and F63's [L, Π²] = 0 for Z-dephasing as structural inputs. " +
                         "F87's dissipator-resonance trichotomy is the spectrum-side reading of the same grading; " +
                         "F112 is the M_anti Π-eigenspace-side reading. Three projections of one Z₂-grading.");
        }
    }
}
