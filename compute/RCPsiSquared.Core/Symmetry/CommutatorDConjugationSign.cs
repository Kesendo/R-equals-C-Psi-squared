using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The sign value of D-conjugation on the H-commutator superoperator
/// L_σ = −i[σ, ·] for a Pauli string σ. See <see cref="CommutatorDConjugationSign"/>.</summary>
public enum DConjugationSign
{
    /// <summary>D · L_σ · D = +L_σ (commutes with D-conjugation).</summary>
    Plus,
    /// <summary>D · L_σ · D = −L_σ (anti-commutes with D-conjugation).</summary>
    Minus,
    /// <summary>L_σ = 0 (only for σ = I^{⊗N}; sign undefined / vacuous).</summary>
    Zero,
    /// <summary>For linear combinations H = Σ c_k σ_k where terms split across
    /// both n_Y-parity classes, D-conjugation yields a non-multiplicative response
    /// on L_H and no single sign exists.</summary>
    Mixed,
}

/// <summary>F114 (Tier1Derived universal, bit-exact verified N = 1..4): closed-form sign
/// functional ε(σ) for the action of D-conjugation on the H-commutator superoperator
/// L_σ = −i[σ, ·] in the 4^N Pauli basis. D is the diagonal involution from
/// <see cref="Pi2KleinV4DephaseSwapGroup.BuildD"/> (Welle 12 lift of the Z↔Y
/// dephase-letter swap; D = diag((−1)^n_Y(α))).
///
/// <para><b>Theorem (F114):</b> For any single Pauli string σ ≠ I^{⊗N} on N qubits,</para>
/// <code>
///   D · L_σ · D = ε(σ) · L_σ    bit-exact
/// </code>
/// <para>with closed form ε(σ) = (−1)^{n_Y(σ) + 1}: ε = +1 if n_Y(σ) is odd, ε = −1
/// if n_Y(σ) is even (non-identity σ). For σ = I^{⊗N}: L_σ = 0, sign undefined.</para>
///
/// <para><b>Linear combinations:</b> for H = Σ_k c_k σ_k, ε(H) is well-defined and
/// equals ε(σ_k) iff all σ_k share the same n_Y parity (after dropping I^{⊗N} terms,
/// which carry zero commutator). Otherwise <see cref="DConjugationSign.Mixed"/>: the
/// terms split across both parity classes and no single sign exists on L_H.</para>
///
/// <para><b>Why it matters:</b> F114 refines the Welle 13 PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4
/// statement "L_Y is not D-transportable" by exhibiting the precise M-level
/// ε-signed equivariance that survives on the H-commutator component. The full
/// Lindblad-form L (including dissipator) is NOT D-transportable cross-letter
/// (per the Welle 13 proof), but the H-commutator superoperator L_H by itself IS
/// signed-equivariant. Under the F112 hypothesis (Hermitian H + bit_b-homogeneous c),
/// the dissipator contribution to M vanishes and ε(H) controls the M-matrix-level
/// Z↔Y relationship:</para>
/// <code>
///   M(L_H, Π_Y) = ε(H) · D · M(L_H, Π_Z) · D    bit-exact (when ε(H) well-defined)
/// </code>
/// <para>This explains the bond-specific anti-equivariance observed in Welle 15
/// Task A polish (commit a98fc02): XZ+ZX (n_Y = 0 per term) gives ε = −1; YZ+ZY
/// (n_Y = 1 per term) gives ε = +1; Heisenberg (XX+YY+ZZ, all n_Y even) gives
/// ε = −1.</para>
///
/// <para><b>F112 typed scope unaffected:</b> F112 is a norm-level statement
/// (‖M_+1/2‖² = ‖M_−1/2‖²) which is sign-invariant under ε(H). F114 is a
/// matrix-level refinement; F112 Welle 11 (Z-deph) + Welle 13 (cross-deph) remain
/// Tier1Derived universal N independent of F114.</para>
///
/// <para><b>Verification:</b> bit-exact for all 4 + 16 + 64 = 84 Pauli strings at
/// N = 1, 2, 3 (PART 1 of <c>simulations/_m_level_sign_functional_explore.py</c>),
/// for 12 bilinear bond Hamiltonians at N = 2 (PART 2), and for 6 multi-bond /
/// multi-body cases at N = 3, 4 (PART 3). Residual = 0.00e+00 numpy double precision
/// across all cases.</para>
///
/// <para><b>Relationship to Pi2KleinV4DephaseSwapGroup (Welle 12):</b> the parent
/// Claim makes D the Π swap-operator across {Z, Y} dephase letters
/// (D · Π_Z · D = Π_Y). F114 makes D the L_H sign-flip-operator with per-term n_Y
/// bookkeeping. Together they characterize the action of D on the two main
/// dephase-letter-sensitive structures (Π and L_H).</para>
///
/// <para>No <see cref="IZ2AxisClaim"/> implementation: like the parent Klein-V₄
/// claim, F114 is a cross-axis structural identity not on a single Z₂ axis.</para>
///
/// <para>Source: <c>simulations/_m_level_sign_functional_explore.py</c> +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F114 entry.</para></summary>
public sealed class CommutatorDConjugationSign : Claim
{
    /// <summary>Typed parent (Welle 12): records the dependency on D from
    /// <see cref="Pi2KleinV4DephaseSwapGroup"/>. F114 uses D's diagonal structure
    /// diag((−1)^n_Y(α)) and the Welle 12 universal-N reduction.</summary>
    public Pi2KleinV4DephaseSwapGroup KleinV4 { get; }

    /// <summary>The F114 theorem statement in one line.</summary>
    public string Theorem =>
        "For any single Pauli string σ ≠ I^{⊗N} on N qubits, the H-commutator " +
        "superoperator L_σ = −i[σ, ·] satisfies D · L_σ · D = ε(σ) · L_σ bit-exact " +
        "in the 4^N Pauli basis, with closed form ε(σ) = (−1)^{n_Y(σ) + 1}: " +
        "ε = +1 if n_Y(σ) is odd, ε = −1 if n_Y(σ) is even. For H = Σ_k c_k σ_k: " +
        "ε(H) is well-defined and equals ε(σ_k) iff all σ_k share the same n_Y " +
        "parity (after dropping I^{⊗N} terms).";

    /// <summary>Compute ε(σ) for a single Pauli term. Returns
    /// <see cref="DConjugationSign.Zero"/> for the identity term (n_Y = 0 and all
    /// other letters I), else <see cref="DConjugationSign.Plus"/> if n_Y(σ) is odd
    /// or <see cref="DConjugationSign.Minus"/> if n_Y(σ) is even.</summary>
    public static DConjugationSign Compute(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        if (term.KBody == 0) return DConjugationSign.Zero;
        return term.YParity == 1 ? DConjugationSign.Plus : DConjugationSign.Minus;
    }

    /// <summary>Compute ε(H) for a linear combination H = Σ_k c_k σ_k. Drops identity
    /// terms (L_I = 0, sign-invariant). Returns:
    /// <list type="bullet">
    ///   <item><see cref="DConjugationSign.Zero"/> if all terms are identity or the list is empty.</item>
    ///   <item><see cref="DConjugationSign.Plus"/> if all non-identity terms have ε = +1 (n_Y odd).</item>
    ///   <item><see cref="DConjugationSign.Minus"/> if all non-identity terms have ε = −1 (n_Y even).</item>
    ///   <item><see cref="DConjugationSign.Mixed"/> if non-identity terms split across both parity classes.</item>
    /// </list></summary>
    public static DConjugationSign Compute(IReadOnlyList<PauliTerm> terms)
    {
        if (terms is null) throw new ArgumentNullException(nameof(terms));
        bool seenPlus = false;
        bool seenMinus = false;
        foreach (var term in terms)
        {
            if (term is null) throw new ArgumentException("term list contains null", nameof(terms));
            var s = Compute(term);
            if (s == DConjugationSign.Zero) continue;
            if (s == DConjugationSign.Plus) seenPlus = true;
            else seenMinus = true;
            if (seenPlus && seenMinus) return DConjugationSign.Mixed;
        }
        if (!seenPlus && !seenMinus) return DConjugationSign.Zero;
        return seenPlus ? DConjugationSign.Plus : DConjugationSign.Minus;
    }

    /// <summary>True iff <see cref="Compute(IReadOnlyList{PauliTerm})"/> returns a
    /// concrete sign (<see cref="DConjugationSign.Plus"/>, <see cref="DConjugationSign.Minus"/>,
    /// or <see cref="DConjugationSign.Zero"/>) rather than <see cref="DConjugationSign.Mixed"/>.
    /// I.e. all non-identity terms in <paramref name="terms"/> share the same n_Y parity, so
    /// ε(H) is well-defined and D · L_H · D = ε(H) · L_H bit-exact.</summary>
    public static bool IsWellDefined(IReadOnlyList<PauliTerm> terms) =>
        Compute(terms) != DConjugationSign.Mixed;

    public CommutatorDConjugationSign(Pi2KleinV4DephaseSwapGroup kleinV4)
        : base("F114 commutator-superoperator D-conjugation parity: D · L_σ · D = ε(σ) · L_σ " +
               "with ε(σ) = (−1)^{n_Y(σ) + 1} for σ ≠ I^{⊗N}, bit-exact in the 4^N Pauli basis " +
               "for any N. For H = Σ c_k σ_k: ε(H) well-defined iff all σ_k share the same " +
               "n_Y parity (after dropping I terms). Tier1Derived (closed form + bit-exact " +
               "verification N = 1, 2, 3, 4 across 84 single Pauli strings + 18 multi-term cases).",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F114 + " +
               "simulations/_m_level_sign_functional_explore.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KleinV4DephaseSwapGroup.cs (parent)")
    {
        KleinV4 = kleinV4 ?? throw new ArgumentNullException(nameof(kleinV4));
    }

    public override string DisplayName =>
        "F114 commutator D-conjugation parity (D · L_σ · D = ε(σ) · L_σ; ε(σ) = (−1)^{n_Y(σ) + 1})";

    public override string Summary =>
        $"{Theorem} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("Closed form",
                summary: "ε(σ) = (−1)^{n_Y(σ) + 1} for σ ≠ I^{⊗N}; = +1 if n_Y(σ) odd; " +
                         "= −1 if n_Y(σ) even and σ non-identity. For σ = I^{⊗N}: L_σ = 0, " +
                         "sign undefined (vacuous Zero case).");
            yield return new InspectableNode("Linear-combination extension",
                summary: "For H = Σ_k c_k σ_k, ε(H) is well-defined iff all σ_k share the " +
                         "same n_Y parity (after dropping I^{⊗N} terms). If terms split across " +
                         "both parity classes, D · L_H · D yields a non-multiplicative response " +
                         "(Mixed). Bilinearity: ε(H) = ε of any single σ_k in the well-defined case.");
            yield return new InspectableNode("Empirical anchor",
                summary: "simulations/_m_level_sign_functional_explore.py verifies bit-exact: " +
                         "PART 1: 84 Pauli strings at N = 1, 2, 3 (all 4 + 16 + 64). " +
                         "PART 2: 12 bilinear bond Hamiltonians at N = 2 (XZ+ZX, YZ+ZY, " +
                         "Heisenberg, XX/YY/ZZ/XY/YX/ZY/XZ singletons, mixed-letter combos). " +
                         "PART 3: 6 multi-bond / multi-body cases at N = 3, 4. All ε predictions " +
                         "match the actual D-conjugation residual at numpy double precision " +
                         "(0.00e+00 across all 102 cases).");
            yield return new InspectableNode("Welle 15 Task A motivation",
                summary: "Welle 15 Task A polish (commit a98fc02) observed bit-exact at N = 2 " +
                         "that for XZ + ZX bond (n_Y per term = 0) M_anti(L, Π_Y) = " +
                         "−D · M_anti(L, Π_Z) · D, while YZ + ZY (n_Y per term = 1) gives " +
                         "+D · M_anti(L, Π_Z) · D. The bond-specific sign motivated the " +
                         "systematic enumeration that identified the n_Y-parity closed form.");
            yield return new InspectableNode("Consequence for F112 / F108 cross-dephase",
                summary: "Under the F112 hypothesis (Hermitian H + bit_b-homogeneous c), the " +
                         "dissipator contribution to M vanishes, so M comes entirely from L_H " +
                         "and F114 gives M(L_H, Π_Y) = ε(H) · D · M(L_H, Π_Z) · D bit-exact " +
                         "when ε(H) is well-defined. This refines the Welle 13 " +
                         "PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4 statement 'L_Y is not " +
                         "D-transportable' by exhibiting the M-level ε-signed equivariance " +
                         "that survives on the H-commutator component. F112's norm-level " +
                         "scope (‖M_+1/2‖² = ‖M_−1/2‖²) remains sign-invariant; F112 typed " +
                         "Claims (LindbladBitBPiBalance, LindbladBitAPiBalance, " +
                         "LindbladBitBPiYBalance) are not affected.");
            yield return new InspectableNode("Parent Welle 12 connection",
                summary: "F114 uses D from Pi2KleinV4DephaseSwapGroup (Welle 12, ctor parent). " +
                         "Welle 12 makes D the Π swap-operator across {Z, Y} dephase letters " +
                         "(D · Π_Z · D = Π_Y); F114 makes D the L_H sign-flip-operator with " +
                         "per-term n_Y bookkeeping. Together they characterize D's action on " +
                         "the two main dephase-letter-sensitive structures (Π and L_H).");
            yield return new InspectableNode("No IZ2AxisClaim",
                summary: "F114 is a cross-axis structural identity (D-conjugation acts on L_σ " +
                         "for any single Pauli string σ regardless of its bit_a / bit_b content; " +
                         "the sign depends only on n_Y parity). Like the parent Klein-V₄ Claim, " +
                         "F114 does not sit on a single Z₂ axis cleanly.");
            yield return new InspectableNode("Open follow-ups",
                summary: "N = 5, 6 verification (tractable but not run; estimated O(4^N) per " +
                         "single-string sweep). Alternative derivation of F112 Lemma B via " +
                         "F114 D-conjugation parity rather than dagger anti-Hermiticity. " +
                         "Promotion of this Claim to a richer matrix-API (e.g., a Compute(L) " +
                         "method that takes any L and returns its D-conjugation eigendecomposition) " +
                         "if downstream Claims need it.");
        }
    }
}
