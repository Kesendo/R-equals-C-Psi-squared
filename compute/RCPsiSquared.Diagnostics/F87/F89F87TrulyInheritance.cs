using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F89 ↔ F87 bridge: the F89 AT-lock Re(λ_n) = −2γ₀ on overlap-subspace
/// F_a modes is a structural consequence of F87 classifying F89's bond Hamiltonian
/// H_B = J·Σ_{(p,q)∈B}(X_pX_q + Y_pY_q) as <see cref="TrichotomyClass.Truly"/>.
///
/// <para>The inheritance chain:</para>
/// <list type="number">
///   <item>F89's per-bond term is the canonical XX+YY Pauli pair — the first
///   <see cref="F87CanonicalWitness"/> in <see cref="F87CanonicalWitness.StandardSet"/>,
///   expected and verified <see cref="TrichotomyClass.Truly"/>.</item>
///   <item>F87-Truly is the operator equation <c>Π·L·Π⁻¹ + L + 2σ·I = 0</c>
///   holding bit-exactly on the full Liouvillian (<see cref="PauliPairTrichotomy"/>
///   palindrome residual M = 0). This forces the spectrum-pairing
///   <c>λ + (−λ − 2σ) = −2σ</c> across the entire spectrum.</item>
///   <item>F87-Truly together with single-letter Z-dephasing γ_l = γ₀ on every site
///   makes the Absorption Theorem hold bit-exactly: <c>Re(λ) = −2γ₀·⟨n_XY⟩</c>
///   for every Liouvillian eigenmode (cf.
///   <see cref="Core.Symmetry.AbsorptionTheoremClaim.HammingComplementPairSum"/>).</item>
///   <item>F89's F_a modes live in the overlap subspace of the (SE, DE) sub-block —
///   coherences |i⟩⟨{j, l_pair}| where the SE index i coincides with one DE index.
///   These coherences have Hamming distance <c>n_diff = 1</c> in the computational
///   basis (one site differs in the bit pattern). By point 3, Re(λ_n) = −2γ₀·1 = −2γ₀
///   exactly — the F89 AT-lock as documented in
///   <see cref="Core.Symmetry.F89PathKAtLockMechanismClaim"/>.</item>
/// </list>
///
/// <para>Equivalently: the F89 AT-lock formula degrades to "Re(λ) ≈ −2γ₀ up to
/// O(‖M‖_F) corrections" when the bond Hamiltonian is not F87-Truly. The bit-exact
/// match is the F87-Truly precondition realized on F89's specific sub-block.</para>
///
/// <para><b>Sibling bond-classes that preserve the F89 AT-lock</b>: any F87-Truly
/// extension that keeps Re(λ) = −2γ₀·n_diff exact also lifts the F89 AT-lock —
/// e.g. Heisenberg H = J·Σ(XX+YY+ZZ) (the second
/// <see cref="F87CanonicalWitness.StandardSet"/> entry, also Truly), since ZZ is
/// Π²-even-truly under Z-dephasing.</para>
///
/// <para><b>What breaks the AT-lock</b>: any bond term that lands in F87-Soft or
/// F87-Hard. For example, replacing one bond's XX+YY with YZ+ZY (F87-Soft per the
/// Marrakesh EQ-030 anchor) mixes the spectrum; F89's F_a modes no longer sit at
/// Re(λ) = −2γ₀ exactly. This is the F87-Soft amplification visible in the
/// hardware-confirmed soft-trichotomy drop on ibm_marrakesh
/// (<see cref="Core.Confirmations.ConfirmationsRegistry"/> entry
/// <c>palindrome_trichotomy</c>).</para>
///
/// <para>Anchors:
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F87 entry,
/// <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c> § "AT-lock",
/// <see cref="F87CanonicalWitness"/> XX+YY witness (Truly),
/// <see cref="Core.Symmetry.F89PathKAtLockMechanismClaim"/> (AT-lock typed claim),
/// <see cref="Core.Symmetry.AbsorptionTheoremClaim"/>.</para></summary>
public sealed class F89F87TrulyInheritance : Claim
{
    /// <summary>The two-bond Pauli-pair specification matching the F89 Hamiltonian
    /// H_b = J·(X_pX_q + Y_pY_q) per bond. Stable across N and across the bond
    /// graph B (orbit-closure invariance).</summary>
    public static readonly IReadOnlyList<PauliPairBondTerm> F89BondTerms = new[]
    {
        new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
        new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
    };

    // Parent-edge markers for Schicht-1 wiring.
    private readonly F87TrichotomyClassification _trichotomy;

    /// <summary>Classify the F89 bond Hamiltonian via F87 on a given chain. Returns
    /// <see cref="TrichotomyClass.Truly"/> for every well-formed N ≥ 2 and uniform γ
    /// (this is the F89 AT-lock precondition); deviation from Truly would indicate
    /// either a chain configuration that breaks F87's classifier or a code regression.</summary>
    public static TrichotomyClass ClassifyF89Hamiltonian(ChainSystem chain) =>
        PauliPairTrichotomy.Classify(chain, F89BondTerms);

    /// <summary>True when F87 classifies F89's bond Hamiltonian as Truly on the
    /// given chain — i.e., the AT-lock precondition is satisfied and Re(λ_n) = −2γ₀
    /// holds exactly on F_a modes.</summary>
    public static bool IsAtLockPreconditionSatisfied(ChainSystem chain) =>
        ClassifyF89Hamiltonian(chain) == TrichotomyClass.Truly;

    public F89F87TrulyInheritance(F87TrichotomyClassification trichotomy)
        : base("F89 ↔ F87 bridge: F89 bond Hamiltonian H_b = J·(XX+YY) is F87-Truly; the AT-lock Re(λ_n) = −2γ₀ on F_a modes is the n_diff=1 instance of F87-Truly's bit-exact Absorption Theorem (Π·L·Π⁻¹ + L + 2σ·I = 0)",
               Tier.Tier1Derived,
               "F87CanonicalWitness 'XX+YY' (Truly anchor); F89PathKAtLockMechanismClaim (Re(λ) = −2γ₀ on F_a); AbsorptionTheoremClaim.HammingComplementPairSum (per-coherence rate 2γ₀·n_diff); docs/ANALYTICAL_FORMULAS.md F87 + F89; docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md § AT-lock")
    {
        _trichotomy = trichotomy ?? throw new ArgumentNullException(nameof(trichotomy));
    }

    public override string DisplayName => "F89 AT-lock inherits F87-Truly (bridge)";

    public override string Summary =>
        $"F89 H_b = J·(XX+YY) classifies as F87-Truly; AT-lock Re(λ_n) = −2γ₀ on F_a modes " +
        $"is the n_diff=1 instance of F87-Truly's bit-exact AT ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F89 bond term",
                summary: "H_b = J·(X_pX_q + Y_pY_q) per bond (p, q) ∈ B");
            yield return new InspectableNode("F87 classification",
                summary: "TrichotomyClass.Truly (palindrome residual M = 0 bit-exactly; F87CanonicalWitness 'XX+YY' anchor)");
            yield return new InspectableNode("Inheritance chain",
                summary: "F87-Truly ⇒ Π·L·Π⁻¹ + L + 2σ·I = 0 bit-exact ⇒ AT exact: Re(λ) = −2γ₀·⟨n_XY⟩ for every eigenmode ⇒ F_a modes (n_diff=1) have Re(λ_n) = −2γ₀");
            yield return new InspectableNode("AT-lock formula",
                summary: "λ_n = −2γ₀ + i·y_n with y_n = 4·cos(πn/(N_block+1)); typed in F89PathKAtLockMechanismClaim");
            yield return new InspectableNode("Bond-class extension",
                summary: "Heisenberg H = J·(XX+YY+ZZ) is also F87-Truly under Z-dephasing (ZZ is Π²-even-truly), preserving the F89 AT-lock formula");
            yield return new InspectableNode("Break case",
                summary: "Replacing any bond's XX+YY with YZ+ZY (F87-Soft, EQ-030 Marrakesh anchor) mixes the spectrum; F_a modes no longer sit at Re(λ) = −2γ₀ exactly");
        }
    }
}
