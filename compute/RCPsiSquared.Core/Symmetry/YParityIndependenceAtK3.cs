using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F102 (Tier1Derived): Y-parity as a term-level Z₂ classifier that
/// collapses to <c>bit_a XOR bit_b</c> at k_body=2 but is independent of the Klein
/// (bit_a, bit_b) signature at k_body≥3.
///
/// <para>For any Pauli string σ = ⊗_l σ_α_l, define y_par(σ) = (Σ_l [σ_α_l = Y]) mod 2.
/// Per-letter contributions are (bit_a, bit_b, y_par) = (0,0,0) for I, (1,0,0) for X,
/// (1,1,1) for Y, (0,1,0) for Z; equivalently bit_a = n_X + n_Y, bit_b = n_Y + n_Z,
/// y_par = n_Y, all mod 2.</para>
///
/// <para>At k_body=2 (every 2-body Pauli bilinear XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ),
/// the algebraic identity y_par = bit_a XOR bit_b holds (mechanical check on the 9
/// bilinears). The identity in general holds iff k_body is even: at k_body=0 trivially
/// (0 = 0), at k_body=2 by the per-letter table; at odd k_body it strictly fails by
/// y_par = (k_body + (bit_a XOR bit_b)) mod 2.</para>
///
/// <para>At k_body=3 the canonical counterexample is σ = XYZ vs σ' = III:
/// XYZ has (bit_a, bit_b, y_par) = (1+1, 1+1, 1) mod 2 = (0, 0, 1);
/// III has (0, 0, 0). Both share Klein (0, 0) but differ on y_par, demonstrating
/// y_par is independent of the Klein signature at k_body≥3.</para>
///
/// <para>Per F34/QUBIT_NECESSITY there is no third independent Π²-operator (Π²_Y
/// collapses to Π²_Z at the operator level). The Y-parity classifier is term-level,
/// not operator-level; this Claim does NOT follow the Pi²-Inheritance naming
/// convention. Implements <see cref="IZ2AxisClaim"/> with
/// <see cref="Z2Axis.YParity"/>; <see cref="BitATwin"/> is always null (Y-axis Claims
/// do not carry a BitA-twin slot semantics).</para>
///
/// <para>Witness: <c>compute/RCPsiSquared.Core.Tests/Pauli/PauliHamiltonianKleinHelpersTests.cs</c>
/// (XYZ_AtK3_IsKleinHomogeneousButZ2HomogeneityRefinesViaYParity).
/// Proof: <c>docs/proofs/PROOF_F102_YPARITY_INDEPENDENCE.md</c>.
/// Documented in <c>docs/ANALYTICAL_FORMULAS.md</c> F102 and
/// <c>hypotheses/THE_POLARITY_LAYER.md</c> Z. 130-134, 313-320.</para></summary>
public sealed class YParityIndependenceAtK3 : Claim, IZ2AxisClaim
{
    /// <summary>YParity axis classification per <see cref="IZ2AxisClaim"/>.</summary>
    public Z2Axis Z2Axis => Z2Axis.YParity;

    /// <summary>YParity-axis Claim; no BitATwin slot semantics.</summary>
    public Claim? BitATwin => null;

    /// <summary>The k_body=2 collapse identity: for every 2-body Pauli bilinear,
    /// y_par(σ) = bit_a(σ) XOR bit_b(σ). The identity equivalently holds iff
    /// k_body is even. Verified by enumerating all 9 bilinears in
    /// <c>YParityIndependenceAtK3Tests.K2CollapseIdentity_HoldsForAllKBody2Bilinears</c>.</summary>
    public string K2CollapseIdentity => "y_par(σ) = bit_a(σ) XOR bit_b(σ)  for all k_body=2 bilinears (equivalently: for all even-k_body terms)";

    /// <summary>Canonical k_body=3 counterexample showing y_par is independent of Klein:
    /// XYZ has Klein (0, 0), y_par = 1; III has Klein (0, 0), y_par = 0.</summary>
    public string K3Counterexample => "XYZ vs III: same Klein (0, 0), different y_par (1 vs 0)";

    /// <summary>Returns true if for the given Pauli term σ (regardless of k-body),
    /// the k_body=2 collapse identity y_par = bit_a XOR bit_b holds. Always true at
    /// even k_body (k_body=0, 2, 4, ...); false at odd k_body (k_body=1, 3, ...).</summary>
    public static bool K2CollapseIdentityHolds(PauliTerm term)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        int bitA = term.TotalBitA & 1;
        int bitB = term.Pi2Parity;
        int yPar = term.YParity;
        return (bitA ^ bitB) == yPar;
    }

    /// <summary>Typed Cubic3 parent: <see cref="KleinEightCellClaim"/>. F102 is
    /// the first YParity-axis Claim and the structural anchor for the 8-cell
    /// Z₂³ decomposition where y_par becomes independent of the Klein
    /// (bit_a, bit_b) signature at k_body ≥ 3. Wired 2026-05-26.</summary>
    public KleinEightCellClaim KleinEightParent { get; }

    public YParityIndependenceAtK3(KleinEightCellClaim klein8)
        : base("F102 Y-parity term-level Z₂ independence at k_body≥3 (k_body=2 collapse to Klein, k_body≥3 independent); typed Cubic3 parent = KleinEightCellClaim",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F102 + " +
               "docs/proofs/PROOF_F102_YPARITY_INDEPENDENCE.md + " +
               "compute/RCPsiSquared.Core.Tests/Pauli/PauliHamiltonianKleinHelpersTests.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (KleinEightCellClaim, typed Cubic3 parent)")
    {
        KleinEightParent = klein8 ?? throw new ArgumentNullException(nameof(klein8));
    }

    public override string DisplayName =>
        "F102 Y-parity term-level Z₂ classifier (collapses at k_body=2, independent at k_body≥3)";

    public override string Summary =>
        $"y_par = (#Y mod 2) per Pauli term; at k_body=2 equals bit_a XOR bit_b (Klein-derived); " +
        $"at k_body≥3 independent (XYZ vs III: same Klein (0,0), different y_par) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Per-letter (bit_a, bit_b, y_par)",
                summary: "I=(0,0,0), X=(1,0,0), Y=(1,1,1), Z=(0,1,0); equivalently bit_a=n_X+n_Y, bit_b=n_Y+n_Z, y_par=n_Y all mod 2");
            yield return new InspectableNode("k_body=2 collapse identity",
                summary: K2CollapseIdentity);
            yield return new InspectableNode("k_body=3 counterexample",
                summary: K3Counterexample);
            yield return new InspectableNode("Why no operator-level third Π² axis",
                summary: "Per F34/QUBIT_NECESSITY only two independent Π² operators exist (Π²_Z = X⊗N, Π²_X = Z⊗N). Y-parity is term-level, surfacing as independent classifier only at k_body≥3.");
            yield return new InspectableNode("Scope (out of scope)",
                summary: "Stage 2a (bit_a twins), 2b (KleinEightCellClaim), 3b (Pi2KleinSpectralView 8-cell), 3c (F87 trichotomy in Z₂³) are separate specs.");
            yield return new InspectableNode("Cubic3 anchor parent",
                summary: $"KleinEightCellClaim ({KleinEightParent.Tier.Label()}): the Z₂³ 8-cell decomposition (bit_a, bit_b, y_par). F102 names the y_par axis where this Claim's third Z₂ classifier becomes independent at k_body ≥ 3.");
        }
    }
}
