using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F63 BitA twin reference (Tier 1 derived, lightweight bit_a-axis sibling
/// of <see cref="F63LCommutesPi2Pi2Inheritance"/>):
///
/// <code>
///   [L, Π²_X] = 0      exactly, for all N
/// </code>
///
/// <para>The BitA-axis statement <c>[L, Π²_X] = 0</c> IS the content of
/// <see cref="F61BitAParityPi2Inheritance"/> already typed in the registry: F61
/// proves <c>[L, Π²_X] = 0</c> with Π²_X = (−1)^{n_XY} = Z⊗N as the bit_a parity
/// operator. This Claim does NOT re-derive F61; it makes F61's role as F63's
/// BitA twin explicit at the typed-Claim level without introducing a constructor
/// cycle.</para>
///
/// <para><b>Why a separate Claim and not a direct edge</b>: F63 currently takes
/// F38 as a typed ctor parent. F61 takes F63 as its typed ctor parent. Wiring
/// F61 as F63's direct BitATwin via a ctor edge would close the cycle
/// F61 → F63 → F61. The cycle is unbreakable by a Lazy<> wrapper on F63's
/// optional twin parameter alone because F61 still needs F63 to be resolved.
/// The pragmatic resolution is: this lightweight reference Claim carries no
/// F61 ctor edge (intentionally; the F61 link lives in the docstring and
/// inspector summary). The typed edge from the BitB side
/// (<see cref="F63LCommutesPi2Pi2Inheritance"/>) is wired via that Claim's
/// optional ctor parameter pointing here.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitA"/>;
/// <see cref="BitATwin"/> is always null (this Claim IS the bit_a side);
/// BitATwinStatus is <see cref="BitATwinClassification.NotApplicableForThisAxis"/>.
/// The canonical bit_a Π² conservation derivation lives in
/// <see cref="F61BitAParityPi2Inheritance"/>; readers seeking the proof should
/// follow that Claim's anchors.</para>
///
/// <para>Tier1Derived: pure reference Claim. F61 is Tier 1 proven
/// (PROOF_PARITY_SELECTION_RULE; 64 configurations verified). F63 is Tier 1 proven
/// (PROOF_BIT_B_PARITY_SYMMETRY). No new derivation; the wiring is composition only.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F63 + F61 +
/// <c>compute/RCPsiSquared.Core/Symmetry/F61BitAParityPi2Inheritance.cs</c>
/// (the canonical bit_a Π² conservation Claim that contains the F63 bit_a sibling
/// content).</para></summary>
public sealed class F63BitAReference : Claim, IZ2AxisClaim
{
    /// <summary>BitA axis: Π²_X = Z⊗N, bit_a parity n_XY = #X + #Y mod 2.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitA;

    /// <summary>BitA-axis Claim: no BitATwin slot semantics.</summary>
    public Claim? BitATwin => null;

    /// <summary>BitA axis Claim: BitATwinStatus is NotApplicableForThisAxis.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.NotApplicableForThisAxis;

    /// <summary>The theorem statement in one line.</summary>
    public string Theorem =>
        "[L, Π²_X] = 0 exactly for all N; identical to F61's content (Π²_X = Z⊗N as bit_a parity operator). " +
        "The canonical derivation lives in F61BitAParityPi2Inheritance; this Claim is the typed bit_a-axis sibling slot for F63.";

    public F63BitAReference()
        : base("F63 BitA twin reference: [L, Π²_X] = 0 (typed bit_a sibling of F63 via F61's derivation; no F61 ctor edge to avoid cycle)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F63 + F61 + " +
               "compute/RCPsiSquared.Core/Symmetry/F61BitAParityPi2Inheritance.cs (canonical bit_a Π² conservation Claim)")
    {
    }

    public override string DisplayName =>
        "F63 BitA twin reference (F61 lives in the canonical bit_a Π² conservation Claim, cycle-free)";

    public override string Summary =>
        $"{Theorem} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("Canonical bit_a Π² conservation Claim",
                summary: "F61BitAParityPi2Inheritance: [L, Π²_X] = 0 exactly all N; 64 configurations verified (N=2..6, 4 topologies, 4 γ profiles); SE-accessibility ceiling on even-n_XY modes. F63BitAReference carries no F61 ctor edge to avoid the F61 → F63 → F61 cycle; the F61 link is documentation-only.");
            yield return new InspectableNode("BitB twin (F63)",
                summary: "F63LCommutesPi2Pi2Inheritance: [L, Π²_Z] = 0 exactly all N; six-line proof in PROOF_BIT_B_PARITY_SYMMETRY. F63BitAReference is the bit_a-axis sibling slot, content delegated to F61.");
            yield return new InspectableNode("Why a separate Claim",
                summary: "F61 → F63 → F38 is the existing ctor chain. Wiring F61 as F63's direct ctor BitATwin would close the cycle F61 → F63 → F61. F63BitAReference sits on the bit_a axis and exposes the bit_a sibling slot; F63 then takes F63BitAReference as an optional ctor parameter via the F1 ↔ F61 pattern. F63BitAReference carries no F61 ctor edge to keep the inheritance graph acyclic.");
        }
    }
}
