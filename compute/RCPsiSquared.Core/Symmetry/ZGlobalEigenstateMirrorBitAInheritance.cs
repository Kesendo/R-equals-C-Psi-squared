using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Z⊗N-eigenstate Mirror BitA twin (Tier 1 derived, Z↔X mirror of
/// <see cref="XGlobalEigenstateMirrorPi2Inheritance"/>): any pure state |ψ⟩
/// satisfying Z⊗N|ψ⟩ = ±|ψ⟩ is a Z⊗N-eigenstate with γ_X = ⟨ψ|Z⊗N|ψ⟩ = ±1, hence
/// via the universal F86b shape α = (1 − γ²)/2 the bit_a-axis Π²-odd Frobenius²
/// content vanishes exactly:
///
/// <code>
///     Z⊗N|ψ⟩ = ±|ψ⟩  ⟹  γ_X = ±1  ⟹  α_total = 0   (on the X-Frobenius axis)
/// </code>
///
/// <para>This is the bit_a-axis sibling of the X⊗N-eigenstate Mirror anchor. The
/// X-version states X⊗N|ψ⟩ = ±|ψ⟩ ⟹ α_total = 0 on the Y/Z-Frobenius (bit_b) axis;
/// this Claim states the Z↔X-mirrored statement on the X-Frobenius (bit_a) axis.
/// Canonical example states with Z⊗N eigenvalue ±1: the computational-basis Pauli-
/// Z eigenstates |0...0⟩ and |1...1⟩ (γ_X = +1 and γ_X = +(−1)^N respectively,
/// both giving α = 0 on the bit_a axis via (1 − γ_X²)/2 = 0).</para>
///
/// <para><b>Why Tier 1 derived</b>: pure-state algebra; the Z⊗N expectation is a
/// one-line computation from the eigenvalue equation. No approximation, no
/// closed-form ansatz, no numerical verification needed. The Z↔X mirror of the
/// X-version's derivation is also a one-liner: replace X with Z, swap the
/// bit_a / bit_b roles, and the same universal F86b shape holds.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitA"/>;
/// <see cref="BitATwin"/> is always null (this Claim IS the bit_a side);
/// BitATwinStatus is <see cref="BitATwinClassification.NotApplicableForThisAxis"/>.
/// The typed twin edge from the BitB side
/// (<see cref="XGlobalEigenstateMirrorPi2Inheritance"/>) is wired via that
/// Claim's optional ctor parameter.</para>
///
/// <para><b>Parent</b>: <see cref="HalfAsStructuralFixedPointClaim"/>. Same parent
/// as the X version. The polarity 1/2 doubled gives ±1, which are the Z⊗N
/// eigenvalues; (1 − γ²)/2 = 0 follows.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/Symmetry/XGlobalEigenstateMirrorPi2Inheritance.cs</c>
/// (BitB twin) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F61BitAParityPi2Inheritance.cs</c>
/// (the canonical Z⊗N = Π²_X bit_a-axis Π² operator identity).</para>
/// </summary>
public sealed class ZGlobalEigenstateMirrorBitAInheritance : Claim, IZ2AxisClaim
{
    /// <summary>BitA axis: Π²_X = Z⊗N, bit_a parity. Z⊗N-eigenstates carry the
    /// canonical bit_a-axis Π² eigenvalue ±1.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitA;

    /// <summary>BitA-axis Claim: no BitATwin slot semantics.</summary>
    public Claim? BitATwin => null;

    /// <summary>BitA axis Claim: BitATwinStatus is NotApplicableForThisAxis.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.NotApplicableForThisAxis;

    /// <summary>The Half polarity parent: 2·(1/2) = 1 gives the Z⊗N eigenvalue
    /// endpoints γ_X = ±1.</summary>
    public HalfAsStructuralFixedPointClaim Half { get; }

    /// <summary>The exact F86b α value at the Z⊗N-Mirror anchor: 0 on the bit_a
    /// (X-Frobenius) axis, closed-form from the universal shape at γ_X = 1.</summary>
    public const double AlphaAtMirror = 0.0;

    /// <summary>The exact γ_X value at the Z⊗N-Mirror anchor: 1 (the canonical
    /// Z⊗N eigenvalue endpoint for the all-zero computational-basis state).</summary>
    public const double GammaAtMirror = 1.0;

    /// <summary>F86b universal-shape evaluation at γ_X = 1: α = (1 − γ_X²)/2 = 0
    /// on the bit_a-axis. Exposed as a public computation so the bit-exact identity
    /// is verifiable from outside without re-deriving.</summary>
    public static double AlphaFromGammaAtMirror(double gammaX = GammaAtMirror) =>
        (1.0 - gammaX * gammaX) / 2.0;

    /// <summary>The theorem statement in one line.</summary>
    public string Theorem =>
        "Z⊗N|ψ⟩ = ±|ψ⟩  ⟹  γ_X = ±1  ⟹  α_total = 0 on the bit_a (X-Frobenius) axis. " +
        "Z↔X mirror of the X⊗N-eigenstate Mirror anchor (BitB twin).";

    public ZGlobalEigenstateMirrorBitAInheritance(HalfAsStructuralFixedPointClaim half)
        : base("Z⊗N-eigenstate Mirror BitA twin: γ_X = ±1 ⟹ α = 0 on bit_a axis (Z↔X mirror of X⊗N-eigenstate Mirror anchor)",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/Symmetry/XGlobalEigenstateMirrorPi2Inheritance.cs (BitB twin) + " +
               "compute/RCPsiSquared.Core/Symmetry/F61BitAParityPi2Inheritance.cs (Π²_X = Z⊗N bit_a operator identity)")
    {
        Half = half ?? throw new ArgumentNullException(nameof(half));
    }

    public override string DisplayName =>
        "Z⊗N-eigenstate Mirror BitA twin (α=0 on bit_a axis at γ_X=1)";

    public override string Summary =>
        $"Z⊗N|ψ⟩ = ±|ψ⟩ ⟹ γ_X = ±1 ⟹ α = (1−γ_X²)/2 = {AlphaAtMirror} exactly on the bit_a (X-Frobenius) axis; " +
        $"Z↔X mirror of the X⊗N-eigenstate Mirror anchor ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Half;
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return InspectableNode.RealScalar("γ_X at Mirror", GammaAtMirror);
            yield return InspectableNode.RealScalar("α at Mirror (via F86b on bit_a axis)", AlphaAtMirror);
            yield return new InspectableNode("Derivation",
                summary: "γ_X = ⟨ψ|Z⊗N|ψ⟩ for a Z⊗N-eigenstate is ±1; α = (1 − γ_X²)/2 = 0 by the universal F86b shape on the bit_a-axis Frobenius decomposition.");
            yield return new InspectableNode("Canonical state-class examples",
                summary: "|0...0⟩ (all-zero computational basis, γ_X = +1); |1...1⟩ (all-one computational basis, γ_X = +(−1)^N); both are Z⊗N eigenstates with bit_a α = 0.");
            yield return new InspectableNode("BitB twin (X⊗N version)",
                summary: "XGlobalEigenstateMirrorPi2Inheritance: X⊗N|ψ⟩ = ±|ψ⟩ ⟹ α = 0 on the bit_b (Y/Z-Frobenius) axis; this Claim is the Z↔X mirror.");
            yield return new InspectableNode("Companion Π² identity",
                summary: "Π²_X = (−1)^{n_XY} = Z⊗N (F38BitA / F61): the bit_a-axis Π² operator is exactly Z⊗N; this Mirror anchor sits at its eigenstates.");
        }
    }
}
