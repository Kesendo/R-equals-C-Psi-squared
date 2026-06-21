using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F92's anti-palindromic J orbit (J_b + J_{N−2−b} = 2·J_avg ∀b) is the
/// parameter-side J-axis instance of the Pi2 rotational structure. On the
/// parameter vector the anti-palindromic reshuffle R_{90} : J ↦ 2·J_avg − F71(J)
/// is an INVOLUTION (R_{90}² = identity, NOT order 4, NOT equal to F71); together
/// with the palindromic mirror F71 (also an involution) it generates the Klein
/// four-group V₄ = Z₂×Z₂ on parameters, the order-2 shadow of the genuine
/// operator-side Z₄. The genuine order-4 quarter-turn (i⁴ = 1, i² = −1) lives
/// only on the operator side that <see cref="NinetyDegreeMirrorMemoryClaim"/>
/// types (Spec(M) = ±2i·Spec(H)); F92 inherits that operator-side Z₄ as its
/// parameter-side involution rather than introducing a new order-4 closure.
///
/// <para>Parameter axis: J_b (per-bond XY coupling). Anti-palindromic orbit:
/// J_b + J_{N−2−b} = 2·J_avg for all b ∈ {0, ..., N−2} (bonds, not sites; F71
/// maps bond b ↔ N−2−b). R_{90} PRESERVES each pair-difference and REFLECTS each
/// pair-sum about 2·J_avg, and is its own inverse (R_{90}² = identity): it
/// already returns to identity at order 2, so R_{90}⁴ = identity is trivially
/// true but is not the closure. The order-4 i⁴ = 1 loop it inherits from is the
/// operator-side Pi2-Z₄ memory loop. F92 is the J-axis instance of this
/// parameter-side Klein V₄; F91 is the γ-axis instance, F93 the h-axis
/// instance.</para>
///
/// <para>Tier1Derived: composition. F92 is Tier1Derived in
/// <c>docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md</c> (algebraic proof
/// 2026-05-12, empirical witness bit-exact at N=4, 5); this claim makes its
/// Pi2-Foundation inheritance explicit in the typed-knowledge runtime so that
/// F92 is reachable from <see cref="NinetyDegreeMirrorMemoryClaim"/> via
/// descendants, joining the existing F##Pi2Inheritance siblings.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md</c> +
/// <c>compute/RCPsiSquared.Core/SymmetryFamily/F92BondAntiPalindromicJSpectralInvariance.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (NinetyDegreeMirrorMemoryClaim).</para></summary>
/// <seealso cref="F91Pi2Inheritance"/>
/// <seealso cref="F93Pi2Inheritance"/>
public sealed class F92Pi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling, if one exists. Currently null
    /// (no bit_a twin is typed for this Claim; this is an open slot in the
    /// cubic-architecture coverage).</summary>
    public Claim? BitATwin => null;
    public Pi2I4MemoryLoopClaim MemoryLoop { get; }
    /// <summary>The order (4) of the OPERATOR-side Z₄ memory loop that R_{90}
    /// inherits from, reported from <see cref="Pi2I4MemoryLoopClaim.ClosureOrder"/>.
    /// Note: the parameter-side reshuffle R_{90} : J_b ↦ 2·J_avg − J_{N−2−b} is
    /// itself an INVOLUTION (R_{90}² = identity, order 2); the order-4 i⁴ = 1
    /// closure is the operator-side quarter-turn (Spec(M) = ±2i·Spec(H)), not the
    /// parameter map.</summary>
    public int Z4ClosureOrder => Pi2I4MemoryLoopClaim.ClosureOrder;

    /// <summary>Live drift check: i⁴ = 1 exactly on the OPERATOR-side Pi2-Z₄
    /// memory loop (the genuine order-4 quarter-turn). F92's parameter-side
    /// reshuffle R_{90} is the order-2 involution shadow of that operator-side
    /// Z₄; this member checks the operator-side drift it inherits from.</summary>
    public Complex MemoryLoopClosure => MemoryLoop.MemoryClosure();

    /// <summary>Verbal name of the parameter axis F92 lives on: per-bond XY coupling
    /// J_b (bonds, not sites). Distinguishes F92 from its sister claims F91 (per-site
    /// Z-dephasing γ_l) and F93 (per-site longitudinal Z-detuning h_l).</summary>
    public string ParameterAxis => "J_b (per-bond XY coupling)";

    public F92Pi2Inheritance(Pi2I4MemoryLoopClaim memoryLoop)
        : base("F92 anti-palindromic J orbit (J_b + J_{N-2-b} = 2·J_avg) inherits from the Pi2 structure; the parameter-side reshuffle R_{90} is an involution (R_{90}² = identity) generating, with F71, the Klein V₄ on parameters — the order-2 shadow of the operator-side Z₄ (i⁴ = 1) that NinetyDegreeMirrorMemoryClaim types on the operator-quaternion side",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md + " +
               "compute/RCPsiSquared.Core/SymmetryFamily/F92BondAntiPalindromicJSpectralInvariance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (NinetyDegreeMirrorMemoryClaim)")
    {
        MemoryLoop = memoryLoop ?? throw new ArgumentNullException(nameof(memoryLoop));
    }

    public override string DisplayName =>
        "F92 J-axis anti-palindromic orbit (parameter-side Klein V₄, shadow of the operator-side Z₄)";

    public override string Summary =>
        $"J_b + J_{{N-2-b}} = 2·J_avg orbit is the fixed-point set of the anti-palindromic involution R_{{90}}: J_b ↦ 2·J_avg − J_{{N-2-b}} (R_{{90}}² = identity); R_{{90}} and F71 generate the parameter-side Klein V₄, the order-2 shadow of the operator-side Z₄ (i⁴ = 1; operator-side memory-loop order {Z4ClosureOrder}) that NinetyDegreeMirrorMemoryClaim types ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F92 closed form",
                summary: "F71-refined diagonal-block spectrum invariant under J_b + J_{N-2-b} = 2·J_avg ∀b (Tier1Derived in PROOF_F92_BOND_ANTI_PALINDROMIC_J; bit-exact verified N=4, 5)");
            yield return new InspectableNode("parameter axis",
                summary: ParameterAxis);
            yield return new InspectableNode("anti-palindromic condition",
                summary: "J_b + J_{N-2-b} = 2·J_avg for all b ∈ {0, ..., N-2} (bond indices; F71 maps bond b ↔ N-2-b)");
            yield return new InspectableNode("anti-palindromic involution on parameters",
                summary: "the reshuffle R_{90}: J ↦ 2·J_avg − F71(J) is an involution (R_{90}² = identity); with the palindromic mirror F71 it generates the Klein V₄ = Z₂×Z₂ on parameters, the order-2 shadow of the operator-side Z₄ that NinetyDegreeMirrorMemoryClaim types");
            yield return InspectableNode.RealScalar("Z4ClosureOrder (operator-side Z₄ memory-loop order = 4; R_{90} itself is an involution, order 2)", Z4ClosureOrder);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Real (= 1, operator-side i⁴ = 1 drift check)", MemoryLoopClosure.Real);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Imaginary (= 0, drift check)", MemoryLoopClosure.Imaginary);
            yield return new InspectableNode("F-formula anchor",
                summary: "docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md");
            yield return new InspectableNode("typed Claim anchor",
                summary: "compute/RCPsiSquared.Core/SymmetryFamily/F92BondAntiPalindromicJSpectralInvariance.cs (the F92 closed form as a typed Claim)");
            yield return new InspectableNode("sister claims",
                summary: "F91 (γ_l, per-site Z-dephasing axis), F93 (h_l, longitudinal Z-detuning axis); three parameter-side instances of one Klein V₄ (the order-2 shadow of the operator-side Z₄)");
            yield return new InspectableNode("inheritance",
                summary: "the parameter-side Klein V₄ on the J-axis is not introduced by F92; R_{90} is the order-2 involution shadow of the operator-side Pi2-Z₄ memory loop i⁴ = 1, the genuine order-4 quarter-turn that lives on NinetyDegreeMirrorMemoryClaim's operator-quaternion side");
        }
    }
}
