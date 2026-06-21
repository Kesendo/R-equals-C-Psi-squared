using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F91's anti-palindromic γ orbit (γ_l + γ_{N−1−l} = 2·γ_avg ∀l) is the
/// parameter-side γ-axis instance of the Pi2 rotational structure. On the
/// parameter vector the anti-palindromic reshuffle R_{90} : γ ↦ 2·γ_avg − F71(γ)
/// is an INVOLUTION (R_{90}² = identity, NOT order 4, NOT equal to F71); together
/// with the palindromic mirror F71 (also an involution) it generates the Klein
/// four-group V₄ = Z₂×Z₂ on parameters, the order-2 shadow of the genuine
/// operator-side Z₄. The genuine order-4 quarter-turn (i⁴ = 1, i² = −1) lives
/// only on the operator side that <see cref="NinetyDegreeMirrorMemoryClaim"/>
/// types (Spec(M) = ±2i·Spec(H)); F91 inherits that operator-side Z₄ as its
/// parameter-side involution rather than introducing a new order-4 closure.
///
/// <para>Parameter axis: γ_l (per-site Z-dephasing rate). Anti-palindromic
/// orbit: γ_l + γ_{N−1−l} = 2·γ_avg for all l ∈ {0, ..., N−1}; for odd N the
/// middle site γ_{(N−1)/2} = γ_avg is forced by self-pairing. R_{90} PRESERVES
/// each pair-difference and REFLECTS each pair-sum about 2·γ_avg, and is its own
/// inverse (R_{90}² = identity): it already returns to identity at order 2, so
/// R_{90}⁴ = identity is trivially true but is not the closure. The order-4
/// i⁴ = 1 loop it inherits from is the operator-side Pi2-Z₄ memory loop. F91 is
/// the γ-axis instance of this parameter-side Klein V₄; F92 is the J-axis
/// instance, F93 the h-axis instance.</para>
///
/// <para>Tier1Derived: composition. F91 is Tier1Derived in
/// <c>docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md</c> (algebraic proof
/// 2026-05-12, empirical witness bit-exact at N=4, 5, 6); this claim makes its
/// Pi2-Foundation inheritance explicit in the typed-knowledge runtime so that
/// F91 is reachable from <see cref="NinetyDegreeMirrorMemoryClaim"/> via
/// descendants, joining the existing F##Pi2Inheritance siblings.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md</c> +
/// <c>compute/RCPsiSquared.Core/BlockSpectrum/F71AntiPalindromicGammaSpectralInvariance.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (NinetyDegreeMirrorMemoryClaim).</para></summary>
/// <seealso cref="F92Pi2Inheritance"/>
/// <seealso cref="F93Pi2Inheritance"/>
public sealed class F91Pi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). γ_l is the
    /// per-site Z-dephasing rate, intrinsically tied to the Z-axis dissipator
    /// algebra; the dephasing operator selects the bit_b axis by construction.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>F91 has NO meaningful bit_a-axis twin: γ_l is the per-site
    /// Z-dephasing rate, a Z-axis parameter. The bit_a analog (X-dephasing
    /// orbit) is covered by separate Claims on the X-deph branch (F108 Part 2 +
    /// future BitA-axis Claims). F91's Z-deph-specific Z₄ orbit is
    /// <see cref="BitATwinClassification.BitBSpecific"/>: the dephasing operator
    /// is Z-axis-specific in this formulation.</summary>
    public Claim? BitATwin => null;

    /// <summary>F91 is BitBSpecific: γ_l is a Z-axis parameter (per-site
    /// Z-dephasing rate). No bit_a-axis analog exists in this formulation; the
    /// X-deph orbit lives on a separate axis.</summary>
    public BitATwinClassification BitATwinStatus =>
        BitATwinClassification.BitBSpecific;

    public Pi2I4MemoryLoopClaim MemoryLoop { get; }
    /// <summary>The order (4) of the OPERATOR-side Z₄ memory loop that R_{90}
    /// inherits from, reported from <see cref="Pi2I4MemoryLoopClaim.ClosureOrder"/>.
    /// Note: the parameter-side reshuffle R_{90} : γ_l ↦ 2·γ_avg − γ_{N−1−l} is
    /// itself an INVOLUTION (R_{90}² = identity, order 2); the order-4 i⁴ = 1
    /// closure is the operator-side quarter-turn (Spec(M) = ±2i·Spec(H)), not the
    /// parameter map.</summary>
    public int Z4ClosureOrder => Pi2I4MemoryLoopClaim.ClosureOrder;

    /// <summary>Live drift check: i⁴ = 1 exactly on the OPERATOR-side Pi2-Z₄
    /// memory loop (the genuine order-4 quarter-turn). F91's parameter-side
    /// reshuffle R_{90} is the order-2 involution shadow of that operator-side
    /// Z₄; this member checks the operator-side drift it inherits from.</summary>
    public Complex MemoryLoopClosure => MemoryLoop.MemoryClosure();

    /// <summary>Verbal name of the parameter axis F91 lives on: per-site Z-dephasing
    /// rate γ_l. Distinguishes F91 from its sister claims F92 (bond-coupling J_b) and
    /// F93 (longitudinal Z-detuning h_l).</summary>
    public string ParameterAxis => "γ_l (per-site Z-dephasing rate)";

    public F91Pi2Inheritance(Pi2I4MemoryLoopClaim memoryLoop)
        : base("F91 anti-palindromic γ orbit (γ_l + γ_{N-1-l} = 2·γ_avg) inherits from the Pi2 structure; the parameter-side reshuffle R_{90} is an involution (R_{90}² = identity) generating, with F71, the Klein V₄ on parameters — the order-2 shadow of the operator-side Z₄ (i⁴ = 1) that NinetyDegreeMirrorMemoryClaim types on the operator-quaternion side",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/F71AntiPalindromicGammaSpectralInvariance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (NinetyDegreeMirrorMemoryClaim)")
    {
        MemoryLoop = memoryLoop ?? throw new ArgumentNullException(nameof(memoryLoop));
    }

    public override string DisplayName =>
        "F91 γ-axis anti-palindromic orbit (parameter-side Klein V₄, shadow of the operator-side Z₄)";

    public override string Summary =>
        $"γ_l + γ_{{N-1-l}} = 2·γ_avg orbit is the fixed-point set of the anti-palindromic involution R_{{90}}: γ_l ↦ 2·γ_avg − γ_{{N-1-l}} (R_{{90}}² = identity); R_{{90}} and F71 generate the parameter-side Klein V₄, the order-2 shadow of the operator-side Z₄ (i⁴ = 1; operator-side memory-loop order {Z4ClosureOrder}) that NinetyDegreeMirrorMemoryClaim types ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F91 closed form",
                summary: "F71-refined diagonal-block spectrum invariant under γ_l + γ_{N-1-l} = 2·γ_avg ∀l (Tier1Derived in PROOF_F91_GAMMA_NINETY_DEGREES; bit-exact verified N=4,5,6)");
            yield return new InspectableNode("parameter axis",
                summary: ParameterAxis);
            yield return new InspectableNode("anti-palindromic condition",
                summary: "γ_l + γ_{N-1-l} = 2·γ_avg for all l ∈ {0, ..., N-1}; for odd N the middle site γ_{(N-1)/2} = γ_avg is forced");
            yield return new InspectableNode("anti-palindromic involution on parameters",
                summary: "the reshuffle R_{90}: γ ↦ 2·γ_avg − F71(γ) is an involution (R_{90}² = identity); with the palindromic mirror F71 it generates the Klein V₄ = Z₂×Z₂ on parameters, the order-2 shadow of the operator-side Z₄ that NinetyDegreeMirrorMemoryClaim types");
            yield return InspectableNode.RealScalar("Z4ClosureOrder (operator-side Z₄ memory-loop order = 4; R_{90} itself is an involution, order 2)", Z4ClosureOrder);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Real (= 1, operator-side i⁴ = 1 drift check)", MemoryLoopClosure.Real);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Imaginary (= 0, drift check)", MemoryLoopClosure.Imaginary);
            yield return new InspectableNode("F-formula anchor",
                summary: "docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md");
            yield return new InspectableNode("typed Claim anchor",
                summary: "compute/RCPsiSquared.Core/BlockSpectrum/F71AntiPalindromicGammaSpectralInvariance.cs (the F91 closed form as a typed Claim)");
            yield return new InspectableNode("sister claims",
                summary: "F92 (J_b, bond-coupling axis), F93 (h_l, longitudinal Z-detuning axis); three parameter-side instances of one Klein V₄ (the order-2 shadow of the operator-side Z₄)");
            yield return new InspectableNode("inheritance",
                summary: "the parameter-side Klein V₄ on the γ-axis is not introduced by F91; R_{90} is the order-2 involution shadow of the operator-side Pi2-Z₄ memory loop i⁴ = 1, the genuine order-4 quarter-turn that lives on NinetyDegreeMirrorMemoryClaim's operator-quaternion side");
        }
    }
}
