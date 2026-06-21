using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F93's anti-palindromic h orbit (h_l + h_{N−1−l} = 2·h_avg ∀l) is the
/// parameter-side h-axis instance of the Pi2 rotational structure. On the
/// parameter vector the anti-palindromic reshuffle R_{90} : h ↦ 2·h_avg − F71(h)
/// is an INVOLUTION (R_{90}² = identity, NOT order 4, NOT equal to F71); together
/// with the palindromic mirror F71 (also an involution) it generates the Klein
/// four-group V₄ = Z₂×Z₂ on parameters, the order-2 shadow of the genuine
/// operator-side Z₄. The genuine order-4 quarter-turn (i⁴ = 1, i² = −1) lives
/// only on the operator side that <see cref="NinetyDegreeMirrorMemoryClaim"/>
/// types (Spec(M) = ±2i·Spec(H)); F93 inherits that operator-side Z₄ as its
/// parameter-side involution rather than introducing a new order-4 closure.
///
/// <para>Parameter axis: h_l (per-site longitudinal Z-detuning). Anti-palindromic
/// orbit: h_l + h_{N−1−l} = 2·h_avg for all l ∈ {0, ..., N−1}; for odd N the
/// middle site h_{(N−1)/2} = h_avg is forced by self-pairing. R_{90} PRESERVES
/// each pair-difference and REFLECTS each pair-sum about 2·h_avg, and is its own
/// inverse (R_{90}² = identity): it already returns to identity at order 2, so
/// R_{90}⁴ = identity is trivially true but is not the closure. The order-4
/// i⁴ = 1 loop it inherits from is the operator-side Pi2-Z₄ memory loop. F93 is
/// the h-axis instance of this parameter-side Klein V₄; F91 is the γ-axis
/// instance, F92 the J-axis instance.</para>
///
/// <para>Scope note: only longitudinal h_l Z_l is in scope. Transverse h_l X_l or
/// h_l Y_l would flip popcount and break joint-popcount conservation, taking F93
/// out of the BlockSpectrum framework entirely.</para>
///
/// <para>Tier1Derived: composition. F93 is Tier1Derived in
/// <c>docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md</c> (algebraic proof
/// 2026-05-12, empirical witness bit-exact at N=4, 5); this claim makes its
/// Pi2-Foundation inheritance explicit in the typed-knowledge runtime so that
/// F93 is reachable from <see cref="NinetyDegreeMirrorMemoryClaim"/> via
/// descendants, joining the existing F##Pi2Inheritance siblings.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md</c> +
/// <c>compute/RCPsiSquared.Core/SymmetryFamily/F93DetuningAntiPalindromicSpectralInvariance.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (NinetyDegreeMirrorMemoryClaim).</para></summary>
/// <seealso cref="F91Pi2Inheritance"/>
/// <seealso cref="F92Pi2Inheritance"/>
public sealed class F93Pi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). h_l is the
    /// per-site longitudinal Z-detuning. The Claim docstring explicitly excludes
    /// transverse fields h_l X_l or h_l Y_l because they break joint-popcount
    /// conservation; the F93 formulation is Z-axis only.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>F93 has NO meaningful bit_a-axis twin: h_l Z_l is the
    /// longitudinal Z-detuning, intrinsically a Z-axis parameter. The Claim
    /// already documents that transverse h_l X_l is out of scope (breaks
    /// joint-popcount per F61's BreakConditions). Classified as
    /// <see cref="BitATwinClassification.BitBSpecific"/>.</summary>
    public Claim? BitATwin => null;

    /// <summary>F93 is BitBSpecific: h_l Z_l is a Z-axis parameter
    /// (longitudinal Z-detuning). Transverse h_l X_l is explicitly out of scope
    /// per the Claim docstring (breaks joint-popcount conservation, takes F93
    /// out of the BlockSpectrum framework entirely).</summary>
    public BitATwinClassification BitATwinStatus =>
        BitATwinClassification.BitBSpecific;

    public Pi2I4MemoryLoopClaim MemoryLoop { get; }
    /// <summary>The order (4) of the OPERATOR-side Z₄ memory loop that R_{90}
    /// inherits from, reported from <see cref="Pi2I4MemoryLoopClaim.ClosureOrder"/>.
    /// Note: the parameter-side reshuffle R_{90} : h_l ↦ 2·h_avg − h_{N−1−l} is
    /// itself an INVOLUTION (R_{90}² = identity, order 2); the order-4 i⁴ = 1
    /// closure is the operator-side quarter-turn (Spec(M) = ±2i·Spec(H)), not the
    /// parameter map.</summary>
    public int Z4ClosureOrder => Pi2I4MemoryLoopClaim.ClosureOrder;

    /// <summary>Live drift check: i⁴ = 1 exactly on the OPERATOR-side Pi2-Z₄
    /// memory loop (the genuine order-4 quarter-turn). F93's parameter-side
    /// reshuffle R_{90} is the order-2 involution shadow of that operator-side
    /// Z₄; this member checks the operator-side drift it inherits from.</summary>
    public Complex MemoryLoopClosure => MemoryLoop.MemoryClosure();

    /// <summary>Verbal name of the parameter axis F93 lives on: per-site longitudinal
    /// Z-detuning h_l. Distinguishes F93 from its sister claims F91 (per-site
    /// Z-dephasing γ_l) and F92 (per-bond XY coupling J_b).</summary>
    public string ParameterAxis => "h_l (per-site longitudinal Z-detuning)";

    public F93Pi2Inheritance(Pi2I4MemoryLoopClaim memoryLoop)
        : base("F93 anti-palindromic h orbit (h_l + h_{N-1-l} = 2·h_avg) inherits from the Pi2 structure; the parameter-side reshuffle R_{90} is an involution (R_{90}² = identity) generating, with F71, the Klein V₄ on parameters — the order-2 shadow of the operator-side Z₄ (i⁴ = 1) that NinetyDegreeMirrorMemoryClaim types on the operator-quaternion side; longitudinal h_l Z_l only (transverse breaks joint-popcount, out of scope)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md + " +
               "compute/RCPsiSquared.Core/SymmetryFamily/F93DetuningAntiPalindromicSpectralInvariance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (NinetyDegreeMirrorMemoryClaim)")
    {
        MemoryLoop = memoryLoop ?? throw new ArgumentNullException(nameof(memoryLoop));
    }

    public override string DisplayName =>
        "F93 h-axis anti-palindromic orbit (parameter-side Klein V₄, shadow of the operator-side Z₄)";

    public override string Summary =>
        $"h_l + h_{{N-1-l}} = 2·h_avg orbit (longitudinal Z-detuning only) is the fixed-point set of the anti-palindromic involution R_{{90}}: h_l ↦ 2·h_avg − h_{{N-1-l}} (R_{{90}}² = identity); R_{{90}} and F71 generate the parameter-side Klein V₄, the order-2 shadow of the operator-side Z₄ (i⁴ = 1; operator-side memory-loop order {Z4ClosureOrder}) that NinetyDegreeMirrorMemoryClaim types ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F93 closed form",
                summary: "F71-refined diagonal-block spectrum invariant under h_l + h_{N-1-l} = 2·h_avg ∀l (Tier1Derived in PROOF_F93_DETUNING_ANTI_PALINDROMIC; bit-exact verified N=4, 5)");
            yield return new InspectableNode("parameter axis",
                summary: ParameterAxis);
            yield return new InspectableNode("scope",
                summary: "longitudinal h_l Z_l only; transverse h_l X_l or h_l Y_l flips popcount and breaks joint-popcount conservation, out of BlockSpectrum scope");
            yield return new InspectableNode("anti-palindromic condition",
                summary: "h_l + h_{N-1-l} = 2·h_avg for all l ∈ {0, ..., N-1}; for odd N the middle site h_{(N-1)/2} = h_avg is forced");
            yield return new InspectableNode("anti-palindromic involution on parameters",
                summary: "the reshuffle R_{90}: h ↦ 2·h_avg − F71(h) is an involution (R_{90}² = identity); with the palindromic mirror F71 it generates the Klein V₄ = Z₂×Z₂ on parameters, the order-2 shadow of the operator-side Z₄ that NinetyDegreeMirrorMemoryClaim types");
            yield return InspectableNode.RealScalar("Z4ClosureOrder (operator-side Z₄ memory-loop order = 4; R_{90} itself is an involution, order 2)", Z4ClosureOrder);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Real (= 1, operator-side i⁴ = 1 drift check)", MemoryLoopClosure.Real);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Imaginary (= 0, drift check)", MemoryLoopClosure.Imaginary);
            yield return new InspectableNode("F-formula anchor",
                summary: "docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md");
            yield return new InspectableNode("typed Claim anchor",
                summary: "compute/RCPsiSquared.Core/SymmetryFamily/F93DetuningAntiPalindromicSpectralInvariance.cs (the F93 closed form as a typed Claim)");
            yield return new InspectableNode("sister claims",
                summary: "F91 (γ_l, per-site Z-dephasing axis), F92 (J_b, per-bond XY coupling axis); three parameter-side instances of one Klein V₄ (the order-2 shadow of the operator-side Z₄)");
            yield return new InspectableNode("inheritance",
                summary: "the parameter-side Klein V₄ on the h-axis is not introduced by F93; R_{90} is the order-2 involution shadow of the operator-side Pi2-Z₄ memory loop i⁴ = 1, the genuine order-4 quarter-turn that lives on NinetyDegreeMirrorMemoryClaim's operator-quaternion side");
        }
    }
}
