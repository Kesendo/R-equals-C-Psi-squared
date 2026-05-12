using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F93's anti-palindromic h orbit (h_l + h_{N−1−l} = 2·h_avg ∀l) is the
/// parameter-side h-axis instance of the Pi2-Z₄ rotational structure. The
/// 90°-rotation h ↦ 2·h_avg − F71(h) is one of the four Z₄ elements; closing
/// under i⁴ = 1 returns to identity. The Z₄ generator on parameters is
/// structurally the same Z₄ that
/// <see cref="NinetyDegreeMirrorMemoryClaim"/> types on the operator-quaternion
/// side; F93 inherits that algebra rather than introducing a new orbit closure.
///
/// <para>Parameter axis: h_l (per-site longitudinal Z-detuning). Anti-palindromic
/// orbit: h_l + h_{N−1−l} = 2·h_avg for all l ∈ {0, ..., N−1}; for odd N the
/// middle site h_{(N−1)/2} = h_avg is forced by self-pairing. Closure under the
/// 90°-rotation: applied four times R_{90}^4 = identity, which is exactly the
/// i⁴ = 1 closure on the Pi2-Z₄ memory loop. F93 is the h-axis instance of one
/// rotational structure; F91 is the γ-axis instance, F92 the J-axis instance.</para>
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
public sealed class F93Pi2Inheritance : Claim
{
    private readonly Pi2I4MemoryLoopClaim _memoryLoop;

    /// <summary>The Z₄ closure order of the 90°-rotation R_{90}: h_l ↦ 2·h_avg − h_{N−1−l}.
    /// Four applications return to identity, matching <see cref="Pi2I4MemoryLoopClaim.ClosureOrder"/>.</summary>
    public int Z4ClosureOrder => Pi2I4MemoryLoopClaim.ClosureOrder;

    /// <summary>Live drift check: i^4 = 1 exactly on the parent Pi2-Z₄ memory loop.
    /// F93's parameter-side 90°-rotation closes at the same order as the operator-side
    /// quaternion algebra.</summary>
    public Complex MemoryLoopClosure => _memoryLoop.MemoryClosure();

    /// <summary>Verbal name of the parameter axis F93 lives on: per-site longitudinal
    /// Z-detuning h_l. Distinguishes F93 from its sister claims F91 (per-site
    /// Z-dephasing γ_l) and F92 (per-bond XY coupling J_b).</summary>
    public string ParameterAxis => "h_l (per-site longitudinal Z-detuning)";

    public F93Pi2Inheritance(Pi2I4MemoryLoopClaim memoryLoop)
        : base("F93 anti-palindromic h orbit (h_l + h_{N-1-l} = 2·h_avg) inherits from Pi2-Z₄ structure (i⁴ = 1 closure); parameter-side h-axis instance of the same Z₄ that NinetyDegreeMirrorMemoryClaim types on the operator-quaternion side; longitudinal h_l Z_l only (transverse breaks joint-popcount, out of scope)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md + " +
               "compute/RCPsiSquared.Core/SymmetryFamily/F93DetuningAntiPalindromicSpectralInvariance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (NinetyDegreeMirrorMemoryClaim)")
    {
        _memoryLoop = memoryLoop ?? throw new ArgumentNullException(nameof(memoryLoop));
    }

    public override string DisplayName =>
        "F93 h-Z₄ anti-palindromic orbit as Pi2-Foundation parameter-side instance";

    public override string Summary =>
        $"h_l + h_{{N-1-l}} = 2·h_avg orbit (longitudinal Z-detuning only) closes under 90°-rotation R_{{90}}: h_l ↦ 2·h_avg − h_{{N-1-l}} at order {Z4ClosureOrder} (i⁴ = 1); parameter-side h-axis instance of the Pi2-Z₄ operator-quaternion structure ({Tier.Label()})";

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
            yield return new InspectableNode("Z₄ generator on parameters",
                summary: "the 90°-rotation h ↦ 2·h_avg − F71(h) is one of four Z₄ elements; structurally the same Z₄ that NinetyDegreeMirrorMemoryClaim types on the operator-quaternion side");
            yield return InspectableNode.RealScalar("Z4ClosureOrder (= 4)", Z4ClosureOrder);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Real (= 1, drift check on i⁴ = 1)", MemoryLoopClosure.Real);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Imaginary (= 0, drift check)", MemoryLoopClosure.Imaginary);
            yield return new InspectableNode("F-formula anchor",
                summary: "docs/proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md");
            yield return new InspectableNode("typed Claim anchor",
                summary: "compute/RCPsiSquared.Core/SymmetryFamily/F93DetuningAntiPalindromicSpectralInvariance.cs (the F93 closed form as a typed Claim)");
            yield return new InspectableNode("sister claims",
                summary: "F91 (γ_l, per-site Z-dephasing axis), F92 (J_b, per-bond XY coupling axis); three parameter-side instances of one Z₄ rotational structure");
            yield return new InspectableNode("inheritance",
                summary: "the Z₄ closure on the h-axis is not introduced by F93; it inherits from the Pi2-Z₄ memory loop i⁴ = 1, which is the same closure that lives on NinetyDegreeMirrorMemoryClaim's operator-quaternion side");
        }
    }
}
