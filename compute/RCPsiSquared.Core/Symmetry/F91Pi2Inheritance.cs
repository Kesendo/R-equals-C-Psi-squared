using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F91's anti-palindromic γ orbit (γ_l + γ_{N−1−l} = 2·γ_avg ∀l) is the
/// parameter-side γ-axis instance of the Pi2-Z₄ rotational structure. The
/// 90°-rotation γ ↦ 2·γ_avg − F71(γ) is one of the four Z₄ elements; closing
/// under i⁴ = 1 returns to identity. The Z₄ generator on parameters is
/// structurally the same Z₄ that
/// <see cref="NinetyDegreeMirrorMemoryClaim"/> types on the operator-quaternion
/// side; F91 inherits that algebra rather than introducing a new orbit closure.
///
/// <para>Parameter axis: γ_l (per-site Z-dephasing rate). Anti-palindromic
/// orbit: γ_l + γ_{N−1−l} = 2·γ_avg for all l ∈ {0, ..., N−1}; for odd N the
/// middle site γ_{(N−1)/2} = γ_avg is forced by self-pairing. Closure under the
/// 90°-rotation: applied four times R_{90}^4 = identity, which is exactly the
/// i⁴ = 1 closure on the Pi2-Z₄ memory loop. F91 is the γ-axis instance of one
/// rotational structure; F92 is the J-axis instance, F93 the h-axis instance.</para>
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
public sealed class F91Pi2Inheritance : Claim
{
    private readonly Pi2I4MemoryLoopClaim _memoryLoop;

    /// <summary>The Z₄ closure order of the 90°-rotation R_{90}: γ_l ↦ 2·γ_avg − γ_{N−1−l}.
    /// Four applications return to identity, matching <see cref="Pi2I4MemoryLoopClaim.ClosureOrder"/>.</summary>
    public int Z4ClosureOrder => Pi2I4MemoryLoopClaim.ClosureOrder;

    /// <summary>Live drift check: i^4 = 1 exactly on the parent Pi2-Z₄ memory loop.
    /// F91's parameter-side 90°-rotation closes at the same order as the operator-side
    /// quaternion algebra.</summary>
    public Complex MemoryLoopClosure => _memoryLoop.MemoryClosure();

    /// <summary>Verbal name of the parameter axis F91 lives on: per-site Z-dephasing
    /// rate γ_l. Distinguishes F91 from its sister claims F92 (bond-coupling J_b) and
    /// F93 (longitudinal Z-detuning h_l).</summary>
    public string ParameterAxis => "γ_l (per-site Z-dephasing rate)";

    public F91Pi2Inheritance(Pi2I4MemoryLoopClaim memoryLoop)
        : base("F91 anti-palindromic γ orbit (γ_l + γ_{N-1-l} = 2·γ_avg) inherits from Pi2-Z₄ structure (i⁴ = 1 closure); parameter-side γ-axis instance of the same Z₄ that NinetyDegreeMirrorMemoryClaim types on the operator-quaternion side",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/F71AntiPalindromicGammaSpectralInvariance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (NinetyDegreeMirrorMemoryClaim)")
    {
        _memoryLoop = memoryLoop ?? throw new ArgumentNullException(nameof(memoryLoop));
    }

    public override string DisplayName =>
        "F91 γ-Z₄ anti-palindromic orbit as Pi2-Foundation parameter-side instance";

    public override string Summary =>
        $"γ_l + γ_{{N-1-l}} = 2·γ_avg orbit closes under 90°-rotation R_{{90}}: γ_l ↦ 2·γ_avg − γ_{{N-1-l}} at order {Z4ClosureOrder} (i⁴ = 1); parameter-side γ-axis instance of the Pi2-Z₄ operator-quaternion structure ({Tier.Label()})";

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
            yield return new InspectableNode("Z₄ generator on parameters",
                summary: "the 90°-rotation γ ↦ 2·γ_avg − F71(γ) is one of four Z₄ elements; structurally the same Z₄ that NinetyDegreeMirrorMemoryClaim types on the operator-quaternion side");
            yield return InspectableNode.RealScalar("Z4ClosureOrder (= 4)", Z4ClosureOrder);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Real (= 1, drift check on i⁴ = 1)", MemoryLoopClosure.Real);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Imaginary (= 0, drift check)", MemoryLoopClosure.Imaginary);
            yield return new InspectableNode("F-formula anchor",
                summary: "docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md");
            yield return new InspectableNode("typed Claim anchor",
                summary: "compute/RCPsiSquared.Core/BlockSpectrum/F71AntiPalindromicGammaSpectralInvariance.cs (the F91 closed form as a typed Claim)");
            yield return new InspectableNode("sister claims",
                summary: "F92 (J_b, bond-coupling axis), F93 (h_l, longitudinal Z-detuning axis); three parameter-side instances of one Z₄ rotational structure");
            yield return new InspectableNode("inheritance",
                summary: "the Z₄ closure on the γ-axis is not introduced by F91; it inherits from the Pi2-Z₄ memory loop i⁴ = 1, which is the same closure that lives on NinetyDegreeMirrorMemoryClaim's operator-quaternion side");
        }
    }
}
