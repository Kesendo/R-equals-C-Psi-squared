using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F92's anti-palindromic J orbit (J_b + J_{N−2−b} = 2·J_avg ∀b) is the
/// parameter-side J-axis instance of the Pi2-Z₄ rotational structure. The
/// 90°-rotation J ↦ 2·J_avg − F71(J) is one of the four Z₄ elements; closing
/// under i⁴ = 1 returns to identity. The Z₄ generator on parameters is
/// structurally the same Z₄ that
/// <see cref="NinetyDegreeMirrorMemoryClaim"/> types on the operator-quaternion
/// side; F92 inherits that algebra rather than introducing a new orbit closure.
///
/// <para>Parameter axis: J_b (per-bond XY coupling). Anti-palindromic orbit:
/// J_b + J_{N−2−b} = 2·J_avg for all b ∈ {0, ..., N−2} (bonds, not sites; F71
/// maps bond b ↔ N−2−b). Closure under the 90°-rotation: applied four times
/// R_{90}^4 = identity, which is exactly the i⁴ = 1 closure on the Pi2-Z₄
/// memory loop. F92 is the J-axis instance of one rotational structure; F91 is
/// the γ-axis instance, F93 the h-axis instance.</para>
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
public sealed class F92Pi2Inheritance : Claim
{
    private readonly Pi2I4MemoryLoopClaim _memoryLoop;

    /// <summary>The Z₄ closure order of the 90°-rotation R_{90}: J_b ↦ 2·J_avg − J_{N−2−b}.
    /// Four applications return to identity, matching <see cref="Pi2I4MemoryLoopClaim.ClosureOrder"/>.</summary>
    public int Z4ClosureOrder => Pi2I4MemoryLoopClaim.ClosureOrder;

    /// <summary>Live drift check: i^4 = 1 exactly on the parent Pi2-Z₄ memory loop.
    /// F92's parameter-side 90°-rotation closes at the same order as the operator-side
    /// quaternion algebra.</summary>
    public Complex MemoryLoopClosure => _memoryLoop.MemoryClosure();

    /// <summary>Verbal name of the parameter axis F92 lives on: per-bond XY coupling
    /// J_b (bonds, not sites). Distinguishes F92 from its sister claims F91 (per-site
    /// Z-dephasing γ_l) and F93 (per-site longitudinal Z-detuning h_l).</summary>
    public string ParameterAxis => "J_b (per-bond XY coupling)";

    public F92Pi2Inheritance(Pi2I4MemoryLoopClaim memoryLoop)
        : base("F92 anti-palindromic J orbit (J_b + J_{N-2-b} = 2·J_avg) inherits from Pi2-Z₄ structure (i⁴ = 1 closure); parameter-side J-axis instance of the same Z₄ that NinetyDegreeMirrorMemoryClaim types on the operator-quaternion side",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md + " +
               "compute/RCPsiSquared.Core/SymmetryFamily/F92BondAntiPalindromicJSpectralInvariance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (NinetyDegreeMirrorMemoryClaim)")
    {
        _memoryLoop = memoryLoop ?? throw new ArgumentNullException(nameof(memoryLoop));
    }

    public override string DisplayName =>
        "F92 J-Z₄ anti-palindromic orbit as Pi2-Foundation parameter-side instance";

    public override string Summary =>
        $"J_b + J_{{N-2-b}} = 2·J_avg orbit closes under 90°-rotation R_{{90}}: J_b ↦ 2·J_avg − J_{{N-2-b}} at order {Z4ClosureOrder} (i⁴ = 1); parameter-side J-axis instance of the Pi2-Z₄ operator-quaternion structure ({Tier.Label()})";

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
            yield return new InspectableNode("Z₄ generator on parameters",
                summary: "the 90°-rotation J ↦ 2·J_avg − F71(J) is one of four Z₄ elements; structurally the same Z₄ that NinetyDegreeMirrorMemoryClaim types on the operator-quaternion side");
            yield return InspectableNode.RealScalar("Z4ClosureOrder (= 4)", Z4ClosureOrder);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Real (= 1, drift check on i⁴ = 1)", MemoryLoopClosure.Real);
            yield return InspectableNode.RealScalar("MemoryLoopClosure.Imaginary (= 0, drift check)", MemoryLoopClosure.Imaginary);
            yield return new InspectableNode("F-formula anchor",
                summary: "docs/proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md");
            yield return new InspectableNode("typed Claim anchor",
                summary: "compute/RCPsiSquared.Core/SymmetryFamily/F92BondAntiPalindromicJSpectralInvariance.cs (the F92 closed form as a typed Claim)");
            yield return new InspectableNode("sister claims",
                summary: "F91 (γ_l, per-site Z-dephasing axis), F93 (h_l, longitudinal Z-detuning axis); three parameter-side instances of one Z₄ rotational structure");
            yield return new InspectableNode("inheritance",
                summary: "the Z₄ closure on the J-axis is not introduced by F92; it inherits from the Pi2-Z₄ memory loop i⁴ = 1, which is the same closure that lives on NinetyDegreeMirrorMemoryClaim's operator-quaternion side");
        }
    }
}
