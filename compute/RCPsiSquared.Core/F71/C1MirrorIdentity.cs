using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F71;

/// <summary>F71 original statement (kinematic, Tier 1 derived): the closure-breaking
/// coefficient c₁ on a uniform N-qubit chain with reflection-symmetric initial state is
/// mirror-symmetric across bonds.
///
/// <code>
///     c₁(N, b, ρ₀)  =  c₁(N, N−2−b, ρ₀)         for all b ∈ {0, …, N−2}
/// </code>
///
/// <para>Independent components: at most ⌈(N−1)/2⌉ instead of N−1. If N is even, the central
/// bond b = (N−2)/2 is self-paired (mirror image is itself); if N is odd, all N−1 bonds pair
/// up in (N−1)/2 disjoint pairs.</para>
///
/// <para>Proof skeleton (cf. <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c>): R commutes with
/// the uniform Liouvillian, [L_A, R_sup] = 0; under R, bond b maps to bond N−2−b via
/// R · T_b · R = T_{N−2−b}; per-site purity is quadratic in ρ so any phase squares away;
/// hence P_B(b, i, t) = P_B(N−2−b, N−1−i, t), and summing ln(α_i) yields the bond-mirror
/// identity for c₁.</para>
///
/// <para>Validity: any [H, R] = 0 (uniform J on a symmetric graph), any [D, R_sup] = 0
/// (uniform or R-symmetric dephasing), any reflection-symmetric initial state. Breaks for
/// non-uniform J_b ≠ J_{N−2−b}, non-uniform γ_i ≠ γ_{N−1−i}, asymmetric ρ₀.</para>
/// </summary>
public sealed class C1MirrorIdentity : Claim
{
    public C1MirrorIdentity()
        : base("F71 c₁ mirror identity: c₁(b) = c₁(N−2−b)",
               Tier.Tier1Derived,
               $"{F71Anchors.Formulas} + {F71Anchors.Proof}")
    { }

    public override string DisplayName => "F71: c₁(N, b, ρ₀) = c₁(N, N−2−b, ρ₀)";

    public override string Summary =>
        "closure-breaking coefficient c₁ is bond-mirror symmetric on uniform reflection-symmetric chains; ⌈(N−1)/2⌉ independent components";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("kinematic origin",
                summary: "[L_A, R_sup] = 0 + R · T_b · R = T_{N−2−b} + per-site purity quadratic in ρ");
            yield return new InspectableNode("validity",
                summary: "uniform J on R-symmetric graph; uniform or R-symmetric dephasing; ρ₀ R-symmetric in per-site purities");
            yield return new InspectableNode("breaks for",
                summary: "non-uniform J_b ≠ J_{N−2−b}; non-uniform γ_i ≠ γ_{N−1−i}; asymmetric ρ₀");
            yield return new InspectableNode("verified",
                summary: "N=3..6 for ψ_1+vac and ψ_2+vac; residuals < 10⁻⁹ (eq021_obc_sine_basis.py)");
        }
    }
}
