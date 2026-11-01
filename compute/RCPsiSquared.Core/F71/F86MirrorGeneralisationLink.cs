using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F71;

/// <summary>F71 → F86 generalisation cross-reference (Tier 1 derived 2026-05-03): the same
/// kinematic mirror argument that gives c₁(b) = c₁(N−2−b) (closure-breaking coefficient)
/// extends to the F86 K_CC_pr per-bond observable on the (n, n+1) popcount coherence block.
///
/// <code>
///     Q_peak(b)  =  Q_peak(N−2−b)              (bit-exactly, all c, N)
/// </code>
///
/// <para>Mechanism: the F86 observable <c>K_b(Q, t) = 2·Re ⟨ρ(t)| S_kernel | ∂ρ/∂J_b ⟩</c>
/// is invariant under R. The Z-dephasing dissipator L_D, the uniform-J Hamiltonian H_xy, the
/// Dicke probe, and the spatial-sum kernel S are R-symmetric; the bond-flip transforms as
/// <c>R · ∂L/∂J_b · R⁻¹ = ∂L/∂J_{N−2−b}</c>. Hence K_b(Q, t) = K_{N−2−b}(Q, t) as functions
/// of (Q, t), and their argmax-Q values coincide.</para>
///
/// <para>Per-F71-orbit substructure: even within the bond-mirror grouping, central
/// (self-paired) bonds give a different Q_peak from flanking bonds. F71 gives the pairing,
/// not the value: e.g. c=2 N=6 central b=2 → 1.436 vs flanking b=1, b=3 → 1.616.</para>
///
/// <para>This is the link between F71 (kinematic, c₁ vac+SE) and F86 (EP-rotation universality,
/// K_CC_pr on (n, n+1)); cf. <see cref="C1MirrorIdentity"/> and the F86 typed knowledge base.</para>
/// </summary>
public sealed class F86MirrorGeneralisationLink : Claim
{
    public F86MirrorGeneralisationLink()
        : base("F71 → F86 mirror generalisation: Q_peak(b) = Q_peak(N−2−b)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F71 generalisation paragraph + docs/proofs/PROOF_F86C_F71_MIRROR.md Statement 3")
    { }

    public override string DisplayName => "F71 → F86: Q_peak(b) = Q_peak(N−2−b) bit-exactly";

    public override string Summary =>
        "F71 spatial mirror lifts from c₁ (vac+SE) to F86 K_CC_pr per-bond Q_peak on (n, n+1) coherence block; bit-exact across c=2 N=5..7 and c=3 N=5..6";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F86 observable",
                summary: "K_b(Q, t) = 2·Re ⟨ρ(t)| S_kernel | ∂ρ/∂J_b ⟩");
            yield return new InspectableNode("R-invariant components",
                summary: "L_D (Z-dephasing), H_xy (uniform J), Dicke probe, S_kernel, all commute with R");
            yield return new InspectableNode("bond-flip transformation",
                summary: "R · ∂L/∂J_b · R⁻¹ = ∂L/∂J_{N−2−b} → K_b(Q, t) = K_{N−2−b}(Q, t)");
            yield return new InspectableNode("substructure refinement",
                summary: "F71 pairs bonds; per-orbit Q_peak values differ (central self-paired ≠ flanking)");
            yield return new InspectableNode("relation to F71 base",
                summary: "same kinematic argument as C1MirrorIdentity (c₁ on vac+SE), lifted to F86 on (n, n+1)");
        }
    }
}
