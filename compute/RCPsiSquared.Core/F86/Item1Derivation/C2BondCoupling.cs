using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum): per-bond projected matrix V_b = B† · M_H_per_bond[b] · B
/// in the 4-mode basis B = [|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩] (PROOF_F86_QPEAK Item 1).
///
/// <para>Houses the three V_b sub-blocks of the 4×4 projection. <b>This task (B1)</b> lands
/// the <i>probe-block</i> (top-left 2×2): entries V_b[α, β] for α, β ∈ {0, 1}, i.e. the
/// projection onto span{|c_1⟩, |c_3⟩} where the Dicke probe lives. B2 will extend this to
/// the cross-block (probe ↔ SVD-top), and B3 will extend to the SVD-block plus an
/// anti-Hermiticity guard.</para>
///
/// <para><b>Two structural facts:</b></para>
/// <list type="bullet">
///   <item><b>F73 sum-rule (generalised):</b> M_H_total = Σ_b M_H_per_bond[b] is diagonal in
///   the channel-uniform basis ⇒ Σ_b V_b[α, β] = 0 for α ≠ β. Each individual V_b[α, β] is
///   generically off-diagonal; the sum vanishes by the bond-uniform-J identity.</item>
///   <item><b>Per-bond entries:</b> Computed by direct projection
///   <c>V_b[α, β] = ⟨c_α | M_H_per_bond[b] | c_β⟩</c> using A2's cached |c_1⟩, |c_3⟩
///   vectors via <see cref="C2ChannelUniformAnalytical"/>. The bond-flip Hamiltonian on
///   (P_1, P_2) gives M_H_per_bond[b] entries with explicit structure — every component is
///   already derived (M_H_per_bond from <see cref="BlockLDecomposition"/>'s explicit Pauli
///   algebra, |c_α⟩ closed-form per A2), so the projection inherits Tier 1 derived.</item>
/// </list>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2).</para>
/// </summary>
public sealed class C2BondCoupling : Claim
{
    public CoherenceBlock Block { get; }

    /// <summary>Composition: the cached |c_1⟩, |c_3⟩ analytical vectors. Reusing A2 keeps
    /// the projection formula honest — the closed forms live in this object's vectors,
    /// V_b[α, β] is just ⟨c_α | M_H_per_bond[b] | c_β⟩ on top.</summary>
    public C2ChannelUniformAnalytical ChannelUniform { get; }

    /// <summary>Public factory: validates c=2, builds the channel-uniform analytical
    /// composition, then constructs the instance. Mirrors the static-factory pattern from
    /// A3 (commit <c>2f62331</c>) so future B2/B3 additions and Tier promotion stay
    /// type-safe.</summary>
    public static C2BondCoupling Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BondCoupling applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var channelUniform = new C2ChannelUniformAnalytical(block);
        return new C2BondCoupling(block, channelUniform);
    }

    private C2BondCoupling(CoherenceBlock block, C2ChannelUniformAnalytical channelUniform)
        : base("c=2 bond coupling V_b in 4-mode basis",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1 (c=2)")
    {
        Block = block;
        ChannelUniform = channelUniform;
    }

    /// <summary>V_b[α, β] for α, β ∈ {0, 1} = probe-block of the 4-mode projected matrix.
    ///
    /// <para>Direct projection <c>⟨c_α | M_H_per_bond[b] | c_β⟩</c> using A2's cached
    /// |c_1⟩, |c_3⟩ vectors. Index 0 → |c_1⟩ (HD=1 channel-uniform), index 1 → |c_3⟩
    /// (HD=3 channel-uniform). Per-entry Tier 1 derived: M_H_per_bond[b] is built by
    /// <see cref="BlockLDecomposition"/> from explicit Pauli algebra; |c_α⟩ are closed-form
    /// per A2.</para>
    /// </summary>
    public Complex ProbeBlockEntry(int bond, int alpha, int beta)
    {
        if (bond < 0 || bond >= Block.NumBonds)
            throw new ArgumentOutOfRangeException(nameof(bond),
                $"bond must be in [0, {Block.NumBonds - 1}]; got {bond}.");
        if (alpha < 0 || alpha > 1)
            throw new ArgumentOutOfRangeException(nameof(alpha),
                $"alpha (probe-block row) must be 0 or 1; got {alpha}.");
        if (beta < 0 || beta > 1)
            throw new ArgumentOutOfRangeException(nameof(beta),
                $"beta (probe-block column) must be 0 or 1; got {beta}.");

        var cAlpha = ProbeVector(alpha);
        var cBeta = ProbeVector(beta);
        var mh = Block.Decomposition.MhPerBond[bond];
        // ⟨c_α | M | c_β⟩ = c_α.Conjugate() · (M · c_β)
        return cAlpha.Conjugate() * (mh * cBeta);
    }

    private ComplexVector ProbeVector(int index) => index switch
    {
        0 => ChannelUniform.C1Vector,
        1 => ChannelUniform.C3Vector,
        _ => throw new ArgumentOutOfRangeException(nameof(index)),
    };

    public override string DisplayName =>
        $"c=2 bond coupling V_b probe-block (N={Block.N}, bonds={Block.NumBonds})";

    public override string Summary =>
        $"probe-block 2×2 entries V_b[α,β] for α,β∈{{0,1}}, F73 sum-rule Σ_b V_b[α,β]=0 for α≠β " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Block.NumBonds.ToString());
            yield return ChannelUniform;
        }
    }
}
