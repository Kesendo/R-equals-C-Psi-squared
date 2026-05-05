using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum): per-bond projected matrix V_b = B† · M_H_per_bond[b] · B
/// in the 4-mode basis B = [|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩] (PROOF_F86_QPEAK Item 1).
///
/// <para>Houses the three V_b sub-blocks of the 4×4 projection.</para>
/// <list type="bullet">
///   <item><b>B1 (probe-block, top-left 2×2):</b> entries V_b[α, β] for α, β ∈ {0, 1} via
///   <see cref="ProbeBlockEntry"/>. Structurally Tier 1 derived (every component already
///   has a closed form: M_H_per_bond[b] from explicit Pauli algebra; |c_α⟩ analytical per A2).</item>
///   <item><b>B2 (cross-block, off-diagonal 2×2):</b> entries V_b[α, j] for α ∈ {0, 1},
///   j ∈ {2, 3} via <see cref="CrossBlockEntry"/>, plus the
///   <see cref="CrossBlockWitnesses"/> collection tagged by <see cref="BondClass"/>.
///   This is the bond-position-dependent block: its Frobenius norm splits Endpoint vs
///   Interior and is the algebraic seed of the F86 universal-shape statement.</item>
///   <item><b>B3 (SVD-block, bottom-right 2×2):</b> reserved.</item>
/// </list>
///
/// <para><b>Class-level Tier: <c>Tier2Verified</c></b> — reflects the weakest link. The
/// probe-block (B1) is structurally Tier 1 derived in isolation, but the cross-block (B2)
/// inherits A3's Tier 2 obstruction: |u_0⟩ and |v_0⟩ are taken from
/// <see cref="C2InterChannelAnalytical"/>, and at even N the top singular value σ_0 of
/// V_inter is exactly degenerate (deg = 2 at N=6, N=8). Inside that 2D top eigenspace the
/// "top SVD vector" is library-dependent (MathNet vs numpy disagree), so the cross-block
/// entries match <c>FourModeEffective.MhPerBondEff[b][α, j]</c> at 1e-12 (library-internal
/// consistency) but are not closed-form-in-N. The Tier 1 promotion path noted in A3's
/// PendingDerivationNote is to lift to projector-overlap onto the 2D top eigenspace, which
/// IS unique. Until then, the per-(N, b, BondClass) numerical witnesses in
/// <see cref="CrossBlockWitnesses"/> are the OOP carry of the partial finding.</para>
///
/// <para><b>Two structural facts:</b></para>
/// <list type="bullet">
///   <item><b>F73 sum-rule (generalised):</b> M_H_total = Σ_b M_H_per_bond[b] is diagonal in
///   the channel-uniform basis ⇒ Σ_b V_b[α, β] = 0 for α ≠ β. Each individual V_b[α, β] is
///   generically off-diagonal; the sum vanishes by the bond-uniform-J identity.</item>
///   <item><b>Per-bond entries:</b> Probe-block via direct projection
///   <c>V_b[α, β] = ⟨c_α | M_H_per_bond[b] | c_β⟩</c>; cross-block via
///   <c>V_b[α, j] = ⟨c_α | M_H_per_bond[b] | x_j⟩</c> where x_2 = |u_0⟩, x_3 = |v_0⟩.</item>
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

    /// <summary>Composition: the SVD-top |u_0⟩, |v_0⟩ inter-channel vectors. Reusing A3
    /// makes the cross-block computation flow through a single source: the A3 Resolution.
    /// Class-level Tier inherits A3's <see cref="C2InterChannelAnalytical.Tier"/> (currently
    /// Tier2Verified due to even-N σ_0 degeneracy — see A3's PendingDerivationNote).</summary>
    public C2InterChannelAnalytical InterChannel { get; }

    /// <summary>One <see cref="CrossBlockWitness"/> per bond, in bond-index order, tagged
    /// with <see cref="BondClass"/>. Captures the cross-block 2×2 entries V_b[α, j]
    /// (α ∈ {0, 1}, j ∈ {2, 3}) plus the Frobenius norm — the magnitude indicator that
    /// drives the Endpoint vs Interior split.</summary>
    public IReadOnlyList<CrossBlockWitness> CrossBlockWitnesses { get; }

    /// <summary>Public factory: validates c=2, builds the channel-uniform + SVD-top
    /// compositions, computes per-bond cross-block witnesses, then constructs the instance.
    /// Mirrors the static-factory pattern from A3 (commit <c>2f62331</c>) so future
    /// promotion to Tier 1 stays type-safe.</summary>
    public static C2BondCoupling Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BondCoupling applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var channelUniform = new C2ChannelUniformAnalytical(block);
        var interChannel = C2InterChannelAnalytical.Build(block);
        var witnesses = BuildCrossBlockWitnesses(block, channelUniform, interChannel);
        return new C2BondCoupling(block, channelUniform, interChannel, witnesses);
    }

    private C2BondCoupling(
        CoherenceBlock block,
        C2ChannelUniformAnalytical channelUniform,
        C2InterChannelAnalytical interChannel,
        IReadOnlyList<CrossBlockWitness> witnesses)
        : base("c=2 bond coupling V_b in 4-mode basis",
               // Class-level Tier reflects the weakest link: the cross-block (B2) inherits
               // A3's Tier (Tier2Verified) because |u_0⟩, |v_0⟩ are not yet closed-form-in-N.
               // The probe-block (B1) is structurally Tier1Derived in isolation; documented
               // in the class-level XML doc above.
               interChannel.Tier,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1 (c=2)")
    {
        Block = block;
        ChannelUniform = channelUniform;
        InterChannel = interChannel;
        CrossBlockWitnesses = witnesses;
    }

    /// <summary>V_b[α, β] for α, β ∈ {0, 1} = probe-block of the 4-mode projected matrix.
    ///
    /// <para>Direct projection <c>⟨c_α | M_H_per_bond[b] | c_β⟩</c> using A2's cached
    /// |c_1⟩, |c_3⟩ vectors. Index 0 → |c_1⟩ (HD=1 channel-uniform), index 1 → |c_3⟩
    /// (HD=3 channel-uniform). Per-entry structurally Tier 1 derived: M_H_per_bond[b] is
    /// built by <see cref="BlockLDecomposition"/> from explicit Pauli algebra; |c_α⟩ are
    /// closed-form per A2.</para>
    /// </summary>
    public Complex ProbeBlockEntry(int bond, int alpha, int beta)
    {
        ValidateBond(bond);
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

    /// <summary>V_b[α, j] for α ∈ {0, 1}, j ∈ {2, 3} = cross-block (probe ↔ EP-partner).
    ///
    /// <para>Direct projection <c>⟨c_α | M_H_per_bond[b] | x_j⟩</c> where x_2 = |u_0⟩
    /// (HD=1 SVD-top, lifted to full block) and x_3 = |v_0⟩ (HD=3 SVD-top). Index 0 → |c_1⟩,
    /// index 1 → |c_3⟩.</para>
    ///
    /// <para><b>Per-entry Tier:</b> matches <c>FourModeEffective.MhPerBondEff[b][α, j]</c>
    /// at 1e-12 because both compositions feed through the same MathNet SVD output. This is
    /// library-internal consistency, not a Tier 1 closed form. The class-level Tier
    /// inherits A3's <see cref="C2InterChannelAnalytical.Tier"/>; until A3 is promoted to
    /// projector-overlap, these per-(N, b) numerics are the OOP carry.</para>
    /// </summary>
    public Complex CrossBlockEntry(int bond, int alpha, int j)
    {
        ValidateBond(bond);
        if (alpha < 0 || alpha > 1)
            throw new ArgumentOutOfRangeException(nameof(alpha),
                $"alpha (cross-block row, probe side) must be 0 or 1; got {alpha}.");
        if (j < 2 || j > 3)
            throw new ArgumentOutOfRangeException(nameof(j),
                $"j (cross-block column, SVD-top side) must be 2 or 3; got {j}.");

        var cAlpha = ProbeVector(alpha);
        var xJ = SvdTopVector(j);
        var mh = Block.Decomposition.MhPerBond[bond];
        // ⟨c_α | M | x_j⟩ = c_α.Conjugate() · (M · x_j)
        return cAlpha.Conjugate() * (mh * xJ);
    }

    private void ValidateBond(int bond)
    {
        if (bond < 0 || bond >= Block.NumBonds)
            throw new ArgumentOutOfRangeException(nameof(bond),
                $"bond must be in [0, {Block.NumBonds - 1}]; got {bond}.");
    }

    private ComplexVector ProbeVector(int index) => index switch
    {
        0 => ChannelUniform.C1Vector,
        1 => ChannelUniform.C3Vector,
        _ => throw new ArgumentOutOfRangeException(nameof(index)),
    };

    private ComplexVector SvdTopVector(int index) => index switch
    {
        2 => InterChannel.U0,
        3 => InterChannel.V0,
        _ => throw new ArgumentOutOfRangeException(nameof(index)),
    };

    private static IReadOnlyList<CrossBlockWitness> BuildCrossBlockWitnesses(
        CoherenceBlock block,
        C2ChannelUniformAnalytical channelUniform,
        C2InterChannelAnalytical interChannel)
    {
        int numBonds = block.NumBonds;
        var list = new CrossBlockWitness[numBonds];
        for (int b = 0; b < numBonds; b++)
        {
            var mh = block.Decomposition.MhPerBond[b];
            // V_b[α, j] = ⟨c_α | M_h_per_bond[b] | x_j⟩
            var c1 = channelUniform.C1Vector;
            var c3 = channelUniform.C3Vector;
            var u0 = interChannel.U0;
            var v0 = interChannel.V0;

            var mhU0 = mh * u0;
            var mhV0 = mh * v0;
            var entryC1U0 = c1.Conjugate() * mhU0;  // V_b[0, 2]
            var entryC1V0 = c1.Conjugate() * mhV0;  // V_b[0, 3]
            var entryC3U0 = c3.Conjugate() * mhU0;  // V_b[1, 2]
            var entryC3V0 = c3.Conjugate() * mhV0;  // V_b[1, 3]

            // Frobenius norm of the 2×2 cross-block as a magnitude indicator.
            double frob = Math.Sqrt(
                entryC1U0.Magnitude * entryC1U0.Magnitude +
                entryC1V0.Magnitude * entryC1V0.Magnitude +
                entryC3U0.Magnitude * entryC3U0.Magnitude +
                entryC3V0.Magnitude * entryC3V0.Magnitude);

            var bondClass = (b == 0 || b == numBonds - 1)
                ? BondClass.Endpoint
                : BondClass.Interior;

            list[b] = new CrossBlockWitness(
                Bond: b,
                BondClass: bondClass,
                EntryC1U0: entryC1U0,
                EntryC1V0: entryC1V0,
                EntryC3U0: entryC3U0,
                EntryC3V0: entryC3V0,
                FrobeniusNorm: frob);
        }
        return list;
    }

    public override string DisplayName =>
        $"c=2 bond coupling V_b probe + cross block (N={Block.N}, bonds={Block.NumBonds})";

    public override string Summary =>
        $"probe-block 2×2 V_b[α,β] (α,β∈{{0,1}}) + cross-block 2×2 V_b[α,j] (α∈{{0,1}}, j∈{{2,3}}); " +
        $"F73 sum-rule Σ_b V_b[α,β]=0 for α≠β; cross-block Frobenius BondClass-tagged " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Block.NumBonds.ToString());
            yield return ChannelUniform;
            yield return InterChannel;
            yield return InspectableNode.Group(
                "CrossBlockWitnesses",
                CrossBlockWitnesses,
                CrossBlockWitnesses.Count);
        }
    }
}

/// <summary>OOP carry of B2's per-(N, bond) cross-block numerical findings, tagged with
/// <see cref="BondClass"/>. One entry per bond; visible in the inspection tree under
/// <see cref="C2BondCoupling.CrossBlockWitnesses"/>.
///
/// <para><b>Why a witness, not a closed form:</b> the cross-block entries
/// V_b[α, j] = ⟨c_α | M_H_per_bond[b] | x_j⟩ depend on |u_0⟩, |v_0⟩, which are A3
/// Tier 2 (numerical SVD, library-tiebreaker-dependent at even N). The Frobenius norm
/// over the 2×2 cross-block is a magnitude indicator that splits Endpoint vs Interior
/// at c=2 (empirical anchor: HWHM_left/Q_peak Endpoint=0.7728 > Interior=0.7506,
/// EQ-022 (b1) 2026-05-02). Future Tier 1 promotion (per A3 PendingDerivationNote: lift
/// to projector-overlap onto the 2D top eigenspace) would replace these per-(N, b)
/// witnesses with a closed-form ratio.</para>
/// </summary>
public sealed record CrossBlockWitness(
    int Bond,
    BondClass BondClass,
    Complex EntryC1U0,    // V_b[0, 2] = ⟨c_1 | M | u_0⟩
    Complex EntryC1V0,    // V_b[0, 3] = ⟨c_1 | M | v_0⟩
    Complex EntryC3U0,    // V_b[1, 2] = ⟨c_3 | M | u_0⟩
    Complex EntryC3V0,    // V_b[1, 3] = ⟨c_3 | M | v_0⟩
    double FrobeniusNorm  // ‖cross-block 2×2‖_F
) : IInspectable
{
    public string DisplayName => $"cross-block witness b={Bond} ({BondClass})";

    public string Summary => $"‖V_b cross‖_F = {FrobeniusNorm:F4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("‖V_b cross‖_F", FrobeniusNorm);
            yield return new InspectableNode("V_b[0, 2] (⟨c_1 | M | u_0⟩)", summary: EntryC1U0.ToString("G6"));
            yield return new InspectableNode("V_b[0, 3] (⟨c_1 | M | v_0⟩)", summary: EntryC1V0.ToString("G6"));
            yield return new InspectableNode("V_b[1, 2] (⟨c_3 | M | u_0⟩)", summary: EntryC3U0.ToString("G6"));
            yield return new InspectableNode("V_b[1, 3] (⟨c_3 | M | v_0⟩)", summary: EntryC3V0.ToString("G6"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
