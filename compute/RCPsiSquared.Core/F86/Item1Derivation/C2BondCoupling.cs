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
///   <item><b>B3 (SVD-block, bottom-right 2×2):</b> entries V_b[j, k] for j, k ∈ {2, 3}
///   via <see cref="SvdBlockEntry"/>; same Tier 2 inheritance from A3 as the cross-block.
///   Plus <see cref="AsMatrix"/> assembling the full 4×4 V_b from the three sub-blocks
///   (top-left probe, top-right cross, bottom-right SVD, bottom-left = -conj(top-right)
///   by anti-Hermiticity of M_H_per_bond). The anti-Hermiticity guard test
///   (<c>‖V_b + V_b†‖_F &lt; 1e-10</c> across all bonds, all N=5..8) catches sign drift
///   within each Hermitian-paired sub-block (probe-block, SVD-block); cross-block
///   consistency requires the parallel <c>AsMatrix_FullVb_MatchesFourModeEffective</c>
///   check against <c>B† M B</c>.</item>
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
///   <c>V_b[α, j] = ⟨c_α | M_H_per_bond[b] | x_j⟩</c>; SVD-block via
///   <c>V_b[j, k] = ⟨x_j | M_H_per_bond[b] | x_k⟩</c> where x_2 = |u_0⟩, x_3 = |v_0⟩.</item>
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

    /// <summary>V_b[j, k] for j, k ∈ {2, 3} = SVD-block (bottom-right 2×2) of the 4-mode
    /// projected matrix.
    ///
    /// <para>Direct projection <c>⟨x_j | M_H_per_bond[b] | x_k⟩</c> where x_2 = |u_0⟩
    /// (HD=1 SVD-top, lifted to full block) and x_3 = |v_0⟩ (HD=3 SVD-top). Same
    /// projection-formula pattern as <see cref="ProbeBlockEntry"/> and
    /// <see cref="CrossBlockEntry"/>.</para>
    ///
    /// <para><b>Per-entry Tier:</b> inherits Tier 2 from A3's |u_0⟩, |v_0⟩ obstruction —
    /// same source as the cross-block. Matches
    /// <c>FourModeEffective.MhPerBondEff[b][j, k]</c> at 1e-12 (library-internal
    /// consistency).</para>
    /// </summary>
    public Complex SvdBlockEntry(int bond, int j, int k)
    {
        ValidateBond(bond);
        if (j < 2 || j > 3)
            throw new ArgumentOutOfRangeException(nameof(j),
                $"j (SVD-block row, SVD-top side) must be 2 or 3; got {j}.");
        if (k < 2 || k > 3)
            throw new ArgumentOutOfRangeException(nameof(k),
                $"k (SVD-block column, SVD-top side) must be 2 or 3; got {k}.");

        var xJ = SvdTopVector(j);
        var xK = SvdTopVector(k);
        var mh = Block.Decomposition.MhPerBond[bond];
        // ⟨x_j | M | x_k⟩ = x_j.Conjugate() · (M · x_k)
        return xJ.Conjugate() * (mh * xK);
    }

    /// <summary>D_eff = B† · D · B in the 4-mode basis at c=2:
    /// <c>diag(−2γ₀, −6γ₀, −2γ₀, −6γ₀)</c>.
    ///
    /// <para>Off-diagonal entries are exactly zero because each of the four basis vectors
    /// <c>|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩</c> lives entirely inside one HD subspace (HD=1 for
    /// indices 0 and 2; HD=3 for indices 1 and 3), and D is HD-diagonal with entries
    /// <c>D[i, i] = −2γ₀ · HD(p, q)</c>. Hence
    /// <c>D_eff[i, j] = ⟨b_i | D | b_j⟩ = (−2γ₀ · HD_i) · ⟨b_i | b_j⟩</c>, which is
    /// <c>−2γ₀ · HD_i · δ_{ij}</c> by orthonormality of B (verified by
    /// <see cref="FourModeBasis.OffOrthonormalityResidual"/>).</para>
    ///
    /// <para><b>Tier 1 derived (structural):</b> closed form from the F73 generalisation
    /// (D HD-diagonal) plus the basis-construction rule (one HD subspace per index). Unlike
    /// the cross/SVD blocks which inherit A3's Tier 2 obstruction, D_eff has no such
    /// dependence — the SVD direction inside the HD subspace is irrelevant because D is
    /// proportional to the identity inside each HD subspace. The class-level Tier remains
    /// Tier 2 verified (governed by the weakest link, B2 cross-block); D_eff being clean
    /// does not promote the class.</para>
    /// </summary>
    public Matrix<Complex> DEffDiagonal()
    {
        var m = Matrix<Complex>.Build.Dense(4, 4);
        double gamma0 = Block.GammaZero;
        // |c_1⟩ (index 0) and |u_0⟩ (index 2) are in the HD=1 subspace ⇒ −2γ₀ · 1 = −2γ₀
        // |c_3⟩ (index 1) and |v_0⟩ (index 3) are in the HD=3 subspace ⇒ −2γ₀ · 3 = −6γ₀
        m[0, 0] = new Complex(-2.0 * gamma0, 0.0);
        m[1, 1] = new Complex(-6.0 * gamma0, 0.0);
        m[2, 2] = new Complex(-2.0 * gamma0, 0.0);
        m[3, 3] = new Complex(-6.0 * gamma0, 0.0);
        return m;
    }

    /// <summary>Diagonal entry <c>D_eff[i, i]</c> for i ∈ {0, 1, 2, 3}: closed form in γ₀.
    /// <c>−2γ₀</c> for indices 0, 2 (HD=1); <c>−6γ₀</c> for indices 1, 3 (HD=3).</summary>
    public double DEffDiagonalEntry(int i)
    {
        if (i < 0 || i > 3)
            throw new ArgumentOutOfRangeException(nameof(i),
                $"D_eff index must be in [0, 3]; got {i}.");
        int hd = (i == 0 || i == 2) ? 1 : 3;
        return -2.0 * Block.GammaZero * hd;
    }

    /// <summary>The full 4×4 V_b matrix at the given bond, assembled from the three
    /// sub-block accessors <see cref="ProbeBlockEntry"/>, <see cref="CrossBlockEntry"/>,
    /// <see cref="SvdBlockEntry"/>, with the bottom-left 2×2 filled in via the
    /// anti-Hermitian relation <c>V_b[j, α] = -conj(V_b[α, j])</c>.
    ///
    /// <para>M_H_per_bond[b] = -i [H_b, ·] is anti-Hermitian as a superoperator, so its
    /// 4×4 projection is also anti-Hermitian: <c>V_b + V_b† = 0</c>. The
    /// <c>Vb_IsAntiHermitian_AcrossAllBondsAndEntries</c> test verifies this at 1e-10
    /// across all bonds and N=5..8. Scope: this guard catches sign drift within each
    /// Hermitian-paired sub-block (probe-block 2×2 and SVD-block 2×2 are each populated
    /// from independent entry computations). It does NOT catch a sign error in
    /// <see cref="CrossBlockEntry"/> on its own, because the bottom-left 2×2 is built
    /// from <c>-conj(top-right)</c> by construction — a buggy <see cref="CrossBlockEntry"/>
    /// yields a buggy mirror, and the cross sector of <c>V_b + V_b†</c> stays zero.
    /// Cross-block consistency is established by the parallel
    /// <c>AsMatrix_FullVb_MatchesFourModeEffective</c> check, which compares entry-by-entry
    /// against <c>B† M B</c> (computed independently by <see cref="FourModeEffective"/>).
    /// Together the two tests cover every sign convention.</para>
    ///
    /// <para>Useful for downstream consumers (Stage C eigenvalue analysis, Stage D
    /// Duhamel evaluation) and for verifying entry-by-entry agreement with
    /// <c>FourModeEffective.MhPerBondEff[b]</c>.</para>
    /// </summary>
    public Matrix<Complex> AsMatrix(int bond)
    {
        ValidateBond(bond);
        var m = Matrix<Complex>.Build.Dense(4, 4);

        // Top-left 2×2: probe-block V_b[α, β] for α, β ∈ {0, 1}.
        for (int alpha = 0; alpha < 2; alpha++)
            for (int beta = 0; beta < 2; beta++)
                m[alpha, beta] = ProbeBlockEntry(bond, alpha, beta);

        // Top-right 2×2: cross-block V_b[α, j] for α ∈ {0, 1}, j ∈ {2, 3}.
        // Bottom-left 2×2: filled by anti-Hermiticity V_b[j, α] = -conj(V_b[α, j]).
        for (int alpha = 0; alpha < 2; alpha++)
            for (int j = 2; j < 4; j++)
            {
                var entry = CrossBlockEntry(bond, alpha, j);
                m[alpha, j] = entry;
                m[j, alpha] = -Complex.Conjugate(entry);
            }

        // Bottom-right 2×2: SVD-block V_b[j, k] for j, k ∈ {2, 3}.
        for (int j = 2; j < 4; j++)
            for (int k = 2; k < 4; k++)
                m[j, k] = SvdBlockEntry(bond, j, k);

        return m;
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
        $"c=2 bond coupling V_b 4×4 (N={Block.N}, bonds={Block.NumBonds})";

    public override string Summary =>
        $"full 4×4 V_b: probe (α,β∈{{0,1}}) + cross (α∈{{0,1}}, j∈{{2,3}}) + SVD (j,k∈{{2,3}}); " +
        $"F73 sum-rule Σ_b V_b[α,β]=0 for α≠β; cross-block Frobenius BondClass-tagged; " +
        $"anti-Hermitian V_b + V_b† = 0 (1e-10) ({Tier.Label()})";

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
