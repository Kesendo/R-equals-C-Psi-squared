using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum), <b>Direction (a'') foundation</b>: per-bond SVD-block
/// magnitude squared lifted to the Riesz projector of <see cref="C2InterChannelProjector"/>.
///
/// <para>For each bond b, the library-invariant magnitude is</para>
/// <code>
///     |V_b|²_proj := Tr(P_top^L · M_H_per_bond[b] · P_top^R · M_H_per_bond[b]†)
///                   = ‖P_top^L · M_H_per_bond[b] · P_top^R‖_F²
/// </code>
///
/// <para>This replaces the library-dependent single-entry <c>|V_b[2, 3]|² = |⟨u_0|M_b|v_0⟩|²</c>
/// from <see cref="C2BondCoupling.SvdBlockEntry"/>. At odd N (where σ_0 is non-degenerate,
/// RankTop = 1) the two agree to machine precision. At even N (where σ_0 is exactly
/// degenerate with RankTop = 2) the projector lift sums over BOTH directions of the
/// 2D top eigenspace; the single-entry version is the library tiebreaker, the projected
/// version is unique (Riesz spectral projection theorem).</para>
///
/// <para><b>Tier outcome: Tier1Derived.</b> Inherits Tier1Derived directly from
/// <see cref="C2InterChannelProjector"/>: once the projectors P_top^L and P_top^R are
/// library-invariant (verified at construction by the runtime stability witness), every
/// trace built from them is library-invariant by basic linear algebra. The class-level Tier
/// thus tracks the projector's Tier; if a future even-N case fails the library-stability
/// witness in C2InterChannelProjector, this primitive's Tier drops to Tier2Verified
/// automatically through the Tier inheritance.</para>
///
/// <para><b>What this unlocks (Direction (a'') downstream).</b> The per-bond magnitude
/// |V_b|²_proj is the SVD-block contribution to the K-resonance K_b(Q, t) restricted to
/// the EP-partner subspace. Direction (a'') asks whether HWHM_left/Q_peak − 0.671535
/// (the "lift" above the bare-doubled-PTF Floor in <see cref="C2HwhmRatio"/>) maps to
/// |V_b|²_proj/σ_0² (the dimensionless SVD-block ratio) via a single bond-class-blind
/// function. With this primitive in place that scan is a direct numerical operation; the
/// previous A3 single-vector formulation made the scan unstable at even N because V_b[2, 3]
/// flipped under library tiebreaker rotation. The scan itself is downstream work
/// (<c>C2BondMagnitudeLiftCollapse</c> or similar), not in this primitive.</para>
///
/// <para><b>Sanity invariants checked at construction:</b></para>
/// <list type="bullet">
///   <item>F73 sum-rule consistency: <c>‖Σ_b P_top^L · M_h_per_bond[b] · P_top^R‖_F² ≈
///   RankTop · σ_0²</c>. Equivalently: the projection of <c>M_H_total</c> on the top
///   eigenspace has Frobenius-square equal to the sum of σ_0² over the degenerate
///   directions.</item>
///   <item>Sign of magnitudes: every per-bond magnitude is real and non-negative
///   (Frobenius-square is by definition).</item>
/// </list>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (a''), and the
/// <see cref="C2HwhmRatio.PendingDerivationNote"/> ranking that promotes (d'')-then-(a'')
/// as the foundational path.</para>
/// </summary>
public sealed class C2SvdBlockProjectedMagnitude : Claim
{
    public CoherenceBlock Block { get; }

    /// <summary>The (d'') Riesz-lifted projectors used to build the magnitudes; held for
    /// downstream consumers that need the projectors directly without rebuilding.</summary>
    public C2InterChannelProjector Projector { get; }

    /// <summary>One <see cref="BondMagnitudeWitness"/> per bond, in bond-index order, tagged
    /// with <see cref="BondClass"/>.</summary>
    public IReadOnlyList<BondMagnitudeWitness> Witnesses { get; }

    /// <summary>Sum-rule residual <c>|‖Σ_b P_top^L · M_b · P_top^R‖_F² − RankTop · σ_0²|</c>.
    /// Below <see cref="SumRuleTolerance"/> the F73-rule projection is consistent at the
    /// tier guaranteed by the projector.</summary>
    public double SumRuleResidual { get; }

    public const double SumRuleTolerance = 1e-8;

    public static C2SvdBlockProjectedMagnitude Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2SvdBlockProjectedMagnitude applies only to the c=2 stratum; got c={block.C}.",
                nameof(block));

        var projector = C2InterChannelProjector.Build(block);
        var witnesses = BuildWitnesses(block, projector);
        double sumRuleResidual = ComputeSumRuleResidual(block, projector);

        return new C2SvdBlockProjectedMagnitude(block, projector, witnesses, sumRuleResidual);
    }

    private C2SvdBlockProjectedMagnitude(
        CoherenceBlock block,
        C2InterChannelProjector projector,
        IReadOnlyList<BondMagnitudeWitness> witnesses,
        double sumRuleResidual)
        : base("c=2 per-bond SVD-block projected magnitude (Direction (a'') foundation)",
               // Tier inherits the projector's Tier: Tier1Derived if (d'') passes its runtime
               // library-stability witness, Tier2Verified if it falls back.
               projector.Tier,
               Item1Anchors.Root)
    {
        Block = block;
        Projector = projector;
        Witnesses = witnesses;
        SumRuleResidual = sumRuleResidual;
    }

    private static IReadOnlyList<BondMagnitudeWitness> BuildWitnesses(
        CoherenceBlock block, C2InterChannelProjector projector)
    {
        int numBonds = block.NumBonds;
        var witnesses = new BondMagnitudeWitness[numBonds];
        double sigma0Squared = projector.Sigma0 * projector.Sigma0;
        for (int b = 0; b < numBonds; b++)
        {
            var mb = block.Decomposition.MhPerBond[b];
            var lifted = projector.PTopL * mb * projector.PTopR;
            double magnitudeSquared = lifted.FrobeniusNorm();
            magnitudeSquared *= magnitudeSquared;
            double normalisedRatio = sigma0Squared > 0 ? magnitudeSquared / sigma0Squared : 0.0;
            BondClass bondClass = ClassifyBond(b, numBonds);
            witnesses[b] = new BondMagnitudeWitness(
                Bond: b,
                BondClass: bondClass,
                MagnitudeSquared: magnitudeSquared,
                NormalisedRatio: normalisedRatio);
        }
        return witnesses;
    }

    private static double ComputeSumRuleResidual(
        CoherenceBlock block, C2InterChannelProjector projector)
    {
        // Σ_b P_top^L · M_b · P_top^R = P_top^L · (Σ_b M_b) · P_top^R = P_top^L · M_h_total · P_top^R.
        // The Frobenius-square of that should equal RankTop · σ_0² because P_top^{L,R} are
        // projections onto the σ_0-eigenspace of the SVD of M_h_total restricted to the HD-pair.
        var mTotalLifted = projector.PTopL * block.Decomposition.MhTotal * projector.PTopR;
        double frobSquared = mTotalLifted.FrobeniusNorm();
        frobSquared *= frobSquared;
        double expected = projector.RankTop * projector.Sigma0 * projector.Sigma0;
        return Math.Abs(frobSquared - expected);
    }

    private static BondClass ClassifyBond(int bond, int numBonds) =>
        bond == 0 || bond == numBonds - 1 ? BondClass.Endpoint : BondClass.Interior;

    /// <summary>Bond-class-mean magnitude over <see cref="Witnesses"/>: average of
    /// <see cref="BondMagnitudeWitness.NormalisedRatio"/> within the chosen
    /// <see cref="BondClass"/>. Returns NaN if the class has no bonds at this N (e.g.
    /// Interior at N = 4 has 1 bond; Endpoint always has 2 for N ≥ 3).</summary>
    public double NormalisedRatioMean(BondClass bondClass)
    {
        var hits = Witnesses.Where(w => w.BondClass == bondClass).Select(w => w.NormalisedRatio).ToArray();
        return hits.Length == 0 ? double.NaN : hits.Average();
    }

    public override string DisplayName =>
        $"c=2 per-bond SVD-block projected magnitude (N={Block.N})";

    public override string Summary
    {
        get
        {
            string endpoint, interior;
            try { endpoint = $"Endpoint mean={NormalisedRatioMean(BondClass.Endpoint):F4}"; }
            catch { endpoint = "Endpoint=n/a"; }
            try { interior = $"Interior mean={NormalisedRatioMean(BondClass.Interior):F4}"; }
            catch { interior = "Interior=n/a"; }
            return $"{endpoint}, {interior}, sum-rule residual={SumRuleResidual:G3} ({Tier.Label()})";
        }
    }

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("RankTop (from projector)", summary: Projector.RankTop.ToString());
            yield return InspectableNode.RealScalar("σ_0", Projector.Sigma0);
            yield return InspectableNode.RealScalar("SumRuleResidual", SumRuleResidual, "G3");
            yield return InspectableNode.RealScalar(
                "NormalisedRatioMean[Endpoint]", NormalisedRatioMean(BondClass.Endpoint), "F4");
            yield return InspectableNode.RealScalar(
                "NormalisedRatioMean[Interior]", NormalisedRatioMean(BondClass.Interior), "F4");
            yield return InspectableNode.Group("Witnesses (per bond)",
                Witnesses.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>Per-bond witness for <see cref="C2SvdBlockProjectedMagnitude"/>:
/// library-invariant SVD-block magnitude lifted via the Riesz projector.</summary>
/// <param name="Bond">Bond index in [0, NumBonds-1].</param>
/// <param name="BondClass">Endpoint (b ∈ {0, NumBonds−1}) or Interior.</param>
/// <param name="MagnitudeSquared">‖P_top^L · M_h_per_bond[b] · P_top^R‖_F² (real, non-negative).</param>
/// <param name="NormalisedRatio"><c>MagnitudeSquared / σ_0²</c>, the dimensionless ratio that
/// would appear in a Direction-(a'') closed-form ansatz <c>Δ(r)</c>.</param>
public sealed record BondMagnitudeWitness(
    int Bond,
    BondClass BondClass,
    double MagnitudeSquared,
    double NormalisedRatio
) : IInspectable
{
    public string DisplayName => $"bond {Bond} [{BondClass}]";

    public string Summary =>
        $"|V_b|²_proj = {MagnitudeSquared:G6}, |V_b|²_proj/σ_0² = {NormalisedRatio:F6}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond index", summary: Bond.ToString());
            yield return new InspectableNode("class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("MagnitudeSquared", MagnitudeSquared, "G6");
            yield return InspectableNode.RealScalar("NormalisedRatio", NormalisedRatio, "F6");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"|V_b|²_proj/σ_0² (bond {Bond}, {BondClass})", NormalisedRatio, "F6");
}
