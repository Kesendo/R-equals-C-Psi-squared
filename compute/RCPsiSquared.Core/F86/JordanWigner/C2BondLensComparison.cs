using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T7 step C: per-bond cross-lens comparison
/// between the EP-resonance lens (<see cref="C2BondLQPeakScan"/>: Q_peak/HWHM of
/// ‖xB(Q)‖_F) and the K-resonance lens (<see cref="C2HwhmRatio"/>: Q_peak/HWHM of the
/// time-evolved K_b(Q, t_peak) from <see cref="ResonanceScan"/>).
///
/// <para>The two lenses measure different aspects of the bond's dynamics:</para>
/// <list type="bullet">
///   <item><b>EP-lens</b>: detects the L(Q)-eigenbasis singularity at Q ≈ Q_EP_b, where
///   R(Q) becomes degenerate. Sharp peak (HWHM small relative to Q_peak), bond-position-
///   distinct Q_peak (= bond's effective Q_EP).</item>
///   <item><b>K-lens</b>: detects the time-averaged K-resonance K_b(Q, t_peak)
///   ≈ 2·Re⟨ρ(t_peak) | S·∂_J ρ(t_peak)⟩, smoothed by the Frechet-derivative integral over
///   the eigenmode lifetimes. Broad peak (HWHM ≈ 0.75·Q_peak), bond-relatively-uniform
///   Q_peak (≈ BareDoubledPtfXPeak ≈ 2.197) but bond-distinct K_max and HWHM.</item>
/// </list>
///
/// <para>Empirical observation: the two Q_peaks differ. EP-Q_peak = Q_EP_b ≈ 1/g_eff_b
/// (bond-position-distinct via OBC sine-mode amplitude); K-Q_peak ≈ 2.197 · Q_EP_b
/// (BareDoubledPtfXPeak universality). The K-Q_peak is MORE bond-uniform than the
/// EP-Q_peak because the Frechet-derivative integral smooths bond-position dependence.
/// HWHM_left/Q_peak gap of factor ~10 between lenses (sharp vs broad).</para>
///
/// <para><b>Class-level Tier: Tier2Verified.</b> Numerical cross-lens comparison; the
/// closed-form relation between the two lenses is open. Tier1 promotion path: derive the
/// K-resonance closed-form via Frechet integral over the EP-projected MhPerBond[b]
/// matrix elements, yielding K-Q_peak / K-HWHM as analytical functions of EP-Q_peak +
/// dispersion ε_k.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c>.</para>
/// </summary>
public sealed class C2BondLensComparison : Claim
{
    public CoherenceBlock Block { get; }
    public C2BondLQPeakScan EpLens { get; }
    public C2HwhmRatio KLens { get; }
    public IReadOnlyList<BondLensComparisonWitness> Bonds { get; }

    public static C2BondLensComparison Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BondLensComparison applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var ep = C2BondLQPeakScan.Build(block);
        var k = C2HwhmRatio.Build(block);

        if (ep.Bonds.Count != k.Witnesses.Count)
            throw new InvalidOperationException(
                $"Lens mismatch: EP-lens has {ep.Bonds.Count} bonds, K-lens has {k.Witnesses.Count}");

        var witnesses = new BondLensComparisonWitness[ep.Bonds.Count];
        for (int b = 0; b < ep.Bonds.Count; b++)
        {
            var epBond = ep.Bonds[b];
            var kBond = k.Witnesses[b];
            witnesses[b] = new BondLensComparisonWitness(
                Bond: b,
                BondClass: epBond.BondClass,
                EpQPeak: epBond.QPeak,
                EpKMax: epBond.XbNormAtPeak,
                EpHwhmLeftOverQPeak: epBond.HwhmLeftOverQPeak,
                KQPeak: kBond.QPeak,
                KKMax: kBond.KMax,
                KHwhmLeftOverQPeak: kBond.HwhmLeft is double hl ? hl / kBond.QPeak : (double?)null);
        }

        return new C2BondLensComparison(block, ep, k, witnesses);
    }

    private C2BondLensComparison(
        CoherenceBlock block,
        C2BondLQPeakScan ep,
        C2HwhmRatio k,
        IReadOnlyList<BondLensComparisonWitness> bonds)
        : base("c=2 per-bond cross-lens comparison (EP-resonance vs K-resonance)",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        Block = block;
        EpLens = ep;
        KLens = k;
        Bonds = bonds;
    }

    public override string DisplayName =>
        $"c=2 BondLensComparison (N={Block.N}, {Bonds.Count} bonds, EP vs K)";

    public override string Summary =>
        $"per-bond Q_peak / HWHM cross-lens diagnostic at c=2 N={Block.N} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Bonds.Count.ToString());
            yield return EpLens;
            yield return KLens;
            yield return InspectableNode.Group("Bonds", Bonds.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One bond's cross-lens witness: pairs EP-lens (sharp ‖xB(Q)‖_F resonance) with
/// K-lens (broad K_b(Q, t_peak) resonance) for direct comparison. Both lenses are
/// F71-mirror-invariant individually; this witness exposes them side-by-side.</summary>
public sealed record BondLensComparisonWitness(
    int Bond,
    BondClass BondClass,
    double EpQPeak,
    double EpKMax,
    double? EpHwhmLeftOverQPeak,
    double KQPeak,
    double KKMax,
    double? KHwhmLeftOverQPeak
) : IInspectable
{
    public double QPeakRatio_KOverEp => KQPeak / EpQPeak;
    public double? HwhmRatioGap =>
        EpHwhmLeftOverQPeak is double e && KHwhmLeftOverQPeak is double k ? k - e : null;

    public string DisplayName => $"bond {Bond} cross-lens ({BondClass})";

    public string Summary =>
        $"EP Q_peak = {EpQPeak:F4} ‖xB‖_max = {EpKMax:F4}; " +
        $"K Q_peak = {KQPeak:F4} K_max = {KKMax:F4}; ratio K/EP = {QPeakRatio_KOverEp:F4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("EP Q_peak", EpQPeak, "F4");
            yield return InspectableNode.RealScalar("EP K_max (‖xB‖_max)", EpKMax, "F4");
            if (EpHwhmLeftOverQPeak is not null)
                yield return InspectableNode.RealScalar("EP HWHM_left/Q_peak", EpHwhmLeftOverQPeak.Value, "F4");
            yield return InspectableNode.RealScalar("K Q_peak", KQPeak, "F4");
            yield return InspectableNode.RealScalar("K K_max", KKMax, "F4");
            if (KHwhmLeftOverQPeak is not null)
                yield return InspectableNode.RealScalar("K HWHM_left/Q_peak", KHwhmLeftOverQPeak.Value, "F4");
            yield return InspectableNode.RealScalar("Q_peak ratio (K / EP)", QPeakRatio_KOverEp, "F4");
            if (HwhmRatioGap is not null)
                yield return InspectableNode.RealScalar("HWHM_left/Q_peak gap (K − EP)", HwhmRatioGap.Value, "F4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
