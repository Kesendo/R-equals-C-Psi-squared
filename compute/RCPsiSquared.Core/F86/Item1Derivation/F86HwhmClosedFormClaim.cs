using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 c=2 HWHM_left/Q_peak prediction per BondSubClass.
///
/// <para><b>Tier 2 empirical</b>: the form is `floor + alpha_subclass · g_eff + beta_subclass`
/// with the floor the bare doubled-PTF envelope ratio 0.7264892771 (derived,
/// <see cref="C2BareDoubledPtfClosedForm.EnvelopeHwhmLeftOverXPeak"/>) and the (alpha, beta) per
/// sub-class read live from <see cref="F86HwhmAlphaExtraction"/>, a line through the N=5..8
/// envelope anchors. No analytical reduction derives them from F89 cyclotomic Phi_{N+1} / F90
/// bridge identity structure.</para>
///
/// <para><b>What the anchors carry:</b> the canonical peaks (first local maximum) on the default
/// Q grid, within 3.5·10⁻⁵ of the continuous envelope's ratio. Endpoint and Flanking have four
/// chain lengths each and miss a line by up to 1.28·10⁻³ and 4.2·10⁻⁴, 12 to 36 times that
/// offset, so the linear-in-g_eff form is not an exact law there; Mid and CentralSelfPaired have two points (the
/// line is exact and tests nothing). Orbit2Escape and CentralEscapeOrbit3 classify a far maximum
/// and have no canonical anchors; <see cref="PredictHwhmRatio"/> throws for them.</para>
///
/// <para><b>What IS derived (Tier 1 sub-results inside this empirical claim):</b>
/// <list type="bullet">
///   <item>the envelope floor 0.7264892771, from <see cref="C2BareDoubledPtfClosedForm"/></item>
///   <item>Sub-class classification (<see cref="BondSubClass"/>): structural via
///         F71 symmetry + Q_peak escape thresholds</item>
///   <item>g_eff = 2 · BareDoubledPtfXPeak / (Q_peak + 2), with the fixed-τ = 1 bare peak
///         2.196910 (the predictor's definition; the bare envelope peaks at 1.7141)</item>
/// </list></para>
///
/// <para><b>What is fitted (open analytical step):</b> the (alpha, beta) per sub-class. The
/// PROOF_F90_F86C2_BRIDGE (Item 1', Direction (b'')) locates the lift above the bare floor in
/// the intra-channel dispersion of the rank-1-bridge + intra-dispersion model (two-dial
/// scout, 2026-06-11; the inter tail only moves Q_peak/g_eff), with the cyclotomic Phi_{N+1}
/// pattern supplying the N-scaling. Deriving the lift from that structure is the open step;
/// the line through the anchors is not that step.</para></summary>
public sealed class F86HwhmClosedFormClaim : Claim
{
    private const double BareFloorValue = F86HwhmAlphaExtraction.BareFloor;

    /// <summary>The g_eff conversion factor: 4.39382 = 2 × C2HwhmRatio.BareDoubledPtfXPeak (2.196910),
    /// the fixed-τ = 1 bare peak. Anchored to PROOF_F86_QPEAK.md.</summary>
    private const double BareDoubledPtfXPeakValue = F86HwhmAlphaExtraction.BareDoubledPtfXPeak;

    public double BareFloor => BareFloorValue;
    public double BareDoubledPtfXPeak => BareDoubledPtfXPeakValue;

    // Per-sub-class (alpha, beta) such that HWHM_ratio = BareFloor + alpha * g_eff + beta, read
    // LIVE from F86HwhmAlphaExtraction (no hand-transcription, so they cannot drift from their
    // source).
    private static readonly Dictionary<BondSubClass, (double Alpha, double Beta)> _subClassParams =
        F86HwhmAlphaExtraction.Extract().Keys.ToDictionary(s => s, F86HwhmAlphaExtraction.FitParams);

    public F86HwhmClosedFormClaim()
        : base("F86 c=2 HWHM_left/Q_peak per BondSubClass (Tier 2 empirical; floor derived, (alpha, beta) a line through the anchors, not an exact law, analytical lift open)",
               Tier.Tier2Empirical,
               "compute/RCPsiSquared.Core/F86/Item1Derivation/BondSubClass.cs + " +
               "compute/RCPsiSquared.Core/F86/Item1Derivation/C2BareDoubledPtfClosedForm.cs (EnvelopeHwhmLeftOverXPeak = 0.7264892771) + " +
               "compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs (PathPolynomial path-3..9) + " +
               "compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs (numerical anchor) + " +
               "docs/proofs/PROOF_F90_F86C2_BRIDGE.md (Item 1', partial closure) + " +
               "docs/superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md + " +
               "simulations/f86_envelope_second_route.py (second route for the anchors) + " +
               "compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmAlphaExtraction.cs (the live line per sub-class)")
    { }

    public double PredictHwhmRatio(int n, int bondIndex, double qPeak)
    {
        var sub = BondSubClassExtensions.Classify(n, bondIndex, qPeak);
        var (alpha, beta) = F86HwhmAlphaExtraction.FitParams(sub);
        double gEff = BareDoubledPtfXPeakValue / (qPeak + 2.0);
        return BareFloorValue + alpha * gEff + beta;
    }

    public override string DisplayName =>
        "F86 c=2 HWHM per BondSubClass (Tier 2 empirical; (alpha, beta) a line through the anchors, analytical lift open)";

    public override string Summary =>
        $"HWHM_ratio = {BareFloorValue} + alpha_subclass * g_eff + beta_subclass, read live from " +
        $"F86HwhmAlphaExtraction on the canonical peaks; {_subClassParams.Count} sub-classes with anchors: Endpoint and " +
        $"Flanking (four chain lengths each) miss a line by up to 1.3e-3 and 4.2e-4, Mid and CentralSelfPaired are " +
        $"two-point lines; analytical derivation of the lift from F89/F90 still open ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("BareFloor (bare doubled-PTF envelope ratio)", BareFloorValue);
            yield return InspectableNode.RealScalar("BareDoubledPtfXPeak (g_eff conversion factor)", BareDoubledPtfXPeakValue);
            // The live witness: each sub-class's line recomputed from the anchors, with how many
            // distinct g_eff values it rests on and how far the anchors sit from it.
            foreach (var ex in F86HwhmAlphaExtraction.Extract().Values)
            {
                yield return new InspectableNode($"{ex.SubClass} [{ex.Verdict}]",
                    summary: $"alpha={ex.FittedAlpha:G6}, beta={ex.FittedBeta:G6}; {ex.DistinctG} distinct g_eff over " +
                             $"{ex.N} bonds{(ex.NEscape > 0 ? $" ({ex.NEscape} past Q = 4)" : "")}, g-lever {ex.GSpan:G3}, " +
                             $"max residual {ex.MaxResidual:G3}",
                    provenance: NodeProvenance.Live);
            }
        }
    }
}
