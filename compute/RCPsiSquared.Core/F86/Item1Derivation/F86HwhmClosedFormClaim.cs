using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 c=2 HWHM_left/Q_peak prediction per BondSubClass.
///
/// <para><b>Tier 1 candidate</b>: the form is `0.671535 + alpha_subclass · g_eff +
/// beta_subclass` with the bare floor 0.671535 analytically derived (from F89 AT-locked
/// F_a/F_b doubled-PTF), but the (alpha, beta) per sub-class are <b>fitted</b> by linear
/// regression on N=5..8 F90 bridge anchors (`simulations/f86_hwhm_closed_form_verification.py`,
/// line 78 `np.polyfit(...deg=1)`). The fit reproduces the 22 anchors within 0.005 residual
/// (verification confirms the form is approximately right), but no analytical reduction
/// derives the (alpha, beta) values from F89 cyclotomic Phi_{N+1} / F90 bridge identity
/// structure.</para>
///
/// <para><b>Honest reshape (f86b2_robust_extraction arc, 2026-06-29):</b> the (alpha, beta) are
/// read live from <see cref="F86HwhmAlphaExtraction"/>, and only Endpoint and Flanking carry a
/// resolved slope (jackknife-tight, though marginal at the 0.005 grid floor). Mid is a flat lift
/// over a microscopic g_eff lever (slope = noise), two classes have a single anchor, and
/// Orbit2Escape sits on grid-edge anchors flagged non-physical; these four collapse to a
/// per-class constant lift (alpha = 0) or a flagged fit, all still reproducing their anchors
/// within 0.005.</para>
///
/// <para><b>What IS derived (Tier 1 sub-results inside this candidate class):</b>
/// <list type="bullet">
///   <item>BareFloor = 0.671535 — from <see cref="C2BareDoubledPtfClosedForm"/></item>
///   <item>Sub-class classification (<see cref="BondSubClass"/>) — structural via
///         F71 symmetry + Q_peak escape thresholds</item>
///   <item>g_eff = 2 · BareDoubledPtfXPeak / (Q_peak + 2) — bare-doubled-PTF coupling</item>
/// </list></para>
///
/// <para><b>What is fitted (open analytical step):</b> the 12 (alpha, beta) values per
/// sub-class. To promote this to Tier 1 derived, these values must derive from F89
/// cyclotomic Phi_{N+1} structure via F90 bridge identity, NOT from polyfit. The
/// PROOF_F90_F86C2_BRIDGE (Item 1', Direction (b'')) locates the lift above the 4-mode
/// floor (0.6715) to 0.7506/0.7728 in the intra-channel dispersion of the rank-1-bridge
/// + intra-dispersion model (two-dial scout, 2026-06-11; the earlier "octic residual of
/// the inter coupling" lift suspicion is refuted, the inter tail only moves Q_peak/g_eff),
/// with the cyclotomic Phi_{N+1} pattern supplying the N-scaling. Deriving (alpha, beta)
/// from that structure is the open analytical step that would close the Tier 1 derivation;
/// the fit IS NOT that step.</para></summary>
public sealed class F86HwhmClosedFormClaim : Claim
{
    private const double BareFloorValue = 0.671535;

    /// <summary>The g_eff conversion factor: 4.39382 = 2 × C2HwhmRatio.BareDoubledPtfXPeak (2.196910).
    /// Matches the Python BARE_DOUBLED_PTF_XPEAK in
    /// simulations/f86_hwhm_subclass_stratification.py used to generate the (alpha, beta) per-sub-class
    /// fit in <see cref="_subClassParams"/>. Anchored to PROOF_F86_QPEAK.md.</summary>
    private const double BareDoubledPtfXPeakValue = 4.39382;

    public double BareFloor => BareFloorValue;
    public double BareDoubledPtfXPeak => BareDoubledPtfXPeakValue;

    // Per-sub-class (alpha, beta) such that HWHM_ratio = BareFloor + alpha * g_eff + beta, read
    // LIVE from F86HwhmAlphaExtraction (no hand-transcription, so they cannot drift from their
    // source). The honest reshape (f86b2_robust_extraction arc, 2026-06-29): only Endpoint and
    // Flanking keep a resolved slope; Mid (noise), the two single-point classes, and the
    // grid-edge Orbit2Escape collapse to a per-class constant lift / flagged fit.
    private static readonly Dictionary<BondSubClass, (double Alpha, double Beta)> _subClassParams =
        Enum.GetValues<BondSubClass>().ToDictionary(s => s, F86HwhmAlphaExtraction.ReshapedParams);

    public F86HwhmClosedFormClaim()
        : base("F86 c=2 HWHM_left/Q_peak per BondSubClass (Tier 1 candidate; form derived, (alpha, beta) fitted via polyfit, analytical lift open)",
               Tier.Tier1Candidate,
               "compute/RCPsiSquared.Core/F86/Item1Derivation/BondSubClass.cs + " +
               "compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs (BareDoubledPtfHwhmRatio = 0.671535) + " +
               "compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs (PathPolynomial path-3..9) + " +
               "compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs (numerical anchor) + " +
               "docs/proofs/PROOF_F90_F86C2_BRIDGE.md (Item 1', partial closure) + " +
               "docs/superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md + " +
               "simulations/f86_hwhm_closed_form_verification.py (polyfit source for (alpha, beta)) + " +
               "compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmAlphaExtraction.cs (the live recompute + honest reshape)")
    { }

    public double PredictHwhmRatio(int n, int bondIndex, double qPeak)
    {
        var sub = BondSubClassExtensions.Classify(n, bondIndex, qPeak);
        var (alpha, beta) = _subClassParams[sub];
        double gEff = BareDoubledPtfXPeakValue / (qPeak + 2.0);
        return BareFloorValue + alpha * gEff + beta;
    }

    public override string DisplayName =>
        "F86 c=2 HWHM per BondSubClass (Tier 1 candidate; (alpha, beta) fitted, analytical lift open)";

    public override string Summary =>
        $"HWHM_ratio = {BareFloorValue} + alpha_subclass * g_eff + beta_subclass, read live from " +
        $"F86HwhmAlphaExtraction; {_subClassParams.Count} sub-classes, but only Endpoint/Flanking carry a " +
        $"resolved slope (jackknife-tight, marginal at the 0.005 grid floor) -- Mid (noise), two single-point " +
        $"classes, and the grid-edge Orbit2Escape are not defensible slopes (f86b2_robust_extraction, 2026-06-29); " +
        $"reproduces 22 anchors within 0.005; analytical (alpha, beta) derivation from F89/F90 still open ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("BareFloor (4-mode F_a/F_b doubled-PTF)", BareFloorValue);
            yield return InspectableNode.RealScalar("BareDoubledPtfXPeak (g_eff conversion factor)", BareDoubledPtfXPeakValue);
            // The live witness: each sub-class's slope recomputed from the anchors with its
            // grid-noise sigma and the honest verdict (resolved / noise / single-point / escape).
            foreach (var ex in F86HwhmAlphaExtraction.Extract().Values)
            {
                var (alpha, beta) = F86HwhmAlphaExtraction.ReshapedParams(ex.SubClass);
                yield return new InspectableNode($"{ex.SubClass} [{ex.Verdict}]",
                    summary: $"used (alpha={alpha:G6}, beta={beta:G6}); fitted slope {ex.FittedAlpha:G4} " +
                             $"+- {ex.SigmaAlpha:G3} (grid sigma_y=0.005), n={ex.N}" +
                             $"{(ex.NEscape > 0 ? $", {ex.NEscape} grid-edge" : "")}, g-lever {ex.GSpan:G3}",
                    provenance: NodeProvenance.Live);
            }
        }
    }
}
