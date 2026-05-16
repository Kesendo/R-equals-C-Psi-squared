using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 c=2 HWHM_left/Q_peak prediction per BondSubClass.
///
/// <para><b>Tier 1 candidate</b>: the form is `0.671535 + alpha_subclass · g_eff +
/// beta_subclass` with the bare floor 0.671535 analytically derived (from F89 AT-locked
/// F_a/F_b doubled-PTF), but the (alpha, beta) per sub-class are <b>fitted</b> by linear
/// regression on N=5..8 F90 bridge anchors (`simulations/_f86_hwhm_closed_form_verification.py`,
/// line 78 `np.polyfit(...deg=1)`). The fit reproduces the 22 anchors within 0.005 residual
/// (verification confirms the form is approximately right), but no analytical reduction
/// derives the (alpha, beta) values from F89 cyclotomic Phi_{N+1} / F90 bridge identity
/// structure.</para>
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
/// PROOF_F90_F86C2_BRIDGE notes (Item 1'-followup) that "Closed-form via F89 AT-locked
/// F_a/F_b (4-mode floor 0.6715) + H_B-mixed octic residual (lift to 0.7506/0.7728) is
/// next analytical step" — this open analytical step is what would close the Tier 1
/// derivation; the fit IS NOT that step.</para></summary>
public sealed class F86HwhmClosedFormClaim : Claim
{
    private const double BareFloorValue = 0.671535;

    /// <summary>The g_eff conversion factor: 4.39382 = 2 × C2HwhmRatio.BareDoubledPtfXPeak (2.196910).
    /// Matches the Python BARE_DOUBLED_PTF_XPEAK in
    /// simulations/_f86_hwhm_subclass_stratification.py used to generate the (alpha, beta) per-sub-class
    /// fit in <see cref="_subClassParams"/>. Anchored to PROOF_F86_QPEAK.md.</summary>
    private const double BareDoubledPtfXPeakValue = 4.39382;

    public double BareFloor => BareFloorValue;
    public double BareDoubledPtfXPeak => BareDoubledPtfXPeakValue;

    // Per-sub-class linear-fit parameters (alpha, beta) such that
    // HWHM_ratio = BareFloor + alpha * g_eff + beta.
    // Values from Phase 4 verification self-fit on F90 bridge anchors N=5..8
    // (simulations/_f86_hwhm_closed_form_verification.py output, 2026-05-13).
    private static readonly Dictionary<BondSubClass, (double Alpha, double Beta)> _subClassParams = new()
    {
        [BondSubClass.Endpoint]            = (-0.129110,  0.227413),
        [BondSubClass.Flanking]            = (-0.094978,  0.193098),
        [BondSubClass.Mid]                 = ( 0.056559,  0.005165),
        [BondSubClass.CentralSelfPaired]   = ( 0.057439,  0.000000),
        [BondSubClass.Orbit2Escape]        = ( 0.698446, -0.086386),
        [BondSubClass.CentralEscapeOrbit3] = (-0.400854,  0.000000),
    };

    public F86HwhmClosedFormClaim()
        : base("F86 c=2 HWHM_left/Q_peak per BondSubClass (Tier 1 candidate; form derived, (alpha, beta) fitted via polyfit, analytical lift open)",
               Tier.Tier1Candidate,
               "compute/RCPsiSquared.Core/F86/Item1Derivation/BondSubClass.cs + " +
               "compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs (BareDoubledPtfHwhmRatio = 0.671535) + " +
               "compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs (PathPolynomial path-3..9) + " +
               "compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs (numerical anchor) + " +
               "docs/proofs/PROOF_F90_F86C2_BRIDGE.md (Item 1', partial closure) + " +
               "docs/superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md + " +
               "simulations/_f86_hwhm_closed_form_verification.py (polyfit source for (alpha, beta))")
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
        $"HWHM_ratio = {BareFloorValue} + alpha_subclass * g_eff + beta_subclass; " +
        $"{_subClassParams.Count} sub-classes; (alpha, beta) fitted via polyfit on N=5..8 anchors; " +
        $"reproduces 22 anchors within 0.005 residual; analytical derivation of (alpha, beta) from F89/F90 structure open ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("BareFloor (4-mode F_a/F_b doubled-PTF)", BareFloorValue);
            yield return InspectableNode.RealScalar("BareDoubledPtfXPeak (g_eff conversion factor)", BareDoubledPtfXPeakValue);
            foreach (var (sub, (alpha, beta)) in _subClassParams)
                yield return new InspectableNode($"{sub} fit",
                    summary: $"alpha = {alpha:G6}, beta = {beta:G6}");
        }
    }
}
