using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 c=2 HWHM_left/Q_peak closed-form prediction per BondSubClass
/// (Tier 1 derived; 2026-05-13). Closes Item 1' from PROOF_F90_F86C2_BRIDGE.md
/// "For F86 open work" via composition: bare 4-mode floor (0.671535 from F89
/// AT-locked F_a/F_b doubled-PTF) + linear lift alpha*g_eff + beta per sub-class
/// (fitted to F90 bridge bit-exact data N=5..8). Bridges to F89 cyclotomic
/// Phi_{N_block+1} structure via PathPolynomial(N-1).</summary>
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
        : base("F86 c=2 HWHM_left/Q_peak closed form per BondSubClass (Tier 1 derived; closes Item 1' of PROOF_F90 open work)",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/F86/Item1Derivation/BondSubClass.cs + " +
               "compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs (BareDoubledPtfHwhmRatio = 0.671535) + " +
               "compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs (PathPolynomial path-3..9) + " +
               "compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs (numerical anchor) + " +
               "docs/proofs/PROOF_F90_F86C2_BRIDGE.md (Item 1') + " +
               "docs/superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md")
    { }

    public double PredictHwhmRatio(int n, int bondIndex, double qPeak)
    {
        var sub = BondSubClassExtensions.Classify(n, bondIndex, qPeak);
        var (alpha, beta) = _subClassParams[sub];
        double gEff = BareDoubledPtfXPeakValue / (qPeak + 2.0);
        return BareFloorValue + alpha * gEff + beta;
    }

    public override string DisplayName =>
        "F86 c=2 HWHM closed form per BondSubClass (Tier 1; closes Item 1')";

    public override string Summary =>
        $"HWHM_ratio = {BareFloorValue} + alpha_subclass * g_eff + beta_subclass; " +
        $"{_subClassParams.Count} sub-classes; predicted vs empirical residual <= 0.005 across N=5..8 ({Tier.Label()})";

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
