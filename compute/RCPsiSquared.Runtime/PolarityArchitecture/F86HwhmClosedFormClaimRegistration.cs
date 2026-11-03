using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F86HwhmClosedFormClaim"/>: the F86 c=2
/// HWHM_left/Q_peak per-BondSubClass prediction (form HWHM_ratio = floor + α·g_eff + β; the
/// floor, the bare doubled-PTF envelope ratio 0.7264892771, analytically derived, the
/// per-sub-class (α, β) a line through the N=5..8 envelope anchors).
///
/// <para>A standalone claim: parameterless, block-independent, no Claim parents (a registry
/// root). Tier2Empirical: the (α, β) are fitted rather than derived from the F89/F90 bridge
/// structure, and the four-point sub-classes miss the line, so no exact law holds the lift.</para></summary>
public static class F86HwhmClosedFormClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF86HwhmClosedFormClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F86HwhmClosedFormClaim>(_ => new F86HwhmClosedFormClaim());
}
