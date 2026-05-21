using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F86HwhmClosedFormClaim"/>: the F86 c=2
/// HWHM_left/Q_peak per-BondSubClass prediction (form HWHM_ratio = 0.671535 + α·g_eff + β;
/// bare floor analytically derived, the per-sub-class (α, β) fitted by polyfit on N=5..8
/// anchors).
///
/// <para>A standalone claim: parameterless, block-independent, no Claim parents (a registry
/// root). Tier1Candidate, not Tier1Derived, because the (α, β) are fitted rather than
/// derived from the F89/F90 bridge structure.</para></summary>
public static class F86HwhmClosedFormClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF86HwhmClosedFormClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F86HwhmClosedFormClaim>(_ => new F86HwhmClosedFormClaim());
}
