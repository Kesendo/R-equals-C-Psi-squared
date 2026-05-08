using RCPsiSquared.Core.F86;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring for <see cref="LocalGlobalEpLink"/> (Locus 5 — F86 ↔
/// FRAGILE_BRIDGE shared exceptional-point structure under AIII chiral algebra).
///
/// <para>Parallel to <see cref="F86PolarityLinkRegistration"/> (Locus 6, polarity-side
/// closure). Together they bracket the F86 c=2 derivation with EP-side and symmetry-side
/// parent-claim references — exactly the inheritance graph reading in
/// <c>project_algebra_is_inheritance.md</c> Locus 5.</para>
///
/// <para>Edge declared: <see cref="LocalGlobalEpLink"/> ← <see cref="ChiralAiiiClassification"/>.
/// Both are F86 main-family claims; ChiralAiiiClassification is the Tier1Derived algebraic
/// label that F86 and FRAGILE_BRIDGE share, so the EP-link inherits from it.</para>
///
/// <para>Tier consistency: LocalGlobalEpLink is Tier2Verified; ChiralAiiiClassification is
/// Tier1Derived. The TierStrength inheritance check (parent at least as strong as child)
/// passes (5 ≥ 3).</para>
///
/// <para>The complex-γ analytic continuation that would promote LocalGlobalEpLink to
/// Tier1Derived is documented as the explicit gap in
/// <see cref="LocalGlobalEpLink.PendingDerivationNote"/>. Until then, the Schicht-1 wiring
/// makes the c=2 N=5..8 PetermannSpikeWitness table queryable through the registry; future
/// drift checks against extended-N or extended-c sweeps land here.</para></summary>
public static class F86LocalGlobalEpLinkRegistration
{
    public static ClaimRegistryBuilder RegisterF86LocalGlobalEpLink(this ClaimRegistryBuilder builder) =>
        builder.Register<LocalGlobalEpLink>(b =>
        {
            _ = b.Get<ChiralAiiiClassification>();
            return LocalGlobalEpLink.Build();
        });
}
