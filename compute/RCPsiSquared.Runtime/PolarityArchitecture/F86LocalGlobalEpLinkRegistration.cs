using RCPsiSquared.Core.F86;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring for <see cref="LocalGlobalEpLink"/> (Locus 5 — F86 ↔
/// FRAGILE_BRIDGE axis-departure comparison; FRAGILE_BRIDGE EP character remains OPEN.
/// The secured shared substrate is the shifted-generator sectorwise P-type anticommutation,
/// while the full SRP class remains open).
///
/// <para>Parallel to <see cref="F86PolarityLinkRegistration"/> (Locus 6, polarity-side
/// closure). Together they bracket the F86 c=2 derivation with comparison-side and symmetry-side
/// parent-claim references — exactly the inheritance graph reading in
/// <c>project_algebra_is_inheritance.md</c> Locus 5.</para>
///
/// <para>Edge declared: <see cref="LocalGlobalEpLink"/> ← <see cref="ShiftedGeneratorSectorwisePClaim"/>.
/// The parent is the Tier1Derived operator relation; it is not a class label.</para>
///
/// <para>Tier consistency: LocalGlobalEpLink is <see cref="Knowledge.Tier.OpenQuestion"/>
/// (demoted from Tier2Verified by the F86a-retraction 2026-06-21); ShiftedGeneratorSectorwisePClaim
/// is Tier1Derived. The TierStrength inheritance check (parent at least as strong as child)
/// passes (5 ≥ 1), so the parent edge survives the demotion.</para>
///
/// <para>The retraction (F86a, 2026-06-21): the full Σγ=N·γ₀ block is genuinely non-normal
/// on the real Q axis, Petermann large but finite. (The retraction's "no real-axis EP,
/// eigenvalues simple" was itself corrected 2026-07-07: F89 locates a real-axis defective seed
/// on this block, see <see cref="LocalGlobalEpLink"/> Summary / PROOF_F86A section The real-axis EP.)
/// The DISTINCT off-axis complex-Q question stays open, documented in
/// <see cref="LocalGlobalEpLink.PendingDerivationNote"/>. The Schicht-1 wiring keeps the c=2
/// N=5..8 PetermannSpikeWitness cautionary non-normality table queryable through the registry;
/// future drift checks against extended-N or extended-c sweeps land here.</para></summary>
public static class F86LocalGlobalEpLinkRegistration
{
    public static ClaimRegistryBuilder RegisterF86LocalGlobalEpLink(this ClaimRegistryBuilder builder) =>
        builder.Register<LocalGlobalEpLink>(b =>
        {
            _ = b.Get<ShiftedGeneratorSectorwisePClaim>();
            return LocalGlobalEpLink.Build();
        });
}
