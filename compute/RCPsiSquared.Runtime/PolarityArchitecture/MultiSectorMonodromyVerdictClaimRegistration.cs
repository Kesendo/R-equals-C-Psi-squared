using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="MultiSectorMonodromyVerdictClaim"/>: the N-dependent multi-sector
/// monodromy verdict (the S₈ braid is confined to the (1,2) orbit at N=4, spreads to a 12-sector joint-popcount
/// diamond at N=5). Three typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89OcticMonodromyClaim"/>: the S₈ braid the (1,2) octic carries, which this census tests for
///         localization vs spread across the joint-popcount sectors.</item>
///   <item><see cref="F89CrossFoldSimilarityClaim"/>: F89d, the exact antiunitary similarity λ ↦ −λ̄ − 2N that
///         pairs Family B with Family A, and whose N=4 self-fold (the cross-fold as a self-map on the |Δ|=1 orbit)
///         IS the confinement.</item>
///   <item><see cref="SpectatorIntertwinerClaim"/>: Theorem B of PROOF_CODIM1_BY_ADDITIVITY, the exact W
///         intertwiner that transports the EP up the diamond (the byte-identity mechanism) and derives the
///         CONTAINMENT half of general-N membership (added 2026-07-02).</item>
/// </list>
///
/// <para>All parents Tier 1 derived and registered above; the builder topo-resolves, so registration order is
/// free. Live: <c>inspect --root sectorbraid</c>.</para></summary>
public static class MultiSectorMonodromyVerdictClaimRegistration
{
    public static ClaimRegistryBuilder RegisterMultiSectorMonodromyVerdictClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<MultiSectorMonodromyVerdictClaim>(b =>
            new MultiSectorMonodromyVerdictClaim(
                b.Get<F89OcticMonodromyClaim>(),
                b.Get<F89CrossFoldSimilarityClaim>(),
                b.Get<SpectatorIntertwinerClaim>()));
}
