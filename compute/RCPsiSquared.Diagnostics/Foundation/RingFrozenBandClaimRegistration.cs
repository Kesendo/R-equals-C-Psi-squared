using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Schicht-1 wiring of <see cref="RingFrozenBandClaim"/> (F170).
///
/// <para>Two typed parents, both resolved from the builder: <see cref="FrozenDivisorClaim"/> supplies the
/// corner floor ⌊N/2⌋ (typed on the open chains; that its proof covers both fixed-flux rings is
/// <see cref="RingFrozenBandClaim"/>'s argument), and <see cref="SpectatorIntertwinerClaim"/> supplies
/// the η ladder whose reverse is the lowest-weight map. The builder resolves topologically, so the order of
/// the Register calls in <c>KnowledgeRegistryFactory</c> is free.</para></summary>
public static class RingFrozenBandClaimRegistration
{
    public static ClaimRegistryBuilder RegisterRingFrozenBandClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<RingFrozenBandClaim>(b =>
        {
            var frozenDivisor = b.Get<FrozenDivisorClaim>();
            var spectator = b.Get<SpectatorIntertwinerClaim>();
            return new RingFrozenBandClaim(frozenDivisor, spectator);
        });
}
