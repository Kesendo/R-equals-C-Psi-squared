using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Schicht-1 wiring of <see cref="GreyRoomLadderImagesClaim"/> (F171).
///
/// <para>Five typed parents, all resolved from the builder: <see cref="SpectatorIntertwinerClaim"/> supplies
/// the η ladder Φ with its exact intertwining, <see cref="PinnedBlockFloorClaim"/> the boundary blocks grey as
/// wholes and the window the landing rule is read in, <see cref="ChiralKClaim"/> the spectrum inversion behind
/// the (1, 1) count, <see cref="FrozenDivisorClaim"/> the involution U whose odd half is the distance-2 image in
/// (1, 1), <see cref="RingFrozenBandClaim"/> Lemma 1, the flux the particle number sets on the spin ring. The
/// builder resolves topologically, so the order of the Register calls in <c>KnowledgeRegistryFactory</c> is
/// free.</para></summary>
public static class GreyRoomLadderImagesClaimRegistration
{
    public static ClaimRegistryBuilder RegisterGreyRoomLadderImagesClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<GreyRoomLadderImagesClaim>(b =>
        {
            var spectator = b.Get<SpectatorIntertwinerClaim>();
            var pinned = b.Get<PinnedBlockFloorClaim>();
            var chiralK = b.Get<ChiralKClaim>();
            var frozenDivisor = b.Get<FrozenDivisorClaim>();
            var ringBand = b.Get<RingFrozenBandClaim>();
            return new GreyRoomLadderImagesClaim(spectator, pinned, chiralK, frozenDivisor, ringBand);
        });
}
