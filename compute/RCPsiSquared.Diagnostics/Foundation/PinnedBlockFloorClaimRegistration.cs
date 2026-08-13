using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Schicht-1 wiring of <see cref="PinnedBlockFloorClaim"/> (F153, the pinning criterion:
/// 4N of the (N+1)² joint-popcount blocks sit entirely on their Absorption floor). Two typed parents,
/// and the claim is not derivable from either alone: <see cref="AbsorptionTheoremClaim"/> supplies the
/// floor and the per-site law whose collapse to a COUNT is what pinning buys, and
/// <see cref="JointPopcountSectors"/> supplies the (N+1)² grading the count 4N is taken in. Must register
/// after both <c>RegisterAbsorptionTheoremClaim</c> and <c>RegisterJointPopcountSectors</c>.</summary>
public static class PinnedBlockFloorClaimRegistration
{
    public static ClaimRegistryBuilder RegisterPinnedBlockFloorClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<PinnedBlockFloorClaim>(b =>
        {
            var absorption = b.Get<AbsorptionTheoremClaim>();   // typed parent edge: the floor
            var sectors = b.Get<JointPopcountSectors>();        // typed parent edge: the grading
            return new PinnedBlockFloorClaim(absorption, sectors);
        });
}
