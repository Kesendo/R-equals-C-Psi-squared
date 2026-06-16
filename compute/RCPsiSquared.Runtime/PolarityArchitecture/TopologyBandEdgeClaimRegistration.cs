using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="TopologyBandEdgeClaim"/>: typed parent edges to
/// <see cref="ClockHandLadderClaim"/> (the chain instance it generalizes) and
/// <see cref="AbsorptionTheoremClaim"/> (the Re=−2γ floor). Both resolved by build time.</summary>
public static class TopologyBandEdgeClaimRegistration
{
    public static ClaimRegistryBuilder RegisterTopologyBandEdgeClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<TopologyBandEdgeClaim>(b =>
        {
            var clock = b.Get<ClockHandLadderClaim>();
            var absorption = b.Get<AbsorptionTheoremClaim>();
            return new TopologyBandEdgeClaim(clock, absorption);
        });
}
