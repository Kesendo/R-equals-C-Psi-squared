using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="HandoverFloorClaim"/> (the handover Q = the F50-floor condition;
/// chain = Q*(N), ring = a distinct (2,2) level crossing; Tier1Candidate). Typed parents (all
/// registered above): <see cref="AbsorptionTheoremClaim"/> (the -2g&lt;n_XY&gt; survivor rate),
/// <see cref="F50WeightOneDegeneracyPi2Inheritance"/> (the off-diagonal floor &lt;n_XY&gt;=1), and
/// <see cref="CoherenceHorizonClaim"/> (the chain solution Q*(N)). Must register AFTER all three.</summary>
public static class HandoverFloorClaimRegistration
{
    public static ClaimRegistryBuilder RegisterHandoverFloorClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<HandoverFloorClaim>(b =>
        {
            var survival = b.Get<AbsorptionTheoremClaim>();                  // -2g<n_XY> (the survivor rate)
            var floor = b.Get<F50WeightOneDegeneracyPi2Inheritance>();       // the F50 off-diagonal floor =1
            var chainSolution = b.Get<CoherenceHorizonClaim>();             // the chain handover = Q*(N)
            return new HandoverFloorClaim(survival, floor, chainSolution);
        });
}
