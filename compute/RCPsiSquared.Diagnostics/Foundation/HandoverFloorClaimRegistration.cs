using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="HandoverFloorClaim"/> (the handover Q = the F50-floor condition;
/// chain = the coherence horizon's floor crossing, at Q*(N) for N = 2, 3 and just below it from N = 4, even ring = a
/// distinct (2,2) level crossing; Tier1Candidate). Typed parents (all
/// registered above): <see cref="AbsorptionTheoremClaim"/> (the -2g&lt;n_XY&gt; survivor rate),
/// <see cref="F50WeightOneDegeneracyPi2Inheritance"/> (the off-diagonal floor &lt;n_XY&gt;=1), and
/// <see cref="CoherenceHorizonClaim"/> (the horizon Q*(N) the chain solution meets or falls just below). Must register AFTER all three.</summary>
public static class HandoverFloorClaimRegistration
{
    public static ClaimRegistryBuilder RegisterHandoverFloorClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<HandoverFloorClaim>(b =>
        {
            var survival = b.Get<AbsorptionTheoremClaim>();                  // -2g<n_XY> (the survivor rate)
            var floor = b.Get<F50WeightOneDegeneracyPi2Inheritance>();       // the F50 off-diagonal floor =1
            var chainSolution = b.Get<CoherenceHorizonClaim>();             // the chain handover = the horizon's floor crossing (= Q*(N) at N=2,3 only)
            return new HandoverFloorClaim(survival, floor, chainSolution);
        });
}
