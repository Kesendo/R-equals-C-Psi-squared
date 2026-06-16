using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="SecondClockRegimeClaim"/> (the stitch): two typed parent edges, one per
/// regime it joins — <see cref="CoherenceHorizonClaim"/> (the EP regime, the dispersive √-EP coherence
/// horizon) and <see cref="StructuralCeilingClaim"/> (the CEILING regime, the degenerate-manifold 4/(m+1)
/// ceiling). Requires both registered first. The child is Tier1Candidate: the tier-inheritance invariant
/// caps it at its weaker parent <see cref="CoherenceHorizonClaim"/> (Tier1Candidate); the other parent is
/// Tier1Derived.</summary>
public static class SecondClockRegimeClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSecondClockRegimeClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<SecondClockRegimeClaim>(b =>
            new SecondClockRegimeClaim(
                b.Get<CoherenceHorizonClaim>(),
                b.Get<StructuralCeilingClaim>()));
}
