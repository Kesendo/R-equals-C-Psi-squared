using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="StarFrozenSeamClaim"/> (the star's frozen seam, the survivor-level reading
/// of "the star has no coherence horizon"): two typed parent edges — <see cref="StructuralCeilingClaim"/>
/// (Tier1Derived, the g2 = 4/(N−1) threshold + the commutant mechanism it reuses) and
/// <see cref="SecondClockRegimeClaim"/> (Tier1Candidate, the {0,2}/second-clock regime map whose star/GRADUAL
/// case this sharpens). Requires both registered first (SecondClockRegimeClaim is registered immediately above
/// in the factory). The child is Tier1Candidate: the tier-inheritance invariant caps it at its weaker parent
/// <see cref="SecondClockRegimeClaim"/>.</summary>
public static class StarFrozenSeamClaimRegistration
{
    public static ClaimRegistryBuilder RegisterStarFrozenSeamClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<StarFrozenSeamClaim>(b =>
            new StarFrozenSeamClaim(
                b.Get<StructuralCeilingClaim>(),
                b.Get<SecondClockRegimeClaim>()));
}
