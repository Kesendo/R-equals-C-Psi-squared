using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="BlindSeatSectorFactorisationClaim"/> (F162; Tier1Derived). TWO typed parents,
/// both Tier1Derived, so the tier rule holds with room to spare. <see cref="SeatCutBlindnessClaim"/> (F157) owns
/// the locus polynomial, the comb, the pole rule and the node modulus as DEFINITIONS and states no sign and no
/// leading coefficient, which is exactly what this claim adds; <see cref="CrackedRingExactCurveClaim"/> (F160)
/// owns the Cassini identity the congruence turns on, cited by section (i) and not repeated there. Resolution is
/// topological, so this registration may sit anywhere after the families that supply the two parents.</summary>
public static class BlindSeatSectorFactorisationClaimRegistration
{
    public static ClaimRegistryBuilder RegisterBlindSeatSectorFactorisationClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<BlindSeatSectorFactorisationClaim>(b =>
            new BlindSeatSectorFactorisationClaim(
                b.Get<SeatCutBlindnessClaim>(),
                b.Get<CrackedRingExactCurveClaim>()));
}
