using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="TransportSeventhMomentClaim"/> (Theorem B of
/// PROOF_TRANSPORT_FORMULA_LOSS_AND_GAIN, F137: the seventh centred trace under cooling and heating at second order
/// in the rates). No typed parent: F137 is not typed.</summary>
public static class TransportSeventhMomentClaimRegistration
{
    public static ClaimRegistryBuilder RegisterTransportSeventhMomentClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<TransportSeventhMomentClaim>(_ => new TransportSeventhMomentClaim());
}
