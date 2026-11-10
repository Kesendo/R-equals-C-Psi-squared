using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="SeventhMomentDegreeFourClaim"/> (Theorems 1 and 2 of
/// PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE, F137: the seventh centred trace under cooling and heating at fourth order
/// in the rates, a sum over the triangles) and of <see cref="EqualTrianglePairingClaim"/> (its Theorem 3: flux
/// reversal on the equal triangle).
///
/// <para>Requires: <see cref="FluxAreaFifthMomentClaimRegistration.RegisterFluxAreaFifthMomentClaim"/>,
/// <see cref="TransportSeventhMomentClaimRegistration.RegisterTransportSeventhMomentClaim"/>; for the second also the first,
/// <see cref="RegisterSeventhMomentDegreeFourClaim"/>.</para></summary>
public static class SeventhMomentDegreeFourClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSeventhMomentDegreeFourClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<SeventhMomentDegreeFourClaim>(b =>
            new SeventhMomentDegreeFourClaim(b.Get<FluxAreaFifthMomentClaim>(), b.Get<TransportSeventhMomentClaim>()));

    public static ClaimRegistryBuilder RegisterEqualTrianglePairingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<EqualTrianglePairingClaim>(b =>
            new EqualTrianglePairingClaim(b.Get<FluxAreaFifthMomentClaim>(), b.Get<SeventhMomentDegreeFourClaim>()));
}
