using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="FluxAreaFifthMomentClaim"/> (F137, the fifth centred trace under cooling and heating
/// as flux times the area of the rate triangle). No typed parent: F137 is not typed.</summary>
public static class FluxAreaFifthMomentClaimRegistration
{
    public static ClaimRegistryBuilder RegisterFluxAreaFifthMomentClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<FluxAreaFifthMomentClaim>(_ => new FluxAreaFifthMomentClaim());
}
