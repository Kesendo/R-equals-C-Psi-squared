using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="AntiFractionTwoBlindSpotsClaim"/>. Two typed parents, both in
/// Diagnostics/F87 and both Tier1Derived: <see cref="F89F87BreakPredictionFromF83"/> (the 0 end) and
/// <see cref="AntiFractionObstructionOrthogonalityClaim"/> (the 1/2 end). Connector is Tier1Derived
/// (5 >= 5). Requires both parents' registrations to have run.</summary>
public static class AntiFractionTwoBlindSpotsRegistration
{
    public static ClaimRegistryBuilder RegisterAntiFractionTwoBlindSpots(
        this ClaimRegistryBuilder builder) =>
        builder.Register<AntiFractionTwoBlindSpotsClaim>(b =>
            new AntiFractionTwoBlindSpotsClaim(
                b.Get<F89F87BreakPredictionFromF83>(),
                b.Get<AntiFractionObstructionOrthogonalityClaim>()));
}
