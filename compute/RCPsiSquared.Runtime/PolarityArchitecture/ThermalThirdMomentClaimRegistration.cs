using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="ThermalThirdMomentClaim"/> (F137, the third moment under cooling and heating, for
/// every traceless H) and of its child <see cref="ConditionedFieldPairingClaim"/> (the conditioned field X₀(1 + Z₁), whose
/// converse is the law's x¹³ coefficient). The law has no typed parent: F137 is not typed.</summary>
public static class ThermalThirdMomentClaimRegistration
{
    public static ClaimRegistryBuilder RegisterThermalThirdMomentClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ThermalThirdMomentClaim>(_ => new ThermalThirdMomentClaim());

    /// <summary>Requires <see cref="RegisterThermalThirdMomentClaim"/>.</summary>
    public static ClaimRegistryBuilder RegisterConditionedFieldPairingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ConditionedFieldPairingClaim>(b => new ConditionedFieldPairingClaim(b.Get<ThermalThirdMomentClaim>()));
}
