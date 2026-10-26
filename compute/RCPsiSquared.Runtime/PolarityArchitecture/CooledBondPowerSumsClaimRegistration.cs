using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="CooledBondPowerSumsClaim"/> (WHICH_BONDS_PAIR_UNDER_COOLING, the seventh
/// power sum and the rings). No typed parent: its parent in prose is F137's criterion, which is not
/// typed.</summary>
public static class CooledBondPowerSumsClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCooledBondPowerSumsClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CooledBondPowerSumsClaim>(_ => new CooledBondPowerSumsClaim());
}
