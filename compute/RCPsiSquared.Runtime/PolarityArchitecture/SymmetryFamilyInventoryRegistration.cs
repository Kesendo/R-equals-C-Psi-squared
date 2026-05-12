using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="SymmetryFamilyInventory"/>: top-level
/// aggregator-Claim with no typed parents (foundational inventory).</summary>
public static class SymmetryFamilyInventoryRegistration
{
    public static ClaimRegistryBuilder RegisterSymmetryFamilyInventory(
        this ClaimRegistryBuilder builder) =>
        builder.Register<SymmetryFamilyInventory>(_ => new SymmetryFamilyInventory());
}
