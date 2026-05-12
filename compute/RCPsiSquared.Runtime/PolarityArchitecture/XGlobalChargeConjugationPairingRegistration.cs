using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="XGlobalChargeConjugationPairing"/>; one typed
/// parent (<see cref="SymmetryFamilyInventory"/>).</summary>
public static class XGlobalChargeConjugationPairingRegistration
{
    public static ClaimRegistryBuilder RegisterXGlobalChargeConjugationPairing(
        this ClaimRegistryBuilder builder) =>
        builder.Register<XGlobalChargeConjugationPairing>(b =>
        {
            var inv = b.Get<SymmetryFamilyInventory>();
            return new XGlobalChargeConjugationPairing(inv);
        });
}
