using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F1PalindromeOrbitPairing"/>; one typed
/// parent (<see cref="SymmetryFamilyInventory"/>). Parallels
/// <see cref="XGlobalChargeConjugationPairingRegistration"/>: the Π-orbit pairing subsumes
/// the X⊗N pairing, both hang off the same inventory parent.</summary>
public static class F1PalindromeOrbitPairingRegistration
{
    public static ClaimRegistryBuilder RegisterF1PalindromeOrbitPairing(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F1PalindromeOrbitPairing>(b =>
        {
            var inv = b.Get<SymmetryFamilyInventory>();
            return new F1PalindromeOrbitPairing(inv);
        });
}
