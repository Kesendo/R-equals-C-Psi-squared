using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="HeisenbergColourLawClaim"/> (Tier1Derived, F169). Two typed parents:
/// <see cref="PalindromeTwoEndCountClaim"/> (F158, whose far kernel a colour enters and whose section (h) covers jumps
/// along any axis) and <see cref="PalindromeComplementConnectionClaim"/> (F165, whose Theorems 1 to 3 settle the axes
/// one global turn makes letters). Resolution is topological.</summary>
public static class HeisenbergColourLawClaimRegistration
{
    public static ClaimRegistryBuilder RegisterHeisenbergColourLawClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<HeisenbergColourLawClaim>(b =>
            new HeisenbergColourLawClaim(b.Get<PalindromeTwoEndCountClaim>(), b.Get<PalindromeComplementConnectionClaim>()));
}
