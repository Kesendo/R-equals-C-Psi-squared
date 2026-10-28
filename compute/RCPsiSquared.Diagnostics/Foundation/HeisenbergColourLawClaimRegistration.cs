using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="HeisenbergColourLawClaim"/> (Tier1Derived, F169). One typed parent:
/// <see cref="PalindromeTwoEndCountClaim"/> (F158, whose far kernel a colour enters and whose section (h) covers jumps
/// along any axis). F165 is adjacent, not a parent: its Theorems 1 to 3 give the same answer at the axes one global
/// turn makes letters, a second route the proof does not rest on. Resolution is topological.</summary>
public static class HeisenbergColourLawClaimRegistration
{
    public static ClaimRegistryBuilder RegisterHeisenbergColourLawClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<HeisenbergColourLawClaim>(b =>
            new HeisenbergColourLawClaim(b.Get<PalindromeTwoEndCountClaim>()));
}
