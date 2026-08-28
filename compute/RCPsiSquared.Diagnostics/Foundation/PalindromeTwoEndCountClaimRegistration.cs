using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="PalindromeTwoEndCountClaim"/> (F158; Tier1Derived). One typed
/// parent, <see cref="F1PalindromeIdentity"/>: F1 is the OPERATOR identity Π·L·Π⁻¹ = −L − 2Σγ·I
/// and F158 is the SPECTRAL statement, so F1 implies F158 wherever F1 holds and not conversely.
/// F1's own validity node already draws that distinction and scopes it to F138, which is why the
/// edge runs this way. Resolution is topological, so this registration may sit anywhere in the
/// chain relative to the family that supplies the parent.</summary>
public static class PalindromeTwoEndCountClaimRegistration
{
    public static ClaimRegistryBuilder RegisterPalindromeTwoEndCountClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<PalindromeTwoEndCountClaim>(b =>
            new PalindromeTwoEndCountClaim(b.Get<F1PalindromeIdentity>()));
}
