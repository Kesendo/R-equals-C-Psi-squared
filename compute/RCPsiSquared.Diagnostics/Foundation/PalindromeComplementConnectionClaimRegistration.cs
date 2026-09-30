using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="PalindromeComplementConnectionClaim"/> (Tier1Derived). One typed parent,
/// <see cref="PalindromeTwoEndCountClaim"/> (F158): the complement connection is F158 read on the class
/// where every site carries one single-site Pauli jump, where both ends become diagonal problems on the
/// bitstrings. Resolution is topological, so this registration may sit anywhere relative to the parent's.</summary>
public static class PalindromeComplementConnectionClaimRegistration
{
    public static ClaimRegistryBuilder RegisterPalindromeComplementConnectionClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<PalindromeComplementConnectionClaim>(b =>
            new PalindromeComplementConnectionClaim(b.Get<PalindromeTwoEndCountClaim>()));
}
