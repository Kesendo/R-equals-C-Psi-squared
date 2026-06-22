using RCPsiSquared.Core.F1;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="GoldenRouterClaim"/> (the Tier1Derived F116 golden/metallic
/// router). Two typed parents, both Tier1Derived (the strength-inheritance check is 5 ≥ 5 against each):
/// <see cref="F1PalindromeIdentity"/> (the global palindrome the router realizes locally for the Z-middle
/// ceiling class) and <see cref="WindowedConverseThresholdClaim"/> (the F87 two-reflection chiral spine
/// whose chiral driving F H F = −H the router's two-sided form distributes into P ≠ Q).
///
/// <para>Must register after <see cref="F1Family.F1FamilyRegistration.RegisterF1Family"/> (which provides
/// <see cref="F1PalindromeIdentity"/>) and after
/// <see cref="WindowedConverseThresholdClaimRegistration.RegisterWindowedConverseThresholdClaim"/>; the
/// builder errors with <c>MissingParent</c> otherwise.</para>
///
/// <para>NOTE: the open arc f116_golden_router_typed_claim suggested parenting on
/// <see cref="PalindromeSoftCertifierClaim"/>, but that is Tier1Candidate (4 &lt; 5) so it would violate
/// the parent ≥ child rule, and it is backwards (the certifier USES this router as a helper, so this claim
/// is logically upstream of it). That relationship is a see-cref in the Claim, not a parent edge.</para></summary>
public static class GoldenRouterClaimRegistration
{
    public static ClaimRegistryBuilder RegisterGoldenRouterClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<GoldenRouterClaim>(b =>
        {
            _ = b.Get<F1PalindromeIdentity>();                   // typed parent edge (parent before child)
            _ = b.Get<WindowedConverseThresholdClaim>();         // typed parent edge (parent before child)
            return new GoldenRouterClaim(
                b.Get<F1PalindromeIdentity>(),
                b.Get<WindowedConverseThresholdClaim>());
        });
        return builder;
    }
}
