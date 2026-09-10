using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F163 uses F131 for order sorting and F89d for partner-block transport.</summary>
public static class RouteBN6A2UnfoldingClaimRegistration
{
    public static ClaimRegistryBuilder RegisterRouteBN6A2UnfoldingClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<RouteBN6A2UnfoldingClaim>(b => new RouteBN6A2UnfoldingClaim(
            b.Get<MirrorOrderSortingClaim>(), b.Get<F89CrossFoldSimilarityClaim>()));
}
