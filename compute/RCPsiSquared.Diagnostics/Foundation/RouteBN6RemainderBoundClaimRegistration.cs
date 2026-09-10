using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>One typed parent: F163 supplies the exact rank, semisimplicity, characteristic
/// factorization and unique repeated root as explicit premises of the ball computation.
/// No second parent is invented. The N=4 ball certificates share the contraction design but
/// run real balls on a scalar system, so they are a method precedent rather than a premise,
/// and typing them here would type a technique instead of a result.</summary>
public static class RouteBN6RemainderBoundClaimRegistration
{
    public static ClaimRegistryBuilder RegisterRouteBN6RemainderBoundClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<RouteBN6RemainderBoundClaim>(b => new RouteBN6RemainderBoundClaim(
            b.Get<RouteBN6A2UnfoldingClaim>()));
}
