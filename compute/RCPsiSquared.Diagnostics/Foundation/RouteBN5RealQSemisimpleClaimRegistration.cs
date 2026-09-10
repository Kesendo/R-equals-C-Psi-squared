using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="RouteBN5RealQSemisimpleClaim"/> (F164; Tier1Derived), with zero typed
/// parents. The exact inputs the proof consumes, the irreducibility of A2_O and the layer split, are carried
/// by anchor rather than as a typed edge: they live in a Python certificate
/// (<c>simulations/o2b_gcd_certificate.py</c>) with no claim of their own, and inventing one here would
/// type a producer instead of a result. The N=4 octic and F163 at N=6 are prior-work siblings, not
/// dependencies, so neither becomes an ancestor.</summary>
public static class RouteBN5RealQSemisimpleClaimRegistration
{
    public static ClaimRegistryBuilder RegisterRouteBN5RealQSemisimpleClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<RouteBN5RealQSemisimpleClaim>(_ => new RouteBN5RealQSemisimpleClaim());
}
