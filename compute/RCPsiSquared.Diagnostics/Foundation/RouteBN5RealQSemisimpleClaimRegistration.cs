using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="RouteBN5RealQSemisimpleClaim"/> (F164; Tier1Derived). ONE typed parent,
/// <see cref="F89Path3OcticEpClaim"/> (Tier1Derived): the same (1,2) coherence block one chain shorter,
/// the same character question, the same verdict reached by compression and live diagnostics rather than
/// by exact arithmetic, so the tier rule holds and this claim is that reading's exact-arithmetic form at
/// N=5. The exact inputs the proof consumes, the irreducibility of A2_O and the layer split, are carried
/// by anchor rather than as a typed edge: they live in a Python certificate
/// (<c>simulations/o2b_gcd_certificate.py</c>) with no claim of their own, and inventing one here would
/// type a producer instead of a result. F163 is the sibling at N=6 and is likewise not an edge, having no
/// claim in Core at all. Resolution is topological, so this registration may sit anywhere after the family
/// that supplies the parent.</summary>
public static class RouteBN5RealQSemisimpleClaimRegistration
{
    public static ClaimRegistryBuilder RegisterRouteBN5RealQSemisimpleClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<RouteBN5RealQSemisimpleClaim>(b =>
            new RouteBN5RealQSemisimpleClaim(b.Get<F89Path3OcticEpClaim>()));
}
