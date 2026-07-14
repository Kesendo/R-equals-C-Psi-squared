using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="CrossTripleOrthogonalityClaim"/> (F127; Tier1Candidate, the
/// code-trust caveat named). One typed parent: <see cref="SeedExistenceCountingClaim"/>, which
/// carries the theorem's structural inputs (the (1,2) pencil and its −2/−6 rungs, and the
/// Conway-Jones fusion-resonance criterion whose doubling is the variety V). The twinning
/// claims are CONSUMERS of F127, not parents, so no edge points at them. Resolution is
/// topological, so registration order among siblings is immaterial.</summary>
public static class CrossTripleOrthogonalityClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCrossTripleOrthogonalityClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<CrossTripleOrthogonalityClaim>(b =>
            new CrossTripleOrthogonalityClaim(b.Get<SeedExistenceCountingClaim>()));
}
