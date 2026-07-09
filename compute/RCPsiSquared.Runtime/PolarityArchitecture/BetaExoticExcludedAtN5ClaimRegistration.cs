using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="BetaExoticExcludedAtN5Claim"/> (the β-exotic excluded
/// exactly at N = 5, both R-parities, by the certified disc-multiplicity reading;
/// <c>experiments/F89_SEED_EXISTENCE_REDUCTION.md</c>). One typed parent edge:
///
/// <list type="bullet">
///   <item><see cref="SeedExistenceCountingClaim"/>: it forces the N − 1 real-to-complex count-drops
///         at every odd N. This claim pins their CHARACTER at N = 5: each is a defective EP2, never
///         the β-exotic, which is precisely the open item that claim's own scope note names.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (an exact algebraic certificate over ℚ(i), one-way lift
/// from a prime certified at both ends of the q-axis, three empty reviews). The parent is Tier 1
/// derived.</para>
///
/// <para>Scope the edge does NOT carry: the parent's identity holds for every odd N, this claim only
/// for N = 5. A per-N certificate does not promote the parent's open ink; it retires one N.</para>
///
/// <para>Requires <see cref="SeedExistenceCountingClaimRegistration.RegisterSeedExistenceCountingClaim"/>
/// (the builder topo-resolves, so registration order is free).</para></summary>
public static class BetaExoticExcludedAtN5ClaimRegistration
{
    public static ClaimRegistryBuilder RegisterBetaExoticExcludedAtN5Claim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<BetaExoticExcludedAtN5Claim>(b =>
            new BetaExoticExcludedAtN5Claim(
                b.Get<SeedExistenceCountingClaim>()));
}
