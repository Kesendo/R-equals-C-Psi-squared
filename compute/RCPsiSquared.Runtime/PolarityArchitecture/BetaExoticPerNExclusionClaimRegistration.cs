using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="BetaExoticPerNExclusionClaim"/> (the β-exotic excluded
/// exactly at N = 5, N = 7, and N = 9, both R-parities, by certified disc-multiplicity readings;
/// <c>experiments/F89_SEED_EXISTENCE_REDUCTION.md</c>). One typed parent edge:
///
/// <list type="bullet">
///   <item><see cref="SeedExistenceCountingClaim"/> supplies the N − 1 endpoint-nullity surplus
///         at every odd N. Separately observed literal drops supply the loci; this claim pins their CHARACTER at the certified chain lengths (N = 5 and
///         N = 7, and N = 9): each is a defective EP2, never
///         the β-exotic. Neither the parent nor this edge proves literal drops for every odd N.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (an exact algebraic certificate over ℚ(i), one-way lift
/// from a prime certified at both ends of the q-axis, three empty reviews). The parent is Tier 1
/// derived.</para>
///
/// <para>Scope the edge does NOT carry: the parent's identity holds for every odd N, this claim only at
/// the chain lengths the certificate has been run at. A per-N certificate does not promote the parent's
/// open ink; it retires chain lengths one at a time.</para>
///
/// <para>Requires <see cref="SeedExistenceCountingClaimRegistration.RegisterSeedExistenceCountingClaim"/>
/// (the builder topo-resolves, so registration order is free).</para></summary>
public static class BetaExoticPerNExclusionClaimRegistration
{
    public static ClaimRegistryBuilder RegisterBetaExoticPerNExclusionClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<BetaExoticPerNExclusionClaim>(b =>
            new BetaExoticPerNExclusionClaim(
                b.Get<SeedExistenceCountingClaim>()));
}
