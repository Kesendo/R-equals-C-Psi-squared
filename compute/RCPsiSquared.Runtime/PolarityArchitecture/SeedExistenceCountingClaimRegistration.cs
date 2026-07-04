using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="SeedExistenceCountingClaim"/> (the counting theorem
/// r(0⁺) − r(∞) = N − 1 for odd N on the (1,2) block,
/// <c>experiments/F89_SEED_EXISTENCE_REDUCTION.md</c>). Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="AbsorptionTheoremClaim"/>: the rate law behind the dephasing diagonal
///         A = −2·diag(n_diff), whose two values −2/−6 are the rungs the theorem counts.</item>
///   <item><see cref="ChiralKClaim"/>: the chiral single-particle pairing λ_{N+1−k} = −λ_k,
///         the mirror powering the 3-to-1 triple↔resonance bijection of the (N1′)
///         ordering-sector theorem and the per-block spectral inheritance.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (three proved counting lemmas, two adversarial
/// reviews, gates in the Python verifier and the live witness). Both parents Tier 1
/// derived.</para>
///
/// <para>Requires <see cref="AbsorptionTheoremClaimRegistration.RegisterAbsorptionTheoremClaim"/>
/// and <see cref="ChiralKClaimRegistration.RegisterChiralK"/> (the builder topo-resolves, so
/// registration order is free).</para></summary>
public static class SeedExistenceCountingClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSeedExistenceCountingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<SeedExistenceCountingClaim>(b =>
            new SeedExistenceCountingClaim(
                b.Get<AbsorptionTheoremClaim>(),
                b.Get<ChiralKClaim>()));
}
