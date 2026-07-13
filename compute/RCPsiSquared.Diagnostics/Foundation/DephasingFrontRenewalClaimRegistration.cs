using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="DephasingFrontRenewalClaim"/> (the exact renewal representation of the
/// watched walk; Tier1Derived). Two typed parents, both Tier1Derived and both registered earlier in the
/// BuildDefault chain: <see cref="AbsorptionTheoremClaim"/> (the uniform sector rate Γ = 4γ) +
/// <see cref="F2bXyChainSpectrumPi2Inheritance"/> (the clean single-particle propagator/band). Resolution is
/// topological (deferred construction), so registration order among siblings is immaterial as long as both
/// parents are registered somewhere in the chain.</summary>
public static class DephasingFrontRenewalClaimRegistration
{
    public static ClaimRegistryBuilder RegisterDephasingFrontRenewalClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<DephasingFrontRenewalClaim>(b =>
        {
            var rateLaw = b.Get<AbsorptionTheoremClaim>();                  // Γ = 4γ, the sector rate
            var band = b.Get<F2bXyChainSpectrumPi2Inheritance>();          // the clean propagator / band
            return new DephasingFrontRenewalClaim(rateLaw, band);
        });
}
