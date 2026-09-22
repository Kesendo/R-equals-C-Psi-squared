using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Composes the owned slow cluster, the physical partial-trace bound and the chiral gauge.
/// Blind-seat geometry and light migration remain ancestors through the spectral claim.</summary>
public static class MissingPhaseSlowReadoutClaimRegistration
{
    public static ClaimRegistryBuilder RegisterMissingPhaseSlowReadoutClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<MissingPhaseSlowReadoutClaim>(context => new MissingPhaseSlowReadoutClaim(
            context.Get<MissingPhaseRelaxationScaleClaim>(),
            context.Get<F70DeltaNSelectionRulePi2Inheritance>(),
            context.Get<ChiralKClaim>()));
}
