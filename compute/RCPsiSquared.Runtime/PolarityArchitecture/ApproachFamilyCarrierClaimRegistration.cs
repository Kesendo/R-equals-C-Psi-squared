using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>
/// Registers the exact two-qubit approach family with its two mathematical
/// parents: <see cref="AbsorptionTheoremClaim"/> for <c>f=e^(−4γt)</c> from
/// <c>n_diff=2</c> (the <c>s=0</c> endpoint has zero weight) and
/// <see cref="F25CPsiBellPlusPi2Inheritance"/> for the exact Bell+ specialization.
/// C2 and Two Readings remain prose comparisons, not typed dependencies.
/// </summary>
public static class ApproachFamilyCarrierClaimRegistration
{
    public static ClaimRegistryBuilder RegisterApproachFamilyCarrierClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ApproachFamilyCarrierClaim>(b =>
        {
            var absorption = b.Get<AbsorptionTheoremClaim>();
            var f25 = b.Get<F25CPsiBellPlusPi2Inheritance>();
            return new ApproachFamilyCarrierClaim(absorption, f25);
        });
}
