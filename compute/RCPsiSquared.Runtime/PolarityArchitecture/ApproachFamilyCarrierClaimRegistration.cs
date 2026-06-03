using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="ApproachFamilyCarrierClaim"/>: the cusp-approach family as a
/// carried c=2 decay reading. Four typed parent edges resolved from the registry:
/// <see cref="UniversalCarrierClaim"/> (the shared 4γ₀ carrier), <see cref="C2BareDoubledPtfClosedForm"/>
/// (the c=2 doubled-PTF kinship), <see cref="TwoReadingsClaim"/> (algebra vs dynamics), and
/// <see cref="F25CPsiBellPlusPi2Inheritance"/> (the Bell+ member). Requires all four parent registrations
/// upstream.</summary>
public static class ApproachFamilyCarrierClaimRegistration
{
    public static ClaimRegistryBuilder RegisterApproachFamilyCarrierClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ApproachFamilyCarrierClaim>(b =>
        {
            var carrier = b.Get<UniversalCarrierClaim>();
            var c2Ptf = b.Get<C2BareDoubledPtfClosedForm>();
            var twoReadings = b.Get<TwoReadingsClaim>();
            var f25 = b.Get<F25CPsiBellPlusPi2Inheritance>();
            return new ApproachFamilyCarrierClaim(carrier, c2Ptf, twoReadings, f25);
        });
}
