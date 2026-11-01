using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="UniversalCarrierClaim"/>: γ₀ as the
/// universal-reference rate-parameter, the special-relativity c analog that additionally
/// carries the calibration substrate. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="AbsorptionTheoremClaim"/>: the one-disagreement basis-cell cost 2γ₀.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: the dyadic-ladder anchors a₀, a_{−1}, a₃.</item>
/// </list>
///
/// <para>Tier consistency: all three Tier1Derived.</para>
///
/// <para>Requires upstream registrations:
/// <see cref="AbsorptionTheoremClaimRegistration.RegisterAbsorptionTheoremClaim"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class UniversalCarrierClaimRegistration
{
    public static ClaimRegistryBuilder RegisterUniversalCarrierClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<UniversalCarrierClaim>(b =>
        {
            var absorption = b.Get<AbsorptionTheoremClaim>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new UniversalCarrierClaim(absorption, ladder);
        });
}
