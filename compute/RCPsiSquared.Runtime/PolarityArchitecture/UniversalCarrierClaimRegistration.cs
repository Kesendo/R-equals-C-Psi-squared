using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="UniversalCarrierClaim"/>: γ₀ as the
/// universal-reference rate-parameter, the special-relativity c analog that additionally
/// carries the calibration substrate. Three parent edges:
///
/// <list type="bullet">
///   <item><see cref="AbsorptionTheoremClaim"/>: the one-disagreement basis-cell cost 2γ₀.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: the dyadic-ladder anchors a₀, a_{−1}, a₃.</item>
///   <item><see cref="PolynomialDiscriminantAnchorClaim"/>: the value a_{−1} = 4 of
///         the four in t_peak = 1/(4γ₀), a four that is the rung-2 rate 2γ₀·2.</item>
/// </list>
///
/// <para>Tier consistency: all four Tier1Derived.</para>
///
/// <para>Requires upstream registrations:
/// <see cref="AbsorptionTheoremClaimRegistration.RegisterAbsorptionTheoremClaim"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="PolynomialDiscriminantAnchorRegistration.RegisterPolynomialDiscriminantAnchor"/>.</para></summary>
public static class UniversalCarrierClaimRegistration
{
    public static ClaimRegistryBuilder RegisterUniversalCarrierClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<UniversalCarrierClaim>(b =>
        {
            var absorption = b.Get<AbsorptionTheoremClaim>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var discriminant = b.Get<PolynomialDiscriminantAnchorClaim>();
            return new UniversalCarrierClaim(absorption, ladder, discriminant);
        });
}
