using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="TwoReadingsClaim"/>: the meta-pattern that any
/// object on the doubled operator-space d² = 4^N admits exactly two coordinate readings of
/// one underlying object (bra/ket, number/angle, M and Π·M·Π⁻¹, inside/outside). Single
/// parent edge:
///
/// <list type="bullet">
///   <item><see cref="PolynomialFoundationClaim"/>: the d=2 polynomial trunk d²−2d=0 whose
///         solution pair {0, 2} is the source of every two-readings instance below.</item>
/// </list>
///
/// <para>Tier consistency: both Tier1Derived.</para>
///
/// <para>Requires upstream registration: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers PolynomialFoundationClaim).</para></summary>
public static class TwoReadingsClaimRegistration
{
    public static ClaimRegistryBuilder RegisterTwoReadingsClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<TwoReadingsClaim>(b =>
        {
            var polynomial = b.Get<PolynomialFoundationClaim>();
            return new TwoReadingsClaim(polynomial);
        });
}
