using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="PolynomialDiscriminantAnchorClaim"/>:
/// makes the connection "polynomial discriminant of d²−2d=0 = 4 = a_{−1} =
/// d² for d=2" explicit in the typed-knowledge runtime. Three typed parent
/// edges:
///
/// <list type="bullet">
///   <item><see cref="PolynomialFoundationClaim"/>: the polynomial itself.</item>
///   <item><see cref="QubitDimensionalAnchorClaim"/>: d=2 root of the polynomial.</item>
///   <item><see cref="Pi2DyadicLadderClaim"/>: a_{−1} = 4 (negative-index access
///         enabled by ZERO_IS_THE_MIRROR; Tom 2026-05-09).</item>
/// </list>
///
/// <para>Tier consistency: pure composition; all three parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class PolynomialDiscriminantAnchorRegistration
{
    public static ClaimRegistryBuilder RegisterPolynomialDiscriminantAnchor(
        this ClaimRegistryBuilder builder) =>
        builder.Register<PolynomialDiscriminantAnchorClaim>(b =>
        {
            var polynomial = b.Get<PolynomialFoundationClaim>();
            var qubit = b.Get<QubitDimensionalAnchorClaim>();
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new PolynomialDiscriminantAnchorClaim(polynomial, qubit, ladder);
        });
}
