using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>
/// Registers the parentless F94 named N=4-ring formula. Stable historical
/// identifiers remain; no dyadic or polarity edge is manufactured.
/// </summary>
public static class F94BornDeviationFourThirdsPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF94BornDeviationFourThirdsPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F94BornDeviationFourThirdsPi2Inheritance>(
            _ => new F94BornDeviationFourThirdsPi2Inheritance());
}
