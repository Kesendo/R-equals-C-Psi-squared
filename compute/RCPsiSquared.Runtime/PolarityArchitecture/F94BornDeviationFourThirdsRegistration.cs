using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>
/// Registers the parentless F94 named N=4-ring formula. Stable historical
/// identifiers remain; no dyadic or polarity edge is manufactured.
/// </summary>
public static class F94BornDeviationFourThirdsRegistration
{
    public static ClaimRegistryBuilder RegisterF94BornDeviationFourThirds(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F94BornDeviationFourThirds>(
            _ => new F94BornDeviationFourThirds());
}
