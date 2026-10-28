using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers the period-one cardioid identity under its one typed parent, the bilinear
/// maximum ¼ (<see cref="QuarterAsBilinearMaxvalClaim"/>, from RegisterPi2Family).</summary>
public static class F97CardioidHalfFixedPointRegistration
{
    public static ClaimRegistryBuilder RegisterF97CardioidHalfFixedPoint(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F97CardioidHalfFixedPoint>(
            b => new F97CardioidHalfFixedPoint(b.Get<QuarterAsBilinearMaxvalClaim>()));
}
