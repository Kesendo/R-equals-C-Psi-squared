using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers the period-one cardioid identity under its one typed parent, the bilinear
/// maximum ¼ (<see cref="QuarterAsBilinearMaxvalClaim"/>, from RegisterPi2Family).</summary>
public static class F97CardioidHalfFixedPointPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF97CardioidHalfFixedPointPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F97CardioidHalfFixedPointPi2Inheritance>(
            b => new F97CardioidHalfFixedPointPi2Inheritance(b.Get<QuarterAsBilinearMaxvalClaim>()));
}
