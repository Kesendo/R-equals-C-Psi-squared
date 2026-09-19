using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers the parentless period-one cardioid identity.</summary>
public static class F97CardioidHalfFixedPointPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF97CardioidHalfFixedPointPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F97CardioidHalfFixedPointPi2Inheritance>(
            _ => new F97CardioidHalfFixedPointPi2Inheritance());
}
