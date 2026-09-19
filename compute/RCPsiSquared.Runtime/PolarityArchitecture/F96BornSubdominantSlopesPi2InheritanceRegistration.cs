using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers the parentless F96 N=4-ring M/U slope table.</summary>
public static class F96BornSubdominantSlopesPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF96BornSubdominantSlopesPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F96BornSubdominantSlopesPi2Inheritance>(
            _ => new F96BornSubdominantSlopesPi2Inheritance());
}
