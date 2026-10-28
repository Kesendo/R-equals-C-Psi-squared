using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers the parentless F96 N=4-ring M/U slope table.</summary>
public static class F96BornSubdominantSlopesRegistration
{
    public static ClaimRegistryBuilder RegisterF96BornSubdominantSlopes(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F96BornSubdominantSlopes>(
            _ => new F96BornSubdominantSlopes());
}
