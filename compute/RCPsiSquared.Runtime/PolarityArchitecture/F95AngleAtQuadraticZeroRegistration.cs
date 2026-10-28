using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>
/// Registers the parentless positive-b closed form. Its domain is finite b &gt; 0;
/// currently registered callers use b &gt; 0.
/// </summary>
public static class F95AngleAtQuadraticZeroRegistration
{
    public static ClaimRegistryBuilder RegisterF95AngleAtQuadraticZero(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F95AngleAtQuadraticZero>(
            _ => new F95AngleAtQuadraticZero());
}
