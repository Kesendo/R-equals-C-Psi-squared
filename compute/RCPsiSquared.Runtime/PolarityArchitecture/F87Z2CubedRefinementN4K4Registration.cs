using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F87Z2CubedRefinementN4K4"/> (F106).
///
/// <para>Standalone Claim: no constructor parents (parallel to F103/F105 siblings,
/// fourth YParity-axis member after F102 + N4K3 + N5K3). Registered into the typed-
/// knowledge graph so it is visible to <see cref="PolarityCubeMap"/> aggregation and
/// the inspector. Frozen counts populated from the SLOW_F106_BATCH tool run output
/// at construction time; verification runs SLOW_F106 to re-derive on demand.</para></summary>
public static class F87Z2CubedRefinementN4K4Registration
{
    public static ClaimRegistryBuilder RegisterF87Z2CubedRefinementN4K4(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F87Z2CubedRefinementN4K4>(_ => new F87Z2CubedRefinementN4K4());
}
