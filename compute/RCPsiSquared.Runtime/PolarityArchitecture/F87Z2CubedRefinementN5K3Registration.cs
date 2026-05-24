using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F87Z2CubedRefinementN5K3"/> (F105).
///
/// <para>Standalone Claim: no constructor parents (parallel to F103 N4K3 sibling and
/// to F102 as a third YParity-axis member). Registered into the typed-knowledge
/// graph so it is visible to <see cref="PolarityCubeMap"/> aggregation and the
/// inspector. Frozen counts are populated from the SLOW_F105_BATCH tool run output
/// at construction time; verification runs SLOW_F105 to re-derive on demand.</para></summary>
public static class F87Z2CubedRefinementN5K3Registration
{
    public static ClaimRegistryBuilder RegisterF87Z2CubedRefinementN5K3(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F87Z2CubedRefinementN5K3>(_ => new F87Z2CubedRefinementN5K3());
}
