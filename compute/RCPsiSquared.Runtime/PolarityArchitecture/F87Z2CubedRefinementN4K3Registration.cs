using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="F87Z2CubedRefinementN4K3"/> (F103).
///
/// <para>Standalone Claim: no constructor parents (parallel to F102 as a second
/// YParity-axis seed, this time an empirical-anchor refinement of F87 in Z₂³
/// at N=4 k=3). Registered into the typed-knowledge graph so it is visible to
/// <see cref="PolarityCubeMap"/> aggregation and the inspector.</para></summary>
public static class F87Z2CubedRefinementN4K3Registration
{
    public static ClaimRegistryBuilder RegisterF87Z2CubedRefinementN4K3(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F87Z2CubedRefinementN4K3>(_ => new F87Z2CubedRefinementN4K3());
}
