using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="YParityIndependenceAtK3"/> (F102).
///
/// <para>Standalone Claim: no constructor parents (F102 is the first YParity-axis
/// seed). Registered into the typed-knowledge graph so it is visible to
/// <see cref="PolarityCubeMap"/> aggregation and the inspector.</para></summary>
public static class YParityIndependenceAtK3Registration
{
    public static ClaimRegistryBuilder RegisterYParityIndependenceAtK3(
        this ClaimRegistryBuilder builder) =>
        builder.Register<YParityIndependenceAtK3>(_ => new YParityIndependenceAtK3());
}
