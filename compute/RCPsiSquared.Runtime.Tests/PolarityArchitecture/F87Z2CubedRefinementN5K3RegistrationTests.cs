using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F87Z2CubedRefinementN5K3RegistrationTests
{
    [Fact]
    public void F87Z2CubedRefinementN5K3_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F87Z2CubedRefinementN5K3>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
        Assert.Equal(5, claim.N);
        Assert.Equal(3, claim.K);
        Assert.Equal(294, claim.TotalPairs);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_IncludesF87Z2CubedRefinementN5K3()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.YParityClaims, c => c is F87Z2CubedRefinementN5K3);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_F87Z2CubedRefinementN5K3_AppearsExactlyOnce()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Equal(1, cubeMap.YParityClaims.Count(c => c is F87Z2CubedRefinementN5K3));
    }
}
