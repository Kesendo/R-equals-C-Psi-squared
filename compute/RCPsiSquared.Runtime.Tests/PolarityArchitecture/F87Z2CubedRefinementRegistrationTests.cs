using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F87Z2CubedRefinementRegistrationTests
{
    [Fact]
    public void F87Z2CubedRefinement_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F87Z2CubedRefinement>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_ContainsExactlyTwo()
    {
        // F102 + F103 are the only Z2Axis.YParity members in PolarityCubeMap; expect count == 2.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.NotNull(cubeMap);
        Assert.Equal(2, cubeMap.YParityClaims.Count);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_IncludesF87Z2CubedRefinement()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.YParityClaims, c => c is F87Z2CubedRefinement);
    }
}
