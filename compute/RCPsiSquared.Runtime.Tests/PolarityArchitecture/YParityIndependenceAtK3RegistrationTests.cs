using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class YParityIndependenceAtK3RegistrationTests
{
    [Fact]
    public void YParityIndependenceAtK3_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<YParityIndependenceAtK3>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_ContainsExactlyOne()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.NotNull(cubeMap);
        Assert.Single(cubeMap.YParityClaims);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaim_IsYParityIndependenceAtK3()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        var yParityClaim = Assert.Single(cubeMap.YParityClaims);
        Assert.IsType<YParityIndependenceAtK3>(yParityClaim);
    }
}
