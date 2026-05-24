using System.Linq;
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
    public void PolarityCubeMap_YParityClaims_IncludesYParityIndependenceAtK3()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.NotNull(cubeMap);
        Assert.Contains(cubeMap.YParityClaims, c => c is YParityIndependenceAtK3);
    }

    [Fact]
    public void PolarityCubeMap_YParityIndependenceAtK3_AppearsExactlyOnce()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Equal(1, cubeMap.YParityClaims.Count(c => c is YParityIndependenceAtK3));
    }
}
