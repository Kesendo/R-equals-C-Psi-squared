using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F87Z2CubedRefinementN4K3RegistrationTests
{
    [Fact]
    public void F87Z2CubedRefinementN4K3_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F87Z2CubedRefinementN4K3>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_ContainsExactlyThree()
    {
        // F102 + F103 (N4K3) + F105 (N5K3) are the Z2Axis.YParity members in PolarityCubeMap; expect count == 3.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.NotNull(cubeMap);
        Assert.Equal(3, cubeMap.YParityClaims.Count);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_IncludesF87Z2CubedRefinementN4K3()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.YParityClaims, c => c is F87Z2CubedRefinementN4K3);
    }
}
