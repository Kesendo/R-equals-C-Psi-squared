using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class MotherSoftYParityOnePurityRegistrationTests
{
    [Fact]
    public void MotherSoftYParityOnePurity_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<MotherSoftYParityOnePurity>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_IncludesMotherSoftYParityOnePurity()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.YParityClaims, c => c is MotherSoftYParityOnePurity);
    }
}
