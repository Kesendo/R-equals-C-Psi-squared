using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class TrulyYParityZeroPurityRegistrationTests
{
    [Fact]
    public void TrulyYParityZeroPurity_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<TrulyYParityZeroPurity>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_IncludesTrulyYParityZeroPurity()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.YParityClaims, c => c is TrulyYParityZeroPurity);
    }
}
