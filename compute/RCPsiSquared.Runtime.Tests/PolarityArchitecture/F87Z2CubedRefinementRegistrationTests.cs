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
        // Before F103: YParityClaims.Count == 1 (F102 only).
        // After F103: YParityClaims.Count == 2 (F102 + F103).
        // ZGlobalMirrorRefinement is on the YParity axis per its Z2Axis property, but
        // it is registered as a separate IZ2AxisClaim entry; the count counted here
        // is the count of YParity-axis entries from F-numbered Claims only? No: the
        // count is the full YParityClaims list from PolarityCubeMap which includes
        // every IZ2AxisClaim with Z2Axis == YParity. Verify by re-reading
        // PolarityCubeMapRegistration.cs to confirm exactly which Claims are wired.
        // At base SHA (after F102 and ZGlobalMirrorRefinement landed): YParity list
        // = [YParityIndependenceAtK3]. ZGlobalMirrorRefinement has Z2Axis.BitB per
        // its implementation; it is NOT a YParity-axis member.
        // After this task lands F103: YParity list = [YParityIndependenceAtK3,
        // F87Z2CubedRefinement] → count 2.
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
