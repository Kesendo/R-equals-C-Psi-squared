using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class HardCellYInversionPatternRegistrationTests
{
    [Fact]
    public void HardCellYInversionPattern_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<HardCellYInversionPattern>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_YParityClaims_IncludesHardCellYInversionPattern()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.YParityClaims, c => c is HardCellYInversionPattern);
    }

    [Fact]
    public void HardCellYInversionPattern_IsTier1Candidate()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<HardCellYInversionPattern>();
        Assert.Equal(RCPsiSquared.Core.Knowledge.Tier.Tier1Candidate, claim.Tier);
    }
}
