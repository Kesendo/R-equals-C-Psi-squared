using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class KleinEightCellClaimRegistrationTests
{
    [Fact]
    public void KleinEightCellClaim_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<KleinEightCellClaim>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.Cubic3, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_Cubic3Claims_IncludesKleinEightCellClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.Cubic3Claims, c => c is KleinEightCellClaim);
    }

    [Fact]
    public void PolarityCubeMap_Cubic3Claims_ContainsExactlyOne()
    {
        // KleinEightCellClaim is currently the only Cubic3-axis Claim (first one,
        // Stage 2b). When future Cubic3 Claims land, bump this count.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Single(cubeMap.Cubic3Claims);
    }
}
