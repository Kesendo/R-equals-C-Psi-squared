using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public sealed class HardCellPureDTemplateRegistrationTests
{
    private static ClaimRegistry BuildRegistry() =>
        new ClaimRegistryBuilder()
            // Welle 7: HardCellPureDTemplate now ctor-takes KleinEightCellClaim as
            // typed Cubic3 parent; need to wire KleinFour + KleinEight + Pi2 chain.
            .RegisterPi2Family()
            .RegisterKleinEightCellClaim()
            .RegisterHardCellPureDTemplate()
            .Build();

    [Fact]
    public void Registered_Resolves_HardCellPureDTemplate()
    {
        var registry = BuildRegistry();
        var claim = registry.Get<HardCellPureDTemplate>();
        Assert.NotNull(claim);
        Assert.Equal(Tier.Tier1Candidate, claim.Tier);
    }

    [Fact]
    public void Registered_ClaimAnchorFile_Exists()
    {
        var registry = BuildRegistry();
        var claim = registry.Get<HardCellPureDTemplate>();
        var anchorFiles = claim.Anchor.Split(" + ", StringSplitOptions.RemoveEmptyEntries);
        Assert.NotEmpty(anchorFiles);
        Assert.Contains(anchorFiles, f => f.Contains("PROOF_F111"));
    }

    [Fact]
    public void Registered_YParityAxis_IsCorrect()
    {
        var registry = BuildRegistry();
        var claim = registry.Get<HardCellPureDTemplate>();
        Assert.Equal(Z2Axis.YParity, ((IZ2AxisClaim)claim).Z2Axis);
    }
}
