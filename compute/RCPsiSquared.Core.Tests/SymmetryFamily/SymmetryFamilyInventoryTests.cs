using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.SymmetryFamily;
using Xunit;

namespace RCPsiSquared.Core.Tests.SymmetryFamily;

public sealed class SymmetryFamilyInventoryTests
{
    [Fact]
    public void Tier_IsTier1Derived()
    {
        var inv = new SymmetryFamilyInventory();
        Assert.Equal(Tier.Tier1Derived, inv.Tier);
    }

    [Fact]
    public void DisplayName_NotEmpty()
    {
        var inv = new SymmetryFamilyInventory();
        Assert.False(string.IsNullOrWhiteSpace(inv.DisplayName));
    }

    [Fact]
    public void Summary_MentionsFamilySize()
    {
        var inv = new SymmetryFamilyInventory();
        Assert.Contains("family", inv.Summary, System.StringComparison.OrdinalIgnoreCase);
    }
}
