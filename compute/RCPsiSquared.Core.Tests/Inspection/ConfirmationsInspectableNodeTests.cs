using RCPsiSquared.Core.Confirmations;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.Tests.Inspection;

public class ConfirmationsInspectableNodeTests
{
    [Fact]
    public void Build_HasExactlyOneChildPerRegistryEntry()
    {
        var node = ConfirmationsInspectableNode.Build();
        var children = node.Children.ToList();
        Assert.Equal(ConfirmationsRegistry.All.Count, children.Count);
    }

    [Fact]
    public void Build_ChildDisplayNamesMatchEntryNames()
    {
        var node = ConfirmationsInspectableNode.Build();
        var names = node.Children.Select(c => c.DisplayName).ToHashSet();
        foreach (var entry in ConfirmationsRegistry.All)
            Assert.Contains(entry.Name, names);
    }

    [Fact]
    public void Build_ChildSummaryCarriesDateMachinePredictedAndMeasured()
    {
        var node = ConfirmationsInspectableNode.Build();
        var entry = ConfirmationsRegistry.Lookup("f25_cusp_trajectory")!;
        var child = node.Children.Single(c => c.DisplayName == entry.Name);
        Assert.Equal(
            $"{entry.Date} {entry.Machine} {entry.Observable}: " +
            $"predicted {entry.PredictedValue} vs measured {entry.MeasuredValue}",
            child.Summary);
    }

    [Fact]
    public void Build_RootSummaryReportsCount()
    {
        var node = ConfirmationsInspectableNode.Build();
        Assert.Contains(ConfirmationsRegistry.All.Count.ToString(), node.Summary);
    }
}
