using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.OpenArcs;

namespace RCPsiSquared.Core.Tests.Inspection;

public class OpenArcsInspectableNodeTests
{
    [Fact]
    public void Registry_HasFourteenEntries()
    {
        Assert.Equal(14, OpenArcsRegistry.All.Count);
    }

    [Fact]
    public void Registry_AllEntriesAreCurrentlyOpen()
    {
        Assert.All(OpenArcsRegistry.All, a => Assert.Equal(OpenArcStatus.Open, a.Status));
        Assert.Equal(OpenArcsRegistry.All.Count, OpenArcsRegistry.OpenCount);
    }

    [Fact]
    public void Lookup_RoundtripsByName()
    {
        var entry = OpenArcsRegistry.Lookup("witness_coverage");
        Assert.NotNull(entry);
        Assert.Equal("witness_coverage", entry!.Name);
        Assert.Equal("2026-06-11", entry.Opened);
    }

    [Fact]
    public void Lookup_UnknownReturnsNull()
    {
        Assert.Null(OpenArcsRegistry.Lookup("no_such_arc"));
    }

    [Fact]
    public void Build_HasExactlyOneChildPerRegistryEntry()
    {
        var node = OpenArcsInspectableNode.Build();
        Assert.Equal(OpenArcsRegistry.All.Count, node.Children.Count());
    }

    [Fact]
    public void Build_ChildSummaryCarriesParkedAndNextStep()
    {
        var node = OpenArcsInspectableNode.Build();
        var entry = OpenArcsRegistry.Lookup("birth_canal_surface")!;
        var child = node.Children.Single(c => c.DisplayName == entry.Name);
        Assert.Equal($"{entry.Opened} parked: {entry.ParkedAt} -> next: {entry.NextStep}", child.Summary);
    }

    [Fact]
    public void Build_RootSummaryReportsOpenAndRetiredCounts()
    {
        var node = OpenArcsInspectableNode.Build();
        Assert.Contains($"{OpenArcsRegistry.OpenCount} open arc(s)", node.Summary);
        Assert.Contains("0 retired", node.Summary);
    }
}
