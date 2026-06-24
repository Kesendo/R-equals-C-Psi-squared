using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.OpenArcs;

namespace RCPsiSquared.Core.Tests.Inspection;

public class OpenArcsInspectableNodeTests
{
    [Fact]
    public void Registry_HasTwentySevenEntries()
    {
        Assert.Equal(28, OpenArcsRegistry.All.Count);
    }

    [Fact]
    public void Registry_OpenAndRetiredCountsPartitionTheRegistry()
    {
        int retired = OpenArcsRegistry.All.Count(a => a.Status == OpenArcStatus.Retired);
        Assert.Equal(OpenArcsRegistry.All.Count, OpenArcsRegistry.OpenCount + retired);
    }

    [Fact]
    public void Registry_RetiredArcsCarryAReason()
    {
        Assert.All(
            OpenArcsRegistry.All.Where(a => a.Status == OpenArcStatus.Retired),
            a => Assert.False(string.IsNullOrWhiteSpace(a.RetiredReason)));
    }

    [Fact]
    public void OneDiagonalMirrorGroup_IsRetired_SpunOutToLinearS3()
    {
        var entry = OpenArcsRegistry.Lookup("one_diagonal_mirror_group");
        Assert.NotNull(entry);
        Assert.Equal(OpenArcStatus.Retired, entry!.Status);
        Assert.Contains("linear_s3_mirror_completion", entry.RetiredReason!);
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
        int retired = OpenArcsRegistry.All.Count(a => a.Status == OpenArcStatus.Retired);
        Assert.Contains($"{OpenArcsRegistry.OpenCount} open arc(s)", node.Summary);
        Assert.Contains($"{retired} retired", node.Summary);
    }
}
