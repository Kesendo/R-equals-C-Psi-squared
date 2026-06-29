using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.OpenArcs;

namespace RCPsiSquared.Core.Tests.Inspection;

public class OpenArcsInspectableNodeTests
{
    [Fact]
    public void Registry_HasThirtyEntries()
    {
        Assert.Equal(30, OpenArcsRegistry.All.Count);   // +diabolic_over_higher_n (2026-06-27)
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
    public void StrangerDoor_IsRetired_FifthDoorHung()
    {
        var entry = OpenArcsRegistry.Lookup("stranger_door");
        Assert.NotNull(entry);
        Assert.Equal(OpenArcStatus.Retired, entry!.Status);
        Assert.Contains("provenance badge", entry.RetiredReason!);
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
        // Decoupled from any specific arc name: the exemplar is whichever arc is still Open,
        // so retiring an arc (e.g. birth_canal_surface, 2026-06-29) cannot re-break this test.
        var entry = OpenArcsRegistry.All.First(a => a.Status == OpenArcStatus.Open);
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
