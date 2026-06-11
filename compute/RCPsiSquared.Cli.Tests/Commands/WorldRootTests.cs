using RCPsiSquared.Cli;
using RCPsiSquared.Cli.Commands;
using RCPsiSquared.Core.Confirmations;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.OpenArcs;

namespace RCPsiSquared.Cli.Tests.Commands;

public class WorldRootTests
{
    private static IInspectable BuildWorld()
    {
        var ctx = new InspectRootContext(
            Parser: new ArgParser(new[] { "--N", "2", "--n", "1", "--gamma", "0.05" }),
            N: 2, WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var world = InspectCommand.Catalog.Single(e => e.Name == "world");
        return world.Factory(ctx);
    }

    [Fact]
    public void World_HasFourSectionChildren_RootsClaimsConfirmationsArcs()
    {
        var children = BuildWorld().Children.ToList();
        var names = children.Select(c => c.DisplayName).ToList();
        Assert.Equal(new[] { "roots", "claims", "confirmations", "arcs" }, names);
    }

    [Fact]
    public void World_ArcsSection_HasExactlyRegistryCountChildren()
    {
        var arcs = BuildWorld().Children.Single(c => c.DisplayName == "arcs");
        Assert.Equal(OpenArcsRegistry.All.Count, arcs.Children.Count());
    }

    [Fact]
    public void World_RootsSection_HasOneNodePerCatalogEntryExceptWorldItself()
    {
        var roots = BuildWorld().Children.Single(c => c.DisplayName == "roots");
        int expected = InspectCommand.Catalog.Count(e => e.Name != "world");
        Assert.Equal(expected, roots.Children.Count());
    }

    [Fact]
    public void World_ConfirmationsSection_HasExactlyRegistryCountChildren()
    {
        var confirmations = BuildWorld().Children.Single(c => c.DisplayName == "confirmations");
        Assert.Equal(ConfirmationsRegistry.All.Count, confirmations.Children.Count());
    }

    [Fact]
    public void World_Summary_ReportsRootClaimAndConfirmationCounts()
    {
        var summary = BuildWorld().Summary;
        int roots = InspectCommand.Catalog.Count(e => e.Name != "world");
        Assert.Contains($"{roots} roots", summary);
        Assert.Contains($"{ConfirmationsRegistry.All.Count} confirmations", summary);
        Assert.Contains("claims", summary);
        Assert.Contains($"{OpenArcsRegistry.OpenCount} open arcs", summary);
    }

    [Fact]
    public void World_DepthTwoWalk_DoesNotEnumerateAnyCatalogRootFactory()
    {
        // The world default render depth is 2: world(0) -> sections(1) -> root/tier/confirmation
        // nodes(2). A depth-2 walk must reach the per-root nodes but never their children
        // (depth 3), so no heavy catalog factory fires. We prove it by walking to depth 2 and
        // asserting no FourModeEffective node (the fourmode factory's product, which only
        // appears one level below its lazy root node) is present.
        var depthTwo = BuildWorld().Walk(maxDepth: 2).ToList();

        var rootsSection = depthTwo.Single(n => n.DisplayName == "roots");
        // The lazy per-root nodes are present at depth 2 (their names + descriptions).
        Assert.Contains(rootsSection.Children, c => c.DisplayName == "fourmode");

        // ... but nothing the fourmode factory would produce has been materialized.
        Assert.DoesNotContain(depthTwo, n => n.DisplayName.Contains("FourModeEffective"));
    }
}
