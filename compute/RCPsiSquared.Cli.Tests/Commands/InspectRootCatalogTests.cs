using RCPsiSquared.Cli;
using RCPsiSquared.Cli.Commands;

namespace RCPsiSquared.Cli.Tests.Commands;

public class InspectRootCatalogTests
{
    // The roots that were hardcoded in the InspectCommand switch before the catalog refactor.
    private static readonly string[] PreviouslyHardcodedRoots =
    {
        "f71", "f1", "f87", "pi2", "mirror", "flow", "between",
        "fourmode", "f86", "c2hwhm", "c2cpsi", "c2cpsi-scan",
    };

    [Fact]
    public void Catalog_ContainsEveryPreviouslyHardcodedRoot()
    {
        var names = InspectCommand.Catalog.Select(e => e.Name).ToHashSet();
        foreach (var root in PreviouslyHardcodedRoots)
            Assert.Contains(root, names);
    }

    [Fact]
    public void Catalog_AddsTheWorldRoot()
    {
        var names = InspectCommand.Catalog.Select(e => e.Name).ToHashSet();
        Assert.Contains("world", names);
    }

    [Fact]
    public void Catalog_DefaultRootIsFourmode_WithDepthFour()
    {
        var fourmode = InspectCommand.Catalog.Single(e => e.Name == "fourmode");
        Assert.Equal(4, fourmode.DefaultDepth);
    }

    [Fact]
    public void Catalog_F86KeepsDefaultDepthOne()
    {
        var f86 = InspectCommand.Catalog.Single(e => e.Name == "f86");
        Assert.Equal(1, f86.DefaultDepth);
    }

    [Fact]
    public void Catalog_WorldDefaultDepthIsTwo()
    {
        var world = InspectCommand.Catalog.Single(e => e.Name == "world");
        Assert.Equal(2, world.DefaultDepth);
    }

    [Fact]
    public void Catalog_NamesAreUnique()
    {
        var names = InspectCommand.Catalog.Select(e => e.Name).ToList();
        Assert.Equal(names.Count, names.Distinct().Count());
    }

    [Fact]
    public void Catalog_EveryEntryHasADescription()
    {
        Assert.All(InspectCommand.Catalog, e => Assert.False(string.IsNullOrWhiteSpace(e.Description)));
    }

    [Fact]
    public void Catalog_AddsTheQuditRoot_WhichDoesNotRequireN()
    {
        var qudit = InspectCommand.Catalog.Single(e => e.Name == "qudit");
        Assert.Equal("F121 qudit partial palindrome, recomputed live", qudit.Description);
        Assert.False(qudit.RequiresN);
    }

    [Fact]
    public void Catalog_QuditFactory_BuildsTheLiveWitness_DefaultD3N2()
    {
        var qudit = InspectCommand.Catalog.Single(e => e.Name == "qudit");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = qudit.Factory(ctx);
        Assert.Contains("d=3", root.Summary);
        Assert.Contains("ceiling met", root.Summary);
    }

    [Fact]
    public void Catalog_WorldAndQudit_DoNotRequireN()
    {
        Assert.False(InspectCommand.Catalog.Single(e => e.Name == "world").RequiresN);
        Assert.False(InspectCommand.Catalog.Single(e => e.Name == "qudit").RequiresN);
        Assert.True(InspectCommand.Catalog.Single(e => e.Name == "fourmode").RequiresN);
    }
}
