using System.Linq;
using RCPsiSquared.Cli;
using RCPsiSquared.Cli.Commands;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Cli.Tests.Commands;

public class TrichotomyRootTests
{
    [Fact]
    public void Catalog_HasTrichotomyRoot() =>
        Assert.Contains(InspectCommand.Catalog, e => e.Name == "trichotomy");

    [Fact]
    public void Catalog_TrichotomyRoot_NFree_HonorsOptionalN_DefaultDepthTwo()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "trichotomy");
        Assert.False(entry.RequiresN);
        Assert.True(entry.HonorsOptionalN);
        Assert.Equal(2, entry.DefaultDepth);
    }

    [Fact]
    public void Catalog_TrichotomyFactory_BuildsTheLiveWitness_DefaultN6Q15()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "trichotomy");
        var ctx = new InspectRootContext(
            Parser: new ArgParser(System.Array.Empty<string>()),
            N: 1, WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = entry.Factory(ctx);
        var witness = Assert.IsType<TrichotomyWitness>(root);
        Assert.Equal(6, witness.N);
        Assert.Equal(1.5, witness.Q);
    }
}
