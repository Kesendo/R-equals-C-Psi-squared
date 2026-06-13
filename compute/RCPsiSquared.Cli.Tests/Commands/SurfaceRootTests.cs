using System.Linq;
using RCPsiSquared.Cli;
using RCPsiSquared.Cli.Commands;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Cli.Tests.Commands;

public class SurfaceRootTests
{
    [Fact]
    public void Catalog_HasSurfaceRoot_ThatBuildsTheN5Witness()
    {
        var entry = InspectCommand.Catalog.SingleOrDefault(e => e.Name == "surface");
        Assert.NotNull(entry);

        var ctx = new InspectRootContext(
            Parser: new ArgParser(new[] { "--root", "surface" }),
            N: 5, WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var witness = entry!.Factory(ctx);
        Assert.IsType<BirthCanalSurfaceWitness>(witness);
        Assert.Equal(5, ((BirthCanalSurfaceWitness)witness).N);     // the witness is N=5 only
    }

    [Fact]
    public void SurfaceRoot_HonorsGridOverride()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "surface");
        var ctx = new InspectRootContext(
            Parser: new ArgParser(new[] { "--root", "surface", "--grid", "5" }),
            N: 5, WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var witness = (BirthCanalSurfaceWitness)entry.Factory(ctx);
        Assert.Equal(5, witness.GridK);
    }
}
