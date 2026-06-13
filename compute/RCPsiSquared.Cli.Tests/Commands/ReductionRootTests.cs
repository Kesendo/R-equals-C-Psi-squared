using System.Linq;
using RCPsiSquared.Cli;
using RCPsiSquared.Cli.Commands;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Cli.Tests.Commands;

public class ReductionRootTests
{
    [Fact]
    public void Catalog_HasReductionRoot_ThatBuildsTheWitness()
    {
        var entry = InspectCommand.Catalog.SingleOrDefault(e => e.Name == "reduction");
        Assert.NotNull(entry);
        var ctx = new InspectRootContext(
            Parser: new ArgParser(new[] { "--root", "reduction" }),
            N: 5, WithQSweep: false, WithMeasured: false, QGridPoints: null);
        Assert.IsType<SectorReductionWitness>(entry!.Factory(ctx));
    }
}
