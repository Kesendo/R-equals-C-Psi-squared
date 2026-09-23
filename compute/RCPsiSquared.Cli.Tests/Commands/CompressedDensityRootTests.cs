using RCPsiSquared.Cli.Commands;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Cli.Tests.Commands;

public class CompressedDensityRootTests
{
    [Fact]
    public void CatalogExposesTheFixedN11LiveWitness()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "compresseddensity");
        Assert.False(entry.RequiresN);
        Assert.False(entry.HonorsOptionalN);
        var context = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 5,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var witness = Assert.IsType<CompressedDensityN11Witness>(entry.Factory(context));
        Assert.Contains("N=11", witness.DisplayName);
        Assert.Contains("identity fails", witness.Summary);
    }
}
