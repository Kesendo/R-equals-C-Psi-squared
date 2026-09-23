using RCPsiSquared.Cli;
using RCPsiSquared.Cli.Commands;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Cli.Tests.Commands;

public class NodePairInspectRootTests
{
    [Fact]
    public void Catalog_ExposesFixedExactNodePairWitnessWithScope()
    {
        var entries = InspectCommand.Catalog.Where(entry => entry.Name == "nodepair").ToArray();
        var entry = Assert.Single(entries);

        Assert.False(entry.RequiresN);
        Assert.False(entry.HonorsOptionalN);
        Assert.Contains("PROOF_NODE_PAIR_RESOLVENT", entry.Description);
        Assert.Contains("uniform-centre iff", entry.Description);
        Assert.Contains("off-centre", entry.Description);
        Assert.Contains("baseline blind energy", entry.Description);
        Assert.Contains("new blind energies", entry.Description);

        var context = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 99,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var witness = Assert.IsType<NodePairResolventWitness>(entry.Factory(context));
        Assert.True(witness.Reading.FactorizationMatches);
    }
}
