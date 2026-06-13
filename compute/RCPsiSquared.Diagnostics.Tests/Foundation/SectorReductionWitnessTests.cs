using System;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class SectorReductionWitnessTests
{
    private static readonly SectorReductionWitness W = new(n: 5);

    [Theory]
    [InlineData(1.0, 1.0, 1.0, 1.0, 1.0)]            // uniform
    [InlineData(0.25, 0.75, 3.0, 0.75, 0.25)]        // peaked-V
    [InlineData(0.25, 1.5, 1.5, 1.5, 0.25)]          // flat-bulk-edge (canal)
    public void VacBlock_SlowestRate_EqualsPostEpFlowField_AtN5(
        double g0, double g1, double g2, double g3, double g4)
    {
        var profile = new[] { g0, g1, g2, g3, g4 };
        foreach (double q in new[] { 1.5, 1000.0 })
        {
            double block = SectorReductionWitness.VacBlockSlowest(5, q, profile, TopologyKind.Chain);
            var field = new PostEpFlowField(5, new[] { 1.5, 1000.0 }, new[] { 0.0, 1.0 }, profile);
            double full = field.ReadAssembly(q).SlowestRate;
            Assert.Equal(full, block, 9);   // the (0,1) block reproduces the full global slowest at N=5
        }
    }

    [Fact]
    public void VacBlock_CanalAnchor_MatchesThePinnedRates()
    {
        var canal = new[] { 0.25, 1.5, 1.5, 1.5, 0.25 };
        Assert.Equal(1.2482918643729715, SectorReductionWitness.VacBlockSlowest(5, 1.5, canal, TopologyKind.Chain), 6);
        Assert.Equal(4.0 / 3.0, SectorReductionWitness.VacBlockSlowest(5, 1000.0, canal, TopologyKind.Chain), 5);
    }
}
