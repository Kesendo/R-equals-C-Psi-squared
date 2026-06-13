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

    [Fact]
    public void FlatGammaBlindness_RateIsTwoGamma_AtEveryN()
    {
        // analytic: at uniform gamma, L = -iQ h - 2 gamma I, -iQh anti-Hermitian -> Re=-2gamma all modes.
        foreach (int n in new[] { 5, 6, 8 })
        {
            var uni = System.Linq.Enumerable.Repeat(1.0, n).ToArray();   // gamma_l = 1
            Assert.Equal(2.0, SectorReductionWitness.VacBlockSlowest(n, 1.5, uni, TopologyKind.Chain), 9);
            Assert.Equal(2.0, SectorReductionWitness.VacBlockSlowest(n, 1000.0, uni, TopologyKind.Chain), 9);
        }
    }

    [Fact]
    public void TheVacReductionNode_RendersAndGoesPastN5()
    {
        var w8 = new SectorReductionWitness(n: 8);                       // past the dense N=6 ceiling
        var node = System.Linq.Enumerable.First(w8.Children);
        Assert.Equal("the |1-exc><vac| reduction", node.DisplayName);
        Assert.Contains("N=8", node.Summary + " " + string.Join(" ",
            System.Linq.Enumerable.Select(System.Linq.Enumerable.ToList(node.Children), c => c.Summary)));
    }
}
