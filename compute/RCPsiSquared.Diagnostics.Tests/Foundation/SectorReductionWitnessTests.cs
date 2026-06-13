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

    [Fact]
    public void Junction_AtN6DeepEdge_GlobalSlowestCrossesToTheO2DensityMode()
    {
        // birth_canal_n6_mode_crossing.py in C#: at N=6 deep-edge, Q=1.5 the (2,2) density block's
        // slowest non-kernel rate is LESS than the (0,1) block's -> the {0,2} mode wins (the crossing).
        var deep = new[] { 0.25, 1.375, 1.375, 1.375, 1.375, 0.25 };
        double vac10 = SectorReductionWitness.VacBlockSlowest(6, 1.5, deep, TopologyKind.Chain);    // ~1.471
        double dens22 = SectorReductionWitness.SectorSlowest(6, 1.5, deep, 2, 2, TopologyKind.Chain); // ~1.120
        Assert.True(dens22 < vac10, $"(2,2) {dens22} should be slower than (0,1) {vac10} at N=6 deep-edge Q=1.5");
        // and at high Q the (0,1) mode is back to (or below) the density mode:
        double vac10Hi = SectorReductionWitness.VacBlockSlowest(6, 1000.0, deep, TopologyKind.Chain);
        double dens22Hi = SectorReductionWitness.SectorSlowest(6, 1000.0, deep, 2, 2, TopologyKind.Chain);
        Assert.True(vac10Hi <= dens22Hi + 1e-9, "at Q=1000 the (0,1) odd mode is the global slowest again");
    }

    [Fact]
    public void JunctionNode_RendersAtN6_WithO2WeightedHistogram_NoCrash()
    {
        // Exercises TheJunctionNode -> DensityMode -> NDiffHistogram (the path the rate tests skip).
        var children = System.Linq.Enumerable.ToList(new SectorReductionWitness(6).Children);
        var junction = System.Linq.Enumerable.Single(children, c => c.DisplayName == "the {0,2} junction");
        // at N=6 deep-edge the {0,2} density mode wins, and the n_diff histogram is {0,2}-weighted:
        Assert.Contains("WINS (the crossing)", junction.Summary);
        Assert.Contains("n_diff hist {0:", junction.Summary);
    }
}
