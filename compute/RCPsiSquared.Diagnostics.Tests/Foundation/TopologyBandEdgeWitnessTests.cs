using System;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

/// <summary>The gate-first map: the band edge is the gap mode (ω_mem == J·ρ) only in a topology-specific
/// regime. Ground truth: simulations/topology_band_edge_review.py (full Q-sweep to Q=1000).</summary>
namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class TopologyBandEdgeWitnessTests
{
    private static readonly TopologyBandEdgeWitness W = new(j: 1.0, gamma: 0.05);

    // The DIMENSIONLESS reading ω_mem / band-edge: 1 = protected, 0 = broken, √2 = ring-N4 high-Q.
    private static double R(TopologyKind t, int n, double q) => W.OmegaOverBand(t, n, q);
    private static double Band(TopologyKind t, int n) => W.BandEdge(t, n);  // J=1 ⟹ = ρ

    [Fact]
    public void Law_BandEdgeIsJTimesSpectralRadius()
    {
        Assert.Equal(2.0 * Math.Cos(Math.PI / 5), Band(TopologyKind.Chain, 4), 9);  // φ
        Assert.Equal(Math.Sqrt(3.0), Band(TopologyKind.Star, 4), 9);                 // √3
        Assert.Equal(2.0, Band(TopologyKind.Ring, 4), 9);
    }

    [Fact]
    public void Chain_Protected_AboveTheHorizon()
    {
        foreach (var n in new[] { 3, 4, 5 })
            Assert.Equal(1.0, R(TopologyKind.Chain, n, 20.0), 4);
    }

    [Fact]
    public void Star_HorizonVsCeiling_TheGateFirstProbe()
    {
        Assert.Equal(1.0, R(TopologyKind.Star, 3, 20.0), 4);
        Assert.Equal(1.0, R(TopologyKind.Star, 4, 20.0), 4);
        // N=5 is a Q-HORIZON: broken at Q=20, protects at Q=1000.
        Assert.Equal(0.0, R(TopologyKind.Star, 5, 20.0), 4);
        Assert.Equal(1.0, R(TopologyKind.Star, 5, 1000.0), 3);
        // N=6 is the STRUCTURAL ceiling: ω_mem=0 at BOTH Q=20 and Q=1000.
        Assert.Equal(0.0, R(TopologyKind.Star, 6, 20.0), 4);
        Assert.Equal(0.0, R(TopologyKind.Star, 6, 1000.0), 4);
    }

    [Fact]
    public void Ring_N4_CoOccupiedFloor()
    {
        Assert.Equal(1.0, R(TopologyKind.Ring, 3, 20.0), 4);
        Assert.Equal(1.0, R(TopologyKind.Ring, 5, 20.0), 4);
        Assert.Equal(0.0, R(TopologyKind.Ring, 4, 20.0), 4);
        Assert.Equal(Math.Sqrt(2.0), R(TopologyKind.Ring, 4, 1000.0), 3);   // 2√2·J ÷ 2J = √2
    }

    [Fact]
    public void Witness_Renders_WithoutThrowing()
    {
        var kids = ((RCPsiSquared.Core.Inspection.IInspectable)W).Children;
        foreach (var k in kids) _ = k.Summary;   // touch every lens; must not throw
        Assert.NotEmpty(kids);
    }
}
