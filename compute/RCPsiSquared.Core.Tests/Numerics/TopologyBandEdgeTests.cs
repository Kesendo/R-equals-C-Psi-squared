using System;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>The law: the XY single-excitation band edge = J × the hopping graph's adjacency spectral
/// radius. Closed forms chain 2cos(π/(N+1)), star √(N−1), ring 2. (The Re=−2γ floor is the Absorption
/// Theorem; this primitive is only the Im/frequency side.)</summary>
public class TopologyBandEdgeTests
{
    private static System.Collections.Generic.IReadOnlyList<Bond> Bonds(int n, TopologyKind topo) =>
        new ChainSystem(N: n, J: 1.0, GammaZero: 0.05, HType: HamiltonianType.XY, Topology: topo).Bonds;

    [Theory]
    [InlineData(TopologyKind.Chain, 3)]
    [InlineData(TopologyKind.Chain, 4)]
    [InlineData(TopologyKind.Chain, 5)]
    [InlineData(TopologyKind.Chain, 8)]
    [InlineData(TopologyKind.Star, 3)]
    [InlineData(TopologyKind.Star, 5)]
    [InlineData(TopologyKind.Star, 8)]
    [InlineData(TopologyKind.Ring, 3)]
    [InlineData(TopologyKind.Ring, 6)]
    public void SpectralRadius_MatchesClosedForm(TopologyKind topo, int n)
    {
        double expected = topo switch
        {
            TopologyKind.Chain => 2.0 * Math.Cos(Math.PI / (n + 1)),
            TopologyKind.Star => Math.Sqrt(n - 1),
            TopologyKind.Ring => 2.0,
            _ => throw new ArgumentOutOfRangeException(nameof(topo)),
        };
        Assert.Equal(expected, TopologyBandEdge.SpectralRadius(n, Bonds(n, topo)), 9);
    }

    [Fact]
    public void BandEdge_IsJTimesSpectralRadius()
    {
        var bonds = Bonds(4, TopologyKind.Star);   // ρ = √3
        Assert.Equal(2.5 * Math.Sqrt(3.0), TopologyBandEdge.BandEdge(4, bonds, j: 2.5), 9);
    }
}
