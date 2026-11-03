using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.Tests.Resonance;

public class EnvelopeRefinementTests
{
    /// <summary>f(t) = 1 − (t − t0)², curvature f'' = −2. Golden section compares values, so it
    /// resolves t only to where two values differ: |Δt| ≈ √(2ε/|f''|) = √ε ≈ 1.5e-8 around t0 = 5
    /// (the relative ε of a value near 1). The gate is twice that; the value itself is exact to a
    /// few ε. t0 sits between grid points, so the grid alone would miss it by up to dt/2 = 0.125.</summary>
    [Fact]
    public void Refine_FindsTheContinuousMaximum_BetweenGridPoints()
    {
        const double t0 = 5.1234567;
        var grid = ResonanceScan.DefaultTGrid(0.05);
        Func<double, double> f = t => 1.0 - (t - t0) * (t - t0);
        int iWin = 0;
        for (int i = 1; i < grid.Length; i++) if (f(grid[i]) > f(grid[iWin])) iWin = i;

        var (k, tStar, edge) = EnvelopeRefinement.Refine(grid, iWin, f(grid[iWin]), f);

        Assert.False(edge);
        Assert.True(Math.Abs(tStar - t0) <= 2 * Math.Sqrt(Math.Pow(2, -52)), $"t* = {tStar:R}");
        Assert.True(Math.Abs(k - 1.0) <= 4 * Math.Pow(2, -52), $"k = {k:R}");
        Assert.True(k > f(grid[iWin]));
    }

    [Fact]
    public void Refine_ReportsAWinnerOnTheWindowEdge()
    {
        var grid = ResonanceScan.DefaultTGrid(0.05);
        Func<double, double> rising = t => t;
        var (k, tStar, edge) = EnvelopeRefinement.Refine(grid, grid.Length - 1, grid[^1], rising);
        Assert.True(edge);
        Assert.True(tStar <= grid[^1] && k <= grid[^1]);
    }

    /// <summary>On the full block every maximiser t* lies inside the seed window and no grid winner
    /// sits on its edge (c = 2, N = 5, Q across the resonance).</summary>
    [Fact]
    public void Envelope_NoCellSitsOnTheWindowEdge_AcrossTheResonance()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var scan = new ResonanceScan(block);
        var q = ResonanceScan.LinearQGrid(0.5, 3.5, 13);
        var curve = scan.ComputeKCurve(q);
        for (int b = 0; b < curve.NumBonds; b++)
            for (int iQ = 0; iQ < q.Length; iQ++)
                Assert.False(curve.TStarAtWindowEdge[b, iQ], $"edge at b={b}, Q={q[iQ]}");
    }
}
