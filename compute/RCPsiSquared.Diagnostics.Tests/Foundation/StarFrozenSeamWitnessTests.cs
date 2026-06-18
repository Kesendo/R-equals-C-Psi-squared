using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gate-first (the C# twin of simulations/star_frozen_seam.py): the star's longest-lived coherence
/// never un-freezes for N≥5, and that threshold IS the structural ceiling g2 = 4/(N−1) ≤ 1. The exact
/// commutant ceiling + the |Im|(Q) dynamic on the full 4^N Liouvillian (N≤5).</summary>
public class StarFrozenSeamWitnessTests
{
    private static readonly double[] Sweep = { 1.0, 2.0, 4.0, 8.0, 16.0, 32.0 };

    /// <summary>The star's (1,1) commutant ceiling is g2 = 4/(N−1) (the leaf-manifold structural ceiling,
    /// recomputed live from the commutant rep). The threshold: ≤ 1 (frozen seam) iff N ≥ 5.</summary>
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void Ceiling_IsFourOverNMinusOne(int n)
    {
        Assert.Equal(4.0 / (n - 1), StarFrozenSeamWitness.Ceiling(n), 5);
    }

    /// <summary>N=5 (g2 = 1, marginal, ≤ 1): the survivor is the frozen commutant coherence at EVERY Q —
    /// it never acquires a frequency. The star's frozen seam.</summary>
    [Fact]
    public void Star5_SurvivorFrozenAtEveryQ()
    {
        double mx = StarFrozenSeamWitness.MaxImOverSweep(5, Sweep);
        Assert.True(mx < 1e-6, $"star N=5 survivor should be frozen at every Q; max |Im| = {mx:e}");
    }

    /// <summary>N=4 (g2 = 4/3 > 1): the commutant is brighter than the −2γ floor, an oscillating band-edge
    /// mode is the slowest, and the star UN-freezes — the known (2,2)/K₄ outlier (the gate caught it).</summary>
    [Fact]
    public void Star4_SurvivorUnFreezes_TheKnownOutlier()
    {
        double mx = StarFrozenSeamWitness.MaxImOverSweep(4, Sweep);
        Assert.True(mx > 1e-2, $"star N=4 survivor should un-freeze (g2=4/3>1, the outlier); max |Im| = {mx:e}");
    }
}
