using System;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gate-first: the live commutant computation reproduces the structural-ceiling closed forms
/// (the C# twin of simulations/topology_ceiling_rep_derivation.py). Small sector matrices, fast.</summary>
public class StructuralCeilingWitnessTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void Complete_OneOneSector_Is4OverN(int n)
    {
        double? g = StructuralCeilingWitness.CommutantDarkest("complete", n, 1, 1);
        Assert.True(g.HasValue);
        Assert.Equal(4.0 / n, g!.Value, 9);
    }

    [Theory]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Star_OneOneSector_Is4OverNMinus1(int n)
    {
        double? g = StructuralCeilingWitness.CommutantDarkest("star", n, 1, 1);
        Assert.True(g.HasValue);
        Assert.Equal(4.0 / (n - 1), g!.Value, 9);
    }

    [Fact]
    public void N4Outlier_LivesInTheTwoTwoHalfFillingSector()
    {
        // the (1,1) ladder hits 1.0 at N=4 (= 4/4 = the band edge), so it makes no ceiling there
        double? k4Ladder = StructuralCeilingWitness.CommutantDarkest("complete", 4, 1, 1);
        Assert.True(k4Ladder.HasValue);
        Assert.Equal(1.0, k4Ladder!.Value, 9);

        // the K_4 ceiling is the (2,2) half-filling sector = 2 − 2/√3, below the floor
        double? k4Ceiling = StructuralCeilingWitness.CommutantDarkest("complete", 4, 2, 2);
        Assert.True(k4Ceiling.HasValue);
        Assert.Equal(2.0 - 2.0 / Math.Sqrt(3.0), k4Ceiling!.Value, 7);
        Assert.True(k4Ceiling.Value < 1.0);

        // ring-4 is special in the SAME (2,2) sector, but co-occupies the floor (= 1.0)
        double? ring4Ceiling = StructuralCeilingWitness.CommutantDarkest("ring", 4, 2, 2);
        Assert.True(ring4Ceiling.HasValue);
        Assert.Equal(1.0, ring4Ceiling!.Value, 9);
    }

    [Fact]
    public void Chain_NeverCeilings()
    {
        // the chain (no adjacency degeneracy) keeps its (1,1) ladder above 1: the band edge protects
        for (int n = 4; n <= 6; n++)
        {
            double? g = StructuralCeilingWitness.CommutantDarkest("chain", n, 1, 1);
            if (g.HasValue) Assert.True(g.Value > 1.0 - 1e-9, $"chain N={n} should not ceiling, got {g.Value}");
        }
    }
}
