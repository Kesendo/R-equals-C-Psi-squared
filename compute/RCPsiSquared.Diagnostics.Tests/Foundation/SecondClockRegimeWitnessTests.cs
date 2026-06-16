using System;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gate-first: the live N=4 full-Liouvillian classification reproduces the gate-verified 2D regime
/// map (all three regimes present at N=4), and the high-Q ceiling knob matches the commutant closed forms.
/// The C# twin of simulations/second_clock_regime_axis.py (N=4 slice). N=4 = 256² Liouvillian, fast.</summary>
public class SecondClockRegimeWitnessTests
{
    [Theory]
    [InlineData("chain", "EP")]
    [InlineData("chain_disordered", "EP")]
    [InlineData("ring", "GRADUAL")]
    [InlineData("star", "EP")]
    [InlineData("complete", "CEILING")]
    public void Classify_N4_MatchesTheGateVerifiedMap(string topo, string expected)
    {
        Assert.Equal(expected, SecondClockRegimeWitness.Classify(topo, 4));
    }

    [Fact]
    public void N4_HasAllThreeRegimes()
    {
        // the point of N=4: EP, GRADUAL, and CEILING are all present, so one slice gates the whole map
        var regimes = new[] { "chain", "chain_disordered", "ring", "star", "complete" }
            .Select(t => SecondClockRegimeWitness.Classify(t, 4));
        Assert.Contains("EP", regimes);
        Assert.Contains("GRADUAL", regimes);
        Assert.Contains("CEILING", regimes);
    }

    [Fact]
    public void HighQCeiling_StitchesToTheCommutantClosedForms()
    {
        // complete-4 ceilings at the (2,2) outlier 2 − 2/√3 (below the floor)
        Assert.Equal(2.0 - 2.0 / Math.Sqrt(3.0), SecondClockRegimeWitness.HighQCeiling("complete", 4), 6);
        // complete-5 = 4/5, complete-6 = 2/3 (the 4/N ladder)
        Assert.Equal(4.0 / 5.0, SecondClockRegimeWitness.HighQCeiling("complete", 5), 6);
        Assert.Equal(2.0 / 3.0, SecondClockRegimeWitness.HighQCeiling("complete", 6), 6);
        // the star graduates: reaches the floor (= 1) at N≤5, first ceilings at N=6 (4/5)
        Assert.Equal(1.0, SecondClockRegimeWitness.HighQCeiling("star", 5), 6);
        Assert.Equal(4.0 / 5.0, SecondClockRegimeWitness.HighQCeiling("star", 6), 6);
        // the chain never ceilings
        Assert.Equal(1.0, SecondClockRegimeWitness.HighQCeiling("chain", 6), 9);
    }

    [Fact]
    public void AdjacencyBand_ReadsDegeneracyAndDispersion()
    {
        // complete K_6: one (N−1)-fold degenerate −1 band + a single +5 level → m = 5, 2 distinct
        var (_, mC, distinctC, _) = SecondClockRegimeWitness.AdjacencyBand("complete", 6);
        Assert.Equal(5, mC);
        Assert.Equal(2, distinctC);
        // chain N=4: a dispersive cosine band, all distinct
        var (_, mChain, distinctChain, _) = SecondClockRegimeWitness.AdjacencyBand("chain", 4);
        Assert.Equal(1, mChain);
        Assert.Equal(4, distinctChain);
    }
}
