using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of the sideways spin ladder renders the shared chain walk
/// (<see cref="SidewaysSpinLadderChain"/>, the same construction the gate asserts on). These tests pin
/// the witness's shape and its N=5 numbers; the physics gates live in
/// <see cref="SidewaysSpinLadderGateTests"/>.</summary>
public class SidewaysSpinLadderWitnessTests
{
    [Fact]
    [Trait("Category", "SIDEWAYS")]
    public void N5_Summary_ReportsExactControlAndCgMatch()
    {
        var w = new SidewaysSpinLadderWitness(5);
        Assert.Contains("exactly 0.0", w.Summary);
        Assert.Contains("= CG coefficients", w.Summary);
        Assert.Contains("= the two chain interiors", w.Summary);
    }

    [Fact]
    [Trait("Category", "SIDEWAYS")]
    public void N5_Children_CarryControlShapeTwoChainsAndTheOpenNode()
    {
        var w = new SidewaysSpinLadderWitness(5);
        var children = w.Children.ToList();
        Assert.Equal(5, children.Count);   // CONTROL + SHAPE/COUNT + two chains + open
        Assert.Contains(children, c => c.DisplayName.StartsWith("CONTROL", StringComparison.Ordinal));
        Assert.Contains(children, c => c.DisplayName.StartsWith("SHAPE + COUNT", StringComparison.Ordinal));
        Assert.Contains(children, c => c.DisplayName.Contains("p+q̃ = 4"));
        Assert.Contains(children, c => c.DisplayName.Contains("p+q̃ = 6"));
        Assert.Contains(children, c => c.DisplayName.StartsWith("N=9 confirmed", StringComparison.Ordinal));
    }

    [Fact]
    [Trait("Category", "SIDEWAYS")]
    public void InvalidN_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SidewaysSpinLadderWitness(6));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SidewaysSpinLadderWitness(9));
    }
}
