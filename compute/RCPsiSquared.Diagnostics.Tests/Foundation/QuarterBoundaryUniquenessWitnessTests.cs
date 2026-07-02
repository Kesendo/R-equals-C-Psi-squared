using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness for UNIQUENESS_PROOF. The gate is two-sided so it cannot pass
/// trivially: α=2 must be state-independent (spread 0) AND equal ¼, while every OTHER swept Rényi
/// order must be genuinely state-dependent (spread &gt; 0); and the α=2 discriminant D=1−4CΨ must have
/// its single zero at ¼.</summary>
public class QuarterBoundaryUniquenessWitnessTests
{
    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new QuarterBoundaryUniquenessWitness(psiLow: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new QuarterBoundaryUniquenessWitness(psiHigh: 1.0));
        Assert.Throws<ArgumentException>(() => new QuarterBoundaryUniquenessWitness(psiLow: 0.5, psiHigh: 0.5));
    }

    [Fact]
    public void AlphaTwo_IsStateIndependent_AndEqualsQuarter()   // the forcing (Step 6)
    {
        var w = new QuarterBoundaryUniquenessWitness(psiLow: 0.3, psiHigh: 0.7);
        var two = w.AtTwo;
        Assert.True(two.Spread < 1e-12, $"α=2 threshold must be state-independent, spread {two.Spread:E3}");
        Assert.Equal(0.25, two.ThresholdLow, 12);
        Assert.Equal(0.25, two.ThresholdHigh, 12);
    }

    [Fact]
    public void AlphaTwo_IsTheUniqueStateIndependentOrder()   // uniqueness, two-sided
    {
        var w = new QuarterBoundaryUniquenessWitness(psiLow: 0.3, psiHigh: 0.7);
        var stateIndependent = w.Rows.Where(r => r.Spread < 1e-12).ToList();
        Assert.Single(stateIndependent);
        Assert.Equal(2.0, stateIndependent[0].Alpha);
        // every other swept order is genuinely state-dependent (the other side of the gate)
        foreach (var r in w.Rows.Where(r => r.Alpha != 2.0))
            Assert.True(r.Spread > 1e-3, $"α={r.Alpha} must be state-dependent, spread {r.Spread:E3}");
    }

    [Fact]
    public void Discriminant_HasItsUniqueZeroAtQuarter()
    {
        Assert.True(QuarterBoundaryUniquenessWitness.Discriminant(0.20) > 0, "D>0 below 1/4 (two real fixed points)");
        Assert.Equal(0.0, QuarterBoundaryUniquenessWitness.Discriminant(0.25), 12);          // the fold tangency
        Assert.True(QuarterBoundaryUniquenessWitness.Discriminant(0.30) < 0, "D<0 above 1/4 (no real fixed point)");
    }

    [Fact]
    public void Witness_SurfacesUniquenessAndDiscriminantChildren()
    {
        var labels = ((IInspectable)new QuarterBoundaryUniquenessWitness()).Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("unique state-independent"));
        Assert.Contains(labels, l => l.Contains("discriminant boundary"));
    }
}
