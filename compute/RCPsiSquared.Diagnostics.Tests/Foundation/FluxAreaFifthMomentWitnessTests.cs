using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using B = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of <see cref="FluxAreaFifthMomentClaim"/> (F137, the fifth centred trace under cooling
/// and heating): the generator's Tr((2L + c)⁵) in exact Gaussian integers against the closed form. Every comparison
/// is exact.</summary>
public class FluxAreaFifthMomentWitnessTests
{
    private static readonly FluxAreaFifthMomentWitness W = new();

    [Fact]
    public void EveryRow_Holds()
    {
        foreach (var x in W.Readings)
            Assert.True(x.Holds, $"{x.Row.Label}: generator {x.Generator}, closed form {x.Law}");
    }

    [Fact]
    public void RecordedValues_AreReproduced()
    {
        BigInteger scale = 32 * 120 * 120;
        Assert.Equal(-656 * scale, W.Readings[0].Generator);
        Assert.Equal(656 * scale, W.Readings[1].Generator);
        Assert.Equal(-3936 * scale, W.Readings[2].Generator);
    }

    [Fact]
    public void Rows_AreNotVacuous()
    {
        Assert.True(W.NotVacuous, W.Summary);
        // the silent rows are silent in the generator, not only in the closed form
        Assert.All(W.Readings.Where(x => !x.Row.IsControl && x.Law.IsZero), x => Assert.True(x.Generator.IsZero, x.Row.Label));
    }

    /// <summary>The rows with several triangles can tell a wrong form from the right one: half the constant, and
    /// one triangle walked against its determinant's columns, each miss the generator there (corollaries of the
    /// generator check on those rows, kept as the record that such rows exist).</summary>
    [Fact]
    public void WrongForms_MissTheRowsWithSeveralTriangles()
    {
        var rows = W.Readings.Where(x => !x.Row.IsControl && x.Triangles.Count(t => !t.Product.IsZero) >= 2).ToArray();
        Assert.True(rows.Length >= 2);
        foreach (var x in rows)
        {
            var scale = 1920 * BigInteger.Pow(4, x.Row.N);
            var flippedFirst = x.Triangles.Skip(1).Aggregate(-x.Triangles[0].Product, (acc, t) => acc + t.Product);
            Assert.NotEqual(x.Generator, scale * flippedFirst);
            Assert.NotEqual(x.Generator, 960 * BigInteger.Pow(4, x.Row.N) * x.Triangles.Aggregate(BigInteger.Zero, (acc, t) => acc + t.Product));
        }
    }

    /// <summary>The pendant's own rate point does not enter: moving it leaves the generator's fifth trace
    /// unchanged (the closed form does not read it at all).</summary>
    [Fact]
    public void PendantRatePoint_DoesNotEnter()
    {
        var bonds = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 0, 0), new B(2, 0, 0, 1, 0), new B(3, 0, 1, 0, 1) };
        var a = FluxAreaFifthMomentWitness.Read(new("pendant at (7, 2)", 4, bonds, new long[] { 2, 1, 3, 7 }, new long[] { 1, 3, 1, 2 }));
        var b = FluxAreaFifthMomentWitness.Read(new("pendant at (1, 9)", 4, bonds, new long[] { 2, 1, 3, 1 }, new long[] { 1, 3, 1, 9 }));
        Assert.False(a.Generator.IsZero);
        Assert.Equal(a.Generator, b.Generator);
        Assert.Equal(a.Law, b.Law);
    }

    /// <summary>The bond's direction: (a, b, J, D) and (b, a, J, −D) are the same bond.</summary>
    [Fact]
    public void ReversedBond_IsTheSameBond()
    {
        var cool = new long[] { 2, 1, 3 };
        var heat = new long[] { 1, 3, 1 };
        var one = FluxAreaFifthMomentClaim.FifthCentredTrace(3, new[] { new B(0, 1, 1, 2, 0), new B(1, 2, 2, -1, 0), new B(2, 0, 1, 3, 0) }, cool, heat);
        var two = FluxAreaFifthMomentClaim.FifthCentredTrace(3, new[] { new B(1, 0, 1, -2, 0), new B(1, 2, 2, -1, 0), new B(0, 2, 1, -3, 0) }, cool, heat);
        Assert.False(one.IsZero);
        Assert.Equal(one, two);
        // and the generator agrees with both listings
        var read = FluxAreaFifthMomentWitness.Read(new("reversed listing", 3, new[] { new B(1, 0, 1, -2, 0), new B(1, 2, 2, -1, 0), new B(0, 2, 1, -3, 0) }, cool, heat));
        Assert.Equal(one, read.Generator);
    }

    [Fact]
    public void BadInput_Throws()
    {
        Assert.Throws<System.ArgumentException>(() => FluxAreaFifthMomentClaim.FifthCentredTrace(3, new[] { new B(0, 3, 1, 0, 0) }, new long[] { 1, 1, 1 }, new long[] { 0, 0, 0 }));
        Assert.Throws<System.ArgumentException>(() => FluxAreaFifthMomentClaim.FifthCentredTrace(3, new[] { new B(1, 1, 1, 0, 0) }, new long[] { 1, 1, 1 }, new long[] { 0, 0, 0 }));
        Assert.Throws<System.ArgumentException>(() => FluxAreaFifthMomentClaim.FifthCentredTrace(3, new[] { new B(0, 1, 1, 0, 0) }, new long[] { 1, 1 }, new long[] { 0, 0, 0 }));
    }
}
