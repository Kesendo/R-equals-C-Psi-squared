using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using B = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of <see cref="TransportSeventhMomentClaim"/> (Theorem B of
/// PROOF_TRANSPORT_FORMULA_LOSS_AND_GAIN, F137): the part of degree two in the rates of Tr((2L + c)⁷), read from the
/// generator in exact Gaussian integers as the words with two rate letters, against the closed form. Every
/// comparison is exact.</summary>
public class TransportSeventhMomentWitnessTests
{
    private static readonly TransportSeventhMomentWitness W = new();

    [Fact]
    public void EveryRow_Holds()
    {
        foreach (var x in W.Readings)
            Assert.True(x.Holds, $"{x.Row.Label}: degree two {x.DegreeTwo}, closed form {x.Law}, trace {x.Trace}");
    }

    [Fact]
    public void Rows_AreNotVacuous()
    {
        Assert.True(W.NotVacuous, W.Summary);
        // the silent rows are silent in the generator, not only in the closed form
        Assert.All(W.Readings.Where(x => !x.Row.IsControl && x.Law.IsZero), x => Assert.True(x.DegreeTwo.IsZero, x.Row.Label));
    }

    /// <summary>The recorded values, as numbers written here and not taken from the rows: gate M11's two breaks and
    /// the fifth trace's gate's ring are all of degree two, the temperature gradient and the triangle with unequal
    /// weights all of degree four.</summary>
    [Fact]
    public void RecordedValues_AreReproduced()
    {
        BigInteger s72 = 128 * 72 * 72;
        var byLabel = W.Readings.ToDictionary(x => x.Row.Label);
        TransportSeventhMomentWitness.Reading Row(string start) => byLabel.Single(p => p.Key.StartsWith(start)).Value;

        var pendant = Row("gate M11's triangle");
        Assert.Equal(-44800 * s72, pendant.DegreeTwo);
        Assert.Equal(pendant.DegreeTwo, pendant.Trace);

        // −20380416⁄25 at H and the rates as recorded; here H times 10 (degree five) and the rates times 840 (degree two)
        var ring = Row("gate M11's ring of four");
        Assert.Equal(BigInteger.Parse("-20380416") * 4000 * 840 * 840 * 128, ring.DegreeTwo);
        Assert.Equal(ring.DegreeTwo, ring.Trace);

        var fifthRing = Row("the ring of four of the fifth trace's gate");
        Assert.Equal(-68239360 * s72, fifthRing.DegreeTwo);
        Assert.Equal(fifthRing.DegreeTwo, fifthRing.Trace);

        var bare = Row("the same ring without its ZZ");
        Assert.True(bare.DegreeTwo.IsZero && bare.Trace.IsZero);

        var gradient = Row("a temperature gradient");
        Assert.True(gradient.DegreeTwo.IsZero);
        Assert.Equal((BigInteger)430080 * 128, gradient.Trace);

        var skew = Row("hoppings 1, 1 + i, 2i");
        Assert.True(skew.DegreeTwo.IsZero);
        Assert.Equal((BigInteger)322560 * 128, skew.Trace);
    }

    /// <summary>β_lk = −β_kl on every row's graph, and not all of them zero.</summary>
    [Fact]
    public void Beta_IsAntisymmetric()
    {
        int nonzero = 0;
        foreach (var row in TransportSeventhMomentWitness.Rows.Where(r => !r.IsControl))
            for (int l = 0; l < row.N; l++)
                for (int k = l + 1; k < row.N; k++)
                {
                    var lk = TransportSeventhMomentClaim.ThriceBeta(row.N, row.Bonds, l, k);
                    Assert.Equal(-lk, TransportSeventhMomentClaim.ThriceBeta(row.N, row.Bonds, k, l));
                    if (!lk.IsZero) nonzero++;
                }
        Assert.True(nonzero >= 20, $"{nonzero}");
    }

    /// <summary>On a triangle alone with three equal |t| and one Δ the three weights are equal, and the part of
    /// degree two is the loop times the determinant of the rate triangle, with the weight
    /// 30·3 − 10·4 − 5·2 − 2 + 10·Δ² = 38 + 10Δ² (S = 3, d = 2, |t_lk|² = 1; Q = 10Δ²).</summary>
    [Fact]
    public void EqualTriangle_IsFluxTimesArea()
    {
        var reading = W.Readings.Single(x => x.Row.Label.StartsWith("the unit triangle at flux π/2 with one Δ"));
        var tri = FluxAreaFifthMomentClaim.Triangles(3, reading.Row.Bonds, reading.Row.Cool, reading.Row.Heat).Single();
        // Σ_{l≠k} s_l·a_k·β_lk with β_lk = ±w·Im(loop) around the triangle is −2·w·Im(loop)·det
        BigInteger w = 38 + 10 * 3 * 3;
        Assert.False(tri.Product.IsZero);
        Assert.Equal(-2 * 3 * 3584 * BigInteger.Pow(4, 3) * w * tri.Product, reading.DegreeTwo);
    }

    /// <summary>The bond's direction: (a, b, J, D) and (b, a, J, −D) are the same bond, in the closed form and in the
    /// generator.</summary>
    [Fact]
    public void ReversedBond_IsTheSameBond()
    {
        var cool = new long[] { 2, 1, 3 };
        var heat = new long[] { 1, 3, 1 };
        var listed = new[] { new B(0, 1, 1, 2, 5), new B(1, 2, 3, -1, -2), new B(0, 2, 2, 1, 1) };
        var reversed = new[] { new B(1, 0, 1, -2, 5), new B(1, 2, 3, -1, -2), new B(2, 0, 2, -1, 1) };
        var one = TransportSeventhMomentClaim.SeventhDegreeTwo(3, listed, cool, heat);
        Assert.False(one.IsZero);
        Assert.Equal(one, TransportSeventhMomentClaim.SeventhDegreeTwo(3, reversed, cool, heat));
        Assert.Equal(one, TransportSeventhMomentWitness.Read(new("reversed listing", 3, reversed, cool, heat)).DegreeTwo);
    }

    /// <summary>Five sites, past the witness's generator: the closed form on the bonds of the two rows at N = 5 of
    /// the proof's gate <c>simulations/seventh_moment_transport.py</c> (a ZZ bond off a triangle; a ring of five with
    /// a chord), at the rates γ↓, γ↑ given here. The numbers are what that gate's generator gives at these rates (its
    /// <c>rate_degrees</c> with s = γ↓ + γ↑, a = γ↓ − γ↑), computed once and written here; the gate's own rows read
    /// the same lists as (s, a) and so print other numbers.</summary>
    [Fact]
    public void FiveSites_AreTheGatesValues()
    {
        var offTriangle = new[]
        {
            new B(0, 1, 1, 1, 1), new B(1, 2, 2, 0, -2), new B(2, 0, 0, 1, 3),
            new B(3, 4, 0, 0, 3), new B(0, 3, 0, 0, 2), new B(2, 4, 0, 0, -1),
        };
        Assert.Equal(BigInteger.Parse("40340815872"),
            TransportSeventhMomentClaim.SeventhDegreeTwo(5, offTriangle, new long[] { 2, 1, 3, 1, 4 }, new long[] { 1, 3, 1, 2, 1 }));
        var chord = new[]
        {
            new B(0, 1, 1, 1, 2), new B(1, 2, 2, 0, -1), new B(2, 3, 1, -1, 3), new B(3, 4, 0, 1, 1), new B(4, 0, 1, 0, -2),
            new B(0, 2, 1, 2, 1), new B(1, 3, 0, 0, 2), new B(1, 4, 0, 0, -3),
        };
        Assert.Equal(BigInteger.Parse("-19502465024"),
            TransportSeventhMomentClaim.SeventhDegreeTwo(5, chord, new long[] { 2, 1, 3, 1, 4 }, new long[] { 1, 3, 1, 2, 2 }));
    }

    /// <summary>A ring of five alone: 2520·4^N·Im(loop)·(pentagram − pentagon) in Tr((L + c/2)⁷), the two polygons
    /// of the five rate points counted the way the loop runs, each as the sum of γ↓_x·γ↑_y − γ↓_y·γ↑_x over its
    /// edges; the gate reads the same number off its generator.</summary>
    [Fact]
    public void RingOfFive_IsPentagramMinusPentagon()
    {
        var ring = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 1, 0), new B(2, 3, 1, 0, 0), new B(3, 4, 2, 0, 0), new B(4, 0, 1, 2, 0) };
        long[] d = { 3, 1, 4, 1, 5 }, u = { 2, 7, 1, 8, 2 };
        var t = FluxAreaFifthMomentClaim.HopAmplitudes(5, ring);
        var loop = t[0, 1] * t[1, 2] * t[2, 3] * t[3, 4] * t[4, 0];
        BigInteger Cross(int x, int y) => (BigInteger)d[x] * u[y] - (BigInteger)d[y] * u[x];
        BigInteger pentagon = 0, pentagram = 0;
        for (int x = 0; x < 5; x++)
        {
            pentagon += Cross(x, (x + 1) % 5);
            pentagram += Cross(x, (x + 2) % 5);
        }
        var want = 128 * 2520 * BigInteger.Pow(4, 5) * loop.Im * (pentagram - pentagon);
        Assert.Equal(BigInteger.Parse("41617981440"), want);
        Assert.Equal(want, TransportSeventhMomentClaim.SeventhDegreeTwo(5, ring, d, u));
    }

    [Fact]
    public void BadInput_Throws()
    {
        Assert.Throws<System.ArgumentException>(() => TransportSeventhMomentClaim.SeventhDegreeTwo(3, new[] { new B(0, 3, 1, 0, 0) }, new long[] { 1, 1, 1 }, new long[] { 0, 0, 0 }));
        Assert.Throws<System.ArgumentException>(() => TransportSeventhMomentClaim.SeventhDegreeTwo(3, new[] { new B(0, 1, 1, 0, 0) }, new long[] { 1, 1 }, new long[] { 0, 0, 0 }));
        Assert.Throws<System.ArgumentException>(() => TransportSeventhMomentClaim.ThriceBeta(3, new[] { new B(0, 1, 1, 0, 0) }, 1, 1));
    }
}
