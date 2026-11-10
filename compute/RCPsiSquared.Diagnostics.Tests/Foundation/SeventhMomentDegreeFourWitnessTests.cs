using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using B = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of <see cref="SeventhMomentDegreeFourClaim"/> (Theorems 1 and 2 of
/// PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE, F137): the part of degree four in the rates of Tr((2L + c)⁷), read from
/// the generator in exact Gaussian integers as the words with four rate letters, against the sum over the triangles.
/// Every comparison is exact.</summary>
public class SeventhMomentDegreeFourWitnessTests
{
    private static readonly SeventhMomentDegreeFourWitness W = new();

    [Fact]
    public void EveryRow_Holds()
    {
        foreach (var x in W.Readings)
            Assert.True(x.Holds, $"{x.Row.Label}: degree four {x.DegreeFour}, sum over the triangles {x.Law}, trace {x.Trace}, together {x.WholeLaw}");
    }

    [Fact]
    public void Rows_AreNotVacuous()
    {
        Assert.True(W.NotVacuous, W.Summary);
        // the silent rows are silent in the generator, not only in the closed form
        Assert.All(W.Readings.Where(x => !x.Row.IsControl && x.Law.IsZero), x => Assert.True(x.DegreeFour.IsZero, x.Row.Label));
        // and the speaking rows carry a part of degree two as well somewhere, so "the whole trace" is not degree four alone
        Assert.Contains(W.Readings, x => !x.Row.IsControl && !x.Law.IsZero && x.Trace != x.DegreeFour);
    }

    /// <summary>The recorded values as numbers written here: F137's 430080 and the transport proof's 322560, all of
    /// degree four.</summary>
    [Fact]
    public void RecordedValues_AreReproduced()
    {
        SeventhMomentDegreeFourWitness.Reading Row(string start) => W.Readings.Single(x => x.Row.Label.StartsWith(start));
        var gradient = Row("F137's temperature gradient");
        Assert.Equal((BigInteger)430080 * 128, gradient.DegreeFour);
        Assert.Equal(gradient.DegreeFour, gradient.Trace);
        var skew = Row("hoppings 1, 1 + i, 2i");
        Assert.Equal((BigInteger)322560 * 128, skew.DegreeFour);
        Assert.Equal(skew.DegreeFour, skew.Trace);
    }

    /// <summary>On the triangle the sum over the triangles is the proof's other written form,
    /// 13440·Im·Σ_cyc ω_{l,l+1}·(s_l² + s_{l+1}² + 2·s_{l+2}²), here times 128 with s = S/2:
    /// 430080·Im·Σ_cyc ω_{l,l+1}·(S_l² + S_{l+1}² + 2·S_{l+2}²), written out independently.</summary>
    [Fact]
    public void Triangle_IsTheSecondWrittenForm()
    {
        int checkedRows = 0;
        foreach (var row in SeventhMomentDegreeFourWitness.Rows.Where(r => r.N == 3 && !r.IsControl))
        {
            var t = FluxAreaFifthMomentClaim.HopAmplitudes(3, row.Bonds);
            var im = (t[0, 1] * t[1, 2] * t[2, 0]).Im;
            BigInteger S(int l) => (BigInteger)row.Cool[l] + row.Heat[l];
            BigInteger Om(int a, int b) => (BigInteger)row.Cool[a] * row.Heat[b] - (BigInteger)row.Cool[b] * row.Heat[a];
            var sum = BigInteger.Zero;
            for (int l = 0; l < 3; l++)
            {
                int k = (l + 1) % 3, j = (l + 2) % 3;
                sum += Om(l, k) * (S(l) * S(l) + S(k) * S(k) + 2 * S(j) * S(j));
            }
            Assert.Equal(430080 * im * sum, SeventhMomentDegreeFourClaim.SeventhDegreeFour(3, row.Bonds, row.Cool, row.Heat));
            if (!sum.IsZero && !im.IsZero) checkedRows++;
        }
        Assert.True(checkedRows >= 3, $"{checkedRows}");
    }

    /// <summary>Theorem 2 by hand on the transport proof's triangle (Im = 2, ⟨γ↓, b⟩ = 0 at this point): on γ↑ = γ↓ − 2
    /// (e = −2, m = 1) the whole seventh trace is −e·3360·Im·(1 + m)²·V(γ↓) = 2·3360·2·4·6 = 322560 at γ↓ = (3, 5, 2).</summary>
    [Fact]
    public void LineForm_GivesThe322560()
    {
        var skew = new[] { new B(0, 1, 1, 0, 0), new B(1, 2, 1, 1, 0), new B(2, 0, 0, 2, 0) };
        long e = -2, m = 1, im = 2;
        long[] gd = { 3, 5, 2 };
        long v = (gd[0] - gd[1]) * (gd[1] - gd[2]) * (gd[2] - gd[0]);
        Assert.Equal(6, v);
        Assert.Equal((BigInteger)(-e * 3360 * im * (1 + m) * (1 + m) * v) * 128,
            SeventhMomentDegreeFourClaim.SeventhCentredTrace(3, skew, gd, new long[] { 1, 3, 0 }));
    }

    /// <summary>Theorem 2 with the flow b not zero, against the generator: the proof's triangle t₀₁ = t₁₂ = 2 + i,
    /// t₂₀ = 3 + 4i with ZZ 2, 2 and 3 (Im = 24, b = Im·(473⁄3, 0, −473⁄3) = (3784, 0, −3784)) on the line γ↑ = 1 + γ↓
    /// (e = 1, m = 1) at γ↓ = (1, 4, 2): −e·[3360·Im·(1 + m)²·V(γ↓) + 10752·⟨γ↓, b⟩] with V = −6 and ⟨γ↓, b⟩ = −3784,
    /// times 128, written out here from the proof's numbers.</summary>
    [Fact]
    public void LineForm_WithAFlow()
    {
        var bonds = new[] { new B(0, 1, 2, 1, 2), new B(1, 2, 2, 1, 2), new B(0, 2, 3, -4, 3) };
        long[] gd = { 1, 4, 2 }, gu = { 2, 5, 3 };
        BigInteger e = 1, m = 1, im = 24, v = -6, flow = 3784 * 1 + 0 * 4 - 3784 * 2;
        var want = -e * (3360 * im * (1 + m) * (1 + m) * v + 10752 * flow) * 128;
        var reading = SeventhMomentDegreeFourWitness.Read(new("the proof's triangle with a flow, on a line", 3, bonds, gd, gu));
        Assert.True(reading.Holds);
        Assert.Equal(want, reading.Trace);
        Assert.NotEqual(reading.Trace, reading.DegreeFour);   // the part of degree two speaks here
    }

    /// <summary>The constants can fail: the sum over the triangles with 630, 420 or 210 moved by one misses the
    /// generator's part of degree four on some speaking row (written here as an independent helper, checked to equal the
    /// claim at the right constants).</summary>
    [Fact]
    public void WrongForms_MissARow()
    {
        BigInteger Form(SeventhMomentDegreeFourWitness.Row r, int p, int q, int w)
        {
            var t = FluxAreaFifthMomentClaim.HopAmplitudes(r.N, r.Bonds);
            var s = Enumerable.Range(0, r.N).Select(x => (BigInteger)r.Cool[x] + r.Heat[x]).ToArray();
            var all = s.Aggregate(BigInteger.Zero, (acc, x) => acc + x * x);
            var total = BigInteger.Zero;
            for (int u = 0; u < r.N; u++)
                for (int v = u + 1; v < r.N; v++)
                    for (int x = v + 1; x < r.N; x++)
                    {
                        var im = (t[u, v] * t[v, x] * t[x, u]).Im;
                        int[] c = { u, v, x };
                        for (int i = 0; i < 3; i++)
                        {
                            int l = c[i], k = c[(i + 1) % 3], j = c[(i + 2) % 3];
                            var om = (BigInteger)r.Cool[l] * r.Heat[k] - (BigInteger)r.Cool[k] * r.Heat[l];
                            total += im * om * (p * all - q * (s[l] * s[l] + s[k] * s[k]) - w * s[j] * s[j]);
                        }
                    }
            return 32 * BigInteger.Pow(4, r.N) * total;
        }
        var speaking = W.Readings.Where(x => !x.Row.IsControl && !x.Law.IsZero).ToArray();
        Assert.All(speaking, x => Assert.Equal(x.DegreeFour, Form(x.Row, 630, 420, 210)));
        foreach (var (p, q, w) in new[] { (631, 420, 210), (630, 421, 210), (630, 420, 211) })
            Assert.Contains(speaking, x => x.DegreeFour != Form(x.Row, p, q, w));
    }

    /// <summary>The two controls read differently on the part of degree four: 3·X₀X₁ moves it off the sum over the
    /// triangles, a Z field on one site leaves it there at this row (both move the whole trace).</summary>
    [Fact]
    public void Controls_ReadOnDegreeFour()
    {
        Assert.True(W.Readings.Single(x => x.Row.Label.Contains("X₀X₁")).DegreeFourMisses);
        Assert.False(W.Readings.Single(x => x.Row.Label.Contains("Z field")).DegreeFourMisses);
    }

    /// <summary>The shear δ_l → δ_l + ν·s_l at ν = −1, in the rates γ↓ → γ↓ − S/2 and γ↑ → γ↑ + S/2 with S = γ↓ + γ↑
    /// (even here, and γ↓ ≥ γ↑, so the sheared rates are rates), keeps every S and every ω, so the closed forms of degree
    /// two and four do not move; the generator agrees on the whole seventh trace and on its part of degree four, and its
    /// ninth trace moves.</summary>
    [Fact]
    public void Shear_KeepsTheSeventhTrace()
    {
        var bonds = new[] { new B(0, 1, 2, 1, 3), new B(1, 2, 1, -2, -1), new B(2, 0, 3, 1, 2) };
        long[] cool = { 3, 5, 4 }, heat = { 1, 1, 2 };
        long[] cool2 = { 1, 2, 1 }, heat2 = { 3, 4, 5 };   // γ↓ − S/2, γ↑ + S/2 with S = 4, 6, 6
        var before = SeventhMomentDegreeFourWitness.Read(new("before", 3, bonds, cool, heat));
        var after = SeventhMomentDegreeFourWitness.Read(new("sheared", 3, bonds, cool2, heat2));
        Assert.False(before.Trace.IsZero);
        Assert.False(before.DegreeFour.IsZero);
        Assert.Equal(before.Trace, after.Trace);
        Assert.Equal(before.DegreeFour, after.DegreeFour);
        var terms = FluxAreaFifthMomentClaim.Terms(3, bonds);
        var ninth = (System.Func<long[], long[], BigInteger>)((c, h) =>
            ExactThermalGenerator.OddPowerSums(ExactThermalGenerator.Build(terms, c, h), 9)[9]);
        Assert.NotEqual(ninth(cool, heat), ninth(cool2, heat2));
    }

    /// <summary>The bond's direction: (a, b, J, D) and (b, a, J, −D) are the same bond, in the closed form and in the
    /// generator.</summary>
    [Fact]
    public void ReversedBond_IsTheSameBond()
    {
        var cool = new long[] { 3, 1, 4 };
        var heat = new long[] { 1, 5, 2 };
        var listed = new[] { new B(0, 1, 2, 1, 3), new B(1, 2, 1, -2, -1), new B(2, 0, 3, 1, 2) };
        var reversed = new[] { new B(1, 0, 2, -1, 3), new B(2, 1, 1, 2, -1), new B(2, 0, 3, 1, 2) };
        var one = SeventhMomentDegreeFourClaim.SeventhDegreeFour(3, listed, cool, heat);
        Assert.False(one.IsZero);
        Assert.Equal(one, SeventhMomentDegreeFourClaim.SeventhDegreeFour(3, reversed, cool, heat));
        Assert.Equal(one, SeventhMomentDegreeFourWitness.Read(new("reversed listing", 3, reversed, cool, heat)).DegreeFour);
    }

    /// <summary>Five sites, past the witness's generator: the proof's gate row (K) with two spectators, one joined by
    /// hopping and one by ZZ alone, at γ↓ = (2, 1, 3, 1, 4), γ↑ = (1, 3, 1, 2, 1). The number is what that gate's generator
    /// gives for the part of degree four (its <c>rate_parts</c>, times 128), computed once and written here.</summary>
    [Fact]
    public void FiveSites_IsTheGatesValue()
    {
        var bonds = new[] { new B(0, 1, 1, 1, 0), new B(1, 2, 1, 0, 0), new B(2, 0, 0, 2, 0), new B(3, 1, 1, 0, 0), new B(2, 4, 0, 0, 1), new B(0, 2, 0, 0, 1) };
        Assert.Equal(BigInteger.Parse(FiveSiteValue),
            SeventhMomentDegreeFourClaim.SeventhDegreeFour(5, bonds, new long[] { 2, 1, 3, 1, 4 }, new long[] { 1, 3, 1, 2, 1 }));
    }

    private const string FiveSiteValue = "-3605790720";

    [Fact]
    public void BadInput_Throws()
    {
        Assert.Throws<System.ArgumentException>(() => SeventhMomentDegreeFourClaim.SeventhDegreeFour(3, new[] { new B(0, 3, 1, 0, 0) }, new long[] { 1, 1, 1 }, new long[] { 0, 0, 0 }));
        Assert.Throws<System.ArgumentException>(() => SeventhMomentDegreeFourClaim.SeventhDegreeFour(3, new[] { new B(0, 1, 1, 0, 0) }, new long[] { 1, 1 }, new long[] { 0, 0, 0 }));
    }
}
