using System.Numerics;
using Xunit;

namespace MirrorWorld.Tests;

// The fifth trace read off the triangles (F137, Sextant.FluxAreaFifth): the closed form against the sextant's own
// centred traces modulo both primes, the recorded value, the class fence, the decision, and a setting past the wall.
public class SextantFluxAreaTests
{
    // J(X_aX_b + Y_aY_b) + D(X_aY_b - Y_aX_b) + delta Z_aZ_b
    static IEnumerable<(PauliString, long)> Bond(int n, int a, int b, long j, long d, long delta)
    {
        string Pair(char pa, char pb)
        {
            var s = new string('I', n).ToCharArray();
            s[a] = pa;
            s[b] = pb;
            return new string(s);
        }
        if (j != 0) { yield return (PauliString.Parse(Pair('X', 'X')), j); yield return (PauliString.Parse(Pair('Y', 'Y')), j); }
        if (d != 0) { yield return (PauliString.Parse(Pair('X', 'Y')), d); yield return (PauliString.Parse(Pair('Y', 'X')), -d); }
        if (delta != 0) yield return (PauliString.Parse(Pair('Z', 'Z')), delta);
    }

    static List<(int, LetterCube.Gauss[,], long)> Bath(long[] cool, long[] heat)
    {
        var list = new List<(int, LetterCube.Gauss[,], long)>();
        for (int l = 0; l < cool.Length; l++)
        {
            list.Add((l, LetterCube.Lowering, cool[l]));
            list.Add((l, LetterCube.Raising, heat[l]));
        }
        return list;
    }

    static List<(PauliString, long)> UnitTriangle(int n = 3)
        => Bond(n, 0, 1, 1, 0, 0).Concat(Bond(n, 1, 2, 1, 0, 0)).Concat(Bond(n, 2, 0, 0, 1, 0)).ToList();

    static long Mod(BigInteger x, long p) => (long)(((x % p) + p) % p);

    // Tr((4 dim M)^5) = 32 dim^5 Tr((2L + c)^5): the sextant's fifth centred trace modulo p against the closed form
    static void AssertMatchesTheTraces(int n, List<(PauliString, long)> h, List<(int, LetterCube.Gauss[,], long)> jumps)
    {
        var fifth = Sextant.FluxAreaFifth(n, h, jumps);
        Assert.NotNull(fifth);
        var s = new Sextant.Setting(n, h, jumps);
        BigInteger dim = BigInteger.Pow(4, n);
        foreach (long p in ModP.Primes)
        {
            var traces = Sextant.CentredTracesModP(s, 5, p).ToDictionary(t => t.Power, t => t.Trace);
            Assert.Equal(0, traces[1]);
            Assert.Equal(0, traces[3]);
            Assert.Equal(Mod(32 * BigInteger.Pow(dim, 5) * fifth!.Value, p), traces[5]);
        }
    }

    [Fact]
    public void RecordedValue_OnF137sTriangle()
    {
        // rates times 120: cooling (1/2, 1, 3/2), heating (1/8, 1/3, 1/5); Tr((L + c/2)^5) = -656 at flux pi/2
        var jumps = Bath(new long[] { 60, 120, 180 }, new long[] { 15, 40, 24 });
        Assert.Equal(new BigInteger(-656) * 32 * 120 * 120, Sextant.FluxAreaFifth(3, UnitTriangle(), jumps));
        AssertMatchesTheTraces(3, UnitTriangle(), jumps);
        // loss and gain exchanged: the sign
        Assert.Equal(new BigInteger(656) * 32 * 120 * 120, Sextant.FluxAreaFifth(3, UnitTriangle(), Bath(new long[] { 15, 40, 24 }, new long[] { 60, 120, 180 })));
    }

    [Fact]
    public void SeveralTriangles_MatchTheTraces()
    {
        // the complete graph on four sites, complex hoppings, unequal ZZ
        var k4 = Bond(4, 0, 1, 1, 2, 1).Concat(Bond(4, 0, 2, 2, -1, -2)).Concat(Bond(4, 0, 3, 1, 1, 0))
            .Concat(Bond(4, 1, 2, 3, 1, 2)).Concat(Bond(4, 1, 3, -1, -2, 1)).Concat(Bond(4, 2, 3, 2, 3, -1)).ToList();
        var jumps = Bath(new long[] { 2, 1, 3, 1 }, new long[] { 1, 3, 1, 2 });
        Assert.False(Sextant.FluxAreaFifth(4, k4, jumps)!.Value.IsZero);
        AssertMatchesTheTraces(4, k4, jumps);
        // two plaquettes and a pendant on five sites
        var five = Bond(5, 0, 1, 1, 0, 1).Concat(Bond(5, 1, 2, 0, 1, 0)).Concat(Bond(5, 0, 2, 1, 1, 2))
            .Concat(Bond(5, 2, 3, 1, 0, -1)).Concat(Bond(5, 0, 3, 2, 1, 0)).Concat(Bond(5, 3, 4, 1, -1, 1)).ToList();
        var jumps5 = Bath(new long[] { 3, 1, 2, 5, 1 }, new long[] { 1, 2, 2, 1, 4 });
        Assert.False(Sextant.FluxAreaFifth(5, five, jumps5)!.Value.IsZero);
        AssertMatchesTheTraces(5, five, jumps5);
    }

    [Fact]
    public void SilentSettings_AreZeroInTheTracesToo()
    {
        var m11 = Bond(3, 0, 1, 1, 0, 7).Concat(Bond(3, 1, 2, 2, 0, 7)).Concat(Bond(3, 2, 0, 0, 3, 7)).ToList();
        foreach (var (cool, heat) in new[]
                 {
                     (new long[] { 2, 4, 6 }, new long[] { 1, 2, 3 }),          // one ratio
                     (new long[] { 2, 1, 3 }, new long[] { 0, 0, 0 }),          // cooling alone
                     (new long[] { 60, 120, 180 }, new long[] { 40, 60, 80 }),  // a line that misses the origin
                     (new long[] { 4, 6, 8 }, new long[] { 2, 4, 6 }),          // a temperature gradient at one coupling
                 })
        {
            Assert.Equal(BigInteger.Zero, Sextant.FluxAreaFifth(3, m11, Bath(cool, heat)));
            AssertMatchesTheTraces(3, m11, Bath(cool, heat));
        }
        // a ring of four has no triangle
        var ring = Bond(4, 0, 1, 1, 0, 1).Concat(Bond(4, 1, 2, 2, 0, 1)).Concat(Bond(4, 2, 3, 3, 0, 1)).Concat(Bond(4, 3, 0, 0, 2, 1)).ToList();
        var jumps4 = Bath(new long[] { 2, 1, 3, 1 }, new long[] { 1, 3, 1, 2 });
        Assert.Equal(BigInteger.Zero, Sextant.FluxAreaFifth(4, ring, jumps4));
        AssertMatchesTheTraces(4, ring, jumps4);
    }

    [Fact]
    public void OutsideTheClass_AnswersNull()
    {
        var jumps = Bath(new long[] { 2, 1, 3 }, new long[] { 1, 3, 1 });
        var tri = UnitTriangle();
        Assert.NotNull(Sextant.FluxAreaFifth(3, tri, jumps));
        // a field, a uniform one included
        Assert.Null(Sextant.FluxAreaFifth(3, tri.Append((PauliString.Parse("ZII"), 1L)), jumps));
        Assert.Null(Sextant.FluxAreaFifth(3, tri.Concat(new[] { (PauliString.Parse("ZII"), 1L), (PauliString.Parse("IZI"), 1L), (PauliString.Parse("IIZ"), 1L) }), jumps));
        // XX without YY, XY without -YX, another letter pair, a three-site term
        Assert.Null(Sextant.FluxAreaFifth(3, tri.Append((PauliString.Parse("XXI"), 2L)), jumps));
        Assert.Null(Sextant.FluxAreaFifth(3, tri.Append((PauliString.Parse("IXY"), 1L)), jumps));
        Assert.Null(Sextant.FluxAreaFifth(3, tri.Append((PauliString.Parse("XZI"), 1L)), jumps));
        Assert.Null(Sextant.FluxAreaFifth(3, tri.Append((PauliString.Parse("ZZZ"), 1L)), jumps));
        // dephasing beside the bath
        Assert.Null(Sextant.FluxAreaFifth(3, tri, jumps.Append((0, LetterCube.LetterMatrix('Z'), 1L))));
        // a jump at rate zero is no jump, and a constant does not enter
        Assert.Equal(Sextant.FluxAreaFifth(3, tri, jumps),
            Sextant.FluxAreaFifth(3, tri.Append((PauliString.Parse("III"), 5L)), jumps.Append((0, LetterCube.LetterMatrix('Z'), 0L))));
    }

    [Fact]
    public void TheSameBondWrittenTwice_AddsUp()
    {
        var jumps = Bath(new long[] { 2, 1, 3 }, new long[] { 1, 3, 1 });
        var once = Bond(3, 0, 1, 2, 4, 0).Concat(Bond(3, 1, 2, 1, 0, 0)).Concat(Bond(3, 2, 0, 0, 1, 0)).ToList();
        // (0, 1) in two halves, one of them listed from site 1: (b, a, J, -D) is (a, b, J, D)
        var twice = Bond(3, 0, 1, 1, 2, 0).Concat(Bond(3, 1, 0, 1, -2, 0)).Concat(Bond(3, 1, 2, 1, 0, 0)).Concat(Bond(3, 2, 0, 0, 1, 0)).ToList();
        Assert.False(Sextant.FluxAreaFifth(3, once, jumps)!.Value.IsZero);
        Assert.Equal(Sextant.FluxAreaFifth(3, once, jumps), Sextant.FluxAreaFifth(3, twice, jumps));
        AssertMatchesTheTraces(3, twice, jumps);
    }

    [Fact]
    public void ATermPastTheLastSite_Throws()
    {
        // XXIZ on three sites is not a bond on (0, 1)
        var jumps = Bath(new long[] { 2, 1, 3 }, new long[] { 1, 3, 1 });
        var h = UnitTriangle().Append((PauliString.Parse("XXIZ"), 1L)).Append((PauliString.Parse("YYIZ"), 1L));
        Assert.Throws<ArgumentException>(() => Sextant.FluxAreaFifth(3, h, jumps));
    }

    [Fact]
    public void Decide_SaysObstructedExactly_AndFallsThroughWhereTheFifthIsSilent()
    {
        var speaking = Sextant.Decide(new Sextant.Setting(3, UnitTriangle(), Bath(new long[] { 2, 1, 3 }, new long[] { 1, 3, 1 })));
        Assert.Equal(Sextant.Palindrome.Obstructed, speaking.Kind);
        Assert.Contains("F137", speaking.Detail);
        // the temperature gradient at one coupling: the fifth is silent, the seventh obstructs
        var gradient = new Sextant.Setting(3, UnitTriangle(), Bath(new long[] { 4, 6, 8 }, new long[] { 2, 4, 6 }));
        Assert.Equal(BigInteger.Zero, Sextant.FluxAreaFifth(3, gradient.H, gradient.Jumps));
        var decided = Sextant.Decide(gradient);
        Assert.Equal(Sextant.Palindrome.Obstructed, decided.Kind);
        Assert.Contains("Tr(M^7)", decided.Detail);
        // one ratio pairs: the reading must not stand in the way of a witness
        var oneRatio = Sextant.Decide(new Sextant.Setting(3, UnitTriangle(), Bath(new long[] { 2, 4, 6 }, new long[] { 1, 2, 3 })));
        Assert.NotEqual(Sextant.Palindrome.Obstructed, oneRatio.Kind);
    }

    [Fact]
    public void PastTheWall_TwentySites()
    {
        // a strip of eighteen triangles (l, l+1, l+2), every bond J = 1 except D = 1 on the long edges (l+2, l): each
        // triangle has loop i; no generator is built
        int n = 20;
        var h = new List<(PauliString, long)>();
        for (int l = 0; l + 1 < n; l++) h.AddRange(Bond(n, l, l + 1, 1, 0, 1));
        for (int l = 0; l + 2 < n; l++) h.AddRange(Bond(n, l + 2, l, 0, 1, 0));
        var cool = Enumerable.Range(0, n).Select(l => (long)(l * l + 1)).ToArray();
        var heat = Enumerable.Range(0, n).Select(l => (long)(l + 2)).ToArray();
        BigInteger sum = 0;
        for (int l = 0; l + 2 < n; l++)
            sum += (cool[l + 1] * heat[l + 2] - cool[l + 2] * heat[l + 1]) - (cool[l] * heat[l + 2] - cool[l + 2] * heat[l]) + (cool[l] * heat[l + 1] - cool[l + 1] * heat[l]);
        Assert.False(sum.IsZero);
        Assert.Equal(1920 * BigInteger.Pow(4, n) * sum, Sextant.FluxAreaFifth(n, h, Bath(cool, heat)));
    }
}
