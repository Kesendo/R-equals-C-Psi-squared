using System.Numerics;
using MirrorWorld;
using static MirrorWorldTests.LetterCubeTests;

namespace MirrorWorldTests;

// From-below guard for the two readings of a colour on the circle that TiltedAxes carries beside its draws: the
// half-turn law (a product of circle letters commutes with a bond exactly when the half-turns about its letters leave
// the coupling matrix as it was) and the Heisenberg bond equations (their solutions counted against the colours).
public class TiltedAxesTests
{
    // Conjugating by a·σ turns σ by π about a, so (a·σ)⊗(b·σ) commutes with Σ W_kl σ_k⊗σ_l exactly when
    // R_a(π)·W·R_b(π)ᵀ = W, with |a|²·R_a(π) = 2aaᵀ − |a|²; dense matrices against the integer identity, on random
    // integer W (diagonal and general) and axes, a third of them drawn with b along Wᵀ·a so that both answers occur.
    // The half-turn law is the cross form (Wᵀ·a)×b = 0 = (W·b)×a, a on the first site, for every W; F138 writes it
    // for diagonal W, without the transpose, and without it the form accepts XX + XY under X and X.
    [Fact]
    public void A_Product_Of_Circle_Letters_Commutes_With_A_Bond_Exactly_When_The_Half_Turns_Keep_Its_Matrix()
    {
        var rnd = new Random(20261027);
        var letters = new[] { Letter('X'), Letter('Y'), Letter('Z') };
        Complex[,] Dir(long[] v) { var m = new Complex[2, 2]; for (int k = 0; k < 3; k++) Acc(m, letters[k], v[k]); return m; }
        int yes = 0, no = 0;
        for (int t = 0; t < 600; t++)
        {
            bool diagonal = t % 2 == 0;
            var W = new long[3, 3];
            for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) W[k, l] = diagonal && k != l ? 0 : rnd.Next(-2, 3);
            long[] a = Enumerable.Range(0, 3).Select(_ => (long)rnd.Next(-2, 3)).ToArray();
            long[] b = Enumerable.Range(0, 3).Select(_ => (long)rnd.Next(-2, 3)).ToArray();
            if (t % 3 == 0) { var wa = Apply(Transpose3(W), a); if (wa.Any(x => x != 0)) b = wa; }
            if (a.All(x => x == 0) || b.All(x => x == 0)) continue;

            var bond = new Complex[4, 4];
            for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) Acc(bond, Kron(letters[k], letters[l]), W[k, l]);
            var g = Kron(Dir(a), Dir(b));
            bool commutes = Same(Mul(g, bond), Mul(bond, g));

            long na = a.Sum(x => x * x), nb = b.Sum(x => x * x);
            var turned = Mul3(HalfTurn(a), W, Transpose3(HalfTurn(b)));
            bool halfTurn = true;
            for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) halfTurn &= turned[k, l] == na * nb * W[k, l];
            Assert.True(commutes == halfTurn, $"W {string.Join(",", W.Cast<long>())}, a {string.Join(",", a)}, b {string.Join(",", b)}");
            if (commutes) yes++; else no++;

            bool cross = Cross(Apply(Transpose3(W), a), b).All(x => x == 0) && Cross(Apply(W, b), a).All(x => x == 0);
            Assert.Equal(commutes, cross);
        }
        Assert.True(yes > 20 && no > 20, $"yes {yes}, no {no}");
        // without the transpose the form needs a symmetric W: XX + XY under X and X has W·a along b and W·b along a,
        // yet X⊗X does not commute with X⊗Y, and the transposed form rejects it
        var w0 = new long[,] { { 1, 1, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
        long[] a0 = { 1, 0, 0 }, b0 = { 1, 0, 0 };
        Assert.True(Cross(Apply(w0, a0), b0).All(x => x == 0) && Cross(Apply(w0, b0), a0).All(x => x == 0));
        Assert.False(Cross(Apply(Transpose3(w0), a0), b0).All(x => x == 0));
        var bond0 = new Complex[4, 4];
        Acc(bond0, Kron(letters[0], letters[0]), 1); Acc(bond0, Kron(letters[0], letters[1]), 1);
        var g0 = Kron(Dir(a0), Dir(b0));
        Assert.False(Same(Mul(g0, bond0), Mul(bond0, g0)));
    }

    // Heisenberg bonds under integer axes, every site watched: the bond equations' solutions in the far span are
    // exactly the colours, Sym^N of the common perpendicular (N + 1 when the axes span a line, 1 when a plane, 0
    // otherwise), 30 rows on each of the chain at N = 2 to 5, the ring at N = 3 to 5 (at N = 2 it is one bond twice)
    // and the star and K4 at N = 4, at random signed weights, the axes generic, coplanar, parallel or in clusters of
    // two and three. The GF(p) nullity is an upper bound and the colours
    // a lower bound, so equality is exact. That the colours are the only solutions is measured here, not proven.
    [Fact]
    public void Heisenberg_Bond_Equations_Admit_Exactly_The_Colours()
    {
        var rnd = new Random(20261027);
        var dirs = TiltedAxes.Directions();
        var seen = new HashSet<string>();
        for (int n = 2; n <= 5; n++)
            foreach (var (name, edges) in TiltedAxes.Graphs(n).Where(g => !(n == 2 && g.Name == "ring")))
                for (int t = 0; t < 30; t++)
                {
                    var a = dirs[rnd.Next(dirs.Count)];
                    var b = dirs[rnd.Next(dirs.Count)];
                    var c = dirs[rnd.Next(dirs.Count)];
                    (long, long, long) InPlane() { long p = rnd.Next(1, 3), q = rnd.Next(-2, 3); var v = (p * a.X + q * b.X, p * a.Y + q * b.Y, p * a.Z + q * b.Z); return v == (0, 0, 0) ? a : v; }
                    var axes = (t % 5) switch
                    {
                        0 => Enumerable.Range(0, n).Select(_ => dirs[rnd.Next(dirs.Count)]).ToArray(),
                        1 => Enumerable.Range(0, n).Select(_ => InPlane()).ToArray(),
                        2 => Enumerable.Repeat(a, n).ToArray(),
                        3 => Enumerable.Range(0, n).Select(_ => rnd.Next(2) == 0 ? a : b).ToArray(),
                        _ => Enumerable.Range(0, n).Select(_ => rnd.Next(3) switch { 0 => a, 1 => b, _ => c }).ToArray(),
                    };
                    var weights = edges.Select(_ => (long)(rnd.Next(2) == 0 ? -1 : 1) * rnd.Next(1, 4)).ToArray();
                    int colours = TiltedAxes.ColourDimension(n, axes);
                    int nullity = TiltedAxes.BondNullity(n, edges, weights, axes);
                    Assert.True(nullity == colours, $"N = {n} {name}, axes {string.Join(" ", axes)}: nullity {nullity}, colours {colours}");
                    seen.Add(colours == 0 ? "none" : colours == 1 ? "plane" : "line");
                }
        Assert.Equal(3, seen.Count);

        // pinned: the chain under X, Z, Y (three letters spanning space) has no solution, under three coplanar axes
        // at other angles one, their normal (4, -3, -2) on every site
        Assert.Equal(0, TiltedAxes.BondNullity(3, new[] { (0, 1), (1, 2) }, new[] { 1L, 1 }, new (long, long, long)[] { (1, 0, 0), (0, 0, 1), (0, 1, 0) }));
        Assert.Equal(1, TiltedAxes.BondNullity(3, new[] { (0, 1), (1, 2) }, new[] { 1L, 1 }, new (long, long, long)[] { (2, 2, 1), (1, 2, -1), (1, 0, 2) }));
    }

    static long[] Apply(long[,] m, long[] v) => Enumerable.Range(0, 3).Select(k => m[k, 0] * v[0] + m[k, 1] * v[1] + m[k, 2] * v[2]).ToArray();
    static long[] Cross(long[] u, long[] v) => new[] { u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0] };
    static long[,] HalfTurn(long[] v) { long s = v.Sum(x => x * x); var r = new long[3, 3]; for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) r[k, l] = 2 * v[k] * v[l] - (k == l ? s : 0); return r; }
    static long[,] Transpose3(long[,] m) { var r = new long[3, 3]; for (int k = 0; k < 3; k++) for (int l = 0; l < 3; l++) r[k, l] = m[l, k]; return r; }
    static long[,] Mul3(long[,] a, long[,] b, long[,] c)
    {
        var ab = new long[3, 3]; var r = new long[3, 3];
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) ab[i, j] += a[i, k] * b[k, j];
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) r[i, j] += ab[i, k] * c[k, j];
        return r;
    }
    static bool Same(Complex[,] a, Complex[,] b)
    {
        for (int i = 0; i < a.GetLength(0); i++) for (int k = 0; k < a.GetLength(1); k++) if (a[i, k] != b[i, k]) return false;
        return true;
    }
    static void Acc(Complex[,] m, Complex[,] add, Complex k)
    {
        for (int i = 0; i < m.GetLength(0); i++) for (int j = 0; j < m.GetLength(1); j++) m[i, j] += k * add[i, j];
    }
}
