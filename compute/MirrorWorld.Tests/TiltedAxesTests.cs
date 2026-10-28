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
    // a lower bound, so equality is exact on each row; that the colours are the only solutions on every connected
    // graph is the theorem (docs/proofs/PROOF_HEISENBERG_COLOURS_OFF_THE_LETTERS.md), which these rows check.
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

    // F169's proof where its last step has content. The local step leaves a site alone only where its bonds meet
    // straight or at a right hinge (the planes of the axes meeting at a right angle); the vertex step then has two
    // words to cut, which needs the frames to agree around every cycle. Integer rows of that kind, off the letters,
    // with the axes spanning space and no bond on one axis, pin the theorem there, nullity 0 against no colour: the
    // six-ring right-hinged at every site, and graphs with three-neighbour sites, the theta graph, the five-ring with
    // a chord and K3,3. The five-ring right-hinged at every site is the other side: an odd cycle of right hinges
    // carries no word, which is also why K4 has no row of the first kind (a triangle on three axes that passes the
    // local step is straight at every site or right-hinged at every site). The theta graph and the chorded ring at generic and coplanar axes
    // too; the law holds on every connected graph and is silent on a disconnected one.
    [Fact]
    public void The_Colour_Law_Holds_Where_Its_Last_Step_Has_Content()
    {
        var ring5 = Enumerable.Range(0, 5).Select(i => (i, (i + 1) % 5)).ToArray();
        var ring6 = Enumerable.Range(0, 6).Select(i => (i, (i + 1) % 6)).ToArray();
        var hinge5 = new (long, long, long)[] { (0, 0, 1), (0, 1, -1), (2, -1, 1), (1, -2, -1), (1, 0, -1) };
        var hinge6 = new (long, long, long)[] { (0, 0, 1), (0, 1, -1), (1, -1, 1), (0, 1, 1), (1, 2, 1), (1, 0, 1) };
        long Dot((long, long, long) a, (long, long, long) b) => a.Item1 * b.Item1 + a.Item2 * b.Item2 + a.Item3 * b.Item3;
        (long, long, long) X((long, long, long) a, (long, long, long) b) => (a.Item2 * b.Item3 - a.Item3 * b.Item2, a.Item3 * b.Item1 - a.Item1 * b.Item3, a.Item1 * b.Item2 - a.Item2 * b.Item1);
        foreach (var axes in new[] { hinge5, hinge6 })
        {
            int L = axes.Length;
            for (int i = 0; i < L; i++)
            {
                var (a, b, c) = (axes[(i + L - 1) % L], axes[i], axes[(i + 1) % L]);
                Assert.Equal(0L, Dot(X(a, b), X(b, c)));                                      // a right hinge at i
                Assert.NotEqual((0L, 0L, 0L), X(b, c));                                        // neighbours on different axes
            }
            Assert.Equal(0, TiltedAxes.ColourDimension(L, axes));
            Assert.True(TiltedAxes.ColourLawProven(L, L == 5 ? ring5 : ring6, axes));
            Assert.Equal(0, TiltedAxes.BondNullity(L, L == 5 ? ring5 : ring6, Enumerable.Range(0, L).Select(k => (long)(k % 2 == 0 ? 1 : -2)).ToArray(), axes));
        }

        var theta = new[] { (0, 2), (2, 1), (0, 3), (3, 1), (0, 4), (4, 1) };
        var chord = ring5.Append((0, 2)).ToArray();
        var generic = new (long, long, long)[] { (1, 2, 0), (0, 1, 2), (2, 0, 1), (1, -1, 2), (2, 1, -1) };
        var coplanar = new (long, long, long)[] { (2, 2, 1), (1, 2, -1), (1, 0, 2), (2, 2, 1), (1, 0, 2) };
        foreach (var edges in new[] { theta, chord })
        {
            Assert.True(TiltedAxes.ColourLawProven(5, edges, generic));
            foreach (var axes in new[] { generic, coplanar })
                Assert.Equal(TiltedAxes.ColourDimension(5, axes), TiltedAxes.BondNullity(5, edges, new[] { 1L, -1, 2, 1, -2, 1, 1 }.Take(edges.Length).ToArray(), axes));
        }

        // two words reach the vertex step at three-neighbour sites; each row: straight or right at every site, the
        // frames agreeing around every cycle, no bond on one axis, the axes spanning space
        var k33 = new[] { (0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5) };
        var twoWords = new (int N, (int, int)[] Edges, (long, long, long)[] Axes, long[] Weights)[]
        {
            (5, theta, new (long, long, long)[] { (1, -1, 1), (-2, 2, -2), (1, -1, -1), (2, -2, 1), (0, -2, 1) }, new[] { 1L, 1, 1, -1, 2, -1 }),
            (5, chord, new (long, long, long)[] { (-2, -2, -2), (-1, 2, 2), (2, 0, 0), (2, -1, 1), (1, 0, 0) }, new[] { 1L, -1, -2, 1, 1, 2 }),
            (6, k33, new (long, long, long)[] { (-1, 1, 1), (-1, -2, 2), (1, -2, -2), (2, 0, 0), (1, 0, 0), (-1, 0, 0) }, new[] { 1L, -3, -2, -2, 3, -1, -2, -2, 1 }),
        };
        foreach (var (n, edges, axes, weights) in twoWords)
        {
            foreach (var (i, j) in edges) Assert.NotEqual((0L, 0L, 0L), X(axes[i], axes[j]));          // no bond on one axis
            for (int l = 0; l < n; l++)
            {
                var ms = edges.Where(e => e.Item1 == l || e.Item2 == l).Select(e => X(axes[e.Item1], axes[e.Item2])).ToArray();
                Assert.All(ms, a => Assert.All(ms, b => Assert.True(X(a, b) == (0, 0, 0) || Dot(a, b) == 0)));   // straight or right
            }
            // the two words survive: with frame e_l the first bond's line at l, a bond (i, j) asks x_i ⊕ x_j = t_ij,
            // t_ij whether its line is e at exactly one end; a walk from site 0 meets no contradiction
            var frame = new (long, long, long)?[n];
            foreach (var (i, j) in edges) { frame[i] ??= X(axes[i], axes[j]); frame[j] ??= X(axes[i], axes[j]); }
            bool OnE(int l, (long, long, long) m) => X(m, frame[l]!.Value) == (0, 0, 0);
            var x = new int?[n]; x[0] = 0;
            for (bool grew = true; grew;)
            {
                grew = false;
                foreach (var (i, j) in edges)
                {
                    var m = X(axes[i], axes[j]);
                    int t = OnE(i, m) ^ OnE(j, m) ? 1 : 0;
                    if (x[i] is { } xi && x[j] is { } xj) Assert.Equal(t, xi ^ xj);
                    else if (x[i] is { } a) { x[j] = a ^ t; grew = true; }
                    else if (x[j] is { } b) { x[i] = b ^ t; grew = true; }
                }
            }
            Assert.Contains(Enumerable.Range(0, n), l => edges.Count(e => e.Item1 == l || e.Item2 == l) >= 3);
            Assert.Equal(0, TiltedAxes.ColourDimension(n, axes));
            Assert.True(TiltedAxes.ColourLawProven(n, edges, axes));
            Assert.Equal(0, TiltedAxes.BondNullity(n, edges, weights, axes));
        }

        var k4 = TiltedAxes.Graphs(4).Single(g => g.Name == "K4").Edges;
        Assert.True(TiltedAxes.ColourLawProven(4, k4, generic.Take(4).ToArray()));
        Assert.True(TiltedAxes.ColourLawProven(4, k4, new (long, long, long)[] { (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 0) }));   // letters
        Assert.False(TiltedAxes.ColourLawProven(4, new[] { (0, 1), (2, 3) }, generic.Take(4).ToArray()));                         // disconnected

        // the admitted colour: the plane's normal when every field lies along it, none when one leaves it or the axes
        // span space; under one axis the fields' common direction if it is perpendicular to the axis
        var none = new (long, long, long)[5];
        Assert.Equal((4L, -3L, -2L), TiltedAxes.AdmittedColour(coplanar, none));
        Assert.Equal((4L, -3L, -2L), TiltedAxes.AdmittedColour(coplanar, new (long, long, long)[] { (-8, 6, 4), (0, 0, 0), (4, -3, -2), (0, 0, 0), (0, 0, 0) }));
        Assert.Null(TiltedAxes.AdmittedColour(coplanar, new (long, long, long)[] { (1, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0) }));
        Assert.Null(TiltedAxes.AdmittedColour(generic, none));
        var line = Enumerable.Repeat((1L, 1L, 0L), 3).ToArray();
        Assert.Equal((1L, -1L, 0L), TiltedAxes.AdmittedColour(line, new (long, long, long)[] { (1, -1, 0), (-2, 2, 0), (0, 0, 0) }));
        Assert.Null(TiltedAxes.AdmittedColour(line, new (long, long, long)[] { (1, -1, 0), (0, 0, 1), (0, 0, 0) }));
        Assert.Null(TiltedAxes.AdmittedColour(line, new (long, long, long)[] { (1, 0, 0), (0, 0, 0), (0, 0, 0) }));
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
