using MirrorWorld;
using static MirrorWorld.LetterCube;

namespace MirrorWorldTests;

// From-below guard for the theodolite (built 2026-10-20 beside experiments/THE_CUBE_SITE_BY_SITE.md and its Python twin
// simulations/framework/diagnostics/cube_letters.py). A reading with rational turns is judged by a route that never
// sees the letters: the generator built row-major from the Lindblad formula, times four in Gaussian integers, permuted
// to the per-site order, must satisfy E·L - (T(L) - c)·E = 0 exactly. A reading that needs a root (x^2 = 3, the golden
// ratio) is pinned as equations: consistent with the right one, inconsistent with a wrong one. Silent rows have a reason
// to be silent: nothing site-local exists, or the spectrum does not pair at every N.
public class TheodoliteTests
{
    static readonly World W = new();
    static readonly Gauss[,] Zm = LetterMatrix('Z');
    static readonly Gauss[,] Xm = LetterMatrix('X');

    [Fact]
    public void The_Theodolite_Owns_The_Letters_The_Reading_And_Every_N_And_Inherits_The_Cube()
    {
        var theodolite = new Theodolite(new LetterCube(W));
        Assert.Equal(new[] { "letters", "reading", "every N" }, theodolite.Own);
        Assert.Equal(new[] { "coordinates", "height", "step", "arrow", "moves", "x", "y", "z" }, theodolite.Inherited);
    }

    // ---- the binomial systems: equations, never numbers ----

    static Binomial B(string monomial, long num, long den = 1, long numIm = 0, long denIm = 0)
    {
        var e = new Dictionary<string, int>();
        foreach (var part in monomial.Split('·'))
        {
            var bits = part.Split('^');
            e[bits[0]] = bits.Length > 1 ? int.Parse(bits[1]) : 1;
        }
        return new Binomial(e, GaussBig.Of(num, numIm), GaussBig.Of(den, denIm), "test");
    }

    [Fact]
    public void A_Root_Stays_An_Equation()
    {
        // x/y = -3 and x·y = -1: the elimination leaves y^2 = 1/3, never a square root of 3
        var red = Binomials.Solve(new[] { B("x·y^-1", -3), B("x·y", -1) }, new[] { "x", "y" });
        Assert.NotNull(red);
        var y = red!.Alone("y");
        Assert.NotNull(y);
        Assert.Equal(2, y!.Value.K);
        Assert.Equal(y.Value.Num * GaussBig.Of(3), y.Value.Den);            // y^2 = 1/3
        Assert.Null(Binomials.RationalSolution(red));                          // a root is needed
    }

    [Fact]
    public void A_Cycle_Closes_Only_When_Its_Constants_Multiply_To_One()
    {
        var ring = new[] { B("x·y^-1", 0, 1, 1), B("y·z^-1", 0, 1, 1) };   // x/y = i, y/z = i
        Assert.NotNull(Binomials.Solve(ring.Append(B("z·x^-1", -1)), new[] { "x", "y", "z" }));
        Assert.Null(Binomials.Solve(ring.Append(B("z·x^-1", 1)), new[] { "x", "y", "z" }));
        Assert.Null(Binomials.Solve(new[] { B("x", 2), B("x", 3) }, new[] { "x" }));
    }

    // ---- the independent route: the full generator, row-major, from the Lindblad formula ----

    static Gauss G(GaussBig g) => new((long)g.Re, (long)g.Im);

    static Gauss[,] Factor(Theodolite.Kind k, (GaussBig Num, GaussBig Den) t)
        => k == Theodolite.Kind.Diagonal
            ? new[,] { { G(t.Den), Gauss.Zero }, { Gauss.Zero, G(t.Num) } }
            : new[,] { { Gauss.Zero, G(t.Num) }, { G(t.Den), Gauss.Zero } };

    static Gauss[,] Kron(Gauss[,] a, Gauss[,] b)
    {
        int ra = a.GetLength(0), rb = b.GetLength(0);
        var m = new Gauss[ra * rb, ra * rb];
        for (int i = 0; i < ra; i++) for (int j = 0; j < ra; j++)
            for (int k = 0; k < rb; k++) for (int l = 0; l < rb; l++)
                m[i * rb + k, j * rb + l] = a[i, j] * b[k, l];
        return m;
    }

    static Gauss[,] Mul(Gauss[,] a, Gauss[,] b)
    {
        int n = a.GetLength(0);
        var c = new Gauss[n, n];
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++)
            {
                var aik = a[i, k];
                if (aik.IsZero) continue;
                for (int j = 0; j < n; j++) c[i, j] = c[i, j] + aik * b[k, j];
            }
        return c;
    }

    static Gauss[,] Conjugate(Gauss[,] a)
    {
        int n = a.GetLength(0);
        var c = new Gauss[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) c[i, j] = a[i, j].Conj;
        return c;
    }

    static Gauss[,] Dagger(Gauss[,] a)
    {
        int n = a.GetLength(0);
        var c = new Gauss[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) c[i, j] = a[j, i].Conj;
        return c;
    }

    static Gauss[,] Embed(int n, int site, Gauss[,] op)
    {
        var m = new Gauss[1, 1] { { Gauss.One } };
        for (int l = 0; l < n; l++) m = Kron(m, l == site ? op : LetterMatrix('I'));
        return m;
    }

    // 4L in the per-site order of the theodolite (site 0 first, cells |i><j| per site)
    static Gauss[,] GeneratorTimesFour(Theodolite.Setting s)
    {
        int n = s.N, d = 1 << n, dd = d * d;
        var h = new Gauss[d, d];
        foreach (var (t, c) in s.H)
        {
            var m = new Gauss[1, 1] { { Gauss.One } };
            for (int l = 0; l < n; l++) m = Kron(m, LetterMatrix(t.Letter(l)));
            for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) h[i, j] = h[i, j] + c * m[i, j];
        }
        var L = new Gauss[dd, dd];
        for (int i = 0; i < d; i++) for (int j = 0; j < d; j++)
            for (int k = 0; k < d; k++) for (int l = 0; l < d; l++)
            {
                // -i [H, rho] times four: -4i (H ⊗ 1 - 1 ⊗ H^T)
                var v = Gauss.Zero;
                if (j == l) v = v + h[i, k];
                if (i == k) v = v - h[l, j];
                if (!v.IsZero) L[i * d + j, k * d + l] = L[i * d + j, k * d + l] + new Gauss(0, -4) * v;
            }
        foreach (var (site, jump, rate) in s.Jumps)
        {
            var J = Embed(n, site, jump);
            var jbar = Conjugate(J);
            var jdj = Mul(Dagger(J), J);
            for (int i = 0; i < d; i++) for (int j = 0; j < d; j++)
                for (int k = 0; k < d; k++) for (int l = 0; l < d; l++)
                {
                    // rate (J ⊗ conj J - (J^dagger J ⊗ 1 + 1 ⊗ (J^dagger J)^T)/2), times four
                    var v = 4 * (J[i, k] * jbar[j, l]);
                    if (j == l) v = v - 2 * jdj[i, k];
                    if (i == k) v = v - 2 * jdj[l, j];
                    if (!v.IsZero) L[i * d + j, k * d + l] = L[i * d + j, k * d + l] + rate * v;
                }
        }
        var perm = new int[dd];
        for (int i = 0; i < d; i++) for (int j = 0; j < d; j++)
        {
            int ps = 0;
            for (int l = 0; l < n; l++)
            {
                int bi = (i >> (n - 1 - l)) & 1, bj = (j >> (n - 1 - l)) & 1;
                ps = ps * 4 + 2 * bi + bj;
            }
            perm[i * d + j] = ps;
        }
        var P = new Gauss[dd, dd];
        for (int r = 0; r < dd; r++) for (int c = 0; c < dd; c++) P[perm[r], perm[c]] = L[r, c];
        return P;
    }

    static bool Reflects(Theodolite.Setting s, Theodolite.Relation rel, bool transposed,
        IReadOnlyList<(Theodolite.Kind A, Theodolite.Kind B)> kinds,
        IReadOnlyList<(GaussBig, GaussBig)> x, IReadOnlyList<(GaussBig, GaussBig)> y, IReadOnlyList<int>? relabelling = null)
    {
        int n = s.N;
        var E = new Gauss[1, 1] { { Gauss.One } };
        int[] perm = { 0, 2, 1, 3 };
        for (int l = 0; l < n; l++)
        {
            var a = Factor(kinds[l].A, x[l]);
            var b = Factor(kinds[l].B, y[l]);
            var el = Kron(a, new[,] { { b[0, 0], b[1, 0] }, { b[0, 1], b[1, 1] } });   // rho -> a rho b: kron(a, b^T)
            if (transposed)
            {
                var t = new Gauss[4, 4];
                for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) t[r, c] = el[r, perm[c]];
                el = t;
            }
            E = Kron(E, el);
        }
        if (relabelling is { } pi)
        {
            // E·P_pi, P_pi moving site l's cell (one base-4 digit, site 0 the most significant) to site pi[l]
            int cells = 1 << (2 * n);
            var moved = new Gauss[cells, cells];
            for (int c = 0; c < cells; c++)
            {
                int image = 0;
                for (int l = 0; l < n; l++)
                {
                    int digit = (c >> (2 * (n - 1 - l))) & 3;
                    image |= digit << (2 * (n - 1 - pi[l]));
                }
                for (int r = 0; r < cells; r++) moved[r, c] = E[r, image];
            }
            E = moved;
        }
        return Satisfies(s, rel, E);
    }

    // E·L = (T(L) - c)·E on the full generator, exactly (a copy: E·L = L·E)
    static bool Satisfies(Theodolite.Setting s, Theodolite.Relation rel, Gauss[,] E)
    {
        var L = GeneratorTimesFour(s);
        int dd = L.GetLength(0);
        var tr = Gauss.Zero;
        for (int i = 0; i < dd; i++) tr = tr + L[i, i];
        long den = 2L * (dd / 4);
        Assert.True(tr.Im == 0 && tr.Re % den == 0, $"4c must be an integer, tr(4L) = {tr}");
        long fourC = -tr.Re / den;                                          // c = -tr(L) / (2·4^(n-1))
        var lhs = Mul(E, L);
        if (rel == Theodolite.Relation.Copy)
        {
            var le = Mul(L, E);                                             // a copy: E·L = L·E, no constant
            for (int i = 0; i < dd; i++)
                for (int j = 0; j < dd; j++)
                    if (lhs[i, j] != le[i, j]) return false;
            return true;
        }
        var T = rel == Theodolite.Relation.Pi ? L : Dagger(L);
        var rhsPart = Mul(T, E);
        for (int i = 0; i < dd; i++)
            for (int j = 0; j < dd; j++)
            {
                // E·L - (T(L) - c)·E with T(L) = -L or -L^dagger
                var v = lhs[i, j] + rhsPart[i, j] + fourC * E[i, j];
                if (!v.IsZero) return false;
            }
        return true;
    }

    // The sector reading's element, densely, at its rational point: for each pair of control bits (sigma, tau) of the
    // input cell the controls go to their kept or flipped bits, and every other site carries A_l.sigma and B_l.tau
    // (transposed: A on the bra's bits, B on the ket's), so A = sum_sigma |pi(sigma)><sigma| ⊗ A_sigma as a matrix.
    static bool ReflectsBySectors(Theodolite.Setting s, Theodolite.Relation rel, Theodolite.Reading r, bool unconditioned = false)
    {
        var sol = Binomials.RationalSolution(r.System);
        Assert.NotNull(sol);
        var controls = r.Controls!.ToList();
        int n = s.N, k = controls.Count, cells = 1 << (2 * n);
        int[] perm = { 0, 2, 1, 3 };
        // unconditioned: every sector takes sector 0's turns, the product the conditioning replaces
        string Bits(int v) => new string(Enumerable.Range(0, k).Select(i => (char)('0' + (unconditioned ? 0 : (v >> i) & 1))).ToArray());
        var E = new Gauss[cells, cells];
        for (int sg = 0; sg < 1 << k; sg++)
            for (int tg = 0; tg < 1 << k; tg++)
            {
                int a = r.Transposed ? tg : sg, b = r.Transposed ? sg : tg;
                var part = new Gauss[1, 1] { { Gauss.One } };
                for (int l = 0; l < n; l++)
                {
                    var el = new Gauss[4, 4];
                    int ci = controls.IndexOf(l);
                    if (ci >= 0)
                    {
                        int ks = (sg >> ci) & 1, kt = (tg >> ci) & 1;
                        int fa = r.Kinds[l].A == Theodolite.Kind.AntiDiagonal ? 1 : 0, fb = r.Kinds[l].B == Theodolite.Kind.AntiDiagonal ? 1 : 0;
                        int outKet = (r.Transposed ? kt : ks) ^ fa, outBra = (r.Transposed ? ks : kt) ^ fb;
                        el[2 * outKet + outBra, 2 * ks + kt] = Gauss.One;
                    }
                    else
                    {
                        var fA = Factor(r.Kinds[l].A, sol![$"x{l}.{Bits(a)}"]);
                        var fB = Factor(r.Kinds[l].B, sol![$"y{l}.{Bits(b)}"]);
                        el = Kron(fA, new[,] { { fB[0, 0], fB[1, 0] }, { fB[0, 1], fB[1, 1] } });
                        if (r.Transposed)
                        {
                            var t = new Gauss[4, 4];
                            for (int row = 0; row < 4; row++) for (int c = 0; c < 4; c++) t[row, c] = el[row, perm[c]];
                            el = t;
                        }
                    }
                    part = Kron(part, el);
                }
                for (int i = 0; i < cells; i++) for (int j = 0; j < cells; j++) E[i, j] = E[i, j] + part[i, j];
            }
        return Satisfies(s, rel, E);
    }

    static bool ReflectsAtRationalPoint(Theodolite.Setting s, Theodolite.Relation rel, Theodolite.Reading r)
    {
        var sol = Binomials.RationalSolution(r.System);
        Assert.NotNull(sol);
        var x = Enumerable.Range(0, s.N).Select(l => sol![$"x{l}"]).ToList();
        var y = Enumerable.Range(0, s.N).Select(l => sol![$"y{l}"]).ToList();
        return Reflects(s, rel, r.Transposed, r.Kinds, x, y, r.Relabelling);
    }

    // ---- settings ----

    static IEnumerable<(string, long)> Bond(int n, int a, int b, long j, long dm, long dz)
    {
        string Two(char p, char q)
        {
            var c = Enumerable.Repeat('I', n).ToArray();
            c[a] = p;
            c[b] = q;
            return new string(c);
        }
        if (j != 0) { yield return (Two('X', 'X'), j); yield return (Two('Y', 'Y'), j); }
        if (dm != 0) { yield return (Two('X', 'Y'), dm); yield return (Two('Y', 'X'), -dm); }
        if (dz != 0) yield return (Two('Z', 'Z'), dz);
    }

    static (string, long) Field(int n, int l, char p, long c)
    {
        var s = Enumerable.Repeat('I', n).ToArray();
        s[l] = p;
        return (new string(s), c);
    }

    static IEnumerable<(string, long)> Heisenberg(int n) => Enumerable.Range(0, n - 1).SelectMany(a => Bond(n, a, a + 1, 1, 0, 1));

    static Gauss[,] Combo(long cx, long cy)        // cx·X + cy·Y
        => new[,] { { Gauss.Zero, new Gauss(cx, -cy) }, { new Gauss(cx, cy), Gauss.Zero } };

    // ---- reading settings ----

    [Fact]
    public void It_Reads_Pi_Under_Z_Dephasing()
    {
        var s = Theodolite.Setting.Of(3, Heisenberg(3), new[] { (0, Zm, 7L), (1, Zm, 10L), (2, Zm, 12L) });
        var r = Theodolite.Read(s, Theodolite.Relation.Pi);
        Assert.NotNull(r);
        Assert.True(r!.Transposed);
        Assert.All(r.Kinds, k => Assert.Equal((Theodolite.Kind.Diagonal, Theodolite.Kind.AntiDiagonal), k));
        Assert.True(ReflectsAtRationalPoint(s, Theodolite.Relation.Pi, r));
    }

    [Fact]
    public void The_Thermal_Weight_Is_The_Steady_State_On_Three_Unaligned_Axes()
    {
        // the flux-pi/2 triangle (hoppings 1, 2, 3i), Delta = 7/10, axes (3X+4Y)/5, X, (5X-12Y)/13 at 1/5, 1/4, 1/3, cooling
        // 1 and heating 1/10, 1/16, 0: every steady-state ratio is 1/4. Each site scaled to integers (a jump (3X+4Y)/5 at
        // rate 1/5 is 3X+4Y at 1/125), every bond by ten; each equation is homogeneous, so the readings stay.
        var h = Bond(3, 0, 1, 10, 0, 7).Concat(Bond(3, 1, 2, 20, 0, 7)).Concat(Bond(3, 2, 0, 0, 30, 7));
        var jumps = new (int, Gauss[,], long)[]
        {
            (0, Lowering, 250), (0, Raising, 25), (0, Combo(3, 4), 2),
            (1, Lowering, 16), (1, Raising, 1), (1, Xm, 4),
            (2, Lowering, 507), (2, Combo(5, -12), 1),
        };
        var s = Theodolite.Setting.Of(3, h, jumps);
        var r = Theodolite.Read(s, Theodolite.Relation.Pi);
        Assert.NotNull(r);
        Assert.True(r!.Transposed);
        Assert.All(r.Kinds, k => Assert.Equal((Theodolite.Kind.AntiDiagonal, Theodolite.Kind.AntiDiagonal), k));
        Assert.True(ReflectsAtRationalPoint(s, Theodolite.Relation.Pi, r));
        // the control: heating on site 0 raised by 1/10 (25 -> 50 at the site's scale) leaves the letters silent
        var raised = Theodolite.Setting.Of(3, h, jumps.Select(j => j.Item1 == 0 && j.Item3 == 25 ? (0, Raising, 50L) : j));
        Assert.Null(Theodolite.Read(raised, Theodolite.Relation.Pi));
    }

    [Fact]
    public void With_X_Dephasing_The_Turn_Is_A_Root_Kept_As_An_Equation()
    {
        // the Heisenberg chain, cooling 1 and X dephasing 1/2 (scaled by two): x/y = -3 and x·y = -1, so x^2 = 3
        var s = Theodolite.Setting.Of(3, Heisenberg(3), Enumerable.Range(0, 3).SelectMany(l => new[] { (l, Lowering, 2L), (l, Xm, 1L) }));
        var r = Theodolite.Read(s, Theodolite.Relation.Pi);
        Assert.NotNull(r);
        Assert.True(r!.Transposed);
        var unknowns = r.System.Unknowns;
        for (int l = 0; l < 3; l++)
        {
            Assert.NotNull(Binomials.Solve(r.Equations.Append(B($"x{l}^2", 3)), unknowns));
            Assert.Null(Binomials.Solve(r.Equations.Append(B($"x{l}^2", 2)), unknowns));
            Assert.NotNull(Binomials.Solve(r.Equations.Append(B($"x{l}·y{l}", -1)), unknowns));
        }
    }

    [Fact]
    public void The_Gradient_Is_Silent_For_Pi_And_Reads_S_XY_For_The_Reflection()
    {
        // cooling 1/2, 1, 3/2 and heating 1/8, 1/3, 1/5 (each site scaled)
        var jumps = new (int, Gauss[,], long)[]
        {
            (0, Lowering, 4), (0, Raising, 1), (1, Lowering, 3), (1, Raising, 1), (2, Lowering, 15), (2, Raising, 2),
        };
        var s = Theodolite.Setting.Of(3, Heisenberg(3), jumps);
        Assert.Null(Theodolite.Read(s, Theodolite.Relation.Pi));
        var r = Theodolite.Read(s, Theodolite.Relation.Reflection);
        Assert.NotNull(r);
        Assert.False(r!.Transposed);                                         // S_XY's kinds: plain, anti-diagonal on both sides
        Assert.All(r.Kinds, k => Assert.Equal((Theodolite.Kind.AntiDiagonal, Theodolite.Kind.AntiDiagonal), k));
        Assert.True(ReflectsAtRationalPoint(s, Theodolite.Relation.Reflection, r));
    }

    static Theodolite.Setting FieldTriangle(long delta, long j2, long d2)
    {
        // F137's triangle: unit hopping, real on (0, 1) and (1, 2), J + iD on (2, 0); fields 3/10, -1/5, 1/2; cooling
        // 1/2, 1, 3/2 at one ratio 1/4 (site 0 times 40, site 1 times 20, site 2 times 8; bonds times 10)
        var h = Bond(3, 0, 1, 10, 0, delta).Concat(Bond(3, 1, 2, 10, 0, delta)).Concat(Bond(3, 2, 0, 10 * j2, 10 * d2, delta))
            .Append(Field(3, 0, 'Z', 12)).Append(Field(3, 1, 'Z', -4)).Append(Field(3, 2, 'Z', 4));
        var jumps = new (int, Gauss[,], long)[]
        {
            (0, Lowering, 20), (0, Raising, 5), (1, Lowering, 20), (1, Raising, 5), (2, Lowering, 12), (2, Raising, 3),
        };
        return Theodolite.Setting.Of(3, h, jumps);
    }

    [Fact]
    public void Fields_Read_F137s_Chain_As_One_Map_And_The_Controls_Stay_Silent()
    {
        var s = FieldTriangle(0, 0, 1);
        var r = Theodolite.Read(s, Theodolite.Relation.Pi);
        Assert.NotNull(r);
        Assert.False(r!.Transposed);
        Assert.All(r.Kinds, k => Assert.Equal((Theodolite.Kind.AntiDiagonal, Theodolite.Kind.AntiDiagonal), k));
        for (int l = 0; l < 3; l++)                                          // the spin flip with the ratio of the arrows
            Assert.NotNull(Binomials.Solve(r.Equations.Append(B($"x{l}·y{l}^-1", -4)), r.System.Unknowns));
        Assert.True(ReflectsAtRationalPoint(s, Theodolite.Relation.Pi, r));
        Assert.Null(Theodolite.Read(FieldTriangle(7, 0, 1), Theodolite.Relation.Pi));
        Assert.Null(Theodolite.Read(FieldTriangle(0, 1, 0), Theodolite.Relation.Pi));
    }

    [Fact]
    public void Identity_Letters_Inside_A_Term_Do_Not_Change_The_Reading()
    {
        var jumps = new (int, Gauss[,], long)[] { (0, Zm, 1), (1, Zm, 1) };
        var asFields = Theodolite.Setting.Of(2, new[] { ("XI", 1L), ("YI", 1L) }, jumps);
        foreach (var rel in new[] { Theodolite.Relation.Pi, Theodolite.Relation.Reflection })
        {
            var r = Theodolite.Read(asFields, rel);
            Assert.NotNull(r);
            Assert.True(ReflectsAtRationalPoint(asFields, rel, r!));
        }
    }

    static IEnumerable<(string, long)> XyChain(int n, long dz = 0) => Enumerable.Range(0, n - 1).SelectMany(a => Bond(n, a, a + 1, 1, 0, dz));

    [Fact]
    public void A_Copy_Reads_The_Sublattice_Transpose_On_The_XY_Chain()
    {
        // U(rho) = Sigma·rho^T·Sigma, Sigma the sign on every second site: the involution
        // docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md splits the block (1, 1) by on the XY chain; it needs no ZZ
        long[] rates = { 3, 5, 7, 4 };
        foreach (int n in new[] { 3, 4 })
        {
            var jumps = Enumerable.Range(0, n).Select(l => (l, Zm, rates[l])).ToList();
            var s = Theodolite.Setting.Of(n, XyChain(n), jumps);
            var r = Theodolite.Read(s, Theodolite.Relation.Copy, branch: true);
            Assert.NotNull(r);
            Assert.All(r!.Kinds, k => Assert.Equal((Theodolite.Kind.Diagonal, Theodolite.Kind.Diagonal), k));
            var sol = Binomials.RationalSolution(r.System);
            Assert.NotNull(sol);
            for (int l = 0; l + 1 < n; l++)                                  // the sign alternates, on both factors
            {
                Assert.Equal(sol![$"x{l + 1}"].Num * sol[$"x{l}"].Den, -(sol[$"x{l}"].Num * sol[$"x{l + 1}"].Den));
                Assert.Equal(sol[$"y{l + 1}"].Num * sol[$"y{l}"].Den, -(sol[$"y{l}"].Num * sol[$"y{l + 1}"].Den));
            }
            Assert.True(ReflectsAtRationalPoint(s, Theodolite.Relation.Copy, r));
            var heis = Theodolite.Setting.Of(n, XyChain(n, 1), jumps);
            Assert.Null(Theodolite.Read(heis, Theodolite.Relation.Copy, branch: true));
            var x = Enumerable.Range(0, n).Select(l => sol![$"x{l}"]).ToList();
            var y = Enumerable.Range(0, n).Select(l => sol![$"y{l}"]).ToList();
            Assert.False(Reflects(heis, Theodolite.Relation.Copy, true, r.Kinds, x, y));   // ZZ breaks it on the full generator
        }
    }

    [Fact]
    public void A_Copy_Leaves_Out_The_Gauge_And_Reads_The_Turn_By_X()
    {
        // plain: the all-diagonal copies hold the identity and are left out; what is read is the turn by X^N, a copy
        // wherever H commutes with X^N, the move carrying (1, 1) onto (N - 1, N - 1)
        var s = Theodolite.Setting.Of(3, Heisenberg(3), new[] { (0, Zm, 1L), (1, Zm, 1L), (2, Zm, 1L) });
        var r = Theodolite.Read(s, Theodolite.Relation.Copy, branch: false);
        Assert.NotNull(r);
        Assert.All(r!.Kinds, k => Assert.Equal((Theodolite.Kind.AntiDiagonal, Theodolite.Kind.AntiDiagonal), k));
        Assert.True(ReflectsAtRationalPoint(s, Theodolite.Relation.Copy, r));
        // silent outside the plain all-diagonal group, which is not read: the turns about Z (and Z^N's moves) live in
        // it, so null says no copy of the other kinds, not no symmetry
        var damped = Theodolite.Setting.Of(3, Heisenberg(3), new[] { (0, Lowering, 1L), (1, Lowering, 1L), (2, Lowering, 1L) });
        Assert.Null(Theodolite.Read(damped, Theodolite.Relation.Copy));
    }

    // ---- every N ----

    static Theodolite.Cell Cell(IEnumerable<(Gauss[,], long)> jumps, params (string, int[], long)[] windows)
        => new(jumps, null, windows);

    static readonly (string, int[], long)[] HeisCell = { ("XX", new[] { 0, 1 }, 1), ("YY", new[] { 0, 1 }, 1), ("ZZ", new[] { 0, 1 }, 1) };

    [Fact]
    public void Every_N_Uniform_Rows()
    {
        foreach (var cell in new[]
                 {
                     Cell(new[] { (Zm, 1L) }, HeisCell),
                     Cell(new[] { (Lowering, 2L), (Xm, 1L) }, HeisCell),
                     Cell(new[] { (Lowering, 4L), (Raising, 1L) }, HeisCell),
                 })
        {
            var t = Theodolite.EveryN(cell, Theodolite.Relation.Pi);
            Assert.NotNull(t);
            Assert.Equal(Theodolite.Shape.Uniform, t!.Shape);
        }
    }

    [Fact]
    public void The_Golden_Bond_Turns_By_The_Phase_Of_Its_Pair_Amplitude()
    {
        // XX + XY + YX: the sigma^+ sigma^+ amplitude is 1 + 2i, and the uniform element's turn is its phase:
        // x^2 = (1 + 2i)^2 / 5 = (-3 + 4i)/5, the golden ratio being tan of half that angle
        var t = Theodolite.EveryN(Cell(new[] { (Zm, 1L) }, ("XX", new[] { 0, 1 }, 1), ("XY", new[] { 0, 1 }, 1), ("YX", new[] { 0, 1 }, 1)),
            Theodolite.Relation.Pi);
        Assert.NotNull(t);
        Assert.Equal(Theodolite.Shape.Uniform, t!.Shape);
        var eqs = t.System.Lines().ToList();
        var rows = t.System;
        var withPhase = Binomials.Solve(Pinned(rows).Append(B("X0^2", -3, 5, 4)), rows.Unknowns);
        var withOne = Binomials.Solve(Pinned(rows).Append(B("X0^2", 1)), rows.Unknowns);
        Assert.True(withPhase is not null || Binomials.Solve(Pinned(rows).Append(B("Y0^2", -3, 5, 4)), rows.Unknowns) is not null,
            string.Join("; ", eqs));
        Assert.True(withOne is null || Binomials.Solve(Pinned(rows).Append(B("Y0^2", 1)), rows.Unknowns) is null);
    }

    // the reduced rows of a system, back as binomials
    static IEnumerable<Binomial> Pinned(Binomials.Reduced red)
    {
        for (int r = 0; r < red.Exponents.Count; r++)
        {
            var e = new Dictionary<string, int>();
            for (int j = 0; j < red.Unknowns.Count; j++)
                if (red.Exponents[r][j] != 0) e[red.Unknowns[j]] = red.Exponents[r][j];
            yield return new Binomial(e, red.Constants[r].Num, red.Constants[r].Den, "row");
        }
    }

    [Fact]
    public void The_C0_Window_Routes_With_One_Map_On_Every_Site()
    {
        var cell = Cell(new[] { (Zm, 1L) }, ("XZY", new[] { 0, 1, 2 }, 1), ("YZX", new[] { 0, 1, 2 }, 1));
        var t = Theodolite.EveryN(cell, Theodolite.Relation.Pi);
        Assert.NotNull(t);
        Assert.Equal(Theodolite.Shape.Uniform, t!.Shape);
        var sol = Binomials.RationalSolution(t.System);
        Assert.NotNull(sol);
        var s = cell.Chain(4);
        var x = Enumerable.Repeat(sol!["X0"], 4).ToList();
        var y = Enumerable.Repeat(sol["Y0"], 4).ToList();
        Assert.True(Reflects(s, Theodolite.Relation.Pi, t.Transposed, Enumerable.Repeat(t.Kinds[0], 4).ToList(), x, y));
    }

    [Fact]
    public void Heisenberg_Plus_DM_Turns_By_Arctan_D_Over_J_And_Pure_DM_By_A_Quarter()
    {
        // D/J = 1/2: J = 2, D = 1, Delta = 1 (scaled); the turn z has z^2 = (2 + i)^2/5 = (3 + 4i)/5
        var dm = Cell(new[] { (Zm, 1L) }, ("XX", new[] { 0, 1 }, 2), ("YY", new[] { 0, 1 }, 2), ("ZZ", new[] { 0, 1 }, 2),
            ("XY", new[] { 0, 1 }, 1), ("YX", new[] { 0, 1 }, -1));
        var t = Theodolite.EveryN(dm, Theodolite.Relation.Pi);
        Assert.NotNull(t);
        Assert.Equal(Theodolite.Shape.Turned, t!.Shape);
        Assert.NotNull(t.Turn);
        var (k, num, den) = t.Turn!.Value;
        Assert.Equal(2, k);
        Assert.Equal(num * GaussBig.Of(5), den * GaussBig.Of(3, 4));
        var pure = Cell(new[] { (Zm, 1L) }, ("XY", new[] { 0, 1 }, 1), ("YX", new[] { 0, 1 }, -1));
        var q = Theodolite.EveryN(pure, Theodolite.Relation.Pi);
        Assert.NotNull(q);
        Assert.Equal(Theodolite.Shape.Turned, q!.Shape);
        var (k2, n2, d2) = q.Turn!.Value;
        Assert.Equal(2, k2);
        Assert.Equal(n2, -d2);                                               // z^2 = -1
    }

    [Fact]
    public void The_Hard_Rows_Are_Silent_For_Every_N()
    {
        foreach (var (t1, t2) in new[] { ("XX", "XY"), ("XY", "XZ") })
            Assert.Null(Theodolite.EveryN(Cell(new[] { (Zm, 1L) }, (t1, new[] { 0, 1 }, 1), (t2, new[] { 0, 1 }, 1)), Theodolite.Relation.Pi));
    }

    [Fact]
    public void A_Copy_Turns_The_XY_Chain_By_A_Quarter_From_Site_To_Site()
    {
        // Sigma on every second site: a transposed sandwich's factor steps by z^-2, so the sign is a quarter turn, z^2 = -1
        var xy = Cell(new[] { (Zm, 1L) }, ("XX", new[] { 0, 1 }, 1), ("YY", new[] { 0, 1 }, 1));
        var t = Theodolite.EveryN(xy, Theodolite.Relation.Copy, branch: true);
        Assert.NotNull(t);
        Assert.Equal(Theodolite.Shape.Turned, t!.Shape);
        var (k, num, den) = t.Turn!.Value;
        Assert.Equal(2, k);
        Assert.Equal(num, -den);
        Assert.Equal(Theodolite.Shape.Uniform, Theodolite.EveryN(xy, Theodolite.Relation.Copy)!.Shape);   // plain first: X^N
        Assert.Null(Theodolite.EveryN(Cell(new[] { (Zm, 1L) }, HeisCell), Theodolite.Relation.Copy, branch: true));
    }

    [Fact]
    public void Rate_Free_Rows()
    {
        var thermal = new[] { (Lowering, 4L), (Raising, 1L) };
        var heis = Theodolite.EveryN(Cell(thermal, HeisCell), Theodolite.Relation.Reflection, rateFree: true);
        Assert.Equal(Theodolite.Shape.Uniform, heis!.Shape);
        var dm = Cell(thermal, ("XX", new[] { 0, 1 }, 4), ("YY", new[] { 0, 1 }, 4), ("XY", new[] { 0, 1 }, 3), ("YX", new[] { 0, 1 }, -3));
        var t = Theodolite.EveryN(dm, Theodolite.Relation.Reflection, rateFree: true);
        Assert.Equal(Theodolite.Shape.Turned, t!.Shape);
    }

    [Fact]
    public void The_Golden_Window_Repeats_With_Period_Four()
    {
        var golden = Cell(new[] { (Zm, 1L) }, ("XZX", new[] { 0, 1, 2 }, 1), ("XZY", new[] { 0, 1, 2 }, 1), ("YZX", new[] { 0, 1, 2 }, 1));
        var t = Theodolite.EveryN(golden, Theodolite.Relation.Pi);
        Assert.NotNull(t);
        Assert.Equal(Theodolite.Shape.Periodic, t!.Shape);
        Assert.Equal(4, t.Period);
    }

    // ---- the closed traverse: every ring ----

    // (name, cell, relation, plain moduli, transposed moduli); an empty list: no route closes any ring
    static IEnumerable<(string, Theodolite.Cell, Theodolite.Relation, int[], int[])> Rings()
    {
        var z = new[] { (Zm, 1L) };
        var ratio = new[] { (Lowering, 4L), (Raising, 1L) };
        var field = new[] { ('Z', 1L) };
        int[] bond = { 0, 1 }, window = { 0, 1, 2 };
        Theodolite.Cell On(IEnumerable<(Gauss[,], long)> j, (char, long)[]? f, int[] offsets, params (string, long)[] terms)
            => new(j, f, terms.Select(t => (t.Item1, offsets, t.Item2)));
        int[] none = Array.Empty<int>();
        yield return ("XY, Z dephasing", On(z, null, bond, ("XX", 1), ("YY", 1)), Theodolite.Relation.Pi, new[] { 2 }, new[] { 1 });
        yield return ("Heisenberg, Z dephasing", On(z, null, bond, ("XX", 1), ("YY", 1), ("ZZ", 1)), Theodolite.Relation.Pi, none, new[] { 1 });
        yield return ("pure DM, Z dephasing", On(z, null, bond, ("XY", 1), ("YX", -1)), Theodolite.Relation.Pi, new[] { 2 }, new[] { 2 });
        yield return ("twisted 3 + 4i, Z dephasing", On(z, null, bond, ("XX", 3), ("YY", 3), ("XY", 4), ("YX", -4)), Theodolite.Relation.Pi, none, none);
        yield return ("XY, one ratio and a field", On(ratio, field, bond, ("XX", 1), ("YY", 1)), Theodolite.Relation.Pi, new[] { 2 }, none);
        yield return ("pure DM, one ratio and a field", On(ratio, field, bond, ("XY", 1), ("YX", -1)), Theodolite.Relation.Pi, new[] { 1 }, none);
        yield return ("XY, Z dephasing, a copy", On(z, null, bond, ("XX", 1), ("YY", 1)), Theodolite.Relation.Copy, new[] { 1 }, new[] { 2 });
        yield return ("the golden bond", On(z, null, bond, ("XX", 1), ("XY", 1), ("YX", 1)), Theodolite.Relation.Pi, new[] { 2 }, new[] { 1 });
        yield return ("the hard row XX + XY", On(z, null, bond, ("XX", 1), ("XY", 1)), Theodolite.Relation.Pi, none, none);
        yield return ("the c = 0 window", On(z, null, window, ("XZY", 1), ("YZX", 1)), Theodolite.Relation.Pi, new[] { 1 }, new[] { 1 });
        yield return ("the golden window", On(z, null, window, ("XZX", 1), ("XZY", 1), ("YZX", 1)), Theodolite.Relation.Pi, new[] { 4 }, new[] { 4 });
        yield return ("Heisenberg + DM, D/J = 1/2", On(z, null, bond, ("XX", 2), ("YY", 2), ("ZZ", 2), ("XY", 1), ("YX", -1)), Theodolite.Relation.Pi, none, none);
        // a twist of an eighth turn per bond, (1 + i)^2/2 = i: four rounds before the turns come back
        yield return ("XX + YY + XY - YX, Z dephasing", On(z, null, bond, ("XX", 1), ("YY", 1), ("XY", 1), ("YX", -1)), Theodolite.Relation.Pi, new[] { 4 }, new[] { 4 });
    }

    // The closed traverse against the rings themselves: Read on a real ring of N sites (every site's kind chosen
    // afresh, no step or period assumed) closes exactly where the traverse says, N up to 8. Pure DM under one ratio and
    // a field closes every ring while under Z dephasing only the even ones: F137's e^(2i·Phi) = (-1)^N, met by the
    // frustration the ring's length brings.
    [Fact]
    public void The_Closed_Traverse_Reads_Every_Ring_The_Ring_Itself_Confirms()
    {
        foreach (var (name, cell, rel, plain, transposed) in Rings())
            foreach (var (branch, moduli) in new[] { (false, plain), (true, transposed) })
            {
                var t = Theodolite.EveryN(cell, rel, branch: branch, closed: true);
                if (moduli.Length == 0) Assert.Null(t);
                else
                {
                    Assert.NotNull(t);
                    Assert.Equal(Theodolite.Shape.Closed, t!.Shape);
                    Assert.Equal(moduli, t.Closures!.Select(c => c.Modulus).ToArray());
                }
                for (int n = cell.Normalized().SmallestRing; n <= 8; n++)
                {
                    bool ring = Theodolite.Read(cell.Normalized().Ring(n), rel, branch) is not null;
                    Assert.True(ring == (t is not null && t.Closes(n)), $"{name}, {(branch ? "transposed" : "plain")}, N = {n}: the ring reads {ring}");
                }
            }
    }

    // The bound from the step lattice against the closure tried at every k: a route closes after k periods exactly
    // when its modulus over its period divides k, k up to 40.
    [Fact]
    public void A_Route_Closes_On_The_Multiples_Of_Its_Modulus_And_Nowhere_Else()
    {
        foreach (var (name, cell, rel, _, _) in Rings())
        {
            var t = Theodolite.EveryN(cell, rel, closed: true);
            if (t is null) continue;
            foreach (var c in t.Closures!)
            {
                int periods = c.Modulus / c.Period;
                for (int k = 1; k <= 40; k++)
                {
                    bool closes = Binomials.Solve(Pinned(c.System).Append(B($"r^{k}", 1)).Append(B($"s^{k}", 1)), c.System.Unknowns) is not null;
                    Assert.True(closes == (k % periods == 0), $"{name}, period {c.Period}, k = {k}");
                }
            }
        }
    }

    // The golden window's ring: a period of two sites whose turns change sign from one period to the next, r = -1,
    // so the turns repeat after four sites (F116's router [a, a, b, b] on the ring, docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md
    // §8); r = 1 does not close.
    [Fact]
    public void The_Golden_Window_Closes_Every_Fourth_Ring_By_A_Sign_Each_Period()
    {
        var golden = Cell(new[] { (Zm, 1L) }, ("XZX", new[] { 0, 1, 2 }, 1), ("XZY", new[] { 0, 1, 2 }, 1), ("YZX", new[] { 0, 1, 2 }, 1));
        var t = Theodolite.EveryN(golden, Theodolite.Relation.Pi, closed: true);
        Assert.NotNull(t);
        var c = Assert.Single(t!.Closures!);
        Assert.Equal(4, c.Modulus);
        Assert.Equal(2, c.Period);
        Assert.NotNull(Binomials.Solve(Pinned(c.System).Append(B("r", -1)), c.System.Unknowns));
        Assert.Null(Binomials.Solve(Pinned(c.System).Append(B("r", 1)), c.System.Unknowns));
        Assert.True(t.Closes(4) && t.Closes(8) && t.Closes(400) && !t.Closes(6) && !t.Closes(10));
    }

    // The ring the census reads is the physical ring: on the XY ring of four sites the reading, built densely at its
    // rational point, carries the full generator, the wrap bond included, onto -L - c.
    [Fact]
    public void The_Ring_Is_The_Ring()
    {
        var xy = Cell(new[] { (Zm, 1L) }, ("XX", new[] { 0, 1 }, 1), ("YY", new[] { 0, 1 }, 1));
        var ring = xy.Ring(4);
        Assert.Equal(8, ring.H.Count);
        var r = Theodolite.Read(ring, Theodolite.Relation.Pi, false);
        Assert.NotNull(r);
        Assert.True(ReflectsAtRationalPoint(ring, Theodolite.Relation.Pi, r!));
        Assert.Null(Theodolite.Read(xy.Ring(3), Theodolite.Relation.Pi, false));
        Assert.Throws<ArgumentOutOfRangeException>(() => xy.Ring(2));
        // a period whose chain would pass twelve sites is not read: the periods stop at (13 - Span)/2
        Assert.Equal(new[] { 2 }, Theodolite.EveryN(xy, Theodolite.Relation.Pi, maxPeriod: 9, branch: false, closed: true)!.Closures!.Select(c => c.Modulus));
    }

    // A route's own element, built densely on the ring from its turns, x_l = P_(l mod p)·r^(l div p), at a rational
    // point of its system closed after n/p periods, carries the full generator of the four-site ring onto -L - c: the
    // XY ring (a sign per site), pure DM (a sign on one side) and the eighth twist (a quarter turn per site, r = ±i).
    // With the step dropped (r = s = 1) the same element fails.
    [Fact]
    public void A_Route_Built_On_The_Ring_Carries_The_Generator()
    {
        var z = new[] { (Zm, 1L) };
        foreach (var cell in new[]
                 {
                     Cell(z, ("XX", new[] { 0, 1 }, 1), ("YY", new[] { 0, 1 }, 1)),
                     Cell(z, ("XY", new[] { 0, 1 }, 1), ("YX", new[] { 0, 1 }, -1)),
                     Cell(z, ("XX", new[] { 0, 1 }, 1), ("YY", new[] { 0, 1 }, 1), ("XY", new[] { 0, 1 }, 1), ("YX", new[] { 0, 1 }, -1)),
                 })
        {
            var c = Theodolite.EveryN(cell, Theodolite.Relation.Pi, branch: false, closed: true)!.Closures![0];
            Assert.True(RouteReflects(cell, Theodolite.Relation.Pi, c, 4, stepped: true));
            Assert.False(RouteReflects(cell, Theodolite.Relation.Pi, c, 4, stepped: false));
        }
    }

    static bool RouteReflects(Theodolite.Cell cell, Theodolite.Relation rel, Theodolite.Closure c, int n, bool stepped)
    {
        int k = n / c.Period;
        var red = Binomials.Solve(Pinned(c.System).Append(B($"r^{k}", 1)).Append(B($"s^{k}", 1)), c.System.Unknowns);
        Assert.NotNull(red);
        var sol = Binomials.RationalSolution(red!);
        Assert.NotNull(sol);
        (GaussBig, GaussBig) Turn(string cls, string step, int l)
        {
            var (num, den) = sol![$"{cls}{l % c.Period}"];
            var (sn, sd) = sol[step];
            for (int b = 0; stepped && b < l / c.Period; b++) { num *= sn; den *= sd; }
            return (num, den);
        }
        var x = Enumerable.Range(0, n).Select(l => Turn("P", "r", l)).ToList();
        var y = Enumerable.Range(0, n).Select(l => Turn("Q", "s", l)).ToList();
        var kinds = Enumerable.Range(0, n).Select(l => c.Kinds[l % c.Period]).ToList();
        return Reflects(cell.Ring(n), rel, c.Transposed, kinds, x, y);
    }

    // ---- relabelling the sites ----

    // Under loss alone a bond with ZZ and s, a, w all nonzero pairs exactly at rates the chain's reversal keeps, and the
    // carrier is the reversal times a uniform turn, then S_XY (experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md, "The
    // reversal"). No site-local sandwich reads it; after the reversal one does, and the element E·P_pi, built densely at
    // the rational point of the reading, carries the full generator onto -L^dagger - c (the second route). At rates the
    // reversal moves no relabelling reads, and the reading at the kept rates does not carry the moved setting.
    [Fact]
    public void The_Reversal_Carries_The_Cooled_Chain()
    {
        var J = new[] { ("XX", -2L), ("XY", -2L), ("YX", 1L), ("YY", -2L), ("ZZ", -2L) };   // s = −2, a = −3/2, |w|² = 1/4
        Theodolite.Setting Chain(params long[] rates) => Theodolite.Setting.Of(3,
            J.SelectMany(t => new[] { (t.Item1 + "I", t.Item2), ("I" + t.Item1, t.Item2) }),
            rates.Select((g, l) => (l, Lowering, g)));
        var kept = Chain(2, 7, 2);
        Assert.Null(Theodolite.Read(kept, Theodolite.Relation.Reflection));
        Assert.Null(Theodolite.Read(kept, Theodolite.Relation.Pi));
        Assert.Equal(new[] { new[] { 2, 1, 0 } }, Theodolite.Automorphisms(kept));
        var r = Theodolite.ReadRelabelled(kept, Theodolite.Relation.Reflection);
        Assert.NotNull(r);
        Assert.Equal(new[] { 2, 1, 0 }, r!.Relabelling);
        Assert.True(ReflectsAtRationalPoint(kept, Theodolite.Relation.Reflection, r));
        var moved = Chain(2, 7, 3);
        Assert.Null(Theodolite.ReadRelabelled(moved, Theodolite.Relation.Reflection));
        Assert.Null(Theodolite.ReadRelabelled(moved, Theodolite.Relation.Pi));
        var sol = Binomials.RationalSolution(r.System)!;
        var x = Enumerable.Range(0, 3).Select(l => sol[$"x{l}"]).ToList();
        var y = Enumerable.Range(0, 3).Select(l => sol[$"y{l}"]).ToList();
        Assert.False(Reflects(moved, Theodolite.Relation.Reflection, r.Transposed, r.Kinds, x, y, r.Relabelling));
        Assert.False(Reflects(kept, Theodolite.Relation.Reflection, r.Transposed, r.Kinds, x, y));   // without the reversal
    }

    // Between two settings, F137's exchange of loss and gain: X on every site carries L(H, loss) onto L(X^N·H·X^N, gain),
    // a copy with c = 0 (the traces agree). Read off the letters, and the element checked densely: E·L(from) = L(to)·E.
    // Gain against the same H does not read when H does not commute with X^N (an X field and a ZZ bond here), and the
    // element of the true pair does not carry the false one.
    [Fact]
    public void Between_Two_Settings_It_Reads_Loss_Against_Gain()
    {
        var h = new[] { ("XYI", 2L), ("IYX", -1L), ("ZZI", 3L), ("IZZ", 1L), ("XII", 1L), ("IIY", 2L) };
        long Sign(string t) => t.Count(c => c is 'Y' or 'Z') % 2 == 0 ? 1 : -1;      // X^N·P·X^N = ±P
        var loss = Theodolite.Setting.Of(3, h, new[] { (0, Lowering, 2L), (1, Lowering, 1L), (2, Lowering, 3L) });
        var gain = Theodolite.Setting.Of(3, h.Select(t => (t.Item1, Sign(t.Item1) * t.Item2)),
            new[] { (0, Raising, 2L), (1, Raising, 1L), (2, Raising, 3L) });
        var r = Theodolite.ReadBetween(loss, gain, Theodolite.Relation.Copy);
        Assert.NotNull(r);
        var sol = Binomials.RationalSolution(r!.System);
        Assert.NotNull(sol);
        var x = Enumerable.Range(0, 3).Select(l => sol![$"x{l}"]).ToList();
        var y = Enumerable.Range(0, 3).Select(l => sol![$"y{l}"]).ToList();
        Assert.True(CopiesBetween(loss, gain, r.Transposed, r.Kinds, x, y));
        var sameH = Theodolite.Setting.Of(3, h, new[] { (0, Raising, 2L), (1, Raising, 1L), (2, Raising, 3L) });
        Assert.False(CopiesBetween(loss, sameH, r.Transposed, r.Kinds, x, y));
        Assert.Null(Theodolite.ReadBetween(loss, sameH, Theodolite.Relation.Copy));
    }

    static bool CopiesBetween(Theodolite.Setting from, Theodolite.Setting to, bool transposed,
        IReadOnlyList<(Theodolite.Kind A, Theodolite.Kind B)> kinds,
        IReadOnlyList<(GaussBig, GaussBig)> x, IReadOnlyList<(GaussBig, GaussBig)> y)
    {
        int n = from.N;
        var E = new Gauss[1, 1] { { Gauss.One } };
        int[] perm = { 0, 2, 1, 3 };
        for (int l = 0; l < n; l++)
        {
            var a = Factor(kinds[l].A, x[l]);
            var b = Factor(kinds[l].B, y[l]);
            var el = Kron(a, new[,] { { b[0, 0], b[1, 0] }, { b[0, 1], b[1, 1] } });
            if (transposed)
            {
                var t = new Gauss[4, 4];
                for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) t[r, c] = el[r, perm[c]];
                el = t;
            }
            E = Kron(E, el);
        }
        var lf = GeneratorTimesFour(from);
        var lt = GeneratorTimesFour(to);
        int dd = lf.GetLength(0);
        Gauss Trace(Gauss[,] m) { var t = Gauss.Zero; for (int i = 0; i < dd; i++) t = t + m[i, i]; return t; }
        if (Trace(lf) != Trace(lt)) return false;
        var left = Mul(E, lf);
        var right = Mul(lt, E);
        for (int i = 0; i < dd; i++)
            for (int j = 0; j < dd; j++)
                if (left[i, j] != right[i, j]) return false;
        return true;
    }

    // The odd-N paragraph of F170's proof, PROOF_RING_FROZEN_BAND §2 after Lemma 1 (U·h₊·U = −h₋; F171 §7 writes it
    // Σh_fΣ = −h_{−f}, its sites counted from 1, here from 0), read on the theodolite's own road. Σ = Z on the odd sites
    // turns the sign of every bond but the wrap bond, which at odd N joins two even sites, so the plain Π that closes the
    // even XY ring under Z dephasing onto itself, X^N on one side and Σ on both, carries the odd ring onto the spin ring
    // with its wrap bond reversed, J_w → −J_w (RingFrozenBandWitness.ReversedSpinRing): E·L(J_w) = (−L(−J_w) − c)·E, read
    // between the two settings at N = 3 and 5, checked densely against both full generators and failing on the ring
    // alone. The raw equations are the stagger itself, x_l·x_{l+1}^{-1} = −1 on both sides for the N − 1 bonds of the
    // chain and +1 for the wrap pair, the sign the reversed ring supplies, so the reduced system is the chain's (2(N − 1)
    // rows, two free turns, the rational point x_l = y_l = (−1)^l). Transposed, Σ·ρᵀ·Σ is a COPY between the odd ring and
    // its reversed ring, the chain's U of PROOF_DEPHASING_FRONT_RENEWAL between two settings, the same raw rows, checked
    // densely: the two rings are isospectral, not only mirror images; the transposed reflection reaches the reversed ring
    // too. Behind all six pairs of relation and branch one parity law: an XY bond's row is x_l/x_{l+1} = ε·J_to/J_from
    // (its two rows agree, J_to = ±J_from here) with ε = −1 for the plain Π, the transposed copy and the transposed
    // reflection (the Π's target and the transposition each turn the sign of the whole Hamiltonian part) and +1 for the
    // other three, and round a ring the ratios multiply to one, so the ring itself is reached when ε^N = 1 and the
    // reversed ring when ε^N = −1: at odd N by exactly those three pairs, at even N by none, every pair then closing the
    // ring onto itself (the plain copy with X^N on both sides). All twelve readings at each of N = 3, 4 and 5 are asserted.
    [Fact]
    public void At_Odd_N_The_Stagger_Carries_The_XY_Ring_Onto_The_Ring_With_Its_Wrap_Bond_Reversed()
    {
        var xy = Cell(new[] { (Zm, 1L) }, ("XX", new[] { 0, 1 }, 1), ("YY", new[] { 0, 1 }, 1));
        var pairs = new[]
        {
            (Theodolite.Relation.Pi, false), (Theodolite.Relation.Pi, true), (Theodolite.Relation.Copy, false),
            (Theodolite.Relation.Copy, true), (Theodolite.Relation.Reflection, false), (Theodolite.Relation.Reflection, true),
        };
        foreach (int n in new[] { 3, 4, 5 })
        {
            var ring = xy.Ring(n);
            var reversed = WrapReversed(ring);
            foreach (var (rel, transposed) in pairs)
            {
                bool minus = (rel == Theodolite.Relation.Pi) != transposed;          // epsilon = -1
                Assert.Equal(!minus || n % 2 == 0, Theodolite.Read(ring, rel, transposed) is not null);
                Assert.Equal(minus && n % 2 == 1, Theodolite.ReadBetween(ring, reversed, rel, transposed) is not null);
            }
            var plainCopy = Theodolite.Read(ring, Theodolite.Relation.Copy, false);
            Assert.All(plainCopy!.Kinds, k => Assert.Equal((Theodolite.Kind.AntiDiagonal, Theodolite.Kind.AntiDiagonal), k));
            if (n % 2 == 0) continue;

            // the plain Pi onto the reversed ring: the stagger on both sides, X^N on one
            var r = Theodolite.ReadBetween(ring, reversed, Theodolite.Relation.Pi, false)!;
            bool antiOnB = r.Kinds[0].B == Theodolite.Kind.AntiDiagonal;
            Assert.All(r.Kinds, k => Assert.True((k.A == Theodolite.Kind.AntiDiagonal) != (k.B == Theodolite.Kind.AntiDiagonal)
                                                 && (k.B == Theodolite.Kind.AntiDiagonal) == antiOnB));
            StaggerRows(r.Equations, n);
            Assert.Equal(2 * (n - 1), r.System.Pivots.Count);
            Assert.Equal(2, r.System.Free.Count);
            var sol = Binomials.RationalSolution(r.System);
            Assert.NotNull(sol);
            var x = Enumerable.Range(0, n).Select(l => sol![$"x{l}"]).ToList();
            var y = Enumerable.Range(0, n).Select(l => sol![$"y{l}"]).ToList();
            for (int l = 0; l < n; l++)
            {
                Assert.Equal((GaussBig.Of(l % 2 == 0 ? 1 : -1), GaussBig.Of(1)), x[l]);
                Assert.Equal((GaussBig.Of(l % 2 == 0 ? 1 : -1), GaussBig.Of(1)), y[l]);
            }
            Assert.True(ReflectsBetween(ring, reversed, Theodolite.Relation.Pi, r.Transposed, r.Kinds, x, y));
            Assert.False(ReflectsBetween(ring, ring, Theodolite.Relation.Pi, r.Transposed, r.Kinds, x, y));

            // the transposed copy between the two rings: Sigma on both sides, the same rows
            var copy = Theodolite.ReadBetween(ring, reversed, Theodolite.Relation.Copy, true)!;
            Assert.All(copy.Kinds, k => Assert.Equal((Theodolite.Kind.Diagonal, Theodolite.Kind.Diagonal), k));
            StaggerRows(copy.Equations, n);
            var csol = Binomials.RationalSolution(copy.System);
            Assert.NotNull(csol);
            var cx = Enumerable.Range(0, n).Select(l => csol![$"x{l}"]).ToList();
            var cy = Enumerable.Range(0, n).Select(l => csol![$"y{l}"]).ToList();
            for (int l = 0; l < n; l++)
            {
                Assert.Equal((GaussBig.Of(l % 2 == 0 ? 1 : -1), GaussBig.Of(1)), cx[l]);
                Assert.Equal((GaussBig.Of(l % 2 == 0 ? 1 : -1), GaussBig.Of(1)), cy[l]);
            }
            Assert.True(CopiesBetween(ring, reversed, copy.Transposed, copy.Kinds, cx, cy));
            Assert.False(CopiesBetween(ring, ring, copy.Transposed, copy.Kinds, cx, cy));

            // the transposed reflection between the two rings, checked densely as well
            var refl = Theodolite.ReadBetween(ring, reversed, Theodolite.Relation.Reflection, true)!;
            Assert.All(refl.Kinds, k => Assert.True((k.A == Theodolite.Kind.AntiDiagonal) != (k.B == Theodolite.Kind.AntiDiagonal)));
            StaggerRows(refl.Equations, n);
            var fsol = Binomials.RationalSolution(refl.System);
            Assert.NotNull(fsol);
            var fx = Enumerable.Range(0, n).Select(l => fsol![$"x{l}"]).ToList();
            var fy = Enumerable.Range(0, n).Select(l => fsol![$"y{l}"]).ToList();
            Assert.True(ReflectsBetween(ring, reversed, Theodolite.Relation.Reflection, refl.Transposed, refl.Kinds, fx, fy));
            Assert.False(ReflectsBetween(ring, ring, Theodolite.Relation.Reflection, refl.Transposed, refl.Kinds, fx, fy));
        }
    }

    // the stagger's raw rows on a ring of n sites: -1 on both sides along the n - 1 bonds of the chain, +1 on the wrap pair
    static void StaggerRows(IReadOnlyList<Binomial> equations, int n)
    {
        for (int l = 0; l < n - 1; l++)
        {
            Assert.True(HasRow(equations, $"x{l}", $"x{l + 1}", -1), $"x{l}/x{l + 1} = -1 at N = {n}");
            Assert.True(HasRow(equations, $"y{l}", $"y{l + 1}", -1), $"y{l}/y{l + 1} = -1 at N = {n}");
        }
        Assert.True(HasRow(equations, $"x{n - 1}", "x0", 1), $"the wrap pair x{n - 1}/x0 = +1 at N = {n}");
        Assert.True(HasRow(equations, $"y{n - 1}", "y0", 1), $"the wrap pair y{n - 1}/y0 = +1 at N = {n}");
    }

    // a raw equation on exactly the two unknowns a and b, exponents +1 and -1 either way, equal to ratio (+1 or -1)
    static bool HasRow(IReadOnlyList<Binomial> equations, string a, string b, int ratio)
        => equations.Any(e =>
        {
            var on = e.Exponents.Where(kv => kv.Value != 0).ToDictionary(kv => kv.Key, kv => kv.Value);
            return on.Count == 2 && on.ContainsKey(a) && on.ContainsKey(b) && Math.Abs(on[a]) == 1 && on[a] == -on[b]
                   && (ratio == 1 ? e.Num == e.Den : e.Num == -e.Den);
        });

    // the ring with its wrap bond reversed, J_w -> -J_w: the only terms touching both site 0 and site n - 1
    static Theodolite.Setting WrapReversed(Theodolite.Setting ring)
    {
        int n = ring.N;
        return new Theodolite.Setting(n,
            ring.H.Select(t => (t.Term, t.Term.Letter(0) != 'I' && t.Term.Letter(n - 1) != 'I' ? -t.Coefficient : t.Coefficient)),
            ring.Jumps);
    }

    // E·L(from) = (T(L(to)) - c)·E on the two full generators, exactly, E built from the kinds and turns as in Reflects;
    // the two traces must agree (the same c on both sides), a copy goes to CopiesBetween
    static bool ReflectsBetween(Theodolite.Setting from, Theodolite.Setting to, Theodolite.Relation rel, bool transposed,
        IReadOnlyList<(Theodolite.Kind A, Theodolite.Kind B)> kinds,
        IReadOnlyList<(GaussBig, GaussBig)> x, IReadOnlyList<(GaussBig, GaussBig)> y)
    {
        if (rel == Theodolite.Relation.Copy) throw new ArgumentException("a copy is checked by CopiesBetween", nameof(rel));
        int n = from.N;
        var E = new Gauss[1, 1] { { Gauss.One } };
        int[] perm = { 0, 2, 1, 3 };
        for (int l = 0; l < n; l++)
        {
            var a = Factor(kinds[l].A, x[l]);
            var b = Factor(kinds[l].B, y[l]);
            var el = Kron(a, new[,] { { b[0, 0], b[1, 0] }, { b[0, 1], b[1, 1] } });
            if (transposed)
            {
                var t = new Gauss[4, 4];
                for (int r = 0; r < 4; r++) for (int c = 0; c < 4; c++) t[r, c] = el[r, perm[c]];
                el = t;
            }
            E = Kron(E, el);
        }
        var lf = GeneratorTimesFour(from);
        var lt = GeneratorTimesFour(to);
        int dd = lf.GetLength(0);
        Gauss tr = Gauss.Zero, trFrom = Gauss.Zero;
        for (int i = 0; i < dd; i++) { tr = tr + lt[i, i]; trFrom = trFrom + lf[i, i]; }
        Assert.True(tr == trFrom, $"the two generators must share their trace, {trFrom} against {tr}");
        long den = 2L * (dd / 4);
        Assert.True(tr.Im == 0 && tr.Re % den == 0, $"4c must be an integer, tr(4L) = {tr}");
        long fourC = -tr.Re / den;                                          // c = -tr(L) / (2·4^(n-1))
        var left = Mul(E, lf);
        var right = Mul(rel == Theodolite.Relation.Pi ? lt : Dagger(lt), E);
        for (int i = 0; i < dd; i++)
            for (int j = 0; j < dd; j++)
                if (!(left[i, j] + right[i, j] + fourC * E[i, j]).IsZero) return false;
        return true;
    }

    // Conditioned on an unwatched bit: -3XY with fields -2X under Z dephasing at (0, 3). X is the only letter on the
    // unwatched site 0, so turned onto Z its bit is a number s and site 1 sees the field (-2, -3s), perpendicular to
    // its watching; no product reads it (plain, relabelled or turned), the sandwich conditioned on that bit does, at
    // a rational point, and the element built densely from it reflects the full turned generator, while the same
    // element with sector 0's turns in both sectors does not. The same for a Y control (-2YY + 2Y + 2YX under Z at
    // (0, 1)) and a Z one (-3Z + 3ZY under X at (0, 1), turned on site 1), which no product reads either. A Z field on
    // site 1 breaks the pairing (F158's word, an exact no) and nothing reads; 3ZY - 3ZX under X at (0, 3, 1) pairs
    // (F158), and its sandwich needs a phase coupling the two watched sites, outside the form.
    [Fact]
    public void Conditioned_On_An_Unwatched_Bit_It_Reads_What_No_Product_Reads()
    {
        var z = LetterMatrix('Z');
        var x = LetterMatrix('X');
        var s = Theodolite.Setting.Of(2, new[] { ("XY", -3L), ("XI", -2L), ("IX", -2L) }, new[] { (0, z, 0L), (1, z, 3L) });
        foreach (var rel in new[] { Theodolite.Relation.Reflection, Theodolite.Relation.Pi })
        {
            Assert.Null(Theodolite.Read(s, rel));
            Assert.Null(Theodolite.ReadRotated(s, rel));
            Assert.Null(Theodolite.ReadRelabelled(s, rel));
            var r = Theodolite.ReadSectors(s, rel);
            Assert.NotNull(r);
            Assert.Equal(new[] { 0 }, r!.Controls);
            Assert.Equal("XZ", r.Axes);
            Assert.True(ReflectsBySectors(s.Rotated(r.Axes!), rel, r));
            Assert.False(ReflectsBySectors(s.Rotated(r.Axes!), rel, r, unconditioned: true));
        }
        var others = new[]
        {
            Theodolite.Setting.Of(2, new[] { ("YY", -2L), ("YI", 2L), ("IY", 2L), ("YX", 2L) }, new[] { (0, z, 0L), (1, z, 1L) }),
            Theodolite.Setting.Of(2, new[] { ("ZI", -3L), ("IZ", -3L), ("ZY", 3L) }, new[] { (0, x, 0L), (1, x, 1L) }),
        };
        foreach (var o in others)
        {
            foreach (var q in new[] { Theodolite.Relation.Reflection, Theodolite.Relation.Pi })
                Assert.True(Theodolite.Read(o, q) is null && Theodolite.ReadRotated(o, q) is null && Theodolite.ReadRelabelled(o, q) is null);
            var r = Theodolite.ReadSectors(o, Theodolite.Relation.Reflection) ?? Theodolite.ReadSectors(o, Theodolite.Relation.Pi);
            Assert.NotNull(r);
            var rel = Theodolite.ReadSectors(o, Theodolite.Relation.Reflection) is null ? Theodolite.Relation.Pi : Theodolite.Relation.Reflection;
            Assert.True(ReflectsBySectors(r!.Axes is { } ax ? o.Rotated(ax) : o, rel, r));
        }
        var broken = Theodolite.Setting.Of(2, new[] { ("XY", -3L), ("XI", -2L), ("IX", -2L), ("IZ", 1L) }, new[] { (0, z, 0L), (1, z, 3L) });
        Assert.Null(Theodolite.ReadSectors(broken, Theodolite.Relation.Reflection));
        Assert.Null(Theodolite.ReadSectors(broken, Theodolite.Relation.Pi));
        Assert.Equal(Sextant.Palindrome.Obstructed, Sextant.Decide(new Sextant.Setting(2, broken.H, broken.Jumps)).Kind);
        Assert.False(ReflectsBySectors(broken.Rotated("XZ"), Theodolite.Relation.Reflection,
            Theodolite.ReadSectors(s, Theodolite.Relation.Reflection)!));             // the check can fail
        var watchedPair = Theodolite.Setting.Of(3, new[] { ("ZYI", 3L), ("IZY", 3L), ("ZXI", -3L), ("IZX", -3L) },
            new[] { (0, x, 0L), (1, x, 3L), (2, x, 1L) });
        Assert.Null(Theodolite.ReadSectors(watchedPair, Theodolite.Relation.Reflection));
        Assert.Null(Theodolite.ReadSectors(watchedPair, Theodolite.Relation.Pi));
        Assert.StartsWith("F158", Sextant.Decide(new Sextant.Setting(3, watchedPair.H, watchedPair.Jumps)).Detail);
    }

    [Fact]
    public void In_The_Quarter_Turn_Frame_The_Damping_Arrow_Lands_On_Y()
    {
        // The reality window (experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md §7, Read on the cube): a compiled
        // circuit puts its damping between gates, i.e. in a turned frame. Rotated("Y") carries the site by 1 - iX,
        // which is sqrt(2) times the quarter-turn about X (sx), the jump by M J M^dagger = 2 sx J sx^dagger, so the
        // transfer matrix is four times that of the sx frame at the same rate. The damping's one arrow follows the
        // carried Z axis: I -> Z in the logical frame (real), I -> X in the Hadamard frame (real), I -> -Y in the sx
        // frame, the odd letter, so the trace feeds n_Y-odd weight at the rate of the lights. Z dephasing is a height
        // in every one of these frames and moves nothing. The same numbers as simulations/reality_window_cube_frames.py
        // up to the arrow's sign: Rotated carries J to M J M^dagger, i.e. sx J sx^dagger, the inverse of the frame the
        // flight document writes (Ad_U^dagger P Ad_U), where the carried arrow is I -> +Y.
        var damping = Theodolite.Setting.Of(1, Array.Empty<(string, long)>(), new[] { (0, Lowering, 1L) });
        var logical = TransferTimesFour(Lowering, 1);
        Assert.Equal(new Gauss[] { new(0, 0), new(2, 0), new(2, 0), new(4, 0) }, HeightsTimesFour(logical));
        Assert.Equal(new[] { ('I', 'Z', new Gauss(4, 0)) }, Arrows(logical));
        var sxFrame = damping.Rotated("Y").Jumps[0];
        Assert.Equal(1L, sxFrame.Rate);
        var turned = TransferTimesFour(sxFrame.Jump, sxFrame.Rate);
        Assert.Equal(new Gauss[] { new(0, 0), new(8, 0), new(16, 0), new(8, 0) }, HeightsTimesFour(turned));
        Assert.Equal(new[] { ('I', 'Y', new Gauss(-16, 0)) }, Arrows(turned));
        var hFrame = damping.Rotated("X").Jumps[0];
        var hadamard = TransferTimesFour(hFrame.Jump, hFrame.Rate);
        Assert.Equal(new Gauss[] { new(0, 0), new(16, 0), new(8, 0), new(8, 0) }, HeightsTimesFour(hadamard));
        Assert.Equal(new[] { ('I', 'X', new Gauss(16, 0)) }, Arrows(hadamard));
        // the arrow is the only entry into or out of Y in the sx frame, and the rate of the lights is its rate:
        // the arrow's modulus equals the sum of the heights on X and on Z (the two lights, 8 each)
        Assert.Equal(16L, Math.Abs(Arrows(turned)[0].CoefficientTimesFour.Re));
        var dephasing = Theodolite.Setting.Of(1, Array.Empty<(string, long)>(), new[] { (0, Zm, 1L) });
        foreach (var axes in new[] { "Z", "Y", "X" })
        {
            var j = dephasing.Rotated(axes).Jumps[0];
            var t = TransferTimesFour(j.Jump, j.Rate);
            Assert.Empty(Arrows(t));
            var heights = HeightsTimesFour(t);
            Assert.Equal(new Gauss(0, 0), heights[0]);
            // a height charges the two letters anticommuting with the carried axis, equally; the carried axis is
            // Z, then -Y (sx), then X (Hadamard); Rotated puts every site on one scale (a Z site keeps its jump at
            // four times the rate), so the charge is 4 * 8 in each frame
            long charged = 32;
            char axis = axes == "Z" ? 'Z' : axes == "Y" ? 'Y' : 'X';
            foreach (char letter in "XYZ")
                Assert.Equal(letter == axis ? new Gauss(0, 0) : new Gauss(charged, 0), heights["IXYZ".IndexOf(letter)]);
        }
    }

    // In a turned frame: uncoupled sites with fields -2Y + Z under X dephasing (one site unwatched) pair site by site,
    // and no Z-basis monomial sandwich reads them; with X carried onto Z (the cube's h_zx) on the watched sites they
    // read. Two checks on the second route: the turned setting's generator is exactly 4·Ad_U·L·Ad_U^-1 (dense, with
    // M = √2·U integral, so 4·S·L·S' = 4^n·L_turned for S = ⊗(M ⊗ conj M), S' = ⊗(M^dagger ⊗ M^T)), and the element
    // read there reflects the turned generator at its rational point. A wrong turn (Y onto Z) does not read.
    [Fact]
    public void In_A_Turned_Frame_It_Reads_X_Dephasing()
    {
        var x = LetterMatrix('X');
        var s = Theodolite.Setting.Of(3, new[] { ("YII", -2L), ("IYI", -2L), ("IIY", -2L), ("ZII", 1L), ("IZI", 1L), ("IIZ", 1L) },
            new[] { (0, x, 2L), (1, x, 0L), (2, x, 3L) });
        Assert.Null(Theodolite.Read(s, Theodolite.Relation.Reflection));
        Assert.Null(Theodolite.Read(s, Theodolite.Relation.Pi));
        var r = Theodolite.ReadRotated(s, Theodolite.Relation.Reflection);
        Assert.NotNull(r);
        Assert.Equal('X', r!.Axes![0]);
        Assert.Equal('X', r.Axes[2]);
        var turned = s.Rotated(r.Axes);
        Assert.True(ReflectsAtRationalPoint(turned, Theodolite.Relation.Reflection, r));
        Assert.Null(Theodolite.Read(s.Rotated("YYY"), Theodolite.Relation.Reflection));
        // the turned generator, densely, on this setting and on one with every letter on turned sites, a three-site
        // term and jumps of four kinds (the X rows of the letter map and the signs across a support)
        var rich = Theodolite.Setting.Of(3, new[] { ("XXI", 2L), ("IYX", -1L), ("ZIX", 3L), ("XYZ", 1L), ("YYY", -2L), ("XII", 1L) },
            new[] { (0, Lowering, 2L), (1, Raising, 1L), (2, LetterMatrix('Y'), 3L), (1, LetterMatrix('X'), 1L) });
        foreach (var (setting, axes) in new[] { (s, r.Axes), (s, "YXZ"), (rich, "XYZ"), (rich, "YXX"), (rich, "XYY") })
        {
            var L = GeneratorTimesFour(setting);
            var Lt = GeneratorTimesFour(setting.Rotated(axes));
            var S = new Gauss[1, 1] { { Gauss.One } };
            var Sinv = new Gauss[1, 1] { { Gauss.One } };
            foreach (char a in axes)
            {
                if (a == 'Z') { S = Kron(S, Scaled(2, Identity(4))); Sinv = Kron(Sinv, Scaled(2, Identity(4))); continue; }
                var m = a == 'X'
                    ? new[,] { { Gauss.One, Gauss.One }, { Gauss.One, -Gauss.One } }                    // X + Z
                    : new[,] { { Gauss.One, new Gauss(0, -1) }, { new Gauss(0, -1), Gauss.One } };     // 1 - iX
                S = Kron(S, Kron(m, Conjugate(m)));                 // rho -> M rho M^dagger in cells: kron(M, conj M)
                Sinv = Kron(Sinv, Kron(Dagger(m), Transpose2(m)));  // rho -> M^dagger rho M: kron(M^dagger, M^T)
            }
            var lhs = Scaled(4, Mul(Mul(S, L), Sinv));
            var rhs = Scaled(1L << (2 * 3), Lt);
            int dd = L.GetLength(0);
            for (int i = 0; i < dd; i++)
                for (int j = 0; j < dd; j++)
                    Assert.Equal(rhs[i, j], lhs[i, j]);
        }
    }

    static Gauss[,] Identity(int d)
    {
        var m = new Gauss[d, d];
        for (int i = 0; i < d; i++) m[i, i] = Gauss.One;
        return m;
    }

    static Gauss[,] Scaled(long k, Gauss[,] a)
    {
        int d = a.GetLength(0);
        var m = new Gauss[d, d];
        for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) m[i, j] = k * a[i, j];
        return m;
    }

    static Gauss[,] Transpose2(Gauss[,] a)
    {
        int d = a.GetLength(0);
        var m = new Gauss[d, d];
        for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) m[i, j] = a[j, i];
        return m;
    }

    // The direction of a relabelling, on a 3-cycle (not an involution) and a setting no permutation keeps: the generator
    // of Relabelled(pi) is P_pi·L·P_pi^-1 with P_pi moving site l's cell to site pi[l], and not P_pi^-1·L·P_pi.
    [Fact]
    public void A_Relabelling_Moves_Site_L_To_Pi_Of_L()
    {
        var s = Theodolite.Setting.Of(3, new[] { ("XYI", 1L), ("IZX", 2L), ("YIZ", -1L), ("XII", 3L) },
            new[] { (0, Lowering, 1L), (1, Lowering, 2L), (2, Lowering, 5L) });
        int[] pi = { 1, 2, 0 };
        var L = GeneratorTimesFour(s);
        var moved = GeneratorTimesFour(s.Relabelled(pi));
        int cells = L.GetLength(0);
        int Image(int c, int[] p)
        {
            int image = 0;
            for (int l = 0; l < 3; l++) image |= ((c >> (2 * (2 - l))) & 3) << (2 * (2 - p[l]));
            return image;
        }
        int[] inverse = { 2, 0, 1 };
        bool Conjugates(int[] p)
        {
            for (int r = 0; r < cells; r++)
                for (int c = 0; c < cells; c++)
                    if (moved[Image(r, p), Image(c, p)] != L[r, c]) return false;
            return true;
        }
        Assert.True(Conjugates(pi));
        Assert.False(Conjugates(inverse));
    }
}
