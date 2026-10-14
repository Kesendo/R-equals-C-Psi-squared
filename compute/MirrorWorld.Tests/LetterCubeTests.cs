using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the letter cube (adopted 2026-10-14 from docs/THE_ONE_SQUARE.md §7; the Python gates
// simulations/diagonal_square_gate.py S8, S9(e)). The numbers are judged by second routes: letters counted one by
// one, the tetrahedron's inequalities, PauliString.Multiply's phases, dense 2^N matrices built letter by letter here,
// the letter map of R and MirrorGroup's R; where a check chains two of LetterCube's own readings it says so.
public class LetterCubeTests
{
    static readonly World W = new();

    [Fact]
    public void The_Cube_Owns_Coordinates_Height_And_Step_And_Inherits_The_Frame()
    {
        var cube = new LetterCube(W);
        Assert.Equal(new[] { "coordinates", "height", "step" }, cube.Own);
        Assert.Equal(new[] { "x", "y", "z" }, cube.Inherited);
    }

    // §7: the four letters are four corners of the unit cube, the even tetrahedron: k_Y = k_Z XOR k_X.
    [Fact]
    public void The_Four_Letters_Are_The_Even_Tetrahedron()
    {
        Assert.Equal(new LetterCube.Point(0, 0, 0), LetterCube.Coordinates("I"));
        Assert.Equal(new LetterCube.Point(1, 0, 1), LetterCube.Coordinates("X"));
        Assert.Equal(new LetterCube.Point(1, 1, 0), LetterCube.Coordinates("Y"));
        Assert.Equal(new LetterCube.Point(0, 1, 1), LetterCube.Coordinates("Z"));
        foreach (char a in "IXYZ")
        {
            var p = LetterCube.Coordinates(a.ToString());
            Assert.Equal(p.KZ ^ p.KX, p.KY);
        }
    }

    // The coordinates read off the masks equal the letter counts (n_X + n_Y, n_Y + n_Z, n_X + n_Z), and w is the
    // number of non-identity letters, for every string on three sites.
    [Fact]
    public void The_Masks_Give_The_Letter_Counts()
    {
        foreach (var s in LetterCube.Strings(3))
        {
            string t = s.ToString(3);
            int nx = t.Count(c => c == 'X'), ny = t.Count(c => c == 'Y'), nz = t.Count(c => c == 'Z');
            Assert.Equal(new LetterCube.Point(nx + ny, ny + nz, nx + nz), LetterCube.Coordinates(s));
            Assert.Equal(nx + ny + nz, LetterCube.Weight(s));
            Assert.Equal(nx + ny + nz, LetterCube.Coordinates(s).Weight);
        }
    }

    // §7, gate S8(c): the strings reach exactly the C(N+3, 3) even lattice points of the scaled tetrahedron,
    // and each point (n_X, n_Y, n_Z) holds the multinomial N!/(n_X! n_Y! n_Z! n_I!) strings.
    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)]
    public void The_Strings_Fill_The_Even_Lattice_Points_Of_The_Tetrahedron(int n)
    {
        var reached = LetterCube.Points(n);
        // the tetrahedron itself, from its four faces: a + b + c even, a + b + c <= 2N, each coordinate at most the
        // sum of the other two (the gate's definition, independent of the letter counts)
        var faces = new HashSet<LetterCube.Point>();
        for (int a = 0; a <= n; a++)
            for (int b = 0; b <= n; b++)
                for (int c = 0; c <= n; c++)
                    if ((a + b + c) % 2 == 0 && a + b + c <= 2 * n && a <= b + c && b <= a + c && c <= a + b)
                        faces.Add(new LetterCube.Point(a, b, c));
        Assert.True(faces.SetEquals(reached.Keys));
        Assert.True(faces.SetEquals(LetterCube.TetrahedronPoints(n)));
        Assert.Equal(LetterCube.TetrahedronCount(n), faces.Count);
        Assert.Equal(1L << (2 * n), reached.Values.Sum());
        foreach (var (p, count) in reached)
        {
            Assert.Equal(0, (p.KZ + p.KX + p.KY) % 2);
            int nx = (p.KZ + p.KY - p.KX) / 2, ny = (p.KZ + p.KX - p.KY) / 2, nz = (p.KX + p.KY - p.KZ) / 2;
            Assert.Equal(Multinomial(n, nx, ny, nz), count);
        }
    }

    // §7, gate S8(f): the centre (N/2, N/2, N/2) holds N!/((N/4)!)^4 strings when 4 | N and none otherwise.
    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)] [InlineData(8)]
    public void The_Centre_Holds_Strings_Only_When_Four_Divides_N(int n)
    {
        long counted = 0;
        if (n % 2 == 0)
            LetterCube.Points(n).TryGetValue(new LetterCube.Point(n / 2, n / 2, n / 2), out counted);
        Assert.Equal(LetterCube.CentreCount(n), (BigInteger)counted);
        if (n == 4) Assert.Equal(24, counted);
        if (n == 8) Assert.Equal(2520, counted);
    }

    // §7: the three dephasing diagonals are the three heights, Q_P(sigma) = sum_l P_l sigma P_l = (N - 2 k_P) sigma.
    // Second route: PauliString.Multiply's phases, site by site.
    [Fact]
    public void The_Three_Diagonals_Are_The_Three_Heights()
    {
        const int n = 3;
        foreach (var s in LetterCube.Strings(n))
        {
            var pt = LetterCube.Coordinates(s);
            foreach (var (letter, k) in new[] { ('Z', pt.KZ), ('X', pt.KX), ('Y', pt.KY) })
            {
                int sum = 0;
                for (int l = 0; l < n; l++)
                {
                    var site = PauliString.Parse(new string('I', l) + letter + new string('I', n - 1 - l));
                    var (left, p1) = PauliString.Multiply(site, s);
                    var (back, p2) = PauliString.Multiply(left, site);
                    Assert.Equal(s, back);
                    int phase = (p1 + p2) % 4;                 // P s P = i^phase s, phase 0 or 2
                    Assert.True(phase == 0 || phase == 2);
                    sum += phase == 0 ? 1 : -1;
                }
                Assert.Equal(n - 2 * k, sum);
            }
        }
    }

    // §7: the dissipator is a height. Dephasing with integer rates per site and letter sends every Pauli string to
    // minus its height times itself. Second route: the Lindblad dissipator built densely here,
    // sum_{l,P} gamma_{l,P} (P_l sigma P_l - sigma), on 2^N x 2^N matrices, compared to 0.0.
    [Fact]
    public void Pauli_Dephasing_Charges_Each_String_Its_Height()
    {
        const int n = 3;
        var rates = new (long Z, long X, long Y)[] { (1, 0, 2), (3, 1, 0), (0, 2, 5) };
        foreach (var s in LetterCube.Strings(n))
        {
            var sigma = DenseString(s.ToString(n));
            var lhs = new Complex[1 << n, 1 << n];
            for (int l = 0; l < n; l++)
                foreach (var (letter, g) in new[] { ('Z', rates[l].Z), ('X', rates[l].X), ('Y', rates[l].Y) })
                {
                    if (g == 0) continue;
                    var P = DenseString(new string('I', l) + letter + new string('I', n - 1 - l));
                    var sandwich = Mul(P, Mul(sigma, P));
                    for (int i = 0; i < 1 << n; i++)
                        for (int j = 0; j < 1 << n; j++)
                            lhs[i, j] += g * (sandwich[i, j] - sigma[i, j]);
                }
            long h = LetterCube.Height(s, rates);
            for (int i = 0; i < 1 << n; i++)
                for (int j = 0; j < 1 << n; j++)
                    Assert.True(lhs[i, j] + h * sigma[i, j] == Complex.Zero);
        }
        // the same rates on every site: the height is linear on the cube (chained: the per-site reading, checked
        // densely above, against the point reading)
        var uniform = new (long Z, long X, long Y)[] { (2, 1, 3), (2, 1, 3), (2, 1, 3) };
        foreach (var s in LetterCube.Strings(n))
            Assert.Equal(LetterCube.Height(s, uniform), LetterCube.Height(LetterCube.Coordinates(s), 2, 1, 3));
    }

    // §7, gate S8(e): the half-turn about P's axis pairs the four one-site rates exactly when gamma_P = 0, every
    // rate triple in {0..3}^3; isotropic rates (depolarizing) pair under none.
    [Fact]
    public void A_HalfTurn_Pairs_The_Rates_Exactly_When_Its_Letter_Is_Free()
    {
        for (long gz = 0; gz <= 3; gz++)
            for (long gx = 0; gx <= 3; gx++)
                for (long gy = 0; gy <= 3; gy++)
                {
                    var axes = LetterCube.PairingHalfTurns(gz, gx, gy);
                    Assert.Equal(gz == 0, axes.Contains('Z'));
                    Assert.Equal(gx == 0, axes.Contains('X'));
                    Assert.Equal(gy == 0, axes.Contains('Y'));
                }
        Assert.Empty(LetterCube.PairingHalfTurns(1, 1, 1));       // a sub-case of the sweep: depolarizing
    }

    // The half-turn about P's axis, built geometrically (flip the two other coordinates), sends I to P and the other
    // two letters to each other, and is an involution.
    [Fact]
    public void The_HalfTurn_About_P_Pairs_I_With_P()
    {
        foreach (char p in "ZXY")
        {
            Assert.Equal(p, LetterCube.HalfTurn('I', p));
            Assert.Equal('I', LetterCube.HalfTurn(p, p));
            var others = "XYZ".Where(c => c != p).ToArray();
            Assert.Equal(others[1], LetterCube.HalfTurn(others[0], p));
            foreach (char c in "IXYZ") Assert.Equal(c, LetterCube.HalfTurn(LetterCube.HalfTurn(c, p), p));
        }
    }

    // §8, gate S9(e): R (I <-> X, Y <-> Z) carries (k_Z, k_X, k_Y) to (N - k_Z, k_X, N - k_Y), the populations'
    // edge k_Z = 0 onto the far edge k_Z = N, N = 1..5. Second route: R's letter map applied letter by letter.
    [Fact]
    public void R_Carries_The_Quiet_Edge_Onto_The_Far_Edge()
    {
        var rLetter = new Dictionary<char, char> { ['I'] = 'X', ['X'] = 'I', ['Y'] = 'Z', ['Z'] = 'Y' };
        for (int n = 1; n <= 5; n++)
            foreach (var s in LetterCube.Strings(n))
            {
                var image = LetterCube.RightByXN(s, n);
                var byLetters = new string(s.ToString(n).Select(c => rLetter[c]).ToArray());
                Assert.Equal(byLetters, image.ToString(n));
                Assert.Equal(LetterCube.RAsPoint(LetterCube.Coordinates(s), n), LetterCube.Coordinates(image));
                // the two edges: k_Z = 0 exactly for the populations {I, Z}^N, k_Z = N exactly for {X, Y}^N
                string t = s.ToString(n);
                Assert.Equal(t.All(c => c == 'I' || c == 'Z'), LetterCube.Coordinates(s).KZ == 0);
                Assert.Equal(t.All(c => c == 'X' || c == 'Y'), LetterCube.Coordinates(s).KZ == n);
            }
        // docked: the letter map is MirrorGroup's R, rho -> rho * X^N, site by site up to phase
        foreach (char c in "IXYZ")
            Assert.Equal(rLetter[c], MirrorGroup.R.ApplySite(c).Letter);
    }

    // The guards: R and the per-site height refuse strings with letters beyond the sites they are given.
    [Fact]
    public void The_Guards_Refuse_Letters_Beyond_The_Sites()
    {
        Assert.Throws<ArgumentException>(() => LetterCube.RightByXN(PauliString.Parse("IIX"), 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => LetterCube.RightByXN(PauliString.Parse("X"), 33));
        Assert.Throws<ArgumentException>(() => LetterCube.Height(PauliString.Parse("IIZ"), new (long, long, long)[] { (1, 1, 1), (1, 1, 1) }));
    }

    // §7, gate S8(h), docked on the hardness bloc: the parity cube (bit_a, bit_b, y_par) is this cube's first two
    // coordinates mod 2 and half a sum, y_par = (k_Z + k_X - k_Y)/2 mod 2 = n_Y mod 2.
    [Fact]
    public void The_Parity_Cube_Is_Two_Coordinates_Mod_Two_And_Half_A_Sum()
    {
        foreach (var s in LetterCube.Strings(3))
            Assert.Equal(Hardness.Cube(s.ToString(3).ToCharArray()), LetterCube.Parities(LetterCube.Coordinates(s)));
    }

    // ---- the steps (§8, gate S9) ----

    static readonly LetterCube.Point O = new(0, 0, 0);

    // [T, sigma] is 0 when they commute and 2 i^k Target when they anticommute (T sigma = i^k Target, Target the
    // phase-free product): every placed term against every
    // string on two sites, judged by the dense commutator built letter by letter, compared exactly.
    [Fact]
    public void The_Commutator_Coefficient_Is_The_Dense_Commutator()
    {
        foreach (var t in LetterCube.Strings(2))
            foreach (var s in LetterCube.Strings(2))
            {
                var T = DenseString(t.ToString(2));
                var S = DenseString(s.ToString(2));
                var ts = Mul(T, S);
                var st = Mul(S, T);
                var c = LetterCube.CommutatorCoefficient(t, s);
                var expected = new Complex[4, 4];
                if (c is { } r)
                {
                    var target = DenseString(r.Target.ToString(2));
                    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) expected[i, j] = 2 * Phase(r.PhasePower) * target[i, j];
                }
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        Assert.True(ts[i, j] - st[i, j] == expected[i, j]);
            }
    }

    // §8, gate S9(a): every nonzero move of a two-site term P (x) Q changes w by exactly one; a bond P (x) P steps
    // +-2 along exactly one axis, never P's own, and ZZ moves only k_X and k_Y, XX only k_Z and k_Y, YY only k_Z and
    // k_X; a mixed term is not an axis step (X (x) Y sends Y (x) I to Z (x) Y, a step (0, 1, 1)).
    [Fact]
    public void A_Bond_Steps_To_The_Adjacent_Plane_And_P_P_Along_One_Axis()
    {
        foreach (char p in "XYZ")
            foreach (char q in "XYZ")
            {
                var term = PauliString.Parse($"{p}{q}");
                var axesSeen = new HashSet<char>();
                foreach (var s in LetterCube.Strings(2))
                {
                    if (LetterCube.Step(term, s) is not { } r) continue;
                    Assert.Equal(1, Math.Abs(LetterCube.Weight(r.Target) - LetterCube.Weight(s)));
                    if (p != q) continue;
                    var d = r.Step;
                    var nonzero = new[] { ('Z', d.KZ), ('X', d.KX), ('Y', d.KY) }.Where(a => a.Item2 != 0).ToList();
                    Assert.Single(nonzero);
                    Assert.Equal(2, Math.Abs(nonzero[0].Item2));
                    Assert.NotEqual(p, nonzero[0].Item1);
                    axesSeen.Add(nonzero[0].Item1);
                }
                if (p == q) Assert.Equal("XYZ".Where(c => c != p).ToHashSet(), axesSeen);
            }
        var mixed = LetterCube.Step(PauliString.Parse("XY"), PauliString.Parse("YI"));
        Assert.NotNull(mixed);
        Assert.Equal("ZY", mixed!.Value.Target.ToString(2));
        Assert.Equal(new LetterCube.Point(0, 1, 1), mixed.Value.Step);
    }

    // §8, gate S9(a): a one-site field P keeps w and steps (+-1, -+1) on the face diagonal of the plane k_P = const.
    [Fact]
    public void A_Field_Steps_On_A_Face_Diagonal()
    {
        foreach (char p in "XYZ")
        {
            int moved = 0;
            foreach (var s in LetterCube.Strings(1))
            {
                if (LetterCube.Step(PauliString.Parse(p.ToString()), s) is not { } r) continue;
                moved++;
                Assert.Equal(LetterCube.Weight(s), LetterCube.Weight(r.Target));
                var d = r.Step;
                int own = p switch { 'Z' => d.KZ, 'X' => d.KX, _ => d.KY };
                Assert.Equal(0, own);
                Assert.Equal(0, d.KZ + d.KX + d.KY);
                Assert.Equal(2, Math.Abs(d.KZ) + Math.Abs(d.KX) + Math.Abs(d.KY));
            }
            Assert.Equal(2, moved);                                  // the two letters other than I and P
        }
    }

    // §8, gate S9(b): on the chain at N = 3 every move of bonds of the XYZ family with no field is an axis step. The
    // moves are read term by term, so the coefficients never enter: the term sets {XX, YY, ZZ} (Heisenberg, XXZ, any
    // triple of bond weights), {XX, YY}, {ZZ} and {XX} cover every such chain. An X field adds face diagonals in the
    // plane k_X = const.
    [Fact]
    public void An_XYZ_Chain_Moves_Only_Along_The_Axes_And_A_Field_Adds_Face_Diagonals()
    {
        static bool Axis(LetterCube.Point d)
            => new[] { d.KZ, d.KX, d.KY }.Count(v => v != 0) == 1 && Math.Abs(d.KZ + d.KX + d.KY) == 2;
        foreach (var bonds in new[] { new[] { "XX", "YY", "ZZ" }, new[] { "XX", "YY" }, new[] { "ZZ" }, new[] { "XX" } })
        {
            var moves = LetterCube.Moves(LetterCube.ChainTerms(bonds, 3), 3);
            Assert.NotEmpty(moves);
            Assert.All(moves, d => Assert.True(Axis(d)));
        }
        var withField = LetterCube.Moves(LetterCube.ChainTerms(new[] { "XX", "YY", "ZZ", "X" }, 3), 3);
        Assert.Contains(withField, d => !Axis(d));
        Assert.All(withField.Where(d => !Axis(d)), d => Assert.True(d.KX == 0 && d.KZ == -d.KY && Math.Abs(d.KZ) == 1));
    }

    // §8: the parity rule's 0. The moves of XX + YY that keep k_Z are exactly its steps along k_X and k_Y; every
    // other move changes k_Z by +-2.
    [Fact]
    public void The_Parity_Rules_Zero_Is_A_Step_Along_Another_Axis()
    {
        var moves = LetterCube.Moves(LetterCube.ChainTerms(new[] { "XX", "YY" }, 3), 3);
        Assert.Equal(
            new HashSet<LetterCube.Point> { new(0, 2, 0), new(0, -2, 0), new(0, 0, 2), new(0, 0, -2) },
            moves.Where(d => d.KZ == 0).ToHashSet());
        Assert.Equal(new HashSet<LetterCube.Point> { new(2, 0, 0), new(-2, 0, 0) }, moves.Where(d => d.KZ != 0).ToHashSet());
    }

    // The edges at N = 3, pinned against the Python primitive's cube_edges (12, 12, 8 and 22 for Heisenberg, XY,
    // Ising and Heisenberg with an X field); the identity, which commutes with everything, lies on no edge; and a
    // three-site term XYZ at N = 4 makes the zero step, which Moves keeps and Edges drops (42 edges in Python).
    [Fact]
    public void The_Edges_Are_Axis_Steps_Between_Reached_Points()
    {
        foreach (var (templates, count) in new[] { (new[] { "XX", "YY", "ZZ" }, 12), (new[] { "XX", "YY" }, 12),
                                                   (new[] { "ZZ" }, 8), (new[] { "XX", "YY", "ZZ", "X" }, 22) })
            Assert.Equal(count, LetterCube.Edges(LetterCube.ChainTerms(templates, 3), 3).Count);
        var xyz = LetterCube.ChainTerms(new[] { "XYZ" }, 4);
        Assert.Contains(O, LetterCube.Moves(xyz, 4));
        var xyzEdges = LetterCube.Edges(xyz, 4);
        Assert.Equal(42, xyzEdges.Count);
        Assert.DoesNotContain(xyzEdges, e => e.Item1 == e.Item2);
        var terms = LetterCube.ChainTerms(new[] { "XX", "YY", "ZZ" }, 3);
        var moves = LetterCube.Moves(terms, 3);
        var points = LetterCube.Points(3);
        foreach (var (a, b) in LetterCube.Edges(terms, 3))
        {
            Assert.True(points.ContainsKey(a) && points.ContainsKey(b));
            var d = new LetterCube.Point(b.KZ - a.KZ, b.KX - a.KX, b.KY - a.KY);
            Assert.True(moves.Contains(d) || moves.Contains(new LetterCube.Point(-d.KZ, -d.KX, -d.KY)));
        }
        Assert.DoesNotContain(LetterCube.Edges(terms, 3), e => e.Item1 == O || e.Item2 == O);
    }

    static Complex Phase(int k) => (k % 4) switch { 0 => 1, 1 => Complex.ImaginaryOne, 2 => -1, _ => -Complex.ImaginaryOne };

    // ---- dense helpers, built letter by letter (the second route) ----

    static Complex[,] Letter(char c) => c switch
    {
        'I' => new Complex[,] { { 1, 0 }, { 0, 1 } },
        'X' => new Complex[,] { { 0, 1 }, { 1, 0 } },
        'Y' => new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } },
        'Z' => new Complex[,] { { 1, 0 }, { 0, -1 } },
        _ => throw new ArgumentException(c.ToString()),
    };

    static Complex[,] Kron(Complex[,] a, Complex[,] b)
    {
        int ra = a.GetLength(0), rb = b.GetLength(0);
        var k = new Complex[ra * rb, ra * rb];
        for (int i = 0; i < ra; i++)
            for (int j = 0; j < ra; j++)
                for (int p = 0; p < rb; p++)
                    for (int q = 0; q < rb; q++)
                        k[i * rb + p, j * rb + q] = a[i, j] * b[p, q];
        return k;
    }

    static Complex[,] DenseString(string letters)
    {
        var m = Letter(letters[0]);
        for (int l = 1; l < letters.Length; l++) m = Kron(m, Letter(letters[l]));
        return m;
    }

    static Complex[,] Mul(Complex[,] a, Complex[,] b)
    {
        int d = a.GetLength(0);
        var c = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int k = 0; k < d; k++)
            {
                if (a[i, k] == Complex.Zero) continue;
                for (int j = 0; j < d; j++) c[i, j] += a[i, k] * b[k, j];
            }
        return c;
    }

    static long Multinomial(int n, int nx, int ny, int nz)
    {
        BigInteger f(int m) { BigInteger r = 1; for (int i = 2; i <= m; i++) r *= i; return r; }
        return (long)(f(n) / (f(nx) * f(ny) * f(nz) * f(n - nx - ny - nz)));
    }
}
