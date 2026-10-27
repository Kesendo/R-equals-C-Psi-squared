namespace MirrorWorld;

using static MirrorWorld.LetterCube;

// F138 off the letters (the open arc f138_clause_two_sweep, its tilted half), the draws and the cube's rule the
// sextant reads them against; machinery for the run mode `sextant tilted` and its test, not a GameObject.
// docs/ANALYTICAL_FORMULAS.md F138 states its two clauses for letter axes; F158 holds for dephasing along any axis n·σ
// (docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md §(h)) and the colouring reads the same way with colours on a circle
// (experiments/THE_PALINDROME_AS_A_COLOURING.md, "Colours on the circle"): a site's colours are the m·σ with m
// perpendicular to its watched axes, and a field on it lies along its colour. A bond sum_k w_k σ_k⊗σ_k with coupling
// matrix W commutes with (a·σ)⊗(b·σ) exactly when (W·a)×b = 0 = (W·b)×a, so where W·a ≠ 0 the bond carries the
// colour a of one end to W·a at the other: Heisenberg (W = 1) takes one colour per component, perpendicular to every watched axis and along
// every field, so the axes lie in one plane and the fields along its normal; XXZ at Δ ≠ ±1 (drawn at Δ = 2) and
// XX + YY take that colour along Z or perpendicular to it (at Δ = -1 the bond would flip its Z part, which the rule
// as written does not read); ZZ takes both ends along Z or both perpendicular to it, one colour per site. The rule is
// sufficient (a colour is an element of F158's far kernel) and not necessary (a sum of circle-letter strings can pair a row
// with no colour, as the test's turned XX + YY row does); the draws read the rule against the sextant's exact decision.
public static class TiltedAxes
{
    public sealed record BondClass(string Name, string[] Terms, long[] Coefficients);

    public static readonly IReadOnlyList<BondClass> Classes = new[]
    {
        new BondClass("Heisenberg", new[] { "XX", "YY", "ZZ" }, new[] { 1L, 1, 1 }),
        new BondClass("XXZ", new[] { "XX", "YY", "ZZ" }, new[] { 1L, 1, 2 }),
        new BondClass("XX+YY", new[] { "XX", "YY" }, new[] { 1L, 1 }),
        new BondClass("ZZ", new[] { "ZZ" }, new[] { 1L }),
    };

    public static IReadOnlyList<(string Name, (int, int)[] Edges)> Graphs(int n)
    {
        var g = new List<(string, (int, int)[])>
        {
            ("chain", Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToArray()),
            ("ring", Enumerable.Range(0, n).Select(i => (i, (i + 1) % n)).ToArray()),
        };
        if (n == 4)
        {
            g.Add(("star", new[] { (0, 1), (0, 2), (0, 3) }));
            g.Add(("K4", new[] { (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3) }));
        }
        return g;
    }

    public sealed record Row(int N, (int, int)[] Edges, BondClass Class, (long X, long Y, long Z)[] Axes, long[] Rates,
        (long X, long Y, long Z)[] Fields)
    {
        public Sextant.Setting Setting()
        {
            var h = new List<(PauliString, long)>();
            string Placed(params (int Site, char Letter)[] at)
            {
                var s = new char[N]; Array.Fill(s, 'I');
                foreach (var (site, letter) in at) s[site] = letter;
                return new string(s);
            }
            foreach (var (i, j) in Edges)
                for (int k = 0; k < Class.Terms.Length; k++)
                    h.Add((PauliString.Parse(Placed((i, Class.Terms[k][0]), (j, Class.Terms[k][1]))), Class.Coefficients[k]));
            for (int l = 0; l < N; l++)
                foreach (var (letter, v) in new[] { ('X', Fields[l].X), ('Y', Fields[l].Y), ('Z', Fields[l].Z) })
                    if (v != 0) h.Add((PauliString.Parse(Placed((l, letter))), v));
            return new Sextant.Setting(N, h, Enumerable.Range(0, N).Select(l => (l, Jump(Axes[l]), Rates[l])));
        }

        // the cube's rule (the head comment); every graph drawn is connected, one component, and every site watched (an
        // unwatched site, free to carry the identity, is outside the rule as written)
        public bool Rule()
        {
            if (Rates.Any(r => r <= 0)) throw new InvalidOperationException("the rule is written for a positive rate on every site");
            var z = (0L, 0L, 1L);
            var fields = Fields.Where(f => f != (0, 0, 0)).Select(Primitive).Distinct().ToList();
            bool kept((long, long, long) m) => Class.Name == "Heisenberg" || Primitive(m) == z || Dot(m, z) == 0;
            if (Class.Name == "ZZ")
                return (Axes.All(x => Dot(x, z) == 0) && fields.All(f => f == z))
                    || Enumerable.Range(0, N).All(l => Fields[l] == (0, 0, 0) || (Dot(Fields[l], z) == 0 && Dot(Fields[l], Axes[l]) == 0));
            if (fields.Count > 1) return false;
            if (fields.Count == 1) return Axes.All(x => Dot(fields[0], x) == 0) && kept(fields[0]);
            var normals = Axes.SelectMany(x => Axes.Select(y => Primitive(Cross(x, y)))).Where(v => v != (0, 0, 0)).Distinct().ToList();
            return normals.Count == 0 || (normals.Count == 1 && Axes.All(x => Dot(normals[0], x) == 0) && kept(normals[0]));
        }

        // F138's two clauses read with directions: at most two distinct axes (any number under the one-letter ZZ bond),
        // and the nonzero fields along one common axis perpendicular to every watched axis
        public bool Clauses()
        {
            var fields = Fields.Where(f => f != (0, 0, 0)).Select(Primitive).Distinct().ToList();
            bool clause1 = Axes.Select(Primitive).Distinct().Count() <= 2 || Class.Name == "ZZ";
            bool clause2 = fields.Count == 0 || (fields.Count == 1 && Axes.All(x => Dot(fields[0], x) == 0));
            return clause1 && clause2;
        }

        public override string ToString()
            => $"{Class.Name} on [{string.Join(" ", Edges.Select(e => $"{e.Item1}-{e.Item2}"))}], axes {string.Join(" ", Axes)}, rates {string.Join(",", Rates)}, fields {string.Join(" ", Fields)}";
    }

    // the integer directions with entries from -2 to 2, primitive, first nonzero entry positive
    public static IReadOnlyList<(long X, long Y, long Z)> Directions()
    {
        var dirs = new List<(long, long, long)>();
        for (long x = -2; x <= 2; x++) for (long y = -2; y <= 2; y++) for (long z = -2; z <= 2; z++)
            if ((x, y, z) != (0, 0, 0) && Primitive((x, y, z)) == (x, y, z)) dirs.Add((x, y, z));
        return dirs;
    }

    // one draw: a bond class, two axes a, b and a third (a third of the time in the plane of the two), each site taking
    // one of them, at least one off the letters, rates 1 to 3; no field, a common field (along a direction perpendicular
    // to every axis when one exists, three times in four) with magnitudes -2 to 2, or a field per site
    public static Row Draw(Random rnd, int n, (int, int)[] edges)
    {
        var dirs = Directions();
        while (true)
        {
            var cls = Classes[rnd.Next(Classes.Count)];
            var a = dirs[rnd.Next(dirs.Count)];
            var b = dirs[rnd.Next(dirs.Count)];
            var c = dirs[rnd.Next(dirs.Count)];
            if (rnd.Next(3) == 0)
            {
                long p = rnd.Next(1, 3), q = rnd.Next(-2, 3);
                c = Primitive((p * a.X + q * b.X, p * a.Y + q * b.Y, p * a.Z + q * b.Z));
                if (c == (0, 0, 0)) c = a;
            }
            var axes = Enumerable.Range(0, n).Select(_ => rnd.Next(4) == 0 ? c : rnd.Next(2) == 0 ? a : b).ToArray();
            var rates = Enumerable.Range(0, n).Select(_ => (long)rnd.Next(1, 4)).ToArray();
            var fields = new (long X, long Y, long Z)[n];
            int mode = rnd.Next(3);
            if (mode == 1)
            {
                var perp = dirs.Where(f => axes.All(x => Dot(f, x) == 0)).ToList();
                var f = perp.Count > 0 && rnd.Next(4) != 0 ? perp[rnd.Next(perp.Count)] : dirs[rnd.Next(dirs.Count)];
                for (int l = 0; l < n; l++) { long m = rnd.Next(-2, 3); fields[l] = (m * f.X, m * f.Y, m * f.Z); }
            }
            else if (mode == 2)
                for (int l = 0; l < n; l++) fields[l] = rnd.Next(3) == 0 ? (0, 0, 0) : dirs[rnd.Next(dirs.Count)];
            if (axes.All(IsLetter)) continue;
            return new Row(n, edges, cls, axes, rates, fields);
        }
    }

    static bool IsLetter((long X, long Y, long Z) v) => (v.X != 0 ? 1 : 0) + (v.Y != 0 ? 1 : 0) + (v.Z != 0 ? 1 : 0) == 1;

    static Gauss[,] Jump((long X, long Y, long Z) n)
        => new[,] { { new Gauss(n.Z, 0), new Gauss(n.X, -n.Y) }, { new Gauss(n.X, n.Y), new Gauss(-n.Z, 0) } };

    static long Dot((long X, long Y, long Z) a, (long X, long Y, long Z) b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;

    static (long, long, long) Cross((long X, long Y, long Z) a, (long X, long Y, long Z) b)
        => (a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);

    static (long, long, long) Primitive((long X, long Y, long Z) v)
    {
        static long Gcd(long a, long b) { a = Math.Abs(a); b = Math.Abs(b); while (b != 0) (a, b) = (b, a % b); return a; }
        long g = Gcd(Gcd(v.X, v.Y), v.Z);
        if (g == 0) return v;
        long sign = (v.X != 0 ? v.X : v.Y != 0 ? v.Y : v.Z) < 0 ? -1 : 1;
        return (sign * v.X / g, sign * v.Y / g, sign * v.Z / g);
    }
}
