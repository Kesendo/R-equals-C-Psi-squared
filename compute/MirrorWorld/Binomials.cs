namespace MirrorWorld;

using System.Numerics;

// Binomial equations over the nonzero complex numbers, kept as equations and never as numbers. Built 2026-10-20 for the
// theodolite, beside its Python twin simulations/framework/diagnostics/cube_letters.py (solve_binomials). Machinery,
// like ModP, not a GameObject: it states no physics and appears in no Own or Inherited bucket.
//
// One equation reads  prod_u u^e = Num / Den  with integer exponents e and nonzero Gaussian integers Num and Den. A
// system is decided by Euclid on the exponent columns: every row operation is unimodular (swap two rows, or subtract an
// integer multiple of one row from another), and the constants ride along multiplicatively, so no root is taken and no
// fraction is formed. The system has a solution in nonzero complex numbers exactly when every row the elimination
// empties reads Num = Den: the operations being unimodular, those rows span every integer relation among the exponent
// rows, and the rows that keep a pivot form a triangle that is solved one root at a time from the last unknown up, the
// free unknowns set to anything nonzero. The answer is the reduced system itself, the equations its unknowns satisfy:
// x^2 = 3, not a square root of 3.
//
// Deliberately outside: the numbers (a radical, a principal root, a representative: the Python twin reads those; only
// RationalSolution gives one, where no root is needed), any ordering or positivity, and Groebner bases (a binomial
// system never needs one).
//
// Words, fenced at the door: "reduced" is the triangle left by the elimination, not a reduced density matrix or a
// reduced word of the letter cube; "free" is an unknown without a pivot (a gauge), not a free fermion.
public readonly record struct GaussBig(BigInteger Re, BigInteger Im)
{
    public static readonly GaussBig One = new(BigInteger.One, BigInteger.Zero);
    public static GaussBig Of(long re, long im = 0) => new(re, im);
    public static GaussBig Of(LetterCube.Gauss g) => new(g.Re, g.Im);
    public static GaussBig operator *(GaussBig a, GaussBig b) => new(a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re);
    public static GaussBig operator -(GaussBig a, GaussBig b) => new(a.Re - b.Re, a.Im - b.Im);
    public static GaussBig operator -(GaussBig a) => new(-a.Re, -a.Im);
    public GaussBig Conj => new(Re, -Im);
    public bool IsZero => Re.IsZero && Im.IsZero;
    public BigInteger Norm => Re * Re + Im * Im;

    // the Gaussian integer nearest to a / b (each part rounded half up); the remainder a - q·b has a smaller norm than b
    public static GaussBig NearestQuotient(GaussBig a, GaussBig b)
    {
        var n = b.Norm;
        var p = a * b.Conj;
        return new(RoundDiv(p.Re, n), RoundDiv(p.Im, n));
    }

    static BigInteger RoundDiv(BigInteger x, BigInteger n)
    {
        var t = 2 * x + n;
        var q = BigInteger.DivRem(t, 2 * n, out var r);
        return r.Sign < 0 ? q - 1 : q;                // floor((2x + n) / 2n)
    }

    // a greatest common divisor, Euclid in Z[i] (up to a unit)
    public static GaussBig Gcd(GaussBig a, GaussBig b)
    {
        while (!b.IsZero)
        {
            var r = a - NearestQuotient(a, b) * b;
            a = b;
            b = r;
        }
        return a;
    }

    // a / g, exact (throws when g does not divide a)
    public static GaussBig ExactQuotient(GaussBig a, GaussBig g)
    {
        var n = g.Norm;
        var p = a * g.Conj;
        if (!(p.Re % n).IsZero || !(p.Im % n).IsZero) throw new ArgumentException($"{g} does not divide {a}");
        return new(p.Re / n, p.Im / n);
    }

    public GaussBig Pow(int k)
    {
        if (k < 0) throw new ArgumentOutOfRangeException(nameof(k), k, "a Gaussian integer has no negative power");
        var r = One;
        for (int i = 0; i < k; i++) r *= this;
        return r;
    }

    public override string ToString()
        => Im.IsZero ? $"{Re}" : Re.IsZero ? $"{Im}i" : $"{Re}{(Im.Sign < 0 ? "-" : "+")}{BigInteger.Abs(Im)}i";
}

// prod_u u^Exponents[u] = Num / Den, with the place it was read from
public sealed record Binomial(IReadOnlyDictionary<string, int> Exponents, GaussBig Num, GaussBig Den, string Origin)
{
    public static string Monomial(IEnumerable<KeyValuePair<string, int>> exponents)
    {
        var parts = exponents.Where(kv => kv.Value != 0).OrderBy(kv => kv.Key, StringComparer.Ordinal)
            .Select(kv => kv.Value == 1 ? kv.Key : $"{kv.Key}^{kv.Value}").ToList();
        return parts.Count == 0 ? "1" : string.Join("·", parts);
    }

    public static string Ratio(GaussBig num, GaussBig den) => den == GaussBig.One ? $"{num}" : $"({num})/({den})";

    public override string ToString() => $"{Monomial(Exponents)} = {Ratio(Num, Den)}";
}

public static class Binomials
{
    // The triangle left by the elimination: row r has its pivot at column Pivots[r] and zeros before it; Free are the
    // unknowns without a pivot.
    public sealed record Reduced(IReadOnlyList<string> Unknowns, IReadOnlyList<int[]> Exponents,
        IReadOnlyList<(GaussBig Num, GaussBig Den)> Constants, IReadOnlyList<int> Pivots)
    {
        public IReadOnlyList<string> Free => Unknowns.Where((_, j) => !Pivots.Contains(j)).ToList();

        // The row whose pivot is this unknown, when that row holds this unknown alone: u^k = Num / Den with k > 0.
        public (int K, GaussBig Num, GaussBig Den)? Alone(string unknown)
        {
            int col = IndexOf(unknown);
            int r = Pivots.ToList().IndexOf(col);
            if (r < 0) return null;
            var e = Exponents[r];
            for (int j = 0; j < e.Length; j++)
                if (j != col && e[j] != 0) return null;
            var (num, den) = Constants[r];
            return e[col] > 0 ? (e[col], num, den) : (-e[col], den, num);
        }

        int IndexOf(string unknown)
        {
            for (int j = 0; j < Unknowns.Count; j++)
                if (Unknowns[j] == unknown) return j;
            throw new ArgumentException($"'{unknown}' is not an unknown of this system", nameof(unknown));
        }

        public IEnumerable<string> Lines()
        {
            for (int r = 0; r < Exponents.Count; r++)
            {
                var e = Exponents[r];
                var exps = Unknowns.Select((u, j) => new KeyValuePair<string, int>(u, e[j]));
                yield return $"{Binomial.Monomial(exps)} = {Binomial.Ratio(Constants[r].Num, Constants[r].Den)}";
            }
        }
    }

    // The reduced system, or null when the equations have no solution in nonzero complex numbers.
    public static Reduced? Solve(IEnumerable<Binomial> equations, IReadOnlyList<string> unknowns)
    {
        var index = new Dictionary<string, int>();
        for (int j = 0; j < unknowns.Count; j++) index[unknowns[j]] = j;
        var rows = new List<(int[] E, GaussBig Num, GaussBig Den)>();
        foreach (var b in equations)
        {
            if (b.Num.IsZero || b.Den.IsZero)
                throw new ArgumentException($"a binomial's constant must be nonzero: {b}", nameof(equations));
            var e = new int[unknowns.Count];
            foreach (var (u, k) in b.Exponents)
            {
                if (!index.TryGetValue(u, out int j))
                    throw new ArgumentException($"'{u}' in {b} is not among the unknowns", nameof(equations));
                e[j] += k;
            }
            var (cn, cd) = Canonical(b.Num, b.Den);
            rows.Add((e, cn, cd));
        }
        int r = 0;
        var pivots = new List<int>();
        for (int col = 0; col < unknowns.Count; col++)
        {
            while (true)
            {
                int best = -1;
                for (int i = r; i < rows.Count; i++)
                    if (rows[i].E[col] != 0 && (best < 0 || Math.Abs(rows[i].E[col]) < Math.Abs(rows[best].E[col])))
                        best = i;
                if (best < 0) break;
                (rows[r], rows[best]) = (rows[best], rows[r]);
                bool clean = true;
                for (int i = r + 1; i < rows.Count; i++)
                {
                    int a = rows[i].E[col];
                    if (a == 0) continue;
                    rows[i] = Subtract(rows[i], rows[r], a / rows[r].E[col]);
                    if (rows[i].E[col] != 0) clean = false;
                }
                if (clean) break;
            }
            bool pivot = false;
            for (int i = r; i < rows.Count; i++)
                if (rows[i].E[col] != 0) pivot = true;
            if (pivot)
            {
                pivots.Add(col);
                r++;
            }
        }
        for (int i = r; i < rows.Count; i++)            // an emptied row: its constant must be 1
            if (rows[i].Num != rows[i].Den) return null;
        for (int i = 0; i < r; i++)                     // every pivot exponent positive (a row and its inverse say the same)
            if (rows[i].E[pivots[i]] < 0)
            {
                var (fn, fd) = Canonical(rows[i].Den, rows[i].Num);
                rows[i] = (rows[i].E.Select(v => -v).ToArray(), fn, fd);
            }
        return new Reduced(unknowns, rows.Take(r).Select(t => t.E).ToList(),
            rows.Take(r).Select(t => (t.Num, t.Den)).ToList(), pivots);
    }

    // row_i - q·row_r: the exponents subtract, the constant is divided by the pivot row's constant to the power q
    static (int[] E, GaussBig Num, GaussBig Den) Subtract((int[] E, GaussBig Num, GaussBig Den) i,
        (int[] E, GaussBig Num, GaussBig Den) p, int q)
    {
        var e = new int[i.E.Length];
        for (int j = 0; j < e.Length; j++) e[j] = i.E[j] - q * p.E[j];
        var (num, den) = q >= 0
            ? (i.Num * p.Den.Pow(q), i.Den * p.Num.Pow(q))
            : (i.Num * p.Num.Pow(-q), i.Den * p.Den.Pow(-q));
        var (rn, rd) = Canonical(num, den);
        return (e, rn, rd);
    }

    // Num and Den divided by their common divisor, so the constants stay small through a long elimination
    static (GaussBig, GaussBig) Lowest(GaussBig num, GaussBig den)
    {
        var g = GaussBig.Gcd(num, den);
        return g.Norm == 1 ? (num, den) : (GaussBig.ExactQuotient(num, g), GaussBig.ExactQuotient(den, g));
    }

    // in lowest terms, the denominator turned by a unit into the quadrant Re > 0, Im >= 0: one ratio, one spelling
    static (GaussBig, GaussBig) Canonical(GaussBig num, GaussBig den)
    {
        (num, den) = Lowest(num, den);
        var i = new GaussBig(BigInteger.Zero, BigInteger.One);
        for (int k = 0; k < 4 && !(den.Re.Sign > 0 && den.Im.Sign >= 0); k++)
        {
            num *= i;
            den *= i;
        }
        return (num, den);
    }

    // A solution with every free unknown 1, as Num/Den pairs, when no irrational root is needed: every pivot exponent
    // is +-1, or +-2 with a ratio whose square root is again a ratio of Gaussian integers; otherwise null.
    public static IReadOnlyDictionary<string, (GaussBig Num, GaussBig Den)>? RationalSolution(Reduced reduced)
    {
        var val = new (GaussBig Num, GaussBig Den)[reduced.Unknowns.Count];
        for (int j = 0; j < val.Length; j++) val[j] = (GaussBig.One, GaussBig.One);
        for (int r = reduced.Exponents.Count - 1; r >= 0; r--)
        {
            int col = reduced.Pivots[r];
            var e = reduced.Exponents[r];
            var (num, den) = reduced.Constants[r];
            for (int j = col + 1; j < e.Length; j++)        // divide by the later unknowns' powers
            {
                if (e[j] == 0) continue;
                var (vn, vd) = val[j];
                if (e[j] > 0) { num *= vd.Pow(e[j]); den *= vn.Pow(e[j]); }
                else { num *= vn.Pow(-e[j]); den *= vd.Pow(-e[j]); }
            }
            if (e[col] < 0) (num, den) = (den, num);
            (num, den) = Lowest(num, den);
            switch (Math.Abs(e[col]))
            {
                case 1:
                    val[col] = (num, den);
                    break;
                case 2:                                       // sqrt(num/den) = sqrt(num·den)/den
                    var root = SquareRoot(num * den);
                    if (root is null) return null;
                    val[col] = (root.Value, den);
                    break;
                default:
                    return null;
            }
        }
        var sol = new Dictionary<string, (GaussBig, GaussBig)>();
        for (int j = 0; j < val.Length; j++) sol[reduced.Unknowns[j]] = val[j];
        return sol;
    }

    // An exact square root in the Gaussian integers, or null: (c + di)^2 = a + bi needs c^2 = (|g| + a)/2 and
    // d^2 = (|g| - a)/2 with |g| = sqrt(a^2 + b^2) an integer, and 2cd = b.
    public static GaussBig? SquareRoot(GaussBig g)
    {
        if (g.IsZero) return g;
        var norm = g.Norm;
        var m = ISqrt(norm);
        if (m * m != norm || !((m + g.Re) % 2).IsZero) return null;
        var c2 = (m + g.Re) / 2;
        var d2 = (m - g.Re) / 2;
        var c = ISqrt(c2);
        var d = ISqrt(d2);
        if (c * c != c2 || d * d != d2) return null;
        if (g.Im.Sign < 0) d = -d;
        var root = new GaussBig(c, d);
        return root * root == g ? root : null;
    }

    // floor(sqrt(n)), Newton from above: the iterates fall monotonically to the floor
    static BigInteger ISqrt(BigInteger n)
    {
        if (n.Sign < 0) throw new ArgumentOutOfRangeException(nameof(n), "no real square root");
        if (n < 2) return n;
        var x = BigInteger.One << (int)((n.GetBitLength() + 1) / 2);
        while (true)
        {
            var y = (x + n / x) / 2;
            if (y >= x) return x;
            x = y;
        }
    }
}
