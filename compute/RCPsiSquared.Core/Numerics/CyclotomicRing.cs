namespace RCPsiSquared.Core.Numerics;

/// <summary>Exact arithmetic in ℤ[ζ_m] = ℤ[x]/Φ_m(x): elements are integer coefficient
/// vectors of length φ(m) on the free ℤ-basis {1, x, …, x^{φ(m)−1}}, so vector equality is
/// ring equality is equality of the represented algebraic numbers (no floats, no primes, no
/// sampling). Since Φ_m | x^m − 1, exponents reduce mod m (x^m ≡ 1); a root power ζ^e is
/// added via the precomputed remainder rem(x^{e mod m}, Φ_m), and products are folded by the
/// monic Φ_m. For every modulus used here (m = 2n ≤ 120, odd part &lt; 105 ⟹ at most two
/// distinct odd primes) all Φ_m and remainder-table coefficients lie in {−1, 0, 1}; all
/// arithmetic is checked 64-bit. This is the exact layer under
/// <see cref="LevelCollisionCensus"/> (F129) and <see cref="CollisionDecoupling"/> (F130),
/// the C# counterpart of the Python gates' ℤ[x]/Φ_2n convention
/// (<c>simulations/f129_level_collision_law.py</c>, <c>root_sum_vec</c>).</summary>
public static class CyclotomicRing
{
    private static readonly Dictionary<int, long[]> PhiCache = new();
    private static readonly Dictionary<int, long[][]> TableCache = new();
    private static readonly object Gate = new();

    /// <summary>Coefficients of the m-th cyclotomic polynomial, ascending, monic;
    /// computed exactly as (x^m − 1) / Π_{d|m, d&lt;m} Φ_d.</summary>
    public static long[] Phi(int m)
    {
        lock (Gate)
        {
            if (PhiCache.TryGetValue(m, out var cached)) return cached;
            var f = new long[m + 1];
            f[0] = -1;
            f[m] = 1;
            for (int d = 1; d < m; d++)
                if (m % d == 0)
                    f = DivideExact(f, Phi(d));
            PhiCache[m] = f;
            return f;
        }
    }

    /// <summary>Euler φ(m) = degree of Φ_m.</summary>
    public static int Degree(int m) => Phi(m).Length - 1;

    /// <summary>rem(x^e, Φ_m) for e = 0..m−1, each a vector of length φ(m).</summary>
    private static long[][] Table(int m)
    {
        lock (Gate)
        {
            if (TableCache.TryGetValue(m, out var cached)) return cached;
            var phi = Phi(m);
            int deg = phi.Length - 1;
            var table = new long[m][];
            var cur = new long[deg];
            cur[0] = 1;
            table[0] = (long[])cur.Clone();
            for (int e = 1; e < m; e++)
            {
                // multiply by x, then fold the (possible) degree-deg term by monic Φ_m
                var next = new long[deg + 1];
                Array.Copy(cur, 0, next, 1, deg);
                long lead = next[deg];
                if (lead != 0)
                    checked
                    {
                        for (int j = 0; j < deg; j++) next[j] -= lead * phi[j];
                    }
                cur = new long[deg];
                Array.Copy(next, cur, deg);
                table[e] = (long[])cur.Clone();
            }
            TableCache[m] = table;
            return table;
        }
    }

    public static long[] Zero(int m) => new long[Degree(m)];

    /// <summary>acc += coeff · ζ^exp (exponent taken mod m; negative exponents allowed).</summary>
    public static void AddRootPower(long[] acc, int m, long exp, long coeff)
    {
        int e = (int)(((exp % m) + m) % m);
        var row = Table(m)[e];
        checked
        {
            for (int j = 0; j < acc.Length; j++) acc[j] += coeff * row[j];
        }
    }

    /// <summary>Product of two reduced vectors, reduced again (fold by the monic Φ_m).</summary>
    public static long[] Multiply(long[] a, long[] b, int m)
    {
        var phi = Phi(m);
        int deg = phi.Length - 1;
        var prod = new long[2 * deg - 1];
        checked
        {
            for (int i = 0; i < a.Length; i++)
            {
                if (a[i] == 0) continue;
                for (int j = 0; j < b.Length; j++)
                    prod[i + j] += a[i] * b[j];
            }
            for (int e = prod.Length - 1; e >= deg; e--)
            {
                long lead = prod[e];
                if (lead == 0) continue;
                prod[e] = 0;
                for (int j = 0; j < deg; j++) prod[e - deg + j] -= lead * phi[j];
            }
        }
        var res = new long[deg];
        Array.Copy(prod, res, deg);
        return res;
    }

    public static long[] Add(long[] a, long[] b)
    {
        var r = new long[a.Length];
        checked
        {
            for (int j = 0; j < a.Length; j++) r[j] = a[j] + b[j];
        }
        return r;
    }

    public static long[] Negate(long[] a)
    {
        var r = new long[a.Length];
        checked
        {
            for (int j = 0; j < a.Length; j++) r[j] = -a[j];
        }
        return r;
    }

    public static bool IsZero(long[] v)
    {
        foreach (long c in v) if (c != 0) return false;
        return true;
    }

    public static bool AreEqual(long[] a, long[] b)
    {
        if (a.Length != b.Length) return false;
        for (int j = 0; j < a.Length; j++) if (a[j] != b[j]) return false;
        return true;
    }

    /// <summary>Canonical dictionary key of a reduced vector (used by the census grouping).</summary>
    public static string Key(long[] v) => string.Join(",", v);

    private static long[] DivideExact(long[] num, long[] den)
    {
        int dn = num.Length - 1, dd = den.Length - 1;
        var rem = (long[])num.Clone();
        var quo = new long[dn - dd + 1];
        checked
        {
            for (int k = dn - dd; k >= 0; k--)
            {
                long q = rem[k + dd] / den[dd];
                quo[k] = q;
                if (q == 0) continue;
                for (int j = 0; j <= dd; j++) rem[k + j] -= q * den[j];
            }
        }
        foreach (long c in rem)
            if (c != 0) throw new InvalidOperationException("cyclotomic division not exact");
        return quo;
    }
}
