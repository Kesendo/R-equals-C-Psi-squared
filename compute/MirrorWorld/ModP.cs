namespace MirrorWorld;

/// <summary>The world's exact arithmetic, in one place (adopted 2026-09-06 from the arc
/// mirrorworld_what_is_missing, NextStep (6)).
///
/// Several objects here walk past the wall by refusing floating point: Seed's nullity surplus,
/// Divisor's multiplicity, BlindSeat's blind count, LevelCollision's cyclotomic levels and Crack's
/// characteristic-polynomial identity are all EXACT ranks and residues over GF(p) at two primes,
/// with no eigensolver anywhere. This is the arithmetic all five of them are written in.
///
/// ONE PRIME LIST, and the choice is forced rather than preferred: both primes are 1 mod 4, so -1
/// is a square at each, and the GAUSSIAN ranks (Divisor works over Z[i], embedding i as a square
/// root of -1) can therefore share the list with the plain integer ranks. A list headed by
/// 2^31 - 1, which is 3 mod 4, cannot do that job at both of its primes, and a shared list is
/// exactly a list that must serve at both.
///
/// The trap this file exists to guard in one place: a reduction such as (a % p + p) evaluated in
/// int before it is widened, which for p near 2^30 wraps silently and returns a plausible wrong
/// residue. Everything below takes and returns long, reduces before it multiplies, refuses a
/// non-positive modulus and a ragged row set rather than returning a plausible number for them,
/// and multiplies through UInt128, which the list's own primes never need (both are below 2^30, so
/// a reduced product fits a long) but which is a contract of the primitive and is read where it
/// bites, at a modulus near 2^61 where a long product wraps.
///
/// ModPTests judges each method against a route independent of it: BigInteger for the modular
/// atoms, trial division for primality, a membership question over all primes for factorisation, a
/// rank known by construction for the elimination. The modular atoms are read at the top residues
/// p-1 and p-2, which are sufficient to separate an int path from a long one; they are not the only
/// inputs that do.
///
/// Deliberately NOT a GameObject: it states no physics, prints nothing and appears in no Own or
/// Inherited bucket. It is machinery the objects are written in, the shape Formulas and Topology
/// already have in this world.</summary>
public static class ModP
{
    /// <summary>The two primes, both 1 mod 4 (so SqrtMinusOne exists at each) and both below 2^31,
    /// so a reduced product stays inside long even before the UInt128 route.</summary>
    public static readonly long[] Primes = { 998244353L, 1004535809L };

    /// <summary>x reduced into [0, p), for either sign and for x far outside the range.</summary>
    public static long Mod(long x, long p)
    {
        if (p <= 0) throw new ArgumentOutOfRangeException(nameof(p), p, "a modulus must be positive");
        long r = x % p;
        return r < 0 ? r + p : r;
    }

    /// <summary>a*b mod p, through UInt128 so two top residues cannot wrap.</summary>
    public static long MulMod(long a, long b, long p)
        => (long)((UInt128)(ulong)Mod(a, p) * (ulong)Mod(b, p) % (ulong)p);

    /// <summary>b^e mod p by square-and-multiply.</summary>
    public static long ModPow(long b, long e, long p)
    {
        if (p <= 0) throw new ArgumentOutOfRangeException(nameof(p), p, "a modulus must be positive");
        long r = 1 % p;
        b = Mod(b, p);
        while (e > 0)
        {
            if ((e & 1) == 1) r = MulMod(r, b, p);
            b = MulMod(b, b, p);
            e >>= 1;
        }
        return r;
    }

    /// <summary>The inverse of x mod p by Fermat, p prime. Inputs of either sign.</summary>
    public static long ModInverse(long x, long p)
    {
        long r = Mod(x, p);
        if (r == 0) throw new ArgumentOutOfRangeException(nameof(x), x, $"0 has no inverse mod {p}");
        return ModPow(r, p - 2, p);
    }

    /// <summary>A square root of -1 mod p, which exists exactly when p = 1 mod 4. This is what lets
    /// a Gaussian-integer matrix be ranked over GF(p): a + b*i maps to a + b*r.</summary>
    public static long SqrtMinusOne(long p)
    {
        if (Mod(p, 4) != 1)
            throw new InvalidOperationException($"no square root of -1 mod {p}: p must be 1 mod 4, this is {Mod(p, 4)} mod 4");
        for (long t = 2; t < 200; t++)
        {
            long r = ModPow(t, (p - 1) / 4, p);
            if (MulMod(r, r, p) == p - 1) return r;
        }
        throw new InvalidOperationException($"no square root of -1 mod {p} found below t = 200");
    }

    /// <summary>Deterministic Miller-Rabin: the twelve bases below decide primality exactly for
    /// every m under 3.3 * 10^24, so this is a proof and not a probable-prime test in our range.</summary>
    public static bool IsPrime(long m)
    {
        if (m < 2) return false;
        long[] bases = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
        foreach (long b in bases) { if (m % b == 0) return m == b; }
        long d = m - 1;
        int s = 0;
        while ((d & 1) == 0) { d >>= 1; s++; }
        foreach (long b in bases)
        {
            long x = ModPow(b, d, m);
            if (x == 1 || x == m - 1) continue;
            bool witness = true;
            for (int i = 1; i < s; i++)
            {
                x = MulMod(x, x, m);
                if (x == m - 1) { witness = false; break; }
            }
            if (witness) return false;
        }
        return true;
    }

    /// <summary>The distinct prime factors of m, ascending.</summary>
    public static int[] PrimeFactors(int m)
    {
        var qs = new List<int>();
        for (int q = 2; q * q <= m; q++)
            if (m % q == 0) { qs.Add(q); while (m % q == 0) m /= q; }
        if (m > 1) qs.Add(m);
        return qs.ToArray();
    }

    /// <summary>An element of EXACT multiplicative order `order` mod p (the cyclotomic atom: it is
    /// the zeta the level and resonance counts evaluate their cosines at). Requires order | p-1;
    /// returns 0 if no such element is found below t = 500.</summary>
    public static long RootOfOrder(int order, long p)
    {
        var qs = PrimeFactors(order);
        for (long x = 2; x < 500; x++)
        {
            long z = ModPow(x, (p - 1) / order, p);
            if (z == 1) continue;
            bool full = true;
            foreach (int q in qs)
                if (ModPow(z, order / q, p) == 1) { full = false; break; }
            if (full) return z;
        }
        return 0;
    }

    /// <summary>The rank of the row set over GF(p), by Gaussian elimination. Rows may be longer or
    /// shorter in number than the column count; entries may be of either sign and unreduced.</summary>
    public static int Rank(IReadOnlyList<long[]> rows, long p)
    {
        if (p <= 0) throw new ArgumentOutOfRangeException(nameof(p), p, "a modulus must be positive");
        int n = rows.Count;
        if (n == 0) return 0;
        int cols = rows[0].Length;
        for (int i = 1; i < n; i++)
            if (rows[i].Length != cols)
                throw new ArgumentException($"ragged input: row 0 has {cols} entries, row {i} has {rows[i].Length}", nameof(rows));

        var a = new long[n][];
        for (int i = 0; i < n; i++)
        {
            a[i] = new long[cols];
            for (int j = 0; j < cols; j++) a[i][j] = Mod(rows[i][j], p);
        }

        int rank = 0;
        for (int c = 0; c < cols && rank < n; c++)
        {
            int piv = -1;
            for (int r = rank; r < n; r++)
                if (a[r][c] != 0) { piv = r; break; }
            if (piv < 0) continue;

            (a[rank], a[piv]) = (a[piv], a[rank]);
            long inv = ModInverse(a[rank][c], p);
            for (int j = c; j < cols; j++) a[rank][j] = MulMod(a[rank][j], inv, p);
            for (int r = 0; r < n; r++)
            {
                if (r == rank || a[r][c] == 0) continue;
                long f = a[r][c];
                for (int j = c; j < cols; j++)
                    a[r][j] = Mod(a[r][j] - MulMod(f, a[rank][j], p), p);
            }
            rank++;
        }
        return rank;
    }

    /// <summary>The rank over the integers (or over Z[i] once i has been embedded), as the MAXIMUM
    /// of the GF(p) ranks at the two primes. One-sided and sound: reduction can only drop the rank,
    /// never raise it, so the larger reading is the one closer to the truth and a rank certified
    /// here is a lower bound that no prime can inflate.</summary>
    public static int Rank(IReadOnlyList<long[]> rows) => Primes.Max(p => Rank(rows, p));
}
