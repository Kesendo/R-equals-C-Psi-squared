using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The exact centre-line count behind <see cref="SelfMirrorObject"/>: how many eigenvalues
/// of L = −i[H, ·] + Σ_l γ_l (Z_l ρ Z_l − ρ), counted with algebraic multiplicity, sit on the line
/// Re λ = −σ and how many sit at the point λ = −σ, σ = Σ_l γ_l. No eigenvalue is computed.
///
/// <para>The route. Every double is a dyadic rational, so after one common power-of-two scale the
/// matrix of L + σ·I in the |a⟩⟨b| basis has Gaussian-integer entries, and its eigenvalues are
/// the scaled μ = λ + σ. The matrix splits into the connected blocks of its sparsity pattern (for
/// a number-conserving H they lie inside the joint-popcount sectors and can be finer, for a diagonal H
/// they are single entries),
/// and each block's characteristic polynomial p is computed division-free over Z[i]
/// (<see cref="GaussianMatrixCharpoly"/>). The point count is the order of the root μ = 0 of p.
/// The line count is the number of roots μ = iy with y real: writing p(iy) = A(y) + i·B(y) with
/// A, B integer polynomials, such a y is a common real root of A and B, with the multiplicity it
/// has in gcd(A, B), and those roots are counted with multiplicity by Sturm sequences over Z.
/// The scale changes neither count.</para>
///
/// <para>What the count belongs to: the doubles as given. γ = 0.1 is read as the double nearest
/// 1/10, and the answer is the exact answer for that input. A block larger than
/// <see cref="MaxBlockDimension"/> is not attempted (the division-free charpoly costs about
/// n⁴/4 big-integer products per n×n block); the count is then reported as unresolved. That bound
/// is a cost bound, not a precision one.</para></summary>
public static class CentreLineExactCount
{
    /// <summary>The largest connected block the exact route attempts: 36 = C(4,2)², the
    /// half-filling sector of a number-conserving H at N = 4.</summary>
    public const int MaxBlockDimension = 36;

    public sealed record Result(int LineCount, int PointCount, int LargestBlock, int BlockCount);

    /// <summary>The exact counts, or null when some connected block exceeds
    /// <see cref="MaxBlockDimension"/>.</summary>
    public static Result? TryCount(int n, ComplexMatrix hamiltonian, IReadOnlyList<double> gammas)
    {
        if (hamiltonian is null) throw new ArgumentNullException(nameof(hamiltonian));
        if (gammas is null) throw new ArgumentNullException(nameof(gammas));
        if (gammas.Count != n) throw new ArgumentException($"need one rate per site (N={n})", nameof(gammas));
        int d = 1 << n;
        if (hamiltonian.RowCount != d || hamiltonian.ColumnCount != d)
            throw new ArgumentException($"Hamiltonian must be {d}×{d}", nameof(hamiltonian));

        // One common exponent E: every input is m·2^e with e ≥ E, so m·2^(e−E) is an integer.
        var values = new List<double>();
        for (int a = 0; a < d; a++)
            for (int c = 0; c < d; c++)
            {
                values.Add(hamiltonian[a, c].Real);
                values.Add(hamiltonian[a, c].Imaginary);
            }
        values.AddRange(gammas);
        int minExp = int.MaxValue;
        foreach (double v in values)
        {
            var (m, e) = Dyadic(v);
            if (!m.IsZero) minExp = Math.Min(minExp, e);
        }
        if (minExp == int.MaxValue) minExp = 0;

        BigInteger Scaled(double v)
        {
            var (m, e) = Dyadic(v);
            return m.IsZero ? BigInteger.Zero : m << (e - minExp);
        }

        var hRe = new BigInteger[d, d];
        var hIm = new BigInteger[d, d];
        for (int a = 0; a < d; a++)
            for (int c = 0; c < d; c++)
            {
                hRe[a, c] = Scaled(hamiltonian[a, c].Real);
                hIm[a, c] = Scaled(hamiltonian[a, c].Imaginary);
            }
        var g = gammas.Select(Scaled).ToArray();

        // L + σI on |a⟩⟨b| (row index a·d + b):
        //   −i(Hρ)_ab = Σ_c (−i H_ac) ρ_cb,   +i(ρH)_ab = Σ_c (i H_cb) ρ_ac,
        //   dephasing + σ = Σ_l γ_l·(−1 if a, b disagree at site l else +1); site l is bit N−1−l.
        int dim = d * d;
        var rows = new Dictionary<int, GaussianInteger>[dim];
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
            {
                int r = a * d + b;
                var row = new Dictionary<int, GaussianInteger>();
                for (int c = 0; c < d; c++)
                {
                    if (!hRe[a, c].IsZero || !hIm[a, c].IsZero)
                        Add(row, c * d + b, new GaussianInteger(hIm[a, c], -hRe[a, c]));     // −i·H_ac
                    if (!hRe[c, b].IsZero || !hIm[c, b].IsZero)
                        Add(row, a * d + c, new GaussianInteger(-hIm[c, b], hRe[c, b]));     // +i·H_cb
                }
                BigInteger diag = BigInteger.Zero;
                for (int l = 0; l < n; l++)
                {
                    int bit = n - 1 - l;
                    diag += ((((a ^ b) >> bit) & 1) != 0) ? -g[l] : g[l];
                }
                Add(row, r, new GaussianInteger(diag, BigInteger.Zero));
                foreach (var key in row.Where(kv => kv.Value == GaussianInteger.Zero).Select(kv => kv.Key).ToList())
                    row.Remove(key);
                rows[r] = row;
            }

        // Connected blocks of the sparsity pattern (union-find).
        var parent = Enumerable.Range(0, dim).ToArray();
        int Find(int x)
        {
            while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
            return x;
        }
        for (int r = 0; r < dim; r++)
            foreach (int c in rows[r].Keys)
            {
                int pr = Find(r), pc = Find(c);
                if (pr != pc) parent[pr] = pc;
            }
        var blocks = new Dictionary<int, List<int>>();
        for (int r = 0; r < dim; r++)
        {
            int root = Find(r);
            if (!blocks.TryGetValue(root, out var list)) blocks[root] = list = new List<int>();
            list.Add(r);
        }
        int largest = blocks.Values.Max(b => b.Count);
        if (largest > MaxBlockDimension) return null;

        int line = 0, point = 0;
        foreach (var members in blocks.Values)
        {
            var p = BlockCharpoly(members, rows);
            point += OrderAtZero(p);
            line += ImaginaryAxisRootCount(p);
        }
        return new Result(line, point, largest, blocks.Count);
    }

    private static void Add(Dictionary<int, GaussianInteger> row, int col, GaussianInteger v) =>
        row[col] = row.TryGetValue(col, out var old) ? old + v : v;

    private static GaussianInteger[] BlockCharpoly(List<int> members, Dictionary<int, GaussianInteger>[] rows)
    {
        int m = members.Count;
        if (m == 1)
        {
            var diag = rows[members[0]].TryGetValue(members[0], out var v) ? v : GaussianInteger.Zero;
            return new[] { -diag, GaussianInteger.One };
        }
        var pos = new Dictionary<int, int>(m);
        for (int i = 0; i < m; i++) pos[members[i]] = i;
        var block = new GaussianInteger[m, m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++) block[i, j] = GaussianInteger.Zero;
        for (int i = 0; i < m; i++)
            foreach (var (col, v) in rows[members[i]])
                block[i, pos[col]] = v;
        return GaussianMatrixCharpoly.Characteristic(block);
    }

    /// <summary>The multiplicity of the root 0: the index of the lowest nonzero coefficient.</summary>
    public static int OrderAtZero(GaussianInteger[] lowFirst)
    {
        int k = 0;
        while (k < lowFirst.Length && lowFirst[k] == GaussianInteger.Zero) k++;
        if (k == lowFirst.Length) throw new ArgumentException("the zero polynomial has no root order");
        return k;
    }

    /// <summary>The number of roots of p on the imaginary axis, with multiplicity: the real roots
    /// of gcd(Re p(iy), Im p(iy)), counted with multiplicity.</summary>
    public static int ImaginaryAxisRootCount(GaussianInteger[] lowFirst)
    {
        int len = lowFirst.Length;
        var re = new BigInteger[len];
        var im = new BigInteger[len];
        for (int k = 0; k < len; k++)
        {
            var (a, b) = (lowFirst[k].Re, lowFirst[k].Im);
            // (a + ib)·i^k
            (re[k], im[k]) = (k & 3) switch
            {
                0 => (a, b),
                1 => (-b, a),
                2 => (-a, -b),
                _ => (b, -a),
            };
        }
        var A = Trim(re);
        var B = Trim(im);
        BigInteger[] gcd = A.Length == 0 ? Primitive(B) : B.Length == 0 ? Primitive(A) : Gcd(A, B);
        return RealRootsWithMultiplicity(gcd);
    }

    /// <summary>Real roots of an integer polynomial, counted with multiplicity: Σ_m of the distinct
    /// real roots of h_m, where h_0 = p and h_(m+1) = gcd(h_m, h_m′). A root of multiplicity r lies in
    /// h_0 … h_(r−1) and is therefore counted r times.</summary>
    public static int RealRootsWithMultiplicity(BigInteger[] p)
    {
        int total = 0;
        var h = Trim(p);
        while (h.Length > 1)
        {
            total += DistinctRealRoots(h);
            h = Gcd(h, Derivative(h));
        }
        return total;
    }

    /// <summary>Distinct real roots by Sturm's theorem (valid for a non-squarefree p too, since the
    /// sign counts are read at ±∞, where no root sits). Every remainder is kept as a POSITIVE multiple
    /// of the true remainder, so the sign pattern is the Sturm sequence's own.</summary>
    public static int DistinctRealRoots(BigInteger[] p)
    {
        var chain = new List<BigInteger[]> { Primitive(p), Primitive(Derivative(p)) };
        while (chain[^1].Length > 1)
        {
            var r = PositivePseudoRemainder(chain[^2], chain[^1]);
            if (r.Length == 0) break;
            chain.Add(Primitive(Negate(r)));
        }
        int SignChanges(bool atPlusInfinity)
        {
            int changes = 0, last = 0;
            foreach (var q in chain)
            {
                if (q.Length == 0) continue;
                int s = q[^1].Sign;
                if (!atPlusInfinity && ((q.Length - 1) & 1) == 1) s = -s;
                if (s == 0) continue;
                if (last != 0 && s != last) changes++;
                last = s;
            }
            return changes;
        }
        return SignChanges(atPlusInfinity: false) - SignChanges(atPlusInfinity: true);
    }

    private static BigInteger[] Gcd(BigInteger[] a, BigInteger[] b)
    {
        var x = Primitive(a);
        var y = Primitive(b);
        if (x.Length < y.Length) (x, y) = (y, x);
        while (y.Length > 0)
        {
            var r = Primitive(PositivePseudoRemainder(x, y));
            (x, y) = (y, r);
        }
        if (x.Length > 0 && x[^1].Sign < 0) x = Negate(x);
        return x;
    }

    /// <summary>c·rem(a, b) for some c &gt; 0: each elimination step multiplies by |lc(b)|.</summary>
    private static BigInteger[] PositivePseudoRemainder(BigInteger[] a, BigInteger[] b)
    {
        var r = (BigInteger[])a.Clone();
        int db = b.Length - 1;
        BigInteger lb = b[^1];
        BigInteger absLb = BigInteger.Abs(lb);
        int signLb = lb.Sign;
        r = Trim(r);
        while (r.Length - 1 >= db && r.Length > 0)
        {
            int shift = r.Length - 1 - db;
            BigInteger t = signLb * r[^1];
            var next = new BigInteger[r.Length];
            for (int i = 0; i < r.Length; i++) next[i] = absLb * r[i];
            for (int j = 0; j <= db; j++) next[j + shift] -= t * b[j];
            r = Trim(next);
        }
        return r;
    }

    private static BigInteger[] Derivative(BigInteger[] p)
    {
        if (p.Length <= 1) return Array.Empty<BigInteger>();
        var d = new BigInteger[p.Length - 1];
        for (int k = 1; k < p.Length; k++) d[k - 1] = k * p[k];
        return Trim(d);
    }

    private static BigInteger[] Negate(BigInteger[] p) => p.Select(c => -c).ToArray();

    private static BigInteger[] Primitive(BigInteger[] p)
    {
        var t = Trim(p);
        if (t.Length == 0) return t;
        BigInteger g = BigInteger.Zero;
        foreach (var c in t) g = BigInteger.GreatestCommonDivisor(g, c);
        return g.IsOne ? t : t.Select(c => c / g).ToArray();
    }

    private static BigInteger[] Trim(BigInteger[] p)
    {
        int len = p.Length;
        while (len > 0 && p[len - 1].IsZero) len--;
        return len == p.Length ? p : p.Take(len).ToArray();
    }

    /// <summary>A finite double as m·2^e exactly (m an integer).</summary>
    private static (BigInteger Mantissa, int Exponent) Dyadic(double x)
    {
        if (double.IsNaN(x) || double.IsInfinity(x))
            throw new ArgumentException("the exact count needs finite inputs", nameof(x));
        if (x == 0.0) return (BigInteger.Zero, 0);
        long bits = BitConverter.DoubleToInt64Bits(x);
        bool negative = bits < 0;
        int exponentField = (int)((bits >> 52) & 0x7FF);
        long fraction = bits & 0xFFFFFFFFFFFFFL;
        long mantissa;
        int exponent;
        if (exponentField == 0) { mantissa = fraction; exponent = -1074; }
        else { mantissa = fraction | (1L << 52); exponent = exponentField - 1075; }
        return (negative ? -new BigInteger(mantissa) : new BigInteger(mantissa), exponent);
    }
}
