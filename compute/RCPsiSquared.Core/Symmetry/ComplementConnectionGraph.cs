using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The exact routine behind <see cref="PalindromeComplementConnectionClaim"/>
/// (<c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>): H read on the bitstrings over the Gaussian
/// integers, and hom(A, B) = dim{g : B_xy·g_y = g_x·A_xy} counted as the good components of the pair graph.
/// Theorem 2 is hom(H, H) at the near end and hom(H, H̄) at the far end with every site dephased; Theorem 4
/// sums the same count over pairs of sectors. Every entry is a Gaussian integer and every edge is checked by
/// cross-multiplication, with no division and no float. The caller turns every jump (and every kept letter)
/// to Z before handing the strings in.</summary>
public static class ComplementConnectionGraph
{
    /// <summary>H in the computational basis, column by column: the result's [c] holds the entries
    /// ⟨r|H|c⟩ keyed by r, zeros dropped. A string P(x, z) = i^|x∧z|·X^x·Z^z (bit l = site l) sends |y⟩ to
    /// i^|x∧z|·(−1)^|z∧y|·|y ⊕ x⟩; each coefficient is a real integer.</summary>
    public static Dictionary<int, GaussianInteger>[] Columns(int n, IEnumerable<(PauliMask S, BigInteger C)> terms)
    {
        if (n < 1 || n > 30) throw new ArgumentOutOfRangeException(nameof(n), n, "the graph has 2^n vertices");
        int dim = 1 << n;
        var h = Enumerable.Range(0, dim).Select(_ => new Dictionary<int, GaussianInteger>()).ToArray();
        foreach (var (s, c) in terms)
        {
            if (s.X >> n != 0 || s.Z >> n != 0)
                throw new ArgumentException($"a string has a letter beyond site {n - 1}", nameof(terms));
            int x = (int)s.X, z = (int)s.Z;
            int phase = BitOperations.PopCount((uint)(x & z)) & 3;
            for (int y = 0; y < dim; y++)
            {
                var v = (BitOperations.PopCount((uint)(z & y)) & 1) == 0 ? c : -c;
                GaussianInteger entry = phase switch
                {
                    0 => new(v, BigInteger.Zero),
                    1 => new(BigInteger.Zero, v),
                    2 => new(-v, BigInteger.Zero),
                    _ => new(BigInteger.Zero, -v),
                };
                int target = y ^ x;
                var column = h[y];
                column[target] = (column.TryGetValue(target, out var old) ? old : GaussianInteger.Zero) + entry;
            }
        }
        foreach (var column in h)
            foreach (var key in column.Where(kv => kv.Value == GaussianInteger.Zero).Select(kv => kv.Key).ToList())
                column.Remove(key);
        return h;
    }

    /// <summary>⟨r|H|c⟩ from <see cref="Columns"/>.</summary>
    public static GaussianInteger Entry(Dictionary<int, GaussianInteger>[] h, int r, int c) =>
        h[c].TryGetValue(r, out var v) ? v : GaussianInteger.Zero;

    /// <summary>hom(A, B) on the given vertices: the graph whose edges are the nonzero off-diagonal entries
    /// of A or of B, its components, and the good ones (equal moduli on every edge, A_xx = B_xx, trivial
    /// holonomy of g_y/g_x = A_xy/B_xy). The candidates must contain every neighbour of a vertex; the entries
    /// decide which are edges.</summary>
    public static (int Components, int Good) Hom(IReadOnlyList<int> vertices,
        Func<int, int, GaussianInteger> a, Func<int, int, GaussianInteger> b, Func<int, IEnumerable<int>> candidates)
    {
        var adjacency = vertices.ToDictionary(x => x, x => candidates(x).Distinct()
            .Where(y => y != x && (a(x, y) != GaussianInteger.Zero || b(x, y) != GaussianInteger.Zero)).ToList());
        var seen = new HashSet<int>();
        int components = 0, good = 0;
        foreach (int start in vertices)
        {
            if (!seen.Add(start)) continue;
            components++;
            var component = new List<int> { start };
            var queue = new Queue<int>();
            queue.Enqueue(start);
            while (queue.Count > 0)
                foreach (int y in adjacency[queue.Dequeue()])
                    if (seen.Add(y)) { component.Add(y); queue.Enqueue(y); }
            if (IsGood(component, adjacency, a, b)) good++;
        }
        return (components, good);
    }

    private static BigInteger Norm(GaussianInteger g) => g.Re * g.Re + g.Im * g.Im;

    /// <summary>One component: walk it from its first vertex with g = 1, carrying g as a quotient num/den of
    /// Gaussian integers (g_y = g_x·A_xy/B_xy), then check the diagonal, the moduli and every edge by
    /// cross-multiplication.</summary>
    private static bool IsGood(List<int> component, Dictionary<int, List<int>> adjacency,
        Func<int, int, GaussianInteger> a, Func<int, int, GaussianInteger> b)
    {
        foreach (int x in component)
        {
            if (a(x, x) != b(x, x)) return false;
            foreach (int y in adjacency[x])
                if (Norm(a(x, y)) != Norm(b(x, y))) return false;
        }
        var num = new Dictionary<int, GaussianInteger> { [component[0]] = GaussianInteger.One };
        var den = new Dictionary<int, GaussianInteger> { [component[0]] = GaussianInteger.One };
        var queue = new Queue<int>();
        queue.Enqueue(component[0]);
        while (queue.Count > 0)
        {
            int x = queue.Dequeue();
            foreach (int y in adjacency[x])
            {
                if (num.ContainsKey(y)) continue;
                num[y] = num[x] * a(x, y);
                den[y] = den[x] * b(x, y);
                queue.Enqueue(y);
            }
        }
        // g_y·B_xy == g_x·A_xy  <=>  num_y·den_x·B_xy == num_x·den_y·A_xy
        foreach (int x in component)
            foreach (int y in adjacency[x])
                if (num[y] * den[x] * b(x, y) != num[x] * den[y] * a(x, y)) return false;
        return true;
    }

    /// <summary>Theorem 2 with every site dephased along Z: the near count (hopping components, every one
    /// good in hom(H, H)), the far count (good components of hom(H, H̄)) and the number of components of the
    /// union with the complement image.</summary>
    public static (int Near, int Far, int UnionComponents) EverySiteDephased(int n, Dictionary<int, GaussianInteger>[] h)
    {
        int all = (1 << n) - 1;
        var vertices = Enumerable.Range(0, 1 << n).ToList();
        var (_, near) = Hom(vertices, (x, y) => Entry(h, x, y), (x, y) => Entry(h, x, y), x => h[x].Keys);
        var (union, far) = Hom(vertices, (x, y) => Entry(h, x, y), (x, y) => Entry(h, x ^ all, y ^ all),
            x => h[x].Keys.Concat(h[x ^ all].Keys.Select(k => k ^ all)));
        return (near, far, union);
    }

    /// <summary>Real doubles scaled exactly to integers by one common power of two: every finite double is
    /// m·2^e with an integer m, so the scaled values are the inputs times 2^s with no rounding, and kernel
    /// dimensions do not see the common factor.</summary>
    public static BigInteger[] DyadicIntegers(IReadOnlyList<double> values)
    {
        int shift = 0;
        foreach (double v in values)
        {
            if (!double.IsFinite(v)) throw new ArgumentException($"a coefficient is not finite: {v}", nameof(values));
            if (v == 0) continue;
            long bits = BitConverter.DoubleToInt64Bits(v);
            int exponent = (int)((bits >> 52) & 0x7FF);
            long mantissa = bits & 0xFFFFFFFFFFFFFL;
            if (exponent == 0) exponent = 1; else mantissa |= 1L << 52;
            int e = exponent - 1075;                                // v = ±mantissa·2^e
            int trailing = BitOperations.TrailingZeroCount(mantissa);
            shift = Math.Max(shift, -(e + trailing));
        }
        return values.Select(v =>
        {
            if (v == 0) return BigInteger.Zero;
            long bits = BitConverter.DoubleToInt64Bits(v);
            int exponent = (int)((bits >> 52) & 0x7FF);
            long mantissa = bits & 0xFFFFFFFFFFFFFL;
            if (exponent == 0) exponent = 1; else mantissa |= 1L << 52;
            int e = exponent - 1075 + shift;
            var m = new BigInteger(mantissa) * (v < 0 ? -1 : 1);
            return e >= 0 ? m << e : m >> -e;                       // exact: shift was chosen so no bit is lost
        }).ToArray();
    }
}
