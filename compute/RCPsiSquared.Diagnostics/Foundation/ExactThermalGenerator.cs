using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The centred generator 2L + c under cooling and heating, in exact Gaussian integers, for the witnesses of
/// <see cref="ThermalThirdMomentClaim"/>, <see cref="ConditionedFieldPairingClaim"/>, <see cref="FluxAreaFifthMomentClaim"/>,
/// <see cref="TransportSeventhMomentClaim"/>, <see cref="SeventhMomentDegreeFourClaim"/> and
/// <see cref="EqualTrianglePairingClaim"/>.
///
/// <para>L = −i[H, ·] + Σ_l (γ↓_l·D[σ⁻_l] + γ↑_l·D[σ⁺_l]), σ⁻ = |0⟩⟨1| with Z|0⟩ = |0⟩, c = Σ_l (γ↓_l + γ↑_l). The image of
/// every matrix unit |i⟩⟨j| is written down directly (the commutator with H row by row and column by column; cooling
/// at site l clears bit l of both indices, heating sets it), vectorized column-major, |i⟩⟨j| at i + d·j, site 0 the
/// most significant bit. This is not the road of the registry's gates <c>simulations/thermal_third_moment.py</c>, which
/// builds L row-major by Kronecker products in exact rationals, and <c>simulations/fifth_moment_flux_area.py</c>, the
/// same in floats on integer inputs. Every entry of 2L + c is a Gaussian integer, and the
/// trace of 2L + c is checked to be zero.</para></summary>
public static class ExactThermalGenerator
{
    public const int MaxN = 4;

    /// <summary>2L + c, dense, 4^N × 4^N.</summary>
    public static GaussianInteger[,] Build(IReadOnlyList<ThermalThirdMomentClaim.Term> terms, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        if (terms is null) throw new ArgumentNullException(nameof(terms));
        if (cool is null) throw new ArgumentNullException(nameof(cool));
        if (heat is null) throw new ArgumentNullException(nameof(heat));
        int n = cool.Count;
        if (n < 1 || n > MaxN) throw new ArgumentOutOfRangeException(nameof(cool), n, $"N must lie in 1..{MaxN}");
        if (heat.Count != n) throw new ArgumentException($"{heat.Count} heating rates for {n} sites", nameof(heat));
        if (cool.Any(r => r < 0)) throw new ArgumentOutOfRangeException(nameof(cool), "rates are nonnegative");
        if (heat.Any(r => r < 0)) throw new ArgumentOutOfRangeException(nameof(heat), "rates are nonnegative");
        int d = 1 << n, dim = d * d;
        var h = Hamiltonian(n, terms);
        var m = new GaussianInteger[dim, dim];
        var minusTwoI = new GaussianInteger(0, -2);
        var plusTwoI = new GaussianInteger(0, 2);
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
            {
                int col = i + d * j;
                // −2i·(H|i⟩⟨j| − |i⟩⟨j|H)
                for (int a = 0; a < d; a++) if (h[a, i] != GaussianInteger.Zero) m[a + d * j, col] += minusTwoI * h[a, i];
                for (int b = 0; b < d; b++) if (h[j, b] != GaussianInteger.Zero) m[i + d * b, col] += plusTwoI * h[j, b];
                for (int l = 0; l < n; l++)
                {
                    int mask = 1 << (n - 1 - l);
                    int bi = (i & mask) != 0 ? 1 : 0, bj = (j & mask) != 0 ? 1 : 0;
                    // cooling: 2γ·σ⁻|i⟩⟨j|σ⁺ − γ·(n|i⟩⟨j| + |i⟩⟨j|n) + γ, n = |1⟩⟨1|
                    long g = cool[l];
                    if (g != 0)
                    {
                        if (bi == 1 && bj == 1) m[(i ^ mask) + d * (j ^ mask), col] += new GaussianInteger(2 * g, 0);
                        m[col, col] += new GaussianInteger(g * (1 - bi - bj), 0);
                    }
                    // heating: 2γ'·σ⁺|i⟩⟨j|σ⁻ − γ'·((1 − n)|i⟩⟨j| + |i⟩⟨j|(1 − n)) + γ'
                    long u = heat[l];
                    if (u != 0)
                    {
                        if (bi == 0 && bj == 0) m[(i ^ mask) + d * (j ^ mask), col] += new GaussianInteger(2 * u, 0);
                        m[col, col] += new GaussianInteger(u * (bi + bj - 1), 0);
                    }
                }
            }
        var tr = Trace(m);
        if (tr != GaussianInteger.Zero) throw new InvalidOperationException($"Tr(2L + c) = {tr.Re} + {tr.Im}i, not zero");
        return m;
    }

    /// <summary>H, d × d, from Pauli strings, site 0 first.</summary>
    private static GaussianInteger[,] Hamiltonian(int n, IReadOnlyList<ThermalThirdMomentClaim.Term> terms)
    {
        int d = 1 << n;
        var h = new GaussianInteger[d, d];
        foreach (var t in terms)
        {
            if (t.Letters is null || t.Letters.Length != n) throw new ArgumentException($"'{t.Letters}' is not a string on {n} sites", nameof(terms));
            for (int col = 0; col < d; col++)
            {
                int row = col;
                var f = GaussianInteger.One;
                for (int l = 0; l < n; l++)
                {
                    var (image, factor) = Letter(t.Letters[l], row, n, l);
                    row = image;
                    f *= factor;
                }
                h[row, col] += new GaussianInteger(t.Coefficient, 0) * f;
            }
        }
        return h;
    }

    /// <summary>P on site l applied to the basis state x: the image and ⟨image|P|x⟩.</summary>
    private static (int Image, GaussianInteger Factor) Letter(char p, int x, int n, int l)
    {
        int mask = 1 << (n - 1 - l);
        bool one = (x & mask) != 0;
        return p switch
        {
            'I' => (x, GaussianInteger.One),
            'X' => (x ^ mask, GaussianInteger.One),
            'Y' => (x ^ mask, one ? new GaussianInteger(0, -1) : new GaussianInteger(0, 1)),   // Y|0⟩ = i|1⟩, Y|1⟩ = −i|0⟩
            'Z' => (x, one ? new GaussianInteger(-1, 0) : GaussianInteger.One),
            _ => throw new ArgumentException($"'{p}' is not a Pauli letter", nameof(p)),
        };
    }

    /// <summary>Site 1's ket-bra sector of the matrix unit at column-major index idx, on two sites only (site 1 the low
    /// bit; |i⟩⟨j| sits at i + 4j).</summary>
    public static (int Ket, int Bra) SiteOneSector(int idx)
    {
        if (idx < 0 || idx >= 16) throw new ArgumentOutOfRangeException(nameof(idx), idx, "two sites: 0..15");
        return ((idx % 4) & 1, (idx / 4) & 1);
    }

    public static GaussianInteger Trace(GaussianInteger[,] m)
    {
        var t = GaussianInteger.Zero;
        for (int i = 0; i < m.GetLength(0); i++) t += m[i, i];
        return t;
    }

    public static GaussianInteger[,] Multiply(GaussianInteger[,] a, GaussianInteger[,] b)
    {
        int n = a.GetLength(0), k = b.GetLength(1), inner = b.GetLength(0);
        var r = new GaussianInteger[n, k];
        for (int i = 0; i < n; i++)
            for (int t = 0; t < inner; t++)
            {
                var x = a[i, t];
                if (x == GaussianInteger.Zero) continue;
                for (int j = 0; j < k; j++)
                {
                    var y = b[t, j];
                    if (y != GaussianInteger.Zero) r[i, j] += x * y;
                }
            }
        return r;
    }

    /// <summary>Tr(M^k) for the odd k up to top, exact; an imaginary part throws, the generator being real in the
    /// Hermitian Pauli basis.</summary>
    public static IReadOnlyDictionary<int, BigInteger> OddPowerSums(GaussianInteger[,] m, int top)
    {
        var sums = new SortedDictionary<int, BigInteger>();
        var p = m;
        for (int k = 1; k <= top; k++)
        {
            if (k > 1) p = Multiply(p, m);
            if (k % 2 == 0) continue;
            var t = Trace(p);
            if (t.Im != 0) throw new InvalidOperationException($"Tr(M^{k}) has an imaginary part");
            sums[k] = t.Re;
        }
        return sums;
    }

    /// <summary>The coefficients of det(y − M), index k the coefficient of y^k (the last one 1), by Faddeev and
    /// LeVerrier; every division is checked to be exact.</summary>
    public static GaussianInteger[] CharPoly(GaussianInteger[,] m)
    {
        int n = m.GetLength(0);
        var c = new GaussianInteger[n + 1];
        c[n] = GaussianInteger.One;
        var mk = new GaussianInteger[n, n];
        for (int k = 1; k <= n; k++)
        {
            var am = Multiply(m, mk);
            for (int i = 0; i < n; i++) am[i, i] += c[n - k + 1];
            mk = am;
            var t = TraceOfProduct(m, mk);
            if (!BigInteger.Remainder(t.Re, k).IsZero || !BigInteger.Remainder(t.Im, k).IsZero)
                throw new InvalidOperationException($"Faddeev-LeVerrier: the trace at step {k} is not divisible by {k}");
            c[n - k] = new GaussianInteger(-t.Re / k, -t.Im / k);
        }
        return c;
    }

    /// <summary>Tr(A·B), without forming the product.</summary>
    public static GaussianInteger TraceOfProduct(GaussianInteger[,] a, GaussianInteger[,] b)
    {
        var t = GaussianInteger.Zero;
        int n = a.GetLength(0), k = a.GetLength(1);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < k; j++)
                if (a[i, j] != GaussianInteger.Zero && b[j, i] != GaussianInteger.Zero) t += a[i, j] * b[j, i];
        return t;
    }

    /// <summary>The principal submatrix on the given indices.</summary>
    public static GaussianInteger[,] Extract(GaussianInteger[,] m, IReadOnlyList<int> idx)
    {
        var r = new GaussianInteger[idx.Count, idx.Count];
        for (int a = 0; a < idx.Count; a++)
            for (int b = 0; b < idx.Count; b++) r[a, b] = m[idx[a], idx[b]];
        return r;
    }

    /// <summary>The polynomial with the given coefficients (index k for y^k) at y.</summary>
    public static GaussianInteger Evaluate(GaussianInteger[] coefficients, GaussianInteger y)
    {
        var r = GaussianInteger.Zero;
        for (int k = coefficients.Length - 1; k >= 0; k--) r = r * y + coefficients[k];
        return r;
    }
}
