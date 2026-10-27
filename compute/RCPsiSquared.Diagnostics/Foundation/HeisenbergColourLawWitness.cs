using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The colour law live (claim <see cref="HeisenbergColourLawClaim"/>, F169, proof
/// <c>docs/proofs/PROOF_HEISENBERG_COLOURS_OFF_THE_LETTERS.md</c>). For each row (a graph, integer axes, Heisenberg
/// weights and integer fields) it computes the dimension of F158's far kernel, the operators W in the span of the
/// products ⊗(w_l·σ) with w_l ⊥ n_l that commute with H, and sets it beside the law's prediction.
///
/// <para><b>A different road from the proof and from the gate.</b> The proof and MirrorWorld's
/// <c>TiltedAxes.BondNullity</c> read the bond equations of Lemma 2's split; this witness does not split. It expands
/// every basis product into Pauli strings with integer coefficients, forms [H, W] with the bonds AND the fields in the
/// string basis (an anticommuting pair t, g contributes 2·t·g with its phase) and takes the rank of the columns over
/// GF(p), sending i to a square root of −1 modulo p. The rank modulo p is at most the rank over ℚ(i), so the count is an
/// upper bound on the far kernel's dimension; where the law predicts colours, those colours are elements, a lower
/// bound, and where it predicts none a count of 0 is exact. Agreement on a proven row is a check of the theorem; on a
/// row beyond the proven class (a site keeping three neighbours in the reduced core) it is a measurement.</para>
///
/// <para>The rows: a tree with a leaf on its neighbour's axis (coplanar axes off the letters, fields along their normal
/// and one field off it); the triangle generic and coplanar; rings of five and six sites with a right hinge at every
/// site (the only rings where the proof's last step has content); a four-ring under one axis, without fields and with
/// fields along one perpendicular; a four-chain under one axis with a uniform field along it, where the palindrome
/// breaks and the far kernel still keeps one element, so the law predicts no dimension there; and, measured, K4 generic and coplanar and the theta graph.</para></summary>
public sealed class HeisenbergColourLawWitness : IInspectable
{
    private const long P = 998_244_353;                    // 1 mod 4, so −1 has a square root
    private static readonly long I = ModPow(3, (P - 1) / 4, P);

    public sealed record Row(string Name, int N, (int A, int B)[] Edges, long[] Weights,
        (long X, long Y, long Z)[] Axes, (long X, long Y, long Z)[] Fields);

    /// <summary>Predicted is null where the law decides only the palindrome (a field along the common axis).</summary>
    public sealed record Reading(Row Row, bool Proven, int Cycles, int? Predicted, int Measured)
    {
        public bool Compared => Predicted is not null;
        public bool Agrees => Predicted is null || Predicted == Measured;
    }

    public static IReadOnlyList<Row> DefaultRows()
    {
        var none4 = new (long, long, long)[4];
        var coplanar = new (long, long, long)[] { (2, 2, 1), (1, 2, -1), (1, 0, 2), (1, 0, 2) };
        var tree = new[] { (0, 1), (1, 2), (2, 3) };
        var ring5 = Enumerable.Range(0, 5).Select(i => (i, (i + 1) % 5)).ToArray();
        var ring6 = Enumerable.Range(0, 6).Select(i => (i, (i + 1) % 6)).ToArray();
        var ring4 = Enumerable.Range(0, 4).Select(i => (i, (i + 1) % 4)).ToArray();
        var k4 = new[] { (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3) };
        var theta = new[] { (0, 2), (2, 1), (0, 3), (3, 1), (0, 4), (4, 1) };
        var generic = new (long, long, long)[] { (1, 2, 0), (0, 1, 2), (2, 0, 1), (1, -1, 2), (2, 1, -1) };
        long[] W(int k) => Enumerable.Range(0, k).Select(j => (long)(j % 3 == 0 ? 1 : j % 3 == 1 ? -2 : 3)).ToArray();
        return new[]
        {
            new Row("tree, leaf on its neighbour's axis, fields along the normal", 4, tree, W(3), coplanar,
                new (long, long, long)[] { (4, -3, -2), (0, 0, 0), (-8, 6, 4), (0, 0, 0) }),
            new Row("tree, one field off the normal", 4, tree, W(3), coplanar,
                new (long, long, long)[] { (4, -3, -2), (0, 1, 0), (0, 0, 0), (0, 0, 0) }),
            new Row("triangle, generic axes", 3, new[] { (0, 1), (1, 2), (2, 0) }, W(3), generic.Take(3).ToArray(), new (long, long, long)[3]),
            new Row("triangle, coplanar axes", 3, new[] { (0, 1), (1, 2), (2, 0) }, W(3), coplanar.Take(3).ToArray(), new (long, long, long)[3]),
            new Row("five-ring, a right hinge at every site", 5, ring5, W(5),
                new (long, long, long)[] { (0, 0, 1), (0, 1, -1), (2, -1, 1), (1, -2, -1), (1, 0, -1) }, new (long, long, long)[5]),
            new Row("six-ring, a right hinge at every site", 6, ring6, W(6),
                new (long, long, long)[] { (0, 0, 1), (0, 1, -1), (1, -1, 1), (0, 1, 1), (1, 2, 1), (1, 0, 1) }, new (long, long, long)[6]),
            new Row("four-ring under one axis, no field", 4, ring4, W(4), Enumerable.Repeat((1L, 1L, 0L), 4).ToArray(), none4),
            new Row("four-ring under one axis, fields along one perpendicular", 4, ring4, W(4), Enumerable.Repeat((1L, 1L, 0L), 4).ToArray(),
                new (long, long, long)[] { (1, -1, 0), (0, 0, 0), (-2, 2, 0), (0, 0, 0) }),
            new Row("four-chain under one axis, a uniform field along it (dimension not predicted)", 4, tree, W(3), Enumerable.Repeat((1L, 1L, 0L), 4).ToArray(),
                Enumerable.Repeat((1L, 1L, 0L), 4).ToArray()),
            new Row("K4, generic axes (measured)", 4, k4, W(6), generic.Take(4).ToArray(), none4),
            new Row("K4, coplanar axes (measured)", 4, k4, W(6), coplanar, none4),
            new Row("theta graph, generic axes (measured)", 5, theta, W(6), generic, new (long, long, long)[5]),
        };
    }

    private readonly IReadOnlyList<Row> _rows;

    public HeisenbergColourLawWitness(IReadOnlyList<Row>? rows = null) => _rows = rows ?? DefaultRows();

    public IReadOnlyList<Reading> Read() => _rows.Select(ReadRow).ToList();

    public static Reading ReadRow(Row r) => new(r, HeisenbergColourLawClaim.Proven(r.N, r.Edges, r.Axes),
        HeisenbergColourLawClaim.ReducedCoreCycles(r.N, r.Edges, r.Axes),
        HeisenbergColourLawClaim.PredictedFarDimension(r.Axes, r.Fields), FarDimension(r));

    /// <summary>The far kernel's dimension by the commutator in the string basis, an upper bound (rank over GF(p)).</summary>
    public static int FarDimension(Row r)
    {
        if (r.N < 1 || r.N > 6) throw new ArgumentOutOfRangeException(nameof(r), r.N, "the string basis has 4^N entries; guarded at N in 1..6");
        // H as strings: a string is one letter per site, 0 = I, 1 = X, 2 = Y, 3 = Z
        var h = new List<(int[] S, long C)>();
        for (int e = 0; e < r.Edges.Length; e++)
            for (int k = 1; k <= 3; k++)
            {
                var s = new int[r.N]; s[r.Edges[e].A] = k; s[r.Edges[e].B] = k;
                h.Add((s, r.Weights[e]));
            }
        for (int l = 0; l < r.N; l++)
            foreach (var (k, v) in new[] { (1, r.Fields[l].X), (2, r.Fields[l].Y), (3, r.Fields[l].Z) })
                if (v != 0) { var s = new int[r.N]; s[l] = k; h.Add((s, v)); }
        var perp = r.Axes.Select(PerpBasis).ToArray();
        int cols = 1 << r.N;
        var columns = new List<Dictionary<int, long>>();
        for (int col = 0; col < cols; col++)
        {
            // W = ⊗ w_l·σ, w_l the chosen perpendicular, expanded into strings
            var w = new Dictionary<int, long> { [0] = 1 };
            for (int l = 0; l < r.N; l++)
            {
                var v = perp[l][(col >> l) & 1];
                var next = new Dictionary<int, long>();
                foreach (var (code, c) in w)
                    foreach (var (k, m) in new[] { (1, v.X), (2, v.Y), (3, v.Z) })
                        if (m != 0) next[code + k * Pow4(l)] = Mod(c * Mod(m));
                w = next;
            }
            var commutator = new Dictionary<int, long>();
            foreach (var (t, ct) in h)
                foreach (var (code, cw) in w)
                {
                    // the product t·g site by site: a phase i^q and a string; t and g anticommute when the number of
                    // sites with two different non-identity letters is odd
                    int prod = 0, q = 0, clashes = 0;
                    for (int l = 0; l < r.N; l++)
                    {
                        int a = t[l], b = code / Pow4(l) % 4;
                        var (c, ph) = Letter(a, b);
                        prod += c * Pow4(l); q += ph;
                        if (a != 0 && b != 0 && a != b) clashes++;
                    }
                    if (clashes % 2 == 0) continue;
                    long phase = (q % 4) switch { 0 => 1, 1 => I, 2 => P - 1, _ => P - I };
                    long add = Mod(Mod(2 * Mod(ct)) * Mod(cw * phase % P));
                    commutator[prod] = Mod(commutator.GetValueOrDefault(prod) + add);
                }
            columns.Add(commutator);
        }
        return cols - Rank(columns);
    }

    // σ_a·σ_b = i^ph σ_c for letters 0..3 (I, X, Y, Z)
    private static (int C, int Ph) Letter(int a, int b)
    {
        if (a == 0) return (b, 0);
        if (b == 0 || a == b) return (a == b ? 0 : a, 0);
        int c = 6 - a - b;
        bool cyclic = (a, b) is (1, 2) or (2, 3) or (3, 1);
        return (c, cyclic ? 1 : 3);
    }

    private static int Rank(List<Dictionary<int, long>> columns)
    {
        var pivots = new Dictionary<int, Dictionary<int, long>>();
        int rank = 0;
        foreach (var original in columns)
        {
            var v = new Dictionary<int, long>(original.Where(kv => kv.Value != 0));
            while (v.Count > 0)
            {
                int lead = v.Keys.Min();
                if (!pivots.TryGetValue(lead, out var p)) { pivots[lead] = Normalise(v, lead); rank++; break; }
                long f = v[lead];
                foreach (var (k, x) in p) v[k] = Mod(v.GetValueOrDefault(k) - f * x % P);
                foreach (var k in v.Where(kv => kv.Value == 0).Select(kv => kv.Key).ToList()) v.Remove(k);
            }
        }
        return rank;
    }

    private static Dictionary<int, long> Normalise(Dictionary<int, long> v, int lead)
    {
        long inv = ModPow(v[lead], P - 2, P);
        return v.ToDictionary(kv => kv.Key, kv => kv.Value * inv % P);
    }

    private static int Pow4(int l) => 1 << (2 * l);
    private static long Mod(long x) => ((x % P) + P) % P;
    private static long ModPow(long b, long e, long m)
    {
        long r = 1; b %= m;
        while (e > 0) { if ((e & 1) == 1) r = r * b % m; b = b * b % m; e >>= 1; }
        return r;
    }

    private static (long X, long Y, long Z)[] PerpBasis((long X, long Y, long Z) n)
    {
        var e = Math.Abs(n.X) <= Math.Abs(n.Y) && Math.Abs(n.X) <= Math.Abs(n.Z) ? (1L, 0L, 0L)
            : Math.Abs(n.Y) <= Math.Abs(n.Z) ? (0L, 1L, 0L) : (0L, 0L, 1L);
        (long, long, long) X((long X, long Y, long Z) a, (long X, long Y, long Z) b)
            => (a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X);
        var u = X(n, e);
        return new[] { u, X(n, u) };
    }

    public string DisplayName => "F169, the colour law live: the far kernel by the commutator, against the law";

    public string Summary
    {
        get
        {
            var r = Read();
            return $"{r.Count(x => x.Compared && x.Agrees)}/{r.Count(x => x.Compared)} predicted rows agree ({r.Count(x => x.Proven)} proven, {r.Count(x => !x.Proven)} measured; {r.Count(x => !x.Compared)} with a field along the common axis, the palindrome broken and the dimension not predicted)";
        }
    }

    public IEnumerable<IInspectable> Children =>
        Read().Select(x => (IInspectable)new InspectableNode(x.Row.Name,
            summary: $"{(x.Proven ? "proven" : "measured")}, reduced core with {x.Cycles} cycle(s): law {(x.Predicted is { } p ? p.ToString() : "palindrome broken, dimension not predicted")}, " +
                     $"commutator {x.Measured}{(x.Agrees ? "" : "  DISAGREE")}"));

    public InspectablePayload Payload => InspectablePayload.Empty;
}
