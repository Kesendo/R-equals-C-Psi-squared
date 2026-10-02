using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The triangle's palindromic locus at N = 3, live (claim <see cref="PalindromeComplementConnectionClaim"/>,
/// section "The triangle's palindromic locus at N = 3" of <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>).
/// Heisenberg triangle on (d, u, v) = (0, 1, 2) with couplings J_du, J_dv, J_uv (at most one zero), the jump Z on d
/// alone, a single-letter field h_l·σ_{w_l} on every site. A row pairs exactly when it lies in
///
/// <code>
///     A  a colouring: every nonzero field along one letter p in {X, Y};
///     B  J_du = J_dv and a pi rotation about an axis n ⊥ z with R_n f_d = f_d, R_n f_u = f_v;
///     C  e₂ = J_duJ_dv + J_duJ_uv + J_dvJ_uv = 0, f_d = 0, f_u and f_v ⊥ z, |f_u| = |f_v|.
/// </code>
///
/// <para><b>A road of its own.</b> The proof's necessity is a Python elimination with a stored certificate
/// (<c>simulations/n3_palindrome_covering.py</c>) and a Python gate. This witness shares neither: it builds the
/// two ends of F158 itself, as the real kernels of W ↦ [H, W] on the Pauli strings whose d-letter is in {I, Z}
/// (near) or {X, Y} (far), with <see cref="BigRational"/> entries and an exact nullspace, and evaluates the three
/// classes on the point. It recomputes the rule at the row asked for and on a fixed grid of special loci over all
/// 27 words (equal, opposite and reciprocal couplings, one coupling zero; equal, opposite and zero fields), with
/// targeted points of C. It is evidence that the rule holds where it reads, not the exhaustion, which stays the
/// certificate's.</para>
///
/// <para>Args: <c>--word</c> (three letters X/Y/Z, one per site d, u, v; default ZXY), <c>--J</c> (J_du,J_dv,J_uv,
/// integers or p/q; default 3,6,-2), <c>--h</c> (h_d,h_u,h_v; default 0,1,1).</para></summary>
public sealed class TrianglePalindromeLocusWitness : IInspectable
{
    private readonly string _word;
    private readonly BigRational[] _j, _h;

    public TrianglePalindromeLocusWitness(string? word = null, string? couplings = null, string? fields = null)
    {
        _word = (word ?? "ZXY").Trim().ToUpperInvariant();
        if (_word.Length != 3 || _word.Any(c => c is not ('X' or 'Y' or 'Z')))
            throw new ArgumentException($"--word takes three letters from X, Y, Z (sites d, u, v); got \"{word}\".");
        _j = ParseTriple(couplings ?? "3,6,-2", "J");
        _h = ParseTriple(fields ?? "0,1,1", "h");
        if (_j.Count(x => x.IsZero) > 1)
            throw new ArgumentException("--J leaves the triangle disconnected (two couplings zero); the locus is stated for connected rows.");
    }

    private static BigRational[] ParseTriple(string spec, string arg)
    {
        var parts = spec.Split(',', StringSplitOptions.TrimEntries);
        if (parts.Length != 3) throw new ArgumentException($"--{arg} takes three comma-separated values; got \"{spec}\".");
        return parts.Select(p =>
        {
            var q = p.Split('/');
            var num = BigInteger.Parse(q[0], CultureInfo.InvariantCulture);
            var den = q.Length > 1 ? BigInteger.Parse(q[1], CultureInfo.InvariantCulture) : BigInteger.One;
            return new BigRational(num, den);
        }).ToArray();
    }

    // ---- the two ends, exactly --------------------------------------------------------------------------

    private static (int Phase, int Letter) Mul1(int a, int b)
    {
        if (a == 0) return (0, b);
        if (b == 0) return (0, a);
        if (a == b) return (0, 0);
        bool cyclic = (a, b) is (1, 2) or (2, 3) or (3, 1);
        return (cyclic ? 1 : 3, 6 - a - b);
    }

    private static bool Anticommute(int[] s, int[] t)
    {
        int n = 0;
        for (int i = 0; i < 3; i++) if (s[i] != 0 && t[i] != 0 && s[i] != t[i]) n++;
        return (n & 1) == 1;
    }

    /// <summary>The terms of H: (coefficient, Pauli string over 0 = I, 1 = X, 2 = Y, 3 = Z).</summary>
    private static List<(BigRational C, int[] S)> Terms(string word, BigRational[] j, BigRational[] h)
    {
        var terms = new List<(BigRational, int[])>();
        var bonds = new[] { (0, 1), (0, 2), (1, 2) };
        for (int b = 0; b < 3; b++)
            for (int c = 1; c <= 3; c++)
            {
                var s = new int[3]; s[bonds[b].Item1] = c; s[bonds[b].Item2] = c;
                terms.Add((j[b], s));
            }
        for (int l = 0; l < 3; l++)
        {
            var s = new int[3]; s[l] = "XYZ".IndexOf(word[l]) + 1;
            terms.Add((h[l], s));
        }
        return terms;
    }

    /// <summary>dim of {W real on the strings with d-letter in {I, Z} (near) or {X, Y} (far) : [H, W] = 0}.
    /// [t, s] = 2·i^k·u for anticommuting strings with k odd; divided by 2i every entry is ± a coefficient.</summary>
    public static int End(string word, BigRational[] j, BigRational[] h, bool far)
    {
        var cols = new List<int[]>();
        for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
                for (int c = 0; c < 4; c++)
                    if ((a is 1 or 2) == far) cols.Add(new[] { a, b, c });
        var rows = new Dictionary<int, BigRational[]>();
        foreach (var (coef, t) in Terms(word, j, h))
        {
            if (coef.IsZero) continue;
            for (int col = 0; col < cols.Count; col++)
            {
                var s = cols[col];
                if (!Anticommute(t, s)) continue;
                int phase = 0, key = 0;
                for (int i = 0; i < 3; i++)
                {
                    var (p, letter) = Mul1(t[i], s[i]);
                    phase += p; key = key * 4 + letter;
                }
                if (!rows.TryGetValue(key, out var row)) rows[key] = row = Enumerable.Repeat(BigRational.Zero, cols.Count).ToArray();
                row[col] = row[col] + ((phase & 3) == 1 ? coef : -coef);
            }
        }
        if (rows.Count == 0) return cols.Count;
        var m = new BigRational[rows.Count, cols.Count];
        int r = 0;
        foreach (var row in rows.Values) { for (int c = 0; c < cols.Count; c++) m[r, c] = row[c]; r++; }
        return BigRationalLinearAlgebra.Nullspace(m).Count;
    }

    // ---- the three classes ------------------------------------------------------------------------------

    /// <summary>The class a row lies in (A, B or C, the first that holds), or null. <paramref name="without"/>
    /// leaves one class out, for the planted mutation in the tests.</summary>
    public static char? Class(string word, BigRational[] j, BigRational[] h, char? without = null)
    {
        var f = Enumerable.Range(0, 3).Select(l => Vec(word[l], h[l])).ToArray();
        foreach (int p in new[] { 0, 1 })
            if (without != 'A' && f.All(v => Enumerable.Range(0, 3).All(q => q == p || v[q].IsZero))) return 'A';
        var rotations = new Func<BigRational[], BigRational[]>[]
        {
            v => new[] { v[0], -v[1], -v[2] }, v => new[] { -v[0], v[1], -v[2] },
            v => new[] { v[1], v[0], -v[2] }, v => new[] { -v[1], -v[0], -v[2] },
        };
        if (without != 'B' && j[0] == j[1] && rotations.Any(R => Same(R(f[0]), f[0]) && Same(R(f[1]), f[2]))) return 'B';
        var e2 = j[0] * j[1] + j[0] * j[2] + j[1] * j[2];
        if (without != 'C' && e2.IsZero && f[0].All(x => x.IsZero) && f[1][2].IsZero && f[2][2].IsZero
            && Norm(f[1]) == Norm(f[2])) return 'C';
        return null;
    }

    private static BigRational[] Vec(char letter, BigRational h)
    {
        var v = new[] { BigRational.Zero, BigRational.Zero, BigRational.Zero };
        v["XYZ".IndexOf(letter)] = h;
        return v;
    }

    private static bool Same(BigRational[] a, BigRational[] b) => a.Zip(b).All(p => p.First == p.Second);
    private static BigRational Norm(BigRational[] v) => v.Aggregate(BigRational.Zero, (s, x) => s + x * x);

    // ---- one reading ------------------------------------------------------------------------------------

    public sealed record Reading(int Near, int Far, bool Palindrome, char? InClass, bool Agrees);

    public static Reading ReadRow(string word, BigRational[] j, BigRational[] h, char? without = null)
    {
        int near = End(word, j, h, false), far = End(word, j, h, true);
        bool pal = far >= 1 && near == far;
        var cls = Class(word, j, h, without);
        return new Reading(near, far, pal, cls, pal == (cls is not null));
    }

    private Reading? _reading;
    public Reading Read() => _reading ??= ReadRow(_word, _j, _h);

    // ---- the grid of special loci -----------------------------------------------------------------------

    /// <summary>Every connected row of the grid: the 27 words, couplings equal, opposite, reciprocal (e₂ = 0,
    /// alone and on the swap branch), one zero and generic, fields equal, opposite, zero and generic.</summary>
    public static IEnumerable<(string Word, BigRational[] J, BigRational[] H)> Grid()
    {
        BigRational R(long p, long q = 1) => new(p, q);
        var couplings = new[]
        {
            new[] { R(1), R(2), R(3) }, new[] { R(2), R(2), R(3) }, new[] { R(2), R(-2), R(3) },
            new[] { R(3), R(6), R(-2) }, new[] { R(1), R(-2), R(-2) }, new[] { R(2), R(2), R(-1) },
            new[] { R(0), R(2), R(3) }, new[] { R(2), R(5), R(-10, 7) },
        };
        var fields = new[]
        {
            new[] { R(0), R(1), R(1) }, new[] { R(0), R(1), R(-1) }, new[] { R(1), R(1), R(1) },
            new[] { R(0), R(2), R(3) }, new[] { R(2), R(0), R(0) }, new[] { R(1), R(2), R(2) },
        };
        foreach (var w in from a in "XYZ" from b in "XYZ" from c in "XYZ" select $"{a}{b}{c}")
            foreach (var j in couplings)
                foreach (var h in fields)
                    yield return (w, j, h);
    }

    public sealed record GridReading(int Rows, int Palindromic, int A, int B, int C, int Disagreements);

    public static GridReading ReadGrid(char? without = null)
    {
        int rows = 0, pal = 0, a = 0, b = 0, c = 0, bad = 0;
        foreach (var (w, j, h) in Grid())
        {
            var r = ReadRow(w, j, h, without);
            rows++;
            if (r.Palindrome) pal++;
            if (r.InClass == 'A') a++; else if (r.InClass == 'B') b++; else if (r.InClass == 'C') c++;
            if (!r.Agrees) bad++;
        }
        return new GridReading(rows, pal, a, b, c, bad);
    }

    // ---- IInspectable -----------------------------------------------------------------------------------

    public string DisplayName =>
        $"the triangle's palindromic locus at N = 3 (word {_word}, J = {string.Join(", ", _j)}, h = {string.Join(", ", _h)})";

    public string Summary
    {
        get
        {
            var r = Read();
            return string.Format(CultureInfo.InvariantCulture, "near {0}, far {1}: {2}; {3}; the rule {4}",
                r.Near, r.Far, r.Palindrome ? "palindrome" : "broken",
                r.InClass is { } k ? $"in class {k}" : "in none of A, B, C",
                r.Agrees ? "holds" : "FAILS");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var g = ReadGrid();
            yield return new InspectableNode("the grid of special loci, recomputed",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "{0} rows over the 27 words: {1} palindromic ({2} in A, {3} in B, {4} in C, first class that " +
                    "holds), {5} rows where the rule and the two ends disagree", g.Rows, g.Palindromic, g.A, g.B, g.C,
                    g.Disagreements), provenance: NodeProvenance.Live);
            yield return new InspectableNode("where the exhaustion lives",
                summary: "This witness reads the rule at points; that every palindromic connected row lies in A, B or " +
                         "C is the elimination's, certificate simulations/results/n3_palindrome_covering.json (1,305 " +
                         "far leaves, 0 open) and its independent gate simulations/n3_palindrome_covering_gate.py.");
        }
    }
}
