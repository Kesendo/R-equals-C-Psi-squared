using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the complement connection (claim <see cref="PalindromeComplementConnectionClaim"/>,
/// proof <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>): with one single-site Pauli jump on
/// every site, turn each jump to Z and the far end of F158 becomes the flat sections of
///
/// <code>
///     d_y / d_x = H_xy / H_x̄ȳ        (x̄ the bitwise complement)
/// </code>
///
/// on H's hopping graph over bitstrings. The near count is the number of components of the hopping graph
/// Γ_H, the far count the number of GOOD components of its union Γ with the complement image (equal
/// moduli on every edge, H_xx = H_x̄x̄, trivial holonomy), and the palindrome holds exactly when every
/// component of Γ is good (Theorem 2).
///
/// <para><b>A different road from the gate.</b> The sober base's <c>EndCount.ComplementConnection</c> is
/// the gate; this witness shares none of its code. It turns the letters with one proper rotation per
/// jump letter on <see cref="PauliMask"/>, reads each matrix entry off the mask form
/// P(x, z) = i^|x∧z|·X^x·Z^z (so ⟨y ⊕ x|P|y⟩ = i^|x∧z|·(−1)^|z∧y|), holds the section as a quotient of two
/// <see cref="GaussianInteger"/>s and checks every edge by cross-multiplication, with no division and no
/// float. Beside it, it reads <see cref="PalindromeStringSpanWitness"/> on the same row, which never forms
/// a basis state: two constructions meeting.</para>
///
/// <para><b>The theorems it reads.</b> The Hamiltonian is the one of <c>twoend</c> and
/// <c>twoendstrings</c>: Heisenberg bonds (weight 10) on a connected topology, letter fields of weight 3.
/// So every row lies in Theorem 1's class (one common axis, turned to Z by a global rotation) or in
/// Theorem 3's (mixed axes), and the witness prints the rule's prediction beside the graph's verdict:
/// under one axis a, the palindrome holds exactly when no field lies along a and the fields do not use
/// both other letters; under mixed axes, exactly when two axes occur and every field lies along the third
/// letter, the near end one-dimensional.</para>
///
/// <para>Args: <c>--N</c> (2..<see cref="MaxN"/>, default 3), <c>--deph</c> (a letter X/Y/Z per site,
/// every site dephased, default all Z), <c>--field</c> (X/Y/Z or '.' per site, default none),
/// <c>--topology</c> (chain|ring|complete, default chain).</para></summary>
public sealed class ComplementConnectionWitness : IInspectable
{
    /// <summary>A cost guard: the graph has 2^N vertices.</summary>
    public const int MaxN = 12;

    private const long Bond = 10, Field = 3;

    private readonly int _n;
    private readonly string _deph, _field, _topology;
    private readonly (int A, int B)[] _edges;

    public ComplementConnectionWitness(int n, string? deph = null, string? field = null, string? topology = null)
    {
        if (n < 2 || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"--root complement builds a graph on 2^N bitstrings and is guarded at N in 2..{MaxN}; got {n}. " +
                "The theorems carry no N, only this witness does.");
        _n = n;
        _deph = Letters(deph, n, 'Z', allowNone: false, nameof(deph));
        _field = Letters(field, n, '.', allowNone: true, nameof(field));
        _topology = (topology ?? "chain").ToLowerInvariant();
        _edges = _topology switch
        {
            "chain" => Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToArray(),
            "ring" => Enumerable.Range(0, n).Select(i => (i, (i + 1) % n)).ToArray(),
            "complete" => (from a in Enumerable.Range(0, n) from b in Enumerable.Range(a + 1, n - a - 1) select (a, b)).ToArray(),
            _ => throw new ArgumentException($"--topology takes chain, ring or complete; got \"{topology}\"."),
        };
    }

    private static string Letters(string? spec, int n, char fallback, bool allowNone, string arg)
    {
        if (string.IsNullOrWhiteSpace(spec)) return new string(fallback, n);
        string s = spec.Trim().ToUpperInvariant();
        if (s.Length != n)
            throw new ArgumentException($"--{arg} takes exactly one letter per site (N = {n}); got \"{s}\".");
        foreach (char c in s)
            if (!(c is 'X' or 'Y' or 'Z' || (allowNone && c is '.' or 'I')))
                throw new ArgumentException(allowNone
                    ? $"--{arg} letter '{c}' is not one of X Y Z or '.'."
                    : $"--{arg} letter '{c}' is not one of X Y Z: the complement connection needs one jump on EVERY " +
                      "site (Theorem 2's scope); an undephased site is read by --root twoendstrings.");
        return new string(s.Select(c => c == 'I' ? '.' : c).ToArray());
    }

    public string Name => "complement";

    public string DisplayName =>
        $"the complement connection live (N = {_n}, {_topology}, dephasing {_deph}, fields {_field})";

    public string Summary
    {
        get
        {
            var r = Read();
            return string.Format(CultureInfo.InvariantCulture,
                "hopping components {0} (near), good components {1} of {2} (far): {3}; {4} predicts {5}",
                r.HoppingComponents, r.Good, r.UnionComponents, r.Palindrome ? "palindrome" : "broken",
                r.TheoremClass, r.Predicted ? "palindrome" : "broken");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    /// <summary>What one inspect recomputes.</summary>
    public sealed record Reading(
        int HoppingComponents, int UnionComponents, int Good,
        bool Palindrome, string TheoremClass, bool Predicted);

    private Reading? _reading;

    public Reading Read()
    {
        if (_reading is not null) return _reading;
        var h = Entries();
        int dim = 1 << _n, all = dim - 1;
        var hop = new UnionFind(dim);
        var union = new UnionFind(dim);
        for (int x = 0; x < dim; x++)
            foreach (var y in h[x].Keys)
                if (y != x)
                {
                    hop.Join(x, y);
                    union.Join(x, y);
                    union.Join(x ^ all, y ^ all);
                }
        int good = 0;
        foreach (var component in Enumerable.Range(0, dim).GroupBy(union.Find))
            if (IsGood(component.ToList(), h, all)) good++;
        int hopping = Enumerable.Range(0, dim).Select(hop.Find).Distinct().Count();
        int unionCount = Enumerable.Range(0, dim).Select(union.Find).Distinct().Count();
        var (cls, predicted) = Prediction();
        return _reading = new Reading(hopping, unionCount, good, good == unionCount, cls, predicted);
    }

    // ---- the turn to Z, and the entries ----

    /// <summary>The proper rotation of one site's letters that sends its jump letter to Z, as
    /// (image letter, sign) for X, Y, Z: identity for Z; for X the Hadamard (X ↔ Z, Y → −Y); for Y the
    /// cycle X → Y → Z → X. Both are rotations (determinant +1), so they are conjugations by one-site
    /// unitaries and keep every kernel dimension.</summary>
    public static (char Letter, int Sign) Turn(char jump, char letter) => (jump, letter) switch
    {
        ('Z', _) => (letter, 1),
        ('X', 'X') => ('Z', 1),
        ('X', 'Z') => ('X', 1),
        ('X', 'Y') => ('Y', -1),
        ('Y', 'X') => ('Y', 1),
        ('Y', 'Y') => ('Z', 1),
        ('Y', 'Z') => ('X', 1),
        _ => throw new ArgumentException($"no turn for jump {jump}, letter {letter}"),
    };

    private IEnumerable<(PauliMask S, long C)> TurnedTerms()
    {
        var terms = new List<(char[] L, long C)>();
        foreach (var (a, b) in _edges)
            foreach (char p in "XYZ")
            {
                var l = Enumerable.Repeat('I', _n).ToArray();
                l[a] = l[b] = p;
                terms.Add((l, Bond));
            }
        for (int s = 0; s < _n; s++)
            if (_field[s] != '.')
            {
                var l = Enumerable.Repeat('I', _n).ToArray();
                l[s] = _field[s];
                terms.Add((l, Field));
            }
        foreach (var (l, c) in terms)
        {
            int sign = 1;
            var turned = new PauliLetter[_n];
            for (int s = 0; s < _n; s++)
            {
                if (l[s] == 'I') continue;
                var (t, g) = Turn(_deph[s], l[s]);
                sign *= g;
                turned[s] = PauliLetterExtensions.FromSymbol(t);
            }
            yield return (PauliMask.FromLetters(turned), sign * c);
        }
    }

    /// <summary>H in the turned basis, row by row: row y holds the column entries H_{y', y} keyed by y'.
    /// A string (x, z) sends |y⟩ to i^|x∧z|·(−1)^|z∧y|·|y ⊕ x⟩.</summary>
    private Dictionary<int, GaussianInteger>[] Entries()
    {
        int dim = 1 << _n;
        var h = Enumerable.Range(0, dim).Select(_ => new Dictionary<int, GaussianInteger>()).ToArray();
        foreach (var (s, c) in TurnedTerms())
        {
            int x = (int)s.X, z = (int)s.Z;
            int phase = BitOperations.PopCount((uint)(x & z)) & 3;
            for (int y = 0; y < dim; y++)
            {
                long v = (BitOperations.PopCount((uint)(z & y)) & 1) == 0 ? c : -c;
                GaussianInteger entry = phase switch
                {
                    0 => new(v, 0),
                    1 => new(0, v),
                    2 => new(-v, 0),
                    _ => new(0, -v),
                };
                int target = y ^ x;
                var row = h[y];
                row[target] = (row.TryGetValue(target, out var old) ? old : GaussianInteger.Zero) + entry;
            }
        }
        foreach (var row in h)
            foreach (var key in row.Where(kv => kv.Value == GaussianInteger.Zero).Select(kv => kv.Key).ToList())
                row.Remove(key);
        return h;
    }

    /// <summary>⟨r|H|c⟩: column c of H is stored as h[c], keyed by the row.</summary>
    private static GaussianInteger Entry(Dictionary<int, GaussianInteger>[] h, int r, int c) =>
        h[c].TryGetValue(r, out var v) ? v : GaussianInteger.Zero;

    private static BigInteger Norm(GaussianInteger g) => g.Re * g.Re + g.Im * g.Im;

    /// <summary>One component of Γ: walk it from its first vertex with d = 1, carrying d as a quotient
    /// num/den of Gaussian integers (d_y = d_x·a/b with a = H_xy, b = H_x̄ȳ), then check the diagonal and
    /// every edge by cross-multiplication. The neighbours of x in Γ are the rows of column x and the
    /// complements of the rows of column x̄.</summary>
    private static bool IsGood(List<int> component, Dictionary<int, GaussianInteger>[] h, int all)
    {
        foreach (int x in component)
            if (Entry(h, x, x) != Entry(h, x ^ all, x ^ all)) return false;
        var edges = new List<(int X, int Y, GaussianInteger A, GaussianInteger B)>();
        var adjacency = new Dictionary<int, List<int>>();
        foreach (int x in component)
        {
            var ys = h[x].Keys.Concat(h[x ^ all].Keys.Select(k => k ^ all)).Where(y => y != x).Distinct().ToList();
            foreach (int y in ys)
            {
                var a = Entry(h, x, y);
                var b = Entry(h, x ^ all, y ^ all);
                if (Norm(a) != Norm(b)) return false;
                edges.Add((x, y, a, b));
            }
            adjacency[x] = ys;
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
                num[y] = num[x] * Entry(h, x, y);
                den[y] = den[x] * Entry(h, x ^ all, y ^ all);
                queue.Enqueue(y);
            }
        }
        if (num.Count != component.Count) return false;
        // d_y·b == d_x·a  <=>  num_y·den_x·b == num_x·den_y·a
        foreach (var (x, y, a, b) in edges)
            if (num[y] * den[x] * b != num[x] * den[y] * a) return false;
        return true;
    }

    // ---- the theorems' rule ----

    private (string Class, bool Predicted) Prediction()
    {
        var axes = _deph.Distinct().ToList();
        var fields = _field.Where(c => c != '.').Distinct().ToList();
        if (axes.Count == 1)
        {
            char a = axes[0];
            var others = "XYZ".Where(c => c != a).ToList();
            bool pal = !fields.Contains(a) && !(fields.Contains(others[0]) && fields.Contains(others[1]));
            return ($"Theorem 1 (one common axis {a})", pal);
        }
        var free = "XYZ".Where(c => !axes.Contains(c)).ToList();
        return ($"Theorem 3 ({axes.Count} axes)", free.Count == 1 && fields.All(f => f == free[0]));
    }

    // ---- the children ----

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var r = Read();
            yield return new InspectableNode("the hopping graph, and its union with the complement image",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "Each jump turned to Z by a proper rotation of its site's letters; H read on the 2^N = {0} " +
                    "bitstrings exactly over the Gaussian integers. The hopping graph has {1} components (the " +
                    "near count, dim ker L), its union with the complement image {2}.",
                    1 << _n, r.HoppingComponents, r.UnionComponents));

            yield return new InspectableNode("the flat sections (Theorem 2)",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "{0} of the {1} union components are good: equal moduli on every edge, H_xx = H_x̄x̄, and a " +
                    "section d_y/d_x = H_xy/H_x̄ȳ that closes on every edge, checked by cross-multiplication with " +
                    "no division. That is the far count, dim ker(L + 2σ), and the verdict: {2}.",
                    r.Good, r.UnionComponents, r.Palindrome ? "every component good, the palindrome holds" : "broken"));

            yield return new InspectableNode("the theorem's rule, beside the graph",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "{0}: the rule predicts {1}, the graph reads {2}. {3}{4}",
                    r.TheoremClass, r.Predicted ? "a palindrome" : "a break", r.Palindrome ? "a palindrome" : "a break",
                    r.Predicted == r.Palindrome ? "They meet." : "THEY DO NOT MEET, READ IT.",
                    _deph.Distinct().Count() > 1
                        ? r.HoppingComponents == 1
                            ? " Under mixed axes Lemma A puts the near count at 1, and it is 1."
                            : $" Under mixed axes Lemma A puts the near count at 1; HERE IT IS {r.HoppingComponents}, READ IT."
                        : ""));

            yield return new InspectableNode("the string route, beside this one",
                summary: MeetTheStringRoute(r));

            yield return new InspectableNode("a break beside the canonical row",
                summary: ComparisonRow());

            yield return new InspectableNode("what this witness does NOT decide",
                summary: "Rows with an undephased site (Theorem 2 needs a jump on every site; --root twoendstrings " +
                         "reads them), jumps that are not single-site letters, and Hamiltonians other than the " +
                         "Heisenberg bonds with letter fields of twoend's family: Theorem 2 itself holds for any " +
                         "real combination of Pauli strings, but this witness builds only that family. The guard " +
                         "on N is the cost of 2^N vertices and carries no physics.");
        }
    }

    private string MeetTheStringRoute(Reading r)
    {
        var s = new PalindromeStringSpanWitness(_n, _deph, _field, _topology).Read();
        bool meet = s.NearUpper == r.HoppingComponents && s.FarUpper == r.Good
                    && PalindromeStringSpanWitness.IsPalindrome(s.Verdict) == r.Palindrome;
        return string.Format(CultureInfo.InvariantCulture,
            "twoendstrings counts the two ends as commutators on spans of Pauli strings and reads near {0}, far {1}, " +
            "{2}; the graph reads {3} and {4}. {5} The two share no construction: one ranks commutators on strings " +
            "and never forms a basis state, the other walks basis states and never forms a commutator.",
            s.NearUpper, s.FarUpper, s.Verdict, r.HoppingComponents, r.Good,
            meet ? "They meet." : "THEY DO NOT MEET, READ IT.");
    }

    private string ComparisonRow()
    {
        var canonical = new ComplementConnectionWitness(_n, new string('Z', _n), null, _topology).Read();
        string mixed = _n >= 3 ? "XYZ" + new string('Z', _n - 3) : "XY";
        var threeAxes = new ComplementConnectionWitness(_n, mixed, null, _topology).Read();
        return string.Format(CultureInfo.InvariantCulture,
            "Canonical: Z on every site, no field: hopping {0}, good {1} (the popcount shells, N + 1 = {2}). " +
            "Beside it, dephasing {3} with no field: hopping {4}, good {5}, {6}. {7}",
            canonical.HoppingComponents, canonical.Good, PalindromeTwoEndCountClaim.CanonicalChainCount(_n),
            mixed, threeAxes.HoppingComponents, threeAxes.Good, threeAxes.Palindrome ? "a palindrome" : "broken",
            (_n >= 3, threeAxes.Palindrome, threeAxes.HoppingComponents) switch
            {
                (true, false, 1) => "Three axes leave no letter to colour with, and the far end is empty (Theorem 3).",
                (false, true, 1) => "Two axes and no field: the third letter Z colours both sites (Theorem 3).",
                _ => "THIS IS NOT WHAT THEOREM 3 SAYS, READ IT.",
            });
    }

    private sealed class UnionFind
    {
        private readonly int[] _p;
        public UnionFind(int n) => _p = Enumerable.Range(0, n).ToArray();
        public int Find(int x)
        {
            while (_p[x] != x) x = _p[x] = _p[_p[x]];
            return x;
        }
        public void Join(int a, int b) => _p[Find(a)] = Find(b);
    }
}
