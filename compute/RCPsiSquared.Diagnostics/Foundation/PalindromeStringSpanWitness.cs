using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F158's string route (claim <see cref="PalindromeTwoEndCountClaim"/>,
/// proof <c>docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md</c>, the second corollary to Lemma 1 and
/// section (f5)): the two ends of the spectrum counted as commutators on spans of Pauli strings,
/// without forming the Liouvillian.
///
/// <code>
///     dim ker L          = |dark| − rank ad_H|dark     (strings commuting with every jump)
///     dim ker(L + 2σ)    = |lit|  − rank ad_H|lit      (strings anticommuting with every jump)
/// </code>
///
/// <para><b>Why a second witness.</b> <see cref="PalindromeTwoEndCountWitness"/> ranks L itself on
/// 4^N columns and stops at N = 4. For Pauli-string jumps conjugation by a jump sends each string to
/// plus or minus itself, so each end is the commutator with H on a span of strings: 2^N strings under
/// one dephasing axis per site. Every entry of ad_H/2i there is plus or minus one coefficient of H,
/// since for a fixed string P distinct terms T give distinct products T·P. This witness takes the same
/// arguments and the same H as <c>twoend</c> (bonds XX + YY + ZZ at J = 1, fields of magnitude 3/10,
/// here as the integers 10 and 3, which moves neither kernel), runs to <see cref="MaxN"/>, and at
/// N ≤ 4 reads the dense witness beside itself, two routes that share no construction.</para>
///
/// <para><b>What is exact and what is a bound.</b> A rank mod p can only drop, so each count read
/// here is an UPPER bound on the rational one, the smaller over two primes kept. Single strings give
/// exact LOWER bounds: a dark string commuting with every term lies in ker L, and a lit one, a
/// COLOURING in the sense of <c>experiments/THE_PALINDROME_AS_A_COLOURING.md</c>, lies in
/// ker(L + 2σ). Three readings are exact: a colouring is an invertible element of the far space, so
/// by Lemma 3 it certifies the palindrome; a word in {H, A_1 .. A_m} with one jump letter and a
/// nonzero trace rules it out (section (f5)); and a far upper bound below the near lower bound rules
/// it out. Only between these do the two ranks decide, and the verdict says which reading decided.
/// The word traces are taken in the Pauli algebra over the Gaussian integers, Tr(H^k A) for
/// k ≤ <see cref="MaxWordPower"/>; a traceless budget decides nothing.</para>
///
/// <para>Args: <c>--N</c> (2..<see cref="MaxN"/>, default 3), <c>--deph</c> (a letter per site from
/// I/X/Y/Z or '.', default all Z), <c>--field</c> (same alphabet, default none), <c>--topology</c>
/// (chain|ring|complete, default chain), the alphabet and defaults of <c>twoend</c>.</para></summary>
public sealed class PalindromeStringSpanWitness : IInspectable
{
    /// <summary>A cost guard: the spans are enumerated string by string.</summary>
    public const int MaxN = 12;

    /// <summary>The largest span, as a power of two; undephased sites double a span per site.</summary>
    public const int MaxSpanBits = 16;

    /// <summary>The largest power of H in a word.</summary>
    public const int MaxWordPower = 4;

    private static readonly long[] Primes = { 998244353L, 1004535809L };
    private const long Bond = 10, Field = 3;

    private readonly int _n;
    private readonly string _deph, _field, _topology;
    private readonly (PauliMask S, long C)[] _terms;
    private readonly PauliMask[] _jumps;

    public PalindromeStringSpanWitness(int n, string? deph = null, string? field = null, string? topology = null)
    {
        if (n < 2 || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"--root twoendstrings enumerates spans of Pauli strings and is guarded at N in 2..{MaxN}; got {n}. " +
                "The theorem carries no N, only this witness does.");
        _n = n;
        _deph = ParseLetters(deph, n, 'Z', nameof(deph));
        _field = ParseLetters(field, n, '.', nameof(field));
        _topology = (topology ?? "chain").ToLowerInvariant();
        var edges = _topology switch
        {
            "chain" => Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToArray(),
            "ring" => Enumerable.Range(0, n).Select(i => (i, (i + 1) % n)).ToArray(),
            "complete" => (from a in Enumerable.Range(0, n) from b in Enumerable.Range(a + 1, n - a - 1) select (a, b)).ToArray(),
            _ => throw new ArgumentException($"--topology takes chain, ring or complete; got \"{topology}\"."),
        };
        _terms = Terms(n, edges, _field);
        _jumps = Enumerable.Range(0, n).Where(l => _deph[l] != '.').Select(l => Site(n, l, _deph[l])).ToArray();
        if (_jumps.Length == 0)
            throw new ArgumentException("--deph names no jump; F158 needs at least one.", nameof(deph));
    }

    private static string ParseLetters(string? spec, int n, char defaultLetter, string argName)
    {
        if (string.IsNullOrWhiteSpace(spec)) return new string(defaultLetter, n);
        string s = spec.Trim().ToUpperInvariant();
        if (s.Length != n)
            throw new ArgumentException(
                $"--{argName} takes exactly one letter per site (N = {n}); got \"{s}\" of length {s.Length}. " +
                "Use one of I X Y Z or '.' for none, e.g. \"Z.Z\".");
        return new string(s.Select(c => c switch
        {
            'I' or '.' or '0' => '.',
            'X' or 'Y' or 'Z' => c,
            _ => throw new ArgumentException($"--{argName} letter '{c}' is not one of I X Y Z or '.'."),
        }).ToArray());
    }

    private static PauliMask Site(int n, int site, char letter)
    {
        var l = new PauliLetter[n];
        l[site] = PauliLetterExtensions.FromSymbol(letter);
        return PauliMask.FromLetters(l);
    }

    private static (PauliMask, long)[] Terms(int n, (int A, int B)[] edges, string field)
    {
        var sum = new Dictionary<PauliMask, long>();
        void Add(PauliMask s, long c) => sum[s] = (sum.TryGetValue(s, out long old) ? old : 0) + c;
        foreach (var (a, b) in edges)
            foreach (char p in "XYZ")
            {
                var l = new PauliLetter[n];
                l[a] = l[b] = PauliLetterExtensions.FromSymbol(p);
                Add(PauliMask.FromLetters(l), Bond);
            }
        for (int s = 0; s < n; s++)
            if (field[s] != '.') Add(Site(n, s, field[s]), Field);
        return sum.Where(kv => kv.Value != 0).Select(kv => (kv.Key, kv.Value)).ToArray();
    }

    public string Name => "twoendstrings";

    public string DisplayName =>
        $"F158 live on strings: the two ends as commutators on Pauli-string spans (N = {_n}, {_topology}, dephasing {_deph})";

    public string Summary
    {
        get
        {
            var r = Read();
            return string.Format(CultureInfo.InvariantCulture,
                "near {0}, far {1} (upper bounds; exact lower bounds {2}, {3}); {4}, {5}",
                r.NearUpper, r.FarUpper, r.NearLower, r.FarLower, r.Verdict,
                IsExact(r.Verdict) ? "an exact reading" : "a rank reading");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public enum Decision
    {
        PalindromeByColour,   // exact: a colouring is an invertible element of the far space
        BrokenByWord,         // exact: a word with one jump letter has a nonzero trace
        BrokenByCount,        // exact: the far upper bound is below the near lower bound
        PalindromeByRank,     // the two upper bounds agree
        BrokenByRank,         // the two upper bounds differ
    }

    public static bool IsExact(Decision d) =>
        d is Decision.PalindromeByColour or Decision.BrokenByWord or Decision.BrokenByCount;

    public static bool IsPalindrome(Decision d) => d is Decision.PalindromeByColour or Decision.PalindromeByRank;

    /// <summary>What one inspect recomputes. Every field is a live number.</summary>
    public sealed record Reading(
        int DarkStrings, int LitStrings,
        int NearUpper, int FarUpper, int NearLower, int FarLower,
        IReadOnlyList<string> Colourings,
        string? Word, GaussianInteger WordTrace,
        Decision Verdict);

    private Reading? _reading;

    public Reading Read()
    {
        if (_reading is not null) return _reading;
        var dark = Span(lit: false);
        var lit = Span(lit: true);
        int nearUp = int.MaxValue, farUp = int.MaxValue;
        var dc = Columns(dark);
        var lc = Columns(lit);
        foreach (long p in Primes)
        {
            nearUp = Math.Min(nearUp, dark.Count - SparseRank(dc, p));
            farUp = Math.Min(farUp, lit.Count - SparseRank(lc, p));
        }
        var colourings = lit.Where(CommutesTermwise).ToList();
        int nearLo = dark.Count(CommutesTermwise), farLo = colourings.Count;
        if (nearUp < nearLo || farUp < farLo)
            throw new InvalidOperationException(
                $"an upper bound fell below a lower one (near {nearUp} < {nearLo} or far {farUp} < {farLo})");

        var (word, trace) = farLo > 0 ? (null, GaussianInteger.Zero) : OddWord();
        Decision verdict =
            farLo > 0 ? Decision.PalindromeByColour
            : word is not null ? Decision.BrokenByWord
            : farUp < nearLo ? Decision.BrokenByCount
            : nearUp == farUp ? Decision.PalindromeByRank
            : Decision.BrokenByRank;
        return _reading = new Reading(dark.Count, lit.Count, nearUp, farUp, nearLo, farLo,
            colourings.Select(c => c.ToString(_n)).ToList(), word, trace, verdict);
    }

    private bool CommutesTermwise(PauliMask p) => _terms.All(t => PauliMask.Commute(t.S, p));

    // ---- the spans, as solution sets over GF(2) ----

    // A string is the 2N-bit vector x | z << N; the symplectic form with a jump a is the GF(2)
    // functional whose x part is a.Z and whose z part is a.X. Dark: every functional 0; lit: every 1.
    private List<PauliMask> Span(bool lit)
    {
        int width = 2 * _n;
        var rows = _jumps.Select(a => (F: a.Z | (a.X << _n), T: lit ? 1 : 0)).ToList();
        var pivots = new List<int>();
        int r = 0;
        for (int c = 0; c < width && r < rows.Count; c++)
        {
            int piv = -1;
            for (int i = r; i < rows.Count; i++)
                if (((rows[i].F >> c) & 1) == 1) { piv = i; break; }
            if (piv < 0) continue;
            (rows[r], rows[piv]) = (rows[piv], rows[r]);
            for (int i = 0; i < rows.Count; i++)
                if (i != r && ((rows[i].F >> c) & 1) == 1)
                    rows[i] = (rows[i].F ^ rows[r].F, rows[i].T ^ rows[r].T);
            pivots.Add(c);
            r++;
        }
        for (int i = r; i < rows.Count; i++)
            if (rows[i].T == 1) return new List<PauliMask>();   // no string anticommutes with every jump

        ulong particular = 0;
        for (int i = 0; i < r; i++)
            if (rows[i].T == 1) particular |= 1UL << pivots[i];
        var basis = new List<ulong>();
        for (int c = 0; c < width; c++)
        {
            if (pivots.Contains(c)) continue;
            ulong v = 1UL << c;
            for (int i = 0; i < r; i++)
                if (((rows[i].F >> c) & 1) == 1) v |= 1UL << pivots[i];
            basis.Add(v);
        }
        if (basis.Count > MaxSpanBits)
            throw new InvalidOperationException(
                $"a span of 2^{basis.Count} strings exceeds 2^{MaxSpanBits}: too many undephased sites for N = {_n}");

        ulong mask = (1UL << _n) - 1;
        var result = new List<PauliMask>(1 << basis.Count);
        ulong cur = particular;
        for (int k = 0; k < 1 << basis.Count; k++)
        {
            if (k > 0) cur ^= basis[BitOperations.TrailingZeroCount(k)];   // Gray code
            result.Add(new PauliMask(cur & mask, cur >> _n));
        }
        return result;
    }

    // ---- the commutator on a span, and its rank ----

    // The column of P holds (ad_H P)/2i: +h_T or -h_T at the row of T·P for each term T
    // anticommuting with P ([T,P] = 2 i^k T·P with k odd, and i^(k-1) = ±1).
    private List<Dictionary<int, long>> Columns(IReadOnlyList<PauliMask> span)
    {
        var rowOf = new Dictionary<PauliMask, int>();
        var cols = new List<Dictionary<int, long>>(span.Count);
        foreach (var p in span)
        {
            var col = new Dictionary<int, long>();
            foreach (var (t, h) in _terms)
            {
                if (PauliMask.Commute(t, p)) continue;
                var (prod, k) = PauliMask.Multiply(t, p);
                if (!rowOf.TryGetValue(prod, out int row)) { row = rowOf.Count; rowOf[prod] = row; }
                col[row] = k == 1 ? h : -h;
            }
            cols.Add(col);
        }
        return cols;
    }

    private static long Mod(long x, long p) { long m = x % p; return m < 0 ? m + p : m; }
    private static long MulMod(long a, long b, long p) => (long)((UInt128)(ulong)Mod(a, p) * (ulong)Mod(b, p) % (ulong)p);

    /// <summary>The rank over GF(p) of sparse vectors, each pivoted on its smallest surviving index
    /// against the pivot vectors found so far.</summary>
    private static int SparseRank(IReadOnlyList<Dictionary<int, long>> vectors, long p)
    {
        var pivots = new Dictionary<int, (int[] Idx, long[] Val)>();
        var work = new SortedDictionary<int, long>();
        foreach (var v in vectors)
        {
            work.Clear();
            foreach (var (c, x) in v) { long m = Mod(x, p); if (m != 0) work[c] = m; }
            while (work.Count > 0)
            {
                int c = work.Keys.First();
                long lead = work[c];
                if (!pivots.TryGetValue(c, out var piv))
                {
                    long inv = CrossFormCertificate.Inv(lead, p);
                    pivots[c] = (work.Keys.ToArray(), work.Values.Select(x => MulMod(x, inv, p)).ToArray());
                    break;
                }
                for (int j = 0; j < piv.Idx.Length; j++)
                {
                    int cc = piv.Idx[j];
                    long nv = Mod((work.TryGetValue(cc, out long old) ? old : 0) - MulMod(lead, piv.Val[j], p), p);
                    if (nv == 0) work.Remove(cc); else work[cc] = nv;
                }
            }
        }
        return pivots.Count;
    }

    // ---- the odd word, section (f5) ----

    private static Dictionary<PauliMask, GaussianInteger> Mul(Dictionary<PauliMask, GaussianInteger> a, Dictionary<PauliMask, GaussianInteger> b)
    {
        var o = new Dictionary<PauliMask, GaussianInteger>();
        foreach (var (s, x) in a)
            foreach (var (t, y) in b)
            {
                var (u, k) = PauliMask.Multiply(s, t);
                var c = x * y;
                for (int j = 0; j < k; j++) c *= GaussianInteger.I;
                o[u] = (o.TryGetValue(u, out var old) ? old : GaussianInteger.Zero) + c;
            }
        foreach (var key in o.Where(kv => kv.Value == GaussianInteger.Zero).Select(kv => kv.Key).ToList()) o.Remove(key);
        return o;
    }

    /// <summary>The first of Tr(A_i), Tr(H A_i), ..., Tr(H^MaxWordPower A_i) that is nonzero, in the
    /// order of fw.odd_word_obstruction's one-jump words, the jump named by its string: Tr(H^k A) = 2^N times the coefficient of A in
    /// H^k, since A·A = 1 and no other string times A is the identity.</summary>
    private (string? Word, GaussianInteger Trace) OddWord()
    {
        var h = _terms.ToDictionary(t => t.S, t => (GaussianInteger)(BigInteger)t.C);
        var power = new Dictionary<PauliMask, GaussianInteger> { [PauliMask.Identity] = GaussianInteger.One };
        BigInteger dim = BigInteger.One << _n;
        for (int k = 0; k <= MaxWordPower; k++)
        {
            if (k > 0) power = Mul(power, h);
            for (int i = 0; i < _jumps.Length; i++)
                if (power.TryGetValue(_jumps[i], out var c) && c != GaussianInteger.Zero)
                    return (k > 0 ? $"H^{k}·{_jumps[i].ToString(_n)}" : _jumps[i].ToString(_n), c * (GaussianInteger)dim);
        }
        return (null, GaussianInteger.Zero);
    }

    // ---- the children ----

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var r = Read();
            yield return new InspectableNode("the two ends, as ranks on string spans",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "{0} dark strings (commuting with every jump) and {1} lit ones (anticommuting with every " +
                    "jump), where the Liouvillian has 4^N = {2} columns. The commutator with H on each span has " +
                    "nullity {3} (near end, ker L) and {4} (far end, ker(L + 2 sigma)), the smaller over two " +
                    "primes, each an upper bound on the rational count.",
                    r.DarkStrings, r.LitStrings, BigInteger.Pow(4, _n), r.NearUpper, r.FarUpper));

            yield return new InspectableNode("exact lower bounds, and the colourings",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "{0} dark strings commute with H and {1} lit ones do; single strings commuting with H are " +
                    "exact elements, so these are lower bounds with no prime in them. The lit ones are the " +
                    "colourings: {2}. Multiplying by one colouring carries the first set onto the second, which " +
                    "is why the two bounds are equal whenever a colouring exists.",
                    r.NearLower, r.FarLower, r.Colourings.Count == 0 ? "none" : string.Join(", ", r.Colourings)));

            yield return new InspectableNode("the odd word (section (f5))",
                summary: r.FarLower > 0
                    ? "Not read: a colouring certifies the palindrome, and then every word with an odd number of jump " +
                      "letters is traceless."
                    : r.Word is null
                        ? $"Every Tr(H^k A_i) with k <= {MaxWordPower} vanishes. That decides nothing: a traceless " +
                          "budget is not a certificate of the palindrome."
                        : string.Format(CultureInfo.InvariantCulture,
                            "{0} has trace {1} + {2} i, exactly, over the Gaussian integers (H scaled to bonds {3}, " +
                            "fields {4}). An invertible element of the far space would commute with H and flip the " +
                            "jump, making this trace vanish, so no palindrome exists at any positive rates.",
                            r.Word, r.WordTrace.Re, r.WordTrace.Im, Bond, Field));

            yield return new InspectableNode("the verdict, and what decided it",
                summary: $"{r.Verdict}: " + (r.Verdict switch
                {
                    Decision.PalindromeByColour => "a colouring, an exact certificate.",
                    Decision.BrokenByWord => "an odd word with a nonzero trace, an exact certificate.",
                    Decision.BrokenByCount => "the far upper bound lies below the near lower bound, exact.",
                    Decision.PalindromeByRank => "no exact reading applies, and the two upper bounds agree: modular evidence.",
                    _ => "no exact reading applies, and the two upper bounds differ: modular evidence.",
                }));

            yield return new InspectableNode("the dense witness, beside this one",
                summary: MeetTheDenseWitness(r));

            yield return new InspectableNode("the canonical row, and a break beside it",
                summary: BuildComparisonRow());

            yield return new InspectableNode("what this witness does NOT decide",
                summary: "Palindromic rows beyond the colouring, whose far space holds only sums of strings (an " +
                         "anticommuting sum, a SWAP-type element): there the verdict is a rank reading, and exhibiting an element " +
                         "is a separate construction, done for the census in simulations/anticommuting_sum_census.py " +
                         "and in compute/MirrorWorld/EndCount.cs's CheckElement. Nor jumps that are not Pauli " +
                         "strings, where conjugation does not act diagonally on strings; the dense witness twoend " +
                         "covers those through L itself. The guards on N and on the span are the cost of " +
                         "enumeration and carry no physics.");
        }
    }

    /// <summary>Both verdicts at this N, side by side, so that a witness that had lost the ability to
    /// report a break would be visible rather than silently reassuring.</summary>
    private string BuildComparisonRow()
    {
        var canonical = new PalindromeStringSpanWitness(_n).Read();
        var broken = new PalindromeStringSpanWitness(_n, new string('Z', _n), "Z" + new string('.', _n - 1)).Read();
        return string.Format(CultureInfo.InvariantCulture,
            "Canonical Heisenberg chain, Z on every site, no field: near {0}, far {1} (N + 1 = {2}), {3}. Beside it " +
            "the SAME dephasing with one Z field on site 0: near {4}, far {5}, {6}{7}. The near end does not move " +
            "and the far one empties; both readings come from the same routine, so the second is the falsifier " +
            "the first would otherwise lack.",
            canonical.NearUpper, canonical.FarUpper, PalindromeTwoEndCountClaim.CanonicalChainCount(_n), canonical.Verdict,
            broken.NearUpper, broken.FarUpper, broken.Verdict,
            broken.Word is null ? "" : $" by the word {broken.Word}");
    }

    private string MeetTheDenseWitness(Reading r)
    {
        if (_n > PalindromeTwoEndCountWitness.MaxN)
            return string.Format(CultureInfo.InvariantCulture,
                "Not run: N = {0} is past the dense witness's MaxN = {1}, which is the point of this route.",
                _n, PalindromeTwoEndCountWitness.MaxN);
        var dense = new PalindromeTwoEndCountWitness(_n, _deph, _field, _topology).Read();
        bool meet = dense.NearCount == r.NearUpper && dense.FarCount == r.FarUpper;
        return string.Format(CultureInfo.InvariantCulture,
            "twoend ranks L on 4^N columns and reads dim ker L = {0}, dim ker(L + 2 sigma) = {1}; the string route " +
            "reads {2} and {3}. {4} The two share no construction: one builds L and never a span, the other " +
            "builds spans and never L.",
            dense.NearCount, dense.FarCount, r.NearUpper, r.FarUpper,
            meet ? "They meet." : "THEY DO NOT MEET, READ IT.");
    }
}
