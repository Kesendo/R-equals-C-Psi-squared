using System.Numerics;

namespace MirrorWorld;

// The end count (adopted 2026-09-30 from F158, docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md, whose
// typed carriers are PalindromeTwoEndCountClaim and the live witness `inspect --root twoend`). For
// L(rho) = -i[H, rho] + sum_l gamma_l (A_l rho A_l - rho) with H Hermitian, every jump A_l Hermitian
// and squaring to the identity and every rate strictly positive, the spectrum pairs about -sigma
// exactly when the two ENDS of the spectrum hold the same number of modes,
//
//     dim ker L = dim ker(L + 2 sigma),     ker L = N_,   ker(L + 2 sigma) = W_,
//
//     N_ = {M : [H, M] = 0, A_l M A_l = +M for every l},   W_ = {M : [H, M] = 0, A_l M A_l = -M for every l}
//
// (Lemma 1 of the proof). The world adopts the COUNT and reads it through Pauli strings. When every
// jump is a Pauli string, conjugation by a jump sends each string to plus or minus itself, so the
// fixed space of all the conjugations is spanned by the DARK strings (those commuting with every
// jump) and the negated space by the LIT ones (those anticommuting with every jump), and each end is
// the kernel of the commutator with H on its own span:
//
//     near = dim N_ = |dark| - rank ad_H|dark,     far = dim W_ = |lit| - rank ad_H|lit.
//
// The spans are solution sets of linear equations over GF(2), one per jump (the dark strings a
// subspace, the lit ones a coset of it, empty exactly when no string anticommutes with every jump),
// so under one dephasing axis per site each span has 2^N strings where the Liouvillian has 4^N
// columns. And the commutator is integer: for a term T of H and a string P that anticommute,
// [T, P]/2i is plus or minus the single string T·P, and for a fixed P distinct terms give distinct
// strings, so every entry of the matrix is plus or minus ONE coefficient of H, never a sum. The
// rank is taken over GF(p) at the two primes of ModP, sparse (ModP.SparseRank).
//
// That answers the premise the arc mirrorworld_what_is_missing set before any adoption: read through
// Block (Grading B), the count could only run where F1 already forces the two ends equal, since Block
// needs a number-conserving H with diagonal jumps. The string route has no such premise. It takes any
// H written as Pauli strings with integer coefficients (a rational H scaled by its common
// denominator, which moves neither kernel) and any Pauli-string jumps: off-axis, multi-site, two axes
// on one site, fields along any letter, the families where the criterion has content.
//
// What is exact and what is a bound. A rank mod p can only DROP, so each nullity read here is an
// UPPER bound on the rational one (the convention of Seed, Divisor and BlindSeat). Two exact LOWER
// bounds come from single strings: a dark string that commutes with every term of H lies in N_ (and
// commuting with H is commuting term by term, since the products T·P cannot cancel), and a lit
// string that does lies in W_ (and once one such lit string F exists, multiplying by F matches the
// two sets string for string, so the two bounds are equal). The single lit strings commuting
// with H are the COLOURINGS of experiments/THE_PALINDROME_AS_A_COLOURING.md, and each is invertible,
// so by Lemma 3 a colouring CERTIFIES the palindrome. The palindrome's absence is certified in two
// ways: by the counts, when the far end's upper bound falls below the near end's lower bound (an
// empty lit coset is the extreme case: far is exactly 0 while the identity keeps near >= 1), and by
// an ODD WORD, F158 section (f5): an invertible element of W_ commutes with H and flips every jump,
// so every word in {H, A_1 .. A_m} with an odd number of jump letters is traceless when the
// palindrome holds, and one nonzero trace rules it out at every positive rate. The words are those
// of fw.odd_word_obstruction (simulations/framework/diagnostics/f158_odd_word.py), in the same order,
// exact over the Gaussian integers; the budget is the caller's, and the default reads one jump letter
// (the Python default also reads three, whose count grows as the cube of the jump count). Over EVERY word the condition is also
// sufficient, the traces of all words fixing a finite-dimensional *-representation up to unitary
// equivalence; a finite budget decides only when it fires. Only between these exact readings do the
// two ranks decide, and then the verdict says so.
//
// An element found elsewhere can certify as well: CheckElement asks whether a real combination of
// strings is lit, commutes with H and squares to a nonzero multiple of the identity, all in exact
// integer arithmetic. The anticommuting sums of the colouring page (stage E's a·YYZ + b·ZXX, the
// defect cascade) square that way. The golden router's identity column G of F116 has its coefficients
// in Z[r] and does not fit CheckElement's integers; the tests split it as G0 + r·G1 and ask
// CommutesWithH of each part, r being irrational, with lit strings only and a square that is
// (1 + r^2)^N times the identity site by site.
//
// The rates never enter. F158 section (f3): inside the open orthant of positive rates both ends are
// rate-blind, so the object takes the jumps and no rates, and every jump it is handed is assumed to
// carry a positive rate (a zero rate removes its jump: hand over fewer jumps). A jump's sign is
// irrelevant, and the identity string is accepted as a jump (it is Hermitian and squares to one; its
// lit coset is empty and the object reports the break).
//
// Parent: the frame. F158's class is every Hermitian involution, and its proof needs neither the
// mirror group D4 nor the fold, so no mechanism of another object is consumed. The two it meets are
// DOCKED in the tests rather than inherited: MirrorGroup's R, the map rho -> rho·X^N, is the
// one-sided multiplication by the canonical element X^N of W_, and Router's identity column is an
// element of W_ that is no single string (MirrorGroup keeps the golden router deliberately outside,
// and most of W_ is not monomial, which is why MirrorGroup is not the parent). Hardness reads a
// break on the same H = Z_0 through the odd power sums of M (F117), a trace on another matrix; how
// the two traces relate is not worked out here.
//
// Deliberately outside: the frame theorem of the colouring page (which frame a given H needs), the
// anticommuting-sum grammar and its census, the similarity L ~ -L - 2 sigma of F158 section (f8), and
// the router's two-sided local W. The count decides the multiset; it does not build a palindromizer.
//
// Words, fenced at the door. An END here is an end of the SPECTRUM, the eigenvalues 0 and
// -2 sigma, not the two ends of Crack's road (the combs) nor the two ends of SpookyAction's chain.
// FAR is the left edge -2 sigma; the glossary's dark "far end" is the long-lived end of a mode's
// lifetime, the opposite end of another line, and F158 section (f1) warns that W_ holds the
// non-oscillating modes at the left edge, not "the fastest modes". DARK and LIT keep the world's
// idiom, unwatched and watched: a dark string is watched by no jump, a lit one by every jump (the
// strings, not states; BlindSeat's dark states are states). A COLOUR is the colouring page's (a
// letter per site, lit there, that every term of H accepts); F87's LinearSiteColoring is its
// sibling, a string that ANTIcommutes with H, and Mirror's "two colors" of a chain are its
// bipartite parity, unrelated.
public sealed class EndCount : GameObject
{
    /// <summary>The largest span, as a power of two: a span is enumerated string by string, so
    /// 2^MaxSpanBits strings is the most either end is asked to rank.</summary>
    public const int MaxSpanBits = 16;

    /// <summary>The largest coefficient magnitude. Every matrix entry is plus or minus one coefficient,
    /// so this is only the bound that keeps a negated coefficient inside long.</summary>
    public const long MaxCoefficient = 1L << 62;

    public int N { get; }
    readonly (PauliString S, long C)[] terms;
    readonly PauliString[] jumpStrings;

    List<PauliString>? dark, lit;
    (int Near, int Far)? upper;

    public EndCount(World world, int n, IReadOnlyList<(string Letters, long Coefficient)> hamiltonian,
        IReadOnlyList<string> jumps) : base(world)
    {
        if (n < 1 || n > PauliString.MaxSites)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"1 to {PauliString.MaxSites} sites");
        if (jumps.Count == 0)
            throw new ArgumentException("F158 needs at least one jump", nameof(jumps));
        N = n;
        terms = Combine(hamiltonian, n).ToArray();
        jumpStrings = jumps.Select(j => ParseSized(j, n)).ToArray();
    }

    // left: what the count itself produces.
    public override IReadOnlyList<string> Own => new[] { "near", "far", "colour", "word" };

    static PauliString ParseSized(string letters, int n)
    {
        if (letters.Length != n)
            throw new ArgumentException($"'{letters}' has {letters.Length} letters on {n} sites");
        return PauliString.Parse(letters);
    }

    static IEnumerable<(PauliString S, long C)> Combine(IReadOnlyList<(string Letters, long Coefficient)> terms, int n)
    {
        // summed exactly and bounded AFTER summing: two in-range coefficients on one string can
        // reach long.MinValue, whose negation in Columns would wrap
        var sum = new Dictionary<PauliString, BigInteger>();
        foreach (var (letters, c) in terms)
        {
            var s = ParseSized(letters, n);
            sum[s] = (sum.TryGetValue(s, out var old) ? old : BigInteger.Zero) + c;
        }
        foreach (var v in sum.Values)
            if (BigInteger.Abs(v) > MaxCoefficient)
                throw new ArgumentOutOfRangeException(nameof(terms), v, $"|coefficient| exceeds {MaxCoefficient}");
        return sum.Where(kv => !kv.Value.IsZero).Select(kv => (kv.Key, (long)kv.Value)).ToList();
    }

    // ---- the two spans, as solution sets over GF(2) ----

    /// <summary>The strings commuting with every jump.</summary>
    public IReadOnlyList<PauliString> DarkStrings() => dark ??= Span(lit: false);

    /// <summary>The strings anticommuting with every jump (possibly none).</summary>
    public IReadOnlyList<PauliString> LitStrings() => lit ??= Span(lit: true);

    // A string is the 2N-bit vector v = x | z << N; the symplectic form with a jump a is the GF(2)
    // functional whose x part is a.Z and whose z part is a.X. Dark: every functional 0; lit: every 1.
    List<PauliString> Span(bool lit)
    {
        int width = 2 * N;
        var rows = jumpStrings.Select(a => (F: a.Z | (a.X << N), T: lit ? 1 : 0)).ToList();
        var pivotCols = new List<int>();
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
            pivotCols.Add(c);
            r++;
        }
        for (int i = r; i < rows.Count; i++)
            if (rows[i].T == 1) return new List<PauliString>();     // 0 = 1: no string anticommutes with every jump

        ulong particular = 0;
        for (int i = 0; i < r; i++)
            if (rows[i].T == 1) particular |= 1UL << pivotCols[i];
        var basis = new List<ulong>();
        for (int c = 0; c < width; c++)
        {
            if (pivotCols.Contains(c)) continue;
            ulong v = 1UL << c;
            for (int i = 0; i < r; i++)
                if (((rows[i].F >> c) & 1) == 1) v |= 1UL << pivotCols[i];
            basis.Add(v);
        }
        if (basis.Count > MaxSpanBits)
            throw new InvalidOperationException(
                $"a span of 2^{basis.Count} strings exceeds 2^{MaxSpanBits}: too few jumps for this many sites");

        ulong mask = N == 64 ? ulong.MaxValue : (1UL << N) - 1;
        var result = new List<PauliString>(1 << basis.Count);
        ulong cur = particular;
        for (int k = 0; k < (1 << basis.Count); k++)
        {
            if (k > 0) cur ^= basis[BitOperations.TrailingZeroCount(k)];   // Gray code: one basis vector per step
            result.Add(new PauliString(cur & mask, cur >> N));
        }
        return result;
    }

    // ---- the two counts ----

    // The commutator on a span, column by column: the column of P holds (ad_H P)/2i, one entry
    // +h_T or -h_T per term T anticommuting with P, at the row of the string T·P.
    List<IReadOnlyDictionary<int, long>> Columns(IReadOnlyList<PauliString> span)
    {
        var rowOf = new Dictionary<PauliString, int>();
        var cols = new List<IReadOnlyDictionary<int, long>>(span.Count);
        foreach (var p in span)
        {
            var col = new Dictionary<int, long>();
            foreach (var (t, h) in terms)
            {
                if (PauliString.Commute(t, p)) continue;
                var (prod, k) = PauliString.Multiply(t, p);
                if (!rowOf.TryGetValue(prod, out int row)) { row = rowOf.Count; rowOf[prod] = row; }
                col[row] = k == 1 ? h : -h;          // [T,P] = 2 i^k T·P, k odd; divided by 2i: i^(k-1) = +1 or -1
            }
            cols.Add(col);
        }
        return cols;
    }

    /// <summary>The nullity of the commutator on the dark (near) and lit (far) spans at one prime:
    /// each an upper bound on the rational count.</summary>
    public (int Near, int Far) NullitiesAt(long p)
    {
        var d = DarkStrings();
        var l = LitStrings();
        return (d.Count - ModP.SparseRank(Columns(d), p), l.Count - ModP.SparseRank(Columns(l), p));
    }

    /// <summary>The counts as upper bounds: the smaller nullity over the two primes at each end.</summary>
    public (int Near, int Far) UpperCounts()
    {
        if (upper is { } u) return u;
        var d = DarkStrings();
        var l = LitStrings();
        var dc = Columns(d);
        var lc = Columns(l);
        int near = d.Count - ModP.SparseRank(dc);
        int far = l.Count - ModP.SparseRank(lc);
        upper = (near, far);
        return upper.Value;
    }

    bool CommutesTermwise(PauliString p) => terms.All(t => PauliString.Commute(t.S, p));

    /// <summary>The colourings: the lit strings commuting with H. Each lies in W_ and is invertible.</summary>
    public IReadOnlyList<PauliString> Colourings() => LitStrings().Where(CommutesTermwise).ToList();

    /// <summary>The counts as exact lower bounds, from single strings. Once one colouring F exists,
    /// multiplication by F carries the dark strings commuting with H one to one onto the lit ones
    /// that do, so the two bounds are then equal and Lemma 3's near >= far adds nothing to them.</summary>
    public (int Near, int Far) LowerCounts()
        => (DarkStrings().Count(CommutesTermwise), Colourings().Count);

    // ---- the odd word, F158 section (f5) ----

    public sealed record OddWord(string Word, BigInteger TraceRe, BigInteger TraceIm);

    // An operator as {string: Gaussian integer}.
    static Dictionary<PauliString, (BigInteger Re, BigInteger Im)> Mul(
        Dictionary<PauliString, (BigInteger Re, BigInteger Im)> a, Dictionary<PauliString, (BigInteger Re, BigInteger Im)> b)
    {
        var o = new Dictionary<PauliString, (BigInteger Re, BigInteger Im)>();
        foreach (var (s, (ar, ai)) in a)
            foreach (var (t, (br, bi)) in b)
            {
                var (u, k) = PauliString.Multiply(s, t);
                var (re, im) = TimesIPower(ar * br - ai * bi, ar * bi + ai * br, k);
                var (r0, i0) = o.TryGetValue(u, out var old) ? old : (BigInteger.Zero, BigInteger.Zero);
                o[u] = (r0 + re, i0 + im);
            }
        foreach (var key in o.Where(kv => kv.Value.Re.IsZero && kv.Value.Im.IsZero).Select(kv => kv.Key).ToList())
            o.Remove(key);
        return o;
    }

    static (BigInteger Re, BigInteger Im) TimesIPower(BigInteger re, BigInteger im, int k)
    {
        for (int j = 0; j < k; j++) (re, im) = (-im, re);
        return (re, im);
    }

    /// <summary>The shortest word with an odd number of jump letters whose trace is nonzero, or null
    /// (no conclusion). The words and their order are fw.odd_word_obstruction's: Tr(H^k A_i) for
    /// k &lt;= maxPower, then with maxJumps = 3 the words H^a A_i H^b A_j H^c A_k with a + b + c &lt;=
    /// maxPower over every ordered jump triple, sorted by total power and then by jump count. The
    /// default reads one jump letter; the triple words grow as m^3 and are asked for explicitly.</summary>
    public OddWord? Word(int maxPower = 4, int maxJumps = 1)
    {
        if (maxPower < 0) throw new ArgumentOutOfRangeException(nameof(maxPower));
        if (maxJumps != 1 && maxJumps != 3) throw new ArgumentOutOfRangeException(nameof(maxJumps), "1 or 3");
        var one = new Dictionary<PauliString, (BigInteger, BigInteger)> { [PauliString.Identity] = (1, 0) };
        var h = terms.ToDictionary(t => t.S, t => ((BigInteger)t.C, BigInteger.Zero));
        var powers = new List<Dictionary<PauliString, (BigInteger Re, BigInteger Im)>> { one };
        for (int k = 0; k < maxPower; k++) powers.Add(Mul(powers[^1], h));

        int m = jumpStrings.Length;
        var words = new List<(int Total, int Jumps, string Name, int[] Hp, int[] Js)>();
        for (int k = 0; k <= maxPower; k++)
            for (int i = 0; i < m; i++)
                words.Add((k, 1, k > 0 ? $"H^{k}·A{i}" : $"A{i}", new[] { k }, new[] { i }));
        if (maxJumps >= 3)
            for (int total = 0; total <= maxPower; total++)
                for (int a = 0; a <= total; a++)
                    for (int b = 0; b <= total; b++)
                    {
                        int c = total - a - b;
                        if (c < 0) continue;
                        for (int t0 = 0; t0 < m; t0++)
                            for (int t1 = 0; t1 < m; t1++)
                                for (int t2 = 0; t2 < m; t2++)
                                    words.Add((total, 3, $"H^{a}·A{t0}·H^{b}·A{t1}·H^{c}·A{t2}",
                                        new[] { a, b, c }, new[] { t0, t1, t2 }));
                    }

        BigInteger dim = BigInteger.One << N;
        // the partial products H^a A_t0 and H^a A_t0 H^b A_t1 recur across words, so each is built once
        var memo = new Dictionary<string, Dictionary<PauliString, (BigInteger Re, BigInteger Im)>>();
        Dictionary<PauliString, (BigInteger Re, BigInteger Im)> TimesJump(Dictionary<PauliString, (BigInteger Re, BigInteger Im)> x, int t)
            => Mul(x, new Dictionary<PauliString, (BigInteger, BigInteger)> { [jumpStrings[t]] = (1, 0) });
        foreach (var w in words.OrderBy(w => w.Total).ThenBy(w => w.Jumps))
        {
            // prefix = everything before the last jump letter; Tr(prefix · A) = 2^N · prefix[A],
            // since A·A = I with phase +1 and no other string times A is the identity.
            var prefix = powers[w.Hp[0]];
            string key = $"{w.Hp[0]}";
            for (int j = 0; j < w.Js.Length - 1; j++)
            {
                key += $"|{w.Js[j]}|{w.Hp[j + 1]}";
                if (!memo.TryGetValue(key, out var next))
                    memo[key] = next = Mul(TimesJump(prefix, w.Js[j]), powers[w.Hp[j + 1]]);
                prefix = next;
            }
            if (prefix.TryGetValue(jumpStrings[w.Js[^1]], out var coef) && !(coef.Re.IsZero && coef.Im.IsZero))
                return new OddWord(w.Name, coef.Re * dim, coef.Im * dim);
        }
        return null;
    }

    // ---- the symmetries of (H, jumps) ----

    /// <summary>A symmetry: a permutation of the sites composed with one proper rotation of the
    /// letter frame, applied to every site. Perm[l] is the site letter l moves to; Rotation[a] and
    /// Sign[a] give the image of X (a = 0), Y (1), Z (2) as sign·letter.</summary>
    public sealed record Symmetry(int[] Perm, int[] Rotation, int[] Sign)
    {
        public bool MovesSites => Perm.Select((t, l) => t != l).Any(x => x);
        public bool IsIdentity => !MovesSites && Rotation.SequenceEqual(new[] { 0, 1, 2 }) && Sign.All(s => s == 1);
        public override string ToString() =>
            $"sites [{string.Join(",", Perm)}], X->{(Sign[0] < 0 ? "-" : "")}{"XYZ"[Rotation[0]]} " +
            $"Y->{(Sign[1] < 0 ? "-" : "")}{"XYZ"[Rotation[1]]} Z->{(Sign[2] < 0 ? "-" : "")}{"XYZ"[Rotation[2]]}";
    }

    /// <summary>The 24 proper rotations that permute the three axes up to sign: signed permutation
    /// matrices of determinant +1.</summary>
    public static IReadOnlyList<(int[] Rotation, int[] Sign)> CubeRotations()
    {
        var list = new List<(int[], int[])>();
        foreach (var perm in new[] { new[] { 0, 1, 2 }, new[] { 0, 2, 1 }, new[] { 1, 0, 2 }, new[] { 1, 2, 0 }, new[] { 2, 0, 1 }, new[] { 2, 1, 0 } })
        {
            // the cyclic orders are the even permutations
            int parity = perm[1] == (perm[0] + 1) % 3 ? 1 : -1;
            for (int s = 0; s < 8; s++)
            {
                var sign = new[] { (s & 1) == 0 ? 1 : -1, (s & 2) == 0 ? 1 : -1, (s & 4) == 0 ? 1 : -1 };
                if (parity * sign[0] * sign[1] * sign[2] == 1) list.Add((perm, sign));
            }
        }
        return list;
    }

    /// <summary>The image of a string under a symmetry, with the product of the letter signs.</summary>
    static (PauliString S, int Sign) Apply(Symmetry g, PauliString p)
    {
        // letter index a: X = 0, Y = 1, Z = 2, read from the (x, z) bits of each site
        ulong x = 0, z = 0;
        int sign = 1;
        for (int l = 0; l < g.Perm.Length; l++)
        {
            int xb = (int)((p.X >> l) & 1), zb = (int)((p.Z >> l) & 1);
            if (xb == 0 && zb == 0) continue;
            int a = xb == 1 ? (zb == 1 ? 1 : 0) : 2;
            int b = g.Rotation[a];
            sign *= g.Sign[a];
            int t = g.Perm[l];
            if (b != 2) x |= 1UL << t;
            if (b != 0) z |= 1UL << t;
        }
        return (new PauliString(x, z), sign);
    }

    /// <summary>Every symmetry of (H, jumps) in the group of site permutations times the 24 cube
    /// rotations: H mapped onto itself coefficient by coefficient, the jump set onto itself up to the
    /// sign of each jump (a jump's sign never enters L). Exact, string by string; N ≤ 8.</summary>
    public IReadOnlyList<Symmetry> Symmetries()
    {
        if (N > 8) throw new InvalidOperationException("the site permutations are enumerated; N <= 8");
        var h = terms.ToDictionary(t => t.S, t => t.C);
        var jumpSet = jumpStrings.ToHashSet();
        var found = new List<Symmetry>();
        var rotations = CubeRotations();
        foreach (var perm in Permutations(N))
            foreach (var (rot, sgn) in rotations)
            {
                var g = new Symmetry(perm, rot, sgn);
                bool ok = true;
                foreach (var (t, c) in terms)
                {
                    var (img, s) = Apply(g, t);
                    if (!h.TryGetValue(img, out long c2) || c2 != s * c) { ok = false; break; }
                }
                if (ok)
                    foreach (var a in jumpStrings)
                        if (!jumpSet.Contains(Apply(g, a).S)) { ok = false; break; }
                if (ok) found.Add(g);
            }
        return found;
    }

    static IEnumerable<int[]> Permutations(int n)
    {
        var a = Enumerable.Range(0, n).ToArray();
        IEnumerable<int[]> Rec(int k)
        {
            if (k == n) { yield return (int[])a.Clone(); yield break; }
            for (int i = k; i < n; i++)
            {
                (a[k], a[i]) = (a[i], a[k]);
                foreach (var r in Rec(k + 1)) yield return r;
                (a[k], a[i]) = (a[i], a[k]);
            }
        }
        return Rec(0);
    }

    // ---- an element found elsewhere ----

    public sealed record ElementCheck(bool AllLit, bool CommutesWithH, BigInteger SquareScalar)
    {
        /// <summary>Lit, commuting with H and squaring to a nonzero multiple of the identity: an
        /// invertible element of W_, which by F158's Lemma 3 carries the palindrome.</summary>
        public bool Certifies => AllLit && CommutesWithH && !SquareScalar.IsZero;
    }

    /// <summary>Is the real combination sum c_P P lit, does it commute with H, and what multiple of
    /// the identity is its square (zero when the square is not a multiple of the identity)?</summary>
    public ElementCheck CheckElement(IReadOnlyList<(string Letters, long Coefficient)> element)
    {
        var e = Combine(element, N).ToList();
        if (e.Count == 0) throw new ArgumentException("the element is zero", nameof(element));
        bool allLit = e.All(x => jumpStrings.All(a => !PauliString.Commute(a, x.S)));
        bool commutes = CommutesWithH(e);
        var w = e.ToDictionary(x => x.S, x => ((BigInteger)x.C, BigInteger.Zero));
        var sq = Mul(w, w);
        BigInteger scalar = sq.Count == 1 && sq.TryGetValue(PauliString.Identity, out var s0) && s0.Im.IsZero ? s0.Re : 0;
        return new ElementCheck(allLit, commutes, scalar);
    }

    /// <summary>[H, sum c_P P] = 0, exactly.</summary>
    public bool CommutesWithH(IReadOnlyList<(string Letters, long Coefficient)> element)
        => CommutesWithH(Combine(element, N).ToList());

    bool CommutesWithH(List<(PauliString S, long C)> element)
    {
        var acc = new Dictionary<PauliString, BigInteger>();
        foreach (var (t, h) in terms)
            foreach (var (p, c) in element)
            {
                if (PauliString.Commute(t, p)) continue;
                var (prod, k) = PauliString.Multiply(t, p);
                BigInteger v = (BigInteger)h * c * (k == 1 ? 1 : -1);
                acc[prod] = (acc.TryGetValue(prod, out var old) ? old : 0) + v;
            }
        return acc.Values.All(v => v.IsZero);
    }

    // ---- the verdict ----

    public enum Reading
    {
        PalindromeByColour,   // exact: a colouring is an invertible element of W_
        BrokenByWord,         // exact: an odd word with a nonzero trace
        BrokenByCount,        // exact: the far end's upper bound is below the near end's lower bound
        PalindromeByRank,     // the two upper bounds agree (a rank reading, not a certificate)
        BrokenByRank,         // the two upper bounds differ (a rank reading, not a certificate)
        PalindromeByElement,  // exact: the far end's one kernel vector, lifted to the integers, certifies
        PalindromeByCount,    // exact: both counts met by lifted, exactly checked kernel vectors, and equal
    }

    public static bool IsPalindrome(Reading r) =>
        r is Reading.PalindromeByColour or Reading.PalindromeByRank or Reading.PalindromeByElement or Reading.PalindromeByCount;
    public static bool IsExact(Reading r) =>
        r is Reading.PalindromeByColour or Reading.BrokenByWord or Reading.BrokenByCount or Reading.PalindromeByElement
            or Reading.PalindromeByCount;

    // ---- the kernels, lifted from GF(p) ----

    /// <summary>The largest span whose kernel is lifted; the elimination is dense.</summary>
    public const int MaxLiftColumns = 1024;

    readonly Dictionary<bool, List<List<(string Letters, long Coefficient)>>> lifted = new();

    /// <summary>A basis of the commutator's kernel on the dark (far = false) or lit (far = true) span,
    /// lifted from GF(p) at the first prime: the reduced-echelon basis, one vector per free column,
    /// each entry lifted to the rationals by rational reconstruction and scaled to coprime integers,
    /// and KEPT only if it commutes with H exactly. The kept vectors lie in the span and are independent
    /// (each carries its own free column), so their number is an exact LOWER bound on that end's count,
    /// whatever the prime did: a lift is a guess and the exact check decides it. Empty past
    /// MaxLiftColumns.</summary>
    public IReadOnlyList<IReadOnlyList<(string Letters, long Coefficient)>> LiftedKernel(bool far)
    {
        if (lifted.TryGetValue(far, out var cached)) return cached;
        var span = far ? LitStrings() : DarkStrings();
        var kept = new List<List<(string, long)>>();
        lifted[far] = kept;
        int s = span.Count;
        if (s == 0 || s > MaxLiftColumns) return kept;
        var cols = Columns(span);
        long p = ModP.Primes[0];
        int rows = cols.Max(c => c.Count == 0 ? 0 : c.Keys.Max() + 1);
        var m = new long[rows][];
        for (int r = 0; r < rows; r++) m[r] = new long[s];
        for (int j = 0; j < s; j++)
            foreach (var (r, x) in cols[j]) m[r][j] = ModP.Mod(x, p);

        var pivotOfRow = new List<int>();
        int rank = 0;
        for (int c = 0; c < s && rank < rows; c++)
        {
            int piv = -1;
            for (int r = rank; r < rows; r++) if (m[r][c] != 0) { piv = r; break; }
            if (piv < 0) continue;
            (m[rank], m[piv]) = (m[piv], m[rank]);
            long inv = ModP.ModInverse(m[rank][c], p);
            for (int j = c; j < s; j++) m[rank][j] = ModP.MulMod(m[rank][j], inv, p);
            for (int r = 0; r < rows; r++)
            {
                if (r == rank || m[r][c] == 0) continue;
                long f = m[r][c];
                for (int j = c; j < s; j++) m[r][j] = ModP.Mod(m[r][j] - ModP.MulMod(f, m[rank][j], p), p);
            }
            pivotOfRow.Add(c);
            rank++;
        }
        foreach (int fc in Enumerable.Range(0, s).Except(pivotOfRow))
        {
            var v = new long[s];
            v[fc] = 1;
            for (int i = 0; i < rank; i++) v[pivotOfRow[i]] = ModP.Mod(-m[i][fc], p);
            var fracs = new (BigInteger Num, BigInteger Den)[s];
            bool ok = true;
            for (int j = 0; j < s && ok; j++)
            {
                var q = Reconstruct(v[j], p);
                if (q is null) ok = false; else fracs[j] = q.Value;
            }
            if (!ok) continue;
            BigInteger lcm = fracs.Aggregate(BigInteger.One, (a, f) => a / BigInteger.GreatestCommonDivisor(a, f.Den) * f.Den);
            var ints = fracs.Select(f => f.Num * (lcm / f.Den)).ToArray();
            BigInteger g = ints.Aggregate(BigInteger.Zero, (a, x) => BigInteger.GreatestCommonDivisor(a, x));
            if (ints.Any(x => BigInteger.Abs(x / g) > MaxCoefficient)) continue;
            var vec = new List<(string, long)>();
            for (int j = 0; j < s; j++)
                if (!ints[j].IsZero) vec.Add((span[j].ToString(N), (long)(ints[j] / g)));
            if (CommutesWithH(vec)) kept.Add(vec);
        }
        return kept;
    }

    /// <summary>The exact lower bounds with the lifted kernels counted in: at each end the larger of
    /// the single-string bound and the number of lifted vectors that passed the exact check.</summary>
    public (int Near, int Far) LiftedLowerCounts()
    {
        var (n0, f0) = LowerCounts();
        return (Math.Max(n0, LiftedKernel(far: false).Count), Math.Max(f0, LiftedKernel(far: true).Count));
    }

    /// <summary>When the far end's upper bound is 1, its one lifted vector together with its exact
    /// check (lit, commuting with H, square a multiple of 1): a certifying check makes it an exact
    /// invertible element of W_. With near = 1 any element of W_ squares to a multiple of the identity,
    /// since its square lies in N_ = span{1}. Null when the bound is not 1 or no vector survived.</summary>
    public (IReadOnlyList<(string Letters, long Coefficient)> Element, ElementCheck Check)? FarElement()
    {
        if (UpperCounts().Far != 1) return null;
        var k = LiftedKernel(far: true);
        if (k.Count != 1) return null;
        return (k[0], CheckElement(k[0]));
    }

    // a/b with |a|, b <= sqrt(p/2) and a = b*x mod p, by the extended Euclidean algorithm; null if none
    static (BigInteger Num, BigInteger Den)? Reconstruct(long x, long p)
    {
        BigInteger bound = new BigInteger(Math.Sqrt(p / 2.0));
        BigInteger r0 = p, r1 = x, t0 = 0, t1 = 1;
        while (r1 > bound)
        {
            BigInteger q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (t0, t1) = (t1, t0 - q * t1);
        }
        if (t1.IsZero || BigInteger.Abs(t1) > bound) return null;
        return t1.Sign < 0 ? (-r1, -t1) : (r1, t1);
    }

    /// <summary>The exact readings first (colouring, odd word, counts, the lifted kernels), the ranks
    /// only where none applies. With a colouring no word is read: F158 makes every odd word traceless
    /// then, so the budget would be spent to learn nothing (the tests hold the two against each other
    /// instead). lift = false skips the lifted kernels, leaving the single-string bounds. An upper bound
    /// below a lower one would contradict the one-sidedness, and throws.</summary>
    public Reading Verdict(int maxPower = 4, int maxJumps = 1, bool lift = true)
    {
        var (upNear, upFar) = UpperCounts();
        var (loNear, loFar) = LowerCounts();
        if (upNear < loNear || upFar < loFar)
            throw new InvalidOperationException(
                $"an upper bound fell below a lower one (near {upNear} < {loNear} or far {upFar} < {loFar})");
        if (loFar > 0) return Reading.PalindromeByColour;
        if (Word(maxPower, maxJumps) is not null) return Reading.BrokenByWord;
        if (upFar < loNear) return Reading.BrokenByCount;
        if (!lift) return upNear == upFar ? Reading.PalindromeByRank : Reading.BrokenByRank;
        // the ranks alone would decide from here; the lifted kernels give exact lower bounds first
        var (liftNear, liftFar) = LiftedLowerCounts();
        if (upFar < liftNear) return Reading.BrokenByCount;
        if (upNear != upFar) return Reading.BrokenByRank;
        if (FarElement() is { Check.Certifies: true }) return Reading.PalindromeByElement;
        if (liftNear == upNear && liftFar == upFar) return Reading.PalindromeByCount;
        return Reading.PalindromeByRank;
    }
}
