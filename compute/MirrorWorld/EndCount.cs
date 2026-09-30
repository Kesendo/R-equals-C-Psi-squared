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

    /// <summary>H's terms as the object holds them (merged, zeros dropped), in the order CliffordMap indexes.</summary>
    public IReadOnlyList<(string Letters, long Coefficient)> Terms => terms.Select(t => (t.S.ToString(N), t.C)).ToList();
    readonly (PauliString S, long C)[] terms;
    readonly PauliString[] jumpStrings;
    readonly World world;

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
        this.world = world;
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

    // ---- the far element a symmetry builds ----

    /// <summary>The axis of a π rotation among the cube rotations, as integer components, or null when
    /// the rotation is not a π rotation. The nine π rotations are the three about the letter axes and
    /// the six about the bisectors of two letters.</summary>
    public static (long X, long Y, long Z)? HalfTurnAxis(int[] rotation, int[] sign)
    {
        var m = new long[3, 3];
        for (int a = 0; a < 3; a++) m[rotation[a], a] = sign[a];
        long trace = m[0, 0] + m[1, 1] + m[2, 2];
        if (trace != -1) return null;                          // a π rotation has trace 1 + 2 cos π = -1
        // its axis spans the kernel of M - 1; for a signed permutation of this kind it is a letter or a bisector
        var axis = new long[3];
        for (int a = 0; a < 3; a++)
            if (rotation[a] == a && sign[a] == 1) { axis[a] = 1; return (axis[0], axis[1], axis[2]); }
        for (int a = 0; a < 3; a++)
        {
            int b = rotation[a];
            if (b != a && rotation[b] == a && sign[a] == sign[b]) { axis[a] = 1; axis[b] = sign[a]; return (axis[0], axis[1], axis[2]); }
        }
        return null;
    }

    /// <summary>The far element a symmetry builds, when one does: a symmetry g = (σ, R) of (H, jumps)
    /// whose rotation R is the π rotation about an axis n, whose permutation σ is an involution fixing
    /// every site a jump acts on, and under which every jump goes to minus itself, gives
    ///
    ///     U_g = P_σ · (n·σ)^⊗N,
    ///
    /// with P_σ the permutation of tensor factors. U_g implements g, so Ad_U fixes H and negates every
    /// jump: U_g lies in W_, and it is invertible (a multiple of a unitary), so by F158 it carries the
    /// palindrome. Built exactly as a combination of strings with integer coefficients (2·SWAP_ij =
    /// 1 + X_iX_j + Y_iY_j + Z_iZ_j) and handed to CheckElement, so the construction is checked rather
    /// than trusted. Null when no symmetry of this kind exists.</summary>
    public (Symmetry G, IReadOnlyList<(string Letters, long Coefficient)> Element)? SymmetryElement()
    {
        foreach (var g in Symmetries())
        {
            var axis = HalfTurnAxis(g.Rotation, g.Sign);
            if (axis is null) continue;
            if (Enumerable.Range(0, N).Any(l => g.Perm[g.Perm[l]] != l)) continue;          // σ an involution
            bool negates = jumpStrings.All(a =>
            {
                for (int l = 0; l < N; l++)
                    if (((a.X >> l) & 1 | (a.Z >> l) & 1) != 0 && g.Perm[l] != l) return false;
                var (img, s) = Apply(g, a);
                return img == a && s == -1;
            });
            if (!negates) continue;

            var one = new Dictionary<PauliString, (BigInteger Re, BigInteger Im)> { [PauliString.Identity] = (1, 0) };
            var u = one;
            for (int i = 0; i < N; i++)
            {
                int j = g.Perm[i];
                if (j <= i) continue;
                var swap = new Dictionary<PauliString, (BigInteger Re, BigInteger Im)> { [PauliString.Identity] = (1, 0) };
                foreach (char c in "XYZ")
                {
                    var letters = Enumerable.Repeat('I', N).ToArray();
                    letters[i] = letters[j] = c;
                    swap[PauliString.Parse(new string(letters))] = (1, 0);
                }
                u = Mul(u, swap);
            }
            var (nx, ny, nz) = axis.Value;
            for (int l = 0; l < N; l++)
            {
                var site = new Dictionary<PauliString, (BigInteger Re, BigInteger Im)>();
                foreach (var (c, v) in new[] { ('X', nx), ('Y', ny), ('Z', nz) })
                {
                    if (v == 0) continue;
                    var letters = Enumerable.Repeat('I', N).ToArray();
                    letters[l] = c;
                    site[PauliString.Parse(new string(letters))] = (v, 0);
                }
                u = Mul(u, site);
            }
            // P_σ (σ an involution) is Hermitian and commutes with (n·σ)^⊗N, so no imaginary part survives
            if (u.Values.Any(c => !c.Im.IsZero) || u.Values.Any(c => BigInteger.Abs(c.Re) > MaxCoefficient)) continue;
            var element = u.Select(kv => (kv.Key.ToString(N), (long)kv.Value.Re)).ToList();
            if (CheckElement(element).Certifies) return (g, element);
        }
        return null;
    }

    // ---- every site dephased: the complement connection ----

    /// <summary>The image of a letter under the proper rotation that turns a site's jump letter into Z,
    /// with its sign: Z fixed; for X the Hadamard, X ↔ Z and Y → −Y; for Y the cyclic X → Y → Z → X.
    /// A proper rotation, so conjugation by a single-site unitary: the product table XY = iZ (and its
    /// cycles) is kept.</summary>
    public static (char Letter, int Sign) TurnToZ(char letter, char jump) => (jump, letter) switch
    {
        ('Z', _) => (letter, 1),
        ('X', 'X') => ('Z', 1), ('X', 'Z') => ('X', 1), ('X', 'Y') => ('Y', -1),
        ('Y', 'Y') => ('Z', 1), ('Y', 'Z') => ('X', 1), ('Y', 'X') => ('Y', 1),
        _ => throw new ArgumentException($"letters and jumps are X, Y or Z; got {letter} under {jump}"),
    };

    /// <summary>The reading of docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md when every site carries
    /// exactly one single-site jump: HoppingComponents = dim of the near end (the components of H's
    /// hopping graph on bitstrings in the jumps' eigenbasis), UnionComponents the components of that
    /// graph joined with its image under the complement x -> x̄, Good the union components on which the
    /// complement connection d_y / d_x = H_xy / H_{x̄ȳ} has a flat section (|H_xy| = |H_{x̄ȳ}| on every
    /// edge, H_xx = H_{x̄x̄} on every vertex, trivial holonomy), which is the dimension of the far end;
    /// AllGood, every union component good, is the palindrome. Null when some site is undephased or a
    /// jump acts on more than one site. Exact: the entries are Gaussian integers, the sections Gaussian
    /// rationals. N ≤ 10.</summary>
    public (int HoppingComponents, int UnionComponents, int Good, bool AllGood)? ComplementConnection()
    {
        if (N > 10) throw new InvalidOperationException("the hopping graph has 2^N vertices; N <= 10");
        var axis = new char[N];
        foreach (var a in jumpStrings)
        {
            var sites = Enumerable.Range(0, N).Where(l => Touches(a, l)).ToList();
            if (sites.Count != 1) return null;
            if (axis[sites[0]] != '\0') return null;
            axis[sites[0]] = a.Letter(sites[0]);
        }
        if (axis.Any(c => c == '\0')) return null;

        int d = 1 << N;
        var re = new Dictionary<(int, int), BigInteger>();
        var im = new Dictionary<(int, int), BigInteger>();
        foreach (var (t, c) in terms)
        {
            var letters = new char[N];
            int sign = 1;
            for (int l = 0; l < N; l++)
            {
                char ch = t.Letter(l);
                if (ch == 'I') { letters[l] = 'I'; continue; }
                var (nl, sg) = TurnToZ(ch, axis[l]);
                letters[l] = nl; sign *= sg;
            }
            var q = PauliString.Parse(new string(letters));
            int ny = System.Numerics.BitOperations.PopCount(q.X & q.Z);
            for (int col = 0; col < d; col++)
            {
                int row = col ^ (int)q.X;
                int k = (ny + 2 * System.Numerics.BitOperations.PopCount((ulong)col & q.Z)) & 3;
                BigInteger v = sign * (BigInteger)c;
                var key = (row, col);
                if (k == 0 || k == 2) re[key] = (re.TryGetValue(key, out var o) ? o : 0) + (k == 0 ? v : -v);
                else im[key] = (im.TryGetValue(key, out var o) ? o : 0) + (k == 1 ? v : -v);
            }
        }
        (BigInteger Re, BigInteger Im) H(int r, int c) =>
            (re.TryGetValue((r, c), out var a) ? a : 0, im.TryGetValue((r, c), out var b) ? b : 0);
        bool Zero((BigInteger Re, BigInteger Im) z) => z.Re.IsZero && z.Im.IsZero;
        int Comp(int x) => x ^ (d - 1);

        var nbr = new List<int>[d];
        for (int x = 0; x < d; x++) nbr[x] = new List<int>();
        foreach (var key in re.Keys.Concat(im.Keys).Distinct())
        {
            var (r, c) = key;
            if (r == c || Zero(H(r, c))) continue;
            nbr[r].Add(c); nbr[c].Add(r);                                  // H's hopping edge
            nbr[Comp(r)].Add(Comp(c)); nbr[Comp(c)].Add(Comp(r));          // and its complement image
        }
        int Components(bool union)
        {
            var seen = new bool[d];
            int count = 0;
            for (int s0 = 0; s0 < d; s0++)
            {
                if (seen[s0]) continue;
                count++;
                var stack = new Stack<int>(); stack.Push(s0); seen[s0] = true;
                while (stack.Count > 0)
                {
                    int x = stack.Pop();
                    foreach (int y in nbr[x])
                        if (!seen[y] && (union || !Zero(H(x, y)))) { seen[y] = true; stack.Push(y); }
                }
            }
            return count;
        }
        int hopping = Components(union: false);

        // the flat sections, component by component of the union graph, as Gaussian rationals (re, im) / den
        var sec = new (BigInteger Re, BigInteger Im, BigInteger Den)?[d];
        int unionCount = 0, good = 0;
        for (int s0 = 0; s0 < d; s0++)
        {
            if (sec[s0] is not null) continue;
            unionCount++;
            bool ok = true;
            sec[s0] = (1, 0, 1);
            var stack = new Stack<int>(); stack.Push(s0);
            while (stack.Count > 0)
            {
                int x = stack.Pop();
                var hxx = H(x, x); var hcc = H(Comp(x), Comp(x));
                if (hxx != hcc) ok = false;
                foreach (int y in nbr[x])
                {
                    var a = H(x, y); var b = H(Comp(x), Comp(y));
                    // equal moduli; also implied by the holonomy, since the edge is walked in both
                    // directions and the reverse ratio is the conjugate one
                    if (a.Re * a.Re + a.Im * a.Im != b.Re * b.Re + b.Im * b.Im || Zero(a)) { ok = false; }
                    // d_y = d_x · a / b = d_x · a · conj(b) / |b|^2
                    var (dr, di, dd) = sec[x]!.Value;
                    BigInteger nr = a.Re * b.Re + a.Im * b.Im, ni = a.Im * b.Re - a.Re * b.Im, n2 = b.Re * b.Re + b.Im * b.Im;
                    (BigInteger, BigInteger, BigInteger) v = n2.IsZero ? (0, 0, 1) : (dr * nr - di * ni, dr * ni + di * nr, dd * n2);
                    var gv = BigInteger.GreatestCommonDivisor(BigInteger.GreatestCommonDivisor(v.Item1, v.Item2), v.Item3);
                    if (!gv.IsZero && !gv.IsOne) v = (v.Item1 / gv, v.Item2 / gv, v.Item3 / gv);
                    if (sec[y] is null) { sec[y] = v; stack.Push(y); }
                    else
                    {
                        var (yr, yi, yd) = sec[y]!.Value;
                        if (yr * v.Item3 != v.Item1 * yd || yi * v.Item3 != v.Item2 * yd) ok = false;
                    }
                }
            }
            if (ok) good++;
        }
        return (hopping, unionCount, good, good == unionCount);
    }

    // ---- undephased sites with a conserved letter: the sector connection ----

    /// <summary>The reading of Theorem 4 of docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md: every
    /// dephased site carries exactly one single-site jump, and every undephased site u a letter P_u that H
    /// conserves (every term of H has I or P_u at u; a site H does not touch conserves Z). Turn each jump and
    /// each conserved letter to Z; H is then block diagonal in the undephased bits, H = ⊕_σ H_σ with each
    /// H_σ a matrix on the dephased bits, and with H̄_σ(x, y) = H_σ(x̄, ȳ) (the complement on the dephased
    /// bits only)
    ///
    ///     near = Σ_{σ,τ} hom(H_τ, H_σ),   far = Σ_{σ,τ} hom(H_τ, H̄_σ),
    ///
    /// hom(A, B) = dim{g : B_xy·g_y = g_x·A_xy}, the good components of the graph whose edges are the
    /// nonzero off-diagonal entries of A or B (|A_xy| = |B_xy| on every edge, A_xx = B_xx, trivial holonomy
    /// of g_y/g_x = A_xy/B_xy). Null when a jump is not a single letter on one site, a site carries two
    /// jumps, or an undephased site has no conserved letter. With no undephased site it is Theorem 2.
    /// Cross is the part of the two counts that pairs two different sectors (σ ≠ τ).
    /// Exact over the Gaussian rationals; N ≤ 10.</summary>
    public (int Near, int Far, int Sectors, int Cross)? SectorConnection()
    {
        if (N > 10) throw new InvalidOperationException("the graphs have 2^N vertices; N <= 10");
        var axis = new char[N];
        foreach (var a in jumpStrings)
        {
            var sites = Enumerable.Range(0, N).Where(l => Touches(a, l)).ToList();
            if (sites.Count != 1 || axis[sites[0]] != '\0') return null;
            axis[sites[0]] = a.Letter(sites[0]);
        }
        var free = Enumerable.Range(0, N).Where(l => axis[l] == '\0').ToList();
        foreach (int u in free)
        {
            var used = terms.Select(t => t.S.Letter(u)).Where(c => c != 'I').Distinct().ToList();
            if (used.Count > 1) return null;
            axis[u] = used.Count == 1 ? used[0] : 'Z';
        }
        var dephased = Enumerable.Range(0, N).Where(l => !free.Contains(l)).ToList();

        // H turned, as (row, col) -> Gaussian integer over all N bits
        int d = 1 << N;
        var h = new Dictionary<(int, int), (BigInteger Re, BigInteger Im)>();
        foreach (var (t, c) in terms)
        {
            var letters = new char[N];
            int sign = 1;
            for (int l = 0; l < N; l++)
            {
                char ch = t.Letter(l);
                if (ch == 'I') { letters[l] = 'I'; continue; }
                var (nl, sg) = TurnToZ(ch, axis[l]);
                letters[l] = nl; sign *= sg;
            }
            var q = PauliString.Parse(new string(letters));
            int ny = System.Numerics.BitOperations.PopCount(q.X & q.Z);
            for (int col = 0; col < d; col++)
            {
                int row = col ^ (int)q.X;
                int k = (ny + 2 * System.Numerics.BitOperations.PopCount((ulong)col & q.Z)) & 3;
                BigInteger v = sign * (BigInteger)c;
                var add = k switch { 0 => (v, BigInteger.Zero), 1 => (BigInteger.Zero, v), 2 => (-v, BigInteger.Zero), _ => (BigInteger.Zero, -v) };
                var old = h.TryGetValue((row, col), out var o) ? o : (BigInteger.Zero, BigInteger.Zero);
                h[(row, col)] = (old.Item1 + add.Item1, old.Item2 + add.Item2);
            }
        }

        // split into sectors of the undephased bits; the dephased bits are re-indexed 0 .. 2^|S| - 1
        int m = dephased.Count, fs = free.Count, dim = 1 << m;
        static int Sub(int v, List<int> sites)
        {
            int r = 0;
            for (int j = 0; j < sites.Count; j++) r |= ((v >> sites[j]) & 1) << j;
            return r;
        }
        var blocks = Enumerable.Range(0, 1 << fs).Select(_ => new Dictionary<(int, int), (BigInteger Re, BigInteger Im)>()).ToArray();
        foreach (var ((r, c), e) in h)
        {
            if (e.Re.IsZero && e.Im.IsZero) continue;
            int sr = Sub(r, free), sc = Sub(c, free);
            if (sr != sc) throw new InvalidOperationException("H does not conserve the turned letters of the undephased sites");
            blocks[sr][(Sub(r, dephased), Sub(c, dephased))] = e;
        }
        var bars = blocks.Select(b => b.ToDictionary(kv => (kv.Key.Item1 ^ (dim - 1), kv.Key.Item2 ^ (dim - 1)), kv => kv.Value)).ToArray();

        int near = 0, far = 0, cross = 0;
        for (int s = 0; s < blocks.Length; s++)
            for (int t = 0; t < blocks.Length; t++)
            {
                int hn = Hom(blocks[t], blocks[s], dim), hf = Hom(blocks[t], bars[s], dim);
                near += hn;
                far += hf;
                if (s != t) cross += hn + hf;
            }
        return (near, far, blocks.Length, cross);
    }

    /// <summary>dim{g on 0 .. dim − 1 : B_xy·g_y = g_x·A_xy for all x, y}, for A and B Hermitian: the good
    /// components of the graph whose edges are the nonzero off-diagonal entries of A or B.</summary>
    static int Hom(Dictionary<(int, int), (BigInteger Re, BigInteger Im)> a,
                   Dictionary<(int, int), (BigInteger Re, BigInteger Im)> b, int dim)
    {
        static (BigInteger Re, BigInteger Im) Get(Dictionary<(int, int), (BigInteger Re, BigInteger Im)> m, int r, int c) =>
            m.TryGetValue((r, c), out var v) ? v : (BigInteger.Zero, BigInteger.Zero);
        var nbr = Enumerable.Range(0, dim).Select(_ => new List<int>()).ToArray();
        foreach (var (r, c) in a.Keys.Concat(b.Keys).Distinct())
            if (r != c) { nbr[r].Add(c); nbr[c].Add(r); }
        var sec = new (BigInteger Re, BigInteger Im, BigInteger Den)?[dim];
        int good = 0;
        for (int s0 = 0; s0 < dim; s0++)
        {
            if (sec[s0] is not null) continue;
            bool ok = true;
            sec[s0] = (1, 0, 1);
            var stack = new Stack<int>(); stack.Push(s0);
            while (stack.Count > 0)
            {
                int x = stack.Pop();
                if (Get(a, x, x) != Get(b, x, x)) ok = false;
                foreach (int y in nbr[x])
                {
                    var p = Get(a, x, y); var q = Get(b, x, y);
                    BigInteger n2 = q.Re * q.Re + q.Im * q.Im;
                    if (p.Re * p.Re + p.Im * p.Im != n2 || n2.IsZero) ok = false;
                    // g_y = g_x · p / q = g_x · p · conj(q) / |q|^2
                    var (dr, di, dd) = sec[x]!.Value;
                    BigInteger nr = p.Re * q.Re + p.Im * q.Im, ni = p.Im * q.Re - p.Re * q.Im;
                    (BigInteger, BigInteger, BigInteger) v = n2.IsZero ? (0, 0, 1) : (dr * nr - di * ni, dr * ni + di * nr, dd * n2);
                    var gv = BigInteger.GreatestCommonDivisor(BigInteger.GreatestCommonDivisor(v.Item1, v.Item2), v.Item3);
                    if (!gv.IsZero && !gv.IsOne) v = (v.Item1 / gv, v.Item2 / gv, v.Item3 / gv);
                    if (sec[y] is null) { sec[y] = v; stack.Push(y); }
                    else
                    {
                        var (yr, yi, yd) = sec[y]!.Value;
                        if (yr * v.Item3 != v.Item1 * yd || yi * v.Item3 != v.Item2 * yd) ok = false;
                    }
                }
            }
            if (ok) good++;
        }
        return good;
    }

    // ---- the far element a Clifford symmetry is ----

    /// <summary>A Clifford symmetry of (H, jumps), read as what it does to the terms: term i of H goes to
    /// Sign[i] times term Image[i] (the index into the row's terms, in the order they were given), and
    /// every jump to minus itself.</summary>
    public sealed record CliffordMap(int[] Image, int[] Sign)
    {
        /// <summary>Every term fixed: then the Clifford acts on the terms and jumps as a Pauli string
        /// does, and the row has a colouring.</summary>
        public bool FixesEveryTerm => Image.Select((j, i) => j == i).All(x => x) && Sign.All(t => t == 1);
    }

    /// <summary>The node budget of the Clifford search; past it the search reports itself unfinished.</summary>
    public const long MaxCliffordNodes = 2_000_000;

    /// <summary>A Clifford unitary U with U·H·U† = H and U·A·U† = −A for every jump A, found by what it does
    /// to the terms, or null when none exists (Exhausted true) or the budget ran out (Exhausted false).
    ///
    /// Such a U lies in F158's far space: [H, U] = 0 and A·U·A = −U, and it is invertible, so it carries
    /// the palindrome. Conversely every Clifford element of the far space is of this kind: conjugation by
    /// a Clifford sends each string to plus or minus a string, distinct strings to distinct ones, so it
    /// fixes H exactly when it permutes H's terms, each onto a term whose coefficient it matches with its
    /// sign. The search assigns every term an image among the terms of equal magnitude and every jump the
    /// image minus itself, and keeps an assignment only while it is a partial Clifford map: commutation
    /// between every two assigned sources equals commutation between their images, and every product
    /// relation among the sources, phase included, holds among the images (a GF(2) elimination on the
    /// sources carrying the images along, each step one exact string product with its power of i). A
    /// complete assignment that passes both is an isometry of the subspace the sources span (injective
    /// although only the source relations are checked: the images are the sources themselves in another
    /// order, terms permuted bijectively and jumps fixed, so they span the same space, and a surjection
    /// onto a space of the same finite dimension is a bijection), which extends to the whole symplectic
    /// space by Witt's theorem (valid for alternating forms in characteristic 2 and degenerate subspaces), and the signs, being consistent on the
    /// subgroup, extend with it; so it is realised by a Clifford. The colourings are the assignments
    /// fixing every term, and the symmetries of Symmetries() the ones a site permutation and a rotation
    /// induce. The unitary itself is not built here; where the far end is one-dimensional the lifted far
    /// element is it, up to scale.</summary>
    public (CliffordMap? Map, bool Exhausted) CliffordSymmetry()
    {
        int t = terms.Length;
        var image = new int[t];
        var sign = new int[t];
        var used = new bool[t];
        // the assigned sources and their images, each an operator i^phase · string
        var src = new List<(PauliString S, int Ph)>();
        var img = new List<(PauliString S, int Ph)>();
        // the elimination basis: reduced source, its image, pivot = highest bit of the source's key (-1 for a relation)
        var basis = new List<(PauliString S, int Ph, PauliString IS, int IPh, int Pivot)>();
        long nodes = 0;

        static ulong Key(PauliString p) => (p.X << 32) | p.Z;
        static int High(ulong k) => 63 - System.Numerics.BitOperations.LeadingZeroCount(k);
        static (PauliString S, int Ph) Times((PauliString S, int Ph) a, (PauliString S, int Ph) b)
        {
            var (prod, k) = PauliString.Multiply(a.S, b.S);
            return (prod, (a.Ph + b.Ph + k) & 3);
        }

        // add one assignment source -> target if it keeps the map a partial Clifford map
        bool Push((PauliString S, int Ph) a, (PauliString S, int Ph) b)
        {
            for (int k = 0; k < src.Count; k++)
                if (PauliString.Commute(a.S, src[k].S) != PauliString.Commute(b.S, img[k].S)) return false;
            var r = a;
            var ri = b;
            while (Key(r.S) != 0)
            {
                int hb = High(Key(r.S));
                int e = basis.FindIndex(x => x.Pivot == hb);
                if (e < 0) break;
                r = Times(r, (basis[e].S, basis[e].Ph));
                ri = Times(ri, (basis[e].IS, basis[e].IPh));
            }
            if (Key(r.S) == 0 && (!ri.S.IsIdentity || ri.Ph != r.Ph)) return false;   // a relation the images break
            src.Add(a); img.Add(b);
            basis.Add((r.S, r.Ph, ri.S, ri.Ph, Key(r.S) == 0 ? -1 : High(Key(r.S))));
            return true;
        }
        void Pop() { src.RemoveAt(src.Count - 1); img.RemoveAt(img.Count - 1); basis.RemoveAt(basis.Count - 1); }

        foreach (var a in jumpStrings)
            if (!Push((a, 0), (a, 2))) return (null, true);          // −A = i^2 · A

        bool exhausted = true;
        bool Assign(int i)
        {
            if (i == t) return true;
            for (int j = 0; j < t; j++)
            {
                if (used[j] || Math.Abs(terms[j].C) != Math.Abs(terms[i].C)) continue;
                if (++nodes > MaxCliffordNodes) { exhausted = false; return false; }
                int sg = terms[j].C == terms[i].C ? 1 : -1;
                if (!Push((terms[i].S, 0), (terms[j].S, sg == 1 ? 0 : 2))) continue;
                used[j] = true; image[i] = j; sign[i] = sg;
                if (Assign(i + 1)) return true;
                used[j] = false;
                Pop();
                if (!exhausted) return false;
            }
            return false;
        }
        return Assign(0) ? (new CliffordMap(image, sign), true) : (null, exhausted);
    }

    // ---- the anticommuting-sum grammar ----

    /// <summary>The largest lit span the grammar's clique search is run on.</summary>
    public const int MaxGrammarStrings = 256;

    /// <summary>An element of the far space built by the grammar of experiments/THE_PALINDROME_AS_A_COLOURING.md
    /// ("What a sum is"), ported from simulations/anticommuting_sum_census.py, with its kind:
    ///   "sum"          the lifted kernel of the commutator with H on a maximal set of pairwise
    ///                  anticommuting lit strings (such a combination squares to its length times 1);
    ///   "product"      H splits into components on disjoint sites, and the element is the tensor
    ///                  product of one element per component;
    ///   "conditioned"  an undephased site u where every term of H carries I or one letter P: P_u
    ///                  commutes with H, each sector P_u = s (s = ±1) is a row on the other sites with
    ///                  the letter replaced by s, and the element is the sum over s of (1 + s·P_u)/2
    ///                  times an element of that sector, recursively.
    /// The kinds are tried in that order (product, sum, conditioned), as the census tries them. The
    /// element is checked, not trusted: every string lit, [H, G] = 0 exactly, and G invertible, its
    /// matrix of full rank modulo a prime p = 1 mod 4 with i sent to a square root of −1, which forces
    /// det G ≠ 0 (Invertible). The check guards the construction and is not expected to fire: sums,
    /// products and sector sums of invertible elements are invertible by construction. Null when the
    /// grammar finds nothing; on a broken row it can find nothing, since an invertible element of the
    /// far space is a palindrome. The sum stage is skipped on a lit span past MaxGrammarStrings, and
    /// Invertible needs N ≤ 10; Explanation reports those rows unfinished.</summary>
    public (string Kind, IReadOnlyList<(string Letters, long Coefficient)> Element)? GrammarElement()
    {
        if (grammar.HasValue) return grammar.Value.Found;
        var found = Grammar();
        if (found is { } f && !(CheckElement(f.Element) is { AllLit: true, CommutesWithH: true } && Invertible(f.Element)))
            throw new InvalidOperationException($"the grammar built a {f.Kind} that fails its check");
        grammar = (found, true);
        return found;
    }
    ((string Kind, IReadOnlyList<(string Letters, long Coefficient)> Element)? Found, bool Done)? grammar;

    (string Kind, IReadOnlyList<(string Letters, long Coefficient)> Element)? Grammar()
    {
        // components of H and the jumps on the sites (a site nothing touches is its own component; a jump
        // on several sites joins them, so that no jump is split between two components)
        var parent = Enumerable.Range(0, N).ToArray();
        int Find(int x) { while (parent[x] != x) x = parent[x]; return x; }
        foreach (var t in terms.Select(x => x.S).Concat(jumpStrings))
        {
            var sites = Enumerable.Range(0, N).Where(l => Touches(t, l)).ToList();
            foreach (int q in sites.Skip(1)) parent[Find(q)] = Find(sites[0]);
        }
        var comps = Enumerable.Range(0, N).GroupBy(Find).Select(g => g.ToArray()).ToList();
        if (comps.Count > 1)
        {
            var product = new Dictionary<PauliString, BigInteger> { [PauliString.Identity] = 1 };
            foreach (var c in comps)
            {
                // a component no jump touches contributes the identity (every string is lit there)
                var sub = SubRow(c);
                IReadOnlyList<(string Letters, long Coefficient)> partElement;
                if (sub is null) partElement = new[] { (new string('I', c.Length), 1L) };
                else if (sub.GrammarElement() is { } part) partElement = part.Element;
                else return null;
                var next = new Dictionary<PauliString, BigInteger>();
                foreach (var (p0, c0) in product)
                    foreach (var (letters, c1) in partElement)
                    {
                        var q = Embed(PauliString.Parse(letters), c);
                        var key = new PauliString(p0.X | q.X, p0.Z | q.Z);
                        next[key] = (next.TryGetValue(key, out var o) ? o : 0) + c0 * c1;
                    }
                product = next;
            }
            return Integral("product", product);
        }

        // a maximal set of pairwise anticommuting lit strings whose kernel is not zero
        var litStrings = LitStrings();
        if (litStrings.Count <= MaxGrammarStrings && CliqueSum(litStrings) is { } sum) return ("sum", sum);

        // a sector split over a letter an undephased site keeps
        for (int u = 0; u < N; u++)
        {
            if (jumpStrings.Any(a => Touches(a, u))) continue;
            foreach (char letter in "XYZ")
            {
                if (!terms.All(t => !Touches(t.S, u) || t.S.Letter(u) == letter)) continue;
                var parts = new List<IReadOnlyList<(string Letters, long Coefficient)>>();
                foreach (int s in new[] { 1, -1 })
                {
                    var sector = Sector(u, letter, s);
                    if (sector?.GrammarElement() is not { } part) break;
                    parts.Add(part.Element);
                }
                if (parts.Count < 2) continue;
                // 2·G = sum_s (1 + s·P_u) ⊗ G_s
                var el = new Dictionary<PauliString, BigInteger>();
                var pu = PauliString.Parse(new string(Enumerable.Range(0, N).Select(l => l == u ? letter : 'I').ToArray()));
                for (int k = 0; k < 2; k++)
                {
                    int s = k == 0 ? 1 : -1;
                    foreach (var (letters, c) in parts[k])
                    {
                        var g = Insert(PauliString.Parse(letters), u);
                        el[g] = (el.TryGetValue(g, out var o1) ? o1 : 0) + c;
                        var gp = new PauliString(g.X | pu.X, g.Z | pu.Z);
                        el[gp] = (el.TryGetValue(gp, out var o2) ? o2 : 0) + s * c;
                    }
                }
                if (Integral("conditioned", el) is { } conditioned) return conditioned;
            }
        }
        return null;
    }

    static bool Touches(PauliString p, int site) => (((p.X | p.Z) >> site) & 1) != 0;

    // the kernel on the first maximal anticommuting clique (Bron–Kerbosch with pivot, in the span's
    // order) that has one, lifted and exactly checked
    IReadOnlyList<(string Letters, long Coefficient)>? CliqueSum(IReadOnlyList<PauliString> span)
    {
        int n = span.Count;
        var adj = new ulong[n][];
        int words = (n + 63) / 64;
        for (int i = 0; i < n; i++) adj[i] = new ulong[words];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (i != j && !PauliString.Commute(span[i], span[j])) adj[i][j >> 6] |= 1UL << (j & 63);
        List<(string, long)>? result = null;
        var r = new List<int>();
        void Bk(ulong[] pSet, ulong[] xSet)
        {
            if (result is not null) return;
            bool pEmpty = pSet.All(w => w == 0), xEmpty = xSet.All(w => w == 0);
            if (pEmpty && xEmpty)
            {
                var clique = r.Select(k => span[k]).ToList();
                var kernel = LiftKernelOf(clique, LiftPrimes.Length);
                if (kernel.Count > 0) result = kernel[0];
                return;
            }
            // pivot: the vertex of P ∪ X with the most neighbours in P
            int pivot = -1, best = -1;
            for (int v = 0; v < n; v++)
            {
                if (((pSet[v >> 6] | xSet[v >> 6]) >> (v & 63) & 1) == 0) continue;
                int deg = 0;
                for (int w = 0; w < words; w++) deg += System.Numerics.BitOperations.PopCount(adj[v][w] & pSet[w]);
                if (deg > best) { best = deg; pivot = v; }
            }
            var pp = (ulong[])pSet.Clone();
            var xx = (ulong[])xSet.Clone();
            for (int v = 0; v < n; v++)
            {
                if ((pp[v >> 6] >> (v & 63) & 1) == 0) continue;
                if ((adj[pivot][v >> 6] >> (v & 63) & 1) != 0) continue;
                r.Add(v);
                Bk(pp.Select((w, k) => w & adj[v][k]).ToArray(), xx.Select((w, k) => w & adj[v][k]).ToArray());
                r.RemoveAt(r.Count - 1);
                if (result is not null) return;
                pp[v >> 6] &= ~(1UL << (v & 63));
                xx[v >> 6] |= 1UL << (v & 63);
            }
        }
        var all = new ulong[words];
        for (int v = 0; v < n; v++) all[v >> 6] |= 1UL << (v & 63);
        Bk(all, new ulong[words]);
        return result;
    }

    // the row restricted to the sites of one component (sites renumbered in order), or null with no jump there
    EndCount? SubRow(int[] sites)
    {
        var index = new Dictionary<int, int>();
        for (int k = 0; k < sites.Length; k++) index[sites[k]] = k;
        string Restrict(PauliString p) => new(sites.Select(l => p.Letter(l)).ToArray());
        bool Inside(PauliString p) => Enumerable.Range(0, N).All(l => !Touches(p, l) || index.ContainsKey(l));
        var h = terms.Where(t => Inside(t.S)).Select(t => (Restrict(t.S), t.C)).ToList();
        var jumps = jumpStrings.Where(Inside).Select(Restrict).ToList();
        return jumps.Count == 0 ? null : new EndCount(world, sites.Length, h, jumps);
    }

    // the sector P_u = s: the site u removed, P_u replaced by s in every term
    EndCount? Sector(int u, char letter, int s)
    {
        if (N == 1) return null;
        var sites = Enumerable.Range(0, N).Where(l => l != u).ToArray();
        string Restrict(PauliString p) => new(sites.Select(l => p.Letter(l)).ToArray());
        var h = terms.Select(t => (Restrict(t.S), Touches(t.S, u) ? s * t.C : t.C))
                     .Where(t => t.Item1.Any(c => c != 'I')).ToList();
        var jumps = jumpStrings.Select(Restrict).ToList();
        return new EndCount(world, N - 1, h, jumps);
    }

    // a string on the component's sites placed back on the full row
    PauliString Embed(PauliString q, int[] sites)
    {
        ulong x = 0, z = 0;
        for (int k = 0; k < sites.Length; k++)
        {
            x |= ((q.X >> k) & 1) << sites[k];
            z |= ((q.Z >> k) & 1) << sites[k];
        }
        return new PauliString(x, z);
    }

    // a string on N - 1 sites with an identity inserted at site u
    static PauliString Insert(PauliString q, int u)
    {
        ulong low = (1UL << u) - 1;
        return new PauliString((q.X & low) | ((q.X & ~low) << 1), (q.Z & low) | ((q.Z & ~low) << 1));
    }

    // coprime integers, zeros dropped; null past MaxCoefficient or when nothing is left
    (string Kind, IReadOnlyList<(string Letters, long Coefficient)> Element)? Integral(string kind, Dictionary<PauliString, BigInteger> el)
    {
        var nz = el.Where(kv => !kv.Value.IsZero).ToList();
        if (nz.Count == 0) return null;
        BigInteger g = nz.Aggregate(BigInteger.Zero, (a, kv) => BigInteger.GreatestCommonDivisor(a, kv.Value));
        if (nz.Any(kv => BigInteger.Abs(kv.Value / g) > MaxCoefficient)) return null;
        return (kind, nz.Select(kv => (kv.Key.ToString(N), (long)(kv.Value / g))).ToList());
    }

    /// <summary>Is the combination sum c_P P invertible? Its 2^N × 2^N matrix over Z[i] is reduced modulo
    /// a prime p = 1 mod 4 with i sent to a square root of −1; full rank there means det ≠ 0 mod p, hence
    /// det ≠ 0. Exact, one-sided (a rank below full at one prime decides nothing, and the second prime is
    /// asked). N ≤ 10.</summary>
    public bool Invertible(IReadOnlyList<(string Letters, long Coefficient)> element)
    {
        if (N > 10) throw new InvalidOperationException("the matrix is dense; N <= 10");
        int d = 1 << N;
        foreach (long p in ModP.Primes)
        {
            long iu = ModP.SqrtMinusOne(p);
            var m = new long[d][];
            for (int r = 0; r < d; r++) m[r] = new long[d];
            foreach (var (letters, c) in element)
            {
                var q = PauliString.Parse(letters);
                // P|b> = i^(#Y) (-1)^(popcount(b & Z)) |b xor X>, Y = X·Z·i on each site
                int ny = System.Numerics.BitOperations.PopCount(q.X & q.Z);
                for (int col = 0; col < d; col++)
                {
                    int row = col ^ (int)q.X;
                    int k = ny + 2 * System.Numerics.BitOperations.PopCount((ulong)col & q.Z);
                    long ph = (k & 3) switch { 0 => 1, 1 => iu, 2 => p - 1, _ => p - iu };
                    m[row][col] = ModP.AddMod(m[row][col], ModP.MulMod(ModP.Mod(c, p), ph, p), p);
                }
            }
            if (ModP.Rank(m, p) == d) return true;
        }
        return false;
    }

    /// <summary>For a palindromic row, what explains it: "colouring" (a single lit string commuting with
    /// H), "symmetry" (the far element of a site symmetry, SymmetryElement), "clifford" (a Clifford
    /// element of the far space that neither of those is, CliffordSymmetry), the grammar's kind
    /// ("sum", "product" or "conditioned", GrammarElement) for an element outside the Clifford group,
    /// "unfinished" when the Clifford search ran out of its budget before either answer
    /// (so nothing past the site symmetry is claimed), or null when the palindrome holds by the counts
    /// with none of these behind it. The
    /// kinds nest, each checked after the ones before it: a colouring is a Clifford element fixing
    /// every term, a site symmetry's element a Clifford one a site permutation and a rotation induce.
    /// N ≤ 8, the bound of Symmetries().</summary>
    public string? Explanation()
    {
        if (Colourings().Count > 0) return "colouring";
        if (SymmetryElement() is not null) return "symmetry";
        var (map, exhausted) = CliffordSymmetry();
        if (map is not null) return "clifford";
        if (!exhausted) return "unfinished";
        if (N > 10) return "unfinished";
        if (GrammarElement() is { } g) return g.Kind;
        return LitStrings().Count > MaxGrammarStrings ? "unfinished" : null;
    }

    static bool PairwiseAnticommuting(IReadOnlyList<(string Letters, long Coefficient)> element)
    {
        var p = element.Select(e => PauliString.Parse(e.Letters)).ToList();
        for (int i = 0; i < p.Count; i++)
            for (int j = 0; j < i; j++)
                if (PauliString.Commute(p[i], p[j])) return false;
        return true;
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

    /// <summary>The primes a lift may combine: ModP's two, then the primes just below the first, found by
    /// ModP.IsPrime. Each added prime widens the rational-reconstruction bound by about 2^15; a vector
    /// whose rationals exceed what all ten reach, or whose integers pass MaxCoefficient, stays unlifted
    /// and its row goes to the ranks.</summary>
    public static readonly long[] LiftPrimes = ModP.Primes.Concat(
        Enumerable.Range(0, 8).Aggregate((List: new List<long>(), Next: ModP.Primes[0] - 2), (acc, _) =>
        {
            long q = acc.Next;
            while (!ModP.IsPrime(q)) q -= 2;
            acc.List.Add(q);
            return (acc.List, q - 2);
        }).List).ToArray();

    /// <summary>A basis of the commutator's kernel on the dark (far = false) or lit (far = true) span,
    /// lifted from GF(p): the reduced-echelon basis at the first prime, one vector per free column; the
    /// same basis at further primes of LiftPrimes, combined by the Chinese remainder theorem, entry by
    /// entry lifted to the rationals by rational reconstruction and scaled to coprime integers, and KEPT
    /// only if it commutes with H exactly. A prime of larger rank, or of equal rank and earlier pivot
    /// columns, becomes the anchor; a prime whose pivot columns come later is skipped (either way one of
    /// the two reductions was bad, and the later profile is the bad one). Primes are added only until the vector passes, since the bound the
    /// reconstruction needs grows with the vector's coefficients, not with the span. The kept vectors lie
    /// in the span and are independent (each carries its own free column), so their number is an exact
    /// LOWER bound on that end's count, whatever the primes did: a lift is a guess and the exact check
    /// decides it. Empty past MaxLiftColumns. primes limits how many of LiftPrimes may be combined (all
    /// by default, the reading the verdict takes; a smaller budget is read fresh, not cached).</summary>
    public IReadOnlyList<IReadOnlyList<(string Letters, long Coefficient)>> LiftedKernel(bool far, int primes = 0)
    {
        if (primes < 0 || primes > LiftPrimes.Length) throw new ArgumentOutOfRangeException(nameof(primes));
        int budget = primes == 0 ? LiftPrimes.Length : primes;
        bool full = budget == LiftPrimes.Length;
        if (full && lifted.TryGetValue(far, out var cached)) return cached;
        var span = far ? LitStrings() : DarkStrings();
        var kept = new List<List<(string, long)>>();
        if (full) lifted[far] = kept;
        kept.AddRange(LiftKernelOf(span, budget));
        return kept;
    }

    // the lifted kernel of the commutator with H on any span of strings, at a budget of LiftPrimes
    List<List<(string, long)>> LiftKernelOf(IReadOnlyList<PauliString> span, int budget)
    {
        var kept = new List<List<(string, long)>>();
        int s = span.Count;
        if (s == 0 || s > MaxLiftColumns) return kept;
        var cols = Columns(span);
        int rows = cols.Max(c => c.Count == 0 ? 0 : c.Keys.Max() + 1);

        // the reduced-echelon kernel basis at one prime: pivot columns and, per free column, its vector
        (List<int> Pivots, Dictionary<int, long[]> Basis) Echelon(long p)
        {
            var m = new long[rows][];
            for (int r = 0; r < rows; r++) m[r] = new long[s];
            for (int j = 0; j < s; j++)
                foreach (var (r, x) in cols[j]) m[r][j] = ModP.Mod(x, p);
            var pivots = new List<int>();
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
                pivots.Add(c);
                rank++;
            }
            var basis = new Dictionary<int, long[]>();
            foreach (int fc in Enumerable.Range(0, s).Except(pivots))
            {
                var v = new long[s];
                v[fc] = 1;
                for (int i = 0; i < rank; i++) v[pivots[i]] = ModP.Mod(-m[i][fc], p);
                basis[fc] = v;
            }
            return (pivots, basis);
        }

        // the anchor: the prime whose profile the others must share. A rank mod p never exceeds the
        // rational rank on any prefix of the columns, so a bad prime's pivots sit at or after the
        // rational ones, pivot by pivot, and the rational profile is the lexicographically smallest of
        // maximal rank: a later prime of larger rank, or of equal rank and a smaller profile, shows the
        // anchor's reduction was bad (it can keep the full rank and only move a pivot). The lift
        // then restarts on that prime, dropping what it kept (exact vectors, but on the bad profile's
        // free columns, whose independence from the new ones is not given)
        List<int> profile = null!;
        Dictionary<int, BigInteger[]> residues = null!;
        BigInteger modulus = 0;
        SortedSet<int> open = null!;
        void Anchor(long p, (List<int> Pivots, Dictionary<int, long[]> Basis) e)
        {
            kept.Clear();
            profile = e.Pivots;
            residues = e.Basis.ToDictionary(kv => kv.Key, kv => kv.Value.Select(x => (BigInteger)x).ToArray());
            modulus = p;
            open = new SortedSet<int>(e.Basis.Keys);
        }
        Anchor(LiftPrimes[0], Echelon(LiftPrimes[0]));
        int next = 1;
        while (true)
        {
            foreach (int fc in open.ToList())
                if (TryLift(residues[fc], modulus) is { } vec) { kept.Add(vec); open.Remove(fc); }
            if (open.Count == 0 || next >= budget) break;
            long p = LiftPrimes[next++];
            var e = Echelon(p);
            if (e.Pivots.Count > profile.Count || (e.Pivots.Count == profile.Count && Earlier(e.Pivots, profile)))
            {
                Anchor(p, e);
                continue;
            }
            if (!e.Pivots.SequenceEqual(profile)) continue;
            foreach (int fc in open)
                for (int j = 0; j < s; j++)
                    residues[fc][j] = ModP.Crt(residues[fc][j], modulus, e.Basis[fc][j], p);
            modulus *= p;
        }
        return kept;

        static bool Earlier(List<int> a, List<int> b)
        {
            for (int i = 0; i < a.Count; i++)
                if (a[i] != b[i]) return a[i] < b[i];
            return false;
        }

        List<(string, long)>? TryLift(BigInteger[] v, BigInteger mod)
        {
            var fracs = new (BigInteger Num, BigInteger Den)[s];
            for (int j = 0; j < s; j++)
            {
                if (ModP.RationalReconstruct(v[j], mod) is not { } q) return null;
                fracs[j] = q;
            }
            BigInteger lcm = fracs.Aggregate(BigInteger.One, (a, f) => a / BigInteger.GreatestCommonDivisor(a, f.Den) * f.Den);
            var ints = fracs.Select(f => f.Num * (lcm / f.Den)).ToArray();
            BigInteger g = ints.Aggregate(BigInteger.Zero, (a, x) => BigInteger.GreatestCommonDivisor(a, x));
            if (ints.Any(x => BigInteger.Abs(x / g) > MaxCoefficient)) return null;
            var vec = new List<(string, long)>();
            for (int j = 0; j < s; j++)
                if (!ints[j].IsZero) vec.Add((span[j].ToString(N), (long)(ints[j] / g)));
            return CommutesWithH(vec) ? vec : null;
        }
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
