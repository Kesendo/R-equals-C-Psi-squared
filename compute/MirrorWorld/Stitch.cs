namespace MirrorWorld;

using System.Numerics;

// The stitch (built 2026-11-10 to empty the sieve's first residue; the name is a working one: a phase stitched across
// every bond). Under Z dephasing on every site the far end of F158 is X^⊗N times the diagonal matrices, and Theorem 2
// of docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md reads it as the flat sections d of the complement
// connection, (★) H_x̄ȳ·d_y = d_x·H_xy on every edge. A colouring is a flat section that is a character, the phase of a
// LINEAR form (−1)^(s·x) ("colourings are characters"). The stitch reads the next order: a section that is a product
// of one-site phases and two-site phases on the pairs of sites a cell's windows reach,
//   d_x = ∏_l a(x_l) · ∏_(l, l+o) b_o(x_l, x_(l+o)),  o = 1 .. span − 1,
// the same factors on every site and every pair (the phase of a QUADRATIC form when the factors are ±1 and ±i). Every
// (★) is then binomial: d_y/d_x is a Laurent monomial in the factors (only those touching the flipped bits change) and
// H_xy/H_x̄ȳ a ratio of Gaussian integers, so the theodolite's Binomials decides it, with no number taken. (★) at (x, y)
// and at (y, x) together, H being Hermitian, force |H_xy| = |H_x̄ȳ| on every edge (Theorem 2's condition (i)), and then
// the two are one equation. Where it reads, W = X^⊗N·D is an invertible element of the far end (every factor nonzero),
// so W·(near end) lies in the far end and W⁻¹·(far end) in the near end: the two ends have one dimension and the
// spectrum pairs (F158).
//
// The no it gives on the way. An edge whose two entries differ in modulus (one of them zero among them) breaks
// Theorem 2's condition (i), a diagonal entry with H_xx ≠ H_x̄x̄ its condition (ii); the component holding it is then not
// good, the far end is smaller than the near end, and the spectrum does not pair. That is a no about the spectrum, not
// about the ansatz, exact, and it carries to every N past the window by the same locality as the yes (Breaks). Its
// condition (iii), the holonomy, is read on short cycles (LocalHolonomy): a square of two flip sets or a triangle whose
// ratios do not multiply to 1 breaks the component too, and a cycle in a local region recurs at every larger N, so that
// no holds for every N from the smallest chain or ring that holds it. A cycle round the whole ring is not local and is
// not read here. Where no condition fails on what is read and the binomial system has no solution, the stitch is
// silent: the section may exist and not be of this form.
//
// Every N. An equation at the edge (x, x ⊕ F), F the flip set of a window, changes only the factors that touch F and
// reads only the terms with flip set F, whose Z letters lie within span − 1 of F: it depends on the bits of at most
// 3·span − 2 sites around F and on nothing else. On a ring of N ≥ 3·span sites that region never meets itself, every
// configuration of it occurs, and no other: the set of distinct equations is the same at every such N, so one solve at
// N = 3·span reads every larger ring. On the chain a region within reach of an end carries the end, and from
// N = 2·(3·span − 2) on every end type (each offset from each end) and the bulk occur and no region meets both ends: one
// solve reads every longer chain. Both are
// checked where they are used: the equation set at the window and one site further must agree, or the reading throws.
// The diagonal check reads a whole string's diagonal, not a window, and is N-independent all the same: H_xx = H_x̄x̄ for
// every x says that the odd-weight diagonal part of H vanishes (the characters (−1)^(z·x) being independent), and past
// the window that is a statement about the cell's own coefficients.
// The N below the window are read one by one, exactly, and say nothing past themselves.
//
// The first reading, the sieve's first ring residue: XZ ± ZY (± ZZ) under Z dephasing, d_x = i^(±|x|)·(−1)^(Σ x_l x_(l+1))
// on the ring, W = X^⊗N·S^⊗N·∏ CZ (S† for −ZY), a Clifford that sends X_l Z_(l+1) to Z_(l−1) Y_l and Z_l Y_(l+1) to
// X_(l+1) Z_(l+2): the two halves of the bond swapped and moved on by one site, a permutation of the ring's terms that the
// chain's ends cut. Adjacent: EndCount.CliffordSymmetry searches, at one N, the Cliffords of the far end that permute
// H's terms; this W is one, read here for every N.
//
// Deliberately outside: jumps other than Z on every site (Theorem 2 turns any one Pauli letter per site to Z; the stitch
// reads Z only), factors that change from site to site (a turned or periodic section: the theodolite reads those for
// products), pairs beyond the cell's span, and a silence's meaning.
//
// Parent: the theodolite, whose cell and every-N reading it extends from products to pairs; the end count's Theorem 2
// is adopted, not inherited.
public sealed class Stitch : GameObject
{
    public Stitch(Theodolite theodolite) : base(theodolite) { }

    // left: what the stitch itself produces.
    public override IReadOnlyList<string> Own => new[] { "stitch", "window" };

    // One N: the stitch reads a section (the spectrum pairs), finds condition (i) or (ii) of Theorem 2 broken (it does not
    // pair), or is silent.
    public enum Answer { Reads, Breaks, Silent }

    // the every-N reading of one geometry: Bulk the system at the window (null when it does not read there; it holds at
    // every N from the window on), Small the N below the window, from the smallest, where a fixed-N solve reads;
    // BreaksFrom the N from which the stitch's no holds at every N (a condition of Theorem 2 failing at the window, or a
    // bad local cycle, LocalHolonomy), null when it gives none; SmallBreaks the N below it where a fixed-N equation set
    // breaks condition (i) or (ii)
    public sealed record Reach(bool Ring, int Window, Binomials.Reduced? Bulk, IReadOnlyList<int> Small, int? BreaksFrom,
        IReadOnlyList<int> SmallBreaks)
    {
        public bool Reads(int n) => n >= Window ? Bulk is not null : Small.Contains(n);
        public bool Breaks(int n) => BreaksFrom is int b && n >= b || SmallBreaks.Contains(n);
        public bool BulkBreaks => BreaksFrom is not null;
    }

    /// <summary>The most sites the every-N reading enumerates (2^MaxSites basis states at the window and one past it):
    /// the chain reaches it at span 2, the ring at span 5.</summary>
    public const int MaxSites = 16;

    /// <summary>The most bits of a local region the holonomy enumerates (2^MaxRegionBits configurations per pair of flip
    /// sets): the region is at most 5·span − 2 sites, so up to span 4.</summary>
    public const int MaxRegionBits = 20;

    // ring: 3·span (a region of 3·span − 2 sites never meets itself); chain: 2·(3·span − 2) (every end type and the bulk)
    public static int WindowOf(Theodolite.Cell cell, bool ring)
    {
        int span = cell.Normalized().Span;
        return ring ? 3 * span : Math.Max(4, 2 * (3 * span - 2));
    }

    // The cell on every N at once; null when the cell's jumps are not Z on every site. When the window and one site more
    // pass MaxSites only the local no is read (no yes, no fixed-N answers).
    public static Reach? EveryN(Theodolite.Cell cell, bool ring)
    {
        var norm = cell.Normalized();
        if (!ZOnly(norm)) return null;
        int window = WindowOf(norm, ring);
        int from = ring ? norm.SmallestRing : Math.Max(2, norm.Span);
        // the fast no first: a bad local cycle that holds from the window on decides every N from there, so the window's
        // solves are not needed and only the N below the no's threshold are read one by one
        var holonomy = LocalHolonomy(norm);
        int? local = ring ? holonomy.Ring : holonomy.Chain;
        if (window + 1 > MaxSites)
            return new Reach(ring, window, null, Array.Empty<int>(), local, Array.Empty<int>());
        if (local is int early && early <= window)
        {
            var below = Enumerable.Range(from, Math.Max(0, early - from)).Select(n => (n, Decide(norm, n, ring))).ToList();
            var reads = below.Where(a => a.Item2 == Answer.Reads).Select(a => a.n).ToList();
            return new Reach(ring, window, null, reads, early, below.Where(a => a.Item2 == Answer.Breaks).Select(a => a.n).ToList());
        }
        var at = Equations(norm, window, ring);
        var next = Equations(norm, window + 1, ring);
        if (at.Contradiction != next.Contradiction || !at.Contradiction && !at.Keys.SetEquals(next.Keys))
            throw new InvalidOperationException($"the equations at N = {window} and {window + 1} differ: the window is not reached");
        var bulk = at.Contradiction ? null : Binomials.Solve(at.Equations, Unknowns(norm.Span));
        var answers = Enumerable.Range(from, Math.Max(0, window - from)).Select(n => (n, Decide(norm, n, ring))).ToList();
        int? breaksFrom = at.Contradiction ? window : null;
        if (local is int h) breaksFrom = Math.Min(breaksFrom ?? int.MaxValue, h);
        var small = answers.Where(a => a.Item2 == Answer.Reads).Select(a => a.n).ToList();
        // a yes and a no at one N would be a contradiction in the construction, not a finding about the cell
        if (bulk is not null && breaksFrom is not null || small.Any(n => breaksFrom is int b && n >= b))
            throw new InvalidOperationException("the stitch reads a section where its own no holds");
        return new Reach(ring, window, bulk, small, breaksFrom,
            answers.Where(a => a.Item2 == Answer.Breaks).Select(a => a.n).ToList());
    }

    // The cell on one ring or chain of n sites, the three answers; Silent also when the jumps are not Z on every site.
    public static Answer Decide(Theodolite.Cell cell, int n, bool ring)
    {
        var norm = cell.Normalized();
        if (!ZOnly(norm)) return Answer.Silent;
        var e = Equations(norm, n, ring);
        if (e.Contradiction) return Answer.Breaks;
        return Binomials.Solve(e.Equations, Unknowns(norm.Span)) is null ? Answer.Silent : Answer.Reads;
    }

    // The cell on one ring or chain of n sites: the reduced system of the factors, or null when silent.
    public static Binomials.Reduced? At(Theodolite.Cell cell, int n, bool ring)
    {
        var norm = cell.Normalized();
        if (!ZOnly(norm)) return null;
        var e = Equations(norm, n, ring);
        return e.Contradiction ? null : Binomials.Solve(e.Equations, Unknowns(norm.Span));
    }

    static bool ZOnly(Theodolite.Cell cell)
    {
        var z = LetterCube.LetterMatrix('Z');
        return cell.Jumps.Count > 0 && cell.Jumps.All(j => j.Rate > 0 && Same(j.Jump, z));
    }

    static bool Same(LetterCube.Gauss[,] a, LetterCube.Gauss[,] b)
    {
        for (int i = 0; i < 2; i++)
            for (int k = 0; k < 2; k++)
                if (a[i, k] != b[i, k]) return false;
        return true;
    }

    static List<string> Unknowns(int span)
    {
        var u = new List<string> { "a0", "a1" };
        for (int o = 1; o < span; o++)
            foreach (var p in new[] { "00", "01", "10", "11" }) u.Add($"b{o}_{p}");
        return u;
    }

    sealed record EquationSet(bool Contradiction, List<Binomial> Equations, HashSet<string> Keys);

    // Every (★) of the cell on n sites, deduplicated by their content (exponents and constants), which carries no site
    // index, so two N with the same local equations give the same set.
    static EquationSet Equations(Theodolite.Cell cell, int n, bool ring)
    {
        if (n > 24) throw new ArgumentOutOfRangeException(nameof(n), n, "the stitch enumerates 2^n basis states; at most 24 sites");
        int span = cell.Span;
        if (ring && n < cell.SmallestRing) throw new ArgumentOutOfRangeException(nameof(n), n, $"a ring of this cell has at least {cell.SmallestRing} sites");
        var byMask = ByMask(cell, n, ring);
        var pairs = new List<(int A, int B, int O)>();
        for (int o = 1; o < span; o++)
            for (int l = 0; l < n; l++)
                if (ring || l + o < n) pairs.Add((l, (l + o) % n, o));
        ulong all = (1UL << n) - 1;
        GaussBig Element(ulong col, ulong mask) => ElementOf(byMask, col, mask);
        var names = Unknowns(span);
        int Pair(int o, ulong x, int a, int b) => 2 + 4 * (o - 1) + 2 * (int)((x >> a) & 1) + (int)((x >> b) & 1);
        // the pairs touching each site, so an edge visits only the factors its flip set changes
        var touching = Enumerable.Range(0, n).Select(l => pairs.Where(p => p.A == l || p.B == l).ToArray()).ToArray();

        var eqs = new Dictionary<string, Binomial>();
        var diff = new int[names.Count];
        foreach (var mask in byMask.Keys)
        {
            var flipped = Enumerable.Range(0, n).Where(l => ((mask >> l) & 1) == 1).ToArray();
            var changed = flipped.SelectMany(l => touching[l]).Distinct().ToArray();
            for (ulong x = 0; x <= all; x++)
            {
                ulong y = x ^ mask;
                if (mask == 0)
                {
                    if (Element(x, 0) != Element(all ^ x, 0)) return new EquationSet(true, new(), new());
                    continue;
                }
                if (y < x) continue;                                     // each edge once: with equal moduli (★) at (y, x) is the same equation
                var hxy = Element(y, mask);                              // <x|H|y>
                var hbar = Element(all ^ y, mask);                       // <x̄|H|ȳ>
                if (hxy.IsZero && hbar.IsZero) continue;
                if (hxy.Norm != hbar.Norm) return new EquationSet(true, new(), new());   // condition (i), zero against nonzero included
                Array.Clear(diff);                                       // exponents of d_y / d_x
                foreach (int l in flipped) { diff[(int)((y >> l) & 1)]++; diff[(int)((x >> l) & 1)]--; }
                foreach (var (pa, pb, o) in changed) { diff[Pair(o, y, pa, pb)]++; diff[Pair(o, x, pa, pb)]--; }
                var key = new System.Text.StringBuilder();
                for (int i = 0; i < diff.Length; i++) if (diff[i] != 0) key.Append(names[i]).Append('^').Append(diff[i]).Append(',');
                key.Append($" = ({hxy})/({hbar})");
                string k = key.ToString();
                if (eqs.ContainsKey(k)) continue;
                var clean = new Dictionary<string, int>();
                for (int i = 0; i < diff.Length; i++) if (diff[i] != 0) clean[names[i]] = diff[i];
                eqs[k] = new Binomial(clean, hxy, hbar, $"edge {x}-{y} on {n} sites");   // d_y/d_x = H_xy/H_x̄ȳ
            }
        }
        return new EquationSet(false, eqs.Values.ToList(), eqs.Keys.ToHashSet());
    }

    // the cell's terms on n sites (ring or chain), grouped by flip set
    static Dictionary<ulong, List<(ulong X, ulong Z, long C)>> ByMask(Theodolite.Cell cell, int n, bool ring)
    {
        var terms = new List<(ulong X, ulong Z, long C)>();
        void Add(int[] sites, string letters, long c)
        {
            ulong x = 0, z = 0;
            for (int i = 0; i < sites.Length; i++)
            {
                ulong bit = 1UL << sites[i];
                if (letters[i] is 'X' or 'Y') x |= bit;
                if (letters[i] is 'Z' or 'Y') z |= bit;
            }
            terms.Add((x, z, c));
        }
        for (int l = 0; l < n; l++)
            foreach (var (letter, c) in cell.Fields) Add(new[] { l }, letter.ToString(), c);
        foreach (var (letters, off, c) in cell.Windows)
            for (int s = 0; s < n; s++)
            {
                if (!ring && s + off.Max() >= n) continue;
                Add(off.Select(o => (s + o) % n).ToArray(), letters, c);
            }
        return terms.GroupBy(t => t.X).ToDictionary(g => g.Key, g => g.ToList());
    }

    // <row|P_mask|col> summed over the terms with that flip set: P|col> = i^|x∧z|·(−1)^|col∧z|·|col ⊕ x> (Y = iXZ)
    static GaussBig ElementOf(Dictionary<ulong, List<(ulong X, ulong Z, long C)>> byMask, ulong col, ulong mask)
    {
        BigInteger re = 0, im = 0;
        if (!byMask.TryGetValue(mask, out var list)) return new GaussBig(0, 0);
        foreach (var t in list)
            switch ((BitOperations.PopCount(t.X & t.Z) + 2 * BitOperations.PopCount(col & t.Z)) % 4)
            {
                case 0: re += t.C; break;
                case 1: im += t.C; break;
                case 2: re -= t.C; break;
                default: im -= t.C; break;
            }
        return new GaussBig(re, im);
    }

    // The local holonomy (Theorem 2's condition (iii) on short cycles). On a chain of 10·span sites, every square
    // u -> u⊕F -> u⊕F⊕G -> u⊕G -> u of two flip sets F, G of the cell's terms within reach of each other, and every
    // triangle where F ⊕ G is a flip set too, is walked with the ratios r(v, F) = H_(v, v⊕F) / H_(v̄, v̄⊕F); the bits
    // outside the region R (F ∪ G widened by span − 1 on each side) are set to 0, since no ratio on the cycle reads them.
    // A cycle whose ratios do not multiply to 1, or an edge on it whose two entries differ in modulus, makes its component
    // not good: the spectrum does not pair. The placement recurs on every longer chain (an end placement at the same end,
    // a bulk one anywhere), so the no holds at every N from the smallest chain that holds it; on a ring only the bulk
    // placements count, from the smallest ring that holds R (no window wraps onto it). Returns those
    // two thresholds; null where no bad local cycle is found, which decides nothing (a cycle round the whole ring is not
    // local, and the cycles walked are squares and triangles only).
    public static (int? Chain, int? Ring) LocalHolonomy(Theodolite.Cell cell)
    {
        var norm = cell.Normalized();
        if (!ZOnly(norm)) return (null, null);
        int span = norm.Span, len = 10 * span;
        if (len > 63 || 5 * span - 2 > MaxRegionBits) return (null, null);
        var byMask = ByMask(norm, len, ring: false);
        ulong all = (1UL << len) - 1;
        var masks = byMask.Keys.Where(m => m != 0).ToList();
        // the edge v -- v⊕f of Γ: absent when both entries vanish; Equal false when their moduli differ
        (bool Present, bool Equal, GaussBig Num, GaussBig Den) Edge(ulong v, ulong f)
        {
            var a = ElementOf(byMask, v ^ f, f);                   // <v|H|v⊕f>
            var b = ElementOf(byMask, all ^ (v ^ f), f);           // <v̄|H|v̄⊕f>
            if (a.IsZero && b.IsZero) return (false, true, a, b);
            return (true, a.Norm == b.Norm, a, b);
        }
        int Lo(ulong m) => BitOperations.TrailingZeroCount(m);
        int Hi(ulong m) => 63 - BitOperations.LeadingZeroCount(m);
        int? chainFrom = null, ringFrom = null;
        // a bulk placement is a translate of another, so the left end, one bulk stretch past it and the right end hold
        // every kind of cycle: f is taken within 4·span of either end
        foreach (var f in masks.Where(m => Lo(m) <= 4 * span || Hi(m) >= len - 1 - 4 * span))
            foreach (var g in masks)
            {
                if (f == g) continue;
                ulong both = f | g;
                if (Math.Max(Lo(f), Lo(g)) - Math.Min(Hi(f), Hi(g)) > span) continue;     // out of each other's reach
                int lo = Math.Max(0, Lo(both) - (span - 1)), hi = Math.Min(len - 1, Hi(both) + (span - 1));
                var region = Enumerable.Range(lo, hi - lo + 1).ToArray();
                ulong k = f ^ g;
                bool triangle = k != 0 && byMask.ContainsKey(k);
                bool bad = false;
                for (long cfg = 0; cfg < 1L << region.Length && !bad; cfg++)
                {
                    ulong u = 0;
                    for (int i = 0; i < region.Length; i++) if (((cfg >> i) & 1) == 1) u |= 1UL << region[i];
                    var e1 = Edge(u, f); var e2 = Edge(u ^ f, g); var e3 = Edge(u, g); var e4 = Edge(u ^ g, f);
                    if (e1.Present && e2.Present && e3.Present && e4.Present)
                    {
                        if (!(e1.Equal && e2.Equal && e3.Equal && e4.Equal)) bad = true;
                        else if (e1.Num * e2.Num * e3.Den * e4.Den != e3.Num * e4.Num * e1.Den * e2.Den) bad = true;
                    }
                    if (!bad && triangle)
                    {
                        var e5 = Edge(u, k);
                        if (e1.Present && e2.Present && e5.Present)
                        {
                            if (!(e1.Equal && e2.Equal && e5.Equal)) bad = true;
                            else if (e1.Num * e2.Num * e5.Den != e5.Num * e1.Den * e2.Den) bad = true;
                        }
                    }
                }
                if (!bad) continue;
                // the smallest chain holding this placement. Every window whose flip set is F, G or F ⊕ G lies inside R, so
                // a chain holding R whole holds the cycle: an end placement keeps its distance to its end, a bulk one needs
                // |R| sites. On a ring of N ≥ |R| no window wraps onto R's flip sets (an arc covering a flip set the other
                // way round would need at least N − extent(F) + 2 ≥ span + 1 sites), so the same cycle stands there.
                bool leftEnd = Lo(both) - (span - 1) < 0, rightEnd = Hi(both) + (span - 1) > len - 1;
                int extent = Hi(both) - Lo(both) + 1;
                int bulkN = extent + 2 * (span - 1);
                int chainN = leftEnd ? Hi(both) + span : rightEnd ? len - Lo(both) + span - 1 : bulkN;
                chainFrom = Math.Min(chainFrom ?? int.MaxValue, chainN);
                if (!leftEnd && !rightEnd) ringFrom = Math.Min(ringFrom ?? int.MaxValue, Math.Max(norm.SmallestRing, bulkN));
            }
        return (chainFrom, ringFrom);
    }
}
