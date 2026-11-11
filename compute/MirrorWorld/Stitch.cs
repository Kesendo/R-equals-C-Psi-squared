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
// condition (iii), the holonomy, is read on every cycle inside a local region (LocalHolonomy): a cycle whose ratios do
// not multiply to 1 breaks the component too, and a cycle in a local region recurs at every larger N, so that no holds
// for every N from the smallest chain or ring that holds it. A cycle round the whole ring is not local and is
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

    /// <summary>The most sites of a local region the holonomy walks (2^MaxHolonomySites configurations per region).</summary>
    public const int MaxHolonomySites = 12;

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

    // The local holonomy (Theorem 2's conditions (i) and (iii) on a local region). A region R is a run of r sites of a
    // long chain; its vertices are the configurations of R with every other bit 0, its edges the flip sets F inside R
    // all of whose terms (each window carrying a term with flip set F) lie inside R, with the ratio
    // r(v, F) = H_(v, v⊕F) / H_(v̄, v̄⊕F). Such an edge reads only bits of R, the bits outside it entering the far
    // configuration v̄ as 1 and nothing else. Walking R breadth first sets a potential along a spanning tree of each
    // component, d(v ⊕ F) = d(v)·r(v, F); an edge whose two entries differ in modulus, or a non-tree edge where the
    // potential does not close, is a broken condition (i) or a cycle of Γ with holonomy ≠ 1: the component holding it
    // is not good, and the spectrum does not pair. Every cycle of the region's graph is read, whatever its length (a
    // cycle through a flip set whose terms stick out of R is not, nor a cycle round the whole ring). The potential is
    // carried modulo a prime p ≡ 1 (mod 4), i a square root of −1 there: every entry on an edge is a nonzero Gaussian
    // integer of norm below p (the coefficients' absolute sum s bounds each entry, and 2s² < p is required, else the
    // reading is silent), so reduction is a ring map on the ratios, and a potential that does not close modulo p does
    // not close exactly (the converse may fail; a cycle missed that way is a silence, not a yes).
    // Where it recurs. A region at the left end of the chain is the chain of r sites, with the same terms on every edge
    // read; at the right end the same; a region in the bulk recurs on every chain of at least r sites, and on every
    // ring of N ≥ max(r, span + e − 1) sites, e the largest extent of a flip set read (a window covering F the long way
    // round the ring needs at least N − e + 2 sites). A region is a subgraph of every larger one of its kind, so the
    // first r at which a region breaks is that kind's threshold: the no holds on every chain of N ≥ r and on every ring
    // past the bulk's. Returns the two thresholds; null where no region up to MaxHolonomySites sites breaks, and for a
    // span past 9 or coefficients too large for the prime, which decides nothing.
    public static (int? Chain, int? Ring) LocalHolonomy(Theodolite.Cell cell)
    {
        var norm = cell.Normalized();
        if (!ZOnly(norm)) return (null, null);
        int span = norm.Span, cap = MaxHolonomySites, len = 2 * cap + 4 * span;
        if (len > 63) return (null, null);
        if (!SmallEntries(norm)) return (null, null);
        var byFlip = Placed(norm, len).Where(t => t.X != 0).GroupBy(t => t.X).ToDictionary(g => g.Key, g => g.ToList());
        ulong all = (1UL << len) - 1;
        long p = ModP.Primes[0], iota = ModP.SqrtMinusOne(p);

        long Reduce((long Re, long Im) g) => ModP.AddMod(ModP.Mod(g.Re, p), ModP.MulMod(ModP.Mod(g.Im, p), iota, p), p);

        // (broken, the largest extent of a flip set read) for the region of r sites from lo
        (bool Bad, int Extent) Region(int lo, int r)
        {
            ulong rm = ((1UL << r) - 1) << lo;
            var edges = byFlip.Where(kv => (kv.Key & ~rm) == 0 && kv.Value.All(t => (t.Window & ~rm) == 0)).ToList();
            int extent = edges.Count == 0 ? 0
                : edges.Max(kv => 64 - BitOperations.LeadingZeroCount(kv.Key) - BitOperations.TrailingZeroCount(kv.Key));
            int size = 1 << r;
            var num = new long[size];
            var den = new long[size];
            var seen = new bool[size];
            var queue = new Queue<int>();
            for (int root = 0; root < size; root++)
            {
                if (seen[root]) continue;
                seen[root] = true; num[root] = 1; den[root] = 1;
                queue.Enqueue(root);
                while (queue.Count > 0)
                {
                    int cu = queue.Dequeue();
                    ulong v = (ulong)cu << lo;
                    foreach (var (f, terms) in edges)
                    {
                        var a = Entry(terms, v ^ f);                     // <v|H|v⊕F>
                        var b = Entry(terms, all ^ v ^ f);               // <v̄|H|v̄⊕F>
                        if (a == (0, 0) && b == (0, 0)) continue;
                        long na = a.Re * a.Re + a.Im * a.Im, nb = b.Re * b.Re + b.Im * b.Im;
                        if (na != nb) return (true, extent);             // condition (i), zero against nonzero included
                        long ra = Reduce(a), rb = Reduce(b);
                        int cw = (int)((v ^ f) >> lo);
                        long nw = ModP.MulMod(num[cu], ra, p), dw = ModP.MulMod(den[cu], rb, p);
                        if (!seen[cw]) { seen[cw] = true; num[cw] = nw; den[cw] = dw; queue.Enqueue(cw); }
                        else if (ModP.MulMod(num[cw], dw, p) != ModP.MulMod(nw, den[cw], p)) return (true, extent);
                    }
                }
            }
            return (false, extent);
        }

        int? chain = null, ring = null;
        for (int r = 1; r <= cap && ring is null; r++)
        {
            if (chain is null && (Region(0, r).Bad || Region(len - r, r).Bad)) chain = r;
            var bulk = Region((len - r) / 2, r);
            if (!bulk.Bad) continue;
            chain ??= r;
            ring = Math.Max(norm.SmallestRing, Math.Max(r, span + bulk.Extent - 1));
        }
        return (chain, ring);
    }

    // each entry sums at most one placement of each template and field, so |Re|, |Im| ≤ s and its norm ≤ 2s² < p
    static bool SmallEntries(Theodolite.Cell norm)
        => norm.Windows.Sum(w => Math.Abs(w.Coefficient)) + norm.Fields.Sum(f => Math.Abs(f.Coefficient)) <= 20000;

    // every term of the cell on a chain of len sites, with its window (the run from its first to its last offset)
    static List<(ulong X, ulong Z, long C, ulong Window)> Placed(Theodolite.Cell norm, int len)
    {
        var placed = new List<(ulong X, ulong Z, long C, ulong Window)>();
        void Add(int first, string letters, int[] offsets, long c)
        {
            ulong x = 0, z = 0, w = 0;
            for (int o = offsets.Min(); o <= offsets.Max(); o++) w |= 1UL << (first + o);
            for (int i = 0; i < offsets.Length; i++)
            {
                ulong bit = 1UL << (first + offsets[i]);
                if (letters[i] is 'X' or 'Y') x |= bit;
                if (letters[i] is 'Z' or 'Y') z |= bit;
            }
            placed.Add((x, z, c, w));
        }
        for (int l = 0; l < len; l++)
            foreach (var (letter, c) in norm.Fields) Add(l, letter.ToString(), new[] { 0 }, c);
        foreach (var (letters, off, c) in norm.Windows)
            for (int s = 0; s + off.Max() < len; s++) Add(s, letters, off, c);
        return placed;
    }

    // <row|H|row ⊕ F> summed over the terms with flip set F (row ⊕ F = col), a Gaussian integer of longs
    static (long Re, long Im) Entry(List<(ulong X, ulong Z, long C, ulong Window)> terms, ulong col)
    {
        long re = 0, im = 0;
        foreach (var t in terms)
            switch ((BitOperations.PopCount(t.X & t.Z) + 2 * BitOperations.PopCount(col & t.Z)) % 4)
            {
                case 0: re += t.C; break;
                case 1: im += t.C; break;
                case 2: re -= t.C; break;
                default: im -= t.C; break;
            }
        return (re, im);
    }

    // The free chain (built 2026-11-11, the window family's chain residue X·(X ± Y)·X read on the letter cube and turned
    // into a reader). Where the local check is complete its silence is a yes. Four conditions, each read on a long chain:
    //   one flip shape: every term that flips anything flips a translate of one set f (XXX and XYX share {0, 1, 2});
    //   no edge vanishes: on every flip set, at every configuration of the bits its terms read, both entries of (★) are
    //     nonzero, with equal moduli (condition (i));
    //   condition (ii): no diagonal string of odd weight keeps a coefficient;
    //   every square closes: for two flip sets F, G whose terms read a bit the other flips (as a Z or Y letter), at every
    //     configuration of the bits their terms cover, r(v, F)·r(v ⊕ F, G) = r(v, G)·r(v ⊕ G, F), exactly, in Gaussian
    //     integers (no prime: a closure modulo p is not a closure).
    // Then on the chain the translates of f are independent (each begins at its own site), so every component of Γ is
    // a whole cube of them, every edge present, and its cycles are generated by squares (a step there and back closes,
    // H being Hermitian and the moduli equal). A square of two flip sets that do not interact closes by itself; an
    // interacting pair's terms lie in a run of at most 3·span − 2 − e sites (e = max f − min f), and a pair an end cuts
    // within the first or last 3·span − 3 − e. From N ≥ max(2·span, 3·span − 3) no flip set is cut at both ends and no
    // interacting pair at opposite ends, so every square of the chain of N sites is one read on the long chain, with the
    // same entries: every component flat, the spectrum pairs on every chain of N ≥ From. On the ring the translates span
    // the ideal (f) of GF(2)[x]/(x^N + 1), of dimension N − deg gcd(f, x^N + 1): they are independent exactly when no
    // order of a root of f divides N (the Moduli, the minimal such N). There, from N ≥ 3·span − 1 (no interacting pair
    // wraps onto itself, no window covers a flip set the long way round), every square is one of the long chain's bulk,
    // and the same reading holds; where a modulus divides N a relation runs round the whole ring and its holonomy is not
    // read (silent). The section it proves is no product and no local quadratic form: in the coordinates of the cube
    // (v = v₀ ⊕ Σ c_l·F_l) it is quadratic in c, and c is not local in v (for X·(X ± Y)·X, 1/(1 + x + x²) =
    // (1 + x)/(1 + x³) over GF(2), sums of every third bit from the end).
    public sealed record Free(bool Ring, int From, IReadOnlyList<int> Moduli)
    {
        public bool Reads(int n) => n >= From && !Moduli.Any(d => n % d == 0);
    }

    /// <summary>The most bits a square of the free chain enumerates (2^MaxSquareBits configurations per pair).</summary>
    public const int MaxSquareBits = 16;

    public static Free? FreeEveryN(Theodolite.Cell cell, bool ring)
    {
        var norm = cell.Normalized();
        if (!ZOnly(norm)) return null;
        int span = norm.Span, len = 8 * span;
        if (len > 63) return null;
        var placed = Placed(norm, len);
        var byFlip = placed.Where(t => t.X != 0).GroupBy(t => t.X).ToDictionary(g => g.Key, g => g.ToList());
        if (byFlip.Count == 0) return null;
        static ulong Shape(ulong m) => m >> BitOperations.TrailingZeroCount(m);
        ulong f = Shape(byFlip.Keys.First());
        if (byFlip.Keys.Any(m => Shape(m) != f)) return null;
        // condition (ii): H_xx = H_x̄x̄ for every x, the odd-weight diagonal strings cancelling
        if (placed.Where(t => t.X == 0).GroupBy(t => t.Z).Any(g => g.Sum(t => t.C) != 0 && BitOperations.PopCount(g.Key) % 2 == 1))
            return null;
        ulong all = (1UL << len) - 1;
        int Lo(ulong m) => BitOperations.TrailingZeroCount(m);
        int Hi(ulong m) => 63 - BitOperations.LeadingZeroCount(m);
        // a placement in the bulk: every window that could carry its flip sets fits on the long chain
        bool Bulk(ulong m) => Lo(m) >= span - 1 && Hi(m) <= len - span;
        GaussBig Big((long Re, long Im) g) => new(g.Re, g.Im);
        var reach = byFlip.ToDictionary(kv => kv.Key, kv => kv.Value.Aggregate(0UL, (acc, t) => acc | t.Window));
        var reads = byFlip.ToDictionary(kv => kv.Key, kv => kv.Value.Aggregate(0UL, (acc, t) => acc | t.Z));
        bool chainOk = true, ringOk = true;
        void Fail(ulong m) { chainOk = false; if (Bulk(m)) ringOk = false; }
        IEnumerable<ulong> Configurations(ulong bits)
        {
            var sites = Enumerable.Range(0, len).Where(l => ((bits >> l) & 1) == 1).ToArray();
            for (long cfg = 0; cfg < 1L << sites.Length; cfg++)
            {
                ulong v = 0;
                for (int i = 0; i < sites.Length; i++) if (((cfg >> i) & 1) == 1) v |= 1UL << sites[i];
                yield return v;
            }
        }
        foreach (var (mask, terms) in byFlip)
            foreach (var col in Configurations(reach[mask]))
            {
                var a = Entry(terms, col);
                var b = Entry(terms, all ^ col);
                if (a == (0, 0) || b == (0, 0) || a.Re * a.Re + a.Im * a.Im != b.Re * b.Re + b.Im * b.Im) { Fail(mask); break; }
            }
        if (!chainOk && !ringOk) return null;
        // r(v, F) = <v|H|v⊕F> / <v̄|H|v̄⊕F>, as the pair (numerator, denominator)
        (GaussBig A, GaussBig B) R(ulong v, ulong m) => (Big(Entry(byFlip[m], v ^ m)), Big(Entry(byFlip[m], all ^ v ^ m)));
        var masks = byFlip.Keys.OrderBy(m => m).ToList();
        foreach (var fm in masks)
            foreach (var gm in masks)
            {
                if (gm <= fm || ((gm & reads[fm]) == 0 && (fm & reads[gm]) == 0)) continue;
                ulong both = fm | gm;
                if (!(chainOk || ringOk && Bulk(both))) continue;
                ulong bits = reach[fm] | reach[gm];
                if (BitOperations.PopCount(bits) > MaxSquareBits) return null;
                foreach (var v in Configurations(bits))
                {
                    var (a1, b1) = R(v, fm);
                    var (a2, b2) = R(v ^ fm, gm);
                    var (a3, b3) = R(v, gm);
                    var (a4, b4) = R(v ^ gm, fm);
                    if (a1 * a2 * b3 * b4 != a3 * a4 * b1 * b2) { Fail(both); break; }
                }
            }
        if (!ring) return chainOk ? new Free(false, Math.Max(2 * span, 3 * span - 3), Array.Empty<int>()) : null;
        if (!ringOk) return null;
        // the minimal N with gcd(f, x^N + 1) ≠ 1: every root of f lies in GF(2^m), m ≤ deg f, of order dividing 2^m − 1
        int deg = Hi(f);
        var moduli = new List<int>();
        for (int d = 1; d < 1 << Math.Max(deg, 1); d++)
            if (!moduli.Any(m => d % m == 0) && Gcd2(f, (1UL << d) | 1) != 1) moduli.Add(d);
        if (moduli.Contains(1)) return null;                             // f(1) = 0: the translates sum to 0 on every ring
        return new Free(true, Math.Max(norm.SmallestRing, 3 * span - 1), moduli);
    }

    // the gcd of two polynomials over GF(2), bit i the coefficient of x^i
    static ulong Gcd2(ulong a, ulong b)
    {
        while (b != 0)
        {
            int db = 63 - BitOperations.LeadingZeroCount(b);
            while (a != 0 && 63 - BitOperations.LeadingZeroCount(a) >= db) a ^= b << (63 - BitOperations.LeadingZeroCount(a) - db);
            (a, b) = (b, a);
        }
        return a;
    }

    // ---------------------------------------------------------------- the ring laws (F173, adopted)
    //
    // Adopted from docs/proofs/PROOF_SIGNED_MASK_BALANCE.md, not re-derived: two laws in N for a cell on the ring under
    // Z dephasing on every site, from the guard N ≥ 2·span − 1 on (and the smallest ring), where the placements are
    // distinct flip sets. Theorem D: every flip-carrying template of its own shape (one string per flip set), read through
    // the (1 + x)-valuations v_i of the masks and the signs s_i = (#Y + #Z) mod 2 alone; with v₀ the least and T = 2^(v₂(N))
    // the ring breaks exactly when T ≤ v₀ and some sign is odd, or T > v₀ and two templates of valuation v₀ differ in sign
    // or one of greater valuation has an odd sign. Corollary E2: every flip-carrying template of one shape f, each flip set
    // a constant times a character of x, ω·(−1)^⟨R, x⟩; the ring pairs exactly when the read sets are read symmetrically
    // and h = gcd(f, x^N + 1) is 1, or h divides R and either 1 + x does not divide h or the relation (x^N + 1)/h walks
    // to 1. Both are laws for every N; the code evaluates them at N up to 63. A cell with an odd-weight diagonal, a flip
    // set whose entries differ in modulus or vanish somewhere, a ratio that is no character, or several shapes one of
    // which carries several strings has no law here (the stitch's no reads the first two).
    public sealed record RingLawReading(int Guard, string Text, Func<int, bool> Pairs);

    public static RingLawReading? RingLaw(Theodolite.Cell cell)
    {
        var norm = cell.Normalized();
        if (!ZOnly(norm)) return null;
        var templates = new List<(ulong X, ulong Z, long C)>();
        foreach (var (letters, off, c) in norm.Windows)
        {
            ulong x = 0, z = 0;
            for (int i = 0; i < letters.Length; i++)
            {
                if (letters[i] is 'X' or 'Y') x |= 1UL << off[i];
                if (letters[i] is 'Y' or 'Z') z |= 1UL << off[i];
            }
            templates.Add((x, z, c));
        }
        foreach (var (letter, c) in norm.Fields)
            templates.Add((letter is 'X' or 'Y' ? 1UL : 0, letter is 'Y' or 'Z' ? 1UL : 0, c));
        // equal strings collected and zero coefficients dropped, so that the laws read the Hamiltonian, not its spelling
        templates = templates.GroupBy(t => (t.X, t.Z)).Select(g => (g.Key.X, g.Key.Z, C: g.Sum(t => t.C))).Where(t => t.C != 0).ToList();
        if (templates.Any(t => t.X == 0 && (BitOperations.PopCount(t.Z) & 1) == 1)) return null;
        var flipping = templates.Where(t => t.X != 0).ToList();
        if (flipping.Count == 0) return null;
        int guard = Math.Max(norm.SmallestRing, 2 * norm.Span - 1);
        static ulong ShapeOf(ulong m) => m >> BitOperations.TrailingZeroCount(m);
        int shapes = flipping.Select(t => ShapeOf(t.X)).Distinct().Count();
        if (shapes == flipping.Count)
        {
            // Theorem D
            var read = flipping.Select(t => (V: Valuation(t.X), S: BitOperations.PopCount(t.Z) & 1)).ToList();
            int v0 = read.Min(r => r.V);
            bool low = read.Any(r => r.S == 1);
            bool high = read.Where(r => r.V == v0).Select(r => r.S).Distinct().Count() > 1 || read.Any(r => r.V > v0 && r.S == 1);
            string text = low == high
                ? $"F173 Theorem D: {(high ? "breaks" : "pairs")} on every ring of N ≥ {guard}"
                : $"F173 Theorem D: {(high ? "breaks" : "pairs")} exactly where 2^v₂(N) > {v0}, from N = {guard}";
            return new RingLawReading(guard, text, n =>
            {
                if (n < guard) throw new ArgumentOutOfRangeException(nameof(n), n, $"the law holds from N = {guard}");
                return (n & -n) > v0 ? !high : !low;
            });
        }
        if (shapes != 1) return null;
        // Corollary E2: the ratio's kind and read set are the same on every ring past the guard; read them once
        ulong f = ShapeOf(flipping[0].X);
        int probe = Math.Max(guard, 3 * norm.Span);
        var first = FlipRatios(norm, probe);
        if (first is null) return null;
        return new RingLawReading(guard, $"F173 Corollary E2: one flip shape, read through gcd(f, x^N + 1), f = {Poly(f)}, from N = {guard}", n =>
        {
            if (n < guard || n > 63) throw new ArgumentOutOfRangeException(nameof(n), n, $"the law is evaluated from N = {guard} to 63");
            var ratios = FlipRatios(norm, n) ?? throw new InvalidOperationException($"N = {n}: a ratio left the character form past the guard");
            static int Parity(ulong v) => BitOperations.PopCount(v) & 1;
            foreach (var (pm, (_, pr)) in ratios)
                foreach (var (qm, (_, qr)) in ratios)
                    if (Parity(pr & qm) != Parity(qr & pm)) return false;
            ulong ringPoly = (1UL << n) | 1UL, h = Gcd2(f, ringPoly);
            if (h == 1) return true;
            var (omega0, read0) = ratios[f];
            if (Mod2(read0, h) != 0) return false;
            if ((BitOperations.PopCount(h) & 1) == 1) return true;
            ulong g = Div2(ringPoly, h), at = 0;
            int quarter = 0;
            for (int p = 0; p < n; p++)
            {
                if (((g >> p) & 1) == 0) continue;
                ulong m = Rot(f, p, n);
                var (om, rd) = ratios[m];
                quarter += om + 2 * Parity(rd & at);
                at ^= m;
            }
            return quarter % 4 == 0;
        });
    }

    // every flip set's ratio on the ring of n sites as (quarter, read set): r(x) = H_(x̄⊕m, x̄)/H_(x⊕m, x) = i^quarter·(−1)^⟨R, x⟩,
    // the proof's ratio conjugated (ElementOf reads H_(x⊕m, x)), the same read set, and every verdict the same,
    // read on the bits its strings' Z letters touch; null when some flip set is not of that form (entries of different
    // modulus, a vanishing entry, no character)
    static Dictionary<ulong, (int Quarter, ulong Read)>? FlipRatios(Theodolite.Cell norm, int n)
    {
        var byMask = ByMask(norm, n, ring: true);
        var result = new Dictionary<ulong, (int, ulong)>();
        foreach (var (mask, list) in byMask)
        {
            if (mask == 0) continue;
            ulong u = 0;
            foreach (var t in list) u |= t.Z;
            var bits = Enumerable.Range(0, n).Where(j => ((u >> j) & 1) == 1).ToArray();
            if (bits.Length > 20) return null;
            ulong Spread(int y) { ulong x = 0; for (int k = 0; k < bits.Length; k++) if (((y >> k) & 1) == 1) x |= 1UL << bits[k]; return x; }
            static GaussBig Turn(GaussBig a, int q) { for (int k = 0; k < q; k++) a = a * new GaussBig(0, 1); return a; }
            int count = 1 << bits.Length;
            var a = new GaussBig[count];
            var b = new GaussBig[count];
            for (int y = 0; y < count; y++)
            {
                a[y] = ElementOf(byMask, Spread(y), mask);
                b[y] = ElementOf(byMask, u ^ Spread(y), mask);
                if (a[y].IsZero || b[y].IsZero || a[y].Norm != b[y].Norm) return null;
            }
            int quarter = Enumerable.Range(0, 4).FirstOrDefault(q => b[0] == Turn(a[0], q), -1);
            if (quarter < 0) return null;
            ulong read = 0;
            for (int k = 0; k < bits.Length; k++)
            {
                var want = Turn(a[1 << k], quarter);
                if (b[1 << k] == want) continue;
                if (b[1 << k] == -want) { read |= 1UL << bits[k]; continue; }
                return null;
            }
            for (int y = 0; y < count; y++)
            {
                var want = Turn(a[y], quarter);
                if ((BitOperations.PopCount(Spread(y) & read) & 1) == 1) want = -want;
                if (b[y] != want) return null;
            }
            result[mask] = (quarter, read);
        }
        return result;
    }

    // v_(1+x)(m): how often 1 + x divides the mask read as a polynomial over GF(2)
    static int Valuation(ulong m)
    {
        int v = 0;
        while (m != 0 && (BitOperations.PopCount(m) & 1) == 0) { m = Div2(m, 0b11); v++; }
        return v;
    }

    static ulong Mod2(ulong a, ulong b)
    {
        int db = 63 - BitOperations.LeadingZeroCount(b);
        while (a != 0 && 63 - BitOperations.LeadingZeroCount(a) >= db) a ^= b << (63 - BitOperations.LeadingZeroCount(a) - db);
        return a;
    }

    static ulong Div2(ulong a, ulong b)
    {
        ulong q = 0;
        int db = 63 - BitOperations.LeadingZeroCount(b);
        while (a != 0 && 63 - BitOperations.LeadingZeroCount(a) >= db)
        {
            int s = 63 - BitOperations.LeadingZeroCount(a) - db;
            q |= 1UL << s;
            a ^= b << s;
        }
        if (a != 0) throw new InvalidOperationException("not divisible");
        return q;
    }

    static ulong Rot(ulong m, int p, int n)
    {
        ulong r = 0;
        for (int j = 0; j < n; j++) if (((m >> j) & 1) == 1) r |= 1UL << ((j + p) % n);
        return r;
    }

    static string Poly(ulong f) => string.Join(" + ", Enumerable.Range(0, 64).Where(j => ((f >> j) & 1) == 1)
        .Select(j => j == 0 ? "1" : j == 1 ? "x" : $"x^{j}"));
}
