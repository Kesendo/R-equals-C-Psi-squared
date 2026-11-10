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
// the two are one equation: an edge with unequal moduli (one of them zero among them), or a diagonal entry with
// H_xx ≠ H_x̄x̄, silences it. Where it
// reads, W = X^⊗N·D is an invertible element of the far end (every factor nonzero), so W·(near end) lies in the far end
// and W⁻¹·(far end) in the near end: the two ends have one dimension and the spectrum pairs (F158).
//
// Every N. An equation at the edge (x, x ⊕ F), F the flip set of a window, changes only the factors that touch F and
// reads only the terms with flip set F, whose Z letters lie within span − 1 of F: it depends on the bits of at most
// 3·span − 2 sites around F and on nothing else. On a ring of N ≥ 3·span sites that region never meets itself, every
// configuration of it occurs, and no other: the set of distinct equations is the same at every such N, so one solve at
// N = 3·span reads every larger ring. On the chain the regions within reach of an end carry the end, and from
// N = 6·span on the two ends' regions and a bulk region are apart: one solve reads every longer chain. Both are
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

    // the every-N reading of one geometry: Bulk the system at the window (null when silent there; it holds at every N
    // from the window on), Small the N below the window, from the smallest, where a fixed-N solve reads
    public sealed record Reach(bool Ring, int Window, Binomials.Reduced? Bulk, IReadOnlyList<int> Small)
    {
        public bool Reads(int n) => n >= Window ? Bulk is not null : Small.Contains(n);
    }

    /// <summary>The most sites the every-N reading enumerates (2^MaxSites basis states at the window and one past it):
    /// the chain reaches it at span 2, the ring at span 5.</summary>
    public const int MaxSites = 16;

    public static int WindowOf(Theodolite.Cell cell, bool ring) => (ring ? 3 : 6) * cell.Normalized().Span;

    // The cell on every N at once; null when the cell's jumps are not Z on every site, or when its window and one site
    // more pass MaxSites (not read, which decides nothing).
    public static Reach? EveryN(Theodolite.Cell cell, bool ring)
    {
        var norm = cell.Normalized();
        if (!ZOnly(norm)) return null;
        int window = WindowOf(norm, ring);
        if (window + 1 > MaxSites) return null;
        var at = Equations(norm, window, ring);
        var next = Equations(norm, window + 1, ring);
        if (at.Contradiction != next.Contradiction || !at.Contradiction && !at.Keys.SetEquals(next.Keys))
            throw new InvalidOperationException($"the equations at N = {window} and {window + 1} differ: the window is not reached");
        var bulk = at.Contradiction ? null : Binomials.Solve(at.Equations, Unknowns(norm.Span));
        int from = ring ? norm.SmallestRing : Math.Max(2, norm.Span);
        var small = Enumerable.Range(from, Math.Max(0, window - from)).Where(n => At(norm, n, ring) is not null).ToList();
        return new Reach(ring, window, bulk, small);
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
        var pairs = new List<(int A, int B, int O)>();
        for (int o = 1; o < span; o++)
            for (int l = 0; l < n; l++)
                if (ring || l + o < n) pairs.Add((l, (l + o) % n, o));
        var byMask = terms.GroupBy(t => t.X).ToDictionary(g => g.Key, g => g.ToList());
        ulong all = n == 64 ? ulong.MaxValue : (1UL << n) - 1;

        // <row|P_mask|col> summed over the terms with that flip set: P|col> = i^|x∧z|·(−1)^|col∧z|·|col ⊕ x> (Y = iXZ)
        GaussBig Element(ulong col, ulong mask)
        {
            BigInteger re = 0, im = 0;
            foreach (var t in byMask[mask])
                switch ((BitOperations.PopCount(t.X & t.Z) + 2 * BitOperations.PopCount(col & t.Z)) % 4)
                {
                    case 0: re += t.C; break;
                    case 1: im += t.C; break;
                    case 2: re -= t.C; break;
                    default: im -= t.C; break;
                }
            return new GaussBig(re, im);
        }
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
}
