using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live battery for Theorem F of F173 (<c>docs/proofs/PROOF_SIGNED_MASK_BALANCE.md</c>, typed
/// <see cref="SignedMaskBalanceClaim.TorusBreaks"/> and <see cref="SignedMaskBalanceClaim.TorusBreaksByDifferences"/>):
/// templates placed at every translate of the L × M torus, one string per flip set, Z dephasing on every site, read on the
/// 2-part T₁ × T₂ of the translations (the second form, a rearrangement of the first, computed beside it and no second
/// route). Held against Theorem 2 of the complement connection walked on the bitstrings
/// (<see cref="ComplementConnectionGraph"/>) on tori of up to 16 sites, 4 × 4 the one among them whose 2-part is not
/// cyclic, and against the balance on the whole torus (<see cref="SignedMaskBalanceClaim.Balanced"/>, Theorem A's
/// elimination on the unfolded placements) to 64 sites: the templates of the 2 × 2 window, the ring as the L × 1 torus
/// against Theorem D, the tori of odd L against Theorem D's one-liner read on the masks summed along x, and the worked
/// examples of the page. Three controls must disagree with the balance: Theorem D's one-liner carried over with one number
/// per mask (its least order in the powers of 1 + x and 1 + y) and T = min(T₁, T₂), read on the tori with T₁, T₂ ≥ 2, where A
/// is no chain, the fold with T₁ and T₂ swapped, and the fold onto the parities of L and M in place of
/// their 2-parts.</summary>
public sealed class TorusMaskBalanceWitness : IInspectable
{
    public sealed record FamilyReading(string Name, int Cells, int Readings, int Pairing, int Mismatches, string? FirstMismatch, string Against)
    {
        public bool Agrees => Mismatches == 0;
    }

    public IReadOnlyList<FamilyReading> Families { get; }
    public int OneNumberMismatches { get; }
    /// <summary>The readings on tori with T₁, T₂ ≥ 2, where the one-number control is read.</summary>
    public int OneNumberReadings { get; }
    /// <summary>The readings of the balance family, where the two fold controls are read.</summary>
    public int FoldControlReadings { get; }
    public int SwappedFoldMismatches { get; }
    public int ParityFoldMismatches { get; }
    public bool AllAgree => Families.All(f => f.Agrees);
    public bool ControlsBreak => OneNumberMismatches > 0 && SwappedFoldMismatches > 0 && ParityFoldMismatches > 0;

    /// <summary>The templates of the 2 × 2 window with a nonzero mask that touches the window's first column and first row
    /// (one template per letter pattern up to translation of its mask).</summary>
    public static IReadOnlyList<string> WindowTemplates { get; } =
        (from a in "IXYZ" from b in "IXYZ" from c in "IXYZ" from d in "IXYZ" select $"{a}{b}/{c}{d}")
        .Where(t =>
        {
            var flips = SignedMaskBalanceClaim.TorusLetters(t).Where(x => x.Letter is 'X' or 'Y').ToList();
            return flips.Count > 0 && flips.Any(x => x.I == 0) && flips.Any(x => x.J == 0);
        }).ToArray();

    static readonly (int L, int M)[] Tori =
        { (3, 3), (3, 4), (4, 3), (4, 4), (3, 8), (8, 3), (4, 6), (6, 4), (6, 6), (5, 8), (8, 4), (4, 8), (8, 8), (16, 4), (4, 16), (12, 5) };

    /// <summary>Theorem 2 on the L × M torus (site (i, j) at position i + L·j), each template at every translate, wrapped,
    /// coefficient 1.</summary>
    public static bool PairsOnTorus(int l, int m, params string[] templates)
    {
        var placed = new List<(PauliMask, BigInteger)>();
        foreach (var t in templates)
        {
            var letters = SignedMaskBalanceClaim.TorusLetters(t);
            for (int b = 0; b < m; b++)
                for (int a = 0; a < l; a++)
                {
                    var s = Enumerable.Repeat('I', l * m).ToArray();
                    foreach (var (i, j, c) in letters) s[(i + a) % l + l * ((j + b) % m)] = c;
                    placed.Add((PauliMask.Parse(new string(s)), BigInteger.One));
                }
        }
        var (near, far, _) = ComplementConnectionGraph.EverySiteDephased(l * m, ComplementConnectionGraph.Columns(l * m, placed));
        return near == far;
    }

    static bool Placeable(string[] ts, int l, int m)
    {
        try { SignedMaskBalanceClaim.OnTorus(ts, l, m); return true; }
        catch (ArgumentException) { return false; }
    }

    static bool BalanceBreaks(string[] ts, int l, int m) => !SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnTorus(ts, l, m));

    static int TwoPart(int n) => 1 << BitOperations.TrailingZeroCount((uint)n);

    // the fold onto t1 × t2 read as a linear system, whatever t1 and t2 are (the controls choose them wrongly)
    static bool FoldBreaks(string[] ts, int t1, int t2)
    {
        var rows = new List<(ulong, int)>();
        foreach (var t in ts)
            for (int b = 0; b < t2; b++)
                for (int a = 0; a < t1; a++) rows.Add((SignedMaskBalanceClaim.FoldedMask(t, t1, t2, a, b), SignedMaskBalanceClaim.Sign(t)));
        return !SignedMaskBalanceClaim.Balanced(rows);
    }

    // control: Theorem D's one-liner with one number per mask, the least total degree of the folded mask written in the
    // powers of s = 1 + x and t = 1 + y (Lucas: x^i holds s^a exactly when a ⊆ i), and T = min(T₁, T₂)
    static bool OneNumberBreaks(string[] ts, int l, int m)
    {
        int t1 = TwoPart(l), t2 = TwoPart(m), big = Math.Min(t1, t2);
        var read = ts.Select(t =>
        {
            ulong f = SignedMaskBalanceClaim.FoldedMask(t, t1, t2);
            int order = int.MaxValue;
            for (int a = 0; a < t1; a++)
                for (int b = 0; b < t2; b++)
                {
                    int coefficient = 0;
                    for (int k = 0; k < t1 * t2; k++)
                        if (((f >> k) & 1) == 1 && ((k % t1) & a) == a && ((k / t1) & b) == b) coefficient ^= 1;
                    if (coefficient == 1) order = Math.Min(order, a + b);
                }
            return (V: order, S: SignedMaskBalanceClaim.Sign(t));
        }).ToList();
        int v0 = read.Min(r => r.V);
        bool pairs = new[] { 0, 1 }.Any(sigma => read.All(r => r.S == sigma * (r.V == v0 && v0 < big ? 1 : 0)));
        return !pairs;
    }

    // odd L: Theorem D's one-liner on the ring of M with each mask summed along x, m(1, y), and T = 2^(v₂(M))
    static bool OddSideBreaks(string[] ts, int m)
    {
        int big = TwoPart(m);
        var read = ts.Select(t =>
        {
            ulong p = 0;
            foreach (var (_, j, c) in SignedMaskBalanceClaim.TorusLetters(t)) if (c is 'X' or 'Y') p ^= 1UL << j;
            return (V: p == 0 ? int.MaxValue : SignedMaskBalanceClaim.Valuation(p), S: SignedMaskBalanceClaim.Sign(t));
        }).ToList();
        int v0 = read.Min(r => r.V);
        bool pairs = new[] { 0, 1 }.Any(sigma => read.All(r => r.S == sigma * (r.V == v0 && v0 < big ? 1 : 0)));
        return !pairs;
    }

    public TorusMaskBalanceWitness(int seed = 20261111)
    {
        var rng = new Random(seed);
        var families = new List<FamilyReading>();
        var pool = WindowTemplates;
        int oneNumber = 0, oneNumberReadings = 0, foldReadings = 0, swapped = 0, parity = 0;

        // two window templates flipping the same cells place the same flip sets (anchored masks of different cells are not
        // translates of each other, and on the tori read here, L, M ≥ 3, their placements never coincide; OnTorus refuses any that would)
        static bool SameMask(string a, string b) => SignedMaskBalanceClaim.FoldedMask(a, 2, 2) == SignedMaskBalanceClaim.FoldedMask(b, 2, 2);
        var singles = pool.Select(t => new[] { t }).ToList();
        var pairs = new List<string[]>();
        for (int i = 0; i < pool.Count; i++)
            for (int j = i + 1; j < pool.Count; j++)
                if (!SameMask(pool[i], pool[j])) pairs.Add(new[] { pool[i], pool[j] });
        var triples = new List<string[]>();
        while (triples.Count < 1500)
        {
            var c = pool.OrderBy(_ => rng.Next()).Take(3).ToArray();
            if (!SameMask(c[0], c[1]) && !SameMask(c[0], c[2]) && !SameMask(c[1], c[2])) triples.Add(c);
        }

        // against Theorem 2 on the bitstrings: every single template, 150 pairs and 60 triples (random, seeded) on 3 × 3, 3 × 4
        // and 4 × 3, whose 2-parts are cyclic, and every single template, 40 of those pairs and 20 of those triples on 4 × 4,
        // whose 2-part ℤ₄ × ℤ₄ is not
        {
            int cells = 0, readings = 0, pairing = 0, bad = 0;
            string? first = null;
            var drawnPairs = pairs.OrderBy(_ => rng.Next()).Take(150).ToList();
            var cellsRead = singles.Concat(drawnPairs).Concat(triples.Take(60)).ToList();
            var atFour = new HashSet<string[]>(singles.Concat(drawnPairs.Take(40)).Concat(triples.Take(20)));
            // the readings listed first, Theorem 2 walked on every core (each walk on its own, the graph routine holds no shared
            // state), then counted in list order so the first mismatch named is the same on every run
            var jobs = new List<(string[] Ts, int L, int M)>();
            foreach (var ts in cellsRead)
            {
                bool any = false;
                foreach (var (l, m) in new[] { (3, 3), (3, 4), (4, 3), (4, 4) })
                {
                    if (l == 4 && m == 4 && !atFour.Contains(ts)) continue;
                    if (!Placeable(ts, l, m)) continue;
                    any = true;
                    jobs.Add((ts, l, m));
                }
                if (any) cells++;
            }
            var graphs = new bool[jobs.Count];
            Parallel.For(0, jobs.Count, k => graphs[k] = PairsOnTorus(jobs[k].L, jobs[k].M, jobs[k].Ts));
            for (int k = 0; k < jobs.Count; k++)
            {
                var (ts, l, m) = jobs[k];
                bool graph = graphs[k];
                bool law = !SignedMaskBalanceClaim.TorusBreaks(ts, l, m);
                bool second = !SignedMaskBalanceClaim.TorusBreaksByDifferences(ts, l, m);
                readings++;
                if (graph) pairing++;
                if (graph != law || graph != second) { bad++; first ??= $"{string.Join(" + ", ts)} on {l} × {m}"; }
            }
            families.Add(new FamilyReading("templates of the 2 × 2 window on 3 × 3, 3 × 4, 4 × 3 and 4 × 4 (the second form beside the first)", cells, readings, pairing, bad, first, "Theorem 2"));
        }

        // against the balance on the whole torus: every single, every pair of different masks, 1500 random triples (seeded,
        // repeats possible), sixteen tori to 64 sites
        {
            int cells = 0, readings = 0, pairing = 0, bad = 0;
            string? first = null;
            foreach (var ts in singles.Concat(pairs).Concat(triples))
            {
                bool any = false;
                foreach (var (l, m) in Tori)
                {
                    if (!Placeable(ts, l, m)) continue;
                    any = true;
                    bool balance = !BalanceBreaks(ts, l, m);
                    bool law = !SignedMaskBalanceClaim.TorusBreaks(ts, l, m);
                    bool second = !SignedMaskBalanceClaim.TorusBreaksByDifferences(ts, l, m);
                    readings++;
                    if (balance) pairing++;
                    if (balance != law || balance != second) { bad++; first ??= $"{string.Join(" + ", ts)} on {l} × {m}"; }
                    if (TwoPart(l) >= 2 && TwoPart(m) >= 2) { oneNumberReadings++; if (balance == OneNumberBreaks(ts, l, m)) oneNumber++; }
                    foldReadings++;
                    if (balance == FoldBreaks(ts, TwoPart(m), TwoPart(l))) swapped++;
                    if (balance == FoldBreaks(ts, l % 2 == 0 ? 2 : 1, m % 2 == 0 ? 2 : 1)) parity++;
                }
                if (any) cells++;
            }
            families.Add(new FamilyReading("templates of the 2 × 2 window, every single, every pair of different masks and 1500 triples, sixteen tori to 64 sites (the second form beside the first)", cells, readings, pairing, bad, first, "the balance"));
        }

        // the tori of odd L: Theorem D's one-liner on the masks summed along x, the ring of M
        {
            int cells = 0, readings = 0, pairing = 0, bad = 0;
            string? first = null;
            foreach (var ts in singles.Concat(pairs).Concat(triples))
            {
                bool any = false;
                foreach (var (l, m) in Tori.Where(t => t.L % 2 == 1))
                {
                    if (!Placeable(ts, l, m)) continue;
                    any = true;
                    bool balance = !BalanceBreaks(ts, l, m);
                    readings++;
                    if (balance) pairing++;
                    if (balance == OddSideBreaks(ts, m)) { bad++; first ??= $"{string.Join(" + ", ts)} on {l} × {m}"; }
                }
                if (any) cells++;
            }
            families.Add(new FamilyReading("odd L: Theorem D's one-liner on m(1, y), the ring of M", cells, readings, pairing, bad, first, "the balance"));
        }

        // the ring as the L × 1 torus: Theorem F against Theorem D, two templates of different shapes up to three sites wide
        {
            var ring = new[] { "XX", "XY", "YY", "XZ", "ZY", "YZ" }
                .Concat(from a in "XYZ" from mid in "IXYZ" from b in "XYZ" select $"{a}{mid}{b}")
                .Where(t => SignedMaskBalanceClaim.Mask(t) != 0).Distinct().ToArray();
            int cells = 0, readings = 0, pairing = 0, bad = 0;
            string? first = null;
            for (int i = 0; i < ring.Length; i++)
                for (int j = i + 1; j < ring.Length; j++)
                {
                    ulong ma = SignedMaskBalanceClaim.Mask(ring[i]), mb = SignedMaskBalanceClaim.Mask(ring[j]);
                    if (ma >> BitOperations.TrailingZeroCount(ma) == mb >> BitOperations.TrailingZeroCount(mb)) continue;
                    var ts = new[] { ring[i], ring[j] };
                    cells++;
                    for (int n = 5; n <= 32; n++)
                    {
                        bool d = !SignedMaskBalanceClaim.RingBreaks(ts, n);
                        bool law = !SignedMaskBalanceClaim.TorusBreaks(ts, n, 1);
                        readings++;
                        if (d) pairing++;
                        if (d != law) { bad++; first ??= $"{ring[i]} + {ring[j]} on the ring of {n}"; }
                    }
                }
            families.Add(new FamilyReading("the ring as the L × 1 torus, two templates of different shapes up to three sites wide, L = 5 to 32", cells, readings, pairing, bad, first, "Theorem D"));
        }

        // the worked examples of the page: Theorem 2 to 4 × 4, the law on the tori to 64 sites with L, M ≥ 3
        {
            var small = new[] { (3, 3), (3, 4), (4, 3), (4, 4) };
            var large = (from l in Enumerable.Range(3, 19) from m in Enumerable.Range(3, 19) where l * m <= 64 select (l, m)).ToArray();
            var examples = new (string Name, string[] Ts, Func<int, int, bool> Pairs)[]
            {
                ("XX across beside XY down pairs exactly at even M", new[] { "XX", "X/Y" }, (l, m) => m % 2 == 0),
                ("the plaquette XX/XY pairs exactly at even L and M", new[] { "XX/XY" }, (l, m) => l % 2 == 0 && m % 2 == 0),
                ("XY across beside XY down pairs exactly at even L and M (the checkerboard)", new[] { "XY", "X/Y" }, (l, m) => l % 2 == 0 && m % 2 == 0),
                ("XX across, XY down, XY on the anti-diagonal pairs exactly at even M", new[] { "XX", "X/Y", "IX/YI" }, (l, m) => m % 2 == 0),
                ("XX across, XY down, XX on the anti-diagonal breaks on every torus", new[] { "XX", "X/Y", "IX/XI" }, (l, m) => false),
                ("the unit mask 1 + x + y with an odd sign, XY/XI, alone pairs on every torus", new[] { "XY/XI" }, (l, m) => true),
                ("the unit mask 1 + x + y beside XY across breaks on every torus", new[] { "XX/XI", "XY" }, (l, m) => false),
            };
            int readings = 0, pairing = 0, bad = 0;
            string? first = null;
            foreach (var (name, ts, expect) in examples)
            {
                foreach (var (l, m) in small)
                {
                    bool graph = PairsOnTorus(l, m, ts);
                    readings++;
                    if (graph) pairing++;
                    if (graph != expect(l, m)) { bad++; first ??= $"{name}: Theorem 2 on {l} × {m}"; }
                }
                foreach (var (l, m) in large)
                {
                    bool law = !SignedMaskBalanceClaim.TorusBreaks(ts, l, m);
                    readings++;
                    if (law) pairing++;
                    if (law != expect(l, m) || law == BalanceBreaks(ts, l, m)) { bad++; first ??= $"{name}: the law on {l} × {m}"; }
                }
            }
            families.Add(new FamilyReading("the worked examples (Theorem 2 to 4 × 4, the law and the balance to 64 sites)", examples.Length, readings, pairing, bad, first, "the page"));
        }

        Families = families;
        OneNumberMismatches = oneNumber;
        OneNumberReadings = oneNumberReadings;
        FoldControlReadings = foldReadings;
        SwappedFoldMismatches = swapped;
        ParityFoldMismatches = parity;
    }

    public string DisplayName => "TorusMaskBalanceWitness (F173 Theorem F)";

    public string Summary =>
        $"Theorem F on {Families.Count} families, {Families.Sum(f => f.Readings)} readings, {Families.Sum(f => f.Mismatches)} disagreeing; " +
        $"the controls disagree with the balance: one number per mask {OneNumberMismatches} of {OneNumberReadings} (T₁, T₂ ≥ 2), the fold with " +
        $"T₁ and T₂ swapped {SwappedFoldMismatches} and the fold onto the parities {ParityFoldMismatches} of {FoldControlReadings}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var f in Families)
                yield return new InspectableNode(
                    displayName: f.Name,
                    summary: $"{f.Cells} {(f.Cells == 1 ? "cell" : "cells")}, {f.Readings} readings ({f.Pairing} pairing), " +
                             (f.Agrees ? $"every one agreeing with {f.Against}" : $"{f.Mismatches} DISAGREEING, first {f.FirstMismatch}"),
                    provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the controls, which must disagree with the balance",
                summary: $"Theorem D's one-liner with one number per mask and T = min(T₁, T₂), on the tori with T₁, T₂ ≥ 2: {OneNumberMismatches} of " +
                         $"{OneNumberReadings} readings; the fold with T₁ and T₂ swapped: {SwappedFoldMismatches} of {FoldControlReadings}; the fold onto " +
                         $"the parities of L and M: {ParityFoldMismatches} of {FoldControlReadings}",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
