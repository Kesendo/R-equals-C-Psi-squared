using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live battery for Theorems F and G of F173 (<c>docs/proofs/PROOF_SIGNED_MASK_BALANCE.md</c>, typed
/// <see cref="SignedMaskBalanceClaim.TorusBreaks"/> and <see cref="SignedMaskBalanceClaim.TorusBreaksByDifferences"/>):
/// templates placed at every translate of the L × M torus, one string per flip set, Z dephasing on every site, read on the
/// 2-part T₁ × T₂ of the translations (the second form, a rearrangement of the first, computed beside it and no second
/// route). Held against Theorem 2 of the complement connection walked on the bitstrings
/// (<see cref="ComplementConnectionGraph"/>) on tori of up to 16 sites, 4 × 4 the one among them whose 2-part is not
/// cyclic, and against the balance on the whole torus (<see cref="SignedMaskBalanceClaim.Balanced"/>, Theorem A's
/// elimination on the unfolded placements) to 64 sites: the templates of the 2 × 2 window, the ring as the L × 1 torus
/// against Theorem D, the tori of odd L against Theorem D's one-liner read on the masks summed along x, and the worked
/// examples of the page. The staircase beside it: the corners read on the small tori (<see cref="SignedMaskBalanceClaim.Corners"/>)
/// rearranged against the balance, one template's corners against the exponents of its mask coarsened to bit lengths (<see cref="SignedMaskBalanceClaim.SingleTemplateCorners"/>),
/// two templates by divisibility in A (<see cref="SignedMaskBalanceClaim.TwoTemplatesBreakByDivisibility"/>) against Theorem F,
/// the corner bound (<see cref="SignedMaskBalanceClaim.CornerBound"/>), the honeycomb and the kagome lattice
/// (<see cref="SignedMaskBalanceClaim.LatticeBreaks"/>) and the size bound (<see cref="SignedMaskBalanceClaim.SufficientTorus"/>). Theorem G on the open
/// patch (<see cref="SignedMaskBalanceClaim.PatchBreaks"/>) against Theorem 2 on patches of 4 to 12 sites and against the balance on every patch to
/// 8 × 8, its row M = 1 against Theorem B, three templates' breaking patches closed upwards, and its worked example; its two
/// controls, the masks as their own quotients and the threshold along x alone, must disagree. Six controls of Theorem F must disagree:
/// Theorem D's one-liner carried over with one number per mask (its least order in the powers of 1 + x and 1 + y) and
/// T = min(T₁, T₂), read on the tori with T₁, T₂ ≥ 2, where A is no chain, the fold with T₁ and T₂ swapped, the fold onto the
/// parities of L and M in place of their 2-parts, the staircase read with one corner, each direction on its own, the corner
/// bound one doubling tighter, and the lattice fold with the sites of a cell merged.</summary>
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
    /// <summary>The control of the staircase: each direction read on its own, one corner at the least e and the least f of
    /// the corners, against the balance on the readings of the staircase family.</summary>
    public int OneCornerMismatches { get; }
    public int StaircaseReadings { get; }
    /// <summary>The template sets of the staircase line whose staircase has more than one corner.</summary>
    public int CornerSetsWithSeveral { get; }
    /// <summary>The one-corner control's disagreements on 3 × 3, the one torus of the sixteen with both sides odd.</summary>
    public int OneCornerMismatchesOnThreeByThree { get; }
    public bool AllAgree => Families.All(f => f.Agrees);
    /// <summary>The corner bound one doubling tighter, 2^eⱼ &lt; n_u and 2^fⱼ &lt; n_v, against the read corners of the bound
    /// lines: the corners that break it.</summary>
    public int TighterBoundMismatches { get; }
    /// <summary>The lattice fold with the Q sites of a cell merged into one, against the balance on the whole lattice.</summary>
    public int MergedSitesMismatches { get; }
    public int LatticeReadings { get; }
    /// <summary>The lattice readings whose merged set the setting accepts, where the merged control is read.</summary>
    public int MergedSitesReadings { get; }
    /// <summary>Theorem G's controls: the quotients taken as the masks themselves (no gcd), and the threshold read along x
    /// alone (M only as wide as the templates), against the balance on the readings of its balance family.</summary>
    public int PatchNoGcdMismatches { get; }
    public int PatchOneDirectionMismatches { get; }
    public int PatchControlReadings { get; }
    /// <summary>The pairs of Theorem G's balance family with c = 1, and those whose corner the gcd moves (or whose c it
    /// changes).</summary>
    public int PatchBreakingPairs { get; }
    public int PatchPairsWhereTheGcdMatters { get; }
    /// <summary>The three-template sets whose breaking patches up to 8 × 8 have more than one corner: no single quadrant.</summary>
    public int PatchTriplesWithSeveralCorners { get; }
    public bool ControlsBreak => OneNumberMismatches > 0 && SwappedFoldMismatches > 0 && ParityFoldMismatches > 0 && OneCornerMismatches > 0
        && TighterBoundMismatches > 0 && MergedSitesMismatches > 0 && PatchNoGcdMismatches > 0 && PatchOneDirectionMismatches > 0;

    // the corner bound against a staircase read on the small tori (2^(e + f) ≤ 64): every read corner inside
    // 2^(eⱼ − 1) < n_u, 2^(fⱼ − 1) < n_v, and the verdict unchanged by a doubling of T₂ from T₂ ≥ n_v (of T₁ from T₁ ≥ n_u)
    // wherever both points are read; returns the checks made, the failures, and the corners outside the bound one doubling tighter
    static (int Checks, int Failures, int Tighter) CheckBound(IReadOnlyDictionary<(int, int), bool> grid, IReadOnlyList<(int E, int F)> corners, int nu, int nv)
    {
        int checks = 0, failures = 0, tighter = 0;
        foreach (var (e, f) in corners)
        {
            checks++;
            if ((e >= 1 && (1 << (e - 1)) >= nu) || (f >= 1 && (1 << (f - 1)) >= nv)) failures++;
            if ((1 << e) >= nu && e >= 1 || (1 << f) >= nv && f >= 1) tighter++;
        }
        foreach (var ((e, f), pairs) in grid)
        {
            if ((1 << f) >= nv && grid.TryGetValue((e, f + 1), out var up)) { checks++; if (up != pairs) failures++; }
            if ((1 << e) >= nu && grid.TryGetValue((e + 1, f), out var right)) { checks++; if (right != pairs) failures++; }
        }
        return (checks, failures, tighter);
    }

    /// <summary>The honeycomb (Q = 2) and kagome (Q = 3) templates of the lattice line: two-site bonds between the sites of a
    /// cell and of its neighbours, with the letters XX, XY, YX, YY, XZ or ZY.</summary>
    public static IReadOnlyList<(int Q, SignedMaskBalanceClaim.LatticeLetter[] Template, string Name)> LatticeTemplates { get; } = BuildLatticeTemplates();

    static IReadOnlyList<(int, SignedMaskBalanceClaim.LatticeLetter[], string)> BuildLatticeTemplates()
    {
        var list = new List<(int, SignedMaskBalanceClaim.LatticeLetter[], string)>();
        // honeycomb: A = 0, B = 1; nearest bonds A(0,0)-B(0,0), A(1,0)-B(0,0), A(0,1)-B(0,0), second neighbours A-A along x and B-B along y
        var honeycomb = new (int I1, int J1, int R1, int I2, int J2, int R2)[]
            { (0, 0, 0, 0, 0, 1), (1, 0, 0, 0, 0, 1), (0, 1, 0, 0, 0, 1), (0, 0, 0, 1, 0, 0), (0, 0, 1, 0, 1, 1) };
        // kagome: A = 0, B = 1, C = 2 (A at the cell's corner, B halfway along x, C halfway along y); the up triangle A-B, A-C, B-C
        // in a cell, and the down triangles from the bonds B(0,0)-A(1,0), C(0,0)-A(0,1), B(0,1)-C(1,0), which close
        // {A(1,0), B(0,0), C(1,-1)}: every site in two triangles
        var kagome = new (int I1, int J1, int R1, int I2, int J2, int R2)[]
            { (0, 0, 0, 0, 0, 1), (0, 0, 0, 0, 0, 2), (0, 0, 1, 0, 0, 2), (1, 0, 0, 0, 0, 1), (0, 1, 0, 0, 0, 2), (0, 1, 1, 1, 0, 2) };
        foreach (var (q, bonds) in new[] { (2, honeycomb), (3, kagome) })
            foreach (var bd in bonds)
                foreach (var pair in new[] { "XX", "XY", "YX", "YY", "XZ", "ZY" })
                    list.Add((q, new[] { new SignedMaskBalanceClaim.LatticeLetter(bd.I1, bd.J1, bd.R1, pair[0]), new SignedMaskBalanceClaim.LatticeLetter(bd.I2, bd.J2, bd.R2, pair[1]) },
                        $"{pair[0]}({bd.I1},{bd.J1},{bd.R1}){pair[1]}({bd.I2},{bd.J2},{bd.R2})"));
        return list;
    }

    /// <summary>Theorem 2 on the L × M lattice with Q sites in a cell (site (i, j, r) at position r + Q·(i + L·j)), each
    /// template at every translate of the cells, coefficient 1.</summary>
    public static bool PairsOnLattice(int l, int m, int q, IReadOnlyList<SignedMaskBalanceClaim.LatticeLetter[]> templates)
    {
        int n = l * m * q;
        var placed = new List<(PauliMask, BigInteger)>();
        foreach (var t in templates)
            for (int b = 0; b < m; b++)
                for (int a = 0; a < l; a++)
                {
                    var s = Enumerable.Repeat('I', n).ToArray();
                    foreach (var x in t) s[x.R + q * (((x.I + a) % l + l) % l + l * (((x.J + b) % m + m) % m))] = x.Letter;
                    placed.Add((PauliMask.Parse(new string(s)), BigInteger.One));
                }
        var (near, far, _) = ComplementConnectionGraph.EverySiteDephased(n, ComplementConnectionGraph.Columns(n, placed));
        return near == far;
    }

    static bool LatticePlaceable(IReadOnlyList<SignedMaskBalanceClaim.LatticeLetter[]> ts, int l, int m, int q)
    {
        try { SignedMaskBalanceClaim.OnLattice(ts, l, m, q); return true; }
        catch (ArgumentException) { return false; }
    }

    /// <summary>One template per mask of the 3 × 3 window touching its first column and first row (400 masks), every
    /// flipped letter X, and with <paramref name="odd"/> the first flipped letter (row by row) a Y; the same mask order
    /// either way.</summary>
    public static IReadOnlyList<string> ThreeByThreeTemplates(bool odd)
    {
        var list = new List<string>();
        for (int bits = 1; bits < 512; bits++)
        {
            bool row0 = (bits & 0b111) != 0, col0 = (bits & 0b001_001_001) != 0;
            if (!row0 || !col0) continue;
            var rows = new string[3];
            bool placedY = false;
            for (int j = 0; j < 3; j++)
            {
                var r = new char[3];
                for (int i = 0; i < 3; i++)
                {
                    bool on = ((bits >> (i + 3 * j)) & 1) == 1;
                    r[i] = !on ? 'I' : (odd && !placedY ? 'Y' : 'X');
                    if (on && odd) placedY = true;
                }
                rows[j] = new string(r);
            }
            list.Add(string.Join("/", rows));
        }
        return list;
    }

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

    /// <summary>Theorem 2 on the open L × M patch (site (i, j) at position i + L·j), each template at every translate that fits,
    /// coefficient 1.</summary>
    public static bool PairsOnPatch(int l, int m, params string[] templates)
    {
        var placed = new List<(PauliMask, BigInteger)>();
        foreach (var t in templates)
        {
            var letters = SignedMaskBalanceClaim.TorusLetters(t);
            var (wx, wy) = SignedMaskBalanceClaim.Extent(t);
            for (int b = 0; b + wy <= m; b++)
                for (int a = 0; a + wx <= l; a++)
                {
                    var s = Enumerable.Repeat('I', l * m).ToArray();
                    foreach (var (i, j, c) in letters) s[i + a + l * (j + b)] = c;
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
        int oneNumber = 0, oneNumberReadings = 0, foldReadings = 0, swapped = 0, parity = 0, oneCorner = 0, staircaseReadings = 0, severalCorners = 0, oneCornerThree = 0;
        int tighter = 0, merged = 0, mergedReadings = 0, latticeReadings = 0;
        int patchNoGcd = 0, patchOneDirection = 0, patchControlReadings = 0, patchBreaking = 0, patchGcdMatters = 0, tripleSeveralCorners = 0;
        int boundSets = 0, boundCommon = 0, boundChecks = 0, boundBad = 0;
        string? boundFirst = null;

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

        // the staircase: on every template set of the balance line, the pairing points (e, f) of the small tori (2^(e + f) ≤ 64)
        // closed upwards, and the corners read there against the balance on the sixteen tori; the control reads each direction
        // on its own, one corner at the least e and the least f of the corners
        {
            int cells = 0, readings = 0, pairing = 0, bad = 0, notClosed = 0, several = 0;
            string? first = null;
            foreach (var ts in singles.Concat(pairs).Concat(triples))
            {
                var grid = new Dictionary<(int, int), bool>();
                for (int a = 0; a <= 6; a++)
                    for (int b = 0; a + b <= 6; b++) grid[(a, b)] = !SignedMaskBalanceClaim.BreaksAtTwoParts(ts, 1 << a, 1 << b);
                foreach (var ((a, b), p) in grid)
                    if (p && ((grid.TryGetValue((a + 1, b), out var r) && !r) || (grid.TryGetValue((a, b + 1), out var u) && !u))) notClosed++;
                var corners = SignedMaskBalanceClaim.Corners(ts);
                if (corners.Count > 1) several++;
                cells++;
                // masks of the 2 × 2 window have degree at most 2 in u, v, so a colength without a common factor is at most 4 < 8
                if (SignedMaskBalanceClaim.CornerBound(ts, 8) is var (_, nu, nv))
                {
                    boundSets++;
                    var (checks, failures, tight) = CheckBound(grid, corners, nu, nv);
                    boundChecks += checks; tighter += tight;
                    if (failures > 0) { boundBad += failures; boundFirst ??= string.Join(" + ", ts); }
                }
                else boundCommon++;
                foreach (var (l, m) in Tori)
                {
                    if (!Placeable(ts, l, m)) continue;
                    bool balance = !BalanceBreaks(ts, l, m);
                    bool byCorners = corners.Any(c => l % (1 << c.E) == 0 && m % (1 << c.F) == 0);
                    readings++;
                    if (balance) pairing++;
                    if (balance != byCorners) { bad++; first ??= $"{string.Join(" + ", ts)} on {l} × {m}"; }
                    staircaseReadings++;
                    if (corners.Count > 0)
                    {
                        int e0 = corners.Min(c => c.E), f0 = corners.Min(c => c.F);
                        if (balance != (l % (1 << e0) == 0 && m % (1 << f0) == 0)) { oneCorner++; if (l == 3 && m == 3) oneCornerThree++; }
                    }
                    else if (balance) oneCorner++;
                }
            }
            if (notClosed > 0) { bad += notClosed; first ??= $"{notClosed} pairing points with a breaking point above them"; }
            severalCorners = several;
            families.Add(new FamilyReading("the staircase: on every template set of the balance line, the pairing points of the small tori closed upwards and the corners read there giving back the balance (a rearrangement, not a second route)", cells, readings, pairing, bad, first, "the balance"));
            families.Add(new FamilyReading($"the corner bound on the template sets of the balance line whose read colength at n = 8 stays below 8 ({boundCommon} reaching it left out): the read corners inside 2^(e − 1) < n_u, 2^(f − 1) < n_v, the verdict unchanged by a doubling past n_u, n_v", boundSets, boundChecks, 0, boundBad, boundFirst, "the colength"));
        }

        // the corner bound on random sets of two to four masks of the 3 × 3 window, random signs, read at n = 16
        {
            var pool3 = ThreeByThreeTemplates(odd: true).Zip(ThreeByThreeTemplates(odd: false)).ToArray();
            int cells = 0, checks = 0, bad = 0, common = 0;
            string? first = null;
            for (int k = 0; k < 600; k++)
            {
                int size = 2 + rng.Next(3);
                var ts = Enumerable.Range(0, pool3.Length).OrderBy(_ => rng.Next()).Take(size)
                    .Select(i => rng.Next(2) == 0 ? pool3[i].First : pool3[i].Second).ToArray();
                if (SignedMaskBalanceClaim.CornerBound(ts, 16) is not var (_, nu, nv)) { common++; continue; }
                var grid = new Dictionary<(int, int), bool>();
                for (int e = 0; e <= 6; e++)
                    for (int f = 0; e + f <= 6; f++) grid[(e, f)] = !SignedMaskBalanceClaim.BreaksAtTwoParts(ts, 1 << e, 1 << f);
                var (c, failures, tight) = CheckBound(grid, SignedMaskBalanceClaim.Corners(ts), nu, nv);
                cells++; checks += c; tighter += tight;
                if (failures > 0) { bad += failures; first ??= string.Join(" + ", ts); }
            }
            families.Add(new FamilyReading($"the corner bound on 600 random sets of two to four masks of the 3 × 3 window, read at n = 16 ({common} whose read colength reaches 16 left out)", cells, checks, 0, bad, first, "the colength"));
        }

        // beyond the square lattice: the honeycomb (Q = 2) and the kagome lattice (Q = 3), random sets of one to three bonds of
        // one lattice, the fold onto T₁ × T₂ cells of Q sites against the balance on the whole lattice, and against Theorem 2 on
        // twelve sites; the control merges the Q sites of a cell into one
        {
            int cells = 0, readings = 0, pairing = 0, bad = 0;
            string? first = null;
            var lattices = new (int L, int M)[] { (3, 3), (3, 4), (4, 3), (4, 4), (5, 4), (6, 3), (3, 6), (6, 5), (4, 8), (8, 4), (5, 5), (10, 3) };
            foreach (int q in new[] { 2, 3 })
            {
                var bonds = LatticeTemplates.Where(x => x.Q == q).ToArray();
                for (int k = 0; k < 300; k++)
                {
                    var drawn = bonds.OrderBy(_ => rng.Next()).Take(1 + rng.Next(3)).ToArray();
                    var ts = drawn.Select(x => x.Template).ToArray();
                    string name = string.Join(" + ", drawn.Select(x => x.Name));
                    bool any = false;
                    foreach (var (l, m) in lattices)
                    {
                        if (l * m * q > 64 || !LatticePlaceable(ts, l, m, q)) continue;
                        any = true;
                        bool balance = SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnLattice(ts, l, m, q));
                        bool law = !SignedMaskBalanceClaim.LatticeBreaks(ts, l, m, q);
                        readings++; latticeReadings++;
                        if (balance) pairing++;
                        if (balance != law) { bad++; first ??= $"{name} on {l} × {m} × {q}"; }
                        // the control: every site of a cell read as the one site of a square-lattice torus
                        var mergedTs = ts.Select(t => t.Select(x => x with { R = 0 }).ToArray()).ToArray();
                        // (a merged set the setting refuses, two letters on one site or a flip set placed twice, is not read)
                        if (LatticePlaceable(mergedTs, l, m, 1))
                        {
                            mergedReadings++;
                            if (SignedMaskBalanceClaim.LatticeBreaks(mergedTs, l, m, 1) == balance) merged++;
                        }
                    }
                    if (any) cells++;
                    if (k < 40)
                        foreach (var (l, m) in q == 2 ? new[] { (3, 2), (2, 3) } : new[] { (2, 2) })
                        {
                            if (!LatticePlaceable(ts, l, m, q)) continue;
                            bool graph = PairsOnLattice(l, m, q, ts);
                            bool law = !SignedMaskBalanceClaim.LatticeBreaks(ts, l, m, q);
                            readings++;
                            if (graph) pairing++;
                            if (graph != law) { bad++; first ??= $"{name} on {l} × {m} × {q}, Theorem 2"; }
                        }
                }
            }
            // the honeycomb is bipartite and the kagome lattice is made of triangles: XY on every nearest bond
            var hc = LatticeTemplates.Where(x => x.Q == 2 && x.Name.StartsWith("X(") && x.Name.Contains(")Y(")).Take(3).Select(x => x.Template).ToArray();
            var kg = LatticeTemplates.Where(x => x.Q == 3 && x.Name.StartsWith("X(") && x.Name.Contains(")Y(")).Select(x => x.Template).ToArray();
            foreach (var (l, m) in lattices)
            {
                if (l * m * 2 <= 64 && LatticePlaceable(hc, l, m, 2))
                { readings++; if (SignedMaskBalanceClaim.LatticeBreaks(hc, l, m, 2) || !SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnLattice(hc, l, m, 2))) { bad++; first ??= $"the honeycomb with XY on its three bonds on {l} × {m}"; } else pairing++; }
                if (l * m * 3 <= 64 && LatticePlaceable(kg, l, m, 3))
                { readings++; if (!SignedMaskBalanceClaim.LatticeBreaks(kg, l, m, 3) || SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnLattice(kg, l, m, 3))) { bad++; first ??= $"the kagome lattice with XY on its six bonds on {l} × {m}"; } }
            }
            families.Add(new FamilyReading("the honeycomb and the kagome lattice, 300 random sets of one to three bonds each: the fold onto T₁ × T₂ cells of Q sites against the balance on the whole lattice, on the listed lattices to 64 sites (twelve for the honeycomb, seven for the kagome lattice), and Theorem 2 on 3 × 2 and 2 × 3 honeycomb cells and 2 × 2 kagome cells (twelve sites, the first 40 sets of each); XY on every nearest bond pairs on every honeycomb and breaks on every kagome lattice", cells, readings, pairing, bad, first, "the balance and Theorem 2"));
        }

        // one template of the 3 × 3 window, odd sign: the corners read on the small tori against the exponents of the mask, coarsened
        {
            int readings = 0, pairing = 0, bad = 0;
            string? first = null;
            foreach (var t in ThreeByThreeTemplates(odd: true))
            {
                var read = SignedMaskBalanceClaim.Corners(new[] { t });
                var star = SignedMaskBalanceClaim.SingleTemplateCorners(t);
                readings++;
                if (read.Count == 1 && read[0] == (0, 0)) pairing++;
                if (!read.SequenceEqual(star)) { bad++; first ??= t; }
            }
            families.Add(new FamilyReading("one template of the 3 × 3 window with an odd sign: its corners against the exponents of its mask in u and v, coarsened to bit lengths, and of those the minimal ones (pairing: the corner (0, 0), every torus)", readings, readings, pairing, bad, first, "the exponents"));
        }

        // two templates of the 3 × 3 window, random signs: the divisibility law against Theorem F on every small torus to 64 cells
        {
            var pool3 = ThreeByThreeTemplates(odd: true).Zip(ThreeByThreeTemplates(odd: false)).ToArray();
            int cells = 0, readings = 0, pairing = 0, bad = 0;
            string? first = null;
            for (int k = 0; k < 600; k++)
            {
                int i = rng.Next(pool3.Length), j = rng.Next(pool3.Length - 1);
                if (j >= i) j++;
                string t1 = rng.Next(2) == 0 ? pool3[i].First : pool3[i].Second, t2 = rng.Next(2) == 0 ? pool3[j].First : pool3[j].Second;
                cells++;
                for (int a = 0; a <= 6; a++)
                    for (int b = 0; a + b <= 6; b++)
                    {
                        bool law = !SignedMaskBalanceClaim.BreaksAtTwoParts(new[] { t1, t2 }, 1 << a, 1 << b);
                        bool division = !SignedMaskBalanceClaim.TwoTemplatesBreakByDivisibility(t1, t2, 1 << a, 1 << b);
                        readings++;
                        if (law) pairing++;
                        if (law != division) { bad++; first ??= $"{t1} + {t2} on {1 << a} × {1 << b}"; }
                    }
            }
            families.Add(new FamilyReading("two templates of the 3 × 3 window, random masks and signs: the divisibility law on every small torus to 64 cells", cells, readings, pairing, bad, first, "Theorem F"));
        }

        // the size bound: every single template and every pair of different masks of the 3 × 3 window placed on
        // L = 2w_x − 1, M = 2w_y − 1 without two placements coinciding
        {
            var pool3 = ThreeByThreeTemplates(odd: false);
            int readings = 0, pairing = 0, bad = 0;
            string? first = null;
            for (int i = 0; i < pool3.Count; i++)
                for (int j = i; j < pool3.Count; j++)
                {
                    var ts = i == j ? new[] { pool3[i] } : new[] { pool3[i], pool3[j] };
                    var (l, m) = SignedMaskBalanceClaim.SufficientTorus(ts);
                    readings++;
                    if (!Placeable(ts, l, m)) { bad++; first ??= $"{string.Join(" + ", ts)} on {l} × {m}"; }
                    else if (!SignedMaskBalanceClaim.TorusBreaks(ts, l, m)) pairing++;
                    else { bad++; first ??= $"{string.Join(" + ", ts)} breaks on {l} × {m} with every sign even"; }
                }
            if (Placeable(new[] { "XIX" }, 4, 1)) { bad++; first ??= "XIX on the ring of 4 (2w − 2) placed without coinciding"; }
            families.Add(new FamilyReading("the size bound L ≥ 2w_x − 1, M ≥ 2w_y − 1: every template and pair of different masks of the 3 × 3 window placed there (all signs even, so every one pairs)", readings, readings, pairing, bad, first, "OnTorus"));
        }

        // the worked examples of the page: Theorem 2 to 4 × 4, the law on the tori to 64 sites with L, M ≥ 3
        {
            var small = new[] { (3, 3), (3, 4), (4, 3), (4, 4) };
            var large = (from l in Enumerable.Range(3, 19) from m in Enumerable.Range(3, 19) where l * m <= 64 select (l, m)).ToArray();
            var examples = new (string Name, string[] Ts, Func<int, int, bool> Pairs)[]
            {
                ("XY on the anti-diagonal alone pairs exactly where L or M is even (two corners)", new[] { "IX/YI" }, (l, m) => l % 2 == 0 || m % 2 == 0),
                ("u² + uv + v², odd sign, pairs exactly at 4 | L, 4 | M, or 2 | L and 2 | M (three corners)", new[] { "XXX/XX/Y" }, (l, m) => l % 4 == 0 || m % 4 == 0 || (l % 2 == 0 && m % 2 == 0)),
                ("the parabola, XX on (1, 0), (0, 2) beside XY on (0, 0), (1, 2), pairs exactly at 4 | L and 8 | M", new[] { "IX//X", "X//IY" }, (l, m) => l % 4 == 0 && m % 8 == 0),
                ("u + v³ beside u + v³ + u³ pairs exactly at 4 | L and 16 | M", new[] { "IX/X/X/X", "XIXX/X/X/Y" }, (l, m) => l % 4 == 0 && m % 16 == 0),
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
                    if (!Placeable(ts, l, m)) continue;
                    bool graph = PairsOnTorus(l, m, ts);
                    readings++;
                    if (graph) pairing++;
                    if (graph != expect(l, m)) { bad++; first ??= $"{name}: Theorem 2 on {l} × {m}"; }
                }
                foreach (var (l, m) in large)
                {
                    if (!Placeable(ts, l, m)) continue;
                    bool law = !SignedMaskBalanceClaim.TorusBreaks(ts, l, m);
                    readings++;
                    if (law) pairing++;
                    if (law != expect(l, m) || law == BalanceBreaks(ts, l, m)) { bad++; first ??= $"{name}: the law on {l} × {m}"; }
                }
            }
            // each example's corners, read on the small tori, as the page names them
            var named = new (string[] Ts, (int, int)[] Corners)[]
            {
                (new[] { "IX/YI" }, new[] { (0, 1), (1, 0) }),
                (new[] { "XXX/XX/Y" }, new[] { (0, 2), (1, 1), (2, 0) }),
                (new[] { "IX//X", "X//IY" }, new[] { (2, 3) }),
                (new[] { "IX/X/X/X", "XIXX/X/X/Y" }, new[] { (2, 4) }),
                (new[] { "XIY" }, new[] { (2, 0) }),
                (new[] { "XXXY" }, new[] { (2, 0) }),
            };
            foreach (var (ts, expected) in named)
            {
                readings++;
                if (!SignedMaskBalanceClaim.Corners(ts).SequenceEqual(expected)) { bad++; first ??= $"the corners of {string.Join(" + ", ts)}"; }
            }
            families.Add(new FamilyReading("the worked examples (Theorem 2 to 4 × 4, the law and the balance to 64 sites, six sets of corners)", examples.Length, readings, pairing, bad, first, "the page"));
        }

        // Theorem G, the open patch: two templates, read through the gcd in GF(2)[x, y]; its own seeds, one per draw, so the
        // families above draw what they drew before and no family's draw moves another's
        {
            var prng = new Random(seed + 1);
            var patchPairs = new List<string[]>(pairs);
            var pool3 = ThreeByThreeTemplates(odd: true).Concat(ThreeByThreeTemplates(odd: false)).ToArray();
            static bool Translates(string a, string b)
            {
                try { SignedMaskBalanceClaim.OnPatch(new[] { a, b }, 8, 8); return false; }
                catch (ArgumentException) { return true; }
            }
            var deepPairs = new List<string[]>();
            var drawnKeys = new HashSet<string>();
            while (deepPairs.Count < 1500)
            {
                string a = pool3[prng.Next(pool3.Length)], b = pool3[prng.Next(pool3.Length)];
                if (!Translates(a, b) && drawnKeys.Add(string.CompareOrdinal(a, b) < 0 ? a + "|" + b : b + "|" + a)) deepPairs.Add(new[] { a, b });
            }

            // against Theorem 2 on the bitstrings: 150 pairs of the 2 × 2 window and 50 of the 3 × 3 window, every patch of 4 to 12 sites
            // with both sides at least 2
            {
                var patches = (from l in Enumerable.Range(2, 7) from m in Enumerable.Range(2, 7) where l * m <= 12 select (l, m)).ToArray();
                var drawRng = new Random(seed + 2);
                var drawn = patchPairs.OrderBy(_ => drawRng.Next()).Take(150).Concat(deepPairs.Take(50)).ToList();
                var jobs = (from ts in drawn from p in patches select (Ts: ts, p.l, p.m)).ToList();
                var graphs = new bool[jobs.Count];
                Parallel.For(0, jobs.Count, k => graphs[k] = PairsOnPatch(jobs[k].l, jobs[k].m, jobs[k].Ts));
                int readings = 0, pairing = 0, bad = 0;
                string? first = null;
                for (int k = 0; k < jobs.Count; k++)
                {
                    var (ts, l, m) = jobs[k];
                    bool law = !SignedMaskBalanceClaim.PatchBreaks(ts[0], ts[1], l, m);
                    readings++;
                    if (graphs[k]) pairing++;
                    if (graphs[k] != law) { bad++; first ??= $"{ts[0]} + {ts[1]} on {l} × {m}"; }
                }
                families.Add(new FamilyReading("Theorem G: two templates on the open patch, 150 pairs of the 2 × 2 window and 50 of the 3 × 3 window, every patch of 4 to 12 sites with both sides at least 2", drawn.Count, readings, pairing, bad, first, "Theorem 2"));
            }

            // against the balance: every pair of different masks of the 2 × 2 window and 1500 pairs of the 3 × 3 window (seeded), every
            // patch up to 8 × 8; the controls read the same patches
            {
                var patches = (from l in Enumerable.Range(1, 8) from m in Enumerable.Range(1, 8) select (l, m)).ToArray();
                int readings = 0, pairing = 0, bad = 0;
                string? first = null;
                foreach (var ts in patchPairs.Concat(deepPairs))
                {
                    var (_, corner) = SignedMaskBalanceClaim.Patch(ts[0], ts[1]);
                    var (_, noGcd) = SignedMaskBalanceClaim.Patch(ts[0], ts[1], divideByGcd: false);
                    var (wy1, wy2) = (SignedMaskBalanceClaim.Extent(ts[0]).Wy, SignedMaskBalanceClaim.Extent(ts[1]).Wy);
                    if (corner is not null) patchBreaking++;
                    if (!Equals(corner, noGcd)) patchGcdMatters++;
                    foreach (var (l, m) in patches)
                    {
                        bool balance = SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnPatch(ts, l, m));
                        bool law = !SignedMaskBalanceClaim.PatchBreaks(ts[0], ts[1], l, m);
                        readings++;
                        if (balance) pairing++;
                        if (balance != law) { bad++; first ??= $"{ts[0]} + {ts[1]} on {l} × {m}"; }
                        if (balance == (noGcd is { } n0 && l >= n0.L && m >= n0.M)) patchNoGcd++;
                        if (balance == (corner is { } c1 && l >= c1.L && m >= Math.Max(wy1, wy2))) patchOneDirection++;
                    }
                }
                patchControlReadings = readings;
                families.Add(new FamilyReading("Theorem G: every pair of different masks of the 2 × 2 window and 1500 different pairs of the 3 × 3 window, every patch up to 8 × 8", patchPairs.Count + deepPairs.Count, readings, pairing, bad, first, "the balance"));
            }

            // masks off the axes: 400 pairs of the 2 × 2 window, each template moved by a leading identity column, row, both or
            // neither (its mask then divisible by x, y or xy, the extents counting the identity letters), every patch up to 8 × 8
            {
                var offRng = new Random(seed + 3);
                static string Move(string t, int how) => string.Join("/", ((how & 2) != 0 ? new[] { "" } : Array.Empty<string>())
                    .Concat(t.Split('/').Select(r => (how & 1) != 0 ? "I" + r : r)));
                var cellsRead = new List<string[]>();
                foreach (var ts in patchPairs.OrderBy(_ => offRng.Next()).Take(400))
                    cellsRead.Add(new[] { Move(ts[0], offRng.Next(4)), Move(ts[1], offRng.Next(4)) });
                int readings = 0, pairing = 0, bad = 0;
                string? first = null;
                foreach (var ts in cellsRead)
                {
                    for (int l = 1; l <= 8; l++)
                        for (int m = 1; m <= 8; m++)
                        {
                            bool balance = SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnPatch(ts, l, m));
                            bool law = !SignedMaskBalanceClaim.PatchBreaks(ts[0], ts[1], l, m);
                            readings++;
                            if (balance) pairing++;
                            if (balance != law) { bad++; first ??= $"{ts[0]} + {ts[1]} on {l} × {m}"; }
                        }
                }
                families.Add(new FamilyReading("Theorem G off the axes, 400 pairs of the 2 × 2 window moved by identity rows and columns, every patch up to 8 × 8", cellsRead.Count, readings, pairing, bad, first, "the balance"));
            }

            // the gcd where its pseudo-remainders pass degree 63 in x: 200 planted factors, masks g·u and g·v with g of two rows
            // and two columns and u, v drawn in a 6 × 6 window with a letter in each row's first column (196 of 200 passing 63 in
            // a scout); the gcd read must be a multiple of g, and it divides both masks or Patch would throw
            {
                var plantRng = new Random(seed + 5);
                static ulong Mul(ulong a, ulong b) { ulong r = 0; for (; b != 0; b >>= 1, a <<= 1) if ((b & 1) == 1) r ^= a; return r; }
                static ulong[] Times(ulong[] a, ulong[] b)
                {
                    var r = new ulong[a.Length + b.Length - 1];
                    for (int i = 0; i < a.Length; i++)
                        for (int j = 0; j < b.Length; j++) r[i + j] ^= Mul(a[i], b[j]);
                    return r;
                }
                static string Template(ulong[] rows, bool odd)
                {
                    bool placedY = false;
                    return string.Join("/", rows.Select(r => new string(Enumerable.Range(0, r == 0 ? 0 : 64 - BitOperations.LeadingZeroCount(r))
                        .Select(i => ((r >> i) & 1) == 0 ? 'I' : odd && !placedY && (placedY = true) ? 'Y' : 'X').ToArray())));
                }
                // exact divisibility of f by g in GF(2)[x][y], rows the coefficients of y^j
                static bool Divides(ulong[] g, IReadOnlyList<ulong> f)
                {
                    var rem = f.ToArray();
                    int dg = 63 - BitOperations.LeadingZeroCount(g[^1]);
                    for (int d = rem.Length - 1; d >= g.Length - 1; d--)
                    {
                        while (rem[d] != 0)
                        {
                            int dr = 63 - BitOperations.LeadingZeroCount(rem[d]);
                            if (dr < dg) return false;
                            for (int j = 0; j < g.Length; j++) rem[d - g.Length + 1 + j] ^= g[j] << (dr - dg);
                        }
                    }
                    return rem.All(r => r == 0);
                }
                int cells = 0, pairing = 0, bad = 0;
                string? first = null;
                while (cells < 200)
                {
                    var g = new[] { (ulong)plantRng.Next(1, 4), (ulong)plantRng.Next(1, 4) };
                    ulong[] Draw() => Enumerable.Range(0, 6).Select(_ => (ulong)plantRng.Next(64) | 1).ToArray();
                    var (t1, t2) = (Template(Times(g, Draw()), false), Template(Times(g, Draw()), true));
                    if (Translates(t1, t2)) continue;
                    cells++;
                    var gcd = SignedMaskBalanceClaim.PatchGcd(t1, t2);
                    if (SignedMaskBalanceClaim.Patch(t1, t2).Corner is null) pairing++;
                    if (!Divides(g, gcd)) { bad++; first ??= $"{t1} + {t2}"; }
                }
                families.Add(new FamilyReading("Theorem G's gcd on 200 planted factors g·u, g·v (g two rows and two columns, u and v drawn in a 6 × 6 window): g divides the gcd read", cells, cells, pairing, bad, first, "the planted factor"));
            }

            // its row M = 1 is Theorem B: two templates of one row up to three sites wide with a letter at each end
            {
                var row = new[] { "XX", "XY", "YY", "XZ", "ZY", "YZ" }
                    .Concat(from a in "XYZ" from mid in "IXYZ" from b in "XYZ" select $"{a}{mid}{b}")
                    .Where(t => SignedMaskBalanceClaim.Mask(t) != 0).Distinct().ToArray();
                int cells = 0, readings = 0, pairing = 0, bad = 0;
                string? first = null;
                for (int i = 0; i < row.Length; i++)
                    for (int j = i + 1; j < row.Length; j++)
                    {
                        ulong ma = SignedMaskBalanceClaim.Mask(row[i]), mb = SignedMaskBalanceClaim.Mask(row[j]);
                        if (ma >> BitOperations.TrailingZeroCount(ma) == mb >> BitOperations.TrailingZeroCount(mb)) continue;
                        cells++;
                        var (bBreaks, from) = SignedMaskBalanceClaim.Chain(row[i], row[j]);
                        var (gBreaks, corner) = SignedMaskBalanceClaim.Patch(row[i], row[j]);
                        readings++;
                        if (!bBreaks) pairing++;
                        if (bBreaks != gBreaks || (gBreaks && corner != (from!.Value, 1))) { bad++; first ??= $"{row[i]} + {row[j]}"; }
                    }
                families.Add(new FamilyReading("Theorem G's row M = 1, two templates of one row up to three sites wide: the corner (N*, 1)", cells, readings, pairing, bad, first, "Theorem B"));
            }

            // three templates: no closed form there; the breaking patches are the union of the quadrants of the unbalanced
            // relations, and the sets whose patches to 8 × 8 show more than one corner are counted. Closure upwards holds by
            // construction (a smaller patch's placements map into a larger one's), so its check tests the placement, not a law
            {
                var tripleRng = new Random(seed + 4);
                int cells = 0, readings = 0, pairing = 0, notClosed = 0;
                string? first = null;
                for (int k = 0; k < 500; k++)
                {
                    var ts = new[] { pool3[tripleRng.Next(pool3.Length)], pool3[tripleRng.Next(pool3.Length)], pool3[tripleRng.Next(pool3.Length)] };
                    if (Translates(ts[0], ts[1]) || Translates(ts[0], ts[2]) || Translates(ts[1], ts[2])) continue;
                    cells++;
                    var breaks = new bool[9, 9];
                    for (int l = 1; l <= 8; l++)
                        for (int m = 1; m <= 8; m++) { breaks[l, m] = !SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnPatch(ts, l, m)); readings++; if (!breaks[l, m]) pairing++; }
                    for (int l = 1; l <= 8; l++)
                        for (int m = 1; m <= 8; m++)
                            if (breaks[l, m] && ((l < 8 && !breaks[l + 1, m]) || (m < 8 && !breaks[l, m + 1]))) { notClosed++; first ??= $"{string.Join(" + ", ts)} at {l} × {m}"; }
                    int corners = 0;
                    for (int l = 1; l <= 8; l++)
                        for (int m = 1; m <= 8; m++)
                            if (breaks[l, m] && (l == 1 || !breaks[l - 1, m]) && (m == 1 || !breaks[l, m - 1])) corners++;
                    if (corners > 1) tripleSeveralCorners++;
                }
                families.Add(new FamilyReading("three templates of the 3 × 3 window: the breaking patches up to 8 × 8 (closed upwards, a check of the placement)", cells, readings, pairing, notClosed, first, "the balance"));
            }

            // the worked example and its counterpart: Theorem 2 on the patches to 16 sites, the law and the balance to 8 × 8
            {
                var examples = new (string Name, string[] Ts, Func<int, int, bool> Pairs)[]
                {
                    ("XY across beside XX/X breaks exactly from 3 × 2", new[] { "XY", "XX/X" }, (l, m) => l < 3 || m < 2),
                    ("XX across beside XY down pairs on every patch", new[] { "XX", "X/Y" }, (l, m) => true),
                };
                int readings = 0, pairing = 0, bad = 0;
                string? first = null;
                foreach (var (name, ts, expect) in examples)
                    for (int l = 1; l <= 8; l++)
                        for (int m = 1; m <= 8; m++)
                        {
                            bool law = !SignedMaskBalanceClaim.PatchBreaks(ts[0], ts[1], l, m);
                            bool balance = SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnPatch(ts, l, m));
                            bool? graph = l * m <= 16 ? PairsOnPatch(l, m, ts) : null;
                            readings++;
                            if (law) pairing++;
                            if (law != expect(l, m) || balance != law || (graph is bool g && g != law)) { bad++; first ??= $"{name}: {l} × {m}"; }
                        }
                families.Add(new FamilyReading("Theorem G's worked example and its counterpart (Theorem 2 to 16 sites, the law and the balance to 8 × 8)", examples.Length, readings, pairing, bad, first, "the page"));
            }
        }

        Families = families;
        PatchNoGcdMismatches = patchNoGcd;
        PatchOneDirectionMismatches = patchOneDirection;
        PatchControlReadings = patchControlReadings;
        PatchBreakingPairs = patchBreaking;
        PatchPairsWhereTheGcdMatters = patchGcdMatters;
        PatchTriplesWithSeveralCorners = tripleSeveralCorners;
        OneCornerMismatches = oneCorner;
        StaircaseReadings = staircaseReadings;
        CornerSetsWithSeveral = severalCorners;
        OneCornerMismatchesOnThreeByThree = oneCornerThree;
        TighterBoundMismatches = tighter;
        MergedSitesMismatches = merged;
        LatticeReadings = latticeReadings;
        MergedSitesReadings = mergedReadings;
        OneNumberMismatches = oneNumber;
        OneNumberReadings = oneNumberReadings;
        FoldControlReadings = foldReadings;
        SwappedFoldMismatches = swapped;
        ParityFoldMismatches = parity;
    }

    public string DisplayName => "TorusMaskBalanceWitness (F173 Theorems F and G)";

    public string Summary =>
        $"Theorems F and G on {Families.Count} families, {Families.Sum(f => f.Readings)} readings, {Families.Sum(f => f.Mismatches)} disagreeing; " +
        $"the controls disagree with the balance: one number per mask {OneNumberMismatches} of {OneNumberReadings} (T₁, T₂ ≥ 2), the fold with " +
        $"T₁ and T₂ swapped {SwappedFoldMismatches} and the fold onto the parities {ParityFoldMismatches} of {FoldControlReadings}, each direction " +
        $"on its own (one corner) {OneCornerMismatches} of {StaircaseReadings}, the corner bound one doubling tighter {TighterBoundMismatches} corners, the " +
        $"lattice sites of a cell merged {MergedSitesMismatches} of {MergedSitesReadings}; Theorem G on the open patch without the gcd " +
        $"{PatchNoGcdMismatches} and along x alone {PatchOneDirectionMismatches} of {PatchControlReadings}";

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
                         $"the parities of L and M: {ParityFoldMismatches} of {FoldControlReadings}; each direction on its own, one corner at the least " +
                         $"e and the least f of the corners: {OneCornerMismatches} of {StaircaseReadings} ({OneCornerMismatchesOnThreeByThree} of them on 3 × 3; {CornerSetsWithSeveral} template sets with more than one corner); the corner bound one doubling tighter, 2^e < n_u " +
                         $"and 2^f < n_v: {TighterBoundMismatches} read corners outside it; the lattice fold with the Q sites of a cell merged: {MergedSitesMismatches} of {MergedSitesReadings} (where the merged set can be placed); " +
                         $"Theorem G with the masks as their own quotients (no gcd): {PatchNoGcdMismatches} of {PatchControlReadings}, and with the threshold along x alone: " +
                         $"{PatchOneDirectionMismatches} of {PatchControlReadings} ({PatchBreakingPairs} pairs with c = 1, {PatchPairsWhereTheGcdMatters} whose corner or c the gcd changes); " +
                         $"three templates with more than one corner on the patches to 8 × 8: {PatchTriplesWithSeveralCorners}",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
