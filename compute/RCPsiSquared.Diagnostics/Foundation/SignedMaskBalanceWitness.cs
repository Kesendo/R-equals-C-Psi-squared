using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live battery for F173 (<c>docs/proofs/PROOF_SIGNED_MASK_BALANCE.md</c>, typed
/// <see cref="SignedMaskBalanceClaim"/>): the claim's readings, GF(2) on masks and signs and, where a flip set carries several
/// strings, Gaussian integers on its ratio and read set, held against
/// Theorem 2 of the complement connection walked on the bitstrings (<see cref="ComplementConnectionGraph"/>, every entry
/// a Gaussian integer, every edge cross-multiplied), Z dephasing on every site, templates placed with their identity
/// letters at every start of the open chain or at every start of the ring, wrapped; Theorem D also against the balance
/// (<see cref="SignedMaskBalanceClaim.Balanced"/>, Theorem A's elimination) on rings to N = 32, past the bitstrings' reach.
///
/// <para>Families, each a count of cells and readings and the readings where the two routes disagree (none expected):
/// two templates of width three with masks of different shapes, both signs of the second coefficient, at N* − 1, N*,
/// N* + 1 where c = 1 and at N = 3 to 6 where c = 0; two templates of unequal widths, where N* pairs each width with the
/// other template's quotient; the diagonal-cell masks at k = 3 and 4 against F115's own valuation
/// (<see cref="WindowedObstructionScan.IsHardPair"/>) and against Theorem 2 at N* − 1 and N*; one template on the ring
/// from N = 2w − 1 to the ring cap; XY + ZZ on the ring; Theorem D on two templates against Theorem 2 to the ring cap,
/// and on two templates up to five sites wide and three templates against the balance to N = 32; three templates at one
/// N against Theorem 2; Theorem E on two strings on one flip mask, chain and ring, Corollary E2 on the ring (also against
/// Theorem E's elimination past the bitstrings' reach) and on the products (aX + bY)⊗(cX + dY), a, b, c, d = ±1,
/// Corollary E3 on two flip shapes on the chain; Theorem B′ on three and four templates on the chain against the balance
/// and Theorem 2, and at two templates against Theorem B; Theorem E's ℤ/4 form, the refinement and its carrier against
/// Theorem 2's (★); Corollary E3′ against E3 at two templates and against Theorem 2 at three; the long chain against
/// Theorem D on the ring of 64; the worked examples of the page. Sixteen controls show the
/// battery can fail: the widths swapped in
/// N*, the quotients swapped in c and the power 2^(v₂(N)) dropped from the ring law each disagree with Theorem 2 somewhere,
/// Theorem D with the parity of N for 2^(v₂(N)) or with i₀ at the greatest valuation disagrees with the balance, and
/// Theorem E with any one of its three conditions dropped, Corollary E2 without its read-set test or its ordering term and
/// Corollary E3 without its read-sum test or its ordering term disagree with Theorem 2, and Theorem B′ read off the
/// unreduced basis of the relations disagrees with the balance, the ℤ/4 carrier without its coset term fails (★), and
/// Corollary E3′ with N₁* read away or with hol walked without its ordering term disagrees with Theorem 2.</para></summary>
public sealed class SignedMaskBalanceWitness : IInspectable
{
    /// <summary>The largest ring read (the graph has 2^N vertices).</summary>
    public const int MaxRing = 12;

    public sealed record FamilyReading(string Name, int Cells, int Readings, int Mismatches, string? FirstMismatch, string Against = "Theorem 2")
    {
        public bool Agrees => Mismatches == 0;
    }

    public int RingCap { get; }
    public IReadOnlyList<FamilyReading> Families { get; }
    public IReadOnlyDictionary<int, int> F115FloorsAtK4 { get; }
    public int SwappedWidthsMismatches { get; }
    public int SwappedQuotientsMismatches { get; }
    public int DroppedRingPowerMismatches { get; }
    public int ParityOnlyMismatches { get; }
    public int MostValuedMismatches { get; }
    /// <summary>Theorem E's readings that fall outside its setting (an absent edge or a ratio that is no character).</summary>
    public int CharacterOutside { get; }
    /// <summary>Theorem E with condition (1), (2) or (3) dropped: readings against Theorem 2.</summary>
    public IReadOnlyList<int> DroppedConditionMismatches { get; }
    /// <summary>Corollary E2 without its read-set test: readings against Theorem 2.</summary>
    public int DroppedReadSetMismatches { get; }
    /// <summary>Corollary E2 with hol(g) walked without its ordering term: readings against Theorem 2.</summary>
    public int DroppedOrderingMismatches { get; }
    /// <summary>Corollary E3 without its read-sum test: readings against Theorem 2.</summary>
    public int DroppedReadSumMismatches { get; }
    /// <summary>Corollary E3 with hol(u₀) walked without its ordering term: readings against Theorem 2.</summary>
    public int DroppedE3OrderingMismatches { get; }
    /// <summary>Theorem B′ read off the unreduced basis of the relations: readings against the balance.</summary>
    public int UnreducedBasisMismatches { get; }
    /// <summary>Theorem E's ℤ/4 carrier without its coset term 2·⟨R(w), x₀⟩: paired readings where it fails (★).</summary>
    public int DroppedCosetTermMismatches { get; }
    /// <summary>Corollary E3′ with N₁* read away (hol on all of S): readings against Theorem 2.</summary>
    public int DroppedE3PrimeReadSumMismatches { get; }
    /// <summary>Corollary E3′ with hol walked without its ordering term: readings against Theorem 2.</summary>
    public int DroppedE3PrimeOrderingMismatches { get; }
    /// <summary>How often each branch of Corollary E3′ decided a three-shape chain read against Theorem 2.</summary>
    public IReadOnlyDictionary<SignedMaskBalanceClaim.ChainShapeBranch, int> E3PrimeBranches { get; private set; } = new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>();
    /// <summary>How many of the ℤ/4 readings pair, each carrier checked against (★).</summary>
    public int Z4Paired { get; private set; }
    /// <summary>How many readings of the four-string ℤ/4 family pair, each carrier checked against (★).</summary>
    public int Z4FourStringPaired { get; private set; }
    /// <summary>How often each branch of Corollary E3 decided a chain read against Theorem 2.</summary>
    public IReadOnlyDictionary<SignedMaskBalanceClaim.ChainShapeBranch, int> E3Branches { get; private set; } = new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>();
    /// <summary>How often each branch of Corollary E2 decided a ring read against Theorem 2.</summary>
    public IReadOnlyDictionary<SignedMaskBalanceClaim.RingShapeBranch, int> E2Branches { get; private set; } = new Dictionary<SignedMaskBalanceClaim.RingShapeBranch, int>();
    /// <summary>How many of Theorem B′'s drawn cells break for large N.</summary>
    public int BPrimeBreaking { get; private set; }
    /// <summary>How many of those cells have an unreduced basis whose leading coefficient matrix falls short of rank k − 1
    /// (the rank test can fail).</summary>
    public int UnreducedBasisDeficient { get; private set; }
    public bool AllAgree => Families.All(f => f.Agrees);
    public bool ControlsBreak => SwappedWidthsMismatches > 0 && SwappedQuotientsMismatches > 0 && DroppedRingPowerMismatches > 0
                                 && ParityOnlyMismatches > 0 && MostValuedMismatches > 0
                                 && DroppedConditionMismatches.All(d => d > 0) && DroppedReadSetMismatches > 0 && DroppedOrderingMismatches > 0
                                 && DroppedReadSumMismatches > 0 && DroppedE3OrderingMismatches > 0
                                 && UnreducedBasisMismatches > 0 && DroppedCosetTermMismatches > 0
                                 && DroppedE3PrimeReadSumMismatches > 0 && DroppedE3PrimeOrderingMismatches > 0;

    static readonly string[] Windows = (from a in "XYZ" from m in "IXYZ" from b in "XYZ" select $"{a}{m}{b}").ToArray();

    /// <summary>Theorem 2 on the open chain of n sites, each template placed with its identity letters at every start
    /// 0 … n − w.</summary>
    public static bool PairsOnChain(int n, params (string T, long C)[] terms)
    {
        var placed = new List<(PauliMask, BigInteger)>();
        foreach (var (t, c) in terms)
            for (int p = 0; p + t.Length <= n; p++)
                placed.Add((PauliMask.Parse(new string('I', p) + t + new string('I', n - p - t.Length)), c));
        var (near, far, _) = ComplementConnectionGraph.EverySiteDephased(n, ComplementConnectionGraph.Columns(n, placed));
        return near == far;
    }

    /// <summary>Theorem 2 on the ring of n sites, each template at every start, wrapped (n ≥ its width).</summary>
    public static bool PairsOnRing(int n, params (string T, long C)[] terms)
    {
        var placed = new List<(PauliMask, BigInteger)>();
        foreach (var (t, c) in terms)
        {
            if (t.Length > n) throw new ArgumentException($"{t} does not fit on the ring of {n} sites");
            for (int p = 0; p < n; p++)
            {
                var letters = Enumerable.Repeat('I', n).ToArray();
                for (int i = 0; i < t.Length; i++) letters[(p + i) % n] = t[i];
                placed.Add((PauliMask.Parse(new string(letters)), c));
            }
        }
        var (near, far, _) = ComplementConnectionGraph.EverySiteDephased(n, ComplementConnectionGraph.Columns(n, placed));
        return near == far;
    }

    // Corollary E3's templates of the given widths: one string with a letter at each end, or two different strings on one
    // mask whose letters together reach both ends (coefficients 1 and ±1, each unordered pair once), kept where the flip
    // set's ratio is a constant times a character
    static List<((string T, long C)[] Strings, SignedMaskBalanceClaim.FlipRatio Ratio)> OneShapeTemplates(params int[] widths)
    {
        var result = new List<((string T, long C)[], SignedMaskBalanceClaim.FlipRatio)>();
        foreach (int w in widths)
        {
            for (ulong f = 1; f < 1UL << w; f++)
            {
                var onMask = new List<string> { "" };
                for (int j = 0; j < w; j++) onMask = onMask.SelectMany(s => (((f >> j) & 1) == 1 ? "XY" : "IZ").Select(d => s + d)).ToList();
                var candidates = new List<(string T, long C)[]>();
                foreach (var t in onMask.Where(t => t[0] != 'I' && t[^1] != 'I')) candidates.Add(new[] { (t, 1L) });
                for (int a = 0; a < onMask.Count; a++)
                    for (int b = a + 1; b < onMask.Count; b++)
                    {
                        string s = onMask[a], t = onMask[b];
                        if ((s[0] == 'I' && t[0] == 'I') || (s[^1] == 'I' && t[^1] == 'I')) continue;
                        foreach (long c in new[] { 1L, -1L }) candidates.Add(new[] { (s, 1L), (t, c) });
                    }
                foreach (var cand in candidates)
                {
                    var r = SignedMaskBalanceClaim.ReadRatio(f, SignedMaskBalanceClaim.PlaceStrings(cand, w, ring: false));
                    if (r.Kind == SignedMaskBalanceClaim.RatioKind.Character) result.Add((cand, r));
                }
            }
        }
        return result;
    }

    static bool SameShape(ulong a, ulong b) => a >> BitOperations.TrailingZeroCount(a) == b >> BitOperations.TrailingZeroCount(b);

    public SignedMaskBalanceWitness(int ringCap = 10)
    {
        if (ringCap < 5 || ringCap > MaxRing) throw new ArgumentOutOfRangeException(nameof(ringCap), ringCap, $"the ring cap runs from 5 to {MaxRing}");
        RingCap = ringCap;
        var families = new List<FamilyReading>();
        int swappedWidths = 0, swappedQuotients = 0, droppedPower = 0, parityOnly = 0, mostValued = 0, characterOutside = 0, droppedReadSet = 0, droppedOrdering = 0, droppedReadSum = 0, droppedE3Ordering = 0, unreducedBasis = 0, droppedCosetTerm = 0, droppedE3PrimeReadSum = 0, droppedE3PrimeOrdering = 0;
        var droppedCondition = new int[3];

        // two templates of width three, masks of different shapes, both signs
        {
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            for (int i = 0; i < Windows.Length; i++)
                for (int j = i + 1; j < Windows.Length; j++)
                {
                    string a = Windows[i], b = Windows[j];
                    ulong ma = SignedMaskBalanceClaim.Mask(a), mb = SignedMaskBalanceClaim.Mask(b);
                    if (ma == 0 || mb == 0 || SameShape(ma, mb)) continue;
                    var (breaks, from) = SignedMaskBalanceClaim.Chain(a, b);
                    foreach (long sign in new[] { 1L, -1L })
                    {
                        cells++;
                        var ns = breaks ? new[] { from!.Value - 1, from.Value, from.Value + 1 }.Where(n => n >= 3) : Enumerable.Range(3, 4);
                        foreach (int n in ns)
                        {
                            bool expected = !breaks || n < from;
                            bool pairs = PairsOnChain(n, (a, 1), (b, sign));
                            bool balanced = SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(new[] { a, b }, n));
                            bool character = SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(new[] { (a, 1L), (b, sign) }, n, ring: false))
                                             == SignedMaskBalanceClaim.CharacterVerdict.Pairs;
                            readings++;
                            if (pairs != expected || balanced != expected || character != expected) { mismatches++; first ??= $"{a} {sign:+0;-0}{b} at N = {n}"; }
                            if (pairs != (SwappedQuotientC(a, b) == 0 || n < NStar(a, b))) swappedQuotients++;
                        }
                    }
                }
            families.Add(new FamilyReading("two templates of width three, masks of different shapes, both signs (the balance and Theorem E agreeing)", cells, readings, mismatches, first));
        }

        // unequal widths
        {
            var two = from x in "XYZ" from y in "XYZ" select $"{x}{y}";
            var longer = new[] { "XIXIX", "XXIX", "YZIX", "XYIZY" };
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            foreach (var a in two.Concat(Windows))
                foreach (var b in Windows.Concat(longer).Where(b => b.Length > a.Length))
                {
                    ulong ma = SignedMaskBalanceClaim.Mask(a), mb = SignedMaskBalanceClaim.Mask(b);
                    if (ma == 0 || mb == 0 || SameShape(ma, mb)) continue;
                    var (breaks, from) = SignedMaskBalanceClaim.Chain(a, b);
                    cells++;
                    // read only where both templates are placed, and where c = 0 from just below the first relation to
                    // just past it, so that every reading but the lowest holds a cycle the law could close or break
                    int star = NStar(a, b);
                    var ns = Enumerable.Range(star - 1, 3).Where(n => n >= Math.Max(3, b.Length));
                    foreach (int n in ns)
                    {
                        bool pairs = PairsOnChain(n, (a, 1), (b, 1));
                        readings++;
                        if (pairs != (!breaks || n < from)) { mismatches++; first ??= $"{a} + {b} at N = {n}"; }
                        if (breaks && pairs != (n < SwappedWidthsFrom(a, b))) swappedWidths++;
                    }
                }
            families.Add(new FamilyReading("two templates of unequal widths, placed with their identity letters", cells, readings, mismatches, first));
        }

        // the F115 face: diagonal-cell masks at k = 3, 4
        {
            var floors = new Dictionary<int, int>();
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            foreach (int k in new[] { 3, 4 })
            {
                var templates = new List<string>();
                for (ulong m = 1; m < 1UL << k; m++)
                {
                    if (BitOperations.PopCount(m) % 2 != 0) continue;
                    var letters = Enumerable.Range(0, k).Select(j => ((m >> j) & 1) == 1 ? 'X' : 'Z').ToArray();
                    if (SignedMaskBalanceClaim.Sign(new string(letters)) == 0) letters[BitOperations.TrailingZeroCount(m)] = 'Y';
                    templates.Add(new string(letters));
                }
                foreach (var a in templates)
                    foreach (var b in templates)
                    {
                        ulong ma = SignedMaskBalanceClaim.Mask(a), mb = SignedMaskBalanceClaim.Mask(b);
                        if (string.CompareOrdinal(a, b) >= 0 || SameShape(ma, mb)) continue;
                        cells++;
                        var (breaks, from) = SignedMaskBalanceClaim.Chain(a, b);
                        if (breaks != WindowedObstructionScan.IsHardPair(ma, mb)) { mismatches++; first ??= $"{a} + {b}: valuation"; }
                        if (!breaks) continue;
                        if (k == 4) floors[from!.Value] = floors.GetValueOrDefault(from.Value) + 1;
                        foreach (int n in new[] { from!.Value - 1, from.Value })
                        {
                            readings++;
                            if (PairsOnChain(n, (a, 1), (b, 1)) != (n < from)) { mismatches++; first ??= $"{a} + {b} at N = {n}"; }
                        }
                    }
            }
            F115FloorsAtK4 = floors;
            families.Add(new FamilyReading("the F115 face: diagonal-cell masks at k = 3, 4, the valuation and the floor", cells, readings, mismatches, first));
        }

        // one template on the ring, and XY + ZZ
        {
            var templates = (from a in "XYZ" from b in "XYZ" select $"{a}{b}").Concat(Windows)
                .Concat(new[] { "XYYY", "XIIY", "XXIY", "YIIY", "XIIIY", "YIIIY", "XYIYX" })
                .Where(t => SignedMaskBalanceClaim.Mask(t) != 0).ToList();
            int readings = 0, mismatches = 0;
            string? first = null;
            foreach (var t in templates)
                for (int n = Math.Max(3, 2 * t.Length - 1); n <= RingCap; n++)
                {
                    bool breaks = SignedMaskBalanceClaim.RingBreaks(t, n);
                    bool pairs = PairsOnRing(n, (t, 1));
                    readings++;
                    if (pairs == breaks) { mismatches++; first ??= $"{t} on the ring of {n}"; }
                    bool dropped = SignedMaskBalanceClaim.Sign(t) == 1 && SignedMaskBalanceClaim.Valuation(SignedMaskBalanceClaim.Mask(t)) >= 1;
                    if (pairs == dropped) droppedPower++;
                }
            families.Add(new FamilyReading("one template on the ring, N = 2w − 1 to the cap", templates.Count, readings, mismatches, first));

            int r2 = 0, m2 = 0;
            string? f2 = null;
            for (int n = 3; n <= RingCap; n++)
            {
                r2++;
                if (PairsOnRing(n, ("XY", 1), ("ZZ", 1)) != (n % 2 == 0)) { m2++; f2 ??= $"XY + ZZ on the ring of {n}"; }
            }
            families.Add(new FamilyReading("XY + ZZ on the ring, paired exactly at even N", 1, r2, m2, f2));
        }

        // Theorem D: two templates on the ring, against Theorem 2 to the ring cap and against the balance to N = 32
        {
            var twoSite = (from x in "XYZ" from y in "XYZ" select $"{x}{y}").Where(x => SignedMaskBalanceClaim.Mask(x) != 0);
            var pool = twoSite.Concat(Windows.Where(x => SignedMaskBalanceClaim.Mask(x) != 0)).ToList();
            var wide = pool.Concat(new[] { "XYYY", "XIIY", "XXIY", "YIIY", "XYIZ", "YZZX", "XIIIY", "YIIIY", "XYIYX" }).ToList();
            int cells = 0, readings = 0, mismatches = 0, rcells = 0, rreadings = 0, rmismatches = 0;
            string? first = null, rfirst = null;
            for (int i = 0; i < wide.Count; i++)
                for (int j = i + 1; j < wide.Count; j++)
                {
                    string a = wide[i], b = wide[j];
                    if (SameShape(SignedMaskBalanceClaim.Mask(a), SignedMaskBalanceClaim.Mask(b))) continue;
                    var ab = new[] { a, b };
                    int lo = Math.Max(3, 2 * Math.Max(a.Length, b.Length) - 1);
                    rcells++;
                    for (int n = lo; n <= 32; n++)
                    {
                        bool breaks = SignedMaskBalanceClaim.RingBreaks(ab, n);
                        rreadings++;
                        if (SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnRing(ab, n)) == breaks) { rmismatches++; rfirst ??= $"{a} + {b} on the ring of {n}"; }
                        if (RingParityOnly(ab, n) != breaks) parityOnly++;
                        if (RingAtTheMostValued(ab, n) != breaks) mostValued++;
                    }
                    foreach (long sign in j < pool.Count ? new[] { 1L, -1L } : new[] { 1L })
                    {
                        cells++;
                        for (int n = lo; n <= RingCap; n++)
                        {
                            readings++;
                            if (PairsOnRing(n, (a, 1), (b, sign)) == SignedMaskBalanceClaim.RingBreaks(ab, n)) { mismatches++; first ??= $"{a} {sign:+0;-0}{b} on the ring of {n}"; }
                        }
                    }
                }
            families.Add(new FamilyReading("Theorem D: two templates on the ring, against Theorem 2 to the cap (both signs below four sites wide)", cells, readings, mismatches, first));
            families.Add(new FamilyReading("Theorem D: two templates up to five sites wide on the ring, against the balance to N = 32", rcells, rreadings, rmismatches, rfirst, "the balance"));

            // three templates, valuations 0 to 4, so that T > v₀ and T ≤ v₀ both meet ties and greater valuations
            var sixteen = new[] { "XY", "XX", "XZ", "XXX", "XIY", "XYZ", "YIX", "XZY", "ZXY", "XXIY", "XIIY", "XYYY", "YIIY", "XIIIY", "YIIIY", "XYIYX" };
            int tcells = 0, treadings = 0, tmismatches = 0;
            string? tfirst = null;
            for (int i = 0; i < sixteen.Length; i++)
                for (int j = i + 1; j < sixteen.Length; j++)
                    for (int l = j + 1; l < sixteen.Length; l++)
                    {
                        var abc = new[] { sixteen[i], sixteen[j], sixteen[l] };
                        if (abc.Select(x => SignedMaskBalanceClaim.Mask(x)).Select(m => m >> BitOperations.TrailingZeroCount(m)).Distinct().Count() < 3) continue;
                        tcells++;
                        for (int n = Math.Max(3, 2 * abc.Max(x => x.Length) - 1); n <= 32; n++)
                        {
                            bool breaks = SignedMaskBalanceClaim.RingBreaks(abc, n);
                            treadings++;
                            if (SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnRing(abc, n)) == breaks) { tmismatches++; tfirst ??= $"{string.Join(" + ", abc)} on the ring of {n}"; }
                            if (RingParityOnly(abc, n) != breaks) parityOnly++;
                            if (RingAtTheMostValued(abc, n) != breaks) mostValued++;
                        }
                    }
            families.Add(new FamilyReading("Theorem D: three templates of different shapes on the ring, against the balance to N = 32", tcells, treadings, tmismatches, tfirst, "the balance"));
        }

        // three templates of different shapes at one N, by the balance
        {
            var pool = new[] { "XIY", "XXX", "XIZ", "XXZ", "XYX", "YZX", "XZZ" };
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            for (int i = 0; i < pool.Length; i++)
                for (int j = i + 1; j < pool.Length; j++)
                    for (int l = j + 1; l < pool.Length; l++)
                    {
                        string a = pool[i], b = pool[j], c = pool[l];
                        var shapes = new[] { a, b, c }.Select(x => SignedMaskBalanceClaim.Mask(x)).Select(m => m >> BitOperations.TrailingZeroCount(m));
                        if (shapes.Distinct().Count() < 3) continue;
                        cells++;
                        for (int n = 4; n <= 6; n++)
                        {
                            readings++;
                            if (SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(new[] { a, b, c }, n)) != PairsOnChain(n, (a, 1), (b, 1), (c, 1)))
                            { mismatches++; first ??= $"{a} + {b} + {c} on the chain of {n}"; }
                            if (n < 5) continue;
                            readings++;
                            bool ringPairs = PairsOnRing(n, (a, 1), (b, 1), (c, 1));
                            if (SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnRing(new[] { a, b, c }, n)) != ringPairs
                                || SignedMaskBalanceClaim.RingBreaks(new[] { a, b, c }, n) == ringPairs)
                            { mismatches++; first ??= $"{a} + {b} + {c} on the ring of {n}"; }
                        }
                    }
            families.Add(new FamilyReading("three templates of different shapes at one N, by the balance (on the ring also by Theorem D)", cells, readings, mismatches, first));
        }

        // Theorem B′: three and four templates of different shapes on the chain, drawn (seeded) from the templates two to
        // four sites wide with a letter at each end, against the balance at every N from 1 to max(Σw + max w, N* + 2)
        // and against Theorem 2 to N = 8, the reduced basis's leading coefficient matrix of full rank; the same reading off the
        // unreduced basis, the control; and two templates against Theorem B's closed form
        {
            var pool = (from a in "XYZ" from b in "XYZ" select $"{a}{b}").Concat(Windows)
                .Concat(from a in "XYZ" from m1 in "IXYZ" from m2 in "IXYZ" from b in "XYZ" select $"{a}{m1}{m2}{b}")
                .Where(t => SignedMaskBalanceClaim.Mask(t) != 0).ToList();
            var rng = new Random(173);
            int cells = 0, readings = 0, mismatches = 0, graphReadings = 0, graphMismatches = 0, breaking = 0, unreducedDeficient = 0, longChainReadings = 0, longChainMismatches = 0;
            string? first = null, graphFirst = null, longChainFirst = null;
            foreach (int k in new[] { 3, 4 })
                for (int drawn = 0; drawn < (k == 3 ? 400 : 200);)
                {
                    var ts = Enumerable.Range(0, k).Select(_ => pool[rng.Next(pool.Count)]).ToArray();
                    if (ts.Select(SignedMaskBalanceClaim.Mask).Select(m => m >> BitOperations.TrailingZeroCount(m)).Distinct().Count() < k) continue;
                    drawn++;
                    cells++;
                    var (breaks, from) = SignedMaskBalanceClaim.Chain(ts);
                    longChainReadings++;
                    if (breaks != SignedMaskBalanceClaim.RingBreaks(ts, 64)) { longChainMismatches++; longChainFirst ??= string.Join(" + ", ts); }
                    var (_, unreducedFrom) = SignedMaskBalanceClaim.Chain(ts, reduced: false);
                    if (breaks) breaking++;
                    // reduced: the leading coefficient matrix (row j: the coefficient of x^(sdeg bⱼ − wᵢ) in bⱼᵢ) of full rank k − 1,
                    // the property the predictable degree rests on, read apart from the reduction's own pivots
                    if (LeadingRank(SignedMaskBalanceClaim.ChainRelationBasis(ts), ts) != k - 1)
                    { mismatches++; first ??= $"{string.Join(" + ", ts)}: the reduced basis is not reduced"; }
                    if (LeadingRank(SignedMaskBalanceClaim.ChainRelationBasis(ts, reduced: false), ts) != k - 1) unreducedDeficient++;
                    int widest = ts.Max(t => t.Length), top = Math.Max(ts.Sum(t => t.Length) + widest, (from ?? 0) + 2);
                    for (int n = 1; n <= top; n++)
                    {
                        bool law = from is not int f || n < f;
                        bool balanced = SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(ts, n));
                        readings++;
                        if (balanced != law) { mismatches++; first ??= $"{string.Join(" + ", ts)} on the chain of {n}"; }
                        if (balanced != (unreducedFrom is not int g || n < g)) unreducedBasis++;
                        if (n > 8) continue;
                        graphReadings++;
                        if (PairsOnChain(n, ts.Select(t => (t, 1L)).ToArray()) != law) { graphMismatches++; graphFirst ??= $"{string.Join(" + ", ts)} on the chain of {n}"; }
                    }
                }
            families.Add(new FamilyReading("Theorem B′: three and four templates on the chain, against the balance at every N from 1 to past N*", cells, readings, mismatches, first, "the balance"));
            families.Add(new FamilyReading("Theorem B′: the same cells against Theorem 2 to N = 8", cells, graphReadings, graphMismatches, graphFirst));
            families.Add(new FamilyReading("Theorem B′: whether the long chain breaks, against Theorem D on the ring of 64 (T > v₀)", cells, longChainReadings, longChainMismatches, longChainFirst, "Theorem D"));
            BPrimeBreaking = breaking;
            UnreducedBasisDeficient = unreducedDeficient;

            var two = (from a in "XYZ" from b in "XYZ" select $"{a}{b}").Concat(Windows).Where(t => SignedMaskBalanceClaim.Mask(t) != 0).ToList();
            int pairs = 0, pairMismatches = 0;
            string? pairFirst = null;
            for (int i = 0; i < two.Count; i++)
                for (int j = i + 1; j < two.Count; j++)
                {
                    ulong ma = SignedMaskBalanceClaim.Mask(two[i]), mb = SignedMaskBalanceClaim.Mask(two[j]);
                    if (SameShape(ma, mb)) continue;
                    pairs++;
                    if (SignedMaskBalanceClaim.Chain(new[] { two[i], two[j] }) != SignedMaskBalanceClaim.Chain(two[i], two[j]))
                    { pairMismatches++; pairFirst ??= $"{two[i]} + {two[j]}"; }
                }
            families.Add(new FamilyReading("Theorem B′ at two templates against Theorem B's N*", pairs, pairs, pairMismatches, pairFirst, "Theorem B"));
        }

        // Theorem E: two strings on one flip mask in a window of two or three sites, against Theorem 2
        {
            var pairsE = new List<(string A, string B)>();
            foreach (int w in new[] { 2, 3 })
                for (ulong f = 1; f < 1UL << w; f++)
                {
                    var strs = new List<string> { "" };
                    for (int j = 0; j < w; j++)
                    {
                        var letters = ((f >> j) & 1) == 1 ? "XY" : "IZ";
                        strs = strs.SelectMany(s => letters.Select(c => s + c)).ToList();
                    }
                    for (int i = 0; i < strs.Count; i++)
                        for (int j = i + 1; j < strs.Count; j++) pairsE.Add((strs[i], strs[j]));
                }
            int cells = 0, readings = 0, mismatches = 0, e2Readings = 0, e2Mismatches = 0, e2Cells = 0, e2Far = 0, e2FarMismatches = 0;
            int z4Readings = 0, z4Paired = 0, z4Mismatches = 0;
            string? first = null, e2First = null, e2FarFirst = null, z4First = null;
            var branches = new Dictionary<SignedMaskBalanceClaim.RingShapeBranch, int>();
            foreach (var (a, b) in pairsE)
                foreach (long sign in new[] { 1L, -1L })
                {
                    cells++;
                    var terms = new[] { (a, 1L), (b, sign) };
                    int w = a.Length;
                    foreach (bool ring in new[] { false, true })
                        for (int n = ring ? Math.Max(3, 2 * w - 1) : Math.Max(3, w); n <= (ring ? RingCap : 7); n++)
                        {
                            var placed = SignedMaskBalanceClaim.PlaceStrings(terms, n, ring);
                            var v = SignedMaskBalanceClaim.DecideCharacter(placed);
                            if (v == SignedMaskBalanceClaim.CharacterVerdict.Outside) { characterOutside++; continue; }
                            bool pairs = ring ? PairsOnRing(n, (a, 1), (b, sign)) : PairsOnChain(n, (a, 1), (b, sign));
                            readings++;
                            if (pairs != (v == SignedMaskBalanceClaim.CharacterVerdict.Pairs)) { mismatches++; first ??= $"{a} {sign:+0;-0}{b} on the {(ring ? "ring" : "chain")} of {n}"; }
                            // the ℤ/4 form: the refinement exists exactly where Theorem 2 pairs, is one (on W up to 128 elements), and
                            // its carrier i^(−Q₀ − 2⟨R(w), x₀⟩) satisfies Theorem 2's (★) entry by entry; without the coset term it must fail
                            {
                                var form = SignedMaskBalanceClaim.QuadraticRefinement(placed);
                                z4Readings++;
                                string at = $"{a} {sign:+0;-0}{b} on the {(ring ? "ring" : "chain")} of {n}";
                                bool bad = (form is not null) != pairs;
                                if (form is not null)
                                {
                                    z4Paired++;
                                    if (form.Count <= 128 && !IsRefinement(form)) bad = true;
                                    if (!StarHolds(placed, n, SignedMaskBalanceClaim.CarrierQuarters(form, n))) bad = true;
                                    if (!StarHolds(placed, n, SignedMaskBalanceClaim.CarrierQuarters(form, n, cosetTerm: false))) droppedCosetTerm++;
                                }
                                if (bad) { z4Mismatches++; z4First ??= at; }
                            }
                            foreach (var (flag, k) in new[] { ((false, true, true), 0), ((true, false, true), 1), ((true, true, false), 2) })
                                if (pairs != (SignedMaskBalanceClaim.DecideCharacter(placed, flag.Item1, flag.Item2, flag.Item3) == SignedMaskBalanceClaim.CharacterVerdict.Pairs))
                                    droppedCondition[k]++;
                            if (!ring) continue;
                            var (breaks, branch) = SignedMaskBalanceClaim.RingOneShape(terms, n);
                            if (branch == SignedMaskBalanceClaim.RingShapeBranch.Outside) continue;
                            if (n == Math.Max(3, 2 * w - 1)) e2Cells++;
                            e2Readings++;
                            branches[branch] = branches.GetValueOrDefault(branch) + 1;
                            if (pairs == breaks) { e2Mismatches++; e2First ??= $"{a} {sign:+0;-0}{b} on the ring of {n}"; }
                            if (pairs == SignedMaskBalanceClaim.RingOneShape(terms, n, readSetTest: false).Breaks) droppedReadSet++;
                        }
                    for (int n = RingCap + 1; n <= 40; n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.RingOneShape(terms, n);
                        if (branch == SignedMaskBalanceClaim.RingShapeBranch.Outside) continue;
                        var v = SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(terms, n, ring: true));
                        e2Far++;
                        if (breaks != (v == SignedMaskBalanceClaim.CharacterVerdict.Breaks)) { e2FarMismatches++; e2FarFirst ??= $"{a} {sign:+0;-0}{b} on the ring of {n}"; }
                    }
                }
            // four strings on one flip set: the products (aX + bY)⊗(cX + dY), a, b, c, d = ±1, where hol(g) is the ordering term alone
            int pcells = 0, preadings = 0, pmismatches = 0;
            string? pfirst = null;
            foreach (int a1 in new[] { 1, -1 }) foreach (int b1 in new[] { 1, -1 }) foreach (int c1 in new[] { 1, -1 }) foreach (int d1 in new[] { 1, -1 })
            {
                pcells++;
                var terms = new[] { ("XX", (long)(a1 * c1)), ("XY", (long)(a1 * d1)), ("YX", (long)(b1 * c1)), ("YY", (long)(b1 * d1)) };
                for (int n = 3; n <= RingCap; n++)
                {
                    bool pairs = PairsOnRing(n, terms.Select(t => (t.Item1, t.Item2)).ToArray());
                    var (breaks, branch) = SignedMaskBalanceClaim.RingOneShape(terms, n);
                    var v = SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(terms, n, ring: true));
                    preadings++;
                    if (branch == SignedMaskBalanceClaim.RingShapeBranch.Outside || pairs == breaks || pairs != (v == SignedMaskBalanceClaim.CharacterVerdict.Pairs))
                    { pmismatches++; pfirst ??= $"({a1:+0;-0}X {b1:+0;-0}Y)⊗({c1:+0;-0}X {d1:+0;-0}Y) on the ring of {n}"; }
                    if (pairs == SignedMaskBalanceClaim.RingOneShape(terms, n, orderingTerm: false).Breaks) droppedOrdering++;
                }
            }
            families.Add(new FamilyReading("Theorem E: two strings on one flip mask, windows of two and three sites, chain to 7 and ring to the cap", cells, readings, mismatches, first));
            families.Add(new FamilyReading("Corollary E2 and Theorem E on (aX + bY)⊗(cX + dY), four strings on one flip set, ring to the cap", pcells, preadings, pmismatches, pfirst));
            families.Add(new FamilyReading("Theorem E's ℤ/4 form: the refinement where Theorem 2 pairs, its carrier against (★)", cells, z4Readings, z4Mismatches, z4First));
            Z4Paired = z4Paired;

            // the ℤ/4 form on four strings per flip set, where ω is real beside a nonempty read set: the sixteen products
            // (aX + bY)⊗(cX + dY) on the ring and X ⊗ 2·CZ on chain and ring, N = 3 to 8
            int fcells = 0, freadings = 0, fpaired = 0, fmismatches = 0;
            string? ffirst = null;
            var fours = new List<((string, long)[] Terms, string Name, bool[] Rings)>();
            foreach (int a1 in new[] { 1, -1 }) foreach (int b1 in new[] { 1, -1 }) foreach (int c1 in new[] { 1, -1 }) foreach (int d1 in new[] { 1, -1 })
                fours.Add((new[] { ("XX", (long)(a1 * c1)), ("XY", (long)(a1 * d1)), ("YX", (long)(b1 * c1)), ("YY", (long)(b1 * d1)) },
                           $"({a1:+0;-0}X {b1:+0;-0}Y)⊗({c1:+0;-0}X {d1:+0;-0}Y)", new[] { true }));
            fours.Add((new[] { ("XII", 1L), ("XIZ", 1L), ("XZI", 1L), ("XZZ", -1L) }, "X ⊗ 2·CZ", new[] { false, true }));
            foreach (var (terms, name, rings) in fours)
            {
                fcells++;
                foreach (bool ring in rings)
                    for (int n = 3; n <= 8; n++)
                    {
                        var placed = SignedMaskBalanceClaim.PlaceStrings(terms, n, ring);
                        bool pairs = ring ? PairsOnRing(n, terms) : PairsOnChain(n, terms);
                        var form = SignedMaskBalanceClaim.QuadraticRefinement(placed);
                        freadings++;
                        bool bad = (form is not null) != pairs;
                        if (form is not null)
                        {
                            fpaired++;
                            if (form.Count <= 128 && !IsRefinement(form)) bad = true;
                            if (!StarHolds(placed, n, SignedMaskBalanceClaim.CarrierQuarters(form, n))) bad = true;
                        }
                        if (bad) { fmismatches++; ffirst ??= $"{name} on the {(ring ? "ring" : "chain")} of {n}"; }
                    }
            }
            families.Add(new FamilyReading("Theorem E's ℤ/4 form on four strings per flip set, a real ω beside a nonempty read set", fcells, freadings, fmismatches, ffirst));
            Z4FourStringPaired = fpaired;
            families.Add(new FamilyReading("Corollary E2: one flip shape on the ring, against Theorem 2 to the cap", e2Cells, e2Readings, e2Mismatches, e2First));
            families.Add(new FamilyReading("Corollary E2 against Theorem E's elimination, rings past the cap to N = 40", e2Cells, e2Far, e2FarMismatches, e2FarFirst, "Theorem E"));
            E2Branches = branches;
        }

        // Corollary E3: two templates of different flip shapes on the chain, each one string or two on one mask, with a
        // character ratio; (a) at least one with a read set or ω = ±i and (1) holding at every offset, N = 3 to 8; (b) one
        // template four sites wide, where D ≠ 0 and (1) holds at N*, read from N* − 1 (at least 4) to N* + 1
        {
            var small = OneShapeTemplates(1, 2, 3);
            var wide = OneShapeTemplates(4);
            int acells = 0, areadings = 0, amismatches = 0, bcells = 0, breadings = 0, bmismatches = 0, e3pReadings = 0, e3pMismatches = 0;
            string? afirst = null, bfirst = null, e3pFirst = null;
            var e3Branches = new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>();
            string Name(IReadOnlyList<(string T, long C)> t) => string.Join(" ", t.Select((s, k) => k == 0 ? s.T : $"{s.C:+0;-0}{s.T}"));
            for (int i = 0; i < small.Count; i++)
                for (int j = i + 1; j < small.Count; j++)
                {
                    var (ta, ra) = small[i];
                    var (tb, rb) = small[j];
                    if (SameShape(ra.Mask, rb.Mask)) continue;
                    if (ra.Read == 0 && ra.Quarter % 2 == 0 && rb.Read == 0 && rb.Quarter % 2 == 0) continue;
                    int wa = ta[0].T.Length, wb = tb[0].T.Length;
                    int every = Math.Max(wa + wb - 1, 2 * Math.Max(wa, wb) - 1);
                    if (SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, every).Branch == SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail) continue;
                    acells++;
                    for (int n = Math.Max(3, Math.Max(wa, wb)); n <= 8; n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, n);
                        bool pairs = PairsOnChain(n, ta.Concat(tb).ToArray());
                        areadings++;
                        e3Branches[branch] = e3Branches.GetValueOrDefault(branch) + 1;
                        if (pairs == breaks) { amismatches++; afirst ??= $"{Name(ta)} with {Name(tb)} at N = {n}"; }
                        if (pairs == SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, n, orderingTerm: false).Breaks) droppedE3Ordering++;
                        e3pReadings++;
                        if (SignedMaskBalanceClaim.ChainShapes(new[] { ta, tb }, n) != (breaks, branch)) { e3pMismatches++; e3pFirst ??= $"{Name(ta)} with {Name(tb)} at N = {n}"; }
                    }
                }
            foreach (var (ta, ra) in small)
                foreach (var (tb, rb) in wide)
                {
                    if (SameShape(ra.Mask, rb.Mask)) continue;
                    var (nStar, readSum) = SignedMaskBalanceClaim.ChainTwoShapeData(ta, tb);
                    if (readSum == 0 || SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, nStar).Branch != SignedMaskBalanceClaim.ChainShapeBranch.ReadSumNonzero) continue;
                    bcells++;
                    for (int n = Math.Max(4, nStar - 1); n <= nStar + 1; n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, n);
                        bool pairs = PairsOnChain(n, ta.Concat(tb).ToArray());
                        breadings++;
                        e3Branches[branch] = e3Branches.GetValueOrDefault(branch) + 1;
                        if (pairs == breaks) { bmismatches++; bfirst ??= $"{Name(ta)} with {Name(tb)} at N = {n}"; }
                        if (pairs == SignedMaskBalanceClaim.ChainTwoShapes(ta, tb, n, readSumTest: false).Breaks) droppedReadSum++;
                        e3pReadings++;
                        if (SignedMaskBalanceClaim.ChainShapes(new[] { ta, tb }, n) != (breaks, branch)) { e3pMismatches++; e3pFirst ??= $"{Name(ta)} with {Name(tb)} at N = {n}"; }
                    }
                }
            // the field X ± Y beside the sixteen products (aX + bY)⊗(cX + dY): read sets on both, and an ordering term of +1 in
            // hol(u₀), where on the families above it is −1 wherever a read set is nonzero
            int pcells = 0, preadings = 0, pmismatches = 0;
            string? pfirst = null;
            foreach (long fy in new[] { 1L, -1L })
                foreach (int a1 in new[] { 1, -1 }) foreach (int b1 in new[] { 1, -1 }) foreach (int c1 in new[] { 1, -1 }) foreach (int d1 in new[] { 1, -1 })
                {
                    pcells++;
                    var field = new[] { ("X", 1L), ("Y", fy) };
                    var product = new[] { ("XX", (long)(a1 * c1)), ("XY", (long)(a1 * d1)), ("YX", (long)(b1 * c1)), ("YY", (long)(b1 * d1)) };
                    for (int n = 2; n <= 8; n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.ChainTwoShapes(field, product, n);
                        bool pairs = PairsOnChain(n, field.Concat(product).ToArray());
                        preadings++;
                        e3Branches[branch] = e3Branches.GetValueOrDefault(branch) + 1;
                        if (pairs == breaks) { pmismatches++; pfirst ??= $"X {fy:+0;-0}Y beside ({a1:+0;-0}X {b1:+0;-0}Y)⊗({c1:+0;-0}X {d1:+0;-0}Y) at N = {n}"; }
                        if (pairs == SignedMaskBalanceClaim.ChainTwoShapes(field, product, n, orderingTerm: false).Breaks) droppedE3Ordering++;
                        e3pReadings++;
                        if (SignedMaskBalanceClaim.ChainShapes(new IReadOnlyList<(string, long)>[] { field, product }, n) != (breaks, branch)) { e3pMismatches++; e3pFirst ??= $"X {fy:+0;-0}Y beside the product at N = {n}"; }
                    }
                }
            families.Add(new FamilyReading("Corollary E3: two flip shapes on the chain, a read set or ω = ±i, (1) at every offset, N = 3 to 8", acells, areadings, amismatches, afirst));
            families.Add(new FamilyReading("Corollary E3 where the read sum D ≠ 0 brings the break forward to N*, one template four sites wide", bcells, breadings, bmismatches, bfirst));
            families.Add(new FamilyReading("Corollary E3 on the field X ± Y beside (aX + bY)⊗(cX + dY), an ordering term of +1, N = 2 to 8", pcells, preadings, pmismatches, pfirst));
            E3Branches = e3Branches;
            families.Add(new FamilyReading("Corollary E3′ at two templates against Corollary E3, branch for branch, on the three families above", acells + bcells + pcells, e3pReadings, e3pMismatches, e3pFirst, "Corollary E3"));

            // Corollary E3′ at three templates of different flip shapes, against Theorem 2: (a) drawn (seeded) from the
            // one-shape templates up to three sites wide, N from the widest to 8; (b) drawn from those and the four-site ones,
            // kept where (1) holds at every offset (read at N = 12), N from the widest to 9; the controls read N₁* away and
            // walk hol without its ordering term
            var rng = new Random(1173);
            var all = small.Concat(wide).ToList();
            var tBranches = new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>();
            foreach (bool rich in new[] { false, true })
            {
                var from = rich ? all : small;
                int tcells = 0, treadings = 0, tmismatches = 0;
                string? tfirst = null;
                while (tcells < (rich ? 200 : 500))
                {
                    var ts = Enumerable.Range(0, 3).Select(_ => from[rng.Next(from.Count)]).ToArray();
                    if (ts.Select(t => t.Ratio.Mask >> BitOperations.TrailingZeroCount(t.Ratio.Mask)).Distinct().Count() < 3) continue;
                    var terms = ts.Select(t => (IReadOnlyList<(string T, long C)>)t.Strings).ToArray();
                    if (rich && SignedMaskBalanceClaim.ChainShapes(terms, 12).Branch == SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail) continue;
                    tcells++;
                    int widest = ts.Max(t => t.Strings[0].T.Length);
                    for (int n = Math.Max(3, widest); n <= (rich ? 9 : 8); n++)
                    {
                        var (breaks, branch) = SignedMaskBalanceClaim.ChainShapes(terms, n);
                        bool pairs = PairsOnChain(n, ts.SelectMany(t => t.Strings).ToArray());
                        treadings++;
                        tBranches[branch] = tBranches.GetValueOrDefault(branch) + 1;
                        if (pairs == breaks || branch == SignedMaskBalanceClaim.ChainShapeBranch.Outside) { tmismatches++; tfirst ??= $"{string.Join(" with ", ts.Select(t => Name(t.Strings)))} at N = {n}"; }
                        if (pairs == SignedMaskBalanceClaim.ChainShapes(terms, n, readSumTest: false).Breaks) droppedE3PrimeReadSum++;
                        if (pairs == SignedMaskBalanceClaim.ChainShapes(terms, n, orderingTerm: false).Breaks) droppedE3PrimeOrdering++;
                    }
                }
                families.Add(new FamilyReading(rich
                    ? "Corollary E3′: three flip shapes on the chain, (1) at every offset, templates up to four sites wide, N to 9"
                    : "Corollary E3′: three flip shapes on the chain, up to three sites wide, N to 8", tcells, treadings, tmismatches, tfirst));
            }
            E3PrimeBranches = tBranches;
        }

        // the worked examples
        {
            var checks = new List<(string Name, bool Holds)>
            {
                ("IXZX + XIZX: N* = 6", SignedMaskBalanceClaim.Chain("IXZX", "XIZX") == (true, 6)),
                ("IXZX + XIZX pairs at N = 5", PairsOnChain(5, ("IXZX", 1), ("XIZX", 1))),
                ("IXZX + XIZX breaks at N = 6, 7", !PairsOnChain(6, ("IXZX", 1), ("XIZX", 1)) && !PairsOnChain(7, ("IXZX", 1), ("XIZX", 1))),
                ("XZ + XIXIX: N* = 6", SignedMaskBalanceClaim.Chain("XZ", "XIXIX") == (true, 6)),
                ("XY + XXIX: N* = 5", SignedMaskBalanceClaim.Chain("XY", "XXIX") == (true, 5)),
                ("+XIY + XXX: N* = 5", SignedMaskBalanceClaim.Chain("XIY", "XXX") == (true, 5)),
                ("+XIY + XIZ: N* = 5", SignedMaskBalanceClaim.Chain("XIY", "XIZ") == (true, 5)),
                ("+XIY + XXZ: N* = 4", SignedMaskBalanceClaim.Chain("XIY", "XXZ") == (true, 4)),
                ("+XIX + XXX never breaks", SignedMaskBalanceClaim.Chain("XIX", "XXX") == (false, null)),
                ("XX + XZ + XXX breaks from N = 3, where each of its three pairs pairs (Theorem B′, Theorem 2 at 3)",
                    SignedMaskBalanceClaim.Chain(new[] { "XX", "XZ", "XXX" }) == (true, 3)
                    && SignedMaskBalanceClaim.Chain("XX", "XZ") == (false, null) && SignedMaskBalanceClaim.Chain("XX", "XXX") == (false, null)
                    && SignedMaskBalanceClaim.Chain("XZ", "XXX") == (true, 4)
                    && !PairsOnChain(3, ("XX", 1), ("XZ", 1), ("XXX", 1)) && PairsOnChain(3, ("XZ", 1), ("XXX", 1))
                    && PairsOnChain(3, ("XX", 1), ("XZ", 1)) && PairsOnChain(3, ("XX", 1), ("XXX", 1))),
                ("XZ beside X·I·(X + Y) and X·X·(X + Y) breaks at N = 3 by the holonomy of (x, 1, 1), where each of its pairs pairs (Corollary E3′, Theorem 2 at 3)",
                    SignedMaskBalanceClaim.ChainShapes(new IReadOnlyList<(string, long)>[] { new[] { ("XZ", 1L) }, new[] { ("XIX", 1L), ("XIY", 1L) }, new[] { ("XXX", 1L), ("XXY", 1L) } }, 3)
                        == (true, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus)
                    && !PairsOnChain(3, ("XZ", 1), ("XIX", 1), ("XIY", 1), ("XXX", 1), ("XXY", 1))
                    && PairsOnChain(3, ("XZ", 1), ("XIX", 1), ("XIY", 1)) && PairsOnChain(3, ("XZ", 1), ("XXX", 1), ("XXY", 1))
                    && PairsOnChain(3, ("XIX", 1), ("XIY", 1), ("XXX", 1), ("XXY", 1))),
                ("XX + XZ + XIY breaks from N = 3, which its unreduced basis reads as 4 (Theorem B′ and its control, Theorem 2 at 3)",
                    SignedMaskBalanceClaim.Chain(new[] { "XX", "XZ", "XIY" }) == (true, 3)
                    && SignedMaskBalanceClaim.Chain(new[] { "XX", "XZ", "XIY" }, reduced: false) == (true, 4)
                    && !PairsOnChain(3, ("XX", 1), ("XZ", 1), ("XIY", 1))),
                ("XXI + XIX pairs at N = 4 to 6", Enumerable.Range(4, 3).All(n => PairsOnChain(n, ("XXI", 1), ("XIX", 1)))),
                ("XXI + XIY pairs at 3 and breaks at 4", PairsOnChain(3, ("XXI", 1), ("XIY", 1)) && !PairsOnChain(4, ("XXI", 1), ("XIY", 1))),
                ("XXI + XIX pairs on the ring of five", PairsOnRing(5, ("XXI", 1), ("XIX", 1))),
                ("+XXZ ± XZY and its seven relatives pair exactly at even N (Theorem D, N = 5 to 64)",
                    new[] { ("XXZ", "XZY"), ("XXZ", "YZX"), ("XZY", "YYZ"), ("XZY", "ZXX"), ("XZY", "ZYY"), ("YYZ", "YZX"), ("YZX", "ZXX"), ("YZX", "ZYY") }
                        .All(p => Enumerable.Range(5, 60).All(n => SignedMaskBalanceClaim.RingBreaks(new[] { p.Item1, p.Item2 }, n) == (n % 2 == 1)))),
                ("XX + XIY breaks on every ring (Theorem D, N = 5 to 64)", Enumerable.Range(5, 60).All(n => SignedMaskBalanceClaim.RingBreaks(new[] { "XX", "XIY" }, n))),
                ("XXX + XYX breaks on the ring exactly at 3 | N (Corollary E2, N = 5 to 63; Theorem 2 to 9)",
                    Enumerable.Range(5, 59).All(n => SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("XYX", 1L) }, n).Breaks == (n % 3 == 0))
                    && Enumerable.Range(5, 5).All(n => PairsOnRing(n, ("XXX", 1), ("XYX", 1)) == (n % 3 != 0))),
                ("XXX + XYX pairs on the chain, N = 3 to 7", Enumerable.Range(3, 5).All(n => PairsOnChain(n, ("XXX", 1), ("XYX", 1)))),
                ("XZ + ZY pairs on the ring and breaks on the chain, N = 3 to 8",
                    Enumerable.Range(3, 6).All(n => PairsOnRing(n, ("XZ", 1), ("ZY", 1)) && !PairsOnChain(n, ("XZ", 1), ("ZY", 1)))
                    && Enumerable.Range(3, 6).All(n => SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(new[] { ("XZ", 1L), ("ZY", 1L) }, n, true)) == SignedMaskBalanceClaim.CharacterVerdict.Pairs)),
                ("XXX + YYY pairs on every ring, relations at 3 | N (Corollary E2, N = 5 to 63)",
                    Enumerable.Range(5, 59).All(n => !SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("YYY", 1L) }, n).Breaks)),
                ("XIY + XZX pairs on the ring exactly at 4 | N (Corollary E2, N = 5 to 63; Theorem 2 to 9)",
                    Enumerable.Range(5, 59).All(n => SignedMaskBalanceClaim.RingOneShape(new[] { ("XIY", 1L), ("XZX", 1L) }, n).Breaks == (n % 4 != 0))
                    && Enumerable.Range(5, 5).All(n => PairsOnRing(n, ("XIY", 1), ("XZX", 1)) == (n % 4 == 0))),
                ("XXX + YYY beside the field X + Y on every site pairs on every chain, and XXX − YYY beside it breaks from N = 3 (Corollary E3 to 64, Theorem 2 to 8)",
                    Enumerable.Range(3, 62).All(n => !SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("X", 1L), ("Y", 1L) }, new[] { ("XXX", 1L), ("YYY", 1L) }, n).Breaks
                                                    && SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("X", 1L), ("Y", 1L) }, new[] { ("XXX", 1L), ("YYY", -1L) }, n).Breaks)
                    && Enumerable.Range(3, 6).All(n => PairsOnChain(n, ("X", 1), ("Y", 1), ("XXX", 1), ("YYY", 1)) && !PairsOnChain(n, ("X", 1), ("Y", 1), ("XXX", 1), ("YYY", -1)))),
                ("X·(X + Y)·X beside XIIX pairs at N = 3, where XIIX is not placed yet (Theorem 2), breaks at N* = 4 by the read sum, and from 5 by a square (Corollary E3 to 64, Theorem 2 to 8)",
                    SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XXX", 1L), ("XYX", 1L) }, new[] { ("XIIX", 1L) }, 4).Branch == SignedMaskBalanceClaim.ChainShapeBranch.ReadSumNonzero
                    && Enumerable.Range(5, 60).All(n => SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XXX", 1L), ("XYX", 1L) }, new[] { ("XIIX", 1L) }, n).Branch == SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail)
                    && PairsOnChain(3, ("XXX", 1), ("XYX", 1)) && Enumerable.Range(4, 5).All(n => !PairsOnChain(n, ("XXX", 1), ("XYX", 1), ("XIIX", 1)))),
            };
            int bad = checks.Count(c => !c.Holds);
            families.Add(new FamilyReading("the worked examples of the page", checks.Count, checks.Count, bad, checks.FirstOrDefault(c => !c.Holds).Name, "the page"));
        }

        Families = families;
        SwappedWidthsMismatches = swappedWidths;
        SwappedQuotientsMismatches = swappedQuotients;
        DroppedRingPowerMismatches = droppedPower;
        ParityOnlyMismatches = parityOnly;
        MostValuedMismatches = mostValued;
        CharacterOutside = characterOutside;
        DroppedConditionMismatches = droppedCondition;
        DroppedReadSetMismatches = droppedReadSet;
        DroppedOrderingMismatches = droppedOrdering;
        DroppedReadSumMismatches = droppedReadSum;
        DroppedE3OrderingMismatches = droppedE3Ordering;
        UnreducedBasisMismatches = unreducedBasis;
        DroppedCosetTermMismatches = droppedCosetTerm;
        DroppedE3PrimeReadSumMismatches = droppedE3PrimeReadSum;
        DroppedE3PrimeOrderingMismatches = droppedE3PrimeOrdering;
    }

    // a ℤ/4 quadratic refinement of β(a, b) = ⟨R(a), b⟩ on all of W: Q(a ⊕ b) = Q(a) + Q(b) + 2β(a, b), β symmetric
    static bool IsRefinement(IReadOnlyDictionary<ulong, (int Q, ulong R)> form)
    {
        static int Beta(ulong r, ulong b) => BitOperations.PopCount(r & b) & 1;
        foreach (var (a, (qa, ra)) in form)
            foreach (var (b, (qb, rb)) in form)
                if (Beta(ra, b) != Beta(rb, a) || form[a ^ b].Q != (qa + qb + 2 * Beta(ra, b)) % 4) return false;
        return true;
    }

    // Theorem 2's far condition (★) for X^⊗N·diag(i^(d_x)): H_x̄ȳ·i^(d_y) = i^(d_x)·H_xy at every entry, in Gaussian integers
    static bool StarHolds(IReadOnlyList<SignedMaskBalanceClaim.PlacedString> placed, int n, int[] d)
    {
        var h = ComplementConnectionGraph.Columns(n, placed.Select(s => (new PauliMask(s.X, s.Z), (BigInteger)s.C)));
        int full = (1 << n) - 1;
        static GaussianInteger Turn(GaussianInteger z, int k) { for (int j = 0; j < k; j++) z = z * GaussianInteger.I; return z; }
        for (int y = 0; y < h.Length; y++)
            foreach (var (x, hxy) in h[y])
            {
                if (Turn(ComplementConnectionGraph.Entry(h, x ^ full, y ^ full), d[y]) != Turn(hxy, d[x])) return false;
                if (Turn(hxy, d[y ^ full]) != Turn(ComplementConnectionGraph.Entry(h, x ^ full, y ^ full), d[x ^ full])) return false;
            }
        return true;
    }

    // the rank over GF(2) of a relation basis's leading coefficient matrix for the shift by the template widths
    static int LeadingRank(IReadOnlyList<ulong[]> basis, string[] ts)
    {
        var rows = new List<ulong>();
        foreach (var b in basis)
        {
            int top = Enumerable.Range(0, ts.Length).Where(i => b[i] != 0).Max(i => BitOperations.Log2(b[i]) + ts[i].Length);
            ulong row = 0;
            for (int i = 0; i < ts.Length; i++)
            {
                int e = top - ts[i].Length;
                if (e >= 0 && ((b[i] >> e) & 1) == 1) row |= 1UL << i;
            }
            rows.Add(row);
        }
        int rank = 0;
        for (int bit = 0; bit < ts.Length; bit++)
        {
            int pivot = rows.FindIndex(rank, r => ((r >> bit) & 1) == 1);
            if (pivot < 0) continue;
            (rows[rank], rows[pivot]) = (rows[pivot], rows[rank]);
            for (int r = 0; r < rows.Count; r++) if (r != rank && ((rows[r] >> bit) & 1) == 1) rows[r] ^= rows[rank];
            rank++;
        }
        return rank;
    }

    // the controls of Theorem D, read the proof's way from one template i₀: 2^(v₂(N)) read off the parity of N alone,
    // and i₀ taken at the greatest valuation
    static bool RingParityOnly(string[] ts, int n) => RingRead(ts, n % 2 == 1 ? 1 : 64, least: true);
    static bool RingAtTheMostValued(string[] ts, int n) => RingRead(ts, 1 << BitOperations.TrailingZeroCount((uint)n), least: false);

    static bool RingRead(string[] ts, int big, bool least)
    {
        var read = ts.Select(t => (V: SignedMaskBalanceClaim.Valuation(SignedMaskBalanceClaim.Mask(t)), S: SignedMaskBalanceClaim.Sign(t))).ToList();
        var idx = Enumerable.Range(0, read.Count);
        int i0 = least ? idx.MinBy(i => read[i].V) : idx.MaxBy(i => read[i].V);
        var (v0, s0) = read[i0];
        if (v0 >= big) return read.Any(r => r.S == 1);
        return idx.Any(j => j != i0 && (read[j].S + s0 * (read[j].V == v0 ? 1 : 0)) % 2 == 1);
    }

    // the controls: Theorem B with its pieces swapped, and Theorem C without 2^(v₂(N))
    static int SwappedQuotientC(string a, string b)
    {
        ulong m1 = SignedMaskBalanceClaim.Mask(a), m2 = SignedMaskBalanceClaim.Mask(b);
        ulong g = Gcd(m1, m2);
        return (SignedMaskBalanceClaim.Sign(a) * (BitOperations.PopCount(Div(m1, g)) & 1)
                + SignedMaskBalanceClaim.Sign(b) * (BitOperations.PopCount(Div(m2, g)) & 1)) % 2;
    }

    // N* = max(w₁ + deg q₂, w₂ + deg q₁) whatever c is
    static int NStar(string a, string b)
    {
        ulong m1 = SignedMaskBalanceClaim.Mask(a), m2 = SignedMaskBalanceClaim.Mask(b);
        ulong g = Gcd(m1, m2);
        return Math.Max(a.Length + Deg(Div(m2, g)), b.Length + Deg(Div(m1, g)));
    }

    static int SwappedWidthsFrom(string a, string b)
    {
        ulong m1 = SignedMaskBalanceClaim.Mask(a), m2 = SignedMaskBalanceClaim.Mask(b);
        ulong g = Gcd(m1, m2);
        return Math.Max(a.Length + Deg(Div(m1, g)), b.Length + Deg(Div(m2, g)));
    }

    static int Deg(ulong p) => 63 - BitOperations.LeadingZeroCount(p);

    static ulong Gcd(ulong a, ulong b)
    {
        while (b != 0)
        {
            ulong r = a;
            while (r != 0 && Deg(r) >= Deg(b)) r ^= b << (Deg(r) - Deg(b));
            (a, b) = (b, r);
        }
        return a;
    }

    static ulong Div(ulong a, ulong b)
    {
        ulong q = 0;
        while (a != 0 && Deg(a) >= Deg(b)) { int s = Deg(a) - Deg(b); q |= 1UL << s; a ^= b << s; }
        return q;
    }

    public string DisplayName => $"SignedMaskBalanceWitness (ring cap N = {RingCap})";

    public string Summary =>
        $"F173 on {Families.Count} families: {Families.Where(f => f.Against == "Theorem 2").Sum(f => f.Readings)} readings against Theorem 2 " +
        $"on the bitstrings, {Families.Where(f => f.Against == "the balance").Sum(f => f.Readings)} of Theorems D and B′ against the balance, " +
        $"{Families.Where(f => f.Against == "Theorem B").Sum(f => f.Readings)} of Theorem B′ against Theorem B, " +
        $"{Families.Where(f => f.Against == "Corollary E3").Sum(f => f.Readings)} of E3′ against E3, " +
        $"{Families.Where(f => f.Against == "Theorem D").Sum(f => f.Readings)} of the long chain against Theorem D, " +
        $"{Families.Where(f => f.Against == "the page").Sum(f => f.Readings)} worked checks, {Families.Sum(f => f.Mismatches)} disagreeing; the controls (widths swapped, quotients swapped, 2^(v₂(N)) dropped, " +
        $"the parity of N for 2^(v₂(N)), i₀ at the greatest valuation) disagree {SwappedWidthsMismatches}, {SwappedQuotientsMismatches}, " +
        $"{DroppedRingPowerMismatches}, {ParityOnlyMismatches}, {MostValuedMismatches} times; Theorem E with (1), (2), (3) dropped and E2 " +
        $"without its read-set test or its ordering term disagree {string.Join(", ", DroppedConditionMismatches)}, {DroppedReadSetMismatches}, " +
        $"{DroppedOrderingMismatches} times, E3 without its read-sum test or its ordering term {DroppedReadSumMismatches}, {DroppedE3OrderingMismatches} times, " +
        $"B′ off the unreduced basis {UnreducedBasisMismatches} times, the ℤ/4 carrier without its coset term {DroppedCosetTermMismatches} times, " +
        $"E3′ without N₁* or without its ordering term {DroppedE3PrimeReadSumMismatches}, {DroppedE3PrimeOrderingMismatches} times";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var f in Families)
                yield return new InspectableNode(
                    displayName: f.Name,
                    summary: $"{f.Cells} {(f.Cells == 1 ? "cell" : "cells")}, {f.Readings} readings, " +
                             (f.Agrees ? $"every one agreeing with {f.Against}" : $"{f.Mismatches} DISAGREEING, first {f.FirstMismatch}"),
                    provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the floors of F115's pairs at k = 4",
                summary: string.Join(", ", F115FloorsAtK4.OrderBy(kv => kv.Key).Select(kv => $"{kv.Value} pairs from N = {kv.Key}")) +
                         " (to be read beside THE_PALINDROME_CLASSIFIER's count of the pairs already hard at N = 5)",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the controls, which must disagree",
                summary: $"N* with each width beside its own quotient: {SwappedWidthsMismatches} readings against Theorem 2; c with the " +
                         $"quotients swapped: {SwappedQuotientsMismatches}; the ring law without 2^(v₂(N)): {DroppedRingPowerMismatches}; " +
                         $"Theorem D with the parity of N in place of 2^(v₂(N)): {ParityOnlyMismatches}; Theorem D with i₀ at the greatest " +
                         $"valuation: {MostValuedMismatches} (the last two against the balance); Theorem E without the symmetric reading: " +
                         $"{DroppedConditionMismatches[0]}, without the read sets summed over the relations: {DroppedConditionMismatches[1]}, " +
                         $"without hol: {DroppedConditionMismatches[2]}; Corollary E2 without its read-set test: {DroppedReadSetMismatches}, " +
                         $"with hol(g) walked without its ordering term: {DroppedOrderingMismatches}; Corollary E3 without its read-sum test: " +
                         $"{DroppedReadSumMismatches}, with hol(u₀) walked without its ordering term: {DroppedE3OrderingMismatches}; Theorem B′ read " +
                         $"off the unreduced basis of the relations: {UnreducedBasisMismatches} (against the balance); Theorem E's ℤ/4 carrier " +
                         $"without its coset term: {DroppedCosetTermMismatches} (against (★)); Corollary E3′ without N₁*: " +
                         $"{DroppedE3PrimeReadSumMismatches}, with hol walked without its ordering term: {DroppedE3PrimeOrderingMismatches}",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "Theorem E's reach",
                summary: $"{CharacterOutside} readings outside its setting (an absent edge or a ratio that is no character); Corollary E2's " +
                         $"branches against Theorem 2: {string.Join(", ", E2Branches.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key} {kv.Value}"))}; " +
                         $"Corollary E3's: {string.Join(", ", E3Branches.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key} {kv.Value}"))}; " +
                         $"Corollary E3′'s at three shapes: {string.Join(", ", E3PrimeBranches.OrderBy(kv => kv.Key).Select(kv => $"{kv.Key} {kv.Value}"))}",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
