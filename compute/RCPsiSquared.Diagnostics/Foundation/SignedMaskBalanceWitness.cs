using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live battery for F173 (<c>docs/proofs/PROOF_SIGNED_MASK_BALANCE.md</c>, typed
/// <see cref="SignedMaskBalanceClaim"/>): the claim's three readings, pure GF(2) on masks and signs, held against
/// Theorem 2 of the complement connection walked on the bitstrings (<see cref="ComplementConnectionGraph"/>, every entry
/// a Gaussian integer, every edge cross-multiplied), Z dephasing on every site, templates placed with their identity
/// letters at every start of the open chain or at every start of the ring, wrapped.
///
/// <para>Families, each a count of cells and readings and the readings where the two routes disagree (none expected):
/// two templates of width three with masks of different shapes, both signs of the second coefficient, at N* − 1, N*,
/// N* + 1 where c = 1 and at N = 3 to 6 where c = 0; two templates of unequal widths, where N* pairs each width with the
/// other template's quotient; the diagonal-cell masks at k = 3 and 4 against F115's own valuation
/// (<see cref="WindowedObstructionScan.IsHardPair"/>) and against Theorem 2 at N* − 1 and N*; one template on the ring
/// from N = 2w − 1 to the ring cap; XY + ZZ on the ring; sixteen two-template ring cells by the balance at each N; the
/// worked examples of the page. Three controls show the battery can fail: the widths swapped in N*, the quotients
/// swapped in c, and the power 2^(v₂(N)) dropped from the ring law each disagree with Theorem 2 somewhere.</para></summary>
public sealed class SignedMaskBalanceWitness : IInspectable
{
    /// <summary>The largest ring read (the graph has 2^N vertices).</summary>
    public const int MaxRing = 12;

    public sealed record FamilyReading(string Name, int Cells, int Readings, int Mismatches, string? FirstMismatch)
    {
        public bool Agrees => Mismatches == 0;
    }

    public int RingCap { get; }
    public IReadOnlyList<FamilyReading> Families { get; }
    public IReadOnlyDictionary<int, int> F115FloorsAtK4 { get; }
    public int SwappedWidthsMismatches { get; }
    public int SwappedQuotientsMismatches { get; }
    public int DroppedRingPowerMismatches { get; }
    public bool AllAgree => Families.All(f => f.Agrees);
    public bool ControlsBreak => SwappedWidthsMismatches > 0 && SwappedQuotientsMismatches > 0 && DroppedRingPowerMismatches > 0;

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

    static bool SameShape(ulong a, ulong b) => a >> BitOperations.TrailingZeroCount(a) == b >> BitOperations.TrailingZeroCount(b);

    public SignedMaskBalanceWitness(int ringCap = 10)
    {
        if (ringCap < 5 || ringCap > MaxRing) throw new ArgumentOutOfRangeException(nameof(ringCap), ringCap, $"the ring cap runs from 5 to {MaxRing}");
        RingCap = ringCap;
        var families = new List<FamilyReading>();
        int swappedWidths = 0, swappedQuotients = 0, droppedPower = 0;

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
                            readings++;
                            if (pairs != expected || balanced != expected) { mismatches++; first ??= $"{a} {sign:+0;-0}{b} at N = {n}"; }
                            if (pairs != (SwappedQuotientC(a, b) == 0 || n < NStar(a, b))) swappedQuotients++;
                        }
                    }
                }
            families.Add(new FamilyReading("two templates of width three, masks of different shapes, both signs", cells, readings, mismatches, first));
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

        // two templates on the ring, by the balance at each N
        {
            var pairs = new[] { ("XXZ", "XZY"), ("XXZ", "YZX"), ("XZY", "YYZ"), ("XZY", "ZXX"), ("XZY", "ZYY"), ("YYZ", "YZX"),
                                ("YZX", "ZXX"), ("YZX", "ZYY") };
            int cells = 0, readings = 0, mismatches = 0;
            string? first = null;
            foreach (var (a, b) in pairs)
                foreach (long sign in new[] { 1L, -1L })
                {
                    cells++;
                    for (int n = 5; n <= Math.Min(RingCap, 9); n++)
                    {
                        bool balanced = SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnRing(new[] { a, b }, n));
                        readings++;
                        if (balanced != (n % 2 == 0) || PairsOnRing(n, (a, 1), (b, sign)) != balanced) { mismatches++; first ??= $"{a} {sign:+0;-0}{b} on the ring of {n}"; }
                    }
                }
            families.Add(new FamilyReading("two templates on the ring by the balance, paired exactly at even N", cells, readings, mismatches, first));
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
                            if (SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnRing(new[] { a, b, c }, n)) != PairsOnRing(n, (a, 1), (b, 1), (c, 1)))
                            { mismatches++; first ??= $"{a} + {b} + {c} on the ring of {n}"; }
                        }
                    }
            families.Add(new FamilyReading("three templates of different shapes at one N, by the balance", cells, readings, mismatches, first));
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
                ("XXI + XIX pairs at N = 4 to 6", Enumerable.Range(4, 3).All(n => PairsOnChain(n, ("XXI", 1), ("XIX", 1)))),
                ("XXI + XIY pairs at 3 and breaks at 4", PairsOnChain(3, ("XXI", 1), ("XIY", 1)) && !PairsOnChain(4, ("XXI", 1), ("XIY", 1))),
                ("XXI + XIX pairs on the ring of five", PairsOnRing(5, ("XXI", 1), ("XIX", 1))),
            };
            int bad = checks.Count(c => !c.Holds);
            families.Add(new FamilyReading("the worked examples of the page", checks.Count, checks.Count, bad, checks.FirstOrDefault(c => !c.Holds).Name));
        }

        Families = families;
        SwappedWidthsMismatches = swappedWidths;
        SwappedQuotientsMismatches = swappedQuotients;
        DroppedRingPowerMismatches = droppedPower;
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
        $"F173 against Theorem 2 on the bitstrings: {Families.Sum(f => f.Readings)} readings over {Families.Count} families, " +
        $"{Families.Sum(f => f.Mismatches)} disagreeing; the controls (widths swapped, quotients swapped, 2^(v₂(N)) dropped) " +
        $"disagree {SwappedWidthsMismatches}, {SwappedQuotientsMismatches}, {DroppedRingPowerMismatches} times";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var f in Families)
                yield return new InspectableNode(
                    displayName: f.Name,
                    summary: $"{f.Cells} {(f.Cells == 1 ? "cell" : "cells")}, {f.Readings} readings, " +
                             (f.Agrees ? "every one agreeing with Theorem 2" : $"{f.Mismatches} DISAGREEING, first {f.FirstMismatch}"),
                    provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the floors of F115's pairs at k = 4",
                summary: string.Join(", ", F115FloorsAtK4.OrderBy(kv => kv.Key).Select(kv => $"{kv.Value} pairs from N = {kv.Key}")) +
                         " (to be read beside THE_PALINDROME_CLASSIFIER's count of the pairs already hard at N = 5)",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the controls, which must disagree",
                summary: $"N* with each width beside its own quotient: {SwappedWidthsMismatches} readings against Theorem 2; c with the " +
                         $"quotients swapped: {SwappedQuotientsMismatches}; the ring law without 2^(v₂(N)): {DroppedRingPowerMismatches}",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
