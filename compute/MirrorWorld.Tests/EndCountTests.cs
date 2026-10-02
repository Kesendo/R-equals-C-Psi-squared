using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the end count (adopted 2026-09-30, F158). Every literal is pinned to a source
// the object never reads:
//   - the live witness PalindromeTwoEndCountWitness (`inspect --root twoend`), which ranks the
//     Liouvillian itself on its 4^N columns and never forms a Pauli-string span: seven rows at N = 3
//     (Heisenberg chain, bonds XX+YY+ZZ at J = 1, fields 3/10, read 2026-09-30 with
//     `inspect --root twoend --N 3 --deph <deph> --field <field>`, one run per InlineData row, its two
//     strings the row's first two arguments), here as bonds 100 and fields 30, since neither kernel
//     moves under a common scale;
//   - the registry and the claim: the canonical chain holds N + 1 at both ends (F158 section (f2),
//     PalindromeTwoEndCountClaim.CanonicalChainCount);
//   - experiments/THE_PALINDROME_AS_A_COLOURING.md: stage E's sum, the defect cascade and its Tr(H^2 A)
//     = 8 eps, the SWAP row outside the anticommuting-sum grammar, the non-commuting-jump row, and the
//     golden router's colours;
//   - Router's own frame and MirrorGroup's own R, read off those objects.
// The counts are upper bounds by construction (a rank mod p can only drop); where a pin says a count
// EQUALS a value, the value is the source's.
public class EndCountTests
{
    static readonly World W = new();

    static string One(int n, int site, char letter)
    {
        var c = Enumerable.Repeat('I', n).ToArray();
        c[site] = letter;
        return new string(c);
    }

    static string Two(int n, int i, int j, char letter)
    {
        var c = Enumerable.Repeat('I', n).ToArray();
        c[i] = letter;
        c[j] = letter;
        return new string(c);
    }

    // Heisenberg bonds on the given edges, fields on the sites marked in `field` (a letter or '.'),
    // jumps one per site marked in `deph`: the witness's rows, in its own alphabet.
    static EndCount Row(int n, (int, int)[] edges, string deph, string field, long bond = 100, long mag = 30)
    {
        var h = new List<(string, long)>();
        foreach (var (i, j) in edges)
            foreach (char p in "XYZ") h.Add((Two(n, i, j, p), bond));
        for (int l = 0; l < n; l++)
            if (field[l] != '.') h.Add((One(n, l, field[l]), mag));
        var jumps = Enumerable.Range(0, n).Where(l => deph[l] != '.').Select(l => One(n, l, deph[l])).ToList();
        return new EndCount(W, n, h, jumps);
    }

    static (int, int)[] Chain(int n) => Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToArray();
    static (int, int)[] Ring(int n) => Enumerable.Range(0, n).Select(i => (i, (i + 1) % n)).ToArray();

    // ---- the ontology ----

    [Fact]
    public void The_Count_Owns_Its_Two_Ends_Its_Colours_And_Its_Word_And_Inherits_The_Frame()
    {
        var e = Row(3, Chain(3), "ZZZ", "...");
        Assert.Equal(new[] { "near", "far", "colour", "word" }, e.Own);
        Assert.IsType<World>(e.Parent);
        Assert.Equal(new[] { "x", "y", "z" }, e.Inherited);
    }

    // ---- against the witness, which ranks L itself ----

    // The reading each row takes is derived, not read off: each coloured row has one colour per
    // component (a Heisenberg bond forces it) that every field accepts, X^3 or Y^3 on the ZZZ row, XXX
    // on (ZZZ, X.X) and on (Z.Z, X.X), YYY on (XZZ, ..Y); on the Z-field row Tr(H Z_0) = 2^3 * 30 fires
    // first; on XYZ no letter is lit at all three sites, and every single-jump word vanishes because a
    // global pi rotation about an axis orthogonal to the jump's letter fixes the Heisenberg chain and
    // flips the jump, so the far bound 0 decides against the near end's identity; on the last row H^1 A has no IXI to pick up and H^2 A picks up 2 * (100 * 30) from each
    // of the two bonds, 12000, times 2^3.
    [Theory]
    [InlineData("ZZZ", "...", 4, 4, EndCount.Reading.PalindromeByColour)]
    [InlineData("ZZZ", "Z..", 4, 0, EndCount.Reading.BrokenByWord)]
    [InlineData("ZZZ", "X.X", 1, 1, EndCount.Reading.PalindromeByColour)]
    [InlineData("XZZ", "..Y", 1, 1, EndCount.Reading.PalindromeByColour)]
    [InlineData("Z.Z", "X.X", 1, 1, EndCount.Reading.PalindromeByColour)]
    [InlineData("XYZ", "...", 1, 0, EndCount.Reading.BrokenByCount)]
    [InlineData(".X.", "X.X", 6, 2, EndCount.Reading.BrokenByWord)]
    public void The_String_Route_Meets_The_Witness_Liouvillian_Route(string deph, string field, int near, int far, EndCount.Reading reading)
    {
        var e = Row(3, Chain(3), deph, field);
        Assert.Equal((near, far), e.UpperCounts());
        foreach (long p in ModP.Primes) Assert.Equal((near, far), e.NullitiesAt(p));
        Assert.Equal(reading, e.Verdict());
        Assert.Equal(near == far, EndCount.IsPalindrome(reading));
        // the exact lower bounds sit under the counts, and one colouring makes them equal
        var (loNear, loFar) = e.LowerCounts();
        Assert.InRange(loNear, 1, near);
        Assert.InRange(loFar, 0, far);
        if (loFar > 0) Assert.Equal(loNear, loFar);
    }

    // ---- the canonical chain, F158 section (f2), and past the witness's MaxN = 4 ----

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(6)]
    [InlineData(8)]
    [InlineData(10)]
    public void The_Canonical_Chain_Holds_N_Plus_One_At_Both_Ends(int n)
    {
        var e = Row(n, Chain(n), new string('Z', n), new string('.', n));
        Assert.Equal((n + 1, n + 1), e.UpperCounts());
        Assert.Equal(1 << n, e.DarkStrings().Count);     // one axis per site: 2^N strings where L has 4^N columns
        Assert.Equal(1 << n, e.LitStrings().Count);
        // X^N and Y^N are its colourings (a Heisenberg bond forces one colour per component), Z^N and I
        // its dark strings commuting with H: exact lower bounds, 2 at each end, below the counts.
        Assert.Equal(new[] { new string('X', n), new string('Y', n) }.OrderBy(s => s),
            e.Colourings().Select(s => s.ToString(n)).OrderBy(s => s));
        Assert.Equal((2, 2), e.LowerCounts());
        Assert.Equal(EndCount.Reading.PalindromeByColour, e.Verdict(maxPower: 2, maxJumps: 1));
    }

    // the witness's own test pins the ring at N = 4 to five at both ends (PalindromeTwoEndCountWitnessTests)
    [Fact]
    public void The_Ring_At_Four_Sites_Holds_Five_At_Both_Ends()
        => Assert.Equal((5, 5), Row(4, Ring(4), "ZZZZ", "....").UpperCounts());

    // ---- the exact readings ----

    [Fact]
    public void A_Z_Field_On_A_Dephased_Site_Is_Ruled_Out_By_The_Shortest_Odd_Word()
    {
        // F158 section (f5): Tr(H Z_l) = 2^N delta_l, here 2^3 * 30.
        var e = Row(3, Chain(3), "ZZZ", "Z..");
        var w = e.Word();
        Assert.NotNull(w);
        Assert.Equal("H^1·A0", w!.Word);
        Assert.Equal(new BigInteger(8 * 30), w.TraceRe);
        Assert.True(w.TraceIm.IsZero);
        Assert.Equal(EndCount.Reading.BrokenByWord, e.Verdict());
        // the anti-vacuity partner: the palindromic canonical row fires no word in the same budget
        Assert.Null(Row(3, Chain(3), "ZZZ", "...").Word());
    }

    [Fact]
    public void A_Depolarized_Site_Empties_The_Lit_Coset_And_Is_Broken_Exactly_Twice_Over()
    {
        // X, Y and Z on site 0: no string anticommutes with all three, so far is exactly 0 while the
        // identity keeps near >= 1; and the three jump letters multiply to i, a word of trace 4i.
        var h = new List<(string, long)> { ("XX", 1), ("YY", 1), ("ZZ", 1) };
        var e = new EndCount(W, 2, h, new[] { "XI", "YI", "ZI" });
        Assert.Empty(e.LitStrings());
        Assert.Equal(0, e.UpperCounts().Far);
        var w = e.Word(maxJumps: 3);
        Assert.Equal("H^0·A0·H^0·A1·H^0·A2", w!.Word);
        Assert.Equal((BigInteger.Zero, new BigInteger(4)), (w.TraceRe, w.TraceIm));
        Assert.Equal(EndCount.Reading.BrokenByWord, e.Verdict(maxJumps: 3));
        Assert.Equal(EndCount.Reading.BrokenByCount, e.Verdict());   // with one jump letter per word the counts alone decide it
    }

    // The triple words, with nonzero powers in them. The pin is fw.odd_word_obstruction's own answer on
    // the XYZ row (simulations/framework/diagnostics/f158_odd_word.py, max_power 4, max_jumps 3), read
    // 2026-09-30: H^0 A0 H^1 A2 H^2 A1 with trace 32000000 i, and nothing at max_power 2. So the two
    // defaults reach this row's break by two exact routes: the Python's by this word, this object's
    // one-letter default by the count.
    [Fact]
    public void The_Triple_Words_Match_The_Python_Original_Where_Their_Powers_Are_Nonzero()
    {
        var e = Row(3, Chain(3), "XYZ", "...");
        var w = e.Word(maxPower: 4, maxJumps: 3);
        Assert.Equal("H^0·A0·H^1·A2·H^2·A1", w!.Word);
        Assert.Equal((BigInteger.Zero, new BigInteger(32_000_000)), (w.TraceRe, w.TraceIm));
        Assert.Null(e.Word(maxPower: 2, maxJumps: 3));
        Assert.Equal(EndCount.Reading.BrokenByWord, e.Verdict(maxPower: 4, maxJumps: 3));
        Assert.Equal(EndCount.Reading.BrokenByCount, e.Verdict());
    }

    [Fact]
    public void A_Rank_Reading_Is_Labelled_As_One()
    {
        // The witness's broken row (.X., X.X) with the word budget cut to H^1: no word fires (the row's
        // first word is H^2 A), the far bound 2 is not below the near lower bound 2 (I and XXX are the
        // dark strings commuting with H: the X fields and every Heisenberg bond accept X on all three
        // sites), so with the lifted kernels switched off only the ranks decide, and the verdict must
        // say so. Switched on, the lifted near kernel holds six exactly checked vectors against a far
        // bound of 2, and the counts decide exactly.
        var e = Row(3, Chain(3), ".X.", "X.X");
        Assert.Equal(new[] { "III", "XXX" }, e.DarkStrings().Where(s => e.CommutesWithH(new List<(string, long)> { (s.ToString(3), 1) }))
            .Select(s => s.ToString(3)).OrderBy(s => s));
        Assert.Equal((2, 0), e.LowerCounts());
        Assert.Equal(EndCount.Reading.BrokenByRank, e.Verdict(maxPower: 1, lift: false));
        Assert.Equal(EndCount.Reading.BrokenByCount, e.Verdict(maxPower: 1));
        Assert.Equal((6, 2), e.LiftedLowerCounts());
        Assert.Equal(EndCount.Reading.BrokenByWord, e.Verdict(maxPower: 2));
        Assert.True(EndCount.IsExact(EndCount.Reading.PalindromeByElement));
        Assert.True(EndCount.IsExact(EndCount.Reading.PalindromeByCount));

        Assert.True(EndCount.IsExact(EndCount.Reading.PalindromeByColour));
        Assert.True(EndCount.IsExact(EndCount.Reading.BrokenByWord));
        Assert.True(EndCount.IsExact(EndCount.Reading.BrokenByCount));
        Assert.False(EndCount.IsExact(EndCount.Reading.PalindromeByRank));
        Assert.False(EndCount.IsExact(EndCount.Reading.BrokenByRank));
    }

    // ---- the lifted kernels ----

    // The rank-only rows of the N = 3 family, decided exactly once the kernels are lifted. The palindromic
    // one's far element is PREDICTED rather than read: the Heisenberg chain is fixed by the mirror 0 <-> 2
    // composed with the pi rotation about (0,1,1), which also carries the field Y on site 0 to the field
    // Z on site 2; that operator is V = (II - XX + YZ + ZY)/2 on sites 0, 2, it commutes with the jump X
    // on site 1, and the middle takes the lit circle colour Y + Z, the rotation's own axis. Expanded,
    // 2 V (x) (Y + Z) is the eight strings below.
    [Fact]
    public void The_Lifted_Kernels_Decide_The_Rank_Rows_Of_The_Three_Site_Family()
    {
        var broken = Row(3, Chain(3), ".X.", "YZY");
        Assert.Equal((2, 1), broken.UpperCounts());
        Assert.Equal((1, 0), broken.LowerCounts());
        Assert.Equal((2, 1), broken.LiftedLowerCounts());
        Assert.Equal(EndCount.Reading.BrokenByRank, broken.Verdict(lift: false));
        Assert.Equal(EndCount.Reading.BrokenByCount, broken.Verdict());

        var pairs = Row(3, Chain(3), ".X.", "Y.Z");
        Assert.Equal(EndCount.Reading.PalindromeByElement, pairs.Verdict());
        var predicted = new List<(string, long)>
        {
            ("IYI", 1), ("IZI", 1), ("XYX", -1), ("XZX", -1), ("YYZ", 1), ("YZZ", 1), ("ZYY", 1), ("ZZY", 1),
        };
        Assert.True(pairs.CheckElement(predicted).Certifies);
        var lifted = pairs.FarElement()!.Value.Element.OrderBy(t => t.Item1).ToList();
        long sign = Math.Sign(lifted[0].Coefficient);
        Assert.Equal(predicted.OrderBy(t => t.Item1), lifted.Select(t => (t.Letters, sign * t.Coefficient)));
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(6)]
    public void On_The_Canonical_Chain_The_Lifted_Kernels_Meet_The_Upper_Bounds(int n)
    {
        var e = Row(n, Chain(n), new string('Z', n), new string('.', n));
        Assert.Equal((n + 1, n + 1), e.LiftedLowerCounts());
        Assert.Equal(e.UpperCounts(), e.LiftedLowerCounts());
        Assert.All(e.LiftedKernel(far: false).Concat(e.LiftedKernel(far: true)), v => Assert.True(e.CommutesWithH(v)));
    }

    // ---- the family, swept ----

    static IEnumerable<EndCount> Family(int n, (int, int)[] edges, bool everySiteDephased)
    {
        const string A = ".XYZ";
        string Pattern(int code) => new(Enumerable.Range(0, n).Select(l => A[(code >> (2 * l)) & 3]).ToArray());
        for (int dc = 1; dc < 1 << (2 * n); dc++)
        {
            string deph = Pattern(dc);
            if (everySiteDephased && deph.Contains('.')) continue;
            for (int fc = 0; fc < 1 << (2 * n); fc++)
                yield return Row(n, edges, deph, Pattern(fc), bond: 10, mag: 3);
        }
    }

    // Derived, not measured: with every site dephased along one letter, a Heisenberg bond forces one
    // colour c on the connected graph; c must avoid every site's jump letter (2^N jump patterns) and
    // every field must be absent or c (2^N field patterns), so 4^N rows per colour; two colours work
    // together only on the one pattern that dephases every site along the third letter with no field.
    // Inclusion-exclusion gives 3 (4^N - 1) coloured rows, and the graph does not enter.
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void Every_Site_Dephased_The_Coloured_Rows_Are_Three_Times_Four_To_The_N_Minus_One(int n)
    {
        foreach (var edges in new[] { Chain(n), Ring(n) })
        {
            int coloured = Family(n, edges, everySiteDephased: true).Count(e => e.LowerCounts().Far > 0);
            Assert.Equal(3 * ((1 << (2 * n)) - 1), coloured);
        }
    }

    // Measured with the run mode `endcount sweep 3 chain` (2026-09-30): once the kernels are lifted no
    // row of the three-site family is left to the ranks: 603 coloured, 6 by a lifted element, 252 by the
    // counts, 3171 by a one-jump word, of 4032.
    [Fact]
    public void No_Row_Of_The_Three_Site_Family_Is_Left_To_The_Ranks()
    {
        var readings = Family(3, Chain(3), everySiteDephased: false).Select(e => e.Verdict()).ToList();
        Assert.Equal(4032, readings.Count);
        Assert.All(readings, r => Assert.True(EndCount.IsExact(r), r.ToString()));
        Assert.Equal(603, readings.Count(r => r == EndCount.Reading.PalindromeByColour));
        Assert.Equal(6, readings.Count(r => r == EndCount.Reading.PalindromeByElement));
        Assert.Equal(252, readings.Count(r => r == EndCount.Reading.BrokenByCount));
        Assert.Equal(3171, readings.Count(r => r == EndCount.Reading.BrokenByWord));
    }

    // ---- F138 (b): the coincident-magnitude exceptions, one formula ----

    // F138's own grid (docs/ANALYTICAL_FORMULAS.md, F138 (b)): the full bond XX + YY + ZZ on three sites,
    // per site one dephasing letter or none and a field along +-X, +-Y, +-Z or none, 28^3 rows, with the
    // two END magnitudes equal. F138 counts 78 rows that pair although its clauses forbid them, and names
    // U = SWAP02 Z0 Z1 Z2 on one of the six whose field is anti-invariant under the end swap; the tests below
    // show that all 78 carry the same form with the axis Z replaced by an axis n the row determines up to
    // sign (on the six, any axis in a plane).
    static IEnumerable<(string Deph, string Field, EndCount E)> F138Grid(long m0, long m1, long m2, bool complete = false)
    {
        var edges = complete ? new[] { (0, 1), (1, 2), (0, 2) } : Chain(3);
        var mag = new[] { m0, m1, m2 };
        var fields = new (char L, int S)[] { ('.', 0), ('X', 1), ('X', -1), ('Y', 1), ('Y', -1), ('Z', 1), ('Z', -1) };
        foreach (var d0 in ".XYZ") foreach (var d1 in ".XYZ") foreach (var d2 in ".XYZ")
        {
            string deph = $"{d0}{d1}{d2}";
            var jumps = Enumerable.Range(0, 3).Where(l => deph[l] != '.').Select(l => One(3, l, deph[l])).ToList();
            if (jumps.Count == 0) continue;
            foreach (var f0 in fields) foreach (var f1 in fields) foreach (var f2 in fields)
            {
                var fs = new[] { f0, f1, f2 };
                var h = new List<(string, long)>();
                foreach (var (i, j) in edges) foreach (char p in "XYZ") h.Add((Two(3, i, j, p), 100));
                for (int l = 0; l < 3; l++) if (fs[l].L != '.') h.Add((One(3, l, fs[l].L), fs[l].S * mag[l]));
                string field = string.Concat(fs.Select(f => f.L == '.' ? "." : (f.S > 0 ? "+" : "-") + f.L));
                yield return (deph, field, new EndCount(W, 3, h, jumps));
            }
        }
    }

    // A product of combinations of Pauli strings with Gaussian-integer coefficients (re, im).
    static Dictionary<PauliString, (long Re, long Im)> Times(Dictionary<PauliString, (long Re, long Im)> a, Dictionary<PauliString, (long Re, long Im)> b)
    {
        var o = new Dictionary<PauliString, (long Re, long Im)>();
        foreach (var (s, x) in a) foreach (var (t, y) in b)
        {
            var (u, k) = PauliString.Multiply(s, t);
            long re = x.Re * y.Re - x.Im * y.Im, im = x.Re * y.Im + x.Im * y.Re;
            for (int j = 0; j < k; j++) (re, im) = (-im, re);
            var old = o.TryGetValue(u, out var v) ? v : (0, 0);
            o[u] = (old.Item1 + re, old.Item2 + im);
        }
        return o.Where(kv => kv.Value != (0, 0)).ToDictionary(kv => kv.Key, kv => kv.Value);
    }

    // The candidate U = SWAP_ij (n.sigma) (x) (n.sigma) (x) (n.sigma), for an axis n with integer
    // components, scaled to integers: 2 SWAP_ij = III + X_iX_j + Y_iY_j + Z_iZ_j. The two switches build the
    // controls: without the swap, or with the axis factor left off the site the swap fixes.
    static List<(string, long)> SwapTimesAxis((long X, long Y, long Z) n, int i = 0, int j = 2,
        bool swap = true, bool fixedSiteFactor = true)
    {
        var sw = new Dictionary<PauliString, (long, long)> { [PauliString.Parse("III")] = (1, 0) };
        if (swap)
            foreach (char c in "XYZ") sw[PauliString.Parse(Two(3, i, j, c))] = (1, 0);
        var axis = new Dictionary<PauliString, (long, long)> { [PauliString.Parse("III")] = (1, 0) };
        int k = 3 - i - j;
        for (int site = 0; site < 3; site++)
        {
            if (site == k && !fixedSiteFactor) continue;
            var one = new Dictionary<PauliString, (long, long)>();
            if (n.X != 0) one[PauliString.Parse(One(3, site, 'X'))] = (n.X, 0);
            if (n.Y != 0) one[PauliString.Parse(One(3, site, 'Y'))] = (n.Y, 0);
            if (n.Z != 0) one[PauliString.Parse(One(3, site, 'Z'))] = (n.Z, 0);
            axis = Times(axis, one);
        }
        var u = Times(sw, axis);
        // SWAP_ij commutes with (n.sigma) on every site, so the product is Hermitian: no imaginary part survives
        Assert.All(u.Values, c => Assert.Equal(0, c.Im));
        return u.Select(kv => (kv.Key.ToString(3), kv.Value.Re)).ToList();
    }

    static readonly List<(long X, long Y, long Z)> NineAxes = BuildNineAxes();

    static List<(long, long, long)> BuildNineAxes()
    {
        var axes = new List<(long, long, long)> { (1, 0, 0), (0, 1, 0), (0, 0, 1) };
        foreach (var (u, v) in new[] { (0, 1), (0, 2), (1, 2) })
            foreach (int sg in new[] { 1, -1 })
            {
                var a = new long[3]; a[u] = 1; a[v] = sg;
                axes.Add((a[0], a[1], a[2]));
            }
        return axes;
    }

    // The field of a row's site as a vector ('.' is zero), read from the "+X.-Z" notation.
    static long[] FieldVector(string field, int site)
    {
        var toks = System.Text.RegularExpressions.Regex.Matches(field, @"\.|[+-][XYZ]").Select(m => m.Value).ToArray();
        var v = new long[3];
        if (toks[site] != ".") v["XYZ".IndexOf(toks[site][1])] = toks[site][0] == '+' ? 1 : -1;
        return v;
    }

    // The mechanism, as a rule on the row: the pi rotation about n, R(v) = 2 (n.v) n - |n|^2 v (scaled by
    // |n|^2), must carry the field of site i onto the field of site j (equal magnitudes, so the unit
    // vectors suffice), fix the field of the site k the swap fixes, and send k's jump letter to minus
    // itself, i.e. n orthogonal to it. No jump may sit on i or j.
    static bool RulePredicts((long X, long Y, long Z) n, string deph, string field, int i, int j)
    {
        int k = 3 - i - j;
        if (deph[i] != '.' || deph[j] != '.' || deph[k] == '.') return false;
        var nv = new[] { n.X, n.Y, n.Z };
        long nn = nv.Sum(x => x * x);
        long[] Rot(long[] v) { long d = nv.Zip(v, (a, b) => a * b).Sum(); return nv.Zip(v, (a, b) => 2 * d * a - nn * b).ToArray(); }
        bool Eq(long[] a, long[] b) => a.Zip(b, (x, y) => x == nn * y).All(t => t);
        var jump = new long[3]; jump["XYZ".IndexOf(deph[k])] = 1;
        return nv.Zip(jump, (a, b) => a * b).Sum() == 0
            && Eq(Rot(FieldVector(field, i)), FieldVector(field, j))
            && Eq(Rot(FieldVector(field, k)), FieldVector(field, k));
    }

    [Fact]
    public void Every_Coincident_Magnitude_Exception_Of_F138_Is_A_Swap_Times_One_Axis_Cubed()
    {
        var exceptions = F138Grid(30, 22, 30)
            .Where(r => r.E.Colourings().Count == 0 && EndCount.IsPalindrome(r.E.Verdict())).ToList();
        Assert.Equal(78, exceptions.Count);
        Assert.Equal(72, exceptions.Count(r => r.E.Verdict() == EndCount.Reading.PalindromeByElement));
        Assert.Equal(6, exceptions.Count(r => r.E.Verdict() == EndCount.Reading.PalindromeByCount));
        foreach (var (deph, field, e) in exceptions)
        {
            // the axes that certify are exactly the axes the rule predicts, among the nine candidates
            var certify = NineAxes.Where(n => e.CheckElement(SwapTimesAxis(n)).Certifies).ToList();
            var predicted = NineAxes.Where(n => RulePredicts(n, deph, field, 0, 2)).ToList();
            Assert.NotEmpty(certify);
            Assert.Equal(predicted, certify);
            // and both factors are load-bearing: without the swap, or without the axis on the middle
            // site, no candidate certifies
            Assert.DoesNotContain(NineAxes, n => e.CheckElement(SwapTimesAxis(n, swap: false)).Certifies);
            Assert.DoesNotContain(NineAxes, n => e.CheckElement(SwapTimesAxis(n, fixedSiteFactor: false)).Certifies);
        }
        // F138's named row: jump X on site 1, fields +X and -X on the ends, n = Z: SWAP02 Z0 Z1 Z2; on the
        // six, n runs over the plane orthogonal to the one letter the jump and the end fields share
        var six = exceptions.Single(r => r.Deph == ".X." && r.Field == "+X.-X").E;
        Assert.True(six.CheckElement(SwapTimesAxis((0, 0, 1))).Certifies);
        Assert.True(six.CheckElement(SwapTimesAxis((0, 1, 2))).Certifies);
        // the committed magnitudes leave no exception at all
        Assert.Equal(0, F138Grid(30, 22, 41).Count(r => r.E.Colourings().Count == 0 && EndCount.IsPalindrome(r.E.Verdict())));
    }

    // K3, F138's other graph: 78 per equal pair, and 234 when all three magnitudes coincide; each carries
    // the same form with the swap of the pair whose fields the rotation exchanges, the third site dephased.
    [Fact]
    public void On_The_Triangle_The_Coincident_Exceptions_Carry_The_Same_Form()
    {
        Assert.Equal(78, F138Grid(30, 22, 30, complete: true).Count(r => r.E.Colourings().Count == 0 && EndCount.IsPalindrome(r.E.Verdict())));
        var all = F138Grid(30, 30, 30, complete: true)
            .Where(r => r.E.Colourings().Count == 0 && EndCount.IsPalindrome(r.E.Verdict())).ToList();
        Assert.Equal(234, all.Count);
        var pairs = new[] { (0, 2), (0, 1), (1, 2) };
        foreach (var (deph, field, e) in all)
            Assert.True(pairs.Any(pq => NineAxes.Any(n =>
                RulePredicts(n, deph, field, pq.Item1, pq.Item2) && e.CheckElement(SwapTimesAxis(n, pq.Item1, pq.Item2)).Certifies)),
                $"{deph} {field}");
    }

    // The symmetry group, from below: the Heisenberg chain is fixed by the reversal and by every rotation;
    // Z-dephasing on every site keeps the rotations that send Z to plus or minus Z, the four about Z and
    // the four pi rotations about axes in the XY plane; so 2 x 8 = 16 symmetries, 8 of them moving sites.
    [Fact]
    public void The_Canonical_Three_Site_Chain_Has_Sixteen_Symmetries()
    {
        var syms = Row(3, Chain(3), "ZZZ", "...").Symmetries();
        Assert.Equal(16, syms.Count);
        Assert.Equal(8, syms.Count(g => g.MovesSites));
        Assert.Single(syms, g => g.IsIdentity);
        Assert.Equal(24, EndCount.CubeRotations().Count);
    }

    // ---- the symmetry element ----

    // Derived: among the 24 proper rotations that permute the axes up to sign, the pi rotations are those
    // of trace -1, three about the letter axes and six about the bisectors of two letters, which are the
    // nine candidate axes above.
    [Fact]
    public void The_Half_Turns_Among_The_Cube_Rotations_Are_The_Nine_Axes()
    {
        var axes = EndCount.CubeRotations().Select(r => EndCount.HalfTurnAxis(r.Rotation, r.Sign))
            .Where(a => a is not null).Select(a => a!.Value).ToList();
        Assert.Equal(9, axes.Count);
        var canon = axes.Select(a => a.X < 0 || (a.X == 0 && (a.Y < 0 || (a.Y == 0 && a.Z < 0))) ? (-a.X, -a.Y, -a.Z) : a);
        Assert.Equal(NineAxes.OrderBy(a => a), canon.OrderBy(a => a));
    }

    // The engine finds the symmetry and builds its element itself: on F138's 78 and on K3's 234 every
    // exception is explained by a symmetry element, and every palindrome of the three-site family by a
    // colouring or a symmetry element, none by neither (measured, run modes endcount f138 / sweep).
    [Fact]
    public void Every_Palindrome_Of_The_Three_Site_Grids_Is_A_Colouring_Or_A_Symmetry_Element()
    {
        var p3 = F138Grid(30, 22, 30).Where(r => EndCount.IsPalindrome(r.E.Verdict())).Select(r => r.E.Explanation()).ToList();
        Assert.Equal(2085, p3.Count(x => x == "colouring"));
        Assert.Equal(78, p3.Count(x => x == "symmetry"));
        Assert.DoesNotContain(null, p3);
        var k3 = F138Grid(30, 30, 30, complete: true).Where(r => EndCount.IsPalindrome(r.E.Verdict())).Select(r => r.E.Explanation()).ToList();
        Assert.Equal(234, k3.Count(x => x == "symmetry"));
        Assert.DoesNotContain(null, k3);
        var family = Family(3, Chain(3), everySiteDephased: false).Where(e => EndCount.IsPalindrome(e.Verdict())).Select(e => e.Explanation()).ToList();
        Assert.Equal((603, 6), (family.Count(x => x == "colouring"), family.Count(x => x == "symmetry")));
        Assert.DoesNotContain(null, family);
    }

    // The prediction that came from the understanding: a reflection of the five-site ring fixes one site
    // and swaps two pairs, so its element needs four undephased sites. With the jump X on site 0 and the
    // fields Z on site 2 and Y on site 3, the reflection through site 0 swaps 1 <-> 4 and 2 <-> 3, and the
    // pi rotation about (0,1,1) carries Z onto Y and negates X.
    [Fact]
    public void A_Five_Site_Ring_Row_Is_Carried_By_Its_Reflection()
    {
        var e = Row(5, Ring(5), "X....", "..ZY.", bond: 10, mag: 3);
        Assert.True(EndCount.IsPalindrome(e.Verdict()));
        var found = e.SymmetryElement();
        Assert.NotNull(found);
        Assert.Equal(new[] { 0, 4, 3, 2, 1 }, found!.Value.G.Perm);
        Assert.Equal((0L, 1L, 1L), EndCount.HalfTurnAxis(found.Value.G.Rotation, found.Value.G.Sign) is { } a && a.Y < 0 ? (-a.X, -a.Y, -a.Z) : EndCount.HalfTurnAxis(found.Value.G.Rotation, found.Value.G.Sign)!.Value);
        Assert.Equal("symmetry", e.Explanation());
    }

    // ---- elements found elsewhere: the colouring page ----

    // Stage E: P3, a(XX + YY) and b(XX + YY) on the two bonds, jumps X, Z, Y. No colouring; U = a YYZ +
    // b ZXX lies in W_ and squares to (a^2 + b^2) I; at a = b = 1 it spans W_.
    static EndCount StageE(long a, long b) => new(W, 3,
        new List<(string, long)> { ("XXI", a), ("YYI", a), ("IXX", b), ("IYY", b) },
        new[] { "XII", "IZI", "IIY" });

    [Fact]
    public void Stage_E_Sum_Certifies_The_Palindrome_That_No_Colouring_Reaches()
    {
        var e = StageE(3, 4);
        Assert.Empty(e.Colourings());
        var ok = e.CheckElement(new List<(string, long)> { ("YYZ", 3), ("ZXX", 4) });
        Assert.True(ok.Certifies);
        Assert.Equal(new BigInteger(25), ok.SquareScalar);
        // the sign-flipped sum is lit and squares the same, and does NOT commute: the cancellation
        // between the two commutators is the whole mechanism
        var flipped = e.CheckElement(new List<(string, long)> { ("YYZ", 3), ("ZXX", -4) });
        Assert.True(flipped.AllLit);
        Assert.False(flipped.CommutesWithH);
        Assert.Equal(EndCount.Reading.PalindromeByRank, e.Verdict(maxPower: 2, lift: false));
        // lifted from GF(p), the far end's one vector IS the page's sum, up to sign
        Assert.Equal(EndCount.Reading.PalindromeByElement, e.Verdict(maxPower: 2));
        var lifted = e.FarElement()!.Value.Element.OrderBy(t => t.Letters).ToList();
        long sign = Math.Sign(lifted[0].Coefficient);
        Assert.Equal(new[] { ("YYZ", 3L), ("ZXX", 4L) }, lifted.Select(t => (t.Letters, sign * t.Coefficient)));

        Assert.Equal((1, 1), StageE(1, 1).UpperCounts());
    }

    [Fact]
    public void The_Defect_Cascade_Commutes_And_A_First_Site_X_Field_Ends_It()
    {
        // H = XX + YY + h0 YI + h1 IZ, jump X on the second site; C = ZZ - h0 IY - h0 h1 YZ.
        const long h0 = 2, h1 = 3, eps = 5;
        var h = new List<(string, long)> { ("XX", 1), ("YY", 1), ("YI", h0), ("IZ", h1) };
        var e = new EndCount(W, 2, h, new[] { "IX" });
        var c = e.CheckElement(new List<(string, long)> { ("ZZ", 1), ("IY", -h0), ("YZ", -h0 * h1) });
        Assert.True(c.Certifies);
        Assert.Equal(new BigInteger(1 + h0 * h0 + h0 * h0 * h1 * h1), c.SquareScalar);
        Assert.False(e.CheckElement(new List<(string, long)> { ("ZZ", 1), ("IY", -h0), ("YZ", h0 * h1) }).CommutesWithH);

        // eps X on the first site: Tr(H^2 A) = 8 eps, and W_ is empty
        var broken = new EndCount(W, 2, h.Append(("XI", eps)).ToList(), new[] { "IX" });
        var w = broken.Word();
        Assert.Equal("H^2·A0", w!.Word);
        Assert.Equal(new BigInteger(8 * eps), w.TraceRe);
        Assert.Equal(0, broken.UpperCounts().Far);
        Assert.Equal(EndCount.Reading.BrokenByWord, broken.Verdict());
    }

    [Fact]
    public void The_Swap_Row_Outside_The_Sum_Grammar_Is_Certified_By_Its_Element()
    {
        // P3 Heisenberg, fields +X and -X on the two ends, jump X on the middle site:
        // 2 U = ZZZ - YZY - XZX + IZI with U = SWAP02 Z0Z1Z2, U^2 = I.
        var e = Row(3, Chain(3), ".X.", "X.X", bond: 10, mag: 3);
        var fields = new List<(string, long)>();
        foreach (var (i, j) in Chain(3)) foreach (char p in "XYZ") fields.Add((Two(3, i, j, p), 10));
        fields.Add(("XII", 3));
        fields.Add(("IIX", -3));
        var swap = new EndCount(W, 3, fields, new[] { "IXI" });
        Assert.Empty(swap.Colourings());
        var u = swap.CheckElement(new List<(string, long)> { ("ZZZ", 1), ("YZY", -1), ("XZX", -1), ("IZI", 1) });
        Assert.True(u.Certifies);
        Assert.Equal(new BigInteger(4), u.SquareScalar);
        Assert.True(EndCount.IsPalindrome(swap.Verdict(maxPower: 2)));
        // equal field signs are the witness's broken row (.X., X.X) above: the sign is the whole difference
        Assert.False(EndCount.IsPalindrome(e.Verdict(maxPower: 2)));
    }

    [Fact]
    public void Non_Commuting_Jumps_Go_Through_The_String_Route()
    {
        // The colouring page's Open item: jumps X and Z on site 0 of three, sites 1, 2 undephased,
        // H = 1 (x) diag(1,2,3,4) + X (x) (|0><3| + h.c.) + Z (x) ((|1> + |2>)<3| + h.c.); W_ is spanned by
        // Y (x) diag(1,1,1,-1). H is expanded into Pauli strings through dense matrices, c_P = Tr(P H)/8.
        var m = new Complex[8, 8];
        for (int b = 0; b < 2; b++)
            for (int k = 0; k < 4; k++) m[4 * b + k, 4 * b + k] += k + 1;
        void Add(int r, int c, Complex v) { m[r, c] += v; m[c, r] += v; }
        Add(0, 4 + 3, 1); Add(4 + 0, 3, 1);                        // X (x) |0><3| + h.c.
        foreach (int k in new[] { 1, 2 }) { Add(k, 3, 1); Add(4 + k, 4 + 3, -1); }   // Z (x) (|k><3| + h.c.)
        var h = new List<(string, long)>();
        const string L = "IXYZ";
        for (int s = 0; s < 64; s++)
        {
            string letters = new(new[] { L[s & 3], L[(s >> 2) & 3], L[(s >> 4) & 3] });
            var p = PauliStringTests.Dense(letters);
            Complex tr = 0;
            for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) tr += p[i, j] * m[j, i];
            Assert.Equal(0.0, tr.Imaginary);
            double c8 = tr.Real;                                    // = 8 c_P, an integer here
            Assert.Equal(Math.Round(c8), c8);
            if (c8 != 0) h.Add((letters, (long)c8));
        }
        var e = new EndCount(W, 3, h, new[] { "XII", "ZII" });
        Assert.Empty(e.Colourings());
        Assert.Equal((1, 1), e.UpperCounts());
        var u = e.CheckElement(new List<(string, long)> { ("YII", 1), ("YIZ", 1), ("YZI", 1), ("YZZ", -1) });
        Assert.True(u.Certifies);
        Assert.Equal(new BigInteger(4), u.SquareScalar);
    }

    // ---- the docks ----

    [Fact]
    public void MirrorGroups_R_Is_Right_Multiplication_By_The_Canonical_Element_Of_W()
    {
        // R: rho -> rho * X^N, read off MirrorGroup's own signed permutation; X^N is the canonical chain's
        // colouring, so R is the one-sided map F158's sufficiency builds from that element.
        const int n = 3;
        var xn = PauliString.Parse(new string('X', n));
        const string L = "IXYZ";
        for (int s = 0; s < 64; s++)
        {
            string letters = new(new[] { L[s & 3], L[(s >> 2) & 3], L[(s >> 4) & 3] });
            var (outLetters, phase) = MirrorGroup.R.Apply(letters.ToCharArray());
            var (prod, k) = PauliString.Multiply(PauliString.Parse(letters), xn);
            Assert.Equal(prod.ToString(n), new string(outLetters));
            Assert.Equal(k switch { 0 => Complex.One, 1 => Complex.ImaginaryOne, 2 => -Complex.One, _ => -Complex.ImaginaryOne }, phase);
        }
        Assert.Contains(xn, Row(n, Chain(n), "ZZZ", "...").Colourings());
    }

    // the router's chains: sliding windows of c XZX + XZY + YZX, Z on every site
    static EndCount RouterChain(int n, long c)
    {
        var h = new List<(string, long)>();
        for (int w = 0; w + 2 < n; w++)
        {
            string Win(string t) { var a = Enumerable.Repeat('I', n).ToArray(); for (int k = 0; k < 3; k++) a[w + k] = t[k]; return new string(a); }
            if (c != 0) h.Add((Win("XZX"), c));
            h.Add((Win("XZY"), 1));
            h.Add((Win("YZX"), 1));
        }
        return new EndCount(W, n, h, Enumerable.Range(0, n).Select(l => One(n, l, 'Z')).ToList());
    }

    // a number u + v r in Z[r], r^2 = c r + 1
    static (long U, long V) QMul((long U, long V) a, (long U, long V) b, long c)
        => (a.U * b.U + a.V * b.V, a.U * b.V + a.V * b.U + c * a.V * b.V);

    // G = (x)_l g_l with g_l = alpha X + beta Y, expanded over the 2^N strings of X and Y, split into its
    // rational and its r part: G = G0 + r G1. Since r is irrational at c = 1, 2 (c^2 + 4 = 5, 8),
    // [H, G] = 0 exactly when [H, G0] = 0 and [H, G1] = 0.
    static (List<(string, long)> G0, List<(string, long)> G1) Expand(((long, long) Alpha, (long, long) Beta)[] g, long c)
    {
        int n = g.Length;
        var g0 = new List<(string, long)>();
        var g1 = new List<(string, long)>();
        for (int s = 0; s < (1 << n); s++)
        {
            (long, long) coef = (1, 0);
            var letters = new char[n];
            for (int l = 0; l < n; l++)
            {
                bool y = ((s >> l) & 1) == 1;
                letters[l] = y ? 'Y' : 'X';
                coef = QMul(coef, y ? g[l].Beta : g[l].Alpha, c);
            }
            if (coef.Item1 != 0) g0.Add((new string(letters), coef.Item1));
            if (coef.Item2 != 0) g1.Add((new string(letters), coef.Item2));
        }
        return (g0, g1);
    }

    static readonly ((long, long), (long, long)) FrameA = ((0, 1), (1, 0));    // a = r X + Y
    static readonly ((long, long), (long, long)) FrameB = ((1, 0), (0, -1));   // b = X - r Y

    [Theory]
    [InlineData(3, 1)]
    [InlineData(4, 1)]
    [InlineData(5, 1)]
    [InlineData(4, 2)]
    public void The_Routers_Identity_Column_Lies_In_W_Exactly(int n, long c)
    {
        // The frame is the router's own: its site maps' identity columns are (r, 1), (r, 1), (1, -r),
        // (1, -r), the rhythm [a, a, b, b], read off Router.SiteMaps at the golden and silver means.
        double r = (c + Math.Sqrt(c * c + 4)) / 2;
        var maps = Router.SiteMaps(r);
        var expected = new[] { (r, 1.0), (r, 1.0), (1.0, -r), (1.0, -r) };
        for (int l = 0; l < 4; l++)
            Assert.Equal(expected[l], (maps[l][1, 0].Real, maps[l][2, 0].Real));

        var e = RouterChain(n, c);
        var frame = Enumerable.Range(0, n).Select(l => (l % 4) < 2 ? FrameA : FrameB).ToArray();
        var (g0, g1) = Expand(frame, c);
        Assert.True(e.CommutesWithH(g0));
        Assert.True(e.CommutesWithH(g1));
        var lit = e.LitStrings().Select(s => s.ToString(n)).ToHashSet();
        Assert.All(g0.Concat(g1), t => Assert.Contains(t.Item1, lit));
        Assert.True(e.UpperCounts().Far >= 1);
        Assert.Equal(EndCount.Reading.PalindromeByRank, e.Verdict(maxPower: 2, maxJumps: 1, lift: false));
        // lifted, both ends are met exactly (four vectors each at these (n, c), an independent sympy
        // computation of both nullities in the review round of 2026-09-30), and the count decides
        Assert.Equal((4, 4), e.LiftedLowerCounts());
        Assert.Equal(EndCount.Reading.PalindromeByCount, e.Verdict(maxPower: 2, maxJumps: 1));

        // the rhythms the colouring page names as failing: [a, b, a, b] and [a, a, a, a]
        foreach (var bad in new[] { (Func<int, bool>)(l => l % 2 == 0), l => true })
        {
            var (b0, b1) = Expand(Enumerable.Range(0, n).Select(l => bad(l) ? FrameA : FrameB).ToArray(), c);
            Assert.False(e.CommutesWithH(b0) && e.CommutesWithH(b1));
        }
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Where_c_Is_Nonzero_No_Pauli_Letter_Colours_The_Router_Chain_And_At_c_Zero_Exactly_Four_Do(int n)
    {
        foreach (long c in new long[] { 1, 2 })
            Assert.Empty(RouterChain(n, c).Colourings());
        string Period(char a, char b) => new(Enumerable.Range(0, n).Select(l => l % 2 == 0 ? a : b).ToArray());
        var four = new[] { Period('X', 'X'), Period('X', 'Y'), Period('Y', 'X'), Period('Y', 'Y') }.OrderBy(s => s);
        Assert.Equal(four, RouterChain(n, 0).Colourings().Select(s => s.ToString(n)).OrderBy(s => s));
    }

    [Fact]
    public void An_Element_That_Is_Not_Lit_Or_Does_Not_Square_To_A_Scalar_Certifies_Nothing()
    {
        var e = Row(3, Chain(3), "ZZZ", "...");
        var dark = e.CheckElement(new List<(string, long)> { ("ZZZ", 1) });       // commutes with H, but dark
        Assert.False(dark.AllLit);
        Assert.True(dark.CommutesWithH);
        Assert.False(dark.Certifies);
        // XXX and YYX are both lit and commute with each other, so the square keeps 2 XXX·YYX = -2 ZZI
        var mixed = e.CheckElement(new List<(string, long)> { ("XXX", 1), ("YYX", 1) });
        Assert.True(mixed.AllLit);
        Assert.True(mixed.SquareScalar.IsZero);
        Assert.False(mixed.Certifies);
        Assert.Throws<ArgumentException>(() => e.CheckElement(new List<(string, long)> { ("XXX", 1), ("XXX", -1) }));
    }

    [Fact]
    public void Thirty_Two_Sites_Run_Through_The_Object()
    {
        // X and Z on every site: the only dark string is I, the only lit one Y^32, which commutes with
        // every Heisenberg bond; the far end holds exactly the colouring, the near end the identity.
        const int n = 32;
        var h = Enumerable.Range(0, n - 1).SelectMany(i => "XYZ".Select(p => (Two(n, i, i + 1, p), 1L))).ToList();
        var jumps = Enumerable.Range(0, n).SelectMany(l => new[] { One(n, l, 'X'), One(n, l, 'Z') }).ToList();
        var e = new EndCount(W, n, h, jumps);
        Assert.Equal(new[] { new string('I', n) }, e.DarkStrings().Select(s => s.ToString(n)));
        Assert.Equal(new[] { new string('Y', n) }, e.LitStrings().Select(s => s.ToString(n)));
        Assert.Equal((1, 1), e.UpperCounts());
        Assert.Equal(EndCount.Reading.PalindromeByColour, e.Verdict());
    }

    [Fact]
    public void A_Span_Of_Exactly_Two_To_The_Bound_Is_Accepted()
    {
        var e = new EndCount(W, 9, new List<(string, long)> { ("XXIIIIIII", 1) }, new[] { "ZIIIIIIII", "IZIIIIIII" });
        Assert.Equal(1 << EndCount.MaxSpanBits, e.DarkStrings().Count);
    }

    [Fact]
    public void Coefficients_Are_Bounded_After_They_Are_Summed()
    {
        long big = EndCount.MaxCoefficient;
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new EndCount(W, 2, new List<(string, long)> { ("XX", -big), ("XX", -big) }, new[] { "ZI" }));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new EndCount(W, 2, new List<(string, long)> { ("XX", big), ("XX", big), ("XX", big) }, new[] { "ZI" }));
        // and a sum that returns inside the bound is accepted
        Assert.Equal(8, new EndCount(W, 2, new List<(string, long)> { ("XX", big), ("XX", 5 - big) }, new[] { "ZI" }).DarkStrings().Count);
    }

    // ---- the guards ----

    [Fact]
    public void A_Span_Past_The_Bound_Is_Refused_Rather_Than_Enumerated()
    {
        // one jump on nine sites leaves 2^17 dark strings
        var e = new EndCount(W, 9, new List<(string, long)> { ("XXIIIIIII", 1) }, new[] { "ZIIIIIIII" });
        Assert.Throws<InvalidOperationException>(() => e.DarkStrings());
    }

    [Fact]
    public void Malformed_Inputs_Are_Refused()
    {
        Assert.Throws<ArgumentException>(() => new EndCount(W, 2, new List<(string, long)> { ("XX", 1) }, Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => new EndCount(W, 2, new List<(string, long)> { ("XXX", 1) }, new[] { "ZI" }));
        Assert.Throws<ArgumentException>(() => new EndCount(W, 2, new List<(string, long)> { ("XX", 1) }, new[] { "ZII" }));
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new EndCount(W, 2, new List<(string, long)> { ("XX", long.MaxValue) }, new[] { "ZI" }));
    }

    // ---- the contrast: one- and two-letter bonds at distinct magnitudes ----

    // The anticommuting-sum census of experiments/THE_PALINDROME_AS_A_COLOURING.md ("The census"): bonds
    // of one letter set, bond weight 100, fields at F138's committed magnitudes 30, 22, 41 (one per site,
    // positive), every dephasing and field pattern. The palindromic rows the colouring does not reach are
    // F138's own rejected-palindrome counts on P3, 62 at a one-letter bond and 22 at a two-letter bond.
    static EndCount Letters(int n, (int, int)[] edges, string bondLetters, long[] mags, string deph, string field, long weight = 100)
    {
        var h = new List<(string, long)>();
        foreach (var (i, j) in edges)
            foreach (char p in bondLetters) h.Add((Two(n, i, j, p), weight));
        for (int l = 0; l < n; l++)
            if (field[l] != '.') h.Add((One(n, l, field[l]), mags[l]));
        var jumps = Enumerable.Range(0, n).Where(l => deph[l] != '.').Select(l => One(n, l, deph[l])).ToList();
        return new EndCount(W, n, h, jumps);
    }

    static bool PairwiseAnticommuting(IReadOnlyList<(string Letters, long Coefficient)> element) =>
        element.SelectMany((a, i) => element.Take(i).Select(b => (a, b)))
            .All(x => !PauliString.Commute(PauliString.Parse(x.a.Letters), PauliString.Parse(x.b.Letters)));

    // Distinct magnitudes leave no site-moving symmetry that carries a field onto another, so the
    // symmetry element cannot reach these rows; every one is read exactly and none is a colouring or a
    // site symmetry. At the two-letter bond the far end is one-dimensional on every such row and its
    // lifted element is an anticommuting sum (the census: "all 22 of its two-letter rows" are one sum),
    // named "clifford" where the sum is a Clifford unitary up to scale (eight rows, each two strings of
    // equal coefficient magnitude) and "sum" otherwise (fourteen); at the one-letter bond no Clifford
    // element exists (the facts further down) and the grammar names the census's kinds, 52 sums and 10
    // conditioned (its stage A).
    // At N = 3 one prime lifts every row here; the wider lift is gated by the two facts below.
    [Theory]
    [InlineData("ZZ", 62)]
    [InlineData("XY", 22)]
    public void Beyond_The_Colouring_At_Distinct_Magnitudes_No_Palindrome_Is_A_Symmetry(string bondLetters, int beyond)
    {
        var kinds = new Dictionary<string, int>();
        const string alphabet = ".XYZ";
        long[] mags = { 30, 22, 41 };
        int count = 0;
        for (int dc = 1; dc < 64; dc++)
            for (int fc = 0; fc < 64; fc++)
            {
                string deph = new(Enumerable.Range(0, 3).Select(l => alphabet[(dc >> (2 * l)) & 3]).ToArray());
                string field = new(Enumerable.Range(0, 3).Select(l => alphabet[(fc >> (2 * l)) & 3]).ToArray());
                var e = Letters(3, Chain(3), bondLetters, mags, deph, field);
                var r = e.Verdict();
                Assert.True(EndCount.IsExact(r), $"{deph} {field}: {r}");
                if (!EndCount.IsPalindrome(r) || r == EndCount.Reading.PalindromeByColour) continue;
                count++;
                kinds[e.Explanation() ?? "none"] = kinds.TryGetValue(e.Explanation() ?? "none", out int kc) ? kc + 1 : 1;
                if (bondLetters == "XY")
                {
                    Assert.Equal(EndCount.Reading.PalindromeByElement, r);
                    Assert.True(PairwiseAnticommuting(e.FarElement()!.Value.Element), $"{deph} {field}");
                }
            }
        Assert.Equal(beyond, count);
        var expectedKinds = bondLetters == "ZZ"
            ? new Dictionary<string, int> { ["sum"] = 52, ["conditioned"] = 10 }
            : new Dictionary<string, int> { ["clifford"] = 8, ["sum"] = 14 };
        Assert.Equal(expectedKinds, kinds);
    }

    // A sum whose coefficients one prime cannot reconstruct: the census's N = 4 chain at XX + YY, fields
    // 30, 22, 41, 17, jump X on site 0 and fields Z, Z, Z, Y. With the first prime alone no far vector
    // survives (its entries exceed the sqrt(p/2) a single-prime reconstruction returns, 4,670,000 the
    // largest), and the row would be left to the ranks; through the Chinese remainder theorem over
    // further primes the element is lifted and certifies exactly.
    [Fact]
    public void A_Sum_Past_One_Primes_Bound_Is_Lifted_Through_Further_Primes()
    {
        var e = Letters(4, Chain(4), "XY", new long[] { 30, 22, 41, 17 }, "X...", "ZZZY");
        Assert.Empty(e.LiftedKernel(far: true, primes: 1));
        Assert.Single(e.LiftedKernel(far: true));
        Assert.Equal(EndCount.Reading.PalindromeByElement, e.Verdict());
        var (element, check) = e.FarElement()!.Value;
        Assert.True(check.Certifies);
        Assert.True(PairwiseAnticommuting(element));
        long bound = (long)Math.Sqrt(EndCount.LiftPrimes[0] / 2.0);
        Assert.True(element.Max(t => Math.Abs(t.Coefficient)) > bound);
        Assert.Equal("sum", e.Explanation());
    }
    // A bad first prime: bonds of weight exactly LiftPrimes[0] vanish mod that prime. At Heisenberg
    // bonds its rank drops on some rows; at a ZZ bond it can keep the full rank and only move a pivot
    // to a later column. Either way a later prime takes over as the anchor, and every row of the
    // three-site chain is read exactly, as it is at any other weight; at ZZ the palindromes beyond the
    // colouring are the 62 of the contrast above.
    [Theory]
    [InlineData("XYZ", -1)]
    [InlineData("ZZ", 62)]
    public void A_Bad_First_Prime_Hands_The_Anchor_To_A_Good_One(string bondLetters, int beyond)
    {
        const string alphabet = ".XYZ";
        long p0 = EndCount.LiftPrimes[0];
        long[] mags = bondLetters == "XYZ" ? new long[] { 3, 3, 3 } : new long[] { 30, 22, 41 };
        int count = 0;
        for (int dc = 1; dc < 64; dc++)
            for (int fc = 0; fc < 64; fc++)
            {
                string deph = new(Enumerable.Range(0, 3).Select(l => alphabet[(dc >> (2 * l)) & 3]).ToArray());
                string field = new(Enumerable.Range(0, 3).Select(l => alphabet[(fc >> (2 * l)) & 3]).ToArray());
                var r = Letters(3, Chain(3), bondLetters, mags, deph, field, weight: p0).Verdict();
                Assert.True(EndCount.IsExact(r), $"{deph} {field}: {r}");
                if (EndCount.IsPalindrome(r) && r != EndCount.Reading.PalindromeByColour) count++;
            }
        if (beyond >= 0) Assert.Equal(beyond, count);
    }
    // The same row at field magnitudes 3001, 2203, 4111, 1709: its sum's rationals need four primes,
    // so the lift runs three CRT steps, two of them on a composite modulus. Empty at three primes, the
    // one certifying sum at four.
    [Fact]
    public void A_Sum_Needing_Four_Primes_Is_Lifted_Through_A_Composite_Modulus()
    {
        var e = Letters(4, Chain(4), "XY", new long[] { 3001, 2203, 4111, 1709 }, "X...", "ZZZY");
        Assert.Empty(e.LiftedKernel(far: true, primes: 3));
        Assert.Single(e.LiftedKernel(far: true, primes: 4));
        Assert.Equal(EndCount.Reading.PalindromeByElement, e.Verdict());
        Assert.True(e.FarElement()!.Value.Check.Certifies);
        Assert.True(PairwiseAnticommuting(e.FarElement()!.Value.Element));
    }
    // ---- the Clifford elements of the far space ----

    static IEnumerable<(string Deph, string Field)> ThreeSitePatterns()
    {
        const string alphabet = ".XYZ";
        for (int dc = 1; dc < 64; dc++)
            for (int fc = 0; fc < 64; fc++)
                yield return (new string(Enumerable.Range(0, 3).Select(l => alphabet[(dc >> (2 * l)) & 3]).ToArray()),
                              new string(Enumerable.Range(0, 3).Select(l => alphabet[(fc >> (2 * l)) & 3]).ToArray()));
    }

    // Soundness and nesting against the exact verdict, over every row of the three-site chain for the
    // three bond sets at one and at distinct magnitudes: a Clifford map the search returns is read
    // palindromic (a map on a broken row would contradict F158, so a search that accepted a map it
    // should not, a dropped commutation or phase check, fails here), the search never runs out of
    // budget, and every row with a colouring or a site symmetry's element has a Clifford map too.
    [Theory]
    [InlineData("XYZ", 30, 30, 30)]
    [InlineData("XYZ", 30, 22, 41)]
    [InlineData("XY", 30, 30, 30)]
    [InlineData("XY", 30, 22, 41)]
    [InlineData("ZZ", 30, 30, 30)]
    [InlineData("ZZ", 30, 22, 41)]
    public void Every_Clifford_Map_The_Search_Returns_Carries_A_Palindrome(string bondLetters, long m0, long m1, long m2)
    {
        int maps = 0;
        foreach (var (deph, field) in ThreeSitePatterns())
        {
            var e = Letters(3, Chain(3), bondLetters, new[] { m0, m1, m2 }, deph, field);
            var (map, exhausted) = e.CliffordSymmetry();
            Assert.True(exhausted, $"{deph} {field}: the search ran out of budget");
            var r = e.Verdict();
            if (map is not null) { maps++; Assert.True(EndCount.IsPalindrome(r), $"{deph} {field}: a Clifford map on a row read {r}"); }
            if (e.Colourings().Count > 0 || e.SymmetryElement() is not null)
                Assert.True(map is not null, $"{deph} {field}: coloured or symmetric, but no Clifford map");
        }
        Assert.True(maps > 0);
    }

    // The fourth carrier. On P3 with XX + YY bonds, fields Y on site 1 and Z on site 2 of one magnitude
    // and the jump X on site 0, the far end is one-dimensional and its element
    //   E = -YII - ZXX - YYZ + ZZY
    // is (twice) a Clifford unitary that is no site symmetry: conjugation by it sends the bond XX_01 to
    // YY_12 and YY_12 back, fixes YY_01 and XX_12, swaps the two fields and negates the jump, so the
    // magnitudes must coincide (at distinct ones the row is broken). The search finds exactly that map,
    // and E conjugates every term of H onto its image, checked here by exact string products.
    [Fact]
    public void A_Clifford_That_Permutes_The_Bonds_Carries_A_Row_No_Site_Symmetry_Reaches()
    {
        var e = Letters(3, Chain(3), "XY", new long[] { 30, 30, 30 }, "X..", ".YZ");
        Assert.Equal(EndCount.Reading.PalindromeByElement, e.Verdict());
        Assert.Null(e.SymmetryElement());
        Assert.Equal("clifford", e.Explanation());
        var (map, _) = e.CliffordSymmetry();
        Assert.NotNull(map);
        // the map read through the terms' own strings, whatever order the object keeps them in
        var names = e.Terms.Select(x => x.Letters).ToArray();
        var byName = Enumerable.Range(0, names.Length).ToDictionary(k => names[k], k => names[map!.Image[k]]);
        Assert.Equal(new Dictionary<string, string>
        {
            ["XXI"] = "IYY", ["IYY"] = "XXI", ["YYI"] = "YYI", ["IXX"] = "IXX", ["IYI"] = "IIZ", ["IIZ"] = "IYI",
        }, byName);
        Assert.All(map!.Sign, s => Assert.Equal(1, s));

        var (element, check) = e.FarElement()!.Value;
        Assert.True(check.Certifies);
        var expected = new Dictionary<string, long> { ["YII"] = -1, ["ZXX"] = -1, ["YYZ"] = -1, ["ZZY"] = 1 };
        long scale = element.First(x => x.Letters == "YII").Coefficient / -1;
        Assert.Equal(expected.Count, element.Count);
        Assert.All(element, x => Assert.Equal(expected[x.Letters] * scale, x.Coefficient));
        string[] termStrings = { "XXI", "YYI", "IXX", "IYY", "IYI", "IIZ" };
        foreach (var t in termStrings)
            Assert.Equal(new Dictionary<string, (long, long)> { [byName[t]] = (4, 0) }, Conjugate(expected, t));
        Assert.Equal(new Dictionary<string, (long, long)> { ["XII"] = (-4, 0) }, Conjugate(expected, "XII"));

        var broken = Letters(3, Chain(3), "XY", new long[] { 30, 22, 41 }, "X..", ".YZ");
        Assert.False(EndCount.IsPalindrome(broken.Verdict()));
        Assert.Equal((null, true), (broken.CliffordSymmetry().Map, broken.CliffordSymmetry().Exhausted));
    }

    // At a ZZ bond the palindromes beyond the colouring have no Clifford element at all, at one
    // magnitude and at distinct ones: whatever carries them has coefficients that are not those of a
    // Clifford unitary.
    [Theory]
    [InlineData(30, 30, 30)]
    [InlineData(30, 22, 41)]
    public void At_A_ZZ_Bond_No_Palindrome_Beyond_The_Colouring_Has_A_Clifford_Element(long m0, long m1, long m2)
    {
        int beyond = 0;
        foreach (var (deph, field) in ThreeSitePatterns())
        {
            var e = Letters(3, Chain(3), "ZZ", new[] { m0, m1, m2 }, deph, field);
            var r = e.Verdict();
            if (!EndCount.IsPalindrome(r) || r == EndCount.Reading.PalindromeByColour) continue;
            beyond++;
            Assert.Null(e.CliffordSymmetry().Map);
            Assert.Contains(e.Explanation(), new[] { "sum", "conditioned" });
        }
        Assert.Equal(62, beyond);
    }

    // ---- the anticommuting-sum grammar ----

    // The census's stage A at a ZZ bond, fields 30, 22, 41 against bonds 100: beyond the colouring P3
    // holds 52 single sums and 10 conditioned elements, K3 18 conditioned ones
    // (simulations/results/anticommuting_sum_census.txt; the kinds agree row for row with
    // anticommuting_sum_certificates.json, read once when the grammar was ported). GrammarElement
    // checks what it builds (lit, commuting with H, invertible) and throws otherwise, so every row
    // here is also a check of the construction; on every broken row it finds nothing, the census's
    // control, since an invertible element of the far space would be a palindrome.
    // The bond with an isolated site is the census's product row: H splits into the bond and the lone
    // site, and every element beyond the colouring is a tensor product, 52 at ZZ and 104 at XX + YY.
    [Theory]
    [InlineData("chain", "ZZ", 52, 10, 0)]
    [InlineData("complete", "ZZ", 0, 18, 0)]
    [InlineData("bond+iso", "ZZ", 0, 0, 52)]
    [InlineData("bond+iso", "XY", 0, 0, 104)]
    public void The_Grammar_Names_The_Census_Kinds(string topology, string bondLetters, int sums, int conditioned, int products)
    {
        var edges = topology switch
        {
            "chain" => Chain(3),
            "complete" => new[] { (0, 1), (1, 2), (0, 2) },
            _ => new[] { (0, 1) },
        };
        var kinds = new Dictionary<string, int>();
        foreach (var (deph, field) in ThreeSitePatterns())
        {
            var e = Letters(3, edges, bondLetters, new long[] { 30, 22, 41 }, deph, field);
            var r = e.Verdict();
            var g = e.GrammarElement();
            if (!EndCount.IsPalindrome(r)) { Assert.Null(g); continue; }
            if (r == EndCount.Reading.PalindromeByColour) continue;
            Assert.NotNull(g);
            kinds[g!.Value.Kind] = kinds.TryGetValue(g.Value.Kind, out int k) ? k + 1 : 1;
        }
        var expected = new Dictionary<string, int>();
        if (conditioned > 0) expected["conditioned"] = conditioned;
        if (sums > 0) expected["sum"] = sums;
        if (products > 0) expected["product"] = products;
        Assert.Equal(expected, kinds);
    }

    // Invertible against matrices whose determinant is known: 1 + Z is singular (eigenvalues 2 and 0),
    // X + Z is invertible (it squares to 2), a sector sum (1 + Z)/2 ⊗ X + (1 − Z)/2 ⊗ 3·Y (twice it:
    // IX + ZX + 3·IY − 3·ZY) is invertible, and with the second sector's element set to zero
    // (IX + ZX) it is singular.
    // Why the turn to Z must be proper, for the counts if not for the verdicts: under Z jumps on two
    // sites, H = XI + XZ + IX + ZX + XX − YY has a connected hopping graph (00 joined to 01, 10 and
    // 11), counts (1, 0); its partial transpose on site 0, the same with + YY, moves the edge 00–11 to
    // 01–10 and has counts (2, 0). Both are broken, so the verdict survives; the near count does not.
    [Fact]
    public void A_Partial_Transpose_Keeps_This_Verdict_And_Moves_The_Near_Count()
    {
        string[] jumps = { "ZI", "IZ" };
        var common = new List<(string, long)> { ("XI", 1), ("XZ", 1), ("IX", 1), ("ZX", 1), ("XX", 1) };
        var h = new EndCount(W, 2, common.Append(("YY", -1)).ToList(), jumps);
        var t = new EndCount(W, 2, common.Append(("YY", 1)).ToList(), jumps);
        Assert.Equal((1, 0), (h.ComplementConnection()!.Value.HoppingComponents, h.ComplementConnection()!.Value.Good));
        Assert.Equal((2, 0), (t.ComplementConnection()!.Value.HoppingComponents, t.ComplementConnection()!.Value.Good));
        Assert.Equal((1, 0), h.UpperCounts());
        Assert.Equal((2, 0), t.UpperCounts());
    }

    // A flat section of modulus one need not meet the frame formula d_x = v̄_x̄·v_x: at N = 2, Heisenberg
    // bond, Z on both sites, D = diag(1, i, i, 1) on 00, 01, 10, 11 is constant on the popcount shells
    // (flat), X^⊗2·D lies in the far space, yet d_x·d_x̄ = −1 on {01, 10}, which the formula would force
    // to be 1. In strings, diag(1, i, i, 1) = ((1+i)/2)·II + ((1−i)/2)·ZZ and XX·ZZ = −YY, so
    // 2·X^⊗2·D = (XX − YY) + i·(XX + YY): its real and imaginary parts are each lit and in the far space.
    [Fact]
    public void A_Unitary_Flat_Section_Need_Not_Be_Of_Frame_Form()
    {
        var e = new EndCount(W, 2, new List<(string, long)> { ("XX", 1), ("YY", 1), ("ZZ", 1) }, new[] { "ZI", "IZ" });
        foreach (var part in new[] { new[] { ("XX", 1L), ("YY", -1L) }, new[] { ("XX", 1L), ("YY", 1L) } })
            Assert.True(e.CheckElement(part) is { AllLit: true, CommutesWithH: true });
        // the far end is the three shells, every one good, so each shell carries its own free constant
        Assert.Equal((3, 3, 3, true), e.ComplementConnection()!.Value);
        // d on 00, 01, 10, 11 and its complement products: 1·1 on {00, 11}, i·i = −1 on {01, 10}
        var d = new[] { (1, 0), (0, 1), (0, 1), (1, 0) };
        (int, int) Times((int a, int b) u, (int c, int e2) v) => (u.a * v.c - u.b * v.e2, u.a * v.e2 + u.b * v.c);
        Assert.Equal((1, 0), Times(d[0], d[3]));
        Assert.Equal((-1, 0), Times(d[1], d[2]));
    }

    // Theorem 3 of the complement-connection proof at random weights, with Theorem 1's one-axis rows among them: Heisenberg bonds of random nonzero weights on the four graphs,
    // every site dephased along a random letter, random letter fields of random signed magnitudes.
    // The palindrome holds exactly when some letter c is no jump axis and every field is c, the
    // colouring rule for a Heisenberg component; on every broken row the far end is zero, not merely
    // short of the near end; and where the axes really are mixed (two letters or more) the far end
    // is at most one-dimensional, the near end exactly one (Lemma A). About a third of the rows force one field letter, and some rows
    // have one common axis (Theorem 1's class).
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void Mixed_Axes_Heisenberg_Pairs_Exactly_When_A_Colouring_Exists(int n)
    {
        var rng = new Random(100 + n);
        var graphs = new[]
        {
            Chain(n), Ring(n), Enumerable.Range(1, n - 1).Select(i => (0, i)).ToArray(),
            (from a in Enumerable.Range(0, n) from b in Enumerable.Range(a + 1, n - a - 1) select (a, b)).ToArray(),
        };
        int pal = 0, broken = 0;
        foreach (var edges in graphs)
            for (int trial = 0; trial < 300; trial++)
            {
                var axes = Enumerable.Range(0, n).Select(_ => "XYZ"[rng.Next(3)]).ToArray();
                var fields = Enumerable.Range(0, n).Select(_ => ".XYZ"[rng.Next(4)]).ToArray();
                if (rng.Next(3) == 0) { char c0 = "XYZ"[rng.Next(3)]; fields = fields.Select(f => f == '.' ? '.' : c0).ToArray(); }
                var h = new List<(string, long)>();
                foreach (var (a, b) in edges)
                {
                    long w = rng.Next(1, 200) * (rng.Next(2) == 0 ? 1 : -1);
                    foreach (char p in "XYZ") h.Add((Two(n, a, b, p), w));
                }
                for (int l = 0; l < n; l++)
                    if (fields[l] != '.') h.Add((One(n, l, fields[l]), rng.Next(1, 90) * (rng.Next(2) == 0 ? 1 : -1)));
                var e = new EndCount(W, n, h, Enumerable.Range(0, n).Select(l => One(n, l, axes[l])).ToList());
                var r = e.Verdict();
                Assert.True(EndCount.IsExact(r), $"{new string(axes)} {new string(fields)}: {r}");
                bool colouring = "XYZ".Any(c => !axes.Contains(c) && fields.All(f => f == '.' || f == c));
                Assert.Equal(colouring, EndCount.IsPalindrome(r));
                if (colouring) { pal++; Assert.NotEmpty(e.Colourings()); }
                else { broken++; Assert.Equal(0, e.ComplementConnection()!.Value.Good); }
                if (axes.Distinct().Count() > 1)
                    Assert.True(e.ComplementConnection()!.Value is { HoppingComponents: 1, Good: <= 1 }, $"{new string(axes)} {new string(fields)}");
            }
        Assert.True(pal > 50 && broken > 50, $"palindromic {pal}, broken {broken}");
    }

    // Theorem 3 of docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md at the inputs most likely to
    // break it: every coupling and every field component of magnitude one (the coincident magnitudes
    // that break F138's converse where a site is undephased), one field per subset of letters
    // (X + Y, X − Y + Z, ...) with unit components of random sign, every mixed assignment of axes, at N = 2
    // and on the path and the triangle at N = 3 (where three axes need the proof's triangle argument).
    // Lemma A: the hopping graph is connected, the near end one-dimensional. The theorem: the palindrome
    // holds exactly when two axes are used and every field component lies along the third letter c,
    // carried by c^⊗N alone; on every broken row the far end is zero.
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void Mixed_Axes_At_Unit_Magnitudes_Follow_Theorem_Three(int n)
    {
        var rng = new Random(300 + n);
        var graphs = n == 2 ? new[] { Chain(2) } : new[] { Chain(3), Ring(3) };
        int pal = 0, broken = 0;
        foreach (var edges in graphs)
            foreach (int axesCode in Enumerable.Range(0, (int)Math.Pow(3, n)))
            {
                var axes = Enumerable.Range(0, n).Select(l => "XYZ"[axesCode / (int)Math.Pow(3, l) % 3]).ToArray();
                if (axes.Distinct().Count() < 2) continue;
                foreach (int fieldCode in Enumerable.Range(0, 1 << (3 * n)))
                {
                    var h = new List<(string, long)>();
                    foreach (var (a, b) in edges)
                    {
                        long w = rng.Next(2) == 0 ? 1 : -1;
                        foreach (char p in "XYZ") h.Add((Two(n, a, b, p), w));
                    }
                    var fieldLetters = new HashSet<char>();
                    for (int l = 0; l < n; l++)
                        for (int k = 0; k < 3; k++)
                            if ((fieldCode >> (3 * l + k) & 1) == 1)
                            {
                                h.Add((One(n, l, "XYZ"[k]), rng.Next(2) == 0 ? 1 : -1));
                                fieldLetters.Add("XYZ"[k]);
                            }
                    var e = new EndCount(W, n, h, Enumerable.Range(0, n).Select(l => One(n, l, axes[l])).ToList());
                    var r = e.Verdict();
                    string row = $"{new string(axes)} fields {fieldCode} on {edges.Length} bonds";
                    Assert.True(EndCount.IsExact(r), $"{row}: {r}");
                    var cc = e.ComplementConnection()!.Value;
                    Assert.True(cc.HoppingComponents == 1, row);
                    var free = "XYZ".Where(c => !axes.Contains(c)).ToList();
                    bool predicted = free.Count == 1 && fieldLetters.All(f => f == free[0]);
                    Assert.True(predicted == EndCount.IsPalindrome(r), row);
                    if (predicted)
                    {
                        pal++;
                        Assert.Equal(new[] { new string(free[0], n) }, e.Colourings().Select(s => s.ToString(n)));
                    }
                    else { broken++; Assert.True(cc.Good == 0, row); }
                }
            }
        Assert.True(pal > 0 && broken > 0, $"palindromic {pal}, broken {broken}");
    }

    // Theorem 3 needs a connected graph: two bonds (0, 1) and (2, 3), axes X, Y, Y, Z, Z fields on
    // sites 0 and 1, X fields on 2 and 3, use three axes and still pair, each bond a component with
    // its own third letter, carried by ZZXX.
    [Fact]
    public void Three_Axes_Pair_On_A_Disconnected_Graph()
    {
        var h = new List<(string, long)>();
        foreach (var (a, b) in new[] { (0, 1), (2, 3) })
            foreach (char p in "XYZ") h.Add((Two(4, a, b, p), 1));
        h.AddRange(new[] { ("ZIII", 1L), ("IZII", -1L), ("IIXI", 1L), ("IIIX", 1L) });
        var e = new EndCount(W, 4, h, new[] { "XIII", "IYII", "IIYI", "IIIZ" });
        var r = e.Verdict();
        Assert.True(EndCount.IsExact(r) && EndCount.IsPalindrome(r), $"{r}");
        Assert.Contains("ZZXX", e.Colourings().Select(c => c.ToString(4)));
    }

    // One undephased site under Heisenberg bonds, measured: every dephasing pattern that leaves exactly
    // one site without a jump and every letter-field pattern, on the path and the triangle at N = 3 and
    // on the path at N = 4, at distinct field magnitudes and at equal ones. Every verdict is exact, and
    // the palindrome holds exactly when a colouring exists: a letter p that is no dephased site's axis
    // with every field along p. Counted (derived): the undephased site (N ways), the axes avoiding p
    // (2^(N-1)), the fields in {none, p} (2^N), for each of three letters, less the 3 rows per site where
    // two letters both colour (every axis the third letter, no field): N·(3·2^(N-1)·2^N − 3), so 279 at
    // N = 3 and 1524 at N = 4. With two undephased sites the rule fails (F138 (b), the SWAP rows).
    [Theory]
    [InlineData(3, "chain", 279)]
    [InlineData(3, "ring", 279)]
    [InlineData(4, "chain", 1524)]
    public void One_Undephased_Site_Heisenberg_Pairs_Exactly_When_A_Colouring_Exists(int n, string topology, int palindromes)
    {
        const string A = ".XYZ";
        var edges = topology == "chain" ? Chain(n) : Ring(n);
        foreach (var mags in new[] { new long[] { 30, 22, 41, 17 }, new long[] { 30, 30, 30, 30 } })
        {
            int pal = 0;
            for (int dc = 0; dc < 1 << (2 * n); dc++)
            {
                string deph = new(Enumerable.Range(0, n).Select(l => A[(dc >> (2 * l)) & 3]).ToArray());
                if (deph.Count(c => c == '.') != 1) continue;
                for (int fc = 0; fc < 1 << (2 * n); fc++)
                {
                    string field = new(Enumerable.Range(0, n).Select(l => A[(fc >> (2 * l)) & 3]).ToArray());
                    var e = Letters(n, edges, "XYZ", mags, deph, field, weight: 10);
                    var r = e.Verdict();
                    Assert.True(EndCount.IsExact(r), $"{deph} {field}: {r}");
                    bool colouring = "XYZ".Any(p => !deph.Contains(p) && field.All(f => f == '.' || f == p));
                    Assert.True(colouring == EndCount.IsPalindrome(r), $"{deph} {field}: {r}");
                    if (colouring) pal++;
                }
            }
            Assert.Equal(palindromes, pal);
        }
    }

    // Theorem 5's branches that uniform bonds never reach, site 0 undephased: the windmill triangle
    // (0, 1, 2) at its coincidence J01 = J02 = −J12 and off it; the bowtie, triangles (0, 1, 2) and
    // (0, 3, 4) each at its own coincidence (with three letters and no field its far end is asserted
    // empty, 24 assignments at each of three couplings); and a star on 1, 2, 3 plus the bond (1, 2) at
    // unequal couplings to site 0 (the S3 branch). Every axis assignment of the dephased sites, fields
    // from a fixed random draw of letters and signed magnitudes (none on a third of the rows). Each
    // verdict exact, the palindrome exactly where the rule says, and then p^⊗N a colouring.
    [Fact]
    public void One_Undephased_Site_Heisenberg_Holds_At_Coincident_Couplings()
    {
        var rng = new Random(20261002);
        int rows = 0, pal = 0, bowtieEmpty = 0;
        void Check(int n, (int A, int B, long J)[] bonds)
        {
            foreach (int ac in Enumerable.Range(0, (int)Math.Pow(3, n - 1)))
                for (int rep = 0; rep < 3; rep++)
                {
                    var axes = new char[n];
                    axes[0] = '.';
                    for (int l = 1; l < n; l++) axes[l] = "XYZ"[ac / (int)Math.Pow(3, l - 1) % 3];
                    var field = Enumerable.Range(0, n).Select(_ => rep == 0 ? '.' : ".XYZ"[rng.Next(4)]).ToArray();
                    var h = new List<(string, long)>();
                    foreach (var (a, b, j) in bonds)
                        foreach (char c in "XYZ") h.Add((Two(n, a, b, c), j));
                    for (int l = 0; l < n; l++)
                        if (field[l] != '.') h.Add((One(n, l, field[l]), rng.Next(1, 6) * (rng.Next(2) == 0 ? 1 : -1)));
                    var e = new EndCount(W, n, h, Enumerable.Range(1, n - 1).Select(l => One(n, l, axes[l])).ToList());
                    var r = e.Verdict();
                    string row = $"{new string(axes)} {new string(field)} bonds {string.Join(" ", bonds)}";
                    Assert.True(EndCount.IsExact(r), $"{row}: {r}");
                    var free = "XYZ".Where(p => !axes.Contains(p) && field.All(f => f == '.' || f == p)).ToList();
                    Assert.True(free.Count > 0 == EndCount.IsPalindrome(r), $"{row}: {r}");
                    rows++;
                    if (free.Count > 0)
                    {
                        pal++;
                        Assert.Contains(new string(free[0], n), e.Colourings().Select(c => c.ToString(n)));
                    }
                    // the bowtie at its coincidence, three letters, no field: the far end is empty, both
                    // routes (the proof's coupling-free system, read at several couplings)
                    if (n == 5 && rep == 0 && axes[1] != axes[2] && axes[3] != axes[4] && axes.Skip(1).Distinct().Count() == 3)
                    {
                        Assert.Equal(0, e.UpperCounts().Far);
                        bowtieEmpty++;
                    }
                }
        }
        foreach (long a in new long[] { 1, -2, 3 })
        {
            Check(3, new[] { (0, 1, a), (0, 2, a), (1, 2, -a) });
            Check(3, new[] { (0, 1, a), (0, 2, a), (1, 2, a) });
            long b = rng.Next(1, 5) * (rng.Next(2) == 0 ? 1 : -1);
            Check(5, new[] { (0, 1, a), (0, 2, a), (1, 2, -a), (0, 3, b), (0, 4, b), (3, 4, -b) });
            Check(4, new[] { (0, 1, a), (0, 2, 2 * a + 1), (0, 3, a), (1, 2, 5L) });
        }
        Assert.True(pal > 20 && rows - pal > 100, $"rows {rows}, palindromic {pal}");
        Assert.Equal(3 * 24, bowtieEmpty);
    }

    // Theorem 5 needs Heisenberg bonds: on one XX + YY bond with the first site undephased, a jump X on the
    // second and fields Y on the first and Z on the second (the colouring page's defect cascade, whose
    // carrier is an anticommuting sum), the row pairs with no colouring.
    [Fact]
    public void One_Undephased_Site_Without_Heisenberg_Bonds_Pairs_Without_A_Colouring()
    {
        var e = new EndCount(W, 2, new List<(string, long)> { ("XX", 1), ("YY", 1), ("YI", 2), ("IZ", 3) }, new[] { "IX" });
        var r = e.Verdict();
        Assert.True(EndCount.IsExact(r) && EndCount.IsPalindrome(r), $"{r}");
        Assert.Empty(e.Colourings());
    }

    // Two undephased sites under Heisenberg bonds pair with neither a colouring nor a site symmetry: on the
    // triangle with the jump Z on site 0, sites 1 and 2 free, fields X on 1 and Y on 2, couplings J01, J02,
    // J12 with 1/J01 + 1/J02 + 1/J12 = 0 (here 3, 6, −2). The weighted triangle has no automorphism
    // (J01 ≠ J02), the fields share no letter, and the far end holds W below, lit, commuting with H, W² =
    // 125·1; it is the s = 3, t = 6 member (scaled) of the family X₀⊗(st²·II + s³·XY + st²·YX − st²·ZZ) +
    // Y₀⊗(s²t·II + t³·XY + s²t·YX − s²t·ZZ) on sites 1, 2, W² = (s² + t²)³. Off the locus (J12 = −3) it breaks.
    [Fact]
    public void Two_Undephased_Sites_Pair_On_The_Reciprocal_Coupling_Locus_Without_Colouring_Or_Symmetry()
    {
        static EndCount Triangle(long j01, long j02, long j12)
        {
            var h = new List<(string, long)> { ("IXI", 1), ("IIY", 1) };
            foreach (char p in "XYZ")
            {
                h.Add(($"{p}{p}I", j01));
                h.Add(($"{p}I{p}", j02));
                h.Add(($"I{p}{p}", j12));
            }
            return new EndCount(W, 3, h, new[] { "ZII" });
        }
        var e = Triangle(3, 6, -2);
        var r = e.Verdict();
        Assert.True(EndCount.IsExact(r) && EndCount.IsPalindrome(r), $"{r}");
        Assert.Empty(e.Colourings());
        Assert.Null(e.SymmetryElement());
        var w = new List<(string, long)>
        {
            ("XII", -4), ("XXY", -1), ("XYX", -4), ("XZZ", 4), ("YII", -2), ("YXY", -8), ("YYX", -2), ("YZZ", 2),
        };
        var check = e.CheckElement(w);
        Assert.True(check.Certifies);
        Assert.Equal(new System.Numerics.BigInteger(125), check.SquareScalar);
        Assert.Equal(1, e.UpperCounts().Far);                       // one far element, so every carrier is W
        Assert.False(EndCount.IsPalindrome(Triangle(3, 6, -3).Verdict()));
        Assert.True(EndCount.IsPalindrome(Triangle(-1, 2, 2).Verdict()));   // another point of the locus
    }

    // Theorem 6, the local system: with a spanning tree of every component whose blocks of H and of H̄ are
    // invertible and the two supports equal, both ends read at one base point per component from the transported generators equal the end
    // count's own string-span counts, on random Pauli Hamiltonians at N = 2 to 4 with one or two undephased
    // sites and random single-letter jumps on the rest (fixed seed). Both readings are nullities over GF(p), the
    // smaller over the primes; they share no code (blocks and transports here, commutators on strings there).
    [Fact]
    public void Local_System_Counts_Equal_The_End_Count()
    {
        var rng = new Random(20261002 + 7);
        int rows = 0, pal = 0, broken = 0, twoFree = 0;
        for (int draw = 0; draw < 900 && rows < 250; draw++)
        {
            int n = rng.Next(2, 5);
            int nFree = n > 2 && rng.Next(3) == 0 ? 2 : 1;
            var free = Enumerable.Range(0, n).OrderBy(_ => rng.Next()).Take(nFree).ToHashSet();
            var jumps = Enumerable.Range(0, n).Where(l => !free.Contains(l))
                .Select(l => One(n, l, "XYZ"[rng.Next(3)])).ToList();
            var h = new List<(string, long)>();
            int nt = rng.Next(2, 7);
            for (int k = 0; k < nt; k++)
            {
                var s = new string(Enumerable.Range(0, n).Select(_ => "IIXYZ"[rng.Next(5)]).ToArray());
                if (s.Any(c => c != 'I')) h.Add((s, new long[] { 1, 2, 3, -1, -2, 5 }[rng.Next(6)]));
            }
            if (h.Count == 0) continue;
            var e = new EndCount(W, n, h, jumps);
            var ls = e.LocalSystem();
            if (ls is null) continue;
            var up = e.UpperCounts();
            Assert.True((ls.Value.Near, ls.Value.Far) == (up.Near, up.Far),
                $"N={n} jumps {string.Join(",", jumps)} H {string.Join(" + ", h.Select(t => $"{t.Item2}{t.Item1}"))}: " +
                $"local system {ls.Value.Near}, {ls.Value.Far}; end count {up.Near}, {up.Far}");
            rows++;
            if (up.Near == up.Far) pal++; else broken++;
            if (nFree == 2) twoFree++;
        }
        Assert.True(rows >= 200 && pal > 30 && broken > 30 && twoFree > 20,
            $"rows {rows}, palindromic {pal}, broken {broken}, two undephased {twoFree}");
    }

    // Theorem 6 with one undephased site: traces of words of length <= 3 in the generators decide the palindrome,
    // on every row where the end count's verdict is exact.
    [Fact]
    public void One_Undephased_Site_Words_Of_Length_Three_Decide()
    {
        var rng = new Random(20261002 + 11);
        int rows = 0, pal = 0, broken = 0;
        for (int draw = 0; draw < 900 && rows < 200; draw++)
        {
            int n = rng.Next(2, 5);
            int u = rng.Next(n);
            var jumps = Enumerable.Range(0, n).Where(l => l != u).Select(l => One(n, l, "XYZ"[rng.Next(3)])).ToList();
            var h = new List<(string, long)>();
            int nt = rng.Next(2, 7);
            for (int k = 0; k < nt; k++)
            {
                var s = new string(Enumerable.Range(0, n).Select(_ => "IIXYZ"[rng.Next(5)]).ToArray());
                if (s.Any(c => c != 'I')) h.Add((s, new long[] { 1, 2, 3, -1, -2, 5 }[rng.Next(6)]));
            }
            if (h.Count == 0) continue;
            var e = new EndCount(W, n, h, jumps);
            var agree = e.LocalSystemTracesAgree(3);
            if (agree is null) continue;
            var r = e.Verdict();
            if (!EndCount.IsExact(r)) continue;
            Assert.True(agree.Value == EndCount.IsPalindrome(r),
                $"N={n} jumps {string.Join(",", jumps)} H {string.Join(" + ", h.Select(t => $"{t.Item2}{t.Item1}"))}: traces {agree}, verdict {r}");
            rows++;
            if (agree.Value) pal++; else broken++;
        }
        Assert.True(rows >= 150 && pal > 30 && broken > 30, $"rows {rows}, palindromic {pal}, broken {broken}");
    }

    // Length 3 is needed: on this row every word of length 1 or 2 has equal traces, a word of length 3 does not,
    // and the row is broken (exactly). Jumps Z on site 0 and X on site 2, site 1 undephased.
    [Fact]
    public void Words_Of_Length_Two_Do_Not_Suffice()
    {
        var e = new EndCount(W, 3, new List<(string, long)>
            { ("XYY", 1), ("IYI", 2), ("IIY", 3), ("ZXX", 3), ("IZX", 1) }, new[] { "ZII", "IIX" });
        Assert.True(e.LocalSystemTracesAgree(2));
        Assert.False(e.LocalSystemTracesAgree(3));
        var r = e.Verdict();
        Assert.True(EndCount.IsExact(r) && !EndCount.IsPalindrome(r), $"{r}");
        Assert.Equal((1, 0), e.UpperCounts());
        Assert.Equal((1, 0, 1), e.LocalSystem());
    }

    // The edge positives H_ab·H_ab^† carry weight: one edge, no vertex blocks, a block 2 + X that is invertible and
    // not unitary, so H·H^† = 5 + 4X is not a scalar. Its commutant is 2-dimensional; the forward transports alone
    // would leave all 4 of the 2 x 2 matrices. The local system meets the end count on this row and on its
    // variants with a Y or a Z beside the X.
    [Theory]
    [InlineData("XX")]
    [InlineData("XY")]
    [InlineData("XZ")]
    public void The_Edge_Positive_Is_A_Generator(string second)
    {
        var e = new EndCount(W, 2, new List<(string, long)> { ("XI", 2), (second, 1) }, new[] { "ZI" });
        var ls = e.LocalSystem();
        Assert.NotNull(ls);
        var up = e.UpperCounts();
        Assert.Equal((up.Near, up.Far), (ls!.Value.Near, ls.Value.Far));
        Assert.Equal(2, up.Near);
    }

    // A singular chord does no harm: two dephased sites under Z, one undephased, the flip of both dephased bits
    // carrying the rank-one block X + i(-1)^y Y = 2 tau on the undephased site (XXX + XYY, no scalar part), the
    // flips of one bit invertible (2 + X and 3 + Z), so a four-cycle of invertible pairs spans the component. The tree avoids the singular pair and the counts still meet the end count.
    [Fact]
    public void A_Singular_Chord_Is_Read()
    {
        var e = new EndCount(W, 3, new List<(string, long)>
            { ("XII", 2), ("XIX", 1), ("IXI", 3), ("IXZ", 1), ("XXX", 1), ("XYY", 1) }, new[] { "ZII", "IZI" });
        var ls = e.LocalSystem();
        Assert.NotNull(ls);
        var up = e.UpperCounts();
        Assert.Equal((up.Near, up.Far), (ls!.Value.Near, ls.Value.Far));
        Assert.Equal((1, 0), (up.Near, up.Far));
    }

    // Both directions of every pair are generators. Keeping only the pairs (a, b) with a > b fails the edge-positive
    // rows; keeping only those with a < b reads this row, broken (near 1, far 0), as far 1, a palindrome that is not
    // there. Found by a search over 2,848 random rows, the only one of them on which that half fails.
    [Fact]
    public void Both_Directions_Of_Every_Pair_Are_Needed()
    {
        var e = new EndCount(W, 3, new List<(string, long)>
            { ("XIZ", 5), ("XIY", 3), ("YZZ", 5), ("XII", 2), ("IZY", 3) }, new[] { "YII", "IXI" });
        Assert.Equal((1, 0), e.UpperCounts());
        Assert.Equal((1, 0), (e.LocalSystem()!.Value.Near, e.LocalSystem()!.Value.Far));
    }

    // Theorem 6 (a): a pair in the support of H and not of H̄ rules out an invertible far element, but the far end need
    // not vanish, since such a pair is never a tree edge and a singular chord leaves room. H = (XII + XIZ + XZI + XZZ)
    // + 2 IXI + 2 XXI under Z on the two dephased sites (the proof's row scaled by 2): near 2,
    // far 1, broken; the local system declines it, the supports differing.
    [Fact]
    public void A_Mismatched_Support_Leaves_The_Far_End_Nonzero()
    {
        var e = new EndCount(W, 3, new List<(string, long)>
            { ("XII", 1), ("XIZ", 1), ("XZI", 1), ("XZZ", 1), ("IXI", 2), ("XXI", 2) }, new[] { "ZII", "IZI" });
        Assert.Equal((2, 1), e.UpperCounts());
        Assert.False(EndCount.IsPalindrome(e.Verdict()));
        Assert.Null(e.LocalSystem());
    }

    // Outside Theorem 6's reading: no undephased site (Theorem 2's case), or a component that only a singular block
    // connects (the Heisenberg flip at a site with no transverse field, 2J·tau, of rank one).
    [Fact]
    public void Local_System_Declines_Where_It_Does_Not_Read()
    {
        var all = new EndCount(W, 2, new List<(string, long)> { ("XX", 1), ("YY", 1) }, new[] { "ZI", "IZ" });
        Assert.Null(all.LocalSystem());
        var heis = new EndCount(W, 2, new List<(string, long)> { ("XX", 1), ("YY", 1), ("ZZ", 1) }, new[] { "ZI" });
        Assert.Null(heis.LocalSystem());
    }

    // A jump on two sites that H does not join: the grammar keeps the two sites in one component (a
    // jump split between two components would make a product that anticommutes with nothing), so what
    // it returns passes its check, or it returns nothing. H = 3·ZI + 5·IZ, jump XX.
    [Fact]
    public void A_Jump_Across_Two_Components_Keeps_Them_Together()
    {
        var e = new EndCount(W, 2, new List<(string, long)> { ("ZI", 3), ("IZ", 5) }, new[] { "XX" });
        var g = e.GrammarElement();
        if (g is { } found) Assert.True(e.CheckElement(found.Element) is { AllLit: true, CommutesWithH: true });
    }

    [Fact]
    public void Invertible_Reads_The_Determinant()
    {
        var e1 = Letters(1, Array.Empty<(int, int)>(), "Z", new long[] { 1 }, "X", ".");
        Assert.False(e1.Invertible(new[] { ("I", 1L), ("Z", 1L) }));
        Assert.False(e1.Invertible(new[] { ("I", 1L), ("Y", 1L) }));            // det = 1 − i·(−i) = 0; without the i it would be 2
        Assert.True(e1.Invertible(new[] { ("X", 1L), ("Z", 1L) }));
        var e2 = Letters(2, Array.Empty<(int, int)>(), "Z", new long[] { 1, 1 }, "X.", "..");
        Assert.True(e2.Invertible(new[] { ("IX", 1L), ("ZX", 1L), ("IY", 3L), ("ZY", -3L) }));
        Assert.False(e2.Invertible(new[] { ("IX", 1L), ("ZX", 1L) }));
    }

    // E · P · E for E a combination of strings with integer coefficients, as strings with Gaussian
    // integer coefficients (re, im); a route through PauliString.Multiply alone.
    static Dictionary<string, (long, long)> Conjugate(Dictionary<string, long> e, string p)
    {
        var acc = new Dictionary<PauliString, (long Re, long Im)>();
        var ps = PauliString.Parse(p);
        foreach (var (a, ca) in e)
            foreach (var (b, cb) in e)
            {
                var (ab, k1) = PauliString.Multiply(PauliString.Parse(a), ps);
                var (abc, k2) = PauliString.Multiply(ab, PauliString.Parse(b));
                long c = ca * cb;
                (long re, long im) v = ((k1 + k2) & 3) switch { 0 => (c, 0), 1 => (0, c), 2 => (-c, 0), _ => (0, -c) };
                var old = acc.TryGetValue(abc, out var o) ? o : (0, 0);
                acc[abc] = (old.Re + v.re, old.Im + v.im);
            }
        return acc.Where(kv => kv.Value != (0, 0)).ToDictionary(kv => kv.Key.ToString(p.Length), kv => ((long, long))kv.Value);
    }
    // ---- every site dephased: the complement connection ----

    // docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md, Theorem 2, against the end count's exact
    // verdict: with every site dephased, the palindrome holds exactly when every component of the
    // hopping graph (joined with its complement image) carries a flat section of the complement
    // connection; the hopping components count the near end and the good components the far end.
    // Every dephasing pattern with one letter on every site and every field pattern, three bond sets,
    // two magnitude tuples, the path and the triangle at N = 3, and Heisenberg on the path at N = 4.
    [Theory]
    [InlineData(3, "chain", "XYZ", new long[] { 30, 22, 41 })]
    [InlineData(3, "chain", "XYZ", new long[] { 30, 30, 30 })]
    [InlineData(3, "complete", "XYZ", new long[] { 30, 22, 41 })]
    [InlineData(3, "chain", "XY", new long[] { 30, 22, 41 })]
    [InlineData(3, "chain", "XY", new long[] { 30, 30, 30 })]
    [InlineData(3, "complete", "XY", new long[] { 30, 30, 30 })]
    [InlineData(3, "chain", "ZZ", new long[] { 30, 22, 41 })]
    [InlineData(3, "complete", "ZZ", new long[] { 30, 30, 30 })]
    [InlineData(4, "chain", "XYZ", new long[] { 30, 22, 41, 17 })]
    public void Every_Site_Dephased_The_Palindrome_Is_A_Flat_Complement_Connection(int n, string topology, string bondLetters, long[] mags)
    {
        var edges = topology == "chain" ? Chain(n)
            : (from a in Enumerable.Range(0, n) from b in Enumerable.Range(a + 1, n - a - 1) select (a, b)).ToArray();
        const string letters = "XYZ", fieldAlphabet = ".XYZ";
        int rows = 0, palindromes = 0;
        for (int dc = 0; dc < (int)Math.Pow(3, n); dc++)
            for (int fc = 0; fc < 1 << (2 * n); fc++)
            {
                var deph = new string(Enumerable.Range(0, n).Select(l => letters[dc / (int)Math.Pow(3, l) % 3]).ToArray());
                var field = new string(Enumerable.Range(0, n).Select(l => fieldAlphabet[(fc >> (2 * l)) & 3]).ToArray());
                var e = Letters(n, edges, bondLetters, mags, deph, field, weight: 100);
                var r = e.Verdict();
                Assert.True(EndCount.IsExact(r), $"{deph} {field}: {r}");
                var c = e.ComplementConnection()!.Value;
                Assert.Equal(EndCount.IsPalindrome(r), c.AllGood);
                // both counts certified by the end count's own route: the lifted, exactly checked
                // vectors (a lower bound) meet the modular ranks (an upper bound) at the graph's value
                Assert.Equal((c.HoppingComponents, c.Good), e.UpperCounts());
                Assert.Equal((c.HoppingComponents, c.Good), e.LiftedLowerCounts());
                rows++;
                if (c.AllGood) palindromes++;
            }
        Assert.True(palindromes > 0 && palindromes < rows);
    }

    // The turn to Z is a proper rotation for each jump letter: it sends the jump to Z, and it keeps the
    // product table, P·Q = i·R for (P, Q, R) cyclic in (X, Y, Z), judged by PauliString.Multiply.
    [Theory]
    [InlineData('X')]
    [InlineData('Y')]
    [InlineData('Z')]
    public void The_Turn_To_Z_Is_A_Proper_Rotation(char jump)
    {
        Assert.Equal(('Z', 1), EndCount.TurnToZ(jump, jump));
        foreach (var (p, q, r) in new[] { ('X', 'Y', 'Z'), ('Y', 'Z', 'X'), ('Z', 'X', 'Y') })
        {
            var (tp, sp) = EndCount.TurnToZ(p, jump);
            var (tq, sq) = EndCount.TurnToZ(q, jump);
            var (tr, sr) = EndCount.TurnToZ(r, jump);
            var (prod, k) = PauliString.Multiply(PauliString.Parse(tp.ToString()), PauliString.Parse(tq.ToString()));
            Assert.Equal(tr.ToString(), prod.ToString(1));
            // sp·sq·(Tp·Tq) = sp·sq·i^k·Tr must equal i·sr·Tr
            Assert.True((k == 1 && sp * sq == sr) || (k == 3 && sp * sq == -sr), $"{p}{q}={r} under jump {jump}: k={k}");
        }
    }

    // The same theorem for Pauli Hamiltonians the families never reach: random strings of any length
    // with random signed coefficients, every site dephased along a random letter, N = 2 to 4, a fixed
    // seed. Rows whose verdict the ranks alone gave are skipped (the connection is exact, the ranks are
    // not); the palindromic and the broken rows among the exact ones are both required.
    [Fact]
    public void The_Complement_Connection_Holds_For_Random_Pauli_Hamiltonians()
    {
        var rng = new Random(20260930);
        int exact = 0, pal = 0;
        for (int trial = 0; trial < 600; trial++)
        {
            int n = rng.Next(2, 5);
            var h = new List<(string, long)>();
            int nt = rng.Next(1, 7);
            for (int k = 0; k < nt; k++)
                h.Add((new string(Enumerable.Range(0, n).Select(_ => "IXYZ"[rng.Next(4)]).ToArray()), rng.Next(1, 4) * (rng.Next(2) == 0 ? 1 : -1)));
            if (h.All(t => t.Item1.All(ch => ch == 'I'))) continue;
            var jumps = Enumerable.Range(0, n).Select(l => One(n, l, "XYZ"[rng.Next(3)])).ToList();
            var e = new EndCount(W, n, h, jumps);
            var r = e.Verdict();
            if (!EndCount.IsExact(r)) continue;
            exact++;
            var c = e.ComplementConnection()!.Value;
            Assert.Equal(EndCount.IsPalindrome(r), c.AllGood);
            Assert.Equal((c.HoppingComponents, c.Good), e.UpperCounts());
            Assert.Equal((c.HoppingComponents, c.Good), e.LiftedLowerCounts());
            if (c.AllGood) pal++;
        }
        Assert.True(exact > 400 && pal > 20 && pal < exact, $"exact {exact}, palindromic {pal}");
    }

    // Theorem 4: undephased sites whose letter H conserves. Random Pauli Hamiltonians at N = 2 to 5, one or
    // two undephased sites, each with a random letter P_u that every term carries or leaves alone at u, the
    // other sites dephased along random letters; the sector connection's two counts against the end
    // count's, both of its routes (modular ranks from above, lifted exactly checked vectors from below),
    // on every row where the two routes meet. Rows with a nonzero cross-sector term are required, and
    // with no undephased site the sector connection is Theorem 2's reading.
    [Fact]
    public void The_Sector_Connection_Counts_Both_Ends_With_Conserved_Letters_At_Undephased_Sites()
    {
        var rng = new Random(20261001);
        int certified = 0, pal = 0, twoFree = 0, crossRows = 0;
        for (int trial = 0; trial < 700; trial++)
        {
            int n = rng.Next(2, 6);
            int nFree = n >= 3 && rng.Next(3) == 0 ? 2 : 1;
            var free = Enumerable.Range(0, n).OrderBy(_ => rng.Next()).Take(nFree).ToHashSet();
            var kept = Enumerable.Range(0, n).ToDictionary(l => l, _ => "XYZ"[rng.Next(3)]);
            var h = new List<(string, long)>();
            int nt = rng.Next(1, 8);
            for (int k = 0; k < nt; k++)
                h.Add((new string(Enumerable.Range(0, n).Select(l => free.Contains(l)
                        ? (rng.Next(2) == 0 ? 'I' : kept[l])
                        : "IXYZ"[rng.Next(4)]).ToArray()),
                    rng.Next(1, 4) * (rng.Next(2) == 0 ? 1 : -1)));
            if (h.All(t => t.Item1.All(ch => ch == 'I'))) continue;
            var jumps = Enumerable.Range(0, n).Where(l => !free.Contains(l)).Select(l => One(n, l, kept[l])).ToList();
            var e = new EndCount(W, n, h, jumps);
            var c = e.SectorConnection()!.Value;
            Assert.Equal(1 << nFree, c.Sectors);
            var up = e.UpperCounts();
            if (up != e.LiftedLowerCounts()) continue;
            certified++;
            Assert.True((c.Near, c.Far) == up, $"{string.Join(" + ", h)} jumps {string.Join(",", jumps)}: sectors {c}, counts {up}");
            if (c.Near == c.Far) pal++;
            if (nFree == 2) twoFree++;
            if (c.Cross > 0) crossRows++;
        }
        Assert.True(certified > 450 && pal > 30 && pal < certified && twoFree > 50 && crossRows > 30,
            $"certified {certified}, palindromic {pal}, two undephased {twoFree}, cross-sector {crossRows}");
    }

    // Theorem 4 on the census's own one-letter family: ZZ bonds (weight 100) with fields 30, 22, 41 on the
    // path at N = 3, every dephasing pattern leaving at least one site undephased and at least one
    // dephased, every field pattern.
    // The sector connection applies exactly where every undephased site's field is absent or Z (derived:
    // one undephased site, 27 dephasing patterns times 2·16 field patterns, plus two, 9 times 4·4, so
    // 864 + 144 = 1008 rows); on each its counts are the end count's upper counts and its verdict the
    // exact one. The palindromic rows among them are 340 (measured). At XX + YY or XX + YY + ZZ bonds
    // every site of the path carries two bond letters or more, so no undephased site keeps one whatever
    // the fields, and the theorem does not reach (pinned on one row each: the bonds alone decide it).
    [Fact]
    public void The_Sector_Connection_Reads_The_Census_Rows_With_An_Undephased_Site()
    {
        const string bonds = "ZZ";
        const string A = ".XYZ";
        int applied = 0, pal = 0, rows = 0;
        for (int dc = 1; dc < 64; dc++)
        {
            string deph = new(Enumerable.Range(0, 3).Select(l => A[(dc >> (2 * l)) & 3]).ToArray());
            if (!deph.Contains('.')) continue;
            for (int fc = 0; fc < 64; fc++)
            {
                string field = new(Enumerable.Range(0, 3).Select(l => A[(fc >> (2 * l)) & 3]).ToArray());
                var e = Letters(3, Chain(3), bonds, new long[] { 30, 22, 41 }, deph, field);
                rows++;
                if (e.SectorConnection() is not { } c) continue;
                applied++;
                Assert.True((c.Near, c.Far) == e.UpperCounts(), $"{bonds} {deph} {field}: {c} vs {e.UpperCounts()}");
                var r = e.Verdict();
                if (EndCount.IsExact(r)) Assert.Equal(EndCount.IsPalindrome(r), c.Near == c.Far);
                if (c.Near == c.Far) pal++;
            }
        }
        Assert.Equal((2304, 1008, 340), (rows, applied, pal));
        foreach (string other in new[] { "XY", "XYZ" })
            Assert.Null(Letters(3, Chain(3), other, new long[] { 30, 22, 41 }, "Z.Z", "...").SectorConnection());
    }

    [Fact]
    public void The_Sector_Connection_With_No_Undephased_Site_Is_Theorem_Two()
    {
        var rng = new Random(7);
        for (int trial = 0; trial < 200; trial++)
        {
            int n = rng.Next(2, 5);
            var h = Enumerable.Range(0, rng.Next(1, 6)).Select(_ =>
                (new string(Enumerable.Range(0, n).Select(_ => "IXYZ"[rng.Next(4)]).ToArray()), (long)rng.Next(1, 4))).ToList();
            var e = new EndCount(W, n, h, Enumerable.Range(0, n).Select(l => One(n, l, "XYZ"[rng.Next(3)])).ToList());
            var s = e.SectorConnection()!.Value;
            var c = e.ComplementConnection()!.Value;
            Assert.Equal((c.HoppingComponents, c.Good, 1, 0), (s.Near, s.Far, s.Sectors, s.Cross));
        }
    }

    // An undephased site carrying two letters conserves none: the sector connection declines (null), and
    // so does a jump on two sites.
    [Fact]
    public void The_Sector_Connection_Declines_Where_No_Letter_Is_Conserved()
    {
        Assert.Null(new EndCount(W, 2, new List<(string, long)> { ("XX", 1), ("ZI", 1) }, new[] { "IZ" }).SectorConnection());
        Assert.Null(new EndCount(W, 2, new List<(string, long)> { ("XX", 1) }, new[] { "ZZ" }).SectorConnection());
        Assert.NotNull(new EndCount(W, 2, new List<(string, long)> { ("XX", 1), ("XI", 1) }, new[] { "IZ" }).SectorConnection());
    }

    // Where Γ_H and its union with the complement image part: H = X_0 + X_0 Z_1 under Z jumps joins
    // |00⟩ to |10⟩ (the two terms add) but not |01⟩ to |11⟩ (they cancel), so the hopping graph has 3
    // components and the union, which adds the complement image |11⟩–|01⟩, has 2. The near end counts
    // the 3, and neither union component is good (the edge |01⟩–|11⟩ has H_xy = 0 against its
    // complement's 2), so the far end is 0: broken.
    [Fact]
    public void The_Near_End_Counts_The_Hopping_Graph_Not_Its_Union_With_The_Complement()
    {
        var e = new EndCount(W, 2, new List<(string, long)> { ("XI", 1), ("XZ", 1) }, new[] { "ZI", "IZ" });
        var c = e.ComplementConnection()!.Value;
        Assert.Equal((3, 2, 0, false), (c.HoppingComponents, c.UnionComponents, c.Good, c.AllGood));
        Assert.Equal((3, 0), e.UpperCounts());
    }

    // Connectedness is load-bearing in Theorem 1: on two disjoint bonds with an X field on one and a Y
    // field on the other, under Z on every site, the palindrome holds, carried by the colouring XXYY
    // (X on the first bond, Y on the second), one colour per component.
    [Fact]
    public void Two_Components_May_Carry_Different_Transverse_Letters()
    {
        var h = new List<(string, long)>();
        foreach (char p in "XYZ") { h.Add((Two(4, 0, 1, p), 100)); h.Add((Two(4, 2, 3, p), 100)); }
        h.Add(("XIII", 30)); h.Add(("IIYI", 22));
        var e = new EndCount(W, 4, h, new[] { "ZIII", "IZII", "IIZI", "IIIZ" });
        Assert.True(e.ComplementConnection()!.Value.AllGood);
        Assert.Equal(EndCount.Reading.PalindromeByColour, e.Verdict());
        Assert.Contains("XXYY", e.Colourings().Select(x => x.ToString(4)));
    }

    // Theorem 1: bonds J(XX + YY) + Δ·ZZ with every J nonzero and every Δ free (zero included, the XY
    // case) on a connected graph, every site dephased along Z, fields along X, Y or Z of ANY
    // magnitudes and signs: the palindrome holds exactly when no field lies along Z and the
    // others use one letter, and then (with a field) its carrier is the colouring X^N or Y^N. Random
    // nonzero magnitudes, N = 3 to 5 on the four graphs and N = 6 on the chain, every pattern of none,
    // X, Y or Z per site.
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Hopping_Bonds_Under_Z_Dephasing_Pair_Exactly_When_The_Fields_Share_A_Transverse_Letter(int n)
    {
        var rng = new Random(n);
        var graphs = new[]
        {
            Chain(n), Ring(n), Enumerable.Range(1, n - 1).Select(i => (0, i)).ToArray(),
            (from a in Enumerable.Range(0, n) from b in Enumerable.Range(a + 1, n - a - 1) select (a, b)).ToArray(),
        };
        foreach (var edges in n < 6 ? graphs : graphs.Take(1))          // N = 6 on the chain only, for time
            for (int fc = 0; fc < 1 << (2 * n); fc++)
            {
                var pattern = Enumerable.Range(0, n).Select(l => ".XYZ"[(fc >> (2 * l)) & 3]).ToArray();
                var h = new List<(string, long)>();
                foreach (var (i, j) in edges)
                {
                    long w = rng.Next(1, 200) * (rng.Next(2) == 0 ? 1 : -1);
                    h.Add((Two(n, i, j, 'X'), w)); h.Add((Two(n, i, j, 'Y'), w));
                    long delta = rng.Next(3) == 0 ? 0 : rng.Next(-200, 200);
                    if (delta != 0) h.Add((Two(n, i, j, 'Z'), delta));
                }
                for (int l = 0; l < n; l++)
                    if (pattern[l] != '.') h.Add((One(n, l, pattern[l]), rng.Next(1, 90) * (rng.Next(2) == 0 ? 1 : -1)));
                var e = new EndCount(W, n, h, Enumerable.Range(0, n).Select(l => One(n, l, 'Z')).ToList());
                bool oneLetter = !pattern.Contains('Z') && !(pattern.Contains('X') && pattern.Contains('Y'));
                Assert.Equal(oneLetter, e.ComplementConnection()!.Value.AllGood);
                var verdict = e.Verdict();
                Assert.True(EndCount.IsExact(verdict), $"{new string(pattern)}: {verdict}");
                Assert.Equal(oneLetter, EndCount.IsPalindrome(verdict));
                if (oneLetter && pattern.Any(ch => ch != '.'))
                {
                    char c = pattern.First(ch => ch != '.');
                    Assert.Equal(new[] { new string(c, n) }, e.Colourings().Select(x => x.ToString(n)).ToArray());
                }
            }
    }
}
