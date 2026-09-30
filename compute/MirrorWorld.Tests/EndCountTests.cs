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
        // sites), so only the ranks decide, and the verdict must say so.
        var e = Row(3, Chain(3), ".X.", "X.X");
        Assert.Equal(new[] { "III", "XXX" }, e.DarkStrings().Where(s => e.CommutesWithH(new List<(string, long)> { (s.ToString(3), 1) }))
            .Select(s => s.ToString(3)).OrderBy(s => s));
        Assert.Equal((2, 0), e.LowerCounts());
        Assert.Equal(EndCount.Reading.BrokenByRank, e.Verdict(maxPower: 1));
        Assert.Equal(EndCount.Reading.BrokenByWord, e.Verdict(maxPower: 2));

        Assert.True(EndCount.IsExact(EndCount.Reading.PalindromeByColour));
        Assert.True(EndCount.IsExact(EndCount.Reading.BrokenByWord));
        Assert.True(EndCount.IsExact(EndCount.Reading.BrokenByCount));
        Assert.False(EndCount.IsExact(EndCount.Reading.PalindromeByRank));
        Assert.False(EndCount.IsExact(EndCount.Reading.BrokenByRank));
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
        Assert.Equal(EndCount.Reading.PalindromeByRank, e.Verdict(maxPower: 2));

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
        Assert.Equal(EndCount.Reading.PalindromeByRank, e.Verdict(maxPower: 2, maxJumps: 1));

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
}
