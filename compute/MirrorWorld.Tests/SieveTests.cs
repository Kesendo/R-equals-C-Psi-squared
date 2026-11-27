using MirrorWorld;
using static MirrorWorld.LetterCube;

namespace MirrorWorldTests;

// From-below guard for the sieve (built 2026-11-10). The theodolite's every-N reading passes first and the end count
// (F158) is asked only on the rest; the guard asks the end count on what the sieve skipped as well, so a yes the
// reading gives is checked by a route that never sees the letters. The first census and its residue are pinned: on
// the chain nothing lies between the reading and the end count, on the ring the residue is the XZ ± ZY family and its
// twin, paired by a far element that is no product.
public class SieveTests
{
    static readonly Gauss[,] Zm = LetterMatrix('Z');

    // the first census, N up to 6, the theodolite alone (about ten seconds), shared by the tests that read it
    static readonly Lazy<List<Sieve.Entry>> Census = new(() =>
        Sieve.BondFamily(3, "Z").Select(c => Sieve.Run(c.Name, c.Cell, 6, stitch: false)).ToList());

    static readonly string[] RingResidue =
    {
        "+XZ+ZY", "+XZ-ZY", "+XZ+ZY+ZZ", "+XZ-ZY+ZZ", "+XZ+ZY-ZZ", "+XZ-ZY-ZZ",
        "+YZ+ZX", "+YZ-ZX", "+YZ+ZX+ZZ", "+YZ-ZX+ZZ", "+YZ+ZX-ZZ", "+YZ-ZX-ZZ",
    };

    [Fact]
    public void The_Sieve_Owns_The_Sieve_And_The_Residue_And_Inherits_The_Theodolite()
    {
        var sieve = new Sieve(new Theodolite(new LetterCube(new World())));
        Assert.Equal(new[] { "sieve", "residue" }, sieve.Own);
        Assert.Equal(new[] { "letters", "reading", "every N", "coordinates", "height", "step", "arrow", "moves", "x", "y", "z" },
            sieve.Inherited);
    }

    // the end count's verdicts on the residue are exact, read directly and not through the sieve's letters; and a jump
    // that is no Pauli letter is refused, not silently replaced
    [Fact]
    public void The_End_Count_Is_Exact_On_The_Residue_And_Only_Pauli_Jumps_Are_Read()
    {
        var family = Sieve.BondFamily(3, "Z").ToDictionary(c => c.Name, c => c.Cell);
        foreach (var name in RingResidue)
        {
            Assert.True(EndCount.IsExact(Sieve.Count(family[name].Normalized().Ring(4)).Verdict()), name);
            Assert.True(EndCount.IsExact(Sieve.Count(family[name].Normalized().Chain(4)).Verdict()), name);
        }
        var lossy = new Theodolite.Cell(new[] { (Lowering, 1L) }, null, new[] { ("XX", new[] { 0, 1 }, 1L) });
        Assert.Throws<ArgumentException>(() => Sieve.Run("loss", lossy, 4, stitch: false));
    }

    [Fact]
    public void The_Bond_Family_Holds_Every_Signed_Cell_Of_One_To_Three_Terms_Once()
    {
        var names = Sieve.BondFamily(3, "Z").Select(c => c.Name).ToList();
        Assert.Equal(9 + 36 * 2 + 84 * 4, names.Count);
        Assert.Equal(names.Count, names.Distinct().Count());
    }

    [Fact]
    public void The_Field_And_Window_Families_Hold_Every_Cell_Once()
    {
        var fields = Sieve.FieldFamily(2, "Z").Select(c => c.Name).ToList();
        Assert.Equal(81 * 27, fields.Count);
        Assert.Equal(fields.Count, fields.Distinct().Count());
        var windows = Sieve.WindowFamily(2, "Z").Select(c => c.Name).ToList();
        Assert.Equal(36 + 630 * 2, windows.Count);
        Assert.Equal(windows.Count, windows.Distinct().Count());
    }

    // The window family's chain residue, X·(X ± Y)·X and its relatives (the outer letters X or Y, the middle a letter of
    // X plus or minus Y), read by the free chain from four sites on (both templates of width three with a flipping letter
    // at each end, so the free chain reads from three, where the stitch's yes comes first): Read, and with it the window
    // family's chain holds no residue.
    [Fact]
    public void The_Window_Residue_On_The_Chain_Is_Read_By_The_Free_Chain()
    {
        var family = Sieve.WindowFamily(2, "Z").ToDictionary(c => c.Name, c => c.Cell);
        foreach (var name in new[] { "+XXX+XYX", "+XXX-XYX", "+XXY+XYY", "+XXY-XYY", "+YXX+YYX", "+YXX-YYX", "+YXY+YYY", "+YXY-YYY" })
        {
            var e = Sieve.Run(name, family[name], 6);
            Assert.Equal((Sieve.Verdict.Read, "sFFF"), (e.Chain.Verdict, e.Chain.Letters));
        }
    }

    // A one-shape cell of span five whose two templates are five wide with a flipping letter at each end: the free chain
    // reads its chain from five sites on, so the sieve leaves no length between the end count and the free chain (a guard
    // of max(2·span, 3·span − 3) = 12 would leave N = 11 unread), and the ring breaks exactly at the multiples of five
    // from ten sites (F173's remark "Which ring sizes break"), the end count breaking at ten.
    [Fact]
    public void A_Span_Five_Cell_Is_Read_On_The_Chain_From_Five_Sites()
    {
        var cell = Cell('Z', ("XXXXX", 1), ("XXYXX", 1));
        Assert.Equal(5, Stitch.FreeEveryN(cell.Normalized(), ring: false)!.From);
        for (int n = 5; n <= 10; n++)                                   // the free chain's new lengths, by the end count
        {
            var v = Sieve.Count(cell.Normalized().Chain(n)).Verdict();
            Assert.True(EndCount.IsPalindrome(v) && EndCount.IsExact(v), $"chain N = {n}: {v}");
        }
        var e = Sieve.Run("+XXXXX+XXYXX", cell, 10);
        Assert.Equal(Sieve.Verdict.Read, e.Chain.Verdict);
        Assert.Equal('B', e.Ring.Letters[^1]);                         // N = 10, the end count
        var law = Stitch.RingLaw(cell)!;
        Assert.All(Enumerable.Range(law.Guard, 63 - law.Guard + 1), n => Assert.Equal(n % 5 != 0, law.Pairs(n)));
    }

    static Theodolite.Cell Cell(char light, params (string Letters, long C)[] terms)
        => new(new[] { (LetterMatrix(light), 1L) }, null, terms.Select(t => (t.Letters, Enumerable.Range(0, t.Letters.Length).ToArray(), t.C)));

    // A ring law that exists decides every N of a pass whose letters below its guard are exact, whatever the verdict:
    // XY on the ring is Read (the theodolite on even rings, broken at the odd N read) and Theorem D decides every odd
    // ring as well; XIY needs its guard, N = 5, among the N read or one past them; XIX − XXZ + ZXX, Mixed at the N read,
    // is decided by Corollary A′ through its solid XIX from N = 6; a cell without a law is decided by none, and the chain
    // has no ring law.
    [Fact]
    public void A_Ring_Law_Decides_The_Lengths_A_Read_Cell_Leaves()
    {
        var xy = Sieve.Run("XY", Cell('Z', ("XY", 1)), 6);
        Assert.Equal((Sieve.Verdict.Read, "BEBE", true), (xy.Ring.Verdict, xy.Ring.Letters, xy.Ring.DecidedByLaw));
        Assert.Contains("Theorem D", xy.Ring.Law!.Text);
        Assert.False(xy.Chain.DecidedByLaw);
        var xiy = Sieve.Run("XIY", Cell('Z', ("XIY", 1)), 9);
        Assert.Equal((Sieve.Verdict.Read, "BEBBBEB", true), (xiy.Ring.Verdict, xiy.Ring.Letters, xiy.Ring.DecidedByLaw));
        Assert.False(Sieve.Run("XIY", Cell('Z', ("XIY", 1)), 3).Ring.DecidedByLaw);
        Assert.True(Sieve.Run("XIY", Cell('Z', ("XIY", 1)), 4).Ring.DecidedByLaw);   // the guard one past the N read
        // a cell whose lone excitation's hops cancel beside a solid shape: the law from N = 6, the N below it read exactly
        var solid = Sieve.Run("XIX-XXZ+ZXX", Cell('Z', ("XIX", 1), ("XXZ", -1), ("ZXX", 1)), 6);
        Assert.Equal((Sieve.Verdict.Law, "sEBE", true), (solid.Ring.Verdict, solid.Ring.Letters, solid.Ring.DecidedByLaw));
        Assert.Equal(6, solid.Ring.Law!.Guard);
        Assert.False(Sieve.Run("XIX-XXZ+ZXX", Cell('Z', ("XIX", 1), ("XXZ", -1), ("ZXX", 1)), 4).Ring.DecidedByLaw);
        // a cell no law reads (XXZ + YYZ + XYZ − YXZ: a hop with a twist beside a Z, strings of both signs): decided only at the N read
        var none = Sieve.Run("XXZ+YYZ+XYZ-YXZ", Cell('Z', ("XXZ", 1), ("YYZ", 1), ("XYZ", 1), ("YXZ", -1)), 6);
        Assert.Null(none.Ring.Law);
        Assert.False(none.Ring.DecidedByLaw);
    }

    // On the first census (with the stitch) the read ring cells that carry a broken letter at an uncovered N are 13, and a
    // ring law decides every N of each: one string per flip set by F173's Theorem D; XX ± XY − YY and XX ± YX − YY by
    // Theorem 9 of the complement connection, the pairing bond; the pure DM cells XY − YX, ± ZZ, by its Theorem 8.
    [Fact]
    public void On_The_First_Census_A_Ring_Law_Decides_All_Thirteen_Read_Cells_With_Broken_Letters()
    {
        var c = Sieve.BondFamily(3, "Z").Select(x => Sieve.Run(x.Name, x.Cell, 6)).ToList();
        var broken = c.Where(e => e.Ring.Verdict == Sieve.Verdict.Read && e.Ring.BrokenFrom is null
            && e.Ring.Letters.Any(ch => ch is 'B' or 'b' or 'C')).ToList();
        Assert.Equal(13, broken.Count);
        Assert.All(broken, e => Assert.True(e.Ring.DecidedByLaw, e.Name));
        Assert.All(broken.Where(e => e.Name.StartsWith("+XX")), e => Assert.Contains("Theorem 9", e.Ring.Law!.Text));
        Assert.Equal(new[] { "+XY-YX", "+XY-YX+ZZ", "+XY-YX-ZZ" },
            broken.Where(e => e.Ring.Law!.Text.Contains("Theorem 8")).Select(e => e.Name).OrderBy(n => n, StringComparer.Ordinal));
    }

    // Theorem 9 of the complement connection (the pairing bond), adopted in Stitch.RingLaw, held against the end count's
    // exact verdict on every bond J_XX·XX + J_YY·YY + J_XY·XY + J_YX·YX + J_ZZ·ZZ with coefficients in {−1, 0, 1, 2}: where
    // the pairing amplitude h is nonzero a ring law exists (F173's Theorem D or Corollary E2 where they read, Theorem 9
    // otherwise), and every one of them agrees with the end count on the ring of N = 3 to 7 sites; the theorem's chain
    // half (g real or imaginary) agrees with it on the open chain of 3 to 6 sites; where h = 0 the bond is Theorem 8's
    // and Theorem 9 is not asked. All three branches of Theorem 9 occur. Negative controls: a Z beside a bond whose
    // pairing vanishes, and a field beside the bond, are outside its class.
    [Fact]
    public void The_Pairing_Bond_Law_Agrees_With_The_End_Count_On_Chain_And_Ring()
    {
        var branches = new HashSet<string>();
        int byTheorem9 = 0, byOthers = 0;
        long[] values = { -1, 0, 1, 2 };
        foreach (long xx in values) foreach (long yy in values) foreach (long xy in values) foreach (long yx in values) foreach (long zz in new long[] { 0, 1 })
        {
            var terms = new[] { ("XX", xx), ("YY", yy), ("XY", xy), ("YX", yx), ("ZZ", zz) }.Where(t => t.Item2 != 0).ToArray();
            if (terms.All(t => t.Item1 == "ZZ")) continue;
            var cell = Cell('Z', terms);
            var law = Stitch.RingLaw(cell);
            string name = string.Join("+", terms);
            if (xx == yy && xy == -yx)                                                   // h = 0: Theorem 8's
            {
                Assert.True(law is null || !law.Text.Contains("Theorem 9"), $"{name}: h = 0 read by Theorem 9");
                continue;
            }
            Assert.NotNull(law);
            bool real = xy == yx, imaginary = xx == -yy;
            if (law!.Text.Contains("Theorem 9")) { byTheorem9++; branches.Add(real ? "real" : imaginary ? "imaginary" : "neither"); }
            else byOthers++;
            for (int n = 3; n <= 7; n++)
            {
                var v = Sieve.Count(cell.Ring(n)).Verdict();
                Assert.True(EndCount.IsExact(v), $"{name} ring {n}: not exact");
                Assert.True(law.Pairs(n) == EndCount.IsPalindrome(v), $"{name} ring {n}: {law.Text}");
                Assert.True((real || imaginary && n % 2 == 0) == EndCount.IsPalindrome(v), $"{name} ring {n}: the theorem");
            }
            for (int n = 3; n <= 6; n++)
            {
                var v = Sieve.Count(cell.Chain(n)).Verdict();
                Assert.True(EndCount.IsExact(v), $"{name} chain {n}: not exact");
                Assert.True((real || imaginary) == EndCount.IsPalindrome(v), $"{name} chain {n}");
            }
        }
        Assert.Equal(new[] { "imaginary", "neither", "real" }, branches.OrderBy(b => b, StringComparer.Ordinal));
        Assert.True(byTheorem9 > 100 && byOthers > 10, $"{byTheorem9} cells read by Theorem 9, {byOthers} by F173's laws");
        // outside the class: a Z beside a bond without a pairing (h = 0), a field beside it
        Theodolite.Cell Window(params (string, long)[] t) => new(new[] { (Zm, 1L) }, null, t.Select(x => (x.Item1, new[] { 0, 1, 2 }, x.Item2)));
        foreach (var outside in new[]
        {
            Window(("XXZ", 1), ("YYZ", 1)), Window(("ZXX", 1), ("ZYY", 1)),
            new Theodolite.Cell(new[] { (Zm, 1L) }, new[] { ('X', 1L) }, new[] { ("XX", new[] { 0, 1 }, 1L), ("YY", new[] { 0, 1 }, -1L) }),
        })
            Assert.True(Stitch.RingLaw(outside) is not { } l || !l.Text.Contains("Theorem 9"));
    }

    // Theorem 9 on bonds d sites apart: J_XX·XIX + J_YY·YIY + J_XY·XIY + J_YX·YIX (and a ZIZ part 0 or 1) with
    // coefficients in {−1, 0, 1}, held against the end count's exact verdict. Where h ≠ 0 a ring law exists, and it agrees
    // with the end count on the rings of 5 to 8 sites, where m = N/gcd(N, 2) takes the values 5, 3, 7 and 4 (so that the
    // six-site ring, two triangles, separates m even from N even); the theorem's chain half agrees on chains of 3 to 7
    // sites (pairing whatever g at N ≤ 2d, where the flips are independent; g real or imaginary from N = 2d + 1). At d = 3
    // the bonds with g = −i, g = 1 and g neither on rings of 7 to 10 sites: m = 10 at N = 10 separates m even from 4 dividing m,
    // and the three triangles at N = 9 break the imaginary bond and keep the real one (g = 1, a hop beside the pairing). Negative controls: two distances
    // in one cell, a Z between the two sites of a bond without a pairing, a field beside the bond.
    [Fact]
    public void The_Pairing_Bond_Law_At_Distance_D_Agrees_With_The_End_Count()
    {
        var branches = new HashSet<string>();
        int byTheorem9 = 0;
        long[] values = { -1, 0, 1 };
        foreach (long xx in values) foreach (long yy in values) foreach (long xy in values) foreach (long yx in values) foreach (long zz in new long[] { 0, 1 })
        {
            if (xx == yy && xy == -yx) continue;                                         // h = 0: Theorem 8's
            var terms = new[] { ("XIX", xx), ("YIY", yy), ("XIY", xy), ("YIX", yx), ("ZIZ", zz) }.Where(t => t.Item2 != 0).ToArray();
            var cell = Cell('Z', terms);
            var law = Stitch.RingLaw(cell);
            string name = string.Join("+", terms);
            Assert.NotNull(law);
            bool real = xy == yx, imaginary = xx == -yy;
            if (law!.Text.Contains("Theorem 9"))
            {
                byTheorem9++;
                branches.Add(real ? "real" : imaginary ? "imaginary" : "neither");
                Assert.Equal(5, law.Guard);
            }
            for (int n = 5; n <= 8; n++)
            {
                var v = Sieve.Count(cell.Ring(n)).Verdict();
                Assert.True(EndCount.IsExact(v), $"{name} ring {n}: not exact");
                int m = n / (n % 2 == 0 ? 2 : 1);
                Assert.True(law.Pairs(n) == EndCount.IsPalindrome(v), $"{name} ring {n}: {law.Text}");
                Assert.True((real || imaginary && m % 2 == 0) == EndCount.IsPalindrome(v), $"{name} ring {n}: the theorem");
            }
            for (int n = 3; n <= 7; n++)
            {
                var v = Sieve.Count(cell.Chain(n)).Verdict();
                Assert.True(EndCount.IsExact(v), $"{name} chain {n}: not exact");
                Assert.True((n <= 4 || real || imaginary) == EndCount.IsPalindrome(v), $"{name} chain {n}");
            }
        }
        Assert.Equal(new[] { "imaginary", "neither", "real" }, branches.OrderBy(b => b, StringComparer.Ordinal));
        Assert.True(byTheorem9 > 30, $"{byTheorem9} cells read by Theorem 9");
        foreach (var (terms, real, imaginary) in new[]
        {
            (new[] { ("XIIX", 1L), ("XIIY", 1L), ("YIIY", -1L) }, false, true),
            (new[] { ("XIIX", 1L), ("XIIY", 1L), ("YIIX", 1L) }, true, false),
            (new[] { ("XIIX", 1L), ("XIIY", 1L), ("YIIY", 1L) }, false, false),
        })
        {
            var cell = Cell('Z', terms);
            var law = Stitch.RingLaw(cell);
            Assert.True(law is { } l && l.Text.Contains("Theorem 9") && l.Guard == 7, $"{string.Join("+", terms)}: {law?.Text}");
            for (int n = 7; n <= 10; n++)
            {
                var v = Sieve.Count(cell.Ring(n)).Verdict();
                Assert.True(EndCount.IsExact(v), $"d = 3, ring {n}: not exact");
                int m = n / (int)System.Numerics.BigInteger.GreatestCommonDivisor(n, 3);
                Assert.True(law!.Pairs(n) == EndCount.IsPalindrome(v), $"{string.Join("+", terms)} ring {n}: {law.Text}");
                Assert.True((real || imaginary && m % 2 == 0) == EndCount.IsPalindrome(v), $"{string.Join("+", terms)} ring {n}: the theorem");
            }
        }
        foreach (var outside in new[]
        {
            Cell('Z', ("XX", 1), ("XIX", 1), ("YIY", -1)),
            Cell('Z', ("XZX", 1), ("YZY", 1)),
            new Theodolite.Cell(new[] { (Zm, 1L) }, new[] { ('X', 1L) }, new[] { ("XIX", new[] { 0, 1, 2 }, 1L), ("YIY", new[] { 0, 1, 2 }, -1L), ("XIY", new[] { 0, 1, 2 }, 1L) }),
        })
            Assert.True(Stitch.RingLaw(outside) is not { } l || !l.Text.Contains("Theorem 9"), Stitch.RingLaw(outside)?.Text);
    }

    // Theorem 9 with Z letters on other sites, the same letters at the same offsets on every term: a Z beside the bond on
    // the right (XXZ, …) or on the left (ZXX, …), one on each side (ZXXZ, …, an even count), and one between two sites two
    // apart (XZX, …). Every bond with J_XX, J_YY, J_XY, J_YX in {−1, 0, 1}, g ≠ 0 and h ≠ 0, has a ring law (Theorem 9,
    // or F173's laws where they read first), held against the end count's exact verdict on the rings from its guard to
    // N = 8, and the rule read here beside it: ε = (−1)^(number of Z letters), ρ = 1 for g real and −1 for g imaginary, the
    // ring pairing exactly when ερ = 1, or ερ = −1 and N/gcd(N, d) is even. Every branch occurs. Between two sites two apart
    // with J_YY = 0 the bond is F116's golden window t₁·XZX + t₂·XZY + t₃·YZX, whose ring soft set PROOF_CEILING_GOLDEN_ROUTER
    // reads by another route (where 4 | N it pairs iff t₁ = 0 or t₂ = t₃, otherwise iff t₁ = 0, every real triple, N ≥ 5),
    // asserted on every such triple, g = 0 and h = 0 included. With g = 0 (XXZ − YYZ among them) the ring pairs at every
    // N, carried by ((X + Y)/√2)^⊗N or X^⊗N after one common turn, and the law says so.
    [Fact]
    public void The_Pairing_Bond_Law_With_Z_Letters_Agrees_With_The_End_Count()
    {
        var branches = new HashSet<string>();
        int byTheorem9 = 0, golden = 0, pureZero = 0;
        long[] values = { -1, 0, 1 };
        var placements = new (string Name, Func<char, char, string> Make, int Zs, int D)[]
        {
            ("right", (a, b) => $"{a}{b}Z", 1, 1), ("left", (a, b) => $"Z{a}{b}", 1, 1),
            ("both", (a, b) => $"Z{a}{b}Z", 2, 1), ("between", (a, b) => $"{a}Z{b}", 1, 2),
        };
        foreach (var (where, make, zs, d) in placements)
            foreach (long xx in values) foreach (long yy in values) foreach (long xy in values) foreach (long yx in values)
            {
                var terms = new[] { (make('X', 'X'), xx), (make('Y', 'Y'), yy), (make('X', 'Y'), xy), (make('Y', 'X'), yx) }
                    .Where(t => t.Item2 != 0).ToArray();
                if (terms.Length == 0) continue;
                var cell = Cell('Z', terms);
                string name = $"{where} " + string.Join("+", terms);
                var law = Stitch.RingLaw(cell);
                bool hZero = xx == yy && xy == -yx, gZero = xx == -yy && xy == yx;
                bool real = xy == yx, imaginary = xx == -yy;
                if (where == "between" && yy == 0)
                {
                    golden++;
                    for (int n = 5; n <= 8; n++)
                    {
                        var v = Sieve.Count(cell.Ring(n)).Verdict();
                        Assert.True(EndCount.IsExact(v), $"{name} ring {n}: not exact");
                        bool soft = n % 4 == 0 ? xx == 0 || xy == yx : xx == 0;
                        Assert.True(soft == EndCount.IsPalindrome(v), $"{name} ring {n}: the golden router's ring soft set");
                    }
                }
                if (hZero) continue;
                Assert.NotNull(law);
                if (gZero)
                {
                    Assert.True(law!.Text.Contains("Theorem 9") && law.Text.Contains("pairs on every ring"), $"{name}: {law.Text}");
                    pureZero++;
                    for (int n = law.Guard; n <= 8; n++)
                    {
                        var v = Sieve.Count(cell.Ring(n)).Verdict();
                        Assert.True(EndCount.IsExact(v) && EndCount.IsPalindrome(v), $"{name} ring {n}: g = 0 pairs");
                    }
                    continue;
                }
                int eps = zs % 2 == 0 ? 1 : -1, rho = real ? 1 : imaginary ? -1 : 0;
                if (law!.Text.Contains("Theorem 9")) { byTheorem9++; branches.Add($"ε={eps} {(rho == 1 ? "real" : rho == -1 ? "imaginary" : "neither")}"); }
                for (int n = law.Guard; n <= 8; n++)
                {
                    var v = Sieve.Count(cell.Ring(n)).Verdict();
                    Assert.True(EndCount.IsExact(v), $"{name} ring {n}: not exact");
                    int m = n / (int)System.Numerics.BigInteger.GreatestCommonDivisor(n, d);
                    Assert.True(law.Pairs(n) == EndCount.IsPalindrome(v), $"{name} ring {n}: {law.Text}");
                    bool rule = rho != 0 && (eps * rho == 1 || m % 2 == 0);
                    Assert.True(rule == EndCount.IsPalindrome(v), $"{name} ring {n}: the rule (ερ)^m");
                }
            }
        Assert.Equal(new[] { "ε=-1 imaginary", "ε=-1 neither", "ε=-1 real", "ε=1 imaginary", "ε=1 neither", "ε=1 real" },
            branches.OrderBy(b => b, StringComparer.Ordinal));
        Assert.True(byTheorem9 > 60 && golden == 26 && pureZero > 10, $"{byTheorem9} cells read by Theorem 9, {golden} golden triples, {pureZero} with g = 0");
    }

    // Theorem 10 of the complement connection (a pairing bond under the control 1 + κ·Z_(l−1)Z_(l+2) of its outer
    // neighbours), adopted in Stitch.RingLaw: every bond J_XX·XX + J_YY·YY + J_XY·XY + J_YX·YX with coefficients in {−1,
    // 0, 1} and h ≠ 0, under κ = +1 and κ = −1, gets the law from N = 7, and the law and the rule (ρ⁴ = 1 and ρ^(wN) = 1,
    // w = 1 under κ = +1, 2 under κ = −1) agree with Theorem 2 read whole (EndCount.ComplementConnection, exact both ways)
    // on the rings of 7 to 10 sites and with the end count's exact verdict on the rings of 7 and 8; on the open chain,
    // where the bond and its control sit only where four sites fit, the rule ρ⁴ = 1 agrees with Theorem 2 on the chains of
    // 6 to 9 sites. Every order of ρ occurs under both controls; a hop alone under the control (h = 0) is not read by
    // Theorem 10.
    [Fact]
    public void The_Controlled_Pairing_Bond_Law_Agrees_With_Theorem_2()
    {
        var branches = new HashSet<string>();
        int read = 0;
        long[] values = { -1, 0, 1 };
        foreach (int kappa in new[] { 1, -1 })
            foreach (long xx in values) foreach (long yy in values) foreach (long xy in values) foreach (long yx in values)
            {
                if (xx == yy && xy == -yx) continue;                                         // h = 0: outside this law
                var bond = new[] { ("XX", xx), ("YY", yy), ("XY", xy), ("YX", yx) }.Where(t => t.Item2 != 0).ToArray();
                var cell = Cell('Z', bond.Concat(bond.Select(t => ("Z" + t.Item1 + "Z", kappa * t.Item2))).ToArray());
                string name = $"κ = {kappa}, " + string.Join("+", bond);
                var law = Stitch.RingLaw(cell);
                Assert.True(law is { Guard: 7 } && law.Text.Contains("Theorem 10"), $"{name}: {law?.Text}");
                read++;
                long re = xx + yy, im = yx - xy;
                int order = im == 0 ? 1 : re == 0 ? 2 : Math.Abs(re) == Math.Abs(im) ? 4 : 0, w = kappa > 0 ? 1 : 2;
                branches.Add($"κ = {kappa}, order {order}");
                for (int n = 7; n <= 10; n++)
                {
                    var count = Sieve.Count(cell.Ring(n));
                    var cc = count.ComplementConnection();
                    Assert.True(cc is not null, $"{name} ring {n}: no reading");
                    Assert.True(law!.Pairs(n) == cc!.Value.AllGood, $"{name} ring {n}: {law.Text}");
                    Assert.True((order != 0 && (w * n) % order == 0) == cc.Value.AllGood, $"{name} ring {n}: the rule");
                    if (n > 8) continue;
                    var v = count.Verdict();
                    Assert.True(EndCount.IsExact(v) && EndCount.IsPalindrome(v) == cc.Value.AllGood, $"{name} ring {n}: the end count");
                }
                for (int n = 6; n <= 9; n++)
                {
                    var cc = ControlledChain(n, bond, kappa).ComplementConnection();
                    Assert.True(cc is not null && (order != 0) == cc.Value.AllGood, $"{name} chain {n}: ρ⁴ = 1");
                }
            }
        Assert.Equal(8, branches.Count);
        Assert.Equal(144, read);
        // h = 0 (a hop alone under the control, XX + YY: the folded XXZ chain) is outside this law
        Assert.True(Stitch.RingLaw(Cell('Z', ("XX", 1), ("YY", 1), ("ZXXZ", 1), ("ZYYZ", 1))) is not { } folded || !folded.Text.Contains("Theorem 10"));
    }

    // The controlled bond on the open chain of n sites under Z on every site: the bond on the sites l, l + 1 and the same
    // bond with Z on l − 1 and l + 2, times κ, for l = 1, …, n − 3, and each diagonal string at every window where it fits.
    static EndCount ControlledChain(int n, (string L, long C)[] bond, int kappa, (string L, long C)[]? diagonal = null)
        => ControlledChain(n, bond, bond.Select(b => (b.L, kappa * b.C)).ToArray(), diagonal);

    // the same with the strings under the control carrying their own bond: bond on the sites l, l + 1 and
    // controlledBond there with Z on l − 1 and l + 2, for l = 1, …, n − 3
    static EndCount ControlledChain(int n, (string L, long C)[] bond, (string L, long C)[] controlledBond, (string L, long C)[]? diagonal = null)
    {
        var h = new List<(string, long)>();
        for (int l = 1; l <= n - 3; l++)
        {
            foreach (var (letters, c) in bond)
            {
                var plain = Enumerable.Repeat('I', n).ToArray();
                plain[l] = letters[0]; plain[l + 1] = letters[1];
                h.Add((new string(plain), c));
            }
            foreach (var (letters, c) in controlledBond)
            {
                var controlled = Enumerable.Repeat('I', n).ToArray();
                controlled[l] = letters[0]; controlled[l + 1] = letters[1];
                controlled[l - 1] = 'Z'; controlled[l + 2] = 'Z';
                h.Add((new string(controlled), c));
            }
        }
        foreach (var (letters, c) in diagonal ?? Array.Empty<(string, long)>())
            for (int s = 0; s + letters.Length <= n; s++)
            {
                var p = Enumerable.Repeat('I', n).ToArray();
                for (int i = 0; i < letters.Length; i++) p[s + i] = letters[i];
                h.Add((new string(p), c));
            }
        var jumps = Enumerable.Range(0, n).Select(j => new string('I', j) + "Z" + new string('I', n - j - 1)).ToList();
        return new EndCount(new World(), n, h, jumps);
    }

    // Theorem 10 beside a diagonal part. The theorem allows a real diagonal part B with B_x = B_x̄, a sum of strings with
    // an even number of Z letters, which Stitch.RingLaw leaves out of the bond it reads. Such a part commutes with X^⊗N
    // and with every diagonal operator, so it moves neither end of the count; the verdicts with the part are checked
    // against the same rule as without it. Each of the 144 cells of the test above (72 bonds under each of the two
    // controls) gets one of four such parts, ZZ, 2·ZIZ, −ZIIZ (on the control's two sites) and ZZZZ, all within four sites
    // so that the guard stays 7, taken in turn within each order of ρ, so that every part meets every order under both
    // controls: the law (Theorem 10 from N = 7) and the rule agree with Theorem 2 read whole on the rings of 7 to 10 sites
    // and with the end count's exact verdict on the rings of 7 and 8; on the chains of 6 to 9 sites, with the part at
    // every window where it fits, the rule ρ⁴ = 1 agrees with Theorem 2. A part with an odd number of Z letters (Z, ZZZ),
    // where B_x̄ = −B_x, is refused by Stitch.RingLaw before any law is asked, and beside XX under κ = +1 and XY under
    // κ = −1, which pair without it, Theorem 2 finds a break on the rings of 7 and 8 and the chains of 6 to 9, and so does
    // the end count on the rings of 7 and 8.
    [Fact]
    public void The_Controlled_Pairing_Bond_Law_Holds_Beside_A_Diagonal_Part()
    {
        var parts = new[] { ("ZZ", 1L), ("ZIZ", 2L), ("ZIIZ", -1L), ("ZZZZ", 1L) };
        var turn = new Dictionary<string, int>();
        var met = new HashSet<string>();
        int read = 0;
        long[] values = { -1, 0, 1 };
        foreach (int kappa in new[] { 1, -1 })
            foreach (long xx in values) foreach (long yy in values) foreach (long xy in values) foreach (long yx in values)
            {
                if (xx == yy && xy == -yx) continue;                                         // h = 0: outside this law
                var bond = new[] { ("XX", xx), ("YY", yy), ("XY", xy), ("YX", yx) }.Where(b => b.Item2 != 0).ToArray();
                long re = xx + yy, im = yx - xy;
                int order = im == 0 ? 1 : re == 0 ? 2 : Math.Abs(re) == Math.Abs(im) ? 4 : 0, w = kappa > 0 ? 1 : 2;
                string branch = $"κ = {kappa}, order {order}";
                int k = turn.TryGetValue(branch, out var seen) ? seen : 0;
                turn[branch] = k + 1;
                var part = parts[k % parts.Length];
                met.Add($"{branch}, {part.Item1}");
                var cell = Cell('Z', bond.Concat(bond.Select(b => ("Z" + b.Item1 + "Z", kappa * b.Item2))).Append(part).ToArray());
                string name = $"κ = {kappa}, " + string.Join("+", bond) + $" + {part}";
                var law = Stitch.RingLaw(cell);
                Assert.True(law is { Guard: 7 } && law.Text.Contains("Theorem 10"), $"{name}: {law?.Text}");
                read++;
                for (int n = 7; n <= 10; n++)
                {
                    var count = Sieve.Count(cell.Ring(n));
                    var cc = count.ComplementConnection();
                    Assert.True(cc is not null, $"{name} ring {n}: no reading");
                    Assert.True(law!.Pairs(n) == cc!.Value.AllGood, $"{name} ring {n}: {law.Text}");
                    Assert.True((order != 0 && (w * n) % order == 0) == cc.Value.AllGood, $"{name} ring {n}: the rule");
                    if (n > 8) continue;
                    var v = count.Verdict();
                    Assert.True(EndCount.IsExact(v) && EndCount.IsPalindrome(v) == cc.Value.AllGood, $"{name} ring {n}: the end count");
                }
                for (int n = 6; n <= 9; n++)
                {
                    var cc = ControlledChain(n, bond, kappa, new[] { part }).ComplementConnection();
                    Assert.True(cc is not null && (order != 0) == cc.Value.AllGood, $"{name} chain {n}: ρ⁴ = 1");
                }
            }
        Assert.Equal(144, read);
        Assert.Equal(32, met.Count);
        foreach (var odd in new[] { ("Z", 1L), ("ZZZ", 1L) })
            foreach (var (bond, kappa) in new[] { (new[] { ("XX", 1L) }, 1), (new[] { ("XY", 1L) }, -1) })
            {
                var cell = Cell('Z', bond.Concat(bond.Select(b => ("Z" + b.Item1 + "Z", kappa * b.Item2))).Append(odd).ToArray());
                string name = $"κ = {kappa}, {bond[0]} + {odd}";
                Assert.Null(Stitch.RingLaw(cell));
                for (int n = 7; n <= 8; n++)
                {
                    var count = Sieve.Count(cell.Ring(n));
                    var v = count.Verdict();
                    Assert.True(count.ComplementConnection() is { AllGood: false }, $"{name} ring {n}: Theorem 2");
                    Assert.True(EndCount.IsExact(v) && !EndCount.IsPalindrome(v), $"{name} ring {n}: the end count");
                }
                for (int n = 6; n <= 9; n++)
                    Assert.True(ControlledChain(n, bond, kappa, new[] { odd }).ComplementConnection() is { AllGood: false }, $"{name} chain {n}: Theorem 2");
            }
    }

    // Theorem 10's hop and pairing under opposite controls, adopted in Stitch.RingLaw (OppositeControlsBond, asked after
    // Corollary A′): the hop t·(XX + YY) + D·(XY − YX) under 1 + κ·Z_(l−1)Z_(l+2) and the pairing a·(XX − YY) + b·(XY + YX)
    // under 1 − κ·Z_(l−1)Z_(l+2), with t, D, a, b in {−1, 0, 1, 2}, hop and pairing nonzero, under κ = +1 and κ = −1
    // (450 cells; κ = +1, D = 0, b = 0 is the folded XYZ chain of Zadnik and Fagotti's Eq. (22)), so g = 2(t − iD). Every
    // cell gets a ring law from N = 7, and the law and the rule (ρ² = 1, and under κ = +1 also ρ^N = 1) agree with Theorem 2
    // read whole on the rings of 7 to 10 sites; the rule agrees with it on the rings of 5 and 6, below the sieve's guard,
    // and with the end count's exact verdict on the rings of 7 and 8 for the coefficients in {−1, 0, 1}; on the chains of
    // 5 to 9 sites, the bond and its controls only where four sites fit, the rule ρ² = 1 agrees with Theorem 2, and the
    // chain of four sites pairs at every g. Real, imaginary and neither occur under both controls, and both Corollary A′
    // and this law read cells; the hop imaginary under κ = +1 breaks at the odd N read, N = 9 among them.
    [Fact]
    public void The_Opposite_Controls_Law_Agrees_With_Theorem_2()
    {
        var branches = new HashSet<string>();
        int read = 0, byThisLaw = 0, byCorollary = 0;
        long[] values = { -1, 0, 1, 2 };
        foreach (int kappa in new[] { 1, -1 })
            foreach (long t in values) foreach (long d in values) foreach (long a in values) foreach (long b in values)
            {
                if ((t, d) == (0, 0) || (a, b) == (0, 0)) continue;
                var plain = new[] { ("XX", t + a), ("YY", t - a), ("XY", d + b), ("YX", -d + b) }.Where(x => x.Item2 != 0).ToArray();
                var controlled = new[] { ("XX", kappa * (t - a)), ("YY", kappa * (t + a)), ("XY", kappa * (d - b)), ("YX", kappa * (-d - b)) }
                    .Where(x => x.Item2 != 0).ToArray();
                var cell = Cell('Z', plain.Concat(controlled.Select(x => ("Z" + x.Item1 + "Z", x.Item2))).ToArray());
                string name = $"κ = {kappa}, t = {t}, D = {d}, a = {a}, b = {b}";
                bool real = d == 0, imaginary = t == 0;
                bool Rule(int n) => real || (imaginary && (kappa < 0 || n % 2 == 0));
                branches.Add($"κ = {kappa}, {(real ? "real" : imaginary ? "imaginary" : "neither")}");
                var law = Stitch.RingLaw(cell);
                Assert.True(law is { Guard: 7 }, $"{name}: {law?.Text}");
                if (law!.Text.Contains("the hop under")) byThisLaw++;
                else if (law.Text.Contains("Corollary A′")) byCorollary++;
                else Assert.Fail($"{name}: {law.Text}");
                read++;
                bool small = Math.Max(Math.Max(Math.Abs(t), Math.Abs(d)), Math.Max(Math.Abs(a), Math.Abs(b))) <= 1;
                for (int n = 5; n <= 10; n++)
                {
                    var count = Sieve.Count(cell.Ring(n));
                    var cc = count.ComplementConnection();
                    Assert.True(cc is not null, $"{name} ring {n}: no reading");
                    Assert.True(Rule(n) == cc!.Value.AllGood, $"{name} ring {n}: the rule");
                    if (n < 7) continue;
                    Assert.True(law.Pairs(n) == cc.Value.AllGood, $"{name} ring {n}: {law.Text}");
                    if (n > 8 || !small) continue;
                    var v = count.Verdict();
                    Assert.True(EndCount.IsExact(v) && EndCount.IsPalindrome(v) == cc.Value.AllGood, $"{name} ring {n}: the end count");
                }
                for (int n = 4; n <= 9; n++)
                {
                    var cc = ControlledChain(n, plain, controlled).ComplementConnection();
                    Assert.True(cc is not null && (n == 4 || real || imaginary) == cc.Value.AllGood, $"{name} chain {n}: ρ² = 1");
                }
            }
        Assert.Equal(6, branches.Count);
        Assert.Equal(450, read);
        Assert.True(byThisLaw > 0 && byCorollary > 0, $"{byThisLaw} cells read by this law, {byCorollary} by Corollary A′");
        var oddRings = Stitch.RingLaw(Cell('Z', ("XX", 1), ("YY", -1), ("XY", 1), ("YX", -1), ("ZXXZ", -1), ("ZYYZ", 1), ("ZXYZ", 1), ("ZYXZ", -1)));
        Assert.True(oddRings is not null && oddRings.Text.Contains("the hop under 1 + Z·Z") && oddRings.Pairs(8) && !oddRings.Pairs(9), oddRings?.Text);
    }

    // Theorem 8 of the complement connection read on a uniform cell (Stitch.RingLaw's HoppingRing), held against the end
    // count's exact verdict on the ring: hopping bonds J·(XX + YY) + D·(XY − YX) with J, D in {−1, 0, 1, 2}, a ZZ part 0
    // or 1, and the uniform fields 0, X, Y, X + Y, X − Y, on rings of 3 to 7 sites; and the same bonds two sites apart
    // (XIX, YIY, XIY, YIX) without fields on rings of 5 to 8 sites, where the cycles have length N/gcd(N, 2). Every order
    // of t/t̄ occurs (1, 2, 4 and none), and fields beside a twisted and an untwisted hop.
    [Fact]
    public void The_Hopping_Ring_Law_Agrees_With_The_End_Count()
    {
        var kinds = new HashSet<string>();
        long[] values = { -1, 0, 1, 2 };
        foreach (int d in new[] { 1, 2 })
            foreach (long j in values) foreach (long dm in values) foreach (long zz in new long[] { 0, 1 })
                foreach (var (a, b) in d == 1 ? new[] { (0L, 0L), (1L, 0L), (0L, 1L), (1L, 1L), (1L, -1L) } : new[] { (0L, 0L) })
                {
                    if (j == 0 && dm == 0) continue;
                    string mid = d == 1 ? "" : "I";
                    var windows = new[] { ("X" + mid + "X", j), ("Y" + mid + "Y", j), ("X" + mid + "Y", dm), ("Y" + mid + "X", -dm), ("Z" + mid + "Z", zz) }
                        .Where(t => t.Item2 != 0).Select(t => (t.Item1, Enumerable.Range(0, t.Item1.Length).ToArray(), t.Item2)).ToList();
                    var fields = new[] { ('X', a), ('Y', b) }.Where(f => f.Item2 != 0).ToList();
                    var cell = new Theodolite.Cell(new[] { (Zm, 1L) }, fields, windows);
                    string name = $"d={d} J={j} D={dm} ZZ={zz} f=({a},{b})";
                    var law = Stitch.RingLaw(cell);
                    Assert.True(law is not null && law.Text.Contains("Theorem 8"), $"{name}: {law?.Text}");
                    kinds.Add(fields.Count > 0 ? (dm == 0 ? "field, real hop" : "field, twisted hop")
                        : dm == 0 ? "order 1" : j == 0 ? "order 2" : Math.Abs(j) == Math.Abs(dm) ? "order 4" : "no order");
                    for (int n = Math.Max(law!.Guard, 3); n <= (d == 1 ? 7 : 8); n++)
                    {
                        var v = Sieve.Count(cell.Ring(n)).Verdict();
                        Assert.True(EndCount.IsExact(v), $"{name} ring {n}: not exact");
                        Assert.True(law.Pairs(n) == EndCount.IsPalindrome(v), $"{name} ring {n}: {law.Text}");
                    }
                }
        Assert.Equal(new[] { "field, real hop", "field, twisted hop", "no order", "order 1", "order 2", "order 4" },
            kinds.OrderBy(k => k, StringComparer.Ordinal));
        // outside the class: a pairing beside the hop (Theorem 9's), a Z beside the bond, two distances at once
        Theodolite.Cell Bonds(params (string, long)[] t) => new(new[] { (Zm, 1L) }, null, t.Select(x => (x.Item1, Enumerable.Range(0, x.Item1.Length).ToArray(), x.Item2)));
        foreach (var outside in new[] { Bonds(("XX", 2), ("YY", 1)), Bonds(("XXZ", 1), ("YYZ", 1)), Bonds(("XX", 1), ("YY", 1), ("XIX", 1), ("YIY", 1)) })
            Assert.True(Stitch.RingLaw(outside) is not { } l || !l.Text.Contains("Theorem 8"));
    }

    // F173's Corollary A′ (one sign per flip set, entries allowed to vanish), adopted in Stitch.RingLaw in both its cases, held
    // against the end count's exact verdict on the ring: every cell of two or three templates, coefficients ±1, from the
    // two-site flips with a Z beside or between (XXZ, YYZ, ZXX, ZYY, XZX, YZY, and XYZ, ZXY of the other sign) and the plain ones (XX, YY, XY, YX, XIX,
    // YIY, XIY, YIX) beside an even diagonal (ZZ, ZIZ), wherever RingLaw returns the corollary, on rings from its guard to
    // N = 8. All four outcomes occur. A cell whose lone excitation's hops cancel is read through a solid flip shape where
    // one exists (XIX − XXZ + ZXX, XIX solid, from N = 6 = the width 4 of XXZ and ZXX on their flip set plus the degree 2 of
    // 1 + x²), and not where none does (XXZ − YYZ, unbalanced at odd N and paired at every N, which Theorem 9 with Z
    // letters reads). Outside it: a flip set
    // carrying strings of both signs (XXZ + XZY share no flip set, so XX + XY's flip set is the control).
    [Fact]
    public void The_One_Sign_Law_Agrees_With_The_End_Count()
    {
        string[] pool = { "XXZ", "YYZ", "ZXX", "ZYY", "XZX", "YZY", "XYZ", "ZXY", "XX", "YY", "XY", "YX", "XIX", "YIY", "XIY", "YIX", "ZZ", "ZIZ" };
        var outcomes = new HashSet<string>();
        int read = 0;
        foreach (int k in new[] { 2, 3 })
            foreach (var combo in Combos(pool.Length, k))
                for (int signs = 0; signs < 1 << (k - 1); signs++)
                {
                    var terms = combo.Select((t, i) => (pool[t], i == 0 || ((signs >> (i - 1)) & 1) == 0 ? 1L : -1L)).ToArray();
                    var cell = Cell('Z', terms);
                    var law = Stitch.RingLaw(cell);
                    if (law is null || !law.Text.Contains("Corollary A′")) continue;
                    read++;
                    string pattern = "";
                    for (int n = law.Guard; n <= 8; n++)
                    {
                        var v = Sieve.Count(cell.Normalized().Ring(n)).Verdict();
                        string name = string.Join("+", terms);
                        Assert.True(EndCount.IsExact(v), $"{name} ring {n}: not exact");
                        Assert.True(law.Pairs(n) == EndCount.IsPalindrome(v), $"{name} ring {n}: {law.Text}");
                    }
                    outcomes.Add(law.Text.Contains("pairs on every") ? "every" : law.Text.Contains("> 1") ? "even" : law.Text.Contains("> 2") ? "four"
                        : law.Text.Contains("breaks on every") ? "never" : law.Text);
                }
        Assert.Equal(new[] { "even", "every", "four", "never" }, outcomes.OrderBy(o => o, StringComparer.Ordinal));
        Assert.True(read > 100, $"{read} cells read by Corollary A′");
        var solidRead = Stitch.RingLaw(Cell('Z', ("XIX", 1), ("XXZ", -1), ("ZXX", 1)));
        Assert.True(solidRead is { Guard: 6 } && solidRead.Text.Contains("solid flip shape") && solidRead.Text.Contains("> 1"), solidRead?.Text);
        // a flip set carrying strings of both signs (XX and XY), beside a second shape so that the corollary is asked
        Assert.True(Stitch.RingLaw(Cell('Z', ("XX", 1), ("XY", 1), ("XIX", 1))) is not { } mixed || !mixed.Text.Contains("Corollary A′"));
        // the pure pairing bond with a Z beside it: one sign, hops vanishing, unbalanced at odd N and paired at every N
        Assert.True(Stitch.RingLaw(Cell('Z', ("XXZ", 1), ("YYZ", -1))) is { } pure && pure.Text.Contains("Theorem 9") && pure.Text.Contains("pairs on every ring"));
    }

    // Corollary A′ through a solid flip shape does not need two-site flip sets: every cell of two or three templates,
    // coefficients ±1, from the three-site flips XXX, XYY, YXY, YYX (sign 0) and XXY, XYX, YXX, YYY (sign 1) beside the
    // two-site ones (XX, YY, XY, YX, XIX, YIY, XZX, YZY, XXZ, YYZ, ZXX, ZYY), wherever RingLaw reads it through a solid
    // shape, agrees with the end count's exact verdict on the rings from its guard to N = 8, cells whose flip shapes have
    // three sites among them, and so do cells of four and five templates with one or two shapes whose entries vanish somewhere beside
    // a solid one, one of them a solid shape of two strings (2·XXZ + ZXX, its entry ±2 ± 1). Two such shapes with no solid one (XXX − XYY beside XXZ − YYZ) are not read. The guard a + δ is
    // needed: XXX + YZZ + ZZY (XXX solid, a = 5 for the one-site flip set Y·ZZ from both sides, δ = 2) is unbalanced at
    // every N and pairs on the ring of six, where the one site off the run of a = 5 reaches two of the four cosets of 1 + x + x²; the
    // law reads it from N = 7, breaking there.
    [Fact]
    public void A_Solid_Flip_Shape_Reads_Vanishing_Flip_Sets_Of_Any_Width()
    {
        string[] pool = { "XXX", "XYY", "YXY", "YYX", "XXY", "XYX", "YXX", "YYY", "XX", "YY", "XY", "YX", "XIX", "YIY", "XZX", "YZY", "XXZ", "YYZ", "ZXX", "ZYY" };
        var outcomes = new HashSet<string>();
        int read = 0, wide = 0;
        void Check(params (string, long)[] terms)
        {
            var cell = Cell('Z', terms);
            var law = Stitch.RingLaw(cell);
            if (law is null || !law.Text.Contains("solid flip shape")) return;
            read++;
            if (terms.Any(t => t.Item1.Count(c => c is 'X' or 'Y') == 3)) wide++;
            string name = string.Join("+", terms);
            for (int n = law.Guard; n <= 8; n++)
            {
                var v = Sieve.Count(cell.Normalized().Ring(n)).Verdict();
                Assert.True(EndCount.IsExact(v), $"{name} ring {n}: not exact");
                Assert.True(law.Pairs(n) == EndCount.IsPalindrome(v), $"{name} ring {n}: {law.Text}");
            }
            outcomes.Add(law.Text.Contains("pairs on every") ? "every" : law.Text.Contains("> 1") ? "even" : law.Text.Contains("> 2") ? "four"
                : law.Text.Contains("breaks on every") ? "never" : law.Text);
        }
        foreach (int k in new[] { 2, 3 })
            foreach (var combo in Combos(pool.Length, k))
                for (int signs = 0; signs < 1 << (k - 1); signs++)
                    Check(combo.Select((t, i) => (pool[t], i == 0 || ((signs >> (i - 1)) & 1) == 0 ? 1L : -1L)).ToArray());
        int before = read;
        Check(("XIX", 1), ("YIY", -1), ("XXZ", 1), ("ZZ", 1));
        Check(("XXX", 1), ("XYY", -1), ("XIX", 1), ("YIY", -1), ("XX", 1));
        Check(("XXZ", 1), ("YYZ", -1), ("XIX", 1), ("YIY", -1), ("XXX", 1));
        Check(("XXZ", 1), ("ZXX", -1), ("XXX", 1), ("ZZ", 1));
        Check(("XXX", 1), ("YZZ", 1), ("ZZY", 1));
        Check(("XXZ", 2), ("ZXX", 1), ("XIX", 1), ("YIY", -1));                   // a solid shape of two strings
        Assert.Equal(6, read - before);
        var needsGuard = Cell('Z', ("XXX", 1), ("YZZ", 1), ("ZZY", 1));
        Assert.True(Stitch.RingLaw(needsGuard) is { Guard: 7 } g7 && g7.Text.Contains("breaks on every"), Stitch.RingLaw(needsGuard)?.Text);
        var six = Sieve.Count(needsGuard.Normalized().Ring(6)).Verdict();
        Assert.True(EndCount.IsExact(six) && EndCount.IsPalindrome(six), "XXX + YZZ + ZZY on the ring of six pairs, below the guard");
        // past N = 8, read whole by Theorem 2: at N = 9 the gcd with x⁹ + 1 is 1 + x + x² again, inside the condition
        // (N − a = 4), and the ring breaks; XIX − XXZ + ZXX breaks at 9 and pairs at 10
        Assert.True(Sieve.Count(needsGuard.Normalized().Ring(9)).ComplementConnection() is { AllGood: false }, "XXX + YZZ + ZZY ring 9");
        var hopless = Cell('Z', ("XIX", 1), ("XXZ", -1), ("ZXX", 1));
        Assert.True(Sieve.Count(hopless.Normalized().Ring(9)).ComplementConnection() is { AllGood: false }, "XIX − XXZ + ZXX ring 9");
        Assert.True(Sieve.Count(hopless.Normalized().Ring(10)).ComplementConnection() is { AllGood: true }, "XIX − XXZ + ZXX ring 10");
        Assert.Equal((1254, 1156), (read, wide));                                    // 1248 and 1152 of them from the pool
        Assert.Contains("every", outcomes);
        Assert.Contains("even", outcomes);
        Assert.Contains("never", outcomes);
        Assert.True(Stitch.RingLaw(Cell('Z', ("XXX", 1), ("XYY", -1), ("XXZ", 1), ("YYZ", -1))) is not { } none
            || !none.Text.Contains("Corollary A′"), "two shapes whose entries vanish somewhere and no solid one");
    }

    static IEnumerable<int[]> Combos(int n, int k)
    {
        if (k == 0) { yield return Array.Empty<int>(); yield break; }
        for (int first = 0; first <= n - k; first++)
            foreach (var rest in Combos(n - first - 1, k - 1))
                yield return new[] { first }.Concat(rest.Select(r => r + first + 1)).ToArray();
    }

    // Light along X or Y on every site is turned to Z by a proper rotation of the letters before the stitch's readers are
    // asked (X ↔ Z with Y → −Y; Y ↔ Z with X → −X): XZ + ZY under X light is ZX − XY under Z, broken for every N by the
    // stitch's no, where without the turn only the N read were broken. Every bond cell of one or two terms under X or Y
    // light reads, letter for letter, as its turned cell under Z.
    [Fact]
    public void Light_Along_X_Or_Y_Is_Turned_To_Z_For_The_Stitch()
    {
        var xl = Cell('X', ("XZ", 1), ("ZY", 1));
        var turned = Sieve.TurnedToZ(xl)!;
        Assert.Equal(new[] { ("ZX", 1L), ("XY", -1L) }, turned.Windows.Select(w => (w.Letters, w.Coefficient)));
        Assert.Equal(new[] { ("XY", -1L), ("YZ", 1L) }, Sieve.TurnedToZ(Cell('Y', ("XZ", 1), ("ZY", 1)))!.Windows.Select(w => (w.Letters, w.Coefficient)));
        Assert.Null(Sieve.TurnedToZ(Cell('Z', ("XZ", 1))));
        // fields turn with their signs: under X light Y → −Y and Z → X, under Y light X → −X and Z → Y
        Theodolite.Cell Fielded(char light, long y, long z) => new(new[] { (LetterMatrix(light), 1L) },
            new[] { ('Y', y), ('Z', z) }, new[] { ("XY", new[] { 0, 1 }, 1L), ("XZ", new[] { 0, 1 }, 1L) });
        Assert.Equal(new[] { ('Y', -1L), ('X', 1L) }, Sieve.TurnedToZ(Fielded('X', 1, 1))!.Fields.Select(f => (f.Letter, f.Coefficient)));
        Assert.Equal(new[] { ('Z', 1L), ('Y', 1L) }, Sieve.TurnedToZ(new Theodolite.Cell(new[] { (LetterMatrix('Y'), 1L) },
            new[] { ('Y', 1L), ('Z', 1L) }))!.Fields.Select(f => (f.Letter, f.Coefficient)));
        Assert.Equal(new[] { ('X', -1L) }, Sieve.TurnedToZ(new Theodolite.Cell(new[] { (LetterMatrix('Y'), 1L) },
            new[] { ('X', 1L) }))!.Fields.Select(f => (f.Letter, f.Coefficient)));
        // two field cells where a dropped field sign would give a wrong every-N answer each way: the turned stitch
        // agrees at every N read with the end count asked on the cell as it is (stitch: false never turns)
        static string Yes(string letters) => new(letters.Select(ch => ch is 'B' or 'b' or 'C' ? 'n' : 'y').ToArray());
        foreach (var (y, z, broken) in new[] { (1L, 1L, true), (-1L, 1L, false) })
        {
            var cell = Fielded('X', y, z);
            var turnedRun = Sieve.Run("XY+XZ|Y,Z", cell, 5);
            var bareRun = Sieve.Run("XY+XZ|Y,Z", cell, 5, stitch: false);
            Assert.Equal((Yes(bareRun.Chain.Letters), Yes(bareRun.Ring.Letters)), (Yes(turnedRun.Chain.Letters), Yes(turnedRun.Ring.Letters)));
            Assert.Equal(broken, turnedRun.Chain.BrokenEveryN && turnedRun.Ring.BrokenEveryN);
            Assert.Equal(!broken, turnedRun.Chain.Verdict == Sieve.Verdict.Read && turnedRun.Ring.Verdict == Sieve.Verdict.Read);
        }
        Assert.Null(Sieve.TurnedToZ(new Theodolite.Cell(new[] { (LetterMatrix('X'), 1L), (LetterMatrix('Y'), 1L) }, null,
            new[] { ("XZ", new[] { 0, 1 }, 1L) })));
        var e = Sieve.Run("XZ+ZY", xl, 6);
        Assert.True(e.Chain.BrokenEveryN);
        Assert.True(e.Ring.BrokenEveryN);
        Assert.Equal("BBBB", Sieve.Run("XZ+ZY", xl, 6, stitch: false).Ring.Letters);
        // the answers agree at every N read; the letters may differ in who answered, since the theodolite reads the cell
        // as it is and not turned: ten cells it reads only under Z, where the sieve answers through the stitch on the turned cell
        static string Pal(string letters) => new(letters.Select(ch => ch is 'B' or 'b' or 'C' ? 'n' : 'y').ToArray());
        var whoDiffers = new List<string>();
        foreach (char light in "XY")
            foreach (var (name, cell) in Sieve.BondFamily(2, light.ToString()))
            {
                var a = Sieve.Run(name, cell, 5);
                var b = Sieve.Run(name, Sieve.TurnedToZ(cell)!, 5);
                Assert.Equal((b.Chain.Verdict, Pal(b.Chain.Letters), b.Chain.BrokenFrom, b.Ring.Verdict, Pal(b.Ring.Letters), b.Ring.BrokenFrom, b.Ring.DecidedByLaw),
                    (a.Chain.Verdict, Pal(a.Chain.Letters), a.Chain.BrokenFrom, a.Ring.Verdict, Pal(a.Ring.Letters), a.Ring.BrokenFrom, a.Ring.DecidedByLaw));
                if (a.Chain.Letters != b.Chain.Letters || a.Ring.Letters != b.Ring.Letters)
                    whoDiffers.Add($"{light}:{name} {a.Chain.Letters}/{a.Ring.Letters} vs {b.Chain.Letters}/{b.Ring.Letters}");
            }
        Assert.Equal(new[]
        {
            "X:+XY+XZ", "X:+XY-XZ", "X:+YX+ZX", "X:+YX-ZX", "X:+YZ+ZY", "Y:+XY+ZY", "Y:+XY-ZY", "Y:+XZ+ZX", "Y:+YX+YZ", "Y:+YX-YZ",
        }, whoDiffers.Select(w => w.Split(' ')[0]));
        // where the theodolite reads the turned cell and not the cell under its own light, the stitch reads it; without
        // the turn those ten were left to the end count, a residue on some geometry
        foreach (var w in whoDiffers.Select(w => w.Split(' ')[0]))
        {
            var cell = Sieve.BondFamily(2, w[0].ToString()).Single(c => c.Name == w[2..]).Cell;
            var bare = Sieve.Run(w, cell, 5, stitch: false);
            Assert.Contains(Sieve.Verdict.Residue, new[] { bare.Chain.Verdict, bare.Ring.Verdict });
            var read = Sieve.Run(w, cell, 5);
            Assert.Equal((Sieve.Verdict.Read, Sieve.Verdict.Read), (read.Chain.Verdict, read.Ring.Verdict));
            // the theodolite covers more lengths on the turned cell than on the cell under its own light, on one geometry
            // at least (on the ring a closure can exist under both lights with different moduli, as for YZ + ZY)
            var turnedCell = Sieve.Run(w, Sieve.TurnedToZ(cell)!, 5);
            static int Es(string letters) => letters.Count(ch => ch == 'E');
            Assert.True(Es(turnedCell.Chain.Letters) > Es(read.Chain.Letters) || Es(turnedCell.Ring.Letters) > Es(read.Ring.Letters));
        }
    }

    // The window family's 24 ring cells that pair by a law in N, decided for every N by F173's ring laws: the eight relatives
    // of X·(X ± Y)·X by Corollary E2 (breaking at 3 | N from N = 6, the ring of three pairing), the sixteen of one string per flip set by Theorem D
    // (pairing exactly at even N); every exact letter from the guard on agrees with the law, or the sieve throws.
    [Fact]
    public void The_Ring_Laws_Decide_The_Window_Familys_Law_Cells()
    {
        var family = Sieve.WindowFamily(2, "Z").ToDictionary(c => c.Name, c => c.Cell);
        var e = Sieve.Run("+XXX+XYX", family["+XXX+XYX"], 9);
        Assert.Equal((Sieve.Verdict.Law, "PssBPFB"), (e.Ring.Verdict, e.Ring.Letters));
        Assert.Contains("Corollary E2", e.Ring.Law!.Text);
        e = Sieve.Run("+XXZ+XZY", family["+XXZ+XZY"], 9);
        Assert.Equal((Sieve.Verdict.Law, "BEBPBEB"), (e.Ring.Verdict, e.Ring.Letters));
        Assert.Contains("Theorem D", e.Ring.Law!.Text);
        // a tail the stitch's no already decides keeps its verdict
        e = Sieve.Run("+XIX+XIY", family["+XIX+XIY"], 7);
        Assert.Equal(Sieve.Verdict.Mixed, e.Ring.Verdict);
        Assert.NotNull(e.Ring.BrokenFrom);
    }

    // A Mixed with a tail: X·I·(X + Y) pairs on the chains of three and four sites (the stitch's yes at those N alone)
    // and breaks on every chain from five sites on, by the local no; XX·(X + Y) pairs at three and breaks from four. The
    // residue X·(X ± Y)·X has no tail on the chain.
    [Fact]
    public void A_Mixed_Cell_Is_Broken_For_Every_N_Past_Its_Last_Pairing()
    {
        var family = Sieve.WindowFamily(2, "Z").ToDictionary(c => c.Name, c => c.Cell);
        var e = Sieve.Run("+XIX+XIY", family["+XIX+XIY"], 6);
        Assert.Equal((Sieve.Verdict.Mixed, "ssCC", (int?)5), (e.Chain.Verdict, e.Chain.Letters, e.Chain.BrokenFrom));
        e = Sieve.Run("+XXX+XXY", family["+XXX+XXY"], 6);
        Assert.Equal((Sieve.Verdict.Mixed, "sCCC", (int?)4), (e.Chain.Verdict, e.Chain.Letters, e.Chain.BrokenFrom));
        Assert.False(e.Chain.BrokenEveryN);
        Assert.Null(Sieve.Run("+XXX+XYX", family["+XXX+XYX"], 6).Chain.BrokenFrom);
    }

    [Fact]
    public void On_The_Chain_Nothing_Lies_Between_The_Reading_And_The_End_Count()
    {
        var c = Census.Value;
        Assert.Equal(149, c.Count(e => e.Chain.Verdict == Sieve.Verdict.Read));
        Assert.Equal(268, c.Count(e => e.Chain.Verdict == Sieve.Verdict.Broken));
        Assert.Equal(0, c.Count(e => e.Chain.Verdict is Sieve.Verdict.Residue or Sieve.Verdict.Mixed));
        // every verdict the end count gave was exact, none a rank reading
        Assert.All(c, e => Assert.DoesNotContain(e.Chain.Letters + e.Ring.Letters, ch => ch is 'p' or 'b'));
    }

    [Fact]
    public void On_The_Ring_The_Residue_Is_The_XZ_ZY_Family_And_Its_Twin()
    {
        var c = Census.Value;
        Assert.Equal(149, c.Count(e => e.Ring.Verdict == Sieve.Verdict.Read));
        Assert.Equal(224, c.Count(e => e.Ring.Verdict == Sieve.Verdict.Broken));
        Assert.Equal(32, c.Count(e => e.Ring.Verdict == Sieve.Verdict.Mixed));
        Assert.Equal(RingResidue.OrderBy(s => s), c.Where(e => e.Ring.Verdict == Sieve.Verdict.Residue).Select(e => e.Name).OrderBy(s => s));
        // the residue breaks on every chain read, exactly
        Assert.All(c.Where(e => e.Ring.Verdict == Sieve.Verdict.Residue), e => Assert.Equal((Sieve.Verdict.Broken, "BBBB"), (e.Chain.Verdict, e.Chain.Letters)));
        // the mixed ones pair on the triangle alone
        Assert.All(c.Where(e => e.Ring.Verdict == Sieve.Verdict.Mixed), e => Assert.Equal("PBBB", e.Ring.Letters));
    }

    // The second pass: with the stitch asked where the theodolite is silent, the ring residue is read (S) at every N
    // read and the residue is empty; what stays mixed pairs on the triangle alone (a yes at one N, of the end count or of
    // the stitch, is no every-N reading); the chain's verdicts do not move, but its broken cells are now broken for every N.
    [Fact]
    public void With_The_Stitch_The_Ring_Residue_Is_Read_And_The_Residue_Is_Empty()
    {
        var c = Sieve.BondFamily(3, "Z").Select(x => Sieve.Run(x.Name, x.Cell, 6)).ToList();
        Assert.Equal((149, 268, 0, 0), (c.Count(e => e.Chain.Verdict == Sieve.Verdict.Read), c.Count(e => e.Chain.Verdict == Sieve.Verdict.Broken),
            c.Count(e => e.Chain.Verdict == Sieve.Verdict.Residue), c.Count(e => e.Chain.Verdict == Sieve.Verdict.Mixed)));
        Assert.Equal((161, 224, 0, 32), (c.Count(e => e.Ring.Verdict == Sieve.Verdict.Read), c.Count(e => e.Ring.Verdict == Sieve.Verdict.Broken),
            c.Count(e => e.Ring.Verdict == Sieve.Verdict.Residue), c.Count(e => e.Ring.Verdict == Sieve.Verdict.Mixed)));
        Assert.All(c.Where(e => RingResidue.Contains(e.Name)), e => Assert.Equal("SSSS", e.Ring.Letters));
        // the stitch's no: every broken cell, chain and ring, is broken for every N, by a bad local cycle from its
        // threshold on and by the end count, exactly, at every N read below it
        Assert.Equal(268, c.Count(e => e.Chain.BrokenEveryN));
        Assert.Equal(224, c.Count(e => e.Ring.BrokenEveryN));
        Assert.All(c.Where(e => RingResidue.Contains(e.Name)), e => Assert.Equal("CCCC", e.Chain.Letters));
        // the mixed ones pair on the triangle alone, at N = 3 read by the end count (P) or by the stitch at that one N (s),
        // and break for every N from 4 on by the stitch's no
        Assert.Equal(16, c.Count(e => e.Ring.Letters == "PCCC"));
        Assert.Equal(16, c.Count(e => e.Ring.Letters == "sCCC"));
    }

    // A proven second route: Theorem 2's corollary Theorem 1 (hopping bonds with a complement-invariant diagonal part,
    // letter fields, Z on every site, docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md) and F138's clause 2 (a field
    // along X + Y, orthogonal to the dephasing, which Theorem 1 leaves outside). In a uniform cell every site carries the
    // same field, so they read: the spectrum pairs exactly when the field has no Z component. The sieve decides all 64
    // cells (ZZ weight 0, ±1, 2; fields from {X, Y, Z} at magnitude 1 or 2) for every N, chain and ring, and agrees.
    [Fact]
    public void On_Hopping_Bonds_With_Fields_The_Sieve_Decides_Every_N_As_Theorem_1_Does()
    {
        int agree = 0;
        foreach (long zz in new[] { 0L, 1L, -1L, 2L })
            for (int mask = 0; mask < 8; mask++)
                foreach (long fc in new[] { 1L, 2L })
                {
                    var fields = new List<(char, long)>();
                    if ((mask & 1) != 0) fields.Add(('X', fc));
                    if ((mask & 2) != 0) fields.Add(('Y', fc));
                    if ((mask & 4) != 0) fields.Add(('Z', fc));
                    var windows = new List<(string, int[], long)> { ("XX", new[] { 0, 1 }, 1), ("YY", new[] { 0, 1 }, 1) };
                    if (zz != 0) windows.Add(("ZZ", new[] { 0, 1 }, zz));
                    var e = Sieve.Run($"ZZ {zz}, fields {mask} at {fc}", new Theodolite.Cell(new[] { (Zm, 1L) }, fields, windows), 6);
                    bool pairs = (mask & 4) == 0;
                    foreach (var p in new[] { e.Chain, e.Ring })
                    {
                        if (pairs) Assert.True(p.Verdict == Sieve.Verdict.Read && p.Letters.All(c => c is 'E' or 'S'), $"{e.Name}: {p.Letters}");
                        else Assert.True(p.Verdict == Sieve.Verdict.Broken && p.BrokenEveryN, $"{e.Name}: {p.Letters}");
                        agree++;
                    }
                }
        Assert.Equal(128, agree);
    }

    // What the sieve skipped, asked anyway at N = 3 and 4: every cell the reading certifies is palindromic there by the
    // end count, exactly, on the chain, and on the ring at those of the two N a route closes.
    [Fact]
    public void Every_Yes_Of_The_Reading_Is_A_Yes_Of_The_End_Count()
    {
        int asked = 0;
        foreach (var e in Census.Value)
        {
            var norm = e.Cell.Normalized();
            if (e.Chain.Verdict == Sieve.Verdict.Read)
                foreach (int n in new[] { 3, 4 })
                {
                    var v = Sieve.Count(norm.Chain(n)).Verdict();
                    Assert.True(EndCount.IsPalindrome(v) && EndCount.IsExact(v), $"{e.Name}, chain N = {n}: {v}");
                    asked++;
                }
            if (e.Ring.Traverse is { } t)
                foreach (int n in new[] { 3, 4 })
                    if (t.Closes(n))
                    {
                        var v = Sieve.Count(norm.Ring(n)).Verdict();
                        Assert.True(EndCount.IsPalindrome(v) && EndCount.IsExact(v), $"{e.Name}, ring N = {n}: {v}");
                        asked++;
                    }
        }
        Assert.True(asked > 2 * 149, $"only {asked} asked");
    }

    // Two cells that pair on a small chain and break on a longer one, and the reading stays silent on both:
    // SOFTNESS_IS_N_DEPENDENT's crossing pairs on four and five sites and breaks on six; its three-body control pairs on
    // three sites only.
    [Fact]
    public void A_Cell_That_Breaks_Past_A_Small_Chain_Is_Not_Read()
    {
        var crossing = new Theodolite.Cell(new[] { (Zm, 1L) }, null,
            new[] { ("XXXX", new[] { 0, 1, 2, 3 }, 1L), ("XYYY", new[] { 0, 1, 2, 3 }, 1L), ("YYYX", new[] { 0, 1, 2, 3 }, 1L) });
        var e = Sieve.Run("crossing", crossing, 6, stitch: false);
        Assert.Null(e.Chain.Traverse);
        Assert.Equal((Sieve.Verdict.Mixed, 4, "PPB"), (e.Chain.Verdict, e.Chain.From, e.Chain.Letters));
        var control = new Theodolite.Cell(new[] { (Zm, 1L) }, null,
            new[] { ("XXX", new[] { 0, 1, 2 }, 1L), ("XXY", new[] { 0, 1, 2 }, 1L), ("YXX", new[] { 0, 1, 2 }, 1L) });
        var c = Sieve.Run("control", control, 6, stitch: false);
        Assert.Equal((Sieve.Verdict.Mixed, "PBBB"), (c.Chain.Verdict, c.Chain.Letters));
    }

    // The residue read from below: for +XZ+ZY on every ring N = 3..6 no reader of the theodolite (site-local,
    // relabelled, turned frame, sectors; the Pi-reflection and the reflection) reads it at that N, the end count's far
    // element certifies the palindrome, and that element is no product of site operators: the far end is
    // one-dimensional, so every far element is a multiple of this one, and its coefficients, arranged by the letter on
    // site 0 against the rest, have rank two where a product has rank one. On the chain the same terms break, exactly,
    // by the count. The other eleven residue cells: the same element test at N = 3..5.
    [Fact]
    public void The_Ring_Residue_Pairs_By_An_Element_That_Is_No_Product()
    {
        var cell = new Theodolite.Cell(new[] { (Zm, 1L) }, null, new[] { ("XZ", new[] { 0, 1 }, 1L), ("ZY", new[] { 0, 1 }, 1L) });
        var family = Sieve.BondFamily(3, "Z").ToDictionary(c => c.Name, c => c.Cell);
        foreach (var name in RingResidue)
            for (int n = 3; n <= 5; n++)
            {
                var count = Sieve.Count(family[name].Normalized().Ring(n));
                Assert.Equal(EndCount.Reading.PalindromeByElement, count.Verdict());
                Assert.Equal((1, 1), count.UpperCounts());
                Assert.Equal(2, SchmidtRankAtSiteZero(count.FarElement()!.Value.Element));
            }
        var readers = new Func<Theodolite.Setting, Theodolite.Relation, Theodolite.Reading?>[]
            { (s, r) => Theodolite.Read(s, r), (s, r) => Theodolite.ReadRelabelled(s, r), (s, r) => Theodolite.ReadRotated(s, r),
              (s, r) => Theodolite.ReadSectors(s, r) };
        for (int n = 3; n <= 6; n++)
        {
            var ring = cell.Ring(n);
            foreach (var read in readers)
                foreach (var rel in new[] { Theodolite.Relation.Pi, Theodolite.Relation.Reflection })
                    Assert.Null(read(ring, rel));
            var count = Sieve.Count(ring);
            Assert.Equal(EndCount.Reading.PalindromeByElement, count.Verdict());
            var far = count.FarElement();
            Assert.NotNull(far);
            Assert.True(far!.Value.Check.Certifies);
            Assert.Equal(2, SchmidtRankAtSiteZero(far.Value.Element));
            Assert.Equal(EndCount.Reading.BrokenByCount, Sieve.Count(cell.Chain(n)).Verdict());
        }
        // at three sites the element is XXX + YYY, up to its scale
        var three = Sieve.Count(cell.Ring(3)).FarElement()!.Value.Element.OrderBy(t => t.Letters).ToList();
        Assert.Equal(new[] { "XXX", "YYY" }, three.Select(t => t.Letters));
        Assert.Equal(three[0].Coefficient, three[1].Coefficient);
    }

    // rank of the coefficient matrix (letter on site 0) x (letters on the other sites), over the rationals: at most
    // two rows here (X, Y), so rank two exactly when some 2x2 minor is nonzero
    static int SchmidtRankAtSiteZero(IReadOnlyList<(string Letters, long Coefficient)> element)
    {
        var rows = element.GroupBy(t => t.Letters[0]).ToDictionary(g => g.Key, g => g.ToDictionary(t => t.Letters[1..], t => (System.Numerics.BigInteger)t.Coefficient));
        if (rows.Count == 0) return 0;
        if (rows.Count == 1) return 1;
        var keys = rows.Keys.ToList();
        Assert.Equal(2, keys.Count);
        var cols = rows.Values.SelectMany(r => r.Keys).Distinct().ToList();
        System.Numerics.BigInteger At(char r, string c) => rows[r].TryGetValue(c, out var v) ? v : 0;
        for (int i = 0; i < cols.Count; i++)
            for (int j = i + 1; j < cols.Count; j++)
                if (At(keys[0], cols[i]) * At(keys[1], cols[j]) - At(keys[0], cols[j]) * At(keys[1], cols[i]) != 0) return 2;
        return 1;
    }

    // The twins are one object up to the sign of H: CZ on every bond of the ring is diagonal, so it commutes with the Z
    // dephasing, and it carries the ring Hamiltonian of XZ + s·ZY + u·ZZ onto that of ZX + s·YZ + u·ZZ term for term,
    // exactly, for s = ±1, u in {0, ±1}, at every N = 3..8 (the family's +YZ - ZX is minus ZX - YZ).
    [Fact]
    public void CZ_On_Every_Ring_Bond_Carries_XZ_ZY_Onto_ZX_YZ()
    {
        foreach (long s in new[] { 1L, -1L })
        foreach (long u in new[] { 0L, 1L, -1L })
        for (int n = 3; n <= 8; n++)
        {
            var from = new Theodolite.Cell(new[] { (Zm, 1L) }, null, new[] { ("XZ", new[] { 0, 1 }, 1L), ("ZY", new[] { 0, 1 }, s), ("ZZ", new[] { 0, 1 }, u) });
            var to = new Theodolite.Cell(new[] { (Zm, 1L) }, null, new[] { ("ZX", new[] { 0, 1 }, 1L), ("YZ", new[] { 0, 1 }, s), ("ZZ", new[] { 0, 1 }, u) });
            var image = from.Ring(n).H.Select(t => ConjugateByRingCz(t.Term, t.Coefficient, n))
                .ToDictionary(t => t.Item1, t => t.Item2);
            var target = to.Ring(n).H.ToDictionary(t => t.Term, t => t.Coefficient);
            Assert.Equal(target.Count, image.Count);
            foreach (var (str, c) in target) Assert.Equal(c, image[str]);
        }
    }

    // CZ on every ring bond: X_a -> Z_(a-1)·X_a·Z_(a+1), Y_a likewise, Z_a fixed; the string's image is the product of
    // its letters' images in site order, the phase tracked exactly (it must come out real)
    static (PauliString, long) ConjugateByRingCz(PauliString p, long coefficient, int n)
    {
        var acc = PauliString.Identity;
        int phase = 0;
        for (int a = 0; a < n; a++)
        {
            char letter = p.Letter(a);
            if (letter == 'I') continue;
            var letters = Enumerable.Repeat('I', n).ToArray();
            letters[a] = letter;
            var image = PauliString.Parse(new string(letters));
            if (letter is 'X' or 'Y')
                foreach (int b in new[] { (a + n - 1) % n, (a + 1) % n })
                {
                    var z = Enumerable.Repeat('I', n).ToArray();
                    z[b] = 'Z';
                    var (prod, ph) = PauliString.Multiply(image, PauliString.Parse(new string(z)));
                    image = prod;
                    phase += ph;
                }
            var (next, phNext) = PauliString.Multiply(acc, image);
            acc = next;
            phase += phNext;
        }
        phase = ((phase % 4) + 4) % 4;
        Assert.True(phase is 0 or 2, $"a complex phase i^{phase}");
        return (acc, phase == 0 ? coefficient : -coefficient);
    }
}
