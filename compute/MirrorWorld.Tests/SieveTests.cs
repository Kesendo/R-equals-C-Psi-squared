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
    // X plus or minus Y), read by the free chain from six sites on, the stitch at three and the end count at four and
    // five below it: Read, and with it the window family's chain holds no residue.
    [Fact]
    public void The_Window_Residue_On_The_Chain_Is_Read_By_The_Free_Chain()
    {
        var family = Sieve.WindowFamily(2, "Z").ToDictionary(c => c.Name, c => c.Cell);
        foreach (var name in new[] { "+XXX+XYX", "+XXX-XYX", "+XXY+XYY", "+XXY-XYY", "+YXX+YYX", "+YXX-YYX", "+YXY+YYY", "+YXY-YYY" })
        {
            var e = Sieve.Run(name, family[name], 6);
            Assert.Equal((Sieve.Verdict.Read, "sPPF"), (e.Chain.Verdict, e.Chain.Letters));
        }
    }

    static Theodolite.Cell Cell(char light, params (string Letters, long C)[] terms)
        => new(new[] { (LetterMatrix(light), 1L) }, null, terms.Select(t => (t.Letters, Enumerable.Range(0, t.Letters.Length).ToArray(), t.C)));

    // A ring law that exists decides every N of a pass whose letters below its guard are exact, whatever the verdict:
    // XY on the ring is Read (the theodolite on even rings, broken at the odd N read) and Theorem D decides every odd
    // ring as well; XIY needs its guard, N = 5, among the N read or one past them; a cell without a law is decided by
    // none, and the chain has no ring law.
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
        var none = Sieve.Run("XY-YX", Cell('Z', ("XY", 1), ("YX", -1)), 6);          // pure DM: Theorem 8's, no ring law here
        Assert.Equal((Sieve.Verdict.Read, "BEBE", false), (none.Ring.Verdict, none.Ring.Letters, none.Ring.DecidedByLaw));
        Assert.Null(none.Ring.Law);
    }

    // On the first census (with the stitch) the read ring cells that carry a broken letter at an uncovered N are 13; a
    // ring law decides every N of 10 of them (one string per flip set by Theorem D; XX ± XY − YY and XX ± YX − YY by
    // Theorem 9, the pairing bond), and 3 stay decided only at the N read (the pure DM cells XY − YX, ± ZZ, Theorem 8's,
    // where RingLaw returns none).
    [Fact]
    public void On_The_First_Census_A_Ring_Law_Decides_Ten_Of_The_Thirteen_Read_Cells_With_Broken_Letters()
    {
        var c = Sieve.BondFamily(3, "Z").Select(x => Sieve.Run(x.Name, x.Cell, 6)).ToList();
        var broken = c.Where(e => e.Ring.Verdict == Sieve.Verdict.Read && e.Ring.BrokenFrom is null
            && e.Ring.Letters.Any(ch => ch is 'B' or 'b' or 'C')).ToList();
        Assert.Equal(13, broken.Count);
        var decided = broken.Where(e => e.Ring.DecidedByLaw).ToList();
        Assert.Equal(new[] { "+XX+XY-YY", "+XX+YX-YY", "+XX-XY-YY", "+XX-YX-YY", "+XY", "+XY+ZZ", "+XY-ZZ", "+YX", "+YX+ZZ", "+YX-ZZ" },
            decided.Select(e => e.Name).OrderBy(n => n, StringComparer.Ordinal));
        Assert.All(decided.Where(e => e.Name.StartsWith("+XX")), e => Assert.Contains("Theorem 9", e.Ring.Law!.Text));
        Assert.Equal(new[] { "+XY-YX", "+XY-YX+ZZ", "+XY-YX-ZZ" },
            broken.Where(e => !e.Ring.DecidedByLaw).Select(e => e.Name).OrderBy(n => n, StringComparer.Ordinal));
        Assert.All(broken.Where(e => !e.Ring.DecidedByLaw), e => Assert.Null(e.Ring.Law));
    }

    // Theorem 9 of the complement connection (the pairing bond), adopted in Stitch.RingLaw, held against the end count's
    // exact verdict on every bond J_XX·XX + J_YY·YY + J_XY·XY + J_YX·YX + J_ZZ·ZZ with coefficients in {−1, 0, 1, 2}: where
    // the pairing amplitude h is nonzero a ring law exists (F173's Theorem D or Corollary E2 where they read, Theorem 9
    // otherwise), and every one of them agrees with the end count on the ring of N = 3 to 7 sites; the theorem's chain
    // half (g real or imaginary) agrees with it on the open chain of 3 to 6 sites; where h = 0 the bond is Theorem 8's
    // and Theorem 9 is not asked. All three branches of Theorem 9 occur. Negative controls: bonds with a Z beside them, a
    // field beside the bond, are outside its class.
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
        // outside the class: a Z beside the bond, a field beside it
        Theodolite.Cell Window(params (string, long)[] t) => new(new[] { (Zm, 1L) }, null, t.Select(x => (x.Item1, new[] { 0, 1, 2 }, x.Item2)));
        foreach (var outside in new[]
        {
            Window(("XXZ", 1), ("YYZ", -1)), Window(("XXZ", 1), ("YYZ", 1)), Window(("ZXX", 1), ("ZYY", -1)),
            new Theodolite.Cell(new[] { (Zm, 1L) }, new[] { ('X', 1L) }, new[] { ("XX", new[] { 0, 1 }, 1L), ("YY", new[] { 0, 1 }, -1L) }),
        })
            Assert.True(Stitch.RingLaw(outside) is not { } l || !l.Text.Contains("Theorem 9"));
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
    // of X·(X ± Y)·X by Corollary E2 (breaking exactly at 3 | N), the sixteen of one string per flip set by Theorem D
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
