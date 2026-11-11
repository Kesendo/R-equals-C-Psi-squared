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
