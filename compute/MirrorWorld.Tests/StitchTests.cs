using System.Numerics;
using MirrorWorld;
using static MirrorWorld.LetterCube;

namespace MirrorWorldTests;

// From-below guard for the stitch (built 2026-11-10). A reading is judged by a route that never sees its equations:
// H built densely as Kronecker products of the Pauli matrices in Gaussian integers, and (★) H_x̄ȳ·d_y = d_x·H_xy checked
// entry by entry, cross-multiplied, at the reading's own rational point. Every yes it gives on the first census is asked
// of the end count; its every-N reading is held against the fixed-N solve past the window; the explicit Clifford of the
// first residue is checked as a matrix.
public class StitchTests
{
    static readonly Gauss[,] Zm = LetterMatrix('Z');

    static Theodolite.Cell Bond(params (string Letters, long C)[] terms)
        => new(new[] { (Zm, 1L) }, null, terms.Select(t => (t.Letters, new[] { 0, 1 }, t.C)));

    [Fact]
    public void The_Stitch_Owns_The_Stitch_And_The_Window_And_Inherits_The_Theodolite()
    {
        var stitch = new Stitch(new Theodolite(new LetterCube(new World())));
        Assert.Equal(new[] { "stitch", "window" }, stitch.Own);
        Assert.Equal(new[] { "letters", "reading", "every N", "coordinates", "height", "step", "arrow", "moves", "x", "y", "z" },
            stitch.Inherited);
    }

    // The sieve's first residue read for every N: the ring from three sites on, the chain at none of N = 2..5 and at no N
    // from the window on; and the reading's rational point satisfies (★) on the dense H of every ring N = 3..8.
    [Fact]
    public void The_Ring_Residue_Is_Read_For_Every_N_And_The_Chain_Is_Not()
    {
        foreach (var cell in new[] { Bond(("XZ", 1), ("ZY", 1)), Bond(("XZ", 1), ("ZY", -1), ("ZZ", 1)), Bond(("YZ", 1), ("ZX", 1)) })
        {
            var ring = Stitch.EveryN(cell, ring: true)!;
            Assert.NotNull(ring.Bulk);
            Assert.Equal(6, ring.Window);
            Assert.Equal(new[] { 3, 4, 5 }, ring.Small);
            var chain = Stitch.EveryN(cell, ring: false)!;
            Assert.Null(chain.Bulk);
            Assert.DoesNotContain(chain.Small, n => n >= 3);
            for (int n = 3; n <= 8; n++)
            {
                var red = Stitch.At(cell, n, ring: true);
                Assert.NotNull(red);
                var sol = Binomials.RationalSolution(red!);
                Assert.NotNull(sol);
                Assert.True(SatisfiesStar(cell.Normalized().Ring(n), sol!), $"N = {n}");
            }
        }
    }

    // The explicit carrier of the residue, W = X^⊗N·D with d_x = (s·i)^|x|·(−1)^(Σ x_l x_(l+1)) round the ring (S^⊗N·∏ CZ,
    // S† for s = −1), commutes with the dense H of XZ + s·ZY + u·ZZ on every ring N = 3..8, exactly; on the chain the same
    // D (the wrap pair left out) does not.
    [Fact]
    public void The_Residue_Is_Carried_By_X_Times_S_Times_CZ_On_Every_Ring_Bond()
    {
        foreach (long s in new[] { 1L, -1L })
            foreach (long u in new[] { 0L, 1L, -1L })
            {
                var cell = Bond(("XZ", 1), ("ZY", s), ("ZZ", u));
                for (int n = 3; n <= 8; n++)
                {
                    var unit = new Gauss(0, s);
                    Gauss D(ulong x, bool ring)
                    {
                        var d = Gauss.One;
                        for (int l = 0; l < n; l++) if (((x >> l) & 1) == 1) d = d * unit;
                        for (int l = 0; l < (ring ? n : n - 1); l++)
                            if (((x >> l) & 1) == 1 && ((x >> ((l + 1) % n)) & 1) == 1) d = -d;
                        return d;
                    }
                    Assert.True(CommutesWithXTimesDiagonal(cell.Normalized().Ring(n), x => D(x, true)), $"ring, s = {s}, u = {u}, N = {n}");
                    Assert.False(CommutesWithXTimesDiagonal(cell.Normalized().Chain(n), x => D(x, false)), $"chain, s = {s}, u = {u}, N = {n}");
                }
            }
    }

    // Every yes of the stitch on the first census, ring and chain, N = 3..5, is a yes of the end count, exactly.
    [Fact]
    public void Every_Yes_Of_The_Stitch_Is_A_Yes_Of_The_End_Count()
    {
        int yes = 0;
        foreach (var (name, cell) in Sieve.BondFamily(3, "Z"))
        {
            var norm = cell.Normalized();
            for (int n = 3; n <= 5; n++)
                foreach (bool ring in new[] { true, false })
                {
                    if (Stitch.At(norm, n, ring) is null) continue;
                    var v = Sieve.Count(ring ? norm.Ring(n) : norm.Chain(n)).Verdict();
                    Assert.True(EndCount.IsPalindrome(v) && EndCount.IsExact(v), $"{name}, {(ring ? "ring" : "chain")} N = {n}: {v}");
                    yes++;
                }
        }
        Assert.Equal(881, yes);                                          // measured: ring and chain, N = 3..5
    }

    // Past the window the every-N reading and the fixed-N solve agree, on cells it reads and cells it does not; the window
    // itself and one site past it are what EveryN solved and checked, so window + 2 is the new evidence.
    [Fact]
    public void Past_The_Window_The_Every_N_Reading_Is_The_Fixed_N_Solve()
    {
        foreach (var cell in new[] { Bond(("XZ", 1), ("ZY", 1)), Bond(("XX", 1), ("YY", 1), ("ZZ", 1)), Bond(("XX", 1), ("XY", 1)),
                                     Bond(("XZ", 1), ("YZ", 1), ("ZX", 1)) })
            foreach (bool ring in new[] { true, false })
            {
                var reach = Stitch.EveryN(cell, ring)!;
                for (int n = reach.Window; n <= reach.Window + 2; n++)
                    Assert.Equal(reach.Reads(n), Stitch.At(cell, n, ring) is not null);
            }
    }

    // Unequal moduli silence it: an X field of 3 with XZ + ZX gives single-flip ratios H_xy/H_x̄ȳ of 5, 1 and 1/5, which a
    // product of factors can absorb in one direction of (★) while the other direction fails; the cell is broken by the
    // end count (near 1, far 0) on rings N = 4..6, and the stitch is silent there and on every ring past its window.
    [Fact]
    public void Unequal_Moduli_On_An_Edge_Silence_The_Stitch()
    {
        var cell = new Theodolite.Cell(new[] { (Zm, 1L) }, new[] { ('X', 3L) },
            new[] { ("XZ", new[] { 0, 1 }, 1L), ("ZX", new[] { 0, 1 }, 1L) });
        for (int n = 4; n <= 6; n++)
        {
            Assert.Null(Stitch.At(cell, n, ring: true));
            var v = Sieve.Count(cell.Normalized().Ring(n)).Verdict();
            Assert.True(!EndCount.IsPalindrome(v) && EndCount.IsExact(v), $"N = {n}: {v}");
        }
        Assert.Null(Stitch.EveryN(cell, ring: true)!.Bulk);
    }

    [Fact]
    public void Jumps_Other_Than_Z_On_Every_Site_And_Windows_Past_The_Cap_Are_Not_Read()
    {
        var wide = new Theodolite.Cell(new[] { (Zm, 1L) }, null, new[] { ("XZX", new[] { 0, 1, 2 }, 1L) });
        Assert.Null(Stitch.EveryN(wide, ring: false));                  // chain window 18 > 16
        Assert.NotNull(Stitch.EveryN(wide, ring: true));                // ring window 9
        var lossy = new Theodolite.Cell(new[] { (Lowering, 1L) }, null, new[] { ("XX", new[] { 0, 1 }, 1L) });
        Assert.Null(Stitch.EveryN(lossy, ring: true));
        var xJump = new Theodolite.Cell(new[] { (LetterMatrix('X'), 1L) }, null, new[] { ("XZ", new[] { 0, 1 }, 1L), ("ZY", new[] { 0, 1 }, 1L) });
        Assert.Null(Stitch.At(xJump, 4, ring: true));
    }

    // ---- the dense route: H as Kronecker products of Pauli matrices, site 0 the lowest bit ----

    static Gauss[,] Dense(Theodolite.Setting s)
    {
        int dim = 1 << s.N;
        var h = new Gauss[dim, dim];
        foreach (var (term, c) in s.H)
            for (int col = 0; col < dim; col++)
            {
                // the column |col> through each site's letter matrix
                int row = 0;
                var amp = Gauss.One;
                for (int l = 0; l < s.N; l++)
                {
                    var m = LetterMatrix(term.Letter(l));
                    int b = (col >> l) & 1;
                    int r = m[0, b].IsZero ? 1 : 0;
                    amp = amp * m[r, b];
                    row |= r << l;
                }
                h[row, col] = h[row, col] + c * amp;
            }
        return h;
    }

    // (★) H_x̄ȳ·d_y = d_x·H_xy for all x, y, with d_x the product of the solution's factors, cross-multiplied; the pair
    // factors at offset 1 only, so for cells of span 2
    static bool SatisfiesStar(Theodolite.Setting s, IReadOnlyDictionary<string, (GaussBig Num, GaussBig Den)> sol)
    {
        int n = s.N, dim = 1 << n, all = dim - 1;
        var h = Dense(s);
        (GaussBig Num, GaussBig Den) D(int x)
        {
            var num = GaussBig.One; var den = GaussBig.One;
            void Times(string u) { var (a, b) = sol[u]; num *= a; den *= b; }
            for (int l = 0; l < n; l++) Times($"a{(x >> l) & 1}");
            for (int l = 0; l < n; l++) Times($"b1_{(x >> l) & 1}{(x >> ((l + 1) % n)) & 1}");
            return (num, den);
        }
        var d = Enumerable.Range(0, dim).Select(D).ToArray();
        for (int x = 0; x < dim; x++)
            for (int y = 0; y < dim; y++)
            {
                var left = GaussBig.Of(h[all ^ x, all ^ y]) * d[y].Num * d[x].Den;
                var right = GaussBig.Of(h[x, y]) * d[x].Num * d[y].Den;
                if (left != right) return false;
            }
        return true;
    }

    // W = X^⊗N·D against the dense H: W·H = H·W entry by entry
    static bool CommutesWithXTimesDiagonal(Theodolite.Setting s, Func<ulong, Gauss> d)
    {
        int dim = 1 << s.N, all = dim - 1;
        var h = Dense(s);
        // W = X^N·D sends |y> to d_y·|ȳ>, so (W·H)_(x̄, y) = d_x·H_xy and (H·W)_(x̄, y) = H_(x̄, ȳ)·d_y
        for (int x = 0; x < dim; x++)
            for (int y = 0; y < dim; y++)
                if (d((ulong)x) * h[x, y] != h[all ^ x, all ^ y] * d((ulong)y)) return false;
        return true;
    }
}
