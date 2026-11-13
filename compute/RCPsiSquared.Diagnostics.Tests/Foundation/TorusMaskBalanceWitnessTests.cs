using System.Linq;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Theorem F of F173: the torus read on the 2-part of its translations agrees, in both forms, with Theorem 2 walked on
/// the bitstrings, with the balance on the whole torus, and with Theorem D on the L × 1 torus and on the tori of odd L; the
/// staircase's corners, one template's exponents coarsened to bit lengths, two templates' divisibility, the corner bound, the
/// honeycomb and kagome lattices and the size bound agree; the six controls (one number per mask, the fold with T₁ and T₂
/// swapped, the fold onto the parities, one corner, the bound one doubling tighter, the sites of a cell merged) disagree, so
/// the battery can fail. Theorem G on the open patch agrees with Theorem 2, the balance and Theorem B, and its two controls (no
/// gcd, the threshold along x alone) disagree, as do the threshold's two (over the near edge, the unreduced basis) and the one-torus
/// reading's three (the bound without the common factor's order, the side one doubling early, the binomial colouring one degree short).</summary>
public class TorusMaskBalanceWitnessTests
{
    [Fact]
    public void EveryFamilyAgrees_AndTheControlsDisagree()
    {
        var w = new TorusMaskBalanceWitness();
        Assert.All(w.Families, f => Assert.True(f.Agrees, $"{f.Name}: {f.FirstMismatch}"));
        Assert.True(w.ControlsBreak, $"controls {w.OneNumberMismatches}, {w.SwappedFoldMismatches}, {w.ParityFoldMismatches}, {w.OneCornerMismatches}, {w.TighterBoundMismatches}, {w.MergedSitesMismatches}, Theorem G {w.PatchNoGcdMismatches}, {w.PatchOneDirectionMismatches}, threshold {w.PatchThresholdOverhangMismatches}, {w.PatchVectorUnreducedMismatches}, one torus {w.BreakBoundNoOrderMismatches}, {w.BreakBoundOneEarlyMismatches}, {w.BreakBoundLowerDegreeMismatches}");
        Assert.Equal(160, TorusMaskBalanceWitness.WindowTemplates.Count);
        Assert.Equal(new[] { (370, 1330, 851), (13180, 210880, 107725), (13180, 52720, 24181), (591, 16548, 8288), (13180, 210880, 107725),
                (12588, 501933, 0), (596, 23866, 0), (522, 5047, 4211), (400, 400, 200), (600, 16800, 8935), (80200, 80200, 80200), (11, 1006, 342),
                (200, 2400, 2126), (13020, 833280, 644262), (400, 25600, 20502), (200, 200, 103), (591, 591, 306), (498, 31872, 23661), (300, 87000, 77672), (300, 3000, 2096), (300, 30000, 17416), (200, 62000, 58662), (200, 2000, 1634), (300, 300, 70), (600, 20928, 210), (2, 128, 86) },
            w.Families.Select(f => (f.Cells, f.Readings, f.Pairing)).ToArray());   // measured
        Assert.Equal((8436, 118620, 4866, 2328, 210880), (w.OneNumberMismatches, w.OneNumberReadings, w.SwappedFoldMismatches, w.ParityFoldMismatches, w.FoldControlReadings));   // measured, quoted on the page
        Assert.Equal((210, 210, 210880, 210), (w.OneCornerMismatches, w.OneCornerMismatchesOnThreeByThree, w.StaircaseReadings, w.CornerSetsWithSeveral));   // measured, quoted on the page
        Assert.Equal((1435, 454, 2216, 4949), (w.TighterBoundMismatches, w.MergedSitesMismatches, w.MergedSitesReadings, w.LatticeReadings));   // measured, quoted on the page
        Assert.Equal((11612, 28909, 833280, 5286, 274, 301), (w.PatchNoGcdMismatches, w.PatchOneDirectionMismatches, w.PatchControlReadings, w.PatchBreakingPairs, w.PatchPairsWhereTheGcdMatters, w.PatchTriplesWithSeveralCorners));   // measured, quoted on the page
        Assert.Equal((10185, 2096), (w.PatchThresholdOverhangMismatches, w.PatchThresholdNever));   // measured, quoted on the page
        Assert.Equal((554, 1200, 3000, 1919), (w.PatchVectorUnreducedMismatches, w.PatchVectorRunsCompared, w.PatchVectorBeyondRuns, w.PatchVectorMovedRunsCompared));   // measured, quoted on the page
        Assert.Equal((30, 43, 201, 8, 55, 108), (w.BreakBoundNoOrderMismatches, w.BreakBoundOneEarlyMismatches, w.BreakBoundSetsWithOrder, w.BreakBoundLargest, w.BreakBoundLowerDegreeMismatches, w.BreakBoundLowerDegreeReadings));   // measured, quoted on the page
    }

    // The page's examples read by the law far past the bitstrings, and the setting's refusals.
    [Fact]
    public void TheLawReadsThePagesExamples()
    {
        for (int l = 3; l <= 16; l++)
            for (int m = 3; l * m <= 64; m++)
            {
                Assert.Equal(m % 2 == 1, SignedMaskBalanceClaim.TorusBreaks(new[] { "XX", "X/Y" }, l, m));
                Assert.Equal(l % 2 == 1 || m % 2 == 1, SignedMaskBalanceClaim.TorusBreaks(new[] { "XX/XY" }, l, m));
                Assert.True(SignedMaskBalanceClaim.TorusBreaks(new[] { "XX", "X/Y", "IX/XI" }, l, m));
            }
        // the staircase's examples: two corners from one template, three, and two templates whose corner outruns their widths
        Assert.Equal(new[] { (0, 1), (1, 0) }, SignedMaskBalanceClaim.Corners(new[] { "IX/YI" }));
        Assert.Equal(new[] { (0, 2), (1, 1), (2, 0) }, SignedMaskBalanceClaim.Corners(new[] { "XXX/XX/Y" }));
        Assert.Equal(SignedMaskBalanceClaim.SingleTemplateCorners("XXX/XX/Y"), SignedMaskBalanceClaim.Corners(new[] { "XXX/XX/Y" }));
        Assert.Equal(new[] { (2, 3) }, SignedMaskBalanceClaim.Corners(new[] { "IX//X", "X//IY" }));
        Assert.Equal(new[] { (2, 4) }, SignedMaskBalanceClaim.Corners(new[] { "IX/X/X/X", "XIXX/X/X/Y" }));
        Assert.Equal(new[] { (2, 0) }, SignedMaskBalanceClaim.Corners(new[] { "XIY" }));
        Assert.Equal(new[] { (2, 0) }, SignedMaskBalanceClaim.Corners(new[] { "XXXY" }));   // exponent 3, bit length 2
        Assert.Equal((3, 5), SignedMaskBalanceClaim.SufficientTorus(new[] { "IX//X", "X//IY" }));
        // the corner bound: the colength and the least powers of u and v in 𝔪·I, met by the corners in both directions
        Assert.Equal((4, 3, 5), SignedMaskBalanceClaim.CornerBound(new[] { "IX//X", "X//IY" }));
        Assert.Equal((9, 4, 10), SignedMaskBalanceClaim.CornerBound(new[] { "IX/X/X/X", "XIXX/X/X/Y" }));
        Assert.Equal((1, 2, 2), SignedMaskBalanceClaim.CornerBound(new[] { "XX", "X/Y" }));
        Assert.Null(SignedMaskBalanceClaim.CornerBound(new[] { "XX/XY" }));   // one even mask: a principal ideal in 𝔪, no finite colength
        Assert.Equal((0, 1, 1), SignedMaskBalanceClaim.CornerBound(new[] { "XY/X" }));   // an odd mask is a unit: I = R
        Assert.Equal((7, 2, 8), SignedMaskBalanceClaim.CornerBound(new[] { "XX", "X/X/X/X/X/X/X/X" }, 8));   // u and v⁷: n_v = ℓ + 1 = n
        // the ring of N as the N × 1 torus: XIY breaks unless 4 | N (Theorem C)
        for (int n = 5; n <= 64; n++) Assert.Equal(n % 4 != 0, SignedMaskBalanceClaim.TorusBreaks(new[] { "XIY" }, n, 1));
        // the mask 1 + x on the 2 × M torus places the same flip set twice
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnTorus(new[] { "XX" }, 2, 3));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnTorus(new[] { "ZZ" }, 3, 3));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnTorus(new[] { "XX", "YX" }, 3, 3));
    }

    // Theorem G on the page's patches: the worked example's corner, a pair whose gcd moves it, the row M = 1, and the refusals.
    [Fact]
    public void TheoremGReadsThePagesPatches()
    {
        Assert.Equal((true, ((int, int)?)(3, 2)), SignedMaskBalanceClaim.Patch("XY", "XX/X"));
        Assert.Equal((false, ((int, int)?)null), SignedMaskBalanceClaim.Patch("XX", "X/Y"));
        // XY across and the plaquette XX/XY share the factor 1 + x: q₁ = 1, q₂ = 1 + y, c = 1·0 + 1·1 = 1, corner (2, 2);
        // without the gcd c = 1·0 + 1·0 = 0 and the control says it never breaks
        Assert.Equal((true, ((int, int)?)(2, 2)), SignedMaskBalanceClaim.Patch("XY", "XX/XY"));
        Assert.Equal((false, ((int, int)?)null), SignedMaskBalanceClaim.Patch("XY", "XX/XY", divideByGcd: false));
        Assert.True(SignedMaskBalanceClaim.PatchBreaks("XY", "XX/XY", 2, 2));
        Assert.False(TorusMaskBalanceWitness.PairsOnPatch(2, 2, "XY", "XX/XY"));
        // Theorem B is the row: +XIY + XXX breaks from N* = 5
        Assert.Equal((true, ((int, int)?)(5, 1)), SignedMaskBalanceClaim.Patch("XIY", "XXX"));
        // Theorem B's width counts identity letters at the end, the patch's extent does not
        Assert.Equal((true, ((int, int)?)(2, 1)), SignedMaskBalanceClaim.Patch("XI", "XY"));
        Assert.Equal((true, (int?)3), SignedMaskBalanceClaim.Chain("XI", "XY"));
        // three templates: two incomparable triangles, X ⊕ X = XY's mask along x and along y, the corners (2, 1) and (1, 2)
        for (int l = 1; l <= 4; l++)
            for (int m = 1; m <= 4; m++)
                Assert.Equal(l >= 2 || m >= 2, !SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnPatch(new[] { "X", "XY", "Y/X" }, l, m)));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.Patch("XX", "IYX"));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnPatch(new[] { "X/X", "IY/IX" }, 4, 4));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnPatch(new[] { "ZZ" }, 4, 4));
    }

    // Whether the patch breaks at some size, on the page's examples: the bound n₀ = ord(g) + ℓ′ and the one small torus
    // T₀ × T₀, T₀ the least power of 2 at least n₀ + 1. The parabola and its sibling (no common factor, n₀ = ℓ′ = 4 and 9, both pairing at every size: their
    // staircases have a corner); XX across beside XY down (ℓ′ = 1, pairs); Theorem G's worked example (breaks); XY across
    // beside the plaquette, the factor 1 + x shared (ord(g) = 1, the quotient 1 a unit, n₀ = 1, breaks); the three masks every two
    // of which share a factor, (1 + x)(1 + y), (1 + y)(x + y), (1 + x)(x + y) (ℓ′ = 3; their quadratic parts are independent, so
    // every relation has its coefficients in 𝔪 and never breaks, the torus 2 × 2 breaking and 4 × 4 pairing); one template, XIY,
    // mask (1 + x)², ord(g) = 2: never breaking on a patch, 4 × 4 pairing and the 2 × 2 of the early control breaking; and the
    // mask u⁵ + v, n₀ = 1, whose corner (3, 0), T₁ = 8, lies past the diagonal's T₀ = 2.
    [Fact]
    public void WhetherThePatchBreaksAtSomeSize()
    {
        Assert.Equal((0, 4, 4), SignedMaskBalanceClaim.PatchBreakBound(new[] { "IX//X", "X//IY" }));
        Assert.Equal((0, 9, 9), SignedMaskBalanceClaim.PatchBreakBound(new[] { "IX/X/X/X", "XIXX/X/X/Y" }));
        Assert.False(SignedMaskBalanceClaim.PatchBreaksAtSomeSize(new[] { "IX//X", "X//IY" }));
        Assert.False(SignedMaskBalanceClaim.PatchBreaksAtSomeSize(new[] { "IX/X/X/X", "XIXX/X/X/Y" }));
        Assert.Equal((0, 1, 1), SignedMaskBalanceClaim.PatchBreakBound(new[] { "XX", "X/Y" }));
        Assert.False(SignedMaskBalanceClaim.PatchBreaksAtSomeSize(new[] { "XX", "X/Y" }));
        Assert.True(SignedMaskBalanceClaim.PatchBreaksAtSomeSize(new[] { "XY", "XX/X" }));
        Assert.Equal((1, 0, 1), SignedMaskBalanceClaim.PatchBreakBound(new[] { "XY", "XX/XY" }));
        Assert.True(SignedMaskBalanceClaim.PatchBreaksAtSomeSize(new[] { "XY", "XX/XY" }));
        Assert.Equal((0, 3, 3), SignedMaskBalanceClaim.PatchBreakBound(new[] { "XX/XY", "IX/XX/X", "IXX/XX" }));
        Assert.False(SignedMaskBalanceClaim.PatchBreaksAtSomeSize(new[] { "XX/XY", "IX/XX/X", "IXX/XX" }));
        Assert.True(SignedMaskBalanceClaim.BreaksAtTwoParts(new[] { "XX/XY", "IX/XX/X", "IXX/XX" }, 2, 2));   // the bound met: 2 × 2 breaks, 4 × 4, the side n₀ + 1, pairs
        Assert.False(SignedMaskBalanceClaim.BreaksAtTwoParts(new[] { "XX/XY", "IX/XX/X", "IXX/XX" }, 4, 4));
        Assert.Equal((2, 0, 2), SignedMaskBalanceClaim.PatchBreakBound(new[] { "XIY" }));
        Assert.False(SignedMaskBalanceClaim.PatchBreaksAtSomeSize(new[] { "XIY" }));
        Assert.True(SignedMaskBalanceClaim.PatchBreaksAtSomeSize(new[] { "XIY" }, oneEarly: true));
        Assert.Equal((1, 0, 1), SignedMaskBalanceClaim.PatchBreakBound(new[] { "IXIIXY/X" }));
        Assert.Equal(new[] { (0, 1), (3, 0) }, SignedMaskBalanceClaim.Corners(new[] { "IXIIXY/X" }));
        // the binomial colourings: XX across beside XY down, the rows alternating, C(q, 1); XIY, XXYY along x, C(p, 2); none where
        // the patch breaks
        var rowsAlternate = SignedMaskBalanceClaim.BinomialColouring(new[] { "XX", "X/Y" })!;
        for (int p = 0; p < 4; p++)
            for (int q = 0; q < 4; q++) Assert.Equal(q & 1, SignedMaskBalanceClaim.BinomialColour(rowsAlternate, p, q));
        var xxyy = SignedMaskBalanceClaim.BinomialColouring(new[] { "XIY" })!;
        Assert.Equal(new[] { 0, 0, 1, 1, 0, 0, 1, 1 }, Enumerable.Range(0, 8).Select(p => SignedMaskBalanceClaim.BinomialColour(xxyy, p, 0)).ToArray());
        Assert.Null(SignedMaskBalanceClaim.BinomialColouring(new[] { "XY", "XX/X" }));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.PatchBreakBound(new[] { "XX", "IYX" }));
    }

    // The threshold at each height, every length at once, on the page's patches: Theorem G's quadrant read row by row and
    // column by column, the counterpart pairing at every height, the three-template set with two corners, Theorem B's row,
    // and the wrong boundary (templates hanging over the near edge onto φ = 0) moving a threshold.
    [Fact]
    public void ThePatchThresholdReadsThePagesPatches()
    {
        int? T(int m, bool alongY, params string[] ts) => SignedMaskBalanceClaim.PatchThreshold(ts, m, alongY);
        Assert.Null(T(1, false, "XY", "XX/X"));
        foreach (int m in new[] { 2, 3, 6 }) Assert.Equal(3, T(m, false, "XY", "XX/X"));
        Assert.Null(T(2, true, "XY", "XX/X"));
        foreach (int l in new[] { 3, 4, 7 }) Assert.Equal(2, T(l, true, "XY", "XX/X"));
        for (int m = 1; m <= 6; m++) { Assert.Null(T(m, false, "XX", "X/Y")); Assert.Null(T(m, true, "XX", "X/Y")); }
        Assert.Equal(2, T(1, false, "X", "XY", "Y/X"));
        Assert.Equal(1, T(2, false, "X", "XY", "Y/X"));
        Assert.Equal(2, T(1, true, "X", "XY", "Y/X"));
        Assert.Equal(1, T(2, true, "X", "XY", "Y/X"));
        Assert.Equal(5, T(1, false, "XIY", "XXX"));
        Assert.NotEqual(T(2, false, "XY", "XX/X"), SignedMaskBalanceClaim.PatchThreshold(new[] { "XY", "XX/X" }, 2, leftOverhang: true));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.PatchThreshold(new[] { "XX", "YX" }, 2));
    }

    // The same thresholds read off a reduced basis of the vector masks (Theorem B′ at height M), and at heights where W·M
    // exceeds the runs' 22 bits:
    // XY across beside XX/X at height 20 still breaks from length 3, along y from height 2 at width 20; the three-template
    // set's two corners; XIY + XXX at height 1, Theorem B′'s 5.
    [Fact]
    public void TheVectorMaskBasisReadsThePagesPatches()
    {
        int? B(int m, bool alongY, params string[] ts) => SignedMaskBalanceClaim.PatchThresholdByBasis(ts, m, alongY);
        Assert.Null(B(1, false, "XY", "XX/X"));
        foreach (int m in new[] { 2, 3, 6, 20 }) Assert.Equal(3, B(m, false, "XY", "XX/X"));
        Assert.Null(B(2, true, "XY", "XX/X"));
        foreach (int l in new[] { 3, 4, 7, 20 }) Assert.Equal(2, B(l, true, "XY", "XX/X"));
        for (int m = 1; m <= 20; m++) { Assert.Null(B(m, false, "XX", "X/Y")); Assert.Null(B(m, true, "XX", "X/Y")); }
        Assert.Equal((2, 1, 2, 1), (B(1, false, "X", "XY", "Y/X"), B(2, false, "X", "XY", "Y/X"), B(1, true, "X", "XY", "Y/X"), B(2, true, "X", "XY", "Y/X")));
        Assert.Equal(5, B(1, false, "XIY", "XXX"));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.PatchThresholdByBasis(new[] { "XX", "YX" }, 2));
    }
}
