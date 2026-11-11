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
/// the battery can fail.</summary>
public class TorusMaskBalanceWitnessTests
{
    [Fact]
    public void EveryFamilyAgrees_AndTheControlsDisagree()
    {
        var w = new TorusMaskBalanceWitness();
        Assert.All(w.Families, f => Assert.True(f.Agrees, $"{f.Name}: {f.FirstMismatch}"));
        Assert.True(w.ControlsBreak, $"controls {w.OneNumberMismatches}, {w.SwappedFoldMismatches}, {w.ParityFoldMismatches}, {w.OneCornerMismatches}, {w.TighterBoundMismatches}, {w.MergedSitesMismatches}");
        Assert.Equal(160, TorusMaskBalanceWitness.WindowTemplates.Count);
        Assert.Equal(new[] { (370, 1330, 851), (13180, 210880, 107725), (13180, 52720, 24181), (591, 16548, 8288), (13180, 210880, 107725),
                (12588, 501933, 0), (596, 23866, 0), (522, 5047, 4211), (400, 400, 200), (600, 16800, 8935), (80200, 80200, 80200), (11, 1006, 342) },
            w.Families.Select(f => (f.Cells, f.Readings, f.Pairing)).ToArray());   // measured
        Assert.Equal((8436, 118620, 4866, 2328, 210880), (w.OneNumberMismatches, w.OneNumberReadings, w.SwappedFoldMismatches, w.ParityFoldMismatches, w.FoldControlReadings));   // measured, quoted on the page
        Assert.Equal((210, 210, 210880, 210), (w.OneCornerMismatches, w.OneCornerMismatchesOnThreeByThree, w.StaircaseReadings, w.CornerSetsWithSeveral));   // measured, quoted on the page
        Assert.Equal((1435, 454, 2216, 4949), (w.TighterBoundMismatches, w.MergedSitesMismatches, w.MergedSitesReadings, w.LatticeReadings));   // measured, quoted on the page
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
}
