using System.Linq;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Theorem F of F173: the torus read on the 2-part of its translations agrees, in both forms, with Theorem 2 walked on
/// the bitstrings, with the balance on the whole torus, and with Theorem D on the L × 1 torus and on the tori of odd L; the
/// three controls (one number per mask, the fold with T₁ and T₂ swapped, the fold onto the parities) disagree, so the battery
/// can fail.</summary>
public class TorusMaskBalanceWitnessTests
{
    [Fact]
    public void EveryFamilyAgrees_AndTheControlsDisagree()
    {
        var w = new TorusMaskBalanceWitness();
        Assert.All(w.Families, f => Assert.True(f.Agrees, $"{f.Name}: {f.FirstMismatch}"));
        Assert.True(w.ControlsBreak, $"controls {w.OneNumberMismatches}, {w.SwappedFoldMismatches}, {w.ParityFoldMismatches}");
        Assert.Equal(160, TorusMaskBalanceWitness.WindowTemplates.Count);
        Assert.Equal(new[] { (370, 1330, 851), (13180, 210880, 107725), (13180, 52720, 24181), (591, 16548, 8288), (7, 672, 222) },
            w.Families.Select(f => (f.Cells, f.Readings, f.Pairing)).ToArray());   // measured
        Assert.Equal((8436, 118620, 4866, 2328, 210880), (w.OneNumberMismatches, w.OneNumberReadings, w.SwappedFoldMismatches, w.ParityFoldMismatches, w.FoldControlReadings));   // measured, quoted on the page
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
        // the ring of N as the N × 1 torus: XIY breaks unless 4 | N (Theorem C)
        for (int n = 5; n <= 64; n++) Assert.Equal(n % 4 != 0, SignedMaskBalanceClaim.TorusBreaks(new[] { "XIY" }, n, 1));
        // the mask 1 + x on the 2 × M torus places the same flip set twice
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnTorus(new[] { "XX" }, 2, 3));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnTorus(new[] { "ZZ" }, 3, 3));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnTorus(new[] { "XX", "YX" }, 3, 3));
    }
}
