using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live F121 witness: the assertions hold against the SPECTRUM the witness builds
/// at construction (not against the closed forms it then compares to). Each (d, N) case
/// materialises the dissipator and counts paired modes; the closed form is the cross-check.</summary>
public class QuditPartialPalindromeWitnessTests
{
    private static System.Collections.Generic.List<IInspectable> Children(QuditPartialPalindromeWitness w) =>
        ((IInspectable)w).Children.ToList();

    [Fact]
    public void Children_AreAllLive_TheTemplateRecomputesEveryRead()
    {
        var w = new QuditPartialPalindromeWitness(d: 2, n: 2);
        // The cited honesty template: every child reads the recomputed spectrum/census (the rate
        // law it checks against is referenced in prose, not surfaced as a banked node), so each
        // child is value-origin Live, not a frozen carrier's default Stored.
        Assert.All(Children(w), c => Assert.Equal(NodeProvenance.Live, c.Provenance));
    }

    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new QuditPartialPalindromeWitness(d: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new QuditPartialPalindromeWitness(n: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new QuditPartialPalindromeWitness(gamma: 0.0));
    }

    [Fact]
    public void Guard_RejectsOversizedGrid()
    {
        // d=4, N=3: d^(2N) = 4096 > MaxDim (1024).
        var ex = Assert.Throws<ArgumentOutOfRangeException>(() => new QuditPartialPalindromeWitness(d: 4, n: 3));
        Assert.Contains("4096", ex.Message);
    }

    [Fact]
    public void Qubit_d2_N2_FullPairing_16of16()
    {
        var w = new QuditPartialPalindromeWitness(d: 2, n: 2);
        Assert.Equal(16, w.Dim);
        Assert.Equal(16, w.PairedCount());                       // live count
        Assert.Equal(16L, QuditPartialPalindromeCeiling.Ceiling(2, 2)); // closed form == live
        Assert.Equal(0L, QuditProductMirrorCap.NonProductPart(2, 2));   // qubit: cap is full
    }

    [Fact]
    public void Qutrit_d3_N2_CeilingMet_54of81_Cap36_NonProduct18()
    {
        var w = new QuditPartialPalindromeWitness(d: 3, n: 2);
        Assert.Equal(81, w.Dim);
        Assert.Equal(54, w.PairedCount());                       // live spectrum pairs 54
        Assert.Equal(54L, QuditPartialPalindromeCeiling.Ceiling(3, 2));
        Assert.Equal(36L, QuditProductMirrorCap.ProductCap(3, 2));
        Assert.Equal(18L, QuditProductMirrorCap.NonProductPart(3, 2));
    }

    [Fact]
    public void Ququart_d4_N2_CeilingMet_128of256_Cap64_NonProduct64()
    {
        var w = new QuditPartialPalindromeWitness(d: 4, n: 2);
        Assert.Equal(256, w.Dim);
        Assert.Equal(128, w.PairedCount());
        Assert.Equal(128L, QuditPartialPalindromeCeiling.Ceiling(4, 2));
        Assert.Equal(64L, QuditProductMirrorCap.ProductCap(4, 2));
        Assert.Equal(64L, QuditProductMirrorCap.NonProductPart(4, 2));
    }

    [Fact]
    public void Qutrit_d3_N3_CeilingMet_378of729()
    {
        var w = new QuditPartialPalindromeWitness(d: 3, n: 3);
        Assert.Equal(729, w.Dim);
        Assert.Equal(378, w.PairedCount());
        Assert.Equal(378L, QuditPartialPalindromeCeiling.Ceiling(3, 3));
        Assert.Equal(216L, QuditProductMirrorCap.ProductCap(3, 3));
        Assert.Equal(162L, QuditProductMirrorCap.NonProductPart(3, 3));
    }

    [Theory]
    [InlineData(2, 2)]
    [InlineData(2, 3)]
    [InlineData(3, 2)]
    [InlineData(4, 2)]
    [InlineData(3, 3)]
    public void LivePairedCount_EqualsClosedFormCeiling_OnEveryGuardedGrid(int d, int n)
    {
        var w = new QuditPartialPalindromeWitness(d, n);
        long ceiling = QuditPartialPalindromeCeiling.Ceiling(d, n);
        Assert.Equal((int)ceiling, w.PairedCount());
        // d = 2 is the unique fully-paired column.
        Assert.Equal(d == 2, w.PairedCount() == w.Dim);
    }

    [Theory]
    [InlineData(2, 2)]
    [InlineData(3, 2)]
    [InlineData(4, 2)]
    [InlineData(3, 3)]
    public void LiveRungCounts_MatchTheMultiplicityFormula(int d, int n)
    {
        var w = new QuditPartialPalindromeWitness(d, n);
        var live = w.LiveRungCounts();
        for (int k = 0; k <= n; k++)
            Assert.Equal(QuditPartialPalindromeCeiling.Multiplicity(d, n, k), live[k]);
        Assert.Equal(QuditPartialPalindromeCeiling.Total(d, n), live.Sum());
    }

    [Fact]
    public void Spectrum_IsRealAndCentredOnTheRungLadder()
    {
        var w = new QuditPartialPalindromeWitness(d: 3, n: 2, gamma: 0.05);
        Assert.Equal(-2 * 0.05, w.Center, 12);  // -N*gamma = -0.1
        // The dissipator-only spectrum is exactly {0, -2γ, -4γ} with rung multiplicities.
        var diag = QuditPartialPalindromeWitness.BuildDissipator(3, 2, 0.05, out var ham);
        var rungs = ham.GroupBy(h => h).ToDictionary(g => g.Key, g => g.Count());
        Assert.Equal(9, rungs[0]);
        Assert.Equal(36, rungs[1]);
        Assert.Equal(36, rungs[2]);
        for (int x = 0; x < diag.RowCount; x++)
            Assert.Equal(0.0, diag[x, x].Imaginary, 12);
    }

    [Fact]
    public void Witness_SurfacesRungAndCapChildren()
    {
        var labels = Children(new QuditPartialPalindromeWitness(d: 3, n: 2)).Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("live count vs the closed form"));
        Assert.Contains(labels, l => l.Contains("rung k=0"));
        Assert.Contains(labels, l => l.Contains("rung k=2"));
        Assert.Contains(labels, l => l.Contains("cap split"));
        Assert.Contains(labels, l => l.Contains("full iff d=2"));
        Assert.Contains(labels, l => l.Contains("live spectrum"));
    }

    [Fact]
    public void Witness_RendersToJson()
    {
        var json = InspectionJsonExporter.ToJson(new QuditPartialPalindromeWitness(d: 3, n: 2));
        Assert.Contains("ceiling met", json);
        Assert.Contains("non-product", json);
    }
}
