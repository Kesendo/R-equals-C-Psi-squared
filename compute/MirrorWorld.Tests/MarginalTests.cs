using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// The local page (Marginal): the partial trace as an atom, lifted from the two inline kernels the F72
// and F73 pins already carried (SmokeTests). The page obeys F70 (PROOF_DELTA_N_SELECTION_RULE, Tier 1):
// a |S|-site page carries only |Delta popcount| <= |S| content. Pins are from-below on a generic cloud,
// with the discriminating contrast the F72 lesson demands: the content is really there, and really
// invisible to the smaller probe while the larger probe sees it.
public class MarginalTests
{
    const double G = 0.05;
    static readonly World W = new();

    // a deterministic generic Field cloud at N=3: every cell set, sign-mixed, all popcount sectors hit.
    static Field GenericCloud()
    {
        var f = new Field(W, 3, G);
        for (int i = 0; i < 8; i++)
            for (int j = i; j < 8; j++)
                f[i, j] = ((i * 3 + j * 7) % 11 - 5) / 5.0;
        return f;
    }

    [Fact]
    public void The_Full_Keep_Set_Is_The_Identity_Page_And_Keep_Order_Sets_The_Page_Bits()
    {
        var f = GenericCloud();
        var page = new Marginal(f, new[] { 0, 1, 2 });     // trace nothing: the identity member
        for (int a = 0; a < 8; a++)
            for (int b = 0; b < 8; b++)
                Assert.Equal(f[a, b], page[a, b].Real, 12);
        // the GIVEN keep order defines the page bits: swapping keep relabels the page. Page bit 0 of
        // the swapped page is site 1, so page index 1 (= bit 0 set) reads the full index |010> = 2.
        var swapped = new Marginal(f, new[] { 1, 0, 2 });
        Assert.Equal(f[2, 2], swapped[1, 1].Real, 12);
        Assert.Equal(f[1, 1], swapped[2, 2].Real, 12);
    }

    [Fact]
    public void F70_From_Below_One_Site_Blind_Past_One_Step_Two_Sites_Not()
    {
        var f = GenericCloud();
        int Dpc(int u, int v) => BitOperations.PopCount((uint)u) - BitOperations.PopCount((uint)v);

        // the contrast needs the |Dpc| >= 2 content to genuinely be there in the cloud
        double high = 0;
        for (int u = 0; u < 8; u++)
            for (int v = 0; v < 8; v++)
                if (Math.Abs(Dpc(u, v)) >= 2) high += Math.Abs(f[u, v]);
        Assert.True(high > 0.5, "the pin needs a cloud with real |Dpc| >= 2 content");

        // every one-site page IS the |Dpc| <= 1 content: the >= 2 sector traces to exactly zero (F70)
        for (int site = 0; site < 3; site++)
        {
            var page = new Marginal(f, new[] { site });
            for (int a = 0; a < 2; a++)
                for (int b = 0; b < 2; b++)
                {
                    var low = Complex.Zero;      // the <= 1 projected cloud, hand-traced onto (a, b)
                    for (int u = 0; u < 8; u++)
                        for (int v = 0; v < 8; v++)
                            if ((u & ~(1 << site)) == (v & ~(1 << site))
                                && ((u >> site) & 1) == a && ((v >> site) & 1) == b
                                && Math.Abs(Dpc(u, v)) <= 1)
                                low += f[u, v];
                    Assert.Equal(0.0, (page[a, b] - low).Magnitude, 15);
                }
        }

        // the two-site page over {0,1} DOES carry |Dpc| = 2 content: its (00, 11) corner is fed only
        // by pairs differing at both kept sites (Dpc = -2 there, forced), and that feed is nonzero.
        var pair = new Marginal(f, new[] { 0, 1 });
        var d2 = Complex.Zero;
        for (int u = 0; u < 8; u++)
            for (int v = 0; v < 8; v++)
                if ((u & ~3) == (v & ~3) && (u & 3) == 0b00 && (v & 3) == 0b11)
                {
                    Assert.Equal(-2, Dpc(u, v));               // the corner is pure |Dpc| = 2 sector
                    d2 += f[u, v];
                }
        Assert.True(d2.Magnitude > 0.05, "|Dpc| = 2 content must reach the two-site page");
        Assert.Equal(0.0, (pair[0b00, 0b11] - d2).Magnitude, 15);
    }

    [Fact]
    public void The_Live_Page_Over_The_Living_World_Matches_The_Hand_Trace()
    {
        // the Restless overload, pinned against the F73-style hand loop and the trace: the page is a
        // LIVE read, so it must follow every Step of the cloud it reduces.
        var w = new Restless(W, 4, 1.0, G);
        w.Seed(0b0000, 0.4);
        w.Seed(0b0011, 0.6);
        w.SeedCoherence(0b0000, 0b0001, 0.3);
        w.SeedCoherence(0b0000, 0b0011, 0.2);
        var page = new Marginal(w, new[] { 1 });
        var pair = new Marginal(w, new[] { 0, 1 });
        for (int tick = 0; tick <= 300; tick++)
        {
            if (tick % 100 == 0)
            {
                var m01 = Complex.Zero;          // (rho_site1)_{01} = sum_{s: bit 1 clear} rho[s, s | 2]
                for (int s = 0; s < 16; s++)
                    if (((s >> 1) & 1) == 0) m01 += w[s, s | 2];
                Assert.Equal(0.0, (page[0, 1] - m01).Magnitude, 12);
                // the partial trace preserves the trace, for any keep-set, at any time
                Assert.Equal(w.Structure, page.Structure, 12);
                Assert.Equal(w.Structure, pair.Structure, 12);
            }
            if (tick < 300) w.Step(0.01);
        }
    }

    [Fact]
    public void The_Keep_Set_Must_Be_Sane_And_Own_Reads_The_Page()
    {
        var f = GenericCloud();
        Assert.Throws<ArgumentException>(() => new Marginal(f, Array.Empty<int>()));
        Assert.Throws<ArgumentException>(() => new Marginal(f, new[] { 0, 0 }));
        Assert.Throws<ArgumentException>(() => new Marginal(f, new[] { 3 }));
        var page = new Marginal(f, new[] { 0 });
        Assert.Equal(2, page.PageDim);
        Assert.Equal(new[] { "page", "structure", "novelty" }, page.Own);
        // the page is a child of the cloud: it inherits the cloud's own, then the frame
        Assert.Equal(new[] { "structure", "novelty", "x", "y", "z" }, page.Inherited);
    }
}
