using MirrorWorld;

namespace MirrorWorldTests;

// From-below test for Resolution: the individuation reading of the watching. The rate -2*gamma*k, read
// as a spectral LINEWIDTH (the boundary a line has), and the one verdict no other object produces --
// one object or two -- pinned to the already-adopted coherence horizon Q*(N) (Survivor's handover, the
// EP where two damped modes coalesce). Sober: the meaning (gamma as the object-maker, "boundaries make
// objects") lives in the docs; here only the numbers, recomputed and cross-checked against Pair.
public class ResolutionTests
{
    const double G = 0.5;
    static readonly World W = new();

    [Fact]
    public void Linewidth_Of_The_Diagonal_Is_Zero()
    {
        // k=0 (agreement, the diagonal): zero width, infinitely sharp -- the immortal, perfectly
        // individuated object. It is the |Re lambda| the Pair already carries at k=0.
        var r = new Resolution(W, n: 4, j: 1.0, gamma: G);
        Assert.Equal(0.0, r.Linewidth(0), 12);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void Linewidth_Is_The_Pairs_Rate_Read_As_A_Width(int k)
    {
        // the linewidth IS |Pair.Rate| = 2*gamma*k: the same number read as a boundary, not a price.
        var r = new Resolution(W, n: 5, j: 1.0, gamma: G);
        Assert.Equal(2.0 * G * k, r.Linewidth(k), 12);

        // cross-check against the Pair itself: a pair whose disagreement is exactly k.
        int j = (1 << k) - 1;                                   // k low bits set, popcount = k
        var p = new Pair(W, 0, j, G);
        Assert.Equal(k, p.Disagreement);
        Assert.Equal(Math.Abs(p.Rate), r.Linewidth(k), 12);
    }

    [Fact]
    public void Resolved_Flips_At_The_Coherence_Horizon()
    {
        // Q*(2) = 1.0 exactly: Survivor's handover, the EP where two damped modes coalesce. Below it the
        // mode pair is over-damped -- one merged object; above it, oscillating -- two distinct objects.
        Assert.Equal(Formulas.Qstar(2), new Resolution(W, 2, 1.0, G).Qstar, 12);

        var below = new Resolution(W, n: 2, j: 0.4, gamma: G);  // Q = 0.8 < 1.0
        var above = new Resolution(W, n: 2, j: 0.6, gamma: G);  // Q = 1.2 > 1.0
        Assert.False(below.Resolved);
        Assert.True(above.Resolved);
    }

    [Fact]
    public void Own_Is_Linewidth_And_Resolved_And_Inherits_The_Frame()
    {
        var r = new Resolution(W, 4, 1.0, G);
        Assert.Equal(new[] { "linewidth", "resolved" }, r.Own);
        Assert.Equal(new[] { "x", "y", "z" }, r.Inherited);     // inherited from the World
    }
}
