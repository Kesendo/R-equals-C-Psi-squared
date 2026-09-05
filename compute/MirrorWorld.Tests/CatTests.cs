using MirrorWorld;

namespace MirrorWorldTests;

// Schrodinger's cat as a named composition (Field at N): the k=N sighting of the same law the double slit
// is the k=1 sighting of. The two definite branches |0..0>,|1..1> are the immortal diagonal ("dead" and
// "alive"); the coherence between them |0..0><1..1> is the maximal-disagreement between (k=N) paying
// -2*gamma*N = -2*Sigma_gamma. Nothing new computed. Meaning in docs/quantum/SCHRODINGERS_CAT_TRANSLATED.md.
public class CatTests
{
    const double G = 0.05;
    static readonly World W = new();

    [Fact]
    public void Branches_Are_Immortal_The_Cat_Coherence_Dies()
    {
        var cat = new Cat(W, 4, G);
        double branches0 = cat.Branches, coh0 = cat.CatCoherence;
        Assert.True(branches0 == 2.0);             // two definite poles: |0000><0000| + |1111><1111| (sums of literal 1.0s)
        Assert.True(coh0 == 2.0);                  // the one cat coherence |0000><1111|, with its mirror twin
        for (int s = 0; s < 100; s++) cat.Watch(0.05);
        Assert.True(cat.Branches == branches0);    // the poles never move: the immortal diagonal (k=0) is never stepped
        Assert.True(cat.CatCoherence < coh0);      // the "both at once" dies
    }

    [Fact]
    public void The_Cat_Coherence_Decays_At_Minus_2NGamma()
    {
        var cat = new Cat(W, 4, G);
        Assert.True(cat.CoherenceRate == -2.0 * G * 4, "k=N=4: rate -2*gamma*N = -0.4, the product Pair forms");
        // Field.Step is forward Euler, so the cat coherence has the EXACT closed form (1 + rate*dt)^n under the
        // scheme; compare exactly (see DoubleSlitTests for the argument), then read the exp law as the
        // first-order limit: the deviation from e^{-2N gamma t} halves when dt halves (t = 2.5, one 1/e time;
        // e^{-1} * rate^2 t/2 * dt ~ 7.4e-5 at dt = 1e-3).
        double ratio1 = WatchToT(cat, dt: 1e-3, t: 2.5, out double euler1);
        Assert.True(ratio1 == euler1, $"the Euler product is exact: {ratio1:R} vs {euler1:R}");
        double exact = Math.Exp(-2.0 * G * 4 * 2.5);
        double ratio2 = WatchToT(new Cat(W, 4, G), dt: 5e-4, t: 2.5, out _);
        double dev1 = Math.Abs(ratio1 - exact), dev2 = Math.Abs(ratio2 - exact);
        double pred1 = exact * (2.0 * G * 4) * (2.0 * G * 4) * 2.5 / 2.0 * 1e-3;   // e^{rt} r^2 t dt / 2
        Assert.InRange(dev1 / pred1, 0.995, 1.005);  // the closed form, to the O(r dt) it leaves out
        Assert.InRange(dev1 / dev2, 1.95, 2.05);
    }

    static double WatchToT(Cat cat, double dt, double t, out double eulerProduct)
    {
        int n = (int)Math.Round(t / dt);
        double c0 = cat.CatCoherence, factor = 1.0 + cat.CoherenceRate * dt;
        eulerProduct = 1.0;
        for (int s = 0; s < n; s++) { cat.Watch(dt); eulerProduct *= factor; }
        return cat.CatCoherence / c0;
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(5)]
    public void The_Cat_Dies_N_Times_Faster_Than_The_Slit(int n)
    {
        // the k=N sighting is exactly N times the k=1 sighting: the bigger the superposition, the faster
        // it decoheres -- the same object read at the two ends of the one law.
        var cat = new Cat(W, n, G);
        double slitRate = new DoubleSlit(W, G).BetweenRate;   // -2*gamma (k=1)
        Assert.True(cat.CoherenceRate == n * slitRate);       // -2*gamma*N = N * (-2*gamma): the same doubles, exact
    }

    [Fact]
    public void Own_Reads_The_Phenomenon_And_Inherits_The_Frame()
    {
        var cat = new Cat(W, 4, G);
        Assert.Equal(new[] { "branches", "catCoherence" }, cat.Own);
        Assert.Equal(new[] { "x", "y", "z" }, cat.Inherited);
    }
}
