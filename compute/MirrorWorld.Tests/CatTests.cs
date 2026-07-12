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
        Assert.Equal(2.0, branches0, 12);          // two definite poles: |0000><0000| + |1111><1111|
        Assert.Equal(2.0, coh0, 12);               // the one cat coherence |0000><1111|, with its mirror twin
        for (int s = 0; s < 100; s++) cat.Watch(0.05);
        Assert.Equal(branches0, cat.Branches, 12); // the poles never move: the immortal diagonal (k=0)
        Assert.True(cat.CatCoherence < coh0);      // the "both at once" dies
    }

    [Fact]
    public void The_Cat_Coherence_Decays_At_Minus_2NGamma()
    {
        var cat = new Cat(W, 4, G);
        Assert.Equal(-2.0 * G * 4, cat.CoherenceRate, 12);   // k=N=4: rate -2*gamma*N = -0.4
        double coh0 = cat.CatCoherence;
        const double dt = 1e-3; const int n = 2500;          // t = 2.5, the 1/e time 1/(2N*gamma)
        for (int s = 0; s < n; s++) cat.Watch(dt);
        Assert.Equal(Math.Exp(-2.0 * G * 4 * 2.5), cat.CatCoherence / coh0, 2);  // -> e^{-2N*gamma*t} = 1/e
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
        Assert.Equal(n * slitRate, cat.CoherenceRate, 12);    // -2*gamma*N = N * (-2*gamma)
    }

    [Fact]
    public void Own_Reads_The_Phenomenon_And_Inherits_The_Frame()
    {
        var cat = new Cat(W, 4, G);
        Assert.Equal(new[] { "branches", "catCoherence" }, cat.Own);
        Assert.Equal(new[] { "x", "y", "z" }, cat.Inherited);
    }
}
