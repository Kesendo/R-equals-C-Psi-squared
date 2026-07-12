using MirrorWorld;

namespace MirrorWorldTests;

// The double slit as a NAMED COMPOSITION of the atoms (Pair + Field at N=1): the access layer that makes
// the phenomenon recognizable where the scattered atoms were not (the missing "function/class that uses it
// composed"). From-below, nothing new computed: the humps are the immortal diagonal (constant), the fringe
// is the between |L><R| (the k=1 coherence) paying -2gamma and fading toward e^{-2gamma t}. Meaning in
// docs/quantum/DOUBLE_SLIT_TRANSLATED.md, not here.
public class DoubleSlitTests
{
    const double G = 0.05;                          // the canonical hardware-anchored watching rate
    static readonly World W = new();

    [Fact]
    public void Humps_Are_Immortal_The_Fringe_Pays()
    {
        var ds = new DoubleSlit(W, G);
        double humps0 = ds.Humps, fringe0 = ds.Fringe;
        Assert.Equal(2.0, humps0, 12);              // two places |L>,|R>: the two diagonal humps
        Assert.Equal(2.0, fringe0, 12);             // the between |L><R|, counted with its mirror twin
        for (int s = 0; s < 200; s++) ds.Watch(0.05);   // t = 10, one coherence 1/e time
        Assert.Equal(humps0, ds.Humps, 12);         // the humps never move: the immortal diagonal (k=0)
        Assert.True(ds.Fringe < fringe0);           // the fringe fades: the watching's price
    }

    [Fact]
    public void The_Between_Decays_At_Minus_2Gamma_Toward_The_Exp_Law()
    {
        var ds = new DoubleSlit(W, G);
        Assert.Equal(-2.0 * G, ds.BetweenRate, 12);         // k=1: rate -2gamma, the generator of the law
        double fringe0 = ds.Fringe;
        const double dt = 1e-3; const int n = 10000;        // t = 10 (the 1/e time), a clean 1/e fall
        for (int s = 0; s < n; s++) ds.Watch(dt);
        double ratio = ds.Fringe / fringe0;
        Assert.Equal(Math.Exp(-2.0 * G * 10.0), ratio, 2);  // |rho_LR(t)|/|rho_LR(0)| -> e^{-2gamma t} = 1/e
    }

    [Fact]
    public void Visibility_Is_Unit_At_The_Balanced_Seed_And_Own_Reads_The_Phenomenon()
    {
        var ds = new DoubleSlit(W, G);
        Assert.Equal(1.0, ds.Visibility, 12);                // V = 2|rho_LR|/(rho_LL+rho_RR) = 1 at the balanced seed
        for (int s = 0; s < 200; s++) ds.Watch(0.05);        // watch a while (t = 10)
        Assert.True(ds.Visibility < 1.0);                    // the watching lowers V below its ceiling
        Assert.Equal(ds.Fringe / ds.Humps, ds.Visibility, 12);  // V = fringe / the two humps
        Assert.Equal(new[] { "humps", "fringe", "visibility" }, ds.Own);
        Assert.Equal(new[] { "x", "y", "z" }, ds.Inherited); // the frame inherited from the World
    }
}
