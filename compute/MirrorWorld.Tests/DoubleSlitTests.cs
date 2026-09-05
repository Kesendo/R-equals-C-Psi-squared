using MirrorWorld;

namespace MirrorWorldTests;

// The double slit as a NAMED COMPOSITION of the atoms (Pair + Field at N=1): the access layer that makes
// the phenomenon recognizable where the scattered atoms were not (the missing "function/class that uses it
// composed"). From-below, nothing new computed: the humps are the immortal diagonal (constant), the fringe
// is the between |L><R| (the k=1 coherence) paying -2gamma and fading toward e^{-2gamma t}. Meaning in
// docs/quantum/DOUBLE_SLIT_TRANSLATED.md, not here.
public class DoubleSlitTests
{
    const double G = 0.05;                          // the code convention gamma_0 (a round number standing in for the hardware's ~5e3 Hz)
    static readonly World W = new();

    [Fact]
    public void Humps_Are_Immortal_The_Fringe_Pays()
    {
        var ds = new DoubleSlit(W, G);
        double humps0 = ds.Humps, fringe0 = ds.Fringe;
        Assert.True(humps0 == 2.0);                 // two places |L>,|R>: the two diagonal humps (sums of literal 1.0s)
        Assert.True(fringe0 == 2.0);                // the between |L><R|, counted with its mirror twin
        for (int s = 0; s < 200; s++) ds.Watch(0.05);   // t = 10, one coherence 1/e time
        Assert.True(ds.Humps == humps0);            // the humps never move: the immortal diagonal (k=0) is never stepped
        Assert.True(ds.Fringe < fringe0);           // the fringe fades: the watching's price
    }

    [Fact]
    public void The_Between_Decays_At_Minus_2Gamma_Toward_The_Exp_Law()
    {
        var ds = new DoubleSlit(W, G);
        Assert.True(ds.BetweenRate == -2.0 * G, "k=1: rate -2gamma, the product Pair forms, bit for bit");
        // Field.Step is forward Euler, w *= 1 + rate*dt, so under the scheme the fringe has an EXACT closed
        // form, (1 + rate*dt)^n, and fringe/fringe0 is that product bit for bit (2|w| over 2.0). Compare
        // exactly: a nonzero residual is a finding about Field.Step, never a tolerance to widen.
        double ratio1 = WatchToT(ds, dt: 1e-3, t: 10.0, out double euler1);
        Assert.True(ratio1 == euler1, $"the Euler product is exact: {ratio1:R} vs {euler1:R}");
        // e^{-2 gamma t} is the dt -> 0 limit, reached at FIRST order: the deviation from it shrinks by the
        // same factor as dt (e^{-1} * rate^2 t/2 * dt ~ 1.8e-5 at dt = 1e-3, t = 10). Read at dt and dt/2:
        // the ratio of the two deviations is 2. A wrong rate leaves both deviations O(1) (ratio ~1); a
        // second-order step would read 4; only the first-order Euler law on the right exponent reads 2.
        double exact = Math.Exp(-2.0 * G * 10.0);
        double ratio2 = WatchToT(new DoubleSlit(W, G), dt: 5e-4, t: 10.0, out _);
        double dev1 = Math.Abs(ratio1 - exact), dev2 = Math.Abs(ratio2 - exact);
        double pred1 = exact * (2.0 * G) * (2.0 * G) * 10.0 / 2.0 * 1e-3;   // the first-order term, e^{rt} r^2 t dt / 2
        Assert.InRange(dev1 / pred1, 0.995, 1.005);  // the closed form, to the O(r dt) it leaves out
        Assert.InRange(dev1 / dev2, 1.95, 2.05);
    }

    // watch to time t in steps of dt; also form the scheme's own closed form (1 + rate*dt)^n by the same
    // repeated multiplication Field.Step performs, so the comparison can be exact.
    static double WatchToT(DoubleSlit ds, double dt, double t, out double eulerProduct)
    {
        int n = (int)Math.Round(t / dt);
        double f0 = ds.Fringe, factor = 1.0 + ds.BetweenRate * dt;
        eulerProduct = 1.0;
        for (int s = 0; s < n; s++) { ds.Watch(dt); eulerProduct *= factor; }
        return ds.Fringe / f0;
    }

    [Fact]
    public void Visibility_Is_Unit_At_The_Balanced_Seed_And_Own_Reads_The_Phenomenon()
    {
        var ds = new DoubleSlit(W, G);
        Assert.True(ds.Visibility == 1.0);                   // V = 2|rho_LR|/(rho_LL+rho_RR) = 2/2 at the balanced seed, exact
        double factor = 1.0 + ds.BetweenRate * 0.05, product = 1.0;
        for (int s = 0; s < 200; s++) { ds.Watch(0.05); product *= factor; }   // t = 10
        Assert.True(ds.Visibility < 1.0);                    // the watching lowers V below its ceiling
        // the humps never move, so V(t) is the fringe's own Euler product, bit for bit (2w/2 = w)
        Assert.True(ds.Visibility == product, $"V(t) = {ds.Visibility:R} vs the Euler product {product:R}");
        Assert.Equal(new[] { "humps", "fringe", "visibility" }, ds.Own);
        Assert.Equal(new[] { "x", "y", "z" }, ds.Inherited); // the frame inherited from the World
    }
}
