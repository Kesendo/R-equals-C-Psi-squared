using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the renewal cut (F126, adopted 2026-07-13 from
// docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md): the watched walk is the unwatched wave repeatedly
// caught and released, P_n(t) = e^{-Gamma t}(|G_n0|^2 + Gamma int sum |G_nm|^2 S_m). Renewal computes
// the WATCHED populations from purely UNWATCHED propagation plus bookkeeping: the dissipator is never
// stepped, only accounted. The cut is faithful exactly because F126 is proven; these tests pin it
// against the Cone engine the same way ConeTests pins Cone against Restless.
public class RenewalTests
{
    static readonly World W = new();

    // the cut is faithful: Renewal reproduces Cone's watched populations without ever stepping D.
    [Fact]
    public void Renewal_Agrees_With_Cone()
    {
        const int n = 15; const double j = 1.0, g = 0.3, tMax = 3.0;
        var renewal = new Renewal(W, n, j, g, seed: 7, dt: 0.01);
        var cone = new Cone(W, n, j, g);
        cone.Seed(7);
        for (int s = 0; s < 300; s++) cone.Step(0.01);
        var p = renewal.Populations(tMax);
        for (int a = 0; a < n; a++)
            Assert.True(Math.Abs(p[a] - cone.Population(a)) < 1e-3,
                $"site {a}: renewal {p[a]} vs cone {cone.Population(a)}");
    }

    // Gamma = 0: no catching, the renewal is the bare clean wave.
    [Fact]
    public void Renewal_Gamma0_Is_The_Clean_Wave()
    {
        const int n = 12; const double j = 1.0, tMax = 2.0;
        var renewal = new Renewal(W, n, j, 0.0, seed: 5, dt: 0.01);
        var cone = new Cone(W, n, j, 0.0);
        cone.Seed(5);
        for (int s = 0; s < 200; s++) cone.Step(0.01);
        var p = renewal.Populations(tMax);
        for (int a = 0; a < n; a++)
            Assert.True(Math.Abs(p[a] - cone.Population(a)) < 1e-6,
                $"site {a}: renewal {p[a]} vs cone {cone.Population(a)}");
    }

    // the watching moves weight around but loses none: the accounted populations stay a distribution.
    // The continuum ladder conserves exactly (the p = 0 pole of F126); the trapezoid grid conserves to
    // O(dt^2), so the pin is grid-limited and halving dt must shrink the drift by ~4 (checked below).
    [Fact]
    public void Renewal_Conserves_The_Excitation()
    {
        double Drift(double dt)
        {
            var renewal = new Renewal(W, 10, 1.0, 0.5, seed: 0, dt: dt);
            var p = renewal.Populations(2.5);
            double sum = 0; foreach (double x in p) sum += x;
            return Math.Abs(sum - 1.0);
        }
        double coarse = Drift(0.01), fine = Drift(0.005);
        Assert.True(coarse < 5e-4, $"coarse-grid drift {coarse}");
        Assert.True(fine < coarse / 3.0, $"drift not O(dt^2): {coarse} -> {fine}");
    }
}
