using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the renewal cut (F126, adopted 2026-07-13 from
// docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md): the walk in the light is the clean wave repeatedly
// caught and released, P_n(t) = e^{-Gamma t}(|G_n0|^2 + Gamma int sum |G_nm|^2 S_m). Renewal computes
// the populations IN THE LIGHT from purely CLEAN propagation plus bookkeeping: the dissipator is never
// stepped, only accounted. The cut is faithful exactly because F126 is proven; these tests pin it
// against the Cone engine the same way ConeTests pins Cone against Restless.
public class RenewalTests
{
    static readonly World W = new();

    // the cut is faithful: Renewal reproduces Cone's populations in the light without ever stepping D.
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

    // At gamma=0 the two-site chain is analytic: site 1 has sin^2(Jt), including off-grid t.
    [Fact]
    public void Renewal_Uses_The_Requested_Time_When_It_Is_Off_Grid()
    {
        const double tMax = 0.06;
        var renewal = new Renewal(W, n: 2, j: 1.0, gamma: 0.0, seed: 0, dt: 0.1);

        var p = renewal.Populations(tMax);
        double expected = Math.Pow(Math.Sin(tMax), 2);

        Assert.True(Math.Abs(p[1] - expected) < 1e-7,
            $"site 1 at t={tMax}: renewal {p[1]} vs exact {expected}");
    }

    // With J=0 the seed never moves; a dose that makes the implicit trapezoid denominator
    // near-singular must be refused instead of producing a population above one.
    [Fact]
    public void Renewal_Rejects_An_Underresolved_Refill_Dose()
    {
        var renewal = new Renewal(W, n: 2, j: 0.0, gamma: 0.5, seed: 0, dt: 0.9);

        Assert.Throws<ArgumentOutOfRangeException>(() => renewal.Populations(0.9));
    }

    // At J=0 the exact population stays one. The trapezoid instead multiplies its total
    // mass each step by exp(-dose)*(1+dose/2)/(1-dose/2); small local drift accumulates.
    [Fact]
    public void Renewal_Rejects_Accumulated_Refill_Drift()
    {
        var renewal = new Renewal(W, n: 2, j: 0.0, gamma: 0.5, seed: 0, dt: 0.2);

        Assert.Throws<ArgumentOutOfRangeException>(() => renewal.Populations(100.0));
    }

    // For the two-site clean chain, an RK4 step at J*h=3 has amplitudes -1/8 and +3i/2.
    // Their squared magnitudes sum to 2.265625, so the step is outside RK4 stability.
    [Fact]
    public void Renewal_Rejects_An_Unstable_Clean_Propagator_Step()
    {
        var renewal = new Renewal(W, n: 2, j: 1.0, gamma: 0.0, seed: 0, dt: 3.0);

        Assert.Throws<ArgumentOutOfRangeException>(() => renewal.Populations(3.0));
    }

    // RK4 stability prevents growth, not a large norm deficit: at J*h=2.8 the two-site
    // clean propagator retains only about 0.866 of the excitation after one step.
    [Fact]
    public void Renewal_Rejects_A_Large_Final_Mass_Deficit()
    {
        var renewal = new Renewal(W, n: 2, j: 1.0, gamma: 0.0, seed: 0, dt: 2.8);

        Assert.Throws<ArgumentOutOfRangeException>(() => renewal.Populations(2.8));
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-0.1)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void Renewal_Rejects_Nonpositive_Or_Nonfinite_Time_Step(double dt)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            new Renewal(W, n: 2, j: 1.0, gamma: 0.0, seed: 0, dt: dt));
    }

    [Theory]
    [InlineData(-0.1)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void Renewal_Rejects_Negative_Or_Nonfinite_Requested_Time(double tMax)
    {
        var renewal = new Renewal(W, n: 2, j: 1.0, gamma: 0.0, seed: 0, dt: 0.1);

        Assert.Throws<ArgumentOutOfRangeException>(() => renewal.Populations(tMax));
    }

    [Fact]
    public void Renewal_At_Zero_Time_Is_The_Seed_Population()
    {
        var renewal = new Renewal(W, n: 2, j: 1.0, gamma: 0.5, seed: 1, dt: 0.1);

        Assert.Equal(new[] { 0.0, 1.0 }, renewal.Populations(0.0));
    }

    // the light moves weight around but loses none: the accounted populations stay a distribution.
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
