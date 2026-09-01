using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the warble (experiments/THE_CRACKED_BELL.md, 2026-08-31): crack the ring's wrap
// bond to J' = J(1-delta) and a launched traveling wave reverses its sense of circulation, fully, at
// T_rev = pi*N/(4*delta*J) + O(delta) -- the crack's O(delta) reflection accumulated over round trips.
// The circulation I(t) = sum_a Im rho[a, a+1 mod N]; R(t) = I(t)/I(0), 1 launched, -1 reversed.
// Pins computed by exact superoperator expm on the same Haken-Strobl block. Every one of them is
// now reproduced by the committed gates simulations/cracked_bell_gate.py (stages D, D5) and
// simulations/ring_renewal_current.py (2026-09-01): N=8, m=1, delta=0.15 -- gamma=0: T_zero=20.295, T_rev=40.6861 (the refined minimiser; a dt=0.02 grid reads 40.68),
// depth=-0.999949;
// gamma=0.01: T_zero=19.347, depth=-0.2752; gamma=0.05: T_zero=16.477. On THIS block the clock is
// gamma-DRESSED (the zero crossing advances) and the reversal outlives the naive scalar model
// e^(-4*gamma*t)*R_{gamma=0}(t) at its own best point (deepest-R ratio 1.251 at gamma=0.01, 3.84 at
// 0.05): the diagonal pays no dephasing and feeds the current back. The (0,1) page of the same crack
// is exactly gamma-free; which page is read decides what the watching does to the clock.
public class WarbleTests
{
    static readonly World W = new();

    // the expensive fixtures, one per gamma (N=8, m=1, delta=0.15, dt=0.02).
    static readonly Lazy<double[]> G000 = new(() => Warble.CurrentSeries(W, 8, 1.0, 0.00, m: 1, delta: 0.15, dt: 0.02, tMax: 58.0));
    static readonly Lazy<double[]> G001 = new(() => Warble.CurrentSeries(W, 8, 1.0, 0.01, m: 1, delta: 0.15, dt: 0.02, tMax: 58.0));
    static readonly Lazy<double[]> G005 = new(() => Warble.CurrentSeries(W, 8, 1.0, 0.05, m: 1, delta: 0.15, dt: 0.02, tMax: 58.0));

    // gamma = 0: the wave comes back whole. Committed: depth -0.999949 at T_rev = 40.6861. The first
    // order pi*N/(4*delta*J) is 41.888, so this reads -2.87% against it; against the EXACT pair split
    // (pi/Delta_E = 40.5932) it reads +0.23%, the (1,1) page's own O(delta^2) hand at delta = 0.15.
    // T_zero = 20.295.
    [Fact]
    public void The_Crack_Reverses_The_Wave_Fully()
    {
        double depth = Warble.ReversalDepth(G000.Value, out int kRev);
        Assert.InRange(depth, -1.0005, -0.995);
        Assert.InRange(kRev * 0.02, 40.29, 41.09);            // committed 40.6861 +/- 0.4
        Assert.InRange(Warble.ZeroCrossing(G000.Value, 0.02), 20.10, 20.50);  // committed 20.295
    }

    // the control: no crack, no warble -- the plane wave is an eigenstate, R(t) stays 1.
    [Fact]
    public void No_Crack_No_Reversal()
    {
        var series = Warble.CurrentSeries(W, 8, 1.0, 0.0, m: 1, delta: 0.0, dt: 0.02, tMax: 58.0);
        for (int k = 0; k < series.Length; k++)
            Assert.True(Math.Abs(series[k] / series[0] - 1.0) < 1e-6, $"R moved at k={k}: {series[k] / series[0]}");
    }

    // the watching DRESSES this clock (unlike the walk-time step): at gamma = 0.05 the zero crossing
    // has advanced from 20.295 to 16.477 -- a shift of 3.82, more than six times the +-0.30 window
    // asserted below (twelve times its half-width).
    [Fact]
    public void The_Watching_Dresses_The_Clock()
    {
        double tz0 = Warble.ZeroCrossing(G000.Value, 0.02);
        double tz5 = Warble.ZeroCrossing(G005.Value, 0.02);
        Assert.InRange(tz5, 16.18, 16.78);                    // committed 16.477
        Assert.True(tz0 - tz5 > 3.0, $"the dressing vanished: {tz0} vs {tz5}");
    }

    // the diagonal feeds the current back: at gamma = 0.01 the measured depth 0.2752 beats the naive
    // scalar model's own deepest point 0.2200 (ratio 1.251, exact-superoperator committed numbers);
    // the cheap in-test floor asserted here is the envelope at the gamma=0 clock, e^(-4*gamma*40.6861).
    [Fact]
    public void The_Diagonal_Feeds_The_Current_Back()
    {
        double depth = Math.Abs(Warble.ReversalDepth(G001.Value, out _));
        Assert.InRange(depth, 0.26, 0.29);                    // committed 0.2752
        Assert.True(depth > Math.Exp(-4.0 * 0.01 * 40.6861), "depth fell under the naive envelope");
    }

    // the crack is flat in mode space at first order, and the ratio the two clocks actually keep is
    // the EXACT split's, Split_1/Split_2 = 0.971754 off the quantization curve (THE_CRACKED_BELL
    // section E, gate E6b), the rest being this page's own O(delta^2) crossing deviation.
    // Committed (N=12, delta=0.1, gamma=0): T_zero 46.806 (m=1) vs 45.425 (m=2), ratio 0.9705.
    [Fact]
    public void Every_Mode_Hears_The_Same_Crack()
    {
        var s1 = Warble.CurrentSeries(W, 12, 1.0, 0.0, m: 1, delta: 0.10, dt: 0.02, tMax: 60.0);
        var s2 = Warble.CurrentSeries(W, 12, 1.0, 0.0, m: 2, delta: 0.10, dt: 0.02, tMax: 60.0);
        double ratio = Warble.ZeroCrossing(s2, 0.02) / Warble.ZeroCrossing(s1, 0.02);
        Assert.InRange(ratio, 0.95, 0.99);                    // committed 0.9705
    }

    // the knob: cracking the wrap bond keeps the trace.
    [Fact]
    public void The_Crack_Keeps_The_Trace()
    {
        var cone = new Cone(W, 8, 1.0, 0.05, Topology.Ring(8));
        cone.SetBond(7, 0, 1.0 * (1 - 0.15));
        cone.SeedPure(Warble.TravelingWave(8, 1));
        for (int k = 0; k < 50; k++) cone.Step(0.02);
        Assert.Equal(1.0, cone.Structure, 9);
    }
}
