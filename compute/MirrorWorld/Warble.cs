using System.Numerics;

namespace MirrorWorld;

// The warble (built 2026-08-31 beside experiments/THE_CRACKED_BELL.md, the bell-founder's word for the
// beat of a cracked bell): crack the ring's wrap bond to J' = J(1-delta) (the walk-time knob strengthens
// its bond, J(1+delta); the crack weakens -- mind the sign when porting delta) and the uniform ring's m <-> N-m
// traveling-wave pairs, degenerate before the crack, split by the SAME 4*delta*J/N at every m -- a point
// defect is flat in mode space. Launch the perfect ring's traveling wave on the cracked ring and the
// split partners beat: the wave's sense of circulation dies at T_zero ~ pi*N/(8*delta*J) and is fully
// REVERSED at T_rev ~ pi*N/(4*delta*J) (first order). This is the walk-time step's discarded O(delta)
// reflection, resonantly accumulated by the closed ring until it is the whole signal.
// The reading runs the Cone on Topology.Ring with the wrap bond cracked and reads the circulation
// I(t) = sum_a Im rho[a, a+1 mod N] (the sign is the sense of travel; R(t) = I(t)/I(0)).
// On THIS block (the (1,1) Haken-Strobl one) the clock is gamma-DRESSED: the zero crossing advances
// under the watching and the depth beats the naive e^(-4*gamma*t) envelope, because the diagonal pays
// no dephasing and feeds the current back. The (0,1) block of the same crack keeps the clock exactly
// gamma-free (there the dissipator is the scalar -2*gamma). What decides is not which block is
// read but whether the dissipator on it is a SCALAR, and then only for reads whose times are
// ZEROS: the (0,1) block's own beat deficit peaks at t = 1/(2*gamma), a time gamma sets
// entirely. Why the walk-time step's clock survives on THIS block is not settled; see
// experiments/THE_CRACKED_BELL.md, which withdrew the timescale answer it first gave.
public static class Warble
{
    // the perfect ring's traveling wave, momentum index m: psi(j) = exp(2*pi*i*m*j/N)/sqrt(N).
    public static Complex[] TravelingWave(int n, int m)
    {
        var psi = new Complex[n];
        for (int j = 0; j < n; j++)
            psi[j] = Complex.FromPolarCoordinates(1.0 / Math.Sqrt(n), 2.0 * Math.PI * m * j / n);
        return psi;
    }

    // record the circulation I(t) on the cracked ring: series[k] at t = k*dt.
    public static double[] CurrentSeries(World w, int n, double j, double gamma, int m, double delta,
                                         double dt, double tMax)
    {
        var cone = new Cone(w, n, j, gamma, Topology.Ring(n));
        cone.SetBond(n - 1, 0, j * (1.0 - delta));
        cone.SeedPure(TravelingWave(n, m));
        int steps = (int)Math.Round(tMax / dt);
        var series = new double[steps + 1];
        series[0] = Circulation(cone);
        for (int k = 1; k <= steps; k++) { cone.Step(dt); series[k] = Circulation(cone); }
        return series;
    }

    // the sense of travel: I = sum_a Im rho[a, a+1 mod N].
    public static double Circulation(Cone cone)
    {
        double s = 0;
        for (int a = 0; a < cone.Sites; a++) s += cone[a, (a + 1) % cone.Sites].Imaginary;
        return s;
    }

    // first zero crossing of I(t), linearly interpolated: the moment the circulation dies before turning.
    public static double ZeroCrossing(double[] series, double dt)
    {
        for (int k = 1; k < series.Length; k++)
        {
            double r0 = series[k - 1] / series[0], r1 = series[k] / series[0];
            if (r0 > 0 && r1 <= 0) return (k - 1) * dt + dt * r0 / (r0 - r1);
        }
        return double.NaN;
    }

    // the deepest reversal: min over t of R(t) = I(t)/I(0) (-1 = the wave came back whole).
    public static double ReversalDepth(double[] series, out int kRev)
    {
        double best = double.PositiveInfinity; kRev = 0;
        for (int k = 0; k < series.Length; k++)
        {
            double r = series[k] / series[0];
            if (r < best) { best = r; kRev = k; }
        }
        return best;
    }
}
