namespace MirrorWorld;

// The walk-time reading (adopted 2026-07-13 from experiments/COUPLING_DEFECT_WALK_TIME_STEP.md): distance
// in this world is walk-time converted by J, and the conversion is locally additive. Deform one bond to
// J' = J(1+delta) and precisely that bond's walk-time changes, for every site beyond it, by the same amount:
// the arrival-time delay profile is a step, zero upstream, near -delta/(2J) downstream (first order). The
// reading runs two Cones, a defect-free reference and the defect run, and reads each site's arrival as the
// first threshold crossing, theta times the site's own defect-free peak (the relative threshold: the front
// peak decays ~ n^(-2/3) along the chain, an absolute threshold collides with that decay).
public static class WalkTime
{
    // record a run: series[k][site] is the population at t = k*dt, defect-free unless a defect is given.
    public static double[][] Reference(World w, int n, double j, double gamma, int seedSite, double dt,
                                       double tMax, (int a, int b)? defect = null, double delta = 0.0)
    {
        var cone = new Cone(w, n, j, gamma);
        if (defect is (int a, int b)) cone.SetBond(a, b, j * (1 + delta));
        cone.Seed(seedSite);
        int steps = (int)Math.Round(tMax / dt);
        var series = new double[steps + 1][];
        series[0] = Snapshot(cone);
        for (int k = 1; k <= steps; k++) { cone.Step(dt); series[k] = Snapshot(cone); }
        return series;
    }

    // the arrival-time delay profile: per site, defect arrival minus defect-free arrival, both read against
    // the defect-free threshold theta * max_t P_i. NaN where a site never crosses; 0 at the seed.
    public static double[] Profile(World w, int n, double j, double gamma, (int a, int b) defect, double delta,
                                   double[][] clean, int seedSite, double dt, double tMax, double theta = 0.2)
    {
        var def = Reference(w, n, j, gamma, seedSite, dt, tMax, defect, delta);
        var profile = new double[n];
        for (int i = 0; i < n; i++)
        {
            if (i == seedSite) { profile[i] = 0.0; continue; }
            double threshold = theta * Max(clean, i);
            double tClean = Crossing(clean, i, threshold, dt);
            double tDefect = Crossing(def, i, threshold, dt);
            profile[i] = double.IsNaN(tClean) || double.IsNaN(tDefect) ? double.NaN : tDefect - tClean;
        }
        return profile;
    }

    public static double Mean(double[] profile, int from, int to)
    {
        double s = 0; int c = 0;
        for (int i = from; i <= to; i++) if (!double.IsNaN(profile[i])) { s += profile[i]; c++; }
        return s / c;
    }

    static double[] Snapshot(Cone cone)
    {
        var p = new double[cone.Sites];
        for (int i = 0; i < cone.Sites; i++) p[i] = cone.Population(i);
        return p;
    }

    static double Max(double[][] series, int site)
    {
        double m = 0;
        foreach (var row in series) if (row[site] > m) m = row[site];
        return m;
    }

    // first crossing of the threshold, linearly interpolated on the dt grid.
    static double Crossing(double[][] series, int site, double threshold, double dt)
    {
        for (int k = 1; k < series.Length; k++)
        {
            double p0 = series[k - 1][site], p1 = series[k][site];
            if (p1 >= threshold && p0 < threshold)
                return (k - 1) * dt + dt * (threshold - p0) / (p1 - p0);
        }
        return double.NaN;
    }
}
