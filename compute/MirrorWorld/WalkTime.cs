namespace MirrorWorld;

// The walk-time reading (adopted 2026-07-13 from experiments/COUPLING_DEFECT_WALK_TIME_STEP.md): distance
// in the measured ballistic-front regime is walk-time converted by J, and the conversion is locally
// additive. Deform one bond to J' = J(1+delta) and the observed downstream arrival-time delay is a step,
// near -delta/(2J) to first order. At strong watching the front is diffusive and this law is not claimed. The
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

    // ---- the inverse (2026-10-30, experiments/THE_DECODER_AND_THE_FRONT.md) ----
    // Read backwards, the step is a dictionary: one letter per bond, letter[b][i] = d(arrival at site i)/d(delta_b),
    // the central difference at +-delta of the profile above (the seed site reads 0). Near -1/(2J) beyond the bond
    // and near 0 before it, with a positive bump on the site just before the bond that grows with the bond's distance
    // from the seed: nested steps, of full rank, so several defects at once are read back by least squares over all
    // bonds (Decode). At theta = 0.2 no pair is anti-collinear through N = 10; past it the worst signed cosine drifts
    // down (about -0.54 at N = 20) and it depends on theta.
    public static double[][] Letters(World w, int n, double j, double gamma, int seedSite, double dt, double tMax,
                                     double delta = 0.02, double theta = 0.2)
    {
        var clean = Reference(w, n, j, gamma, seedSite, dt, tMax);
        var letters = new double[n - 1][];
        for (int b = 0; b < n - 1; b++)
        {
            var up = Profile(w, n, j, gamma, (b, b + 1), +delta, clean, seedSite, dt, tMax, theta);
            var down = Profile(w, n, j, gamma, (b, b + 1), -delta, clean, seedSite, dt, tMax, theta);
            letters[b] = new double[n];
            for (int i = 0; i < n; i++) letters[b][i] = (up[i] - down[i]) / (2 * delta);
        }
        return letters;
    }

    // the delay profile of several defect bonds at once, deltas[b] on bond (b, b+1), against the same clean run.
    public static double[] Profile(World w, int n, double j, double gamma, double[] deltas, double[][] clean,
                                   int seedSite, double dt, double tMax, double theta = 0.2)
    {
        if (deltas.Length != n - 1) throw new ArgumentException("one delta per bond of the chain", nameof(deltas));
        var cone = new Cone(w, n, j, gamma);
        for (int b = 0; b < n - 1; b++) if (deltas[b] != 0.0) cone.SetBond(b, b + 1, j * (1 + deltas[b]));
        cone.Seed(seedSite);
        int steps = (int)Math.Round(tMax / dt);
        var def = new double[steps + 1][];
        def[0] = Snapshot(cone);
        for (int k = 1; k <= steps; k++) { cone.Step(dt); def[k] = Snapshot(cone); }
        var profile = new double[n];
        for (int i = 0; i < n; i++)
        {
            if (i == seedSite) continue;
            double threshold = theta * Max(clean, i);
            double tClean = Crossing(clean, i, threshold, dt), tDefect = Crossing(def, i, threshold, dt);
            profile[i] = double.IsNaN(tClean) || double.IsNaN(tDefect) ? double.NaN : tDefect - tClean;
        }
        return profile;
    }

    // read every bond at once: the least-squares deltas x with sum_b x_b * letters[b] closest to the profile
    // (the normal equations, solved by Gauss-Jordan elimination with partial pivoting; NaN sites are skipped).
    public static double[] Decode(double[][] letters, double[] profile)
    {
        Check(letters);
        int m = letters.Length;
        var g = new double[m, m + 1];
        for (int a = 0; a < m; a++)
        {
            for (int b = 0; b < m; b++) g[a, b] = Dot(letters[a], letters[b], profile);
            g[a, m] = Dot(letters[a], profile, profile);
        }
        for (int c = 0; c < m; c++)
        {
            int p = c;
            for (int r = c + 1; r < m; r++) if (Math.Abs(g[r, c]) > Math.Abs(g[p, c])) p = r;
            for (int k = 0; k <= m; k++) (g[c, k], g[p, k]) = (g[p, k], g[c, k]);
            for (int r = 0; r < m; r++)
            {
                if (r == c) continue;
                double f = g[r, c] / g[c, c];
                for (int k = c; k <= m; k++) g[r, k] -= f * g[c, k];
            }
        }
        var x = new double[m];
        for (int a = 0; a < m; a++) x[a] = g[a, m] / g[a, a];
        return x;
    }

    // the DefectDecoder's rule on these letters: project the profile on each letter alone, the smallest squared
    // residual wins; the runner-up's residual over the winner's is the confidence (the decoder flags below 3).
    public static (int Bond, double Strength, int RunnerUp, double Ratio) DecodeOne(double[][] letters, double[] profile)
    {
        Check(letters);
        var scored = new List<(double Residual, int Bond, double Strength)>();
        for (int b = 0; b < letters.Length; b++)
        {
            double s = Dot(letters[b], profile, profile) / Dot(letters[b], letters[b], profile);
            double r = 0;
            for (int i = 0; i < profile.Length; i++)
                if (!double.IsNaN(profile[i])) { double e = profile[i] - s * letters[b][i]; r += e * e; }
            scored.Add((r, b, s));
        }
        scored.Sort((u, v) => u.Residual.CompareTo(v.Residual));
        return (scored[0].Bond, scored[0].Strength, scored[1].Bond, scored[1].Residual / scored[0].Residual);
    }

    // the most anti-collinear pair of letters: the minimum signed cosine (-1 is the sign-location confusion). It does
    // not see two letters that are nearly the same (the location blindness); Conditioning does. Throws on a letter
    // that is not finite (a site that never crossed) or zero, so a silent pair cannot pass for a clean one.
    public static (double Cos, int A, int B) WorstSignedCos(double[][] letters)
    {
        var c = Cosines(letters);
        (double Cos, int A, int B) worst = (double.PositiveInfinity, -1, -1);
        for (int a = 0; a < letters.Length; a++)
            for (int b = a + 1; b < letters.Length; b++)
                if (c[a, b] < worst.Cos) worst = (c[a, b], a, b);
        return worst;
    }

    // the smallest singular value of the dictionary with every letter normalised: 0 when some combination of letters
    // vanishes (a direction no reading sees), 1 for orthogonal letters. The square root of the smallest eigenvalue of
    // the cosine matrix, by cyclic Jacobi rotations (the matrix is (N-1) x (N-1)).
    public static double Conditioning(double[][] letters) => Math.Sqrt(Math.Max(MinEigenvalue(Cosines(letters)), 0.0));

    // the smallest eigenvalue of a real symmetric matrix, by cyclic Jacobi rotations (the matrix is copied).
    public static double MinEigenvalue(double[,] a)
    {
        int m = a.GetLength(0);
        var g = (double[,])a.Clone();
        for (int sweep = 0; sweep < 100; sweep++)
        {
            double off = 0;
            for (int p = 0; p < m; p++) for (int q = p + 1; q < m; q++) off += g[p, q] * g[p, q];
            if (off < 1e-30) break;
            for (int p = 0; p < m; p++)
                for (int q = p + 1; q < m; q++)
                {
                    if (g[p, q] == 0.0) continue;
                    double theta = (g[q, q] - g[p, p]) / (2 * g[p, q]);
                    double t = Math.Sign(theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
                    if (theta == 0.0) t = 1.0;
                    double cs = 1 / Math.Sqrt(t * t + 1), sn = t * cs;
                    for (int k = 0; k < m; k++)
                    {
                        double gkp = g[k, p], gkq = g[k, q];
                        g[k, p] = cs * gkp - sn * gkq; g[k, q] = sn * gkp + cs * gkq;
                    }
                    for (int k = 0; k < m; k++)
                    {
                        double gpk = g[p, k], gqk = g[q, k];
                        g[p, k] = cs * gpk - sn * gqk; g[q, k] = sn * gpk + cs * gqk;
                    }
                }
        }
        double min = double.PositiveInfinity;
        for (int k = 0; k < m; k++) min = Math.Min(min, g[k, k]);
        return min;
    }

    // The same readings at equal cost (experiments/THE_DECODER_AND_THE_FRONT.md §5): one shot is one preparation of
    // the carrier and one Z measurement of every site at one time, and it returns the excitation's position j with
    // probability p_j(t). The per-shot Fisher matrix of every bond's relative delta, averaged over the schedule's
    // times (equal shots per time), is F_ab = mean_t sum_j dp_j/da dp_j/db / p_j, the central difference at
    // +-delta for the slopes. The carrier and the cone's step are as in TrajectoryLetters; every time must be a
    // multiple of dt.
    public static double[,] PositionFisher(World w, int n, double j, double gamma,
                                           (double Weight, double[] Wave)[] carrier,
                                           double[] times, double dt, double delta = 1e-5)
    {
        var steps = times.Select(t => (int)Math.Round(t / dt)).ToArray();
        for (int k = 0; k < times.Length; k++)
            if (Math.Abs(steps[k] * dt - times[k]) > 1e-9 * Math.Max(1.0, times[k]) || (k > 0 && steps[k] <= steps[k - 1]))
                throw new ArgumentException("times must rise and be multiples of dt", nameof(times));
        double[][] Run(int bond, double d)
        {
            var p = new double[times.Length][];
            for (int k = 0; k < times.Length; k++) p[k] = new double[n];
            foreach (var (weight, wave) in carrier)
            {
                var cone = new Cone(w, n, j, gamma);
                if (bond >= 0) cone.SetBond(bond, bond + 1, j * (1 + d));
                cone.SeedPure(wave);
                int done = 0;
                for (int k = 0; k < times.Length; k++)
                {
                    for (; done < steps[k]; done++) cone.Step(dt);
                    for (int i = 0; i < n; i++) p[k][i] += weight * cone.Population(i);
                }
            }
            return p;
        }
        var p0 = Run(-1, 0.0);
        var slope = new double[n - 1][][];
        for (int b = 0; b < n - 1; b++)
        {
            var up = Run(b, +delta);
            var down = Run(b, -delta);
            slope[b] = new double[times.Length][];
            for (int k = 0; k < times.Length; k++)
            {
                slope[b][k] = new double[n];
                for (int i = 0; i < n; i++) slope[b][k][i] = (up[k][i] - down[k][i]) / (2 * delta);
            }
        }
        var f = new double[n - 1, n - 1];
        for (int a = 0; a < n - 1; a++)
            for (int b = 0; b < n - 1; b++)
            {
                double s = 0;
                for (int k = 0; k < times.Length; k++)
                    for (int i = 0; i < n; i++)
                        if (p0[k][i] > 1e-13) s += slope[a][k][i] * slope[b][k][i] / p0[k][i];
                f[a, b] = s / times.Length;
            }
        return f;
    }

    static double[,] Cosines(double[][] letters)
    {
        Check(letters);
        int m = letters.Length;
        var c = new double[m, m];
        for (int a = 0; a < m; a++)
            for (int b = 0; b < m; b++)
                c[a, b] = Dot(letters[a], letters[b], null)
                          / Math.Sqrt(Dot(letters[a], letters[a], null) * Dot(letters[b], letters[b], null));
        return c;
    }

    // The same letters read on the whole population trajectory, for any single-excitation carrier: letter[b] is
    // d p_i(t)/d delta_b over every site and every sample time (t = sample, 2*sample, .., tMax), the central
    // difference at +-delta. The carrier is a list of weighted pure waves, rho = sum_c weight_c |psi_c><psi_c|
    // (one wave for a pure carrier, two weighted site seeds for the end mixture); the cone runs each wave and adds.
    // This is the observable held fixed while the carrier changes: the bonding mode's letters are anti-collinear and
    // poorly conditioned, a localized carrier's (the site seed, the end mixture) are neither.
    public static double[][] TrajectoryLetters(World w, int n, double j, double gamma,
                                               (double Weight, double[] Wave)[] carrier,
                                               double dt, double sample, double tMax, double delta = 0.01)
    {
        int every = (int)Math.Round(sample / dt), samples = (int)Math.Round(tMax / sample);
        double[] Run(int bond, double d)
        {
            var row = new double[samples * n];
            foreach (var (weight, wave) in carrier)
            {
                var cone = new Cone(w, n, j, gamma);
                if (bond >= 0) cone.SetBond(bond, bond + 1, j * (1 + d));
                cone.SeedPure(wave);
                for (int s = 0; s < samples; s++)
                {
                    for (int k = 0; k < every; k++) cone.Step(dt);
                    for (int i = 0; i < n; i++) row[s * n + i] += weight * cone.Population(i);
                }
            }
            return row;
        }
        var letters = new double[n - 1][];
        for (int b = 0; b < n - 1; b++)
        {
            var up = Run(b, +delta);
            var down = Run(b, -delta);
            letters[b] = new double[up.Length];
            for (int k = 0; k < up.Length; k++) letters[b][k] = (up[k] - down[k]) / (2 * delta);
        }
        return letters;
    }

    // the clean chain's mode psi_k(l) = sqrt(2/(N+1)) sin(pi k (l+1)/(N+1)), k = 1..N; psi_1 is the bonding mode.
    public static double[] Mode(int n, int k)
    {
        var m = new double[n];
        for (int l = 0; l < n; l++) m[l] = Math.Sqrt(2.0 / (n + 1)) * Math.Sin(Math.PI * k * (l + 1) / (n + 1));
        return m;
    }

    // every letter finite (a site that never crossed reads NaN) and not zero (a bond the reading cannot see).
    static void Check(double[][] letters)
    {
        foreach (var l in letters)
        {
            if (l.Any(x => !double.IsFinite(x))) throw new ArgumentException("a letter is not finite (a site never crossed)", nameof(letters));
            if (l.All(x => x == 0.0)) throw new ArgumentException("a letter is zero (a bond the reading cannot see)", nameof(letters));
        }
    }

    static double Dot(double[] a, double[] b, double[]? mask)
    {
        double s = 0;
        for (int i = 0; i < a.Length; i++) if (mask is null || !double.IsNaN(mask[i])) s += a[i] * b[i];
        return s;
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
