namespace MirrorWorld;

// The flavor split of experiments/FLAVOR_RESOLVED_T2_INHERITANCE.md's neural reading, read as F36's pairing in
// closed form (adopted 2026-10-08 from that page's neural section and simulations/neural/flavor_split_is_f36.py).
//
// The generator is J = [[-gE I + c A, -h I], [h I, -gI I - c A]] on an adjacency A scaled to spectral radius one
// (the scaling is the producer's; here A arrives already scaled, and its eigenvalues a arrive as inputs, so no
// eigensolver runs). With Q the per-node E <-> I swap and s = (gE + gI)/2 it meets F36 identically, Q J Q = -J - 2s I:
// the leak condition pairs every E seat with its own I seat, the off-diagonal condition holds because the coupling
// enters as +cA on E and -cA on I and the h block is antisymmetric.
//
// Per eigenvalue a of A with eigenvector w, the pair (w, 0), (0, w) spans the invariant block
//     -s I + [[x, -h], [h, -x]],   x = Delta + c a,   Delta = (gI - gE)/2,
// with eigenvalues -s +- sqrt(x^2 - h^2). The block is OVERDAMPED (two real rates mirrored about s) iff |x| > h, its
// slow mode excitation-dominant iff x > h and inhibition-dominant iff x < -h; an underdamped block (|x| < h) sits on
// Re = -s with E weight exactly one half; |x| = h is a Jordan block at -s. So the split of lifetimes by flavor is
// the mirror pairing read on the overdamped blocks, and the page's reported ratio is the slowest E rate (the block
// a = 1) over the faster rate s + sqrt(.) of the weakest overdamped block. Which flavor is slow is the sign of x,
// not F36. The page's classifier (neural_flavor_rule.classify_mode) calls a mode E-dominant at I weight below 0.35
// and I-dominant above 0.65 and reports the slowest stable mode of each; ReportedRatio implements that rule on the
// closed form's modes, so it agrees with the producer wherever the producer's eigensolver does. The coupling c is
// taken non-negative (the producer's sign). Own: the generator, the swap, the block, the discriminant, the rates,
// the mode weights, the flavor, the reported ratio and the h = 0 form (gI + c a_min)/(gE - c a_max). Parent: none
// (a static companion of NeuralPalindrome, which owns the identity itself). Brains, baths, temperature and the
// Absorption split are outside: this generator has no density matrix and no dissipator, only leaks.
public static class FlavorSplit
{
    private static void RequireLeak(double gamma, string name)
    {
        if (!double.IsFinite(gamma) || gamma <= 0)
            throw new ArgumentOutOfRangeException(name, "A leak rate must be finite and positive.");
    }
    private static void RequireFinite(double value, string name)
    {
        if (!double.IsFinite(value)) throw new ArgumentOutOfRangeException(name, "Must be finite.");
    }
    private static void RequireCoupling(double coupling)
    {
        if (!double.IsFinite(coupling) || coupling < 0)
            throw new ArgumentOutOfRangeException(nameof(coupling), "The graph coupling must be finite and non-negative (the producer's sign; a_max then carries the slowest excitatory rate).");
    }

    // s = (gE + gI)/2, the mirror centre's offset (F36's s for this generator).
    public static double Centre(double gammaE, double gammaI)
    {
        RequireLeak(gammaE, nameof(gammaE)); RequireLeak(gammaI, nameof(gammaI));
        return 0.5 * (gammaE + gammaI);
    }

    // Delta = (gI - gE)/2, the half leak difference.
    public static double HalfLeakDifference(double gammaE, double gammaI)
    {
        RequireLeak(gammaE, nameof(gammaE)); RequireLeak(gammaI, nameof(gammaI));
        return 0.5 * (gammaI - gammaE);
    }

    // The 2N x 2N generator on a scaled adjacency (rows 0..N-1 excitatory, N..2N-1 inhibitory).
    public static double[,] Jacobian(double[,] scaledAdjacency, double h, double coupling, double gammaE, double gammaI)
    {
        ArgumentNullException.ThrowIfNull(scaledAdjacency);
        int n = scaledAdjacency.GetLength(0);
        if (n == 0 || scaledAdjacency.GetLength(1) != n)
            throw new ArgumentException("The adjacency must be square and non-empty.", nameof(scaledAdjacency));
        RequireFinite(h, nameof(h)); RequireCoupling(coupling);
        RequireLeak(gammaE, nameof(gammaE)); RequireLeak(gammaI, nameof(gammaI));
        var j = new double[2 * n, 2 * n];
        for (int r = 0; r < n; r++)
            for (int k = 0; k < n; k++)
            {
                double a = scaledAdjacency[r, k];
                if (!double.IsFinite(a)) throw new ArgumentException("Adjacency entries must be finite.", nameof(scaledAdjacency));
                j[r, k] = coupling * a;
                j[n + r, n + k] = -coupling * a;
            }
        for (int r = 0; r < n; r++)
        {
            j[r, r] -= gammaE;
            j[n + r, n + r] -= gammaI;
            j[r, n + r] = -h;
            j[n + r, r] = h;
        }
        return j;
    }

    // Q: seat r <-> seat n + r.
    public static int[] SwapPermutation(int n)
    {
        if (n <= 0) throw new ArgumentOutOfRangeException(nameof(n));
        var q = new int[2 * n];
        for (int r = 0; r < n; r++) { q[r] = n + r; q[n + r] = r; }
        return q;
    }

    // Max |Q J Q + J + 2 s I| with Q the swap and s the centre: F36 on this generator, summed in the order that makes
    // it exact for EVERY input: off the diagonal every pair is (c a) + (-c a) or (-h) + h, exactly 0; on the diagonal
    // the two leaks are summed first, (-gE) + (-gI), and 2 s = gE + gI is the same rounded sum with the opposite sign
    // (halving and doubling are exact), so the residual is exactly 0.0 at generic leaks. NeuralPalindrome.MaxResidual
    // centres each diagonal entry before adding and carries one rounding at generic leaks (a summation-order
    // residual, CLAUDE.md's case 3); it is read beside this route in the tests, not gated.
    public static double F36Residual(double[,] scaledAdjacency, double h, double coupling, double gammaE, double gammaI)
    {
        var j = Jacobian(scaledAdjacency, h, coupling, gammaE, gammaI);
        int n = scaledAdjacency.GetLength(0); int[] q = SwapPermutation(n); double twoS = gammaE + gammaI;
        double worst = 0;
        for (int i = 0; i < 2 * n; i++)
            for (int k = 0; k < 2 * n; k++)
            {
                double r = i == k ? (j[i, i] + j[q[i], q[i]]) + twoS : j[i, k] + j[q[i], q[k]];
                worst = Math.Max(worst, Math.Abs(r));
            }
        return worst;
    }

    // The host's route to the same residual (diagonals centred first): read, not gated, at generic leaks.
    public static double HostRouteResidual(double[,] scaledAdjacency, double h, double coupling, double gammaE, double gammaI)
    {
        var j = Jacobian(scaledAdjacency, h, coupling, gammaE, gammaI);
        return NeuralPalindrome.MaxResidual(j, SwapPermutation(scaledAdjacency.GetLength(0)), Centre(gammaE, gammaI));
    }

    // x = Delta + c a, the block's off-centre diagonal.
    public static double BlockDiagonal(double a, double coupling, double gammaE, double gammaI)
    {
        RequireFinite(a, nameof(a)); RequireCoupling(coupling);
        return HalfLeakDifference(gammaE, gammaI) + coupling * a;
    }

    // x^2 - h^2: positive overdamped, negative underdamped, zero a Jordan block. The zero is decided by float
    // equality, exact where x is computed exactly (dyadic inputs, the ring's a = 0 pin); an exceptional point
    // reached through generic inputs rounds off zero and reads as one of the two sides.
    public static double Discriminant(double a, double h, double coupling, double gammaE, double gammaI)
    {
        RequireFinite(h, nameof(h));
        double x = BlockDiagonal(a, coupling, gammaE, gammaI);
        return x * x - h * h;
    }

    public static bool IsOverdamped(double a, double h, double coupling, double gammaE, double gammaI)
        => Discriminant(a, h, coupling, gammaE, gammaI) > 0;

    public static bool IsDefective(double a, double h, double coupling, double gammaE, double gammaI)
        => Discriminant(a, h, coupling, gammaE, gammaI) == 0;

    // The two real rates s -+ sqrt(x^2 - h^2) of an overdamped block; throws otherwise.
    public static (double Slow, double Fast) OverdampedRates(double a, double h, double coupling, double gammaE, double gammaI)
    {
        double disc = Discriminant(a, h, coupling, gammaE, gammaI);
        if (disc <= 0) throw new InvalidOperationException("The block is not overdamped; its rates are not two reals.");
        double s = Centre(gammaE, gammaI); double root = Math.Sqrt(disc);
        return (s - root, s + root);
    }

    // Which flavor the block's SLOW mode carries in the weight-above-one-half sense: "E" if x > h, "I" if x < -h,
    // "mixed" when |x| <= h (weight exactly one half in the underdamped case, a Jordan block at |x| = h). The page's
    // classifier is stricter (0.35 / 0.65); ModeInhibitionWeight gives the number it thresholds.
    public static string SlowFlavor(double a, double h, double coupling, double gammaE, double gammaI)
    {
        RequireFinite(h, nameof(h));
        double x = BlockDiagonal(a, coupling, gammaE, gammaI);
        if (x > Math.Abs(h)) return "E";
        if (x < -Math.Abs(h)) return "I";
        return "mixed";
    }

    // The inhibition weight of the block mode with relative eigenvalue mu (lambda = -s + mu, rate s - mu): the
    // eigenvector solves (x - mu) alpha = h beta, so the weight is (x - mu)^2 / (h^2 + (x - mu)^2); at h = 0 the
    // modes are the pure seats, weight 0 at mu = x and 1 at mu = -x.
    public static double ModeInhibitionWeight(double x, double h, double mu)
    {
        RequireFinite(x, nameof(x)); RequireFinite(h, nameof(h)); RequireFinite(mu, nameof(mu));
        if (h == 0) return mu == x ? 0 : 1;
        double b = x - mu; return b * b / (h * h + b * b);
    }

    // Max |J (alpha w, beta w) - ((x alpha - h beta) w, (h alpha - x beta) w) + s (alpha w, beta w)|: the block action on
    // an adjacency eigenvector w with eigenvalue a. Exactly 0 on dyadic inputs; the invariance of the block.
    public static double BlockActionResidual(double[,] scaledAdjacency, double[] w, double a, double alpha, double beta,
                                             double h, double coupling, double gammaE, double gammaI)
    {
        ArgumentNullException.ThrowIfNull(w);
        int n = scaledAdjacency.GetLength(0);
        if (w.Length != n) throw new ArgumentException("The eigenvector must match the adjacency.", nameof(w));
        var j = Jacobian(scaledAdjacency, h, coupling, gammaE, gammaI);
        double s = Centre(gammaE, gammaI); double x = BlockDiagonal(a, coupling, gammaE, gammaI);
        var v = new double[2 * n];
        for (int r = 0; r < n; r++) { v[r] = alpha * w[r]; v[n + r] = beta * w[r]; }
        double worst = 0;
        for (int r = 0; r < 2 * n; r++)
        {
            double jv = 0;
            for (int k = 0; k < 2 * n; k++) jv += j[r, k] * v[k];
            double expected = r < n ? (x * alpha - h * beta) * w[r] - s * v[r] : (h * alpha - x * beta) * w[r - n] - s * v[r];
            worst = Math.Max(worst, Math.Abs(jv - expected));
        }
        return worst;
    }

    // The ratio the page reports, by the producer's rule on the closed form's modes: over every block's two modes
    // (relative eigenvalues +-sqrt(disc) when overdamped; an underdamped or Jordan block has weight one half and is
    // "mixed"), keep the stable ones (rate s - mu > 0), call a mode E-dominant at inhibition weight below 0.35 and
    // I-dominant above 0.65, and return the slowest I-dominant rate over the slowest E-dominant rate. Throws when
    // either flavor is missing, the page's "classifier-failed". On the page's own rows every overdamped block has
    // x > h, so this is the slowest E rate (block a = 1) over the faster rate of the weakest overdamped block; off
    // those rows (a block with x < -h, or a weakest block whose fast mode carries less than 0.65) the rule and that
    // shortcut part, and the rule is what the page's numbers are.
    public static double ReportedRatio(double[] adjacencyEigenvalues, double h, double coupling, double gammaE, double gammaI)
    {
        ArgumentNullException.ThrowIfNull(adjacencyEigenvalues);
        RequireFinite(h, nameof(h));
        double s = Centre(gammaE, gammaI);
        double slowE = double.NaN, slowI = double.NaN;
        foreach (double a in adjacencyEigenvalues)
        {
            double x = BlockDiagonal(a, coupling, gammaE, gammaI); double disc = x * x - h * h;
            if (disc <= 0) continue;                            // underdamped or Jordan: weight one half, mixed
            double root = Math.Sqrt(disc);
            foreach (double mu in new[] { root, -root })
            {
                double rate = s - mu;
                if (rate <= 0) continue;                        // unstable, skipped by the producer too
                double wI = ModeInhibitionWeight(x, h, mu);
                if (wI < 0.35 && (double.IsNaN(slowE) || rate < slowE)) slowE = rate;
                if (wI > 0.65 && (double.IsNaN(slowI) || rate < slowI)) slowI = rate;
            }
        }
        if (double.IsNaN(slowE) || double.IsNaN(slowI))
            throw new InvalidOperationException("The classifier finds no stable mode of one flavor: the page's 'classifier-failed'.");
        return slowI / slowE;
    }

    // At h = 0 every block with x != 0 is overdamped at rates gE - c a (E) and gI + c a (I), so the reported ratio is
    // (gI + c a_min)/(gE - c a_max).
    public static double ZeroCrossCouplingRatio(double aMin, double aMax, double coupling, double gammaE, double gammaI)
    {
        RequireFinite(aMin, nameof(aMin)); RequireFinite(aMax, nameof(aMax)); RequireCoupling(coupling);
        RequireLeak(gammaE, nameof(gammaE)); RequireLeak(gammaI, nameof(gammaI));
        double slowE = gammaE - coupling * aMax; double slowI = gammaI + coupling * aMin;
        if (slowE <= 0 || slowI <= 0) throw new InvalidOperationException("A slowest rate is not positive: the network is not stable there.");
        return slowI / slowE;
    }
}
