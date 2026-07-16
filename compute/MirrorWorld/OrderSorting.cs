namespace MirrorWorld;

// The mirror's order-sorting law, Theorem A (F131, adopted 2026-07-16 from
// docs/proofs/PROOF_MIRROR_ORDER_SORTING.md; minted in docs/ANALYTICAL_FORMULAS.md): a mirror does
// not throw information away, it sorts it. Let R be the F71 site reversal and split the chain's
// parameters into an R-even base plus t times an R-odd direction. The conjugation identity
//
//     (R x R) . L(base + t*dir) . (R x R) = L(base - t*dir)          (exact, entry-wise)
//
// is the pencil face (ParameterKlein.MirrorConjugationResidual, machine zero on all three axes);
// this object runs the TRAJECTORY face: for an OPERATOR-R-even preparation (R rho0 R = rho0) and a
// readout O of definite R-parity q (R O R = q O), every evolution time obeys
//
//     <O>(t; time) = q * <O>(-t; time),
//
// so the response orders sort by q * sigma_eff into four cells: (+,-) EVEN response, (-,-) ODD
// response, (-,+) IDENTICALLY ZERO (a pure selection rule, not the Pi-protected cluster
// cancellation), (+,+) generic. Here sigma_eff = sigma_op * chi_M = -1 for an R-odd direction
// under the linear mirror R (chi_R = +1); the antiunitary column chi_M = -1 (Theorem B, the zeta^2
// anti-protection law on the Floquet step) stays in the main repo with its tracking hypotheses,
// and the spectral-evenness corollary stays out too (an eigensolver has no place here). The
// hypotheses are the physics: a preparation that is merely expectation-even leaks, and the leak is
// EXACTLY affine in the odd admixture eps (the master equation is linear in rho0) -- halving eps
// halves the forbidden channel, ratio 2.000.
//
// Twin RK4 runs on the gamma axis (Restless.siteGammas; the site-resolved watching is all the law
// needs -- every cell is reachable without per-bond J). No eigensolver anywhere.
public sealed class OrderSorting : GameObject
{
    public int N { get; }
    public double J { get; }

    public OrderSorting(World world, int n, double j = 1.0) : base(world)
    {
        N = n;
        J = j;
    }

    // left: what the order-sorting law itself produces.
    public override IReadOnlyList<string> Own => new[] { "conjugation", "cells", "leak" };

    public sealed record SortingReport(
        double ConjugationResidual,   // T0, gamma axis: the pencil identity under the twin scan
        double EvenCellResidual,      // (+,-): worst |<O_even>(+t) - <O_even>(-t)| over probe times
        double OddCellResidual,       // (-,-): worst |<O_odd>(+t) + <O_odd>(-t)| over probe times
        double OddMagnitude,          // (-,-) non-vacuity: max |<O_odd>| (the allowed orders are O(1))
        double ZeroCellWorst,         // (-,+): worst |<O_odd>| under the R-even generator, EVERY tick
        double GenericGap,            // (+,+): |<O_even>(+t) - <O_even>(-t)| under an R-even direction, O(1)
        double LeakSlope,             // the forbidden odd-in-t channel opened by eps * rho_odd
        double LeakRatio);            // slope(eps) / slope(eps/2) = 2.000 exactly (affine)

    // ---- the profiles: an R-even base (bitwise symmetric by construction) and the two scan
    // directions -- R-odd (dir_{N-1-l} = -dir_l bitwise, sigma_op = -1) and R-even (sigma_op = +1).
    public static double[] SymmetricBase(int n)
    {
        var r = new double[n];
        for (int l = 0; 2 * l <= n - 1; l++)
        {
            double v = 0.05 + 0.03 * l;
            r[l] = v;
            r[n - 1 - l] = v;
        }
        return r;
    }

    public static double[] OddDirection(int n)
    {
        var r = new double[n];
        for (int l = 0; 2 * l < n - 1; l++)
        {
            double d = 0.1 * (2.0 * l - (n - 1)) / (n - 1);
            r[l] = d;
            r[n - 1 - l] = -d;
        }
        return r;                                            // middle site (odd n) stays 0
    }

    public static double[] EvenDirection(int n)
    {
        var r = new double[n];
        r[0] = 0.1;
        r[n - 1] = 0.1;
        return r;
    }

    // the F71 site reversal on a basis configuration.
    public static int Reverse(int c, int n)
    {
        int r = 0;
        for (int s = 0; s < n; s++)
            if (((c >> s) & 1) == 1) r |= 1 << (n - 1 - s);
        return r;
    }

    // ---- the preparations (single-excitation, block (1,1)): psi = (|e_0> + |e_{N-1}>)/sqrt(2),
    // rho_even = |psi><psi| with R rho R = rho as an OPERATOR identity; the leak adds
    // eps * rho_odd with rho_odd = (|e_0><e_0| - |e_{N-1}><e_{N-1}|)/2, operator-R-ODD and purely
    // diagonal (psi and its R-odd partner share the diagonal; the cross terms cancel). ----
    public static void SeedEvenPrep(Restless r, int n)
    {
        int u = 1, v = 1 << (n - 1);
        r.Seed(u, 0.5);
        r.Seed(v, 0.5);
        r.SeedCoherence(u, v, 0.5);
    }

    public static void SeedLeakPrep(Restless r, int n, double eps)
    {
        int u = 1, v = 1 << (n - 1);
        r.Seed(u, 0.5 + eps / 2);
        r.Seed(v, 0.5 - eps / 2);
        r.SeedCoherence(u, v, 0.5);
    }

    // ---- the readouts: <Z_s> off the diagonal of rho; O_even = Z_0 + Z_{N-1} (q = +1),
    // O_odd = Z_0 - Z_{N-1} (q = -1). Populations only -- the law needs no generic expectation. ----
    static double ExpectZ(Restless r, int site)
    {
        double e = 0;
        for (int i = 0; i < r.Dim; i++) e += (1 - 2 * ((i >> site) & 1)) * r[i, i].Real;
        return e;
    }

    static double EvenReadout(Restless r, int n) => ExpectZ(r, 0) + ExpectZ(r, n - 1);
    static double OddReadout(Restless r, int n) => ExpectZ(r, 0) - ExpectZ(r, n - 1);

    Restless NewRun(double[] gammas, Action<Restless, int> seed)
    {
        var run = new Restless((World)Parent!, N, J, gammas.Average(), siteGammas: gammas);
        seed(run, N);
        return run;
    }

    double[] Scan(double[] baseP, double[] dir, double t)
    {
        var r = new double[N];
        for (int l = 0; l < N; l++) r[l] = baseP[l] + t * dir[l];
        return r;
    }

    // ---- the four-cell run: twin RK4 at +t and -t (R-odd direction) reads the two response
    // cells; the R-even generator run reads the zero cell at EVERY tick; the R-even-direction twin
    // reads the generic cell from below; the leak twins read the affine hypothesis violation. ----
    public SortingReport Run(double t = 0.3, double eps = 0.02, double dt = 0.02,
                             int ticksPerProbe = 35, int probes = 3)
    {
        var gBase = SymmetricBase(N);
        var gOdd = OddDirection(N);
        var gEven = EvenDirection(N);

        var klein = new ParameterKlein((World)Parent!, N);
        var jU = Enumerable.Repeat(1.0, N - 1).ToArray();
        var hZ = new double[N];
        double conj = klein.MirrorConjugationResidual(
            Scan(gBase, gOdd, +t), jU, hZ, Scan(gBase, gOdd, -t), jU, hZ);

        var plus = NewRun(Scan(gBase, gOdd, +t), SeedEvenPrep);
        var minus = NewRun(Scan(gBase, gOdd, -t), SeedEvenPrep);
        var zero = NewRun(gBase, SeedEvenPrep);
        var genP = NewRun(Scan(gBase, gEven, +t), SeedEvenPrep);
        var genM = NewRun(Scan(gBase, gEven, -t), SeedEvenPrep);

        double evenRes = 0, oddRes = 0, oddMag = 0, genericGap = 0;
        double zeroWorst = Math.Abs(OddReadout(zero, N));
        for (int probe = 0; probe < probes; probe++)
        {
            for (int k = 0; k < ticksPerProbe; k++)
            {
                plus.Step(dt);
                minus.Step(dt);
                genP.Step(dt);
                genM.Step(dt);
                zero.Step(dt);
                zeroWorst = Math.Max(zeroWorst, Math.Abs(OddReadout(zero, N)));
            }
            evenRes = Math.Max(evenRes, Math.Abs(EvenReadout(plus, N) - EvenReadout(minus, N)));
            oddRes = Math.Max(oddRes, Math.Abs(OddReadout(plus, N) + OddReadout(minus, N)));
            oddMag = Math.Max(oddMag, Math.Abs(OddReadout(plus, N)));
            genericGap = Math.Max(genericGap, Math.Abs(EvenReadout(genP, N) - EvenReadout(genM, N)));
        }

        double s1 = LeakSlope(gBase, gOdd, eps, dt, ticksPerProbe);
        double s2 = LeakSlope(gBase, gOdd, eps / 2, dt, ticksPerProbe);
        return new SortingReport(conj, evenRes, oddRes, oddMag, zeroWorst, genericGap, s1, s1 / s2);
    }

    // the forbidden channel under the broken hypothesis: the odd-in-t component of the EVEN
    // readout, [<O_even>(+t) - <O_even>(-t)] / (2t), at a small scan t and one probe time.
    double LeakSlope(double[] gBase, double[] gOdd, double eps, double dt, int ticks)
    {
        const double tScan = 0.05;
        var lp = NewRun(Scan(gBase, gOdd, +tScan), (r, n) => SeedLeakPrep(r, n, eps));
        var lm = NewRun(Scan(gBase, gOdd, -tScan), (r, n) => SeedLeakPrep(r, n, eps));
        for (int k = 0; k < ticks; k++)
        {
            lp.Step(dt);
            lm.Step(dt);
        }
        return (EvenReadout(lp, N) - EvenReadout(lm, N)) / (2 * tScan);
    }
}
