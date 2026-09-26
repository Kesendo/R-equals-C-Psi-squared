using System.Numerics;

namespace MirrorWorld;

// The renewal cut (F126, adopted 2026-07-13 from docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md): the
// walk in the light is the clean wave, repeatedly caught and released. The single-excitation populations
// in the light obey exactly
//
//     P_n(t) = e^{-Gamma t} S_n(t),   S_n(t) = |G_{n,seed}(t)|^2 + Gamma int_0^t ds sum_m |G_{nm}(t-s)|^2 S_m(s),
//
// with G the CLEAN propagator and Gamma = 4 gamma (every refill order carries the same universal decay;
// the never-caught term is the coherent front, everything once caught is re-born and runs again). So the
// world in the light is computable from purely clean propagation plus bookkeeping: this object never steps
// the dissipator, it only accounts for it. A knower's cut in the exact sense of the README: what F126
// proves, the engine does not simulate. The pin against Cone (which DOES step the dissipator) is the
// faithfulness guard, the same pattern as Cone against Restless.
public sealed class Renewal : GameObject
{
    public int N { get; }
    public double J { get; }
    public double Gamma { get; }        // the site rate gamma; the sector coherence rate is 4 gamma
    readonly int seed;
    readonly double dt;

    // This implementation builds the chain's clean propagator. F126's renewal equation is
    // topology-blind; its Graf closed form uses the infinite chain.
    public Renewal(World world, int n, double j, double gamma, int seed, double dt) : base(world)
    {
        if (!double.IsFinite(dt) || dt <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(dt), dt, "The maximum time step must be finite and positive.");
        N = n; J = j; Gamma = gamma;
        this.seed = seed;
        this.dt = dt;
    }

    // P_n(tMax): the populations in the light, from clean propagation + the refill ladder. No dissipator step.
    public double[] Populations(double tMax)
    {
        if (!double.IsFinite(tMax) || tMax < 0.0)
            throw new ArgumentOutOfRangeException(nameof(tMax), tMax, "The requested time must be finite and nonnegative.");
        if (tMax == 0.0)
        {
            var initial = new double[N];
            initial[seed] = 1.0;
            return initial;
        }

        double stepCount = Math.Ceiling(tMax / dt);
        if (!double.IsFinite(stepCount) || stepCount >= int.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(tMax), tMax, "The requested time needs too many steps.");
        int steps = Math.Max(1, (int)stepCount);
        double h = tMax / steps; // dt is an upper bound; this uniform mesh ends at the requested time.
        // The open-chain hopping matrix has spectral radius 2|J|cos(pi/(N+1)). On an
        // imaginary eigenvalue i*x, RK4 has |R(i*x)|^2 = 1 + x^6(x^2-8)/576, so a clean
        // step outside |x| <= 2sqrt(2) amplifies norm. This is stability, not accuracy.
        double maxCleanPhase = N <= 1 ? 0.0 : 2.0 * Math.Abs(J) * Math.Cos(Math.PI / (N + 1.0)) * h;
        if (!double.IsFinite(maxCleanPhase) || maxCleanPhase > 2.0 * Math.Sqrt(2.0))
            throw new ArgumentOutOfRangeException(nameof(dt), dt, "The clean RK4 step is outside its imaginary-axis stability interval; reduce the maximum time step.");
        // Stability controls norm, not phase. RK4 is the fourth Taylor polynomial of exp(-iHh):
        // its one-step operator error is at most x^5/120 for x = h*||H||. Both the exact step and
        // the stable RK4 step have norm <= 1, so telescoping bounds the error after any k <= steps
        // by steps*x^5/120. A measurement distribution differs in L1 by at most twice that norm.
        const double maxCleanPopulationError = 1e-3;
        double cleanPopulationErrorBound = 2.0 * steps * Math.Pow(maxCleanPhase, 5) / 120.0;
        if (cleanPopulationErrorBound > maxCleanPopulationError)
            throw new ArgumentOutOfRangeException(nameof(dt), dt,
                "The clean RK4 population-error bound exceeds 0.1%; reduce the maximum time step.");
        double gPhi = 4.0 * Gamma;
        double dose = gPhi * h;
        const double maxMassDrift = 1e-3;
        // At dose = 0.5, even the exactly stationary J=0 control drifts by 1.1% in one
        // trapezoid step: exp(-dose) * (1 + dose/2) / (1 - dose/2). The pole is at dose = 2.
        if (!double.IsFinite(dose) || dose >= 0.5)
            throw new ArgumentOutOfRangeException(nameof(dt), dt, "The renewal dose 4*gamma*h must be below 0.5; reduce the maximum time step.");
        // For J=0, the excess log mass per step is 2*atanh(dose/2)-dose. Its positive
        // series is bounded by dose^3/[12*(1-dose^2/4)] at the allowed positive doses.
        // Bound the accumulated inflation before building the O(steps^2) convolution.
        if (dose > 0.0)
        {
            double logInflationBound = steps * dose * dose * dose / (12.0 * (1.0 - dose * dose / 4.0));
            if (logInflationBound > Math.Log(1.0 + maxMassDrift))
                throw new ArgumentOutOfRangeException(nameof(dt), dt, "The accumulated refill-grid error bound exceeds the 0.1% budget; reduce the maximum time step.");
        }

        double[] coarse = ComputePopulations(steps, tMax, gPhi);
        CheckMass(coarse, maxMassDrift);
        if (gPhi == 0.0 || J == 0.0 || N <= 1) return coarse;

        // Total mass sees only the column sums of the clean kernel. It cannot detect a coarse
        // refill grid that puts the right mass on the wrong sites. Compare with a doubled grid;
        // this is a convergence check, not an absolute error certificate for the continuum law.
        if (steps > int.MaxValue / 2)
            throw new ArgumentOutOfRangeException(nameof(dt), dt, "The refill grid cannot be doubled for a spatial convergence check.");
        double[] fine = ComputePopulations(2 * steps, tMax, gPhi);
        CheckMass(fine, maxMassDrift);
        double discrepancy = 0.0;
        for (int n = 0; n < N; n++) discrepancy = Math.Max(discrepancy, Math.Abs(coarse[n] - fine[n]));
        if (discrepancy > 1e-3)
            throw new ArgumentOutOfRangeException(nameof(dt), dt,
                "The refill populations change by more than 0.1% when the grid is halved; reduce the maximum time step.");
        return fine;
    }

    void CheckMass(double[] populations, double maxMassDrift)
    {
        double totalMass = populations.Sum();
        if (!double.IsFinite(totalMass) || Math.Abs(totalMass - 1.0) > maxMassDrift)
            throw new ArgumentOutOfRangeException(nameof(dt), dt, "The computed populations violate the 0.1% mass-conservation budget; reduce the maximum time step.");
    }

    double[] ComputePopulations(int steps, double tMax, double gPhi)
    {
        double h = tMax / steps;
        double dose = gPhi * h;
        // the clean kernel K[k][m,n] = |<n| e^{-i H (k h)} |m>|^2: evolve U-dot = -i H U from U(0) = I.
        var u = new Complex[N, N];
        for (int m = 0; m < N; m++) u[m, m] = Complex.One;
        var kernel = new double[steps + 1][,];
        kernel[0] = Squares(u);
        for (int k = 1; k <= steps; k++)
        {
            var k1 = Rhs(u);
            var k2 = Rhs(Axpy(u, k1, h / 2));
            var k3 = Rhs(Axpy(u, k2, h / 2));
            var k4 = Rhs(Axpy(u, k3, h));
            for (int a = 0; a < N; a++)
                for (int b = 0; b < N; b++)
                    u[a, b] += (h / 6) * (k1[a, b] + 2 * k2[a, b] + 2 * k3[a, b] + k4[a, b]);
            kernel[k] = Squares(u);
        }

        // the Volterra refill ladder: trapezoid in s, the s = k self-term (kernel[0] = identity) implicit.
        var S = new double[steps + 1][];
        S[0] = new double[N]; S[0][seed] = 1.0;
        double denom = 1.0 - 0.5 * dose;
        for (int k = 1; k <= steps; k++)
        {
            var row = new double[N];
            for (int n = 0; n < N; n++) row[n] = kernel[k][seed, n];          // the never-caught term
            for (int s = 0; s < k; s++)
            {
                double w = (s == 0) ? 0.5 : 1.0;
                var ks = kernel[k - s];
                var Ss = S[s];
                for (int m = 0; m < N; m++)
                {
                    if (Ss[m] == 0.0) continue;
                    double c = w * dose * Ss[m];
                    for (int n = 0; n < N; n++) row[n] += c * ks[m, n];
                }
            }
            for (int n = 0; n < N; n++) row[n] /= denom;
            S[k] = row;
        }

        var p = new double[N];
        double damp = Math.Exp(-gPhi * tMax);
        for (int n = 0; n < N; n++)
        {
            p[n] = damp * S[steps][n];
        }
        return p;
    }

    Complex[,] Rhs(Complex[,] x)
    {
        int n = x.GetLength(0);
        var r = new Complex[n, n];
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
            {
                Complex hx = Complex.Zero;
                if (a > 0) hx += J * x[a - 1, b];
                if (a < n - 1) hx += J * x[a + 1, b];
                r[a, b] = -Complex.ImaginaryOne * hx;
            }
        return r;
    }

    static Complex[,] Axpy(Complex[,] u, Complex[,] v, double s)
    {
        int n = u.GetLength(0);
        var r = new Complex[n, n];
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++)
                r[a, b] = u[a, b] + s * v[a, b];
        return r;
    }

    static double[,] Squares(Complex[,] u)
    {
        int n = u.GetLength(0);
        var k = new double[n, n];
        for (int m = 0; m < n; m++)
            for (int a = 0; a < n; a++)
                k[m, a] = u[a, m].Real * u[a, m].Real + u[a, m].Imaginary * u[a, m].Imaginary;
        return k;
    }

    public override IReadOnlyList<string> Own => new[] { "populations" };
}
