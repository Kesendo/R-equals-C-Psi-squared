using System.Numerics;

namespace MirrorWorld;

// The renewal cut (F126, adopted 2026-07-13 from docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md): the
// watched walk is the unwatched wave, repeatedly caught and released. The single-excitation populations
// under the watching obey exactly
//
//     P_n(t) = e^{-Gamma t} S_n(t),   S_n(t) = |G_{n,seed}(t)|^2 + Gamma int_0^t ds sum_m |G_{nm}(t-s)|^2 S_m(s),
//
// with G the CLEAN propagator and Gamma = 4 gamma (every refill order carries the same universal decay;
// the never-caught term is the coherent front, everything once caught is re-born and runs again). So the
// watched world is computable from purely unwatched propagation plus bookkeeping: this object never steps
// the dissipator, it only accounts for it. A knower's cut in the exact sense of the README: what F126
// proves, the engine does not simulate. The pin against Cone (which DOES step the watching) is the
// faithfulness guard, the same pattern as Cone against Restless.
public sealed class Renewal : GameObject
{
    public int N { get; }
    public double J { get; }
    public double Gamma { get; }        // the site rate gamma; the sector coherence rate is 4 gamma
    readonly int seed;
    readonly double dt;

    // chain topology only: F126 is proven on the chain (the clean propagator G and its Graf closure).
    public Renewal(World world, int n, double j, double gamma, int seed, double dt) : base(world)
    {
        N = n; J = j; Gamma = gamma;
        this.seed = seed;
        this.dt = dt;
    }

    // P_n(tMax): the watched populations, from clean propagation + the refill ladder. No dissipator step.
    public double[] Populations(double tMax)
    {
        int steps = (int)Math.Round(tMax / dt);
        double gPhi = 4.0 * Gamma;

        // the clean kernel K[k][m,n] = |<n| e^{-ih (k dt)} |m>|^2: evolve U-dot = -i h U from U(0) = I.
        var u = new Complex[N, N];
        for (int m = 0; m < N; m++) u[m, m] = Complex.One;
        var kernel = new double[steps + 1][,];
        kernel[0] = Squares(u);
        for (int k = 1; k <= steps; k++)
        {
            var k1 = Rhs(u);
            var k2 = Rhs(Axpy(u, k1, dt / 2));
            var k3 = Rhs(Axpy(u, k2, dt / 2));
            var k4 = Rhs(Axpy(u, k3, dt));
            for (int a = 0; a < N; a++)
                for (int b = 0; b < N; b++)
                    u[a, b] += (dt / 6) * (k1[a, b] + 2 * k2[a, b] + 2 * k3[a, b] + k4[a, b]);
            kernel[k] = Squares(u);
        }

        // the Volterra refill ladder: trapezoid in s, the s = k self-term (kernel[0] = identity) implicit.
        var S = new double[steps + 1][];
        S[0] = new double[N]; S[0][seed] = 1.0;
        double denom = 1.0 - 0.5 * gPhi * dt;
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
                    double c = w * gPhi * dt * Ss[m];
                    for (int n = 0; n < N; n++) row[n] += c * ks[m, n];
                }
            }
            for (int n = 0; n < N; n++) row[n] /= denom;
            S[k] = row;
        }

        var p = new double[N];
        double damp = Math.Exp(-gPhi * steps * dt);
        for (int n = 0; n < N; n++) p[n] = damp * S[steps][n];
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
