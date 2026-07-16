using System.Numerics;

namespace MirrorWorld;

// The bridged lattice of worlds (2026-07-16, the engine beat the fourth play offered): the Klein
// V4 of watchings, run dynamically. The anti-turn (2026-07-03) showed one second world; the
// double turn showed the way home. Here all four vertices run side by side as their own worlds:
//
//     e  = rho(t)                the normal world      rule -2*gamma*k       (disagreement watched)
//     L  = X^N * rho(t)          the ket-side reading  rule -2*gamma*(N-k)   (agreement watched)
//     R  = rho(t) * X^N          the bra-side reading  rule -2*gamma*(N-k)   (the 2026-07-03 anti-world)
//     LR = X^N * rho(t) * X^N    the double turn       rule -2*gamma*k       (home again)
//
// Every edge is an exact always-open identity, because [H, X^N] = 0 (the handshake s+s- + s-s+
// is invariant when every spin flips; the ZZ bond too, since X flips the sign of every Z and the
// signs cancel pairwise) and k(~i, j) = k(i, ~j) = N - k(i, j), k(~i, ~j) = k(i, j). Nothing is
// interpreted here: L and R are not new states (they are not Hermitian), they are the SAME state
// seen through one complemented leg, and the loop runs the reading as its own world under the
// turned rule. The bridges checked at every tick:
//
//     L(t)[i,j]  = e(t)[~i, j]        R(t)[i,j] = e(t)[i, ~j]        LR(t)[i,j] = e(t)[~i, ~j]
//     LR(t)[i,j] = L(t)[i, ~j] = R(t)[~i, j]                         (the V4 composes: L o R = LR)
//     L(t) = R(t)-dagger                                             (the dagger exchanges the legs)
//
// The conservation law is never broken, it MOVES: the carried unit sits on the trace at e and LR
// and on the anti-trace at L and R; the immortal set moves with it (the diagonal at e/LR has rate
// 0 under the normal rule, the ANTI-diagonal at L/R has rate -2*gamma*(N-N) = 0 under the turned
// one). The dagger edge is the antilinear triangle's dagger read at the lattice level (F119: the
// two one-sided readings are each other's adjoint; both pay the same turned rule).
public sealed class Lattice : GameObject
{
    public int N { get; }
    public double J { get; }
    public double Gamma { get; }
    readonly double zz;

    public Lattice(World world, int n, double j = 1.0, double gamma = 0.5, double zz = 0.0) : base(world)
    {
        N = n;
        J = j;
        Gamma = gamma;
        this.zz = zz;
    }

    // left: what the lattice itself produces.
    public override IReadOnlyList<string> Own => new[] { "vertices", "bridges", "carried-unit" };

    public sealed record LatticeReport(
        double WorstBridgeL,        // worst |L[i,j] - e[~i,j]| over every probed tick
        double WorstBridgeR,        // worst |R[i,j] - e[i,~j]|
        double WorstBridgeLR,       // worst |LR[i,j] - e[~i,~j]|
        double WorstComposition,    // worst |LR[i,j] - L[i,~j]| and |LR[i,j] - R[~i,j]| (L o R = LR = R o L)
        double WorstDagger,         // worst |L[i,j] - conj(R[j,i])| (the dagger exchanges the legs)
        double WorstUnitE,          // worst |1 - trace(e)|
        double WorstUnitLR,         // worst |1 - trace(LR)|
        double WorstUnitL,          // worst |1 - antitrace(L)|
        double WorstUnitR,          // worst |1 - antitrace(R)|
        double VertexSeparation,    // min over vertex pairs of the max entry difference at the final tick
        double[] HistogramE,        // final k-histogram of e (weight per disagreement k)
        double[] HistogramL);       // final k-histogram of L (the e one, read backward)

    // ---- run the four worlds side by side and check every edge at every tick. ----
    public LatticeReport Run(int seed, double dt, int ticks)
    {
        var w = (World)Parent!;
        int dim = 1 << N, bar = dim - 1;
        int sbar = bar - seed;

        var e = new Restless(w, N, J, Gamma, zz: zz);
        e.Seed(seed);
        var l = new Restless(w, N, J, Gamma, antiWatching: true, zz: zz);
        l.SeedRaw(sbar, seed, 1.0);                                  // X^N |s><s|  = |~s><s|
        var r = new Restless(w, N, J, Gamma, antiWatching: true, zz: zz);
        r.SeedRaw(seed, sbar, 1.0);                                  // |s><s| X^N  = |s><~s|
        var lr = new Restless(w, N, J, Gamma, zz: zz);
        lr.Seed(sbar);                                               // X^N |s><s| X^N = |~s><~s|

        double bL = 0, bR = 0, bLR = 0, comp = 0, dag = 0;
        double uE = 0, uLR = 0, uL = 0, uR = 0;
        for (int tick = 0; tick <= ticks; tick++)
        {
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                {
                    bL = Math.Max(bL, (l[i, j] - e[bar - i, j]).Magnitude);
                    bR = Math.Max(bR, (r[i, j] - e[i, bar - j]).Magnitude);
                    bLR = Math.Max(bLR, (lr[i, j] - e[bar - i, bar - j]).Magnitude);
                    comp = Math.Max(comp, (lr[i, j] - l[i, bar - j]).Magnitude);
                    comp = Math.Max(comp, (lr[i, j] - r[bar - i, j]).Magnitude);
                    dag = Math.Max(dag, (l[i, j] - Complex.Conjugate(r[j, i])).Magnitude);
                }
            uE = Math.Max(uE, Math.Abs(1.0 - e.Structure));
            uLR = Math.Max(uLR, Math.Abs(1.0 - lr.Structure));
            uL = Math.Max(uL, Math.Abs(1.0 - l.AntiStructure));
            uR = Math.Max(uR, Math.Abs(1.0 - r.AntiStructure));
            if (tick == ticks) break;
            e.Step(dt); l.Step(dt); r.Step(dt); lr.Step(dt);
        }

        // the vertices are four different matrices (the bridges are relabelings, not identities
        // of the vertex contents): min over pairs of the max entry difference at the final tick.
        var worlds = new[] { e, l, r, lr };
        double separation = double.MaxValue;
        for (int a = 0; a < 4; a++)
            for (int b = a + 1; b < 4; b++)
            {
                double diff = 0;
                for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++)
                        diff = Math.Max(diff, (worlds[a][i, j] - worlds[b][i, j]).Magnitude);
                separation = Math.Min(separation, diff);
            }

        return new LatticeReport(bL, bR, bLR, comp, dag, uE, uLR, uL, uR, separation,
            e.WeightByDisagreement(), l.WeightByDisagreement());
    }

    // ---- the opening law (experiments/LATTICE_OPENING_LAW.md, found playing 2026-07-16): on the
    // cat pair psi(theta) = cos|0..0> + sin|1..1> the entry-wise distance between e and its
    // one-sided reading L has the closed form
    //     opening(t) = max(cos^2, sin^2) - cos*sin * e^(-2*G*t),   G = sum of site rates = N*gamma here,
    // "the heavier sock's weight minus the LIVING spook". Exact because the cat sector is H-dead
    // (a hop needs an excitation beside a hole; |0..0> has none, |1..1> has no hole; the ZZ term
    // gives both ends the same energy), so the e trajectory is pure dephasing and the bridge
    // relabels it entry-wise. Returns the worst |measured - closed form| over the run. ----
    public double OpeningLawDeviation(double theta, double dt, int ticks)
    {
        var w = (World)Parent!;
        int dim = 1 << N, full = dim - 1;
        double c = Math.Cos(theta), s = Math.Sin(theta);

        var e = new Restless(w, N, J, Gamma, zz: zz);
        e.Seed(0, c * c);
        e.Seed(full, s * s);
        e.SeedCoherence(0, full, c * s);
        var l = new Restless(w, N, J, Gamma, antiWatching: true, zz: zz);
        l.SeedRaw(full, 0, c * c);                                   // X^N applied to the rows
        l.SeedRaw(0, full, s * s);
        l.SeedRaw(0, 0, c * s);
        l.SeedRaw(full, full, c * s);

        double bigGamma = N * Gamma;
        double worst = 0;
        for (int tick = 0; tick <= ticks; tick++)
        {
            double opening = 0;
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    opening = Math.Max(opening, (e[i, j] - l[i, j]).Magnitude);
            double predicted = Math.Max(c * c, s * s) - c * s * Math.Exp(-2.0 * bigGamma * (dt * tick));
            worst = Math.Max(worst, Math.Abs(opening - predicted));
            if (tick == ticks) break;
            e.Step(dt); l.Step(dt);
        }
        return worst;
    }

    // ---- the discriminator: run the R reading under the NORMAL rule instead of the turned one;
    // the bridge must break at O(1) (the watching assignment is load-bearing, not decoration). ----
    public double BrokenBridgeR(int seed, double dt, int ticks)
    {
        var w = (World)Parent!;
        int dim = 1 << N, bar = dim - 1;
        var e = new Restless(w, N, J, Gamma, zz: zz);
        e.Seed(seed);
        var wrong = new Restless(w, N, J, Gamma, zz: zz);            // normal rule: the WRONG watching
        wrong.SeedRaw(seed, bar - seed, 1.0);
        double worst = 0;
        for (int tick = 0; tick < ticks; tick++)
        {
            e.Step(dt); wrong.Step(dt);
        }
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                worst = Math.Max(worst, (wrong[i, j] - e[i, bar - j]).Magnitude);
        return worst;
    }
}
