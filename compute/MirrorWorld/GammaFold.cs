using System.Numerics;

namespace MirrorWorld;

// The pair of mirrors on the gamma axis (adopted 2026-07-21, the F134/F139 arc's home-side move):
// the two involutions of the watching parameter and the identity that chains them.
//
//     s   : gamma_l -> -gamma_l          the gain turn (reflection through the unwatched zero;
//                                        negative rates amplify instead of cull)
//     s0  : the anti-watch turn          agreement watched instead of disagreement
//                                        (Restless antiWatching, the Lattice's turned rule)
//
// THE IDENTITY (the generator level, entry-wise, no eigensolver): the Hamiltonian leg never sees
// gamma, and on a cell |i><j| the turned rate is plain arithmetic,
//
//     -2 * sum_{l agrees} gamma_l  =  +2 * sum_{l differs} gamma_l - 2*sigma,   sigma = sum_l gamma_l,
//
// so   L_anti(gamma) = L(-gamma) - 2*sigma*Id   exactly, for ANY site profile, and the trajectory
// wears the shift as a scalar veil:
//
//     rho_anti(t) = e^(-2*sigma*t) * rho_gain(t)     (gain = the same seed run at -gamma).
//
// The anti-watched world IS the gain world in the price veil. On the rate functions r the two turns
// read s: r -> -r and s0: r -> -r - 2*sigma; their composition is the translation r -> r + 2*sigma,
// so <s, s0> is the infinite dihedral group with the full price 2*sigma as its translation step --
// the same two-mirrors-make-a-translation shape F134 carries on the character side (s: mu -> -mu,
// s0: mu -> 22 - mu, step 22). The fixed locus of s0 is r = -sigma: the palindrome center. This is
// a DIFFERENT object from the F1 fold Pi L Pi^-1 = -L - 2*sigma (that one flips the sign of L_H;
// this one keeps H untouched and flips only gamma), and different from the Lattice bridges (those
// relabel the NORMAL trajectory through X^N; here a genuine gain trajectory is propagated and the
// scalar e^(-2*sigma*t) is the whole bridge). Cross-dock: composing this veil with the Lattice's
// one-sided X^N reading gives "the watched world read through the complement = the gain world in
// the price veil" (ReadThroughVeil below).
public sealed class GammaFold : GameObject
{
    public int N { get; }
    public double J { get; }
    readonly double[] gammas;      // the site profile gamma_l (non-uniform by default: per-site is load-bearing)
    readonly double zz;

    public double Sigma => gammas.Sum();

    public GammaFold(World world, int n, double j = 1.0, double[]? siteGammas = null, double zz = 0.0) : base(world)
    {
        N = n;
        J = j;
        this.zz = zz;
        gammas = siteGammas ?? Enumerable.Range(0, n).Select(l => 0.2 + 0.1 * l).ToArray();
    }

    // left: what the fold itself produces.
    public override IReadOnlyList<string> Own => new[] { "mask-identity", "veil", "dihedral" };

    // ---- the mask level: the three exact laws of the rate arithmetic. ----
    // rate masks as plain sums; nothing dynamical, the identity is bit arithmetic per cell.
    double RateNormal(int i, int j)
    {
        double r = 0; int diff = i ^ j;
        for (int l = 0; l < N; l++) if (((diff >> l) & 1) == 1) r -= 2.0 * gammas[l];
        return r;
    }
    double RateAnti(int i, int j)
    {
        double r = 0; int diff = i ^ j;
        for (int l = 0; l < N; l++) if (((diff >> l) & 1) == 0) r -= 2.0 * gammas[l];
        return r;
    }

    public sealed record MaskLawsReport(
        double WorstIdentity,       // worst |r_anti(cell) - (r at -gamma - 2*sigma)| over all cells
        double WorstInvolution,     // worst |s0(s0(r)) - r|: the turn is its own inverse
        double WorstTranslation,    // worst |(s . s0)(r) - (r + 2*sigma)|: two mirrors make the translation
        double Step);               // the translation step 2*sigma (non-vacuity: must be O(1))

    public MaskLawsReport MaskLaws()
    {
        int dim = 1 << N;
        double sigma = Sigma;
        double id = 0, inv = 0, tr = 0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                double r = RateNormal(i, j);
                double ra = RateAnti(i, j);
                id = Math.Max(id, Math.Abs(ra - (-r - 2.0 * sigma)));       // s0 = (gain, then shift by -2*sigma)
                double s0r = -r - 2.0 * sigma;
                inv = Math.Max(inv, Math.Abs((-s0r - 2.0 * sigma) - r));    // s0(s0(r)) = r
                tr = Math.Max(tr, Math.Abs(-s0r - (r + 2.0 * sigma)));      // s(s0(r)) = r + 2*sigma
            }
        return new MaskLawsReport(id, inv, tr, 2.0 * sigma);
    }

    // ---- the trajectory level: the veil law rho_anti(t) = e^(-2*sigma*t) * rho_gain(t). ----
    // Twin RK4: the anti-watched world and the gain world (-gamma, normal watching) from the same
    // seed; the whole bridge is the scalar veil. The discriminator runs the veil against the
    // NORMAL (+gamma) world instead: the gain flip is load-bearing, not decoration.
    public sealed record GammaFoldReport(
        double WorstVeil,           // worst |anti[i,j] - e^(-2*sigma*t) * gain[i,j]| over every probed tick
        double NoveltyRatio,        // gain novelty / anti novelty at the final tick (= e^(+2*sigma*t): the amplification)
        double VertexSeparation,    // max |anti - gain| entry at the final tick (the veil is not vacuous)
        double BrokenFlip);         // worst |anti - e^(-2*sigma*t) * normal| at the final tick (must be O(1))

    public GammaFoldReport Run(int seed, double dt, int ticks)
    {
        var w = (World)Parent!;
        int dim = 1 << N;
        int seed2 = seed ^ 1;                                        // one flipped bit: a live coherence
        var negated = gammas.Select(g => -g).ToArray();

        var anti = new Restless(w, N, J, 0.0, siteGammas: gammas, antiWatching: true, zz: zz);
        var gain = new Restless(w, N, J, 0.0, siteGammas: negated, zz: zz);
        var norm = new Restless(w, N, J, 0.0, siteGammas: gammas, zz: zz);
        foreach (var r in new[] { anti, gain, norm })
        {
            r.Seed(seed, 0.5);
            r.Seed(seed2, 0.5);
            r.SeedCoherence(seed, seed2, 0.5);
        }

        double sigma = Sigma;
        double veil = 0, broken = 0, separation = 0;
        for (int tick = 0; tick <= ticks; tick++)
        {
            double scale = Math.Exp(-2.0 * sigma * (dt * tick));
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    veil = Math.Max(veil, (anti[i, j] - scale * gain[i, j]).Magnitude);
            if (tick == ticks)
            {
                for (int i = 0; i < dim; i++)
                    for (int j = 0; j < dim; j++)
                    {
                        separation = Math.Max(separation, (anti[i, j] - gain[i, j]).Magnitude);
                        broken = Math.Max(broken, (anti[i, j] - scale * norm[i, j]).Magnitude);
                    }
                break;
            }
            anti.Step(dt); gain.Step(dt); norm.Step(dt);
        }
        return new GammaFoldReport(veil, gain.Novelty / anti.Novelty, separation, broken);
    }

    // ---- the cross-dock with the Lattice (the scout's third find): the one-sided X^N reading of
    // the NORMAL world equals the gain-evolved reading in the price veil,
    //     X^N * rho(t)  =  e^(-2*sigma*t) * [gain evolution of X^N * rho(0)].
    // Left side: the committed Lattice bridge L(t)[i,j] = e(t)[~i,j]. Right side: the veil law
    // applied to the one-sided seed. Returns the worst residual over the probed ticks. ----
    public double ReadThroughVeil(int seed, double dt, int ticks)
    {
        var w = (World)Parent!;
        int dim = 1 << N, bar = dim - 1;
        var negated = gammas.Select(g => -g).ToArray();

        var e = new Restless(w, N, J, 0.0, siteGammas: gammas, zz: zz);
        e.Seed(seed);
        var gainRead = new Restless(w, N, J, 0.0, siteGammas: negated, zz: zz);
        gainRead.SeedRaw(bar - seed, seed, 1.0);                     // X^N |s><s| = |~s><s|

        double sigma = Sigma;
        double worst = 0;
        for (int tick = 0; tick <= ticks; tick++)
        {
            double scale = Math.Exp(-2.0 * sigma * (dt * tick));
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    worst = Math.Max(worst, (e[bar - i, j] - scale * gainRead[i, j]).Magnitude);
            if (tick == ticks) break;
            e.Step(dt); gainRead.Step(dt);
        }
        return worst;
    }
}
