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
    public override IReadOnlyList<string> Own => new[] { "mask-identity", "veil", "dihedral", "site-turn" };

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

    // ---- the PER-SITE turns (built here 2026-08-20, and DERIVED here rather than adopted). ----
    // Every other object in this world adopts a result proven elsewhere and names it in its first
    // line. This half names none, because none exists: the arc site_resolved_vacuum_block asked for
    // a per-site sign involution and whether the group it makes with s0 closes, and the answer was
    // worked out here. It stays in genre for the reason README gives (genre, not topic): it is
    // entry-wise rate arithmetic, exact and cell by cell, with no eigensolver and no path object.
    //
    // GammaFold's two mirrors turn the WHOLE profile at once. Turning one site, s_l, is a different
    // animal, and the arithmetic says so before any experiment does: on a cell |i><j| the rate is
    // -2*sum_{l differs} gamma_l, so s_l
    //
    //     leaves every cell where site l AGREES untouched, and
    //     moves every cell where site l DIFFERS by exactly +4*gamma_l.
    //
    // Exactly half the cells (half of the 4^N pairs differ at any fixed bit, at every N), all by
    // the same step, while s moves EVERY cell by twice its own rate, which is what makes s the
    // reflection r -> -r. The composite of the N single turns is s; the turns commute and each is
    // an involution, so ON THE PROFILE they generate (Z/2)^|support|, since turning an unwatched
    // site returns the same profile.
    //
    // ON THE PROFILE THE TURN IS UNCONDITIONAL. It is exact, it is an involution, it asks nothing
    // of the profile. What is conditional is its SHADOW ON THE RATE AXIS, where GammaFold's two
    // mirrors live: a rate is -2 times a subset sum, so s_l is a function of the rate exactly when
    // the rate pins down whether site l disagrees, and that holds iff the profile RESTRICTED TO ITS
    // NONZERO SITES is dissociated (all subset sums of the support distinct). Dissociation of the
    // whole profile is sufficient and NOT necessary: an unwatched site moves by 4*0 = 0, so it may
    // collide freely. Uniform gamma is the extreme failure, and it is the case the repo ran for
    // months: the rate then sees only HOW MANY sites disagree, the spectrum collapses to N+1
    // values, and for N >= 2 no single-site shadow exists.
    //
    // AND THE GROUP DOES NOT CLOSE, which is the arc's second question. On the rate axis s and s0
    // are reflections, slope -1 apiece. Where s_l descends it is PIECEWISE the identity and a
    // translation: identity on the half of the spectrum where site l agrees, a shift by 4*gamma_l
    // on the other half. So it is affine on each piece and on neither the whole, it is no
    // reflection, and (s_l o s0) squared is not even TOTAL on the spectrum, its orbit leaving the
    // set. The infinite dihedral <s, s0> therefore does not grow into a larger reflection group; a
    // partial, piecewise action sits beside it instead.
    //
    // The subset-sum test is EXACT rather than a float comparison, and the reason is worth
    // stating precisely, because the tempting reason is wrong. It is NOT that 1/10 + 2/10 = 3/10
    // slips through: a double array cannot hold those rationals at all, and the values it does hold
    // for 0.1, 0.2 and 0.3 have distinct sums, so both tests agree there. The real gap runs the
    // other way, as a FALSE COLLISION: distinct exact sums can round to the same double, so a float
    // test reports a collision that is not there. {1.0, 1e-20} is the smallest witness, its four
    // subset sums exactly distinct while 1 + 1e-20 evaluates to 1. (The opposite error cannot
    // happen: if two exact sums coincide, the coincident value is representable and the float sum
    // hits it.) Every finite double IS a binary rational, so the sums are compared as integers
    // after a common power-of-two scaling and the question is decided rather than estimated.
    public double[] Turn(IEnumerable<int> sites) => Turn(gammas, sites);

    public static double[] Turn(double[] profile, IEnumerable<int> sites)
    {
        var set = sites.ToHashSet();
        // 0.0 turned stays 0.0 and never -0.0: two equal profiles must compare equal
        // (the repo's group_closure_negative_zero_key trap).
        return profile.Select((g, l) => set.Contains(l) ? (g == 0.0 ? 0.0 : -g) : g).ToArray();
    }

    static double RateOf(double[] g, int i, int j)
    {
        double r = 0; int diff = i ^ j;
        for (int l = 0; l < g.Length; l++) if (((diff >> l) & 1) == 1) r -= 2.0 * g[l];
        return r;
    }

    // every finite double is m * 2^e exactly; scaled to a common exponent the subset sums are
    // integers, so their distinctness is decided with no rounding anywhere.
    static System.Numerics.BigInteger[] ExactScaled(double[] g)
    {
        var mant = new System.Numerics.BigInteger[g.Length];
        var exp = new int[g.Length];
        for (int l = 0; l < g.Length; l++)
        {
            double v = g[l];
            int e = 0;
            while (v != Math.Floor(v)) { v *= 2.0; e++; }        // terminates: doubles are dyadic
            mant[l] = new System.Numerics.BigInteger(v);
            exp[l] = e;
        }
        int max = exp.Length == 0 ? 0 : exp.Max();
        for (int l = 0; l < g.Length; l++)
            mant[l] *= System.Numerics.BigInteger.Pow(2, max - exp[l]);
        return mant;
    }

    // the rate of a cell over the INTEGERS: -2 * sum of the scaled rates where the two indices
    // differ. The float route is the object's own arithmetic and is what the trajectory uses; this
    // one exists because the DESCENT question is about whether two cells SHARE a rate, and float
    // sums answer that wrongly. The default profile is the witness: 0.2 + 0.5 and 0.3 + 0.4 are the
    // same double and different integers, so the float route reported no shadow where the
    // mathematics has one.
    static System.Numerics.BigInteger RateExact(System.Numerics.BigInteger[] scaled, int i, int j)
    {
        System.Numerics.BigInteger r = 0;
        int diff = i ^ j;
        for (int l = 0; l < scaled.Length; l++) if (((diff >> l) & 1) == 1) r -= 2 * scaled[l];
        return r;
    }

    // the common power of two the scaling used, so a scaled integer can be read back as a double
    static int ExactScaleExponent(double[] g)
    {
        int max = 0;
        foreach (var x in g)
        {
            double v = x; int e = 0;
            while (v != Math.Floor(v)) { v *= 2.0; e++; }
            if (e > max) max = e;
        }
        return max;
    }

    static bool SumsDistinct(double[] g)
    {
        var scaled = ExactScaled(g);
        var seen = new HashSet<System.Numerics.BigInteger>();
        for (int m = 0; m < (1 << g.Length); m++)
        {
            System.Numerics.BigInteger t = 0;
            for (int l = 0; l < g.Length; l++) if (((m >> l) & 1) == 1) t += scaled[l];
            if (!seen.Add(t)) return false;
        }
        return true;
    }

    public sealed record SiteTurnReport(
        int CellsWhereSiteDiffers,  // the cells the turn acts on: exactly half of 4^N, at every N
        int CellsWhereSiteAgrees,
        double Step,                // the common step 4*gamma_l (zero at an unwatched site)
        double WorstStepResidual);  // worst |delta - step| where it differs, |delta| where it agrees

    public SiteTurnReport SiteTurn(int site)
    {
        int dim = 1 << N;
        var turned = Turn(new[] { site });
        double step = 4.0 * gammas[site];
        int differs = 0, agrees = 0;
        double worst = 0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                double d = RateOf(turned, i, j) - RateOf(gammas, i, j);
                if (((i ^ j) >> site & 1) == 1) { differs++; worst = Math.Max(worst, Math.Abs(d - step)); }
                else { agrees++; worst = Math.Max(worst, Math.Abs(d)); }
            }
        return new SiteTurnReport(differs, agrees, step, worst);
    }

    public sealed record SiteTurnGroupReport(
        int OrbitSize,                  // distinct profiles reached: 2^|support|
        int Support,                    // how many rates are nonzero
        bool AllInvolutions,            // s_l o s_l = id, by COMPOSING the turn with itself
        bool AllCommute,                // s_a o s_b = s_b o s_a, likewise composed
        double WorstCompositeResidual,  // worst |r(s_{N-1} o ... o s_0) + r|: the composite IS s
        double WorstSigmaResidual,      // worst |sigma(s_S g) - (sigma - 2*sum_S g)|
        bool SupportSumsDistinct,       // the criterion, decided over the integers
        bool WholeProfileSumsDistinct,  // the stronger version: sufficient, not necessary
        bool DescendsToRateAxis,        // s_l is a function of the rate value, at every site
        bool DescendedIsInvolution,     // the map back from the turned profile undoes the map out
        int PieceCount,                 // pieces of the descended map: 2 (identity, and a shift)
        double PieceShift,              // the shift on the moving piece: 4*gamma_l
        bool SquaredWithAntiWatchIsTotal, // is (s_l o s0)^2 total on the spectrum? it is not
        int DistinctRateValues);        // the spectrum (exactly N+1 when gamma is uniform)

    public SiteTurnGroupReport SiteTurnGroup()
    {
        int dim = 1 << N, subsets = 1 << N;
        double sigma = Sigma;

        var orbit = new HashSet<string>();
        double sigWorst = 0;
        for (int m = 0; m < subsets; m++)
        {
            var S = Enumerable.Range(0, N).Where(l => ((m >> l) & 1) == 1).ToArray();
            var g = Turn(S);
            orbit.Add(string.Join(",", g.Select(x => x.ToString("R"))));
            sigWorst = Math.Max(sigWorst, Math.Abs(g.Sum() - (sigma - 2.0 * S.Sum(l => gammas[l]))));
        }

        // COMPOSITION, actually composed. Turn takes a profile now, so s_a o s_b is a composite;
        // an earlier version turned the stored profile twice and re-negated a coordinate by hand,
        // which tested a HashSet and nothing else.
        bool invol = true, comm = true;
        for (int a = 0; a < N; a++)
        {
            var twice = Turn(Turn(gammas, new[] { a }), new[] { a });
            for (int l = 0; l < N; l++) invol &= twice[l] == gammas[l];
            for (int b = 0; b < N; b++)
            {
                var ab = Turn(Turn(gammas, new[] { b }), new[] { a });
                var ba = Turn(Turn(gammas, new[] { a }), new[] { b });
                for (int l = 0; l < N; l++) comm &= ab[l] == ba[l];
            }
        }

        var composed = gammas;
        for (int l = 0; l < N; l++) composed = Turn(composed, new[] { l });
        double compWorst = 0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                compWorst = Math.Max(compWorst,
                    Math.Abs(RateOf(composed, i, j) + RateOf(gammas, i, j)));

        bool suppDistinct = SumsDistinct(gammas.Where(x => x != 0.0).ToArray());
        bool wholeDistinct = SumsDistinct(gammas);

        var scaled = ExactScaled(gammas);
        var spectrum = new HashSet<System.Numerics.BigInteger>();
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                spectrum.Add(RateExact(scaled, i, j));

        bool descends = true, dInvol = true, total = true;
        int pieces = 0;
        System.Numerics.BigInteger shift = 0;
        for (int site = 0; site < N; site++)
        {
            var turnedScaled = ExactScaled(Turn(new[] { site }));
            var map = new Dictionary<System.Numerics.BigInteger, System.Numerics.BigInteger>();
            bool ok = true;
            for (int i = 0; i < dim && ok; i++)
                for (int j = 0; j < dim; j++)
                {
                    var r = RateExact(scaled, i, j);
                    var img = RateExact(turnedScaled, i, j);
                    if (map.TryGetValue(r, out var seen)) { if (seen != img) { ok = false; break; } }
                    else map[r] = img;
                }
            if (!ok) { descends = false; continue; }

            // the map runs from the spectrum of gamma to the spectrum of the turned profile, two
            // DIFFERENT sets, so composing it with itself is not defined; what is an involution is
            // the pair. A collision in the back map is a failure, not a silent overwrite.
            var back = new Dictionary<System.Numerics.BigInteger, System.Numerics.BigInteger>();
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                {
                    var r = RateExact(turnedScaled, i, j);
                    var img = RateExact(scaled, i, j);
                    if (back.TryGetValue(r, out var seen)) { if (seen != img) dInvol = false; }
                    else back[r] = img;
                }
            foreach (var (r, img) in map)
                if (!back.TryGetValue(img, out var b2) || b2 != r) dInvol = false;

            // PIECEWISE, which is the sharp description: the descended map shifts by one of at
            // most two values, 0 and 4*gamma_l, so it is the identity on one piece and a
            // translation on the other. Neither piece is a reflection, and that is why the
            // dihedral cannot absorb it.
            var shifts = map.Select(kv => kv.Value - kv.Key).ToHashSet();
            pieces = Math.Max(pieces, shifts.Count);
            foreach (var sh in shifts) if (sh > shift) shift = sh;

            // (s_l o s0)^2 with s0: r -> -r - 2*sigma, all of it over the integers. Total on the
            // spectrum, or does the orbit leave it?
            var twoSigma = 2 * scaled.Aggregate(System.Numerics.BigInteger.Zero, (a, b) => a + b);
            foreach (var r in spectrum)
            {
                var a1 = -r - twoSigma;
                if (!map.TryGetValue(a1, out var b3)) { total = false; break; }
                var c = -b3 - twoSigma;
                if (!map.ContainsKey(c)) { total = false; break; }
            }
        }
        if (!descends) { dInvol = false; pieces = 0; shift = 0; total = false; }

        // the shift back in the object's own units: the scaling is one common power of two,
        // so dividing by it is exact for every value the object can hold.
        double pieceShift = (double)shift / Math.Pow(2.0, ExactScaleExponent(gammas));

        return new SiteTurnGroupReport(orbit.Count, gammas.Count(x => x != 0.0), invol, comm,
            compWorst, sigWorst, suppDistinct, wholeDistinct, descends, dInvol, pieces, pieceShift,
            total, spectrum.Count);
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
