using System.Numerics;

namespace MirrorWorld;

// The atmosphere (adopted 2026-10-11 from experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md, typed in the
// main repo as CorrelatedZDephasingPriceClaim, Tier 1 derived, live at `inspect --root sharedfield`).
// The Mirror watches every site at one rate gamma and pays its fold one number, 2*N*gamma. That number
// is a featureless atmosphere's: the main repo's gamma_0, which "carries no information ITSELF" and which
// the mirror identity expels to a scalar (docs/THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md, Mechanism B).
// Let the sites' Z watchings CO-VARY, a Kossakowski matrix a = gamma*I + C over the Z_k with C the
// covariances, and the price stops being a number:
//
//   rate(ket,bra) = -1/2 u^T a u = -2*gamma*|D| - 4 * sum_{k<l in D} s_k s_l C_kl     (u = z(ket) - z(bra))
//   price         = K(ket) + K(bra),   K(i) = z(i)^T a z(i)   (the Ising energy; the folds pay it per cell)
//
// D the sites where ket and bra disagree, s_k the ket's Z sign there: the covariance reads the SIGNS of
// the disagreement, and a disagreement at one site alone is blind to it. The rates of a cell and of its
// bra flip sum to -(K(ket) + K(bra)) (the source page's Abstract: the two rates sum to the Ising energy of
// the two configurations); with the Mirror's hopping identity that puts the Mirror's own leg maps (t, f_P,
// f_Q, Klein, with the bipartite gauge) back to exact on the dressed block once each fold pays its cell's
// own price instead of the one number. That step is worked out here, from that sum and the hopping
// identity, at the Mirror's Delta = 0 (the page's gauged linear fold breaks at Delta = 2/5). At the one
// number 2*N*gamma the folds FAIL on every block where K(ket) + K(bra) is not 2*N*gamma, and a test pins it.
// Of the Mirror's price only the MEAN survives: K averages to N*gamma, so 2*N*gamma is the mean cell price
// in every weather, but it is not the deepest rate any more: a cell of full disagreement, u = 2s, pays
// 2*N*gamma + 2*s^T C s (the GHZ cell 2*N*gamma + 4*sum_{k<l} C_kl), and since the s^T C s sum to zero over
// the sign vectors s, every nonzero C has such a cell paying more. The orbit's "one fold: spectrum paid"
// holds entry-wise at a block's own single price where K is one value on both its popcounts, and the
// self-folded trace law at the block's mean price (both worked out here from the same sum).
//
// THE PARENT IS THE MIRROR. What changes is the mirror's price, so the leg maps, the orbit and the price
// as a mean are inherited, and the atmosphere's own outputs are the covariance-dressed rate and price and
// what follows from them. The VARIANCE is the Mirror's too: every site is watched at the Mirror's gamma,
// the Pair rate Mirror.BuildBlock already writes, and the atmosphere's own input is the covariance alone.
// Unequal rates together with a covariance stay in the main repo (GammaFold, Divisor and Cone hold unequal
// local rates).
//
// INTEGER INPUTS, EXACT COMPARISONS. The Mirror's gamma and J must be integers and C an integer matrix, each
// of size at most 10^9, so every block entry is an integer or +-iJ and every residual is compared with
// == 0.0: the same quantity computed twice, a nonzero residual a finding about the construction.
//
// Words, fenced at the door. "Atmosphere" here is the covariance structure of the Z watchings at one rate:
// it generalises the main repo's gamma_0 along the covariances, not along a gamma profile, and at C = 0
// only the inherited rate gamma_0 remains. "Dark" is EndCount's word, a cell commuting with every jump
// (cells a third carrier beside EndCount's strings and BlindSeat's states), which for a positive
// semidefinite a is a cell of rate zero (for an indefinite a, which no jumps realise, the count is of
// zero-rate cells only); a dark CLASS is a set of configurations the noise cannot tell apart, the whole a
// and not the covariance alone (at C = 0 and gamma != 0 every class is one configuration), dark to the
// noise and not thereby to a Hamiltonian. It is not called a code, which would claim protection under the
// whole dynamics, and it is none of the other classes here: not ParameterKlein's anti-palindromic class,
// not Router's {I,Z} <-> {X,Y} letter classes, not Hardness's classes, not EndCount's F158 class, not the
// Formulas' Pi^2-class, and not Witness's Dark record family. A GHZ or Neel "pair" is a pair of
// configurations, whose cross cells are Pairs; the Neel pair is the Z-basis 0101.. and 1010.., the pair the
// main repo's LiouvillianBlockSpectrum calls "both Neel states", not the X-basis Neel state of its
// ZGlobalMirror; independent "pairs" are pairs of sites. The four SINGLE-DISAGREEMENT blocks (0,1), (1,0), (N-1,N), (N,N-1), which a
// covariance cannot touch (nor the population corners (0,0) and (N,N)), are named so here (the source page's
// "four edge blocks"; "edge" is also Crack's band edge
// and EndCount's, "corner" Divisor's and F140's, and Mirror calls its 1x1 blocks one-cell). The lower end
// of the covariance range at one common covariance, eta = -1/(N-1), where every row of a sums to zero, is
// the anti-correlated END, not F158's far end nor CollisionGap's chain or ring end; "centres" are block
// centres, not BlindSeat's centre seat or Divisor's centre room. "Footprint" is the covariance's footprint in the mirror's
// residual, the part that does not collapse to Mechanism B's scalar 2*sigma. "rate" is the Pair's rate
// dressed by the covariance. kappa is a common covariance, not Crack's k = i*kappa; eta = kappa/gamma is
// not CollisionGap's reflection sign. "Centres" holds only the sufficient direction: a block whose two
// popcounts each see one value of K pairs with its image about -(K_p + K_q)/2; that no other block pairs is
// read on the source page, not proved, and is not held here. Which popcounts see one value is the page's
// four-shapes law, adopted beside the enumeration: 0 and N always, 1 and N-1 exactly when every site carries
// the same summed covariance, N/2 exactly when the covariances are additive, C_kl = y_k + y_l, every other
// popcount exactly when they are one number, and there K is the shared field's value at the mean covariance.
// "Additive" is that form, not PROOF_CODIM1_BY_ADDITIVITY's additivity; the page's "cut" is a graph's cut,
// not the Cone's memory cut.
//
// WHAT STAYS OUTSIDE: the global break (F137's one candidate centre), the two ladders, longitudinal fields,
// the XXZ chain and complex or DM Hamiltonians, unequal rates with a covariance, and every eigen-question
// (spectra, Jordan structure, the centres' converse).
public sealed class Atmosphere : GameObject
{
    public Mirror Mirror { get; }
    public int N => Mirror.N;
    public long Gamma { get; }        // the variance, the Mirror's rate: inherited, not owned

    readonly long[,] c;               // the covariances: the atmosphere's own input

    public Atmosphere(Mirror mirror, long[,] covariance) : base(mirror)
    {
        ArgumentNullException.ThrowIfNull(mirror);
        ArgumentNullException.ThrowIfNull(covariance);
        int n = mirror.N;
        if (n > 8)
            throw new ArgumentOutOfRangeException(nameof(mirror), n, "the cell sums run over 4^N cells; N <= 8");
        if (!IsSmallInteger(mirror.Gamma) || !IsSmallInteger(mirror.J))
            throw new ArgumentException("the Mirror's gamma and J must be integers of size at most 10^9, so every comparison is exact", nameof(mirror));
        if (covariance.GetLength(0) != n || covariance.GetLength(1) != n)
            throw new ArgumentException($"the covariance needs N x N = {n} x {n} entries", nameof(covariance));
        for (int k = 0; k < n; k++)
        {
            if (covariance[k, k] != 0)
                throw new ArgumentException("the diagonal is the Mirror's variance; the covariance's diagonal must be 0", nameof(covariance));
            for (int l = 0; l < n; l++)
            {
                if (covariance[k, l] != covariance[l, k])
                    throw new ArgumentException("a Kossakowski matrix is symmetric", nameof(covariance));
                if (covariance[k, l] > 1_000_000_000L || covariance[k, l] < -1_000_000_000L)
                    throw new ArgumentOutOfRangeException(nameof(covariance), covariance[k, l],
                        "|C_kl| <= 10^9 keeps every entry a small exact integer");
            }
        }
        Mirror = mirror;
        Gamma = (long)mirror.Gamma;
        c = (long[,])covariance.Clone();
    }

    static bool IsSmallInteger(double x) => x == Math.Floor(x) && Math.Abs(x) <= 1e9;

    // left: what the atmosphere itself produces.
    public override IReadOnlyList<string> Own => new[] { "rate", "ising", "footprint", "centres", "dark" };

    public long Covariance(int k, int l) => c[k, l];

    // one covariance kappa common to every pair of sites: the shared field
    public static long[,] Uniform(int n, long kappa)
    {
        var m = new long[n, n];
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
                if (k != l) m[k, l] = kappa;
        return m;
    }

    // the Z sign of site k in configuration i: +1 for an empty site, -1 for an excited one (bit k)
    static int Z(int i, int k) => 1 - 2 * ((i >> k) & 1);

    // the covariance part of a cell's rate: -4 * sum over pairs of disagreement sites of s_k s_l C_kl
    long CovarianceRate(int ket, int bra)
    {
        int diff = ket ^ bra;
        long sum = 0;
        for (int k = 0; k < N; k++)
        {
            if (((diff >> k) & 1) == 0) continue;
            for (int l = k + 1; l < N; l++)
                if (((diff >> l) & 1) == 1)
                    sum += Z(ket, k) * Z(ket, l) * c[k, l];
        }
        return -4 * sum;
    }

    // ---- rate: what a cell pays, two routes ----

    // Re lambda of |ket><bra| with no Hamiltonian: the Pair's -2*gamma*k plus the covariance part
    public long Rate(int ket, int bra)
        => -2 * Gamma * BitOperations.PopCount((uint)(ket ^ bra)) + CovarianceRate(ket, bra);

    // the same rate read from the Kossakowski form, -1/2 u^T a u with u = z(ket) - z(bra) (u even, so exact)
    public long RateFromKossakowski(int ket, int bra)
    {
        long quad = 0;
        for (int k = 0; k < N; k++)
            for (int l = 0; l < N; l++)
            {
                long uk = Z(ket, k) - Z(bra, k), ul = Z(ket, l) - Z(bra, l);
                quad += uk * (k == l ? Gamma : c[k, l]) * ul;
            }
        return -quad / 2;
    }

    // ---- ising: the fold's price, cell by cell ----

    // K(i) = z(i)^T a z(i), the Ising energy of configuration i under the atmosphere
    public long Ising(int config)
    {
        long e = 0;
        for (int k = 0; k < N; k++)
            for (int l = 0; l < N; l++)
                e += Z(config, k) * (k == l ? Gamma : c[k, l]) * Z(config, l);
        return e;
    }

    // what the fold legs pay on |ket><bra|: K(ket) + K(bra); the Mirror's 2*N*gamma when C = 0
    public long FoldPrice(int ket, int bra) => Ising(ket) + Ising(bra);

    // ---- the block, dressed: the Mirror's block plus the covariance on the diagonal ----
    public Complex[,] BuildBlock(int p, int q, out List<int> kets, out List<int> bras)
    {
        var l = Mirror.BuildBlock(p, q, out kets, out bras);
        for (int ka = 0; ka < kets.Count; ka++)
            for (int bb = 0; bb < bras.Count; bb++)
            {
                int cell = ka * bras.Count + bb;
                l[cell, cell] += CovarianceRate(kets[ka], bras[bb]);
            }
        return l;
    }

    // the sum of K over all 2^N configurations: 2^N * N * gamma in every atmosphere (the covariance terms cancel)
    public long IsingSum()
    {
        long s = 0;
        for (int i = 0; i < (1 << N); i++) s += Ising(i);
        return s;
    }

    // ---- the four legs on the dressed block: the folds at the cell's own price, transpose and Klein free ----
    public double TransposeResidual(int p, int q) => Mirror.TransposeResidual(BuildBlock, p, q);
    public double BraFoldResidual(int p, int q) => Mirror.BraFoldResidual(BuildBlock, p, q, (a, b) => FoldPrice(a, b));
    public double KetFoldResidual(int p, int q) => Mirror.KetFoldResidual(BuildBlock, p, q, (a, b) => FoldPrice(a, b));
    public double KleinResidual(int p, int q) => Mirror.KleinResidual(BuildBlock, p, q);

    // the bra fold with the price replaced by one number, the Mirror's own form: exact on a block whose two
    // popcounts each see one value of K, at the price K_p + K_q
    public double BraFoldResidualAtScalarPrice(int p, int q, long price)
        => Mirror.BraFoldResidual(BuildBlock, p, q, (_, _) => price);

    // ---- footprint: what no longer cancels ----

    // sum over all 4^N cells of (K(i) + K(j) - 2*N*gamma)^2: the squared Frobenius norm of the dissipator's
    // residual -{K - sigma, .}, read on the matrix units
    public long FootprintByCells()
    {
        int d = 1 << N;
        long twoSigma = 2 * N * Gamma, s = 0;
        var k = new long[d];
        for (int i = 0; i < d; i++) k[i] = Ising(i);
        checked
        {
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                {
                    long r = k[i] + k[j] - twoSigma;
                    s += r * r;
                }
        }
        return s;
    }

    // the closed form 8 * 4^N * sum_{k<l} C_kl^2 (the source page; ResidualNormSquared in the main repo)
    public long FootprintClosedForm()
    {
        long sq = 0;
        for (int k = 0; k < N; k++)
            for (int l = k + 1; l < N; l++)
                sq = checked(sq + c[k, l] * c[k, l]);
        return checked(8L * (1L << (2 * N)) * sq);
    }

    // ---- centres: where the price is one number per block ----

    // the popcounts on which K takes a single value (the empty and the full configuration always do)
    public int[] QualifyingPopcounts()
    {
        var res = new List<int>();
        for (int p = 0; p <= N; p++)
        {
            long? first = null;
            bool one = true;
            for (int i = 0; i < (1 << N) && one; i++)
            {
                if (BitOperations.PopCount((uint)i) != p) continue;
                long e = Ising(i);
                if (first is null) first = e;
                else if (e != first) one = false;
            }
            if (one) res.Add(p);
        }
        return res.ToArray();
    }

    // the same popcounts read from the covariances alone, the page's four-shapes law: no configuration is
    // enumerated. Integer tests: equal row sums; additive, 2*C_kl = y2_k + y2_l with y2 the doubled site values;
    // one number. The tests pin this against QualifyingPopcounts, two routes meeting.
    public int[] QualifyingPopcountsByLaw()
    {
        if (N <= 2) return Enumerable.Range(0, N + 1).ToArray();
        bool constant = true, equalRows = true, additive = true;
        for (int k = 0; k < N; k++)
            for (int l = k + 1; l < N; l++)
                if (c[k, l] != c[0, 1]) constant = false;
        long RowSum(int k) { long r = 0; for (int l = 0; l < N; l++) if (l != k) r += c[k, l]; return r; }
        for (int k = 1; k < N; k++) if (RowSum(k) != RowSum(0)) equalRows = false;
        var y2 = new long[N];
        y2[0] = c[0, 1] + c[0, 2] - c[1, 2];
        for (int k = 1; k < N; k++) y2[k] = 2 * c[0, k] - y2[0];
        for (int k = 0; k < N && additive; k++)
            for (int l = k + 1; l < N; l++)
                if (2 * c[k, l] != y2[k] + y2[l]) { additive = false; break; }
        var res = new List<int>();
        for (int p = 0; p <= N; p++)
            if (p == 0 || p == N || constant || ((p == 1 || p == N - 1) && equalRows) || (2 * p == N && additive))
                res.Add(p);
        return res.ToArray();
    }

    // the summed covariance T = sum_{k<l} C_kl: on every popcount that sees one value, N(N-1)/2 * K_p =
    // N(N-1)/2 * N*gamma - N*T + T*(N-2p)^2, the shared field's value at the mean covariance T/(N(N-1)/2)
    public long TotalCovariance()
    {
        long t = 0;
        for (int k = 0; k < N; k++)
            for (int l = k + 1; l < N; l++) t += c[k, l];
        return t;
    }

    // K_p + K_q, the one fold price of a block whose two popcounts each see one value of K
    public long BlockPrice(int p, int q)
    {
        var ok = QualifyingPopcounts();
        if (Array.IndexOf(ok, p) < 0 || Array.IndexOf(ok, q) < 0)
            throw new InvalidOperationException($"block ({p},{q}): K is not one number on both popcounts");
        return IsingOnPopcount(p) + IsingOnPopcount(q);
    }

    long IsingOnPopcount(int p)
    {
        for (int i = 0; i < (1 << N); i++)
            if (BitOperations.PopCount((uint)i) == p) return Ising(i);
        throw new ArgumentOutOfRangeException(nameof(p), p, "a popcount lies in [0,N]");
    }

    // sigma(p,q) under one common covariance kappa: N*gamma - N*kappa + kappa*(d^2 + m^2), d = p - q,
    // m = p + q - N; block (p,q) and its image pair about -sigma(p,q) (the main repo's BlockSigma)
    public static long SharedBlockSigma(int n, long gamma, long kappa, int p, int q)
    {
        long d = p - q, m = p + q - n;
        return n * gamma - n * kappa + kappa * (d * d + m * m);
    }

    // ---- dark: what the atmosphere does not see ----

    // the cells of rate zero; for a positive semidefinite a exactly the cells commuting with every jump
    public int DarkCount()
    {
        int d = 1 << N, n = 0;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                if (Rate(i, j) == 0) n++;
        return n;
    }

    // the dark classes: the configurations the noise cannot tell apart, grouped by a*z(i), largest first.
    // A cell lies inside a class exactly when a u = 0, so for a positive semidefinite a the dark cells are the
    // cells inside the classes and DarkCount is the sum of the squared class sizes: on the configurations of
    // one class every coherence commutes with every jump. Dark to the noise, not to the dynamics: a Hamiltonian
    // may still move a class (the Mirror's XY chain leaves the GHZ cells still, each alone in its block (0,N)
    // or (N,0), and carries out of the dark cells the Neel pair, which a positive covariance between ring
    // neighbours at half the rate, gamma = 2*kappa, leaves dark at even N >= 4). For an indefinite a, which no
    // jumps realise, a rate can vanish off the classes, and where one does the zero-rate cells outnumber the
    // cells inside them.
    public IReadOnlyList<int[]> DarkClasses()
    {
        var groups = new Dictionary<string, List<int>>();
        for (int i = 0; i < (1 << N); i++)
        {
            var az = new long[N];
            for (int k = 0; k < N; k++)
                for (int l = 0; l < N; l++)
                    az[k] += (k == l ? Gamma : c[k, l]) * Z(i, l);
            string key = string.Join(",", az);
            if (!groups.TryGetValue(key, out var members)) groups[key] = members = new List<int>();
            members.Add(i);
        }
        return groups.Values.OrderByDescending(m => m.Count).ThenBy(m => m[0]).Select(m => m.ToArray()).ToList();
    }
}
