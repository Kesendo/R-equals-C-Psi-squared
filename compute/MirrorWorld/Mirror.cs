using System.Numerics;

namespace MirrorWorld;

// The first mirror IN the world of mirrors (adopted 2026-07-03 from the fold-lattice lemma,
// docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md section 7). Until now this world held only mirror
// SIGNATURES (the Pi^2 parity signs, the Klein cells); the mirror itself lived elsewhere. Here it is,
// as an object -- and it is not an eigen-story but a REARRANGEMENT: on the block lattice (p,q) the
// three legs
//
//   t   : swap bra and ket        (p,q) -> (q,p)      spectrum unchanged
//   f_P : flip the bra            (p,q) -> (p,N-q)    spectrum folded: lambda -> -lambda - price
//   f_Q : flip the ket            (p,q) -> (N-p,q)    spectrum folded: same price
//
// generate a group of eight, and EVERY leg is an exact entry-wise identity between the two block
// matrices at the SAME coupling -- checkable cell by cell, machine zero, no eigensolver:
//
//   L(q,p)  [t(i),t(j)] =  e_i e_j * L(p,q)[i,j]
//   L(p,N-q)[m(i),m(j)] = -e_i e_j * L(p,q)[i,j] - price * delta_ij     (price = 2*N*gamma)
//
// with e the bipartite gauge sign (excitations on odd sites) and m/t the index relabelings. The fold
// legs carry the affine cocycle s(lambda) = -lambda - price: the mirror shows the PARTNER block
// running BACKWARD, and matching it costs the deepest rate. A block the fold maps to itself (even N,
// q = N/2) must pay the price out of its own pocket: trace L = -(price/2)*dim, no computation needed.
// States and their mirrors live here; the paths (the braid, the monodromy) stay outside -- a catalog
// cannot hold a way, only what the way leaves behind.
public sealed class Mirror : GameObject
{
    public int N { get; }
    public double J { get; }
    public double Gamma { get; }

    readonly World world;
    readonly (int a, int b)[] bonds;                    // the gauge needs a bipartite geometry: the chain

    public Mirror(World world, int n, double j, double gamma) : base(world)
    {
        this.world = world;
        N = n;
        J = j;
        Gamma = gamma;
        bonds = Topology.Chain(n);
    }

    // the fold's price: the deepest rate in the world, 2*N*gamma (every site watched at once).
    public double Price => 2.0 * N * Gamma;

    // left: what the mirror itself produces.
    public override IReadOnlyList<string> Own => new[] { "legs", "price", "orbit" };

    // ---- the block pencil, built from the atoms this world already owns ----
    // basis = (ket config a of popcount p) x (bra config b of popcount q); diagonal = the Pair's own
    // rate -2*gamma*k(a,b); ket hop -iJ, bra hop +iJ per bond (the handshake read on one side each,
    // exactly Restless's -i[H,rho] restricted to the block).
    public Complex[,] BuildBlock(int p, int q, out List<int> kets, out List<int> bras)
    {
        kets = Configs(N, p);
        bras = Configs(N, q);
        int dim = kets.Count * bras.Count;
        var braIndex = new Dictionary<int, int>();
        for (int i = 0; i < bras.Count; i++) braIndex[bras[i]] = i;
        var ketIndex = new Dictionary<int, int>();
        for (int i = 0; i < kets.Count; i++) ketIndex[kets[i]] = i;

        var l = new Complex[dim, dim];
        for (int ka = 0; ka < kets.Count; ka++)
            for (int bb = 0; bb < bras.Count; bb++)
            {
                int col = ka * bras.Count + bb;
                l[col, col] = new Pair(world, kets[ka], bras[bb], Gamma).Rate;     // the atom, read directly
                foreach (int a2 in Hops(kets[ka]))
                    l[ketIndex[a2] * bras.Count + bb, col] += -Complex.ImaginaryOne * J;
                foreach (int b2 in Hops(bras[bb]))
                    l[ka * bras.Count + braIndex[b2], col] += Complex.ImaginaryOne * J;
            }
        return l;
    }

    // ---- the four legs, each a worst entry-wise residual over the two blocks (0 = the identity holds) ----

    // t: (a,b) -> (b,a). L(q,p)[t(i),t(j)] = e_i e_j L(p,q)[i,j].
    public double TransposeResidual(int p, int q)
        => LegResidual(p, q, q, p, (a, b) => (b, a), sign: +1, shift: 0.0);

    // f_P: (a,b) -> (a, complement of b). L(p,N-q)[m(i),m(j)] = -e_i e_j L(p,q)[i,j] - price*delta.
    public double BraFoldResidual(int p, int q)
        => LegResidual(p, q, p, N - q, (a, b) => (a, Complement(b)), sign: -1, shift: Price);

    // f_Q: the ket-side twin.
    public double KetFoldResidual(int p, int q)
        => LegResidual(p, q, N - p, q, (a, b) => (Complement(a), b), sign: -1, shift: Price);

    // f_P f_Q = the Klein full flip: both complemented, same q-couplings, no gauge, no price.
    public double KleinResidual(int p, int q)
    {
        var src = BuildBlock(p, q, out var kets, out var bras);
        var tgt = BuildBlock(N - p, N - q, out var tk, out var tb);
        var map = BuildMap(kets, bras, tk, tb, (a, b) => (Complement(a), Complement(b)));
        double worst = 0;
        int dim = kets.Count * bras.Count;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                worst = Math.Max(worst, (tgt[map[i], map[j]] - src[i, j]).Magnitude);
        return worst;
    }

    double LegResidual(int p, int q, int tp, int tq, Func<int, int, (int, int)> relabel, int sign, double shift)
    {
        var src = BuildBlock(p, q, out var kets, out var bras);
        var tgt = BuildBlock(tp, tq, out var tk, out var tb);
        var map = BuildMap(kets, bras, tk, tb, relabel);
        int dim = kets.Count * bras.Count;
        var gauge = new int[dim];
        for (int ka = 0; ka < kets.Count; ka++)
            for (int bb = 0; bb < bras.Count; bb++)
                gauge[ka * bras.Count + bb] = GaugeSign(kets[ka]) * GaugeSign(bras[bb]);
        double worst = 0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                var expected = sign * gauge[i] * gauge[j] * src[i, j] - (i == j ? shift : 0.0);
                worst = Math.Max(worst, (tgt[map[i], map[j]] - expected).Magnitude);
            }
        return worst;
    }

    int[] BuildMap(List<int> kets, List<int> bras, List<int> tk, List<int> tb, Func<int, int, (int, int)> relabel)
    {
        var tkIndex = new Dictionary<int, int>();
        for (int i = 0; i < tk.Count; i++) tkIndex[tk[i]] = i;
        var tbIndex = new Dictionary<int, int>();
        for (int i = 0; i < tb.Count; i++) tbIndex[tb[i]] = i;
        var map = new int[kets.Count * bras.Count];
        for (int ka = 0; ka < kets.Count; ka++)
            for (int bb = 0; bb < bras.Count; bb++)
            {
                var (a2, b2) = relabel(kets[ka], bras[bb]);
                map[ka * bras.Count + bb] = tkIndex[a2] * tb.Count + tbIndex[b2];
            }
        return map;
    }

    // ---- the orbit: the eight images of a block, with the fold parity (how many fold legs it took) ----
    public static (int P, int Q, int FoldParity)[] OrbitImages(int n, int p, int q) => new[]
    {
        (p, q, 0), (q, p, 0), (n - p, n - q, 0), (n - q, n - p, 0),                 // no fold: spectrum kept
        (p, n - q, 1), (n - p, q, 1), (n - q, p, 1), (q, n - p, 1),                 // one fold: spectrum paid
    };

    // the block's representative: the lexicographic least of its orbit (the fundamental domain), and
    // the fold parity of the path from the representative to this block (a self-folded block reads 0:
    // its two parities coincide because its spectrum is s-symmetric).
    public (int P, int Q, int FoldParity) Representative(int p, int q)
    {
        var rep = OrbitImages(N, p, q).OrderBy(x => x.P).ThenBy(x => x.Q).First();
        int parity = OrbitImages(N, rep.P, rep.Q).First(x => x.P == p && x.Q == q).FoldParity;
        return (rep.P, rep.Q, parity);
    }

    public int OrbitCount()
    {
        var reps = new HashSet<(int, int)>();
        for (int p = 0; p <= N; p++)
            for (int q = 0; q <= N; q++)
            {
                var img = OrbitImages(N, p, q).OrderBy(x => x.P).ThenBy(x => x.Q).First();
                reps.Add((img.P, img.Q));
            }
        return reps.Count;
    }

    // ---- the self-folded price (even N, q = N/2): the block is its own fold image, so it pays the
    // price out of its own trace: L + (its mirror image) = -price*I, hence trace L = -(price/2)*dim. ----
    public (double Trace, double Law) SelfFoldedTrace(int p)
    {
        if (N % 2 != 0) throw new InvalidOperationException("self-folded blocks need even N (q = N/2).");
        var l = BuildBlock(p, N / 2, out var kets, out var bras);
        double tr = 0;
        for (int i = 0; i < kets.Count * bras.Count; i++) tr += l[i, i].Real;
        return (tr, -(Price / 2.0) * kets.Count * bras.Count);
    }

    // ---- the trajectory fold: the mirror shows the partner world running backward, at the price ----
    // If x(t) solves x-dot = L(p,q) x, then w(t) = exp(price*t) * (fold of x(t)) solves
    // w-dot = -L(p,N-q) w: the PARTNER block, time reversed, the growth paid for by exp(-price*t).
    // Both sides are run independently (RK4) and compared tick by tick.
    public (double[] T, double[] NormX, double[] NormW, double WorstResidual) TrajectoryFold(
        int p, int q, double dt, int ticks)
    {
        var l = BuildBlock(p, q, out var kets, out var bras);
        var lp = BuildBlock(p, N - q, out var tk, out var tb);
        var map = BuildMap(kets, bras, tk, tb, (a, b) => (a, Complement(b)));
        int dim = kets.Count * bras.Count;
        var gauge = new int[dim];
        for (int ka = 0; ka < kets.Count; ka++)
            for (int bb = 0; bb < bras.Count; bb++)
                gauge[ka * bras.Count + bb] = GaugeSign(kets[ka]) * GaugeSign(bras[bb]);

        // seed: an even spread over the block (any vector does; the identity is exact).
        var x = new Complex[dim];
        for (int i = 0; i < dim; i++) x[i] = 1.0 / Math.Sqrt(dim);
        var w = FoldVector(x, map, gauge, dim);          // w(0) = fold of x(0), price not yet accrued

        var minusLp = new Complex[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                minusLp[i, j] = -lp[i, j];

        var ts = new double[ticks + 1];
        var nx = new double[ticks + 1];
        var nw = new double[ticks + 1];
        double worst = 0;
        for (int tick = 0; tick <= ticks; tick++)
        {
            double t = tick * dt;
            ts[tick] = t;
            nx[tick] = Norm(x);
            nw[tick] = Norm(w);
            var predicted = FoldVector(x, map, gauge, dim);
            double scale = Math.Exp(Price * t);
            double res = 0;
            for (int i = 0; i < dim; i++) res = Math.Max(res, (w[i] - scale * predicted[i]).Magnitude);
            worst = Math.Max(worst, res / Math.Max(Norm(w), 1e-300));
            if (tick == ticks) break;
            x = Rk4(l, x, dt, dim);
            w = Rk4(minusLp, w, dt, dim);
        }
        return (ts, nx, nw, worst);
    }

    static Complex[] FoldVector(Complex[] x, int[] map, int[] gauge, int dim)
    {
        var y = new Complex[dim];
        for (int i = 0; i < dim; i++) y[map[i]] = gauge[i] * x[i];
        return y;
    }

    static Complex[] Rk4(Complex[,] m, Complex[] x, double dt, int dim)
    {
        var k1 = Mul(m, x, dim);
        var k2 = Mul(m, Axpy(x, k1, dt / 2, dim), dim);
        var k3 = Mul(m, Axpy(x, k2, dt / 2, dim), dim);
        var k4 = Mul(m, Axpy(x, k3, dt, dim), dim);
        var r = new Complex[dim];
        for (int i = 0; i < dim; i++) r[i] = x[i] + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        return r;
    }

    static Complex[] Mul(Complex[,] m, Complex[] x, int dim)
    {
        var r = new Complex[dim];
        for (int i = 0; i < dim; i++)
        {
            Complex s = Complex.Zero;
            for (int j = 0; j < dim; j++) s += m[i, j] * x[j];
            r[i] = s;
        }
        return r;
    }

    static Complex[] Axpy(Complex[] a, Complex[] b, double s, int dim)
    {
        var r = new Complex[dim];
        for (int i = 0; i < dim; i++) r[i] = a[i] + s * b[i];
        return r;
    }

    static double Norm(Complex[] x)
    {
        double s = 0;
        foreach (var v in x) s += v.Real * v.Real + v.Imaginary * v.Imaginary;
        return Math.Sqrt(s);
    }

    // ---- small shared pieces ----
    int Complement(int c) => ((1 << N) - 1) ^ c;

    // the bipartite gauge sign: parity of the excitations sitting on odd sites (the chain's two colors).
    static int GaugeSign(int c)
    {
        int sign = 1;
        for (int l = 1; l < 31; l += 2)
            if (((c >> l) & 1) == 1) sign = -sign;
        return sign;
    }

    int[] Hops(int c)
    {
        var moves = new List<int>();
        foreach (var (a, b) in bonds)
        {
            bool ea = ((c >> a) & 1) == 1, eb = ((c >> b) & 1) == 1;
            if (ea != eb) moves.Add(c ^ (1 << a) ^ (1 << b));
        }
        return moves.ToArray();
    }

    static List<int> Configs(int n, int w)
    {
        var res = new List<int>();
        for (int m = 0; m < (1 << n); m++)
            if (BitOperations.PopCount((uint)m) == w) res.Add(m);
        return res;
    }
}
