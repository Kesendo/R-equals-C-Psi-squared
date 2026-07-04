using System.Numerics;

namespace MirrorWorld;

// The parameter-side Klein V4 (adopted 2026-07-04 from F91 + F92 + F93,
// docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md and its J / h twins): on a parameter axis --
// gamma per site (F91), J per bond (F92), longitudinal detuning h per site (F93) -- the
// palindromic mirror F71 (reverse) and the anti-palindromic reshuffle R90 (x_l -> 2 avg - x_{L-1-l})
// are two commuting involutions generating the Klein four-group V4 = Z2 x Z2, the order-2 shadow
// of the operator-side Z4. The anti-palindromic class (every pair-sum x_l + x_{L-1-l} = 2 avg) is
// exactly the FIXED-POINT set of R90.
//
// The adopted law, in its sharper entry-wise form (the proofs' actual conclusion, and the reason
// this object fits MirrorWorld's genre -- no eigensolver anywhere): the F71-refined DIAGONAL-block
// matrix elements of L = -i[H,.] + D depend only on the F71 pair-sums (S_l = gamma_l + gamma_{N-1-l},
// T_b = J_b + J_{N-2-b}, h_l + h_{N-1-l}), never on individual values or pair-differences. So any
// two profiles in the anti-palindromic orbit give IDENTICAL diagonal blocks, cell for cell -- the
// spectral invariance (the registry's F91-F93 statement) is the trivial corollary. The
// pair-differences live only in the F71 cross-blocks: the breaking sits in the eigenvectors.
public sealed class ParameterKlein : GameObject
{
    public int N { get; }

    public ParameterKlein(World world, int n) : base(world) => N = n;

    // left: what the parameter Klein itself produces.
    public override IReadOnlyList<string> Own => new[] { "v4", "orbit", "blocks" };

    // ---- the V4 on a parameter axis (any length: N sites for gamma/h, N-1 bonds for J) ----
    public static double[] Mirror(double[] x)                       // F71: the palindromic reverse
    {
        var r = new double[x.Length];
        for (int i = 0; i < x.Length; i++) r[i] = x[x.Length - 1 - i];
        return r;
    }

    public static double[] Reshuffle(double[] x)                    // R90: 2 avg - reverse
    {
        double m = 2.0 * x.Average();
        var r = new double[x.Length];
        for (int i = 0; i < x.Length; i++) r[i] = m - x[x.Length - 1 - i];
        return r;
    }

    public static double[] MeanReflection(double[] x)               // F71 o R90: per-site 2 avg - x
    {
        double m = 2.0 * x.Average();
        return x.Select(v => m - v).ToArray();
    }

    // ---- the entry-wise law: worst |difference| between the F71-refined diagonal blocks of two
    // parameter sets, over every joint-popcount block (p,q). Zero iff the two sets agree on all
    // refined diagonal entries -- for two members of the anti-palindromic orbit, machine zero. ----
    public double DiagonalBlocksResidual(double[] gammasA, double[] jsA, double[] hsA,
                                         double[] gammasB, double[] jsB, double[] hsB)
    {
        double worst = 0;
        for (int p = 0; p <= N; p++)
            for (int q = 0; q <= N; q++)
            {
                var (la, perm) = BuildBlock(p, q, gammasA, jsA, hsA);
                var (lb, _) = BuildBlock(p, q, gammasB, jsB, hsB);
                var reps = Reps(perm);
                foreach (int i in reps)
                    foreach (int j in reps)
                    {
                        worst = Math.Max(worst, (EvenElement(la, i, j, perm) - EvenElement(lb, i, j, perm)).Magnitude);
                        if (perm[i] > i && perm[j] > j)
                            worst = Math.Max(worst, (OddElement(la, i, j, perm) - OddElement(lb, i, j, perm)).Magnitude);
                    }
            }
        return worst;
    }

    // ---- where the breaking lives: the worst F71 cross-block entry (even <-> odd). Zero iff F71
    // is an exact L-symmetry (the palindromic case); O(1) on a non-palindromic anti-palindromic
    // profile, where the diagonal blocks still match. ----
    public double CrossBlockNorm(double[] gammas, double[] js, double[] hs)
    {
        double worst = 0;
        for (int p = 0; p <= N; p++)
            for (int q = 0; q <= N; q++)
            {
                var (l, perm) = BuildBlock(p, q, gammas, js, hs);
                var reps = Reps(perm);
                foreach (int i in reps)
                    foreach (int j in reps)
                    {
                        if (perm[j] > j)
                            worst = Math.Max(worst, CrossElement(l, i, j, perm).Magnitude);
                        if (perm[i] > i)
                            worst = Math.Max(worst, CrossElementT(l, i, j, perm).Magnitude);
                    }
            }
        return worst;
    }

    // ---- the block pencil with site-resolved parameters (the same atoms Mirror builds from, plus
    // per-bond J and per-site h): diagonal = -2 sum_l gamma_l [a_l != b_l] - i (E_h(a) - E_h(b)),
    // ket hop -i J_b, bra hop +i J_b on the chain bond b = (l, l+1). ----
    (Complex[,] L, int[] Perm) BuildBlock(int p, int q, double[] gammas, double[] js, double[] hs)
    {
        var kets = Configs(N, p);
        var bras = Configs(N, q);
        var ketIndex = new Dictionary<int, int>();
        for (int i = 0; i < kets.Count; i++) ketIndex[kets[i]] = i;
        var braIndex = new Dictionary<int, int>();
        for (int i = 0; i < bras.Count; i++) braIndex[bras[i]] = i;

        int dim = kets.Count * bras.Count;
        var l = new Complex[dim, dim];
        var perm = new int[dim];
        for (int ka = 0; ka < kets.Count; ka++)
            for (int bb = 0; bb < bras.Count; bb++)
            {
                int a = kets[ka], b = bras[bb], col = ka * bras.Count + bb;
                double watching = 0, detune = 0;
                for (int s = 0; s < N; s++)
                {
                    if (((a >> s) & 1) != ((b >> s) & 1)) watching += 2.0 * gammas[s];
                    detune += hs[s] * ((1 - 2 * ((a >> s) & 1)) - (1 - 2 * ((b >> s) & 1)));
                }
                l[col, col] = new Complex(-watching, -detune);
                for (int bond = 0; bond < N - 1; bond++)
                {
                    int mask = (1 << bond) | (1 << (bond + 1));
                    if (((a >> bond) & 1) != ((a >> (bond + 1)) & 1))
                        l[ketIndex[a ^ mask] * bras.Count + bb, col] += -Complex.ImaginaryOne * js[bond];
                    if (((b >> bond) & 1) != ((b >> (bond + 1)) & 1))
                        l[ka * bras.Count + braIndex[b ^ mask], col] += Complex.ImaginaryOne * js[bond];
                }
                perm[col] = ketIndex[Reverse(a)] * bras.Count + braIndex[Reverse(b)];
            }
        return (l, perm);
    }

    // ---- the F71-refined elements (the bond mirror as bit reversal; orbit reps i <= perm[i]):
    // even vector (d_i + d_pi)/sqrt2 (or d_i when self-mirrored), odd vector (d_i - d_pi)/sqrt2. ----
    static List<int> Reps(int[] perm)
    {
        var reps = new List<int>();
        for (int i = 0; i < perm.Length; i++) if (perm[i] >= i) reps.Add(i);
        return reps;
    }

    static Complex EvenElement(Complex[,] m, int i, int j, int[] perm)
    {
        int pi = perm[i], pj = perm[j];
        if (pi == i && pj == j) return m[i, j];
        if (pi == i) return (m[i, j] + m[i, pj]) / Math.Sqrt(2.0);
        if (pj == j) return (m[i, j] + m[pi, j]) / Math.Sqrt(2.0);
        return 0.5 * (m[i, j] + m[i, pj] + m[pi, j] + m[pi, pj]);
    }

    static Complex OddElement(Complex[,] m, int i, int j, int[] perm)
        => 0.5 * (m[i, j] - m[i, perm[j]] - m[perm[i], j] + m[perm[i], perm[j]]);

    static Complex CrossElement(Complex[,] m, int i, int j, int[] perm)          // even i, odd j
    {
        int pi = perm[i], pj = perm[j];
        if (pi == i) return (m[i, j] - m[i, pj]) / Math.Sqrt(2.0);
        return 0.5 * (m[i, j] - m[i, pj] + m[pi, j] - m[pi, pj]);
    }

    static Complex CrossElementT(Complex[,] m, int i, int j, int[] perm)         // odd i, even j
    {
        int pi = perm[i], pj = perm[j];
        if (pj == j) return (m[i, j] - m[pi, j]) / Math.Sqrt(2.0);
        return 0.5 * (m[i, j] + m[i, pj] - m[pi, j] - m[pi, pj]);
    }

    int Reverse(int c)
    {
        int r = 0;
        for (int s = 0; s < N; s++)
            if (((c >> s) & 1) == 1) r |= 1 << (N - 1 - s);
        return r;
    }

    static List<int> Configs(int n, int w)
    {
        var res = new List<int>();
        for (int m = 0; m < (1 << n); m++)
            if (BitOperations.PopCount((uint)m) == w) res.Add(m);
        return res;
    }
}
