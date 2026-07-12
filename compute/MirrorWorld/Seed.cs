using System.Numerics;

namespace MirrorWorld;

// The within-block self-dual seed (adopted 2026-07-07 from the F89 seed-existence reduction,
// experiments/F89_SEED_EXISTENCE_REDUCTION.md + docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md).
//
// Mirror.cs holds the BETWEEN-block folds (t / f_P / f_Q map one block onto another, exactly). This
// holds the WITHIN-block self-duality the folds leave untouched: the point where a state meets the
// mirror's null -- self-orthogonal under the transpose bilinear form, v^T v = 0 -- a DEFECTIVE seed.
// That seed is the static source the shadow (a large projector norm, non-normality) and the i^4=1
// eigenvector holonomy leave behind in the main repo; those readings are EIGEN-stories and PATHS, so
// they stay outside. What stays inside is the seed's EXISTENCE, held as a COUNT with no eigensolver.
//
// The (1,2) block is the affine pencil L(q) = A + q*C: A = diag(-2*gamma*n_diff) (the Pair rate, two
// rungs n_diff in {1,3}); C = the hop part (ket hop -i, bra hop +i per bond, so C/i is a 0/+-1 integer
// matrix). Count the real strands FORCED to leave the axis as q grows (F89):
//
//   r(inf) = nullity(C)                            (the strands already merged at q = infinity)
//   r(0+)  = sum over rungs of nullity(P_r C P_r)  (the strands born real at q = 0+, per A-rung)
//   surplus = r(0+) - r(inf) = the seed count.
//
// Proven surplus = N-1 for every ODD N: the chain's unmirrorable middle site forces each reflection
// half to face itself, making eigenvalues real, which must then collide (Kato: a simple discriminant
// zero is defective). The COUNT is F89's theorem; that each collision is defective rather than the
// non-generic semisimple beta-exotic rests on Kato's simple zero modulo the still-open codim-2
// genericity (census-clean through N=11). This object holds the count, never the locus. For EVEN N there is no unmirrorable seat, every seat mirrors to another, and the
// surplus is 0 -- no real strand is FORCED off the axis (the count is a forced lower bound; whether an
// accidental seed exists off it is a separate question this count does not settle). Ranks over GF(p) by exact integer
// Gaussian elimination (the genre of Hardness's GF(2)[x]); a rank mod p can only DROP at a bad prime, so
// the max over two primes pins the true rank (exact for the N run here, verified against F89's table).
//
// SCOPE (do not naively extend; verified from below 2026-07-08). This surplus = defective-seed count is
// F89's theorem for the OPEN CHAIN (1,2) block ONLY. The same nullity construction on other blocks or
// topologies is NOT the seed count: on the RING it is an artifact (surplus -4 at N=4 where the actual
// real-eigenvalue count-change is 0 -- the cyclic structure breaks the rung-nullity = real-count
// correspondence the chain proof rests on); on the STAR the real strands do drain (12 -> 0 at N=4, even N,
// the hub being its own kind of unmirrorable seat), but whether those transitions are DEFECTIVE needs the
// Riesz EP-character instrument, not this count. A topology/block map of seeds is a separate careful arc.
public sealed class Seed : GameObject
{
    public int N { get; }
    public double Gamma { get; }
    readonly (int a, int b)[] bonds;

    // two large primes; the rank over Q is the MAX of the ranks mod p (a prime is "bad" only when it
    // divides a maximal minor, vanishingly rare), so max over two pins the exact integer rank.
    static readonly long[] Primes = { 2147483647L, 999999937L };

    public Seed(World world, int n, double gamma) : base(world)
    {
        N = n;
        Gamma = gamma;
        bonds = Topology.Chain(n);
    }

    // left: what the seed object produces itself.
    public override IReadOnlyList<string> Own => new[] { "count", "rungs" };

    // the seed count and its two halves, plus the per-rung (n_diff, dim, nullity) breakdown.
    public (int RInf, int R0, int Surplus, (int NDiff, int Dim, int Nullity)[] Parts) Count()
    {
        var (rung, M, dim) = BuildPencil();
        var all = Enumerable.Range(0, dim).ToArray();
        int rInf = dim - Rank(M, all);

        var byRung = new SortedDictionary<int, List<int>>();
        for (int i = 0; i < dim; i++)
        {
            if (!byRung.TryGetValue(rung[i], out var list)) { list = new List<int>(); byRung[rung[i]] = list; }
            list.Add(i);
        }
        int r0 = 0;
        var parts = new List<(int, int, int)>();
        foreach (var (nd, idx) in byRung)
        {
            int nl = idx.Count - Rank(M, idx.ToArray());
            parts.Add((nd, idx.Count, nl));
            r0 += nl;
        }
        return (rInf, r0, r0 - rInf, parts.ToArray());
    }

    // ---- the fusion-resonance count, closed (adopted 2026-07-12; same doc, Piece 3 + the 2026-07-10c
    // resonance criterion). PROVED for every N, even N included: nullity(C) = 3*Z3, where Z3 counts the
    // unordered mode triples a < b < c with lam_a + lam_b + lam_c = 0 (Step 4's cyclotomic integrality
    // kills every degenerate relation, so r(inf) is exactly divisible by 3). For ODD N the Conway-Jones
    // closed form: Z3 = (N-1)/2 + [3|n]*(n/3 - 2) + 2*[15|n], n = N+1 (criterion and TRIV/PENT counts
    // proved; the ROT3 multiplicity n/3 - 2 verified, not derived). RESONANT -- the kernel exceeding the
    // always-present TRIV floor -- iff 3 | N+1 and N >= 11; the next resonant N after 17 is 23, not 29.
    // The even rows carry no closed form here (n odd, no TRIV family); their Z3 is still forced by the
    // divisibility (N = 8: r(inf) = 6, Z3 = 2, the 2cos(pi/9) family -- the smallest even laboratory of
    // the twinning arc, run 2026-07-12).
    public static int Z3FromRInf(int rInf) =>
        rInf % 3 == 0 ? rInf / 3
                      : throw new InvalidOperationException($"r(inf) = {rInf} not divisible by 3: the 3*Z3 theorem broke");

    public static int Z3ClosedFormOdd(int nSites)
    {
        if (nSites % 2 == 0) throw new ArgumentException("the closed form is the odd-N Conway-Jones count", nameof(nSites));
        int n = nSites + 1;
        return (nSites - 1) / 2 + (n % 3 == 0 ? n / 3 - 2 : 0) + (n % 15 == 0 ? 2 : 0);
    }

    // the criterion is the ODD-N statement (the twinning arc's "kernel exceeds the TRIV floor"); at
    // even N there is no TRIV floor and the notion does not apply, so the predicate is false there.
    public static bool IsResonant(int nSites) => nSites % 2 == 1 && nSites >= 11 && (nSites + 1) % 3 == 0;

    // ---- the (1,2) pencil, from the atoms: rung n_diff per basis index, and the integer hop M = C/i ----
    (int[] rung, int[,] M, int dim) BuildPencil()
    {
        var kets = Configs(N, 1);
        var bras = Configs(N, 2);
        var ki = new Dictionary<int, int>();
        for (int i = 0; i < kets.Count; i++) ki[kets[i]] = i;
        var bi = new Dictionary<int, int>();
        for (int i = 0; i < bras.Count; i++) bi[bras[i]] = i;
        int nb = bras.Count, dim = kets.Count * nb;

        var rung = new int[dim];
        var M = new int[dim, dim];
        for (int ka = 0; ka < kets.Count; ka++)
            for (int bb = 0; bb < nb; bb++)
            {
                int col = ka * nb + bb;
                rung[col] = BitOperations.PopCount((uint)(kets[ka] ^ bras[bb]));   // n_diff, the Pair's read
                foreach (int a2 in Hops(kets[ka]))
                    M[ki[a2] * nb + bb, col] += -1;                                // ket hop: C/i = -1
                foreach (int b2 in Hops(bras[bb]))
                    M[ka * nb + bi[b2], col] += +1;                                // bra hop: C/i = +1
            }
        return (rung, M, dim);
    }

    // rank of the submatrix M[idx, idx] over Q, as max over the primes (exact).
    static int Rank(int[,] M, int[] idx)
    {
        int r = 0;
        foreach (long p in Primes) r = Math.Max(r, RankModP(M, idx, p));
        return r;
    }

    // Gaussian elimination of the submatrix M[idx, idx] over GF(p); returns its rank.
    static int RankModP(int[,] M, int[] idx, long p)
    {
        int d = idx.Length;
        var a = new long[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                a[i, j] = ((M[idx[i], idx[j]] % p) + p) % p;

        int rank = 0;
        for (int col = 0; col < d && rank < d; col++)
        {
            int piv = -1;
            for (int row = rank; row < d; row++)
                if (a[row, col] != 0) { piv = row; break; }
            if (piv < 0) continue;
            for (int j = 0; j < d; j++) (a[rank, j], a[piv, j]) = (a[piv, j], a[rank, j]);
            long inv = ModInverse(a[rank, col], p);
            for (int j = 0; j < d; j++) a[rank, j] = a[rank, j] * inv % p;
            for (int row = 0; row < d; row++)
            {
                if (row == rank || a[row, col] == 0) continue;
                long f = a[row, col];
                for (int j = 0; j < d; j++)
                    a[row, j] = ((a[row, j] - f * a[rank, j]) % p + p) % p;
            }
            rank++;
        }
        return rank;
    }

    static long ModInverse(long x, long p) => ModPow(((x % p) + p) % p, p - 2, p);   // Fermat, p prime

    static long ModPow(long b, long e, long p)
    {
        long r = 1; b %= p;
        while (e > 0)
        {
            if ((e & 1) == 1) r = r * b % p;
            b = b * b % p;
            e >>= 1;
        }
        return r;
    }

    // ---- shared atoms (the chain hop and the popcount basis, as everywhere in this world) ----
    int[] Hops(int c)
    {
        var moves = new List<int>();
        foreach (var (a, b) in bonds)
            if (((c >> a) & 1) != ((c >> b) & 1)) moves.Add(c ^ (1 << a) ^ (1 << b));
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
