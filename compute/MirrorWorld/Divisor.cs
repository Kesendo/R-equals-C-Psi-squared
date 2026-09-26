using System.Numerics;

namespace MirrorWorld;

// The frozen divisor (adopted 2026-07-25 from docs/proofs/PROOF_R90_FROZEN_DIVISOR.md, registry F140).
//
// Mirror holds the BETWEEN-block folds; Seed holds the WITHIN-block self-duality they leave untouched.
// This holds a third thing the same mirror leaves behind, and it is neither a fold nor a seed: a VALUE
// the Hamiltonian cannot move. On the R90 locus (every reflection pair of site rates carrying the same
// total, gamma_l + gamma_{R(l)} = 2*gbar) the corner block carries lambda = -4*gbar with multiplicity
// AT LEAST floor(N/2) for EVERY coupling J, and, when gbar != 0, exactly floor(N/2) away from finitely
// many couplings. One frozen mode per balanced pair. Three distinctions are visible from here:
// at J = 0 a generic locus profile doubles the count, at real exceptional couplings the root goes defective without
// the rank moving, and at gbar = 0 the diagonal cells stop charging Rooms() and start paying it, so
// the proved kernel lower bound is N at every J; equality of the geometric count was
// exact-checked at selected nonzero J for N=3..7
// (at J = 0 it is N + 2*floor(N/2) for a generic locus profile; that stratum has defective couplings of its own,
// with a size-3 Jordan block at N = 3, which is the proof's Section 9).
//
// THE PARENT IS THE MIRROR, not the frame. (Marginal got to a non-frame parent first, on 2026-07-12,
// hanging on the running cloud it reads; this is the second, and the first among the closed-form
// objects.) The choice is the content: what does the freezing is a ROOM SHORTAGE, and the shortage is
// not the divisor's. The cell mirror
//
//   tauQ : (a,b) -> (R(b), R(a))        (Mirror's transpose leg t, dressed with the site reversal R)
//
// fixes exactly the cells (a, R(a)) -- the coherence between a site and its own mirror partner, one per
// site -- and an involution's even rooms outnumber its odd rooms by exactly its fixed count (every
// 2-cycle contributes one of each and cancels; only a room that is its OWN mirror image stands alone).
// So dim O+ - dim O- = 2*floor(N/2), which the proof doc states in its Section 1 as a dimension count
// "by inspection". An operator ODD under that mirror must send the bigger half into the smaller, and
// what will not fit must freeze. The divisor's OWN possessions are therefore only what it adds on top:
// the FROZEN count (the subtraction that survives the diagonal cells' even defect), the root those
// rooms sit at, and the ladder that says how long each frozen mode holds. Both counts that go into the
// subtraction -- the surplus and the tax -- are tauQ fixed-counts, so both are the mirror's.
//
// INTEGER MATRIX, NOT A FLOATING-POINT RANK. Everything here is a count, an entry-wise residual,
// or a rank over GF(p) (Seed's genre: a rank, never an eigensolver). Two-prime modular nullity is
// an upper bound on rational nullity when both primes are bad; equality with a proved room lower
// bound certifies the rational count at that input. The rank matters most: the frozen root is an exact
// eigenvalue at EVERY coupling, and a floating-point rank silently miscounts it once the coupling is
// small and the chain long, where the other eigenvalues can crowd the root no earlier than
// order J^(2d). So the
// inputs are integers over one common denominator, the block is scaled to Gaussian integers, and the
// rank is taken mod two primes p = 1 (mod 4), where i is a genuine square root of -1.
//
// WHAT STAYS OUTSIDE. The proof's analytic half (the cofactor determinant, the DERIVATION of the two
// boundary clocks, the Schur complements, the defective exceptional couplings) is eigen-work and stays
// in the main repo; the boundary-clock MODULUS is adopted here as a number, the way Mirror adopts its
// price.
public sealed class Divisor : GameObject
{
    public int N { get; }

    readonly long jNum;                 // J = jNum / den
    readonly long[] gNum;               // gamma_l = gNum[l] / den, on the R90 locus
    readonly long den;
    readonly (int a, int b)[] bonds;
    readonly bool zz;

    // the locus profile as integers over the common denominator: any half-profile, mirrored
    // antisymmetrically about the mean, so that every reflection pair sums to 2*gbar exactly.
    public static long[] Locus(int n, long gbarNum, params long[] halfNum)
    {
        var g = new long[n];
        for (int i = 0; i < n; i++) g[i] = gbarNum;
        for (int i = 0; i < n / 2 && i < halfNum.Length; i++)
        {
            g[i] = checked(gbarNum + halfNum[i]);
            g[n - 1 - i] = checked(gbarNum - halfNum[i]);
        }
        return g;
    }

    public Divisor(Mirror mirror, int n, long jNum, long[] gammaNum, long den, bool zz = false)
        : base(mirror)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), "the corner block needs at least two sites");
        ArgumentNullException.ThrowIfNull(gammaNum);
        if (gammaNum.Length != n)
            throw new ArgumentException("the rate profile needs one value per site", nameof(gammaNum));
        if (den <= 0)
            throw new ArgumentOutOfRangeException(nameof(den), "the common denominator must be positive");
        N = n;
        this.jNum = jNum;
        gNum = (long[])gammaNum.Clone();
        this.den = den;
        bonds = Topology.Chain(n);
        this.zz = zz;
    }

    public double J => (double)jNum / den;
    BigInteger SigmaNumerator => gNum.Aggregate(BigInteger.Zero, static (sum, gamma) => sum + gamma);
    public bool IsOnLocus
    {
        get
        {
            BigInteger twiceSigma = 2 * SigmaNumerator;
            for (int site = 0; site < N; site++)
                if ((BigInteger)N * ((BigInteger)gNum[site] + gNum[N - 1 - site]) != twiceSigma)
                    return false;
            return true;
        }
    }

    void RequireLocus()
    {
        if (!IsOnLocus)
            throw new InvalidOperationException("the frozen room and ladder claims require the R90 rate locus");
    }
    public double Sigma => (double)SigmaNumerator / den;
    public double GBar => Sigma / N;
    public double Root => -4.0 * GBar;                      // the frozen value
    public double FoldRoot => 4.0 * GBar - 2.0 * Sigma;     // its GammaFold partner, r -> -r - 2*sigma

    // left: what the divisor produces itself. NOT the rooms and NOT the tax -- both are tauQ counts,
    // and tauQ is the mirror's. What is the divisor's is the subtraction and what sits in the result.
    public override IReadOnlyList<string> Own => new[] { "frozen", "root", "ladder" };

    // ---- the room count: pure counting, no matrix at all ----
    public (int Fixed, int DimOPlus, int DimOMinus, int Surplus, int Tax, int Frozen) Rooms()
    {
        RequireLocus();
        int fixedCells = 0, o = 0;
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
            {
                if (a == b) continue;
                o++;
                if (b == N - 1 - a) fixedCells++;
            }
        int dimOPlus = (o - fixedCells) / 2 + fixedCells;
        int dimOMinus = (o - fixedCells) / 2;
        int surplus = dimOPlus - dimOMinus;

        // The diagonal cells decide the rest, and WHICH WAY they decide it is gbar's doing.
        // With the even defect 4*gbar*P_D present they are even where everything else is odd, so
        // D- cuts floor(N/2) conditions off the surplus: that is the tax, and the count is one
        // frozen mode per balanced PAIR. At gbar = 0 there is no defect to be even, the whole
        // block is odd, and D pays a surplus of its own instead (ceil(N/2) - floor(N/2), the
        // odd-N centre's room), so the count is N: one frozen mode per SITE. Note it is NOT
        // 2*floor(N/2): at odd N the centre room is what makes up the difference. The locus
        // reaches gbar = 0 because it never asks the site rates to be positive.
        // Tax is what the diagonal cells take OUT of the surplus, so frozen = surplus - tax on
        // both strata; on the untaxed one they put in instead, and the tax is MINUS what they pay,
        // -(ceil(N/2) - floor(N/2)). That is the odd-N centre's own room, so it is -1 at odd N and
        // 0 at even N, where there is no centre cell to pay it. One caveat on the names: Surplus
        // stays the O-side count 2*floor(N/2) on both strata, so on the untaxed one it is no
        // longer the whole index that produces Frozen; the missing piece is exactly this Tax.
        bool taxed = SigmaNumerator != 0;
        int tax = taxed ? N / 2 : -((N + 1) / 2 - N / 2);
        int frozen = surplus - tax;
        return (fixedCells, dimOPlus, dimOMinus, surplus, tax, frozen);
    }

    // ---- the hypothesis, entry by entry: is the generator odd under the cell mirror? ----
    // tauQ K tauQ = -K always (h real symmetric and R-invariant); tauQ (2Gamma - 4 gbar) tauQ =
    // -(2Gamma - 4 gbar) on O EXACTLY on the locus. Off the locus the SECOND one is what breaks, and
    // with it the whole freezing. Both worst-case residuals, no eigensolver.
    public (double Hop, double Rate) OddnessResidual()
    {
        int n2 = N * N;
        var k = BuildK();
        var tq = TauQ();
        double hop = 0;
        for (int i = 0; i < n2; i++)
            for (int j = 0; j < n2; j++)
                hop = Math.Max(hop, (k[tq[i], tq[j]] + k[i, j]).Magnitude);

        double rate = 0, gbar = GBar;
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
            {
                if (a == b) continue;                        // the defect lives on D; that is the tax
                int ra = N - 1 - a, rb = N - 1 - b;
                double here = 2.0 * (double)((BigInteger)gNum[a] + gNum[b]) / den - 4 * gbar;
                double there = 2.0 * (double)((BigInteger)gNum[rb] + gNum[ra]) / den - 4 * gbar;
                rate = Math.Max(rate, Math.Abs(there + here));
            }
        return (hop, rate);
    }

    // ---- the mechanism on the untaxed stratum, entry by entry ----
    // OddnessResidual above deliberately skips the diagonal cells, because on the taxed stratum
    // they are EXACTLY where the oddness fails: that failure is the tax. So it cannot tell the
    // two strata apart. This one reads the WHOLE recentered block, D included, which is odd only
    // when the even defect 4*gbar*P_D is absent. Worst-case residual, no eigensolver.
    public double WholeBlockOddnessResidual()
    {
        int n2 = N * N;
        var k = BuildK();
        var tq = TauQ();
        double gbar = GBar, worst = 0;
        double Entry(int r, int c)
        {
            // M~ = J*K - 2*Gamma + 4*gbar on the cell basis; only the diagonal carries the rest.
            double v = 0;
            if (r == c)
            {
                int a = r / N, b = r % N;
                v = 4 * gbar - (a == b ? 0 : 2.0 * (double)((BigInteger)gNum[a] + gNum[b]) / den);
            }
            return v;
        }
        for (int i = 0; i < n2; i++)
            for (int j = 0; j < n2; j++)
            {
                Complex here = J * k[i, j] + new Complex(Entry(i, j), 0);
                Complex there = J * k[tq[i], tq[j]] + new Complex(Entry(tq[i], tq[j]), 0);
                worst = Math.Max(worst, (there + here).Magnitude);
            }
        return worst;
    }

    // ---- a one-sided multiplicity reading, by modular rank ----
    // The block is scaled by N*den into Gaussian integers and ranked mod two primes = 1 (mod 4).
    // Modular nullity can exceed rational nullity at bad primes; matching the proved room lower
    // bound confirms the rational count for an ON-LOCUS input. The matrix is built before reduction with
    // unbounded integers, so a large but valid long profile cannot silently become another matrix.
    public int KernelDimensionUpperBound()
    {
        int n2 = N * N;
        var (re, im) = BuildScaledBlock();
        int rank = 0;
        foreach (long p in ModP.Primes) rank = Math.Max(rank, RankModP(re, im, n2, p));
        return n2 - rank;
    }

    // ---- the four corners, and which root each one carries ----
    // Each one-sided X^N bridge (GammaFold, already in this world) sends a rate r to -r - 2*sigma and
    // flips one popcount index to N - it. So the four blocks with p, q in {1, N-1} carry the root
    // chosen by how many folds separate them from (1,1): an even number returns to -4*gbar, an odd one
    // lands on 4*gbar - 2*sigma. That is what this object owns, and it holds on both chains.
    //
    // What it does NOT say, and used to: that no other block carries. That confinement is the proof
    // document's census (its Section 5) and it is a HEISENBERG statement -- what confines the divisor
    // is the ZZ term as a QUARTIC term (it removes the ladder that carries the corner up the diagonal),
    // not the diagonal it puts on h. This object's default is zz = false, the XY chain, where the same root is
    // carried at the same multiplicity by many more blocks. So the census is not adopted here at all;
    // only the fold parity is, which is the part that travels (found 2026-07-25, gate G2c).
    public (int P, int Q, int Folds, double Root)[] Corners()
    {
        RequireLocus();
        var outp = new List<(int, int, int, double)>();
        foreach (int p in new[] { 1, N - 1 })
            foreach (int q in new[] { 1, N - 1 })
            {
                int folds = (p == N - 1 ? 1 : 0) + (q == N - 1 ? 1 : 0);
                outp.Add((p, q, folds, folds % 2 == 0 ? Root : FoldRoot));
            }
        return outp.ToArray();
    }

    // ---- the ladder: how long each frozen mode holds ----
    // The pair through site c (1-based) and its mirror R(c) sit d_c = N + 1 - 2c apart. The proof
    // bounds its departure order BELOW by 2d_c, and the total valuation BELOW by 2*sum_c d_c.
    // Equality depends on a nonvanishing still open at all N; it was exact-checked at N=3..8.
    // On the zero-mean stratum the room minimum is N, versus floor(N/2) pair entries;
    // these entries do not assign one departing mode to each pair there.
    public (int[] Distances, int[] PerPairLowerBounds, long TotalLowerBound) Ladder()
    {
        RequireLocus();
        int m = N / 2;
        var d = new int[m];
        var per = new int[m];
        for (int c = 1; c <= m; c++)
        {
            d[c - 1] = checked((int)((long)N + 1 - 2L * c));
            per[c - 1] = 2 * d[c - 1];
        }
        return (d, per, 2L * ((long)N * N / 4));
    }

    // the boundary clock the chain carries, adopted as a number (the derivation stays outside):
    // N + 1 for the hopping-only chain, N for the chain that also carries the ZZ diagonal.
    public int ClockModulus => zz ? N : N + 1;

    // ================= the atoms =================

    int[] TauQ()
    {
        var t = new int[N * N];
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
                t[a * N + b] = (N - 1 - b) * N + (N - 1 - a);
        return t;
    }

    long[,] H()
    {
        var h = new long[N, N];
        // XX+YY gives a two-unit single-excitation hop when the ZZ bond is present.
        // The XY normalization used by this object's default is one unit.
        long hop = zz ? 2 : 1;
        foreach (var (a, b) in bonds) { h[a, b] = hop; h[b, a] = hop; }
        if (zz)
            for (int a = 0; a < N; a++)
            {
                long dsum = 0;
                foreach (var (x, y) in bonds) dsum += (a == x || a == y) ? -1 : 1;
                h[a, a] = dsum;
            }
        return h;
    }

    // K at J = 1, the same handshake Mirror reads, restricted to the (1,1) block (floating point, used
    // only for the oddness residual, which is a comparison and not a rank).
    Complex[,] BuildK()
    {
        int n2 = N * N;
        var h = H();
        var k = new Complex[n2, n2];
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
            {
                int col = a * N + b;
                for (int c = 0; c < N; c++)
                {
                    if (h[a, c] != 0) k[c * N + b, col] += -Complex.ImaginaryOne * h[a, c];
                    if (h[c, b] != 0) k[a * N + c, col] += Complex.ImaginaryOne * h[c, b];
                }
            }
        return k;
    }

    // (N*den) * (L_block(J) + 4*gbar), as Gaussian integers: the hop part is purely imaginary
    // (N*jNum times the integer h), the rate part purely real.
    (BigInteger[,] Re, BigInteger[,] Im) BuildScaledBlock()
    {
        int n2 = N * N;
        var h = H();
        var re = new BigInteger[n2, n2];
        var im = new BigInteger[n2, n2];
        BigInteger sigmaNum = SigmaNumerator;                 // 4*gbar * (N*den) = 4*sigmaNum
        BigInteger hopScale = (BigInteger)N * jNum;
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
            {
                int col = a * N + b;
                for (int c = 0; c < N; c++)
                {
                    if (h[a, c] != 0) im[c * N + b, col] -= hopScale * h[a, c];
                    if (h[c, b] != 0) im[a * N + c, col] += hopScale * h[c, b];
                }
                re[col, col] += 4 * sigmaNum;
                if (a != b) re[col, col] -= 2 * (BigInteger)N * ((BigInteger)gNum[a] + gNum[b]);
            }
        return (re, im);
    }

    // rank over F_p with i realized as a square root of -1 (p = 1 mod 4): a + b*i  ->  a + b*r.
    // The embedding is Divisor's own step; the root and the elimination are ModP's, and the world's
    // one prime list is 1 mod 4 precisely so this map exists at both of its primes.
    static int RankModP(BigInteger[,] re, BigInteger[,] im, int d, long p)
    {
        long r = ModP.SqrtMinusOne(p);
        var rows = new long[d][];
        for (int i = 0; i < d; i++)
        {
            rows[i] = new long[d];
            for (int j = 0; j < d; j++)
            {
                long real = ModP.Mod((long)(re[i, j] % p), p);
                long imaginary = ModP.Mod((long)(im[i, j] % p), p);
                rows[i][j] = ModP.Mod(real + ModP.MulMod(r, imaginary, p), p);
            }
        }
        return ModP.Rank(rows, p);
    }
}
