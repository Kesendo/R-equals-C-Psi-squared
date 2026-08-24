namespace MirrorWorld;

// The blind seat (adopted 2026-08-24, the seat-cut blindness arc: experiments/THE_SEAT_THAT_CUTS.md +
// experiments/THE_BLIND_SITE.md; no F-number yet, adopted on the genre gate the way GammaFold's per-site
// turns and Cyclotomy were). Put the watching on ONE seat of a chain and count what it cannot touch: the
// blind dimension is N minus the number of single-excitation eigenspaces the seat's unit vector MEETS,
// the dimension of the largest H-invariant subspace carrying no amplitude at the seat, held here in the
// world's genre as a COUNT with no eigensolver anywhere,
//
//     blind(seat) = N - rank of the Krylov matrix [e, He, H^2 e, ...] on the seat,
//
// an exact GF(p) rank on integer inputs at two primes (a rank mod p can only DROP at a bad prime, so the
// max over two bounds it from the one side that matters, and the reported blind dimension can only ever
// be too large, never too small -- Seed's convention, restated per the world's per-object
// self-containment). Since 2026-08-24 the count equals the main repo's gcd criterion by a Cramer theorem
// with no hypothesis beyond real symmetry; the world adopts the COUNT and leaves the gcd phrasing, the
// blind-projector corner and every path object outside (genre, not topic). On the UNIFORM chain the count
// closes to integer arithmetic, one form per book:
//
//     ZZ on  (Heisenberg):  blind(seat) = (gcd(2*seat+1, N) - 1) / 2
//     ZZ off (XY):          blind(seat) = gcd(seat+1, N+1) - 1
//
// and on the XY book a third, PARITY-forced kind exists beside the mirror-forced and the met: a
// zero-diagonal Jacobi block of odd size is singular outright, so every odd seat of every odd chain is
// blind at every zero-free profile, however irregular. The stationary manifold of the watched sector is
// one larger, span = 1 + blind, on the zero-free chain; the identity can break only where the spectrum
// is degenerate (simplicity is SUFFICIENT for it and not necessary), and the zero bond is NOT the
// discriminator: [1,1,0,1,1] breaks it at seats 1 and 4 while the zero-bond [1,0,1] holds at every seat.
//
// Words, fenced at the door: a SEAT here is the site the single watching sits on, not Seed's unmirrorable
// reflection seat; BLIND keeps the world's own idiom (the watching is blind to these states); the seat
// leaves two principal submatrices behind and this file never calls that a cut, because the world's cuts
// (memory, knower's, block, renewal) are all licences not to compute; and the main repo's "divisor law"
// name stays outside, Divisor being F140's frozen object with an unrelated meaning here.
public sealed class BlindSeat : GameObject
{
    static readonly long[] Primes = { 2147483647L, 999999937L };

    public int N { get; }
    readonly long[] bonds;   // integer couplings on the chain's N-1 bonds; zero allowed (and fenced where it matters)
    readonly bool zz;        // the book: true = Heisenberg (ZZ diagonal on), false = XY (zero diagonal)

    // parent = the frame. The single-excitation sector is the Cone's memory cut REBUILT here in integer
    // arithmetic rather than borrowed (the world's per-object self-containment), so no mechanism of
    // another object is consumed and Divisor's parent rule sends this one to the frame. The count still
    // PREDICTS what a per-site-lit Cone does -- on the XY book, the only book the Cone runs, and the
    // guard pins that prediction dynamically -- but a prediction is not an inheritance.
    public BlindSeat(World world, int n, long[] bondCouplings, bool heisenberg = true) : base(world)
    {
        N = n;
        if (bondCouplings.Length != N - 1)
            throw new ArgumentException($"a chain on {N} sites carries {N - 1} bonds.");
        bonds = bondCouplings;
        zz = heisenberg;
    }

    // The site-indexed integer single-excitation H, in the main repo's integer convention
    // (simulations/seat_cut_blindness.py se_hamiltonian_int): hop 2J on a bond, and on the ZZ book a
    // diagonal that pays -J at a bond's own two ends and +J elsewhere, bond by bond.
    long[,] H()
    {
        var h = new long[N, N];
        for (int b = 0; b < N - 1; b++)
        {
            h[b, b + 1] += 2 * bonds[b];
            h[b + 1, b] += 2 * bonds[b];
        }
        if (zz)
            for (int s = 0; s < N; s++)
                for (int b = 0; b < N - 1; b++)
                    h[s, s] += (s == b || s == b + 1) ? -bonds[b] : bonds[b];
        return h;
    }

    // blind(seat) = N - rank of the Krylov matrix the seat generates, exact over GF(p) at two primes.
    // The Krylov columns are taken mod p from the start (reduction commutes with the matrix product), so
    // nothing overflows at any N; the rank mod p can only drop, so the max over the primes pins it and the
    // reported blind dimension can only ever be too large, never too small.
    public int Blind(int seat)
    {
        var h = H();
        int rank = 0;
        foreach (long p in Primes)
            rank = Math.Max(rank, KrylovRankModP(h, seat, p));
        return N - rank;
    }

    int KrylovRankModP(long[,] h, int seat, long p)
    {
        var cols = new List<long[]>();
        var vec = new long[N];
        vec[seat] = 1;
        for (int k = 0; k <= N; k++)
        {
            cols.Add((long[])vec.Clone());
            var next = new long[N];
            for (int a = 0; a < N; a++)
            {
                long s = 0;
                for (int b = 0; b < N; b++)
                    s = (s + Mod(h[a, b], p) * vec[b]) % p;
                next[a] = s;
            }
            vec = next;
        }
        return RankModP(cols, p);
    }

    // The stationary span of the watched sector: dim of { X : [H, X] = 0 and X carries nothing on the
    // seat's off-diagonal row and column }, exact over GF(p) at two primes; the nullity mod p can only be
    // too LARGE, so the min over the primes errs the same way Blind does. The complex dimension equals
    // this rational nullity because the constraint matrix is integer and rank is field-invariant in
    // characteristic 0. On the zero-free chain the span is 1 + blind (the +1 is the direction that does
    // sit on the seat); the identity can break only at a degenerate spectrum, and the zero bond does not
    // decide it, which is why span and count live as separate readings here. Unlike the count, the span
    // is an N^2-sized elimination: it is read at small N and does not walk past the wall.
    public int Span(int seat)
    {
        var h = H();
        int dim = int.MaxValue;
        foreach (long p in Primes)
            dim = Math.Min(dim, SpanModP(h, seat, p));
        return dim;
    }

    int SpanModP(long[,] h, int seat, long p)
    {
        // unknowns: X[a,b] with a,b both != seat, plus the surviving X[seat,seat].
        var index = new int[N, N];
        for (int a = 0; a < N; a++) for (int b = 0; b < N; b++) index[a, b] = -1;
        int unknowns = 0;
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
                if ((a == seat) == (b == seat)) index[a, b] = unknowns++;
        var rows = new List<long[]>();
        for (int a = 0; a < N; a++)
            for (int b = 0; b < N; b++)
            {
                var row = new long[unknowns];
                bool any = false;
                for (int c = 0; c < N; c++)
                {
                    if (index[c, b] >= 0 && h[a, c] != 0)
                    { row[index[c, b]] = (row[index[c, b]] + Mod(h[a, c], p)) % p; any = true; }
                    if (index[a, c] >= 0 && h[c, b] != 0)
                    { row[index[a, c]] = (row[index[a, c]] - Mod(h[c, b], p) % p + p) % p; any = true; }
                }
                if (any) rows.Add(row);
            }
        return unknowns - RankModP(rows, p);
    }

    // The closed form on the UNIFORM chain, one per book: integer arithmetic over the public Cyclotomy gcd.
    public int UniformLaw(int seat) => zz
        ? (Cyclotomy.Gcd(2 * seat + 1, N) - 1) / 2
        : Cyclotomy.Gcd(seat + 1, N + 1) - 1;

    // The parity forcing, XY book only: both leftover blocks have odd size exactly at an odd seat of an
    // odd chain, a zero-diagonal Jacobi block of odd size is singular (det T_m = -b^2 det T_{m-2} with b
    // the block's own last bond and det T_1 = 0), so the two blocks share the root 0 at every zero-free
    // profile and the seat is blind, however irregular the couplings.
    public bool ParityForced(int seat) => !zz && N % 2 == 1 && seat % 2 == 1;

    // --- the GF(p) atoms, restated per the world's per-object self-containment ---

    static long Mod(long x, long p) { long r = x % p; return r < 0 ? r + p : r; }

    static int RankModP(List<long[]> rows, long p)
    {
        int cols = rows.Count == 0 ? 0 : rows[0].Length;
        int rank = 0;
        var work = rows.Select(r => (long[])r.Clone()).ToList();
        for (int c = 0; c < cols && rank < work.Count; c++)
        {
            int piv = -1;
            for (int r = rank; r < work.Count; r++)
                if (work[r][c] % p != 0) { piv = r; break; }
            if (piv < 0) continue;
            (work[rank], work[piv]) = (work[piv], work[rank]);
            long inv = ModInverse(Mod(work[rank][c], p), p);
            for (int k = c; k < cols; k++) work[rank][k] = Mod(work[rank][k], p) * inv % p;
            for (int r = 0; r < work.Count; r++)
            {
                if (r == rank) continue;
                long f = Mod(work[r][c], p);
                if (f == 0) continue;
                for (int k = c; k < cols; k++)
                    work[r][k] = (Mod(work[r][k], p) - f * work[rank][k] % p + p * p) % p;
            }
            rank++;
        }
        return rank;
    }

    static long ModInverse(long a, long p) => ModPow(a, p - 2, p);

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

    public override IReadOnlyList<string> Own => new[] { "blind", "law", "span" };
}
