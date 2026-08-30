namespace MirrorWorld;

// The blind seat (adopted 2026-08-24, the seat-cut blindness arc: experiments/THE_SEAT_THAT_CUTS.md +
// experiments/THE_BLIND_SITE.md; F157 since 2026-08-24, adopted on the genre gate the way
// GammaFold's per-site turns and Cyclotomy were). Put the watching on ONE seat of a chain and count what it cannot touch: the
// blind dimension is N minus the number of single-excitation eigenspaces the seat's unit vector MEETS,
// the dimension of the largest H-invariant subspace carrying no amplitude at the seat, held here in the
// world's genre as a COUNT with no eigensolver anywhere,
//
//     blind(seat) = N - rank of the Krylov matrix [e, He, H^2 e, ...] on the seat,
//
// an exact GF(p) rank on integer inputs at two primes (a rank mod p can only DROP at a bad prime, so the
// max over two bounds it from the one side that matters, and the reported blind dimension can only ever
// be too large, never too small -- Seed's convention, restated per the world's per-object
// self-containment). That one-sidedness has two premises and both are measured rather than assumed. It
// is a statement about an EXACT H, which is what MaxCoupling below buys: a wrapped H is a different
// matrix and its mod-p ranks bound nothing about this one. And two primes make a bad reduction a
// CONSTRUCTIBLE coincidence rather than an improbable one, not an impossible one: bonds built to meet
// both moduli report too large (XY [2147483647, 999999937] at N = 3 gives 1 where the count is 0), while
// 15682 random profiles at |J| up to 10^9 disagree nowhere. So the COUNT is scale-free and exactly free
// of the coupling's sign, and this ROUTE to it is scale-free only away from the ranking primes.
// Since 2026-08-24 the count equals the main repo's gcd criterion by a Cramer theorem
// with no hypothesis beyond real symmetry; the world adopts the COUNT and leaves the gcd phrasing, the
// blind-projector corner and every path object outside (genre, not topic).
//
// One seat parts the COUNT from the gloss above it, and the fence is F157's Breaks-for. Where the seat's
// OWN ray is H-invariant (an isolated seat, or its incident bonds detuned to zero, both of which this
// constructor accepts) that ray is dark as well, and the dark states are then a UNION of two subspaces
// rather than one: measured at N = 3, XY, bonds [0,1], seat 0, where |seat> stays pure under the watching,
// a state of the Krylov complement stays pure, and their superposition falls to purity 1/2. The count
// reports the larger of the two, and no single number reports a union. So the exact reading is the
// H-invariant one above, which holds at every seat; "count what it cannot touch" is the physical gloss
// and it is one short exactly there.
//
// On the UNIFORM chain the count closes to integer arithmetic, one form per book:
//
//     ZZ on  (Heisenberg):  blind(seat) = (gcd(2*seat+1, N) - 1) / 2
//     ZZ off (XY):          blind(seat) = gcd(seat+1, N+1) - 1
//
// The ZZ book is the ISOTROPIC point and its law is fine-tuned to it, and the fine-tuning is a ZERO SET
// rather than a single point. Carry an anisotropy Delta on the ZZ term alone. The mirror-forced CENTRE
// seat keeps its blindness at EVERY Delta, and that half is a PROOF rather than a sample: at the
// reflection-fixed seat the two principal submatrices the seat leaves behind are conjugate by the chain
// reflection, so they carry the same characteristic polynomial and share every root, at any Delta and any
// palindromic profile (their resultant is identically zero). Every OTHER seat the gcd law names is blind
// exactly on the zero set of ONE equation in Delta: at N = 9 seat 1 that resultant factors as
// 128*D*(D-1)*(D+1)*(D^2-3), so the seat is blind at Delta = 0, +-1 and +-sqrt(3) and nowhere else, and
// the sqrt(3) row reproduces the isotropic row [0,1,0,0,4,0,0,1,0] exactly. So the MET kind is not an
// accident of the isotropic point, it is an accident of a codimension-one LOCUS that happens to contain
// it; sampling Delta can only ever show the seats are not GENERICALLY blind (they are gone at 1/2, 3 and
// 10^-6 alike, so |Delta| = 1 is an isolated point of that locus and not a limit). The first draft of
// this paragraph said "an accident of the isotropic point" and that was a sample promoted to a rule.
// The two signs agree for F152's reason (at Delta = -1 the object is the signless Laplacian, cospectral
// on a bipartite graph). The two booleans built here are Delta = 1 and Delta = 0, BOTH on the locus;
// nothing between them is built.
//
// and on the XY book a third, PARITY-forced kind exists beside the mirror-forced and the met: a
// zero-diagonal Jacobi block of odd size is singular outright, so every odd seat of every odd chain is
// blind at every zero-free profile, however irregular. The stationary manifold of the watched sector is
// one larger, span = 1 + blind, on the zero-free chain, and since 2026-08-25 that is a commutant identity
// with a sharp criterion rather than a pattern with a fence (PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA,
// Theorem A and Corollary C): dim ker = 1 + dim commutant(H restricted to the seat's KRYLOV COMPLEMENT
// W), which is 1 + blind exactly when H|_W has a simple spectrum. The condition sits on W, and NOT on H
// and NOT on the blind subspace: watch [1,1,0] at its isolated site and the whole spectrum is degenerate,
// {-4, 0, 2, 2}, while H|_W is simple and the identity holds at 4 = 1 + 3. A zero-free chain always
// supplies a simple H|_W, which is why it holds there; [1,1,0,1,1] breaks it at seats 1 and 4 while the
// zero-bond [1,0,1] holds at every seat, so the zero bond was never the discriminator.
//
// Words, fenced at the door: a SEAT here is the site the single watching sits on, not Seed's unmirrorable
// reflection seat; BLIND keeps the world's own idiom (the watching is blind to these states); the seat
// leaves two principal submatrices behind and this file never calls that a cut, because the world's cuts
// (memory, knower's, block, renewal) are all licences not to compute; and the main repo's "divisor law"
// name stays outside, Divisor being F140's frozen object with an unrelated meaning here.
public sealed class BlindSeat : GameObject
{
    static readonly long[] Primes = { 2147483647L, 999999937L };

    // The largest |J| a bond may carry, the sibling witness's guard adopted with its reason
    // (SeatCutBlindnessWitness): H() doubles a coupling for the hop and SUMS N-1 of them on the ZZ
    // diagonal, so past this a chain would wrap int64 silently, and a wrapped H is a DIFFERENT matrix
    // whose mod-p ranks say nothing about this one -- the whole "never too small" argument below is an
    // argument about an exact H and is void without it. Measured, before the guard existed: N = 6 on the
    // ZZ book at |J| = 4378862956477877167 reported blind 0 at seats 1 and 4 where the count is 1.
    public const long MaxCoupling = 1L << 40;

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
        if (N < 2)
            throw new ArgumentException($"a chain needs at least two sites; got {N}.");
        if (bondCouplings.Length != N - 1)
            throw new ArgumentException($"a chain on {N} sites carries {N - 1} bonds.");
        foreach (long b in bondCouplings)
            if (Math.Abs(b) > MaxCoupling)
                throw new ArgumentOutOfRangeException(nameof(bondCouplings),
                    $"|J| = {Math.Abs(b)} exceeds {MaxCoupling}: the ZZ diagonal sums up to N-1 couplings and the hop " +
                    "doubles one, so a larger magnitude wraps int64 silently and the count loses its one-sidedness.");
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
    // nothing overflows in THIS loop at any N; the rank mod p can only drop, so the max over the primes
    // pins it and the reported blind dimension can only ever be too large, never too small. The place
    // that can overflow is H() before it, which is built in full precision, and MaxCoupling is what keeps
    // the sentence above true rather than the loop's own arithmetic.
    public int Blind(int seat)
    {
        Seat(seat);
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
    // this rational nullity because the constraint matrix has integer entries and the rank of a FIXED
    // matrix is unchanged by a field EXTENSION; what can change it is the reduction mod p, which is not an
    // extension, and that is the one residual this min bounds. Theorem A says what the dimension IS at
    // every profile, 1 + dim commutant(H|_W), and Corollary C when it is 1 + blind: exactly when H|_W is
    // simple, which a zero-free chain always supplies (the +1 is the direction that does sit on the seat).
    // Read the two as two READINGS and never as one gate: span and count come off the SAME H at the SAME
    // two primes, so a bad reduction inflates both and the identity survives it unbroken. The proof's own
    // G6 takes both sides over the rationals for exactly that reason, and the guard here is the committed
    // literal beside the identity, not the identity alone. Unlike the count, the span is an N^2-sized
    // elimination: it is read at small N and does not walk past the wall.
    public int Span(int seat)
    {
        Seat(seat);
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
    public int UniformLaw(int seat) => Seat(seat) && zz
        ? (Cyclotomy.Gcd(2 * seat + 1, N) - 1) / 2
        : Cyclotomy.Gcd(seat + 1, N + 1) - 1;

    // The parity forcing, XY book only: both leftover blocks have odd size exactly at an odd seat of an
    // odd chain, a zero-diagonal Jacobi block of odd size is singular (det T_m = -b^2 det T_{m-2} with b
    // the block's own last bond and det T_1 = 0), so the two blocks share the root 0 at every zero-free
    // profile and the seat is blind, however irregular the couplings.
    public bool ParityForced(int seat) => Seat(seat) && !zz && N % 2 == 1 && seat % 2 == 1;

    // Every public reading takes a seat, and a seat off the chain is a caller error rather than a
    // number: UniformLaw would otherwise return -1 through a gcd that takes no absolute value.
    bool Seat(int seat) => seat >= 0 && seat < N
        ? true
        : throw new ArgumentOutOfRangeException(nameof(seat), $"seat {seat} is off a chain of {N} sites.");

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
