namespace MirrorWorld;

// The grey room (adopted 2026-11-06 from docs/THE_GREY_GLASS.md and docs/proofs/PROOF_GREY_GLASS_LADDER_IMAGES.md,
// F171; BlindSeat's sibling on the same primitive). Under equal light, Z dephasing at one rate on every site, the
// dissipator charges a cell |a><b| its Hamming distance k = popcount(a xor b) and moves nothing, so a cell's
// price is its HEIGHT on the letter cube and the dissipator is one number on each height plane. A grey room is
// a subspace the Hamiltonian's commutator keeps on which the dissipator is one number: there the motion comes
// through whole, e^{tL} = e^{mu t} e^{t L_H}, and the room's modes do not move with gamma. The largest grey room
// at height k inside a joint-popcount block (p, q) is the largest ad_H-invariant subspace inside that block's
// plane of Hamming distance k, which control theory calls the unobservable subspace of the pair (ad_H, the
// projection off the plane). Held here in the world's genre as a COUNT with no eigensolver anywhere,
//
//     room(p, q, k) = dim(block (p, q)) - rank of the Krylov rows [e_c, e_c A, e_c A^2, ...]
//                     over every cell c of the block OFF the plane, A = ad_H on the block
//                     (the off-plane unit rows are among them, so the rank is at least dim(block) - dim(plane)),
//
// an exact GF(p) rank on integer inputs at two primes: a rank mod p can only DROP at a bad prime, so the
// reading can only ever be too LARGE, never too small, and the min over the primes is the tighter of two
// UPPER bounds on the rational room (the convention BlindSeat reads from the same primitive, ModP, where
// the blind count errs the same way). The rate gamma never enters: the planes are the dissipator's eigenspaces
// at every gamma > 0, so the room is a property of H and the Hamming geometry alone.
//
// What the world adopts as the LOWER bound is F171's theorem, not a second computation: under equal light the
// eta ladder Phi (ket and bra one excitation up) and the spin ladder S+ (ket up, bra down) of
// PROOF_FROZEN_BAND_SO4 commute with ad_H on the XY chain at any real couplings and with the dissipator, so
// S-^s Phi^r carries the boundary block (k, 0), whole and one-to-one, into every block (p, q) whose Hamming
// window |p - q| <= k <= min(p + q, 2N - p - q), in steps of two, holds k: a grey room of dimension C(N, k) at
// height k in every such block. On the ring the wrap bond's sign follows the particle number (Jordan-Wigner),
// and the image is kept only where the fluxes the two sides see serve the ladders: on the EVEN ring at every
// even k, and at odd k, with the wrap bond nonzero, only where the image fills a corner, a block with an empty or a
// full side (proved at every N for the keeping, the losses exact at N = 4, 6 and 8 and for (5, 2) at N = 10); on the ODD
// ring, with the wrap bond nonzero, exactly when the image reaches
// its block from a corner by the ladder its parity of k allows, s in {0, k} at even k and r in {0, N - k} at
// odd k (proved for the keeping, the losses measured at N = 3, 5 and 7). In the block (1, 1) of the UNIFORM
// open XY chain the two ladder images are all there is: C(N, 2) at height 2 and the identity at height 0,
// proved at every N. Where the upper bound meets the ladder bound the rational room is certified at that
// input, the way Divisor certifies its count against the mirror's lower bound.
//
// Beyond the ladders, at the height N/2 strictly inside a block's window, the chain's reflection builds grey
// operators more at even N and at couplings the reflection keeps: the products of the filled mirror orbitals
// with the staggered mirror pairs, |Psi_e,o><(Lambda-dagger)^j 0| in the block (N/2, 2j) and their transposes in
// (2j, N/2), 1 <= j <= N/2 - 1, two per block and four in the self-transposed middle block (N/2, N/2) at
// 4 | N, where a product and its transpose land together. Proved grey; what they add to the ladder images is
// measured and not proved (2 at N = 4 and 6, where at N = 4 the middle block's transposes fall into the
// ladder span, 2 in (4, 2) and 4 in (4, 4) at N = 8), so Products reports the construction's count beside the
// ladder bound and is never added to it, and the products' own ladder images in other blocks are not
// counted. The ring's bond-axis and antipodal products stay outside, with the counts to N = 12, on the page.
//
// Words, fenced at the door: a ROOM here is a subspace the motion never leaves, not Divisor's even and odd
// rooms; HEIGHT is the letter cube's k_Z, the Hamming distance of a cell; a CORNER is a block with an empty or
// a full side, NOT Divisor's corner block (1, 1), which is the block both objects read; the LADDERS are
// PROOF_FROZEN_BAND_SO4's and are not built here, only their theorem is adopted, so Own says "images" where
// Divisor's "ladder" is how long a frozen mode holds and CollisionGap's the comb rung. Adjacent inside the
// world and not docked: on the Heisenberg book the (1, 1) room is the frozen modes of the uniform point, the
// count Divisor reads in its corner block on the R90 locus (PROOF_R90_FROZEN_DIVISOR, the same source).
// The elimination below is this file's own and not ModP.Rank's: the closure must know when the Krylov space
// stops growing, which an incremental insert says and a batch rank cannot; the same routine stands in
// RCPsiSquared.Diagnostics' FreeFermionFloorWitness (KrylovRowRankModP), which reads the row k = 1 of the open
// chain's blocks (n, n + 1), the k = 1 case of this count, certified there against its own lower bound.
// The Heisenberg book is carried as a reading for the block (1, 1) alone (its ZZ term loses the spin ladder,
// and its room there is the floor(N/2) frozen modes on the open chain, measured at N = 3 to 8): the ladder
// bound is an XY theorem and reports 0 on that book.
public sealed class GreyRoom : GameObject
{
    // BlindSeat's guard, adopted with its reason: the hop doubles a coupling and the ZZ diagonal sums N of them,
    // so past this an H wraps int64 silently, and a wrapped H is a different matrix whose mod-p ranks say
    // nothing about this one.
    public const long MaxCoupling = 1L << 40;

    // The rows are dense, C(N, p)·C(N, q) longs each, and the middle block at N = 10 holds 63 504 cells: the
    // closure would want tens of gigabytes per prime. The page's own readings end at N = 8 ((4, 2) and its
    // images), and that is the wall this object stands at.
    public const int MaxSites = 8;

    public int N { get; }
    public bool Ring { get; }
    public bool Heisenberg { get; }
    readonly long[] bonds;   // N - 1 couplings on the chain, N on the ring (the last one the wrap bond)

    // parent = the frame: the block and its planes are rebuilt here in integer arithmetic (the joint-popcount
    // block is Block's object, the height LetterCube's; neither mechanism is consumed), so Divisor's parent
    // rule sends this one to the frame, as it sent BlindSeat.
    public GreyRoom(World world, int n, long[] bondCouplings, bool ring = false, bool heisenberg = false) : base(world)
    {
        N = n;
        if (N < 2 || (ring && N < 3))
            throw new ArgumentException($"a {(ring ? "ring needs at least three" : "chain needs at least two")} sites; got {N}.");
        if (N > MaxSites)
            throw new ArgumentOutOfRangeException(nameof(n), $"N = {N} exceeds {MaxSites}: the dense closure would exhaust memory.");
        int expected = ring ? N : N - 1;
        if (bondCouplings.Length != expected)
            throw new ArgumentException($"a {(ring ? "ring" : "chain")} on {N} sites carries {expected} bonds.");
        foreach (long b in bondCouplings)
            if (Math.Abs(b) > MaxCoupling)
                throw new ArgumentOutOfRangeException(nameof(bondCouplings),
                    $"|J| = {Math.Abs(b)} exceeds {MaxCoupling}: a larger magnitude wraps int64 silently.");
        bonds = (long[])bondCouplings.Clone();
        Ring = ring;
        Heisenberg = heisenberg;
    }

    // ---- the Hamming window and the ladders' theorem (F171), adopted ----

    public static int Choose(int n, int k) => k < 0 || k > n ? 0 : checked((int)Block.Binomial(n, k));

    // The Hamming distances a block (p, q) holds: |p - q| <= k <= min(p + q, 2N - p - q) in steps of two.
    public IReadOnlyList<int> Window(int p, int q)
    {
        CheckBlock(p, q);
        var w = new List<int>();
        int top = Math.Min(p + q, 2 * N - p - q);
        for (int k = Math.Abs(p - q); k <= top; k += 2) w.Add(k);
        return w;
    }

    public bool IsCorner(int p, int q) => CheckBlock(p, q) && (Math.Min(p, q) == 0 || Math.Max(p, q) == N);

    // Whether F171 proves the ladder image of (k, 0) in the block (p, q) kept by the Hamiltonian part: always
    // on the chain; on the even ring at even k, and at odd k only in a corner; on the odd ring exactly from a
    // corner by the ladder its parity of k allows; a ring whose wrap bond is zero is the open chain and keeps
    // them all. False means NOT PROVED KEPT: a loss where F171 measures one (the even ring to N = 8 and (5, 2)
    // at N = 10, the odd ring to N = 7), beyond that an open row. The XY theorem; on the Heisenberg book the
    // spin ladder is lost and nothing is claimed.
    public bool Kept(int p, int q, int k)
    {
        if (!Window(p, q).Contains(k) || Heisenberg) return false;
        if (!Ring || bonds[N - 1] == 0) return true;
        int s = (k - p + q) / 2, r = q - s;     // steps down the spin ladder and up the eta ladder
        if (N % 2 == 0) return k % 2 == 0 || IsCorner(p, q);
        return k % 2 == 0 ? (s == 0 || s == k) : (r == 0 || r == N - k);
    }

    // The proved lower bound on the room at height k in (p, q): C(N, k) where the ladder image is kept.
    public int LadderBound(int p, int q, int k) => Kept(p, q, k) ? Choose(N, k) : 0;

    // The mirror products at the height N/2 on the open chain at even N and couplings the reflection keeps:
    // two in each block (N/2, 2j) and (2j, N/2), 1 <= j <= N/2 - 1, four in the self-transposed middle block
    // (N/2, N/2) at 4 | N; proved grey, what they add to the ladder images measured and never added to the
    // bound. Zero elsewhere, and zero off the window.
    public int Products(int p, int q, int k)
    {
        if (!Window(p, q).Contains(k)) return 0;
        if (Ring || Heisenberg || N % 2 != 0 || k != N / 2 || !ReflectionKeepsCouplings()) return 0;
        int m = N / 2;
        bool ketSide = p == m && q % 2 == 0 && q >= 2 && q <= N - 2;
        bool braSide = q == m && p % 2 == 0 && p >= 2 && p <= N - 2;
        return (ketSide ? 2 : 0) + (braSide ? 2 : 0);
    }

    public bool ReflectionKeepsCouplings()
    {
        for (int b = 0; b < N - 1; b++)
            if (bonds[b] != bonds[N - 2 - b]) return false;
        return true;
    }

    // ---- the room, read from below the wall ----

    // The mod-p reading, min over the two primes: an upper bound on the rational room, exact where it meets
    // the ladder bound.
    public int Room(int p, int q, int k)
    {
        if (!Window(p, q).Contains(k)) return 0;
        int room = int.MaxValue;
        foreach (long prime in ModP.Primes)
            room = Math.Min(room, RoomModP(p, q, k, prime));
        return room;
    }

    public bool Certified(int p, int q, int k) => Room(p, q, k) == LadderBound(p, q, k);

    // The cells of the block (p, q): every (a, b) with popcount(a) = p and popcount(b) = q.
    (int[] Ket, int[] Bra, Dictionary<(int, int), int> Index) Cells(int p, int q)
    {
        var kets = new List<int>(); var bras = new List<int>();
        var index = new Dictionary<(int, int), int>();
        var withP = Enumerable.Range(0, 1 << N).Where(a => System.Numerics.BitOperations.PopCount((uint)a) == p).ToArray();
        var withQ = Enumerable.Range(0, 1 << N).Where(b => System.Numerics.BitOperations.PopCount((uint)b) == q).ToArray();
        foreach (int a in withP)
            foreach (int b in withQ)
            {
                index[(a, b)] = kets.Count;
                kets.Add(a); bras.Add(b);
            }
        return (kets.ToArray(), bras.ToArray(), index);
    }

    // H on a basis state, in BlindSeat's integer convention: the hop 2J on a bond (XX + YY), and on the
    // Heisenberg book the ZZ diagonal +J where a bond's two sites agree and -J where they differ.
    IEnumerable<(int State, long Coefficient)> HAction(int a)
    {
        int nb = bonds.Length;
        long diag = 0;
        for (int b = 0; b < nb; b++)
        {
            int i = b, j = (b + 1) % N;
            bool ai = ((a >> i) & 1) == 1, aj = ((a >> j) & 1) == 1;
            if (ai != aj) yield return (a ^ (1 << i) ^ (1 << j), 2 * bonds[b]);
            if (Heisenberg) diag += ai == aj ? bonds[b] : -bonds[b];
        }
        if (Heisenberg && diag != 0) yield return (a, diag);
    }

    // ad_H on a cell: H|a><b| - |a><b|H = sum_a' H[a',a] |a'><b| - sum_b' H[b',b] |a><b'|, H real symmetric.
    List<(int Cell, long Coefficient)>[] Commutator(int[] ket, int[] bra, Dictionary<(int, int), int> index)
    {
        var image = new List<(int, long)>[ket.Length];
        for (int c = 0; c < ket.Length; c++)
        {
            var img = new List<(int, long)>();
            foreach (var (a2, coef) in HAction(ket[c])) img.Add((index[(a2, bra[c])], coef));
            foreach (var (b2, coef) in HAction(bra[c])) img.Add((index[(ket[c], b2)], -coef));
            image[c] = img;
        }
        return image;
    }

    int RoomModP(int p, int q, int k, long prime)
    {
        var (ket, bra, index) = Cells(p, q);
        int d = ket.Length;
        var image = Commutator(ket, bra, index);
        // the observable rows: start from every cell off the plane and close under w -> w A
        var basis = new List<(int Pivot, long[] Row)>();
        var frontier = new List<long[]>();
        int plane = 0;
        for (int c = 0; c < d; c++)
        {
            if (System.Numerics.BitOperations.PopCount((uint)(ket[c] ^ bra[c])) == k) { plane++; continue; }
            var e = new long[d]; e[c] = 1;
            if (Insert(basis, e, prime)) frontier.Add(e);
        }
        while (frontier.Count > 0)
        {
            var next = new List<long[]>();
            foreach (var w in frontier)
            {
                // (w A)[c'] = sum over cells c of w[c] * A[c, c'], A[c, c'] the coefficient of c in ad_H(c')
                var wa = new long[d];
                for (int c2 = 0; c2 < d; c2++)
                {
                    long s = 0;
                    foreach (var (c, coef) in image[c2])
                        if (w[c] != 0) s = ModP.Mod(s + ModP.MulMod(w[c], ModP.Mod(coef, prime), prime), prime);
                    wa[c2] = s;
                }
                if (Insert(basis, wa, prime)) next.Add(wa);
            }
            frontier = next;
        }
        // every off-plane unit row is in the basis, so the rank splits as (d - plane) + (the observed part of
        // the plane); the room is what the plane keeps unobserved.
        return d - basis.Count;
    }

    // Reduce v against the basis in insertion order (each row carries zeros at every earlier pivot), then
    // insert its first nonzero as a new pivot, normalised to 1. Returns false when v reduces to zero.
    static bool Insert(List<(int Pivot, long[] Row)> basis, long[] v, long prime)
    {
        int d = v.Length;
        for (int i = 0; i < d; i++) v[i] = ModP.Mod(v[i], prime);
        foreach (var (pivot, row) in basis)
        {
            long f = v[pivot];
            if (f == 0) continue;
            for (int i = 0; i < d; i++)
                if (row[i] != 0) v[i] = ModP.Mod(v[i] - ModP.MulMod(f, row[i], prime), prime);
        }
        int lead = Array.FindIndex(v, x => x != 0);
        if (lead < 0) return false;
        long inv = ModP.ModInverse(v[lead], prime);
        for (int i = 0; i < d; i++) v[i] = ModP.MulMod(v[i], inv, prime);
        basis.Add((lead, v));
        return true;
    }

    bool CheckBlock(int p, int q) => p >= 0 && p <= N && q >= 0 && q <= N
        ? true
        : throw new ArgumentOutOfRangeException($"({p}, {q}) is not a block of {N} sites.");

    public override IReadOnlyList<string> Own => new[] { "room", "images", "kept", "products" };
}
