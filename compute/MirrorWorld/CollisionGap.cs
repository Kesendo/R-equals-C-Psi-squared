using System.Numerics;

namespace MirrorWorld;

// The collision gap (F161), adopted 2026-09-08 from docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md
// + experiments/THE_COMB_ON_THE_ROAD.md + simulations/collision_gap_odd_orders.py: what happens to
// two levels that agree, when the road starts to move under them.
//
// THE LAW. Write a chain level's motion in the wrap bond as E_k(u) = 2cos(theta_k) + sum_m d_m u^m.
// The proof computes five coefficients and puts FOUR of them in multiplier form; those four come home.
// Each is a signed combination of NEIGHBOURING evaluations of ONE comb under an integer multiplier, a
// difference of two through third order and three cosines at fifth. ETA IS NOT IN THEM: Lemma B
// gives d_m = eta^m * D_m with D_m free of eta, which is why the odd orders carry the reflection sign
// and the even ones do not, and the closed forms below are the D_m. A reader who maps
// Prefactor(n, 1) = (2, n) onto d_1 is off by that sign at every odd order. The odd
// orders read
//     X_2j(tau) = sum_{k in tau} (-1)^(k+1) cos(2jk pi / n) = -M_{n+2j}(tau)
// and the even orders read M_{2j+1}(tau) = sum_{k in tau} cos((2j+1)k pi / n). Differencing a colliding
// pair gives the GAP, and which of its orders vanish is decided by a gcd: the multiplier is a Galois
// automorphism of Q(zeta_2n) exactly when gcd(n + 2j, 2n) = 1, which at ODD n reduces to gcd(j, n) = 1.
// A collision is DeltaM_1 = 0, so an automorphism carries it to the whole rung it reaches, and the rung
// dies. At every odd n and for EVERY pair, standing or separating, DeltaX_2 = DeltaX_4 = 0 and hence
// c_3 = 0. The rung j = 0 is never an automorphism, gcd(n, 2n) = n, and that is the rung which leaves
// the first order standing. Since F129 fires only at 3|n or 10|n, and 10|n forces n even, every odd
// firing modulus has 3|n, so j = 3 is the first surviving rung and X_6 first enters at FIFTH order.
//
// At EVEN n the collision never reaches the X ladder by that route: n + 2j is even, so no rung is
// killed this way. What carries the vanishing instead is the SHAPE of the triples, the Conway-Jones
// ROT3 form read on the doubled labels, for which X_2j = 0 is FORCED whenever 3 does not divide j.
// The forcing has exactly one direction: at 3|j the coset collapses to a specific cosine, and n = 24
// with (1, 7, 9) shows it nonzero at j = 3 and zero at j = 6, so a version reading "zero iff 3 does not
// divide j" is false. At an even n parity-uniformity comes free with the shape (6|n makes every coset
// parity-homogeneous); at an odd n it is a real and load-bearing hypothesis, and (1, 2, 4) at n = 9
// breaks it.
//
// THE TWO ROUTES ARE ONE MECHANISM READ AT TWO RESOLUTIONS, which is Corollary G and not an aside:
// Theorem D decomposes the PAIR's vanishing sum in Z[zeta_2n] and concludes DeltaM_m = 0; the ROT3 rung
// lemma decomposes a SINGLE triple's doubled-label sum in Z[zeta_n] and concludes X_2j(tau) = 0, which
// the collision hypothesis does not supply and the shape does. The argument is the same word for word
// with 2n replaced by n: Theorem D is the case where the multiplier is invertible on every piece at
// once, the rung lemma the case of pieces of ratio-order 3, where "3 does not divide j" is exactly
// "the multiplier is coprime to their ratio-order". They name the same first surviving rung, three.
//
// Exactness, LevelCollision's own convention: the combs are read in GF(p) at two independent primes
// p = 1 (mod 2n) with zeta of exact order 2n, since 2cos(m pi / n) = zeta^m + zeta^-m. The reduction is
// a ring homomorphism, so NONZERO at one prime is exactly nonzero and every SEPARATES verdict here is
// proof grade. THE OTHER DIRECTION IS NOT FREE AND THIS OBJECT DOES NOT PRETEND OTHERWISE: STANDING is
// a vanishing read at both primes, and so is the second order's death. What the theorems give is the
// EXPLANATION of SOME of those zeros, never a decision procedure this file runs, and the first thing
// to say is which zeros are left over. Of the 223 standing pairs of the census, 200 stand at an ODD
// comb, and there the theorems give c_3 = 0 and nothing else: the first order runs on the rung j = 0,
// which no automorphism reaches, so those 200 stand because their two triples happen to carry the SAME
// odd-label count, and nothing here forces that. What the two theorems do explain is the vanishing of
// the LADDER rungs, and the two orders take DIFFERENT splits, which is worth keeping straight because
// they look alike:
//   * the FIRST order lives on the X ladder and splits by the PARITY of n. At odd n the multiplier
//     reaches the ladder and kills the rung; at even n it never does, and where 3|n the shape carries
//     the vanishing OF c_3. The shape does not touch c_1: fifty-eight pairs of this census satisfy the
//     ROT3 hypothesis in full and separate anyway, all with an odd-label difference of three. What
//     makes the twelve even-n non-mirror pairs stand is the shape for c_3 AND a parity match for c_1,
//     two facts, and the proof is emphatic that this is not a technicality.
//   * the SECOND order lives on the M ladder and splits by 3|n. Its rung is m = 3, so gcd(3, 2n) = 1
//     exactly when 3 does not divide n, and there the gcd kills c_2 for EVERY pair. Among the firing
//     moduli that is the 10|n family with 3 not dividing n, smallest member n = 20; at every 3|n the
//     gcd fails, and only Corollary G's LOCAL piece criterion can speak, which this file does not run.
// So the 20 deaths at n = 20 are the gcd's and the 40 at n = 30 are not explained here at all.
// No eigensolver runs and nothing at u != 0 is ever diagonalised. The u = 0 levels ARE computed, in
// LevelCollision, as zeta^k + zeta^-k; they are that object's and they are how the pairs are found.
//
// PARENT the Crack, not the frame. The series expands the road's own polynomial at u = 0, so what this
// object consumes is the road, and the combs arrive through it from the Cyclotomy: both ends of the
// road are that object's, and this one owns neither. It owns the SERIES (the four carried coefficients
// in multiplier form), the LADDER (which comb each order reads, and the gcd that kills a rung) and the
// GAP (what a colliding pair does as the road moves). Marginal is the precedent that READING an object
// is enough to hang on it, and this one reads LevelCollision's census; the edge is the Crack's anyway,
// because what fixes it is the SERIES, which is a single level's motion and exists with no collision in
// it at all. The census is read here and not rebuilt, and the pairs stay LevelCollision's.
//
// The expansion point is u = 0, the CHAIN END. F160's Theorem E lives at the ring end, in delta = 1 - u,
// and is a pair SPLIT at modulus N; nothing transfers between them and this file claims nothing there.
//
// What stays outside. The multipliers' PARITY is a theorem at every order; their RANGE is not, so the
// general shape past the fifth is open and no order above five is offered here. d_4 is computed in the
// proof and deliberately not carried: it has no eta and enters no odd order, but it IS needed to reach
// d_5 and the proof's script keeps it. Corollary G's piece criterion is sharper than the gcd and is
// decomposition-local; deciding it needs the label and orbit engines that F129's own inventory leaves
// in its proof, so what comes home is the LOWER BOUND those pieces force, which is nothing but
// LevelCollision's families C and L. Two things about that bound, because it is easy to credit it with
// more than it has. Inside the n <= 30 census family L never occurs at all, its door being 70, so the
// equality at n = 20 and n = 30 is family C alone and no c_2 reading here exercises L. What L does have
// from below, in this same assembly, is its contribution to a census TOTAL: LevelCollisionTests ties
// CensusOf(70) = 140 against FamilyCount(C, 70) + FamilyCount(L, 70) = 120 + 20. So the arc's note that
// L has no from-below check is about its MEMBERSHIP and its c_2 consequence, and this adoption does not
// narrow that gap. And the converse
// gated at all 2558 pairs of the nine moduli (L13) is a different statement from this bound: it is
// Corollary G's tiling criterion, c_2 = 0 exactly when some minimal tiling is 3-free, which needs the
// decomposition this file does not carry. Whether an even-n standing pair MUST be ROT3, and whether a
// mixed-parity pair can stand, are open.
//
// Words, fenced at the door. A COMB here is F129's mode ladder at modulus n and the letter n is that
// modulus, never the site count; the world's combs proper are the Cyclotomy's two turn-fraction
// families, which arrive as inheritance. A RUNG is the index j on the multiplier ladder: it is NOT
// Seed's rungs, which are an Own output name in this same assembly (n_diff on the (1,2) pencil), and it
// is not the main repo's rate-ladder rung (the 2 gamma rung, F65's k = 1 rung). The cosine indices the
// closed forms are written on are called that and not rungs, because they are 2j on the odd orders and
// 2j + 1 on the even ones and confusing the two silently is the trap this naming exists to close. The
// GAP is the difference of two levels along the road, never a spectral gap (D6_Gap, Q*_gap). ROT3 is
// the Conway-Jones shape on the DOUBLED labels: Seed's TripleFamily.Rot3 is the same shape on the
// narrower zero-sum domain and the two do not coincide here. And THEOREM E is ambiguous between this
// object's proof, where it IS the ROT3 rung lemma and is what the letter means unqualified, and F160's,
// where it is the ring-end split. PROOF_COLLISION_GAP_ODD_ORDERS handles the overlap by qualifying
// Theorem A by file wherever both are in one sentence; this file goes further and calls the F161 one
// the ROT3 rung lemma throughout, so that the bare letter never has to be read at all. Where F160's is
// meant it is written as F160's Theorem E, and where Theorem A is named the file is named with it.
// One more, and it is the sharpest in the assembly: LADDER is an Own output name here AND in Divisor,
// where it is J^(2d), how long a frozen mode holds. Duplicate Own names are this world's practice
// (orbit, cube, vertices, conjugation are each shared by two objects) and no bucket is confused, since
// Divisor is not in this chain; it is named because a fence that lists Seed's rungs and skips the
// identical collision one row down would be claiming to be exhaustive and not be. Rung also carries a
// third in-assembly sense, Block's and Redistribution's disagreement rung -2*gamma*k, and Gap a second,
// OrderSorting's GenericGap.
public sealed class CollisionGap : GameObject
{
    /// <summary>The road this series expands, at its chain end. Only its N is read: the road as an
    /// object is the whole road and its u is a location on it, while the expansion point is this
    /// object's own statement and is always nought, so a parent built at u != 0 is accepted and its u
    /// ignored.</summary>
    public Crack Road { get; }

    /// <summary>The comb modulus n = N_sites + 1, never the site count.</summary>
    public int Ncomb { get; }

    public CollisionGap(Crack road) : base(road)
    {
        Road = road;
        Ncomb = road.N + 1;
        CheckComb(Ncomb);
    }

    // left: the series (the four carried coefficients in multiplier form), the ladder (which comb each
    // order reads, and the gcd that kills a rung), and the gap (what a colliding pair does as the road
    // moves). NOT the road and NOT the combs: the first is the Crack's, the second the Cyclotomy's, and
    // neither is this object's to own. Nor the collisions: those are LevelCollision's census, read here.
    public override IReadOnlyList<string> Own => new[] { "series", "ladder", "gap" };

    static void CheckComb(int n)
    {
        if (n < 5) throw new ArgumentOutOfRangeException(nameof(n), n, "a comb needs n >= 5");
    }

    // ------------------------------------------------------------------ the ladder

    /// <summary>The multiplier an ODD order reads its comb under: X_2j = -M_{n+2j}.</summary>
    public static int OddOrderMultiplier(int n, int j)
    {
        CheckComb(n);
        return n + 2 * j;
    }

    /// <summary>The multiplier an EVEN order reads its comb under: M_{2j+1}.</summary>
    public static int EvenOrderMultiplier(int j) => 2 * j + 1;

    /// <summary>Whether reading the comb under this multiplier is a Galois automorphism of
    /// Q(zeta_2n), which is what carries a collision onto the rung and kills it. The modulus is 2n and
    /// not n: the comb lives at 2n-th roots, since 2cos(m pi / n) = zeta^m + zeta^-m.</summary>
    public static bool IsAutomorphism(int n, int multiplier)
    {
        CheckComb(n);
        return Cyclotomy.Gcd(multiplier, 2 * n) == 1;
    }

    /// <summary>Whether the collision reaches the odd rung j and kills it. False at every even n, where
    /// n + 2j is even and no rung is reachable this way: there the shape carries it.</summary>
    public static bool GaloisKillsOddRung(int n, int j) => IsAutomorphism(n, OddOrderMultiplier(n, j));

    /// <summary>The odd-n reduction: gcd(n + 2j, 2n) = 1 and gcd(j, n) = 1 are the same condition. It
    /// has no content at an even comb, where the left side is false at every j, so this refuses one
    /// rather than returning a true that would mean nothing.
    ///
    /// It is a THEOREM at an odd comb and therefore cannot return false: n + 2j is odd, so
    /// gcd(n + 2j, 2n) = gcd(n + 2j, n) = gcd(2j, n) = gcd(j, n), the last step because n is odd. This
    /// is said rather than gated, because a gate on a predicate no input can falsify is a gate that
    /// cannot fail: a body reading `return true` would pass every assertion about THIS method. What a
    /// test can hold is that the two SIDES are live over the range, both taking both values, which is
    /// what makes the equality a statement rather than a tautology between constants.</summary>
    public static bool GcdIdentityHolds(int n, int j)
    {
        CheckComb(n);
        RequireOddComb(n, "the reduction");
        return GaloisKillsOddRung(n, j) == (Cyclotomy.Gcd(j, n) == 1);
    }

    /// <summary>The first rung j >= 1 that the collision does NOT kill, by the gcd route. Three at every
    /// odd firing modulus, since those all carry 3|n. It refuses an EVEN comb for the same reason
    /// GcdIdentityHolds does: there this route kills nothing at all, so its "first survivor" would be
    /// the first rung asked about rather than a fact. At an even comb with 3|n the shape answers
    /// instead; at an even comb WITHOUT it, n = 20 being the census member, neither route speaks, and
    /// none is needed there since no pair stands. (The rung j = 0 is never killed at any n, which is why
    /// the search starts at one: it is the rung that leaves the first order standing, not a survivor of
    /// anything.)</summary>
    public static int FirstSurvivingOddRungByGcd(int n)
    {
        CheckComb(n);
        RequireOddComb(n, "the first surviving rung by the gcd");
        for (int j = 1; ; j++)
            if (!GaloisKillsOddRung(n, j)) return j;
    }

    /// <summary>The first rung the ROT3 rung lemma does not force to vanish, read off the forcing rather
    /// than stated: three, since the forcing is exactly 3 not dividing j. The same number the gcd route
    /// reaches at an odd firing modulus, by the same argument one resolution finer. The lemma needs 3|n
    /// to have any triple to speak about, so this is the route's answer WHERE IT APPLIES and not a fact
    /// about every comb; at n = 20 there is no ROT3 triple at all.</summary>
    public static int FirstSurvivingRungByShape()
        => Enumerable.Range(1, 16).First(j => !Rot3ForcesVanishing(j));

    /// <summary>The ROT3 rung lemma's forced direction: on a parity-uniform doubled-label ROT3 triple,
    /// X_2j vanishes whenever 3 does not divide j. One direction only; at 3|j the coset collapses to a
    /// cosine which may or may not be zero.</summary>
    public static bool Rot3ForcesVanishing(int j) => j % 3 != 0;

    static void RequireOddComb(int n, string what)
    {
        if (n % 2 == 0)
            throw new ArgumentException($"{what} is the ODD comb's; n = {n} is even, where no multiplier is an automorphism", nameof(n));
    }

    // ------------------------------------------------------------------ the series, in multiplier form

    /// <summary>The cosine indices r the order is written on, cos(r theta), consecutive in every carried
    /// case. These are the COMB subscripts, 2j on the odd orders and 2j + 1 on the even ones, and not the
    /// rung index j. Order 4 is computed in the proof and not carried; past the fifth the shape is
    /// open.</summary>
    public static IReadOnlyList<int> CosineIndices(int order) => order switch
    {
        1 => new[] { 0, 2 },
        2 => new[] { 1, 3 },
        3 => new[] { 2, 4 },
        5 => new[] { 0, 2, 4, 6 },
        4 => throw new ArgumentOutOfRangeException(nameof(order), order,
                 "d_4 is computed in the proof and deliberately not carried: it has no eta and enters no odd order, though it is needed to reach d_5"),
        _ => throw new ArgumentOutOfRangeException(nameof(order), order,
                 "PROOF_COLLISION_GAP_ODD_ORDERS Theorem A reaches order 5; the ladder's general shape past it is open"),
    };

    /// <summary>The order's overall rational prefactor, as the formula writes it and unreduced:
    /// 2/n, (n-3)/n^2, 2(n-2)(n-4)/(3n^3), and F = 4(n-6)(n-2)/(15 n^5) at fifth order.</summary>
    public static (BigInteger Num, BigInteger Den) Prefactor(int n, int order)
    {
        CheckComb(n);
        _ = CosineIndices(order);
        BigInteger b = n;
        return order switch
        {
            1 => (2, b),
            2 => (b - 3, b * b),
            3 => (2 * (b - 2) * (b - 4), 3 * b * b * b),
            _ => (4 * (b - 6) * (b - 2), 15 * BigInteger.Pow(b, 5)),
        };
    }

    /// <summary>The common denominator the order's weights are written over: one through third order,
    /// four at fifth, where the weights are 5(n-1)/4, (3n^2-16n+16)/2 and -3(2n-3)(n-3)/4.</summary>
    public static BigInteger WeightDenominator(int order)
    {
        _ = CosineIndices(order);
        return order == 5 ? 4 : 1;
    }

    /// <summary>The order's weight on cos(r theta), as a numerator over WeightDenominator. It REFUSES a
    /// cosine index the order is not written on rather than returning zero, because the two indices in
    /// play differ by a factor of two and a silent zero would tell a caller who passed the rung j that
    /// the coefficient vanishes, which at j = 3 is the exact opposite of the truth. The constant index
    /// at fifth order is carried and weighs zero: the order spans the ladder from nought without
    /// touching it, and the zero is what the sum-to-zero identity needs.</summary>
    public static BigInteger WeightNumerator(int n, int order, int cosineIndex)
    {
        CheckComb(n);
        var indices = CosineIndices(order);
        if (!indices.Contains(cosineIndex))
            throw new ArgumentOutOfRangeException(nameof(cosineIndex), cosineIndex,
                $"order {order} is written on cos({string.Join(", ", indices)} theta); mind that a rung j reads the cosine index 2j at the odd orders and 2j+1 at the even ones");
        BigInteger b = n;
        if (order != 5) return cosineIndex == indices[0] ? BigInteger.One : BigInteger.MinusOne;
        return cosineIndex switch
        {
            0 => BigInteger.Zero,
            2 => 5 * (b - 1),
            4 => 2 * (3 * b * b - 16 * b + 16),
            _ => -3 * (2 * b - 3) * (b - 3),
        };
    }

    // ------------------------------------------------------------------ the combs, exactly

    /// <summary>Two independent fields the comb can be read in: p = 1 (mod 2n) with zeta of exact order
    /// 2n. Two, for the same reason the level census takes two.</summary>
    public static (long P, long Zeta)[] Settings(int n)
    {
        CheckComb(n);
        var first = ModP.CyclotomicPrime(2 * n, 0);
        var second = ModP.CyclotomicPrime(2 * n, first.P);
        return new[] { first, second };
    }

    /// <summary>2cos(e pi / n) in GF(p), which is zeta^e + zeta^-e once zeta has order 2n. The exponent
    /// is taken mod 2n, so a multiplier may run past the comb without folding by hand.</summary>
    public static long TwoCos(int n, long e, long p, long zeta)
    {
        long order = 2L * n;
        long r = ModP.Mod(e, order);
        return (ModP.ModPow(zeta, r, p) + ModP.ModPow(zeta, (order - r) % order, p)) % p;
    }

    /// <summary>2 X_2j(tau), the odd ladder's rung, in GF(p). The sign (-1)^(k+1) is the mirror sign the
    /// endpoint amplitude carries; it is what makes the odd orders odd.</summary>
    public static long TwoX(int n, int j, (int K1, int K2, int K3) t, long p, long zeta)
    {
        long s = 0;
        foreach (int k in new[] { t.K1, t.K2, t.K3 })
        {
            long term = TwoCos(n, 2L * j * k, p, zeta);
            s += k % 2 == 1 ? term : p - term;
        }
        return s % p;
    }

    /// <summary>2 M_m(tau), the even ladder's rung, in GF(p).</summary>
    public static long TwoM(int n, int m, (int K1, int K2, int K3) t, long p, long zeta)
    {
        long s = 0;
        foreach (int k in new[] { t.K1, t.K2, t.K3 })
            s += TwoCos(n, (long)m * k, p, zeta);
        return s % p;
    }

    /// <summary>2 DeltaX_2j for a pair, in GF(p).</summary>
    public static long DeltaTwoX(int n, int j, (int, int, int) a, (int, int, int) b, long p, long zeta)
        => ModP.Mod(TwoX(n, j, a, p, zeta) - TwoX(n, j, b, p, zeta), p);

    /// <summary>2 DeltaM_m for a pair, in GF(p).</summary>
    public static long DeltaTwoM(int n, int m, (int, int, int) a, (int, int, int) b, long p, long zeta)
        => ModP.Mod(TwoM(n, m, a, p, zeta) - TwoM(n, m, b, p, zeta), p);

    // ------------------------------------------------------------------ the gap of a pair

    /// <summary>The triple's odd-label count, which is the whole of X_0: X_0 = 2o - 3.</summary>
    public static int OddLabelCount((int K1, int K2, int K3) t)
        => (t.K1 % 2) + (t.K2 % 2) + (t.K3 % 2);

    /// <summary>DeltaX_0 = 2(o_tau - o_sigma), an integer count on labels and no cosine at all.</summary>
    public static int DeltaXZero((int, int, int) a, (int, int, int) b)
        => 2 * (OddLabelCount(a) - OddLabelCount(b));

    /// <summary>Whether the pair separates at first order: c_1 is proportional to DeltaX_0 - DeltaX_2,
    /// nonzero at either prime being an exact proof of separation. At an odd comb DeltaX_2 is dead by
    /// Theorem D and this reduces to the label count.</summary>
    public static bool FirstOrderSeparates(int n, (int, int, int) a, (int, int, int) b)
        => Settings(n).Any(s => FirstOrderSeparates(n, a, b, s.P, s.Zeta));

    static bool FirstOrderSeparates(int n, (int, int, int) a, (int, int, int) b, long p, long zeta)
        => ModP.Mod(2L * DeltaXZero(a, b) - DeltaTwoX(n, 1, a, b, p, zeta), p) != 0;

    /// <summary>Whether the pair is a theta-mirror, sigma = n - tau. At EVEN n those are carried term by
    /// term and need no shape at all, which matters because one of them is not ROT3: eta_{n-k} = eta_k
    /// and cos(2j(n-k) pi / n) = cos(2jk pi / n), so X_2j(sigma) = X_2j(tau) for EVERY j. At ODD n the
    /// same reflection flips every eta instead, so nothing is carried, and an odd-n theta-mirror always
    /// SEPARATES: o_sigma = 3 - o_tau cannot equal o_tau. The predicate is the reflection at any n; only
    /// the carrying is the even comb's.</summary>
    public static bool IsThetaMirror(int n, (int K1, int K2, int K3) a, (int K1, int K2, int K3) b)
    {
        CheckComb(n);
        var mirror = new[] { n - a.K1, n - a.K2, n - a.K3 };
        Array.Sort(mirror);
        return mirror[0] == b.K1 && mirror[1] == b.K2 && mirror[2] == b.K3;
    }

    /// <summary>All three labels of one parity, so eta is constant over the triple. The ROT3 rung
    /// lemma's second hypothesis: free at an even comb, load-bearing at an odd one.</summary>
    public static bool IsParityUniform((int K1, int K2, int K3) t)
        => t.K1 % 2 == t.K2 % 2 && t.K2 % 2 == t.K3 % 2;

    /// <summary>The Conway-Jones ROT3 shape read on the DOUBLED labels: the plus/minus label set
    /// {+-k mod n} has six elements and is the union of two cosets of the order-3 subgroup of Z/n. Needs
    /// 3|n. Seed's TripleFamily.Rot3 is the same shape on the narrower zero-sum domain.</summary>
    public static bool IsDoubledLabelRot3(int n, (int K1, int K2, int K3) t)
    {
        CheckComb(n);
        // a shortcut and not a second condition: at 3 not dividing n the step floor(n/3) has no order 3
        // (n would have to divide 3*floor(n/3) = n - (n mod 3), impossible for 0 < n mod 3 < n), so the
        // closure below already refuses every triple there. Kept because the reader should not have to
        // run that argument, and stated because a gate on it would be a gate that cannot fail.
        if (n % 3 != 0) return false;
        var p = new HashSet<int>();
        foreach (int k in new[] { t.K1, t.K2, t.K3 })
        {
            p.Add(((k % n) + n) % n);
            p.Add(((-k % n) + n) % n);
        }
        // zero is unreachable from a comb label 1..n-1 and reachable from a caller's k = 0 or k = n,
        // which is the only reason the guard is here.
        if (p.Count != 6 || p.Contains(0)) return false;
        int step = n / 3;
        // closure under one step suffices on a finite set: it forces the orbit of every element, and
        // six elements closed under an order-3 shift are exactly two cosets.
        foreach (int e in p)
            if (!p.Contains((e + step) % n)) return false;
        return true;
    }

    // ------------------------------------------------------------------ the census the law explains

    /// <summary>One comb's gap census: how the colliding pairs of LevelCollision's census behave as the
    /// road moves. Standing means the first order dies, read at both primes; the two halves of the
    /// standing set are EXPLAINED by different arguments, the theta-mirrors term by term and, at an even
    /// comb, the rest by the ROT3 shape.</summary>
    public sealed record GapCensus(
        int N,
        int Pairs,
        int Separating,
        int Standing,
        int StandingMirror,
        int StandingNonMirror,
        int SecondOrderZero);

    public GapCensus Census() => CensusOf(Ncomb);

    public static GapCensus CensusOf(int n)
    {
        CheckComb(n);
        var settings = Settings(n);
        int pairs = 0, separating = 0, mirror = 0, nonMirror = 0, secondZero = 0;
        foreach (var (a, b) in LevelCollision.CollidingPairs(n))
        {
            pairs++;
            if (settings.Any(s => FirstOrderSeparates(n, a, b, s.P, s.Zeta))) separating++;
            else if (IsThetaMirror(n, a, b)) mirror++;
            else nonMirror++;

            // c_2 = -((n-3)/n^2) DeltaM_3, and the prefactor vanishes only at n = 3, which CheckComb
            // refuses, so the second order dies exactly where the third comb reading does
            if (settings.All(s => DeltaTwoM(n, 3, a, b, s.P, s.Zeta) == 0)) secondZero++;
        }
        return new GapCensus(n, pairs, separating, mirror + nonMirror, mirror, nonMirror, secondZero);
    }

    /// <summary>What the pieces force at second order, past the census where the pair list is out of
    /// reach: family C (door 10, zero mode plus two R_5 pieces) and family L (door 70). It is a BOUND
    /// and not a count, and it is nothing but the inventory's own two closed forms.</summary>
    public static long SecondOrderZeroLowerBound(int n)
    {
        CheckComb(n);
        return LevelCollision.FamilyCount(LevelCollision.CollisionFamily.C, n)
             + LevelCollision.FamilyCount(LevelCollision.CollisionFamily.L, n);
    }
}
