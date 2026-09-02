using System.Numerics;

namespace MirrorWorld;

// The crack: the ring's wrap bond detuned to u*J, SOLVED rather than sampled (built 2026-09-01 beside
// experiments/THE_CRACKED_BELL.md, its section "The crack is exactly solvable", whose gate is stage E of
// simulations/cracked_bell_gate.py;
// registry F160, minted the same evening as an INDEX of that section and of the sibling paragraph in
// experiments/COUPLING_DEFECT_WALK_TIME_STEP.md, not as a finding of its own). Since 2026-09-02 the law
// has a proof, docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md (the determinant by Laplace, the curve, its
// SIMPLICITY at every u >= 0 except u = 1 by the factorization G = 2AB into the two reflection sectors,
// which is what licenses the root reading and the near-double-pair guard below, the join, c_m, the
// departures by Weyl, and the chain-end velocity dE_k/du = (-1)^(k+1) alpha_k/gamma_0, F65's rates
// signed), and a typed claim in the main repo, CrackedRingExactCurveClaim (RCPsiSquared.Core), which
// carries the closed forms and NOT this elimination; the elimination stays here.
//
// THE PARENT IS THE CYCLOTOMY, not the frame and not the Warble. The choice is the content: u = J'/J
// on one bond is a BOUNDARY-CONDITION parameter, and its two ends are exactly Cyclotomy's two Own
// outputs. u = 1 gives cos(Nk) = 1, the ring comb k = 2*pi*m/N (Cyclotomy.RingOrders); u = 0 gives
// sin((N+1)k) = 0, the open N-site chain comb k = pi*m/(N+1), which is F2b (Cyclotomy.PathOrders). The
// world held that pair as a SWITCH three times over (Cyclotomy.Own, Divisor's clock modulus, BlindSeat's
// two gcd laws) and as a road zero times. So this object owns the ROAD between the two combs and owns
// neither comb: both arrive in the right bucket. It is NOT a road between two topologies: for every
// u > 0 the graph is still a ring, only the endpoint u = 0 is a chain, so this says boundary condition
// and legislates nothing about the word "topology", which the repo spends twice already.
//
// Because Cyclotomy.Parent is null, Inherited here is the two combs and NOTHING ELSE. Cyclotomy itself
// inherits nothing (its parent is deliberately open), and Divisor and Marginal hang on parents that hang
// on the frame; this is the first object whose inherited bucket is non-empty and has no frame in it.
// That is the two-bucket ontology reporting where the crack hangs (the crack's arithmetic sits on the
// cyclotomic lattice, and the lattice sits on nothing of ours), not an omission.
//
// NO TIME IS OWNED HERE, and that is deliberate. The objects that hold a time (Warble's T_zero and
// T_rev, WalkTime's arrival step, a few Formulas faces) are static classes, not GameObjects; no
// GameObject owns one, and this one does not start. The road is gamma-free (gamma does not appear in G
// at all) and its outputs are a polynomial, a count and energies; the beat on top of the split, its
// clocks and the gamma-dressing are Warble's and stay there.
//
// EXACT, NOT FLOATING POINT, Seed's and Divisor's genre. The single-excitation block H_se = J*A of the
// cracked ring (adjacency, the XY book of PROOF_RING_GAP_DOMINANCE; NOT the Heisenberg Laplacian book)
// has the characteristic polynomial, in units of J and with u = p/q rational,
//
//   det(x*I - H) = U_N(x/2) - u^2 * U_{N-2}(x/2) - 2u,        U_n the Chebyshev polynomials of the
//                                                              second kind, U_n(x/2) monic with
//                                                              integer coefficients,
//
// which at x = 2 cos k is the determinant identity det(2J cos k I - H) = J^N G(k)/sin k of section E,
// G the quantization curve
//
//   G(k) = (1 - u^2) sin(Nk) cos k + [(1 + u^2) cos(Nk) - 2u] sin k
//        = sin((N+1)k) - u^2 sin((N-1)k) - 2u sin k.
//
// Section E meets the identity through a determinant to 4.1e-13 relative (gate E1c). Here it is met EXACTLY:
// the left side is the characteristic polynomial computed from the matrix by Faddeev-LeVerrier over
// the integers (the matrix scaled by q, every division exact), the right side is the Chebyshev
// recursion, and the two integer coefficient lists are compared cell for cell. Past the exact route's
// own cost (N^4 loop iterations, N^3 big-integer products since the matrix has three nonzeros per row;
// this is not the 4^N wall of the spectrum, which the polynomial
// never meets), the same identity is met mod two primes at N+1 points by a bordered elimination that
// is O(N log p) per point (one Fermat inverse per pivot), Seed's two-prime convention.
//
// THE DEPARTURE COUNT is read off that polynomial without a root: H is real symmetric, so the
// polynomial is real-rooted, and for a real-rooted polynomial Descartes' rule is EXACT (the number of
// sign variations of P(x + 2) is the number of roots above the band, of P(-x - 2) the number below).
// Section E's law, proven there from the band-bottom expansion and repaired once already (the first
// version forgot the parity of N, docs/CAUGHT_ERRORS.md 2026-08-31): at u <= 1 nothing leaves; past
// u = 1 the top level leaves at once for every N, and the bottom one at once for even N but only past
// u = (N+1)/(N-1) for odd N, because P(-2) = ((N-1)u - (N+1))(u+1) at odd N and -((N-1)u + (N+1))(u-1)
// at even N. AT the odd threshold the level sits ON the band edge, P(-2) = 0 exactly, and is counted
// as not departed; the edge is reported separately rather than folded into a tolerance. The COUNT is
// exact at every rational u.
//
// The roots themselves are a READING, not a pin: the in-band momenta are the sign changes of G on the
// OPEN interval (0, pi) refined by bisection (open because G vanishes at both ends for free), the
// departed levels the same curve continued to k = i*kappa and pi + i*kappa. The scan is guarded by the
// exact count: in-band roots + departures + edge levels must be N, and a scan too coarse to separate a
// near-double pair throws instead of returning the smaller number.
//
// ONE HONEST LIMIT OF THE READING (the counts have none). Just past u = 1 the departed level's distance
// from the band edge goes as E - 2 ~ (u-1)^2, and it is read no better than the representability floor
// 2 eps/(E - 2): 1e-9 to 1e-8 relative at u - 1 = 1e-6 (N = 4, 12, 30), 1e-5 at 1e-8, and below
// u - 1 ~ 1e-9 the distance falls under one ulp of 2, E returns 2.0 exactly, and DepartedLevels throws
// rather than returning a level the count says has left. (A first sentence here put that loss at 1e-15,
// six decades too late; the fourth review round measured it.) There is no limit at large u: the continued curve is written in the GROUPED form 2 e^(-N kappa) G, which has no
// cancellation and no overflow at any (N, u). A first draft wrote it as (1-u^2) sinh(N kappa) cosh kappa
// + ..., whose two terms are each ~u^2 e^(N kappa) near the root kappa = ln u and cancel to nothing
// (e^kappa - u^2 e^(-kappa) is what is left, algebraically), so it lost every digit past u ~ 1e3 at
// N = 12 and returned 1e5.0107 for an exact 1e5.00001; a guard on N ln u then protected the wrong
// variable, and a test row certified the wrong number as fine. The third review round's mathematician
// caught it with the exact polynomial's own Newton step, which the neighbouring test already carried.
// And the modular route
// proves equality mod p, not over Z: what it cannot exclude is a residual divisible by both primes,
// Seed's accepted risk; the exact route (CharacteristicPolynomial) has no such gap and is the one the
// tests run to N = 21.
public sealed class Crack : GameObject
{
    public int N { get; }
    public long UNum { get; }                // u = UNum / UDen, reduced; u < 1 weakens (the crack proper), u > 1 strengthens
    public long UDen { get; }
    public double U => (double)UNum / UDen;

    // two primes below 2^31, so that a product of two residues fits a long; Seed's convention that a
    // coincidence mod one prime does not survive the second.
    static readonly long[] Primes = { 2147483647L, 999999937L };

    public Crack(Cyclotomy comb, int n, long uNum, long uDen = 1) : base(comb)
    {
        if (n < 3) throw new ArgumentOutOfRangeException(nameof(n), "a ring needs N >= 3");
        if (uDen <= 0) throw new ArgumentOutOfRangeException(nameof(uDen), "u = p/q needs q > 0");
        if (uNum < 0) throw new ArgumentOutOfRangeException(nameof(uNum), "u = J'/J is a ratio of couplings, u >= 0");
        long g = Gcd(uNum, uDen);
        N = n; UNum = uNum / g; UDen = uDen / g;
    }

    // left: the road between the two combs, and what the road costs, the departures: the count of
    // levels the road pushes out of the band past u = 1 (a departure here is that, not the departure from
    // normality of the F2b corollary and F89; a first version filed that word under F86, which never uses
    // it). NOT the combs: both ends are the Cyclotomy's, and neither is this object's to own.
    public override IReadOnlyList<string> Own => new[] { "road", "departures" };

    static long Gcd(long a, long b) { a = Math.Abs(a); b = Math.Abs(b); while (b != 0) (a, b) = (b, a % b); return a == 0 ? 1 : a; }

    // ------------------------------------------------------------------ the curve, both written forms

    public static double G(double k, int n, double u)
        => (1 - u * u) * Math.Sin(n * k) * Math.Cos(k) + ((1 + u * u) * Math.Cos(n * k) - 2 * u) * Math.Sin(k);

    // the second printed line of section E; kept so the two lines are checked against each other
    public static double GSecondForm(double k, int n, double u)
        => Math.Sin((n + 1) * k) - u * u * Math.Sin((n - 1) * k) - 2 * u * Math.Sin(k);

    // ------------------------------------------------------------------ the exact polynomials

    // U_n(x/2), ascending coefficients: U_0 = 1, U_1 = x, U_n = x U_{n-1} - U_{n-2}. Monic, integer.
    // The same polynomial is docs/proofs/PROOF_F139_SEAM_IDENTITY.md's S_m (the 2cos normalization); it
    // is named U here and not S because Renewal.cs already spends S_n for the survival function.
    public static BigInteger[] ChebyshevSecondKind(int n)
    {
        if (n < 0) return new[] { BigInteger.Zero };
        var prev = new[] { BigInteger.One };
        if (n == 0) return prev;
        var cur = new[] { BigInteger.Zero, BigInteger.One };
        for (int m = 2; m <= n; m++)
        {
            var next = new BigInteger[m + 1];
            for (int i = 0; i < cur.Length; i++) next[i + 1] += cur[i];
            for (int i = 0; i < prev.Length; i++) next[i] -= prev[i];
            prev = cur; cur = next;
        }
        return cur;
    }

    // 2 T_n(x/2) - 2, the perfect ring's characteristic polynomial prod (x - 2 cos(2 pi m/N)):
    // D_0 = 2, D_1 = x, D_n = x D_{n-1} - D_{n-2}, then minus 2. A second recursion, independent of U.
    public static BigInteger[] RingPolynomial(int n)
    {
        if (n == 0) return new[] { BigInteger.Zero };           // 2 T_0 - 2 = 0
        var prev = new[] { new BigInteger(2) };
        var cur = new[] { BigInteger.Zero, BigInteger.One };
        for (int m = 2; m <= n; m++)
        {
            var next = new BigInteger[m + 1];
            for (int i = 0; i < cur.Length; i++) next[i + 1] += cur[i];
            for (int i = 0; i < prev.Length; i++) next[i] -= prev[i];
            prev = cur; cur = next;
        }
        cur[0] -= 2;
        return cur;
    }

    // q^N * det(x I - H) = q^N U_N - p^2 q^(N-2) U_{N-2} - 2 p q^(N-1), ascending, integer. The road.
    BigInteger[]? roadCache;

    public BigInteger[] RoadPolynomial()
    {
        roadCache ??= BuildRoad();
        return (BigInteger[])roadCache.Clone();
    }

    BigInteger[] BuildRoad()
    {
        BigInteger p = UNum, q = UDen;
        var uN = ChebyshevSecondKind(N);
        var uN2 = ChebyshevSecondKind(N - 2);
        var qN = BigInteger.Pow(q, N);
        var road = new BigInteger[N + 1];
        for (int i = 0; i <= N; i++) road[i] = qN * uN[i];
        var w = p * p * BigInteger.Pow(q, N - 2);
        for (int i = 0; i < uN2.Length; i++) road[i] -= w * uN2[i];
        road[0] -= 2 * p * BigInteger.Pow(q, N - 1);
        return road;
    }

    // q^N * det(x I - H) from the MATRIX, by Faddeev-LeVerrier on the integer matrix q*H (chain bonds q,
    // wrap bond p): every division by k is exact for an integer matrix. Ascending coefficients.
    // N^4 loop iterations (N^3 big-integer products, three nonzeros per row): the exact route up to a few
    // dozen sites; past that, IdentityModP.
    public BigInteger[] CharacteristicPolynomial()
    {
        int n = N;
        var a = new BigInteger[n, n];
        for (int i = 0; i < n - 1; i++) { a[i, i + 1] = UDen; a[i + 1, i] = UDen; }
        a[0, n - 1] += UNum; a[n - 1, 0] += UNum;
        var c = new BigInteger[n + 1];
        c[n] = BigInteger.One;
        var m = new BigInteger[n, n];
        for (int k = 1; k <= n; k++)
        {
            // M_k = A M_{k-1} + c_{n-k+1} I
            var am = new BigInteger[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    BigInteger s = BigInteger.Zero;
                    for (int l = 0; l < n; l++) if (!a[i, l].IsZero) s += a[i, l] * m[l, j];
                    am[i, j] = s;
                }
            for (int i = 0; i < n; i++) am[i, i] += c[n - k + 1];
            m = am;
            // c_{n-k} = -tr(A M_k) / k
            BigInteger tr = BigInteger.Zero;
            for (int i = 0; i < n; i++)
                for (int l = 0; l < n; l++) if (!a[i, l].IsZero) tr += a[i, l] * m[l, i];
            BigInteger num = -tr;
            if (!BigInteger.Remainder(num, k).IsZero) throw new InvalidOperationException("Faddeev-LeVerrier division is not exact: the matrix is not integer");
            c[n - k] = num / k;
        }
        // det(x I - qH) = sum c_j x^j; det(x I - H) = q^-N det(q x I - qH), so q^N det(x I - H) has coefficient c_j q^j.
        var road = new BigInteger[n + 1];
        for (int j = 0; j <= n; j++) road[j] = c[j] * BigInteger.Pow(UDen, j);
        return road;
    }

    // the identity, exactly: the coefficient list of the matrix route minus the Chebyshev route.
    public BigInteger[] IdentityResidual()
    {
        var chi = CharacteristicPolynomial();
        var road = RoadPolynomial();
        var r = new BigInteger[N + 1];
        for (int i = 0; i <= N; i++) r[i] = chi[i] - road[i];
        return r;
    }

    public bool IdentityHolds() => IdentityResidual().All(x => x.IsZero);

    // Past the wall: the identity mod a prime, at N+1 integer points x0 = 1, 2, ... (x0 = 0 is always a
    // dead pivot, a = 0, and any other dead pivot is skipped the same way).
    // Left side det(q x0 I - qH) mod prime by a bordered elimination of the cyclic tridiagonal, O(N log p)
    // matrix (a Gaussian elimination that knows the sparsity and nothing else); right side the road
    // polynomial mod prime by Horner. Two monic-scaled polynomials of degree N that agree at N+1 points
    // mod prime agree mod prime. Returns the number of mismatching points; 0 is the pass.
    public int IdentityMismatchesModP(long prime) => MismatchesAgainst(RoadPolynomial(), prime);

    // the same route with the road supplied from outside, so a WRONG road can be fed to it (the control).
    public int MismatchesAgainst(BigInteger[] road, long prime)
    {
        if (road.Length != N + 1) throw new ArgumentException("a road for this ring has N + 1 coefficients", nameof(road));
        var roadMod = new long[N + 1];
        for (int i = 0; i <= N; i++) roadMod[i] = (long)(((road[i] % prime) + prime) % prime);
        int agreed = 0, mismatches = 0;
        for (long x0 = 1; agreed < N + 1; x0++)
        {
            if (x0 > 4L * N + 8) throw new InvalidOperationException("too many dead pivots mod " + prime);
            long? det = CyclicTridiagonalDeterminantModP(x0, prime);
            if (det is null) continue;                       // a zero pivot: this point is skipped, not counted
            long rhs = 0;
            for (int i = N; i >= 0; i--) rhs = (MulMod(rhs, x0 % prime, prime) + roadMod[i]) % prime;
            agreed++;
            if (det.Value != rhs) mismatches++;
        }
        return mismatches;
    }

    public int IdentityMismatchesModTwoPrimes() => Primes.Sum(IdentityMismatchesModP);

    static long MulMod(long a, long b, long p) => (long)((UInt128)(ulong)a * (ulong)b % (ulong)p);

    static long PowMod(long b, long e, long p) { long r = 1; b %= p; while (e > 0) { if ((e & 1) == 1) r = MulMod(r, b, p); b = MulMod(b, b, p); e >>= 1; } return r; }

    // det(M) mod p for M = q x0 I - qH: diagonal a = q x0, chain bonds -q, wrap bond -p, N >= 3.
    // Rows 0..N-2 are eliminated in order; the fill-in lives only in the last column and the last row.
    long? CyclicTridiagonalDeterminantModP(long x0, long p)
    {
        long a = MulMod(((UDen % p) + p) % p, x0 % p, p);
        long bond = (p - UDen % p) % p;                       // -q
        long corner = (p - UNum % p) % p;                     // -p_u
        int n = N;
        long det = 1;
        long d = a;                                           // current row's pivot
        long last = corner;                                   // current row's entry in the last column
        long r = corner;                                      // last row's entry in the current column
        long s = a;                                           // last row's diagonal
        for (int i = 0; i <= n - 2; i++)
        {
            if (i == n - 2) last = (last + bond) % p;         // row N-2's right neighbour IS the last column
            if (d == 0) return null;
            det = MulMod(det, d, p);
            long inv = PowMod(d, p - 2, p);
            // the last row loses its column-i entry
            long g = MulMod(r, inv, p);
            s = (s - MulMod(g, last, p) + p) % p;
            if (i == n - 2) break;
            // the next row's left entry (bond) is eliminated by row i
            long f = MulMod(bond, inv, p);
            long dNext = (a - MulMod(f, bond, p) + p) % p;
            long lastNext = (p - MulMod(f, last, p)) % p;     // next row's last-column entry, 0 before fill-in
            // the last row's column-(i+1) entry: -g * bond, plus its own entry there (bond when i+1 == n-2)
            long rNext = (p - MulMod(g, bond, p)) % p;
            if (i + 1 == n - 2) rNext = (rNext + bond) % p;
            d = dNext; last = lastNext; r = rNext;
        }
        return MulMod(det, s, p);
    }

    // ------------------------------------------------------------------ the departure count

    // Taylor shift P(x + c), ascending coefficients, by repeated synthetic division. Exact.
    static BigInteger[] Shift(BigInteger[] poly, BigInteger c)
    {
        var a = (BigInteger[])poly.Clone();
        int n = a.Length;
        for (int i = 0; i < n; i++)
            for (int j = n - 2; j >= i; j--) a[j] += c * a[j + 1];
        return a;
    }

    // sign variations of the nonzero coefficients: for a real-rooted polynomial, EXACTLY the number of
    // positive roots (a zero root, a vanishing constant term, is not positive and drops out).
    static int Descartes(BigInteger[] poly)
    {
        int v = 0, prev = 0;
        foreach (var c in poly)
        {
            int sgn = c.Sign;
            if (sgn == 0) continue;
            if (prev != 0 && sgn != prev) v++;
            prev = sgn;
        }
        return v;
    }

    static BigInteger Eval(BigInteger[] poly, BigInteger x)
    {
        BigInteger r = BigInteger.Zero;
        for (int i = poly.Length - 1; i >= 0; i--) r = r * x + poly[i];
        return r;
    }

    // levels strictly above the band, read off the road: roots of P(x + 2) with x > 0.
    public int DeparturesAbove()
    {
        var road = RoadPolynomial();
        return Descartes(Shift(road, 2));
    }

    // levels strictly below the band: roots of P(-x - 2) with x > 0 (shift by -2, then x -> -x).
    public int DeparturesBelow()
    {
        var shifted = Shift(RoadPolynomial(), -2);
        for (int i = 1; i < shifted.Length; i += 2) shifted[i] = -shifted[i];
        return Descartes(shifted);
    }

    public int Departures => DeparturesAbove() + DeparturesBelow();

    // levels sitting ON the band edge, P(+-2) = 0 exactly: u = 1 (top, every N; bottom, even N) and the
    // odd threshold u = (N+1)/(N-1) (bottom). Reported, never tolerated into the count.
    public (bool top, bool bottom) OnTheEdge()
    {
        var road = RoadPolynomial();
        return (Eval(road, 2).IsZero, Eval(road, -2).IsZero);
    }

    public int EdgeLevels() { var (t, b) = OnTheEdge(); return (t ? 1 : 0) + (b ? 1 : 0); }

    // section E's law as a predicate over the integers: 0 at u <= 1; past u = 1 the top level for every
    // N, the bottom one at once for even N and only PAST u = (N+1)/(N-1) for odd N (at the threshold it
    // is on the edge). Compared exactly as p(N-1) > q(N+1).
    public static int DepartureLaw(int n, long uNum, long uDen)
    {
        if (n < 3) throw new ArgumentOutOfRangeException(nameof(n), "a ring needs N >= 3");
        if (uNum <= uDen) return 0;
        if (n % 2 == 0) return 2;
        return (BigInteger)uNum * (n - 1) > (BigInteger)uDen * (n + 1) ? 2 : 1;   // BigInteger: no silent overflow at large u
    }

    public int LawCount => DepartureLaw(N, UNum, UDen);

    // the odd ring's threshold u = (N+1)/(N-1), reduced; in the sibling's strengthening convention
    // u = 1 + delta it reads delta = 2/(N-1). Meaningless at even N, where the prefactor (1-u) vanishes.
    public static (long num, long den) OddThreshold(int n)
    {
        if (n % 2 == 0) throw new ArgumentException("the threshold exists at odd N only; the even ring sheds its level at every u > 1", nameof(n));
        long g = Gcd(n + 1, n - 1);
        return ((n + 1) / g, (n - 1) / g);
    }

    // ------------------------------------------------------------------ the reading: roots of the curve

    // sign changes of G on the OPEN interval (0, pi), each refined by bisection to the last double.
    // At u = 1 the roots are double and a sign count is blind to them, so the ring comb is returned
    // (that end of the road is the parent's, not scanned). The count is tied to the exact one: in-band
    // roots + departures + edge levels must be N, else the scan was too coarse and this THROWS.
    double[]? momentaCache;              // the default-resolution scan, kept so the pair table is one scan

    public double[] InBandMomenta(int points = 0)
    {
        if (points <= 0 && momentaCache is not null) return (double[])momentaCache.Clone();
        var ks = ScanInBandMomenta(points);
        if (points <= 0) momentaCache = (double[])ks.Clone();
        return ks;
    }

    double[] ScanInBandMomenta(int points)
    {
        if (UNum == UDen)
        {
            // the interior of the ring comb: k = 2 pi m/N with m = 1..N-1, minus k = pi at even N (m = N/2),
            // which is the band edge E = -2 and belongs to OnTheEdge, as k = 0 (m = 0) does at every N.
            var comb = new List<double>();
            // listed with multiplicity: m and N-m share E = 2cos(2 pi m/N), i.e. the one k in (0, pi)
            for (int m = 1; m <= N - 1; m++) if (2 * m != N) comb.Add(2.0 * Math.PI * Math.Min(m, N - m) / N);
            return comb.ToArray();
        }
        double u = U;
        if (points <= 0)
        {
            double delta = Math.Abs(1.0 - u);
            points = (int)Math.Min(20_000_000, Math.Max(200_001, Math.Ceiling(64.0 * N / delta)));
        }
        var roots = new List<double>();
        double lo = Math.PI / points, glo = G(lo, N, u);
        for (int i = 2; i < points; i++)
        {
            double hi = Math.PI * i / points, ghi = G(hi, N, u);
            if (glo == 0.0) roots.Add(lo);
            else if (glo * ghi < 0) roots.Add(Bisect(k => G(k, N, u), lo, hi, glo));
            lo = hi; glo = ghi;
        }
        int expected = N - Departures - EdgeLevels();
        if (roots.Count != expected)
            throw new InvalidOperationException($"N={N}, u={UNum}/{UDen}: the scan found {roots.Count} in-band roots, the exact count says {expected}; the scan is too coarse for a near-double pair");
        return roots.ToArray();
    }

    static double Bisect(Func<double, double> f, double lo, double hi, double flo)
    {
        for (int it = 0; it < 200; it++)
        {
            double mid = 0.5 * (lo + hi);
            if (mid <= lo || mid >= hi) break;
            double fm = f(mid);
            if (fm == 0.0) return mid;
            if ((fm < 0) == (flo < 0)) { lo = mid; flo = fm; } else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    // G continued above the band, k = i*kappa, E = +2 cosh kappa: 2 e^(-N kappa) * G(i kappa)/i, real, the
    // positive factor keeping the zeros. Written in the GROUPED form
    //   (e^kappa - u^2 e^(-kappa)) + e^(-2N kappa) (u^2 e^kappa - e^(-kappa)) - 4u e^(-N kappa) sinh kappa,
    // which is (1-u^2) sinh(N kappa) cosh kappa + ((1+u^2) cosh(N kappa) - 2u) sinh kappa times 2 e^(-N kappa)
    // (sympy, symbolic in u), with no cancellation near the root kappa ~ ln u and no overflow at any (N, u).
    public static double GAbove(double kappa, int n, double u)
        => (Math.Exp(kappa) - u * u * Math.Exp(-kappa))
         + Math.Exp(-2.0 * n * kappa) * (u * u * Math.Exp(kappa) - Math.Exp(-kappa))
         - 4 * u * Math.Exp(-n * kappa) * Math.Sinh(kappa);

    // G continued below the band, k = pi + i*kappa, E = -2 cosh kappa: 2 e^(-N kappa) * G/(-i), real, the same
    // grouped form with (-1)^N on the two bracketed terms and not on the 2u term; that parity is what makes
    // the odd ring's second departure a threshold and the even ring's immediate.
    public static double GBelow(double kappa, int n, double u)
    {
        double s = n % 2 == 0 ? 1.0 : -1.0;
        return s * ((Math.Exp(kappa) - u * u * Math.Exp(-kappa)) + Math.Exp(-2.0 * n * kappa) * (u * u * Math.Exp(kappa) - Math.Exp(-kappa)))
             - 4 * u * Math.Exp(-n * kappa) * Math.Sinh(kappa);
    }

    // the departed levels, in units of J, bisected on kappa > 0 up to the Gershgorin edge (1 + u).
    public double[] DepartedLevels()
    {
        int above = DeparturesAbove(), below = DeparturesBelow();
        if (above > 1 || below > 1) throw new InvalidOperationException("more than one level per side left the band; a rank-two update of signature (1,1) cannot do that");
        var levels = new List<double>();
        double u = U;
        double kappaMax = Math.Acosh(Math.Max(1.0, (1.0 + u) / 2.0)) + 1e-9;
        if (above == 1) levels.Add(2.0 * Math.Cosh(RootOnKappa(k => GAbove(k, N, u), kappaMax)));
        if (below == 1) levels.Add(-2.0 * Math.Cosh(RootOnKappa(k => GBelow(k, N, u), kappaMax)));
        // the one place a reading could return a number the count contradicts: just past u = 1 the distance
        // from the edge falls under one ulp of 2 and 2 cosh kappa rounds to 2.0; refuse, as the scan does
        if (levels.Any(e => Math.Abs(e) <= 2.0))
            throw new InvalidOperationException($"N={N}, u={UNum}/{UDen}: a departed level's distance from the band edge is below one ulp of 2 (E - 2 ~ (u-1)^2); the count says it left, the reading cannot resolve it");
        return levels.ToArray();
    }

    static double RootOnKappa(Func<double, double> f, double kappaMax)
    {
        const int pts = 200_001;
        double lo = kappaMax / pts, flo = f(lo);
        for (int i = 2; i <= pts; i++)
        {
            double hi = kappaMax * i / pts, fhi = f(hi);
            if (flo * fhi < 0) return Bisect(f, lo, hi, flo);
            lo = hi; flo = fhi;
        }
        throw new InvalidOperationException("the departed level was counted but not bracketed on kappa");
    }

    // every level in units of J, descending: the in-band ones E = 2 cos k plus the departed ones.
    public double[] Levels(int points = 0)
    {
        var e = new List<double>();
        if (UNum == UDen)
        {
            for (int m = 0; m < N; m++) e.Add(2.0 * Math.Cos(2.0 * Math.PI * m / N));
        }
        else
        {
            foreach (double k in InBandMomenta(points)) e.Add(2.0 * Math.Cos(k));
            var (top, bottom) = OnTheEdge();
            if (top) e.Add(2.0);
            if (bottom) e.Add(-2.0);
            e.AddRange(DepartedLevels());
        }
        e.Sort((x, y) => y.CompareTo(x));
        return e.ToArray();
    }

    // the m-th degenerate pair's split, read off the curve alone: the two in-band levels nearest the
    // perfect ring's E_m = 2 cos(2 pi m/N). Valid for 1 <= m < N/2 (m = 0 and m = N/2 are unpaired and
    // shift instead), and for a crack small enough that the pair is still the nearest two.
    public double SplitOfPair(int m, int points = 0)
    {
        if (m < 1 || 2 * m >= N) throw new ArgumentOutOfRangeException(nameof(m), "a paired mode has 1 <= m < N/2");
        if (UNum > UDen) throw new InvalidOperationException("the split is read on the crack (u <= 1); past u = 1 a departed level sits beside the pair and the nearest-two rule is not the pair");
        if (UNum == UDen) return 0.0;      // the perfect ring: the pair is degenerate, the split is zero by construction
        double em = 2.0 * Math.Cos(2.0 * Math.PI * m / N);
        var byDistance = InBandMomenta(points).Select(k => 2.0 * Math.Cos(k)).OrderBy(e => Math.Abs(e - em)).ToArray();
        double a = byDistance[0], b = byDistance[1];
        // the fence, enforced rather than described: the partners shift by -(2 delta J/N)(cos theta_m +- 1),
        // one up and one down, so the pair STRADDLES E_m, and the next level must sit farther from E_m
        // than the pair's own width; a crack too deep for that (N = 25, u = 1/2 fails both) would
        // otherwise return a wrong number with no signal
        bool straddle = (a - em) * (b - em) < 0;
        bool isolated = byDistance.Length < 3 || Math.Abs(byDistance[2] - em) > Math.Abs(a - b);
        if (!straddle || !isolated)
            throw new InvalidOperationException($"N={N}, u={UNum}/{UDen}, m={m}: the two levels nearest E_m are not an isolated straddling pair; the crack is too deep for the pair reading");
        return Math.Abs(a - b);
    }

}
