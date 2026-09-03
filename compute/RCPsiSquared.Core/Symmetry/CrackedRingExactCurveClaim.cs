using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F160, the cracked ring is exactly solvable: the XY ring of N sites, every bond J except the
/// wrap bond N-1 &lt;-&gt; 0, which carries u*J (u = J'/J), has the single-excitation characteristic polynomial
///
/// <code>
///     det(x*I - H) = U_N(x/2) - u^2 * U_{N-2}(x/2) - 2u              (units of J; U_n Chebyshev, second kind)
///
///     at x = 2 cos k:   det(2J cos k * I - H) = J^N * G(k) / sin k,
///     G(k) = (1 - u^2) sin(Nk) cos k + [(1 + u^2) cos(Nk) - 2u] sin k = sin((N+1)k) - u^2 sin((N-1)k) - 2u sin k,
/// </code>
///
/// so for 0 &lt;= u &lt; 1 the whole spectrum is the zero set of G on the open interval (0, pi), multiplicities
/// included, and past u = 1 the departed levels ride the same curve continued to complex k. The curve
/// factors into the two reflection sectors, G = 2AB with A = cos(k(N+1)/2) - u cos(k(N-1)/2) and
/// B = sin(k(N+1)/2) + u sin(k(N-1)/2), each a nonvanishing prefactor times an unreduced Jacobi block's
/// characteristic polynomial, and the Bezout identity (cos a + u cos b)A + (sin a - u sin b)B = 1 - u^2
/// forbids a common zero unless u^2 = 1: the spectrum is SIMPLE at every u &gt;= 0 with u != 1, a theorem
/// where the experiment page had a 1050-point sweep.
///
/// <para><b>The road, and what u is.</b> u = 1 gives cos(Nk) = 1, the perfect ring's comb k = 2 pi m/N; u = 0
/// gives sin((N+1)k) = 0, the OPEN N-site chain's comb k = pi m/(N+1), which is F2b read as the matching
/// condition psi_{-1} = u psi_{N-1}, psi_N = u psi_0 with u = 0. So u interpolates the ring's modulus N and the
/// chain's modulus N+1 and is a BOUNDARY-CONDITION parameter: for every u &gt; 0 the graph is still a ring, only
/// the endpoint is a chain. The word topology is not used for u here (the repo spends it twice already).</para>
///
/// <para><b>The join.</b> 1/t(k) = (1 + u^2)/(2u) + i (1 - u^2) cot k/(2u) is the inverse of the transmission
/// amplitude the walk-time step carries for the same bond on the infinite chain, and
/// Re[e^(-iNk)/t(k)] - 1 = G(k)/(2u sin k): the ring's quantization condition is that amplitude with one round
/// trip of phase, the chain's scatterer closed into a loop.</para>
///
/// <para><b>The split's next order.</b> Near k_m = 2 pi m/N the curve truncates to a quadratic in x = N(k - k_m)
/// whose roots are the committed branch shifts of the flat split 4 delta J/N, delta = 1 - u, and the exact law
/// adds the order the first-order reading could only measure:
/// Delta E_m = (4 delta J/N) [1 + delta c_m + O(delta^2)], c_m = 1/2 - 1/(N sin^2 k_m). c_m is the delta -&gt; 0
/// form; at a finite delta the split's own next order has already moved its zero.</para>
///
/// <para><b>The departures, without a root.</b> U_n(+-1) = (+-1)^n (n+1) makes the band-edge values two linear
/// factors: P(+2) = -((N-1)u + (N+1))(u - 1) at every N; P(-2) = ((N-1)u - (N+1))(u + 1) at odd N and
/// -((N-1)u + (N+1))(u - 1) at even N. H is real symmetric, so P is real-rooted and a sign at the edge counts
/// an odd number of roots beyond it; the crack is the rank-two update V = u(|0&gt;&lt;N-1| + |N-1&gt;&lt;0|) with
/// eigenvalues (u, 0, ..., -u), so by Weyl's inequality lambda_2(H) &lt;= lambda_1(H_path) + 0 &lt; 2 and at
/// most one level passes each edge. Hence: at u &lt;= 1 nothing leaves (Perron-Frobenius puts every level strictly
/// inside for u &lt; 1); past u = 1 the top level leaves at every N, the bottom one at every u &gt; 1 for even N and
/// only past u = (N+1)/(N-1) for odd N, where P(-2) = 0 exactly and the level sits ON the edge, reported as an
/// edge level and not folded into a tolerance. The first draft of this count forgot the parity of N
/// (docs/CAUGHT_ERRORS.md 2026-08-31: a count is a claim).</para>
///
/// <para><b>The velocity at the chain end, where two routes meet.</b> Differentiating the polynomial at u = 0
/// gives dE_k/du = 2/P'(x_k) = (-1)^(k+1) (4/(N+1)) sin^2(k pi/(N+1)) for the chain level k = 1..N; the
/// eigenvector route gives the same number as 2 psi_k(0) psi_k(N-1). That is F65's endpoint rate
/// alpha_k/gamma_0 carrying an alternating sign, so the crack's first-order motion of the chain comb IS the
/// dephased endpoint's rates, read as a comb. It is what decides which of F129's comb collisions survive the road to first
/// order (experiments/THE_COMB_ON_THE_ROAD.md).</para>
///
/// <para><b>What this class is NOT.</b> The XY adjacency book only (PROOF_RING_GAP_DOMINANCE's block, not the
/// Heisenberg Laplacian of D10 and not the (1,1) Haken-Strobl block); uniform bonds off the crack; no gamma
/// anywhere in G (the Liouvillian rate of these coherences is the Absorption Theorem's -2 gamma and is not this
/// claim's); no time (the beat on top of the split, its clocks and the gamma-dressing are the experiment's and
/// MirrorWorld's Warble); not the blind seat (THE_SEAT_THAT_CUTS's detuned-bond item is a different object);
/// what does cross is Corollary B's fold, which PROOF_BLIND_SEAT_TWO_AXES spends on the
/// blind seat's locus and, in its section (g), on its count.
/// The knob's sign: the crack weakens, u = 1 - delta; the walk-time step strengthens, u = 1 + delta; the
/// statement in u is the one to port.</para>
///
/// <para>Proof <c>docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md</c>, whose nine symbolic gates are
/// <c>simulations/cracked_ring_exact_curve_proof.py</c> (exact, no floating point); gate
/// <c>simulations/cracked_bell_gate.py</c> stage E (the curve through an eigensolver, 27 gates). Typed sibling
/// in the sober base: <c>compute/MirrorWorld/Crack.cs</c>, which meets the identity EXACTLY over the integers
/// (Faddeev-LeVerrier against the Chebyshev recursion, and mod two primes to N = 1001) and owns the road and
/// the departures; this claim carries the closed forms and deliberately not the elimination.</para></summary>
public sealed class CrackedRingExactCurveClaim : Claim
{
    /// <summary>Parent: the topology band edge, whose ring row is this curve's u = 1 end and whose §Scope fence
    /// in PROOF_RING_GAP_DOMINANCE is this curve's Perron root (the committed 1.9999500 at N = 4, delta = 1e-4
    /// is 2 - delta/2 + delta^2/8 - delta^4/128 + ... read off the same polynomial). The edge runs this way
    /// because that claim states the band edge of the UNIFORM ring and this one refines the ring's spectrum
    /// along one bond; the parent knows the crack only as the deformation that lifts its N = 4 exception.</summary>
    public TopologyBandEdgeClaim BandEdge { get; }

    /// <summary>Parent: F2b, the open chain's spectrum, which is this curve's u = 0 end read as a matching
    /// condition with the two virtual sites set to zero. The chain comb is the one object the road leaves and
    /// returns to; the claim stands on F2b for the endpoint and adds the road.</summary>
    public F2bXyChainSpectrumPi2Inheritance ChainEnd { get; }

    public CrackedRingExactCurveClaim(TopologyBandEdgeClaim bandEdge, F2bXyChainSpectrumPi2Inheritance chainEnd)
        : base("F160 the cracked ring is exactly solvable: the XY ring of N sites with every bond J except the wrap " +
               "bond, which carries u*J (u = J'/J), has the single-excitation characteristic polynomial " +
               "det(x I - H) = U_N(x/2) - u^2 U_{N-2}(x/2) - 2u in units of J, so at x = 2 cos k the spectrum is the " +
               "zero set on (0, pi) of G(k) = (1 - u^2) sin(Nk) cos k + [(1 + u^2) cos(Nk) - 2u] sin k, multiplicities " +
               "included, for 0 <= u < 1, and simple for every u >= 0 except u = 1 (G = 2AB, the two reflection sectors); " +
               "u is a boundary-condition parameter whose two ends are the ring comb (u = 1) and " +
               "the open chain comb (u = 0, F2b); the condition is the walk-time step's transmission amplitude with one " +
               "round trip of phase, Re[e^(-iNk)/t(k)] = 1; the flat split 4 delta J/N of every pair m <-> N-m gets its " +
               "next order Delta E_m = (4 delta J/N)[1 + delta (1/2 - 1/(N sin^2 k_m)) + O(delta^2)]; past u = 1 the " +
               "number of levels leaving the band is a parity law read off two linear factors, P(+2) = " +
               "-((N-1)u + (N+1))(u - 1) at every N and P(-2) = ((N-1)u - (N+1))(u + 1) at odd N, " +
               "-((N-1)u + (N+1))(u - 1) at even N, at most one per side by Weyl's inequality for the rank-two " +
               "update with eigenvalues (u, 0, ..., -u): the top " +
               "level leaves at every u > 1, the bottom one at every u > 1 for even N and only past u = (N+1)/(N-1) " +
               "for odd N, where it sits exactly on the edge; and the road's first-order velocity at the chain end is " +
               "F65's endpoint rate comb with alternating sign, dE_k/du = (-1)^(k+1) (4/(N+1)) sin^2(k pi/(N+1)). " +
               "XY adjacency book only, uniform bonds off the crack, gamma absent from G, no time owned here",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md (primary: the determinant, the curve, the two ends, the " +
               "join, the split's next order, the departure count, the chain-end velocity) + " +
               "experiments/THE_CRACKED_BELL.md (the section 'The crack is exactly solvable', where the law was derived " +
               "and gated 2026-08-31, stage E of simulations/cracked_bell_gate.py) + " +
               "experiments/COUPLING_DEFECT_WALK_TIME_STEP.md (the transmission amplitude since 2026-07-12 and the " +
               "paragraph 'The same amplitude quantizes a ring with one deformed bond') + " +
               "experiments/THE_COMB_ON_THE_ROAD.md (the chain-end velocity as F65's comb, and F129's collisions on " +
               "the road) + " +
               "docs/ANALYTICAL_FORMULAS.md (F160; F2b the u = 0 end, F65 the rate comb, F129 the comb law, F122 the " +
               "ring degeneracy it lifts) + " +
               "docs/proofs/PROOF_RING_GAP_DOMINANCE.md (the uniform ring and the Scope fence this curve's Perron " +
               "root sits in) + " +
               "docs/proofs/PROOF_F139_SEAM_IDENTITY.md (the Chebyshev polynomial S_m in the 2cos normalization, " +
               "cited not re-derived) + " +
               "compute/MirrorWorld/Crack.cs (the adopted object: the identity met exactly over the integers, the " +
               "Descartes count, the root reading, run mode crack N [u])")
    {
        BandEdge = bandEdge ?? throw new ArgumentNullException(nameof(bandEdge));
        ChainEnd = chainEnd ?? throw new ArgumentNullException(nameof(chainEnd));
    }

    /// <summary>The smallest N this claim speaks about: a ring. N = 2 would make the wrap bond a second bond
    /// between the same two sites (the identity still holds there, P = x^2 - (1+u)^2), and N = 1 has no bond.</summary>
    public const int MinN = 3;

    /// <summary>U_n(x/2) as ascending integer coefficients: p_0 = 1, p_1 = x, p_n = x p_{n-1} - p_{n-2}. Monic of
    /// degree n; at x = 2 cos k it is sin((n+1)k)/sin k, the S_n of PROOF_F139_SEAM_IDENTITY.</summary>
    public static BigInteger[] ChebyshevSecondKindMonic(int n)
    {
        if (n < 0) throw new ArgumentOutOfRangeException(nameof(n), n, "n must be non-negative");
        var prev = new BigInteger[] { BigInteger.One };
        if (n == 0) return prev;
        var cur = new BigInteger[] { BigInteger.Zero, BigInteger.One };
        for (int m = 2; m <= n; m++)
        {
            var next = new BigInteger[m + 1];
            for (int i = 0; i < cur.Length; i++) next[i + 1] += cur[i];
            for (int i = 0; i < prev.Length; i++) next[i] -= prev[i];
            prev = cur;
            cur = next;
        }
        return cur;
    }

    /// <summary>uDen^2 * det(x I - H) as ascending integer coefficients, u = uNum/uDen with uDen &gt; 0 (any
    /// representation; the result is uDen^2 P(x) for each): uDen^2 U_N(x/2) - uNum^2 U_{N-2}(x/2) - 2 uNum uDen.
    /// The scaling by uDen^2 is what keeps the polynomial over the integers; at uDen = 1 it is det(x I - H) itself.</summary>
    public static BigInteger[] RoadPolynomialScaled(int n, long uNum, long uDen)
    {
        CheckRing(n);
        CheckU(uNum, uDen);
        var uN = ChebyshevSecondKindMonic(n);
        var uN2 = ChebyshevSecondKindMonic(n - 2);
        var p = new BigInteger[n + 1];
        BigInteger den2 = (BigInteger)uDen * uDen, num2 = (BigInteger)uNum * uNum;
        for (int i = 0; i < uN.Length; i++) p[i] += den2 * uN[i];
        for (int i = 0; i < uN2.Length; i++) p[i] -= num2 * uN2[i];
        p[0] -= 2 * (BigInteger)uNum * uDen;
        return p;
    }

    /// <summary>The curve G(k) in its sin(Nk) form, (1 - u^2) sin(Nk) cos k + [(1 + u^2) cos(Nk) - 2u] sin k, real k
    /// (the proof names the three forms sum / sin(Nk) / sector rather than numbering them).</summary>
    public static double Curve(double k, int n, double u)
    {
        CheckRing(n);
        return (1 - u * u) * Math.Sin(n * k) * Math.Cos(k) + ((1 + u * u) * Math.Cos(n * k) - 2 * u) * Math.Sin(k);
    }

    /// <summary>uDen^2 * P(+2) as the closed factor form -((N-1)u + (N+1))(u - 1), scaled: -((N-1) uNum + (N+1) uDen)(uNum - uDen).
    /// U_n(1) = n + 1 is the whole derivation.</summary>
    public static BigInteger BandTopValueScaled(int n, long uNum, long uDen)
    {
        CheckRing(n);
        CheckU(uNum, uDen);
        return -((n - 1) * (BigInteger)uNum + (n + 1) * (BigInteger)uDen) * ((BigInteger)uNum - uDen);
    }

    /// <summary>uDen^2 * P(-2) as the closed factor form: ((N-1)u - (N+1))(u + 1) at odd N, -((N-1)u + (N+1))(u - 1)
    /// at even N; scaled likewise. U_n(-1) = (-1)^n (n + 1) is the whole derivation, the parity of N entering
    /// through that sign and nowhere else.</summary>
    public static BigInteger BandBottomValueScaled(int n, long uNum, long uDen)
    {
        CheckRing(n);
        CheckU(uNum, uDen);
        return n % 2 == 1
            ? ((n - 1) * (BigInteger)uNum - (n + 1) * (BigInteger)uDen) * ((BigInteger)uNum + uDen)
            : -((n - 1) * (BigInteger)uNum + (n + 1) * (BigInteger)uDen) * ((BigInteger)uNum - uDen);
    }

    /// <summary>Levels above the band, 0 or 1. P is monic and real-rooted, so P(2) &lt; 0 means an odd number of
    /// roots above 2, and Weyl's inequality for the rank-two update with eigenvalues (u, 0, ..., -u) makes that
    /// number 1. P(2) = 0 is an edge level, counted as not departed (see <see cref="OnTheEdge"/>).</summary>
    public static int DeparturesAbove(int n, long uNum, long uDen) =>
        BandTopValueScaled(n, uNum, uDen).Sign < 0 ? 1 : 0;

    /// <summary>Levels below the band, 0 or 1: the sign of P at -infinity is (-1)^N, so (-1)^N P(-2) &lt; 0 means an odd
    /// number of roots below -2, and Weyl's inequality makes it 1.</summary>
    public static int DeparturesBelow(int n, long uNum, long uDen)
    {
        var v = BandBottomValueScaled(n, uNum, uDen);
        int signAtMinusInfinity = n % 2 == 0 ? 1 : -1;
        return v.Sign * signAtMinusInfinity < 0 ? 1 : 0;
    }

    /// <summary>Whether a level sits exactly ON an edge of the band, P(+2) = 0 or P(-2) = 0, compared as zero
    /// over the integers. The top edge is met at u = 1 for every N; the bottom edge at u = 1 for even N and at
    /// u = (N+1)/(N-1) for odd N.</summary>
    public static (bool Top, bool Bottom) OnTheEdge(int n, long uNum, long uDen) =>
        (BandTopValueScaled(n, uNum, uDen).IsZero, BandBottomValueScaled(n, uNum, uDen).IsZero);

    /// <summary>The odd-N threshold u = (N+1)/(N-1) in lowest terms; for even N there is none (the bottom level
    /// leaves at every u &gt; 1) and the call throws rather than returning a number.</summary>
    public static (long Num, long Den) OddThreshold(int n)
    {
        CheckRing(n);
        if (n % 2 == 0) throw new ArgumentException("even N has no bottom threshold: the alternating level leaves at every u > 1", nameof(n));
        long g = Gcd(n + 1, n - 1);
        return ((n + 1) / g, (n - 1) / g);
    }

    /// <summary>c_m = 1/2 - 1/(N sin^2 k_m), k_m = 2 pi m/N, the coefficient of delta in
    /// Delta E_m/(4 delta J/N) = 1 + delta c_m + O(delta^2); 1 &lt;= m &lt; N/2 (a paired mode). Positive at the band
    /// centre, negative at the edge, zero at N sin^2 k_m = 2; the delta -&gt; 0 form and nothing else.</summary>
    public static double SplitCorrection(int n, int m)
    {
        CheckRing(n);
        if (m < 1 || 2 * m >= n) throw new ArgumentOutOfRangeException(nameof(m), m, "a paired mode has 1 <= m < N/2");
        double s = Math.Sin(2 * Math.PI * m / n);
        return 0.5 - 1.0 / (n * s * s);
    }

    /// <summary>dE_k/du at u = 0 for the chain level k = 1..N (E_k = 2 cos(k pi/(N+1)) there), in units of J:
    /// (-1)^(k+1) (4/(N+1)) sin^2(k pi/(N+1)). From the polynomial, -(dP/du)/(dP/dx) = 2/P'(x_k); from the
    /// eigenvectors, 2 psi_k(0) psi_k(N-1). Its magnitude is F65's alpha_k/gamma_0.</summary>
    public static double ChainEndVelocity(int n, int k)
    {
        CheckRing(n);
        if (k < 1 || k > n) throw new ArgumentOutOfRangeException(nameof(k), k, "the chain comb has k = 1..N");
        double s = Math.Sin(k * Math.PI / (n + 1));
        return (k % 2 == 1 ? 1.0 : -1.0) * 4.0 / (n + 1) * s * s;
    }

    private static void CheckRing(int n)
    {
        if (n < MinN) throw new ArgumentOutOfRangeException(nameof(n), n, "the cracked ring needs N >= 3");
    }

    private static void CheckU(long uNum, long uDen)
    {
        if (uDen <= 0) throw new ArgumentOutOfRangeException(nameof(uDen), uDen, "u = uNum/uDen needs uDen > 0");
        if (uNum < 0) throw new ArgumentOutOfRangeException(nameof(uNum), uNum, "u = J'/J is a non-negative coupling ratio here");
    }

    private static long Gcd(long a, long b)
    {
        while (b != 0) (a, b) = (b, a % b);
        return a;
    }

    public override string DisplayName =>
        "F160: the cracked ring is exactly solvable (det(x I - H) = U_N(x/2) - u^2 U_{N-2}(x/2) - 2u)";

    public override string Summary =>
        "one bond of the XY ring detuned to u*J: the whole single-excitation spectrum is the zero set of one " +
        "curve in k, the chain's transmission amplitude closed into a loop; u walks from the chain comb to the ring " +
        $"comb, and past u = 1 the departures are a parity law read off two linear factors ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the polynomial, and why it is exact",
                summary: "Expand det(x I - H) along the two corner entries -u: the path determinant U_N(x/2), minus " +
                         "u^2 times the path determinant of the N-2 inner sites, and the two cyclic terms, each " +
                         "(-1)^(N+1) times a product of N-1 entries -1 and one entry -u, which together give -2u. " +
                         "U_n(x/2) is monic with integer coefficients, so at rational u the identity lives over Q " +
                         "and is met coefficient for coefficient; MirrorWorld's Crack does that by Faddeev-LeVerrier " +
                         $"against the recursion. N = 3, u = 1: {Show(RoadPolynomialScaled(3, 1, 1))} (x^3 - 3x - 2 = " +
                         "(x - 2)(x + 1)^2, the triangle's 2, -1, -1).");

            yield return new InspectableNode("the curve is the polynomial at x = 2 cos k",
                summary: "U_n(cos k) = sin((n+1)k)/sin k turns P(2 cos k) into G(k)/sin k with G = sin((N+1)k) - " +
                         "u^2 sin((N-1)k) - 2u sin k, and the sum formulas for sin((N+-1)k) give the first printed " +
                         "form. G vanishes at k = 0 and pi for free, which is why the spectrum is the zero set on " +
                         "the OPEN interval and why a level can hide there only when P(+-2) = 0, the edge cases below. " +
                         "For u < 1 Perron-Frobenius (weighted row sums 2J except at the crack) keeps every level " +
                         "strictly inside. And the spectrum is SIMPLE at every u >= 0 with u != 1: G = 2AB, " +
                         "A = cos(k(N+1)/2) - u cos(k(N-1)/2), B = sin(k(N+1)/2) + u sin(k(N-1)/2), the two " +
                         "reflection sectors through the crack, each a nonvanishing prefactor times the " +
                         "characteristic polynomial of an unreduced Jacobi block (simple spectrum), and the " +
                         "Bezout identity (cos a + u cos b)A + (sin a - u sin b)B = 1 - u^2 forbids a common zero " +
                         "unless u^2 = 1. So for 0 <= u < 1 G has exactly N simple zeros on (0, pi); the experiment " +
                         "page had swept this over 1050 (N, delta) points (gate E2d(i)) and fenced it as unproven.");

            yield return new InspectableNode("the two ends, and what u is",
                summary: "u = 1: P = 2T_N(x/2) - 2, cos(Nk) = 1, the ring comb 2 pi m/N with its m <-> N-m pairs. " +
                         "u = 0: P = U_N(x/2), sin((N+1)k) = 0, the open chain comb pi m/(N+1), F2b, reached as the " +
                         "matching condition psi_{-1} = u psi_{N-1}, psi_N = u psi_0 with u set to zero. So u " +
                         "interpolates the moduli N and N+1 (the two turn-fraction families MirrorWorld's Cyclotomy " +
                         "owns) and is a boundary condition, not a graph: at every u > 0 the graph is a ring.");

            yield return new InspectableNode("the join with the walk-time step",
                summary: "t(k) = -2iu sin k/(e^(-ik) - u^2 e^(ik)) is the exact transmission amplitude of the same " +
                         "bond on the infinite chain (COUPLING_DEFECT_WALK_TIME_STEP, since 2026-07-12), and " +
                         "1/t = (1 + u^2)/(2u) + i(1 - u^2) cot k/(2u) is G's coefficient pair over 2u, so " +
                         "Re[e^(-iNk)/t(k)] - 1 = G(k)/(2u sin k). The upstream O(delta) reflection the sibling " +
                         "sets aside as 'not read as signal' is on the closed ring the whole signal. Two chains " +
                         "meet here and stay apart: the infinite one whose leads define t, and the open N-site " +
                         "one the same equation reaches at u = 0.");

            yield return new InspectableNode("the split, and its next order",
                summary: "Near k_m = 2 pi m/N, with sin(Nk) = sin x and cos(Nk) = cos x for x = N(k - k_m), the " +
                         "curve is delta(2 - delta) sin x cos(k_m + x/N) + sin(k_m + x/N)[2(1 - delta)(cos x - 1) + " +
                         "delta^2 cos x] = 0, whose small-(delta, x) truncation sin k_m x^2 - 2 delta cos k_m x - " +
                         "delta^2 sin k_m = 0 has roots x = delta(cos k_m +- 1)/sin k_m, the flat split 4 delta J/N. " +
                         "Carrying the series one order further (the proof, checked symbolically): " +
                         $"Delta E_m/(4 delta J/N) = 1 + delta c_m + O(delta^2), c_m = 1/2 - 1/(N sin^2 k_m); N = 12, m = 1: c = {SplitCorrection(12, 1):F6}, " +
                         $"m = 2: c = {SplitCorrection(12, 2):F6}. The sign of c_m changes at N sin^2 k_m = 2 (for m = 1 " +
                         "between N = 19 and 20). c_m is the delta -> 0 form; at a finite delta the split's own next " +
                         "order has moved that zero, and the (1,1) block's zero-crossing reading sits elsewhere again.");

            yield return new InspectableNode("the departures: two linear factors and a parity",
                summary: "U_n(+-1) = (+-1)^n (n + 1) gives P(+2) = -((N-1)u + (N+1))(u - 1) at every N and P(-2) = " +
                         "((N-1)u - (N+1))(u + 1) at odd N, -((N-1)u + (N+1))(u - 1) at even N. P is real-rooted " +
                         "(H real symmetric) and monic, so a negative P(2) means an odd number of roots above the " +
                         "band; the crack is the rank-two update with eigenvalues (u, 0, ..., -u), so Weyl's " +
                         "inequality lets at most one level pass each edge. Hence the count (u >= 0): top leaves at every u > 1; bottom " +
                         $"at every u > 1 for even N, past u = (N+1)/(N-1) for odd N. N = 7: threshold u = {OddThreshold(7).Num}/{OddThreshold(7).Den}, " +
                         $"departures at u = 5/4: {DeparturesAbove(7, 5, 4) + DeparturesBelow(7, 5, 4)}, at u = 4/3 exactly: " +
                         $"{DeparturesAbove(7, 4, 3) + DeparturesBelow(7, 4, 3)} with the bottom level on the edge " +
                         $"({OnTheEdge(7, 4, 3).Bottom}), at u = 3/2: {DeparturesAbove(7, 3, 2) + DeparturesBelow(7, 3, 2)}; " +
                         $"N = 8 at u = 101/100: {DeparturesAbove(8, 101, 100) + DeparturesBelow(8, 101, 100)}. The first " +
                         "draft of this count forgot the parity of N (CAUGHT_ERRORS 2026-08-31).");

            yield return new InspectableNode("the velocity at the chain end is F65's comb, signed",
                summary: "At u = 0 the roots are x_k = 2 cos(k pi/(N+1)) and dx_k/du = -(dP/du)/(dP/dx) = 2/P'(x_k) " +
                         "= (-1)^(k+1) (4/(N+1)) sin^2(k pi/(N+1)), because (d/dx) U_N(x/2) at x_k is " +
                         "-(N+1)(-1)^k/(2 sin^2 theta_k). The eigenvector route gives the same number as " +
                         "2 psi_k(0) psi_k(N-1), the chain reflection supplying the sign (F75's mirror sign). Its " +
                         "magnitude is exactly F65's alpha_k/gamma_0 = (4/(N+1)) sin^2(k pi/(N+1)), the rate of the " +
                         "level under one dephased endpoint, twice the Absorption Theorem's light there; velocity " +
                         "means dE/du along the road, not F2b's group velocity dE/dk. " +
                         $"N = 5, k = 1..5: {string.Join(", ", Enumerable.Range(1, 5).Select(k => ChainEndVelocity(5, k).ToString("F6")))} " +
                         "(F65's 1/6, 1/2, 2/3, 1/2, 1/6 with alternating sign). So F129's collisions on the chain comb " +
                         "move at first order by the signed F65 sums, and which survive is decided exactly there " +
                         "(experiments/THE_COMB_ON_THE_ROAD.md).");

            yield return new InspectableNode("scope, and the fences that do not lift",
                summary: "The XY adjacency book (H_se = J * adjacency), uniform bonds off the crack, u >= 0 rational " +
                         "for the exact statements; the roots themselves are a reading in MirrorWorld, held to the " +
                         "exact count. gamma does not appear in G: the Liouvillian rate of these coherences is the " +
                         "Absorption Theorem's, and the beat, the reversal and the visibility wall of the experiment " +
                         "are first-order statements this claim does not carry. Not a statement about the Heisenberg " +
                         "Laplacian, not about the blind seat, and 'departure' here is a level leaving the band, not " +
                         "the departure from normality of the F2b corollary and F89.");

            yield return new InspectableNode("where it is met from below",
                summary: "compute/MirrorWorld/Crack.cs (run mode crack N [u]; 74 CrackTests): the identity over the " +
                         "integers at twelve (N, u) and mod two primes to N = 1001, the Descartes count against the law " +
                         "over N = 3..20 with the odd threshold met at the exact rational, the 0.971754 cross-dock with " +
                         "WarbleTests, the c_m decade law. simulations/cracked_bell_gate.py stage E: the curve through " +
                         "an eigensolver (Newton step 1.8e-14 over 48 (N, delta) points, exactly N sign changes, the " +
                         "determinant identity to 4.1e-13, the join, the c_m law, the road past u = 1). This claim's " +
                         "own tests meet the polynomial against a Bareiss determinant over the integers and the " +
                         "factor forms against the unfactored polynomial at x = +-2, two routes each.");
        }
    }

    private static string Show(BigInteger[] ascending)
    {
        var parts = new List<string>();
        for (int i = ascending.Length - 1; i >= 0; i--)
        {
            if (ascending[i].IsZero) continue;
            string c = ascending[i].ToString();
            string coef = (c == "1" && i > 0) ? "" : (c == "-1" && i > 0) ? "-" : c + " ";
            parts.Add(i == 0 ? c : i == 1 ? $"{coef}x" : $"{coef}x^{i}");
        }
        return string.Join(" + ", parts).Replace("+ -", "- ");
    }
}
