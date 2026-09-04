using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Which of the two reflection sectors a halves-resultant belongs to. The names are the
/// reflection's, R-even and R-odd; they are not a parity of the knob and not a parity of N.</summary>
public enum ReflectionSector
{
    /// <summary>R-even. Its comb carries the ODD node indices.</summary>
    Even,

    /// <summary>R-odd. Its comb carries the EVEN node indices.</summary>
    Odd,
}

/// <summary>F162, the blind seat's sector factorisation: on the end-pair-anisotropic open chain
/// A + t·D, the two reflection sectors' halves-resultants factor F157's own locus polynomial
/// exactly, sign and leading coefficient included, and the sign reads the FOLD COORDINATE alone.
///
/// <code>
///     p = min(j, N−1−j)      the fold coordinate          n = |N − 1 − 2j|   the node modulus
///     S_m                    sin(mθ)/sin θ, monic in x = 2cos θ  (S_0 = 0, S_1 = 1, S_{−1} = −1)
///     α_p = S_{p+1} − t·S_p  the p-site chain carrying the knob at coordinate 0, monic of degree p
///
///     S_p · α_{N−1−p}  ≡  −S_n           (mod α_p)                            Lemma 9
///     Res(α_p, S_p)    =  (−1)^C(p,2)     with no t in it                      Lemma 10
///     Res(α_p, α_{N−1−p}) = (−1)^C(p+1,2) · Q_E · Q_O                          Corollary 10a
///     Q_S = c_S · ∏ (t − Δ_k) over the NON-POLE roots of β_S                   Corollary 10b
///     c_S = (−1)^(p·r_S + n_S) · Res(h_S, S_p) · Res(g_S, S_{p+1})
/// </code>
///
/// <para><b>Every resultant here is taken FOLD HALF FIRST</b>, outside and inside the sectors
/// alike. That is not tidiness: a resultant is antisymmetric up to (−1)^(deg f · deg g), so a sign
/// law is a statement about a NAMED argument order and about nothing else. Reading the halves in
/// seat order instead of fold order does move the verdict, on both sides: <see cref="OuterResultant(int, int, bool)"/>
/// and <see cref="ReadingHoldsFor(int, BigInteger[], bool)"/> take the order as an argument so a control can
/// feed the wrong one, and where it breaks is pinned there rather than counted in prose.</para>
///
/// <para><b>The mechanism, in one line.</b> α_p is monic in x, so reduction modulo it is a
/// division. Writing q = N−1−p, so that q − p = n, the Chebyshev addition formula at (n, p+1) and
/// at (n, p) turns α_q into S_n·α_{p+1} − S_{n−1}·α_p; the second summand dies mod
/// α_p, and multiplying the first by S_p and using t·S_p ≡ S_{p+1} leaves S_n·(S_p·S_{p+2} −
/// S_{p+1}²), which the CASSINI identity S_{p+1}² − S_p·S_{p+2} = 1 collapses to −S_n. The Cassini
/// step is F160's, cited and not repeated; the addition formula appears in no proof of this
/// repository and is written out in §(i).</para>
///
/// <para><b>The pole split, which is what makes the constant an integer.</b> The middle route the
/// seat leaves standing carries the modulus n, and the polynomial S_n whose roots are its nodes
/// 2cos(kπ/n), k = 1..n−1, which has degree n−1 and is the characteristic polynomial of an
/// (n−1)-site path, splits by the reflection into β_E·β_O. A node index k with S_p(x_k) = 0 is a POLE: Δ_k does not exist
/// there and that root contributes the constant S_{p+1}(x_k) instead of a linear factor. Splitting
/// β_S = g_S·h_S along g_S = gcd(β_S, S_p) makes this exact,
/// <c>Q_S = (−1)^(p·r_S)·Res(g_S, S_{p+1})·Res(h_S, α_p)</c>, so deg_t Q_S = n_S = deg h_S and
/// lc_t Q_S = c_S. The constant is nonzero because S_n is SQUAREFREE, which makes h_S coprime to
/// S_p; that hypothesis is load-bearing and not a convenience.</para>
///
/// <para><b>What this class recomputes, and on which road.</b> Everything below is exact integer
/// and rational arithmetic from (N, j) alone: no eigensolver, no float, no tolerance. The road is
/// deliberately NOT the gate's. The gate (<c>simulations/blind_seat_two_axes_proof.py</c>, block W)
/// builds every law-carrying resultant as a SYLVESTER DETERMINANT in sympy over ℤ[t][x], because
/// sympy's own routine is the object under suspicion there. This class instead
/// (1) computes resultants by the EUCLIDEAN remainder sequence over ℚ, never forming a Sylvester
/// matrix; (2) works one integer t at a time and recovers the ℤ[t] answer by Lagrange
/// interpolation, with two sample points held back so an underestimated t-degree reddens rather
/// than passing silently; and (3) decides the sector parity as an exact REMAINDER,
/// S_{n+1} + 1 ≡ 0 (mod β_E) and S_{n+1} − 1 ≡ 0 (mod β_O), where the gate reads it at five hundred
/// points of the node field. Three different roads to the same three statements.</para>
///
/// <para><b>What is READ and not derived here.</b> That the two sector LEFT halves are literally
/// α_p and the two right halves are the combs β_E, β_O is the gate's W5, read off the
/// reflection-adapted sector blocks; this class takes it as the definition of Q_S and does not
/// rebuild those blocks. That the multiplicity of a shared Δ in Q_E is §(h) Corollary 8a's b_E is a
/// composition of Lemma 8 with the parity, and is measured against no b_E by any gate, this one
/// included. F157's four committed rows carry a primitive-and-positive normalisation whose
/// discarded sign is exactly what this claim computes; the pinned per-row scales stay in the
/// gate.</para>
///
/// <para><b>What this class is NOT.</b> The fences are §(k)'s, restated for the words this file actually
/// uses; §(k) fences more words than appear here and is the fuller list. A <i>pair</i> is the site pair
/// {p, N−1−p} the reflection joins, never <c>MirrorWorld.Pair</c>'s bare coherence |i⟩⟨j|; the <i>end
/// pair</i> is F157's {0, N−1}, which the road page fences against F140's corner block. A <i>block</i>
/// here is the matrix of H restricted to one reflection sector, never <c>MirrorWorld.Block</c>'s
/// joint-popcount block (p, q); "block W" is the one place a block here is a lettered group of gates
/// instead. A <i>sector</i> is always a reflection sector, the whole object sitting inside one excitation
/// number, so no popcount grading is in play and neither
/// <c>RCPsiSquared.Core.BlockSpectrum.JointPopcountSectors</c> nor its witnesses are meant. A <i>comb</i>
/// is F157's node comb and, in the pole split, the middle route's own comb of x-values; never
/// <c>Cyclotomy</c>'s two combs, which are turn fractions and not Δ values. A <i>pole</i> is a k at which
/// S_p(x_k) = 0, so Δ_k is undefined and the degree drops; it carries no complex-analytic sense and is not
/// <c>Cat</c>'s two immortal poles. The <i>fold coordinate</i> is min(j, N−1−j), not <c>Mirror</c>'s folds
/// f_P and f_Q, which pay λ → −λ − 2Nγ; no spectrum is folded here. A <i>mirror seat</i> is a seat past
/// the reflection's fixed point and a <i>mirror image</i> a site's partner under R, never <c>Mirror</c>'s
/// block-lattice group of eight. Lemma 10's object is called THE COMMON FACTOR and never a divisor,
/// because F139 already uses "Chebyshev divisor" in the same polynomial-algebra sense; and the DIVISOR LAW
/// this file inherits from its parent is F157's gcd law on the SEAT INDEX, never <c>Divisor</c>'s frozen
/// divisor of F140 and never this assembly's own <c>FrozenDivisorClaim</c>, which live on the R₉₀ locus
/// and have no seat in them. No γ, no Liouvillian and no dephasing appears anywhere here: the whole
/// statement lives in the single-excitation XY block of a real symmetric chain.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md</c> §(i), gate block W of
/// <c>simulations/blind_seat_two_axes_proof.py</c>, registry entry F162.</para>
/// </summary>
public sealed class BlindSeatSectorFactorisationClaim : Claim
{
    /// <summary>Parent: F157, which owns the locus polynomial P_j, the comb Δ_k, the pole rule and
    /// the node modulus as DEFINITIONS, states no sign and no leading coefficient, and prints its
    /// worked rows under the normalisation that discards precisely the number computed here.</summary>
    public SeatCutBlindnessClaim BlindSeat { get; }

    /// <summary>Parent: F160, whose proof owns the Cassini identity S_{p+1}² − S_p·S_{p+2} = 1 that
    /// Lemma 9 turns on, and the U_m(x/2) normalisation the sine quotient is written in.</summary>
    public CrackedRingExactCurveClaim CrackedRing { get; }

    /// <summary>The smallest chain that has an interior seat the reflection does not fix.</summary>
    public const int MinChain = 4;

    /// <summary>The largest chain the gate sweeps. The statements are algebraic in p and N and hold
    /// past it; this bound is what has been READ, and the survey below stops here for that reason.</summary>
    public const int GateMaxChain = 14;

    public BlindSeatSectorFactorisationClaim(
        SeatCutBlindnessClaim blindSeat,
        CrackedRingExactCurveClaim crackedRing)
        : base("F162 the blind seat's sector factorisation: on the end-pair-anisotropic open chain A + t*D, with " +
               "S_m the sine quotient monic in x = 2cos(theta), alpha_p = S_{p+1} - t*S_p, p = min(j, N-1-j) the " +
               "fold coordinate and n = |N-1-2j| the node modulus, the congruence S_p*alpha_{N-1-p} = -S_n " +
               "(mod alpha_p) is the Chebyshev addition formula twice closed by the Cassini identity F160's proof " +
               "owns, and it carries the two reflection sectors' halves-resultants onto the seat's own locus " +
               "polynomial up to the single common factor Res(alpha_p, S_p) = (-1)^C(p,2), so that with every " +
               "resultant taken FOLD HALF FIRST, in the sectors as well as outside, " +
               "Res(alpha_p, alpha_{N-1-p}) = (-1)^C(p+1,2) * Q_E * Q_O, an exponent reading the fold coordinate " +
               "and NOT N; that object is F157's own definition route Res(S_n, Delta*S_j - S_{j+1}) times (-1)^e " +
               "with e = (n-1)(p+1) + p + C(p,2) + [j > N-1-j]*C(n,2), the last term the price of carrying F157's " +
               "SEAT index onto the fold coordinate through the node identity S_{p+n}(x_k) = (-1)^k*S_p(x_k); and " +
               "each factor is Q_S = c_S * prod (t - Delta_k) over the NON-POLE roots of beta_S, one factor per " +
               "non-pole root and a pole root contributing a constant instead, with " +
               "c_S = (-1)^(p*r_S + n_S)*Res(h_S, S_p)*Res(g_S, S_{p+1}) a nonzero integer off the pole split " +
               "g_S = gcd(beta_S, S_p), h_S = beta_S/g_S, so that deg_t Q_S = n_S and lc_t Q_S = c_S and the two " +
               "constants compose as lc_t(P_j) = (-1)^C(p+1,2)*c_E*c_O with deg_t P_j = n_E + n_O. Squarefreeness " +
               "of beta_S is a hypothesis and not a convenience. A sign law is a statement about a NAMED argument " +
               "order: a resultant is antisymmetric up to (-1)^(deg f * deg g), and the +-1 this closes had no " +
               "convention-free answer to give",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md (primary: section (i), Lemmas 9 to 11 and " +
               "Corollaries 10a, 10b, 11a, 11b; sections (g) and (h) are the two open constants it closes) + " +
               "docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md (F160's section (b), the Cassini identity Lemma 9 " +
               "turns on, cited and not repeated) + " +
               "docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md (F157's Lemma J, the node route the locus " +
               "polynomial is defined on) + " +
               "experiments/THE_BLIND_SEAT_ON_THE_ROAD.md (the road page whose open items (g) and (h) left) + " +
               "docs/ANALYTICAL_FORMULAS.md (F162; F157 the locus it factorises, F160 the Cassini step) + " +
               "docs/CAUGHT_ERRORS.md (2026-09-03, thirteen items: the repair of a check that cannot fail is " +
               "where the next one is written) + " +
               "simulations/blind_seat_two_axes_proof.py (the gate, block W, twenty-three checks)")
    {
        BlindSeat = blindSeat ?? throw new ArgumentNullException(nameof(blindSeat));
        CrackedRing = crackedRing ?? throw new ArgumentNullException(nameof(crackedRing));
    }

    // ---------------------------------------------------------------------------------------
    // the two integers a seat is, and the scope
    // ---------------------------------------------------------------------------------------

    /// <summary>p = min(j, N−1−j), the fold coordinate: the size of the shorter of the two pieces
    /// the seat cuts the chain into.</summary>
    public static int FoldCoordinate(int chain, int seat)
    {
        RequireSeat(chain, seat);
        return Math.Min(seat, chain - 1 - seat);
    }

    /// <summary>n = |N − 1 − 2j|, F157's node modulus: the difference in size between the two
    /// pieces, and the length of the middle route left standing between them.</summary>
    public static int NodeModulus(int chain, int seat)
    {
        RequireSeat(chain, seat);
        return Math.Abs(chain - 1 - 2 * seat);
    }

    /// <summary>Whether the seat is one the statement speaks about. It is NOT at the two end seats
    /// (p = 0: α_p is the empty determinant, both sides are 1, so the statement is true and says
    /// nothing) and NOT at the reflection-fixed centre seat of an odd chain (n = 0: the two halves
    /// share a characteristic polynomial, the outer resultant vanishes identically in t, and there
    /// is no ratio left to carry a sign).</summary>
    public static bool IsInScope(int chain, int seat) =>
        FoldCoordinate(chain, seat) >= 1 && NodeModulus(chain, seat) >= 1;

    /// <summary>The seats of one chain the statement speaks about, in increasing order.</summary>
    public static IReadOnlyList<int> SeatsInScope(int chain)
    {
        RequireChain(chain);
        var seats = new List<int>();
        for (int j = 0; j < chain; j++) if (IsInScope(chain, j)) seats.Add(j);
        return seats;
    }

    // ---------------------------------------------------------------------------------------
    // the polynomials: the sine quotient, the knob chain, the two sector combs
    // ---------------------------------------------------------------------------------------

    /// <summary>S_m as a monic integer polynomial in x = 2cos θ, ascending coefficients:
    /// S_0 = 0, S_1 = 1, S_{m+1} = x·S_m − S_{m−1}, and S_{−1} = −1 by the same recursion run back.
    /// This is F157's U_{m−1}(x/2), written in the book §(i) uses.</summary>
    public static BigInteger[] SineQuotient(int index)
    {
        if (index < -1) throw new ArgumentOutOfRangeException(nameof(index), index, "S_m is carried from m = -1 up.");
        if (index == -1) return new[] { BigInteger.MinusOne };
        var prev = new[] { BigInteger.Zero };                 // S_0
        if (index == 0) return prev;
        var cur = new[] { BigInteger.One };                   // S_1
        for (int m = 1; m < index; m++)
        {
            var next = new BigInteger[cur.Length + 1];
            for (int k = 0; k < cur.Length; k++) next[k + 1] += cur[k];   // x·S_m
            for (int k = 0; k < prev.Length; k++) next[k] -= prev[k];     // − S_{m−1}
            prev = cur;
            cur = Trim(next);
        }
        return cur;
    }

    /// <summary>α_p = S_{p+1} − t·S_p at one integer value of the knob, ascending in x. Monic of
    /// degree p, because S_{p+1} is monic of degree p and t·S_p has degree p−1.</summary>
    public static BigInteger[] KnobChainAt(int fold, BigInteger knob)
    {
        if (fold < 0) throw new ArgumentOutOfRangeException(nameof(fold), fold, "the fold coordinate is a size.");
        var top = SineQuotient(fold + 1);
        var low = SineQuotient(fold);
        var res = new BigInteger[Math.Max(top.Length, low.Length)];
        for (int k = 0; k < top.Length; k++) res[k] += top[k];
        for (int k = 0; k < low.Length; k++) res[k] -= knob * low[k];
        return Trim(res);
    }

    /// <summary>The two sector combs of the middle route, β_E·β_O = S_n. For even n they are
    /// S_n/S_{n/2} and S_{n/2}; for odd n with r = (n−1)/2 they are S_{r+1} − S_r and
    /// S_{r+1} + S_r. Both monic over ℤ.</summary>
    public static BigInteger[] SectorComb(int nodeModulus, ReflectionSector sector)
    {
        if (nodeModulus < 1)
            throw new ArgumentOutOfRangeException(nameof(nodeModulus), nodeModulus,
                "the node modulus is at least 1 wherever the statement speaks.");
        if (nodeModulus % 2 == 0)
        {
            var half = SineQuotient(nodeModulus / 2);
            if (sector == ReflectionSector.Odd) return half;
            var (quotient, remainder) = DivMod(SineQuotient(nodeModulus), half);
            if (!IsZero(remainder))
                throw new InvalidOperationException(
                    "S_{n/2} must divide S_n at even n; the recursion above is wrong rather than the arithmetic.");
            return quotient;
        }
        int r = (nodeModulus - 1) / 2;
        var up = SineQuotient(r + 1);
        var down = SineQuotient(r);
        return sector == ReflectionSector.Even ? Sub(up, down) : Add(up, down);
    }

    // ---------------------------------------------------------------------------------------
    // Lemma 9, the congruence
    // ---------------------------------------------------------------------------------------

    /// <summary>Lemma 9: S_p·α_{N−1−p} + S_n is divisible by α_p, exactly.
    ///
    /// <para>Checked at N + 3 integer knobs, which is enough to force the zero polynomial rather than to
    /// sample it: the probe is degree 1 in the knob, each of the N−1−p reduction steps by the knob-linear
    /// α_p raises the remainder's knob-degree by at most one, so the remainder's coefficients have degree at
    /// most N − p, and N + 3 &gt; N − p at every p ≥ 1. A nonzero remainder at ANY sample is a break.</para></summary>
    public static bool CongruenceIsExact(int chain, int seat)
    {
        RequireInScope(chain, seat);
        int p = FoldCoordinate(chain, seat), n = NodeModulus(chain, seat);
        var middle = SineQuotient(n);
        for (int knob = 0; knob <= chain + 2; knob++)
        {
            var alphaFold = KnobChainAt(p, knob);
            var alphaFar = KnobChainAt(chain - 1 - p, knob);
            var probe = Add(Mul(SineQuotient(p), alphaFar), middle);
            if (!IsZero(DivMod(probe, alphaFold).Remainder)) return false;
        }
        return true;
    }

    // ---------------------------------------------------------------------------------------
    // Lemma 10 and Corollary 10a, the common factor and the outer sign
    // ---------------------------------------------------------------------------------------

    /// <summary>C(m, 2) = m(m−1)/2, the binomial the two sign laws are written in.</summary>
    public static int Pairs(int m) => m * (m - 1) / 2;

    /// <summary>Lemma 10's common factor Res(α_p, S_p) as a polynomial in the knob. The content of
    /// the lemma is that it is CONSTANT, so this returns the whole polynomial and the caller reads
    /// both halves of the statement off it.</summary>
    public static BigInteger[] CommonFactor(int fold)
    {
        // deg S_p = p-1 rows of the Sylvester matrix carry the knob, so the knob-degree is at most p-1. The
        // bound is read off the matrix and not off Lemma 10, which is the thing under test here.
        if (fold < 1) throw new ArgumentOutOfRangeException(nameof(fold), fold, "Lemma 10 speaks from p = 1 up.");
        return InKnob(fold - 1, knob => Resultant(KnobChainAt(fold, knob), ToRational(SineQuotient(fold))));
    }

    /// <summary>Q_S = Res(α_p, β_S), the sector halves-resultant, FOLD HALF FIRST, as a polynomial
    /// in the knob.</summary>
    public static BigInteger[] SectorHalvesResultant(int chain, int seat, ReflectionSector sector)
    {
        RequireInScope(chain, seat);
        int p = FoldCoordinate(chain, seat);
        var comb = SectorComb(NodeModulus(chain, seat), sector);
        return HalvesResultant(p, comb);
    }

    /// <summary>Res(α_p, α_{N−1−p}), FOLD HALF FIRST, as a polynomial in the knob. This is the
    /// object the composition calls P_j.</summary>
    public static BigInteger[] OuterResultant(int chain, int seat) => OuterResultant(chain, seat, true);

    /// <summary>The same, with the argument ORDER handed in. <paramref name="foldHalfFirst"/> false reads the
    /// halves in seat order, the shorter one second at a mirror seat, and the parameter exists so a control
    /// can feed that: a resultant is antisymmetric up to (−1)^(deg f·deg g), so the sign law is a statement
    /// about this argument and about nothing else.</summary>
    public static BigInteger[] OuterResultant(int chain, int seat, bool foldHalfFirst)
    {
        RequireInScope(chain, seat);
        int p = FoldCoordinate(chain, seat), far = chain - 1 - p;
        if (!foldHalfFirst && seat > chain - 1 - seat) (p, far) = (far, p);
        // Both halves carry the knob, so all p + far = N - 1 rows of the Sylvester matrix can: the bound is
        // N - 1, read off the matrix and not off the composition, whose degree n_E + n_O is under test.
        return InKnob(chain - 1, knob => Resultant(KnobChainAt(p, knob), ToRational(KnobChainAt(far, knob))));
    }

    /// <summary>Corollary 10a: Res(α_p, α_{N−1−p}) = (−1)^C(p+1,2)·Q_E·Q_O, an exponent reading the
    /// fold coordinate and not N.</summary>
    public static bool OuterSignLawHolds(int chain, int seat)
    {
        int p = FoldCoordinate(chain, seat);
        var product = Mul(SectorHalvesResultant(chain, seat, ReflectionSector.Even),
                          SectorHalvesResultant(chain, seat, ReflectionSector.Odd));
        return Same(OuterResultant(chain, seat), Scale(product, Sign(Pairs(p + 1))));
    }

    // ---------------------------------------------------------------------------------------
    // Lemma 11 and the sector parity, decided as a remainder
    // ---------------------------------------------------------------------------------------

    /// <summary>The node identity S_{p+n} ≡ S_p·S_{n+1} (mod S_n), which is Lemma 11 before the
    /// value of S_{n+1} at a node is read: the addition formula at (p, n) with S_n ≡ 0.</summary>
    public static bool NodeIdentityIsExact(int nodeModulus, int fold)
    {
        if (nodeModulus < 1) throw new ArgumentOutOfRangeException(nameof(nodeModulus));
        if (fold < 0) throw new ArgumentOutOfRangeException(nameof(fold));
        var middle = SineQuotient(nodeModulus);
        if (middle.Length <= 1) return true;                  // n = 1: the middle route is a point
        var lhs = DivMod(SineQuotient(fold + nodeModulus), middle).Remainder;
        var rhs = DivMod(Mul(SineQuotient(fold), SineQuotient(nodeModulus + 1)), middle).Remainder;
        return Same(lhs, rhs);
    }

    /// <summary>The other half of Lemma 11, and §(h)'s Lemma 8 with it: S_{n+1}(x_k) = (−1)^k at
    /// every node, read as an exact REMAINDER rather than at points of the node field. β_E carries
    /// the odd indices exactly when S_{n+1} + 1 dies modulo it, β_O the even ones exactly when
    /// S_{n+1} − 1 does. Both are asserted here, so a swapped pair of combs breaks it.</summary>
    public static bool SectorParityIsExact(int nodeModulus) =>
        SectorParityIsExact(nodeModulus, ReflectionSector.Even);

    /// <summary>The same, with the ASSIGNMENT handed in: <paramref name="carriesOddIndices"/> names the
    /// sector claimed to carry the odd node indices. Lemma 8 says that is the R-even one, so a control can
    /// pass the R-odd one and require the reading to fail.</summary>
    public static bool SectorParityIsExact(int nodeModulus, ReflectionSector carriesOddIndices)
    {
        var top = SineQuotient(nodeModulus + 1);
        var odds = SectorComb(nodeModulus, carriesOddIndices);
        var evens = SectorComb(nodeModulus,
                               carriesOddIndices == ReflectionSector.Even
                                   ? ReflectionSector.Odd : ReflectionSector.Even);
        bool oddSide = odds.Length <= 1 || IsZero(DivMod(Add(top, new BigInteger[] { 1 }), odds).Remainder);
        bool evenSide = evens.Length <= 1 || IsZero(DivMod(Sub(top, new BigInteger[] { 1 }), evens).Remainder);
        return oddSide && evenSide;
    }

    // ---------------------------------------------------------------------------------------
    // Corollary 10b, the pole split and the constant
    // ---------------------------------------------------------------------------------------

    /// <summary>The pole split of one sector comb: g_S = gcd(β_S, S_p) monic, h_S = β_S/g_S, and
    /// the two degrees r_S = deg β_S, n_S = deg h_S the constant is written in.</summary>
    public static (BigInteger[] Poles, BigInteger[] Free, int CombDegree, int FreeDegree)
        PoleSplit(int chain, int seat, ReflectionSector sector)
    {
        RequireInScope(chain, seat);
        return PoleSplitOf(FoldCoordinate(chain, seat), SectorComb(NodeModulus(chain, seat), sector));
    }

    /// <summary>The same split for ANY monic integer comb. Corollary 10b's identity is stated for a monic
    /// SQUAREFREE β, and that hypothesis is load-bearing rather than a convenience, so the general entry
    /// exists in order that a control can feed a comb which is not squarefree and require the reading to
    /// break. Nothing in the sector path calls it with anything but β_E or β_O.</summary>
    public static (BigInteger[] Poles, BigInteger[] Free, int CombDegree, int FreeDegree)
        PoleSplitOf(int fold, BigInteger[] comb)
    {
        RequireComb(fold, comb);
        var poles = MonicGcd(comb, SineQuotient(fold));
        var free = ExactQuotient(comb, poles);
        return (poles, free, comb.Length - 1, free.Length - 1);
    }

    /// <summary>Res(α_p, β) fold half first, for any monic integer comb, as a polynomial in the knob.</summary>
    public static BigInteger[] HalvesResultant(int fold, BigInteger[] comb) => HalvesResultant(fold, comb, true);

    /// <summary>The same, with the argument ORDER handed in, so a control can read the halves right half
    /// first and require the law to break.</summary>
    public static BigInteger[] HalvesResultant(int fold, BigInteger[] comb, bool foldHalfFirst)
    {
        RequireComb(fold, comb);
        // Only the deg β rows built from alpha_p carry the knob, so the knob-degree is at most deg β. The
        // bound is the Sylvester row count and not Corollary 10b's n_S, which is what is under test.
        var target = ToRational(comb);
        return foldHalfFirst
            ? InKnob(comb.Length - 1, knob => Resultant(KnobChainAt(fold, knob), target))
            : InKnob(comb.Length - 1, knob => Resultant(target, ToRational(KnobChainAt(fold, knob))));
    }

    /// <summary>Corollary 10b read on any monic integer comb: whether the pole split reproduces the
    /// halves-resultant, and whether its knob-degree and leading knob-coefficient are n_S and c_S. Both are
    /// theorems for a SQUAREFREE comb only.</summary>
    public static (bool ProductHolds, bool DegreeAndLeadHold) ReadingHoldsFor(int fold, BigInteger[] comb) =>
        ReadingHoldsFor(fold, comb, true);

    /// <inheritdoc cref="ReadingHoldsFor(int, BigInteger[])"/>
    public static (bool ProductHolds, bool DegreeAndLeadHold)
        ReadingHoldsFor(int fold, BigInteger[] comb, bool foldHalfFirst)
    {
        var (poles, free, combDegree, freeDegree) = PoleSplitOf(fold, comb);
        var halves = HalvesResultant(fold, comb, foldHalfFirst);

        var polePart = ResultantOfIntegerPolynomials(poles, SineQuotient(fold + 1));
        var freePart = InKnob(free.Length - 1,
                              knob => Resultant(ToRational(free), ToRational(KnobChainAt(fold, knob))));
        bool product = Same(halves, Scale(freePart, Sign(fold * combDegree) * polePart));

        var constant = Sign(fold * combDegree + freeDegree)
                       * ResultantOfIntegerPolynomials(free, SineQuotient(fold))
                       * polePart;
        bool reading = halves.Length - 1 == freeDegree && halves[^1] == constant;
        return (product, reading);
    }

    private static void RequireComb(int fold, BigInteger[] comb)
    {
        if (fold < 1) throw new ArgumentOutOfRangeException(nameof(fold), fold, "the fold coordinate is at least 1.");
        if (comb is null) throw new ArgumentNullException(nameof(comb));
        if (comb.Length == 0 || comb[^1] != BigInteger.One)
            throw new ArgumentOutOfRangeException(nameof(comb), "the comb must be monic over the integers.");
    }

    /// <summary>Corollary 10b's constant c_S = (−1)^(p·r_S + n_S)·Res(h_S, S_p)·Res(g_S, S_{p+1}),
    /// a nonzero integer.</summary>
    public static BigInteger SectorConstant(int chain, int seat, ReflectionSector sector)
    {
        int p = FoldCoordinate(chain, seat);
        var (poles, free, combDegree, freeDegree) = PoleSplit(chain, seat, sector);
        var withFold = ResultantOfIntegerPolynomials(free, SineQuotient(p));
        var withNext = ResultantOfIntegerPolynomials(poles, SineQuotient(p + 1));
        return Sign(p * combDegree + freeDegree) * withFold * withNext;
    }

    /// <summary>Corollary 10b's product line: a POLE root contributes a constant and not a linear factor,
    /// so Q_S = (−1)^(p·r_S)·Res(g_S, S_{p+1})·Res(h_S, α_p) exactly.
    ///
    /// <para>What it reads, stated narrowly because the wide reading would be wrong: dividing out the
    /// identity Res(α_p, β_S) = (−1)^(p·r_S)·Res(g_S, α_p)·Res(h_S, α_p), which is multiplicativity and
    /// holds whatever the split is, what is left is Res(g_S, α_p) = Res(g_S, S_{p+1}), that is, g_S | S_p.
    /// Where g_S = 1 that is 1 = 1, and g_S = 1 at 124 of the 144 sector readings over N = 4..14, the
    /// complement being exactly the 20 pole shortfalls the survey counts. So this line carries the pole
    /// rule at 20 readings and the argument ORDER everywhere. The falsifiable reading of the constant is
    /// <see cref="DegreeAndLeadingCoefficientHold"/>, where Q_S arrives by Euclid and interpolation while
    /// n_S and c_S come off the pole split, two roads that meet rather than one.</para></summary>
    public static bool PoleSplitIdentityHolds(int chain, int seat, ReflectionSector sector)
    {
        RequireInScope(chain, seat);
        return ReadingHoldsFor(FoldCoordinate(chain, seat),
                               SectorComb(NodeModulus(chain, seat), sector)).ProductHolds;
    }

    /// <summary>Corollary 10b's two readings on Q_S itself: its knob-degree is n_S and its leading
    /// knob-coefficient is c_S. One factor per non-pole root, and a pole drops the degree.</summary>
    public static bool DegreeAndLeadingCoefficientHold(int chain, int seat, ReflectionSector sector)
    {
        RequireInScope(chain, seat);
        return ReadingHoldsFor(FoldCoordinate(chain, seat),
                               SectorComb(NodeModulus(chain, seat), sector)).DegreeAndLeadHold;
    }

    // ---------------------------------------------------------------------------------------
    // Corollary 11b and the composition
    // ---------------------------------------------------------------------------------------

    /// <summary>Corollary 11b's exponent e = (n−1)(p+1) + p + C(p,2) + [j &gt; N−1−j]·C(n,2). The
    /// last term is the price of carrying F157's SEAT index onto the fold coordinate, and it is
    /// paid only at the mirror seats.</summary>
    public static int GeneratorTieExponent(int chain, int seat)
    {
        RequireInScope(chain, seat);
        int p = FoldCoordinate(chain, seat), n = NodeModulus(chain, seat);
        int mirrored = seat > chain - 1 - seat ? Pairs(n) : 0;
        return (n - 1) * (p + 1) + p + Pairs(p) + mirrored;
    }

    /// <summary>F157's own definition route G_j = Res(S_n, Δ·S_j − S_{j+1}) = Res(S_n, −α_j), as a
    /// polynomial in the knob. Note the SEAT index j, not the fold coordinate.</summary>
    public static BigInteger[] SeatGenerator(int chain, int seat)
    {
        RequireInScope(chain, seat);
        int n = NodeModulus(chain, seat);
        var middle = ToRational(SineQuotient(n));
        // Only the deg S_n = n - 1 rows built from alpha_j carry the knob.
        return InKnob(n - 1, knob => Resultant(middle, ToRational(Negate(KnobChainAt(seat, knob)))));
    }

    /// <summary>Corollary 11b: Res(α_p, α_{N−1−p}) = (−1)^e·G_j, the fold-coordinate object and
    /// F157's seat-indexed generator differing by one closed-form sign.</summary>
    public static bool GeneratorTieHolds(int chain, int seat) =>
        Same(OuterResultant(chain, seat),
             Scale(SeatGenerator(chain, seat), Sign(GeneratorTieExponent(chain, seat))));

    /// <summary>The composition of the two constants: deg_t P_j = n_E + n_O and
    /// lc_t P_j = (−1)^C(p+1,2)·c_E·c_O, where P_j is the outer resultant.</summary>
    public static bool CompositionHolds(int chain, int seat)
    {
        int p = FoldCoordinate(chain, seat);
        var (_, _, _, freeEven) = PoleSplit(chain, seat, ReflectionSector.Even);
        var (_, _, _, freeOdd) = PoleSplit(chain, seat, ReflectionSector.Odd);
        var outer = OuterResultant(chain, seat);
        var lead = Sign(Pairs(p + 1))
                   * SectorConstant(chain, seat, ReflectionSector.Even)
                   * SectorConstant(chain, seat, ReflectionSector.Odd);
        return outer.Length - 1 == freeEven + freeOdd && outer[^1] == lead;
    }

    // ---------------------------------------------------------------------------------------
    // the survey, recomputed
    // ---------------------------------------------------------------------------------------

    /// <summary>What one sweep over the seats in scope found. Every field is a count of seats, and
    /// the last three are populations rather than verdicts: <paramref name="SeatsByFold"/> so a
    /// narrowed sweep reddens rather than quietly thinning, and <paramref name="RepeatedFactorSeats"/>
    /// so the per-non-pole-root wording has something that could contradict "every factor simple".</summary>
    public sealed record FactorisationSurvey(
        int MaxChain,
        int Seats,
        int CongruenceHolds,
        int CommonFactorHolds,
        int OuterSignLawHolds,
        int GeneratorTieHolds,
        int PoleSplitHolds,
        int DegreeAndLeadHolds,
        int CompositionHolds,
        int NodeIdentityHolds,
        int SectorParityHolds,
        int PoleShortfallReadings,
        IReadOnlyDictionary<int, int> SeatsByFold,
        IReadOnlyList<(int Chain, int Seat, ReflectionSector Sector, int Multiplicity)> RepeatedFactorSeats,
        int MirrorSeats);

    /// <summary>Runs every statement of §(i) this class carries over the seats in scope of chains
    /// 4..<paramref name="maxChain"/>, the two node-level ones, Lemma 11 and the sector parity, included at
    /// each seat's own (n, p) and counted there even where the modulus makes them EMPTY rather than true:
    /// at n = 1 both are vacuous, which is 12 of the 72 seats over N = 4..14, and at n = 2 one half of the
    /// parity is. Nothing is stored; the counts below are recomputed at call time. What it does
    /// NOT run is what the class does not carry: Corollary 11a's seat-versus-fold reading of Δ_k, §(i)'s
    /// convention reading of K1B_SIGNS, and the resultant-routine checks W0 to W0d, all of which stay in the
    /// gate.</summary>
    public static FactorisationSurvey Survey(int maxChain = GateMaxChain)
    {
        RequireChain(maxChain);
        int seats = 0, congruence = 0, common = 0, outer = 0, tie = 0, split = 0, degree = 0, composed = 0;
        int mirror = 0, shortfall = 0, nodeIdentity = 0, parity = 0;
        var byFold = new SortedDictionary<int, int>();
        var repeated = new List<(int, int, ReflectionSector, int)>();

        for (int chain = MinChain; chain <= maxChain; chain++)
        foreach (int seat in SeatsInScope(chain))
        {
            seats++;
            int p = FoldCoordinate(chain, seat);
            byFold[p] = byFold.TryGetValue(p, out int had) ? had + 1 : 1;
            if (seat > chain - 1 - seat) mirror++;

            if (CongruenceIsExact(chain, seat)) congruence++;
            var factor = CommonFactor(p);
            if (factor.Length == 1 && factor[0] == Sign(Pairs(p))) common++;
            if (OuterSignLawHolds(chain, seat)) outer++;
            if (GeneratorTieHolds(chain, seat)) tie++;
            if (CompositionHolds(chain, seat)) composed++;
            if (NodeIdentityIsExact(NodeModulus(chain, seat), p)) nodeIdentity++;
            if (SectorParityIsExact(NodeModulus(chain, seat))) parity++;

            bool splitOk = true, degreeOk = true;
            foreach (var sector in new[] { ReflectionSector.Even, ReflectionSector.Odd })
            {
                splitOk &= PoleSplitIdentityHolds(chain, seat, sector);
                degreeOk &= DegreeAndLeadingCoefficientHold(chain, seat, sector);
                var (_, _, combDegree, freeDegree) = PoleSplit(chain, seat, sector);
                if (freeDegree < combDegree) shortfall++;
                int multiplicity = LargestRepeatedFactorMultiplicity(chain, seat, sector);
                if (multiplicity > 1) repeated.Add((chain, seat, sector, multiplicity));
            }
            if (splitOk) split++;
            if (degreeOk) degree++;
        }

        return new FactorisationSurvey(maxChain, seats, congruence, common, outer, tie, split, degree, composed,
                                       nodeIdentity, parity, shortfall, byFold, repeated, mirror);
    }

    /// <summary>The largest multiplicity with which one factor repeats in Q_S, read off the gcd with the
    /// derivative and never off a root. A constant Q_S has no factor at all and reports 0; 1 means every
    /// factor is simple THERE, which is a reading and not the statement, since §(i) says one factor per
    /// non-pole root and a shared Δ repeats.
    ///
    /// <para>The gcd here is made PRIMITIVE over ℤ and not monic. Q_S is monic in x but not in the knob:
    /// its leading knob-coefficient is c_S, which is 4 at N = 19 seat 3, where Q_E = 4t⁶ − 12t⁴ + 9t² =
    /// t²·(2t² − 3)². A monic normalisation would ask for t³ − (3/2)t there and leave the integers, and the
    /// exception it raised would have named the construction rather than the normalisation.</para></summary>
    public static int LargestRepeatedFactorMultiplicity(int chain, int seat, ReflectionSector sector)
    {
        var q = SectorHalvesResultant(chain, seat, sector);
        if (q.Length <= 2) return q.Length - 1 <= 0 ? 0 : 1;
        int multiplicity = 1;
        var repeated = PrimitiveGcd(q, Derivative(q));
        while (repeated.Length > 1)
        {
            multiplicity++;
            repeated = PrimitiveGcd(repeated, Derivative(repeated));
        }
        return multiplicity;
    }

    // ---------------------------------------------------------------------------------------
    // exact arithmetic: the Euclidean resultant, the knob sweep, integer polynomial helpers
    // ---------------------------------------------------------------------------------------

    /// <summary>Res(f, g) by the EUCLIDEAN remainder sequence, never by a Sylvester matrix:
    /// Res(f, g) = (−1)^(deg f · deg g)·lc(g)^(deg f − deg r)·Res(g, r) with r = f mod g. The
    /// argument ORDER is the statement's; nothing here symmetrises it.</summary>
    public static BigRational Resultant(BigInteger[] first, RationalPolynomial second) =>
        Resultant(ToRational(first), second);

    /// <inheritdoc cref="Resultant(BigInteger[], RationalPolynomial)"/>
    public static BigRational Resultant(RationalPolynomial first, RationalPolynomial second)
    {
        if (first.IsZero || second.IsZero) return BigRational.Zero;
        int df = first.Degree, dg = second.Degree;
        if (df == 0 && dg == 0) return BigRational.One;
        if (dg == 0) return Power(second[0], df);
        if (df == 0) return Power(first[0], dg);

        var rest = first.Mod(second);
        if (rest.IsZero) return BigRational.Zero;
        var sign = (df * dg) % 2 == 0 ? BigRational.One : -BigRational.One;
        return sign * Power(second[dg], df - rest.Degree) * Resultant(second, rest);
    }

    /// <summary>Res(f, g) for two integer polynomials, as the integer it must be.</summary>
    public static BigInteger ResultantOfIntegerPolynomials(BigInteger[] first, BigInteger[] second) =>
        RequireInteger(Resultant(ToRational(first), ToRational(second)),
                       "a resultant of integer polynomials");

    /// <summary>Evaluates <paramref name="atKnob"/> at the integer knobs 0..<paramref name="degreeBound"/>+2 and
    /// returns the integer polynomial in the knob through them. TWO points more than the
    /// interpolation needs are held back and checked against the interpolant, so an underestimated
    /// knob-degree throws here instead of passing as a shorter answer, PROVIDED the underestimate is by one
    /// or two: the difference between the true polynomial and the interpolant already vanishes at the
    /// bound+1 sample points, so from bound+3 upward it can absorb both held-back points as well. Every
    /// caller passes a bound read off the Sylvester matrix's row count and never off a statement under test,
    /// and every one of them is a valid upper bound, so the window is not being relied on; it is stated
    /// because the method is public, so that a control can hand it a bound that is too small and require the
    /// throw.</summary>
    public static BigInteger[] InKnob(int degreeBound, Func<BigInteger, BigRational> atKnob)
    {
        int samples = degreeBound + 1;
        var values = new BigRational[samples + 2];
        for (int k = 0; k < values.Length; k++) values[k] = atKnob(k);

        var interpolant = Lagrange(values.Take(samples).ToArray());
        for (int k = samples; k < values.Length; k++)
            if (Evaluate(interpolant, k) != values[k])
                throw new InvalidOperationException(
                    $"the knob-degree bound {degreeBound} is too small: the interpolant through {samples} points " +
                    $"misses the held-back sample at t = {k}. The bound above is wrong, not the arithmetic.");

        var coefficients = new BigInteger[Math.Max(interpolant.Degree + 1, 1)];
        for (int k = 0; k < coefficients.Length; k++)
            coefficients[k] = RequireInteger(interpolant[k], "an interpolated resultant coefficient");
        return Trim(coefficients);
    }

    private static RationalPolynomial Lagrange(BigRational[] values)
    {
        var acc = RationalPolynomial.Zero;
        for (int i = 0; i < values.Length; i++)
        {
            if (values[i].IsZero) continue;
            var term = new RationalPolynomial(BigRational.One);
            var weight = BigRational.One;
            for (int m = 0; m < values.Length; m++)
            {
                if (m == i) continue;
                term *= RationalPolynomial.LinearFactor(new BigRational(m));
                weight *= new BigRational(i - m);
            }
            acc += (values[i] / weight) * term;
        }
        return acc;
    }

    private static BigRational Evaluate(RationalPolynomial poly, BigInteger at)
    {
        var acc = BigRational.Zero;
        for (int k = poly.Degree; k >= 0; k--) acc = acc * new BigRational(at) + poly[k];
        return acc;
    }

    private static BigRational Power(BigRational value, int exponent)
    {
        if (exponent < 0) throw new ArgumentOutOfRangeException(nameof(exponent));
        var acc = BigRational.One;
        for (int k = 0; k < exponent; k++) acc *= value;
        return acc;
    }

    private static BigInteger RequireInteger(BigRational value, string what)
    {
        if (!value.IsInteger)
            throw new InvalidOperationException($"{what} came out non-integral ({value}), which cannot happen " +
                                                "for monic integer arguments; the construction above is wrong.");
        return value.Numerator;
    }

    private static RationalPolynomial ToRational(BigInteger[] poly) =>
        new(poly.Select(c => new BigRational(c)));

    /// <summary>gcd over ℚ, returned monic over ℤ (which it is, both arguments being monic).</summary>
    private static BigInteger[] MonicGcd(BigInteger[] first, BigInteger[] second)
    {
        var a = ToRational(first);
        var b = ToRational(second);
        while (!b.IsZero) (a, b) = (b, a.Mod(b));
        if (a.IsZero) return new BigInteger[] { 0 };
        var lead = a[a.Degree];
        var monic = new BigInteger[a.Degree + 1];
        for (int k = 0; k <= a.Degree; k++)
            monic[k] = RequireInteger(a[k] / lead, "a monic gcd coefficient");
        return monic;
    }

    /// <summary>gcd over ℚ, returned PRIMITIVE over ℤ: the coefficients are cleared to integers and divided
    /// by their content, with the sign of the leading coefficient made positive. Unlike
    /// <see cref="MonicGcd"/> this asks nothing of the arguments beyond integrality.</summary>
    private static BigInteger[] PrimitiveGcd(BigInteger[] first, BigInteger[] second)
    {
        var a = ToRational(first);
        var b = ToRational(second);
        while (!b.IsZero) (a, b) = (b, a.Mod(b));
        if (a.IsZero) return new BigInteger[] { 0 };

        var lcm = a.DenominatorLcm();
        var cleared = new BigInteger[a.Degree + 1];
        for (int k = 0; k <= a.Degree; k++)
        {
            var scaled = a[k] * new BigRational(lcm);
            cleared[k] = RequireInteger(scaled, "a cleared gcd coefficient");
        }
        var content = BigInteger.Zero;
        foreach (var c in cleared) content = BigInteger.GreatestCommonDivisor(content, c);
        if (content.IsZero) return new BigInteger[] { 0 };
        if (cleared[^1].Sign < 0) content = -content;
        return Trim(cleared.Select(c => c / content).ToArray());
    }

    private static BigInteger[] ExactQuotient(BigInteger[] numerator, BigInteger[] denominator)
    {
        var (quotient, remainder) = DivMod(numerator, denominator);
        if (!IsZero(remainder))
            throw new InvalidOperationException("the pole split must divide the comb exactly; it did not.");
        return quotient;
    }

    private static (BigInteger[] Quotient, BigInteger[] Remainder) DivMod(BigInteger[] numerator, BigInteger[] denominator)
    {
        var (q, r) = ToRational(numerator).DivMod(ToRational(denominator));
        return (FromRational(q), FromRational(r));
    }

    private static BigInteger[] FromRational(RationalPolynomial poly)
    {
        if (poly.IsZero) return new BigInteger[] { 0 };
        var coefficients = new BigInteger[poly.Degree + 1];
        for (int k = 0; k <= poly.Degree; k++)
            coefficients[k] = RequireInteger(poly[k], "a quotient or remainder coefficient");
        return coefficients;
    }

    private static BigInteger[] Add(BigInteger[] a, BigInteger[] b) => Combine(a, b, 1);
    private static BigInteger[] Sub(BigInteger[] a, BigInteger[] b) => Combine(a, b, -1);

    private static BigInteger[] Combine(BigInteger[] a, BigInteger[] b, int sign)
    {
        var res = new BigInteger[Math.Max(a.Length, b.Length)];
        for (int k = 0; k < a.Length; k++) res[k] += a[k];
        for (int k = 0; k < b.Length; k++) res[k] += sign * b[k];
        return Trim(res);
    }

    private static BigInteger[] Mul(BigInteger[] a, BigInteger[] b)
    {
        var res = new BigInteger[a.Length + b.Length - 1];
        for (int i = 0; i < a.Length; i++)
        for (int k = 0; k < b.Length; k++)
            res[i + k] += a[i] * b[k];
        return Trim(res);
    }

    private static BigInteger[] Negate(BigInteger[] a) => a.Select(c => -c).ToArray();

    private static BigInteger[] Scale(BigInteger[] a, BigInteger by) => Trim(a.Select(c => c * by).ToArray());

    private static BigInteger[] Derivative(BigInteger[] a)
    {
        if (a.Length <= 1) return new BigInteger[] { 0 };
        var res = new BigInteger[a.Length - 1];
        for (int k = 1; k < a.Length; k++) res[k - 1] = k * a[k];
        return Trim(res);
    }

    private static BigInteger[] Trim(BigInteger[] a)
    {
        int last = a.Length - 1;
        while (last > 0 && a[last].IsZero) last--;
        return a.Take(last + 1).ToArray();
    }

    private static bool IsZero(BigInteger[] a) => a.Length == 1 && a[0].IsZero;

    private static bool Same(BigInteger[] a, BigInteger[] b) => Trim(a).SequenceEqual(Trim(b));

    private static BigInteger Sign(int exponent) => exponent % 2 == 0 ? BigInteger.One : BigInteger.MinusOne;

    private static void RequireChain(int chain)
    {
        if (chain < MinChain)
            throw new ArgumentOutOfRangeException(nameof(chain), chain,
                $"the statement needs an interior seat the reflection does not fix, so N >= {MinChain}.");
    }

    private static void RequireSeat(int chain, int seat)
    {
        RequireChain(chain);
        if (seat < 0 || seat >= chain)
            throw new ArgumentOutOfRangeException(nameof(seat), seat, "the seat is a site of the chain.");
    }

    private static void RequireInScope(int chain, int seat)
    {
        if (!IsInScope(chain, seat))
            throw new ArgumentOutOfRangeException(nameof(seat), seat,
                "the statement is empty at an end seat (p = 0) and fenced out at the reflection-fixed centre " +
                "seat of an odd chain (n = 0), where the outer resultant vanishes identically.");
    }

    // ---------------------------------------------------------------------------------------
    // the claim's own face
    // ---------------------------------------------------------------------------------------

    public override string DisplayName =>
        "F162: the blind seat's two sector halves-resultants factor F157's locus, and the sign reads the fold";

    public override string Summary =>
        "with every resultant taken FOLD HALF FIRST, Res(alpha_p, alpha_{N-1-p}) = (-1)^C(p+1,2)*Q_E*Q_O at " +
        "p = min(j, N-1-j), an exponent reading the fold coordinate and not N; each factor is " +
        "c_S * prod (t - Delta_k) over the NON-POLE roots of its sector comb, a pole contributing a constant " +
        $"instead, and the two constants compose into the locus polynomial's own leading coefficient ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode(
                "the two integers a seat is",
                summary: "p = min(j, N-1-j) is the fold coordinate and n = |N-1-2j| the node modulus, with " +
                         $"N-1 = 2p + n. At N = 11 seat 2: p = {FoldCoordinate(11, 2)}, n = {NodeModulus(11, 2)}; " +
                         $"at its mirror seat 8: p = {FoldCoordinate(11, 8)}, n = {NodeModulus(11, 8)}, the same " +
                         "pair, which is why the seat index has a price to pay to reach the fold coordinate",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                "Lemma 9, the congruence, recomputed as a division",
                summary: "S_p*alpha_{N-1-p} + S_n is divisible by the monic alpha_p, exactly. " +
                         $"At N = 11 seat 2: {CongruenceIsExact(11, 2)}; at N = 14 seat 3: {CongruenceIsExact(14, 3)}. " +
                         "The Chebyshev addition formula twice, closed by F160's Cassini identity",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                "Lemma 10, the common factor, which carries no knob",
                summary: $"Res(alpha_p, S_p) = (-1)^C(p,2): at p = 2 it is {Format(CommonFactor(2))}, " +
                         $"at p = 3 it is {Format(CommonFactor(3))}, at p = 5 it is {Format(CommonFactor(5))}. " +
                         "A polynomial is returned rather than a number so that BOTH halves of the lemma, the " +
                         "value and the absence of t, can break",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                "Corollary 10a, the outer sign law",
                summary: $"Res(alpha_p, alpha_{{N-1-p}}) = (-1)^C(p+1,2)*Q_E*Q_O, fold half first. " +
                         $"N = 11 seat 2: {OuterSignLawHolds(11, 2)}; N = 12 seat 4: {OuterSignLawHolds(12, 4)}; " +
                         $"N = 13 seat 9 (a mirror seat): {OuterSignLawHolds(13, 9)}. The exponent reads p alone",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                "Lemma 11 and the sector parity, decided as a remainder",
                summary: $"S_{{p+n}} = S_p*S_{{n+1}} (mod S_n) at (n, p) = (5, 3): {NodeIdentityIsExact(5, 3)}; " +
                         $"and S_{{n+1}} + 1 dies mod beta_E while S_{{n+1}} - 1 dies mod beta_O at n = 5: " +
                         $"{SectorParityIsExact(5)}, at n = 8: {SectorParityIsExact(8)}. That is beta_E carrying " +
                         "the ODD node indices, read exactly rather than at points of the node field",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                "Corollary 10b, the pole split and its constant",
                summary: $"at N = 10 seat 2, sector E: comb degree {PoleSplit(10, 2, ReflectionSector.Even).CombDegree}, " +
                         $"non-pole degree {PoleSplit(10, 2, ReflectionSector.Even).FreeDegree}, " +
                         $"c_E = {SectorConstant(10, 2, ReflectionSector.Even)}, and Q_E is " +
                         $"{Format(SectorHalvesResultant(10, 2, ReflectionSector.Even))}, a repeated factor from two " +
                         "odd node indices and not a coincidence of the product",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                "Corollary 11b, the tie to F157's own generator",
                summary: $"Res(alpha_p, alpha_{{N-1-p}}) = (-1)^e * Res(S_n, Delta*S_j - S_{{j+1}}) with " +
                         "e = (n-1)(p+1) + p + C(p,2) + [mirror seat]*C(n,2). " +
                         $"N = 9 seat 2 (e = {GeneratorTieExponent(9, 2)}): {GeneratorTieHolds(9, 2)}; " +
                         $"N = 9 seat 6, its mirror (e = {GeneratorTieExponent(9, 6)}): {GeneratorTieHolds(9, 6)}",
                provenance: NodeProvenance.Live);

            yield return new InspectableNode("the sweep, recomputed here", summary: Describe(Survey()),
                                             provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                "scope, and what is read rather than derived",
                summary: "the end-pair-anisotropic open chain A + t*D in the single-excitation XY block, every " +
                         "N >= 4 and every interior seat the reflection does not fix. EMPTY at the end seats " +
                         "(p = 0, both sides are 1) and OUT at the odd chain's centre seat (n = 0, the outer " +
                         "resultant vanishes identically and there is no ratio to carry a sign). Says nothing " +
                         "about the ring, where the same identity stops being a constant at 20 of 32 seats over " +
                         "N = 4..10, a different population from the 20 pole shortfalls above (the proof's gate " +
                         "K1b3, 2026-09-04, says what stands there and why those 12 stand, the divisibility " +
                         "(N-1-2jr) | jr, no claim of its own yet), and " +
                         "nothing about F157's hop-2 normalisation, so MirrorWorld's un-normalised 128 at N = 9 " +
                         "seat 1 is not identified here. READ and not derived: that the sector left halves ARE " +
                         "alpha_p and the right halves the two combs (the gate's W5); and that the multiplicity " +
                         "of a shared Delta is section (h) Corollary 8a's b_E, which no gate measures against a b_E");
        }
    }

    private static string Describe(FactorisationSurvey s) =>
        $"chains 4..{s.MaxChain}, {s.Seats} seats in scope ({s.MirrorSeats} of them mirror seats), by fold " +
        $"coordinate {string.Join(", ", s.SeatsByFold.Select(kv => $"{kv.Key}:{kv.Value}"))}; " +
        $"congruence {s.CongruenceHolds}, common factor {s.CommonFactorHolds}, outer sign law {s.OuterSignLawHolds}, " +
        $"generator tie {s.GeneratorTieHolds}, pole split {s.PoleSplitHolds}, degree and leading coefficient " +
        $"{s.DegreeAndLeadHolds}, composition {s.CompositionHolds}, node identity {s.NodeIdentityHolds}, " +
        $"sector parity {s.SectorParityHolds}; a pole drops the knob-degree at " +
        $"{s.PoleShortfallReadings} of the {2 * s.Seats} sector readings; " +
        (s.RepeatedFactorSeats.Count == 0
            ? "no repeated factor anywhere in the range"
            : $"{s.RepeatedFactorSeats.Count} sector readings carry a repeated factor, the first being " +
              $"N = {s.RepeatedFactorSeats[0].Chain} seat {s.RepeatedFactorSeats[0].Seat} " +
              $"sector {s.RepeatedFactorSeats[0].Sector} at multiplicity {s.RepeatedFactorSeats[0].Multiplicity}");

    private static string Format(BigInteger[] poly)
    {
        if (poly.Length == 1) return poly[0].ToString();
        var parts = new List<string>();
        for (int k = poly.Length - 1; k >= 0; k--)
        {
            if (poly[k].IsZero) continue;
            var magnitude = BigInteger.Abs(poly[k]);
            var body = k == 0 ? magnitude.ToString()
                     : magnitude.IsOne ? (k == 1 ? "t" : $"t^{k}")
                     : k == 1 ? $"{magnitude}*t" : $"{magnitude}*t^{k}";
            parts.Add(parts.Count == 0
                          ? (poly[k].Sign < 0 ? "-" + body : body)
                          : (poly[k].Sign < 0 ? " - " + body : " + " + body));
        }
        return parts.Count == 0 ? "0" : string.Concat(parts);
    }
}
