using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the F157 anisotropy locus: which Δ leave a seat blind (claim
/// <see cref="SeatCutBlindnessClaim"/>, sibling witness <see cref="SeatCutBlindnessWitness"/>).
///
/// <para><b>The question this object exists for.</b> F157's two closed forms are stated at two
/// POINTS of one axis: <c>(gcd(2j+1, N) − 1)/2</c> at the isotropic point Δ = 1 and
/// <c>gcd(j+1, N+1) − 1</c> at the XY point Δ = 0, where the chain carries
/// <c>H = Σ_b J_b (XX + YY + Δ·ZZ)</c>. Nothing between them was ever built. This witness carries
/// the knob and reports the whole locus.</para>
///
/// <code>
///     N_node = |N − 1 − 2j|                                  the single controlling integer
///     N_node = 0     ->  blind at EVERY Δ, blind = (N−1)/2   the reflection-fixed centre seat
///     N_node >= 1    ->  blind exactly at Δ_k, k = 1..N_node−1   (the locus, below)
///     j = 0 or N−1   ->  blind = 0 at every Δ                an end seat, one half is empty
///
///     Δ_k = sin((j+1)kπ/N_node) / sin(jkπ/N_node)
/// </code>
///
/// <para><b>Why there is no radical and no float in this file.</b> Writing x = cos θ: the hop is
/// 2J on every bond, so the single-excitation block is 2A + diagonal and its EIGENVALUE is
/// 4cos θ, four times the Chebyshev variable, not twice it (measured: the N = 9, Δ = 0 spectrum
/// is exactly 4cos(kπ/10)). With that, the nodes θ_k = kπ/N_node are exactly the roots of U_{N_node−1},
/// and Δ_k = U_j(cos θ_k)/U_{j−1}(cos θ_k). So the locus is the real root set of
///
/// <code>
///     P_j(Δ)  =  Res_x( U_{N_node−1}(x),  Δ·U_{j−1}(x) − U_j(x) )      in ℤ[Δ]
/// </code>
///
/// an INTEGER polynomial. <see cref="LocusPolynomial"/> returns its primitive integer coefficients,
/// computed by exact Sylvester determinants over <see cref="BigInteger"/> (fraction-free Bareiss)
/// at N_node integer points followed by an exact interpolation. A trig closed form is therefore reported
/// as an integer identity: the irrational members of the locus (√3 at N = 9 seat 1, √2/2 at N = 9
/// seat 2, √2 and √(2±√2) at N = 11 seat 1, 2√3/3 at N = 11 seat 2) are NAMED by a polynomial
/// rather than printed as doubles. That list is an illustration and not the content: the content is
/// the polynomial, which is why the list can be checked at all.</para>
///
/// <para><b>What is recomputed against what, and what each route can actually certify.</b> THREE
/// routes meet here and they do not carry the same weight, so they are named separately.
/// (1) <see cref="LocusPolynomial"/>, Chebyshev arithmetic, which PRESUPPOSES that the shared
/// eigenvalues sit at the roots of U_{N_node−1}. (2) <see cref="DefinitionPolynomial(int)"/>,
/// Res_λ(χ_L, χ_R) built from the two submatrices by the Jacobi three-term recursion with no
/// Chebyshev identity used, a second DERIVATION rather than a second implementation; its end shift
/// is a parameter precisely so a gate can feed a wrong chain and require the two to part.
/// (3) the Krylov rank, over ℚ in <see cref="BlindAtRational(int, long, long)"/> and over the
/// whole field in <see cref="LocusOverFp"/>.
///
/// The rank over ℚ decides only RATIONAL Δ, and the locus meets ℚ at 0 and ±1 alone, which are
/// the two committed endpoints and the sign partner, so <see cref="RationalRootAgreements"/>
/// certifies NOTHING this witness is the first to name. Saying it closed a loop was an overclaim,
/// made here before it was measured. <see cref="LocusOverFp"/> is what closes it: over GF(p) every
/// Δ is rational, the irrational members' images included, and sweeping every residue makes the
/// comparison two-sided by construction.</para>
///
/// <para><b>Falsifiers are printed beside the counts.</b> Every locus is reported beside a Δ that
/// is NOT in it and is shown to give blind = 0 by the same routine; the end seats are reported as
/// the zeros they must be; and the ONE input that cannot break the law is named so it is never
/// counted as a probe: Δ = −1, because ΣH(Δ)Σ = −H(−Δ) with Σ the staggering makes blind(−Δ) =
/// blind(Δ) identically (F152's bipartite cospectrality reason). Prime N is NOT a second one, and
/// the arc note that said so was wrong: it is a weak probe of the ZZ book alone, where
/// gcd(2j+1, N) forces every seat but the centre to zero, and at the XY endpoint prime N = 11
/// blinds seven seats.</para>
///
/// <para><b>Scope, and it is the honest part.</b> The repo already owns the trigonometric node
/// route at BOTH endpoints, exactly and past the wall (<c>simulations/seat_cut_blindness.py</c>
/// "Route 3, EXACT" to N = 200, <c>simulations/blind_site.py</c>'s <c>mode_set</c>). What is new
/// here is the general FORM at every (N, j), and N_node as the single controlling integer; the word
/// "interpolation" is deliberately not used for it, both because this object is not between the
/// endpoints (see the neither-endpoint seats below) and because the word is already spoken for in
/// this file by the Lagrange interpolation that builds the polynomials. The MECHANISM was
/// Δ-general before this: Lemma J (J4) of PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md gives
/// blind = deg gcd(χ_L, χ_R) for any unreduced Jacobi matrix, and its §(f) records that the ZZ term
/// "changes only the diagonal, which Lemma J never constrains". The locus is derived from the two
/// principal submatrices, which is F157's FENCED phrasing: it needs a chain with no zero bond, and
/// this witness runs on the UNIFORM chain only, where that fence cannot bite.</para>
///
/// <para>Children: the per-seat locus table (live), the polynomial route vs the count route (live),
/// the two committed endpoints (live), the falsifiers (live), and the scope fences. Root:
/// <c>inspect --root blindlocus</c>.</para></summary>
public sealed class SeatBlindnessDeltaLocusWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Two primes for the rank. A rank over GF(p) can only DROP relative to the rank over
    /// ℚ, so the LARGER rank is kept and the reported count can only ever be too LARGE.
    ///
    /// <para><b>Which failure mode this actually defends against, which is NOT the sibling's.</b>
    /// <see cref="SeatCutBlindnessWitness"/> takes two primes because a coupling divisible by a
    /// ranking prime would fake blindness. Here that cannot happen at all: <see cref="MaxDeltaTerm"/>
    /// holds every matrix entry strictly below both primes, so no entry is a nonzero multiple of
    /// one. What remains is the residual mode, a Krylov MINOR that happens to vanish mod p while
    /// being nonzero over ℚ. A minor is a determinant and is not bounded by the entries, so this is
    /// not excluded by any guard; two primes make it need a simultaneous accident. It is rare rather
    /// than impossible, and the direction is pinned by
    /// <see cref="BlindAtRationalModP"/> so that taking the smaller rank could not pass unnoticed.
    /// </para></summary>
    private static readonly long[] PrimeSet = { 2147483647L, 999999937L };

    /// <summary>The ranking primes, exposed so a gate can measure the entry bound against the REAL
    /// array instead of a hand copy of it, and so the one-sidedness POLICY can be tested: see
    /// <see cref="BlindAtRational(int, long, long, long[])"/>.</summary>
    public static IReadOnlyList<long> Primes => PrimeSet;

    /// <summary>The largest chain the COUNT side will run on. An N × N elimination, so the guard is
    /// generous: the congruence law is checked past the 2^N wall.</summary>
    public const int MaxN = 200;

    /// <summary>The largest chain the LOCUS POLYNOMIAL will be built on. The polynomial costs N_node
    /// Sylvester determinants of size about (N_node + j) over <see cref="BigInteger"/>, so it is capped
    /// far below <see cref="MaxN"/>. Cost, not physics: the congruence counts hold to
    /// <see cref="MaxN"/> whether or not the polynomial is built.</summary>
    public const int MaxLocusN = 24;

    /// <summary>The largest |numerator| or |denominator| a rational Δ may carry, and the number is
    /// pinned by an exact bound rather than chosen for comfort.
    ///
    /// <para>Every entry of the cleared matrix den·H(Δ) is bounded by max(2·den, (N−1)·|num|), which
    /// at <see cref="MaxN"/> = 200 and this cap is max(2²¹, 199·2²⁰) = 208 666 624, strictly BELOW
    /// both entries of <see cref="Primes"/>. So no entry can be a nonzero multiple of a ranking
    /// prime, and the sibling witness's failure mode, a coupling divisible by a prime reporting the
    /// whole space blind, is structurally impossible here rather than merely defended against. That
    /// is what this constant buys; int64 headroom is not the binding constraint and would permit a
    /// cap 2²⁵ times larger.</para></summary>
    public const long MaxDeltaTerm = 1L << 20;

    /// <summary>The largest prime <see cref="LocusOverFp"/> will sweep. The sweep costs p Krylov
    /// eliminations, so this is a COST bound and not physics; it also keeps the arithmetic honest,
    /// since <c>KrylovRankModP</c> multiplies two reduced residues and 2003² is 4.01·10⁶, ten orders
    /// of magnitude below int64, while the ranking primes near 2³¹ leave only a factor of two.</summary>
    public const long MaxSweepPrime = 2003;

    /// <summary>The largest prime any rank route will accept. <c>KrylovRankModP</c> forms
    /// <c>Mod(h, p) · vec</c> with both factors below p, so the product must stay inside int64:
    /// p ≤ 3 037 000 499 keeps p² ≤ 9.22·10¹⁸. Both entries of <see cref="Primes"/> sit below it,
    /// the larger by a factor of about 1.41. Without this the two public rank entries would take
    /// any long and wrap in silence.</summary>
    public const long MaxRankingPrime = 3037000499L;

    public int N { get; }

    public SeatBlindnessDeltaLocusWitness(int n = 9)
    {
        if (n < 2)
            throw new ArgumentOutOfRangeException(nameof(n), $"a chain needs at least 2 sites; got N = {n}.");
        if (n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"N = {n} exceeds the live-count guard {MaxN}; the count is an N x N elimination and the guard " +
                "is cost, not physics.");
        N = n;
    }

    // -----------------------------------------------------------------------------------------
    // the controlling integer, and the two congruence counts that live on it
    // -----------------------------------------------------------------------------------------

    /// <summary>N_node = |N − 1 − 2j|, the single integer the whole locus hangs on. N_node = 0 exactly at the
    /// reflection-fixed centre seat of an odd chain.</summary>
    public int NodeModulus(int seat)
    {
        RequireSeat(seat);
        return Math.Abs(N - 1 - 2 * seat);
    }

    /// <summary>Whether the seat is blind at EVERY Δ, which happens exactly at N_node = 0. The two
    /// principal submatrices are then conjugate by the chain reflection, so they carry the same
    /// characteristic polynomial and share every root; the resultant vanishes identically. This is
    /// a proof, not a sample.</summary>
    public bool BlindAtEveryDelta(int seat) => NodeModulus(seat) == 0;

    /// <summary>Whether this seat's locus is a nonempty FINITE set, i.e. an interior seat with
    /// N_node ≥ 1 that is not the forced centre and whose locus does not collapse to poles.
    ///
    /// <para>The criterion is <b>empty ⇔ N_node divides j</b>, so nonempty is <c>j % nodeModulus != 0</c>. The
    /// first version of <see cref="Summary"/> wrote <c>nodeModulus % j != 0</c>, the reverse divisibility,
    /// which is wrong at every N and silently drops seat 1 always, since nodeModulus % 1 == 0. It reported 4
    /// interior seats at N = 9 where there are 6, and the two it dropped were seats 1 and 2, the
    /// pair the registry entry writes out as its worked examples. Nothing tested
    /// <see cref="Summary"/>, so the tree below it printed the truth while the headline did not.
    /// </para></summary>
    public bool HasNonemptyLocus(int seat)
    {
        RequireSeat(seat);
        if (seat == 0 || seat == N - 1) return false;
        int nodeModulus = NodeModulus(seat);
        return nodeModulus >= 1 && seat % nodeModulus != 0;
    }

    /// <summary>The count the closed form predicts at Δ = 1, as an integer congruence:
    /// #{k ∈ 1..N_node−1 : (2j+1)k ≡ N_node (mod 2·N_node)}. Null at an end seat and at the forced centre, where
    /// the k-count is not the governing form. Compare against
    /// <see cref="SeatCutBlindnessClaim.BlindHeisenberg"/>.</summary>
    public int? CountAtIsotropic(int seat)
    {
        RequireSeat(seat);
        int nodeModulus = NodeModulus(seat);
        if (nodeModulus == 0 || seat == 0 || seat == N - 1) return null;
        int c = 0;
        for (int k = 1; k < nodeModulus; k++)
            if ((2L * seat + 1) * k % (2L * nodeModulus) == nodeModulus % (2L * nodeModulus)) c++;
        return c;
    }

    /// <summary>The count the closed form predicts at Δ = 0, as an integer congruence:
    /// #{k ∈ 1..N_node−1 : N_node | (j+1)k}. Compare against <see cref="SeatCutBlindnessClaim.BlindXy"/>.
    ///
    /// <para>The pole exclusion N_node ∤ jk is IMPLIED and is kept below only to make the derivation
    /// legible: if N_node divided both (j+1)k and jk it would divide their difference k, which is
    /// impossible for 0 &lt; k &lt; N_node. Measured: over N_node = 1..399 and j = 0..399 the second clause
    /// rejects 0 of 462 504 qualifying k. It is not doing work, and no sentence anywhere should
    /// present it as if it were.</para></summary>
    public int? CountAtXy(int seat)
    {
        RequireSeat(seat);
        int nodeModulus = NodeModulus(seat);
        if (nodeModulus == 0 || seat == 0 || seat == N - 1) return null;
        int c = 0;
        for (int k = 1; k < nodeModulus; k++)
            if ((seat + 1L) * k % nodeModulus == 0 && (long)seat * k % nodeModulus != 0) c++;
        return c;
    }

    // -----------------------------------------------------------------------------------------
    // the locus polynomial, in exact integer arithmetic
    // -----------------------------------------------------------------------------------------

    /// <summary>The primitive integer coefficients of P_j(Δ) = Res_x(U_{N_node−1}, Δ·U_{j−1} − U_j),
    /// ASCENDING (index = power of Δ), whose real roots are exactly the Δ that leave the seat blind.
    ///
    /// <para>Returns <c>{1}</c> (the constant polynomial, no roots) when the locus is empty, which
    /// happens exactly when N_node divides j: every node k is then a pole of the ratio. Returns null at
    /// the forced centre, where the locus is the whole axis rather than a root set, and at an end
    /// seat, where it is empty for a different reason (one half of the chain is not there).</para>
    ///
    /// <para>P may carry repeated roots; the LOCUS is its set of DISTINCT real roots. No squarefree
    /// part is taken here because it would not change the set and would add a polynomial gcd this
    /// file does not otherwise need.</para></summary>
    public BigInteger[]? LocusPolynomial(int seat)
    {
        RequireSeat(seat);
        if (N > MaxLocusN)
            throw new InvalidOperationException(
                $"N = {N} exceeds the locus-polynomial guard {MaxLocusN}; the polynomial costs N_node Sylvester " +
                "determinants over BigInteger. The congruence counts are unaffected and still run to " +
                $"{MaxN}. The guard is cost, not physics.");
        int nodeModulus = NodeModulus(seat);
        if (nodeModulus == 0 || seat == 0 || seat == N - 1) return null;
        if (nodeModulus < 2) return new[] { BigInteger.One };

        BigInteger[] a = ChebyshevU(nodeModulus - 1);
        BigInteger[] uPrev = ChebyshevU(seat - 1);
        BigInteger[] uThis = ChebyshevU(seat);

        // deg_Δ P <= nodeModulus-1, because Res is homogeneous of degree deg(U_{nodeModulus-1}) = nodeModulus-1 in the
        // coefficients of the second argument and those are linear in Δ. So nodeModulus points determine it.
        int pts = nodeModulus;
        var ys = new BigInteger[pts];
        for (int t = 0; t < pts; t++)
        {
            BigInteger[] b = SubScaled(uPrev, t, uThis);   // t*U_{j-1} - U_j
            ys[t] = ResultantInt(a, b);
        }
        return Primitive(InterpolateAtZeroToN(ys));
    }

    /// <summary>Every rational root of an integer polynomial, exactly, by the rational-root
    /// theorem: p/q in lowest terms with p | a₀ and q | aₙ, each candidate decided by an exact
    /// integer evaluation. Not a search and not a sample. Separated from
    /// <see cref="RationalLocusPoints"/> so that the enumeration can be pinned on polynomials whose
    /// roots are chosen, rather than only on the ones this object happens to produce: over the whole
    /// guarded range the locus's rational points are only ever 0 and ±1, so nothing here that
    /// handles a denominator would otherwise ever run.</summary>
    public static IReadOnlyList<(long Num, long Den)> RationalRootsOf(BigInteger[] poly)
    {
        var found = new List<(long, long)>();
        BigInteger[] p = Trim(poly);
        if (p.Length < 2 || p.All(c => c.IsZero)) return found;

        int shift = 0;
        while (shift < p.Length && p[shift].IsZero) shift++;
        if (shift > 0)
        {
            found.Add((0L, 1L));
            p = p.Skip(shift).ToArray();          // a0 and an are BOTH re-read from here on
        }
        if (p.Length < 2) return found;

        BigInteger a0 = p[0], an = p[^1];
        foreach (BigInteger num in Divisors(BigInteger.Abs(a0)))
        foreach (BigInteger den in Divisors(BigInteger.Abs(an)))
        {
            if (BigInteger.GreatestCommonDivisor(num, den) != BigInteger.One) continue;
            foreach (int sign in new[] { 1, -1 })
            {
                BigInteger n2 = sign * num;
                if (!EvaluatesToZero(p, n2, den)) continue;
                if (BigInteger.Abs(n2) > long.MaxValue || den > long.MaxValue)
                    throw new InvalidOperationException(
                        $"a rational root {n2}/{den} does not fit in a (long, long) pair; the caller's " +
                        "polynomial is far outside anything this witness builds.");
                var pair = ((long)n2, (long)den);
                if (!found.Contains(pair)) found.Add(pair);
            }
        }
        // Ordered EXACTLY, by cross-multiplication: an exact route exists, so no double appears.
        found.Sort((x, y) => ((BigInteger)x.Item1 * y.Item2).CompareTo((BigInteger)y.Item1 * x.Item2));
        return found;
    }

    /// <summary>The locus polynomial reached the OTHER way: Res_λ(χ_L, χ_R) over ℤ[Δ], formed from
    /// the two principal submatrices the definition names, with no Chebyshev identity used
    /// anywhere. Primitive integer coefficients, ascending; null exactly where
    /// <see cref="LocusPolynomial"/> is.
    ///
    /// <para><b>Why this exists, and it is the difference between a cross-check and a claim.</b>
    /// <see cref="LocusPolynomial"/> builds Res_x(U_{N_node−1}, Δ·U_{j−1} − U_j), which already ENCODES
    /// the answer: it presupposes that the shared eigenvalues sit at the roots of U_{N_node−1}, which is
    /// the very thing being claimed. Agreeing with a second implementation of that same formula
    /// would only pin the arithmetic. This method instead forms the object F157's own criterion
    /// names, the resultant of the two submatrices' characteristic polynomials in λ, by the
    /// three-term Jacobi recursion and the same Bareiss/interpolation machinery. The two agreeing
    /// is the closed form being CORRECT, not merely computed twice.</para>
    ///
    /// <para>Both sides are compared as PRIMITIVE integer polynomials, which is the right
    /// comparison: the two resultants are formed from different matrices and differ by a nonzero
    /// rational scale, never in their roots.</para></summary>
    public BigInteger[]? DefinitionPolynomial(int seat) => DefinitionPolynomial(seat, N - 3);

    /// <summary>The definition route with its END SHIFT exposed. The physical value is N − 3, which
    /// is what <see cref="DefinitionPolynomial(int)"/> supplies: the ZZ term gives every interior
    /// site Δ(N−5) and the two chain ends Δ(N−3), a defect of 2Δ at the outer end of each half.
    ///
    /// <para><b>Why the parameter exists.</b> It is the seam a mutation needs. The claim that this
    /// route and <see cref="LocusPolynomial"/> are two derivations rather than one routine cannot be
    /// tested by comparing their outputs at different seats, which is what an earlier gate here did:
    /// replacing this method's body with <c>return LocusPolynomial(seat);</c> passed that gate
    /// unchanged. With the shift as a parameter a test can feed a physically WRONG chain and require
    /// the two routes to part, which the disguised version cannot do. The shift is not a knob on the
    /// physics and no caller should pass anything but N − 3.</para></summary>
    public BigInteger[]? DefinitionPolynomial(int seat, int endShift)
    {
        RequireSeat(seat);
        if (N > MaxLocusN)
            throw new InvalidOperationException(
                $"N = {N} exceeds the locus-polynomial guard {MaxLocusN}; the guard is cost, not physics.");
        int nodeModulus = NodeModulus(seat);
        if (nodeModulus == 0 || seat == 0 || seat == N - 1) return null;

        int left = seat, right = N - 1 - seat;
        // deg_Δ of the Sylvester determinant is bounded by right*left + left*right; loose is safe,
        // an over-sampled interpolation recovers the same polynomial and Trim drops the zero head.
        int pts = 2 * left * right + 1;
        var ys = new BigInteger[pts];
        for (int t = 0; t < pts; t++)
            ys[t] = ResultantInt(CharPolyOfPath(left, t, false, endShift), CharPolyOfPath(right, t, true, endShift));
        return Primitive(InterpolateAtZeroToN(ys));
    }

    /// <summary>The characteristic polynomial, ascending in λ, of the principal submatrix on one
    /// side of the seat at integer Δ = t, by the three-term recursion for a tridiagonal matrix:
    /// p_i = (x − a_i)·p_{i−1} − b²·p_{i−2}, with b = 2 the hop.
    ///
    /// <para>The diagonal is NOT constant, and that is the step the derivation turns on: an end site
    /// of the WHOLE chain touches one bond and pays t·(N−3), an interior site touches two and pays
    /// t·(N−5). Each submatrix inherits exactly one such end, the outer one, never the cut end.
    /// <paramref name="outer"/> says which end of this block that is, and it is kept because it makes
    /// the construction legible, NOT because the answer depends on it: every off-diagonal here is the
    /// same hop 2, so the block is conjugate to its own reversal by the flip permutation and the
    /// characteristic polynomial cannot see which end carries the defect. Swapping the branch changes
    /// no output and no gate could distinguish it, which is worth saying rather than leaving a reader
    /// to assume a parameter this visible must be load-bearing. What IS load-bearing in this routine
    /// is the interior shift t·(N−5) and the hop, and both are caught by the two-route agreement.</para></summary>
    private BigInteger[] CharPolyOfPath(int size, long t, bool outer, int endShift)
    {
        var diag = new BigInteger[size];
        for (int i = 0; i < size; i++) diag[i] = t * (N - 5);
        if (size > 0) diag[outer ? size - 1 : 0] = t * endShift;

        // p_{-1} = 1, p_0 = 1; here p is built in ascending-λ coefficient arrays.
        var prev2 = new[] { BigInteger.Zero };
        var prev1 = new[] { BigInteger.One };
        for (int i = 0; i < size; i++)
        {
            var cur = new BigInteger[i + 2];
            for (int k = 0; k < prev1.Length; k++)
            {
                cur[k + 1] += prev1[k];              // x * p_{i-1}
                cur[k] -= diag[i] * prev1[k];        // - a_i * p_{i-1}
            }
            for (int k = 0; k < prev2.Length; k++)
                cur[k] -= 4 * prev2[k];              // - b^2 * p_{i-2}, b = 2
            prev2 = prev1;
            prev1 = cur;
        }
        return Trim(prev1);
    }

    /// <summary>The RATIONAL members of this seat's locus. Empty at an end seat and at the forced
    /// centre, where the locus is not a root set.
    ///
    /// <para>Measured over the whole guarded range N = 4..<see cref="MaxLocusN"/>, every interior
    /// seat: the rational points are only ever <b>0 and ±1</b>, never a proper fraction. Δ = 0 and
    /// Δ = 1 are the two committed endpoints and Δ = −1 is the sign partner forced by
    /// Σ H(Δ) Σ = −H(−Δ), so the locus meets ℚ exactly where F157 was already standing and is
    /// irrational everywhere else. That is a measurement about this family, not an assumption:
    /// <see cref="RationalRootsOf"/> is fully general and is pinned as such.</para></summary>
    public IReadOnlyList<(long Num, long Den)> RationalLocusPoints(int seat)
    {
        BigInteger[]? p = LocusPolynomial(seat);
        return p is null ? Array.Empty<(long, long)>() : RationalRootsOf(p);
    }

    // -----------------------------------------------------------------------------------------
    // the independent COUNT route: an exact rank at a rational Δ
    // -----------------------------------------------------------------------------------------

    /// <summary>blind(seat) at the rational Δ = num/den, recomputed from the DEFINITION: N minus the
    /// rank of the seat's Krylov matrix.
    ///
    /// <para>Clearing the denominator is exact and costs nothing: den·H(Δ) has integer entries and a
    /// Krylov rank is invariant under a nonzero scaling of H, so this is the same count the sibling
    /// witness takes, with the same two primes and the same one-sidedness (a rank mod p can only
    /// drop, so the count can only ever be too LARGE).</para></summary>
    public int BlindAtRational(int seat, long num, long den) =>
        BlindAtRational(seat, num, den, PrimeSet);

    /// <summary>The count over a GIVEN set of ranking primes. The shipped set is
    /// <see cref="Primes"/> and no caller should pass another one; the parameter exists so a gate
    /// can reach the ONE-SIDEDNESS POLICY. A rank mod p can only be too small, so a nullity mod p
    /// can only be too large, so the tighter answer is the LARGER rank, and this method keeps it.
    /// With the shipped primes that choice never does any work, since the two agree on every row
    /// measured; feed a set containing a deliberately bad prime and the two policies part.</summary>
    public int BlindAtRational(int seat, long num, long den, long[] primes)
    {
        RequireSeat(seat);
        if (primes is null || primes.Length == 0)
            throw new ArgumentException("at least one ranking prime is needed.", nameof(primes));
        if (den == 0) throw new ArgumentOutOfRangeException(nameof(den), "a rational needs a nonzero denominator.");
        if (Math.Abs(num) > MaxDeltaTerm || Math.Abs(den) > MaxDeltaTerm)
            throw new ArgumentOutOfRangeException(nameof(num),
                $"|num| and |den| must stay within {MaxDeltaTerm}: clearing the denominator scales the whole " +
                "matrix and the ZZ diagonal already sums N-1 couplings, so a larger value could wrap int64. " +
                "The count is scale-free, so reduce the fraction instead.");
        if (den < 0) { num = -num; den = -den; }
        long[,] h = BuildCleared(num, den);

        int rank = 0;
        foreach (long p in primes) rank = Math.Max(rank, KrylovRankModP(h, seat, p));
        return N - rank;
    }

    /// <summary>The same count read at ONE named prime, exposed so that the one-sidedness argument
    /// above is falsifiable rather than merely asserted.
    ///
    /// <para>A rank over GF(p) can only drop, so this can only ever be too LARGE, and
    /// <see cref="BlindAtRational"/> is by construction the MINIMUM of this over
    /// <see cref="Primes"/>. Handing in a small prime makes a rank drop easy to produce, which is
    /// what pins the direction: if the maximum rank were not the one kept, a deliberately bad prime
    /// would change the public answer. Nothing in the object uses this; it exists to be tested
    /// against, and its guard is the same as <see cref="BlindAtRational"/>'s.</para></summary>
    public int BlindAtRationalModP(int seat, long num, long den, long prime)
    {
        RequireSeat(seat);
        if (prime < 2) throw new ArgumentOutOfRangeException(nameof(prime), "a ranking prime must be at least 2.");
        if (prime > MaxRankingPrime)
            throw new ArgumentOutOfRangeException(nameof(prime),
                $"a ranking prime must stay at or below {MaxRankingPrime}: the elimination multiplies two " +
                "reduced residues in int64, so a larger prime wraps silently and the rank it returns is " +
                "arithmetic noise rather than a bound. The shipped primes sit just inside it.");
        if (den == 0) throw new ArgumentOutOfRangeException(nameof(den), "a rational needs a nonzero denominator.");
        if (Math.Abs(num) > MaxDeltaTerm || Math.Abs(den) > MaxDeltaTerm)
            throw new ArgumentOutOfRangeException(nameof(num),
                $"|num| and |den| must stay within {MaxDeltaTerm}; see BlindAtRational for the bound's reason.");
        if (den < 0) { num = -num; den = -den; }
        return N - KrylovRankModP(BuildCleared(num, den), seat, prime);
    }

    /// <summary>The loop closed: (rational roots that ARE blind, rational roots total, deliberate
    /// non-roots that are correctly NOT blind, non-roots tried). The polynomial route and the
    /// Krylov route are computed from different arithmetic and must agree on both sides.</summary>
    public (int RootsBlind, int RootsTotal, int NonRootsClean, int NonRootsTried) RationalRootAgreements(int seat)
    {
        var roots = RationalLocusPoints(seat);
        int rootsBlind = roots.Count(r => BlindAtRational(seat, r.Num, r.Den) > 0);

        // Deliberate non-roots: small rationals that the polynomial does NOT vanish at. Without
        // these the agreement would be one-sided and a polynomial naming too MANY points would pass.
        var nonRoots = new List<(long, long)>();
        foreach ((long num, long den) in new[] { (1L, 3L), (2L, 5L), (5L, 2L), (-3L, 7L), (7L, 3L), (4L, 1L) })
            if (!roots.Contains((num, den))) nonRoots.Add((num, den));
        int clean = nonRoots.Count(r => BlindAtRational(seat, r.Item1, r.Item2) == 0);
        return (rootsBlind, roots.Count, clean, nonRoots.Count);
    }

    /// <summary>The blind count at an ALGEBRAIC Δ, exactly, with no eigensolver and no float: the
    /// MULTIPLICITY of Δ's minimal polynomial as a factor of P_j.
    ///
    /// <para><b>Why that is the count and not merely a membership test.</b> The locus is
    /// Δ_k = sin((j+1)kπ/N_node)/sin(jkπ/N_node) over k = 1..N_node−1, and distinct k give distinct
    /// θ_k in (0, π), hence distinct shared eigenvalues, so blind(Δ) = #{k : Δ_k = Δ}. P_j vanishes
    /// once per such k, so the multiplicity of a root IS the count there. Two k landing on one Δ is
    /// not exotic: at N = 11 seat 2 the polynomial is 3Δ⁴ − 4Δ² and 0 is a double root, where the XY
    /// law needs gcd(3, 12) − 1 = 2. Counting DISTINCT roots would give 1 and the law would fail.</para>
    ///
    /// <para><b>What this closes.</b> Every other route here decides rational Δ or membership only,
    /// so the rows this witness exists for, the ones at irrational Δ, had no producer in committed
    /// code: they came from a local scout. This one reaches them by integer polynomial division, so
    /// √3 and √2/2 are counted rather than cited.</para>
    ///
    /// <para>Pass the minimal polynomial as primitive integer coefficients, ASCENDING
    /// (Δ² − 3 is <c>{-3, 0, 1}</c>). It must be irreducible over ℚ, which is not checked: a
    /// reducible input would count the multiplicity of the whole product and silently under-report.
    /// End seats give 0 and the forced centre gives (N−1)/2 at every Δ, algebraic or not.</para>
    /// </summary>
    public int BlindAtAlgebraic(int seat, BigInteger[] minimalPolynomial)
    {
        RequireSeat(seat);
        if (minimalPolynomial is null)
            throw new ArgumentNullException(nameof(minimalPolynomial));
        if (N > MaxLocusN)
            throw new InvalidOperationException(
                $"N = {N} exceeds the locus-polynomial guard {MaxLocusN}; the guard is cost, not physics.");
        if (seat == 0 || seat == N - 1) return 0;
        if (NodeModulus(seat) == 0) return (N - 1) / 2;

        BigInteger[] q = Trim(LocusPolynomial(seat)!);
        BigInteger[] mu = Primitive(Trim(minimalPolynomial));
        // The degree is checked AFTER normalising, not before. An earlier version tested the raw
        // array's Length, which {5, 0} passes: it trims to a CONSTANT, every division by a constant
        // succeeds without lowering the degree, and the counting loop below never terminates. A
        // guard on the caller's formatting is not a guard on the object.
        RequireIrreducible(mu);
        if (mu.Length < 2)
            throw new ArgumentException(
                "a minimal polynomial must have degree at least 1 after trailing zeros are dropped; " +
                $"this one reduces to the constant [{string.Join(", ", mu)}], and dividing by a constant " +
                "never lowers a degree, so the multiplicity count would not terminate.",
                nameof(minimalPolynomial));
        int count = 0;
        while (q.Length >= mu.Length)
        {
            BigInteger[]? next = TryDivideExact(q, mu);
            if (next is null) break;
            q = next; count++;
        }
        return count;
    }

    /// <summary>The whole blind row at an algebraic Δ, seat by seat, so a committed row can be
    /// COMPARED rather than quoted. See <see cref="BlindAtAlgebraic"/>.</summary>
    public int[] BlindRowAtAlgebraic(BigInteger[] minimalPolynomial) =>
        Enumerable.Range(0, N).Select(s => BlindAtAlgebraic(s, minimalPolynomial)).ToArray();

    /// <summary>Refuses a minimal polynomial that is not irreducible over ℚ, because the count below
    /// is only meaningful for one Δ and a reducible input names several unrelated ones.
    ///
    /// <para><b>Why this is a guard and not a comment.</b> The local prototype this method ports,
    /// the number-field rank scout, raises on a reducible μ and tells its porter in as many words to
    /// "guard it and refuse rather than return a wrong rank". An earlier version of this file
    /// documented the hazard in prose instead, which is the shape where the repo answered and the
    /// answer was read as advice.</para>
    ///
    /// <para><b>What it decides, exactly.</b> Squarefree by gcd(μ, μ′), and no rational root by the
    /// rational-root theorem, which together settle every degree up to 3. At degree 4 and 5 a
    /// reducible polynomial with no rational root must have an integer quadratic factor, and that is
    /// searched exhaustively over the divisors of the constant and leading terms, so the decision is
    /// COMPLETE up to degree 5. At degree 6 and above it refuses rather than guessing: pass a factor.
    /// The realistic mistake, passing the whole locus polynomial (at N = 9 seat 1 that is
    /// Δ⁵ − 4Δ³ + 3Δ = Δ(Δ−1)(Δ+1)(Δ²−3)), is caught by the rational root 0.</para></summary>
    private static void RequireIrreducible(BigInteger[] mu)
    {
        int deg = mu.Length - 1;
        if (deg <= 1) return;
        if (deg >= 6)
            throw new ArgumentException(
                $"irreducibility is not decided here above degree 5, and this μ has degree {deg}. " +
                "Pass an irreducible factor of the locus polynomial rather than a product; the count " +
                "is the multiplicity of ONE Δ's minimal polynomial and a product names several.",
                nameof(mu));

        // squarefree: gcd(mu, mu') must be constant
        var d = new BigInteger[deg];
        for (int i = 1; i <= deg; i++) d[i - 1] = mu[i] * i;
        if (PolyGcdDegree(mu, Trim(d)) > 0)
            throw new ArgumentException(
                "this μ is not squarefree, so it is reducible and names no single Δ.", nameof(mu));

        // a rational root p/q needs p | mu[0] and q | mu[deg]
        foreach ((long num, long den) in RationalRootsOf(mu))
            throw new ArgumentException(
                $"this μ has the rational root {num}/{den}, so (({den}·Δ) − {num}) is a proper factor " +
                "and μ is reducible. The whole locus polynomial is the usual way to hit this: at " +
                "N = 9 seat 1 it is Δ⁵ − 4Δ³ + 3Δ, which has the rational root 0.", nameof(mu));

        if (deg < 4) return;                       // degree 2 and 3 are settled by the root test

        // degree 4 or 5 with no rational root: reducible iff it has an integer quadratic factor
        foreach (BigInteger a in DivisorsOf(mu[deg]))
            foreach (BigInteger c in DivisorsOf(mu[0]))
                for (BigInteger b = -MaxQuadraticMiddle; b <= MaxQuadraticMiddle; b++)
                {
                    var q = Trim(new[] { c, b, a });
                    if (q.Length < 3) continue;
                    if (TryDivideExact(mu, q) is not null)
                        throw new ArgumentException(
                            $"this μ factors: {a}Δ² + {b}Δ + {c} divides it, so it is reducible.",
                            nameof(mu));
                }
    }

    /// <summary>The search width for the middle coefficient of a quadratic factor. A factor
    /// aΔ² + bΔ + c of a degree ≤ 5 integer polynomial has |b| bounded by the sum of |coefficients|,
    /// which for every locus polynomial in the guarded range is far inside this.</summary>
    private const int MaxQuadraticMiddle = 64;

    private static IEnumerable<BigInteger> DivisorsOf(BigInteger n)
    {
        n = BigInteger.Abs(n);
        if (n.IsZero) { yield return BigInteger.One; yield break; }
        for (BigInteger k = 1; k * k <= n; k++)
            if ((n % k).IsZero)
            {
                yield return k; yield return -k;
                if (k * k != n) { yield return n / k; yield return -(n / k); }
            }
    }

    /// <summary>Degree of gcd over ℚ, by a fraction-free Euclid on primitive parts.</summary>
    private static int PolyGcdDegree(BigInteger[] a, BigInteger[] b)
    {
        a = Primitive(Trim(a)); b = Primitive(Trim(b));
        while (b.Length > 1 || !b[0].IsZero)
        {
            BigInteger[] r = PseudoRemainder(a, b);
            a = b; b = Primitive(Trim(r));
        }
        return a.Length - 1;
    }

    private static BigInteger[] PseudoRemainder(BigInteger[] a, BigInteger[] b)
    {
        int da = a.Length - 1, db = b.Length - 1;
        if (da < db) return a;
        var r = (BigInteger[])a.Clone();
        BigInteger lb = b[db];
        for (int k = da - db; k >= 0; k--)
        {
            BigInteger top = r[k + db];
            if (top.IsZero) continue;
            for (int i = 0; i < r.Length; i++) r[i] *= lb;
            top = r[k + db];
            BigInteger f = top / lb;
            for (int i = 0; i <= db; i++) r[k + i] -= f * b[i];
        }
        return Trim(r);
    }

    /// <summary>Exact division of integer polynomials, ascending coefficients, or null when the
    /// divisor does not divide over ℚ. Pseudo-division by the divisor's leading coefficient keeps
    /// everything in ℤ; the quotient is returned primitive, which is what the multiplicity count
    /// needs and what makes a repeated division terminate.</summary>
    private static BigInteger[]? TryDivideExact(BigInteger[] num, BigInteger[] den)
    {
        int dn = num.Length - 1, dd = den.Length - 1;
        if (dd < 0 || dn < dd) return null;
        BigInteger lead = den[dd];
        var rem = (BigInteger[])num.Clone();
        var quo = new BigInteger[dn - dd + 1];
        for (int k = dn - dd; k >= 0; k--)
        {
            BigInteger top = rem[k + dd];
            if (top.IsZero) { quo[k] = BigInteger.Zero; continue; }
            if (!BigInteger.Remainder(top, lead).IsZero)
            {
                // scale the whole remainder so the division stays exact in Z
                BigInteger g = BigInteger.GreatestCommonDivisor(BigInteger.Abs(top), BigInteger.Abs(lead));
                BigInteger f = lead / g;
                for (int i = 0; i < rem.Length; i++) rem[i] *= f;
                for (int i = 0; i < quo.Length; i++) quo[i] *= f;
                top = rem[k + dd];
            }
            BigInteger c = top / lead;
            quo[k] = c;
            for (int i = 0; i <= dd; i++) rem[k + i] -= c * den[i];
        }
        for (int i = 0; i < dd; i++) if (!rem[i].IsZero) return null;
        return Primitive(Trim(quo));
    }

    /// <summary>THE ROUTE THE POLYNOMIAL NEVER TOUCHES: sweep the WHOLE field 𝔽_p and ask the
    /// Krylov rank, at every Δ in it, whether the seat is blind. Returns the two sets, sorted:
    /// where the RANK says blind, and where the POLYNOMIAL vanishes.
    ///
    /// <para><b>Why this is the load-bearing cross-check and the rational one is not.</b>
    /// <see cref="RationalRootAgreements"/> can only ask about RATIONAL Δ, and the locus meets ℚ
    /// at 0 and ±1 alone (see <see cref="RationalLocusPoints"/>), which are the two committed
    /// endpoints and the sign partner. So that check never reaches a point this witness is the
    /// first to name. Over 𝔽_p every Δ is rational, the irrational members included: √2 has an
    /// image whenever 2 is a quadratic residue mod p, and the sweep meets it like any other
    /// residue. The check is also TWO-SIDED for free, since every Δ in the field is tried, so a
    /// polynomial naming too many or too few points fails on one side or the other.</para>
    ///
    /// <para><b>What each side is made of.</b> The rank side runs <see cref="BuildCleared"/> at
    /// den = 1 and eliminates over 𝔽_p: no Chebyshev polynomial, no Sylvester matrix, no Bareiss
    /// division, no interpolation. The polynomial side is <see cref="LocusPolynomial"/> evaluated
    /// by Horner. The two share the physical convention (hop 2, end shift 2Δ) and nothing else,
    /// which is why replacing either implementation by the other reddens this comparison while
    /// the rational one stays green.</para>
    ///
    /// <para><b>The one thing this does NOT do.</b> A rank mod p can only be too SMALL, so a
    /// nullity mod p can only be too LARGE: the sweep is an upper bound on blindness at each Δ,
    /// not a proof of it. What it certifies is that the polynomial's ROOT SET and the rank's
    /// BLIND SET are the same set, at every prime tried. Blindness itself at an irrational Δ is
    /// carried by Lemma J through <see cref="DefinitionPolynomial"/>, not by this.</para>
    ///
    /// <para>Returns null at an end seat and at the forced centre, where there is no finite locus
    /// to compare, matching <see cref="LocusPolynomial"/>.</para></summary>
    public (long[] RankSaysBlind, long[] PolynomialSaysBlind)? LocusOverFp(int seat, long prime)
    {
        RequireSeat(seat);
        if (!IsPrime(prime))
            throw new ArgumentOutOfRangeException(nameof(prime),
                $"{prime} is not prime. The elimination inverts by Fermat, m^(p-2), which is not an " +
                "inverse in a ring with zero divisors, so a composite modulus does not fail loudly: " +
                "both sides come back empty and the sweep reports agreement. Measured at 9, 15, 25 " +
                "and 2001.");
        if (prime < 3 || prime > MaxSweepPrime)
            throw new ArgumentOutOfRangeException(nameof(prime),
                $"the sweep prime must lie in 3..{MaxSweepPrime}: the sweep is p Krylov eliminations, " +
                "and the entries stay well inside int64 only while p is small. It is a cost bound, not physics.");
        if (NodeModulus(seat) == 0 || seat == 0 || seat == N - 1) return null;

        // No "p divides every coefficient" guard here, and the reason is that it could never fire:
        // LocusPolynomial returns a PRIMITIVE polynomial, so its content is 1 by construction and no
        // prime divides all of it. A guard that cannot fire reads as a defence and is a decoration.
        BigInteger[] poly = LocusPolynomial(seat)!;

        var byRank = new List<long>();
        var byPoly = new List<long>();
        for (long delta = 0; delta < prime; delta++)
        {
            if (N - KrylovRankModP(BuildCleared(delta, 1), seat, prime) > 0) byRank.Add(delta);

            long acc = 0;                                          // Horner, mod p
            for (int i = poly.Length - 1; i >= 0; i--)
                acc = (acc * delta + (long)Mod(poly[i], prime)) % prime;
            if (acc == 0) byPoly.Add(delta);
        }
        return (byRank.ToArray(), byPoly.ToArray());
    }

    private static BigInteger Mod(BigInteger v, long p)
    {
        BigInteger r = v % p;
        return r.Sign < 0 ? r + p : r;
    }

    // -----------------------------------------------------------------------------------------
    // exact integer helpers: Chebyshev, Sylvester/Bareiss, interpolation
    // -----------------------------------------------------------------------------------------

    /// <summary>U_n(x), ascending integer coefficients. U_0 = 1, U_1 = 2x, U_n = 2x·U_{n−1} − U_{n−2}.
    /// U_{−1} is the zero polynomial, which is what the j = 0 end seat would need and never asks for.
    /// </summary>
    public static BigInteger[] ChebyshevU(int n)
    {
        if (n < 0) return new[] { BigInteger.Zero };
        var prev = new[] { BigInteger.One };                       // U_0
        if (n == 0) return prev;
        var cur = new[] { BigInteger.Zero, new BigInteger(2) };    // U_1
        for (int i = 2; i <= n; i++)
        {
            var next = new BigInteger[i + 1];
            for (int k = 0; k < cur.Length; k++) next[k + 1] += 2 * cur[k];
            for (int k = 0; k < prev.Length; k++) next[k] -= prev[k];
            prev = cur;
            cur = next;
        }
        return cur;
    }

    /// <summary>t·p − q, ascending coefficients.</summary>
    private static BigInteger[] SubScaled(BigInteger[] p, long t, BigInteger[] q)
    {
        int len = Math.Max(p.Length, q.Length);
        var r = new BigInteger[len];
        for (int i = 0; i < p.Length; i++) r[i] += t * p[i];
        for (int i = 0; i < q.Length; i++) r[i] -= q[i];
        return Trim(r);
    }

    private static BigInteger[] Trim(BigInteger[] p)
    {
        int last = p.Length - 1;
        while (last > 0 && p[last].IsZero) last--;
        return p.Take(last + 1).ToArray();
    }

    /// <summary>Res(a, b) for integer polynomials given ascending, as the Sylvester determinant taken
    /// by fraction-free Bareiss elimination. Exact: every division in Bareiss is exact by
    /// construction, and the assertion below is not decoration, it is the invariant that makes the
    /// routine integer.</summary>
    public static BigInteger ResultantInt(BigInteger[] a, BigInteger[] b)
    {
        a = Trim(a); b = Trim(b);
        int m = a.Length - 1, n = b.Length - 1;
        if (m <= 0 && n <= 0) return BigInteger.One;
        if (m < 0 || n < 0) return BigInteger.Zero;
        if (n == 0) return BigInteger.Pow(b[0], m);
        if (m == 0) return BigInteger.Pow(a[0], n);

        int size = m + n;
        var s = new BigInteger[size, size];
        for (int r = 0; r < n; r++)
            for (int k = 0; k <= m; k++)
                s[r, r + k] = a[m - k];                 // descending coefficients of a
        for (int r = 0; r < m; r++)
            for (int k = 0; k <= n; k++)
                s[n + r, r + k] = b[n - k];             // descending coefficients of b
        return Bareiss(s, size);
    }

    private static BigInteger Bareiss(BigInteger[,] mat, int n)
    {
        BigInteger prev = BigInteger.One;
        int sign = 1;
        for (int k = 0; k < n - 1; k++)
        {
            if (mat[k, k].IsZero)
            {
                int swap = -1;
                for (int r = k + 1; r < n; r++)
                    if (!mat[r, k].IsZero) { swap = r; break; }
                if (swap < 0) return BigInteger.Zero;
                for (int c = 0; c < n; c++) (mat[k, c], mat[swap, c]) = (mat[swap, c], mat[k, c]);
                sign = -sign;
            }
            for (int i = k + 1; i < n; i++)
            for (int j = k + 1; j < n; j++)
            {
                BigInteger num = mat[i, j] * mat[k, k] - mat[i, k] * mat[k, j];
                mat[i, j] = BigInteger.Divide(num, prev);
                if (mat[i, j] * prev != num)
                    throw new InvalidOperationException(
                        "Bareiss elimination hit an inexact division, which cannot happen for an integer " +
                        "Sylvester matrix; the construction above is wrong rather than the arithmetic.");
            }
            prev = mat[k, k];
        }
        return sign * mat[n - 1, n - 1];
    }

    /// <summary>The unique polynomial of degree &lt; ys.Length through (0, ys[0]) .. (n, ys[n]),
    /// ascending. Lagrange over a common denominator n!, so every step is integer and the final
    /// division is exact by construction.</summary>
    public static BigInteger[] InterpolateAtZeroToN(BigInteger[] ys)
    {
        int n = ys.Length - 1;
        if (n == 0) return new[] { ys[0] };
        var acc = new BigInteger[n + 1];
        for (int i = 0; i <= n; i++)
        {
            if (ys[i].IsZero) continue;
            // prod_{m != i} (t - m), expanded
            var poly = new BigInteger[] { BigInteger.One };
            for (int m = 0; m <= n; m++)
            {
                if (m == i) continue;
                var next = new BigInteger[poly.Length + 1];
                for (int k = 0; k < poly.Length; k++)
                {
                    next[k + 1] += poly[k];
                    next[k] -= m * poly[k];
                }
                poly = next;
            }
            BigInteger weight = ys[i] * Binomial(n, i) * ((n - i) % 2 == 0 ? 1 : -1);
            for (int k = 0; k < poly.Length; k++) acc[k] += weight * poly[k];
        }
        BigInteger fact = BigInteger.One;
        for (int i = 2; i <= n; i++) fact *= i;
        var res = new BigInteger[n + 1];
        for (int k = 0; k <= n; k++)
        {
            res[k] = BigInteger.Divide(acc[k], fact);
            if (res[k] * fact != acc[k])
                throw new InvalidOperationException(
                    "the interpolated coefficient is not integral, which cannot happen for a resultant of " +
                    "integer polynomials; the sample points or the weights above are wrong.");
        }
        return Trim(res);
    }

    private static BigInteger Binomial(int n, int k)
    {
        BigInteger r = BigInteger.One;
        for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1);
        return r;
    }

    /// <summary>Divide out the content and make the leading coefficient positive.</summary>
    private static BigInteger[] Primitive(BigInteger[] p)
    {
        p = Trim(p);
        BigInteger g = BigInteger.Zero;
        foreach (BigInteger c in p) g = BigInteger.GreatestCommonDivisor(g, BigInteger.Abs(c));
        if (g.IsZero) return new[] { BigInteger.Zero };
        var q = p.Select(c => c / g).ToArray();
        if (q[^1].Sign < 0) q = q.Select(c => -c).ToArray();
        return q;
    }

    private static IEnumerable<BigInteger> Divisors(BigInteger v)
    {
        if (v.IsZero) { yield return BigInteger.One; yield break; }
        for (BigInteger i = BigInteger.One; i * i <= v; i++)
        {
            if (!(v % i).IsZero) continue;
            yield return i;
            if (i * i != v) yield return v / i;
        }
    }

    /// <summary>Whether p(num/den) = 0, evaluated as the integer Σ c_k·num^k·den^(deg−k).</summary>
    private static bool EvaluatesToZero(BigInteger[] p, BigInteger num, BigInteger den)
    {
        int deg = p.Length - 1;
        BigInteger acc = BigInteger.Zero;
        for (int k = 0; k <= deg; k++)
            acc += p[k] * BigInteger.Pow(num, k) * BigInteger.Pow(den, deg - k);
        return acc.IsZero;
    }

    /// <summary>den·H(num/den), integral by construction: the hop pays 2·den on each bond and the
    /// ZZ diagonal pays num per bond, −num at a bond's own two ends. A Krylov rank is invariant
    /// under a nonzero scaling of H, so clearing the denominator changes no count.</summary>
    private long[,] BuildCleared(long num, long den)
    {
        var h = new long[N, N];
        for (int b = 0; b < N - 1; b++)
        {
            h[b, b + 1] += 2 * den;
            h[b + 1, b] += 2 * den;
        }
        for (int s = 0; s < N; s++)
            for (int b = 0; b < N - 1; b++)
                h[s, s] += (s == b || s == b + 1) ? -num : num;
        return h;
    }

    private static int KrylovRankModP(long[,] h, int seat, long p)
    {
        int n = h.GetLength(0);
        var rows = new List<long[]>();
        var vec = new long[n];
        vec[seat] = 1;
        for (int k = 0; k <= n; k++)
        {
            rows.Add((long[])vec.Clone());
            var next = new long[n];
            for (int a = 0; a < n; a++)
            {
                long s = 0;
                for (int b = 0; b < n; b++) s = (s + Mod(h[a, b], p) * vec[b]) % p;
                next[a] = s;
            }
            vec = next;
        }
        return RankModP(rows, n, p);
    }

    private static long Mod(long v, long p) { long r = v % p; return r < 0 ? r + p : r; }

    /// <summary>Trial division, which is ample below <see cref="MaxRankingPrime"/> and is only ever
    /// asked once per call.</summary>
    private static bool IsPrime(long n)
    {
        if (n < 2) return false;
        if (n % 2 == 0) return n == 2;
        for (long f = 3; f <= n / f; f += 2) if (n % f == 0) return false;
        return true;
    }

    private static int RankModP(List<long[]> rows, int cols, long p)
    {
        var m = rows.Select(r => r.Select(v => Mod(v, p)).ToArray()).ToList();
        int rank = 0;
        for (int c = 0; c < cols && rank < m.Count; c++)
        {
            int pivot = -1;
            for (int r = rank; r < m.Count; r++) if (m[r][c] != 0) { pivot = r; break; }
            if (pivot < 0) continue;
            (m[rank], m[pivot]) = (m[pivot], m[rank]);
            long inv = ModPow(m[rank][c], p - 2, p);
            for (int r = rank + 1; r < m.Count; r++)
            {
                if (m[r][c] == 0) continue;
                long f = m[r][c] * inv % p;
                for (int k = c; k < cols; k++) m[r][k] = (m[r][k] - f * m[rank][k]) % p;
                for (int k = c; k < cols; k++) if (m[r][k] < 0) m[r][k] += p;
            }
            rank++;
        }
        return rank;
    }

    private static long ModPow(long b, long e, long p)
    {
        long r = 1; b = Mod(b, p);
        while (e > 0)
        {
            if ((e & 1) == 1) r = (long)((System.Int128)r * b % p);
            b = (long)((System.Int128)b * b % p);
            e >>= 1;
        }
        return r;
    }

    private void RequireSeat(int seat)
    {
        if (seat < 0 || seat >= N)
            throw new ArgumentOutOfRangeException(nameof(seat), $"seat must lie in 0..{N - 1}; got {seat}.");
    }

    internal static string Format(BigInteger[] p)
    {
        var parts = new List<string>();
        for (int k = p.Length - 1; k >= 0; k--)
        {
            if (p[k].IsZero) continue;
            string c = BigInteger.Abs(p[k]) == BigInteger.One && k > 0
                ? "" : BigInteger.Abs(p[k]).ToString(Inv);
            string x = k == 0 ? "" : k == 1 ? "D" : $"D^{k}";
            parts.Add($"{(p[k].Sign < 0 ? "-" : parts.Count == 0 ? "" : "+")} {c}{x}".Trim());
        }
        return parts.Count == 0 ? "0" : string.Join(" ", parts);
    }

    // -----------------------------------------------------------------------------------------
    // the inspect tree
    // -----------------------------------------------------------------------------------------

    public string DisplayName => $"SeatBlindnessDeltaLocusWitness (F157 on the anisotropy axis, N={N})";

    public string Summary
    {
        get
        {
            var forced = Enumerable.Range(0, N).Where(BlindAtEveryDelta).ToArray();
            var met = Enumerable.Range(1, Math.Max(0, N - 2)).Where(HasNonemptyLocus).ToArray();
            return $"N={N}: N_node = |N-1-2j| governs. Forced (blind at every Delta): " +
                   $"{(forced.Length == 0 ? "none, N is even" : string.Join(", ", forced))}. " +
                   $"Interior seats with a nonempty finite locus: {met.Length}. " +
                   $"Locus polynomial {(N <= MaxLocusN ? "built live" : $"not built (N > {MaxLocusN})")}.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. Per seat: nodeModulus, the locus, and the two endpoint counts.
            var rows = new List<IInspectable>();
            for (int seat = 0; seat < N && seat < 24; seat++)
            {
                int nodeModulus = NodeModulus(seat);
                string body;
                if (seat == 0 || seat == N - 1)
                    body = "an END seat: one principal submatrix is empty, so there is no shared root to have; " +
                           $"blind = {BlindAtRational(seat, 1, 1)} at Delta = 1 and " +
                           $"{BlindAtRational(seat, 0, 1)} at Delta = 0, and 0 at every Delta.";
                else if (nodeModulus == 0)
                    body = $"the reflection-fixed CENTRE seat, N_node = 0: blind at EVERY Delta, and the count is " +
                           $"(N-1)/2 = {(N - 1) / 2}. Live at Delta = 1: {BlindAtRational(seat, 1, 1)}; " +
                           $"at Delta = 3: {BlindAtRational(seat, 3, 1)}. The two halves are conjugate by the " +
                           "chain reflection, so their resultant vanishes identically.";
                else
                {
                    string poly = N <= MaxLocusN ? Format(LocusPolynomial(seat)!) : "not built";
                    var rat = N <= MaxLocusN ? RationalLocusPoints(seat) : Array.Empty<(long, long)>();
                    string ratTxt = rat.Count == 0 ? "none" :
                        string.Join(", ", rat.Select(r => r.Den == 1
                            ? r.Num.ToString(Inv)
                            : $"{r.Num.ToString(Inv)}/{r.Den.ToString(Inv)}"));
                    body = $"N_node = {nodeModulus}; locus polynomial P(D) = {poly}; its RATIONAL points: {ratTxt}; " +
                           $"counts from the congruences: {CountAtIsotropic(seat)} at Delta = 1, " +
                           $"{CountAtXy(seat)} at Delta = 0.";
                }
                rows.Add(new InspectableNode($"seat {seat}", summary: body, provenance: NodeProvenance.Live));
            }
            yield return new InspectableNode("per seat: N_node, the locus polynomial, the endpoint counts",
                summary: N > 24 ? $"the first 24 of {N} seats" : $"all {N} seats, recomputed at inspect time",
                children: rows, provenance: NodeProvenance.Live);

            // 2. The two routes meeting, which is the whole point.
            if (N <= MaxLocusN)
            {
                int rb = 0, rt = 0, nc = 0, nt = 0;
                for (int seat = 1; seat < N - 1; seat++)
                {
                    if (NodeModulus(seat) == 0) continue;
                    var (a, b, c, e) = RationalRootAgreements(seat);
                    rb += a; rt += b; nc += c; nt += e;
                }
                yield return new InspectableNode("the rational check, and what it does NOT reach",
                    summary: $"every RATIONAL root of the locus polynomial must leave the seat blind, and a " +
                             $"deliberate non-root must not. Roots blind: {rb} of {rt}. Non-roots correctly " +
                             $"clean: {nc} of {nt}. " +
                             (rb == rt && nc == nt
                                ? "They agree, on both sides. "
                                : "THEY DISAGREE, which on a uniform chain is a finding about the construction. ") +
                             "BUT READ THE SCOPE: the locus meets Q only at 0 and +-1, which are the two " +
                             "COMMITTED endpoints and the sign partner, so this check certifies nothing this " +
                             "witness is the first to name. It is a regression guard, not the evidence.",
                    provenance: NodeProvenance.Live);

                // 2b. The two POLYNOMIAL routes, which the registry says are computed beside each
                //     other and which nothing in this tree used to show.
                int same = 0, both = 0;
                for (int seat = 1; seat < N - 1; seat++)
                {
                    if (NodeModulus(seat) == 0) continue;
                    both++;
                    if (LocusPolynomial(seat)!.SequenceEqual(DefinitionPolynomial(seat)!)) same++;
                }
                yield return new InspectableNode("the two polynomial routes, live",
                    summary: $"the Chebyshev resultant Res_x(U_(N_node-1), D*U_(j-1) - U_j) against the " +
                             $"definition route Res_lambda(chi_L, chi_R), built from the two submatrices by the " +
                             $"Jacobi three-term recursion with no Chebyshev identity used: {same} of {both} " +
                             "interior non-centre seats agree as primitive integer polynomials. The second is a " +
                             "DERIVATION and not a second implementation, which is why its end shift is a " +
                             "parameter: feed it the physically wrong shift and the two must part.",
                    provenance: NodeProvenance.Live);

                // 2c. The field sweep: the only route whose arithmetic the polynomial never touches.
                const long SweepPrime = 101;
                int agree = 0, tried = 0, withLocus = 0;
                for (int seat = 1; seat < N - 1; seat++)
                {
                    var got = LocusOverFp(seat, SweepPrime);
                    if (got is null) continue;
                    tried++;
                    if (got.Value.RankSaysBlind.SequenceEqual(got.Value.PolynomialSaysBlind)) agree++;
                    if (got.Value.PolynomialSaysBlind.Length > 0) withLocus++;
                }
                yield return new InspectableNode("the field sweep, which is the load-bearing one",
                    summary: $"at every one of the {SweepPrime} residues of GF({SweepPrime}), the Krylov rank " +
                             $"is asked whether the seat is blind, and the blind SET is compared to the " +
                             $"polynomial's ROOT SET: {agree} of {tried} interior non-centre seats agree, " +
                             (withLocus > 0
                                ? $"{withLocus} of them with a nonempty locus, so the comparison is not 0 == 0. "
                                : "and NONE of them has a nonempty locus at this N and this prime, so the " +
                                  "comparison here IS 0 == 0 and carries nothing. Read it at a larger N. ") +
                             "No Chebyshev polynomial, Sylvester matrix, Bareiss division or interpolation " +
                             "takes part on the rank side, and every Delta in the field is tried, so the check " +
                             "is two-sided by construction. This is also where the IRRATIONAL members are " +
                             "reached: over GF(p) they have images like any other residue. What it does NOT " +
                             "do is prove blindness, since a rank mod p can only be too small; it certifies " +
                             "that the two SETS agree, and blindness itself is Lemma J through the resultant.",
                    provenance: NodeProvenance.Live);
            }

            // 3. The two committed endpoints, reproduced from nodeModulus alone.
            int agreeH = 0, agreeX = 0, seats = 0;
            for (int seat = 1; seat < N - 1; seat++)
            {
                if (NodeModulus(seat) == 0) continue;
                seats++;
                if (CountAtIsotropic(seat) == SeatCutBlindnessClaim.BlindHeisenberg(N, seat)) agreeH++;
                if (CountAtXy(seat) == SeatCutBlindnessClaim.BlindXy(N, seat)) agreeX++;
            }
            yield return new InspectableNode("the two committed endpoints, recovered from N_node",
                summary: $"the k-congruences on N_node reproduce (gcd(2j+1, N) - 1)/2 at {agreeH} of {seats} interior " +
                         $"non-centre seats and gcd(j+1, N+1) - 1 at {agreeX} of {seats}. The two committed gcd " +
                         "laws are therefore the Delta = 1 and Delta = 0 sections of one locus, and N_node is what " +
                         "carries both.",
                provenance: NodeProvenance.Live);

            // 4. Falsifiers, and the one input that cannot break the law.
            yield return new InspectableNode("the falsifiers, and the one non-probe",
                summary: $"end seats, live: blind(0) = {BlindAtRational(0, 1, 1)} and " +
                         $"blind({N - 1}) = {BlindAtRational(N - 1, 1, 1)} at Delta = 1, both required 0. " +
                         $"A generic Delta = 1/3 leaves blind = " +
                         $"[{string.Join(", ", Enumerable.Range(0, Math.Min(N, 12)).Select(s => BlindAtRational(s, 1, 3)))}], " +
                         "which must be zero away from the forced centre. ONE INPUT CANNOT BREAK THE LAW " +
                         "and is never to be counted as a probe: Delta = -1, because with Sigma the staggering " +
                         "Sigma*H(Delta)*Sigma = -H(-Delta) makes blind(-Delta) = blind(Delta) identically " +
                         "(F152's bipartite cospectrality reason). AND ONE THAT IS OFTEN MISTAKEN FOR ONE: prime " +
                         "N is a weak probe of the ZZ book ONLY, where gcd(2j+1, N) forces every seat but the " +
                         "centre to zero. It is NOT a non-probe of this law. At the XY endpoint gcd(j+1, N+1) " +
                         "blinds seven seats at N = 11 and three at N = 7, and the closed form predicts them " +
                         "rather than merely permitting them: the N = 7 seat 1 polynomial D^3 - 2D has root 0. " +
                         "(That row is printed above only when this witness is built at N = 7.)",
                provenance: NodeProvenance.Live);

            // 5. Scope, written down rather than recomputed.
            yield return new InspectableNode("scope and fences",
                summary: "The UNIFORM chain only. The locus is derived from the two principal submatrices the " +
                         "seat leaves behind, which is F157's FENCED phrasing: SeatCutBlindnessClaim records it " +
                         "as wrong on all 1682 zero-bond (profile, seat) pairs at N = 3..6 ON THE HEISENBERG " +
                         "BOOK and on most, not all, of the XY ones; the 60 that still come out right are counted in " +
                         "experiments/THE_BLIND_SITE.md and in F4's seat bullet, not in the claim. So the " +
                         "derivation " +
                         "needs a chain with no zero bond, which a uniform chain is. The repo already owns the " +
                         "trigonometric node route at BOTH endpoints, exactly and to N = 200 " +
                         "(simulations/seat_cut_blindness.py 'Route 3, EXACT', simulations/blind_site.py " +
                         "mode_set); what is new here is the general FORM at every (N, j), with N_node as the " +
                         "single controlling integer; the MECHANISM was Delta-general already: Lemma J (J4) of " +
                         "PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md gives blind = deg gcd(chi_L, chi_R) for any " +
                         "unreduced Jacobi matrix, and that proof's section (f) already records that the ZZ " +
                         "term 'changes only the diagonal, which Lemma J never constrains'. What is added is " +
                         "the EVALUATION of that intersection as Delta moves. TWO LIMITS THAT ARE REAL. (1) N_node " +
                         "is a statement about block SIZES, while THE_SEAT_THAT_CUTS records that with the ZZ " +
                         "term the two submatrices are not free-standing subchains, carrying a boundary term " +
                         "at the cut and a different shift on each side; that size alone governs is the " +
                         "surprising part, not a step. (2) The independent Krylov route here decides only " +
                         "RATIONAL Delta, and the locus meets Q only at 0 and +-1 (MEASURED over N = 4..24 at every "
                         + "interior seat, a statement about this family and not a theorem), the two committed endpoints " +
                         "and the sign partner, so every IRRATIONAL member is carried by the polynomial route " +
                         "alone. The ceiling blind(j) <= min(j, N-1-j) belongs to the two-halves criterion and " +
                         "not to the fence-free count; it is respected here at every Delta. A Delta is NOT the " +
                         "detuned bond that THE_SEAT_THAT_CUTS leaves open; do not report one as the other.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
