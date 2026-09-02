using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F161, each computed order of a collision gap reads ONE comb under an integer multiplier, and which
/// orders vanish is decided by a gcd. Switch on the cracked ring's wrap bond (F160's road, u = J'/J) and expand a
/// chain level in u, E_k(u) = 2cos θ_k + Σ_m d_m·u^m with θ_k = kπ/n. Then through fifth order
///
/// <code>
///     d_1 = η·(2/n)·(1 − cos 2θ)                              η = (−1)^(k+1)
///     d_2 = ((n−3)/n^2)·(cos θ − cos 3θ)
///     d_3 = η·(2(n−2)(n−4)/(3n^3))·(cos 2θ − cos 4θ)
///     d_5 = η·F·[(5/4)(n−1)·cos 2θ + (3n^2/2 − 8n + 8)·cos 4θ − (3/4)(2n−3)(n−3)·cos 6θ],  F = 4(n−6)(n−2)/(15 n^5)
/// </code>
///
/// so each coefficient is a signed combination of neighbouring evaluations of one comb: the EVEN orders on
/// M_j(τ) = Σ_{k∈τ} cos(jkπ/n) with ODD multipliers, the ODD orders on X_2j(τ) = Σ_{k∈τ} η_k·cos(2jkπ/n) with even
/// ones, and X is the same comb SHIFTED, X_2j(τ) = −M_{n+2j}(τ). That split is a theorem at every order (two
/// symmetries of the angle equation, §(d) of the proof), not a pattern in the five computed.
///
/// <para><b>The gcd.</b> k ↦ k·m mod 2n induces a Galois automorphism of ℚ(ζ_2n) exactly when gcd(m, 2n) = 1. For a
/// COLLISION pair (two clean triples with equal level, F129) the collision IS ΔM_1 = 0, so any automorphism sends it
/// to zero again. On the X ladder the multiplier is n + 2j, and at ODD n that reduces to gcd(j, n) = 1, which j = 1
/// and j = 2 always satisfy: hence <b>ΔX_2 = ΔX_4 = 0 and c_3 = 0 at every odd n, for every collision pair, standing
/// or separating</b>. The rung j = 0 is never an automorphism (gcd(n, 2n) = n) and it is exactly the rung that leaves
/// the first order standing, c_1 = (4/n)(o_τ − o_σ) at odd n, the difference of odd-label counts. At EVEN n no rung of
/// the X ladder is reached from the collision at all, and the twelve non-mirror standing pairs of the n ≤ 30 census
/// are carried instead by the Conway-Jones ROT3 shape of their doubled labels, which forces X_2j = 0 at every rung
/// with 3 ∤ j; at 3 | j the coset collapses to 3η·cos(2πja/n), a cosine that need not vanish (n = 24, τ = (1,7,9) is
/// the break-input, so the lemma has ONE forced direction).</para>
///
/// <para><b>Where the ladder stops is F129's own divisor.</b> F129 fires only at 3|n ≥ 9 or 10|n ≥ 20, and 10|n forces
/// n even, so every ODD firing modulus has 3|n. At odd n the first rung that is not an automorphism is therefore
/// j = 3; at even n the first rung that collapses a ROT3 coset is also j = 3. Either way X_6 is the first rung that
/// can survive, and by the multiplier form X_6 first enters at FIFTH order. So for a pair that stands at first order,
/// g(u) = c_2u^2 + c_4u^4 + c_5u^5 + O(u^6) with c_5 = −(3/4)(2n−3)(n−3)·F·ΔX_6: the same prime that lets the
/// coincidence exist at an odd modulus is the prime that ends its protection.</para>
///
/// <para><b>The criterion is local to the PIECES, and that is the sharp form.</b> gcd(m, 2n) = 1 is sufficient and too
/// crude. A collision reduces to a vanishing sum of roots of unity, which decomposes into minimal vanishing pieces;
/// what governs a piece is its RATIO-ORDER o (all its exponents sit in one coset of ⟨2n/o⟩, so the piece is a unit
/// times a vanishing element of ℤ[ζ_o]), NOT a prime-coset structure. Reading under the multiplier m keeps a piece
/// zero exactly while gcd(m, o) = 1, so the statement is: if SOME minimal tiling has every ratio-order coprime to m,
/// the reading stays zero. Global invertibility is the special case where every tiling qualifies. Under it Theorem D
/// (the pair's sum in ℤ[ζ_2n]) and Theorem E (a single triple's doubled-label sum in ℤ[ζ_n]) are one criterion at its
/// two ends. At m = 3 the census adds the converse exactly: c_2 = 0 on exactly the pairs admitting some entirely
/// 3-free minimal tiling, which inside the census are F129's family C, and family L is 3-free too, so past the census
/// <see cref="SecondOrderZeroPairLowerBound"/> = Count(C, n) + Count(L, n) is a LOWER BOUND and not the count.</para>
///
/// <para><b>What this class is NOT.</b> No hardware claim. The letter n here is the COMB modulus n = N + 1, never the
/// site count (the crack code spends the same letter the other way; every member below says <c>nComb</c>). The letter
/// c_m is F160's split correction ½ − 1/(N·sin²k_m) on a different expansion point, at the RING end and about a pair
/// split; nothing here transfers to it, which is why the coefficients below are named by their order and not by that
/// letter. A <i>rung</i> here is one value of j and not the glossary's rate rung; a <i>pair</i> is two triples and not
/// MirrorWorld's coherence <c>Pair</c>; a <i>ladder</i> is the sequence of evaluations of this one comb; and
/// <i>divisor</i> in the heading above is the integer 3 dividing n, not F140's frozen <c>Divisor</c> and not F157's
/// divisor law. LEMMA B's reflection is the chain reversal R with R·ψ_k = η_k·ψ_k (F71's mirror sign, carried in
/// PROOF_K_PARTNERSHIP's closing section), NOT <see cref="ChiralKClaim"/>'s sublattice K = diag((−1)^ℓ), which SWAPS
/// modes; the two share a proof file and are different operators. That fence is about Lemma B alone: the chiral K
/// does appear here, at <see cref="IsThetaMirrorPair"/>, where K·H(u)·K = −H(−u) is the all-orders route for the
/// eleven Θ-mirror pairs, and it is cited there rather than used, this claim reaching those pairs through the ladder
/// and only to fifth order. The parity split above is proved from the angle equation and not from either operator, a
/// route a review lens supplied after the first draft went through the chiral K and left a degree bound open; that is
/// also why ChiralKClaim is not a typed parent although it would be tier-legal.</para>
///
/// <para><b>What is READ and not derived</b> (the proof's own fence): that all twelve even-n non-mirror standing pairs
/// of the n ≤ 30 census have both triples of the ROT3 shape AND share a parity class. The consequence is a theorem;
/// the hypothesis is a census observation. Whether a standing pair at even n must be ROT3, and whether a mixed-parity
/// pair can stand, are open (the arc <c>the_forced_and_the_met</c>).</para>
///
/// <para>Proof <c>docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md</c>, gate
/// <c>simulations/collision_gap_odd_orders.py</c> (fifteen blocks L1..L14, every one exact except L4, an error-model
/// law on an eigensolver). The census this proof explains is <c>experiments/THE_COMB_ON_THE_ROAD.md</c>. Everything
/// below recomputes at call time, in ℤ[ζ_2n] over <see cref="CyclotomicRing"/> or in exact rationals; no count is
/// stored.</para></summary>
public sealed class CollisionGapOddOrdersClaim : Claim
{
    /// <summary>Parent: F160, whose Theorem A polynomial det(x·I − H(u)) = U_N(x/2) − u²·U_{N−2}(x/2) − 2u is the
    /// whole input to the series, and whose Theorem G is this series read at m = 1. F129 is NOT a typed parent: its
    /// typed home <see cref="CrossTripleOrthogonalityClaim"/> is Tier1Candidate (its code-trust caveat), and the tier
    /// rule forbids a Tier1Candidate parent for a Tier1Derived child, so F129's firing condition and family inventory
    /// are carried by anchor and by the executable <see cref="LevelCollisionCensus.Fires"/> and
    /// <see cref="CollisionFamilyInventory.Count"/> instead. F2b and F65 are ancestors through this edge and are not
    /// re-added.</summary>
    public CrackedRingExactCurveClaim CrackedRing { get; }

    public CollisionGapOddOrdersClaim(CrackedRingExactCurveClaim crackedRing)
        : base("F161 each computed order of a collision gap reads one comb under an integer multiplier, and which " +
               "orders vanish is a gcd: expanding a chain level of the cracked ring in the wrap bond, " +
               "E_k(u) = 2cos(k pi/n) + sum_m d_m u^m with n = N+1 the COMB modulus and eta = (-1)^(k+1), the five " +
               "computed coefficients each sit on a multiplier ladder, the even orders on M_j = sum cos(j k pi/n) " +
               "with odd multipliers and the odd orders on X_2j = sum eta_k cos(2j k pi/n) with even ones, " +
               "X_2j = -M_{n+2j}; the PARITY of the multiplier is a theorem at every order, from two symmetries of " +
               "the angle equation, and the RANGE is not proved. k -> k*m mod 2n is a Galois automorphism of " +
               "Q(zeta_2n) exactly when gcd(m, 2n) = 1, which on the X ladder at ODD n is gcd(j, n) = 1, so for " +
               "EVERY collision pair at odd n (all 627 of the n <= 30 census, standing or separating) " +
               "dX_2 = dX_4 = 0 and hence c_3 = 0, while the rung j = 0 is never an automorphism and is the rung " +
               "that leaves c_1 = (4/n)(o_tau - o_sigma) standing; at EVEN n no rung is reached from the collision " +
               "and the twelve non-mirror standing pairs are carried instead by the doubled-label ROT3 shape, which " +
               "forces X_2j = 0 whenever 3 does not divide j and collapses to 3 eta cos(2 pi j a/n) when it does, so " +
               "the lemma has ONE forced direction (n = 24 breaks the converse); since F129 fires only at 3|n or " +
               "10|n and 10|n forces n even, every odd firing modulus has 3|n, so j = 3 is the first surviving rung " +
               "on both sides and X_6 first enters at FIFTH order, c_5 = -(3/4)(2n-3)(n-3)*F*dX_6 with " +
               "F = 4(n-6)(n-2)/(15 n^5). The sharp criterion is local to the minimal PIECES of the vanishing sum, " +
               "not the whole ring: some minimal tiling with every RATIO-ORDER coprime to m keeps the reading zero, " +
               "of which global invertibility is the special case, and at m = 3 the census adds the converse " +
               "exactly, the 60 pairs with c_2 = 0 being F129's family C inside it, with family L 3-free as well, " +
               "so 2(n-10)*[10|n] + 20*[70|n] is a LOWER BOUND past n = 30. Every count is at n <= 30; the ROT3 " +
               "shape of the twelve is READ and not derived",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md (primary: Theorem A the five coefficients, Lemma B the " +
               "rank-one sector fold, Corollary C the multiplier form, Theorem D the Galois kill and its converse, " +
               "Theorem E the ROT3 rung lemma, Theorem F where the ladder stops, Corollary G the second order and " +
               "the piece criterion) + " +
               "experiments/THE_COMB_ON_THE_ROAD.md (the census this proof explains: 2558 pairs at the nine firing " +
               "n <= 30, 2335 separating at first order and 223 standing) + " +
               "docs/proofs/PROOF_CRACKED_RING_EXACT_CURVE.md (F160 Theorem A, the polynomial that is the whole " +
               "input, and Theorem G, this series at first order) + " +
               "docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md (the firing condition 3|n or 10|n, used and not " +
               "reproved) + " +
               "docs/proofs/PROOF_F129_FAMILY_INVENTORY_COUNTS.md (families C and L, their doors and closed-form " +
               "counts, which are the second order's bound) + " +
               "experiments/F89_SEED_EXISTENCE_REDUCTION.md (Conway-Jones TRIV/ROT3/PENT and the general Galois " +
               "criterion Theorem D cites) + " +
               "docs/proofs/PROOF_K_PARTNERSHIP.md (R psi_k = (-1)^(k+1) psi_k, Lemma B's only input beyond " +
               "arithmetic; the sign itself is F71's) + " +
               "docs/proofs/PROOF_ZETA2_ANTI_PROTECTION.md (the all-orders route for the eleven Theta-mirror pairs) " +
               "+ " +
               "docs/proofs/PROOF_MIRROR_ORDER_SORTING.md (F131, the repo's home for which orders may enter a " +
               "response; its carrier is an involution and grades by Z/2 at most, which is why it reaches the " +
               "eleven and not the 212) + " +
               "docs/ANALYTICAL_FORMULAS.md (F161; F160 the road, F129 the law under test, F65 the rate comb that " +
               "is the first order, F2b the chain comb) + " +
               "simulations/collision_gap_odd_orders.py (the gate, fifteen blocks L1 to L14)")
    {
        CrackedRing = crackedRing ?? throw new ArgumentNullException(nameof(crackedRing));
    }

    /// <summary>The smallest COMB modulus this claim speaks about: F129's smallest firing modulus. Below it there is
    /// no collision to expand. The series itself (<see cref="LevelCoefficient"/>) holds at every n ≥ 4.</summary>
    public const int MinFiringComb = 9;

    /// <summary>The census range of every count on the proof's page, and of <see cref="Census"/>'s default.</summary>
    public const int CensusMaxComb = 30;

    // ---------------------------------------------------------------- the two ladders, and the gcd

    /// <summary>The multiplier the ODD orders read the comb under at rung j: X_2j = −M_{n+2j}, so the multiplier is
    /// n + 2j. The shift by n is the reflection sign itself, η_k = −cos(k·nπ/n).</summary>
    public static int OddOrderMultiplier(int nComb, int rungJ)
    {
        CheckComb(nComb);
        CheckRung(rungJ);
        return nComb + 2 * rungJ;
    }

    /// <summary>The multiplier the EVEN orders read the comb under at rung j: M_{2j+1}, odd already.</summary>
    public static int EvenOrderMultiplier(int rungJ)
    {
        CheckRung(rungJ);
        return 2 * rungJ + 1;
    }

    /// <summary>Whether k ↦ k·multiplier mod 2n induces a Galois automorphism of ℚ(ζ_2n), i.e. gcd(m, 2n) = 1. This
    /// is the GLOBAL criterion: sufficient for the kill and, by Corollary G, too crude (the sharp one is per
    /// piece).</summary>
    public static bool IsGaloisMultiplier(int nComb, int multiplier)
    {
        CheckComb(nComb);
        return Gcd(multiplier, 2L * nComb) == 1;
    }

    /// <summary>Whether the ODD-order rung j is a Galois kill, gcd(n + 2j, 2n) = 1. False at j = 0 for every n
    /// (gcd(n, 2n) = n), which is why the first order survives; false at every rung when n is even.</summary>
    public static bool OddOrderRungIsGaloisKill(int nComb, int rungJ) =>
        IsGaloisMultiplier(nComb, OddOrderMultiplier(nComb, rungJ));

    /// <summary>Whether the EVEN-order rung j is a Galois kill, gcd(2j + 1, 2n) = gcd(2j + 1, n) = 1.</summary>
    public static bool EvenOrderRungIsGaloisKill(int nComb, int rungJ) =>
        IsGaloisMultiplier(nComb, EvenOrderMultiplier(rungJ));

    /// <summary>The odd-n reduction of the X ladder's criterion, returned as both sides so a caller can compare them:
    /// (gcd(n + 2j, 2n), gcd(j, n)). They are EQUAL at every odd n, because n + 2j is then odd so the factor 2 never
    /// contributes and gcd(n + 2j, n) = gcd(2j, n) = gcd(j, n). At even n they differ in general.</summary>
    public static (long Full, long Reduced) OddCombGcdIdentity(int nComb, int rungJ) =>
        (Gcd(OddOrderMultiplier(nComb, rungJ), 2L * nComb), Gcd(rungJ, nComb));

    /// <summary>The first rung j ≥ 1 of the ODD-order ladder that is NOT a Galois kill, so the first that can leave a
    /// coefficient standing. At every odd firing modulus this is 3, because F129 forces 3 | n there; at even n it is
    /// 1, the ladder being unreachable from the collision altogether (which is what Theorem E has to replace).</summary>
    public static int FirstSurvivingOddOrderRung(int nComb)
    {
        CheckComb(nComb);
        for (int j = 1; ; j++)
            if (!OddOrderRungIsGaloisKill(nComb, j)) return j;
    }

    // ---------------------------------------------------------------- the comb readings, exact in Z[zeta_2n]

    /// <summary>2·M_m(τ) = Σ_{k∈τ} (ζ^{mk} + ζ^{−mk}) ∈ ℤ[ζ_2n], the comb read under the multiplier m. Doubled to
    /// stay on the integer free basis; vector equality IS equality over ℚ(ζ), so no tolerance arises. At m = 1 this is
    /// <see cref="LevelCollisionCensus.LevelVector"/>.</summary>
    public static long[] TwoTimesM(int nComb, (int K1, int K2, int K3) triple, int multiplier)
    {
        CheckTriple(nComb, triple);
        int m = 2 * nComb;
        var acc = CyclotomicRing.Zero(m);
        foreach (int k in stackalloc[] { triple.K1, triple.K2, triple.K3 })
        {
            CyclotomicRing.AddRootPower(acc, m, (long)multiplier * k, 1);
            CyclotomicRing.AddRootPower(acc, m, -(long)multiplier * k, 1);
        }
        return acc;
    }

    /// <summary>2·X_2j(τ) = Σ_{k∈τ} η_k·(ζ^{2jk} + ζ^{−2jk}) ∈ ℤ[ζ_2n], η_k = (−1)^{k+1}: the comb read with the
    /// reflection sign attached. Equal to −2·M_{n+2j}(τ) as a vector, which is Corollary C's identity.</summary>
    public static long[] TwoTimesX(int nComb, (int K1, int K2, int K3) triple, int rungJ)
    {
        CheckTriple(nComb, triple);
        CheckRung(rungJ);
        int m = 2 * nComb;
        var acc = CyclotomicRing.Zero(m);
        foreach (int k in stackalloc[] { triple.K1, triple.K2, triple.K3 })
        {
            long sign = k % 2 == 1 ? 1 : -1;
            CyclotomicRing.AddRootPower(acc, m, 2L * rungJ * k, sign);
            CyclotomicRing.AddRootPower(acc, m, -2L * rungJ * k, sign);
        }
        return acc;
    }

    /// <summary>2·ΔX_2j = 2·(X_2j(τ) − X_2j(σ)), the odd orders' input.</summary>
    public static long[] DeltaTwoTimesX(int nComb, (int K1, int K2, int K3) tau, (int K1, int K2, int K3) sigma,
                                        int rungJ) =>
        CyclotomicRing.Add(TwoTimesX(nComb, tau, rungJ), CyclotomicRing.Negate(TwoTimesX(nComb, sigma, rungJ)));

    /// <summary>2·ΔM_m = 2·(M_m(τ) − M_m(σ)), the even orders' input. ΔM_1 = 0 IS the collision.</summary>
    public static long[] DeltaTwoTimesM(int nComb, (int K1, int K2, int K3) tau, (int K1, int K2, int K3) sigma,
                                        int multiplier) =>
        CyclotomicRing.Add(TwoTimesM(nComb, tau, multiplier), CyclotomicRing.Negate(TwoTimesM(nComb, sigma, multiplier)));

    /// <summary>Whether the two triples collide, decided exactly in ℤ[ζ_2n] (F129's level equality, ΔM_1 = 0).</summary>
    public static bool Collide(int nComb, (int K1, int K2, int K3) tau, (int K1, int K2, int K3) sigma)
    {
        CheckTriple(nComb, tau);
        CheckTriple(nComb, sigma);
        return CyclotomicRing.AreEqual(
            LevelCollisionCensus.LevelVector(nComb, tau.K1, tau.K2, tau.K3),
            LevelCollisionCensus.LevelVector(nComb, sigma.K1, sigma.K2, sigma.K3));
    }

    /// <summary>Every collision pair at the comb modulus n: distinct CLEAN triples with equal level vector, each
    /// unordered pair yielded once, in ascending triple order. Recomputed exactly, never stored. Cost is
    /// <see cref="LevelCollisionCensus.CleanTriples"/>'s Θ(n³) plus the pairs within each bucket.</summary>
    public static IEnumerable<((int K1, int K2, int K3) Tau, (int K1, int K2, int K3) Sigma)> CollisionPairs(int nComb)
    {
        CheckComb(nComb);
        var buckets = new Dictionary<string, List<(int K1, int K2, int K3)>>();
        foreach (var t in LevelCollisionCensus.CleanTriples(nComb))
        {
            string key = CyclotomicRing.Key(LevelCollisionCensus.LevelVector(nComb, t.K1, t.K2, t.K3));
            if (!buckets.TryGetValue(key, out var list)) buckets[key] = list = new List<(int, int, int)>();
            list.Add(t);
        }
        foreach (var list in buckets.Values)
            for (int i = 0; i < list.Count; i++)
                for (int j = i + 1; j < list.Count; j++)
                    yield return (list[i], list[j]);
    }

    // ---------------------------------------------------------------- what each order does, per pair

    /// <summary>The number of ODD labels in the triple. ΔX_0 = 2(o_τ − o_σ), because Σ_τ η = 2o − 3.</summary>
    public static int OddLabelCount((int K1, int K2, int K3) triple) =>
        (triple.K1 & 1) + (triple.K2 & 1) + (triple.K3 & 1);

    /// <summary>Whether ΔX_2 = 0, which is what makes <see cref="FirstOrderCoefficient"/> the odd-label count. True
    /// for every collision pair at odd n (Theorem D) and for the ROT3 and Θ-mirror pairs at even n.</summary>
    public static bool FirstOrderIsRational(int nComb, (int K1, int K2, int K3) tau, (int K1, int K2, int K3) sigma) =>
        CyclotomicRing.IsZero(DeltaTwoTimesX(nComb, tau, sigma, 1));

    /// <summary>c_1 = (2/n)(ΔX_0 − ΔX_2), returned as an exact rational ONLY where ΔX_2 vanishes, where it reads
    /// (4/n)(o_τ − o_σ). Elsewhere ΔX_2 is a cyclotomic integer with no rational value and the call throws rather
    /// than rounding one; check first with <see cref="FirstOrderIsRational"/> and otherwise read
    /// <see cref="DeltaTwoTimesX"/> directly.</summary>
    public static BigRational FirstOrderCoefficient(int nComb, (int K1, int K2, int K3) tau,
                                                    (int K1, int K2, int K3) sigma)
    {
        if (!FirstOrderIsRational(nComb, tau, sigma))
            throw new ArgumentException(
                $"DeltaX_2 does not vanish at n = {nComb} for these triples, so c_1 is not (4/n)(o_tau - o_sigma) " +
                "but (2/n)(DeltaX_0 - DeltaX_2) with an irrational second term. Read DeltaTwoTimesX instead.",
                nameof(tau));
        return new BigRational(4L * (OddLabelCount(tau) - OddLabelCount(sigma)), nComb);
    }

    /// <summary>Whether the pair STANDS at first order, c_1 = 0. At odd n that is o_τ = o_σ; in general it is
    /// ΔX_0 = ΔX_2 as ring elements, which is what is decided here.</summary>
    public static bool StandsAtFirstOrder(int nComb, (int K1, int K2, int K3) tau, (int K1, int K2, int K3) sigma) =>
        CyclotomicRing.AreEqual(DeltaTwoTimesX(nComb, tau, sigma, 0), DeltaTwoTimesX(nComb, tau, sigma, 1));

    /// <summary>Whether c_2 = −((n−3)/n²)·ΔM_3 vanishes, decided exactly. Theorem D gives it for free when 3 ∤ n;
    /// Corollary G sharpens that to a per-piece criterion and the census identifies the 60 as family C.</summary>
    public static bool SecondOrderVanishes(int nComb, (int K1, int K2, int K3) tau, (int K1, int K2, int K3) sigma) =>
        CyclotomicRing.IsZero(DeltaTwoTimesM(nComb, tau, sigma, 3));

    /// <summary>Whether c_3 ∝ ΔX_2 − ΔX_4 vanishes. A THEOREM at every odd n (Theorem D kills both rungs
    /// separately); at even n it holds for the ROT3 and Θ-mirror pairs and is decided here rather than
    /// assumed.</summary>
    public static bool ThirdOrderVanishes(int nComb, (int K1, int K2, int K3) tau, (int K1, int K2, int K3) sigma) =>
        CyclotomicRing.AreEqual(DeltaTwoTimesX(nComb, tau, sigma, 1), DeltaTwoTimesX(nComb, tau, sigma, 2));

    /// <summary>Whether the ODD part of the gap begins exactly at fifth order. The condition is the proof's, §(g):
    /// the pair stands at first order AND ΔX_2 = ΔX_4 = 0 AND ΔX_6 ≠ 0. Both middle clauses are needed, and
    /// <see cref="ThirdOrderVanishes"/> is not enough for them: c_5 carries ΔX_2 and ΔX_4 with their own weights, so
    /// only their VANISHING (not their equality, which is all c_3 = 0 says) reduces c_5 to
    /// −(3/4)(2n−3)(n−3)·F·ΔX_6. Every standing pair of the n ≤ 30 census satisfies the stronger form, so the
    /// distinction changes no count there; it is the difference between the predicate and the statement it licenses.
    /// False for the Θ-mirror pairs, whose odd part vanishes at every order.</summary>
    public static bool OddPartStartsAtFifthOrder(int nComb, (int K1, int K2, int K3) tau,
                                                  (int K1, int K2, int K3) sigma) =>
        StandsAtFirstOrder(nComb, tau, sigma)
        && CyclotomicRing.IsZero(DeltaTwoTimesX(nComb, tau, sigma, 1))
        && CyclotomicRing.IsZero(DeltaTwoTimesX(nComb, tau, sigma, 2))
        && !CyclotomicRing.IsZero(DeltaTwoTimesX(nComb, tau, sigma, 3));

    /// <summary>Whether σ = n − τ, the Θ-mirror shape, whose gap is even in u at ALL orders (K·H(u)·K = −H(−u) at
    /// even n; the ladder alone reaches only through the fifth).</summary>
    public static bool IsThetaMirrorPair(int nComb, (int K1, int K2, int K3) tau, (int K1, int K2, int K3) sigma)
    {
        CheckTriple(nComb, tau);
        CheckTriple(nComb, sigma);
        var mirrored = new[] { nComb - tau.K1, nComb - tau.K2, nComb - tau.K3 };
        Array.Sort(mirrored);
        return mirrored[0] == sigma.K1 && mirrored[1] == sigma.K2 && mirrored[2] == sigma.K3;
    }

    // ---------------------------------------------------------------- the ROT3 shape

    /// <summary>Whether every label of the triple has the same parity, so η is constant on it. At EVEN n this comes
    /// free with the ROT3 shape (6 | n makes every coset parity-homogeneous); at odd n it is a real
    /// hypothesis.</summary>
    public static bool IsParityUniform((int K1, int K2, int K3) triple) =>
        ((triple.K1 ^ triple.K2) & 1) == 0 && ((triple.K1 ^ triple.K3) & 1) == 0;

    /// <summary>The triple's DOUBLED-LABEL ROT3 cosets, or null if it does not have the shape: the ±label set
    /// P = {±k mod n} must have six elements and be the union of two cosets of the order-3 subgroup ⟨n/3⟩ ≤ ℤ/n
    /// (which are then necessarily C and −C). Conway-Jones' ROT3 read on the doubled labels; NOT
    /// <c>MirrorWorld.Seed</c>'s <c>TripleFamily.Rot3</c>, which is the same shape on the narrower domain of triples
    /// whose own level sum vanishes. Each coset ascending, the two ordered by first element.</summary>
    public static (int[] First, int[] Second)? DoubledLabelRot3Cosets(int nComb, (int K1, int K2, int K3) triple)
    {
        CheckTriple(nComb, triple);
        if (nComb % 3 != 0) return null;
        int step = nComb / 3;
        var labels = new SortedSet<int>();
        foreach (int k in stackalloc[] { triple.K1, triple.K2, triple.K3 })
        {
            labels.Add(Mod(k, nComb));
            labels.Add(Mod(-k, nComb));
        }
        if (labels.Count != 6) return null;
        var seen = new HashSet<string>();
        var kept = new List<int[]>();
        foreach (int a in labels)
        {
            var coset = new[] { Mod(a, nComb), Mod(a + step, nComb), Mod(a + 2 * step, nComb) };
            Array.Sort(coset);
            if (seen.Add(string.Join(",", coset))) kept.Add(coset);
        }
        if (kept.Count != 2) return null;
        foreach (var coset in kept)
            foreach (int e in coset)
                if (!labels.Contains(e)) return null;
        kept.Sort((a, b) => a[0].CompareTo(b[0]));
        return (kept[0], kept[1]);
    }

    /// <summary>Whether the ROT3 rung lemma FORCES X_2j(τ) = 0: the triple is parity-uniform and doubled-label ROT3
    /// and 3 ∤ j. The converse is false: at 3 | j the coset collapses to a cosine that may or may not vanish, and
    /// n = 24 with τ = (1,7,9) exhibits both outcomes (j = 3 gives cos(π/4) ≠ 0, j = 6 gives cos(π/2) = 0).</summary>
    public static bool Rot3ForcesVanishing(int nComb, (int K1, int K2, int K3) triple, int rungJ)
    {
        CheckRung(rungJ);
        if (rungJ % 3 == 0) return false;
        return IsParityUniform(triple) && DoubledLabelRot3Cosets(nComb, triple) is not null;
    }

    /// <summary>The collapse value at 3 | j, as 2·X_2j(τ) in ℤ[ζ_2n]: each coset collapses to one point taken three
    /// times, so X_2j(τ) = 3·η_τ·cos(2πja/n) with a the coset label, and twice it is
    /// 3η_τ(ζ_2n^{2ja} + ζ_2n^{−2ja}). Built from that closed form alone, so comparing it with
    /// <see cref="TwoTimesX"/> is two routes to one number.</summary>
    public static long[] Rot3CollapseValue(int nComb, (int K1, int K2, int K3) triple, int rungJ)
    {
        CheckRung(rungJ);
        if (rungJ % 3 != 0)
            throw new ArgumentOutOfRangeException(nameof(rungJ), rungJ, "the coset collapses only at 3 | j");
        var cosets = DoubledLabelRot3Cosets(nComb, triple)
                     ?? throw new ArgumentException("the triple is not doubled-label ROT3", nameof(triple));
        if (!IsParityUniform(triple))
            throw new ArgumentException("the collapse formula needs a parity-uniform triple", nameof(triple));
        long eta = triple.K1 % 2 == 1 ? 1 : -1;
        int m = 2 * nComb;
        var acc = CyclotomicRing.Zero(m);
        // One representative per coset: at 3 | j every element of a coset has the same image, and the two cosets are
        // C and -C, so the two terms are the conjugate pair 3 eta (zeta^{2ja} + zeta^{-2ja}).
        foreach (int a in new[] { cosets.First[0], cosets.Second[0] })
            CyclotomicRing.AddRootPower(acc, m, 2L * rungJ * a, 3 * eta);
        return acc;
    }

    // ---------------------------------------------------------------- the coefficients, in closed form

    /// <summary>The orders Theorem A computes and this claim carries: 1, 2, 3 and 5. d_4 exists in the proof, carries
    /// no η and plays no part in an odd order; it is deliberately not carried here, and asking for it throws rather
    /// than returning a number this claim does not hold.</summary>
    public static IReadOnlyList<int> CarriedExpansionOrders { get; } = new[] { 1, 2, 3, 5 };

    /// <summary>The multipliers r that the order m reads the comb under, so that
    /// d_m = η^m · <see cref="MultiplierFormPrefactor"/>(n, m) · Σ_r <see cref="MultiplierFormWeight"/>(n, m, r)·cos rθ.
    /// This is Corollary C, and it is the claim's actual content: a difference of two NEIGHBOURING rungs through third
    /// order (r = 0,2 then 1,3 then 2,4) and three rungs at fifth (r = 2,4,6, with the constant rung r = 0 present but
    /// weighted zero). Every r satisfies r ≡ m + 1 (mod 2), which is the parity theorem of §(d): the EVEN orders read
    /// odd multipliers (the M ladder) and the ODD orders even ones (the X ladder, r = 2j).</summary>
    public static IReadOnlyList<int> MultiplierFormRungs(int expansionOrder) => expansionOrder switch
    {
        1 => new[] { 0, 2 },
        2 => new[] { 1, 3 },
        3 => new[] { 2, 4 },
        5 => new[] { 0, 2, 4, 6 },
        4 => throw new ArgumentOutOfRangeException(nameof(expansionOrder), expansionOrder,
                 "d_4 is computed in the proof and deliberately not carried: it has no eta and enters no odd order"),
        _ => throw new ArgumentOutOfRangeException(nameof(expansionOrder), expansionOrder,
                 "Theorem A reaches order 5; the ladder's general shape past it is open"),
    };

    /// <summary>The order's overall rational prefactor: 2/n, (n−3)/n², 2(n−2)(n−4)/(3n³) and, at fifth order,
    /// F = 4(n−6)(n−2)/(15n^5), which vanishes at n = 6, outside every firing modulus.</summary>
    public static BigRational MultiplierFormPrefactor(int nComb, int expansionOrder)
    {
        CheckComb(nComb);
        _ = MultiplierFormRungs(expansionOrder);
        long n = nComb;
        return expansionOrder switch
        {
            1 => new BigRational(2, n),
            2 => new BigRational(n - 3, n * n),
            3 => new BigRational(2 * (n - 2) * (n - 4), 3 * n * n * n),
            _ => new BigRational(4 * (n - 6) * (n - 2), 15 * BigInteger.Pow(nComb, 5)),
        };
    }

    /// <summary>The order's weight on the rung cos rθ, exactly. Through third order the two weights are +1 and −1;
    /// at fifth they are 0, (5/4)(n−1), (3n²/2 − 8n + 8) and −(3/4)(2n−3)(n−3) on r = 0, 2, 4, 6. In every carried
    /// order the weights SUM TO ZERO, which is d_m → 0 as θ → 0, where the mode has no endpoint amplitude; and the
    /// fifth order's r = 6 weight is nonzero at every firing modulus (its factors are positive for n ≥ 9), which is
    /// why X_6 first enters there. A rung the order does not carry weighs zero.</summary>
    public static BigRational MultiplierFormWeight(int nComb, int expansionOrder, int rung)
    {
        CheckComb(nComb);
        var rungs = MultiplierFormRungs(expansionOrder);
        if (!rungs.Contains(rung)) return BigRational.Zero;
        long n = nComb;
        if (expansionOrder != 5) return rung == rungs[0] ? BigRational.One : -BigRational.One;
        return rung switch
        {
            0 => BigRational.Zero,
            2 => new BigRational(5 * (n - 1), 4),
            4 => new BigRational(3 * n * n - 16 * n + 16, 2),
            _ => new BigRational(-3 * (2 * n - 3) * (n - 3), 4),
        };
    }

    /// <summary>d_m, the u^m coefficient of the level E_k(u) at the chain end, in units of J: assembled from the
    /// multiplier form above, so the closed form lives in exactly one place. n = nComb = N + 1, k the mode label
    /// 1..n−1. d_1 is F160's Theorem G, and its magnitude is F65's endpoint rate comb. Doubles, because cos(rkπ/n) is
    /// not rational in general; where it is, <see cref="LevelCoefficientExact"/> gives the same number over ℚ.</summary>
    public static double LevelCoefficient(int nComb, int k, int expansionOrder)
    {
        CheckComb(nComb);
        if (k < 1 || k >= nComb) throw new ArgumentOutOfRangeException(nameof(k), k, "the comb has k = 1..n-1");
        double theta = k * Math.PI / nComb;
        double sum = 0;
        foreach (int r in MultiplierFormRungs(expansionOrder))
            sum += ToDouble(MultiplierFormWeight(nComb, expansionOrder, r)) * Math.Cos(r * theta);
        double eta = (expansionOrder % 2 == 1 && k % 2 == 0) ? -1.0 : 1.0;
        return eta * ToDouble(MultiplierFormPrefactor(nComb, expansionOrder)) * sum;
    }

    /// <summary>d_m over ℚ, for the comb points where cos θ_k is rational (k = n/3, n/2, 2n/3 give ½, 0, −½) and for
    /// any other exact cosine a caller can supply. cos rθ is reached from cos θ by the Chebyshev recursion, so the
    /// whole evaluation stays in ℚ and can be compared for equality rather than tolerated.</summary>
    public static BigRational LevelCoefficientExact(int nComb, int k, int expansionOrder, BigRational cosTheta)
    {
        CheckComb(nComb);
        if (k < 1 || k >= nComb) throw new ArgumentOutOfRangeException(nameof(k), k, "the comb has k = 1..n-1");
        var sum = BigRational.Zero;
        foreach (int r in MultiplierFormRungs(expansionOrder))
            sum += MultiplierFormWeight(nComb, expansionOrder, r) * ChebyshevFirstKind(r, cosTheta);
        var eta = (expansionOrder % 2 == 1 && k % 2 == 0) ? -BigRational.One : BigRational.One;
        return eta * MultiplierFormPrefactor(nComb, expansionOrder) * sum;
    }

    /// <summary>T_r(c) = cos(r·arccos c), by the recursion T_0 = 1, T_1 = c, T_r = 2c·T_{r−1} − T_{r−2}.</summary>
    private static BigRational ChebyshevFirstKind(int r, BigRational c)
    {
        var prev = BigRational.One;
        if (r == 0) return prev;
        var cur = c;
        for (int i = 2; i <= r; i++) (prev, cur) = (cur, 2 * c * cur - prev);
        return cur;
    }

    private static double ToDouble(BigRational r) => (double)r.Numerator / (double)r.Denominator;

    // ---------------------------------------------------------------- the second order past the census

    /// <summary>A LOWER BOUND on the number of collision pairs with c_2 = 0 at the comb modulus n, and the census's
    /// exact count at n ≤ 30: the two families of F129's inventory whose every piece has a RATIO-ORDER coprime to 3,
    /// C (zero mode plus two R₅ pieces, door 10 | n) and L (zero mode plus one piece of ratio-order 70, door 70 | n).
    /// Read straight off <see cref="CollisionFamilyInventory.Count"/>, doors included, rather than restated: that is
    /// 2(n−10)·[10|n] + 20·[70|n], giving 20 at n = 20, 40 at n = 30, 140 at n = 70 and 420 at n = 210. Only the
    /// PROVED direction is in it; the converse, which would make it the count, is checked at n ≤ 30 and nowhere
    /// else.</summary>
    public static long SecondOrderZeroPairLowerBound(int nComb)
    {
        CheckComb(nComb);
        return CollisionFamilyInventory.Count(CollisionFamilyInventory.Family.C, nComb)
             + CollisionFamilyInventory.Count(CollisionFamilyInventory.Family.L, nComb);
    }

    /// <summary>The families of F129's inventory that are 3-free, so that Corollary G's criterion applies to them
    /// unchanged: C and L, and no others, read off the committed table by its own ratio-orders.</summary>
    public static IReadOnlyList<CollisionFamilyInventory.Family> ThreeFreeFamilies { get; } =
        new[] { CollisionFamilyInventory.Family.C, CollisionFamilyInventory.Family.L };

    // ---------------------------------------------------------------- the census, recomputed

    /// <summary>What <see cref="Census"/> recomputes. Every field is a count of pairs, every one decided in ℤ[ζ_2n];
    /// nothing is stored or read back.</summary>
    public sealed record CensusReport(
        int MaxComb,
        int Pairs,
        int Separating,
        int Standing,
        int OddCombPairs,
        int OddCombThirdOrderVanishes,
        int ThetaMirrorStanding,
        int NonMirrorStanding,
        int NonMirrorStandingWithSixthRungAlive,
        int EvenCombStandingBothRot3AndParityMatched,
        int SecondOrderVanishes);

    /// <summary>The census the proof reports, recomputed from the triples: every collision pair at every firing comb
    /// modulus up to <paramref name="maxComb"/>, classified by which orders vanish. At the default bound this is the
    /// nine firing moduli 9..30.</summary>
    public static CensusReport Census(int maxComb = CensusMaxComb)
    {
        if (maxComb < MinFiringComb)
            throw new ArgumentOutOfRangeException(nameof(maxComb), maxComb,
                $"the first firing modulus is {MinFiringComb}");
        int pairs = 0, standing = 0, oddPairs = 0, oddThird = 0, mirrorStanding = 0,
            nonMirrorStanding = 0, sixthAlive = 0, rot3Matched = 0, secondZero = 0;
        for (int n = MinFiringComb; n <= maxComb; n++)
        {
            if (!LevelCollisionCensus.Fires(n)) continue;
            bool odd = n % 2 == 1;
            foreach (var (tau, sigma) in CollisionPairs(n))
            {
                pairs++;
                if (odd)
                {
                    oddPairs++;
                    if (ThirdOrderVanishes(n, tau, sigma)) oddThird++;
                }
                if (SecondOrderVanishes(n, tau, sigma)) secondZero++;
                if (!StandsAtFirstOrder(n, tau, sigma)) continue;
                standing++;
                if (IsThetaMirrorPair(n, tau, sigma)) { mirrorStanding++; continue; }
                nonMirrorStanding++;
                if (OddPartStartsAtFifthOrder(n, tau, sigma)) sixthAlive++;
                if (!odd
                    && IsParityUniform(tau) && IsParityUniform(sigma) && (tau.K1 % 2) == (sigma.K1 % 2)
                    && DoubledLabelRot3Cosets(n, tau) is not null && DoubledLabelRot3Cosets(n, sigma) is not null)
                    rot3Matched++;
            }
        }
        return new CensusReport(maxComb, pairs, pairs - standing, standing, oddPairs, oddThird,
            mirrorStanding, nonMirrorStanding, sixthAlive, rot3Matched, secondZero);
    }

    // ---------------------------------------------------------------- guards

    private static void CheckComb(int nComb)
    {
        if (nComb < 4) throw new ArgumentOutOfRangeException(nameof(nComb), nComb,
            "the comb modulus is n = N + 1 with N >= 3 sites");
    }

    private static void CheckRung(int rungJ)
    {
        if (rungJ < 0) throw new ArgumentOutOfRangeException(nameof(rungJ), rungJ, "a rung is j >= 0");
    }

    private static void CheckTriple(int nComb, (int K1, int K2, int K3) t)
    {
        CheckComb(nComb);
        if (!(0 < t.K1 && t.K1 < t.K2 && t.K2 < t.K3 && t.K3 < nComb))
            throw new ArgumentOutOfRangeException(nameof(t), t, "a triple is 1 <= k1 < k2 < k3 <= n-1");
    }

    private static int Mod(int a, int m) => ((a % m) + m) % m;

    private static long Gcd(long a, long b)
    {
        a = a < 0 ? -a : a;
        b = b < 0 ? -b : b;
        while (b != 0) (a, b) = (b, a % b);
        return a;
    }

    public override string DisplayName =>
        "F161: each order of a collision gap reads one comb under a multiplier, and which orders vanish is a gcd";

    public override string Summary =>
        "the cracked ring's road expands a chain level in u; the even orders read the comb at odd multipliers and " +
        "the odd orders at n + 2j, so a collision, being a rational identity among roots of unity, is invisible to " +
        "every multiplier that merely permutes the comb: at odd n that empties the third order for every pair, and " +
        "the first multiplier that fails is 3, the prime F129 needed for the coincidence to exist at all " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the five coefficients, and the two ladders",
                summary: "Solving det(x I - H(u)) = 0 order by order at a comb root gives d_1 = eta(4/n)s^2, " +
                         "d_2 = (4(n-3)/n^2)s^2 c, d_3 = eta(4(n-2)(n-4)/(3n^3))s^2(4c^2-1) and d_5 = eta F[...], " +
                         "with d_4 computed and not carried. In the cosine basis each is a signed combination of " +
                         "NEIGHBOURING evaluations of one comb, a difference of two through third order and three " +
                         "rungs at fifth. The even orders sit on M with odd multipliers, the odd orders on X with " +
                         "even ones, and X is the comb shifted by n because the reflection sign is itself a " +
                         "harmonic, eta_k = -cos(k n pi/n). n = 12, order by order as (rungs r of cos r theta; " +
                         "weights; prefactor): " +
                         string.Join("   ", CarriedExpansionOrders.Select(m =>
                             $"d_{m} ({string.Join(",", MultiplierFormRungs(m))}; " +
                             $"{string.Join(",", MultiplierFormRungs(m).Select(r => MultiplierFormWeight(12, m, r)))}; " +
                             $"{MultiplierFormPrefactor(12, m)})")) +
                         ". Every order's weights sum to zero, which is d_m -> 0 as theta -> 0, and every rung r " +
                         "satisfies r = m + 1 mod 2, which is the parity theorem.");

            yield return new InspectableNode("the gcd that decides a rung",
                summary: "k -> k m mod 2n is a Galois automorphism of Q(zeta_2n) exactly when gcd(m, 2n) = 1, and " +
                         "an automorphism sends the collision's Delta M_1 = 0 to zero again. On the X ladder " +
                         "m = n + 2j, and at ODD n that is odd, so the factor 2 never contributes and the criterion " +
                         "collapses to gcd(j, n) = 1. n = 9, as (gcd(n+2j, 2n), gcd(j, n)): rung 1 " +
                         $"{OddCombGcdIdentity(9, 1)}, rung 2 {OddCombGcdIdentity(9, 2)}, rung 3 " +
                         $"{OddCombGcdIdentity(9, 3)}; the first surviving rung is j = " +
                         $"{FirstSurvivingOddOrderRung(9)}. At n = 20, where 3 does not divide n, the EVEN ladder's " +
                         $"rung 1 (multiplier 3) is a kill: {EvenOrderRungIsGaloisKill(20, 1)}, and that is the " +
                         "second order gone for every collision pair there. The rung j = 0 is never an automorphism, " +
                         "gcd(n, 2n) = n, and it is exactly the rung that leaves the first order standing.");

            yield return new InspectableNode("where the ladder stops is F129's own divisor",
                summary: "F129 fires only at 3|n >= 9 or 10|n >= 20; 10|n forces n even, so every ODD firing " +
                         "modulus has 3|n. Hence at odd n the first rung that is not an automorphism is 3, and at " +
                         "even n the first rung that collapses a ROT3 coset is 3 as well. Either way X_6 is the " +
                         "first rung that can survive, and by the multiplier form it first enters at FIFTH order. " +
                         "The same arithmetic that lets the coincidence exist at an odd modulus is the arithmetic " +
                         "that ends its protection. First surviving rung at the odd firing moduli 9, 15, 21, 27: " +
                         $"{string.Join(", ", new[] { 9, 15, 21, 27 }.Select(FirstSurvivingOddOrderRung))}.");

            yield return new InspectableNode("the ROT3 shape, and its one forced direction",
                summary: "At even n the collision reaches no rung of the X ladder at all, and what carries the " +
                         "twelve non-mirror standing pairs is the SHAPE of the triples: a parity-uniform triple " +
                         "whose doubled-label set is a union of two order-3 cosets has X_2j = 0 forced at every rung " +
                         "with 3 not dividing j, a full coset of a nontrivial cyclic group summing to zero. At 3 | j " +
                         "the coset collapses to 3 eta cos(2 pi j a/n), which need not vanish. BREAK-INPUT, n = 24, " +
                         $"tau = (1,7,9): ROT3 = {DoubledLabelRot3Cosets(24, (1, 7, 9)) is not null}, forced at " +
                         $"j = 1: {Rot3ForcesVanishing(24, (1, 7, 9), 1)}, at j = 3: " +
                         $"{Rot3ForcesVanishing(24, (1, 7, 9), 3)} with the collapse NONZERO there " +
                         $"({!CyclotomicRing.IsZero(TwoTimesX(24, (1, 7, 9), 3))}) and the same collapse vanishing " +
                         $"at j = 6 ({CyclotomicRing.IsZero(TwoTimesX(24, (1, 7, 9), 6))}). A reading of the lemma " +
                         "as an equivalence is wrong, and n = 30 alone could not have shown it.");

            yield return new InspectableNode("the second order, and the piece criterion",
                summary: "c_2 = -((n-3)/n^2) Delta M_3, so it vanishes with Delta M_3. The global criterion " +
                         "gcd(3, 2n) = 1 gives it for the 10|n family with 3 not dividing n, smallest member " +
                         "n = 20. It is too crude: what a multiplier must be invertible on is every minimal PIECE of " +
                         "the vanishing sum and its RATIO-ORDER, not the whole ring, and under that reading the 40 " +
                         "further pairs at n = 30 are as forced as the 20 at n = 20. Inside the census they are " +
                         "F129's family C; family L is 3-free as well, so the bound past the census is " +
                         "Count(C, n) + Count(L, n) off the committed inventory: n = 20 gives " +
                         $"{SecondOrderZeroPairLowerBound(20)}, n = 30 gives {SecondOrderZeroPairLowerBound(30)}, " +
                         $"n = 70 gives {SecondOrderZeroPairLowerBound(70)} and n = 210 gives " +
                         $"{SecondOrderZeroPairLowerBound(210)}. The proved direction only; the converse is checked " +
                         "at n <= 30 and nowhere else.");

            yield return new InspectableNode("the census, recomputed here", summary: Describe(Census()));

            yield return new InspectableNode("scope, and what is read rather than derived",
                summary: "Five orders, not a general form: the PARITY of the multipliers is a theorem at every " +
                         "order (two symmetries of the angle equation), the RANGE is not. Theorem E has one forced " +
                         "direction. The ROT3 shape of the twelve even-n standing pairs is READ at n <= 30; the " +
                         "consequence is a theorem, the hypothesis a census observation, and whether a standing pair " +
                         "at even n must be ROT3 stays open. Every count is at n <= 30 except those read off the " +
                         "committed family table. The expansion point is u = 0, the CHAIN end, unlike F160's " +
                         "Theorem E, which expands about the ring end and is about a pair split on a different comb. " +
                         "No hardware claim: Confirmation 24 is the n = 9 pair seen on another axis, and nothing " +
                         "here is confirmed by it.");
        }
    }

    private static string Describe(CensusReport r) =>
        $"firing moduli up to n = {r.MaxComb}: {r.Pairs} collision pairs, {r.Separating} separating at first order " +
        $"and {r.Standing} standing. Of the {r.OddCombPairs} at ODD n, {r.OddCombThirdOrderVanishes} have a " +
        $"vanishing third order, which by Theorem D is all of them. Of the standing pairs " +
        $"{r.ThetaMirrorStanding} are Theta-mirror (odd part zero at every order) and {r.NonMirrorStanding} are not, " +
        $"of which {r.NonMirrorStandingWithSixthRungAlive} have the sixth rung alive, so their odd part begins " +
        $"exactly at fifth order; {r.EvenCombStandingBothRot3AndParityMatched} of the even-n ones have both triples " +
        $"doubled-label ROT3 and parity-matched, which is Theorem E's hypothesis met. {r.SecondOrderVanishes} pairs " +
        "in all have a vanishing SECOND order.";
}
