using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below pins for the adopted collision gap (F161), the ladder the crack's road carries at its
// chain end. Sources: docs/proofs/PROOF_COLLISION_GAP_ODD_ORDERS.md (Theorems A, D, E, Corollaries
// C, G) and the census it explains, experiments/THE_COMB_ON_THE_ROAD.md.
//
// The literals below are that census's committed table and the proof's own gate populations, asserted
// BESIDE the closed forms, so a gate built only out of the formulas cannot pass by construction: the
// nine-column table comes from the experiment page, the ladder counts from the proof's L blocks, and
// the second-order bound is checked against LevelCollision's two family counts rather than against a
// copy of itself.
//
// The one-sidedness is LevelCollision's, and it points the useful way here: a comb reading that is
// NONZERO mod one prime is exactly nonzero, so every "separates" verdict is proof grade, and the
// vanishing is never measured but decided by a gcd (Theorem D) or by the ROT3 shape (Theorem E).
public class CollisionGapTests
{
    static CollisionGap Gap(int ncomb) => new(new Crack(new Cyclotomy(), ncomb - 1, 0));

    // the nine firing moduli of the census, and its committed table
    static readonly int[] Moduli = { 9, 12, 15, 18, 20, 21, 24, 27, 30 };
    static readonly int[] TablePairs = { 1, 25, 127, 162, 20, 255, 411, 244, 1313 };
    static readonly int[] TableSeparate = { 1, 24, 107, 160, 20, 195, 404, 124, 1300 };
    static readonly int[] TableStand = { 0, 1, 20, 2, 0, 60, 7, 120, 13 };

    // ---- the ladder: which comb each order reads, under which integer multiplier ----

    [Fact]
    public void TheTwoLaddersAreTheOddAndEvenMultipliers()
    {
        // odd orders read X_2j = -M_{n+2j}; even orders read M_{2j+1}. The letter n is the COMB
        // modulus, never the site count.
        Assert.Equal(11, CollisionGap.OddOrderMultiplier(9, 1));
        Assert.Equal(15, CollisionGap.OddOrderMultiplier(9, 3));
        Assert.Equal(9, CollisionGap.OddOrderMultiplier(9, 0));
        Assert.Equal(1, CollisionGap.EvenOrderMultiplier(0));
        Assert.Equal(3, CollisionGap.EvenOrderMultiplier(1));
    }

    // The criterion itself, asked with FREE multipliers. Its only caller passes n + 2j, and at odd n
    // gcd(n + 2j, n) and gcd(n + 2j, 2n) agree, so without these rows the modulus could be n rather
    // than 2n and every other gate in this file would stay green. The comb lives at 2n-th roots.
    [Fact]
    public void TheAutomorphismCriterionReadsTheTwoNthRootsAndNotTheNth()
    {
        Assert.False(CollisionGap.IsAutomorphism(9, 2));    // gcd(2, 18) = 2, while gcd(2, 9) = 1
        Assert.False(CollisionGap.IsAutomorphism(9, 4));    // gcd(4, 18) = 2, while gcd(4, 9) = 1
        Assert.False(CollisionGap.IsAutomorphism(9, 8));    // gcd(8, 18) = 2, while gcd(8, 9) = 1
        Assert.False(CollisionGap.IsAutomorphism(15, 8));   // gcd(8, 30) = 2, while gcd(8, 15) = 1
        Assert.False(CollisionGap.IsAutomorphism(21, 4));   // gcd(4, 42) = 2, while gcd(4, 21) = 1
        Assert.True(CollisionGap.IsAutomorphism(9, 5));     // gcd(5, 18) = 1
        Assert.True(CollisionGap.IsAutomorphism(9, 7));
        Assert.False(CollisionGap.IsAutomorphism(9, 9));    // gcd(9, 18) = 9: the j = 0 rung
        Assert.False(CollisionGap.IsAutomorphism(9, 6));
        Assert.False(CollisionGap.IsAutomorphism(12, 2));   // an even comb: nothing even is coprime to 24
        Assert.True(CollisionGap.IsAutomorphism(12, 5));
    }

    // Theorem D's criterion and its odd-n reduction, on the adopted range. The reduction is the
    // load-bearing half: it is what makes j = 3 the first survivor at every odd firing modulus.
    [Fact]
    public void TheGcdIdentityHoldsAtEveryOddComb()
    {
        bool sawKill = false, sawLive = false;
        for (int n = 9; n <= 59; n += 2)
            for (int j = 0; j <= 12; j++)
            {
                Assert.True(CollisionGap.GcdIdentityHolds(n, j), $"n={n}, j={j}");
                if (CollisionGap.GaloisKillsOddRung(n, j)) sawKill = true; else sawLive = true;
            }
        // the identity is a theorem, so no input can make it false and a `return true` would pass the
        // loop above on its own: what carries content is that BOTH sides take both values over the
        // range, so the equality is between two live predicates and not between two constants
        Assert.True(sawKill && sawLive);
        Assert.True(CollisionGap.GaloisKillsOddRung(9, 1));
        Assert.False(CollisionGap.GaloisKillsOddRung(9, 3));
        Assert.Equal(1, Cyclotomy.Gcd(1, 9));
        Assert.Equal(3, Cyclotomy.Gcd(3, 9));
    }

    // The rung j = 0 is never an automorphism, gcd(n, 2n) = n, and that is the rung which leaves the
    // first order standing. Without it the whole ladder would collapse and every gap would vanish.
    [Fact]
    public void TheZerothRungIsNeverAnAutomorphism()
    {
        for (int n = 9; n <= 40; n++)
            Assert.False(CollisionGap.GaloisKillsOddRung(n, 0), $"n={n}");
    }

    // The literal list, so the filter cannot assert itself: three at every odd firing modulus, and at
    // an even one the Galois route kills nothing at all, which is why Theorem E has to carry it.
    [Fact]
    public void ThreeIsTheFirstSurvivingRungAtEveryOddFiringModulus()
    {
        int[] oddFiring = { 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 81, 87, 93, 99 };
        foreach (int n in oddFiring)
        {
            Assert.True(LevelCollision.Fires(n), $"{n} is not a firing modulus");
            Assert.Equal(3, CollisionGap.FirstSurvivingOddRungByGcd(n));
        }
        // the shape reaches the same three, and it is the route that means something at an even comb,
        // where the gcd kills nothing and therefore refuses rather than naming a survivor of nothing
        Assert.Equal(3, CollisionGap.FirstSurvivingRungByShape());
        foreach (int n in new[] { 12, 18, 20, 24, 30 })
            Assert.Throws<ArgumentException>(() => CollisionGap.FirstSurvivingOddRungByGcd(n));
        // three is a property of the FIRING moduli (they all carry 3|n), not of the method: at an odd
        // comb the gcd route names the smallest prime factor of n, so a suite asking only firing moduli
        // cannot tell the search from a constant three
        Assert.Equal(5, CollisionGap.FirstSurvivingOddRungByGcd(25));
        Assert.Equal(7, CollisionGap.FirstSurvivingOddRungByGcd(49));
        Assert.Equal(11, CollisionGap.FirstSurvivingOddRungByGcd(121));
        Assert.Equal(5, CollisionGap.FirstSurvivingOddRungByGcd(55));
    }

    // The ROT3 route reaches three as well, by a different mechanism: forced whenever 3 does not
    // divide j. Both routes name the same rung and neither is the other.
    [Fact]
    public void TheRot3RouteAlsoFirstSurvivesAtThree()
    {
        for (int j = 1; j <= 30; j++)
            Assert.Equal(j % 3 != 0, CollisionGap.Rot3ForcesVanishing(j));
        // and the route's answer is READ OFF the forcing, not stated beside it: this dies with the
        // predicate above rather than agreeing with a literal three of its own
        Assert.Equal(3, CollisionGap.FirstSurvivingRungByShape());
        Assert.False(CollisionGap.Rot3ForcesVanishing(CollisionGap.FirstSurvivingRungByShape()));
        for (int j = 1; j < CollisionGap.FirstSurvivingRungByShape(); j++)
            Assert.True(CollisionGap.Rot3ForcesVanishing(j));
    }

    // ---- the closed forms, in multiplier form (Corollary C) ----

    [Fact]
    public void TheCosineIndicesAreConsecutiveAndTheWeightsSumToZero()
    {
        // consecutive on the order's own step: two apart on the odd orders, two on the even ones
        foreach (int order in new[] { 1, 2, 3, 5 })
        {
            var idx = CollisionGap.CosineIndices(order);
            for (int i = 1; i < idx.Count; i++)
                Assert.Equal(2, idx[i] - idx[i - 1]);
        }

        // d_m -> 0 as theta -> 0, where the mode has no endpoint amplitude: that is the identity the
        // weights carry. Below fifth order it is structurally true, the weights being +-1 by
        // construction and independent of n, so the content of this loop is the fifth-order column,
        // where three polynomials in n cancel identically. The magnitudes below fifth are pinned
        // separately, because this gate cannot see them.
        foreach (int order in new[] { 1, 2, 3, 5 })
            for (int n = 9; n <= 60; n++)
            {
                // the sum is formed HERE and not read off a boolean of the object's. The object used to
                // carry one; it was deleted, because the identity is a theorem and no input could make
                // such a predicate return false, so every gate on it was a gate that could not fail.
                BigInteger sum = BigInteger.Zero;
                foreach (int r in CollisionGap.CosineIndices(order))
                    sum += CollisionGap.WeightNumerator(n, order, r);
                Assert.Equal(BigInteger.Zero, sum);
            }
        // and the fifth-order column is the one with content: drop any one of its three weights and the
        // sum is a nonzero polynomial in n
        Assert.NotEqual(BigInteger.Zero,
            CollisionGap.WeightNumerator(30, 5, 2) + CollisionGap.WeightNumerator(30, 5, 4));
    }

    [Fact]
    public void TheFifthOrderSpansTheLadderFromZeroButDoesNotTouchTheConstantRung()
    {
        Assert.Equal(new[] { 0, 2, 4, 6 }, CollisionGap.CosineIndices(5));
        for (int n = 9; n <= 30; n++)
            Assert.Equal(BigInteger.Zero, CollisionGap.WeightNumerator(n, 5, 0));
    }

    // Why X_6 first enters at fifth order: the r = 6 weight is nonzero at every firing modulus, its
    // factors being positive for n >= 9.
    [Fact]
    public void TheSixthRungCarriesWeightAtEveryFiringModulus()
    {
        foreach (int n in Moduli)
            Assert.NotEqual(BigInteger.Zero, CollisionGap.WeightNumerator(n, 5, 6));
    }

    [Fact]
    public void TheClosedFormsAreTheProofsOwnAtTheSmallestFiringComb()
    {
        // d_1 = (2/n)(1 - cos 2t); d_2 = ((n-3)/n^2)(cos t - cos 3t);
        // d_3 = (2(n-2)(n-4)/(3n^3))(cos 2t - cos 4t); F = 4(n-6)(n-2)/(15 n^5).
        // three moduli, not one: at n = 9 alone (n-3) and (2n-12) both read 6, (n-2)(n-4) and 35*1 both
        // read 35, and (n-6)(n-2) and 21*1 both read 21, so a single comb pins none of the factorisations
        foreach (int n in new[] { 9, 12, 15 })
        {
            BigInteger b = n;
            Assert.Equal((new BigInteger(2), b), CollisionGap.Prefactor(n, 1));
            Assert.Equal((b - 3, b * b), CollisionGap.Prefactor(n, 2));
            Assert.Equal((2 * (b - 2) * (b - 4), 3 * b * b * b), CollisionGap.Prefactor(n, 3));
            Assert.Equal((4 * (b - 6) * (b - 2), 15 * BigInteger.Pow(b, 5)), CollisionGap.Prefactor(n, 5));
        }
        // and the arithmetic values, so the shapes above and the numbers cannot drift together
        Assert.Equal((new BigInteger(2), new BigInteger(9)), CollisionGap.Prefactor(9, 1));
        Assert.Equal((new BigInteger(6), new BigInteger(81)), CollisionGap.Prefactor(9, 2));
        Assert.Equal((new BigInteger(70), new BigInteger(2187)), CollisionGap.Prefactor(9, 3));
        Assert.Equal((new BigInteger(84), new BigInteger(885735)), CollisionGap.Prefactor(9, 5));
        Assert.Equal((new BigInteger(160), new BigInteger(5184)), CollisionGap.Prefactor(12, 3));   // 2*10*8 over 3*12^3

        // the denominator the weights are written over, which nothing else in the world reads
        Assert.Equal(BigInteger.One, CollisionGap.WeightDenominator(1));
        Assert.Equal(BigInteger.One, CollisionGap.WeightDenominator(2));
        Assert.Equal(BigInteger.One, CollisionGap.WeightDenominator(3));
        Assert.Equal(new BigInteger(4), CollisionGap.WeightDenominator(5));

        // the fifth-order weights as the proof writes them, at two moduli
        foreach (int n in new[] { 9, 30 })
        {
            BigInteger b = n;
            Assert.Equal(5 * (b - 1), CollisionGap.WeightNumerator(n, 5, 2));
            Assert.Equal(2 * (3 * b * b - 16 * b + 16), CollisionGap.WeightNumerator(n, 5, 4));
            Assert.Equal(-3 * (2 * b - 3) * (b - 3), CollisionGap.WeightNumerator(n, 5, 6));
        }

        // through third order the two weights are +1 and -1 and not merely antisymmetric, which is all
        // the sum-to-zero identity can see
        foreach (int order in new[] { 1, 2, 3 })
        {
            var idx = CollisionGap.CosineIndices(order);
            Assert.Equal(BigInteger.One, CollisionGap.WeightNumerator(21, order, idx[0]));
            Assert.Equal(BigInteger.MinusOne, CollisionGap.WeightNumerator(21, order, idx[1]));
        }

        Assert.Equal(new[] { 0, 2 }, CollisionGap.CosineIndices(1));
        Assert.Equal(new[] { 1, 3 }, CollisionGap.CosineIndices(2));
        Assert.Equal(new[] { 2, 4 }, CollisionGap.CosineIndices(3));
        Assert.Equal(BigInteger.One, CollisionGap.WeightNumerator(9, 1, 0));
        Assert.Equal(BigInteger.MinusOne, CollisionGap.WeightNumerator(9, 1, 2));
    }

    // d_4 is computed in the proof and deliberately not carried; the general shape past the fifth is
    // open. Both refusals are the fence, not a gap in the implementation.
    [Fact]
    public void TheLadderRefusesTheOrdersItDoesNotCarry()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.CosineIndices(4));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.CosineIndices(6));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.CosineIndices(7));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.Prefactor(9, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.WeightDenominator(4));
        // and a cosine index the order is not written on is REFUSED, not weighed zero: the rung j = 3
        // reads the cosine index 6, and a silent zero there would say the surviving rung vanishes
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.WeightNumerator(9, 5, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.WeightNumerator(9, 1, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.WeightNumerator(9, 3, 0));
        Assert.Equal(BigInteger.Zero, CollisionGap.WeightNumerator(9, 5, 0));   // carried, and zero
        // the object is not laxer than its Formulas face
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.OddOrderMultiplier(4, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.IsAutomorphism(4, 1));
    }

    // ---- the combs themselves, exactly, in GF(p) at a 2n-th root ----

    // X_2j = -M_{n+2j}, the identity that makes one ladder out of two, on every clean triple at every
    // firing modulus and every carried rung.
    [Fact]
    public void TheXLadderIsTheMLadderShiftedByTheModulus()
    {
        foreach (int n in Moduli)
        {
            var (p, zeta) = ModP.CyclotomicPrime(2 * n, 0);
            foreach (var t in LevelCollision.CleanTriples(n))
                for (int j = 0; j <= 4; j++)
                {
                    long x = CollisionGap.TwoX(n, j, t, p, zeta);
                    long m = CollisionGap.TwoM(n, CollisionGap.OddOrderMultiplier(n, j), t, p, zeta);
                    Assert.Equal(x, ModP.Mod(-m, p));
                }
        }
    }

    // ---- Theorem E: the ROT3 rung lemma ----

    [Fact]
    public void ParityUniformRot3TriplesHaveEveryUnthirdedRungForced()
    {
        int seen = 0;
        foreach (int n in new[] { 9, 12, 15, 18, 21, 24, 27, 30 })
        {
            var (p, zeta) = ModP.CyclotomicPrime(2 * n, 0);
            int here = 0;
            foreach (var t in LevelCollision.CleanTriples(n))
            {
                if (!CollisionGap.IsDoubledLabelRot3(n, t) || !CollisionGap.IsParityUniform(t)) continue;
                here++;
                for (int j = 1; j <= 9; j++)
                    if (CollisionGap.Rot3ForcesVanishing(j))
                        Assert.Equal(0L, CollisionGap.TwoX(n, j, t, p, zeta));
            }
            // per modulus, not as one global threshold: the four even combs alone supply 80, so a
            // single "more than eighty" would stay green with every odd comb gone silent
            Assert.True(here > 0, $"no triple met the hypothesis at n = {n}");
            seen += here;
        }
        Assert.Equal(100, seen);
    }

    // The break-input, and it is not hypothetical: at 3 | j the coset collapses to a specific cosine,
    // which at n = 24, (1,7,9) is cos(pi/4) at j = 3 and cos(pi/2) at j = 6. A lemma reading "zero iff
    // 3 does not divide j" would be false here.
    [Fact]
    public void TheCollapseAtThirdedRungsIsACosineAndNotAForcedZero()
    {
        const int n = 24;
        var t = (1, 7, 9);
        Assert.True(CollisionGap.IsDoubledLabelRot3(n, t));
        Assert.True(CollisionGap.IsParityUniform(t));
        var (p, zeta) = ModP.CyclotomicPrime(2 * n, 0);
        Assert.NotEqual(0L, CollisionGap.TwoX(n, 3, t, p, zeta));
        Assert.Equal(0L, CollisionGap.TwoX(n, 6, t, p, zeta));
    }

    // At an even n parity-uniformity comes with the shape and is not a second assumption: 6 | n makes
    // every coset parity-homogeneous. Eighty ROT3 triples at the four even moduli, all of them.
    [Fact]
    public void AtEvenCombsTheShapeCarriesTheParityForFree()
    {
        int rot3 = 0;
        foreach (int n in new[] { 12, 18, 24, 30 })
            foreach (var t in LevelCollision.CleanTriples(n))
                if (CollisionGap.IsDoubledLabelRot3(n, t))
                {
                    rot3++;
                    Assert.True(CollisionGap.IsParityUniform(t), $"n={n}, {t}");
                }
        Assert.Equal(80, rot3);
    }

    // At an odd n the hypothesis is real and load-bearing, and the proof names the witness.
    [Fact]
    public void AtOddCombsTheParityHypothesisIsLoadBearing()
    {
        const int n = 9;
        var t = (1, 2, 4);
        Assert.Contains(t, LevelCollision.CleanTriples(n));
        Assert.True(CollisionGap.IsDoubledLabelRot3(n, t));
        Assert.False(CollisionGap.IsParityUniform(t));
        var (p, zeta) = ModP.CyclotomicPrime(2 * n, 0);
        Assert.NotEqual(0L, CollisionGap.TwoX(n, 1, t, p, zeta));
    }

    // ---- the census the law explains ----

    [Fact]
    public void TheCensusReproducesTheCommittedTableColumnForColumn()
    {
        int pairs = 0, separate = 0, stand = 0;
        for (int i = 0; i < Moduli.Length; i++)
        {
            var c = CollisionGap.CensusOf(Moduli[i]);
            Assert.Equal(TablePairs[i], c.Pairs);
            Assert.Equal(TableSeparate[i], c.Separating);
            Assert.Equal(TableStand[i], c.Standing);
            Assert.Equal(c.Pairs, c.Separating + c.Standing);
            pairs += c.Pairs; separate += c.Separating; stand += c.Standing;
        }
        Assert.Equal(2558, pairs);
        Assert.Equal(2335, separate);
        Assert.Equal(223, stand);
    }

    // The pair count is LevelCollision's own, so the two objects are reading one census and not two.
    [Fact]
    public void TheGapCensusPairsAreTheLevelCensusPairs()
    {
        for (int i = 0; i < Moduli.Length; i++)
        {
            Assert.Equal(TablePairs[i], LevelCollision.CensusOf(Moduli[i]).CollidingPairs);
            Assert.Equal(TablePairs[i], (int)LevelCollision.InventoryOf(Moduli[i]).Total);
        }
    }

    // Theorem D's consequence, read on the census: at every odd modulus and for EVERY pair, standing
    // or separating, the second and fourth rungs are dead, and hence c_3 = 0.
    [Fact]
    public void EveryOddCombPairHasBothOddRungsDead()
    {
        int odd = 0, separating = 0;
        foreach (int n in new[] { 9, 15, 21, 27 })
        {
            var (p, zeta) = ModP.CyclotomicPrime(2 * n, 0);
            foreach (var (a, b) in LevelCollision.CollidingPairs(n))
            {
                odd++;
                Assert.Equal(0L, CollisionGap.DeltaTwoX(n, 1, a, b, p, zeta));
                Assert.Equal(0L, CollisionGap.DeltaTwoX(n, 2, a, b, p, zeta));
                if (CollisionGap.FirstOrderSeparates(n, a, b)) separating++;
            }
        }
        Assert.Equal(627, odd);
        Assert.Equal(427, separating);
    }

    // At an odd comb the first order is an integer count on labels: c_1 = (4/n)(o_tau - o_sigma), so
    // the speed takes only four values. No cosine survives to carry anything else.
    [Fact]
    public void AtOddCombsTheFirstOrderIsALabelCount()
    {
        foreach (int n in new[] { 9, 15, 21, 27 })
            foreach (var (a, b) in LevelCollision.CollidingPairs(n))
            {
                int d = CollisionGap.OddLabelCount(a) - CollisionGap.OddLabelCount(b);
                Assert.InRange(Math.Abs(d), 0, 3);
                Assert.Equal(d != 0, CollisionGap.FirstOrderSeparates(n, a, b));
            }
    }

    // The twenty-three even-n standing pairs split eleven to twelve, and the two halves are carried by
    // DIFFERENT arguments: the eleven term by term as theta-mirrors, the twelve by the ROT3 shape.
    [Fact]
    public void TheEvenStandingPairsSplitElevenToTwelve()
    {
        int standing = 0, mirror = 0, nonMirror = 0;
        var nonMirrorByComb = new Dictionary<int, int>();
        foreach (int n in new[] { 12, 18, 20, 24, 30 })
        {
            var c = CollisionGap.CensusOf(n);
            standing += c.Standing;
            mirror += c.StandingMirror;
            nonMirror += c.StandingNonMirror;
            nonMirrorByComb[n] = c.StandingNonMirror;
        }
        Assert.Equal(23, standing);
        Assert.Equal(11, mirror);
        Assert.Equal(12, nonMirror);
        Assert.Equal(4, nonMirrorByComb[24]);
        Assert.Equal(8, nonMirrorByComb[30]);
    }

    // Read, not derived, and the proof says so: all twelve have both triples ROT3 and share a parity
    // class. Nothing here says a standing pair must.
    [Fact]
    public void TheTwelveAreRot3AndShareAParityClass()
    {
        int seen = 0;
        foreach (int n in new[] { 12, 18, 20, 24, 30 })
            foreach (var (a, b) in LevelCollision.CollidingPairs(n))
            {
                if (CollisionGap.FirstOrderSeparates(n, a, b)) continue;
                if (CollisionGap.IsThetaMirror(n, a, b)) continue;
                seen++;
                Assert.True(CollisionGap.IsDoubledLabelRot3(n, a) && CollisionGap.IsDoubledLabelRot3(n, b), $"n={n}");
                Assert.Equal(CollisionGap.OddLabelCount(a), CollisionGap.OddLabelCount(b));
            }
        Assert.Equal(12, seen);
    }

    // The eleven theta-mirrors need NO shape, which matters because one of them is not ROT3 at all.
    [Fact]
    public void TheThetaMirrorsAreCarriedTermByTermAndOneIsNotRot3()
    {
        var a = (6, 18, 20);
        var b = (10, 12, 24);
        Assert.True(CollisionGap.IsThetaMirror(30, a, b));
        Assert.False(CollisionGap.IsDoubledLabelRot3(30, a) && CollisionGap.IsDoubledLabelRot3(30, b));

        var (p, zeta) = ModP.CyclotomicPrime(60, 0);
        for (int j = 0; j <= 3; j++)
            Assert.Equal(0L, CollisionGap.DeltaTwoX(30, j, a, b, p, zeta));
    }

    // The carrying is the EVEN comb's half. At an odd n the same reflection flips every eta, so nothing
    // is carried, and the census contains the witness: the single pair at n = 9 IS a theta-mirror and it
    // separates. So no odd-n theta-mirror can ever reach the mirror column, and that column being all
    // even is a theorem rather than an accident of these nine moduli.
    [Fact]
    public void AnOddCombThetaMirrorAlwaysSeparates()
    {
        Assert.True(CollisionGap.IsThetaMirror(9, (1, 5, 7), (2, 4, 8)));
        Assert.True(CollisionGap.FirstOrderSeparates(9, (1, 5, 7), (2, 4, 8)));
        foreach (int n in new[] { 9, 15, 21, 27 })
            foreach (var (a, b) in LevelCollision.CollidingPairs(n))
                if (CollisionGap.IsThetaMirror(n, a, b))
                {
                    // o_sigma = 3 - o_tau, so DeltaX_0 = 2(2 o_tau - 3) is never zero
                    Assert.Equal(3 - CollisionGap.OddLabelCount(a), CollisionGap.OddLabelCount(b));
                    Assert.True(CollisionGap.FirstOrderSeparates(n, a, b));
                }
    }

    // The shape gives c_3, the parity match gives c_1, and neither gives the other: fifty-eight pairs
    // satisfy the ROT3 hypothesis in full and still separate at first order, on both sides of the parity.
    [Fact]
    public void TheShapeAloneDoesNotStopTheFirstOrder()
    {
        int both = 0, odd = 0, even = 0;
        foreach (int n in Moduli)
            foreach (var (a, b) in LevelCollision.CollidingPairs(n))
            {
                bool shaped = CollisionGap.IsDoubledLabelRot3(n, a) && CollisionGap.IsParityUniform(a)
                           && CollisionGap.IsDoubledLabelRot3(n, b) && CollisionGap.IsParityUniform(b);
                if (!shaped || !CollisionGap.FirstOrderSeparates(n, a, b)) continue;
                both++;
                if (n % 2 == 1) odd++; else even++;
                Assert.Equal(3, Math.Abs(CollisionGap.OddLabelCount(a) - CollisionGap.OddLabelCount(b)));
            }
        Assert.Equal(58, both);
        Assert.Equal(30, odd);
        Assert.Equal(28, even);
    }

    // The smallest modulus that carries one is n = 9, and it is the pair Confirmation 24 flew, on
    // another axis. Nothing here is confirmed by that flight.
    [Fact]
    public void TheSmallestShapedSeparatingPairIsTheNineOne()
    {
        var a = (1, 5, 7);
        var b = (2, 4, 8);
        Assert.True(CollisionGap.IsDoubledLabelRot3(9, a) && CollisionGap.IsDoubledLabelRot3(9, b));
        Assert.Equal(6, CollisionGap.DeltaXZero(a, b));
        Assert.True(CollisionGap.FirstOrderSeparates(9, a, b));
    }

    // DeltaX_0 is an integer count, and the ladder's zeroth rung is the same number read in the field.
    // Tying them is what keeps the count honest: without this the factor in X_0 = 2o - 3 is held by one
    // literal at one modulus, and a mutation that halves it survives the whole census.
    [Fact]
    public void TheZerothRungIsTheLabelCountReadInTheField()
    {
        foreach (int n in Moduli)
        {
            var (p, zeta) = ModP.CyclotomicPrime(2 * n, 0);
            foreach (var (a, b) in LevelCollision.CollidingPairs(n))
                Assert.Equal(
                    ModP.Mod(2L * CollisionGap.DeltaXZero(a, b), p),
                    CollisionGap.DeltaTwoX(n, 0, a, b, p, zeta));
        }
    }

    // The two primes are the census's insurance and not its discriminator, and this gate says so rather
    // than dressing it up. Both readings run through the same TwoCos/TwoX/TwoM and are exact evaluations
    // of ONE element of Z[zeta_2n], so under any mutation of the physics they move together and this
    // stays green; what it can catch is a corrupted SECOND setting. So the property that actually
    // matters is asserted directly here as well: the second zeta has exact multiplicative order 2n,
    // walked power by power rather than taken from RootOfOrder's own test. The agreement over all 2558
    // pairs is kept beside it as a MEASUREMENT of this census, which is what the README cites.
    [Fact]
    public void TheSecondSettingIsARealSecondFieldAndTheTwoNeverDisagreeHere()
    {
        foreach (int n in Moduli)
        {
            var settings = CollisionGap.Settings(n);
            Assert.NotEqual(settings[0].P, settings[1].P);
            foreach (var (p, zeta) in settings)
            {
                Assert.Equal(1L, p % (2 * n));
                System.Numerics.BigInteger acc = 1;
                int first = 0;
                for (int e = 1; e <= 2 * n; e++)
                {
                    acc = acc * zeta % p;
                    if (acc == 1) { first = e; break; }
                }
                Assert.Equal(2 * n, first);
            }
            foreach (var (a, b) in LevelCollision.CollidingPairs(n))
            {
                var first = settings.Select(s =>
                    ModP.Mod(2L * CollisionGap.DeltaXZero(a, b) - CollisionGap.DeltaTwoX(n, 1, a, b, s.P, s.Zeta), s.P) != 0).ToArray();
                Assert.Equal(first[0], first[1]);
                var third = settings.Select(s => CollisionGap.DeltaTwoM(n, 3, a, b, s.P, s.Zeta) != 0).ToArray();
                Assert.Equal(third[0], third[1]);
            }
        }
    }

    // ---- the second order ----

    [Fact]
    public void TheSecondOrderDiesAtExactlySixtyPairsOfTheCensus()
    {
        int zero = 0;
        var byComb = new Dictionary<int, int>();
        foreach (int n in Moduli)
        {
            var c = CollisionGap.CensusOf(n);
            zero += c.SecondOrderZero;
            byComb[n] = c.SecondOrderZero;
        }
        Assert.Equal(60, zero);
        Assert.Equal(20, byComb[20]);
        Assert.Equal(40, byComb[30]);
        Assert.Equal(0, byComb[9]);
        Assert.Equal(0, byComb[21]);
    }

    // The second order splits by 3|n and NOT by parity, which is the trap this gate exists to hold: its
    // rung is m = 3, so the gcd kills c_2 exactly where 3 does not divide n. Among the nine firing
    // moduli that is n = 20 alone, and there it kills EVERY pair. At every 3|n modulus the gcd fails and
    // only Corollary G's local criterion can speak, which this object does not run: so the 40 at n = 30
    // are not explained here, and nothing in this file should be read as explaining them.
    [Fact]
    public void TheSecondOrderSplitsByThreeDividingNAndNotByParity()
    {
        foreach (int n in Moduli)
        {
            bool gcdKills = Cyclotomy.Gcd(3, 2 * n) == 1;
            Assert.Equal(n % 3 != 0, gcdKills);
            var c = CollisionGap.CensusOf(n);
            if (gcdKills) Assert.Equal(c.Pairs, c.SecondOrderZero);   // every pair, at n = 20
        }
        // n = 20 is the only census modulus with 3 not dividing n, and it is EVEN: so the parity split
        // that governs the first order gets this one exactly backwards
        Assert.Equal(new[] { 20 }, Moduli.Where(n => n % 3 != 0).ToArray());
        Assert.Equal(0, 20 % 2);
        // and it carries no ROT3 triple at all, so the shape route has nothing to speak about there
        Assert.DoesNotContain(LevelCollision.CleanTriples(20), t => CollisionGap.IsDoubledLabelRot3(20, t));

        // "neither route, and none is needed" is asserted in four places and was gated in none. The
        // second half is the load-bearing one: at n = 20 there is no odd-order vanishing to explain,
        // every pair separating at first order AND no pair having a dead second or fourth rung.
        var (p, zeta) = ModP.CyclotomicPrime(40, 0);
        foreach (var (a, b) in LevelCollision.CollidingPairs(20))
        {
            Assert.True(CollisionGap.FirstOrderSeparates(20, a, b));
            Assert.NotEqual(0L, CollisionGap.DeltaTwoX(20, 1, a, b, p, zeta));
        }
        Assert.Equal(0, CollisionGap.CensusOf(20).Standing);
    }

    // The lower bound past the census is not a new formula: it IS LevelCollision's families C and L,
    // read off the inventory. The bound counts only what family membership forces.
    [Fact]
    public void TheSecondOrderBoundIsTheInventorysOwnTwoFamilies()
    {
        foreach (int n in new[] { 9, 20, 21, 30, 70, 210 })
            Assert.Equal(
                LevelCollision.FamilyCount(LevelCollision.CollisionFamily.C, n)
                    + LevelCollision.FamilyCount(LevelCollision.CollisionFamily.L, n),
                CollisionGap.SecondOrderZeroLowerBound(n));

        Assert.Equal(0L, CollisionGap.SecondOrderZeroLowerBound(9));
        Assert.Equal(0L, CollisionGap.SecondOrderZeroLowerBound(21));
        Assert.Equal(20L, CollisionGap.SecondOrderZeroLowerBound(20));
        Assert.Equal(40L, CollisionGap.SecondOrderZeroLowerBound(30));
        Assert.Equal(140L, CollisionGap.SecondOrderZeroLowerBound(70));
        Assert.Equal(420L, CollisionGap.SecondOrderZeroLowerBound(210));

        // and the half of the bound this census cannot see: family L needs 70 | n, so inside n <= 30 it
        // is zero everywhere and the equality at n = 20 and n = 30 is family C ALONE. Nothing here
        // exercises L, which is what the arc records about it too.
        foreach (int n in Moduli)
            Assert.Equal(0L, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.L, n));
        Assert.Equal(20L, LevelCollision.FamilyCount(LevelCollision.CollisionFamily.L, 70));
    }

    // A bound is a bound: at the two moduli inside the census where it fires it is MET, and that is a
    // reading of this census, not a theorem about every modulus.
    [Fact]
    public void TheBoundIsMetInsideTheCensusAndIsNotClaimedBeyondIt()
    {
        Assert.Equal(CollisionGap.SecondOrderZeroLowerBound(20), CollisionGap.CensusOf(20).SecondOrderZero);
        Assert.Equal(CollisionGap.SecondOrderZeroLowerBound(30), CollisionGap.CensusOf(30).SecondOrderZero);
    }

    // ---- the object ----

    [Fact]
    public void TheInstanceCensusIsTheCombItWasBuiltOn()
    {
        var g = Gap(21);
        Assert.Equal(21, g.Census().N);
        Assert.Equal(CollisionGap.CensusOf(21), g.Census());
        Assert.NotEqual(CollisionGap.CensusOf(27), g.Census());
    }

    [Fact]
    public void TheObjectOwnsTheSeriesAndInheritsTheRoad()
    {
        var g = Gap(9);
        Assert.Equal(new[] { "series", "ladder", "gap" }, g.Own);
        Assert.Equal(9, g.Ncomb);
        Assert.IsType<Crack>(g.Parent);
        // the road arrives as inheritance, and the combs arrive through it from the Cyclotomy
        Assert.Contains("road", g.Inherited);
        Assert.Contains(g.Inherited, s => s.Contains("turn fractions"));
        // the frame is not in the chain: Cyclotomy is a second parentless root
        Assert.DoesNotContain("x", g.Inherited);
        // and the comb is not this object's to own
        Assert.DoesNotContain("combs", g.Own);
    }

    [Fact]
    public void TheObjectRefusesWhatItCannotMean()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Gap(4));
        // and the contract really is n >= 5 and not n >= 9: below the smallest FIRING comb the object
        // still answers, it just answers about a comb with no collision on it
        foreach (int n in new[] { 5, 6, 7, 8 })
        {
            Assert.Equal(n, Gap(n).Ncomb);
            Assert.Equal(0, CollisionGap.CensusOf(n).Pairs);
            Assert.Equal(n + 2, CollisionGap.OddOrderMultiplier(n, 1));
            Assert.False(LevelCollision.Fires(n));
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.CensusOf(4));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.FirstSurvivingOddRungByGcd(4));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.Settings(4));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.IsThetaMirror(4, (1, 2, 3), (1, 2, 3)));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGap.IsDoubledLabelRot3(4, (1, 2, 3)));
        // the odd-n reduction has no content at an even comb and says so rather than returning true
        Assert.Throws<ArgumentException>(() => CollisionGap.GcdIdentityHolds(12, 1));
    }
}
