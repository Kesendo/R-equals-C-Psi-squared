using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>From-below guard for <see cref="CollisionGapOddOrdersClaim"/> (F161). Two routes meet everywhere they
/// can: the closed coefficients against an EXACT power-series solution of F160's polynomial over ℚ (possible at the
/// comb points where cos θ_k is rational, which is every k = n/3, n/2, 2n/3), the multiplier identity against the
/// committed level vector, the ROT3 collapse against the ring sum, and every census count against the triples. No
/// tolerance appears below except one: the cyclotomic route compares integer vectors and the series route compares
/// rationals, and where F160's ChainEndVelocity is met the comparison crosses into floating point, where no exact
/// route exists and the bound is an ulp-scale error model rather than a round number.</summary>
public class CollisionGapOddOrdersClaimTests
{
    private static readonly int[] OddFiringCombs = { 9, 15, 21, 27 };
    private static readonly int[] FiringCombsToThirty = { 9, 12, 15, 18, 20, 21, 24, 27, 30 };

    // ---------------------------------------------------------------- Theorem A, against an exact series over Q

    /// <summary>The comb points with a RATIONAL cosine, (n, k, cos θ_k): k = n/3 gives 1/2, k = n/2 gives 0,
    /// k = 2n/3 gives −1/2. Both signs of η and all three cosines appear, and the c = 0 rows are the only ones where
    /// d_3 is not identically zero (4c² − 1 = 0 at c = ±1/2), so the set exercises every coefficient.</summary>
    public static TheoryData<int, int, int, int> RationalCosineCombPoints()
    {
        var data = new TheoryData<int, int, int, int>();
        foreach (var (n, k, num, den) in new[]
                 {
                     (6, 2, 1, 2), (6, 3, 0, 1), (6, 4, -1, 2),
                     (9, 3, 1, 2), (9, 6, -1, 2),
                     (12, 4, 1, 2), (12, 6, 0, 1), (12, 8, -1, 2),
                     (15, 5, 1, 2), (15, 10, -1, 2),
                     (18, 6, 1, 2), (18, 9, 0, 1),
                     (20, 10, 0, 1),
                     (21, 7, 1, 2),
                     (30, 10, 1, 2), (30, 15, 0, 1), (30, 20, -1, 2),
                 })
            data.Add(n, k, num, den);
        return data;
    }

    [Theory]
    [MemberData(nameof(RationalCosineCombPoints))]
    public void TheoremA_TheClosedCoefficients_AreTheExactSeriesOfF160sPolynomial(int nComb, int k, int cosNum, int cosDen)
    {
        // The claim's own exact evaluation against a power-series solution of F160's polynomial. Both live in ℚ, so
        // this is an equality and not a threshold; a nonzero difference would be a finding about the closed form.
        var cos = new BigRational(cosNum, cosDen);
        var d = ExactRootSeries(nComb, 2 * cos);

        foreach (int order in CollisionGapOddOrdersClaim.CarriedExpansionOrders)
            Assert.Equal(CollisionGapOddOrdersClaim.LevelCoefficientExact(nComb, k, order, cos), d[order]);
    }

    [Fact]
    public void CorollaryC_EveryRungCarriesTheParityOfItsOrder_AndEveryOrdersWeightsSumToZero()
    {
        // r = m + 1 mod 2 is §(d)'s theorem: the even orders read odd multipliers and the odd orders even ones.
        for (int n = 9; n <= 60; n++)
            foreach (int m in CollisionGapOddOrdersClaim.CarriedExpansionOrders)
            {
                var rungs = CollisionGapOddOrdersClaim.MultiplierFormRungs(m);
                Assert.All(rungs, r => Assert.Equal((m + 1) % 2, r % 2));

                var total = BigRational.Zero;
                foreach (int r in rungs) total += CollisionGapOddOrdersClaim.MultiplierFormWeight(n, m, r);
                Assert.True(total.IsZero, $"n = {n}, d_{m}: the rung weights sum to {total}, not 0");

                Assert.True(CollisionGapOddOrdersClaim.MultiplierFormPrefactor(n, m).Sign > 0);
                Assert.True(CollisionGapOddOrdersClaim.MultiplierFormWeight(n, m, rungs[^1]).Sign < 0);
            }

        Assert.True(CollisionGapOddOrdersClaim.MultiplierFormWeight(12, 5, 0).IsZero);   // the absent constant rung
        Assert.True(CollisionGapOddOrdersClaim.MultiplierFormWeight(12, 1, 4).IsZero);   // a rung the order misses
    }

    [Fact]
    public void FirstOrder_IsF160sChainEndVelocity_ByTheIndependentSeriesRoute()
    {
        // F160's Theorem G and this file's d_1 are the same number by two derivations; the series is the third route
        // and is what ties them, so the comparison is not one formula against itself.
        foreach (var (nComb, k, num, den) in new[] { (6, 3, 0, 1), (12, 6, 0, 1), (9, 3, 1, 2), (30, 20, -1, 2) })
        {
            var d1 = ExactRootSeries(nComb, 2 * new BigRational(num, den))[1];
            double velocity = CrackedRingExactCurveClaim.ChainEndVelocity(nComb - 1, k);
            // F160's member is a floating-point evaluation of sin^2 and a division, so there is no exact route on
            // that side. At THESE theta the relative error is a few ulp of the value, and 4 * eps * |value| is what
            // is gated: not a round number that happens to pass. The bound is stated for these rows only, since
            // sin's relative error carries the argument amplification |theta * cot theta|, which grows with n as
            // k -> n - 1; the four rows here sit at theta = pi/2, pi/3 and 2pi/3, where that factor is below 1.3.
            Assert.True(Math.Abs(velocity - ToDouble(d1)) <= 4 * 2.220446049250313e-16 * Math.Abs(ToDouble(d1)),
                $"n = {nComb}, k = {k}: series {d1}, F160 velocity {velocity}");
        }
    }

    [Fact]
    public void FourthOrder_IsRefused_RatherThanReturned()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGapOddOrdersClaim.LevelCoefficient(12, 5, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGapOddOrdersClaim.LevelCoefficient(12, 5, 6));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGapOddOrdersClaim.MultiplierFormRungs(4));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGapOddOrdersClaim.MultiplierFormPrefactor(12, 4));
    }

    [Theory]
    [MemberData(nameof(RationalCosineCombPoints))]
    public void TheDoublePath_AssemblesTheSameNumberAsTheExactOne(int nComb, int k, int cosNum, int cosDen)
    {
        // LevelCoefficient and LevelCoefficientExact share the rungs, the weights and the prefactor but each writes
        // the eta line and the cos r*theta step for itself, so this pins that pair of copies. It is arithmetic and
        // not physics: the closed form itself is under test one Theory up, against the series.
        var cos = new BigRational(cosNum, cosDen);
        foreach (int order in CollisionGapOddOrdersClaim.CarriedExpansionOrders)
        {
            double exact = ToDouble(CollisionGapOddOrdersClaim.LevelCoefficientExact(nComb, k, order, cos));
            double read = CollisionGapOddOrdersClaim.LevelCoefficient(nComb, k, order);
            // Every order's weights sum to zero, so the double path is a CANCELLING sum and its absolute error is set
            // by the size of the terms, not of the result: at c = +-1/2 the exact d_3 is 0 while the float route
            // returns ~1e-17. The model is therefore eps times the sum of |weight| times the prefactor, which is the
            // largest quantity the evaluation actually handles.
            double terms = ToDouble(CollisionGapOddOrdersClaim.MultiplierFormPrefactor(nComb, order))
                         * CollisionGapOddOrdersClaim.MultiplierFormRungs(order)
                               .Sum(r => Math.Abs(ToDouble(CollisionGapOddOrdersClaim.MultiplierFormWeight(nComb, order, r))));
            Assert.True(Math.Abs(read - exact) <= 8 * 2.220446049250313e-16 * terms,
                $"d_{order} at n = {nComb}, k = {k}: exact {exact}, double path {read}, term scale {terms}");
        }
    }


    // ---------------------------------------------------------------- Corollary C, the multiplier identity

    [Fact]
    public void CorollaryC_XIsTheCombShiftedByN_AndMAtMultiplierOneIsTheCommittedLevelVector()
    {
        foreach (int n in FiringCombsToThirty)
            foreach (var t in LevelCollisionCensus.CleanTriples(n))
            {
                Assert.Equal(
                    CyclotomicRing.Key(LevelCollisionCensus.LevelVector(n, t.K1, t.K2, t.K3)),
                    CyclotomicRing.Key(CollisionGapOddOrdersClaim.TwoTimesM(n, t, 1)));

                for (int j = 0; j <= 4; j++)
                {
                    var x = CollisionGapOddOrdersClaim.TwoTimesX(n, t, j);
                    var shifted = CyclotomicRing.Negate(
                        CollisionGapOddOrdersClaim.TwoTimesM(n, t, CollisionGapOddOrdersClaim.OddOrderMultiplier(n, j)));
                    Assert.True(CyclotomicRing.AreEqual(x, shifted),
                        $"n = {n}, tau = {t}, j = {j}: X_2j != -M_(n+2j)");
                }
            }
    }

    [Fact]
    public void XAtRungZero_IsTwiceTheOddLabelSurplus()
    {
        // Sum over tau of eta = 2o - 3, so 2*X_0 = 2(2o - 3) as a rational integer in the ring.
        foreach (int n in FiringCombsToThirty)
            foreach (var t in LevelCollisionCensus.CleanTriples(n))
            {
                var expected = CyclotomicRing.Zero(2 * n);
                CyclotomicRing.AddRootPower(expected, 2 * n, 0, 2L * (2 * CollisionGapOddOrdersClaim.OddLabelCount(t) - 3));
                Assert.True(CyclotomicRing.AreEqual(CollisionGapOddOrdersClaim.TwoTimesX(n, t, 0), expected));
            }
    }

    // ---------------------------------------------------------------- Theorem D, the Galois kill

    [Fact]
    public void TheoremD_TheOddCombGcdIdentity_HoldsAtEveryOddNAndEveryRung()
    {
        for (int n = 9; n <= 59; n += 2)
            for (int j = 0; j <= 12; j++)
            {
                var (full, reduced) = CollisionGapOddOrdersClaim.OddCombGcdIdentity(n, j);
                Assert.Equal(reduced, full);
            }
    }

    [Fact]
    public void TheoremD_RungZeroIsNeverAnAutomorphism_AndTheEvenCombLadderIsNeverReached()
    {
        for (int n = 4; n <= 40; n++)
        {
            Assert.False(CollisionGapOddOrdersClaim.OddOrderRungIsGaloisKill(n, 0));
            if (n % 2 == 0)
                for (int j = 0; j <= 10; j++)
                    Assert.False(CollisionGapOddOrdersClaim.OddOrderRungIsGaloisKill(n, j));
        }
    }

    [Fact]
    public void TheoremF_TheFirstSurvivingRungIsThreeAtEveryOddFiringModulus()
    {
        // The moduli are a LITERAL list, not the firing predicate: at odd n the predicate's second clause (10|n) is
        // unreachable, so filtering on it and then asserting 3|n would be the filter asserting itself. That exact
        // gate was removed from the Python gate the same day (docs/CAUGHT_ERRORS.md, 2026-09-02, L11).
        int[] oddFiring = { 9, 15, 21, 27, 33, 39, 45, 51, 57, 63, 69, 75, 81, 87, 93, 99 };
        Assert.Equal(oddFiring, Enumerable.Range(0, 46).Select(i => 9 + 2 * i)   // 9, 11, .., 99
                                          .Where(LevelCollisionCensus.Fires).ToArray());
        foreach (int n in oddFiring)
            Assert.Equal(3, CollisionGapOddOrdersClaim.FirstSurvivingOddOrderRung(n));
    }

    [Fact]
    public void TheoremD_AtOddN_EveryCollisionPairHasAVanishingThirdOrder()
    {
        int pairs = 0;
        foreach (int n in OddFiringCombs)
            foreach (var (tau, sigma) in CollisionGapOddOrdersClaim.CollisionPairs(n))
            {
                pairs++;
                Assert.True(CyclotomicRing.IsZero(CollisionGapOddOrdersClaim.DeltaTwoTimesX(n, tau, sigma, 1)));
                Assert.True(CyclotomicRing.IsZero(CollisionGapOddOrdersClaim.DeltaTwoTimesX(n, tau, sigma, 2)));

                // The first order is then the odd-label count. The JUDGE is the ring: with DeltaX_2 = 0 just
                // established, c_1 = (2/n)(DeltaX_0 - DeltaX_2) = (2/n)*DeltaX_0, and DeltaTwoTimesX at rung 0 is
                // the rational integer 2*DeltaX_0, so c_1 = that integer over n. Comparing instead against
                // 4(o_tau - o_sigma)/n would be the method's own body retyped.
                var c1 = CollisionGapOddOrdersClaim.FirstOrderCoefficient(n, tau, sigma);
                Assert.Equal(new BigRational(RationalIntegerValue(n, tau, sigma), n), c1);
                Assert.Equal(c1.IsZero, CollisionGapOddOrdersClaim.StandsAtFirstOrder(n, tau, sigma));
            }
        Assert.Equal(627, pairs);
    }

    [Fact]
    public void TheoremD_AtEveryCoprimeRung_TheKillFires()
    {
        int checks = 0;
        foreach (int n in OddFiringCombs)
            foreach (var (tau, sigma) in CollisionGapOddOrdersClaim.CollisionPairs(n).Take(30))
                for (int j = 1; j <= 8; j++)
                {
                    if (CollisionGapOddOrdersClaim.OddCombGcdIdentity(n, j).Reduced != 1) continue;
                    checks++;
                    Assert.True(CyclotomicRing.IsZero(CollisionGapOddOrdersClaim.DeltaTwoTimesX(n, tau, sigma, j)),
                        $"n = {n}, j = {j}: a coprime rung left DeltaX nonzero");
                }
        // Pinned, because a Take + a continue is exactly the shape that goes vacuous under an unrelated edit.
        Assert.Equal(486, checks);
    }

    [Fact]
    public void TheoremD_NegativeControl_DeltaX2VanishesOnTheCollisionsAndOnNothingElse()
    {
        // The population is EVERY clean-triple pair, not the collisions, so a law that fired everywhere would not
        // pass. What it does NOT do is distinguish a false law: at odd n this is the converse of Theorem D, proved
        // in the same section, so it is a regression check on the implementation. The committed gate says exactly
        // that of its own L8 and the wording is carried here rather than upgraded.
        foreach (int n in new[] { 9, 15 })
        {
            var triples = LevelCollisionCensus.CleanTriples(n).ToList();
            int checkedPairs = 0, collisions = 0;
            for (int i = 0; i < triples.Count; i++)
                for (int j = i + 1; j < triples.Count; j++)
                {
                    checkedPairs++;
                    bool collides = CollisionGapOddOrdersClaim.Collide(n, triples[i], triples[j]);
                    bool x2Zero = CyclotomicRing.IsZero(
                        CollisionGapOddOrdersClaim.DeltaTwoTimesX(n, triples[i], triples[j], 1));
                    Assert.Equal(collides, x2Zero);
                    if (collides) collisions++;
                }
            Assert.True(checkedPairs > 0 && collisions > 0);
        }
    }

    // ---------------------------------------------------------------- Theorem E, the ROT3 rung lemma

    [Fact]
    public void TheoremE_TheForcedDirection_HoldsOnEveryParityUniformRot3TripleAtEveryRungNotDivisibleByThree()
    {
        int triplesSeen = 0;
        foreach (int n in FiringCombsToThirty)
        {
            if (n % 3 != 0) continue;
            foreach (var t in LevelCollisionCensus.CleanTriples(n))
            {
                if (CollisionGapOddOrdersClaim.DoubledLabelRot3Cosets(n, t) is null) continue;
                if (!CollisionGapOddOrdersClaim.IsParityUniform(t)) continue;
                triplesSeen++;
                for (int j = 1; j <= 9; j++)
                {
                    // Rot3ForcesVanishing is the conjunction of the three filters above and so cannot fail here;
                    // the content is that the RING SUM really vanishes where the lemma says it must.
                    if (j % 3 == 0) continue;
                    Assert.True(CyclotomicRing.IsZero(CollisionGapOddOrdersClaim.TwoTimesX(n, t, j)),
                        $"n = {n}, tau = {t}, j = {j}: the forced direction failed");
                }
            }
        }
        // 2 + 8 + 4 + 16 + 6 + 24 + 8 + 32 at n = 9..30, so a narrowing of either predicate would show here.
        Assert.Equal(100, triplesSeen);
    }

    [Fact]
    public void TheoremE_AtEvenN_ParityUniformityComesFreeWithTheShape()
    {
        // 3|n and 2|n give 6|n, so n/3 is even and every coset is parity-homogeneous. Read as a count too: the
        // eighty ROT3 triples at n = 12, 18, 24, 30.
        int rot3 = 0;
        foreach (int n in new[] { 12, 18, 24, 30 })
            foreach (var t in LevelCollisionCensus.CleanTriples(n))
                if (CollisionGapOddOrdersClaim.DoubledLabelRot3Cosets(n, t) is not null)
                {
                    rot3++;
                    Assert.True(CollisionGapOddOrdersClaim.IsParityUniform(t), $"n = {n}, tau = {t}");
                }
        Assert.Equal(80, rot3);
    }

    [Fact]
    public void TheoremE_AtOddN_ParityUniformityIsALoadBearingHypothesis()
    {
        var mixed = (1, 2, 4);
        Assert.NotNull(CollisionGapOddOrdersClaim.DoubledLabelRot3Cosets(9, mixed));
        Assert.False(CollisionGapOddOrdersClaim.IsParityUniform(mixed));
        Assert.False(CyclotomicRing.IsZero(CollisionGapOddOrdersClaim.TwoTimesX(9, mixed, 1)));
        Assert.False(CollisionGapOddOrdersClaim.Rot3ForcesVanishing(9, mixed, 1));
    }

    [Fact]
    public void TheoremE_TheConverseIsFalse_AndTheBreakInputExhibitsBothOutcomes()
    {
        var tau = (1, 7, 9);
        Assert.NotNull(CollisionGapOddOrdersClaim.DoubledLabelRot3Cosets(24, tau));
        Assert.True(CollisionGapOddOrdersClaim.IsParityUniform(tau));

        // j = 1: the lemma forces the vanishing, which is the TRUE branch and is pinned here on a triple whose
        // shape is established two lines up rather than by restating the method's own conjunction.
        Assert.True(CollisionGapOddOrdersClaim.Rot3ForcesVanishing(24, tau, 1));
        // j = 3: the coset collapses and the cosine does NOT vanish (cos(pi/4)).
        Assert.False(CollisionGapOddOrdersClaim.Rot3ForcesVanishing(24, tau, 3));
        Assert.False(CyclotomicRing.IsZero(CollisionGapOddOrdersClaim.TwoTimesX(24, tau, 3)));
        // j = 6: the same collapse, and this time it does (cos(pi/2)).
        Assert.True(CyclotomicRing.IsZero(CollisionGapOddOrdersClaim.TwoTimesX(24, tau, 6)));

        // n = 30 alone could not have broken it: there the collapse never vanishes.
        foreach (var t in LevelCollisionCensus.CleanTriples(30))
        {
            if (CollisionGapOddOrdersClaim.DoubledLabelRot3Cosets(30, t) is null) continue;
            if (!CollisionGapOddOrdersClaim.IsParityUniform(t)) continue;
            Assert.False(CyclotomicRing.IsZero(CollisionGapOddOrdersClaim.TwoTimesX(30, t, 3)));
        }
    }

    [Fact]
    public void TheoremE_TheCollapseClosedForm_MeetsTheRingSumExactly()
    {
        int seen = 0;
        foreach (int n in new[] { 12, 18, 24, 30 })
            foreach (var t in LevelCollisionCensus.CleanTriples(n))
            {
                if (CollisionGapOddOrdersClaim.DoubledLabelRot3Cosets(n, t) is null) continue;
                if (!CollisionGapOddOrdersClaim.IsParityUniform(t)) continue;
                foreach (int j in new[] { 3, 6, 9 })
                {
                    seen++;
                    Assert.True(CyclotomicRing.AreEqual(
                            CollisionGapOddOrdersClaim.Rot3CollapseValue(n, t, j),
                            CollisionGapOddOrdersClaim.TwoTimesX(n, t, j)),
                        $"n = {n}, tau = {t}, j = {j}: the collapse form and the ring sum differ");
                }
            }
        Assert.Equal(240, seen);
        // The three refusals: a rung the coset does not collapse at, a triple without the shape, and a ROT3 triple
        // whose labels are of mixed parity, where eta is not a single sign and the closed form has no meaning.
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGapOddOrdersClaim.Rot3CollapseValue(24, (1, 7, 9), 4));
        Assert.Throws<ArgumentException>(() => CollisionGapOddOrdersClaim.Rot3CollapseValue(24, (1, 2, 4), 3));
        Assert.Throws<ArgumentException>(() => CollisionGapOddOrdersClaim.Rot3CollapseValue(9, (1, 2, 4), 3));
    }

    // ---------------------------------------------------------------- the census, and the second order

    [Fact]
    public void TheCensus_RecomputesEveryCountTheProofReports()
    {
        var c = CollisionGapOddOrdersClaim.Census();

        Assert.Equal(30, c.MaxComb);
        Assert.Equal(2558, c.Pairs);
        Assert.Equal(2335, c.Separating);
        Assert.Equal(223, c.Standing);
        Assert.Equal(627, c.OddCombPairs);
        Assert.Equal(627, c.OddCombThirdOrderVanishes);          // Theorem D, on the whole odd population
        Assert.Equal(11, c.ThetaMirrorStanding);
        Assert.Equal(212, c.NonMirrorStanding);
        Assert.Equal(212, c.NonMirrorStandingWithSixthRungAlive); // the odd part begins exactly at fifth order
        Assert.Equal(12, c.EvenCombStandingBothRot3AndParityMatched);
        Assert.Equal(60, c.SecondOrderVanishes);
    }

    [Fact]
    public void TheSecondOrder_VanishesOnNoStandingPair_AndOnExactlyTheFamilyCPairsPerModulus()
    {
        foreach (int n in FiringCombsToThirty)
        {
            int zeros = 0;
            foreach (var (tau, sigma) in CollisionGapOddOrdersClaim.CollisionPairs(n))
            {
                bool second = CollisionGapOddOrdersClaim.SecondOrderVanishes(n, tau, sigma);
                if (second) zeros++;
                // "every one of the 223 leaves at second order": standing and c_2 = 0 never meet.
                Assert.False(second && CollisionGapOddOrdersClaim.StandsAtFirstOrder(n, tau, sigma),
                    $"n = {n}: {tau} ~ {sigma} stands at first order AND has a vanishing second order");
            }
            Assert.Equal(CollisionGapOddOrdersClaim.SecondOrderZeroPairLowerBound(n), zeros);
        }
    }

    [Fact]
    public void TheSecondOrderBound_IsTheCommittedInventorysFamiliesCAndL_DoorsIncluded()
    {
        Assert.Equal(new[] { CollisionFamilyInventory.Family.C, CollisionFamilyInventory.Family.L },
            CollisionGapOddOrdersClaim.ThreeFreeFamilies);

        Assert.Equal(0, CollisionGapOddOrdersClaim.SecondOrderZeroPairLowerBound(9));   // the door, not -2
        Assert.Equal(0, CollisionGapOddOrdersClaim.SecondOrderZeroPairLowerBound(21));
        Assert.Equal(20, CollisionGapOddOrdersClaim.SecondOrderZeroPairLowerBound(20));
        Assert.Equal(40, CollisionGapOddOrdersClaim.SecondOrderZeroPairLowerBound(30));
        Assert.Equal(140, CollisionGapOddOrdersClaim.SecondOrderZeroPairLowerBound(70));
        Assert.Equal(420, CollisionGapOddOrdersClaim.SecondOrderZeroPairLowerBound(210));
    }

    [Fact]
    public void TheThetaMirrorPairs_AreTheOnesWhoseSixthRungIsAlsoEmpty()
    {
        int mirror = 0;
        foreach (int n in FiringCombsToThirty)
            foreach (var (tau, sigma) in CollisionGapOddOrdersClaim.CollisionPairs(n))
            {
                if (!CollisionGapOddOrdersClaim.StandsAtFirstOrder(n, tau, sigma)) continue;
                bool sixthEmpty = CyclotomicRing.IsZero(CollisionGapOddOrdersClaim.DeltaTwoTimesX(n, tau, sigma, 3));
                Assert.Equal(CollisionGapOddOrdersClaim.IsThetaMirrorPair(n, tau, sigma), sixthEmpty);
                if (sixthEmpty) mirror++;
            }
        Assert.Equal(11, mirror);
    }

    // ---------------------------------------------------------------- guards and the claim's own surface

    [Fact]
    public void TheFirstOrderCoefficient_RefusesToBeRationalWhereDeltaX2DoesNotVanish()
    {
        // A mixed-parity pair at an even modulus: the claim must throw rather than return the odd-label form.
        var tau = (1, 2, 4);
        var sigma = (1, 2, 5);
        Assert.False(CollisionGapOddOrdersClaim.FirstOrderIsRational(12, tau, sigma));
        Assert.Throws<ArgumentException>(() => CollisionGapOddOrdersClaim.FirstOrderCoefficient(12, tau, sigma));
    }

    [Fact]
    public void TheGuards_RefuseAMalformedTripleAndANegativeRung()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGapOddOrdersClaim.TwoTimesX(9, (3, 2, 1), 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGapOddOrdersClaim.TwoTimesX(9, (1, 2, 9), 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGapOddOrdersClaim.TwoTimesX(9, (1, 2, 4), -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => CollisionGapOddOrdersClaim.Census(8));
    }

    [Fact]
    public void TheClaimRefusesANullParent() =>
        Assert.Throws<ArgumentNullException>(() => new CollisionGapOddOrdersClaim(null!));

    // The tier, the statement's earned fences and the parent edge are pinned on the registry instance, in
    // RCPsiSquared.Runtime.Tests/F1Family/CollisionGapOddOrdersClaimRegistrationTests: F160's own parents are four
    // claims deep, and hand-building that chain here would test the chain rather than this claim.

    // ---------------------------------------------------------------- the exact route, local to this file

    /// <summary>d_m for m = 0..6 of the branch of det(x·I − H(u)) = 0 through x0, solved by Newton in the truncated
    /// power series ring ℚ[u]/(u^7). Exact: F160's polynomial has integer coefficients and x0 is rational at every
    /// comb point used above, so nothing here leaves ℚ. Independent of the closed forms under test.</summary>
    private static BigRational[] ExactRootSeries(int nComb, BigRational x0)
    {
        int sites = nComb - 1;
        var a = TaylorAt(CrackedRingExactCurveClaim.ChebyshevSecondKindMonic(sites), x0);
        var b = TaylorAt(CrackedRingExactCurveClaim.ChebyshevSecondKindMonic(sites - 2), x0);

        var u = Constant(0); u[1] = BigRational.One;
        var u2 = Mul(u, u);

        // P(x0 + d, u) = sum_i (a_i - u^2 b_i) d^i - 2u. delta has no constant term (the root is exactly x0 at u = 0),
        // so d^i vanishes mod u^Trunc for i >= Trunc and P's sum below is the whole polynomial, not a truncation of
        // it. DP is NOT exact in the same way: it carries d^(i-1), so its i = Trunc term is dropped. That is harmless
        // and deliberate. Newton needs the derivative only to converge, and what certifies the answer is the exact
        // residual P(delta) == 0 together with delta[0] == 0, whose solution is unique because dP/ddelta has the unit
        // constant term a_1 != 0 (the branch is simple, asserted in Inv).
        BigRational[] P(BigRational[] d)
        {
            var acc = Mul(Constant(-2), u);
            var power = Constant(1);
            for (int i = 0; i < Trunc; i++)
            {
                acc = Add(acc, Mul(Coefficient(i, 1), power));
                power = Mul(power, d);
            }
            return acc;
        }

        BigRational[] DP(BigRational[] d)
        {
            var acc = Constant(0);
            var power = Constant(1);                 // d^(i-1)
            for (int i = 1; i < Trunc; i++)
            {
                acc = Add(acc, Mul(Coefficient(i, i), power));
                power = Mul(power, d);
            }
            return acc;
        }

        BigRational[] Coefficient(int i, long factor) =>
            Sub(Constant(i < a.Length ? a[i] * factor : BigRational.Zero),
                Mul(u2, Constant(i < b.Length ? b[i] * factor : BigRational.Zero)));

        var delta = Constant(0);
        for (int iteration = 0; iteration < 10; iteration++)
            delta = Sub(delta, Mul(P(delta), Inv(DP(delta))));

        // The residual is exact, so the check is an equality and not a threshold.
        Assert.All(P(delta), v => Assert.True(v.IsZero, $"Newton left a residual: {v}"));
        Assert.True(delta[0].IsZero, "the branch starts at x0");
        return delta;
    }

    /// <summary>The value of 2·ΔX_0 read out of the ring: X_0 is a rational integer, so its vector is that integer
    /// times the basis element 1 and every other coordinate is zero, which is asserted rather than assumed.</summary>
    private static long RationalIntegerValue(int nComb, (int, int, int) tau, (int, int, int) sigma)
    {
        var v = CollisionGapOddOrdersClaim.DeltaTwoTimesX(nComb, tau, sigma, 0);
        for (int i = 1; i < v.Length; i++) Assert.Equal(0, v[i]);
        return v[0];
    }

    // -------- truncated power series over Q, and one polynomial helper

    private const int Trunc = 7;

    private static BigRational[] Constant(BigRational v)
    {
        var r = new BigRational[Trunc];
        for (int i = 0; i < Trunc; i++) r[i] = BigRational.Zero;
        r[0] = v;
        return r;
    }

    private static BigRational[] Constant(long v) => Constant(new BigRational(v));

    private static BigRational[] Add(BigRational[] a, BigRational[] b)
    {
        var r = new BigRational[Trunc];
        for (int i = 0; i < Trunc; i++) r[i] = a[i] + b[i];
        return r;
    }

    private static BigRational[] Sub(BigRational[] a, BigRational[] b)
    {
        var r = new BigRational[Trunc];
        for (int i = 0; i < Trunc; i++) r[i] = a[i] - b[i];
        return r;
    }

    private static BigRational[] Mul(BigRational[] a, BigRational[] b)
    {
        var r = Constant(0);
        for (int i = 0; i < Trunc; i++)
        {
            if (a[i].IsZero) continue;
            for (int j = 0; i + j < Trunc; j++)
                if (!b[j].IsZero) r[i + j] += a[i] * b[j];
        }
        return r;
    }

    private static BigRational[] Inv(BigRational[] a)
    {
        Assert.False(a[0].IsZero, "the branch is simple, so dP/ddelta has a unit constant term");
        var r = Constant(0);
        r[0] = BigRational.One / a[0];
        for (int k = 1; k < Trunc; k++)
        {
            var s = BigRational.Zero;
            for (int i = 1; i <= k; i++) s += a[i] * r[k - i];
            r[k] = -s / a[0];
        }
        return r;
    }

    /// <summary>Taylor coefficients of an integer polynomial about x0, by repeated synthetic division.</summary>
    private static BigRational[] TaylorAt(BigInteger[] ascending, BigRational x0)
    {
        var c = ascending.Select(v => new BigRational(v)).ToArray();
        var outp = new List<BigRational>();
        for (int step = 0; step <= Trunc; step++)
        {
            if (c.Length == 0) { outp.Add(BigRational.Zero); c = Array.Empty<BigRational>(); continue; }
            var q = new BigRational[Math.Max(c.Length - 1, 0)];
            var acc = BigRational.Zero;
            for (int i = c.Length - 1; i >= 0; i--)
            {
                acc = c[i] + acc * x0;
                if (i > 0) q[i - 1] = acc;
            }
            outp.Add(acc);
            c = q;
        }
        return outp.ToArray();
    }

    private static double ToDouble(BigRational r) => (double)r.Numerator / (double)r.Denominator;
}
