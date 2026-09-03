using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>From-below pins for <see cref="BlindSeatSectorFactorisationClaim"/> (F162, proof
/// <c>docs/proofs/PROOF_BLIND_SEAT_TWO_AXES.md</c> §(i)). Every statement is exact integer or rational
/// arithmetic, so every comparison below is <c>==</c> against a literal and there is no tolerance anywhere.
///
/// <para>Five of these tests are CONTROLS that assert a break, and each predicts the COUNT of that break from
/// a stated reason before running it, rather than reading the count off the run. A resultant is antisymmetric
/// up to (−1)^(deg f·deg g), so reading the halves in seat order instead of fold order must part exactly where
/// that exponent is odd, at 6 of 50 seats outside and 24 of 100 sector readings inside; dropping Corollary
/// 11b's mirror term must part exactly where C(n,2) is odd at a mirror seat, at 12 of 50; the pole split needs
/// a comb with no repeated root, so it must break exactly where S_m and S_p share one, at 16 of 42 draws; the
/// sector parity is an assignment, so the swapped assignment must fail; and an interpolation given too small a
/// bound must throw rather than return a shorter polynomial. A control that could pass for a healthy reason is
/// not a control, and neither is one that asserts whatever it happens to find.</para>
///
/// <para>The two lowest tests do not go through the claim's laws at all: they pin the sine quotient against its
/// own recursion and against the Chebyshev addition formula, which is the one identity §(i) writes out because
/// no proof of this repository carried it.</para></summary>
public class BlindSeatSectorFactorisationClaimTests
{
    private static BlindSeatSectorFactorisationClaim BuildClaim() =>
        new(BuildBlindSeat(), BuildCrackedRing());

    private static SeatCutBlindnessClaim BuildBlindSeat() => new(new F4KernelDimensionByComponentsClaim());

    private static CrackedRingExactCurveClaim BuildCrackedRing()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var qubitAnchor = new QubitDimensionalAnchorClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubitAnchor);
        var f65 = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        var f2b = new F2bXyChainSpectrumPi2Inheritance(ladder, f65);
        var absorption = new AbsorptionTheoremClaim(ladder);
        var discriminant = new PolynomialDiscriminantAnchorClaim(new PolynomialFoundationClaim(), qubitAnchor, ladder);
        var carrier = new UniversalCarrierClaim(absorption, ladder, discriminant);
        var clock = new ClockHandLadderClaim(f2b, absorption, carrier);
        return new CrackedRingExactCurveClaim(new TopologyBandEdgeClaim(clock, absorption), f2b);
    }

    private const int SweepMax = 12;

    private static IEnumerable<(int Chain, int Seat)> SeatsInScope(int maxChain = SweepMax)
    {
        for (int chain = BlindSeatSectorFactorisationClaim.MinChain; chain <= maxChain; chain++)
        foreach (int seat in BlindSeatSectorFactorisationClaim.SeatsInScope(chain))
            yield return (chain, seat);
    }

    // -------------------------------------------------------------------------------------------
    // the polynomials themselves, below every law
    // -------------------------------------------------------------------------------------------

    [Fact]
    public void TheSineQuotient_IsItsOwnRecursion_AndCarriesTheBackwardStep()
    {
        Assert.Equal(new BigInteger[] { 0 }, BlindSeatSectorFactorisationClaim.SineQuotient(0));
        Assert.Equal(new BigInteger[] { 1 }, BlindSeatSectorFactorisationClaim.SineQuotient(1));
        Assert.Equal(new BigInteger[] { 0, 1 }, BlindSeatSectorFactorisationClaim.SineQuotient(2));
        Assert.Equal(new BigInteger[] { -1, 0, 1 }, BlindSeatSectorFactorisationClaim.SineQuotient(3));
        Assert.Equal(new BigInteger[] { 0, -2, 0, 1 }, BlindSeatSectorFactorisationClaim.SineQuotient(4));
        Assert.Equal(new BigInteger[] { -1 }, BlindSeatSectorFactorisationClaim.SineQuotient(-1));

        for (int m = 1; m <= 12; m++)
        {
            var up = BlindSeatSectorFactorisationClaim.SineQuotient(m + 1);
            var here = BlindSeatSectorFactorisationClaim.SineQuotient(m);
            var down = BlindSeatSectorFactorisationClaim.SineQuotient(m - 1);
            for (int x = -3; x <= 3; x++)
                Assert.Equal(Evaluate(up, x), x * Evaluate(here, x) - Evaluate(down, x));
            Assert.Equal(BigInteger.One, up[^1]);            // monic
            Assert.Equal(m, up.Length - 1);                  // of degree m
        }
    }

    [Fact]
    public void TheAdditionFormula_TheOneIdentityTheProofWritesOut_HoldsAtEveryPairInRange()
    {
        // S_{a+b} = S_a*S_{b+1} - S_{a-1}*S_b, the step Lemma 9 uses twice.
        for (int a = 1; a <= 8; a++)
        for (int b = 0; b <= 8; b++)
        for (int x = -3; x <= 3; x++)
        {
            var lhs = Evaluate(BlindSeatSectorFactorisationClaim.SineQuotient(a + b), x);
            var rhs = Evaluate(BlindSeatSectorFactorisationClaim.SineQuotient(a), x)
                        * Evaluate(BlindSeatSectorFactorisationClaim.SineQuotient(b + 1), x)
                      - Evaluate(BlindSeatSectorFactorisationClaim.SineQuotient(a - 1), x)
                        * Evaluate(BlindSeatSectorFactorisationClaim.SineQuotient(b), x);
            Assert.Equal(lhs, rhs);
        }
    }

    [Fact]
    public void TheTwoSectorCombs_MultiplyBackToTheMiddleRoute()
    {
        for (int n = 1; n <= 14; n++)
        {
            var even = BlindSeatSectorFactorisationClaim.SectorComb(n, ReflectionSector.Even);
            var odd = BlindSeatSectorFactorisationClaim.SectorComb(n, ReflectionSector.Odd);
            Assert.Equal(BlindSeatSectorFactorisationClaim.SineQuotient(n), Multiply(even, odd));
            Assert.Equal(BigInteger.One, even[^1]);
            Assert.Equal(BigInteger.One, odd[^1]);
        }
    }

    // -------------------------------------------------------------------------------------------
    // the two integers a seat is, and the fence
    // -------------------------------------------------------------------------------------------

    [Fact]
    public void TheFoldCoordinateAndTheNodeModulus_AddBackToTheChain_AndAMirrorSeatCarriesTheSamePair()
    {
        for (int chain = 4; chain <= SweepMax; chain++)
        for (int seat = 0; seat < chain; seat++)
        {
            int p = BlindSeatSectorFactorisationClaim.FoldCoordinate(chain, seat);
            int n = BlindSeatSectorFactorisationClaim.NodeModulus(chain, seat);
            Assert.Equal(chain - 1, 2 * p + n);
            Assert.Equal(p, BlindSeatSectorFactorisationClaim.FoldCoordinate(chain, chain - 1 - seat));
            Assert.Equal(n, BlindSeatSectorFactorisationClaim.NodeModulus(chain, chain - 1 - seat));
        }
    }

    [Fact]
    public void TheFencedSeats_AreRefused_AndTheyAreTheEndsAndTheOddChainsCentre()
    {
        Assert.False(BlindSeatSectorFactorisationClaim.IsInScope(9, 0));
        Assert.False(BlindSeatSectorFactorisationClaim.IsInScope(9, 8));
        Assert.False(BlindSeatSectorFactorisationClaim.IsInScope(9, 4));      // the reflection-fixed centre
        Assert.True(BlindSeatSectorFactorisationClaim.IsInScope(9, 3));
        Assert.DoesNotContain(4, BlindSeatSectorFactorisationClaim.SeatsInScope(9));
        Assert.Equal(new[] { 1, 2, 3, 5, 6, 7 }, BlindSeatSectorFactorisationClaim.SeatsInScope(9));
        Assert.Equal(new[] { 1, 2, 3, 4, 5, 6 }, BlindSeatSectorFactorisationClaim.SeatsInScope(8));

        Assert.Throws<ArgumentOutOfRangeException>(() => BlindSeatSectorFactorisationClaim.CongruenceIsExact(9, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => BlindSeatSectorFactorisationClaim.OuterResultant(9, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => BlindSeatSectorFactorisationClaim.FoldCoordinate(3, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => BlindSeatSectorFactorisationClaim.FoldCoordinate(9, 9));
    }

    // -------------------------------------------------------------------------------------------
    // Lemmas 9 to 11 and the corollaries, at every seat in scope
    // -------------------------------------------------------------------------------------------

    [Fact]
    public void Lemma9_TheCongruence_IsExactAtEverySeatInScope()
    {
        foreach (var (chain, seat) in SeatsInScope())
            Assert.True(BlindSeatSectorFactorisationClaim.CongruenceIsExact(chain, seat),
                        $"the congruence broke at N = {chain}, seat {seat}");
    }

    [Fact]
    public void Lemma10_TheCommonFactor_IsTheSignAndCarriesNoKnob()
    {
        for (int p = 1; p <= 7; p++)
        {
            var factor = BlindSeatSectorFactorisationClaim.CommonFactor(p);
            Assert.Single(factor);                                                    // no t in it
            Assert.Equal(BlindSeatSectorFactorisationClaim.Pairs(p) % 2 == 0 ? BigInteger.One : BigInteger.MinusOne,
                         factor[0]);
        }
    }

    [Fact]
    public void Corollary10a_TheOuterSignLaw_HoldsAtEverySeatInScope()
    {
        foreach (var (chain, seat) in SeatsInScope())
            Assert.True(BlindSeatSectorFactorisationClaim.OuterSignLawHolds(chain, seat),
                        $"the outer sign law broke at N = {chain}, seat {seat}");
    }

    [Fact]
    public void Corollary10a_TheExponentReadsTheFoldCoordinate_AndNotTheChain()
    {
        // Two seats of DIFFERENT chains sharing a fold coordinate must carry the same sign, and two seats of
        // one chain with fold coordinates of different parity in C(p+1,2) must not.
        var byFold = new Dictionary<int, List<(int Chain, int Seat)>>();
        foreach (var (chain, seat) in SeatsInScope())
        {
            int p = BlindSeatSectorFactorisationClaim.FoldCoordinate(chain, seat);
            if (!byFold.TryGetValue(p, out var list)) byFold[p] = list = new List<(int, int)>();
            list.Add((chain, seat));
        }

        foreach (var (p, seats) in byFold)
        {
            var expected = BlindSeatSectorFactorisationClaim.Pairs(p + 1) % 2 == 0 ? 1 : -1;
            foreach (var (chain, seat) in seats)
                Assert.Equal(expected, SignOfRatio(chain, seat));
        }

        // The reading is only about N where one fold coordinate is carried by several chains. The largest fold
        // coordinates of a bounded sweep sit on one chain each, so the demand is made where it can be met.
        foreach (int p in new[] { 1, 2, 3, 4 })
            Assert.True(byFold[p].Select(s => s.Chain).Distinct().Count() > 1,
                        $"fold coordinate {p} is carried by one chain only, so it would not be reading N");

        Assert.Contains(byFold.Keys, p => BlindSeatSectorFactorisationClaim.Pairs(p + 1) % 2 == 0);
        Assert.Contains(byFold.Keys, p => BlindSeatSectorFactorisationClaim.Pairs(p + 1) % 2 != 0);
    }

    [Fact]
    public void Lemma11_TheNodeIdentity_AndTheSectorParity_AreExact()
    {
        for (int n = 1; n <= 14; n++)
        {
            for (int p = 0; p <= 8; p++)
                Assert.True(BlindSeatSectorFactorisationClaim.NodeIdentityIsExact(n, p),
                            $"the node identity broke at n = {n}, p = {p}");
            Assert.True(BlindSeatSectorFactorisationClaim.SectorParityIsExact(n),
                        $"the sector parity broke at n = {n}");
        }
    }

    [Fact]
    public void Corollary10b_ThePoleSplit_TheDegreeAndTheLeadingCoefficient_HoldAtEverySeatInScope()
    {
        foreach (var (chain, seat) in SeatsInScope())
        foreach (var sector in new[] { ReflectionSector.Even, ReflectionSector.Odd })
        {
            Assert.True(BlindSeatSectorFactorisationClaim.PoleSplitIdentityHolds(chain, seat, sector),
                        $"the pole split broke at N = {chain}, seat {seat}, sector {sector}");
            Assert.True(BlindSeatSectorFactorisationClaim.DegreeAndLeadingCoefficientHold(chain, seat, sector),
                        $"the degree or the leading coefficient broke at N = {chain}, seat {seat}, sector {sector}");
            // c_S is pinned on ONE factor and not on the product c_E*c_O, because a sign error common to
            // both sectors cancels there: the composition reads only their product, so it is blind to
            // exactly the mutation that drops c_S's own (-1)^(p*r_S + n_S). This line is where that sign
            // has to be right.
            var constant = BlindSeatSectorFactorisationClaim.SectorConstant(chain, seat, sector);
            Assert.NotEqual(BigInteger.Zero, constant);
            Assert.Equal(BlindSeatSectorFactorisationClaim.SectorHalvesResultant(chain, seat, sector)[^1],
                         constant);
        }
    }

    [Fact]
    public void Corollary11b_TheTieToF157sOwnGenerator_HoldsAtEverySeatInScope()
    {
        foreach (var (chain, seat) in SeatsInScope())
            Assert.True(BlindSeatSectorFactorisationClaim.GeneratorTieHolds(chain, seat),
                        $"the generator tie broke at N = {chain}, seat {seat}");
    }

    [Fact]
    public void TheComposition_OfTheTwoConstants_HoldsAtEverySeatInScope()
    {
        foreach (var (chain, seat) in SeatsInScope())
            Assert.True(BlindSeatSectorFactorisationClaim.CompositionHolds(chain, seat),
                        $"the composition broke at N = {chain}, seat {seat}");
    }

    // -------------------------------------------------------------------------------------------
    // the repeated factor, which is what "one factor per root" says and "every factor simple" would not
    // -------------------------------------------------------------------------------------------

    [Fact]
    public void ARepeatedFactor_IsReal_AndTheFirstOneIsWhereTheProofSaysItIs()
    {
        // N = 10 seat 2: Q_E = -(t-1)^2, from the two odd node indices k = 1 and k = 3.
        var qEven = BlindSeatSectorFactorisationClaim.SectorHalvesResultant(10, 2, ReflectionSector.Even);
        Assert.Equal(new BigInteger[] { -1, 2, -1 }, qEven);
        Assert.Equal(2, BlindSeatSectorFactorisationClaim.LargestRepeatedFactorMultiplicity(10, 2, ReflectionSector.Even));

        // N = 14 seat 3 carries a multiplicity of three.
        Assert.Equal(3, BlindSeatSectorFactorisationClaim.LargestRepeatedFactorMultiplicity(14, 3, ReflectionSector.Even));

        // Below N = 10 nothing in the block repeats, which is why a narrowed sweep could not see any of this.
        foreach (var (chain, seat) in SeatsInScope(9))
        foreach (var sector in new[] { ReflectionSector.Even, ReflectionSector.Odd })
            Assert.True(BlindSeatSectorFactorisationClaim.LargestRepeatedFactorMultiplicity(chain, seat, sector) <= 1,
                        $"N = {chain}, seat {seat}, sector {sector} repeats below N = 10");
    }

    [Fact]
    public void PastTheGatesRange_TheMultiplicityIsStillRead_WhereTheLeadingCoefficientIsNotPlusOrMinusOne()
    {
        // At N = 19 seat 3 the even sector gives Q_E = 4t^6 - 12t^4 + 9t^2 = t^2*(2t^2 - 3)^2, whose leading
        // KNOB-coefficient is 4. Reading the repeated factor through a monic gcd asks for t^3 - (3/2)t there
        // and leaves the integers; the primitive gcd does not. The statements are algebraic in p and N, so
        // the class must answer past the range the gate swept, and this is the first seat that shows it.
        Assert.Equal(new BigInteger[] { 0, 0, 9, 0, -12, 0, 4 },
                     BlindSeatSectorFactorisationClaim.SectorHalvesResultant(19, 3, ReflectionSector.Even));
        Assert.Equal(new BigInteger(4),
                     BlindSeatSectorFactorisationClaim.SectorConstant(19, 3, ReflectionSector.Even));
        Assert.Equal(2, BlindSeatSectorFactorisationClaim.LargestRepeatedFactorMultiplicity(
                            19, 3, ReflectionSector.Even));
        Assert.True(BlindSeatSectorFactorisationClaim.OuterSignLawHolds(19, 3));
        Assert.True(BlindSeatSectorFactorisationClaim.GeneratorTieHolds(19, 3));
        Assert.True(BlindSeatSectorFactorisationClaim.CompositionHolds(19, 3));
    }

    // -------------------------------------------------------------------------------------------
    // controls: each one names why the break is forced
    // -------------------------------------------------------------------------------------------

    [Fact]
    public void Control_ReadingTheHalvesInSeatOrderInsteadOfFoldOrder_MustBreakBothLaws_AtPinnedCounts()
    {
        // "Fold half first" is the headline fence, so it is read on the LAWS and not on the antisymmetry of
        // the resultant routine, which that routine satisfies in one recursion step by construction and
        // which no statement of section (i) depends on. Both counts are predicted from the exponent
        // (-1)^(deg f * deg g) and not read off the run: the outer resultant differs only at a mirror seat,
        // and there only when p*(N-1-p) is odd, which is 6 of the 50 seats over N = 4..12; a sector reading
        // differs when p*deg(beta_S) is odd, which is 24 of the 100 sector readings over the same range.
        int outerBroken = 0, outerStood = 0, sectorBroken = 0, sectorStood = 0;
        foreach (var (chain, seat) in SeatsInScope())
        {
            int p = BlindSeatSectorFactorisationClaim.FoldCoordinate(chain, seat);
            var product = Multiply(
                BlindSeatSectorFactorisationClaim.SectorHalvesResultant(chain, seat, ReflectionSector.Even),
                BlindSeatSectorFactorisationClaim.SectorHalvesResultant(chain, seat, ReflectionSector.Odd));
            var signed = product
                .Select(c => c * (BlindSeatSectorFactorisationClaim.Pairs(p + 1) % 2 == 0
                                      ? BigInteger.One : BigInteger.MinusOne))
                .ToArray();
            var seatOrder = BlindSeatSectorFactorisationClaim.OuterResultant(chain, seat, foldHalfFirst: false);
            if (seatOrder.SequenceEqual(signed)) outerStood++; else outerBroken++;

            foreach (var sector in new[] { ReflectionSector.Even, ReflectionSector.Odd })
            {
                var comb = BlindSeatSectorFactorisationClaim.SectorComb(
                    BlindSeatSectorFactorisationClaim.NodeModulus(chain, seat), sector);
                var swapped = BlindSeatSectorFactorisationClaim.ReadingHoldsFor(p, comb, foldHalfFirst: false);
                if (swapped.DegreeAndLeadHold) sectorStood++; else sectorBroken++;
            }
        }
        Assert.Equal(6, outerBroken);
        Assert.Equal(44, outerStood);
        Assert.Equal(24, sectorBroken);
        Assert.Equal(76, sectorStood);

        // and the routine's own antisymmetry, which is what makes an argument order a statement at all
        Assert.Equal(1, BlindSeatSectorFactorisationClaim.FoldCoordinate(9, 1));
        var foldFirst = BlindSeatSectorFactorisationClaim.Resultant(
            BlindSeatSectorFactorisationClaim.KnobChainAt(1, 3),
            ToPoly(BlindSeatSectorFactorisationClaim.KnobChainAt(7, 3)));
        var swappedOrder = BlindSeatSectorFactorisationClaim.Resultant(
            BlindSeatSectorFactorisationClaim.KnobChainAt(7, 3),
            ToPoly(BlindSeatSectorFactorisationClaim.KnobChainAt(1, 3)));
        Assert.NotEqual(BigInteger.Zero, foldFirst.Numerator);
        Assert.Equal(foldFirst, -swappedOrder);
    }

    [Fact]
    public void Control_DroppingCorollary11bsMirrorTerm_MustBreakTheTie_AtAPinnedCountOfSeats()
    {
        // Corollary 11b's last summand [j > N-1-j]*C(n,2) is the price of carrying F157's SEAT index onto the
        // fold coordinate. This control recomputes the tie with that summand deleted and requires the two
        // outcomes the parity predicts: it must break at the mirror seats where C(n,2) is odd, 12 of the 50
        // seats over N = 4..12, and stand at the other 38. A term with no effect reddens the first literal, a
        // term that broke the tie everywhere reddens the second. Counting parities alone, which an earlier
        // version of this control did, asserted neither and survived deleting the very term it names.
        int broken = 0, survived = 0;
        foreach (var (chain, seat) in SeatsInScope())
        {
            int n = BlindSeatSectorFactorisationClaim.NodeModulus(chain, seat);
            int without = BlindSeatSectorFactorisationClaim.GeneratorTieExponent(chain, seat)
                          - (seat > chain - 1 - seat ? BlindSeatSectorFactorisationClaim.Pairs(n) : 0);
            if (TieHoldsAtExponent(chain, seat, without)) survived++; else broken++;
        }
        Assert.Equal(12, broken);
        Assert.Equal(38, survived);
    }

    [Fact]
    public void Control_TheFactorisationNeedsASquarefreeComb_AndBreaksWithoutOne()
    {
        // Corollary 10b is stated for a monic SQUAREFREE comb and the word is load-bearing: with beta = S_m
        // squared, the pole split leaves the shared root in h as well as in g, so Res(h, S_p) vanishes and
        // the constant with it, while the halves-resultant does not. That happens exactly when S_m and S_p
        // share a root, which is when gcd(m, p) > 1, and over m = 2..8 and p = 1..6 that is 16 of the 42
        // draws. The count is predicted from the shared-root condition, not read off the run. The squarefree
        // half must not break at all, so a reading that broke everywhere reddens here rather than passing as
        // a strong control.
        int squarefreeDraws = 0, squarefreeFailures = 0, squaredDraws = 0, squaredFailures = 0;
        for (int fold = 1; fold <= 6; fold++)
        for (int m = 2; m <= 8; m++)
        {
            var comb = BlindSeatSectorFactorisationClaim.SineQuotient(m);
            squarefreeDraws++;
            var plain = BlindSeatSectorFactorisationClaim.ReadingHoldsFor(fold, comb);
            if (!plain.ProductHolds || !plain.DegreeAndLeadHold) squarefreeFailures++;

            squaredDraws++;
            var squared = BlindSeatSectorFactorisationClaim.ReadingHoldsFor(fold, Multiply(comb, comb));
            if (!squared.DegreeAndLeadHold) squaredFailures++;
        }
        Assert.Equal(42, squarefreeDraws);
        Assert.Equal(0, squarefreeFailures);
        Assert.Equal(42, squaredDraws);
        Assert.Equal(16, squaredFailures);
    }

    [Fact]
    public void Control_TheSwappedSectorAssignment_MustFail()
    {
        // beta_E carries the ODD node indices. If the assignment were the other way round, S_{n+1} + 1 would die
        // modulo beta_O; at n = 3 the two combs are t - 1 and t + 1 and it plainly does not.
        var top = BlindSeatSectorFactorisationClaim.SineQuotient(4);                 // S_4 = x^3 - 2x
        var oddIndices = BlindSeatSectorFactorisationClaim.SectorComb(3, ReflectionSector.Even);
        var evenIndices = BlindSeatSectorFactorisationClaim.SectorComb(3, ReflectionSector.Odd);
        Assert.Equal(new BigInteger[] { -1, 1 }, oddIndices);
        Assert.Equal(new BigInteger[] { 1, 1 }, evenIndices);
        Assert.Equal(BigInteger.MinusOne, Evaluate(top, 1));      // S_{n+1}(x_1) = -1, k = 1 odd
        Assert.Equal(BigInteger.One, Evaluate(top, -1));          // S_{n+1}(x_2) = +1, k = 2 even

        // and the break itself, which the two lines above only make plausible: handing the law the OTHER
        // assignment must fail, at every modulus where the two combs differ. They coincide only at n = 1,
        // where both are the constant 1 and the statement is empty, so the swapped reading survives there
        // and nowhere else, and that split is asserted rather than the failures alone.
        int failed = 0, survived = 0;
        for (int n = 1; n <= 14; n++)
        {
            Assert.True(BlindSeatSectorFactorisationClaim.SectorParityIsExact(n, ReflectionSector.Even));
            if (BlindSeatSectorFactorisationClaim.SectorParityIsExact(n, ReflectionSector.Odd)) survived++;
            else failed++;
        }
        Assert.Equal(13, failed);
        Assert.Equal(1, survived);
    }

    [Fact]
    public void Control_AnUnderestimatedKnobDegree_MustThrowRatherThanReturnAShorterAnswer()
    {
        // Every resultant here is recovered from evaluations at integer knobs, so a bound below the true
        // knob-degree would return a shorter polynomial that agrees at the sample points and nowhere else.
        // The interpolation holds two points back to stop that, and this hands it a bound that is too small
        // and requires the throw. Res(alpha_1, alpha_4) at N = 6 seat 1 has knob-degree 2.
        var honest = BlindSeatSectorFactorisationClaim.OuterResultant(6, 1);
        Assert.Equal(2, honest.Length - 1);

        BigRational AtKnob(BigInteger knob) => BlindSeatSectorFactorisationClaim.Resultant(
            BlindSeatSectorFactorisationClaim.KnobChainAt(1, knob),
            ToPoly(BlindSeatSectorFactorisationClaim.KnobChainAt(4, knob)));

        Assert.Equal(honest, BlindSeatSectorFactorisationClaim.InKnob(2, AtKnob));
        Assert.Throws<InvalidOperationException>(() => BlindSeatSectorFactorisationClaim.InKnob(1, AtKnob));
        Assert.Throws<InvalidOperationException>(() => BlindSeatSectorFactorisationClaim.InKnob(0, AtKnob));
    }

    // -------------------------------------------------------------------------------------------
    // the claim's own face
    // -------------------------------------------------------------------------------------------

    [Fact]
    public void TheSurvey_IsRecomputed_AndEveryLawHoldsAtEverySeatItCounts()
    {
        var survey = BlindSeatSectorFactorisationClaim.Survey(SweepMax);
        Assert.Equal(SweepMax, survey.MaxChain);
        Assert.True(survey.Seats > 0);
        Assert.Equal(survey.Seats, survey.CongruenceHolds);
        Assert.Equal(survey.Seats, survey.CommonFactorHolds);
        Assert.Equal(survey.Seats, survey.OuterSignLawHolds);
        Assert.Equal(survey.Seats, survey.GeneratorTieHolds);
        Assert.Equal(survey.Seats, survey.PoleSplitHolds);
        Assert.Equal(survey.Seats, survey.DegreeAndLeadHolds);
        Assert.Equal(survey.Seats, survey.CompositionHolds);
        Assert.Equal(survey.Seats, survey.SeatsByFold.Values.Sum());
        Assert.True(survey.MirrorSeats > 0 && survey.MirrorSeats < survey.Seats);
    }

    [Fact]
    public void TheSweepReproducesTheGatesOwnPinnedPopulations_AndCarriesTheLawsPastTheSeatBySeatRange()
    {
        // The seat-by-seat tests above stop at N = 12. The survey runs to N = 14, which is where the gate's
        // populations are pinned, so these counts are the only place N = 13 and N = 14 are asserted at all.
        // The seat count and the fold profile are bookkeeping on min(j, N-1-j) and touch no resultant on
        // either road: they fence the RANGE. The 20 and the 16 are read through the resultants, and there the
        // two roads are genuinely different, Euclidean remainder sequences at integer knobs here against
        // Sylvester determinants in Z[t][x] in the gate.
        var survey = BlindSeatSectorFactorisationClaim.Survey(BlindSeatSectorFactorisationClaim.GateMaxChain);
        Assert.Equal(72, survey.Seats);                                              // W1's seat count
        Assert.Equal(72, survey.CongruenceHolds);
        Assert.Equal(72, survey.CommonFactorHolds);
        Assert.Equal(72, survey.OuterSignLawHolds);                                  // W3
        Assert.Equal(72, survey.GeneratorTieHolds);                                  // W7c
        Assert.Equal(72, survey.PoleSplitHolds);
        Assert.Equal(72, survey.DegreeAndLeadHolds);
        Assert.Equal(72, survey.CompositionHolds);                                   // W7
        Assert.Equal(72, survey.NodeIdentityHolds);
        Assert.Equal(72, survey.SectorParityHolds);
        Assert.Equal(new Dictionary<int, int> { [1] = 22, [2] = 18, [3] = 14, [4] = 10, [5] = 6, [6] = 2 },
                     survey.SeatsByFold);                                            // W3b's pinned profile
        Assert.Equal(20, survey.PoleShortfallReadings);                              // W6's pinned shortfall

        // W6b pins a SET and not a count, so a sector mislabel that permuted it would survive a count alone.
        // The gate's list is in its own sector letters; what is convention-free between the two roads is the
        // (N, seat) multiset, and that is what is pinned.
        Assert.Equal(16, survey.RepeatedFactorSeats.Count);
        Assert.Equal(new[] { (10, 2), (10, 2), (10, 7), (10, 7), (11, 2), (11, 3), (11, 7), (11, 8),
                             (14, 3), (14, 3), (14, 4), (14, 4), (14, 9), (14, 9), (14, 10), (14, 10) },
                     survey.RepeatedFactorSeats.Select(r => (r.Chain, r.Seat))
                           .OrderBy(r => r.Chain).ThenBy(r => r.Seat).ToArray());
    }

    [Fact]
    public void TheSeatBySeatRangeIsFenced_SoItCannotBeNarrowedInSilence()
    {
        // Every test above sweeps N = 4..SweepMax. Without these literals the range could be narrowed to 11
        // with nothing going red, because the only range-sensitive assertion up there, that each of the fold
        // coordinates 1 to 4 is carried by more than one chain, still holds at 11.
        Assert.Equal(12, SweepMax);
        Assert.Equal(50, BlindSeatSectorFactorisationClaim.Survey(SweepMax).Seats);
        Assert.Equal(40, BlindSeatSectorFactorisationClaim.Survey(11).Seats);
    }

    [Fact]
    public void TheClaim_IsTier1Derived_KeepsItsFences_AndRefusesANullParent()
    {
        // The Name and Anchor fences are pinned once, in BlindSeatSectorFactorisationClaimRegistrationTests,
        // off the registry rather than off a hand-built claim. A second copy here is a second place for a
        // repair to miss.
        var claim = BuildClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);

        Assert.Throws<ArgumentNullException>(() =>
            new BlindSeatSectorFactorisationClaim(null!, BuildCrackedRing()));
        Assert.Throws<ArgumentNullException>(() =>
            new BlindSeatSectorFactorisationClaim(BuildBlindSeat(), null!));
    }

    [Fact]
    public void TheChildren_ReadLive_AndNameWhatIsNotDerivedHere()
    {
        var children = BuildClaim().Children.ToList();
        Assert.Contains(children, c => c.DisplayName == "tier");
        Assert.Contains(children, c => c.DisplayName == "anchor");
        var scope = Assert.Single(children, c => c.DisplayName.StartsWith("scope"));
        Assert.Contains("READ and not derived", scope.Summary);
        Assert.Contains("hop-2", scope.Summary);
        Assert.Contains("ring", scope.Summary);

        // The children interpolate ten live readings, two of which, N = 13 seat 9 and N = 14 seat 3, sit
        // outside the seat-by-seat range. They are asserted here directly rather than sniffed for the string
        // "False", which would also fire on the word appearing in prose and would stop protecting anything
        // the moment a summary stopped printing a boolean.
        Assert.True(BlindSeatSectorFactorisationClaim.CongruenceIsExact(11, 2));
        Assert.True(BlindSeatSectorFactorisationClaim.CongruenceIsExact(14, 3));
        Assert.True(BlindSeatSectorFactorisationClaim.OuterSignLawHolds(11, 2));
        Assert.True(BlindSeatSectorFactorisationClaim.OuterSignLawHolds(12, 4));
        Assert.True(BlindSeatSectorFactorisationClaim.OuterSignLawHolds(13, 9));
        Assert.True(BlindSeatSectorFactorisationClaim.NodeIdentityIsExact(5, 3));
        Assert.True(BlindSeatSectorFactorisationClaim.SectorParityIsExact(5));
        Assert.True(BlindSeatSectorFactorisationClaim.SectorParityIsExact(8));
        Assert.True(BlindSeatSectorFactorisationClaim.GeneratorTieHolds(9, 2));
        Assert.True(BlindSeatSectorFactorisationClaim.GeneratorTieHolds(9, 6));
        Assert.DoesNotContain(children, c => c.Summary.Contains("False"));
    }

    // -------------------------------------------------------------------------------------------

    /// <summary>Corollary 11b's tie read at an exponent handed in, so a control can hand in a wrong one.</summary>
    private static bool TieHoldsAtExponent(int chain, int seat, int exponent)
    {
        var outer = BlindSeatSectorFactorisationClaim.OuterResultant(chain, seat);
        var generator = BlindSeatSectorFactorisationClaim.SeatGenerator(chain, seat);
        var sign = exponent % 2 == 0 ? BigInteger.One : BigInteger.MinusOne;
        return outer.SequenceEqual(generator.Select(c => c * sign).ToArray());
    }

    private static int SignOfRatio(int chain, int seat)
    {
        var outer = BlindSeatSectorFactorisationClaim.OuterResultant(chain, seat);
        var product = Multiply(BlindSeatSectorFactorisationClaim.SectorHalvesResultant(chain, seat, ReflectionSector.Even),
                               BlindSeatSectorFactorisationClaim.SectorHalvesResultant(chain, seat, ReflectionSector.Odd));
        if (outer.SequenceEqual(product)) return 1;
        if (outer.SequenceEqual(product.Select(c => -c).ToArray())) return -1;
        throw new InvalidOperationException(
            $"at N = {chain}, seat {seat} the outer resultant is not the product up to a sign at all.");
    }

    private static BigInteger Evaluate(BigInteger[] poly, BigInteger at)
    {
        var acc = BigInteger.Zero;
        for (int k = poly.Length - 1; k >= 0; k--) acc = acc * at + poly[k];
        return acc;
    }

    private static BigInteger[] Multiply(BigInteger[] a, BigInteger[] b)
    {
        var res = new BigInteger[a.Length + b.Length - 1];
        for (int i = 0; i < a.Length; i++)
        for (int k = 0; k < b.Length; k++)
            res[i + k] += a[i] * b[k];
        int last = res.Length - 1;
        while (last > 0 && res[last].IsZero) last--;
        return res.Take(last + 1).ToArray();
    }

    private static Core.Numerics.RationalPolynomial ToPoly(BigInteger[] poly) =>
        new(poly.Select(c => new Core.Numerics.BigRational(c)));
}
