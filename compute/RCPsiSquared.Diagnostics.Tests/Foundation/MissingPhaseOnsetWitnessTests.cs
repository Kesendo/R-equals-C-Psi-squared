using System;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The exact-arithmetic checks §9 of <c>experiments/THE_MOTION_AND_THE_MISSING_PHASE.md</c>
/// names, run here rather than by hand.
/// <para>How independent the oracle is, exactly, since the answer differs by reading. The PROPAGATION
/// is independent: §8 counts hops on the page and the recurrence iterates the generator here. Of the
/// three readout maps, only the leakage is independent too, computing Tr(QδA) off the coefficient
/// matrix where §8.4 runs a Dyson integral; <c>TwoSiteDistance</c> and <c>DecodedDistanceSquared</c>
/// evaluate the reductions §8.2 and §8.3 derive, so on those two the comparison confirms the
/// propagation and the coefficients rather than the reduction. The reduction has its own gates: the
/// two-block form against the full Lindbladian, and the readout maps against hand-built inputs that
/// make their unexercised channels nonzero.</para>
/// <para>Every assertion is exact equality in ℚ or ℚ(i); there is no tolerance in this file, because
/// both sides are the same rational.</para></summary>
public class MissingPhaseOnsetWitnessTests
{
    private static readonly BigRational Eps = new(1, 3);
    private static readonly BigRational Gam = new(2, 5);

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(9)]
    [InlineData(11)]
    [InlineData(13)]
    public void EveryReadingMeetsThePageClosedFormExactly(int n)
    {
        var r = MissingPhaseOnsetWitness.Read(n, Eps, Gam);
        int m = r.M;

        Assert.Equal(m + 1, r.DOutOrder);
        Assert.True(r.DOutCoefficientSquared == MissingPhaseOnsetWitness.PageDOutSquared(m, Eps, Gam),
            $"N={n}: d_out² coefficient {r.DOutCoefficientSquared}, page gives " +
            $"{MissingPhaseOnsetWitness.PageDOutSquared(m, Eps, Gam)}.");

        Assert.Equal(m % 2 == 0 ? m + 1 : m + 2, r.DTwoOrder);
        Assert.True(r.DTwoCoefficient == MissingPhaseOnsetWitness.PageDTwo(m, Eps, Gam),
            $"N={n}: d_2 coefficient {r.DTwoCoefficient}, page gives " +
            $"{MissingPhaseOnsetWitness.PageDTwo(m, Eps, Gam)}.");

        Assert.Equal(2 * m + 1, r.LeakageOrder);
        Assert.True(r.LeakageCoefficient == MissingPhaseOnsetWitness.PageLeakage(m, Eps, Gam),
            $"N={n}: leakage coefficient {r.LeakageCoefficient}, page gives " +
            $"{MissingPhaseOnsetWitness.PageLeakage(m, Eps, Gam)}.");

        Assert.True(r.CentreEndA == MissingPhaseOnsetWitness.PageCentreEndA(m, Eps, Gam),
            $"N={n}: dA[c,0] = {r.CentreEndA}, page gives {MissingPhaseOnsetWitness.PageCentreEndA(m, Eps, Gam)}.");
        Assert.True(r.CentreEndB == MissingPhaseOnsetWitness.PageCentreEndB(m, Eps, Gam),
            $"N={n}: dB[c,0] = {r.CentreEndB}, page gives {MissingPhaseOnsetWitness.PageCentreEndB(m, Eps, Gam)}.");
    }

    [Theory]
    [InlineData(5, +1)]
    [InlineData(7, -1)]
    [InlineData(9, +1)]
    [InlineData(11, -1)]
    [InlineData(13, +1)]
    public void TheLeakageDifferenceAlternatesInSignWithM(int n, int expectedSign)
    {
        var r = MissingPhaseOnsetWitness.Read(n, Eps, Gam);
        Assert.Equal(expectedSign, r.LeakageCoefficient.Sign);
    }

    [Fact]
    public void DTwoLosesExactlyOneOrderToDOutAtOddM()
    {
        // The page's §8.2 claim, read as a difference rather than as two separate numbers: at even m the
        // two distances start together, at odd m the two-site readout starts one order later.
        foreach (var n in MissingPhaseOnsetWitness.PageSites)
        {
            var r = MissingPhaseOnsetWitness.Read(n, Eps, Gam);
            Assert.Equal(r.M % 2 == 0 ? 0 : 1, r.DTwoOrder - r.DOutOrder);
        }
    }

    [Theory]
    [InlineData(-1, 3)]   // ε = −1/3
    [InlineData(-1, 1)]   // ε = −1, the severed end bond
    [InlineData(-3, 1)]   // ε = −3, the sign-reversed end bond
    public void TheOtherCouplingsSection9NamesAlsoMeetThePage(int numerator, int denominator)
    {
        var epsilon = new BigRational(numerator, denominator);
        var r = MissingPhaseOnsetWitness.Read(7, epsilon, Gam);
        Assert.True(MissingPhaseOnsetWitness.Agrees(r, epsilon, Gam),
            $"N=7 at eps={epsilon}: orders (d_out {r.DOutOrder}, d_2 {r.DTwoOrder}, leak {r.LeakageOrder}); " +
            $"d_2 = {r.DTwoCoefficient} against page {MissingPhaseOnsetWitness.PageDTwo(r.M, epsilon, Gam)}.");
    }

    [Fact]
    public void TheOddMBranchActuallyUsesItsStretchFactor()
    {
        // The odd-m branch carries max(1, |1+ε|), which is not a constant dressed up as a formula. ε = +3
        // and ε = −3 have the SAME |ε|, so every other factor of the branch cancels and the ratio of the
        // two coefficients is the ratio of the stretch alone: max(1,4)/max(1,2) = 2, exactly.
        var atPlusThree = MissingPhaseOnsetWitness.Read(7, new BigRational(3), Gam);
        var atMinusThree = MissingPhaseOnsetWitness.Read(7, new BigRational(-3), Gam);
        Assert.True(atPlusThree.DTwoCoefficient == atMinusThree.DTwoCoefficient * new BigRational(2),
            $"stretch not seen: at eps=+3 {atPlusThree.DTwoCoefficient}, at eps=-3 " +
            $"{atMinusThree.DTwoCoefficient}; the ratio should be exactly 2.");
    }

    [Fact]
    public void AtZeroDefectEveryDifferenceVanishesAtEveryOrder()
    {
        // The page's control: at ε = 0 the encoding is invariant, so the light never reaches the motion.
        // The three −1 onsets alone would stay green under any mutation that reports no onset at all, so
        // the positive statement is asserted cell by cell instead, and it is not that everything
        // vanishes. δA does: the reflection-odd state keeps zero centre amplitude and the dephasing
        // charges it nothing. δB does NOT, because B wears the veil e^(−2γt) whatever the chain does, so
        // δB = (e^(−2γt) − 1)·A₀ ≠ 0. What vanishes in B is its CENTRE ROW, the only part of it that
        // reaches the decoded output, and that is why d_out has no onset here.
        foreach (var n in MissingPhaseOnsetWitness.PageSites)
        {
            int m = (n - 1) / 2, maxOrder = 2 * m + 2, c = m;
            var (a, b, reference) = MissingPhaseOnsetWitness.Coefficients(n, BigRational.Zero, Gam, maxOrder);
            var bMovesSomewhere = false;
            for (int k = 0; k <= maxOrder; k++)
                for (int p = 0; p < n; p++)
                    for (int q = 0; q < n; q++)
                    {
                        Assert.True(a[k][p, q] == reference[k][p, q],
                            $"N={n}, order {k}, cell ({p},{q}): δA = {a[k][p, q] - reference[k][p, q]} at ε = 0.");
                        if (p == c || q == c)
                            Assert.True(b[k][p, q] == reference[k][p, q],
                                $"N={n}, order {k}, centre cell ({p},{q}): δB = " +
                                $"{b[k][p, q] - reference[k][p, q]} at ε = 0.");
                        else if (b[k][p, q] != reference[k][p, q]) bMovesSomewhere = true;
                    }
            Assert.True(bMovesSomewhere,
                $"N={n}: δB vanished everywhere at ε = 0, so B is not wearing the veil and the centre-row " +
                "assertion above is passing on an empty B rather than on the physics.");

            var r = MissingPhaseOnsetWitness.Read(n, BigRational.Zero, Gam);
            Assert.Equal(-1, r.DOutOrder);
            Assert.Equal(-1, r.DTwoOrder);
            Assert.Equal(-1, r.LeakageOrder);
        }
    }

    [Fact]
    public void TheMatchedReferenceIsTheSameTrajectoryAtEveryRate()
    {
        // The reference arm must be the γ = 0 run of the SAME chain, so it may not depend on γ at all.
        // Asserting instead that the readings vanish at γ = 0 would assert nothing: there A, B and the
        // reference come out of the same arithmetic on the same seed and the difference is zero however
        // the generator is wired. This comparison can fail, and does the moment γ leaks into the
        // reference.
        foreach (var n in MissingPhaseOnsetWitness.PageSites)
        {
            var (_, _, atGamma) = MissingPhaseOnsetWitness.Coefficients(n, Eps, Gam, n);
            var (_, _, atOther) = MissingPhaseOnsetWitness.Coefficients(n, Eps, new BigRational(1, 7), n);
            for (int k = 0; k <= n; k++)
                for (int p = 0; p < n; p++)
                    for (int q = 0; q < n; q++)
                        Assert.True(atGamma[k][p, q] == atOther[k][p, q],
                            $"N={n}, order {k}, cell ({p},{q}): the reference moved with γ, " +
                            $"{atGamma[k][p, q]} against {atOther[k][p, q]}.");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(9)]
    public void TheBBlockWearsGammaFoldsVeilCoefficientByCoefficient(int n)
    {
        // B_γ(t) = e^(−2γt)·A_{−γ}(t) is GammaFold's price veil at Σ_l γ_l = γ. Read on Taylor
        // coefficients it says B_k = Σ_j (−2γ)^j/j! · A^{(−γ)}_{k−j}.
        // What this gates, exactly: the two dissipative factors are −γ(z_i z_j + 1) and, at −γ,
        // +γ(1 − z_i z_j), whose difference is −2γ in EVERY cell, so the veil is an algebraic consequence
        // of the construction rather than an independent confirmation of the physics. It is kept because
        // that scalar offset is the whole content of the identity and a wiring slip breaks it, not because
        // it corroborates §8; the independent route to §8 is the closed-form comparison above.
        int maxOrder = n + 2;
        var (_, b, _) = MissingPhaseOnsetWitness.Coefficients(n, Eps, Gam, maxOrder);
        var (aGain, _, _) = MissingPhaseOnsetWitness.Coefficients(n, Eps, -Gam, maxOrder);

        for (int k = 0; k <= maxOrder; k++)
        {
            for (int p = 0; p < n; p++)
                for (int q = 0; q < n; q++)
                {
                    var sum = GaussianRational.Zero;
                    var power = BigRational.One;
                    var fact = BigRational.One;
                    for (int j = 0; j <= k; j++)
                    {
                        if (j > 0) { power *= new BigRational(-2) * Gam; fact *= new BigRational(j); }
                        sum += new GaussianRational(power / fact, BigRational.Zero) * aGain[k - j][p, q];
                    }
                    Assert.True(sum == b[k][p, q],
                        $"N={n}, order {k}, cell ({p},{q}): veil gives {sum}, recurrence gives {b[k][p, q]}.");
                }
        }
    }

    [Fact]
    public void MovingTheLightOffTheCentreBreaksTheAgreement()
    {
        // The mutation control. Everything above would also pass if the readings were insensitive to
        // where the light sits; §8's coefficients are specifically the centre seat's. Off the centre the
        // decoded difference opens EARLIER, at t³ rather than t⁴, because the light now reaches the outer
        // block directly instead of having to walk m hops to the centre row. Both onsets are real
        // readings: a sentinel would make this test pass on an absence rather than on a difference, so
        // both are asserted to exist first.
        var offCentre = MissingPhaseOnsetWitness.ReadWithLightOnSeat(7, Eps, Gam, seat: 2);
        var centre = MissingPhaseOnsetWitness.Read(7, Eps, Gam);

        Assert.True(offCentre.DOutOrder >= 0, "the off-centre run reported no onset at all.");
        Assert.True(centre.DOutOrder >= 0, "the centre run reported no onset at all.");
        Assert.Equal(4, centre.DOutOrder);
        Assert.Equal(3, offCentre.DOutOrder);
    }

    [Fact]
    public void TheTwoBlockReductionIsThePhysicalLindbladian()
    {
        // §8.1's actual claim, and the one nothing else here gates: that ½·W[[A,B],[B,A]]W† is the real
        // 2^N × 2^N state under L(ρ) = −i[H,ρ] + γ(Z_c ρ Z_c − ρ). Comparing the reduced generator with a
        // doubled generator built from the same h and z would only check the split; this evolves the full
        // space in the computational basis, applying the action of the XX+YY bond terms and of Z_c ρ Z_c
        // per cell rather than building Pauli products, and asks whether the lifted blocks reproduce it
        // cell for cell.
        // The seed matters here: PhysicalSeed is built from the state ½|Φ⟩⟨Φ| and not through Embed, so
        // order 0 tests the embedding instead of restating it. With a seed taken through Embed, a wrong
        // embedding sits on both sides of the comparison and passes.
        const int n = 5, maxOrder = 6;
        var full = MissingPhaseOnsetWitness.FullLindbladCoefficients(n, Eps, Gam, maxOrder);
        var (a, b, _) = MissingPhaseOnsetWitness.Coefficients(n, Eps, Gam, maxOrder);
        int dim = 1 << n;

        for (int k = 0; k <= maxOrder; k++)
        {
            var lifted = MissingPhaseOnsetWitness.Embed(a[k], b[k], n);
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    Assert.True(lifted[i, j] == full[k][i, j],
                        $"order {k}, cell ({i},{j}): the reduction gives {lifted[i, j]}, the physical " +
                        $"Lindbladian gives {full[k][i, j]}.");
        }
    }

    [Theory]
    [InlineData("A and B swapped")]
    [InlineData("B negated")]
    public void TheReductionGateRejectsAWrongTwoBlockForm(string mutation)
    {
        // The controls for the gate above. Both mutations are applied to the lifted side only, exactly
        // as a wrong reduction would be, and both are caught. What a comparison like this has to avoid
        // is a seed DERIVED FROM WHAT IS UNDER TEST, and the reason is structural: H conserves the
        // excitation number, so the four index families the embedding writes into are each invariant,
        // and a sign put on B at the seed would simply travel down B's own family and meet the same sign
        // on the lifted side at every order. PhysicalSeed is built from the state's four amplitudes
        // instead, so a slip in Embed or InitialBlock shows at order 0 rather than cancelling. It
        // produces the SAME matrix as Embed(A₀, B₀), by a different route; that the two routes meet is
        // asserted separately below.
        const int n = 5, maxOrder = 3;
        var full = MissingPhaseOnsetWitness.FullLindbladCoefficients(n, Eps, Gam, maxOrder);
        var (a, b, _) = MissingPhaseOnsetWitness.Coefficients(n, Eps, Gam, maxOrder);
        int dim = 1 << n;

        var mismatches = 0;
        for (int k = 0; k <= maxOrder; k++)
        {
            var second = mutation == "A and B swapped" ? a[k] : Negate(b[k], n);
            var wrong = mutation == "A and B swapped"
                ? MissingPhaseOnsetWitness.Embed(b[k], second, n)
                : MissingPhaseOnsetWitness.Embed(a[k], second, n);
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    if (wrong[i, j] != full[k][i, j]) mismatches++;
        }
        Assert.True(mismatches > 0,
            $"the reduction gate accepted a wrong two-block form ({mutation}) at every cell of every " +
            "order, so it is not gating the two-block form of the state, only the generator.");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void TheStateBuiltSeedAndTheEmbeddedSeedAreTheSameMatrix(int n)
    {
        // Two routes to ρ(0): from the state ½|Φ⟩⟨Φ| with |Φ⟩ = |d₀⟩ + F|d₀⟩, and through the block
        // embedding of A₀ = B₀ = |d₀⟩⟨d₀|. They must agree cell for cell, and the reason to have both is
        // that the reduction gate needs a seed that does not come from the code it is testing. If this
        // ever parts, one of the two is wrong and the gate above is comparing a state with something
        // that is not the same state.
        var fromState = MissingPhaseOnsetWitness.PhysicalSeed(n);
        var fromBlocks = MissingPhaseOnsetWitness.Embed(
            MissingPhaseOnsetWitness.InitialBlock(n), MissingPhaseOnsetWitness.InitialBlock(n), n);
        int dim = 1 << n;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                Assert.True(fromState[i, j] == fromBlocks[i, j],
                    $"N={n}, cell ({i},{j}): from the state {fromState[i, j]}, from the blocks " +
                    $"{fromBlocks[i, j]}.");
    }

    private static GaussianRational[,] Negate(GaussianRational[,] x, int n)
    {
        var y = new GaussianRational[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                y[i, j] = -x[i, j];
        return y;
    }

    // ---- the readout maps themselves ---------------------------------------------------------------
    // The page's own system never exercises parts of these maps: the diagonal channel s of §8.2 is zero
    // at every onset (a diagonal change needs 2m hops plus a dephasing action), the trace of δA is zero
    // by conservation, and B_cc never reaches the decoded output. Those branches are gated here on
    // hand-built inputs instead, where they can be made nonzero.

    private static GaussianRational[,] Cells(int n, params (int Row, int Col, int Re, int Im)[] entries)
    {
        var m = new GaussianRational[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                m[i, j] = GaussianRational.Zero;
        foreach (var (row, col, re, im) in entries)
            m[row, col] = new GaussianRational(new BigRational(re), new BigRational(im));
        return m;
    }

    [Fact]
    public void TheTwoSiteDistanceReadsItsDiagonalChannel()
    {
        // s = δA_pp + δA_qq with no off-diagonal at all: the reading is |s|/2 + max(|s|/2, 0) = |s|.
        // A map that only ever looked at Re δA_pq would return 0 here. The input is TRACELESS, because
        // §8.2's corners are −s/2 only when the sites outside the pair carry −s between them; a
        // trace-4 input would exercise the written expression on a matrix it does not describe.
        var deltaA = Cells(5, (1, 1, 3, 0), (2, 2, 1, 0), (0, 0, -4, 0));
        Assert.True(MissingPhaseOnsetWitness.TwoSiteDistance(deltaA, 5) == new BigRational(4),
            "the diagonal channel of §8.2 is not being read.");
    }

    [Fact]
    public void TheTwoSiteDistanceTakesTheLargerOfItsTwoChannels()
    {
        // |v| = 1 against |s|/2 = 3: the formula is |s|/2 + max(|s|/2, |v|) = 6, not |s|/2 + |v| = 4.
        // Dropping the max is invisible on the page's own system, where s is always zero at the onset.
        var deltaA = Cells(5, (1, 1, 4, 0), (2, 2, 2, 0), (1, 2, 1, 0), (2, 1, 1, 0), (0, 0, -6, 0));
        Assert.True(MissingPhaseOnsetWitness.TwoSiteDistance(deltaA, 5) == new BigRational(6),
            "the two-site reading is not taking the larger of its two channels.");
    }

    [Fact]
    public void TheLeakageReadsBothHalvesOfItsProjector()
    {
        // Q = (I + R)/2 has an identity half and a reflection half. The physical coefficients are
        // traceless, so the identity half is never exercised by the readings themselves.
        Assert.True(MissingPhaseOnsetWitness.LeakageDifference(Cells(5, (0, 0, 2, 0)), 5) == BigRational.One,
            "the identity half of Q is not being read.");
        Assert.True(MissingPhaseOnsetWitness.LeakageDifference(Cells(5, (4, 0, 2, 0)), 5) == BigRational.One,
            "the reflection half of Q is not being read.");
    }

    [Fact]
    public void TheLeakageRefusesAnImaginaryPartInsteadOfDiscardingIt()
    {
        // Tr(Q·δA) is real because both factors are Hermitian. A δA that breaks that is a finding about
        // the recurrence, so the map must not quietly project it away.
        Assert.Throws<InvalidOperationException>(
            () => MissingPhaseOnsetWitness.LeakageDifference(Cells(5, (0, 0, 0, 1)), 5));
    }

    [Fact]
    public void TheDecodedOutputIgnoresExactlyTheEntriesThatDoNotReachIt()
    {
        // §8.3: ρ_out carries δA off the centre row and column, δA_cc, and the couplings δB_cj for
        // j ≠ c. The centre row of A and the entry B_cc do NOT reach it. Both directions are asserted,
        // so a map inspecting too much and a map inspecting too little each go red.
        const int n = 5, c = 2;
        var zero = Cells(n);
        Assert.True(
            MissingPhaseOnsetWitness.DecodedDifferenceVanishes(Cells(n, (c, 1, 5, 0), (1, c, 5, 0)), zero, n, c),
            "the centre row of A does not reach the decoded output, but the reading claims it does.");
        Assert.True(MissingPhaseOnsetWitness.DecodedDifferenceVanishes(zero, Cells(n, (c, c, 5, 0)), n, c),
            "B_cc does not reach the decoded output, but the reading claims it does.");
        Assert.False(MissingPhaseOnsetWitness.DecodedDifferenceVanishes(zero, Cells(n, (c, 0, 5, 0)), n, c),
            "the coupling B_c0 does reach the decoded output and must not be ignored.");
        Assert.False(MissingPhaseOnsetWitness.DecodedDifferenceVanishes(Cells(n, (0, 1, 5, 0)), zero, n, c),
            "the non-centre block of A does reach the decoded output and must not be ignored.");
    }

    [Fact]
    public void TheDecodedDistanceReadsTheCentreRowOfBAndNotItsColumn()
    {
        // δB is Hermitian on the page's system, so B[c,j] and B[j,c] agree there and the choice is
        // invisible. Here they are made to differ.
        const int n = 5, c = 2;
        var rowOnly = Cells(n, (c, 0, 3, 0));
        Assert.True(MissingPhaseOnsetWitness.DecodedDistanceSquared(Cells(n), rowOnly, n, c) == new BigRational(9),
            "the centre ROW of B is the coupling |f⟩⟨j|, and it is not being read.");
        var columnOnly = Cells(n, (0, c, 3, 0));
        Assert.True(MissingPhaseOnsetWitness.DecodedDistanceSquared(Cells(n), columnOnly, n, c).IsZero,
            "the reading is picking up B[j,c] where §8.3 puts B[c,j].");
        var withCentre = Cells(n, (c, 0, 3, 0), (c, c, 7, 0));
        Assert.True(MissingPhaseOnsetWitness.DecodedDistanceSquared(Cells(n), withCentre, n, c) == new BigRational(9),
            "B_cc is being counted into a norm it does not enter.");
    }

    [Fact]
    public void TheRankTwoConditionSeesTheCentrePopulation()
    {
        const int n = 5, c = 2;
        Assert.False(MissingPhaseOnsetWitness.DecodedIsRankTwo(Cells(n, (c, c, 1, 0)), c, n),
            "δA_cc is the |f⟩⟨f| population and takes the difference out of the rank-two case.");
        Assert.True(MissingPhaseOnsetWitness.DecodedIsRankTwo(Cells(n, (c, 1, 1, 0)), c, n),
            "the centre row of A never reaches the decoded output and must not disqualify the closed form.");
        Assert.Throws<InvalidOperationException>(
            () => MissingPhaseOnsetWitness.DecodedDistanceSquared(Cells(n, (c, c, 1, 0)), Cells(n), n, c));
    }

    [Fact]
    public void OffTheCentreSeatTheCoefficientIsReportedAsAbsentRatherThanGuessed()
    {
        // The one live path that returns a null coefficient: off the centre the A part does not vanish
        // at the onset, so no exact trace norm is returned, while the ONSET stays exact.
        var offCentre = MissingPhaseOnsetWitness.ReadWithLightOnSeat(7, Eps, Gam, seat: 2);
        Assert.Equal(3, offCentre.DOutOrder);
        Assert.False(offCentre.DOutCoefficientSquared.HasValue,
            "a coefficient was returned for a difference that is not the rank-two case of §8.3.");
    }

    [Fact]
    public void ThePageForbidsANonPositiveRateAndSoDoesTheWitness()
    {
        // §8's forms carry |ε|γ, so at γ ≤ 0 the oracle changes sign while the distance does not.
        Assert.Throws<ArgumentOutOfRangeException>(() => new MissingPhaseOnsetWitness(Eps, BigRational.Zero));
        Assert.Throws<ArgumentOutOfRangeException>(() => new MissingPhaseOnsetWitness(Eps, -Gam));
    }

    [Fact]
    public void TheGuardRejectsEvenAndTooShortChains()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => MissingPhaseOnsetWitness.Read(6, Eps, Gam));
        Assert.Throws<ArgumentOutOfRangeException>(() => MissingPhaseOnsetWitness.Read(3, Eps, Gam));
        Assert.Throws<ArgumentOutOfRangeException>(
            () => MissingPhaseOnsetWitness.Read(MissingPhaseOnsetWitness.MaxSites + 2, Eps, Gam));
    }
}
