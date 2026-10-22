using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness for F167 (<c>experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md</c>): one
/// spin 3/2 under quadrupolar relaxation secular in the Larmor frequency, the generator built from the spin matrices and
/// read exactly over ℚ(i) (the √3 of I_± factored out of every T₂,ₘ with m ≠ 0, the √6 of T₂₀ taken
/// into the integer parameter q). Every gate is two-sided: the statement holds and its control
/// breaks; at every exact point the kernel count (F158) and the characteristic polynomial agree
/// with each other and with the claim's rule.</summary>
public class SodiumTwoMirrorsWitnessTests
{
    private static readonly SodiumTwoMirrorsWitness Default = new();

    [Fact]
    public void SqrtThree_FactorsOutOfEveryOffDiagonalTensor_ReadFromTheSpinMatrices()
    {
        Assert.True(Default.SqrtThreeFactorisationHolds);
        Assert.True(Default.AdjointRelationsHold);       // T₂,₁† = −T₂,₋₁, T₂,₂† = T₂,₋₂
        Assert.True(Default.SplittingParameterHolds);    // ω_Q·T₂₀ = q·diag(1, −1, −1, 1), q = ω_Q·√6/2
    }

    [Fact]
    public void FiveComponents_AreDiracMatrices_SquaringToThePagesConstants_AsTheFivePauliStrings()
    {
        Assert.True(Default.DiracAlgebraHolds);
        var expected = new[] { new BigRational(3, 2), 3, 3, 3, 3 }.Select(b => (GaussianRational)b).ToArray();
        Assert.Equal(expected, Default.SquareConstants);
        Assert.True(Default.PauliStringsMatch);
        Assert.Equal(new[] { "ZZ", "ZX", "ZY", "XI", "YI" }, Default.PauliLabels);
    }

    [Fact]
    public void DiracAlgebra_HoldsAtSpinThreeHalvesOnly_TheControlSpinsBreakIt()
    {
        var byTwiceSpin = Default.SpinReadings.ToDictionary(r => r.TwiceSpin);
        Assert.Equal(new[] { 2, 3, 4, 5 }, byTwiceSpin.Keys.OrderBy(k => k));
        Assert.True(byTwiceSpin[3].IsDirac);
        Assert.False(byTwiceSpin[2].IsDirac);
        Assert.False(byTwiceSpin[4].IsDirac);
        Assert.False(byTwiceSpin[5].IsDirac);
        // the similarity route at 3/2 reads the same square constants as the direct √3 route
        Assert.Equal(Default.SquareConstants, byTwiceSpin[3].Squares);
    }

    [Fact]
    public void OneMap_MirrorIdentityHoldsExactly_AtTheDefaultAndAtEveryPoint()
    {
        Assert.True(Default.MirrorIdentityHolds);
        Assert.True(Default.InverseHolds);
        Assert.All(Default.Points, p => Assert.True(p.MirrorIdentityHolds, $"S·L·S⁻¹ identity fails at {p}"));
    }

    [Fact]
    public void OneMap_ControlsBreakIt()
    {
        Assert.True(Default.TransposeDroppedBreaks);
        Assert.True(Default.UnitInPlaceOfGammaBreaks);
    }

    [Theory]
    [InlineData(0, 0, 0, 0, 0)]
    [InlineData(-3, 5, -2, 7, -4)]
    [InlineData(2, 0, 1, 3, 6)]
    public void OneMap_HoldsAtUserPoints_AnySign(long j0, long j1, long j2, long q, long omega)
    {
        var w = new SodiumTwoMirrorsWitness(j0, j1, j2, q, omega);
        Assert.True(w.MirrorIdentityHolds);
        Assert.True(w.SectorMirrorsHold);
    }

    [Fact]
    public void BrightProjector_CommutesWithLAndS_CarriesTheJ0Part_AndBothSectorsPair()
    {
        Assert.True(Default.ProjectorCommutesWithL);
        Assert.True(Default.ProjectorCommutesWithS);
        Assert.True(Default.J0PartIsMinusThreeJ0P);
        Assert.True(Default.DarkSectorPairs);
        Assert.True(Default.BrightSectorPairs);
        Assert.All(Default.Points, p => Assert.True(p.SectorMirrorsHold, $"sector mirrors fail at {p}"));
    }

    [Fact]
    public void GlobalCentre_BothRoutesAgreeWithEachOtherAndWithTheRule_AtEveryPoint()
    {
        Assert.All(Default.Points, p => Assert.True(p.Agrees, $"routes or rule disagree at {p}"));
        Assert.True(Default.TheoremHolds);
    }

    [Fact]
    public void GlobalCentre_PointsCoverBothOutcomes_AndBothRoutesWhereTheyApply()
    {
        Assert.Contains(Default.Points, p => p.Rule);
        Assert.Contains(Default.Points, p => !p.Rule);
        // route A only for non-negative densities, not all zero (at 2s = 0 its two ends are one matrix);
        // at least one negative-density point read by route B alone
        Assert.All(Default.Points, p => Assert.Equal(p.J0 >= 0 && p.J1 >= 0 && p.J2 >= 0 && p.TwoS != 0, p.RouteAApplies));
        Assert.Contains(Default.Points, p => p.J0 < 0 && !p.Rule);
        Assert.Contains(Default.Points, p => !p.RouteAApplies);
        // the required cases: J0 = 0 with both frequencies on (pairs); J0 > 0 with J1·J2 = 0 away from
        // ω_Q = ω = 0 (breaks); every density positive (breaks)
        Assert.Contains(Default.Points, p => p.J0 == 0 && p.Q != 0 && p.Omega != 0 && p.CharpolySymmetric);
        Assert.Contains(Default.Points, p => p.J0 > 0 && p.J1 * p.J2 == 0 && (p.Q != 0 || p.Omega != 0) && !p.CharpolySymmetric);
        Assert.Contains(Default.Points, p => p.J0 > 0 && p.J1 > 0 && p.J2 > 0 && !p.CharpolySymmetric);
    }

    [Fact]
    public void CountingPoint_ReadsThePagesKernels_FourAgainstTwo()
    {
        // J1 = J2 = 0, J0 = 1, ω_Q = √6 (q = 3), ω = √6·ω_Q = 6 = 2q: the satellite sits at −2s, but 0 is fourfold and −3 twofold
        var p = Default.Points.Single(x => x.J0 == 1 && x.J1 == 0 && x.J2 == 0 && x.Q == 3 && x.Omega == 6);
        Assert.Equal(4, p.NearKernel);
        Assert.Equal(2, p.FarKernel);
        Assert.False(p.CharpolySymmetric);
        Assert.False(p.Rule);
    }

    [Fact]
    public void AllDensitiesPositive_ReadOneAgainstNone()
    {
        var p = Default.Points.First(x => x.J0 > 0 && x.J1 > 0 && x.J2 > 0);
        Assert.Equal(1, p.NearKernel);
        Assert.Equal(0, p.FarKernel);
    }

    [Fact]
    public void Rule_IsTheClaimsRule()
    {
        Assert.All(Default.Points, p => Assert.Equal(
            Spin32QuadrupolarMirrorsClaim.HasGlobalCentre(p.J0, p.J1, p.J2, p.Q, p.Omega), p.Rule));
        Assert.True(Spin32QuadrupolarMirrorsClaim.HasGlobalCentre(0, 4, 5, 2, 3));
        Assert.True(Spin32QuadrupolarMirrorsClaim.HasGlobalCentre(2, 0, 5, 0, 0));
        Assert.False(Spin32QuadrupolarMirrorsClaim.HasGlobalCentre(2, 0, 5, 0, 1));
        Assert.False(Spin32QuadrupolarMirrorsClaim.HasGlobalCentre(2, 1, 5, 0, 0));
        Assert.Equal(15, Spin32QuadrupolarMirrorsClaim.DarkCentre(2, 3));
        Assert.Equal(18, Spin32QuadrupolarMirrorsClaim.BrightCentre(1, 2, 3));
        Assert.Equal(33, Spin32QuadrupolarMirrorsClaim.TwiceGlobalCentre(1, 2, 3));
    }

    [Fact]
    public void UserPoint_IsReadByBothRoutes()
    {
        var w = new SodiumTwoMirrorsWitness(0, 2, 1, 4, 5);
        var user = w.UserPoint;
        Assert.True(user.Rule);
        Assert.True(user.Agrees);
        var neg = new SodiumTwoMirrorsWitness(1, -2, 3, 0, 0);
        Assert.False(neg.UserPoint.RouteAApplies);
        Assert.False(neg.UserPoint.CharpolySymmetric);
        Assert.True(neg.UserPoint.Agrees);
    }

    [Fact]
    public void Handshake_TheTextbookRelaxation_IsStringsAndCounting()
    {
        Assert.True(Default.LongitudinalSplitHolds);            // I_z = ZI + ½·IZ, ZI at 6J₂, IZ at 6J₁
        Assert.Equal(new GaussianRational[] { new BigRational(4, 5), new BigRational(1, 5) }, Default.LongitudinalWeights);
        Assert.True(Default.TransverseSplitHolds);              // at ω_Q = 0: CT at c₁, satellites at 3(J₀ + J₁)
        Assert.Equal(new GaussianRational[] { new BigRational(3, 5), new BigRational(2, 5) }, Default.TransverseWeights);
        Assert.True(Default.SplittingDetunesSatellites);        // the controls: a splitting, and swapped rate assignments
        Assert.True(Default.SwappedAssignmentBreaks);
        Assert.True(Default.SplittingFeedsPartner);             // −i[Γ̂₀, satellites] = −2i·partner
        Assert.Equal(new GaussianRational[] { 3, 4, 3 }, Default.TripletIntensities);
        Assert.True(Default.DarkLinesBlindToJ0AndSplitting);
        Assert.True(Default.TripleQuantumOnDarkCentre);
        Assert.True(Default.PopulationsPairAboutDarkCentre);
        Assert.True(Default.ExtremeNarrowingOneRate);
        Assert.True(Default.SatellitePartnerOnBrightMirror);      // −S(satellites)† = T₂,₁ up to a factor, at 2c₂ − 3(J₀ + J₁)
        Assert.True(Default.HandshakeHolds);
        Assert.Contains("the handshake with NMR holds", Default.Summary);
        Assert.Contains(Default.Children, c => c.DisplayName == "the handshake with NMR");
    }

    [Theory]
    [InlineData(0, 5, 1, 3, 7)]
    [InlineData(4, 1, 9, 0, 0)]
    [InlineData(-2, 3, 5, 2, -1)]
    public void Handshake_HoldsAtOtherPoints(long j0, long j1, long j2, long q, long omega)
    {
        var w = new SodiumTwoMirrorsWitness(j0, j1, j2, q, omega);
        Assert.True(w.HandshakeHolds);
    }

    [Fact]
    public void DynamicShift_KeepsBothMirrors_MovesEachPairAsAWhole_AndKeepsTheJordanBlock()
    {
        Assert.True(Default.ShiftClosedForm);          // built from the spin matrices: −3K₂·ZI − 3K₁·IZ
        Assert.True(Default.ShiftKeepsMirror);
        Assert.True(Default.ShiftMovesPairsWhole);      // satellites by 6i·K₁, double quantum by 6i·K₂, CT by 6i(K₂ − K₁)
        Assert.True(Default.ShiftKeepsJordan);          // ranks 15, 14, 14 at λ* + 6i·K₁, rank 16 one step off in q
        Assert.True(Default.DynamicShiftHolds);
    }

    [Fact]
    public void DynamicShift_GlobalCentre_FollowsTheShiftedRule_BothRoutesAgainstLiterals()
    {
        // J₁ = 0 at ω = 6K₁ pairs, at ω = 0 not; J₂ = 0 at ω = 3K₂ pairs, at ω = 0 not; the same with both K on, and at
        // an ω a rule mixing the K would accept, not; with a splitting, not; J₀ = 0 pairs; all densities on, not
        var expected = new[] { true, false, true, false, true, false, true, false, false, true, false };
        Assert.Equal(expected, Default.ShiftCentreReadings.Select(r => r.CharpolySymmetric).ToArray());
        Assert.Equal(expected, Default.ShiftCentreReadings.Select(r => r.Rule).ToArray());
        Assert.Equal(expected, Default.ShiftCentreReadings.Select(r => r.NearKernel == r.FarKernel).ToArray());
    }

    [Theory]
    [InlineData(1, 0, 3, 0, 0)]
    [InlineData(1, 2, 0, 0, 0)]
    [InlineData(1, 2, 3, 0, 0)]
    [InlineData(0, 2, 3, 1, 1)]
    [InlineData(1, 0, 3, 1, 0)]
    [InlineData(1, 0, 0, 0, 2)]
    public void DynamicShift_RuleAtZeroShift_IsTheRuleWithoutIt(long j0, long j1, long j2, long q, long omega) =>
        Assert.Equal(Spin32QuadrupolarMirrorsClaim.HasGlobalCentre(j0, j1, j2, q, omega),
                     Spin32QuadrupolarMirrorsClaim.HasGlobalCentreWithShift(j0, j1, j2, q, omega, 0, 0));

    [Fact]
    public void ExceptionalPoint_BrightPairs_AreClosedBlocks_WithTheSumRuleAndTheirDiscriminants()
    {
        Assert.True(Default.BrightPairsClosed);              // the control set {satellite, central transition} is not
        Assert.True(Default.BrightPairTracesHold);           // tr = −2c₂ − 2ip·ω: a pair's two rates sum to 2c₂
        Assert.True(Default.BrightPairDiscriminantsHold);    // 4(9J₂² − 4q²), 4(9J₁² − 4q²); the swapped assignment fails
    }

    [Fact]
    public void ExceptionalPoint_SatellitesAndDoubleQuantum_AreJordanPairsOfTheWholeGenerator()
    {
        foreach (var e in new[] { Default.SatelliteEp, Default.DoubleQuantumEp })
        {
            Assert.True(e.PartNonZero);
            Assert.True(e.PartSquaresToZero);
            Assert.Equal((15, 14, 14), (e.Rank1, e.Rank2, e.Rank3));
            Assert.Equal((1, 2), (e.Geometric, e.Algebraic));
            Assert.Equal((16, 16), (e.RankStepDown, e.RankStepUp));
            Assert.True(e.IsJordanPair);
        }
        Assert.Equal(new GaussianRational(-21, -1), Default.SatelliteEp.Lambda);      // −c₂ − iω at (1, 4, 2, 3, 1)
        Assert.Equal(new GaussianRational(-21, -2), Default.DoubleQuantumEp.Lambda);  // −c₂ − 2iω at (1, 2, 4, 3, 1)
    }

    [Fact]
    public void ExceptionalPoint_TheJordanReading_FailsOffTheExceptionalPoint_AndOnTheScalarControl()
    {
        var below = Default.ReadEpAt((1, 4, 2, 2, 1), 1);           // one step below the satellites' EP
        Assert.Equal(16, below.Rank1);                                // λ* is no eigenvalue at all
        Assert.False(below.IsJordanPair);
        Assert.False(Default.ReadEpAt((1, 4, 2, 3, 1), 2).IsJordanPair);   // the double quantum at the satellites' EP
        var scalar = Default.ReadEpAt(SodiumTwoMirrorsWitness.DiabolicPoint, 1);
        Assert.False(scalar.PartNonZero);                             // the scalar block: N = 0
        Assert.Equal((2, 2), (scalar.Geometric, scalar.Algebraic));
        Assert.False(scalar.IsJordanPair);
        Assert.True(Default.DiabolicControlSemisimple);
    }

    [Fact]
    public void ExceptionalPoint_TheCoalescedVector_IsSelfOrthogonal()
    {
        Assert.True(Default.EpVectorSelfOrthogonal);
        Assert.True(Default.GeneratorComplexSymmetric);
    }

    [Fact]
    public void ExceptionalPoint_TheSatelliteSignal_InitialRate_BelowAndAboveTheExceptionalPoint()
    {
        Assert.True(Default.InitialSlopeSplittingFree);              // R₂′ = 3(J₀ + J₁) at every splitting, one cell alone not
        Assert.Equal(9, Default.BelowEpOmega);                       // pin of the closed form √(9·25 − 4·36) at (1, 3, 5, 6, 1)
        Assert.Equal(16, Default.AboveEpOmega);                      // pin of √(4·100 − 9·16) at (1, 3, 4, 10, 1); the ranks test it on L
        Assert.True(Default.BelowEpTwoRates);
        Assert.Equal(new GaussianRational[] { new BigRational(4, 3), new BigRational(-1, 3) }, Default.BelowEpWeights);
        Assert.True(Default.AboveEpFrequenciesPulled);               // modes at ±16, neither at ±20; D_ω² + R₁(4/5)² = D²
        Assert.True(Default.EpAtTheLongitudinalRates);               // 4q = R₁(4/5) for the satellites, R₁(1/5) for the double quantum
        Assert.True(Default.ExceptionalPointHolds);
        Assert.Contains("the bright pairs' exceptional point holds", Default.Summary);
        Assert.Contains(Default.Children, c => c.DisplayName == "the bright pairs' exceptional point");
    }

    [Theory]
    [InlineData(1, 4, 2, 3, 1, "on its exceptional point")]
    [InlineData(1, 2, 3, 1, 1, "below its exceptional point")]
    [InlineData(1, 3, 4, 10, 1, "past its exceptional point")]
    [InlineData(1, 3, 0, 2, 1, "uncoupled")]
    [InlineData(1, 4, 0, 0, 1, "a scalar")]
    public void ExceptionalPoint_TheSatelliteRegime_IsReadOffTheGenerator(long j0, long j1, long j2, long q, long omega, string expected)
    {
        var w = new SodiumTwoMirrorsWitness(j0, j1, j2, q, omega);
        Assert.StartsWith(expected, w.SatelliteRegime);
        Assert.True(w.BrightPairsClosed && w.BrightPairTracesHold && w.BrightPairDiscriminantsHold);
    }

    [Fact]
    public void ExceptionalPoint_TheRepositoryEpInstrument_ReadsDefectiveOnBothPairs_DiabolicAtTheControl()
    {
        // a third route, numerical and artifact-free (EpCharacter: a Riesz projector, the departure from normality,
        // geometric against algebraic multiplicity) on the same exact generator in doubles: the exact gate
        // certifies, this classifies
        static double ToDouble(BigRational r) => (double)r.Numerator / (double)r.Denominator;
        static Matrix<Complex> ToComplex(GaussianRational[,] m) =>
            Matrix<Complex>.Build.Dense(16, 16, (i, j) => new Complex(ToDouble(m[i, j].Re), ToDouble(m[i, j].Im)));
        // at two isolating radii (the nearest other eigenvalue is at least 2 away at every point read)
        foreach (double radius in new[] { 0.25, 0.5 })
        {
            EpCharacter.Reading At((long J0, long J1, long J2, long Q, long Omega) p, Complex centre) =>
                EpCharacter.Characterize(ToComplex(Default.GeneratorMatrix(p.J0, p.J1, p.J2, p.Q, p.Omega)), centre, radius);

            var sat = At(SodiumTwoMirrorsWitness.SatelliteEpPoint, new Complex(-21, -1));
            Assert.Equal(EpCharacter.EpKind.Defective, sat.Kind);
            Assert.Equal((2, 1), (sat.Algebraic, sat.Geometric));
            var dq = At(SodiumTwoMirrorsWitness.DoubleQuantumEpPoint, new Complex(-21, -2));
            Assert.Equal(EpCharacter.EpKind.Defective, dq.Kind);
            Assert.Equal((2, 1), (dq.Algebraic, dq.Geometric));
            var control = At(SodiumTwoMirrorsWitness.DiabolicPoint, new Complex(-15, -1));
            Assert.Equal(EpCharacter.EpKind.Diabolic, control.Kind);
            Assert.Equal((2, 2), (control.Algebraic, control.Geometric));
        }
    }

    [Fact]
    public void Summary_AndFirstChild_ShowTheGate()
    {
        Assert.Contains("hold", Default.Summary);
        Assert.DoesNotContain("FAIL", Default.Summary);
        var first = Default.Children.First();
        Assert.Contains("gate", first.DisplayName);
        Assert.True(Default.Children.Count() >= 6);
    }

    [Fact]
    public void Ctor_RejectsOversizedInputs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SodiumTwoMirrorsWitness(SodiumTwoMirrorsWitness.MaxMagnitude + 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SodiumTwoMirrorsWitness(1, 1, 1, 1, -SodiumTwoMirrorsWitness.MaxMagnitude - 1));
    }
}
