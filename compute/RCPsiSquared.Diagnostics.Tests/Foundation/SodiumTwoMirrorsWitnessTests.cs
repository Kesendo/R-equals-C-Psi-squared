using System;
using System.Linq;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness for F167 (<c>experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md</c>): one
/// spin 3/2 under secular quadrupolar relaxation, the generator built from the spin matrices and
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
    public void Summary_AndFirstChild_ShowTheGate()
    {
        Assert.Contains("hold", Default.Summary);
        Assert.DoesNotContain("FAIL", Default.Summary);
        var first = Default.Children.First();
        Assert.Contains("gate", first.DisplayName);
        Assert.True(Default.Children.Count() >= 5);
    }

    [Fact]
    public void Ctor_RejectsOversizedInputs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SodiumTwoMirrorsWitness(SodiumTwoMirrorsWitness.MaxMagnitude + 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SodiumTwoMirrorsWitness(1, 1, 1, 1, -SodiumTwoMirrorsWitness.MaxMagnitude - 1));
    }
}
