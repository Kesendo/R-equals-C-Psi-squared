using System;
using System.Linq;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness for the chain part of F167 (<c>experiments/THE_SODIUM_CHAIN.md</c>) at
/// two nuclei: the 256 × 256 generator exact over ℚ(i), the coupling read from the spin-3/2
/// matrices, the √3 of the full dipolar coupling carried by a diagonal similarity that makes the
/// generator rational. Every gate is two-sided: the statement holds and its control breaks. The
/// sharp-subspace dimensions are readings at the witness's parameters, pinned against the page's
/// bounds over two primes.</summary>
public class SodiumChainWitnessTests
{
    private static readonly SodiumChainWitness Default = new();

    [Fact]
    public void Dictionary_SpinOperatorsAreTheStringForms_ReadFromTheSpinMatrices()
    {
        Assert.True(Default.IzMatches);      // I_z = ZI + ½·IZ
        Assert.True(Default.IxMatches);      // I_x = (√3/2)·IX + ½(XX + YY)
        Assert.True(Default.IyMatches);      // I_y = (√3/2)·IY + ½(YX − XY)
        Assert.True(Default.T20Matches);     // T₂₀ = (3/√6)·ZZ
        Assert.True(Default.SqrtThreePartNeeded);
    }

    [Fact]
    public void DipolarCoupling_22Strings_TheEightSqrtThreeOnesExactlyThoseAnticommutingWithQ()
    {
        Assert.Equal(22, Default.DipolarStringCount);
        Assert.Equal(new[] { "IXXX", "IXYY", "IYXY", "IYYX", "XXIX", "XYIY", "YXIY", "YYIX" }, Default.SqrtThreeStrings);
        Assert.True(Default.PartsDisjoint);
        Assert.True(Default.SqrtThreeAreTheQAnticommuting);  // and NOT those anticommuting with ZZ₁ alone
        Assert.Equal(new[] { "IXIX", "IYIY" }, Default.SatelliteSatelliteStrings);
        Assert.True(Default.SatelliteSatelliteAnticommuteWithEach);
        Assert.True(Default.CouplingsReal);
    }

    [Fact]
    public void Truncations_AreTheShellProjections_AndDiffer()
    {
        Assert.True(Default.UnequalTruncationIsShellProjection);
        Assert.True(Default.EqualTruncationIsShellProjection);
        Assert.True(Default.TruncationsDiffer);
        Assert.True(Default.DictionaryHolds);
    }

    [Fact]
    public void Lemma_HoldsEntryByEntry_ForIsingAndBothTruncations()
    {
        Assert.True(Default.LemmaIsingPerNucleus);
        Assert.True(Default.LemmaIsingIdentical);
        Assert.True(Default.LemmaUnequal);
        Assert.True(Default.LemmaEqual);
        Assert.True(Default.LemmaHolds);
    }

    [Fact]
    public void Lemma_TheSqrtThreeStringsBreakIt()
    {
        Assert.True(Default.SqrtThreeStringsBreakLemma);
    }

    [Fact]
    public void Ladder_IsingAndUnequalKeepBrightness_TheEqualTruncationMovesIt()
    {
        Assert.True(Default.IsingKeepsBrightness);
        Assert.True(Default.UnequalKeepsBrightness);
        Assert.False(Default.EqualKeepsBrightness);
    }

    [Fact]
    public void Ladder_EveryBrightSetBlockPairsAboutItsOwnRung_AndAboutNoNeighbour()
    {
        var perNucleus = new[] { (Default.J0, Default.J1, Default.J2), SodiumChainWitness.SecondNucleusRates };
        var identical = new[] { (Default.J0, Default.J1, Default.J2), (Default.J0, Default.J1, Default.J2) };
        foreach (var (blocks, rates) in new[] { (Default.IsingBlocks, perNucleus), (Default.UnequalBlocks, identical) })
        {
            Assert.Equal(4, blocks.Count);
            Assert.All(blocks, b =>
            {
                Assert.Equal(64, b.Size);
                Assert.Equal(SodiumChainMirrorLadderClaim.TwiceCentre(rates, b.BrightMask), b.TwoS);
                Assert.True(b.Pairs, $"block {b.Label} does not pair about its rung");
                Assert.All(b.Neighbours, n => Assert.False(n.Pairs, $"block {b.Label} pairs about the neighbour 2s = {n.TwoS}"));
            });
        }
        // rates per nucleus: four distinct rungs, up to 2ⁿ
        Assert.Equal(4, Default.IsingBlocks.Select(b => b.TwoS).Distinct().Count());
        Assert.True(Default.LadderHolds);
    }

    [Fact]
    public void FullCoupling_FarSpanAtJ0Zero_IsQAlone_AndEmptyWithTheJ0Jumps()
    {
        Assert.Equal(new[] { "ZZZZ" }, Default.FarSpanAtJ0Zero);
        Assert.Equal(0, Default.FarSpanWithJ0Jumps);
        Assert.True(Default.RestCommutesWithQ);
        Assert.True(Default.SqrtThreeAnticommutesWithQ);
    }

    [Fact]
    public void FullCoupling_GlobalMirrorBreaks_ByBothRoutes_AtSqrtThreeAndAtTheStandIn()
    {
        Assert.True(Default.SimilarityIsRational);
        Assert.Equal(72, Default.TwoSigma);
        Assert.Equal(1, Default.NearKernel);
        Assert.Equal(0, Default.FarKernel);
        Assert.Equal(1, Default.NearKernelAtC);
        Assert.Equal(0, Default.FarKernelAtC);
        Assert.False(Default.CharpolyPairsAboutSigma);
    }

    [Fact]
    public void FullCoupling_Control_TheEqualTruncationAtJ0Zero_PairsByBothRoutes()
    {
        Assert.Equal(Default.ControlNearKernel, Default.ControlFarKernel);
        Assert.True(Default.ControlCharpolyPairs);
        Assert.True(Default.GlobalBreakHolds);
    }

    [Fact]
    public void Readings_SharpSubspaces_MeetThePagesBounds_AndPairAsThePageReads()
    {
        Assert.Equal(new[] { 24, 128, 24 }, Default.EqualSharp.Select(h => h.Dimension));
        Assert.All(Default.EqualSharp, h => Assert.True(h.Pairs));
        Assert.Equal(new[] { 10, 12, 1 }, Default.FullSharp.Select(h => h.Dimension));
        Assert.Equal(new bool?[] { false, true, true }, Default.FullSharp.Select(h => h.Pairs));   // K₀ fails, proved on the page
        Assert.Equal(4, Default.SqrtThreeRank);
        Assert.Equal(new[] { 5, 12, 1 }, Default.FullKernelHull.Select(h => h.Dimension));
        Assert.All(Default.FullKernelHull, h => Assert.True(h.Pairs));                             // the proposition, every rung
        Assert.Equal(new[] { 64, 128, 64 }, Default.EqualSharp.Select(h => h.CellCount));
        Assert.Equal(new long[] { 72, 108, 144 }, Default.EqualSharp.Select(h => h.TwoS));
        Assert.True(Default.K0BeyondU0MatchesPage);
        Assert.True(Default.HullsInvariant);
    }

    [Fact]
    public void AnotherPoint_TheGatesHold_TheStandInBreaksForAnyCoefficient_TheReadingsRepeat()
    {
        // the page's second rate point (8, 6, 2); another coupling sign and size; a negative stand-in
        var w = new SodiumChainWitness(8, 6, 2, -3, -2);
        Assert.True(w.DictionaryHolds);
        Assert.True(w.LemmaHolds);
        Assert.True(w.SqrtThreeStringsBreakLemma);
        Assert.True(w.LadderHolds);
        Assert.True(w.GlobalBreakHolds);
        Assert.Equal(w.TwoSigma, SodiumChainMirrorLadderClaim.TwiceRung(2, 0, 0, 6, 2));
        Assert.Equal(new[] { 24, 128, 24 }, w.EqualSharp.Select(h => h.Dimension));
        Assert.Null(w.K0BeyondU0MatchesPage);   // the page states the quotient at its own point only
        Assert.True(w.HullsInvariant);
    }

    [Fact]
    public void Summary_AndFirstChild_ShowTheGate()
    {
        Assert.Contains("holds", Default.Summary);
        Assert.DoesNotContain("FAIL", Default.Summary);
        Assert.DoesNotContain("NOT", Default.Summary);
        var first = Default.Children.First();
        Assert.Contains("gate", first.DisplayName);
        Assert.Equal(5, Default.Children.Count());
    }

    [Fact]
    public void Ctor_RejectsInputsOutsideItsScope()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SodiumChainWitness(SodiumChainWitness.MaxMagnitude + 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SodiumChainWitness(0));          // J₀ = 0 is gate 4's own point
        Assert.Throws<ArgumentOutOfRangeException>(() => new SodiumChainWitness(6, 0, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SodiumChainWitness(6, 4, 2, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SodiumChainWitness(6, 4, 2, 4, 0));
    }
}
