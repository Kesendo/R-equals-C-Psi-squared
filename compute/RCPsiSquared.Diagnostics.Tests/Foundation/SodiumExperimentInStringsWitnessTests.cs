using System.Linq;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The experiment in strings, live (<c>experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md</c>, F167): the twelve
/// statements of the page's verifier read exactly over ℚ(ζ₂₄), most gates carrying a control inside the witness; the
/// readings the gates hold are pinned here as numbers.</summary>
public class SodiumExperimentInStringsWitnessTests
{
    private static readonly SodiumExperimentInStringsWitness W = new();

    [Fact]
    public void Field_RanksThroughTheRegularRepresentation_ReadTheControls() => Assert.True(W.FieldRankControl);

    [Fact]
    public void Labelling_IsACnotAndASwap_EightStringsMoveInKz()
    {
        Assert.True(W.LabellingReproducesIz);
        Assert.True(W.LabellingSignedImages);
        Assert.Equal(8, W.StringsMovedInKz);
        Assert.True(W.WrongAssignmentFails);
    }

    [Fact]
    public void CoherenceOrder_IsTwiceTheFirstBitsOrderPlusTheSecond_InSouzaBitsOnly()
    {
        Assert.True(W.OrderIsWeightedBitOrders);
        Assert.True(W.SwappedWeightsFail);
        Assert.True(W.PageBitsDoNotFixOrder);
    }

    [Fact]
    public void Generator_IsUnitalAndLinear_TheThermalSourceIsF120sForm()
    {
        Assert.True(W.GeneratorLinear);
        Assert.True(W.UnitalAndTracePreserving);
        Assert.True(W.AmplitudeDampingMovesIdentity);
        Assert.True(W.ThermalSourceHolds);
    }

    [Fact]
    public void InversionRecovery_FourFifthsAgainstTheSaturationsTwoFifths()
    {
        Assert.True(W.LongitudinalEigenoperators);
        Assert.Equal(new[] { BigRational.Zero, new BigRational(4, 5), new BigRational(-4, 5) }, W.InversionAmplitudes);
        Assert.Equal(new[] { BigRational.Zero, new BigRational(2, 5), new BigRational(-2, 5) }, W.SaturationAmplitudes);
        Assert.True(W.InversionRecoveryHolds);
    }

    [Fact]
    public void TurnAboutZ_FactorsInSouzaBits_InThePagesBitsAProductOnlyAtHalfTurns()
    {
        Assert.True(W.TurnFactorsInSouzaBits);
        Assert.True(W.TurnFormInPageBits);
        Assert.Equal(24, W.TurnRanksInPageBits.Count);
        Assert.Equal(new[] { 0, 6, 12, 18 }, W.TurnRanksInPageBits.Select((r, k) => (r, k)).Where(t => t.r == 1).Select(t => t.k));
        Assert.True(W.TurnRanksInPageBits.All(r => r is 1 or 2));
        Assert.True(W.TurnHolds);
    }

    [Fact]
    public void PhaseCycle_KeepsTheResidue_ThePathwayPhase_EightStepsOnTheSpin()
    {
        Assert.True(W.CycleKeepsResidue);
        Assert.True(W.PathwayPhase);
        Assert.True(W.EightStepsOnTheSpin);
    }

    [Fact]
    public void Pulse180_IsOneAntiDiagonalFactorPerBit_KeepsKz_SignedAtQuarterTurns()
    {
        Assert.True(W.Pulse180IsAProduct);
        Assert.True(W.Pulse180KeepsKz);
        Assert.True(W.Pulse180SignedAtQuarterTurns);
        Assert.True(W.Pulse90MovesKz);
    }

    [Fact]
    public void Echo_ReversesOffsetKeepsSplitting_NoPlainSandwichReversesBoth()
    {
        Assert.True(W.EchoReversesOffsetKeepsSplitting);
        Assert.True(W.NoPlainSandwichReversesBoth);
        Assert.True(W.TransposeReversesBothAndIsNoPulse);
    }

    [Fact]
    public void Pulse180_InThePagesBits_RanksOneOneTwoTwo()
    {
        Assert.Equal(new[] { 1, 1, 2, 2 }, W.Pulse180RanksInPageBits);
        Assert.True(W.Pulse180PageBitsHolds);
    }

    [Fact]
    public void Pulse90_Images_RankFourInSouzaBits_TwoFourFourFourInThePages()
    {
        Assert.True(W.Pulse90Images);
        Assert.Equal(new[] { 4, 4 }, W.Pulse90RanksInSouzaBits);
        Assert.Equal(new[] { 2, 4, 4, 4 }, W.Pulse90RanksInPageBits);
        Assert.True(W.Pulse90Holds);
    }

    [Fact]
    public void Receiver_ReadsOnlyRankOne_RelaxationShiftAndSplittingBringRankThreeIntoView()
    {
        Assert.True(W.ReceiverCells);
        Assert.True(W.ReceiverReadsRankOne);
        Assert.True(W.RelaxationCarriesRankThree);
        Assert.True(W.SplittingAndShiftInView);
        Assert.Equal(new[] { 2, 3 }, W.ReceiverRanks);
    }

    [Fact]
    public void Pools_TwoOneSiteJumpsOnALabel_TheFractionsAreTheSteadyState()
    {
        Assert.True(W.PoolTransferMatrix);
        Assert.True(W.PoolSteadyState);
    }

    [Fact]
    public void AllTwelveGates_Hold_AndTheNodesSayWhichReadingsStayTheVerifiers()
    {
        Assert.True(W.AllHold);
        Assert.Equal(13, W.Children.Count());
        var text = string.Join(" ", W.Children.Select(c => c.Summary).Prepend(W.Summary));
        Assert.Contains("the Theodolite's reading of the same is the verifier's", text);
        Assert.Contains("what LetterCube reads for a one-site jump (named, not run)", text);
        Assert.DoesNotContain("FAIL", text);
        Assert.DoesNotContain(" NOT ", text);                   // the failure branches' word ("CNOT" is not one)
        Assert.DoesNotContain("(NOT", text);
    }
}
