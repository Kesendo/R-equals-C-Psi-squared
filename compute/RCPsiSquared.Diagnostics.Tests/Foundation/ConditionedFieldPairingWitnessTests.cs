using System.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of <see cref="ConditionedFieldPairingClaim"/> (F137, the conditioned field): the
/// triangularity on the integer pieces, the block identities proven on a grid one wider than their degree, the
/// pairing on the line, the law off it, the controls. Every comparison is exact.</summary>
public class ConditionedFieldPairingWitnessTests
{
    private static readonly ConditionedFieldPairingWitness W = new();

    [Fact]
    public void Generator_IsBlockTriangular_OnEveryPiece()
    {
        Assert.True(W.PiecesExact, W.Summary);
        Assert.Equal(0, W.Leaks);
        Assert.True(W.ArrowEntries > 0, "site 1's jump must show, else the check reads nothing");
    }

    [Fact]
    public void BlockIdentities_HoldOnTheWholeGrid()
    {
        Assert.Equal(3125, W.GridPoints);
        Assert.Equal(0, W.FieldFreeFailures);
        Assert.Equal(0, W.OneSidedFailures);
        Assert.Equal(0, W.FieldBlockFailures);
    }

    [Fact]
    public void Converse_HoldsOnTheWholeGrid()
    {
        Assert.True(W.ConverseProven, W.Summary);
    }

    [Fact]
    public void OnTheLine_EveryOddCoefficientAndOddSumVanishes()
    {
        Assert.True(W.PairsOnTheLine, W.Summary);
        foreach (var x in W.OnLineReadings)
        {
            Assert.Empty(x.NonzeroOddCoefficients);
            Assert.Empty(x.NonzeroOddSums);
        }
    }

    [Fact]
    public void OffTheLine_TheThirteenthIsTheLaw_AndNotZero()
    {
        Assert.True(W.ThirteenthIsTheLaw, W.Summary);
        Assert.True(W.BreaksOffTheLine, W.Summary);
    }

    [Fact]
    public void Controls_CancelTheThird_AndBreakAtTheFifth()
    {
        Assert.True(W.ControlsCancelThenBreak, W.Summary);
        foreach (var x in W.ControlReadings)
            Assert.True(x.CancelsThenBreaks, $"{x.Setting.Label}: third {x.Third}, law {x.Law}, fifth {x.Fifth}");
    }

    [Fact]
    public void ThirteenthCoefficient_ClosedForm()
    {
        Assert.Equal(BigInteger.Zero, ConditionedFieldPairingClaim.ThirteenthCoefficientOfTwoLPlusC(1, 1, 1, 1));
        Assert.Equal(new BigInteger(-64 * 4 * 2), ConditionedFieldPairingClaim.ThirteenthCoefficientOfTwoLPlusC(2, 3, 1, 1));
    }
}
