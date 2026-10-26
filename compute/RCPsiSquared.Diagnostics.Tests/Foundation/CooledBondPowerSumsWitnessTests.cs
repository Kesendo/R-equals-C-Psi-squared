using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of <see cref="CooledBondPowerSumsClaim"/>
/// (<c>experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md</c>): a chain and a ring of one class bond under loss
/// alone, the odd power sums of 2L + c read exactly. Every comparison is exact.</summary>
public class CooledBondPowerSumsWitnessTests
{
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void ChainLowOddSums_VanishExactly(int n)
    {
        var w = new CooledBondPowerSumsWitness(n);
        Assert.True(w.ChainLowVanish, w.Summary);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void Seventh_IsTheClosedForm_AndNotZero(int n)
    {
        var w = new CooledBondPowerSumsWitness(n);
        Assert.True(w.SeventhIsLaw, $"p₇ = {w.ChainSums[7]} against the law {w.Law}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void EqualEndProducts_AndTheControl_AreSilent(int n)
    {
        var w = new CooledBondPowerSumsWitness(n);
        Assert.True(w.EqualEndsSilent, w.Summary);
        Assert.True(w.ControlSilent, w.Summary);
        Assert.Equal(n >= 4, w.EqualEndsMovedByReversal);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void Ring_IsSilentThroughTheSeventh_AndTheNinthBreaks(int n)
    {
        var w = new CooledBondPowerSumsWitness(n);
        Assert.True(w.RingSilentThroughSeventh, w.Summary);
        Assert.True(w.RingNinthBreaks, w.Summary);
    }

    /// <summary>The generator itself at the reversed profile: its seventh is −p₇ (the page's p_k(Pγ) = −p_k(γ)),
    /// and the closed form at the reversed profile agrees.</summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void Reversal_FlipsTheSeventh_OnTheGenerator(int n)
    {
        var w = new CooledBondPowerSumsWitness(n);
        Assert.True(w.ReversalFlipsSeventh, w.Summary);
        var reversed = w.Breaking.Reverse().ToArray();
        Assert.Equal(w.ReversedSeventh, CooledBondPowerSumsClaim.SeventhPowerSum(CooledBondPowerSumsWitness.MainBond, reversed));
    }

    [Fact]
    public void OutOfRange_Throws()
    {
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CooledBondPowerSumsWitness(2));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CooledBondPowerSumsWitness(6));
    }
}
