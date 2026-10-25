using System.Linq;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of F168 (<c>experiments/THE_CAT_UNDER_A_TURNED_FIELD.md</c>): the cat under
/// Z dephasing, a field at cos φ = 3/5 and a coupling commuting with X^⊗N and Σ Z_l, its centred odd
/// moments read exactly in the Pauli basis. Every comparison is exact (rational arithmetic).</summary>
public class CatTurnedFieldWitnessTests
{
    [Theory]
    [InlineData(1, false)]
    [InlineData(2, false)]
    [InlineData(3, false)]
    [InlineData(4, false)]
    [InlineData(3, true)]
    [InlineData(4, true)]
    public void OddMomentsBelowTwoNPlusOne_VanishExactly(int n, bool ring)
    {
        var w = new CatTurnedFieldWitness(n, ring);
        Assert.True(w.VanishingHolds, w.Summary);
    }

    [Theory]
    [InlineData(1, false)]
    [InlineData(2, false)]
    [InlineData(3, false)]
    [InlineData(4, false)]
    [InlineData(3, true)]
    [InlineData(4, true)]
    [InlineData(5, true)]
    public void FirstOddMoment_IsTheClaimsLaw_AndNotZero(int n, bool ring)
    {
        var w = new CatTurnedFieldWitness(n, ring);
        Assert.True(w.LawHolds, $"m_{2 * n + 1} = {w.OddMoments[2 * n + 1]} against the law {w.LawValue}");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void OneStrengthForEveryBit_WouldGiveAnotherNumber(int n) =>
        Assert.True(new CatTurnedFieldWitness(n).StrengthControlDiffers);

    [Theory]
    [InlineData(2)]
    [InlineData(4)]
    [InlineData(3)]
    public void AtTheZeroAngle_EveryComputedOddMomentVanishes(int n)
    {
        var w = new CatTurnedFieldWitness(n);
        Assert.True(w.ZeroAngleHolds, string.Join(", ", w.ZeroAngleOddMoments.Select(kv => $"m_{kv.Key} = {kv.Value}")));
    }

    [Fact]
    public void TheNextOddMoment_IsNotZero_SoTheZeroAngleCheckCanFail()
    {
        var w = new CatTurnedFieldWitness(3);
        Assert.False(w.OddMoments[9].IsZero);
    }

    [Fact]
    public void SinOfMultiple_ReadsTheCircle()
    {
        var c = new BigRational(3, 5);
        var s = new BigRational(4, 5);
        Assert.Equal(new BigRational(24, 25), CatTurnedFieldSelfEncounterClaim.SinOfMultiple(2, c, s));
        Assert.Equal(new BigRational(44, 125), CatTurnedFieldSelfEncounterClaim.SinOfMultiple(3, c, s));
        Assert.Throws<System.ArgumentException>(() =>
            CatTurnedFieldSelfEncounterClaim.SinOfMultiple(2, new BigRational(1, 2), new BigRational(1, 2)));
    }

    [Fact]
    public void Guard_RejectsOutOfRange()
    {
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CatTurnedFieldWitness(0));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CatTurnedFieldWitness(6));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CatTurnedFieldWitness(2, ring: true));
    }
}
