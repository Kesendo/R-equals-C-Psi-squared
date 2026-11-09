using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using T = RCPsiSquared.Core.Symmetry.ThermalThirdMomentClaim.Term;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness of <see cref="ThermalThirdMomentClaim"/> (F137, the third moment under cooling and
/// heating): the generator's Tr((2L + c)³) in exact Gaussian integers against the closed form. Every comparison is
/// exact.</summary>
public class ThermalThirdMomentWitnessTests
{
    private static readonly ThermalThirdMomentWitness W = new();

    [Fact]
    public void EveryRow_GeneratorEqualsTheClosedForm()
    {
        foreach (var x in W.Readings)
            Assert.True(x.Holds, $"{x.Row.Label}: generator {x.Generator}, law {x.Law}");
    }

    [Fact]
    public void RecordedValues_AreReproduced()
    {
        Assert.Equal(new BigInteger(8 * 288), W.Readings[0].Generator);
        Assert.Equal(new BigInteger(8 * -768), W.Readings[1].Generator);
    }

    [Fact]
    public void Rows_AreNotVacuous()
    {
        Assert.True(W.BothPartsNonzero >= 3, W.Summary);
        Assert.True(W.NotVacuous, W.Summary);
    }

    /// <summary>The rows discriminate: the law with the arrow's factor 2 dropped misses every row that carries both parts
    /// (a corollary of the generator check on rows with R ≠ 0, kept as the record that those rows exist).</summary>
    [Fact]
    public void WrongArrowFactor_MissesEveryRowWithBothParts()
    {
        foreach (var x in W.Readings.Where(r => !r.Field.IsZero && !r.Arrow.IsZero))
            Assert.NotEqual(x.Generator, 12 * BigInteger.Pow(4, x.Row.N) * (x.Field - x.Arrow));
    }

    /// <summary>The closed form drops the identity string, the implementation choice; that dropping it is right is the
    /// witness's constant row, generator against law.</summary>
    [Fact]
    public void Constant_IsDropped()
    {
        var cool = new long[] { 1, 2 };
        var heat = new long[] { 0, 1 };
        var with = ThermalThirdMomentClaim.CubedCentredTrace(new[] { new T("II", 7), new T("ZI", 1), new T("XX", 2) }, cool, heat);
        var without = ThermalThirdMomentClaim.CubedCentredTrace(new[] { new T("ZI", 1), new T("XX", 2) }, cool, heat);
        Assert.Equal(without, with);
    }

    [Fact]
    public void BadInput_Throws()
    {
        Assert.Throws<System.ArgumentException>(() => ThermalThirdMomentClaim.Parts(new[] { new T("XA", 1) }, new long[] { 1, 1 }, new long[] { 0, 0 }));
        Assert.Throws<System.ArgumentException>(() => ThermalThirdMomentClaim.Parts(new[] { new T("XI", 1) }, new long[] { 1, 1 }, new long[] { 0 }));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => ThermalThirdMomentClaim.Parts(new[] { new T("XI", 1) }, new long[] { -1, 1 }, new long[] { 0, 0 }));
    }
}
