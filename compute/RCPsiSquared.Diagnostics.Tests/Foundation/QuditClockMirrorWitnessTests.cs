using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live lab of F172: every mirror identity holds exactly in ℤ[ζ_m] on every matrix unit, every control
/// breaks, and the hidden qubit reads −1 at even d while no shift does at odd d.</summary>
public class QuditClockMirrorWitnessTests
{
    [Theory]
    [InlineData(4, 1)]
    [InlineData(4, 2)]
    [InlineData(4, 3)]
    [InlineData(6, 2)]
    [InlineData(8, 2)]
    public void IdentitiesHold_AndControlsBreak(int d, int n)
    {
        var w = new QuditClockMirrorWitness(d, n, 1);
        Assert.True(w.HalfShiftAnticommutes);
        Assert.True(w.HiddenQubitFactorizes);
        Assert.True(w.HamiltoniansHermitian);
        Assert.All(w.OddDimensionReadings, o => Assert.True(o.NoShiftAnticommutes));
        Assert.All(w.MustHold, x => Assert.True(x.Holds, $"{x.Name} at d={d}, N={n}"));
        Assert.All(w.MustBreak, x => Assert.False(x.Holds, $"{x.Name} at d={d}, N={n}"));
        Assert.All(w.AlsoReached, x => Assert.True(x.Holds, $"{x.Name} at d={d}, N={n}"));
        Assert.True(w.RoutesConditionsOnOwnH);
        if (n >= 2) Assert.Empty(w.AlsoReached);   // on two or more sites each route's H is reached by its route alone
        Assert.True(w.FarElementHamiltonianMissesU);
        Assert.True(w.AllHold);
    }

    [Fact]
    public void ARateOtherThanOneKeepsTheIdentities()
    {
        var w = new QuditClockMirrorWitness(4, 2, 3);
        Assert.True(w.AllHold);
    }

    [Fact]
    public void EveryUnitIsRead()
    {
        var w = new QuditClockMirrorWitness(4, 2, 1);
        Assert.All(w.MustHold, x => Assert.Equal(16 * 16, x.UnitsRead));
        Assert.All(w.MustBreak, x => Assert.True(x.FirstFailure.HasValue && x.UnitsRead <= 16 * 16));
    }

    [Fact]
    public void OnOneQuquart_TheImaginaryHamiltoniansAreAlsoReversedByY()
    {
        // d = 4, N = 1: an antisymmetric real matrix commuting with U = X^2 joins only levels of opposite parity, so i times
        // it is reversed by Y = Z^2; route 3 reaches route 4's H there, and the witness reads it holding
        var w = new QuditClockMirrorWitness(4, 1, 1);
        Assert.Contains(w.LoneConjugationControls, c => c.ConditionHolds && c.Reading.Holds && c.Reading.Name.StartsWith("route 3"));
        Assert.True(w.AllHold);
    }

    [Theory]
    [InlineData(3, 1)]
    [InlineData(5, 1)]
    [InlineData(1, 1)]
    [InlineData(2, 2)]
    [InlineData(4, 0)]
    [InlineData(4, 4)]
    public void RejectsOutOfScope(int d, int n) =>
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new QuditClockMirrorWitness(d, n, 1));

    [Fact]
    public void Summary_SaysReadLive()
    {
        var w = new QuditClockMirrorWitness(4, 2, 1);
        Assert.Contains("F172 read live", w.Summary);
        Assert.Equal(4, w.Children.Count());
    }
}
