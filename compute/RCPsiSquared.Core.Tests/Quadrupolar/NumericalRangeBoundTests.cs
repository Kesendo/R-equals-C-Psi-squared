using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using Xunit;

namespace RCPsiSquared.Core.Tests.Quadrupolar;

/// <summary>The bound decides (no Hermitian H gives a centre) exactly where D + 2s·1 is positive
/// definite. At 3/2 that is every density positive (F167's "no" side once more); at J₀ = 0 the string
/// ZZ anticommutes with every jump and the bound is silent, as at J₁·J₂ = 0. At extreme narrowing the
/// closed form gives r_max/2s = 3/4, 4/5, 11/14 at 2I = 2, 3, 4 and below 1 up to 2I = 12.</summary>
public class NumericalRangeBoundTests
{
    private static bool Decides(int twoI, long j0, long j1, long j2) =>
        NumericalRangeBound.ExcludesEveryCentre(new QuadrupolarGenerator(SpinSpace.Of(twoI), new QuadrupolarRates(j0, j1, j2), SecularHamiltonian.None));

    [Theory]
    [InlineData(1, 2, 3, true)]
    [InlineData(1, 1, 1, true)]
    [InlineData(5, 1, 7, true)]
    [InlineData(0, 2, 3, false)]
    [InlineData(1, 2, 0, false)]
    [InlineData(1, 0, 3, false)]
    [InlineData(1, 0, 0, false)]
    public void AtThreeHalves_DecidesExactlyWhenEveryDensityIsPositive(long j0, long j1, long j2, bool expected)
    {
        Assert.Equal(expected, Decides(3, j0, j1, j2));
    }

    /// <summary>At spin 1 on J₀ = 0 the sodium page's two rays with a centre (a tuned splitting),
    /// J₁ = 4J₂/5 and J₁ = 8J₂, are where the bound falls silent: between them it decides, on them
    /// r_max = 2s exactly (a zero pivot), beyond them it is silent (where the page's moment route
    /// excludes a centre by other means).</summary>
    [Theory]
    [InlineData(0, 4, 5, false)]     // the first ray
    [InlineData(0, 8, 1, false)]     // the second ray
    [InlineData(0, 41, 50, true)]    // J₁/J₂ = 0.82, just inside
    [InlineData(0, 79, 10, true)]    // 7.9, just inside
    [InlineData(0, 6, 6, true)]
    [InlineData(0, 39, 50, false)]   // 0.78, just outside
    [InlineData(0, 81, 10, false)]   // 8.1, just outside
    public void AtSpinOne_TheBoundFallsSilentExactlyOnThePagesTwoRays(long j0, long j1, long j2, bool expected)
    {
        Assert.Equal(expected, Decides(2, j0, j1, j2));
    }

    [Fact]
    public void ExtremeNarrowingRatio_ClosedForm_AndTheDecisionAgree()
    {
        Assert.Equal(new BigRational(3, 4), NumericalRangeBound.ExtremeNarrowingRatio(SpinSpace.Of(2)));
        Assert.Equal(new BigRational(4, 5), NumericalRangeBound.ExtremeNarrowingRatio(SpinSpace.Of(3)));
        Assert.Equal(new BigRational(11, 14), NumericalRangeBound.ExtremeNarrowingRatio(SpinSpace.Of(4)));
        for (int twoI = 2; twoI <= 12; twoI++)
        {
            var ratio = NumericalRangeBound.ExtremeNarrowingRatio(SpinSpace.Of(twoI));
            Assert.True((ratio - 1).Sign < 0);
            Assert.True(Decides(twoI, 1, 1, 1));
        }
    }

    [Fact]
    public void Pivots_ControlFails_OnASemidefiniteMatrix()
    {
        var semidefinite = new BigRational[,] { { 1, 1 }, { 1, 1 } };
        var definite = new BigRational[,] { { 2, 1 }, { 1, 2 } };
        Assert.False(NumericalRangeBound.AllPivotsPositive(semidefinite));
        Assert.True(NumericalRangeBound.AllPivotsPositive(definite));
    }
}
