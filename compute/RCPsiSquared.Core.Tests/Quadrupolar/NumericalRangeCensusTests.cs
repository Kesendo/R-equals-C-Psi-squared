using System.Text;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Quadrupolar;

/// <summary>A reading, not a law: per spin 2I = 2 … 12, the integer points J₀ + J₁ + J₂ = 12 (91 points)
/// where the numerical-range bound decides that no Hermitian H gives a centre, printed for the page.
/// The asserts are only what is known exactly: the extreme-narrowing point decides at every spin, and at
/// 3/2 the bound decides exactly where all three densities are positive.</summary>
public class NumericalRangeCensusTests
{
    private readonly ITestOutputHelper _out;
    public NumericalRangeCensusTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Census_OnTheSimplex_J0PlusJ1PlusJ2Is12()
    {
        const int total = 12;
        var sb = new StringBuilder();
        for (int twoI = 2; twoI <= 12; twoI++)
        {
            var s = SpinSpace.Of(twoI);
            int decided = 0, points = 0;
            var silent = new StringBuilder();
            for (int j0 = 0; j0 <= total; j0++)
                for (int j1 = 0; j0 + j1 <= total; j1++)
                {
                    int j2 = total - j0 - j1;
                    points++;
                    bool d = NumericalRangeBound.ExcludesEveryCentre(new QuadrupolarGenerator(s, new QuadrupolarRates(j0, j1, j2), SecularHamiltonian.None));
                    if (d) decided++; else silent.Append($" ({j0},{j1},{j2})");
                    if (j0 == 4 && j1 == 4) Assert.True(d);
                    if (twoI == 3) Assert.Equal(j0 > 0 && j1 > 0 && j2 > 0, d);
                }
            sb.AppendLine($"2I = {twoI}: EN ratio {NumericalRangeBound.ExtremeNarrowingRatio(s)}, decided {decided} of {points}; silent at{silent}");
        }
        _out.WriteLine(sb.ToString());
    }
}
