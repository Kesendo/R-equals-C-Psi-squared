using System.Text;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Quadrupolar;

/// <summary>A reading, not a law: per spin 2I = 2 … 12, the integer points J₀ + J₁ + J₂ = 12 (91 points)
/// where the numerical-range bound decides that no Hermitian H gives a centre, printed for the page.
/// The asserts are only what is known exactly: the extreme-narrowing point decides at every spin, at 3/2 the
/// bound decides exactly where all three densities are positive, and it is silent at (12, 0, 0), (11, 1, 0) and
/// (11, 0, 1) at every spin (the page quotes these and the counts printed here).</summary>
public class NumericalRangeCensusTests
{
    private readonly ITestOutputHelper _out;
    public NumericalRangeCensusTests(ITestOutputHelper output) => _out = output;

    /// <summary>The decided counts at 2I = 2 … 12, pinned because the page quotes them (53 at I = 2, 25 at I = 6).</summary>
    private static readonly int[] DecidedCounts = { 57, 55, 53, 52, 46, 37, 35, 33, 31, 26, 25 };

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
                    if ((j0, j1, j2) is (12, 0, 0) or (11, 1, 0) or (11, 0, 1)) Assert.False(d);   // silent at J₀ alone and next to it
                }
            sb.AppendLine($"2I = {twoI}: EN ratio {NumericalRangeBound.ExtremeNarrowingRatio(s)}, decided {decided} of {points}; silent at{silent}");
            Assert.Equal(DecidedCounts[twoI - 2], decided);
        }
        _out.WriteLine(sb.ToString());
    }
}
