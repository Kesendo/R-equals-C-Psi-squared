using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F110 SLOW enumeration tests: re-run the F103/F105/F106 anchors and
/// assert F110's Aspect A (hard only in diagonal cell), Aspect B (Y-inversion),
/// and Aspect C (k-purity sharpening 42:8 to 228:0) bit-exactly.
///
/// <para>Runtime: ~10s for k=3 N=4, ~30s for k=3 N=5, ~3min PLINQ for k=4 N=4.
/// Skip-by-default in CI; manual re-run via
/// <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" -c Release --filter "Category=SLOW_F110"</c>.</para></summary>
[Trait("Category", "SLOW_F110")]
public class HardCellYInversionPatternEnumerationTests :
    IClassFixture<HardCellYInversionPatternCountsFixture>
{
    private readonly HardCellYInversionPatternCountsFixture _fixture;

    public HardCellYInversionPatternEnumerationTests(HardCellYInversionPatternCountsFixture fixture)
    {
        _fixture = fixture;
    }

    [Fact]
    public void K3N4_HardDiagonalCounts_Match42_8_WithYInversion()
    {
        var counts = _fixture.CountsK3N4;
        // Z-deph + Klein (0,1) -> (42, 8): dominant y_par=0 per F110 Aspect B
        Assert.Equal((42, 8), HardCounts(counts, (0, 1), 'Z'));
        // X-deph + Klein (1,0) -> (42, 8): dominant y_par=0
        Assert.Equal((42, 8), HardCounts(counts, (1, 0), 'X'));
        // Y-deph + Klein (1,1) -> (8, 42): Y-INVERSION dominant y_par=1
        Assert.Equal((8, 42), HardCounts(counts, (1, 1), 'Y'));
    }

    [Fact]
    public void K3N4_HardOnlyInDiagonalCells_AllOtherCellsZeroHard()
    {
        // Aspect A: hard does not appear in any non-diagonal cell.
        var counts = _fixture.CountsK3N4;
        foreach (var dephase in new[] { 'Z', 'X', 'Y' })
        {
            var dephaseLetter = CharToDephase(dephase);
            var diagonalCell = HardCellYInversionPattern.DiagonalKleinCellForDephase(dephaseLetter);
            foreach (var klein in AllKleinCells)
            {
                if (klein == diagonalCell) continue;
                Assert.Equal(0, HardCount(counts, klein, dephase, 0));
                Assert.Equal(0, HardCount(counts, klein, dephase, 1));
            }
        }
    }

    [Fact]
    public void K3N5_HardDiagonalCounts_NStableWithK3N4()
    {
        // F105 N-stability: identical 42:8 at N=5 k=3 as at N=4 k=3
        var counts = _fixture.CountsK3N5;
        Assert.Equal((42, 8), HardCounts(counts, (0, 1), 'Z'));
        Assert.Equal((42, 8), HardCounts(counts, (1, 0), 'X'));
        Assert.Equal((8, 42), HardCounts(counts, (1, 1), 'Y'));
    }

    [Fact]
    public void K4N4_HardDiagonalCounts_Match228_0_PureWithYInversion()
    {
        // F106 k-sharpening: pure 228:0 with Y-inversion preserved
        var counts = _fixture.CountsK4N4;
        Assert.Equal((228, 0), HardCounts(counts, (0, 1), 'Z'));
        Assert.Equal((228, 0), HardCounts(counts, (1, 0), 'X'));
        Assert.Equal((0, 228), HardCounts(counts, (1, 1), 'Y'));
    }

    [Fact]
    public void K4N4_HardOnlyInDiagonalCells_AllOtherCellsZeroHard()
    {
        var counts = _fixture.CountsK4N4;
        foreach (var dephase in new[] { 'Z', 'X', 'Y' })
        {
            var dephaseLetter = CharToDephase(dephase);
            var diagonalCell = HardCellYInversionPattern.DiagonalKleinCellForDephase(dephaseLetter);
            foreach (var klein in AllKleinCells)
            {
                if (klein == diagonalCell) continue;
                Assert.Equal(0, HardCount(counts, klein, dephase, 0));
                Assert.Equal(0, HardCount(counts, klein, dephase, 1));
            }
        }
    }

    [Fact]
    public void DominantYParity_StructuralReading_MatchesEmpiricalDominance()
    {
        // Verify HardCellYInversionPattern.DominantYParityForDephase predicts the
        // empirically dominant y_par bit in the hard cell for each dephase letter,
        // at both k=3 and k=4.
        VerifyDominanceAt(_fixture.CountsK3N4);
        VerifyDominanceAt(_fixture.CountsK4N4);
    }

    private static void VerifyDominanceAt(
        IReadOnlyDictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> counts)
    {
        foreach (var dephase in new[] { PauliLetter.Z, PauliLetter.X, PauliLetter.Y })
        {
            var klein = HardCellYInversionPattern.DiagonalKleinCellForDephase(dephase);
            var expectedDominant = HardCellYInversionPattern.DominantYParityForDephase(dephase);
            int yPar0Hard = HardCount(counts, klein, DephaseChar(dephase), 0);
            int yPar1Hard = HardCount(counts, klein, DephaseChar(dephase), 1);
            int dominantCount = expectedDominant == 0 ? yPar0Hard : yPar1Hard;
            int otherCount = expectedDominant == 0 ? yPar1Hard : yPar0Hard;
            Assert.True(dominantCount > otherCount,
                $"F110 Aspect B violated at klein={klein}, dephase={dephase}: " +
                $"expected dominant y_par={expectedDominant}, got y_par=0={yPar0Hard} vs y_par=1={yPar1Hard}");
        }
    }

    // ---------- Helpers ----------

    private static readonly (int, int)[] AllKleinCells =
        { (0, 0), (0, 1), (1, 0), (1, 1) };

    private static PauliLetter CharToDephase(char c) => c switch
    {
        'Z' => PauliLetter.Z,
        'X' => PauliLetter.X,
        'Y' => PauliLetter.Y,
        _ => throw new ArgumentOutOfRangeException(nameof(c)),
    };

    private static char DephaseChar(PauliLetter l) => l switch
    {
        PauliLetter.Z => 'Z',
        PauliLetter.X => 'X',
        PauliLetter.Y => 'Y',
        _ => throw new ArgumentOutOfRangeException(nameof(l)),
    };

    private static (int yPar0, int yPar1) HardCounts(
        IReadOnlyDictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> counts,
        (int, int) klein, char dephase) =>
        (HardCount(counts, klein, dephase, 0), HardCount(counts, klein, dephase, 1));

    private static int HardCount(
        IReadOnlyDictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> counts,
        (int, int) klein, char dephase, int yPar) =>
        counts.GetValueOrDefault(((klein.Item1, klein.Item2), dephase, yPar, TrichotomyClass.Hard));
}

/// <summary>Shared lazy classification grids at (N=4, k=3), (N=5, k=3), and (N=4, k=4).
/// Each grid is computed once on first access; subsequent test methods reuse the
/// same dictionary. Thread-safe via <see cref="Lazy{T}"/>.</summary>
public sealed class HardCellYInversionPatternCountsFixture
{
    private static readonly PauliLetter[] DephaseLetters =
        { PauliLetter.Z, PauliLetter.X, PauliLetter.Y };

    private readonly Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>> _countsK3N4;
    private readonly Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>> _countsK3N5;
    private readonly Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>> _countsK4N4;

    public HardCellYInversionPatternCountsFixture()
    {
        _countsK3N4 = new Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>>(
            () => ClassifyAndGroup(N: 4, k: 3), System.Threading.LazyThreadSafetyMode.ExecutionAndPublication);
        _countsK3N5 = new Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>>(
            () => ClassifyAndGroup(N: 5, k: 3), System.Threading.LazyThreadSafetyMode.ExecutionAndPublication);
        _countsK4N4 = new Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>>(
            () => ClassifyAndGroup(N: 4, k: 4), System.Threading.LazyThreadSafetyMode.ExecutionAndPublication);
    }

    public IReadOnlyDictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> CountsK3N4 => _countsK3N4.Value;
    public IReadOnlyDictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> CountsK3N5 => _countsK3N5.Value;
    public IReadOnlyDictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> CountsK4N4 => _countsK4N4.Value;

    private static char DephaseChar(PauliLetter l) => l switch
    {
        PauliLetter.Z => 'Z',
        PauliLetter.X => 'X',
        PauliLetter.Y => 'Y',
        _ => throw new ArgumentOutOfRangeException(nameof(l)),
    };

    private static Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>
        ClassifyAndGroup(int N, int k)
    {
        var chain = new ChainSystem(N: N, J: 1.0, GammaZero: 0.05);
        var classifications = Z2HomogeneousKBodyEnumeration.Enumerate(k)
            .SelectMany(item => DephaseLetters.Select(d => (item, dephase: d)))
            .AsParallel()
            .Select(x =>
            {
                var templates = new[] { x.item.Term1, x.item.Term2 };
                var cls = PauliPairTrichotomy.Classify(chain, templates, dephaseLetter: x.dephase);
                return (x.item.Klein, Dephase: DephaseChar(x.dephase), x.item.YPar, Cls: cls);
            })
            .ToList();

        return classifications
            .GroupBy(c => (c.Klein, c.Dephase, c.YPar, c.Cls))
            .ToDictionary(g => g.Key, g => g.Count());
    }
}
