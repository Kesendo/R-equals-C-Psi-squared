using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F106 verification: re-classify the 4248 N=4 k=4 Z₂³-homogeneous pairs in C#
/// (via the <see cref="PauliPairTrichotomy.Classify"/> k≥3 overload at N=4) and assert
/// per-cell counts match <see cref="F87Z2CubedRefinementN4K4"/>'s frozen records
/// bit-exactly. On-demand re-verification mechanism for the F106 anchor; parallel to
/// <see cref="F105KBodyTrichotomyVerificationTestsN5K3"/> at k=3.
///
/// <para>Runtime: ~2-3min PLINQ on 24 cores at N=4 (12744 classifications). The Task 13
/// SLOW_F106_BATCH tool run already supplied the JSON anchor; this verification suite
/// is skip-by-default in CI. Manual re-run via
/// <c>dotnet test "compute\RCPsiSquared.Diagnostics.Tests" -c Release --filter "Category=SLOW_F106"</c>.
/// The full classification grid is computed once via <see cref="F106CountsFixtureN4K4"/>
/// and shared across all 4 cell-asserting tests (xUnit <c>IClassFixture</c>).</para></summary>
[Trait("Category", "SLOW_F106")]
public class F106KBodyTrichotomyVerificationTestsN4K4 : IClassFixture<F106CountsFixtureN4K4>
{
    private readonly F106CountsFixtureN4K4 _fixture;

    public F106KBodyTrichotomyVerificationTestsN4K4(F106CountsFixtureN4K4 fixture)
    {
        _fixture = fixture;
    }

    [Fact]
    public void Verification_Enumerate4248Pairs_TotalCountMatches()
    {
        var items = Z2HomogeneousKBodyEnumeration.Enumerate(4);
        Assert.Equal(4248, items.Count);
    }

    [Fact]
    public void Verification_TrulyCountsMatch_F106_FrozenAnchor()
    {
        var counts = _fixture.Counts;
        var f106 = new F87Z2CubedRefinementN4K4();

        int totalTruly = counts
            .Where(kv => kv.Key.Cls == TrichotomyClass.Truly)
            .Sum(kv => kv.Value);
        int trulyYParOne = counts
            .Where(kv => kv.Key.Cls == TrichotomyClass.Truly && kv.Key.YPar == 1)
            .Sum(kv => kv.Value);

        Assert.Equal(f106.TrulyPurity.TotalTrulyClassifications, totalTruly);
        Assert.Equal(f106.TrulyPurity.YParityOneCount, trulyYParOne);
    }

    [Fact]
    public void Verification_HardDiagonalCounts_MatchF106()
    {
        var counts = _fixture.Counts;
        var f106 = new F87Z2CubedRefinementN4K4();

        Assert.Equal(f106.HardDiagonal.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Hard))));
        Assert.Equal(f106.HardDiagonal.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Hard))));
        Assert.Equal(f106.HardDiagonal.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Hard))));
    }

    [Fact]
    public void Verification_SoftCounts_AllPatternsMatch_F106()
    {
        var counts = _fixture.Counts;
        var f106 = new F87Z2CubedRefinementN4K4();

        // DiagonalSoft 3 cells (Klein matches dephase)
        Assert.Equal(f106.DiagonalSoft.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f106.DiagonalSoft.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f106.DiagonalSoft.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Soft))));

        // MotherSoft Klein (0,0) across 3 letters
        Assert.Equal(f106.MotherSoft.ZDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f106.MotherSoft.XDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f106.MotherSoft.YDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Y', 1, TrichotomyClass.Soft))));

        // Off-diagonal cells: whatever is in OffDiagonalSoft.Cells
        foreach (var (key, expected) in f106.OffDiagonalSoft.Cells)
        {
            var (kA, kB, dephase) = key;
            var actual = (counts.GetValueOrDefault(((kA, kB), dephase, 0, TrichotomyClass.Soft)),
                          counts.GetValueOrDefault(((kA, kB), dephase, 1, TrichotomyClass.Soft)));
            Assert.Equal(expected, actual);
        }
    }

    [Fact]
    public void Verification_F106_FrozenCounts_AreReproducibleFromCSharp()
    {
        // Top-level composite: enumerate, classify, group, and verify every cell in
        // F106's frozen records matches the C# re-classification bit-exactly. Catches
        // any drift not caught by the per-record tests above.
        var counts = _fixture.Counts;
        var f106 = new F87Z2CubedRefinementN4K4();

        // Truly: total + y_par=1 count
        int trulyTotal = counts.Where(kv => kv.Key.Cls == TrichotomyClass.Truly).Sum(kv => kv.Value);
        int trulyYOne = counts.Where(kv => kv.Key.Cls == TrichotomyClass.Truly && kv.Key.YPar == 1).Sum(kv => kv.Value);
        Assert.Equal(f106.TrulyPurity.TotalTrulyClassifications, trulyTotal);
        Assert.Equal(f106.TrulyPurity.YParityOneCount, trulyYOne);

        // Hard diagonals
        Assert.Equal(f106.HardDiagonal.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Hard))));
        Assert.Equal(f106.HardDiagonal.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Hard))));
        Assert.Equal(f106.HardDiagonal.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Hard))));

        // Diagonal soft
        Assert.Equal(f106.DiagonalSoft.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f106.DiagonalSoft.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f106.DiagonalSoft.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Soft))));

        // Mother soft
        Assert.Equal(f106.MotherSoft.ZDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f106.MotherSoft.XDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f106.MotherSoft.YDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Y', 1, TrichotomyClass.Soft))));

        // Off-diagonal soft cells
        foreach (var (key, expected) in f106.OffDiagonalSoft.Cells)
        {
            var (kA, kB, dephase) = key;
            var actual = (counts.GetValueOrDefault(((kA, kB), dephase, 0, TrichotomyClass.Soft)),
                          counts.GetValueOrDefault(((kA, kB), dephase, 1, TrichotomyClass.Soft)));
            Assert.Equal(expected, actual);
        }
    }
}

/// <summary>Shared cache of the F106 classification grid at N=4 (12744 = 4248 pairs × 3
/// dephase letters), reused across all tests in
/// <see cref="F106KBodyTrichotomyVerificationTestsN4K4"/>. Eliminates per-test
/// re-classification: <see cref="PauliPairTrichotomy.Classify"/> is called exactly
/// 12744 times per class run instead of 12744 × (number of tests calling it).
/// Thread-safe via <see cref="Lazy{T}"/> with ExecutionAndPublication mode.</summary>
public sealed class F106CountsFixtureN4K4
{
    private static readonly PauliLetter[] DephaseLetters =
        { PauliLetter.Z, PauliLetter.X, PauliLetter.Y };

    private readonly Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>> _counts;

    public F106CountsFixtureN4K4()
    {
        _counts = new Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>>(
            ClassifyAndGroup, System.Threading.LazyThreadSafetyMode.ExecutionAndPublication);
    }

    /// <summary>Lazily computed on first access; returns the same dictionary instance on
    /// every subsequent call across the test class.</summary>
    public Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> Counts => _counts.Value;

    private static ChainSystem MakeChainN4() => new(N: 4, J: 1.0, GammaZero: 0.05);

    private static char DephaseChar(PauliLetter l) => l switch
    {
        PauliLetter.Z => 'Z',
        PauliLetter.X => 'X',
        PauliLetter.Y => 'Y',
        _ => throw new ArgumentOutOfRangeException(nameof(l)),
    };

    /// <summary>Classify all 4248 enumerated pairs × 3 dephase letters at N=4 and group
    /// by (Klein, dephase letter, y_par, trichotomy class). Returns a fully-populated
    /// count grid for downstream assertion. Enumeration delegated to
    /// <see cref="Z2HomogeneousKBodyEnumeration.Enumerate"/> (shared single source of
    /// truth across F104/F105/F106).</summary>
    private static Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> ClassifyAndGroup()
    {
        // Parallelize over (pair × dephase letter). PauliPairTrichotomy.Classify is pure
        // (constructs fresh H/L/M per call); ChainSystem and PauliTerm are immutable records.
        // 12744 independent classifications saturate all available cores via PLINQ.
        var chain = MakeChainN4();
        var classifications = Z2HomogeneousKBodyEnumeration.Enumerate(4)
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
