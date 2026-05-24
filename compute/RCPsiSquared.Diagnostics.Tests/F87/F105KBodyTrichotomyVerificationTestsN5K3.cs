using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F105 verification: re-classify F105's 294 N=5 k=3 Z₂³-homogeneous pairs
/// in C# (via the <see cref="PauliPairTrichotomy.Classify"/> k≥3 overload at N=5) and
/// assert per-cell counts match <see cref="F87Z2CubedRefinementN5K3"/>'s frozen records
/// bit-exactly. On-demand re-verification mechanism for the F105 anchor; parallel to
/// <see cref="F104KBodyTrichotomyVerificationTests"/> at N=4.
///
/// <para>Runtime: ~3h dense at N=5 (882 classifications × per-call cost). The Task 7
/// SLOW_F105_BATCH tool run already supplied the JSON anchor; this verification suite
/// is skip-by-default in CI. Manual re-run via
/// <c>dotnet test "compute\RCPsiSquared.Diagnostics.Tests" -c Release --filter "Category=SLOW_F105"</c>.
/// The full classification grid is computed once via <see cref="F105CountsFixtureN5K3"/>
/// and shared across all 4 cell-asserting tests (xUnit <c>IClassFixture</c>).</para></summary>
[Trait("Category", "SLOW_F105")]
public class F105KBodyTrichotomyVerificationTestsN5K3 : IClassFixture<F105CountsFixtureN5K3>
{
    private readonly F105CountsFixtureN5K3 _fixture;

    public F105KBodyTrichotomyVerificationTestsN5K3(F105CountsFixtureN5K3 fixture)
    {
        _fixture = fixture;
    }

    [Fact]
    public void Verification_Enumerate294Pairs_TotalCountMatches()
    {
        var items = F105CountsFixtureN5K3.EnumerateZ2HomogeneousK3();
        Assert.Equal(294, items.Count);
    }

    [Fact]
    public void Verification_TrulyCountsMatch_F105_FrozenAnchor()
    {
        var counts = _fixture.Counts;
        var f105 = new F87Z2CubedRefinementN5K3();

        int totalTruly = counts
            .Where(kv => kv.Key.Cls == TrichotomyClass.Truly)
            .Sum(kv => kv.Value);
        int trulyYParOne = counts
            .Where(kv => kv.Key.Cls == TrichotomyClass.Truly && kv.Key.YPar == 1)
            .Sum(kv => kv.Value);

        Assert.Equal(f105.TrulyPurity.TotalTrulyClassifications, totalTruly);
        Assert.Equal(f105.TrulyPurity.YParityOneCount, trulyYParOne);
    }

    [Fact]
    public void Verification_HardDiagonalCounts_MatchF105()
    {
        var counts = _fixture.Counts;
        var f105 = new F87Z2CubedRefinementN5K3();

        Assert.Equal(f105.HardDiagonal.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Hard))));
        Assert.Equal(f105.HardDiagonal.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Hard))));
        Assert.Equal(f105.HardDiagonal.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Hard))));
    }

    [Fact]
    public void Verification_SoftCounts_AllPatternsMatch_F105()
    {
        var counts = _fixture.Counts;
        var f105 = new F87Z2CubedRefinementN5K3();

        // DiagonalSoft 13:13 universal
        Assert.Equal(f105.DiagonalSoft.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f105.DiagonalSoft.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f105.DiagonalSoft.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Soft))));

        // MotherSoft 0:21 across 3 letters
        Assert.Equal(f105.MotherSoft.ZDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f105.MotherSoft.XDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f105.MotherSoft.YDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Y', 1, TrichotomyClass.Soft))));

        // Off-diagonal 6 cells: Pattern B + Pattern C
        foreach (var (key, expected) in f105.OffDiagonalSoft.Cells)
        {
            var (kA, kB, dephase) = key;
            var actual = (counts.GetValueOrDefault(((kA, kB), dephase, 0, TrichotomyClass.Soft)),
                          counts.GetValueOrDefault(((kA, kB), dephase, 1, TrichotomyClass.Soft)));
            Assert.Equal(expected, actual);
        }
    }

    [Fact]
    public void Verification_F105_FrozenCounts_AreReproducibleFromCSharp()
    {
        // Top-level composite: enumerate, classify, group, and verify every cell in
        // F105's frozen records matches the C# re-classification bit-exactly. Catches
        // any drift not caught by the per-record tests above.
        var counts = _fixture.Counts;
        var f105 = new F87Z2CubedRefinementN5K3();

        // Truly: 300 total, 0 y_par=1
        int trulyTotal = counts.Where(kv => kv.Key.Cls == TrichotomyClass.Truly).Sum(kv => kv.Value);
        int trulyYOne = counts.Where(kv => kv.Key.Cls == TrichotomyClass.Truly && kv.Key.YPar == 1).Sum(kv => kv.Value);
        Assert.Equal(f105.TrulyPurity.TotalTrulyClassifications, trulyTotal);
        Assert.Equal(f105.TrulyPurity.YParityOneCount, trulyYOne);

        // Hard diagonals
        Assert.Equal(f105.HardDiagonal.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Hard))));
        Assert.Equal(f105.HardDiagonal.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Hard))));
        Assert.Equal(f105.HardDiagonal.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Hard))));

        // Diagonal soft 13:13
        Assert.Equal(f105.DiagonalSoft.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f105.DiagonalSoft.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f105.DiagonalSoft.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Soft))));

        // Mother soft
        Assert.Equal(f105.MotherSoft.ZDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f105.MotherSoft.XDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f105.MotherSoft.YDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Y', 1, TrichotomyClass.Soft))));

        // Off-diagonal soft cells
        foreach (var (key, expected) in f105.OffDiagonalSoft.Cells)
        {
            var (kA, kB, dephase) = key;
            var actual = (counts.GetValueOrDefault(((kA, kB), dephase, 0, TrichotomyClass.Soft)),
                          counts.GetValueOrDefault(((kA, kB), dephase, 1, TrichotomyClass.Soft)));
            Assert.Equal(expected, actual);
        }
    }
}

/// <summary>Shared cache of the F105 classification grid at N=5 (882 = 294 pairs × 3
/// dephase letters), reused across all tests in
/// <see cref="F105KBodyTrichotomyVerificationTestsN5K3"/>. Eliminates per-test
/// re-classification: <see cref="PauliPairTrichotomy.Classify"/> is called exactly 882
/// times per class run instead of 882 × (number of tests calling it). Thread-safe via
/// <see cref="Lazy{T}"/> with ExecutionAndPublication mode.</summary>
public sealed class F105CountsFixtureN5K3
{
    private static readonly PauliLetter[] AllLetters =
        { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };

    private static readonly PauliLetter[] DephaseLetters =
        { PauliLetter.Z, PauliLetter.X, PauliLetter.Y };

    private readonly Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>> _counts;

    public F105CountsFixtureN5K3()
    {
        _counts = new Lazy<Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int>>(
            ClassifyAndGroup, System.Threading.LazyThreadSafetyMode.ExecutionAndPublication);
    }

    /// <summary>Lazily computed on first access; returns the same dictionary instance on
    /// every subsequent call across the test class.</summary>
    public Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> Counts => _counts.Value;

    private static ChainSystem MakeChainN5() => new(N: 5, J: 1.0, GammaZero: 0.05);

    /// <summary>Enumerate unordered (term1, term2) pairs at k=3 such that both terms
    /// are non-trivial (not all-I), share the same Klein index, and share the same
    /// Y-parity. Same algorithm as F104's <c>EnumerateZ2HomogeneousK3</c>; pair count
    /// is N-independent at k=3 (294 pairs regardless of chain length).</summary>
    public static List<(PauliTerm Term1, PauliTerm Term2, (int A, int B) Klein, int YPar)> EnumerateZ2HomogeneousK3()
    {
        var seen = new HashSet<string>();
        var items = new List<(PauliTerm, PauliTerm, (int, int), int)>();

        foreach (var a in AllLetters)
            foreach (var b in AllLetters)
                foreach (var c in AllLetters)
                {
                    var term1Letters = new[] { a, b, c };
                    if (term1Letters.All(l => l == PauliLetter.I)) continue;
                    var term1 = new PauliTerm(term1Letters, Complex.One);

                    foreach (var d in AllLetters)
                        foreach (var e in AllLetters)
                            foreach (var f in AllLetters)
                            {
                                var term2Letters = new[] { d, e, f };
                                if (term2Letters.All(l => l == PauliLetter.I)) continue;
                                var term2 = new PauliTerm(term2Letters, Complex.One);

                                if (term1.KleinIndex != term2.KleinIndex) continue;
                                if (term1.YParity != term2.YParity) continue;

                                string key1 = LettersKey(term1Letters);
                                string key2 = LettersKey(term2Letters);
                                string key = string.Compare(key1, key2, StringComparison.Ordinal) <= 0
                                    ? key1 + "|" + key2
                                    : key2 + "|" + key1;
                                if (!seen.Add(key)) continue;

                                items.Add((term1, term2, term1.KleinIndex, term1.YParity));
                            }
                }

        return items;
    }

    private static string LettersKey(IReadOnlyList<PauliLetter> letters) =>
        string.Concat(letters.Select(LetterChar));

    private static char LetterChar(PauliLetter l) => l switch
    {
        PauliLetter.I => 'I',
        PauliLetter.X => 'X',
        PauliLetter.Y => 'Y',
        PauliLetter.Z => 'Z',
        _ => throw new ArgumentOutOfRangeException(nameof(l)),
    };

    private static char DephaseChar(PauliLetter l) => l switch
    {
        PauliLetter.Z => 'Z',
        PauliLetter.X => 'X',
        PauliLetter.Y => 'Y',
        _ => throw new ArgumentOutOfRangeException(nameof(l)),
    };

    /// <summary>Classify all 294 enumerated pairs × 3 dephase letters at N=5 and group
    /// by (Klein, dephase letter, y_par, trichotomy class). Returns a fully-populated
    /// count grid for downstream assertion.</summary>
    private static Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> ClassifyAndGroup()
    {
        var chain = MakeChainN5();
        var counts = new Dictionary<((int, int), char, int, TrichotomyClass), int>();

        foreach (var item in EnumerateZ2HomogeneousK3())
        {
            var templates = new[] { item.Term1, item.Term2 };
            foreach (var dephase in DephaseLetters)
            {
                var cls = PauliPairTrichotomy.Classify(chain, templates,
                    dephaseLetter: dephase);
                var key = (item.Klein, DephaseChar(dephase), item.YPar, cls);
                counts.TryGetValue(key, out int existing);
                counts[key] = existing + 1;
            }
        }

        return counts;
    }
}
