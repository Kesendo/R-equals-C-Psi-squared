using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F104 verification: re-classify F103's 294 N=4 k=3 Z₂³-homogeneous pairs
/// in C# (via the new <see cref="PauliPairTrichotomy.Classify"/> k≥3 overload) and
/// assert per-cell counts match <see cref="F87Z2CubedRefinement"/>'s frozen records
/// bit-exactly. Closes F103's explicit out-of-scope item "C# k≥3 classifier lift".
///
/// <para>Runtime: ~30-60s for 882 classifications (294 pairs × 3 dephase letters).
/// Tagged Slow so CI can filter via <c>--filter "Category!=Slow"</c> when needed.</para></summary>
[Trait("Category", "Slow")]
public class F104KBodyTrichotomyVerificationTests
{
    private static readonly PauliLetter[] AllLetters =
        { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };

    private static readonly PauliLetter[] DephaseLetters =
        { PauliLetter.Z, PauliLetter.X, PauliLetter.Y };

    private static ChainSystem MakeChainN4() => new(N: 4, J: 1.0, GammaZero: 0.05);

    /// <summary>Enumerate unordered (term1, term2) pairs at N=4, k=3 such that both
    /// terms are non-trivial (not all-I), share the same Klein index, and share the
    /// same Y-parity. Same algorithm as Python <c>enumerate_z2_homogeneous_k3</c>.</summary>
    private static List<(PauliTerm Term1, PauliTerm Term2, (int A, int B) Klein, int YPar)> EnumerateZ2HomogeneousK3()
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

                                string key = string.Compare(LettersKey(term1Letters), LettersKey(term2Letters), StringComparison.Ordinal) <= 0
                                    ? LettersKey(term1Letters) + "|" + LettersKey(term2Letters)
                                    : LettersKey(term2Letters) + "|" + LettersKey(term1Letters);
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

    /// <summary>Classify all 294 enumerated pairs × 3 dephase letters and group by
    /// (Klein, dephase letter, y_par, trichotomy class). Returns a fully-populated
    /// count grid for downstream assertion.</summary>
    private static Dictionary<((int A, int B) Klein, char Dephase, int YPar, TrichotomyClass Cls), int> ClassifyAndGroup()
    {
        var chain = MakeChainN4();
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

    [Fact]
    public void Verification_Enumerate294Pairs_TotalCountMatches()
    {
        var items = EnumerateZ2HomogeneousK3();
        Assert.Equal(294, items.Count);
    }

    [Fact]
    public void Verification_TrulyCountsMatch_F103_FrozenAnchor()
    {
        var counts = ClassifyAndGroup();
        var f103 = new F87Z2CubedRefinement();

        int totalTruly = counts
            .Where(kv => kv.Key.Cls == TrichotomyClass.Truly)
            .Sum(kv => kv.Value);
        int trulyYParOne = counts
            .Where(kv => kv.Key.Cls == TrichotomyClass.Truly && kv.Key.YPar == 1)
            .Sum(kv => kv.Value);

        Assert.Equal(f103.TrulyPurity.TotalTrulyClassifications, totalTruly);
        Assert.Equal(f103.TrulyPurity.YParityOneCount, trulyYParOne);
    }

    [Fact]
    public void Verification_HardDiagonalCounts_MatchF103()
    {
        var counts = ClassifyAndGroup();
        var f103 = new F87Z2CubedRefinement();

        var zHard = ((0, 1), 'Z');
        var xHard = ((1, 0), 'X');
        var yHard = ((1, 1), 'Y');

        Assert.Equal(f103.HardDiagonal.ZDephKlein01,
            (counts.GetValueOrDefault((zHard.Item1, zHard.Item2, 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault((zHard.Item1, zHard.Item2, 1, TrichotomyClass.Hard))));
        Assert.Equal(f103.HardDiagonal.XDephKlein10,
            (counts.GetValueOrDefault((xHard.Item1, xHard.Item2, 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault((xHard.Item1, xHard.Item2, 1, TrichotomyClass.Hard))));
        Assert.Equal(f103.HardDiagonal.YDephKlein11,
            (counts.GetValueOrDefault((yHard.Item1, yHard.Item2, 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault((yHard.Item1, yHard.Item2, 1, TrichotomyClass.Hard))));
    }

    [Fact]
    public void Verification_SoftCounts_AllPatternsMatch_F103()
    {
        var counts = ClassifyAndGroup();
        var f103 = new F87Z2CubedRefinement();

        // DiagonalSoft 13:13 universal
        Assert.Equal(f103.DiagonalSoft.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f103.DiagonalSoft.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f103.DiagonalSoft.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Soft))));

        // MotherSoft 0:21 across 3 letters
        Assert.Equal(f103.MotherSoft.ZDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f103.MotherSoft.XDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f103.MotherSoft.YDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Y', 1, TrichotomyClass.Soft))));

        // Off-diagonal 6 cells: Pattern B + Pattern C
        foreach (var (key, expected) in f103.OffDiagonalSoft.Cells)
        {
            var (kA, kB, dephase) = key;
            var actual = (counts.GetValueOrDefault(((kA, kB), dephase, 0, TrichotomyClass.Soft)),
                          counts.GetValueOrDefault(((kA, kB), dephase, 1, TrichotomyClass.Soft)));
            Assert.Equal(expected, actual);
        }
    }

    [Fact]
    public void Verification_F103_FrozenCounts_AreReproducibleFromCSharp()
    {
        // Top-level composite: enumerate, classify, group, and verify every cell in
        // F103's frozen records matches the C# re-classification bit-exactly. Catches
        // any drift not caught by the per-record tests above.
        var counts = ClassifyAndGroup();
        var f103 = new F87Z2CubedRefinement();

        // Truly: 300 total, 0 y_par=1
        int trulyTotal = counts.Where(kv => kv.Key.Cls == TrichotomyClass.Truly).Sum(kv => kv.Value);
        int trulyYOne = counts.Where(kv => kv.Key.Cls == TrichotomyClass.Truly && kv.Key.YPar == 1).Sum(kv => kv.Value);
        Assert.Equal(f103.TrulyPurity.TotalTrulyClassifications, trulyTotal);
        Assert.Equal(f103.TrulyPurity.YParityOneCount, trulyYOne);

        // Hard diagonals
        Assert.Equal(f103.HardDiagonal.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Hard))));
        Assert.Equal(f103.HardDiagonal.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Hard))));
        Assert.Equal(f103.HardDiagonal.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Hard)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Hard))));

        // Diagonal soft 13:13
        Assert.Equal(f103.DiagonalSoft.ZDephKlein01,
            (counts.GetValueOrDefault(((0, 1), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 1), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f103.DiagonalSoft.XDephKlein10,
            (counts.GetValueOrDefault(((1, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f103.DiagonalSoft.YDephKlein11,
            (counts.GetValueOrDefault(((1, 1), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((1, 1), 'Y', 1, TrichotomyClass.Soft))));

        // Mother soft
        Assert.Equal(f103.MotherSoft.ZDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Z', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Z', 1, TrichotomyClass.Soft))));
        Assert.Equal(f103.MotherSoft.XDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'X', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'X', 1, TrichotomyClass.Soft))));
        Assert.Equal(f103.MotherSoft.YDephCounts,
            (counts.GetValueOrDefault(((0, 0), 'Y', 0, TrichotomyClass.Soft)),
             counts.GetValueOrDefault(((0, 0), 'Y', 1, TrichotomyClass.Soft))));

        // Off-diagonal soft cells
        foreach (var (key, expected) in f103.OffDiagonalSoft.Cells)
        {
            var (kA, kB, dephase) = key;
            var actual = (counts.GetValueOrDefault(((kA, kB), dephase, 0, TrichotomyClass.Soft)),
                          counts.GetValueOrDefault(((kA, kB), dephase, 1, TrichotomyClass.Soft)));
            Assert.Equal(expected, actual);
        }
    }
}
