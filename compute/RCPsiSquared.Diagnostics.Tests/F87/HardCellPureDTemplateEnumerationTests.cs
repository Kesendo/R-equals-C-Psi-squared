using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F111 SLOW enumeration tests: classify all 4248 Z₂³-homogeneous k=4 N=4
/// pairs via <see cref="PauliPairTrichotomy.Classify"/> × 3 dephase letters
/// (12744 classifications, ~2-3 min PLINQ on 24 cores) and assert the Pure-D
/// Template Rule's hard/non-hard prediction matches the F87 classifier bit-exactly,
/// plus the F110 Aspect B Y-inversion corollary holds on every hard pair.
///
/// <para>Runtime: ~2-3 min for the shared per-pair classification (one-shot via
/// <see cref="F111PerPairClassificationFixture"/>; reused across all 9 tests in
/// this class via xUnit <see cref="IClassFixture{TFixture}"/>). Subsequent test
/// methods are O(1) lookups into the precomputed per-cell pair lists.</para>
///
/// <para>Skip-by-default in CI via <c>[Trait("Category", "SLOW_F111")]</c>; manual
/// re-run via:
/// <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" -c Release --filter "Category=SLOW_F111"</c>.</para></summary>
[Trait("Category", "SLOW_F111")]
public sealed class HardCellPureDTemplateEnumerationTests
    : IClassFixture<F111PerPairClassificationFixture>
{
    private readonly F111PerPairClassificationFixture _fixture;

    public HardCellPureDTemplateEnumerationTests(F111PerPairClassificationFixture fixture)
    {
        _fixture = fixture;
    }

    // ============================================================
    // Test 3: F87-hard pairs match Pure-D Template Rule prediction
    // ============================================================

    [Theory]
    [InlineData(PauliLetter.Z)]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void HardPairsInDiagonalCell_MatchPureDTemplateRule(PauliLetter dephase)
    {
        // Diagonal Klein cell is (D.BitA(), D.BitB()) by construction: Z → (0, 1),
        // X → (1, 0), Y → (1, 1). Deriving here rather than passing explicit InlineData
        // bits avoids the redundant bit-pun bookkeeping (matches Test 5's pattern).
        var (bitA, bitB) = (dephase.BitA(), dephase.BitB());
        var hardPairs = _fixture.GetHardPairsInCell(dephase, bitA, bitB);
        Assert.Equal(228, hardPairs.Count);
        foreach (var (p, q) in hardPairs)
        {
            Assert.True(
                HardCellPureDTemplate.IsPredictedHardAtK4N4(p, q, dephase),
                $"Pair ({p.Label}, {q.Label}) is F87-hard in {dephase}-deph diagonal cell ({bitA},{bitB}) " +
                $"but Pure-D Template Rule predicts NOT hard");
        }
    }

    // ============================================================
    // Test 4: F87-soft and F87-truly pairs are NOT predicted hard
    // ============================================================

    [Theory]
    [InlineData(PauliLetter.Z)]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void NonHardPairsInDiagonalCell_AreNotPredictedHard(PauliLetter dephase)
    {
        var (bitA, bitB) = (dephase.BitA(), dephase.BitB());
        var nonHardPairs = _fixture.GetSoftPairsInCell(dephase, bitA, bitB)
            .Concat(_fixture.GetTrulyPairsInCell(dephase, bitA, bitB))
            .ToList();
        foreach (var (p, q) in nonHardPairs)
        {
            Assert.False(
                HardCellPureDTemplate.IsPredictedHardAtK4N4(p, q, dephase),
                $"Pair ({p.Label}, {q.Label}) is NOT F87-hard in {dephase}-deph diagonal cell ({bitA},{bitB}) " +
                $"but Pure-D Template Rule predicts hard");
        }
    }

    // ============================================================
    // Test 5: Y-inversion corollary holds for all hard pairs (F110 Aspect B at k=4)
    // ============================================================

    [Theory]
    [InlineData(PauliLetter.Z, 0)]
    [InlineData(PauliLetter.X, 0)]
    [InlineData(PauliLetter.Y, 1)]
    public void HardPairsInDiagonalCell_AllSatisfyYInversionCorollary(
        PauliLetter dephase, int expectedYpar)
    {
        var diagonalCellKlein = (dephase.BitA(), dephase.BitB());
        var hardPairs = _fixture.GetHardPairsInCell(
            dephase, diagonalCellKlein.Item1, diagonalCellKlein.Item2);
        Assert.Equal(228, hardPairs.Count);
        foreach (var (p, q) in hardPairs)
        {
            Assert.True(
                HardCellPureDTemplate.VerifyYInversionCorollaryAtK4N4(p, q, dephase),
                $"Y-inversion corollary fails for pair ({p.Label}, {q.Label}) in {dephase}-deph diagonal cell");
            Assert.Equal(expectedYpar, p.YParity);
            Assert.Equal(expectedYpar, q.YParity);
        }
    }
}

/// <summary>F111-dedicated shared fixture: classify all 4248 Z₂³-homogeneous
/// k=4 N=4 pairs × 3 dephase letters once (~2-3 min PLINQ), and expose per-cell
/// per-class pair lists for the F111 enumeration tests. Stores actual
/// <see cref="PauliTerm"/> tuples (not counts) so the Pure-D Template Rule can be
/// checked per pair against the classifier's output.
///
/// <para>Why a separate fixture from <see cref="F106CountsFixtureN4K4"/> and
/// <see cref="HardCellYInversionPatternCountsFixture"/>: those fixtures store only
/// aggregated counts (4-tuple → int) and re-deriving pair tuples from them is
/// impossible. Extending them would change their type signatures and risk breaking
/// existing tests. Adding a parallel fixture is the minimal, side-effect-free
/// change. The PLINQ classification cost is the same (~2-3 min once); we just
/// keep the pair tuples in addition to grouping them by class.</para></summary>
public sealed class F111PerPairClassificationFixture
{
    private static readonly PauliLetter[] DephaseLetters =
        { PauliLetter.Z, PauliLetter.X, PauliLetter.Y };

    private readonly Lazy<Dictionary<(PauliLetter Dephase, (int A, int B) Klein, TrichotomyClass Cls),
        List<(PauliTerm Term1, PauliTerm Term2)>>> _classifiedPairs;

    public F111PerPairClassificationFixture()
    {
        _classifiedPairs = new Lazy<Dictionary<(PauliLetter, (int A, int B), TrichotomyClass),
            List<(PauliTerm, PauliTerm)>>>(
            ClassifyAndGroup, System.Threading.LazyThreadSafetyMode.ExecutionAndPublication);
    }

    /// <summary>All F87-hard pairs whose Klein index matches
    /// (<paramref name="bitA"/>, <paramref name="bitB"/>) under
    /// <paramref name="dephase"/>. Returns empty list if no pairs classified as hard
    /// in that cell. List is computed once via <see cref="Lazy{T}"/>.</summary>
    public IReadOnlyList<(PauliTerm Term1, PauliTerm Term2)> GetHardPairsInCell(
        PauliLetter dephase, int bitA, int bitB) =>
        Get(dephase, bitA, bitB, TrichotomyClass.Hard);

    /// <summary>All F87-soft pairs whose Klein index matches
    /// (<paramref name="bitA"/>, <paramref name="bitB"/>) under
    /// <paramref name="dephase"/>.</summary>
    public IReadOnlyList<(PauliTerm Term1, PauliTerm Term2)> GetSoftPairsInCell(
        PauliLetter dephase, int bitA, int bitB) =>
        Get(dephase, bitA, bitB, TrichotomyClass.Soft);

    /// <summary>All F87-truly pairs whose Klein index matches
    /// (<paramref name="bitA"/>, <paramref name="bitB"/>) under
    /// <paramref name="dephase"/>.</summary>
    public IReadOnlyList<(PauliTerm Term1, PauliTerm Term2)> GetTrulyPairsInCell(
        PauliLetter dephase, int bitA, int bitB) =>
        Get(dephase, bitA, bitB, TrichotomyClass.Truly);

    private IReadOnlyList<(PauliTerm, PauliTerm)> Get(
        PauliLetter dephase, int bitA, int bitB, TrichotomyClass cls) =>
        _classifiedPairs.Value.TryGetValue((dephase, (bitA, bitB), cls), out var list)
            ? list
            : Array.Empty<(PauliTerm, PauliTerm)>();

    private static Dictionary<(PauliLetter Dephase, (int A, int B) Klein, TrichotomyClass Cls),
        List<(PauliTerm, PauliTerm)>>
        ClassifyAndGroup()
    {
        var chain = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);

        // Parallelize over (pair × dephase letter). PauliPairTrichotomy.Classify is
        // pure; ChainSystem and PauliTerm are immutable records. Cap PLINQ degree
        // so MKL's all-core EVD inside Classify() has thread headroom (mirrors the
        // F110 fixture's WithDegreeOfParallelism call; avoids O(P²) oversubscription).
        var classifications = Z2HomogeneousKBodyEnumeration.Enumerate(4)
            .SelectMany(item => DephaseLetters.Select(d => (item, dephase: d)))
            .AsParallel()
            .WithDegreeOfParallelism(Math.Max(1, Environment.ProcessorCount / 4))
            .Select(x =>
            {
                var templates = new[] { x.item.Term1, x.item.Term2 };
                var cls = PauliPairTrichotomy.Classify(chain, templates, dephaseLetter: x.dephase);
                return (x.dephase, x.item.Klein, x.item.Term1, x.item.Term2, Cls: cls);
            })
            .ToList();

        var grouped = new Dictionary<(PauliLetter, (int A, int B), TrichotomyClass),
            List<(PauliTerm, PauliTerm)>>();
        foreach (var c in classifications)
        {
            var key = (c.dephase, c.Klein, c.Cls);
            if (!grouped.TryGetValue(key, out var list))
            {
                list = new List<(PauliTerm, PauliTerm)>();
                grouped[key] = list;
            }
            list.Add((c.Term1, c.Term2));
        }
        return grouped;
    }
}
