using System.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F111 fast combinatorial tests: verify the Pure-D Template Rule's
/// counting bookkeeping (pure-D template count = 8 per diagonal cell at k=N=4,
/// and the 8/24 split → 36/192/300 unordered-pair decomposition) without any
/// F87 classification work. These tests are pure enumeration over 4^4 = 256
/// length-4 Pauli strings; the per-test cost is microseconds.
///
/// <para>Sister to the SLOW <see cref="HardCellPureDTemplateEnumerationTests"/>,
/// which couples these combinatorial buckets to the F87 trichotomy classifier.
/// Split out to keep the fast bookkeeping independent of the ~2-3 min PLINQ
/// fixture so default CI sees the combinatorial sanity checks even when
/// SLOW_F111 is filtered out.</para></summary>
public sealed class HardCellPureDTemplateCombinatorialTests
{
    [Theory]
    [InlineData(PauliLetter.Z)]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void PureDTemplateCount_IsEightPerDiagonalCellAtK4N4(PauliLetter dephase)
    {
        var diagCellTerms = EnumerateDiagonalCellTermsAtK4N4(dephase);
        var pureD = diagCellTerms
            .Where(t => HardCellPureDTemplate.IsPureDTemplate(t, dephase))
            .ToList();

        // 8 = 4 single-D placements (#D=1) + 4 triple-D placements (#D=3).
        // Even-#D placements live in the (0, 0) Klein cell (identity / all-D),
        // not the diagonal cell whose bit_b = 1 selects odd #D for D ∈ {Z, X, Y}.
        Assert.Equal(8, pureD.Count);
    }

    [Theory]
    [InlineData(PauliLetter.Z)]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void DecompositionPerCell_Matches_36_192_300(PauliLetter dephase)
    {
        // The "cell" in the F111 sense is the (Klein, y_par) Z₂³-bucket: pairs
        // in Z₂HomogeneousKBodyEnumeration share both Klein and y_par. The
        // diagonal-cell y_par equals y_par(dephase) by F110 Aspect B; this is
        // the bucket carrying the 228 hard pairs.
        int dominantYpar = dephase.BitA() & dephase.BitB();
        var bucketTerms = EnumerateDiagonalCellTermsAtK4N4(dephase)
            .Where(t => t.YParity == dominantYpar)
            .ToList();

        // 32 = 4^4 / 2^3 (one of eight Z₂³ buckets) on length-4 strings.
        Assert.Equal(32, bucketTerms.Count);

        var pureD = bucketTerms
            .Where(t => HardCellPureDTemplate.IsPureDTemplate(t, dephase))
            .ToList();
        var mixed = bucketTerms
            .Where(t => !HardCellPureDTemplate.IsPureDTemplate(t, dephase))
            .ToList();
        Assert.Equal(8, pureD.Count);
        Assert.Equal(24, mixed.Count);

        // Unordered pair counts including self-pairs: C(n+1, 2) = n(n+1)/2.
        int purePure = pureD.Count * (pureD.Count + 1) / 2;
        int pureMixed = pureD.Count * mixed.Count;
        int mixedMixed = mixed.Count * (mixed.Count + 1) / 2;
        Assert.Equal(36, purePure);
        Assert.Equal(192, pureMixed);
        Assert.Equal(300, mixedMixed);
        Assert.Equal(528, purePure + pureMixed + mixedMixed);
    }

    // ============================================================
    // Helpers
    // ============================================================

    private static readonly PauliLetter[] AllLetters =
        { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };

    /// <summary>Enumerate all 4^4 = 256 length-4 Pauli strings (no filters).</summary>
    internal static IEnumerable<PauliTerm> EnumerateLength4PauliTerms()
    {
        for (int idx = 0; idx < 256; idx++)
        {
            var seq = new PauliLetter[4];
            int x = idx;
            for (int j = 0; j < 4; j++)
            {
                seq[j] = AllLetters[x & 3];
                x >>= 2;
            }
            yield return new PauliTerm(seq, Complex.One);
        }
    }

    /// <summary>Filter to length-4 terms in the diagonal Klein cell of
    /// <paramref name="dephase"/>, i.e., Klein index equals
    /// (<paramref name="dephase"/>.BitA(), <paramref name="dephase"/>.BitB()).</summary>
    private static List<PauliTerm> EnumerateDiagonalCellTermsAtK4N4(PauliLetter dephase)
    {
        var diagKlein = (dephase.BitA(), dephase.BitB());
        return EnumerateLength4PauliTerms()
            .Where(t => t.KleinIndex == diagKlein)
            .ToList();
    }
}
