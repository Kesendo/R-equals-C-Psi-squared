using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;
using Verdict = RCPsiSquared.Diagnostics.F87.PalindromeMaskClassifier.Verdict;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The Liouvillian-free <see cref="PalindromeMaskClassifier"/> gives the correct soft/hard
/// verdict for bit_b-homogeneous hopping Hamiltonians (matching the trichotomy at N = 4), flags the
/// out-of-scope cases (mixed Klein cell / pure diagonal), and scales to large N with no Liouvillian.</summary>
public class PalindromeMaskClassifierTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> H(params string[] labels) => labels.Select(T).ToList();

    [Fact]
    public void Classify_MatchesTrichotomy_OnBitBHomogeneousHopping()
    {
        var chain = new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        var inScope = new[]
        {
            H("XX", "YY"),          // XY model: truly -> not hard
            H("XY", "YX"),          // soft
            H("XZ", "ZX"),          // soft
            H("XXI", "YYI"),        // truly -> not hard
            H("XYI", "YIX"),        // hard
            H("XIY", "IXY"),        // hard
            H("XYI", "YIX", "IXY"), // hard, three bit_b-homogeneous terms
        };
        foreach (var terms in inScope)
        {
            var spectral = PauliPairTrichotomy.Classify(chain, terms, dephaseLetter: PauliLetter.Z);
            var expected = spectral == TrichotomyClass.Hard ? Verdict.Hard : Verdict.Soft;
            Assert.Equal(expected, PalindromeMaskClassifier.Classify(terms, 4));
        }
    }

    [Fact]
    public void Classify_ReturnsOutOfScope_ForMixedCellOrPureDiagonal()
    {
        Assert.Equal(Verdict.OutOfScope, PalindromeMaskClassifier.Classify(H("ZI"), 4));              // pure diagonal
        Assert.Equal(Verdict.OutOfScope, PalindromeMaskClassifier.Classify(H("XX", "YY", "XY"), 4));  // mixed cell
        Assert.Equal(Verdict.OutOfScope, PalindromeMaskClassifier.Classify(H("XX", "XZ"), 4));        // mixed cell
    }

    [Fact]
    public void Classify_ScalesToLargeN_WithNoLiouvillian()
    {
        // The same hard / soft pairs classified at N = 40: pure GF(2) on 40-bit masks, no 2^40 matrix.
        Assert.Equal(Verdict.Hard, PalindromeMaskClassifier.Classify(H("XYI", "YIX"), 40));
        Assert.Equal(Verdict.Soft, PalindromeMaskClassifier.Classify(H("XY", "YX"), 40));
    }
}
