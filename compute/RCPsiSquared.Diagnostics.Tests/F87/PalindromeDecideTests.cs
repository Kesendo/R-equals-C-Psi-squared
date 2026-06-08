using System.Numerics;
using System.Collections.Generic;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class PalindromeDecideTests
{
    private static PauliTerm T(params PauliLetter[] ls) => new(ls, Complex.One);
    private const PauliLetter I = PauliLetter.I, X = PauliLetter.X, Y = PauliLetter.Y, Z = PauliLetter.Z;

    [Fact]
    public void CertifyHard_True_OnDiagonalCellHardOddCyclePair()
    {
        // XXZ + XZX: Klein (0,1), Mixed, the §7 hard-odd-cycle witness (BipartiteChirality StandardSet).
        Assert.True(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(X, X, Z), T(X, Z, X) }));
    }

    [Fact]
    public void CertifyHard_False_OnDiagonalCellSoftPair()
    {
        // XXZ + ZXX: Klein (0,1), Mixed, soft (bipartite); equal valuation => not hard.
        Assert.False(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(X, X, Z), T(Z, X, X) }));
    }

    [Fact]
    public void CertifyHard_False_OutOfScope()
    {
        Assert.False(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(X, X, Z) }));               // 1 term
        Assert.False(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(X, X, X), T(X, X, Y), T(Y, X, X) })); // 3 terms
        Assert.False(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(Z), T(Z) }));               // pure-diagonal (na<2), longitudinal field
    }

    [Fact]
    public void Decide_Hard_OnDiagonalCellHardPair_ExhibitsObstruction()
    {
        var d = PalindromeSoftCertifier.Decide(new[] { T(X, X, Z), T(X, Z, X) }, n: 4);
        Assert.Equal(PalindromeSoftCertifier.Decision.Hard, d.Verdict);
        Assert.Equal(PalindromeSoftCertifier.HardStrategy.DiagonalCellValuation, d.HardStrategy);
        Assert.Contains("valuation", d.Reason);          // the exhibited (1+x)-valuation obstruction
    }

    [Fact]
    public void Decide_Soft_OnXYModel()
    {
        var d = PalindromeSoftCertifier.Decide(new[] { T(X, X), T(Y, Y) }, n: 4);
        Assert.Equal(PalindromeSoftCertifier.Decision.Soft, d.Verdict);
        Assert.NotEqual(PalindromeSoftCertifier.SoftStrategy.None, d.SoftStrategy);
    }

    [Fact]
    public void Decide_Undetermined_OutOfScope()
    {
        // A 3-term frustrated hard set: soft strategies decline AND it is out of the F115 pair scope.
        var d = PalindromeSoftCertifier.Decide(new[] { T(X, X, X), T(X, X, Y), T(Y, X, X) }, n: 4);
        Assert.Equal(PalindromeSoftCertifier.Decision.Undetermined, d.Verdict);
    }

    private static ChainSystem MakeChainN4() => new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);

    [Fact]
    public void Decide_HardVerdict_MatchesEngineAndAuthority_OverDiagonalCell_N4()
    {
        var chain = MakeChainN4();
        int checkedPairs = 0, hardPairs = 0;
        foreach (var t0 in PalindromeHardSweepTests.DiagonalCellMixedTerms(3))
            foreach (var t1 in PalindromeHardSweepTests.DiagonalCellMixedTerms(3))
            {
                if (PalindromeHardSweepTests.YParity(t0) != PalindromeHardSweepTests.YParity(t1)) continue;
                var terms = new[] { t0, t1 };
                bool decideHard = PalindromeSoftCertifier.Decide(terms, 4).Verdict == PalindromeSoftCertifier.Decision.Hard;
                bool engineHard = WindowedObstructionScan.IsHardPair(
                    PalindromeHardSweepTests.XyMask(t0), PalindromeHardSweepTests.XyMask(t1));
                Assert.True(engineHard == decideHard,                                  // Decide-hard == engine
                    $"Decide/engine disagree for [{string.Concat(t0.Letters)},{string.Concat(t1.Letters)}]: engine={engineHard}, Decide={decideHard}");
                if (decideHard)
                {
                    hardPairs++;
                    Assert.Equal(TrichotomyClass.Hard,                                 // == authority hard
                        PauliPairTrichotomy.Classify(chain, terms, dephaseLetter: PauliLetter.Z));
                }
                checkedPairs++;
            }
        Assert.True(checkedPairs > 0, "the consistency loop enumerated no pairs");
        Assert.True(hardPairs > 0, "no Decide-Hard pair occurred, so the authority cross-check never ran (vacuous)");
    }

    [Fact]
    public void Certify_StandardModels_Unchanged()
    {
        // Certify keeps its soft-only contract and verdicts (the doc-table models).
        Assert.True(PalindromeSoftCertifier.Certify(new[] { T(X, X), T(Y, Y) }, 4).Certified);          // XY model -> soft
        Assert.True(PalindromeSoftCertifier.Certify(new[] { T(X, Y), T(Y, X) }, 4).Certified);          // DM XY+YX -> soft
        Assert.False(PalindromeSoftCertifier.Certify(new[] { T(X, X, Z), T(X, Z, X) }, 4).Certified);   // hard pair -> NOT soft-certified
    }
}
