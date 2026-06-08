using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Task-1 soundness gate for the F115 hard-certifier: the (1+x)-valuation engine
/// (WindowedObstructionScan.IsHardPair) must agree with the spectral authority
/// (PauliPairTrichotomy.Classify) on EVERY in-scope diagonal-cell Mixed pair. Soundness = engine-Hard
/// implies authority-Hard (0 false-positives); completeness-in-scope = authority-Hard implies engine-Hard.</summary>
public class PalindromeHardSweepTests
{
    /// <summary>All Klein-(0,1) Mixed strings of body k (na=#(X/Y) even and >=2, nb=#(Y/Z) odd).</summary>
    internal static List<PauliTerm> DiagonalCellMixedTerms(int k)
    {
        var letters = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        var terms = new List<PauliTerm>();
        long total = 1L << (2 * k);                       // 4^k
        for (long code = 0; code < total; code++)
        {
            var ls = new PauliLetter[k];
            int na = 0, nb = 0;
            for (int i = 0; i < k; i++)
            {
                var l = letters[(code >> (2 * i)) & 3];
                ls[i] = l;
                if (l == PauliLetter.X || l == PauliLetter.Y) na++;
                if (l == PauliLetter.Y || l == PauliLetter.Z) nb++;
            }
            if (na % 2 == 0 && na >= 2 && nb % 2 == 1)
                terms.Add(new PauliTerm(ls, System.Numerics.Complex.One));
        }
        return terms;
    }

    internal static ulong XyMask(PauliTerm t)
    {
        ulong m = 0;
        for (int i = 0; i < t.Letters.Count; i++)
            if (t.Letters[i] == PauliLetter.X || t.Letters[i] == PauliLetter.Y) m |= 1UL << i;
        return m;
    }

    internal static int YParity(PauliTerm t)
    {
        int ny = 0;
        for (int i = 0; i < t.Letters.Count; i++) if (t.Letters[i] == PauliLetter.Y) ny++;
        return ny & 1;
    }

    [Theory]
    [InlineData(2, 4)]
    [InlineData(3, 4)]
    [InlineData(2, 5)]
    [InlineData(3, 5)]
    public void Engine_AgreesWithAuthority_OnDiagonalCellPairs(int k, int n)
    {
        var chain = MakeChain(n);
        var cell = DiagonalCellMixedTerms(k);
        int checkedPairs = 0;
        for (int a = 0; a < cell.Count; a++)
            for (int b = a; b < cell.Count; b++)
            {
                if (YParity(cell[a]) != YParity(cell[b])) continue;   // y_par-homogeneous
                var terms = new[] { cell[a], cell[b] };
                bool engineHard = WindowedObstructionScan.IsHardPair(XyMask(cell[a]), XyMask(cell[b]));
                var actual = PauliPairTrichotomy.Classify(chain, terms, dephaseLetter: PauliLetter.Z);
                if (engineHard)
                    Assert.True(actual == TrichotomyClass.Hard,
                        $"false-positive: engine HARD but authority {actual} for [{Label(cell[a])},{Label(cell[b])}] N={n}");
                if (actual == TrichotomyClass.Hard)
                    Assert.True(engineHard,
                        $"in-scope miss: authority HARD but engine SOFT for [{Label(cell[a])},{Label(cell[b])}] N={n}");
                checkedPairs++;
            }
        Assert.True(checkedPairs > 0, "the sweep enumerated no pairs");
    }

    [Fact]
    public void Engine_AgreesWithAuthority_N6_Spotcheck()
    {
        // N=6 stability spot-check on two known diagonal-cell pairs (a full sweep at N=6 is needlessly slow).
        var chain = MakeChain(6);
        PauliTerm T(params PauliLetter[] ls) => new(ls, System.Numerics.Complex.One);
        var hard = new[] { T(PauliLetter.X, PauliLetter.X, PauliLetter.Z), T(PauliLetter.X, PauliLetter.Z, PauliLetter.X) };
        var soft = new[] { T(PauliLetter.X, PauliLetter.X, PauliLetter.Z), T(PauliLetter.Z, PauliLetter.X, PauliLetter.X) };
        Assert.True(WindowedObstructionScan.IsHardPair(XyMask(hard[0]), XyMask(hard[1])));
        Assert.Equal(TrichotomyClass.Hard, PauliPairTrichotomy.Classify(chain, hard, dephaseLetter: PauliLetter.Z));
        Assert.False(WindowedObstructionScan.IsHardPair(XyMask(soft[0]), XyMask(soft[1])));
        Assert.Equal(TrichotomyClass.Soft, PauliPairTrichotomy.Classify(chain, soft, dephaseLetter: PauliLetter.Z));
    }

    private static string Label(PauliTerm t) => string.Concat(t.Letters);

    private static ChainSystem MakeChain(int n) =>
        new ChainSystem(N: n, J: 1.0, GammaZero: 0.05);   // matches F87DiagonalCellBipartiteWitnessTests.MakeChainN4
}
