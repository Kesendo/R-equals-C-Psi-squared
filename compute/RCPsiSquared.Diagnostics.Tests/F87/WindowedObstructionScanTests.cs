using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The shape of the windowed F87 hardness obstruction, by (k, N). The §7.5/§7.6 proof made
/// soft ⟺ bipartite a theorem for any k, leaving the SHAPE of the minimal odd 𝔽₂-relation as the
/// last k-specific piece. <see cref="WindowedObstructionScan"/> computes it in pure GF(2) bit
/// arithmetic (no Hamiltonian), reaching body counts the Liouvillian route cannot.
///
/// <para>Finding: the obstruction is a triangle at the smallest window count N = k+1 (forced, since
/// only 2 windows give |S| ≤ 4 and any odd relation among ≤ 4 masks has size 3), and for k = 3 it
/// stays a triangle even with many windows; but for k ≥ 4 with more windows a size-5 odd cycle
/// appears, so the obstruction shape is a genuine (k, N) family, NOT a uniform triangle. The
/// soft ⟺ bipartite derivation is k-agnostic regardless of the shape.</para></summary>
public class WindowedObstructionScanTests
{
    private readonly ITestOutputHelper _out;
    public WindowedObstructionScanTests(ITestOutputHelper o) => _out = o;

    [Fact]
    public void Scan_K3N4_ReproducesF103Anchor()
    {
        var r = WindowedObstructionScan.Scan(k: 3, n: 4);
        Assert.Equal(42, r.Pairs);
        Assert.Equal(16, r.Hard);
        Assert.Equal(26, r.Soft);
        Assert.True(r.ObstructionIsAlwaysTriangle);
        Assert.Equal(16, r.MinOddCycleSizes[3]);
    }

    // Pair/hard counts cross-checked against the Python scout simulations/_f87_oddcycle_kscaling.py.
    [Theory]
    [InlineData(3, 4, 42, 16)]
    [InlineData(4, 5, 828, 192)]
    [InlineData(5, 6, 14520, 1792)]
    public void Scan_PairAndHardCounts_MatchPythonScout(int k, int n, int pairs, int hard)
    {
        var r = WindowedObstructionScan.Scan(k, n);
        Assert.Equal(pairs, r.Pairs);
        Assert.Equal(hard, r.Hard);
    }

    /// <summary>The pure-GF(2) mask scan agrees with the canonical H-based bipartite classifier
    /// (<see cref="BipartiteChirality"/>) on every windowed k=3, N=4 pair: MinOddCycle == 0 ⟺
    /// bipartite ⟺ soft. This grounds the fast scan against the Hamiltonian route.</summary>
    [Fact]
    public void PureGF2Scan_AgreesWithHBasedBipartite_K3N4()
    {
        const int k = 3, n = 4;
        var chain = new ChainSystem(N: n, J: 1.0, GammaZero: 0.05);
        var withLetters = EnumerateCellTermsWithLetters(k);

        int checked_ = 0;
        for (int a = 0; a < withLetters.Count; a++)
            for (int b = a; b < withLetters.Count; b++)
            {
                var (letters1, mask1, yp1) = withLetters[a];
                var (letters2, mask2, yp2) = withLetters[b];
                if (yp1 != yp2) continue;

                var s = new HashSet<ulong>();
                for (int w = 0; w <= n - k; w++) { s.Add(mask1 << w); s.Add(mask2 << w); }
                bool gf2Bipartite = WindowedObstructionScan.MinOddCycle(new List<ulong>(s)) == 0;

                var terms = new List<PauliTerm> { new(letters1, Complex.One), new(letters2, Complex.One) };
                bool hBipartite = BipartiteChirality.Classify(chain, terms, PauliLetter.Z).IsBipartite;

                Assert.Equal(hBipartite, gf2Bipartite);
                checked_++;
            }
        Assert.Equal(42, checked_);
    }

    /// <summary>At the smallest window count N = k+1, the obstruction is FORCED to be a triangle:
    /// two windows give each term 2 masks, so |S| ≤ 4, and the only odd-cardinality subset that can
    /// XOR to 0 has size 3 (size 1 is impossible, size ≥ 5 needs ≥ 5 masks). Verified k = 3..8.</summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Obstruction_AtSmallestWindow_IsTriangle(int k)
    {
        var r = WindowedObstructionScan.Scan(k, k + 1);
        Assert.True(r.Hard > 0);
        Assert.True(r.ObstructionIsAlwaysTriangle, $"k={k}: sizes {string.Join(",", r.MinOddCycleSizes)}");
    }

    /// <summary>With more windows (N &gt; k+1) the triangle is no longer forced, and for k ≥ 4 a
    /// size-5 odd cycle genuinely appears: the obstruction is a (k, N) family, not a uniform
    /// triangle. (k = 3 happens to stay a triangle even with many windows; that is a k = 3 special
    /// case, not the general rule.)</summary>
    [Fact]
    public void Obstruction_GrowsBeyondTriangle_WithMoreWindows()
    {
        var r = WindowedObstructionScan.Scan(k: 4, n: 7);   // 4 windows, |S| up to 8
        Assert.False(r.ObstructionIsAlwaysTriangle);
        Assert.True(r.MinOddCycleSizes.ContainsKey(5), "expected a size-5 odd cycle at k=4, N=7");
        Assert.True(r.MinOddCycleSizes.Keys.Max() > 3, "obstruction exceeds a triangle");

        // k = 3 is the special case that stays a triangle even with many windows.
        Assert.True(WindowedObstructionScan.Scan(k: 3, n: 10).ObstructionIsAlwaysTriangle);
    }

    /// <summary>The obstruction-size law: the maximal minimal-odd-cycle over the hard pairs equals
    /// min(2W-1, 2k-3), with W = N-k+1 windows. The 2W-1 leg is the window count (|S| ≤ 2W masks, so
    /// an odd subset has ≤ 2W-1 of them); the 2k-3 leg is the body count (in the GF(2)[x] view a mask
    /// is a polynomial, even-popcount so divisible by 1+x, so the gcd-quotients p_i/g have degree ≤ k-2
    /// hence popcount ≤ k-1, and an odd coprime pair maxes at (k-1)+(k-2) = 2k-3). Below saturation the
    /// window leg binds; past W = k-1 the body leg binds. k=3 gives 2k-3 = 3, the always-triangle case.</summary>
    [Theory]
    [InlineData(3, 4)] [InlineData(3, 6)] [InlineData(3, 9)]
    [InlineData(4, 5)] [InlineData(4, 6)] [InlineData(4, 7)] [InlineData(4, 9)]
    [InlineData(5, 6)] [InlineData(5, 7)] [InlineData(5, 8)] [InlineData(5, 10)]
    [InlineData(6, 7)] [InlineData(6, 8)] [InlineData(6, 9)] [InlineData(6, 10)] [InlineData(6, 11)]
    public void SizeLaw_MaxIsMinOf2WMinus1And2kMinus3(int k, int n)
    {
        int w = n - k + 1;
        int expected = Math.Min(2 * w - 1, 2 * k - 3);
        var r = WindowedObstructionScan.Scan(k, n);
        Assert.True(r.Hard > 0);
        Assert.Equal(expected, r.MinOddCycleSizes.Keys.Max());
        Assert.True(r.MinOddCycleSizes.Keys.All(s => s <= 2 * k - 3), "no obstruction exceeds 2k-3");
    }

    /// <summary>At saturation (N = 2k, W = k+1 &gt; k-1) a pair is hard ⟺ its two window-masks have
    /// different (1+x)-adic valuations: hard ⟺ v(p1) ≠ v(p2). This is the non-bipartite criterion in
    /// GF(2)[x] language: an odd relation q_A p1 = q_B p2 exists iff popcount(p1/g)+popcount(p2/g) is
    /// odd, which holds iff the valuations differ (equal valuations force every relation to even size).
    /// Cross-checked here against the actual minimal-odd-cycle search on every saturated pair.</summary>
    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)]
    public void HardnessCriterion_IsValuationDifference_AtSaturation(int k)
    {
        int n = 2 * k;
        var terms = WindowedObstructionScan.CellTerms(k);
        int checkedPairs = 0;
        for (int a = 0; a < terms.Count; a++)
            for (int b = a; b < terms.Count; b++)
            {
                if (terms[a].YParity != terms[b].YParity) continue;
                var set = new HashSet<ulong>();
                for (int w = 0; w <= n - k; w++) { set.Add(terms[a].WindowMask << w); set.Add(terms[b].WindowMask << w); }
                bool hard = WindowedObstructionScan.MinOddCycle(new List<ulong>(set)) > 0;
                int v1 = WindowedObstructionScan.ValuationAtOnePlusX(terms[a].WindowMask);
                int v2 = WindowedObstructionScan.ValuationAtOnePlusX(terms[b].WindowMask);
                Assert.Equal(v1 != v2, hard);
                checkedPairs++;
            }
        Assert.True(checkedPairs > 0);
    }

    /// <summary>The extremal family saturating the body bound: masks p1 = x + x^{k-1} and
    /// p2 = 1 + x^{k-1} (the terms I·X·I…I·Y and X·I…I·Y) have gcd (1+x) and quotients of popcount
    /// k-2 and k-1, so the gcd-formula size is (k-2)+(k-1) = 2k-3. Pure valuation/gcd arithmetic, so
    /// this scales far past the exponential cycle-search range, confirming the bound is achieved for
    /// every k. (Their (1+x)-valuations differ, so the pair is genuinely hard.)</summary>
    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)]
    [InlineData(8)] [InlineData(10)] [InlineData(12)] [InlineData(16)] [InlineData(20)]
    public void ExtremalFamily_AchievesBodyBound(int k)
    {
        ulong p1 = (1UL << 1) | (1UL << (k - 1));   // x + x^{k-1}
        ulong p2 = 1UL | (1UL << (k - 1));          // 1 + x^{k-1}
        Assert.Equal(2 * k - 3, WindowedObstructionScan.GcdFormulaSize(p1, p2));
        Assert.NotEqual(WindowedObstructionScan.ValuationAtOnePlusX(p1),
                        WindowedObstructionScan.ValuationAtOnePlusX(p2));
    }

    /// <summary>Dumps the obstruction-size distribution over a (k, N) grid (run with a detailed
    /// console logger to read it). Documents how the shape grows with the window count.</summary>
    [Fact]
    public void DumpObstructionGrid()
    {
        foreach (int k in new[] { 3, 4, 5, 6 })
            foreach (int n in Enumerable.Range(k + 1, 6))
            {
                var r = WindowedObstructionScan.Scan(k, n);
                string sizes = string.Join(", ", r.MinOddCycleSizes.OrderBy(x => x.Key)
                    .Select(x => $"{x.Key}:{x.Value}"));
                _out.WriteLine($"k={k} N={n} (windows={n - k + 1}): pairs={r.Pairs} hard={r.Hard} " +
                               $"obstruction-sizes [{sizes}]");
            }
    }

    private static List<(PauliLetter[] letters, ulong mask, int yPar)> EnumerateCellTermsWithLetters(int k)
    {
        var outp = new List<(PauliLetter[], ulong, int)>();
        var lut = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        long total = 1L << (2 * k);
        for (long code = 0; code < total; code++)
        {
            int na = 0, nb = 0, ny = 0;
            ulong mask = 0;
            var letters = new PauliLetter[k];
            for (int i = 0; i < k; i++)
            {
                int letter = (int)((code >> (2 * i)) & 3);
                letters[i] = lut[letter];
                bool bitA = letter == 1 || letter == 2;
                bool bitB = letter == 2 || letter == 3;
                if (bitA) { na++; mask |= 1UL << i; }
                if (bitB) nb++;
                if (letter == 2) ny++;
            }
            if (na % 2 == 0 && nb % 2 == 1 && na >= 2)
                outp.Add((letters, mask, ny & 1));
        }
        return outp;
    }
}
