using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F115 closed-form hard-count ports (PROOF_F103 §7.8-§7.9), each asserted bit-exact
/// against direct enumeration over the even-popcount k-bit mask space (the diagonal-cell X/Y flip
/// patterns), using the same <see cref="WindowedObstructionScan"/> helpers the scan itself uses
/// (ValuationAtOnePlusX, GRestDegree, MinOddCycle). The oracles are the C# mirror of the verified
/// scouts simulations/_f87_hardcount_closedform.py, _f87_dlayer_count.py, _f87_size_second_layer.py:
/// the A203241 total, the d=0 base B(k), the d-layered count 2^(d-1)·B(k-d), the size-3 (triangle)
/// sub-count, and the d-layered obstruction-size cap 2k-3-2d.</summary>
public class WindowedHardnessCountClosedFormTests
{
    // ---- independent oracles: direct enumeration over the mask space ----

    private static long OracleTotalHard(int k)
    {
        var masks = WindowedObstructionScan.EvenPopcountMasks(k);
        long n = 0;
        for (int i = 0; i < masks.Count; i++)
            for (int j = i + 1; j < masks.Count; j++)
                if (WindowedObstructionScan.ValuationAtOnePlusX(masks[i])
                    != WindowedObstructionScan.ValuationAtOnePlusX(masks[j])) n++;
        return n;
    }

    private static Dictionary<int, long> OracleHardByD(int k)
    {
        var masks = WindowedObstructionScan.EvenPopcountMasks(k);
        var byD = new Dictionary<int, long>();
        for (int i = 0; i < masks.Count; i++)
            for (int j = i + 1; j < masks.Count; j++)
            {
                if (WindowedObstructionScan.ValuationAtOnePlusX(masks[i])
                    == WindowedObstructionScan.ValuationAtOnePlusX(masks[j])) continue;
                int d = WindowedObstructionScan.GRestDegree(masks[i], masks[j]);
                byD[d] = byD.GetValueOrDefault(d) + 1;
            }
        return byD;
    }

    private static List<ulong> WindowedEdgeSet(ulong p1, ulong p2, int k, int nSites)
    {
        var set = new HashSet<ulong>();
        for (int w = 0; w <= nSites - k; w++) { set.Add(p1 << w); set.Add(p2 << w); }
        return set.ToList();
    }

    private static long OracleTriangleCountAtN2k(int k)
    {
        var masks = WindowedObstructionScan.EvenPopcountMasks(k);
        int nSites = 2 * k;
        long n = 0;
        for (int i = 0; i < masks.Count; i++)
            for (int j = i + 1; j < masks.Count; j++)
            {
                if (WindowedObstructionScan.ValuationAtOnePlusX(masks[i])
                    == WindowedObstructionScan.ValuationAtOnePlusX(masks[j])) continue;
                if (WindowedObstructionScan.MinOddCycle(WindowedEdgeSet(masks[i], masks[j], k, nSites)) == 3) n++;
            }
        return n;
    }

    private static Dictionary<int, int> OracleMaxSizeByDAtN2k(int k)
    {
        var masks = WindowedObstructionScan.EvenPopcountMasks(k);
        int nSites = 2 * k;
        var cap = new Dictionary<int, int>();
        for (int i = 0; i < masks.Count; i++)
            for (int j = i + 1; j < masks.Count; j++)
            {
                if (WindowedObstructionScan.ValuationAtOnePlusX(masks[i])
                    == WindowedObstructionScan.ValuationAtOnePlusX(masks[j])) continue;
                int d = WindowedObstructionScan.GRestDegree(masks[i], masks[j]);
                int s = WindowedObstructionScan.MinOddCycle(WindowedEdgeSet(masks[i], masks[j], k, nSites));
                cap[d] = System.Math.Max(cap.GetValueOrDefault(d), s);
            }
        return cap;
    }

    // ---- size-3 floor (the MacWilliams-kernel floor, mirror of simulations/f87_size_cells.py) ----

    /// <summary>Oracle for c(D,3): reduced size-3 pairs are (monomial x^j, popcount-2 b) that are
    /// coprime, with max(deg a, deg b) = D. Coprime ⟺ j==0 or low-bit(b)==0 (PolyGcd checks it).</summary>
    private static long OracleTriangleReducedPairsByMaxDegree(int D)
    {
        long n = 0;
        for (int j = 0; j <= D; j++)
        {
            ulong a = 1UL << j;
            for (int p = 0; p <= D; p++)
                for (int q = p + 1; q <= D; q++)
                {
                    ulong b = (1UL << p) | (1UL << q);
                    if (System.Math.Max(j, q) != D) continue;
                    if (WindowedObstructionScan.PolyGcd(a, b) != 1) continue;
                    n++;
                }
        }
        return n;
    }

    // ---- helper sanity ----

    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)]
    public void EvenPopcountMasks_HasExpectedCount(int k)
    {
        // nonzero even-popcount k-bit masks: 2^(k-1) - 1 of them, all with even popcount.
        var masks = WindowedObstructionScan.EvenPopcountMasks(k);
        Assert.Equal((1 << (k - 1)) - 1, masks.Count);
        Assert.All(masks, m => Assert.True(System.Numerics.BitOperations.PopCount(m) % 2 == 0 && m != 0));
    }

    [Fact]
    public void GRestDegree_StripsOnePlusXFactors()
    {
        // gcd(x(1+x)=0b110, (1+x)^2=0b101) = (1+x); stripping (1+x) leaves a unit -> degree 0.
        Assert.Equal(0, WindowedObstructionScan.GRestDegree(0b110, 0b101));
        // (1+x)(1+x+x^2)=0b1001 and (1+x+x^2)=0b111 share the irreducible 1+x+x^2 beyond any (1+x):
        // gcd = 0b111 (odd popcount, no (1+x) factor), so g_rest degree = deg(1+x+x^2) = 2.
        Assert.Equal(2, WindowedObstructionScan.GRestDegree(0b1001, 0b111));
    }

    // ---- closed forms vs oracles ----

    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)]
    [InlineData(7)] [InlineData(8)] [InlineData(9)] [InlineData(10)]
    public void TotalHardCount_MatchesA203241(int k)
    {
        Assert.Equal(OracleTotalHard(k), WindowedObstructionScan.HardMaskPairCount(k));
    }

    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)] [InlineData(8)]
    public void BaseCountB_MatchesDZeroLayer(int k)
    {
        Assert.Equal(OracleHardByD(k).GetValueOrDefault(0), WindowedObstructionScan.HardCountBaseB(k));
    }

    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)] [InlineData(8)]
    public void DLayerCount_MatchesEnumeration_AndSumsToTotal(int k)
    {
        var byD = OracleHardByD(k);
        foreach (var (d, count) in byD)
            Assert.Equal(count, WindowedObstructionScan.HardCountByGRestDegree(k, d));
        Assert.Equal(WindowedObstructionScan.HardMaskPairCount(k), byD.Values.Sum());
    }

    [Theory]
    [InlineData(3)] [InlineData(4)] [InlineData(5)] [InlineData(6)] [InlineData(7)]
    public void TriangleSubCount_MatchesClosedForm(int k)
    {
        Assert.Equal(OracleTriangleCountAtN2k(k), WindowedObstructionScan.TriangleHardMaskCount(k));
    }

    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(4)]
    [InlineData(5)] [InlineData(6)] [InlineData(7)]
    public void TriangleReducedPairCountByMaxDegree_Is3DMinus1_AndMatchesOracle(int D)
    {
        Assert.Equal(3L * D - 1, WindowedObstructionScan.TriangleReducedPairCountByMaxDegree(D));
        Assert.Equal(OracleTriangleReducedPairsByMaxDegree(D),
                     WindowedObstructionScan.TriangleReducedPairCountByMaxDegree(D));
    }

    [Theory]
    [InlineData(4)] [InlineData(5)] [InlineData(6)]
    public void DLayeredMaxObstructionSize_Is2kMinus3Minus2d(int k)
    {
        var cap = OracleMaxSizeByDAtN2k(k);
        Assert.NotEmpty(cap);
        foreach (var (d, maxSize) in cap)
            Assert.Equal(WindowedObstructionScan.MaxObstructionSizeForGRestDegree(k, d), maxSize);
    }
}
