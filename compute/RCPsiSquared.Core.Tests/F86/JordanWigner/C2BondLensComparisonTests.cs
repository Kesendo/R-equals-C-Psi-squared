using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86.JordanWigner;

public class C2BondLensComparisonTests
{
    private readonly ITestOutputHelper _out;

    public C2BondLensComparisonTests(ITestOutputHelper output) => _out = output;

    // Both lenses are F71-mirror-invariant individually. The combined witness preserves
    // this; mirror-paired bonds have identical (EP-Q_peak, K-Q_peak, ...) pairs.
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void F71MirrorInvariance_BothLenses_BitExact(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var cmp = C2BondLensComparison.Build(block);
        int numBonds = cmp.Bonds.Count;
        for (int b = 0; b < numBonds; b++)
        {
            int mirror = numBonds - 1 - b;
            if (mirror <= b) continue;
            var bA = cmp.Bonds[b];
            var bB = cmp.Bonds[mirror];
            Assert.Equal(bA.EpQPeak, bB.EpQPeak, precision: 8);
            Assert.Equal(bA.EpKMax, bB.EpKMax, precision: 6);
            Assert.Equal(bA.KQPeak, bB.KQPeak, precision: 6);
            Assert.Equal(bA.KKMax, bB.KKMax, precision: 6);
        }
    }

    // **Cross-lens divergence finding (verified N=5..7):** the EP-Q_peak (sharp ‖xB(Q)‖_F
    // peak) and the K-Q_peak (broad K_b(Q, t_peak) peak) differ both in scale (the EP-Q_peak
    // tracks bond's effective Q_EP_b, K-Q_peak is broader t-averaged structure) and in
    // bond-class ordering. At N=7:
    //   Endpoint: EP-Q_peak = 1.84, K-Q_peak = 2.53 (K > EP)
    //   Innermost: EP-Q_peak = 1.84, K-Q_peak = 1.49 (K < EP)
    //   Flanking b=1: EP-Q_peak = 1.84, K-Q_peak = 4.00 (escaping the Q-grid upper bound!)
    // The flanking-Interior orbit-escape that becomes loud at N≥9 (per memory
    // project_q_peak_ep_structure.md) is already detectable in the K-Lens at N=7. EP-Lens
    // does NOT see this escape because it tracks Q_EP_b which stays bounded.
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    public void CrossLens_KAndEpQPeaks_AreNotProportional(int N)
    {
        var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
        var cmp = C2BondLensComparison.Build(block);
        // Q_peak ratio (K / EP) is NOT a constant across bonds — that's the cross-lens
        // divergence. We check that at least one bond has ratio significantly different
        // from another bond (i.e. spread of ratios > 0.05).
        var ratios = cmp.Bonds.Select(b => b.QPeakRatio_KOverEp).ToArray();
        double spread = ratios.Max() - ratios.Min();
        Assert.True(spread > 0.05,
            $"N={N}: Q_peak ratio spread (K/EP) = {spread:F4} should be > 0.05 (cross-lens divergence)");
    }

    [Fact]
    public void Build_RejectsNonC2Block()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);
        Assert.Throws<ArgumentException>(() => C2BondLensComparison.Build(block));
    }

    [Fact]
    public void Tier_IsTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var cmp = C2BondLensComparison.Build(block);
        Assert.Equal(Tier.Tier2Verified, cmp.Tier);
    }

    [Fact]
    public void Anchor_References_PROOF_F86_QPEAK_AndDirectionB()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var cmp = C2BondLensComparison.Build(block);
        Assert.Contains("PROOF_F86_QPEAK", cmp.Anchor);
        Assert.Contains("Direction (b'')", cmp.Anchor);
    }

    private static BondLensComparisonWitness? ClosestToCenter(IReadOnlyList<BondLensComparisonWitness> bonds, int N)
    {
        double center = (N - 2) / 2.0;
        BondLensComparisonWitness? best = null;
        double bestDist = double.PositiveInfinity;
        foreach (var b in bonds)
        {
            if (b.BondClass != BondClass.Interior) continue;
            double d = Math.Abs(b.Bond - center);
            if (d < bestDist) { bestDist = d; best = b; }
        }
        return best;
    }

    [Fact]
    public void Reconnaissance_EmitsEpVsKLens_AcrossN5To7()
    {
        // Emits both lenses' Q_peak/KMax/HWHM side-by-side. Expensive (C2HwhmRatio
        // builds a 153-pt ResonanceScan per N), so limited to N=5..7 for fast turnaround.
        _out.WriteLine("  N | b | class    | EP Q_pk | EP ‖xB‖ | EP H/Q | K Q_pk | K K_max | K H/Q | K-Q/EP-Q");
        _out.WriteLine("  --|---|----------|---------|---------|--------|---------|---------|--------|---------");
        foreach (int N in new[] { 5, 6, 7 })
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var cmp = C2BondLensComparison.Build(block);
            foreach (var b in cmp.Bonds)
            {
                string epHQ = b.EpHwhmLeftOverQPeak?.ToString("F4") ?? "n/a";
                string kHQ = b.KHwhmLeftOverQPeak?.ToString("F4") ?? "n/a";
                _out.WriteLine(
                    $"  {N} | {b.Bond} | {b.BondClass,-8} | {b.EpQPeak,7:F4} | {b.EpKMax,7:F4} | " +
                    $"{epHQ,6} | {b.KQPeak,7:F4} | {b.KKMax,7:G4} | {kHQ,6} | {b.QPeakRatio_KOverEp,7:F4}");
            }
            _out.WriteLine("");
        }
    }
}
