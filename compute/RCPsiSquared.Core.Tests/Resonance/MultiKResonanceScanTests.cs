using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Resonance;

public class MultiKResonanceScanTests(ITestOutputHelper output)
{
    [Fact]
    public void MultiK_C2_N5_MatchesFourModeScan()
    {
        // Sanity check: for c=2, MultiK should give the same K-curve as FourMode (same span).
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var multiCurve = new MultiKResonanceScan(block).ComputeKCurve();
        var fourCurve = new FourModeResonanceScan(block).ComputeKCurve();

        var multiInt = multiCurve.Peak(BondClass.Interior);
        var fourInt = fourCurve.Peak(BondClass.Interior);
        Assert.True(Math.Abs(multiInt.QPeak - fourInt.QPeak) < 0.01,
            $"c=2 N=5 Interior: MultiK={multiInt.QPeak:F4}, Four={fourInt.QPeak:F4}");
    }

    [Theory]
    [InlineData(5, 2)]  // c=3 N=5
    [InlineData(6, 2)]  // c=3 N=6
    [InlineData(7, 2)]  // c=3 N=7
    public void MultiK_C3_KMaxIsExactlyZero_NaiveExtensionFails(int N, int n)
    {
        // PROOF Item 2 verification: does the 7-mode multi-k effective reproduce the c=3
        // K-curve? Answer: NO — K_max collapses to exactly 0.
        //
        // Structural diagnosis: Gram-Schmidt orthogonalisation of the SVD-top vectors
        // against the channel-uniform vectors (|c_1⟩, |c_3⟩, |c_5⟩ at c=3) pushes them
        // into the CU-complement. Because M_H respects the CU/CU-complement decomposition
        // (channel-uniform-diagonal property of M_H_total, F73 generalisation), the probe
        // (which lives entirely in CU span) is uncoupled from the GS-modified SVD modes.
        // ∂ρ/∂J_b cannot move ρ out of CU → K = 0 identically.
        //
        // This is a real negative result for Item 2: the naïve "one quartet per k" picture
        // doesn't extend to c ≥ 3 with Gram-Schmidt. A correct effective model for c ≥ 3
        // would need a non-orthogonal frame that preserves CU ↔ SVD coupling, or a
        // different choice of the c−1 quartets that maintains coupling under projection.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        Assert.Equal(3, block.C);

        var multiCurve = new MultiKResonanceScan(block).ComputeKCurve();
        var multiInt = multiCurve.Peak(BondClass.Interior);
        var multiEnd = multiCurve.Peak(BondClass.Endpoint);

        output.WriteLine($"c=3 N={N}: multi-k Interior K_max={multiInt.KMax:E2}, Endpoint K_max={multiEnd.KMax:E2}");

        Assert.True(multiInt.KMax < 1e-10,
            $"c=3 N={N} Interior K_max should be ~0 (naive multi-k extension fails); got {multiInt.KMax:E3}");
        Assert.True(multiEnd.KMax < 1e-10,
            $"c=3 N={N} Endpoint K_max should be ~0 (naive multi-k extension fails); got {multiEnd.KMax:E3}");
    }
}
