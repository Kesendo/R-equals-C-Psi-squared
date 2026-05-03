using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.Tests.Resonance;

/// <summary>F86 documentary tests for the 4-mode minimal effective model.
///
/// <para>The structural tests verify that the 4-mode scan produces a well-formed K-curve.
/// The numerical fingerprint tests record the empirical answer to F86's open question
/// (PROOF_F86_QPEAK "What's missing for full Tier 1" item 1):</para>
/// <list type="bullet">
///   <item>4-mode <b>Interior HWHM_left/Q_peak ≈ 0.74</b> across c=2 N=5..8 — within ~2%
///         of the universal 0.756 from the full block-L. The shape symmetry survives.</item>
///   <item>4-mode <b>Interior Q_peak ≈ 1.9–3.0</b> vs full block-L ≈ 1.48–1.60 — Q_peak
///         position is NOT reproduced; the 4-mode model misses ~50% of the eigenstructure
///         that determines absolute peak position.</item>
///   <item>4-mode <b>Endpoint Q_peak</b> exceeds the default Q ≤ 4 grid — Endpoint shape
///         not captured at all by the minimal effective.</item>
///   <item>4-mode <b>K_max</b> is 3–10× smaller than full block-L — the 4-mode subspace
///         carries only a fraction of the observable signal.</item>
/// </list>
/// <para>Conclusion: the 4-mode minimal effective is a structural witness (probe lives in
/// channel-uniform 2D, EP partners in SVD-top 2D, with cross-coupling carrying bond-position
/// information) but NOT sufficient to reproduce the universal shape numerically. More modes
/// are needed for Tier 1 promotion of the K_class(Q) = f_class(Q/Q_EP) statement.</para>
/// </summary>
public class FourModeResonanceScanTests
{
    [Fact]
    public void KCurve_DefaultGrids_ProducesShapeCompatibleWithFullScan()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var curve = new FourModeResonanceScan(block).ComputeKCurve();

        Assert.Equal(block.NumBonds, curve.NumBonds);
        Assert.Equal(ResonanceScan.DefaultQGrid().Length, curve.QGrid.Count);
        Assert.Equal(0.20, curve.QGrid[0], 12);
        Assert.Equal(4.00, curve.QGrid[^1], 12);
    }

    [Theory]
    [InlineData(5, 1, 1.85, 1.99, 0.72, 0.76)]
    [InlineData(6, 1, 2.85, 3.05, 0.73, 0.76)]
    [InlineData(7, 1, 2.95, 3.15, 0.73, 0.76)]
    [InlineData(8, 1, 2.95, 3.10, 0.73, 0.76)]
    public void Interior_FourMode_HasShiftedQPeakButCloseHwhmRatio(
        int N, int n, double qLo, double qHi, double hwhmLo, double hwhmHi)
    {
        // Documented finding: the 4-mode Interior HWHM/Q ratio is within ~2% of the
        // universal 0.756 (full block-L), but Q_peak itself is shifted significantly
        // higher. Together this is the algebraic witness that:
        //   (a) the bond-class shape symmetry partially lives in the 4-mode subspace, but
        //   (b) the absolute peak position does NOT — more modes are needed.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var curve = new FourModeResonanceScan(block).ComputeKCurve();

        var interior = curve.Peak(BondClass.Interior);
        Assert.InRange(interior.QPeak, qLo, qHi);
        Assert.NotNull(interior.HwhmLeftOverQPeak);
        Assert.InRange(interior.HwhmLeftOverQPeak!.Value, hwhmLo, hwhmHi);
    }

    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 1)]
    [InlineData(7, 1)]
    [InlineData(8, 1)]
    public void Endpoint_FourMode_QPeakBeyondDefaultGrid(int N, int n)
    {
        // Documented finding: the 4-mode Endpoint K-curve is still rising at Q = 4
        // (default grid edge). The peak is off-grid and the universal Endpoint
        // HWHM/Q ≈ 0.770 is NOT recovered by the 4-mode reduction. This is part of
        // the structural answer to PROOF_F86_QPEAK item 1.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var curve = new FourModeResonanceScan(block).ComputeKCurve();

        var endpoint = curve.Peak(BondClass.Endpoint);
        Assert.True(endpoint.QPeak >= 3.99,
            $"c=2 N={N} Endpoint 4-mode unexpectedly peaks at Q={endpoint.QPeak:F3} inside default grid; if true, " +
            "the 4-mode reduction may capture more Endpoint structure than measured 2026-05-02 — investigate.");
    }

    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 1)]
    [InlineData(7, 1)]
    [InlineData(8, 1)]
    public void Interior_FourMode_KMax_FractionOfFullBlockL(int N, int n)
    {
        // Documented finding: 4-mode K_max is a fraction of the full block-L K_max —
        // the 4-mode subspace carries only part of the K_CC_pr signal. We assert the
        // fraction is small (< 0.5) to lock in the "more modes needed" conclusion.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var fullCurve = new ResonanceScan(block).ComputeKCurve();
        var fourCurve = new FourModeResonanceScan(block).ComputeKCurve();

        double kFull = fullCurve.Peak(BondClass.Interior).KMax;
        double kFour = fourCurve.Peak(BondClass.Interior).KMax;
        double ratio = kFour / kFull;
        Assert.True(ratio > 0 && ratio < 0.5,
            $"c=2 N={N} Interior K_max ratio (4-mode / full) = {ratio:F3}; expected (0, 0.5).");
    }
}
