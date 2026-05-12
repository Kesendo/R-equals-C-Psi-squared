using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>BlockSpectrumPerformanceWitness (Tier 1 derived; 2026-05-11): empirical witness
/// at N=5, 6 (small enough to run full-L eig in test budget) that
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrum"/> produces multiset-equal spectrum
/// to full <c>L.Evd()</c> AND runs at least <see cref="ExpectedSpeedupAtN5"/>× faster
/// wall-clock time. Recorded reference timings (Windows 11 release-mode, MathNet/MKL
/// provider):
///
/// <list type="bullet">
///   <item>N=5: block ≈ 0.05 s vs full ≈ 0.4 s (~8×)</item>
///   <item>N=6: block ≈ 1.0 s vs full ≈ 10 s (~10×)</item>
/// </list>
///
/// <para>Tests assert: multiset-equal exactly at N=5; block-time/full-time
/// ratio &lt; <c>1/ExpectedSpeedupAtN5</c> (i.e. at least
/// <see cref="ExpectedSpeedupAtN5"/>× speedup) on min-of-3 measurement after warm-up.
/// The conservative <see cref="ExpectedSpeedupAtN5"/> = 2.0 lower bound absorbs CI
/// machine-variance well below the observed ~8× headroom.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/BlockSpectrum/LiouvillianBlockSpectrum.cs</c>
/// (parent claim, structural correctness),
/// <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/BlockSpectrumPerformanceWitnessTests.cs</c>
/// (timing assertion at N=5).</para></summary>
public sealed class BlockSpectrumPerformanceWitness : Claim
{
    private readonly LiouvillianBlockSpectrum _blockSpectrum;

    /// <summary>Conservative lower bound on the N=5 block/full speedup. Observed ~8× on
    /// the reference machine; the test asserts >= 2× to absorb CI variance without false
    /// positives on shared CI hardware.</summary>
    public const double ExpectedSpeedupAtN5 = 2.0;

    public BlockSpectrumPerformanceWitness(LiouvillianBlockSpectrum blockSpectrum)
        : base("BlockSpectrumPerformanceWitness: per-block eig wall-time at N=5,6 is at least 2× faster than full-L eig with multiset-equal spectrum; reference timings 0.05 s vs 0.4 s (N=5) and 1.0 s vs 10 s (N=6).",
               Tier.Tier1Derived,
               "LiouvillianBlockSpectrum (parent, structural correctness) + BlockSpectrumPerformanceWitnessTests (N=5 timing assertion: min-of-3 after warm-up, block/full ratio < 1/ExpectedSpeedupAtN5)")
    {
        _blockSpectrum = blockSpectrum ?? throw new ArgumentNullException(nameof(blockSpectrum));
    }

    public override string DisplayName =>
        "BlockSpectrumPerformanceWitness: per-block eig is ≥ 2× faster than full-L eig at N=5,6 with multiset-equal spectrum";

    public override string Summary =>
        $"reference timings (Win11/MKL release): N=5 ≈ 0.05 s vs 0.4 s, N=6 ≈ 1.0 s vs 10 s; conservative test bound {ExpectedSpeedupAtN5}× ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent",
                summary: "LiouvillianBlockSpectrum (per-block eig structural correctness)");
            yield return new InspectableNode("witness-N5",
                summary: "block ≈ 0.05 s vs full ≈ 0.4 s (~8×) on Win11 release/MKL");
            yield return new InspectableNode("witness-N6",
                summary: "block ≈ 1.0 s vs full ≈ 10 s (~10×) on Win11 release/MKL");
            yield return new InspectableNode("test-bound",
                summary: $"asserted >= {ExpectedSpeedupAtN5}× speedup at N=5, multiset-equal spectrum to 1e-9");
            yield return new InspectableNode("CLI-smoke-N7",
                summary: "block-spectrum --N 7 --refine f71 in tens of seconds, peak < 1 GB");
            yield return new InspectableNode("CLI-smoke-N8",
                summary: "block-spectrum --N 8 --refine f71 (no full-L verify), peak well under 1 GB");
        }
    }
}
