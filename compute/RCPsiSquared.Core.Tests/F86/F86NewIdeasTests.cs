using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Resonance;

using RCPsiSquared.Core.Knowledge;
namespace RCPsiSquared.Core.Tests.F86;

/// <summary>Tests for the two new F86 findings discovered 2026-05-03 in the OOP-driven
/// exploration: σ_0 chromaticity scaling and F71 mirror invariance of per-bond Q_peak.
/// </summary>
public class F86NewIdeasTests
{
    [Fact]
    public void SigmaZeroChromaticityScaling_AsymptoteFormula_IsTwoSqrtTwoCMinusOne()
    {
        Assert.Equal(2.0 * Math.Sqrt(2.0), SigmaZeroChromaticityScaling.Asymptote(c: 2), 12);
        Assert.Equal(4.0, SigmaZeroChromaticityScaling.Asymptote(c: 3), 12);
        Assert.Equal(2.0 * Math.Sqrt(6.0), SigmaZeroChromaticityScaling.Asymptote(c: 4), 12);
        // Q_EP asymptote
        Assert.Equal(1.0 / Math.Sqrt(2.0), SigmaZeroChromaticityScaling.QEpAsymptote(c: 2), 12);
        Assert.Equal(0.5, SigmaZeroChromaticityScaling.QEpAsymptote(c: 3), 12);
        Assert.Equal(1.0 / Math.Sqrt(6.0), SigmaZeroChromaticityScaling.QEpAsymptote(c: 4), 12);
    }

    [Fact]
    public void SigmaZeroChromaticityScaling_C2N7_HitsAsymptoteWithinTenMicro()
    {
        // Sweet-spot observation: c=2 N=7 σ_0 = 2√2 to within 10⁻⁵.
        var w = new SigmaZeroScalingWitness(chromaticity: 2, n: 7, gammaZero: 0.05);
        Assert.InRange(w.NormalisedRatio, 2.0 - 1e-5, 2.0 + 1e-5);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void SigmaZeroChromaticityScaling_NormalisedRatio_GrowsMonotonicallyWithN(int c)
    {
        // Within each chromaticity, σ_0/√(2(c−1)) grows monotonically with N. The 2.0
        // value is a TRAJECTORY CROSSING (sweet-spot at c=2 N=7, bit-exact to ~10⁻¹⁵),
        // NOT an asymptote: σ_0/√(2·1) at c=2 keeps climbing past 2.0 to 2.014 at N=9
        // and ~2.020 at N=11 (verified 2026-05-08 in the σ_0 bridge sweep). So the
        // monotone-growth claim survives, but the upper-bound envelope must now include
        // values above 2.0; see the class-level summary in SigmaZeroChromaticityScaling
        // for the full retraction context (c=2 case is empirically the only one tested
        // past the 2.0 crossing within the current N range; c=3 and c=4 may have their
        // own sweet-spot N_c* > 8).
        int[] Ns = c switch
        {
            2 => new[] { 5, 6, 7, 8, 9 },
            3 => new[] { 5, 6, 7, 8 },
            4 => new[] { 7, 8 },
            _ => throw new ArgumentException(),
        };

        var ratios = Ns
            .Select(N => new SigmaZeroScalingWitness(c, N, 0.05).NormalisedRatio)
            .ToArray();

        for (int i = 0; i < ratios.Length - 1; i++)
            Assert.True(ratios[i] < ratios[i + 1],
                $"σ_0/√(2(c−1)) should grow monotonically with N at c={c}: " +
                $"N={Ns[i]}→{ratios[i]:F4}, N={Ns[i + 1]}→{ratios[i + 1]:F4}");

        // Envelope: c=2 ratio at N=9 is 2.014 (above 2.0, the post-crossing region);
        // c=3 and c=4 still climb from below. Range covers both pre- and post-crossing.
        Assert.All(ratios, r => Assert.InRange(r, 1.65, 2.05));
    }

    [Theory]
    [InlineData(2, 5)]
    [InlineData(2, 6)]
    [InlineData(3, 5)]
    [InlineData(3, 6)]
    public void F71MirrorInvariance_PerBondQPeak_BitExactSymmetricUnderBondMirror(int c, int N)
    {
        // For every bond pair (b, N−2−b), Q_peak(b) = Q_peak(N−2−b) bit-exactly. This is
        // the F71 spatial-mirror algebraic identity applied to F86's per-bond observable.
        int n = c - 1;
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        // Wide grid so Endpoint Q_peak (~2.5) lands inside.
        var qGrid = ResonanceScan.LinearQGrid(0.20, 6.00, 60);
        var curve = new ResonanceScan(block).ComputeKCurve(qGrid);

        double maxDev = F71MirrorInvariance.MaxMirrorDeviation(curve);
        Assert.True(maxDev < 1e-10,
            $"F71 mirror invariance violated at c={c} N={N}: max |Q_peak(b) − Q_peak(N−2−b)| = {maxDev:E3}");
    }

    [Fact]
    public void F86KnowledgeBase_HasNewSigma0ScalingAndF71MirrorClaims()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);

        Assert.NotNull(kb.Sigma0Scaling);
        Assert.NotNull(kb.F71Mirror);
        // Tier was demoted from Tier1Candidate to Tier2Empirical 2026-05-08 after the σ_0
        // bridge sweep falsified the asymptote reading: 2√(2(c−1)) is a trajectory crossing
        // at sweet-spot N (=7 at c=2, bit-exact), not an N → ∞ limit.
        Assert.Equal(Tier.Tier2Empirical, kb.Sigma0Scaling.Tier);
        Assert.Equal(Tier.Tier1Derived, kb.F71Mirror.Tier);
    }

    [Fact]
    public void PerF71OrbitObservation_HasNineWitnessCases_AndFlagsCInversion()
    {
        var obs = new PerF71OrbitObservation();
        Assert.Equal(Tier.Tier2Empirical, obs.Tier);
        Assert.Equal(9, obs.Witnesses.Count);

        // c=2 N=6 has central < flanking
        var c2N6 = obs.Witnesses.First(w => w.Chromaticity == 2 && w.N == 6);
        Assert.Equal(2, c2N6.CentralIndex);
        Assert.True(c2N6.QPeakPerOrbit[2] < c2N6.QPeakPerOrbit[1],
            $"c=2 N=6: central {c2N6.QPeakPerOrbit[2]} should be < flanking {c2N6.QPeakPerOrbit[1]}");

        // c=3 N=6 has central > flanking — the inversion
        var c3N6 = obs.Witnesses.First(w => w.Chromaticity == 3 && w.N == 6);
        Assert.Equal(2, c3N6.CentralIndex);
        Assert.True(c3N6.QPeakPerOrbit[2] > c3N6.QPeakPerOrbit[1],
            $"c=3 N=6: central {c3N6.QPeakPerOrbit[2]} should be > flanking {c3N6.QPeakPerOrbit[1]}");
    }
}
