using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Guards <see cref="SectorEpProbe"/>: the core Multi-Sector Monodromy Census primitive
/// that, for a sector at a coupling q0, builds the RAW (p, q̃) block (Task 1's
/// <see cref="RCPsiSquared.Core.BlockSpectrum.SectorBlock"/>), finds the coalescing eigenvalue pair,
/// and classifies it defective-vs-diabolic via the trusted
/// <see cref="RCPsiSquared.Core.Numerics.EpCharacter"/>.
///
/// <para>The two q-values below are the review's decisive from-below result
/// (<c>_msm_ep_character.py</c>) on the N=4 (p=1, q̃=2) block: q≈0.658983 is the SEMISIMPLE DIABOLIC
/// (a (3q⁴+q²−1) root, SVD nullity 2, the pair merges at λ≈−4+1.318i and stays independent) and
/// q≈0.460212 is a DEFECTIVE EP (a P₁₀ real root, SVD nullity 1, a braid-carrying Jordan block).
/// Reproducing this in C# is the load-bearing validation that the Task-1 convention + EpCharacter
/// agree with the exact math.</para></summary>
public class SectorEpProbeTests
{
    [Theory]
    [InlineData(0.658983, 0.0, "Diabolic")]   // (3q^4+q^2-1) root: semisimple node, silent (linear split)
    [InlineData(0.460212, 0.0, "Defective")]  // P_10 real root: braid-carrying EP (√-split; refined below)
    public void Probe_ClassifiesOcticCoalescence(double qre, double qim, string kind)
    {
        // The DEFECTIVE EP splits as gap ∝ √(q−q*), so the 6-digit reference seed sits ~3e-3 off the true
        // gap floor. Refine that seed to the actual EP on the REAL q-axis before probing (the EP has its
        // λ on the Re=−4 line, so a 1D real-axis refine reaches it): the gap drops to ~2e-6. The DIABOLIC
        // node splits LINEARLY and is already at gap ≈ 4e-7 from its seed, so it needs no refine. The
        // refined q is computed here from our own minimization, never hardcoded.
        double qProbe = kind == "Defective" ? RefineToEpRealAxis(qre - 0.005, qre + 0.005) : qre;

        var r = SectorEpProbe.Probe(N: 4, p: 1, qTilde: 2, q0: new Complex(qProbe, qim));
        Assert.True(r.MinGap < 1e-3, $"expected a coalescence, gap={r.MinGap} at q={qProbe}");
        Assert.Equal(kind, r.Kind.ToString());
    }

    // A size-1 sector (the corners C(N,0)² = 1 at every N) has no eigenvalue pair. The census sweep will
    // hit these, so the probe must NOT index a second eigenvalue: it reports MinGap = +∞ / Normal (no
    // coalescence) instead of throwing. This pins that guard.
    [Fact]
    public void Probe_Size1Sector_ReportsNoCoalescence_WithoutThrowing()
    {
        var r = SectorEpProbe.Probe(N: 4, p: 0, qTilde: 0, q0: new Complex(2.0, 0.0));
        Assert.Equal(double.PositiveInfinity, r.MinGap);
        Assert.Equal(EpCharacter.EpKind.Normal, r.Kind);
    }

    /// <summary>Golden-section minimize the RAW block's min-eigenvalue-gap over real q in [lo, hi],
    /// returning the minimizing q (the actual EP). About a defective EP the coalescence gap is
    /// gap(q) ∝ √|q−q*|, a unimodal V with its cusp at q*, so golden-section converges to it; 50 iters
    /// shrink the 1e-2 bracket to ~4e-13, well past the point where the √-law puts the gap below 1e-3.
    /// The objective is the probe's own <see cref="SectorEpProbe.ProbeReading.MinGap"/> (reuse of the
    /// primitive, not a re-derivation).</summary>
    private static double RefineToEpRealAxis(double lo, double hi)
    {
        const double invPhi = 0.6180339887498949;   // 1/φ = (√5 − 1)/2
        static double Gap(double q) => SectorEpProbe.Probe(4, 1, 2, new Complex(q, 0.0)).MinGap;

        double a = lo, b = hi;
        double c = b - invPhi * (b - a), d = a + invPhi * (b - a);
        double fc = Gap(c), fd = Gap(d);
        for (int it = 0; it < 50; it++)
        {
            if (fc < fd) { b = d; d = c; fd = fc; c = b - invPhi * (b - a); fc = Gap(c); }
            else         { a = c; c = d; fc = fd; d = a + invPhi * (b - a); fd = Gap(d); }
        }
        return 0.5 * (a + b);
    }
}
