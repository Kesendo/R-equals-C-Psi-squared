using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The step-3 adjudication of the 2026-07-03 σ-scout's "NEW N=9 seed" q* = 1.4994 (recorded in
/// the OpenArcs ledger with gap 4.6e-5): it is NOT a defective seed. The AT-aware Riesz probe reads NO
/// defective cluster at the locus, and the near-coalescing pair's gap closes LINEARLY in |q − q*|
/// (ratio ≈ 2 per doubling; a defective EP closes as √). Verdict: a semisimple real-real (near-)crossing
/// in the R-even sector — exactly the count-change census's declared blind-spot class, so the census's
/// N=9 seed list (7 entries, RealDefectiveSeeds) was RIGHT to omit it. The scout's σ_min readings at this
/// locus stay valid as instrument anchors (λ = −4.3807 is an eigenvalue of (1,2) there regardless of
/// Jordan character; containment transports it), but the locus is not census-seed material.</summary>
public class ScoutSeedAdjudicationProbe
{
    [Fact]
    [Trait("Category", "SHELLCENSUS_ADJUDICATE")]
    public void N9_ScoutLocus_1_4994_IsSemisimpleCrossing_NotADefectiveSeed()
    {
        var q0 = new Complex(1.4994, 0);

        // (a) the AT-aware Riesz reader: no defective cluster anywhere in the (1,2) block at the locus
        var reading = SectorEpProbe.ProbeDefectiveAnywhere(9, 1, 2, q0);
        Assert.False(reading.HasDefective);
        Assert.True(reading.MultiMemberClusterCount >= 1);

        // (b) the pair lives in the R-even sector with the scout's recorded gap...
        double GapAt(double q, bool odd)
        {
            var (a, d) = WeightCoherenceBlock.BuildReflectionSectorColumnMajor(9, 1, 2, new Complex(q, 0), odd);
            var m = Matrix<Complex>.Build.Dense(d, d);
            for (int c = 0; c < d; c++)
                for (int r = 0; r < d; r++)
                    m[r, c] = a[(long)c * d + r];
            var ev = m.Evd().EigenValues.Enumerate()
                .Where(z => Math.Abs(z.Real - (-4.3807)) < 0.3).ToList();
            double minGap = double.PositiveInfinity;
            for (int i = 0; i < ev.Count; i++)
                for (int j = i + 1; j < ev.Count; j++)
                    minGap = Math.Min(minGap, (ev[i] - ev[j]).Magnitude);
            return minGap;
        }
        Assert.True(GapAt(1.4994, odd: false) < 1e-4);   // the scout's 4.6e-5
        Assert.True(GapAt(1.4994, odd: true) > 1e-1);    // R-odd far from coalescing

        // (c) LINEAR gap closure (semisimple crossing), not the defective √-law:
        // gap(2δ)/gap(δ) ≈ 2 for linear, ≈ √2 ≈ 1.41 for defective
        double g1 = GapAt(1.4994 + 0.001, odd: false);
        double g2 = GapAt(1.4994 + 0.002, odd: false);
        Assert.InRange(g2 / g1, 1.8, 2.2);
    }
}
