using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class DimensionSweepTests
{
    [Fact]
    public void Crossover_Marks_Invariant_While_Subspace_Rotates()
    {
        var axis = DimensionAxis.Crossover(N: 3, gamma: 0.5, thetaPoints: 13);
        var sweep = DimensionSweep.Compute(axis, slowCount: 16);
        Assert.True(sweep.MaxEigenvalueDriftAcrossTheta < 1e-7,
            $"marks should be invariant (similarity); drift {sweep.MaxEigenvalueDriftAcrossTheta:E2}");
        Assert.True(sweep.SubspaceRotation.Max() > 1e-3,
            $"the in-between (slow subspace) should rotate; max {sweep.SubspaceRotation.Max():E2}");
    }

    [Fact]
    public void Polarity_Follows_F99_Ladder_SinSquaredOverTwo()
    {
        // thetaPoints=13 gives the grid θ_k = (π/2)·k/12, so π/4 lands on k=6 and π/2 on k=12.
        // Polarity[k] = sin²(θ_k)/2 reads the F99 ladder: ¼ at the symmetric crossover (π/4),
        // ½ at the pure-YZ endpoint (π/2).
        var axis = DimensionAxis.Crossover(N: 3, gamma: 0.5, thetaPoints: 13);
        var sweep = DimensionSweep.Compute(axis, slowCount: 16);

        Assert.Equal(0.25, sweep.Polarity[6], 1e-12);   // θ = π/4
        Assert.Equal(0.5, sweep.Polarity[12], 1e-12);   // θ = π/2
        Assert.Equal(0.0, sweep.Polarity[0], 1e-12);    // θ = 0 (pure XZ)
    }

    [Fact]
    public void CumulativeRotation_RisesFromZero_BoundedByHalfPi()
    {
        // The resolving eyepiece: the largest principal angle (radians) between the slow subspace
        // at θ₀ and at θ[p], read from the singular values of Q(θ₀)ᴴ·Q(θ[p]). By construction
        // CumulativeRotation[0] = 0 (the subspace compared to itself), it rises as the in-between
        // turns with R_z(θ), and every principal angle is bounded by π/2.
        //
        // Measured note (N=3 crossover, slowCount=16): the LARGEST principal angle saturates near
        // π/2 already at the first nonzero θ-step (the 16-mode slow manifold is highly degenerate,
        // so some direction in it rotates out almost immediately). It does NOT grow linearly here,
        // it pins near π/2 and then jitters in that band; the spec anticipates this saturation
        // ("if it saturates at π/2 partway, that is fine and expected"). So we assert the robustly
        // true facts: starts at exactly 0, rises at the first step, stays in [0, π/2], ends > 0. We
        // deliberately do NOT assert strict step-to-step monotonicity: near the saturation ceiling
        // the value jitters, and at near-degenerate slow-mode crossings the membership of the
        // "16 smallest |Re λ|" set can swap, which the guard tolerates rather than crashing.
        var axis = DimensionAxis.Crossover(N: 3, gamma: 0.5, thetaPoints: 13);
        var sweep = DimensionSweep.Compute(axis, slowCount: 16);

        var cum = sweep.CumulativeRotation;
        Assert.Equal(13, cum.Count);
        Assert.Equal(0.0, cum[0], 1e-12);

        const double bound = Math.PI / 2.0 + 1e-9;
        foreach (double a in cum)
        {
            Assert.True(a >= -1e-12, $"principal angle must be ≥ 0; got {a:E3}");
            Assert.True(a <= bound, $"principal angle must be bounded by π/2; got {a:E3}");
        }

        Assert.True(cum[1] > 1e-6, $"the in-between should start rotating at the first step; got cum[1]={cum[1]:E3}");
        Assert.True(cum[^1] > 0.0, $"the in-between should have rotated by the end; got cum[^1]={cum[^1]:E3}");

        // The per-θ Q bases are exposed so a cumulative metric can be re-formed downstream.
        Assert.Equal(13, sweep.SlowBasis.Count);
        Assert.True(sweep.SlowBasis[0].ColumnCount > 0, "the slow basis should be non-empty");
    }

    [Fact]
    public void PrincipalAngleSpectrum_Resolves_The_Slow_Manifold_Split()
    {
        // The sharpened eyepiece. The single largest principal angle saturates on the degenerate slow
        // manifold (one direction rotates fully out at once), so it reports a saturated number and
        // hides the structure. The FULL spectrum, all k angles per θ sorted ascending, resolves it:
        // the small angles staying near 0 are the invariant core (the part of the slow manifold the
        // rotation fixes), the large angles that grow are the rotating in-between. The manifold
        // carries its own marks-and-in-between split, and the spectrum makes it visible.
        var axis = DimensionAxis.Crossover(N: 3, gamma: 0.5, thetaPoints: 13);
        var sweep = DimensionSweep.Compute(axis, slowCount: 16);

        var spectrum = sweep.PrincipalAngleSpectrum;
        Assert.Equal(13, spectrum.Count);

        // θ₀ vs itself: every principal angle is 0, one per slow-basis column.
        Assert.Equal(sweep.SlowBasis[0].ColumnCount, spectrum[0].Length);
        foreach (double a in spectrum[0])
            Assert.Equal(0.0, a, 1e-9);

        const double bound = Math.PI / 2.0 + 1e-9;
        for (int p = 0; p < spectrum.Count; p++)
        {
            double[] angles = spectrum[p];
            // Sorted ascending (invariant core first, rotated-out last).
            for (int i = 1; i < angles.Length; i++)
                Assert.True(angles[i] >= angles[i - 1] - 1e-12,
                    $"spectrum[{p}] must be ascending; {angles[i]:E3} < {angles[i - 1]:E3}");
            // Every angle in [0, π/2].
            foreach (double a in angles)
            {
                Assert.True(a >= -1e-12, $"angle must be ≥ 0; got {a:E3}");
                Assert.True(a <= bound, $"angle must be ≤ π/2; got {a:E3}");
            }
            // Consistency: CumulativeRotation keeps the largest angle (the last entry, ascending).
            double largest = angles.Length == 0 ? 0.0 : angles[^1];
            Assert.Equal(sweep.CumulativeRotation[p], largest, 1e-12);
        }

        // The resolving content: at the final θ the spectrum SPREADS (min ≪ max). The largest-angle
        // eyepiece reports only the max; the spread is the structure it hides, the manifold's own
        // invariant-core-vs-rotating split.
        double[] final = spectrum[^1];
        Assert.True(final.Length >= 2, "need at least two angles to resolve a split");
        double lo = final.Min(), hi = final.Max();
        Assert.True(hi - lo > 1e-3,
            $"the spectrum should spread at the final θ (the fan resolves the split); lo={lo:E3} hi={hi:E3}");
    }
}
