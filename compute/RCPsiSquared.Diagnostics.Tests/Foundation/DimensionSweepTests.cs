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
}
