using System.Linq;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The three-ladder hinge (inspect --root ladders): Q is the hinge. The bridge identity
/// P_{m,1} = m·Tr(Q·A^{m−1}) = the girth moments is the gate (port fidelity to
/// simulations/three_ladders_bridge.py); the rung-essential check shows Q is not incidental.</summary>
public class LadderHingeWitnessTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void bridge_identity_holds_at_m3_and_m5(int n)
    {
        // P_{m,1} = m*Tr(Q.A^{m-1}) (rung-weighted girth-walks) == the supertrace girth-moment form.
        var w = new LadderHingeWitness(n);
        Assert.True(w.BridgeIdentityDev(3) < 1e-9, $"bridge broke at m=3, N={n}: dev {w.BridgeIdentityDev(3)}");
        Assert.True(w.BridgeIdentityDev(5) < 1e-9, $"bridge broke at m=5, N={n}: dev {w.BridgeIdentityDev(5)}");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void the_rung_is_essential(int n)
    {
        // Remove the rung weighting (Q -> I): the coefficient is nowhere near the girth-moment form.
        var w = new LadderHingeWitness(n);
        double gap = Math.Abs(w.RungLessCoefficient(3) - w.GirthMomentCoefficient(3));
        Assert.True(gap > 1.0, $"rung-less should differ from the girth-moment form; gap {gap}");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void per_rung_decomposition_sums_to_the_coefficient(int n)
    {
        var w = new LadderHingeWitness(n);
        double sum = w.PerRungContribution(3).Values.Sum();
        Assert.True(Math.Abs(sum - w.RungWeightedCoefficient(3)) < 1e-7,
            $"per-rung sum {sum} != coefficient {w.RungWeightedCoefficient(3)}");
    }

    [Fact]
    public void guards_against_too_large_N()
    {
        // 4^6 = 4096 > MaxDim (1024)
        Assert.Throws<ArgumentOutOfRangeException>(() => new LadderHingeWitness(6));
    }
}
