using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The stone (felt_time_dimensions arc, step B), graduated to a live C# witness. The PTF
/// painter closure Sum_i ln(alpha_i), run on the mode-isolating probe rho_0 = I/d + eps*Herm(mode)
/// through the CANONICAL Symphony FitAlpha, reads the chosen mode's first-order RATE shift: it BREAKS
/// (out of the +-0.05 window) AND is sign-coherent for the soft survivor (Re moves), and HOLDS (in
/// window) for the rigid band edge (Re frozen). Mirrors simulations/stone_survivor_alpha_closure.py
/// (two-lens reviewed 2026-06-19); the review-pinned scope is the probe-state-specific, sign-coherence-
/// certified claim, not a universal trajectory law.</summary>
public class StoneSurvivorClosureWitnessTests
{
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void Survivor_breaks_the_closure_as_a_sign_coherent_rate_shift(int n)
    {
        var w = new StoneSurvivorClosureWitness(n);
        Assert.False(w.Survivor.InWindow);              // OUT of +-0.05: K_decay (Re lambda) moves
        Assert.True(w.Survivor.SignCoherence > 0.8);    // sign-coherent reliable f => a real rate shift
        Assert.True(w.Survivor.IsRateShift);
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void Band_edge_holds_the_closure(int n)
    {
        var w = new StoneSurvivorClosureWitness(n);
        Assert.True(w.BandEdge.InWindow);               // IN window: K_decay frozen (rigid darkness)
        Assert.False(w.BandEdge.IsRateShift);
    }
}
