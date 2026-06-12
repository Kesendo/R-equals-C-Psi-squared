using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class ClockHandLadderWitnessTests
{
    private const double Tol = 1e-9;

    [Theory]
    [InlineData(3, 1.4142135623730951)] // sqrt(2)
    [InlineData(4, 1.6180339887498949)] // golden ratio phi
    [InlineData(5, 1.7320508075688772)] // sqrt(3)
    public void OmegaMem_AtNge3_IsTheF2bBandEdge_2JcosPiOverNPlus1(int n, double expected)
    {
        var w = new ClockHandLadderWitness();
        double omega = w.OmegaMem(n); // live from Symphony.Clock at J=1
        Assert.Equal(expected, omega, 6);
        double bandEdge = 2.0 * 1.0 * Math.Cos(Math.PI / (n + 1));
        Assert.True(Math.Abs(omega - bandEdge) < Tol,
            $"N={n}: live Omega {omega} must equal F2b band edge {bandEdge}");
    }

    [Fact]
    public void GammaProtection_AtN3_OmegaUnmoved_WhileGapTracks2Gamma()
    {
        var weak = new ClockHandLadderWitness(j: 1.0, gamma: 0.1);
        var strong = new ClockHandLadderWitness(j: 1.0, gamma: 0.4);
        // the coherence hand is γ-protected: quadrupling γ leaves it unmoved
        Assert.Equal(weak.OmegaMem(3), strong.OmegaMem(3), 9);
        // the Takt hand is NOT protected: it is exactly 2γ and tracks γ
        Assert.True(Math.Abs(weak.Gap(3) - 2.0 * 0.1) < Tol, $"weak gap {weak.Gap(3)} must be 2γ = 0.2");
        Assert.True(Math.Abs(strong.Gap(3) - 2.0 * 0.4) < Tol, $"strong gap {strong.Gap(3)} must be 2γ = 0.8");
        Assert.True(strong.Gap(3) > weak.Gap(3) + Tol, "the Takt hand moved with γ");
    }

    [Fact]
    public void N2_OmegaIsTheGammaPulledForm_2SqrtJ2MinusGamma2_NotTheBandEdge()
    {
        var w = new ClockHandLadderWitness(j: 1.0, gamma: 0.2); // Q = 5
        double omega = w.OmegaMem(2);
        double pulled = 2.0 * Math.Sqrt(1.0 - 0.04); // 1.95959
        Assert.Equal(pulled, omega, 6);
        // and it is NOT the F2b band edge (which would be 2cos(π/3) = 1.0)
        Assert.True(Math.Abs(omega - w.BandEdge(2)) > 0.5,
            $"N=2 Omega {omega} is the pulled block, distinct from the F2b band edge {w.BandEdge(2)}");
    }

    [Fact]
    public void N2_TheHandStops_AtTheExceptionalPoint_Q1()
    {
        // approach the EP: γ → J makes Q → 1 and the pulled coherence hand 2√(J²−γ²) → 0.
        // The HONEST witness near the EP is the closed form N2PulledOmega(), NOT the live
        // OmegaMem(2): past γ = √3/2·J the pulled coherence (here 0.28213) drops below the
        // {Im = ±J} modes that share the 2γ gap, so the live clock's "max |Im| at the gap"
        // tie-breaks to those ±J modes and reports 1.0, a true mode at the gap, but not the
        // pulled hand. See the N2_OmegaMem_TracksPulledForm_OnlyAboveTheCrossover test for the
        // boundary, and N2PulledOmega's note for the mechanism.
        var nearEp = new ClockHandLadderWitness(j: 1.0, gamma: 0.99);
        double pulled = nearEp.N2PulledOmega(); // 2√(1 − 0.9801) = 0.28213
        Assert.True(pulled < 0.3, $"near Q=1 the pulled hand is nearly stopped: {pulled}");
        Assert.Equal(2.0 * Math.Sqrt(1.0 - 0.99 * 0.99), pulled, 6);
        // at the EP exactly (γ = J, Q = 1) the hand is fully stopped:
        var atEp = new ClockHandLadderWitness(j: 1.0, gamma: 1.0);
        Assert.Equal(0.0, atEp.N2PulledOmega(), 12);
    }

    [Fact]
    public void N2_OmegaMem_TracksPulledForm_OnlyAboveTheCrossover()
    {
        // The live OmegaMem(2) equals the pulled closed form 2√(J²−γ²) while that pulled
        // coherence sits above the equal-rate {Im=±J} modes at the 2γ gap, i.e. for
        // 2√(J²−γ²) > J ⟺ γ < √3/2·J (Q > 2/√3 ≈ 1.1547). Below the crossover the live clock
        // reports the ±J modes (max |Im| at the gap = J), while the closed form keeps falling
        // to the EP. This documents the live/closed-form divergence as the finding, not a bug.
        var above = new ClockHandLadderWitness(j: 1.0, gamma: 0.2); // Q = 5, well above crossover
        Assert.Equal(above.N2PulledOmega(), above.OmegaMem(2), 6);

        var below = new ClockHandLadderWitness(j: 1.0, gamma: 0.99); // Q ≈ 1.01, below crossover
        Assert.Equal(1.0, below.OmegaMem(2), 6);                     // live clock: the ±J gap modes
        Assert.True(below.N2PulledOmega() < below.OmegaMem(2),
            $"below the crossover the pulled hand {below.N2PulledOmega()} has dropped under the " +
            $"live ±J gap modes {below.OmegaMem(2)}");
    }

    [Fact]
    public void Angle_AtN3_DefaultRegime_MatchesArctanOmegaOverGap()
    {
        var w = new ClockHandLadderWitness(); // J=1, γ=0.2
        double expected = Math.Atan(Math.Sqrt(2.0) / (2.0 * 0.2)) * 180.0 / Math.PI;
        Assert.Equal(expected, w.AngleDegrees(3), 6);
    }

    [Fact]
    public void Witness_SurfacesAllFourStories_AndNamesTheClaimAndDocs()
    {
        var w = new ClockHandLadderWitness();
        Assert.Contains("ClockHandLadderClaim", w.Summary);
        var labels = ((IInspectable)w).Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("ladder"));
        Assert.Contains(labels, l => l.Contains("protection"));
        Assert.Contains(labels, l => l.Contains("exceptional point"));
        Assert.Contains(labels, l => l.Contains("dial"));
    }
}
