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
        var w = new ClockHandLadderWitness(); // J=1, default γ
        // robust against the witness default: θ(N=3) = arctan(bandEdge √2 / Gap 2γ)
        double expected = Math.Atan(Math.Sqrt(2.0) / (2.0 * w.Gamma)) * 180.0 / Math.PI;
        Assert.Equal(expected, w.AngleDegrees(3), 6);
    }

    [Fact]
    public void Angle_AtN2_BelowCrossover_UsesHonestPulledHand_NotRawClock()
    {
        // Below the crossover Q=2/√3 the raw clock's max|Im| at the gap is the ±J band line, not the
        // pulled hand. The dial must report the angle of the hand that actually stops at the EP, i.e.
        // θ = arctan(√(Q²−1)) from the closed-form pulled hand, NOT arctan(J/2γ) from the raw clock.
        var w = new ClockHandLadderWitness(j: 1.0, gamma: 0.9); // Q = 1.111, below Q=2/√3 ≈ 1.155
        double q = 1.0 / 0.9;
        double honest = Math.Atan(Math.Sqrt(q * q - 1.0)) * 180.0 / Math.PI;   // 25.84°
        double rawClock = Math.Atan(1.0 / (2.0 * 0.9)) * 180.0 / Math.PI;      // 29.05°, the bug
        Assert.Equal(honest, w.AngleDegrees(2), 4);
        Assert.True(Math.Abs(w.AngleDegrees(2) - rawClock) > 1.0,
            $"the dial must not show the raw-clock angle {rawClock} below the crossover");
    }

    [Fact]
    public void Angle_AtN2_GoesToZero_AtTheExceptionalPoint()
    {
        // At the EP (γ=J, Q=1) the pulled hand stops, so the dial angle is exactly 0.
        var w = new ClockHandLadderWitness(j: 1.0, gamma: 1.0);
        Assert.Equal(0.0, w.AngleDegrees(2), 9);
    }

    [Fact]
    public void Ladder_AtDefaultGamma_IsInTheProtectedRegime_ShowsTheBandEdge()
    {
        var w = new ClockHandLadderWitness(); // γ=0.2, H-competitive
        Assert.True(w.BandEdgeIsTheGapMode(3) && w.BandEdgeIsTheGapMode(4) && w.BandEdgeIsTheGapMode(5));
        var ladder = ((IInspectable)w).Children.First(c => c.DisplayName.Contains("ladder"));
        Assert.Contains("F2b band edge", ladder.Summary);
        Assert.DoesNotContain("out of the protected regime", ladder.Summary);
    }

    [Fact]
    public void Ladder_InStrongDephasing_IsHonest_DoesNotClaimTheBandEdgeAsLive()
    {
        // At γ=0.9 the band-edge mode (rate 2γ) is no longer the slowest: a real overdamped mode takes
        // the gap, so OmegaMem(3) → 0. The ladder node must NOT print √2 as the live value; it must
        // report the regime exit honestly. (Regression guard for the γ=0.9 "N=3 → 0 (√2)" lie.)
        var w = new ClockHandLadderWitness(j: 1.0, gamma: 0.9);
        Assert.False(w.BandEdgeIsTheGapMode(3), "γ=0.9 is out of the protected regime");
        Assert.True(w.OmegaMem(3) < 0.5, "the live coherence hand has dropped (the gap mode is real)");
        Assert.True(w.Gap(3) < 2.0 * 0.9 - 1e-6, "the gap is no longer 2γ at strong dephasing");
        var ladder = ((IInspectable)w).Children.First(c => c.DisplayName.Contains("ladder"));
        Assert.Contains("out of the protected regime", ladder.Summary);
        var protection = ((IInspectable)w).Children.First(c => c.DisplayName.Contains("protection"));
        Assert.Contains("protection is exited", protection.Summary);
    }

    [Fact]
    public void Witness_SurfacesAllFourStories_AndNamesTheClaimAndDocs()
    {
        var w = new ClockHandLadderWitness();
        Assert.Contains("ClockHandLadderClaim", w.Summary);
        Assert.Equal(4, ((IInspectable)w).Children.Count());
        var labels = ((IInspectable)w).Children.Select(c => c.DisplayName).ToList();
        Assert.Contains(labels, l => l.Contains("ladder"));
        Assert.Contains(labels, l => l.Contains("protection"));
        Assert.Contains(labels, l => l.Contains("exceptional point"));
        Assert.Contains(labels, l => l.Contains("dial"));
    }
}
