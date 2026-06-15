using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The envelope_n4_rise boundary chart, gate-first. The CΨ Envelope Theorem is proven for N=2;
/// the full-state envelope GENUINELY RISES at N≥4 strong coupling (CpsiEnvelopeTheoremClaim). These tests
/// chart WHERE that rise turns on, reusing the witness's own detector (EnvelopeTheoremWitness.GlobalReading).
///
/// The finding: the boundary is not a sharp N-step and not a pure J/γ contour — it is BOTH, cleanly
/// factored. (1) The rise count is a pure (N, Q=J/γ) observable: the J-sweep and the γ-sweep collapse to
/// one Q-axis (the (Q,K)-purity gate below — it CAN fail on an absolute-time leak, and passing validates
/// the collapse). (2) An N≥4 FLOOR: N=3 never rises (one internal site cannot inject), N≥4 can. (3) Above
/// the floor a Q threshold Q_c(N) that RISES with N: Q_c(4)≈27, Q_c(5)≈45, the rise strength at fixed Q
/// falling with N. All readings at the witness regime's dose window (γ·tMax = 0.25). See
/// experiments/ENVELOPE_RISE_BOUNDARY.md.</summary>
public class EnvelopeBoundaryTests
{
    private const double Bar = EnvelopeTheoremWitness.GenuinenessBar;

    [Fact]
    public void Purity_GlobalRiseIsAPure_N_Q_Observable_JAndGammaCollapseToOneAxis()
    {
        // The clock movement certifies the dimensionless lenses are pure (Q,K)-observables: under L→rL
        // (J,γ both scaled by r, the dose window K=γt held fixed) the global CΨ(K) curve is invariant. So
        // the rise reading over a fixed dose window depends only on (N, Q=J/γ) — the "J_c sweep" and the
        // "γ_c sweep" are the SAME axis. Two regimes at Q=500 related by r=2 (5,0.01,tMax=25) vs
        // (10,0.02,tMax=12.5), both K_max=0.25, must give the bit-identical reading. A firing here would be
        // an absolute-time leak in the detector, not physics — the gate that earns the J,γ→Q collapse.
        var a = EnvelopeTheoremWitness.GlobalReading(4, 5.0, 0.01, 25.0, 1600);
        var b = EnvelopeTheoremWitness.GlobalReading(4, 10.0, 0.02, 12.5, 1600);
        Assert.Equal(a.RiseCount, b.RiseCount);
        Assert.Equal(a.MaxRiseMagnitude, b.MaxRiseMagnitude, 6);
    }

    [Fact]
    public void NFloor_N3HoldsEvenAtHugeQ_TheRiseNeedsAnInternalPair()
    {
        // N=3 (carrier pair 0-1, ONE internal site) holds non-increasing even at Q=2000 — 100× the N=4
        // threshold. The rise is the Part-6 coherence injection: it needs an internal ≥2-site coherent
        // subsystem to pump CΨ back, which N=3 lacks. So the boundary has an N≥4 FLOOR (Q_c(3)=∞), not a
        // threshold reachable by turning Q up at N=3.
        Assert.Equal(0, EnvelopeTheoremWitness.GlobalRiseCount(3, 20.0, 0.01, 25.0, 1600));   // Q=2000
        Assert.True(EnvelopeTheoremWitness.GlobalReading(4, 5.0, 0.01, 25.0, 1600).MaxRiseMagnitude
                    > Bar);                                                                    // Q=500: N=4 rises
    }

    [Fact]
    public void QContour_N4_RisesAboveThreshold_HoldsWellBelow()
    {
        // Above the floor the rise is a Q=J/γ contour: at N=4 the global envelope rises genuinely and
        // grid-stably for Q≳30 (here Q=40, maxMag≈0.011 ≫ bar) and holds (non-increasing) well below the
        // ~18-28 fuzzy band (here Q=13). Q_c(4)≈27.
        Assert.True(EnvelopeTheoremWitness.GlobalReading(4, 0.40, 0.01, 25.0, 1600).MaxRiseMagnitude
                    > Bar);                                                                    // Q=40: rises
        Assert.Equal(0, EnvelopeTheoremWitness.GlobalRiseCount(4, 0.13, 0.01, 25.0, 1600));    // Q=13: holds
    }

    [Fact]
    public void QContour_ThresholdRisesWithN_N5NeedsMoreQThanN4()
    {
        // The threshold Q_c(N) climbs with N. At Q=40 the N=4 envelope rises genuinely but the N=5 one does
        // not (its whole rise curve is shifted down — injection weakens as the internal bath grows from a
        // pair to a trio). N=5 still rises higher up (Q=500), so it is a shifted contour, not a new floor.
        double n4_q40 = EnvelopeTheoremWitness.GlobalReading(4, 0.40, 0.01, 25.0, 1600).MaxRiseMagnitude;
        double n5_q40 = EnvelopeTheoremWitness.GlobalReading(5, 0.40, 0.01, 25.0, 1600).MaxRiseMagnitude;
        double n5_q500 = EnvelopeTheoremWitness.GlobalReading(5, 5.0, 0.01, 25.0, 1600).MaxRiseMagnitude;
        Assert.True(n4_q40 > Bar);    // N=4 rises at Q=40
        Assert.True(n5_q40 < Bar);    // N=5 holds at Q=40 (Q_c(5) > Q_c(4))
        Assert.True(n5_q500 > Bar);   // but N=5 does rise at Q=500
    }
}
