using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's telescope onto the exceptional point: the birth of the rotation. Sweeps
/// Q across Q_EP = 2/g_eff and reads, via <see cref="ExceptionalPointClock"/>, the in-between between the
/// anchors Q=0, Q_EP, Q_peak. Five readings: the marks (Q_EP, Q_peak, t_peak=1/(4γ₀)); the Takt (the slow
/// mode's decay coalescing the two real channels at −4γ₀ and pinning there); the Rotation (the angle θ
/// lifting off 0 above the EP, the oscillation born, = the F95 angle); the defectiveness (the eigenvector
/// overlap min(x,1/x) → 1 at the EP, the Jordan-block pinch); and the hardware (the IBM Kingston EP-onset,
/// the memory revival lifting off the 1/N floor as Q crosses Q_EP). Closed-form, N-free. The mirror of the
/// interior axis: there the rotation stills at ¼, here it is born at Q_EP (the cusp/EP F95 siblinghood).
/// The post-EP regime it opens is the BirthCanal of <see cref="PostEpFlowField"/>.</summary>
public sealed class EpField : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private readonly double _gamma;
    private readonly double _gEff;
    private readonly double[] _qGrid;

    public EpField(double gamma = 0.5, double gEff = 4.0 / 3.0, double qLo = 0.3, double qHi = 4.0,
        int qPoints = 41)
    {
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        if (gEff <= 0) throw new ArgumentOutOfRangeException(nameof(gEff), $"gEff must be positive; got {gEff}");
        if (qLo < 0 || qHi <= qLo) throw new ArgumentOutOfRangeException(nameof(qHi), $"need 0 ≤ qLo < qHi; got [{qLo}, {qHi}]");
        if (qPoints < 2) throw new ArgumentOutOfRangeException(nameof(qPoints), $"need at least two Q points; got {qPoints}");
        _gamma = gamma;
        _gEff = gEff;
        _qGrid = new double[qPoints];
        for (int i = 0; i < qPoints; i++) _qGrid[i] = qLo + (qHi - qLo) * i / (qPoints - 1);
    }

    private double QEp => ExceptionalPointClock.QEp(_gEff);
    private double QPeak => ExceptionalPointClock.QPeak(_gEff);

    public string DisplayName =>
        $"EpField (the exceptional point, g_eff={_gEff.ToString("0.###", Inv)}, Q_EP={QEp.ToString("0.###", Inv)}, Q∈[{_qGrid[0].ToString("0.##", Inv)}, {_qGrid[^1].ToString("0.##", Inv)}])";

    public string Summary =>
        $"sweeping Q across Q_EP={QEp.ToString("0.###", Inv)}: two real decay channels coalesce defectively at −4γ₀ " +
        $"(the Takt pins, t_peak=1/(4γ₀)), the Rotation angle lifts off 0 (the oscillation born), the eigenvectors " +
        $"collapse (overlap→1). The cusp's mirror: there the rotation stills at ¼, here it is born. Hardware-anchored " +
        "(IBM Kingston EP onset). The post-EP regime it opens is the BirthCanal.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. The marks (the anchors, the contract).
            yield return new InspectableNode(
                displayName: "the marks (the anchors)",
                summary: $"Q=0 (rest), Q_EP={QEp.ToString("0.###", Inv)} (the defective coalescence, the rotation born), " +
                         $"Q_peak={QPeak.ToString("0.###", Inv)} (the K_b resonance peak, past the EP). At the EP the slow mode " +
                         $"pins at 4γ₀ so the Takt period is t_peak=1/(4γ₀). Inert, the contract.");

            // 2. The Takt: the eigenvalue coalescence, the decay pinning at 4γ₀.
            var decay = _qGrid.Select(q => ExceptionalPointClock.Decay(_gamma, q, _gEff)).ToArray();
            yield return new InspectableNode(
                displayName: "the Takt: the two channels coalesce at −4γ₀",
                summary: $"the slow mode's decay −Re(λ) rises from ~2γ₀ to 4γ₀ as Q→Q_EP (the two real rates merging), then pins at 4γ₀ above the EP (the Takt hand stops). The absorption saturates at the ⟨n⟩=2 rung.",
                payload: new InspectablePayload.Curve("decay −Re(λ) vs Q", _qGrid, decay, "Q", "decay [γ₀] (pins at 4γ₀)"));

            // 3. The Rotation: the angle lifting off (the oscillation born).
            var theta = _qGrid.Select(q => ExceptionalPointClock.RotationAngleDegrees(_gamma, q, _gEff)).ToArray();
            yield return new InspectableNode(
                displayName: "the Rotation: the angle lifts off (the rotation born)",
                summary: $"θ=arctan(ω/gap) (the F95 angle): 0 below Q_EP (pure Takt, overdamped), lifting off above it ({theta[^1].ToString("0.#", Inv)}° at the far end). The imaginary part appears: the oscillation, the memory, is born.",
                payload: new InspectablePayload.Curve("Rotation angle θ° vs Q", _qGrid, theta, "Q", "θ° (lifts off at Q_EP)"));

            // 4. The defectiveness: the eigenvector overlap → 1 (the fragile pinch).
            var overlap = _qGrid.Select(q => ExceptionalPointClock.EigenvectorOverlap(q, _gEff)).ToArray();
            yield return new InspectableNode(
                displayName: "the defectiveness (the fragile-bridge pinch)",
                summary: "the eigenvector overlap |⟨v₊|v₋⟩| = min(x,1/x), x=Q/Q_EP: peaks at 1 exactly at the EP (the two eigenvectors coalesce, a Jordan block, Petermann sensitivity → ∞), falling off on either side. The defective signature that makes it an exceptional point, not a mere degeneracy.",
                payload: new InspectablePayload.Curve("eigenvector overlap vs Q", _qGrid, overlap, "Q", "|⟨v₊|v₋⟩| (→1 at the EP)"));

            // 5. The hardware (IBM Kingston EP onset, Part B).
            var hwQ = new double[] { 0.5, 1.0, 1.5, 2.5, 5.0, 20.0 };
            var hwRev = new double[] { 0.30, 0.36, 0.34, 0.49, 0.56, 0.70 };   // data/ibm_ep_onset_may2026 README, job d8drjbfd0j8c73f4mobg
            yield return new InspectableNode(
                displayName: "the hardware (IBM Kingston EP onset)",
                summary: "IBM Kingston 2026-05-31: a single-excitation walk under injected dephasing, Q swept {0.5,1,1.5,2.5,5,20}. The memory revival (max⟨n₀⟩) stays at the 1/N≈0.33 floor up to Q≈Q_EP, then lifts off (0.34→0.49→0.56→0.70): the rotation born on a real chip (ep_onset_may2026, job d8drjbfd0j8c73f4mobg).",
                payload: new InspectablePayload.Curve("Kingston revival vs Q", hwQ, hwRev, "Q", "revival (memory return)"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
