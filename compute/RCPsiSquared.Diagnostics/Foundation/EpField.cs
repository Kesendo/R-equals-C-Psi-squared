using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's telescope onto the <b>toy 2×2 rate-channel exceptional point</b>: the
/// birth of the rotation in the 2-level reduction. Sweeps Q across Q_EP = 2/g_eff and reads, via
/// <see cref="ExceptionalPointClock"/>, the in-between between the anchors Q=0, Q_EP, Q_peak. Five
/// readings: the marks (Q_EP, Q_peak, t_peak=1/(4γ₀)); the Takt (the toy slow mode's decay coalescing the
/// two real channels at −4γ₀ and pinning there); the Rotation (the angle θ lifting off 0 above the toy
/// EP, the oscillation born, = the F95 angle); the defectiveness (the toy eigenvector overlap min(x,1/x)
/// → 1 at the toy EP, the Jordan-block pinch of the 2×2 reduction); plus a separately scoped hardware
/// population handover. The IBM Kingston curve is converted from its coherence-rate labels to the
/// repository's canonical Lindblad Q; it is not evidence for the toy EP or for chain spectral character.
/// Closed-form, N-free, 2-level.
///
/// <para><b>Scope (F86a-retraction, 2026-06-21):</b> every "EP / defective / Jordan-block / overlap→1"
/// reading on this field is a property of the toy 2×2 rate-channel reduction, which IS genuinely
/// defective at Q_EP. The full physical (n,n+1)-coherence block of the chain is genuinely non-normal on
/// the real Q axis (large but finite Petermann), and F89 has separate scattered real-axis defective
/// seeds; it does not inherit the toy's clean Q_EP pinch. Do not read this telescope as the chain
/// block's real-axis behaviour. See
/// <c>RCPsiSquared.Core.F86.LocalGlobalEpLink</c> (OpenQuestion) and
/// <c>docs/proofs/PROOF_F86A_EP_MECHANISM.md</c> §The real-axis EP.</para></summary>
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
        $"the toy 2×2 rate-channel EP, sweeping Q across Q_EP={QEp.ToString("0.###", Inv)}: the toy's two real decay " +
        $"channels coalesce defectively at −4γ₀ (the Takt pins, t_peak=1/(4γ₀)), the Rotation angle lifts off 0 (the " +
        $"oscillation born), the 2×2 eigenvectors collapse (overlap→1). This defectiveness is the toy reduction's clean pinch AT Q_EP, NOT " +
        $"the physical chain block's behaviour there (the block is non-normal at Q_peak; its own real-axis defective EPs are " +
        $"F89's scattered seeds, not a clean Q_EP pinch, the retraction's 'no real-axis EP' corrected 2026-07-07; see LocalGlobalEpLink). " +
        $"The hardware child is a separate population handover, not a calibration or confirmation of this toy EP; " +
        "its spectral character remains open.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. The marks (the anchors, the contract).
            yield return new InspectableNode(
                displayName: "the marks (the anchors)",
                summary: $"Q=0 (rest), Q_EP={QEp.ToString("0.###", Inv)} (the toy 2×2's defective coalescence, the rotation born), " +
                         $"Q_peak={QPeak.ToString("0.###", Inv)} (the K_b resonance peak, past the toy EP). At the toy EP the slow mode " +
                         $"pins at 4γ₀ so the Takt period is t_peak=1/(4γ₀). Inert, the contract.");

            // 2. The Takt: the toy 2×2's eigenvalue coalescence, the decay pinning at 4γ₀.
            var decay = _qGrid.Select(q => ExceptionalPointClock.Decay(_gamma, q, _gEff)).ToArray();
            yield return new InspectableNode(
                displayName: "the Takt: the two channels coalesce at −4γ₀",
                summary: $"the toy 2×2's slow-mode decay −Re(λ) rises from ~2γ₀ to 4γ₀ as Q→Q_EP (the toy's two real rates merging), then pins at 4γ₀ above the toy EP (the Takt hand stops). The absorption saturates at the ⟨n⟩=2 rung.",
                payload: new InspectablePayload.Curve("decay −Re(λ) vs Q", _qGrid, decay, "Q", "decay [γ₀] (pins at 4γ₀)"));

            // 3. The Rotation: the angle lifting off (the oscillation born).
            var theta = _qGrid.Select(q => ExceptionalPointClock.RotationAngleDegrees(_gamma, q, _gEff)).ToArray();
            yield return new InspectableNode(
                displayName: "the Rotation: the angle lifts off (the rotation born)",
                summary: $"θ=arctan(ω/gap) (the F95 angle) of the toy 2×2: 0 below Q_EP (pure Takt, overdamped), lifting off above it ({theta[^1].ToString("0.#", Inv)}° at the far end). The toy's imaginary part appears: the oscillation, the memory, is born.",
                payload: new InspectablePayload.Curve("Rotation angle θ° vs Q", _qGrid, theta, "Q", "θ° (lifts off at Q_EP)"));

            // 4. The defectiveness: the toy 2×2's eigenvector overlap → 1 (the toy EP pinch).
            var overlap = _qGrid.Select(q => ExceptionalPointClock.EigenvectorOverlap(q, _gEff)).ToArray();
            yield return new InspectableNode(
                displayName: "the defectiveness (the toy 2×2 EP pinch)",
                summary: "the toy 2×2's eigenvector overlap |⟨v₊|v₋⟩| = min(x,1/x), x=Q/Q_EP: peaks at 1 exactly at the toy EP " +
                         "(the two 2×2 eigenvectors coalesce, a Jordan block, the toy's Petermann sensitivity diverging), falling off on " +
                         "either side. This is the defective signature of the toy 2×2 rate-channel reduction, a genuine EP OF THE TOY. " +
                         "It is NOT the physical (n,n+1) chain block's behaviour at Q_peak, where the block is non-normal (large but " +
                         "FINITE Petermann). The block's OWN real-axis defective EPs are F89's scattered seeds (the retraction's " +
                         "'no real-axis coalescence' was corrected 2026-07-07); see LocalGlobalEpLink / PROOF_F86A section The real-axis EP.",
                payload: new InspectablePayload.Curve("eigenvector overlap vs Q", _qGrid, overlap, "Q", "|⟨v₊|v₋⟩| (→1 at the EP)"));

            // 5. The hardware (IBM Kingston finite-time population handover, Part B).
            // Registry anchor: ConfirmationsRegistry.Lookup("ibm_ep_onset_may2026") (mirrored in simulations/framework/confirmations.py).
            var hwQLabel = new double[] { 0.5, 1.0, 1.5, 2.5, 5.0, 20.0 };
            var hwQLindblad = hwQLabel.Select(q => 2.0 * q).ToArray();
            // Exact raw-record values; simulations/docs_verify.py checks this array against the JSON populations.
            var hwRev = new double[] { 0.2978515625, 0.3623809814453125, 0.34356689453125,
                0.4898834228515625, 0.560699462890625, 0.70330810546875 };
            yield return new InspectableNode(
                displayName: "the hardware (IBM Kingston single-excitation-walk population handover)",
                summary: "IBM Kingston 2026-05-31 (real data): single-excitation populations under injected random-Z-twirl dephasing. " +
                         "The runner uses coherence-rate labels Q_label=J/Γ: σ²=2Γdt gives exp(−Γt). The repository jump √γ Z gives exp(−2γt), " +
                         "so γ=Γ/2 and Q_Lindblad = 2 Q_label. The sampled population handover Q_label=1.5→2.5 is therefore Q_Lindblad=3→5 " +
                         "(revival 0.344→0.490). The 1/N line is only a reference level for this finite-time scalar, not a lower floor or an equipartition certificate. " +
                         "Populations alone do not locate an EP, critical damping, or mode coalescence; spectral character remains open. " +
                         "The historical ep_onset_may2026 name identifies the run, not its present verdict.",
                payload: new InspectablePayload.Curve("Kingston revival vs canonical Lindblad Q", hwQLindblad, hwRev,
                    "Q_Lindblad = 2 Q_label", "revival (population return)"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
