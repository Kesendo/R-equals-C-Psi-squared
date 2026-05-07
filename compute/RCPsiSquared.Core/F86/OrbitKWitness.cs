using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.F86;

/// <summary>F71-orbit-decorated K-resonance witness: one orbit's measured Q_peak,
/// HWHM_left, HWHM_left/Q_peak, |K|max from the c=2 C2HwhmRatio pipeline. K-prefix
/// distinguishes from the existing per-(c,N) <see cref="OrbitWitness"/> in
/// <see cref="PerF71OrbitObservation"/> (Q_peak-array record, frozen empirical).
///
/// <para>For 2-bond orbits the values are the F71-mirror-pair average (max-deviation
/// guarded by <see cref="PerF71OrbitKTable.Build"/>); for self-paired orbits the
/// values are the single-bond witness directly. No baked-in EscapeFlag — see
/// <see cref="IsEscaped"/> for grid-relative escape classification.</para></summary>
public sealed record OrbitKWitness(
    F71BondOrbit Orbit,
    double QPeak,
    double HwhmLeft,
    double HwhmLeftOverQPeak,
    double KMax
) : IInspectable
{
    /// <summary>Returns true when QPeak sits within one grid spacing of the upper bound
    /// of the given Q-grid — i.e., the peak finder pinned to the grid edge and the
    /// orbit has no detected resonance peak in the scanned range. Grid-relative by
    /// design (no pre-baked EscapeFlag).</summary>
    public bool IsEscaped(IReadOnlyList<double> qGrid)
    {
        if (qGrid is null) throw new ArgumentNullException(nameof(qGrid));
        if (qGrid.Count < 2)
            throw new ArgumentException("qGrid must have at least 2 points", nameof(qGrid));
        // Use the last grid interval as dQ: correct for any monotone grid (uniform
        // or widening toward qMax). For the default uniform grid both endpoints
        // give identical dQ; this form additionally works for log-spaced grids
        // without changing the contract.
        double dQ = qGrid[^1] - qGrid[^2];
        double qMax = qGrid[^1];
        return QPeak >= qMax - dQ;
    }

    public string DisplayName =>
        Orbit.IsSelfPaired
            ? $"OrbitKWitness {{b={Orbit.BondA}}} (self-paired)"
            : $"OrbitKWitness {{b={Orbit.BondA} ↔ b={Orbit.BondB}}}";

    public string Summary =>
        $"Q_peak={QPeak:F4}, HWHM-/Q*={HwhmLeftOverQPeak:F4}, |K|max={KMax:G4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return Orbit;
            yield return InspectableNode.RealScalar("Q_peak", QPeak, "F4");
            yield return InspectableNode.RealScalar("HWHM_left", HwhmLeft, "F4");
            yield return InspectableNode.RealScalar("HWHM_left/Q_peak", HwhmLeftOverQPeak, "F4");
            yield return InspectableNode.RealScalar("|K|max", KMax, "G4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
