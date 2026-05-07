using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.F86;

/// <summary>One F71 orbit's measured K-resonance values from the c=2 C2HwhmRatio pipeline:
/// Q_peak, HWHM_left, |K|max, plus the derived HWHM_left/Q_peak ratio. 2-bond orbits hold
/// the F71-mirror-pair average (max-deviation guarded in <see cref="PerF71OrbitKTable.Build"/>);
/// self-paired orbits hold the single-bond witness.</summary>
public sealed record OrbitKWitness(
    F71BondOrbit Orbit,
    double QPeak,
    double HwhmLeft,
    double KMax
) : IInspectable
{
    /// <summary>HWHM_left/Q_peak ratio, derived from the two stored fields. Computing it
    /// here avoids the mean-of-ratios pitfall when 2-bond orbits average the underlying
    /// witnesses (mean of HwhmLeft / mean of QPeak ≠ mean of HwhmLeft/QPeak in general,
    /// though the F71-mirror guard makes the difference negligible).</summary>
    public double HwhmLeftOverQPeak => HwhmLeft / QPeak;

    /// <summary>True when QPeak sits within one grid spacing of the upper bound of
    /// <paramref name="qGrid"/> — the peak finder pinned to the edge, no detected
    /// resonance in the scanned range. Uses the last grid interval as dQ so log-spaced
    /// grids work without changing the contract.</summary>
    public bool IsEscaped(IReadOnlyList<double> qGrid)
    {
        if (qGrid is null) throw new ArgumentNullException(nameof(qGrid));
        if (qGrid.Count < 2)
            throw new ArgumentException("qGrid must have at least 2 points", nameof(qGrid));
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
