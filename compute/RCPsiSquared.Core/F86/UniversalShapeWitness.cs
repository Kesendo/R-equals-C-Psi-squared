using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>One empirical witness for the F86 universal-shape claim at (c, N, γ₀, BondClass).
/// Self-computing <see cref="IInspectable"/>: the <see cref="HwhmLeftOverQPeak"/>,
/// <see cref="QPeak"/>, and <see cref="KMax"/> values are produced on demand by running a
/// <see cref="ResonanceScan"/> through the shared <see cref="WitnessCache"/> — no hardcoded
/// numbers.
///
/// <para>Together the witness list is what makes <see cref="UniversalShapePrediction"/> a
/// Tier-1-candidate (multi-N stability across c=2..4) rather than a single-anchor
/// coincidence. By computing each witness from the actual physics, the claim's evidence
/// base is always synchronised with the current scan implementation — no drift from
/// historical hardcoded numbers.</para>
///
/// <para>Lazy contract: <see cref="DisplayName"/> and <see cref="Children"/> shape do not
/// trigger computation; only accessing <see cref="Summary"/>, <see cref="Payload"/>, or any
/// of the Peak-derived properties forces the underlying scan. This means a
/// <c>Walk()</c>-with-tier-filter that doesn't materialise summaries stays scan-free.</para>
/// </summary>
public sealed class UniversalShapeWitness : IInspectable
{
    public int Chromaticity { get; }
    public int N { get; }
    public double GammaZero { get; }
    public Resonance.BondClass BondClass { get; }
    public WitnessCache Cache { get; }

    private readonly Lazy<PeakResult> _peak;

    public PeakResult Peak => _peak.Value;
    public double QPeak => Peak.QPeak;
    public double KMax => Peak.KMax;

    /// <summary>The HWHM_left / Q_peak ratio — the F86 universal-shape numerical witness.
    /// Throws if the scan returned no left HWHM (curve doesn't drop below half-max on the left).</summary>
    public double HwhmLeftOverQPeak => Peak.HwhmLeftOverQPeak
        ?? throw new InvalidOperationException(
            $"no left HWHM at c={Chromaticity} N={N} γ₀={GammaZero} {BondClass}");

    public UniversalShapeWitness(int chromaticity, int n, double gammaZero, Resonance.BondClass bondClass,
        WitnessCache? cache = null)
    {
        Chromaticity = chromaticity;
        N = n;
        GammaZero = gammaZero;
        BondClass = bondClass;
        Cache = cache ?? WitnessCache.Default;
        _peak = new Lazy<PeakResult>(() => Cache.GetOrCompute(Chromaticity, N, GammaZero).Peak(BondClass));
    }

    public string DisplayName => $"witness c={Chromaticity} N={N} {BondClass}";

    public string Summary =>
        $"HWHM_left/Q_peak = {HwhmLeftOverQPeak:F4}, Q_peak = {QPeak:F4}, |K|max = {KMax:F4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("chromaticity", Chromaticity);
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("γ₀", GammaZero, "G3");
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("HWHM_left/Q_peak", HwhmLeftOverQPeak, "F4");
            yield return InspectableNode.RealScalar("Q_peak", QPeak, "F4");
            yield return InspectableNode.RealScalar("|K|max", KMax, "F4");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real("HWHM_left/Q_peak", HwhmLeftOverQPeak, "F4");
}
