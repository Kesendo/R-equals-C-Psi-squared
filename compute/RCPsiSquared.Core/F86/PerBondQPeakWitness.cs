using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 per-bond Q_peak witness at (c, N, BondClass) — self-computing
/// <see cref="IInspectable"/> backed by the shared <see cref="WitnessCache"/>. Reads
/// <see cref="PeakResult.QPeak"/> on demand from the cached <see cref="ResonanceScan"/>.
///
/// <para>The values do not match clean closed forms (the csc-formulas were retracted).</para>
/// </summary>
public sealed class PerBondQPeakWitness : IInspectable
{
    public int Chromaticity { get; }
    public int N { get; }
    public BondClass BondClass { get; }
    public double GammaZero { get; }
    public WitnessCache Cache { get; }

    private readonly Lazy<PeakResult> _peak;
    public PeakResult Peak => _peak.Value;
    public double QPeak => Peak.QPeak;

    public PerBondQPeakWitness(int chromaticity, int n, BondClass bondClass,
        double gammaZero = 0.05, WitnessCache? cache = null)
    {
        Chromaticity = chromaticity;
        N = n;
        BondClass = bondClass;
        GammaZero = gammaZero;
        Cache = cache ?? WitnessCache.Default;
        _peak = new Lazy<PeakResult>(() => Cache.GetOrCompute(Chromaticity, N, GammaZero).Peak(BondClass));
    }

    public string DisplayName => $"Q_peak c={Chromaticity} N={N} {BondClass}";
    public string Summary => $"Q_peak = {QPeak:F4} (computed)";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("chromaticity", Chromaticity);
            yield return InspectableNode.RealScalar("N", N);
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("Q_peak", QPeak, "F4");
        }
    }

    public InspectablePayload Payload => new InspectablePayload.Real("Q_peak", QPeak, "F4");
}

/// <summary>Group container of <see cref="PerBondQPeakWitness"/> for one bond class —
/// self-computing typed table form of the per-bond Q_peak rows in F86.</summary>
public sealed class PerBondQPeakWitnessTable : F86Claim
{
    public BondClass BondClass { get; }
    public IReadOnlyList<PerBondQPeakWitness> Witnesses { get; }
    public bool ExpectedToSaturate { get; }
    public double? ExpectedSaturationValue { get; }

    public PerBondQPeakWitnessTable(BondClass bondClass, IReadOnlyList<PerBondQPeakWitness> witnesses,
        bool expectedToSaturate, double? expectedSaturationValue)
        : base($"per-bond Q_peak table ({bondClass})",
               Tier.Tier2Empirical,
               "docs/ANALYTICAL_FORMULAS.md F86 + docs/proofs/PROOF_F86_QPEAK.md Empirical Q_peak data")
    {
        BondClass = bondClass;
        Witnesses = witnesses;
        ExpectedToSaturate = expectedToSaturate;
        ExpectedSaturationValue = expectedSaturationValue;
    }

    public override string DisplayName => $"per-bond Q_peak ({BondClass}, {Witnesses.Count} cases)";

    public override string Summary =>
        ExpectedToSaturate && ExpectedSaturationValue is { } sv
            ? $"expected to saturate near {sv:F2} ({Tier.Label()}, computed live)"
            : $"no expected saturation, {Witnesses.Count} witnesses ({Tier.Label()}, computed live)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            foreach (var w in Witnesses) yield return w;
            if (ExpectedSaturationValue is { } sv)
                yield return InspectableNode.RealScalar("expected saturation value", sv, "F4");
        }
    }

    public static PerBondQPeakWitnessTable BuildEndpoint(double gammaZero = 0.05, WitnessCache? cache = null) =>
        new(BondClass.Endpoint,
            F86StandardLocations.Full.Select(loc =>
                new PerBondQPeakWitness(loc.C, loc.N, BondClass.Endpoint, gammaZero, cache)).ToArray(),
            expectedToSaturate: true,
            expectedSaturationValue: 2.53);

    public static PerBondQPeakWitnessTable BuildInterior(double gammaZero = 0.05, WitnessCache? cache = null) =>
        new(BondClass.Interior,
            F86StandardLocations.Full.Select(loc =>
                new PerBondQPeakWitness(loc.C, loc.N, BondClass.Interior, gammaZero, cache)).ToArray(),
            expectedToSaturate: false,
            expectedSaturationValue: null);
}
