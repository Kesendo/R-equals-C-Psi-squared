using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum), <b>Direction (c'') falsification probe</b>: tests
/// whether the 4-mode K-resonance reproduces the empirical HWHM_left/Q_peak per bond
/// class, or sits at the bare-doubled-PTF floor 0.671535.
///
/// <para>Direction (c'') proposed a three-block superposition K_total = K_pb + K_sv +
/// 2·Re(K_cross) within the 4-mode L_eff basis. The structural concern is that ANY
/// decomposition within the 4-mode subspace cannot exceed the 4-mode K-resonance itself.
/// If the full 4-mode K_b sits at the bare-doubled-PTF floor (≈ 0.671) for both bond
/// classes, then no additive (pb + sv + cross) splitter can produce the empirical
/// 0.7506 / 0.7728 — the +0.08 / +0.10 lift simply does not live in the 4-mode subspace.
/// Direction (c'') would be structurally falsified by the same argument that already
/// disqualified the 4-mode reduction itself.</para>
///
/// <para>This probe runs the falsification in ~5 seconds: per N ∈ {5, 6, 7, 8}, build
/// <see cref="FourModeResonanceScan"/>, compute the K-curve, extract per-class
/// HWHM_left/Q_peak via <see cref="ParabolicPeakFinder"/>, and compare with the empirical
/// values from <see cref="C2HwhmRatio"/>. The verdict labels each N:</para>
/// <list type="bullet">
///   <item><b>FourModeMatchesEmpirical</b>: 4-mode HWHM/Q_peak agrees with empirical
///   within tolerance for both bond classes. (c'') has structural content the diagonal
///   anatomy missed.</item>
///   <item><b>FourModeAtBareFloor</b>: 4-mode HWHM/Q_peak sits at the bare-doubled-PTF
///   floor regardless of bond class. (c'') is structurally trivial; 4-mode insufficient.</item>
///   <item><b>FourModePartial</b>: 4-mode lifts above the floor but does not reach
///   empirical. The gap quantifies how much of the lift requires modes beyond the 4-mode
///   subspace.</item>
/// </list>
///
/// <para><b>Tier outcome: Tier2Verified.</b> Empirical falsification probe; the verdict is
/// derived from the witnesses but the witnesses themselves are numerical measurements.</para>
///
/// <para>Anchor: <c>compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs</c>
/// PendingDerivationNote (c'') ranking.</para>
/// </summary>
public sealed class C2DirectionCFalsificationProbe : Claim
{
    public const double BareFloor = C2HwhmRatio.BareDoubledPtfHwhmRatio;

    /// <summary>Tolerance for "FourModeMatchesEmpirical": gap |4-mode − empirical| below
    /// this threshold for BOTH bond classes counts as a match.</summary>
    public const double MatchTolerance = 0.005;

    /// <summary>Tolerance for "FourModeAtBareFloor": gap |4-mode − BareFloor| below this
    /// threshold for BOTH bond classes counts as floor-locked.</summary>
    public const double FloorTolerance = 0.01;

    public IReadOnlyList<int> NValues { get; }
    public IReadOnlyList<DirectionCProbeWitness> Witnesses { get; }
    public DirectionCVerdict Verdict { get; }

    public static C2DirectionCFalsificationProbe Build(IEnumerable<int>? nValues = null)
    {
        int[] ns = (nValues ?? new[] { 5, 6, 7, 8 }).ToArray();
        var witnesses = new List<DirectionCProbeWitness>(ns.Length);
        foreach (int N in ns)
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: 0.05);
            var scan = new FourModeResonanceScan(block);
            var kCurve = scan.ComputeKCurve();

            double endpointHwhmOverQPeak = SafePeakRatio(kCurve, BondClass.Endpoint);
            double interiorHwhmOverQPeak = SafePeakRatio(kCurve, BondClass.Interior);

            var empiric = C2HwhmRatio.Build(block);
            double endpointEmpirical = SafeEmpiricalMean(empiric, BondClass.Endpoint);
            double interiorEmpirical = SafeEmpiricalMean(empiric, BondClass.Interior);

            witnesses.Add(new DirectionCProbeWitness(
                N: N,
                FourModeEndpointHwhmRatio: endpointHwhmOverQPeak,
                FourModeInteriorHwhmRatio: interiorHwhmOverQPeak,
                EmpiricalEndpointHwhmRatio: endpointEmpirical,
                EmpiricalInteriorHwhmRatio: interiorEmpirical,
                FourModeFloorGapEndpoint: endpointHwhmOverQPeak - BareFloor,
                FourModeFloorGapInterior: interiorHwhmOverQPeak - BareFloor,
                FourModeEmpiricalGapEndpoint: endpointHwhmOverQPeak - endpointEmpirical,
                FourModeEmpiricalGapInterior: interiorHwhmOverQPeak - interiorEmpirical));
        }

        DirectionCVerdict verdict = AggregateVerdict(witnesses);
        return new C2DirectionCFalsificationProbe(ns, witnesses, verdict);
    }

    private C2DirectionCFalsificationProbe(
        IReadOnlyList<int> nValues,
        IReadOnlyList<DirectionCProbeWitness> witnesses,
        DirectionCVerdict verdict)
        : base("c=2 Direction (c'') falsification probe (4-mode HWHM/Q_peak vs empirical)",
               Tier.Tier2Verified,
               Item1Anchors.Root)
    {
        NValues = nValues;
        Witnesses = witnesses;
        Verdict = verdict;
    }

    private static double SafePeakRatio(KCurve curve, BondClass bondClass)
    {
        try
        {
            var peak = curve.Peak(bondClass);
            if (peak.HwhmLeft is null) return double.NaN;
            return peak.HwhmLeft.Value / peak.QPeak;
        }
        catch { return double.NaN; }
    }

    private static double SafeEmpiricalMean(C2HwhmRatio hwhm, BondClass bondClass)
    {
        try { return hwhm.HwhmLeftOverQPeakMean(bondClass); }
        catch { return double.NaN; }
    }

    private static DirectionCVerdict AggregateVerdict(IReadOnlyList<DirectionCProbeWitness> witnesses)
    {
        bool allMatch = witnesses.All(w =>
            !double.IsNaN(w.FourModeEmpiricalGapEndpoint) &&
            !double.IsNaN(w.FourModeEmpiricalGapInterior) &&
            Math.Abs(w.FourModeEmpiricalGapEndpoint) < MatchTolerance &&
            Math.Abs(w.FourModeEmpiricalGapInterior) < MatchTolerance);
        if (allMatch) return DirectionCVerdict.FourModeMatchesEmpirical;

        bool allFloor = witnesses.All(w =>
            !double.IsNaN(w.FourModeFloorGapEndpoint) &&
            !double.IsNaN(w.FourModeFloorGapInterior) &&
            Math.Abs(w.FourModeFloorGapEndpoint) < FloorTolerance &&
            Math.Abs(w.FourModeFloorGapInterior) < FloorTolerance);
        if (allFloor) return DirectionCVerdict.FourModeAtBareFloor;

        return DirectionCVerdict.FourModePartial;
    }

    public override string DisplayName =>
        $"c=2 Direction (c'') falsification probe (N ∈ {{{string.Join(", ", NValues)}}})";

    public override string Summary => $"Verdict: {Verdict}";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("BareFloor", BareFloor, "F6");
            yield return new InspectableNode("Verdict", summary: Verdict.ToString());
            yield return InspectableNode.Group("Per-N witnesses",
                Witnesses.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>Per-N witness for <see cref="C2DirectionCFalsificationProbe"/>.</summary>
public sealed record DirectionCProbeWitness(
    int N,
    double FourModeEndpointHwhmRatio,
    double FourModeInteriorHwhmRatio,
    double EmpiricalEndpointHwhmRatio,
    double EmpiricalInteriorHwhmRatio,
    double FourModeFloorGapEndpoint,    // 4-mode minus BareFloor (positive = 4-mode lifted above floor)
    double FourModeFloorGapInterior,
    double FourModeEmpiricalGapEndpoint, // 4-mode minus empirical (negative = 4-mode below empirical)
    double FourModeEmpiricalGapInterior
) : IInspectable
{
    public string DisplayName => $"N={N}";

    public string Summary =>
        $"4-mode E={FourModeEndpointHwhmRatio:F4} I={FourModeInteriorHwhmRatio:F4}; " +
        $"empirical E={EmpiricalEndpointHwhmRatio:F4} I={EmpiricalInteriorHwhmRatio:F4}; " +
        $"4-mode−empirical: E={FourModeEmpiricalGapEndpoint:+0.0000;-0.0000} " +
        $"I={FourModeEmpiricalGapInterior:+0.0000;-0.0000}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("FourModeEndpointHwhmRatio", FourModeEndpointHwhmRatio, "F4");
            yield return InspectableNode.RealScalar("FourModeInteriorHwhmRatio", FourModeInteriorHwhmRatio, "F4");
            yield return InspectableNode.RealScalar("EmpiricalEndpointHwhmRatio", EmpiricalEndpointHwhmRatio, "F4");
            yield return InspectableNode.RealScalar("EmpiricalInteriorHwhmRatio", EmpiricalInteriorHwhmRatio, "F4");
            yield return InspectableNode.RealScalar("FourModeFloorGapEndpoint", FourModeFloorGapEndpoint, "F4");
            yield return InspectableNode.RealScalar("FourModeFloorGapInterior", FourModeFloorGapInterior, "F4");
            yield return InspectableNode.RealScalar("FourModeEmpiricalGapEndpoint", FourModeEmpiricalGapEndpoint, "F4");
            yield return InspectableNode.RealScalar("FourModeEmpiricalGapInterior", FourModeEmpiricalGapInterior, "F4");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}

/// <summary>Aggregate verdict of <see cref="C2DirectionCFalsificationProbe"/>:
/// FourModeMatchesEmpirical (4-mode reaches the empirical anchor; (c'') has structural
/// content), FourModeAtBareFloor (4-mode locked at BareFloor; (c'') structurally trivial),
/// or FourModePartial (4-mode partially lifted above the floor but not all the way to
/// empirical; the lift requires modes beyond the 4-mode subspace).</summary>
public enum DirectionCVerdict
{
    FourModeMatchesEmpirical,
    FourModeAtBareFloor,
    FourModePartial,
}
