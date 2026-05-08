using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum), <b>Direction (a'') falsification probe</b>: tests
/// whether the empirical HWHM/Q* lift above the bare-doubled-PTF floor 0.671535 maps
/// monotonically to the Riesz-lifted SVD-block per-bond magnitude
/// <c>|V_b|²_proj/σ_0²</c> (<see cref="C2SvdBlockProjectedMagnitude"/>).
///
/// <para>If the bond-class signature on the magnitude has the SAME sign as the bond-class
/// signature on the HWHM-lift, Direction (a'') stays viable and a closed-form
/// <c>Δ(r) = HWHM_ratio − BareFloor</c> ansatz is the next analytical step. If the signs
/// are OPPOSITE — meaning Endpoint has larger Δ but smaller r, or vice versa — there is no
/// monotone single-variable map between the two, and Direction (a'') is structurally
/// falsified at the magnitude-only level (any closed form would need at least one
/// additional variable, e.g. cross-block Frobenius or a phase factor).</para>
///
/// <para><b>Tier outcome: Tier2Verified.</b> The probe is empirical: it runs the C2HwhmRatio
/// pipeline (full Q-scan per bond, which is the empirical anchor) plus the projector lift
/// across N=5..8, then compares per-class means. The verdict (<see cref="DirectionAVerdict"/>)
/// is derived bit-exactly from the witnesses, but the witnesses themselves are empirical
/// measurements. A future analytical proof of monotonicity (or refutation of it) can lift
/// this to Tier1Derived; until then the probe is a falsification gate, not a closed form.
/// </para>
///
/// <para>Construction is expensive (full Q-scan at every N). The typed claim caches the
/// witnesses on first build; downstream consumers should reuse the same instance. The
/// default N grid is {5, 6, 7, 8} matching the C2HwhmRatio anchor envelope.</para>
///
/// <para>Anchor: <c>compute/RCPsiSquared.Core/F86/Item1Derivation/C2HwhmRatio.cs</c>
/// PendingDerivationNote ranking (a'' is HIGH-RISK, this probe is the gate).</para>
/// </summary>
public sealed class C2DirectionAFalsificationProbe : Claim
{
    /// <summary>The bare-doubled-PTF floor against which the HWHM lift is measured.
    /// <c>Δ_b := HwhmLeftOverQPeak(b) − BareDoubledPtfHwhmRatio</c>.</summary>
    public const double BareFloor = C2HwhmRatio.BareDoubledPtfHwhmRatio;

    public IReadOnlyList<int> NValues { get; }
    public double GammaZero { get; }
    public IReadOnlyList<DirectionAProbeWitness> Witnesses { get; }
    public DirectionAVerdict Verdict { get; }

    public static C2DirectionAFalsificationProbe Build(
        IEnumerable<int>? nValues = null, double gammaZero = 0.05)
    {
        int[] ns = (nValues ?? new[] { 5, 6, 7, 8 }).ToArray();
        var witnesses = new List<DirectionAProbeWitness>(ns.Length);
        foreach (int N in ns)
        {
            var block = new CoherenceBlock(N: N, n: 1, gammaZero: gammaZero);
            var hwhm = C2HwhmRatio.Build(block);
            var magnitude = C2SvdBlockProjectedMagnitude.Build(block);

            double endpointHwhm = SafeMean(hwhm, BondClass.Endpoint);
            double interiorHwhm = SafeMean(hwhm, BondClass.Interior);
            double endpointR = magnitude.NormalisedRatioMean(BondClass.Endpoint);
            double interiorR = magnitude.NormalisedRatioMean(BondClass.Interior);

            double dEndpoint = endpointHwhm - BareFloor;
            double dInterior = interiorHwhm - BareFloor;
            double dDelta = dEndpoint - dInterior;
            double dR = endpointR - interiorR;

            SignRelation rel = ClassifySignRelation(dDelta, dR);

            witnesses.Add(new DirectionAProbeWitness(
                N: N,
                EndpointHwhmOverQpeak: endpointHwhm,
                InteriorHwhmOverQpeak: interiorHwhm,
                EndpointDelta: dEndpoint,
                InteriorDelta: dInterior,
                EndpointMagnitudeRatio: endpointR,
                InteriorMagnitudeRatio: interiorR,
                DeltaSignDifference: dDelta,
                RatioSignDifference: dR,
                SignRelation: rel));
        }

        DirectionAVerdict verdict = AggregateVerdict(witnesses);
        return new C2DirectionAFalsificationProbe(ns, gammaZero, witnesses, verdict);
    }

    private C2DirectionAFalsificationProbe(
        IReadOnlyList<int> nValues, double gammaZero,
        IReadOnlyList<DirectionAProbeWitness> witnesses, DirectionAVerdict verdict)
        : base("c=2 Direction (a'') falsification probe (HWHM-lift vs SVD-magnitude sign coherence)",
               Tier.Tier2Verified,
               Item1Anchors.Root)
    {
        NValues = nValues;
        GammaZero = gammaZero;
        Witnesses = witnesses;
        Verdict = verdict;
    }

    private static double SafeMean(C2HwhmRatio hwhm, BondClass bondClass)
    {
        try { return hwhm.HwhmLeftOverQPeakMean(bondClass); }
        catch { return double.NaN; }
    }

    private static SignRelation ClassifySignRelation(double dDelta, double dR)
    {
        const double tied = 1e-6;
        if (Math.Abs(dDelta) < tied || Math.Abs(dR) < tied) return SignRelation.Tied;
        return Math.Sign(dDelta) == Math.Sign(dR) ? SignRelation.Match : SignRelation.Mismatch;
    }

    private static DirectionAVerdict AggregateVerdict(IReadOnlyList<DirectionAProbeWitness> witnesses)
    {
        int matches = witnesses.Count(w => w.SignRelation == SignRelation.Match);
        int mismatches = witnesses.Count(w => w.SignRelation == SignRelation.Mismatch);
        int tied = witnesses.Count(w => w.SignRelation == SignRelation.Tied);

        if (mismatches == 0 && matches > 0) return DirectionAVerdict.StructurallyValid;
        if (matches == 0 && mismatches > 0) return DirectionAVerdict.Falsified;
        return DirectionAVerdict.Ambiguous;
    }

    public override string DisplayName =>
        $"c=2 Direction (a'') falsification probe (N ∈ {{{string.Join(", ", NValues)}}}, γ₀={GammaZero})";

    public override string Summary =>
        $"Verdict: {Verdict}. " +
        $"Witnesses: {Witnesses.Count(w => w.SignRelation == SignRelation.Match)} match / " +
        $"{Witnesses.Count(w => w.SignRelation == SignRelation.Mismatch)} mismatch / " +
        $"{Witnesses.Count(w => w.SignRelation == SignRelation.Tied)} tied across N values.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("BareFloor (BareDoubledPtfHwhmRatio)", BareFloor, "F6");
            yield return new InspectableNode("Verdict", summary: Verdict.ToString());
            yield return InspectableNode.Group("Per-N witnesses (signs of Endpoint − Interior)",
                Witnesses.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>Per-N witness for <see cref="C2DirectionAFalsificationProbe"/>: bond-class
/// means of the HWHM/Q* ratio and the SVD-block magnitude ratio, plus their sign-of-difference
/// relation. The probe asks whether <c>sign(Δ_E − Δ_I) = sign(r_E − r_I)</c> at every N.</summary>
public sealed record DirectionAProbeWitness(
    int N,
    double EndpointHwhmOverQpeak,
    double InteriorHwhmOverQpeak,
    double EndpointDelta,
    double InteriorDelta,
    double EndpointMagnitudeRatio,
    double InteriorMagnitudeRatio,
    double DeltaSignDifference,
    double RatioSignDifference,
    SignRelation SignRelation
) : IInspectable
{
    public string DisplayName => $"N={N} [{SignRelation}]";

    public string Summary =>
        $"Δ_E={EndpointDelta:+0.0000;-0.0000} Δ_I={InteriorDelta:+0.0000;-0.0000} " +
        $"(Δ_E−Δ_I={DeltaSignDifference:+0.0000;-0.0000}); " +
        $"r_E={EndpointMagnitudeRatio:F4} r_I={InteriorMagnitudeRatio:F4} " +
        $"(r_E−r_I={RatioSignDifference:+0.0000;-0.0000}); " +
        $"signs {SignRelation}.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("EndpointHwhmOverQpeak", EndpointHwhmOverQpeak, "F4");
            yield return InspectableNode.RealScalar("InteriorHwhmOverQpeak", InteriorHwhmOverQpeak, "F4");
            yield return InspectableNode.RealScalar("EndpointDelta", EndpointDelta, "F4");
            yield return InspectableNode.RealScalar("InteriorDelta", InteriorDelta, "F4");
            yield return InspectableNode.RealScalar("EndpointMagnitudeRatio", EndpointMagnitudeRatio, "F4");
            yield return InspectableNode.RealScalar("InteriorMagnitudeRatio", InteriorMagnitudeRatio, "F4");
            yield return InspectableNode.RealScalar("DeltaSignDifference (E-I)", DeltaSignDifference, "F4");
            yield return InspectableNode.RealScalar("RatioSignDifference (E-I)", RatioSignDifference, "F4");
            yield return new InspectableNode("SignRelation", summary: SignRelation.ToString());
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}

/// <summary>Sign relation between (Δ_E − Δ_I) and (r_E − r_I) at a fixed N: Match means
/// both have the same sign (HWHM-lift and SVD-magnitude split align), Mismatch means
/// opposite signs (anti-monotone), Tied means at least one of the two differences is below
/// numerical-noise tolerance and the relation is undefined.</summary>
public enum SignRelation { Match, Mismatch, Tied }

/// <summary>Aggregate verdict of <see cref="C2DirectionAFalsificationProbe"/> over the N
/// scan. StructurallyValid: every N is Match (no mismatches). Falsified: every N is
/// Mismatch (consistent anti-monotonicity rules out a single-variable monotone map).
/// Ambiguous: the N values disagree, suggesting at least one extra variable is needed.</summary>
public enum DirectionAVerdict { StructurallyValid, Falsified, Ambiguous }
