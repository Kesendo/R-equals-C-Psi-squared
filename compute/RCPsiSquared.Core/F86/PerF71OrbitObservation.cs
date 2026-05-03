using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier 2 empirical observation refining Statement 3 (F71 mirror invariance).
/// Within the Interior bond class, per-F71-orbit Q_peak values are not uniform — the orbit
/// structure shows non-trivial sub-classification. Three effects observed:
///
/// <list type="bullet">
///   <item><b>F71-pairing (Tier 1):</b> Q_peak(b) = Q_peak(N−2−b) bit-exact (Statement 3).</item>
///   <item><b>N=6 central-vs-flanking inversion between c=2 and c=3:</b> the self-paired
///         central bond gives Q_peak BELOW the flanking inner orbit at c=2 N=6
///         (1.43 vs 1.63), but ABOVE it at c=3 N=6 (1.71 vs 1.66).</item>
///   <item><b>High-Q secondary peaks at c=2 N≥7:</b> inner-non-central bonds (b=1, b=N−3)
///         show Q_peak shifted to high Q (>6 in observed grid), while central bonds keep
///         the canonical Interior peak ~1.5. Likely related to coupling with an
///         Endpoint-like eigenmode at higher Q.</item>
/// </list>
///
/// <para>The simple "Endpoint vs Interior" dichotomy is the leading approximation; the full
/// per-F71-orbit structure is finer-grained and c-dependent. Closed-form classification
/// remains open.</para>
/// </summary>
public sealed class PerF71OrbitObservation : F86Claim
{
    public IReadOnlyList<OrbitWitness> Witnesses { get; }

    public PerF71OrbitObservation()
        : base("per-F71-orbit substructure within Interior class",
               Tier.Tier2Empirical,
               "PROOF_F86_QPEAK Open Element 5 (partially addressed by Statement 3)")
    {
        Witnesses = new[]
        {
            new OrbitWitness(2, 5, new[] { 2.50, 1.49 }),
            new OrbitWitness(2, 6, new[] { 2.57, 1.63, 1.43 }, centralIndex: 2),
            new OrbitWitness(2, 7, new[] { 2.56, 6.00, 1.50 }, note: "inner orbit Q_peak off-grid (>6)"),
            new OrbitWitness(2, 8, new[] { 2.53, 6.00, 1.52, 1.58 }, centralIndex: 3, note: "inner orbit off-grid"),
            new OrbitWitness(3, 5, new[] { 2.39, 1.61 }),
            new OrbitWitness(3, 6, new[] { 2.54, 1.66, 1.71 }, centralIndex: 2,
                note: "central > flanking, opposite of c=2 N=6"),
            new OrbitWitness(3, 7, new[] { 2.59, 1.74, 1.71 }),
            new OrbitWitness(3, 8, new[] { 2.60, 1.76, 1.70, 1.71 }, centralIndex: 3),
            new OrbitWitness(4, 7, new[] { 2.61, 1.75, 1.78 }),
        };
    }

    public override string DisplayName => "Per-F71-orbit Q_peak substructure";

    public override string Summary =>
        $"{Witnesses.Count} (c, N) cases; central/flanking inversion at N=6 between c=2 and c=3 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            foreach (var w in Witnesses) yield return w;
            yield return new InspectableNode("inversion observation",
                summary: "c=2 N=6: central 1.43 < flanking 1.63; c=3 N=6: central 1.71 > flanking 1.66");
            yield return new InspectableNode("high-Q feature",
                summary: "c=2 N≥7: inner-non-central bonds (b=1, b=N−3) have Q_peak >6 (off-grid); central pair stays at canonical ~1.5");
        }
    }
}

/// <summary>One (c, N) per-orbit Q_peak witness. Orbit indices run from 0 (endpoint pair)
/// inward; <see cref="CentralIndex"/> marks the self-paired central orbit if present.</summary>
public sealed class OrbitWitness : IInspectable
{
    public int Chromaticity { get; }
    public int N { get; }
    public IReadOnlyList<double> QPeakPerOrbit { get; }
    public int? CentralIndex { get; }
    public string? Note { get; }

    public OrbitWitness(int c, int n, IReadOnlyList<double> qPeakPerOrbit,
        int? centralIndex = null, string? note = null)
    {
        Chromaticity = c;
        N = n;
        QPeakPerOrbit = qPeakPerOrbit;
        CentralIndex = centralIndex;
        Note = note;
    }

    public string DisplayName => $"c={Chromaticity} N={N}";
    public string Summary
    {
        get
        {
            var parts = QPeakPerOrbit
                .Select((q, i) => i == CentralIndex ? $"orbit{i}*={q:F2}" : $"orbit{i}={q:F2}");
            string body = string.Join(" ", parts);
            return Note is null ? body : $"{body} ({Note})";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("c", Chromaticity);
            yield return InspectableNode.RealScalar("N", N);
            for (int i = 0; i < QPeakPerOrbit.Count; i++)
            {
                string label = i == CentralIndex ? $"orbit {i} (self-paired)" : $"orbit {i}";
                yield return InspectableNode.RealScalar(label, QPeakPerOrbit[i], "F4");
            }
            if (Note is not null) yield return new InspectableNode("note", summary: Note);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
