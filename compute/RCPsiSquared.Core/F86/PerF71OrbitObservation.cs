using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier 2 empirical observation refining Statement 3 (F71 mirror invariance):
/// per-F71-orbit Q_peak within Interior is not uniform; the orbit structure shows
/// non-trivial sub-classification with two regimes (canonical bonds peak-and-decay vs
/// saturation-plateau bonds at large N).
///
/// <para>Full discussion + numeric data: see PROOF_F86_QPEAK Open Element 5. The
/// witnesses below are <b>frozen measurements</b> from extended-grid scans
/// (LinearQGrid(0.2, 6.0, 60) for the 9-case sweep; LinearQGrid(0.2, 20.0, 100) for the
/// c=2 N=7 and N=8 high-Q verifications, using <c>F86NewIdeasProbe</c>) and not derived
/// from a closed form. Closed-form classification of orbit pattern and plateau enhancement
/// remains open.</para>
/// </summary>
public sealed class PerF71OrbitObservation : Claim
{
    public IReadOnlyList<OrbitWitness> Witnesses { get; }

    public PerF71OrbitObservation()
        : base("per-F71-orbit substructure within Interior class",
               Tier.Tier2Empirical,
               "PROOF_F86_QPEAK Open Element 5 (partially addressed by Statement 3)")
    {
        // Frozen empirical witnesses, measured via F86NewIdeasProbe.
        // 9-case sweep: LinearQGrid(0.2, 6.0, 60) at γ₀=0.05.
        // c=2 N=7 and N=8 numbers (orbit 1 Q≈7-8, central Q≈17): LinearQGrid(0.2, 20.0, 100).
        // Re-measure if the underlying ResonanceScan implementation or grid changes.
        Witnesses = new[]
        {
            new OrbitWitness(2, 5, new[] { 2.50, 1.49 }),
            new OrbitWitness(2, 6, new[] { 2.57, 1.63, 1.43 }, centralIndex: 2),
            new OrbitWitness(2, 7, new[] { 2.56, 7.24, 1.50 }, note: "orbit 1 saturation plateau Q≈7-8 K≈0.15 dominant over canonical 1.8 K≈0.11"),
            new OrbitWitness(2, 8, new[] { 2.55, 8.07, 1.53, 16.79 }, centralIndex: 3, note: "orbit 1 plateau Q≈8 K=0.15; central also plateau Q≈17 K=0.10"),
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
