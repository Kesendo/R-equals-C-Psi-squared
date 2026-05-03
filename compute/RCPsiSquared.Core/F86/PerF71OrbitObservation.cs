using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier 2 empirical observation refining Statement 3 (F71 mirror invariance).
/// Within the Interior bond class, per-F71-orbit Q_peak values are not uniform — the orbit
/// structure shows non-trivial sub-classification with two qualitatively distinct regimes:
///
/// <list type="bullet">
///   <item><b>F71-pairing (Tier 1):</b> Q_peak(b) = Q_peak(N−2−b) bit-exact (Statement 3).</item>
///   <item><b>N=6 central-vs-flanking inversion between c=2 and c=3:</b> the self-paired
///         central bond gives Q_peak BELOW the flanking inner orbit at c=2 N=6
///         (1.43 vs 1.63), but ABOVE it at c=3 N=6 (1.71 vs 1.66).</item>
///   <item><b>Two-regime K-curve structure at c=2 N≥7:</b> verified at c=2 N=7 and N=8
///         on extended grid [0.2, 20] × 100 pts. Two distinct regimes coexist:
///         <list type="bullet">
///           <item><i>Canonical bonds</i> (Endpoint, orbit 2): K_b peaks at Q ≈ Q_EP_canonical
///                 (2.5 for Endpoint, 1.5 for Interior orbit 2) and decays afterwards.</item>
///           <item><i>Saturation-plateau bonds</i> (orbit 1 at N≥7, central self-paired at N=8):
///                 K_b grows past Q_EP and forms a BROAD HIGH-Q PLATEAU where K stays elevated.
///                 N=7 orbit 1: plateau around Q ≈ 7-8 with K ≈ 0.15. N=8 orbit 1: plateau at
///                 Q ≈ 7.6-8.6, K ≈ 0.153. N=8 central self-paired b=3: plateau at Q ≈ 15-17,
///                 K ≈ 0.099. The "Q_peak" reading is the argmax of this plateau.</item>
///         </list>
///         The plateau is NOT a near-EP region: Petermann factor probe at c=2 N=7 across
///         Q ∈ [0.3, 20] shows max K_n ≈ 6 at the plateau (Q ≈ 7-8) — eigenvectors are
///         well-separated. Two clear near-EPs sit at the canonical peaks instead: Q ≈ 1
///         (max K_n = 123, near Interior canonical 1.5) and Q ≈ 3 (max K_n = 181, near
///         Endpoint canonical 2.5-2.6). The plateau is therefore a multi-mode coherent
///         enhancement effect at deep post-EP, not an EP-driven Petermann collapse.</item>
/// </list>
///
/// <para>The simple "Endpoint vs Interior" dichotomy is the leading approximation. The
/// full per-F71-orbit structure is finer-grained, c-dependent, and N-dependent: at large
/// N the orbit-1 (and eventually central) bonds transition into a saturation regime where
/// K_b plateaus past the slowest-pair Q_EP rather than decaying. Closed-form classification
/// of both the orbit pattern and the plateau enhancement remains open.</para>
///
/// <para>Connection to <see href="hypotheses/FRAGILE_BRIDGE.md">FRAGILE_BRIDGE</see>:
/// same Petermann diagnostic on a different parameter axis. There K = 403 above γ_crit
/// signals an EP in the complex γ plane; here K = 181 at Q ≈ 3 signals a near-EP on the
/// real Q axis at Endpoint canonical peak. The c=2 N=7 full block-L has structural near-EPs
/// at canonical Q_peak positions, confirming the EP mechanism for those peaks while
/// disconnecting it from the high-Q saturation plateau.</para>
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
