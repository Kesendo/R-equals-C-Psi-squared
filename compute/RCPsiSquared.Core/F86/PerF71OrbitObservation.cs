using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier 2 empirical observation refining Statement 3 (F71 mirror invariance):
/// per-F71-orbit Q_peak within Interior is not uniform, and some orbits carry further maxima
/// far past their first one.
///
/// <para>Full discussion: PROOF_F86_QPEAK Open Element 5. The witnesses are <b>frozen
/// measurements</b> of the envelope max_t |K_b(Q, t)| (the <see cref="ResonanceScan"/> reading) at
/// γ₀ = 0.05, one value per orbit (the lower-index bond of the mirror pair; the pair agrees to the
/// scan's resolution), in two named readings:</para>
/// <list type="bullet">
///   <item><b>canonical</b> (<see cref="OrbitWitness.QPeakPerOrbit"/>): the first local maximum
///         from below (<see cref="KCurve.CanonicalPeakAtBond"/>), read on the default Q grid
///         (153 points on [0.2, 4]). It lies between Q 1.428 and 2.591 for every orbit here and
///         does not move when the window is extended;</item>
///   <item><b>global on a named window</b> (<see cref="OrbitWitness.Plateaus"/>): where the
///         orbit's maximum over a wider window lies above its canonical one, that maximum and the
///         window. Read for c=2 N=7 and N=8 on [0.2, 20] and N=8 also on [0.2, 60]; the curves keep further
///         local maxima (c=2 N=8 orbit 1 has three on [0.2, 60], the endpoint a second one near
///         Q 32 below its canonical height), so this reading is a property of the window.</item>
/// </list>
/// <para>What the numbers are: the reading on the default Q grid (step 0.025), the parabolic
/// vertex through three grid points and a linearly interpolated half-maximum. The second route,
/// simulations/f86_envelope_second_route.py (global and canonical readings, a Q_HI argument for
/// wider windows), reproduces that reading to ~10⁻¹⁰ where it has been run (c=2 N=5..8 per bond; the class values of the
/// fixture): the peak value K is fixed to ~10⁻¹³, and the vertex on a grid of step 0.025 turns
/// that into a few 10⁻¹¹ in Q_peak at c=2 N=5..8 and up to ~10⁻¹⁰ at c=3, the size of the spread
/// between the two bonds of a mirror pair within one route (1.0·10⁻¹⁰ at c=3 N=8, 7·10⁻¹⁰ at
/// c=4 N=8). The reading is not the continuous envelope
/// peak: against a continuous search in Q it sits high by 7·10⁻⁵ to 2.2·10⁻⁴ in Q_peak and off by
/// up to 3.5·10⁻⁵ in HWHM_left/Q_peak (c=2 N=5..8, simulations/f86_envelope_second_route.py --continuous N). The two far maxima of c=2 N=8 beyond Q = 20 (orbit 2 and the
/// central orbit) are read by the second route alone. Not derived from a closed form; the
/// closed-form classification of the orbit pattern and of the far maxima remains open.</para>
/// </summary>
public sealed class PerF71OrbitObservation : Claim
{
    public IReadOnlyList<OrbitWitness> Witnesses { get; }

    public PerF71OrbitObservation()
        : base("per-F71-orbit substructure within Interior class",
               Tier.Tier2Empirical,
               "PROOF_F86_QPEAK open front 2, the F86d Tier-1 promotion (partially addressed by PROOF_F86C_F71_MIRROR Statement 3)")
    {
        // Frozen envelope readings; re-measure if the ResonanceScan envelope or its grid changes.
        Witnesses = new[]
        {
            new OrbitWitness(2, 5, new[] { 2.517071370614458, 1.4831232668717471 }),
            new OrbitWitness(2, 6, new[] { 2.580604998574177, 1.6366804166405813, 1.4283108372718716 }, centralIndex: 2),
            new OrbitWitness(2, 7, new[] { 2.5600734278105572, 1.698919918615899, 1.5010592087425514 },
                plateaus: new[] { new PlateauReading(1, 7.255447920341961, 0.14886801406506361, 20) }),
            new OrbitWitness(2, 8, new[] { 2.543127121124538, 1.727187504810217, 1.5207477539505165, 1.5698435864972127 }, centralIndex: 3,
                plateaus: new[]
                {
                    new PlateauReading(1, 8.224722736060128, 0.15351937853386294, 20),
                    new PlateauReading(2, 38.207793066213945, 0.09285593229315857, 60),
                    new PlateauReading(3, 20.091690478319904, 0.10091516738120127, 60),
                }),
            new OrbitWitness(3, 5, new[] { 2.3926106535196983, 1.5996673298568214 }),
            new OrbitWitness(3, 6, new[] { 2.5481946106587965, 1.6622167263981307, 1.7117377732328158 }, centralIndex: 2,
                note: "central > flanking, opposite of c=2 N=6"),
            new OrbitWitness(3, 7, new[] { 2.5836685556845613, 1.7276621884295194, 1.6965676610188718 }),
            new OrbitWitness(3, 8, new[] { 2.583262803492589, 1.7628004489513558, 1.6929165430288617, 1.705793514967033 }, centralIndex: 3),
            new OrbitWitness(4, 7, new[] { 2.5752121936601067, 1.7422447718487073, 1.7787917896506178 }),
            new OrbitWitness(4, 8, new[] { 2.5908339756199688, 1.7807955537329054, 1.8070132707126416, 1.781729923217086 }, centralIndex: 3),
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
                summary: "c=2 N=6: central 1.428 < flanking 1.637; c=3 N=6: central 1.712 > flanking 1.662");
            yield return new InspectableNode("high-Q feature",
                summary: "c=2 N≥7: past a minimum the curves rise again, and some far maxima exceed the canonical one: orbit 1 at Q 7.26 (N=7) and 8.22 (N=8), at N=8 also orbit 2 at Q 38.2 and the central orbit at Q 20.09, read on [0.2, 60]");
        }
    }
}

/// <summary>The maximum of one orbit over a named window [0.2, <see cref="WindowHigh"/>], where it
/// lies above the orbit's canonical (first) maximum.</summary>
public sealed record PlateauReading(int Orbit, double QPeak, double KMax, double WindowHigh);

/// <summary>One (c, N) per-orbit Q_peak witness. Orbit indices run from 0 (endpoint pair)
/// inward; <see cref="CentralIndex"/> marks the self-paired central orbit if present.</summary>
public sealed class OrbitWitness : IInspectable
{
    public int Chromaticity { get; }
    public int N { get; }
    public IReadOnlyList<double> QPeakPerOrbit { get; }
    public int? CentralIndex { get; }
    public IReadOnlyList<PlateauReading> Plateaus { get; }
    public string? Note { get; }

    /// <param name="qPeakPerOrbit">The canonical Q_peak per orbit (first local maximum).</param>
    /// <param name="plateaus">The orbits whose maximum over a wider window lies higher.</param>
    public OrbitWitness(int c, int n, IReadOnlyList<double> qPeakPerOrbit,
        int? centralIndex = null, IReadOnlyList<PlateauReading>? plateaus = null, string? note = null)
    {
        Chromaticity = c;
        N = n;
        QPeakPerOrbit = qPeakPerOrbit;
        CentralIndex = centralIndex;
        Plateaus = plateaus ?? Array.Empty<PlateauReading>();
        Note = note;
    }

    public string DisplayName => $"c={Chromaticity} N={N}";
    public string Summary
    {
        get
        {
            var parts = QPeakPerOrbit
                .Select((q, i) => i == CentralIndex ? $"orbit{i}*={q:F4}" : $"orbit{i}={q:F4}");
            string body = string.Join(" ", parts);
            if (Plateaus.Count > 0)
                body += "; far maxima " + string.Join(" ", Plateaus.Select(
                    p => $"orbit{p.Orbit}={p.QPeak:F4} (K {p.KMax:F4}, on [0.2, {p.WindowHigh:G}])"));
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
            foreach (var p in Plateaus)
                yield return InspectableNode.RealScalar($"orbit {p.Orbit} maximum on [0.2, {p.WindowHigh:G}]", p.QPeak, "F4");
            if (Note is not null) yield return new InspectableNode("note", summary: Note);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
