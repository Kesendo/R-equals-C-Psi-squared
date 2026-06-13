using System.Globalization;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The dynamic survival probe, live (<c>inspect --root survivor</c>): WHERE the longest-lived
/// dissipative mode lives across the three physically-grounded topologies, and HOW its lifetime scales
/// with N. The thread-(b) companion to the <c>survival_incompleteness_mirror</c> arc; the C# witness for
/// <see cref="RCPsiSquared.Core.Symmetry.SurvivalIncompletenessMirrorClaim"/>.
///
/// <para><b>The labels are physical, not abstract graph types:</b>
/// <list type="bullet">
/// <item>CHAIN = extended <b>dispersive</b> 1D matter: conjugated polyene chains (polyacetylene/SSH, the
/// carbon work), spin chains, the Grotthuss proton wire (water). k_min = pi/(N+1).</item>
/// <item>RING = extended <b>dispersive</b> 1D periodic: aromatic rings (benzene, the Frost circle),
/// light-harvesting macrocycles. k_min = 2pi/N (4x the chain).</item>
/// <item>STAR = hub-spoke <b>NON-dispersive</b> central-spin: NV centre / quantum dot / donor-spin bath
/// (Bortz-Stolze), the optical-cavity point-focus (STAR_CONFOCAL_LIMIT), the mediator.</item>
/// </list></para>
///
/// <para><b>Finding (verified):</b> on the DISPERSIVE substrates (chain, ring) the survivor is the
/// INTERIOR (incompleteness, C=0.5) coherence - chain at the dead-centre half-filling, ring at the
/// off-centre (2,2)/(N-2,N-2); the central-spin STAR is the COUNTEREXAMPLE (boundary (1,1), no central
/// momentum mode). Lifetime &lt;n_XY&gt; ~ Q^2/N^2, ring/chain -&gt; 4 (cyclic-vs-open k_min^2,
/// model-independent). Reuses <see cref="SectorReductionWitness.SectorSlowest"/> (no full 4^N).</para>
///
/// <para>Convention (carbon, XY/free-fermion, no ZZ): J=1, gamma=1/Q. SectorReductionWitness builds
/// H=qh*(XX+YY) (qh=0.5 reproduces J=1) with an absolute gamma profile; rate = -2*gamma*&lt;n_XY&gt;.
/// The python twin is <c>simulations/carbon/incompleteness_survivor.py</c>.</para></summary>
public sealed class IncompletenessSurvivorWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Qh = 0.5;          // H = qh*(XX+YY) reproduces the carbon J=1
    private const double KernelTol = 1e-7;

    public int N { get; }
    public double Q { get; }                 // J/gamma in the carbon convention (J=1, gamma=1/Q)

    public IncompletenessSurvivorWitness(int n = 6, double q = 1.5)
    {
        if (n < 3 || n > 8) throw new ArgumentOutOfRangeException(nameof(n), n, "N in 3..8 for the dense sector survivor");
        if (q <= 0) throw new ArgumentOutOfRangeException(nameof(q), q, "Q must be > 0");
        N = n; Q = q;
    }

    /// <summary>The global survivor (smallest dissipation gap) over the candidate low-light blocks: the
    /// diagonal (p,p) for p=1..N-1 (the interior/incompleteness coherences) plus the (0,1) odd band edge.
    /// Returns (gap, pCol, pRow, n_XY = gap/2gamma). Carbon convention J=1, gamma=1/q.</summary>
    public static (double Gap, int PCol, int PRow, double NXy) Survivor(int n, double q, TopologyKind topology)
    {
        double gamma = 1.0 / q;
        var profile = Enumerable.Repeat(gamma, n).ToArray();
        double best = double.PositiveInfinity;
        int bc = 0, br = 1;
        var cands = new List<(int PCol, int PRow)>();
        for (int p = 1; p < n; p++) cands.Add((p, p));
        cands.Add((0, 1));
        foreach (var (pc, pr) in cands)
        {
            double gap = SectorReductionWitness.SectorSlowest(n, Qh, profile, pc, pr, topology);
            if (gap > KernelTol && gap < best) { best = gap; bc = pc; br = pr; }
        }
        return (best, bc, br, best / (2.0 * gamma));
    }

    private static string Kind(int n, int pc, int pr)
    {
        if (pc == 0 && pr == 1) return "band-edge (above the handover)";
        if (pc == pr && pc >= 2 && pc <= n - 2) return "INTERIOR = incompleteness (the survivor)";
        if (pc == pr && (pc == 1 || pc == n - 1)) return "BOUNDARY (the central-spin counterexample)";
        return "?";
    }

    private InspectableNode TheWhereNode()
    {
        var kids = new List<IInspectable>();
        foreach (var topo in new[] { TopologyKind.Chain, TopologyKind.Ring, TopologyKind.Star })
        {
            var (gap, pc, pr, nxy) = Survivor(N, Q, topo);
            kids.Add(new InspectableNode($"{topo}",
                summary: $"survivor sector ({pc},{pr}), <n_XY>={nxy.ToString("0.000", Inv)}, gap={gap.ToString("0.000", Inv)} " +
                         $"-> {Kind(N, pc, pr)}"));
        }
        return new InspectableNode($"where the survivor lives (N={N}, Q={Q.ToString("0.##", Inv)})",
            summary: "dispersive chain/ring -> the interior incompleteness coherence; the central-spin star -> " +
                     "the boundary hub coherence (no momentum mode). The label is physical, not the graph type.",
            children: kids);
    }

    private InspectableNode TheScalingNode()
    {
        var kids = new List<IInspectable>();
        foreach (int n in new[] { 4, 5, 6 })
        {
            double cNxy = Survivor(n, Q, TopologyKind.Chain).NXy;
            double rNxy = Survivor(n, Q, TopologyKind.Ring).NXy;
            double ratio = cNxy > 1e-9 ? rNxy / cNxy : double.NaN;
            kids.Add(new InspectableNode($"N={n}",
                summary: $"chain <n_XY>={cNxy.ToString("0.0000", Inv)}, ring <n_XY>={rNxy.ToString("0.0000", Inv)}, " +
                         $"ring/chain={ratio.ToString("0.00", Inv)}"));
        }
        return new InspectableNode($"lifetime scaling (Q={Q.ToString("0.##", Inv)})",
            summary: "<n_XY> ~ c*Q^2/N^2 (the magnon-admixture inheritance); ring/chain -> 4 (the cyclic-vs-open " +
                     "k_min^2 ratio, model-independent - the SAME 4x CHAIN_GAP reports for Heisenberg). A separate " +
                     "1/N^2 inheritance from the Pi2 dyadic CONSTANT ladder.",
            children: kids);
    }

    public string DisplayName => $"IncompletenessSurvivorWitness (N={N}, Q={Q.ToString("0.##", Inv)})";

    public string Summary =>
        "the dynamic survival probe: on dispersive extended matter (chain, ring) the longest-lived mode is the " +
        "interior incompleteness (C=0.5) coherence; the hub-localized central-spin star is the boundary " +
        "counterexample. Lifetime <n_XY> ~ Q^2/N^2, ring/chain -> 4. Reuses SectorReductionWitness.SectorSlowest.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheWhereNode();
            yield return TheScalingNode();
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
