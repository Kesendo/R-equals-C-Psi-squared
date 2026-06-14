using System.Globalization;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The dynamic survival probe, live (<c>inspect --root survivor</c>): WHERE the longest-lived
/// dissipative mode lives across the three physically-grounded topologies, and HOW its lifetime scales
/// with N. The thread-(b) companion to the <c>survival_incompleteness_mirror</c> arc; the C# witness for
/// <see cref="SurvivalIncompletenessMirrorClaim"/>.
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
/// INTERIOR (incompleteness, C=0.5) coherence - the open XY chain filling-degenerate; the ring at the even
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

    /// <summary>The slowest non-kernel rate of each diagonal (p,p) filling sector, p=1..N-1 (carbon
    /// XY convention). On the OPEN chain these are degenerate to machine precision (the free-fermion
    /// dressed-magnon rate is filling-blind), so the "survivor sector" is a tie-break artifact; the ZZ
    /// of Heisenberg lifts it to the dead-centre (the CHAIN_GAP result). On the ring they split by
    /// parity (the even fillings are slower).</summary>
    public static double[] SectorRates(int n, double q, TopologyKind topology)
    {
        double gamma = 1.0 / q;
        var profile = Enumerable.Repeat(gamma, n).ToArray();
        var rates = new double[n - 1];
        for (int p = 1; p < n; p++)
            rates[p - 1] = SectorReductionWitness.SectorSlowest(n, Qh, profile, p, p, topology);
        return rates;
    }

    /// <summary>True iff every diagonal (p,p) filling sector shares the slowest rate (max-min &lt; tol):
    /// the open XY chain's free-fermion degeneracy, where NO single filling is the survivor (the
    /// dead-centre is the Heisenberg/ZZ result). The ring splits by parity, so it is not degenerate.</summary>
    public static bool IsFillingDegenerate(int n, double q, TopologyKind topology, double tol = 1e-9)
    {
        var r = SectorRates(n, q, topology);
        return r.Length > 0 && r.Max() - r.Min() < tol;
    }

    /// <summary>The interior survivor's darkness ⟨n_XY⟩ in the (2,2) two-fermion block at this Q
    /// (= rate/2γ). On the RING this is the darkest interior (the V-Effect seam, p*=2); on the open
    /// XY CHAIN every (p,p) ties (filling-degenerate), so (2,2) equals the global interior survivor -
    /// so the (2,2) block gives the handover for BOTH topologies (N≥4), cheaply.</summary>
    public static double Interior22NXy(int n, double q, TopologyKind topology)
    {
        double gamma = 1.0 / q;
        double rate = SectorReductionWitness.SectorSlowest(n, Qh, Enumerable.Repeat(gamma, n).ToArray(), 2, 2, topology);
        return rate / (2.0 * gamma);
    }

    /// <summary>The HANDOVER Q: where the interior incompleteness survivor brightens to the F50
    /// off-diagonal floor ⟨n_XY⟩ = 1 (the (0,1) band edge, Re = −2γ exactly), so the band edge takes
    /// over. Below it the dressed interior mode is darker (the incomplete survives); above it the F50
    /// floor wins. A closed, F50-grounded condition; spectral, γ-independent (depends only on Q=J/γ).
    /// Bisects <see cref="Interior22NXy"/> = 1; NaN if no crossing in [lo,hi]. (CHAIN handover = the
    /// coherence horizon Q*(N) by filling-degeneracy; RING = a distinct (2,2) level crossing that grows
    /// ~linearly. Verifier: simulations/carbon/handover_q.py.)</summary>
    public static double HandoverQ(int n, TopologyKind topology, double lo = 0.5, double hi = 12.0, double tol = 1e-4)
    {
        if (n < 4) return double.NaN;   // (2,2) is interior only for N≥4
        double F(double q) => Interior22NXy(n, q, topology) - 1.0;
        if (F(lo) > 0 || F(hi) < 0) return double.NaN;
        for (int it = 0; it < 60 && hi - lo > tol; it++)
        {
            double m = 0.5 * (lo + hi);
            if (F(m) < 0) lo = m; else hi = m;
        }
        return 0.5 * (lo + hi);
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
            if (IsFillingDegenerate(N, Q, topo))
            {
                var rates = SectorRates(N, Q, topo);
                kids.Add(new InspectableNode($"{topo}",
                    summary: $"FILLING-DEGENERATE: every (p,p) sector shares the slowest rate {rates[0].ToString("0.000", Inv)} " +
                             $"(spread {(rates.Max() - rates.Min()).ToString("E1", Inv)}) -> the free-fermion XY chain has NO " +
                             "unique survivor sector. The dead-centre (N/2,N/2) winner is the HEISENBERG (ZZ) result " +
                             $"(CHAIN_GAP); the ZZ lifts this degeneracy. Shared <n_XY>={nxy.ToString("0.000", Inv)}, gap={gap.ToString("0.000", Inv)}."));
            }
            else
                kids.Add(new InspectableNode($"{topo}",
                    summary: $"survivor sector ({pc},{pr}), <n_XY>={nxy.ToString("0.000", Inv)}, gap={gap.ToString("0.000", Inv)} " +
                             $"-> {Kind(N, pc, pr)}"));
        }
        return new InspectableNode($"where the survivor lives (N={N}, Q={Q.ToString("0.##", Inv)})",
            summary: "the dispersive RING puts the survivor at the even off-centre interior (incompleteness); the open XY " +
                     "CHAIN is filling-degenerate (no unique sector, the ZZ-free accident; Heisenberg's ZZ pins the dead-centre); " +
                     "the hub-spoke STAR sits at the boundary (no momentum mode). The label is physical.",
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

    /// <summary>The handover Q across topologies: the chain handover IS the coherence horizon Q*(N)
    /// (filling-degeneracy), the ring is a distinct (2,2) free-fermion level crossing that GROWS
    /// ~linearly (not saturating). N-independent of the witness's own N (shows the ladder).</summary>
    private InspectableNode TheHandoverNode()
    {
        // Q*(N) canonical (the coherence horizon, F2b corollary / CoherenceHorizonClaim), indexed by N.
        double[] qStar = { 0, 0, 1.0, 1.41421, 1.87874, 2.37367, 2.88925 };
        var kids = new List<IInspectable>();
        foreach (int n in new[] { 4, 5, 6 })
        {
            double qh = HandoverQ(n, TopologyKind.Chain);
            kids.Add(new InspectableNode($"chain N={n} (= Q*(N))",
                summary: $"handover Q={qh.ToString("0.0000", Inv)} = the coherence horizon Q*({n})=" +
                         $"{qStar[n].ToString("0.0000", Inv)} (filling-degenerate; gap {(qStar[n] - qh).ToString("+0.0000;-0.0000", Inv)} " +
                         "= the trace dressing, exact only at the clean-2x2 N=2,3 - a coalescence/EP)"));
        }
        foreach (int n in new[] { 6, 8 })
        {
            double qh = HandoverQ(n, TopologyKind.Ring);
            kids.Add(new InspectableNode($"ring N={n} ((2,2) seam)",
                summary: $"handover Q={qh.ToString("0.0000", Inv)} (the (2,2) free-fermion V-Effect seam, a LEVEL " +
                         $"CROSSING; c_eff=(N/Qh)^2={((n / qh) * (n / qh)).ToString("0.00", Inv)} - flat in N => linear growth)"));
        }
        return new InspectableNode("the handover Q (the incomplete meets the F50 floor)",
            summary: "the diagonal (p,p) survivor brightens with Q until ⟨n_XY⟩ reaches the F50 OFF-diagonal floor =1 " +
                     "(the (0,1) band edge / Uhr 1, Re=-2γ exactly), where the band edge takes over: a closed, F50-grounded " +
                     "condition (spectral, depends only on Q=J/γ). CHAIN: filling-degenerate, so the handover IS the " +
                     "coherence horizon Q*(N) (a coalescence/EP). RING: a distinct (2,2) free-fermion LEVEL CROSSING, growing " +
                     "~linearly Q_h~0.29N (c_eff~12 flat, ~4x chain darkness => ~half the slope), NOT saturating; NOT " +
                     "co-located with the ring SE-EP (curves cross near N~10; benzene's 2.0-vs-1.609 split is small-N). " +
                     "Verifier simulations/carbon/handover_q.py; the F50 floor = F50WeightOneDegeneracyPi2Inheritance.",
            children: kids);
    }

    public string DisplayName => $"IncompletenessSurvivorWitness (N={N}, Q={Q.ToString("0.##", Inv)})";

    public string Summary =>
        "the dynamic survival probe: the dispersive RING puts the longest-lived mode at the even off-centre interior " +
        "(incompleteness, C=0.5) coherence; the open XY CHAIN is filling-degenerate (no unique survivor sector - the " +
        "dead-centre is the Heisenberg/ZZ result, CHAIN_GAP); the central-spin STAR is the boundary counterexample. " +
        "Lifetime <n_XY> ~ Q^2/N^2, ring/chain -> 4. Reuses SectorReductionWitness.SectorSlowest.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheWhereNode();
            yield return TheScalingNode();
            yield return TheHandoverNode();
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
