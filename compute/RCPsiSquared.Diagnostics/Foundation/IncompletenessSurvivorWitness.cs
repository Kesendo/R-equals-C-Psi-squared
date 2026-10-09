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
/// carbon work), spin chains, the Grotthuss proton wire (water). k_min = pi/N, the lowest Neumann cosine of
/// the population diffusion (PROOF_WEIGHT1_DEGENERACY, the Zeno end).</item>
/// <item>RING = extended <b>dispersive</b> 1D periodic: aromatic rings (benzene, the Frost circle),
/// light-harvesting macrocycles. k_min = 2pi/N (2x the chain's k, 4x its k^2).</item>
/// <item>STAR = hub-spoke <b>NON-dispersive</b> central-spin: NV centre / quantum dot / donor-spin bath
/// (Bortz-Stolze), the optical-cavity point-focus (STAR_CONFOCAL_LIMIT), the mediator.</item>
/// </list></para>
///
/// <para><b>Finding (verified):</b> on the even RING the survivor is the INTERIOR (incompleteness)
/// coherence of the even fillings, every one of which carries the (2,2) block's level ((2,2)/(N-2,N-2) at N = 6, the
/// half filling among them when 4 divides N), the wrapped string's split; the odd ring ties every filling, read at every Q
/// tried (PROOF_WEIGHT1_DEGENERACY: a pi flux through it changes no block's spectrum, and at the Zeno end the ladder
/// carries the tie wherever each block's slowest level is the lifted pair, read at N = 5 to 9); the open XY chain is
/// filling-degenerate (every filling ties, measured to rounding here and proved at the Zeno end by
/// PROOF_CODIM1_BY_ADDITIVITY's ladder with Theorem E (c)), the ZZ term picking the half
/// filling, or the central pair at odd N (PROOF_WEIGHT1_DEGENERACY, Theorem E (d) at the Zeno end); the central-spin STAR is the COUNTEREXAMPLE
/// (boundary (1,1); with ZZ its hub hops are detuned by the other arms' imbalance, the Zeno end's detuning form
/// with the hops' own share read on the stars N = 4..7; in XY the boundary wins only at order (J/gamma)^6 at the
/// Zeno end, consistent with the spins' exchange statistics, free fermions on the same star tying every filling). Lifetime &lt;n_XY&gt; ~ Q^2/N^2, ring/chain -&gt; 4 (cyclic-vs-open k_min^2,
/// model-independent). Reuses <see cref="SectorReductionWitness.SectorSlowest"/> (no full 4^N).</para>
///
/// <para>Convention (carbon, XY/free-fermion, no ZZ): J=1, gamma=1/Q. SectorReductionWitness takes
/// canonical Hamiltonian Q=1 and builds H=(1/2)*(XX+YY), reproducing J=1; rate = -2*gamma*&lt;n_XY&gt;.
/// The python twin is <c>simulations/carbon/incompleteness_survivor.py</c>.</para></summary>
public sealed class IncompletenessSurvivorWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double HamiltonianQ = 1.0; // H=(Q/2)*(XX+YY), Q=1 reproduces carbon J=1
    private const double KernelTol = 1e-7;

    /// <summary>The canonical coherence-horizon EP ladder Q*(N) (F2b corollary / <see cref="CoherenceHorizonClaim"/>),
    /// indexed by N: the chain handover sits AT this at N=2,3 and just below it from N=4 (the trace dressing,
    /// <c>CoherenceHorizonWitness</c>). One shared source for the handover node and the <c>HandoverFloorClaim</c>
    /// battery, instead of two literal copies of 1.87874 / 2.88925.</summary>
    public static readonly double[] CoherenceHorizonQStar = { 0, 0, 1.0, 1.41421, 1.87874, 2.37367, 2.88925 };

    // HandoverQ is pure in (N, topology) (gamma-independent by construction, J fixed via Qh), so memoize it:
    // the battery, the witness node, and the Shared static otherwise each re-run the same (heavy, N=8) bisection.
    private static readonly System.Collections.Concurrent.ConcurrentDictionary<(int, TopologyKind), double> _handoverCache = new();

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
            double gap = SectorReductionWitness.SectorSlowest(n, HamiltonianQ, profile, pc, pr, topology);
            if (gap > KernelTol && gap < best) { best = gap; bc = pc; br = pr; }
        }
        return (best, bc, br, best / (2.0 * gamma));
    }

    /// <summary>The slowest non-kernel rate of each diagonal (p,p) filling sector, p=1..N-1 (carbon
    /// XY convention). On the OPEN chain these are degenerate to machine precision (the free-fermion
    /// dressed-magnon rate is filling-blind), so the "survivor sector" is a tie-break artifact; the ZZ
    /// of Heisenberg lifts it to the half filling, the central pair at odd N (the CHAIN_GAP result, proved at the Zeno end by
    /// PROOF_WEIGHT1_DEGENERACY Theorem E (d): a hop that changes the ZZ energy is a slower hop). On the even ring they split by
    /// parity, the even fillings slower at Q = 1 and 2 (N = 4, 6) and the odd ones at Q = 6 (N = 6); the odd ring ties
    /// every filling (read; a pi flux through it changes no block's spectrum).</summary>
    public static double[] SectorRates(int n, double q, TopologyKind topology)
    {
        double gamma = 1.0 / q;
        var profile = Enumerable.Repeat(gamma, n).ToArray();
        var rates = new double[n - 1];
        for (int p = 1; p < n; p++)
            rates[p - 1] = SectorReductionWitness.SectorSlowest(n, HamiltonianQ, profile, p, p, topology);
        return rates;
    }

    /// <summary>True iff every diagonal (p,p) filling sector shares the slowest rate (max-min &lt; tol):
    /// the open XY chain's free-fermion degeneracy, where NO single filling is the survivor (the
    /// half filling, or the central pair at odd N, is the Heisenberg/ZZ result). The even ring splits by parity, at the
    /// order gamma (J/gamma)^N at small Q, and reads as degenerate only where that split is below tol; the odd ring ties
    /// (read; a pi flux through it changes no block's spectrum).</summary>
    public static bool IsFillingDegenerate(int n, double q, TopologyKind topology, double tol = 1e-9)
    {
        var r = SectorRates(n, q, topology);
        return r.Length > 0 && r.Max() - r.Min() < tol;
    }

    /// <summary>The interior survivor's light content ⟨n_XY⟩ in the (2,2) two-fermion block at this Q
    /// (= rate/2γ). On the even RING this is the darkest interior (the V-Effect seam, p*=2; on the odd ring the level every filling shares); on the open
    /// XY CHAIN every (p,p) ties (filling-degenerate), so (2,2) equals the global interior survivor -
    /// so the (2,2) block gives the handover for BOTH topologies (N≥4), cheaply.</summary>
    public static double Interior22NXy(int n, double q, TopologyKind topology)
    {
        double gamma = 1.0 / q;
        double rate = SectorReductionWitness.SectorSlowest(n, HamiltonianQ, Enumerable.Repeat(gamma, n).ToArray(), 2, 2, topology);
        return rate / (2.0 * gamma);
    }

    /// <summary>The HANDOVER Q: where the interior incompleteness survivor brightens to the F50
    /// off-diagonal floor ⟨n_XY⟩ = 1 (the (0,1) band edge, Re = −2γ exactly), so the band edge takes
    /// over. Below it the dressed interior mode is darker (the incomplete survives); above it the F50
    /// floor wins. A closed, F50-grounded condition; spectral, γ-independent (depends only on Q=J/γ).
    /// Bisects <see cref="Interior22NXy"/> = 1; NaN if no crossing in [lo,hi]. (CHAIN handover = the coherence
    /// horizon's floor crossing by filling-degeneracy: the EP Q*(N) at N=2,3, below it from N=4 by the floor
    /// strain ((2w2−1)/c)²; RING at even N = a distinct (2,2) level crossing that grows ~linearly, the odd ring's (2,2) block
    /// holding the level every filling shares. Verifier:
    /// simulations/carbon/handover_q.py.)</summary>
    public static double HandoverQ(int n, TopologyKind topology, double lo = 0.5, double hi = 12.0, double tol = 1e-4)
    {
        if (n < 4) return double.NaN;   // (2,2) is interior only for N≥4
        if (_handoverCache.TryGetValue((n, topology), out var hit)) return hit;
        double F(double q) => Interior22NXy(n, q, topology) - 1.0;
        double result;
        if (F(lo) > 0 || F(hi) < 0)
            result = double.NaN;
        else
        {
            for (int it = 0; it < 60 && hi - lo > tol; it++)
            {
                double m = 0.5 * (lo + hi);
                if (F(m) < 0) lo = m; else hi = m;
            }
            result = 0.5 * (lo + hi);
        }
        return _handoverCache[(n, topology)] = result;
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
            bool partsBelowTol = (topo == TopologyKind.Ring && N % 2 == 0) || topo == TopologyKind.Star;
            if (partsBelowTol && IsFillingDegenerate(N, Q, topo))
            {
                var rates = SectorRates(N, Q, topo);
                string how = topo == TopologyKind.Star
                    ? "the spin XY star's fillings, which part toward the boundary at the order gamma (J/gamma)^6 (read at N = 4 to 7), part"
                    : "the even ring's two parity classes, each tied within itself and parting at the order gamma (J/gamma)^N, part";
                kids.Add(new InspectableNode($"{topo}",
                    summary: $"survivor sector ({pc},{pr}), <n_XY>={nxy.ToString("0.000", Inv)}, gap={gap.ToString("0.000", Inv)}: " +
                             $"{how} here by {(rates.Max() - rates.Min()).ToString("E1", Inv)}, below this reading's tolerance 1e-9, so " +
                             "the sector is the solver's pick among the fillings"));
            }
            else if (IsFillingDegenerate(N, Q, topo))
            {
                var rates = SectorRates(N, Q, topo);
                kids.Add(new InspectableNode($"{topo}",
                    summary: $"FILLING-DEGENERATE: every (p,p) sector shares the slowest rate {rates[0].ToString("0.000", Inv)} " +
                             $"(spread {(rates.Max() - rates.Min()).ToString("E1", Inv)}) -> the free-fermion XY {(topo == TopologyKind.Ring ? "odd ring" : topo.ToString().ToLowerInvariant())} has NO " +
                             "unique survivor sector. The half filling, or the central pair at odd N, is the HEISENBERG (ZZ) result " +
                             "(CHAIN_GAP); the ZZ lifts this degeneracy. " +
                             ((pc, pr) == (0, 1)
                                 ? $"The global survivor is the (0,1) band edge past its handover, <n_XY>={nxy.ToString("0.000", Inv)}, gap={gap.ToString("0.000", Inv)}."
                                 : $"Shared <n_XY>={nxy.ToString("0.000", Inv)}, gap={gap.ToString("0.000", Inv)}.")));
            }
            else
                kids.Add(new InspectableNode($"{topo}",
                    summary: $"survivor sector ({pc},{pr}), <n_XY>={nxy.ToString("0.000", Inv)}, gap={gap.ToString("0.000", Inv)} " +
                             $"-> {Kind(N, pc, pr)}"));
        }
        return new InspectableNode($"where the survivor lives (N={N}, Q={Q.ToString("0.##", Inv)})",
            summary: "the dispersive RING at even N puts the survivor in its even fillings, every one of which carries the (2,2) " +
                     "block's level ((2,2)/(N-2,N-2) at N = 6, the half filling among them when 4 divides N; incompleteness), the " +
                     "wrapped string's split, and at odd N ties every filling (read; a pi flux through it changes no block's spectrum); " +
                     "the open XY CHAIN is filling-degenerate (no unique sector, measured to rounding, at the Zeno end " +
                     "at every order by the ladder of PROOF_CODIM1_BY_ADDITIVITY; a ZZ term picks the half filling, the central pair at odd N, " +
                     "PROOF_WEIGHT1_DEGENERACY Theorem E (d)); the hub-spoke STAR sits at the boundary (with ZZ by the " +
                     "other arms' imbalance, read; in XY at the order (J/gamma)^6, consistent with the exchange statistics). The label is physical.",
            children: kids);
    }

    private InspectableNode TheScalingNode()
    {
        var kids = new List<IInspectable>();
        foreach (int n in new[] { 4, 5, 6 })
        {
            double cNxy = Survivor(n, Q, TopologyKind.Chain).NXy;
            var ringSurvivor = Survivor(n, Q, TopologyKind.Ring);
            double rNxy = ringSurvivor.NXy;
            double ratio = cNxy > 1e-9 ? rNxy / cNxy : double.NaN;
            string ringEdge = (ringSurvivor.PCol, ringSurvivor.PRow) == (0, 1) ? " (the (0,1) band edge, past its handover)" : "";
            kids.Add(new InspectableNode($"N={n}",
                summary: $"chain <n_XY>={cNxy.ToString("0.0000", Inv)}, ring <n_XY>={rNxy.ToString("0.0000", Inv)}{ringEdge}, " +
                         $"ring/chain={ratio.ToString("0.00", Inv)}"));
        }
        return new InspectableNode($"lifetime scaling (Q={Q.ToString("0.##", Inv)})",
            summary: "<n_XY> ~ c*Q^2/N^2 (the magnon-admixture inheritance); ring/chain -> 4 (the cyclic-vs-open " +
                     "k_min^2 ratio, model-independent - the SAME 4x CHAIN_GAP reports for Heisenberg). A separate " +
                     "1/N^2 inheritance from the Pi2 dyadic CONSTANT ladder.",
            children: kids);
    }

    /// <summary>The handover Q across topologies: the chain handover is the coherence horizon's floor crossing
    /// (filling-degeneracy), equal to the EP Q*(N) at N=2,3 and below it from N=4 by the trace dressing; the ring
    /// is a distinct (2,2) free-fermion level crossing that GROWS ~linearly (not saturating). N-independent of
    /// the witness's own N (shows the ladder).</summary>
    private InspectableNode TheHandoverNode()
    {
        var kids = new List<IInspectable>();
        foreach (int n in new[] { 4, 5, 6 })
        {
            double qh = HandoverQ(n, TopologyKind.Chain);
            kids.Add(new InspectableNode($"chain N={n} (below the EP Q*(N) by the trace dressing)",
                summary: $"handover Q={qh.ToString("0.0000", Inv)} vs the EP Q*({n})=" +
                         $"{CoherenceHorizonQStar[n].ToString("0.0000", Inv)} (filling-degenerate; gap {(CoherenceHorizonQStar[n] - qh).ToString("+0.0000;-0.0000", Inv)} " +
                         "= ((2w2-1)/c)^2, the trace dressing: zero only at the clean-2x2 N=2,3 where the EP sits on the floor)"));
        }
        foreach (int n in new[] { 6, 8 })
        {
            double qh = HandoverQ(n, TopologyKind.Ring);
            kids.Add(new InspectableNode($"ring N={n} ((2,2) seam)",
                summary: $"handover Q={qh.ToString("0.0000", Inv)} (the level of the (2,2) block, which every even filling carries, a LEVEL " +
                         $"CROSSING; c_eff=(N/Qh)^2={((n / qh) * (n / qh)).ToString("0.00", Inv)}, climbs toward 4pi^2/3=13.16 => " +
                         $"slope sqrt3/(2pi)~0.276 derived, PROOF_RING_HANDOVER_SLOPE)"));
        }
        return new InspectableNode("the handover Q (the incomplete meets the F50 floor)",
            summary: "the diagonal (p,p) survivor brightens with Q until ⟨n_XY⟩ reaches the F50 OFF-diagonal floor =1 " +
                     "(the (0,1) band edge / Uhr 1, Re=-2γ exactly), where the band edge takes over: a closed, F50-grounded " +
                     "condition (spectral, depends only on Q=J/γ). CHAIN: filling-degenerate, so the handover is the " +
                     "coherence horizon's floor crossing: the EP Q*(N) at N=2,3, below it from N=4 by the trace dressing. RING at even N: a distinct 2-excitation (2,2)/(N-2,N-2) level (NOT half-filling at N = 6; every even filling " +
                     "carries it, the half filling too when 4 divides N; at odd N every filling ties, read) free-fermion LEVEL CROSSING, asymptotic slope sqrt3/(2pi)~0.276 DERIVED (PROOF_RING_HANDOVER_SLOPE, " +
                     "reviewed 2026-07-19; the <n_XY> = 1 sibling of Q*, ratio sqrt3/2; " +
                     "c_eff climbs toward 4pi^2/3=13.16); NOT " +
                     "co-located with the ring SE-EP (curves cross near N~10; benzene's 2.0-vs-1.609 split is small-N). " +
                     "Verifier simulations/carbon/handover_q.py; the F50 floor = F50WeightOneDegeneracyPi2Inheritance.",
            children: kids);
    }

    public string DisplayName => $"IncompletenessSurvivorWitness (N={N}, Q={Q.ToString("0.##", Inv)})";

    public string Summary =>
        "the dynamic survival probe: the dispersive RING at even N puts the longest-lived mode in its even fillings' " +
        "interior coherence ((2,2)/(N-2,N-2) at N = 6; incompleteness); the odd ring (read) and the open XY CHAIN are filling-degenerate " +
        "(no unique survivor sector - the half filling, the central pair at odd N, is the Heisenberg/ZZ result, CHAIN_GAP); the central-spin STAR is the boundary counterexample. " +
        "Lifetime <n_XY> ~ Q^2/N^2, ring/chain -> 4. Reuses SectorReductionWitness.SectorSlowest. " +
        "Sector overview: inspect --root blockspectrum (this zooms the diagonal (p,p) interior sectors — the " +
        "(2,2) two-excitation {0,2}, half-filling at N=4; DISTINCT from horizon's (1,1) single-excitation {0,2}).";

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
