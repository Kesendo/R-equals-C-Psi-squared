using System.Globalization;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The HANDOVER Q has a closed condition: the diagonal incompleteness survivor meets the
/// F50 floor (Tier1Candidate). Resolves the open loose end of the <c>survival_incompleteness_mirror</c>
/// and <c>clock_hand_ladder</c> arcs and the F2b corollary "Open remainder". The live witness is
/// <see cref="IncompletenessSurvivorWitness"/> (<c>inspect --root survivor</c>, the handover node).
///
/// <para><b>THE CONDITION (closed, F50-grounded).</b> The diagonal (p,p) "incompleteness survivor"
/// decays at −2γ·⟨n_XY⟩ with fractional ⟨n_XY⟩ &lt; 1 (the Absorption Theorem, parent
/// <see cref="AbsorptionTheoremClaim"/>), out-surviving the bare band edge - until, as Q rises, its
/// darkness reaches the F50-pinned OFF-diagonal floor ⟨n_XY⟩ = 1 (the (0,1) band edge / Uhr 1, Re = −2γ
/// exactly; parent <see cref="F50WeightOneDegeneracyPi2Inheritance"/>). That meeting is the handover.
/// Spectral (state-independent), depends only on Q = J/γ.</para>
///
/// <para><b>THE CHAIN SOLUTION = the coherence horizon Q*(N).</b> The open XY chain is filling-
/// degenerate (free-fermion/OBC), so the crossing is the single-excitation {0,2}-coherence point =
/// Q*(N) (parent <see cref="CoherenceHorizonClaim"/>), a coalescence/EP. It coincides with Q*(N)
/// exactly only at the clean-2×2 N=2,3 (a tangency) and sits just below it by the trace dressing
/// O((tr−1)²) at N≥4 (gap 0.0002/0.0015/0.0050 at N=4/5/6).</para>
///
/// <para><b>THE RING SOLUTION = a distinct (2,2) free-fermion level crossing, growing ~linearly.</b>
/// The wrap bond breaks filling-degeneracy (Fermi degeneracy at half-filling); the survivor is the
/// (2,2) two-fermion V-Effect seam (in pure XY a free-fermion dephasing mode, NOT a Hamiltonian bound
/// pair). Its handover is a frozen LEVEL CROSSING (a different sector than the SE-EP), growing
/// ~linearly Q_h ≈ N/√c_eff with c_eff ≈ 12 flat in N (≈4× the chain darkness, so ≈half the chain
/// slope, ~0.29N) - NOT saturating, faster than √N. The handover and the ring SE-EP are mechanistically
/// distinct but their values merely CROSS near N≈10 (benzene's 2.0-vs-1.609 split is a small-N feature).</para>
///
/// <para>Anchors: <c>simulations/carbon/handover_q.py</c> (self-validating) +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F2b corollary + <c>IncompletenessSurvivorWitness.cs</c>
/// (the live handover node).</para></summary>
public sealed class HandoverFloorClaim : Claim
{
    /// <summary>Typed parent: the survival law Re = −2γ·⟨n_XY⟩ that sets the survivor's darkness.</summary>
    public AbsorptionTheoremClaim Survival { get; }

    /// <summary>Typed parent: the F50 off-diagonal floor (Re = −2γ, ⟨n_XY⟩ = 1) the survivor rises to.</summary>
    public F50WeightOneDegeneracyPi2Inheritance Floor { get; }

    /// <summary>Typed parent: the chain solution - the handover IS the coherence horizon Q*(N).</summary>
    public CoherenceHorizonClaim ChainSolution { get; }

    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public HandoverFloorClaim(
        AbsorptionTheoremClaim survival,
        F50WeightOneDegeneracyPi2Inheritance floor,
        CoherenceHorizonClaim chainSolution)
        : base("The handover Q (where the incompleteness survivor stops winning and the band edge takes over) has a CLOSED, " +
               "F50-grounded condition: the diagonal (p,p) survivor decays at -2g<n_XY> (Absorption Theorem, <n_XY><1, the " +
               "incomplete out-survives) and brightens with Q until <n_XY> reaches the F50 off-diagonal floor =1 (the (0,1) band " +
               "edge, Re=-2g exactly). CHAIN: filling-degenerate, so the handover IS the coherence horizon Q*(N) (a coalescence/EP; " +
               "= Q* exactly at the clean-2x2 N=2,3, just below by trace dressing O((tr-1)^2) at N>=4). RING: a DISTINCT (2,2) " +
               "free-fermion LEVEL CROSSING, growing ~linearly Q_h~0.29N (c_eff~12 flat, ~4x chain darkness so ~half the slope), " +
               "NOT saturating; NOT co-located with the ring SE-EP (their values cross near N~10, benzene's 2.0-vs-1.609 is small-N). " +
               "Tier1Candidate: the chain identity and the mechanism are derived+live (witness bit-for-bit vs Python), the ring " +
               "growth law and the general handover proof are verified (N<=10), not proven.",
               Tier.Tier1Candidate,
               "simulations/carbon/handover_q.py (self-validating) + " +
               "docs/ANALYTICAL_FORMULAS.md F2b corollary + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/IncompletenessSurvivorWitness.cs (the handover node, inspect --root survivor)")
    {
        Survival = survival ?? throw new ArgumentNullException(nameof(survival));
        Floor = floor ?? throw new ArgumentNullException(nameof(floor));
        ChainSolution = chainSolution ?? throw new ArgumentNullException(nameof(chainSolution));
        Cases = BuildBattery();
    }

    private static string Fmt(double v) => v.ToString("0.####", CultureInfo.InvariantCulture);

    public override string DisplayName =>
        "The handover Q = the F50-floor condition: chain = Q*(N), ring = a distinct (2,2) level crossing (Tier1Candidate)";

    public override string Summary =>
        "the handover (the incompleteness survivor's darkness rises to the F50 floor <n_XY>=1) is closed: the CHAIN handover IS " +
        "the coherence horizon Q*(N), the RING is a distinct (2,2) free-fermion level crossing growing ~0.29N (not saturating); " +
        $"the same band-edge floor governs the XXZ Delta-axis handover too (cross-axis universal, verified); " +
        $"{PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the condition: the survivor meets the F50 floor",
                summary: "the diagonal (p,p) survivor decays at -2g<n_XY> (Absorption Theorem, a_0) with fractional <n_XY><1, " +
                         "out-surviving the bare band edge; the handover is where it brightens to the F50 OFF-diagonal floor " +
                         "<n_XY>=1 (the (0,1) band edge, Re=-2g). Closed, spectral, depends only on Q=J/g.");
            yield return new InspectableNode("chain = Q*(N) (a coalescence/EP)",
                summary: "the open XY chain is filling-degenerate, so the handover is the single-excitation {0,2}-coherence point " +
                         "= the coherence horizon Q*(N); = Q* exactly only at the clean-2x2 N=2,3, just below by trace dressing at N>=4.");
            yield return new InspectableNode("ring = a distinct (2,2) level crossing, growing",
                summary: "the wrap bond breaks filling-degeneracy; the (2,2) two-fermion V-Effect seam (a free-fermion dephasing " +
                         "mode, not a bound pair) hands over by a frozen level crossing, growing ~linearly Q_h~0.29N (c_eff~12 flat), " +
                         "NOT saturating; its values cross the ring SE-EP near N~10 (benzene's 2.0-vs-1.609 split is small-N). " +
                         "Live: inspect --root survivor (the handover node).");
            yield return new InspectableNode("cross-axis universality: the same floor governs the XXZ Delta-axis",
                summary: "CONFIRMED 2026-06-14 (simulations/xxz_handover_unification.py; arc xxz_axis_handover): the band-edge " +
                         "floor (darkness=1 = the Absorption-Theorem 2*gamma) is the handover threshold on BOTH the dephasing (Q) " +
                         "axis AND the Hamiltonian-anisotropy (Delta) axis. Walking H=J(XX+YY)+Delta*ZZ, the Lebensader (the " +
                         "dead-centre I/Z survivor) overtakes the band edge exactly where its darkness crosses 1 (Delta*(4)=1.618=phi, " +
                         "Delta*(5)=1.525), a LEVEL CROSSING (frozen Lebensader meets oscillating band edge, the ring family). The " +
                         "band edge stays at 2*gamma for ALL Delta - model-independent, since |vac><magnon| is an eigenoperator of " +
                         "[H,.] (the ZZ shifts only Im); F50's 2N count breaks for Delta!=1 but the FLOOR persists.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return Survival;       // typed parent: the -2g<n_XY> rate
            yield return Floor;          // typed parent: the F50 floor =1
            yield return ChainSolution;  // typed parent: the chain = Q*(N)
        }
    }

    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        var cases = new List<BatteryCase>();

        // 1. THE FLOOR: the (0,1) band edge sits at <n_XY> = 1 (Re=-2g) exactly (the F50 floor the survivor meets).
        double gamma = 1.0 / 2.0;
        double floor = SectorReductionWitness.SectorSlowest(6, 0.5, Enumerable.Repeat(gamma, 6).ToArray(), 0, 1, TopologyKind.Chain) / (2.0 * gamma);
        cases.Add(new BatteryCase("the F50 floor: the (0,1) band edge <n_XY> = 1 (Re=-2g)",
            $"N=6, Q=2: (0,1) band-edge <n_XY> = {Fmt(floor)}", "1", Math.Abs(floor - 1.0) < 1e-6 ? "1" : Fmt(floor)));

        // 2. CHAIN = Q*(N): the chain handover reproduces the coherence horizon Q*(4)=1.87874.
        double qhChain4 = IncompletenessSurvivorWitness.HandoverQ(4, TopologyKind.Chain);
        cases.Add(new BatteryCase("chain handover = the coherence horizon Q*(4)",
            $"HandoverQ(4, chain) = {Fmt(qhChain4)} vs Q*(4) = 1.87874", "= Q*(4)",
            Math.Abs(qhChain4 - 1.87874) < 0.01 ? "= Q*(4)" : Fmt(qhChain4)));

        // 3. CHAIN below Q* by the trace dressing (N=6): 0 < Q*(6) - handover < 0.02.
        double qhChain6 = IncompletenessSurvivorWitness.HandoverQ(6, TopologyKind.Chain);
        double gap6 = 2.88925 - qhChain6;
        cases.Add(new BatteryCase("chain handover just below Q*(6) (the trace dressing)",
            $"Q*(6) - HandoverQ(6, chain) = {Fmt(gap6)}", "in (0, 0.02)",
            gap6 > 0 && gap6 < 0.02 ? "in (0, 0.02)" : Fmt(gap6)));

        // 4. RING is DISTINCT (the (2,2) V-Effect seam): benzene ring handover != chain handover.
        double qhRing6 = IncompletenessSurvivorWitness.HandoverQ(6, TopologyKind.Ring);
        double sep = Math.Abs(qhRing6 - qhChain6);
        cases.Add(new BatteryCase("ring handover is the distinct (2,2) V-Effect seam (benzene)",
            $"|HandoverQ(6, ring) {Fmt(qhRing6)} - HandoverQ(6, chain) {Fmt(qhChain6)}| = {Fmt(sep)}", "> 0.5",
            sep > 0.5 ? "> 0.5" : Fmt(sep)));

        // 5. RING GROWS (not saturating): HandoverQ(8, ring) - HandoverQ(6, ring) > 0.2.
        double qhRing8 = IncompletenessSurvivorWitness.HandoverQ(8, TopologyKind.Ring);
        double grow = qhRing8 - qhRing6;
        cases.Add(new BatteryCase("ring handover GROWS with N (not saturating)",
            $"HandoverQ(8, ring) {Fmt(qhRing8)} - HandoverQ(6, ring) {Fmt(qhRing6)} = {Fmt(grow)}", "> 0.2",
            grow > 0.2 ? "> 0.2" : Fmt(grow)));

        return cases;
    }

    /// <summary>The shared instance with its three parents built fresh from the Pi2-Foundation roots,
    /// mirroring <see cref="SurvivalIncompletenessMirrorClaim.Shared"/>: lets metadata be read without
    /// the full registry.</summary>
    public static HandoverFloorClaim Shared { get; } =
        new HandoverFloorClaim(
            new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim()),
            new F50WeightOneDegeneracyPi2Inheritance(new Pi2DyadicLadderClaim()),
            CoherenceHorizonClaim.Build());
}
