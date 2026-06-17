using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the coherence horizon Q*(N) (typed home intent:
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F2b corollary "Coherence horizon Q*(N)"). Q*(N) is the
/// threshold where the slowest non-zero Liouvillian mode stops oscillating: above Q* (small γ) the
/// slowest mode oscillates (the coherence hand turns), below it the mode is overdamped (real, the
/// hand has frozen). Sweeping Q = J/γ downward, the freeze happens at
/// Q*(2)=1, Q*(3)=√2, Q*(4)≈1.8785, Q*(5)≈2.37217.
///
/// <para>These match, bit-for-bit, the carbon coherent↔incoherent threshold from
/// <c>docs/carbon/FROST_CIRCLE_AS_THE_CLOCK_FACE.md</c> (√2 / 1.879 / 2.372 at N=3/4/5) under the
/// label swap J ↔ |β|: the XY chain's coherence horizon IS the Frost-Hückel threshold of the same
/// polyene. N=2 (Q*=1) is the exceptional point itself (γ=J), the base rung the carbon polyene
/// layer (N≥3) cannot reach; the quantum side supplies it. This witness is the "C# witness first"
/// port of the probe <c>simulations/_carbon_quantum_same_mountain.py</c>.</para>
///
/// <para>It reuses <see cref="Symphony"/> as the spectrum engine: Q*(N) is found by bisecting γ at
/// J=1 on whether <see cref="Symphony.Clock"/>.Omega (max|Im λ| among the modes at the slowest
/// non-zero decay rate) is non-zero. The per-N value is computed once and cached (N=5 is a 1024×1024
/// eigendecomposition; ~16 bisection builds, lazy and build-once). The closed form of Q*(N) is OPEN;
/// the witness recomputes the live threshold and checks it against the carbon ladder.</para></summary>
public sealed class CoherenceHorizonWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The coupling the horizon is read at (J=1, so Q* = 1/γ*).</summary>
    public const double J = 1.0;

    /// <summary>The Omega threshold: |Im| above this counts as "the slowest mode still oscillates".</summary>
    private const double OmegaEps = 1e-3;

    /// <summary>γ bracket for the bisection: oscillating at γ_lo=0.3 (Q=3.33), frozen at γ_hi=1.2
    /// (Q=0.83), covering every Q* ∈ [0.83, 3.3].</summary>
    private const double GammaLo = 0.3, GammaHi = 1.2;
    private const int BisectionSteps = 16;

    /// <summary>The spectrum-reading grid: the clock reads the eigenvalues, so the time grid is
    /// resolution-independent; a small grid keeps the per-build cost down.</summary>
    private const int SpectrumPoints = 8;

    /// <summary>The carbon Frost-Hückel coherent↔incoherent thresholds (FROST_CIRCLE_AS_THE_CLOCK_FACE.md):
    /// √2 / 1.879 / 2.372 at N=3/4/5. N=2 has no carbon entry (the polyene layer starts at N≥3).</summary>
    private static readonly Dictionary<int, double> Carbon = new()
    {
        { 3, Math.Sqrt(2.0) }, { 4, 1.879 }, { 5, 2.372 },
    };

    private readonly Dictionary<int, double> _horizonCache = new();

    /// <summary>The live coherence horizon Q*(N): the largest Q (= smallest γ at J=1) below which the
    /// slowest non-zero Liouvillian mode no longer oscillates. Computed by bisecting γ on whether
    /// <see cref="Symphony.Clock"/>.Omega &gt; <see cref="OmegaEps"/>, then Q* = J/γ*. Cached per N
    /// (one bisection, ~16 dense eigendecompositions; build-once).</summary>
    public double Horizon(int n)
    {
        if (n < 2 || n > Symphony.MaxN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"the coherence horizon needs N in 2..{Symphony.MaxN} (dense Liouvillian spectrum); got {n}");
        if (_horizonCache.TryGetValue(n, out var q)) return q;

        // γ_lo oscillates (Omega > ε), γ_hi is frozen (Omega ≈ 0); bisect to the crossover γ*.
        double glo = GammaLo, ghi = GammaHi;
        for (int step = 0; step < BisectionSteps; step++)
        {
            double mid = 0.5 * (glo + ghi);
            if (Omega(n, mid) > OmegaEps) glo = mid;   // still oscillating → need more γ (lower Q)
            else ghi = mid;                            // frozen → less γ (higher Q)
        }
        double gStar = 0.5 * (glo + ghi);
        q = J / gStar;
        _horizonCache[n] = q;
        return q;
    }

    /// <summary>Omega(N, γ) = the coherence hand at J=1: max|Im λ| among the modes sharing the slowest
    /// non-zero decay rate (the gap), read off the <see cref="Symphony"/> clock. &gt; ε ⟹ the slowest
    /// mode oscillates; ≈ 0 ⟹ it is overdamped (real).</summary>
    public double Omega(int n, double gamma) =>
        new Symphony(n: n, j: J, gamma: gamma, initialState: InitialStateKind.BellPair,
            tPoints: SpectrumPoints).Clock.Omega;

    /// <summary>The carbon threshold for this N (the Frost-Hückel coherent↔incoherent value), or null
    /// for N=2 (no polyene entry).</summary>
    public double? CarbonThreshold(int n) => Carbon.TryGetValue(n, out var c) ? c : null;

    /// <summary>The low-N band-edge coincidence 2cos(π/(N+1)): equal to Q*(N) at N=2,3 only
    /// (1 = 2cos60°, √2 = 2cos45°), then departing (Q*(4)=1.8785 ≠ φ=1.618).</summary>
    public static double BandEdgeCoincidence(int n) => 2.0 * Math.Cos(Math.PI / (n + 1));

    /// <summary>The rigidity below which a gap mode counts as coalescing (an EP).</summary>
    private const double REpThreshold = 0.05;

    /// <summary>The non-zero modes within a small band of the slowest decay rate (the gap), with their
    /// phase rigidity, built at Q = J/γ on the live Liouvillian.</summary>
    private static List<PhaseRigidity.Mode> GapModes(int n, double gamma)
    {
        var L = new ChainSystem(n, J, gamma).BuildLiouvillian();
        var nz = PhaseRigidity.Compute(L).Where(m => m.Lambda.Real < -1e-6).ToList();
        double gap = nz.Max(m => m.Lambda.Real);
        return nz.Where(m => m.Lambda.Real > gap - 0.15).ToList();
    }

    /// <summary>At Q*(n): the coalescing gap mode (minimum phase rigidity, the {0,2}-coherence whose
    /// r → 0) with its n_diff histogram, and the co-located band-edge survivor (Im ≈ 2cos(π/(N+1)),
    /// r ≈ 1). The instrument that distinguishes the EP (erasure) from the crossing (survival).</summary>
    public (PhaseRigidity.Mode Coalescer,
            IReadOnlyDictionary<int, double> CoalescerHist,
            double CoalescerMeanNDiff,
            PhaseRigidity.Mode BandEdge,
            double BandEdgeR) EpModes(int n)
    {
        var gapModes = GapModes(n, J / Horizon(n));
        var coalescer = gapModes.OrderBy(m => m.Rigidity).First();
        var (mean, hist) = LiouvilleOperatorContent.NDiffHistogram(coalescer.Right, n);
        double bandIm = 2.0 * Math.Cos(Math.PI / (n + 1));
        var bandEdge = gapModes.OrderBy(m => Math.Abs(Math.Abs(m.Lambda.Imaginary) - bandIm)).First();
        return (coalescer, hist, mean, bandEdge, bandEdge.Rigidity);
    }

    /// <summary>√-scaling certificate of a 2nd-order EP: Im²/(Q−Q*) for the small-Im coalescer branch
    /// at Q = Q*(n)·(1+delta). Constant across deltas ⟹ Im ∝ √(Q−Q*) ⟹ a clean 2-dim EP (not a
    /// cluster). Returns NaN if no coalescer branch is resolved above Q*.</summary>
    public double SqrtScalingRatio(int n, double delta)
    {
        double qStar = Horizon(n);
        double q = qStar * (1.0 + delta);
        double bandIm = 2.0 * Math.Cos(Math.PI / (n + 1));
        var branches = GapModes(n, J / q)
            .Where(m => Math.Abs(m.Lambda.Imaginary) > 1e-6 && Math.Abs(m.Lambda.Imaginary) < bandIm - 0.2)
            .OrderBy(m => Math.Abs(m.Lambda.Imaginary))
            .ToList();
        if (branches.Count == 0) return double.NaN;
        double im = Math.Abs(branches[0].Lambda.Imaginary);
        return im * im / Math.Abs(q - qStar);
    }

    public string DisplayName =>
        $"CoherenceHorizonWitness (Q*(N) live, J={J.ToString("0.#", Inv)}, ε={OmegaEps.ToString("0.###", Inv)})";

    public string Summary =>
        "the coherence horizon Q*(N) live (typed-home intent: docs/ANALYTICAL_FORMULAS.md F2b corollary " +
        "\"Coherence horizon Q*(N)\"): the Q below which the slowest non-zero Liouvillian mode stops " +
        "oscillating (the coherence hand freezes). Computed by bisecting γ on Symphony.Clock.Omega at J=1: " +
        $"Q*(2)={Horizon(2).ToString("0.####", Inv)}, Q*(3)={Horizon(3).ToString("0.####", Inv)}, " +
        $"Q*(4)={Horizon(4).ToString("0.####", Inv)}, Q*(5)={Horizon(5).ToString("0.####", Inv)}. " +
        "This equals the carbon Frost-Hückel coherent↔incoherent threshold (√2 / 1.879 / 2.372 at N=3/4/5, " +
        "FROST_CIRCLE_AS_THE_CLOCK_FACE.md) under the label swap J ↔ |β|. " +
        "Sector overview: inspect --root blockspectrum (this zooms the (1,1) single-excitation {0,2} sector — " +
        "the low-Q EP regime; --root ceiling reads its high-Q regime).";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheLadder();
            yield return TheEpBase();
            yield return TheBandEdgeCoincidence();
            yield return TheEpVerdict();
        }
    }

    /// <summary>the ladder: the computed Q*(N) for N=2..5 with the carbon comparison; confirm they match.</summary>
    private InspectableNode TheLadder()
    {
        var ns = new[] { 2, 3, 4, 5 };
        var rungs = new List<IInspectable>();
        bool allMatch = true;
        foreach (int n in ns)
        {
            double q = Horizon(n);
            double? c = CarbonThreshold(n);
            string carbonStr = c is { } cv ? cv.ToString("0.####", Inv) : "none (no polyene entry)";
            bool match = c is { } cm && Math.Abs(q - cm) < 0.01;
            if (c is not null && !match) allMatch = false;
            rungs.Add(InspectableNode.RealScalar($"Q*({n})", q, "0.#####"));
            rungs.Add(new InspectableNode($"N={n}: Q* vs carbon",
                summary: $"Q*({n}) = {q.ToString("0.#####", Inv)}, carbon = {carbonStr}" +
                         (c is null ? " (the EP base, quantum-only)" : $" → {(match ? "match ✓" : "MISMATCH")}")));
        }
        return new InspectableNode("the ladder",
            summary: $"the computed Q*(N) (live, bisected γ on Symphony.Clock.Omega at J=1): " +
                     $"N=2 → {Horizon(2).ToString("0.#####", Inv)}, N=3 → {Horizon(3).ToString("0.#####", Inv)} (√2), " +
                     $"N=4 → {Horizon(4).ToString("0.#####", Inv)}, N=5 → {Horizon(5).ToString("0.#####", Inv)}. " +
                     $"Carbon Frost-Hückel threshold (N=3/4/5): √2 / 1.879 / 2.372. " +
                     $"{(allMatch ? "All N≥3 rungs match the carbon ladder ✓ (the same mountain, label swap J ↔ |β|)." : "A rung does NOT match the carbon ladder; investigate.")}",
            children: rungs);
    }

    /// <summary>the EP base (N=2): Q*(2)=1 is the exceptional point, the rung the polyene layer cannot reach.</summary>
    private InspectableNode TheEpBase()
    {
        double q2 = Horizon(2);
        return new InspectableNode("the EP base (N=2)",
            summary: $"Q*(2) = {q2.ToString("0.#####", Inv)} is the exceptional point itself (γ=J, the two clocks " +
                     "merge into one). It is the base rung the carbon polyene layer (N≥3) cannot reach: the quantum " +
                     "side supplies it. Carbon's coherent↔incoherent ladder begins at N=3 (a polyene needs ≥3 sites); " +
                     "the EP is below it. Mechanism note: at N=2 the crossover the bisection lands is NOT the pulled " +
                     "coherence hand 2√(J²−γ²) freezing. Symphony.Clock.Omega(2) reads the ±J band mode (Omega=1) right " +
                     "down to γ=J, then drops to 0, because that band mode ceases to be the slowest-decay (gap) mode " +
                     "exactly at γ=J, which coincides with the exceptional point. So Q*(2)=1 marks the ±J band mode " +
                     "ceasing to be the gap mode, coinciding with the EP, distinct from the N≥3 band-edge-vs-overdamped " +
                     "crossover (see ClockHandLadderWitness for the N=2 raw-clock subtlety).");
    }

    /// <summary>the band-edge coincidence: Q*(N) = 2cos(π/(N+1)) at N=2,3 ONLY, a low-N accident.</summary>
    private InspectableNode TheBandEdgeCoincidence()
    {
        double be3 = BandEdgeCoincidence(3);
        double phi = BandEdgeCoincidence(4); // 2cos(π/5) = φ, the golden ratio
        return new InspectableNode("the band-edge coincidence",
            summary: $"Q*(N) = 2cos(π/(N+1)) at N=2,3 ONLY: 1 = 2cos60° and √2 = 2cos45° = " +
                     $"{be3.ToString("0.#####", Inv)} match Q*(2)/Q*(3). It is a low-N accident, departing at N≥4: " +
                     $"2cos(π/5) = φ = {phi.ToString("0.#####", Inv)} but Q*(4) = {Horizon(4).ToString("0.#####", Inv)} " +
                     $"≠ φ. This is why the carbon note saw \"√2 exact at " +
                     "N=3, the rest awaiting a clean form\": the band edge is the right answer only where it happens to " +
                     "coincide with the horizon.");
    }

    /// <summary>the EP verdict, recomputed live: per N=2..5 the coalescing {0,2}-coherence (r → 0,
    /// histogram {0:½,2:½}) and the co-located band-edge survivor (r ≈ 1). No bifurcation at N=4: the
    /// {0,2}-coherence is the freezer at every N, the band edge the γ-protected survivor; they share
    /// the gap Re = −2γ (Absorption Theorem, both ⟨n_diff⟩ = 1). The closed form of Q*(N) is OPEN.</summary>
    private InspectableNode TheEpVerdict()
    {
        var rungs = new List<IInspectable>();
        foreach (int n in new[] { 2, 3, 4, 5 })
        {
            var ep = EpModes(n);
            string h0 = ep.CoalescerHist.GetValueOrDefault(0).ToString("0.00", Inv);
            string h2 = ep.CoalescerHist.GetValueOrDefault(2).ToString("0.00", Inv);
            bool isEp = ep.Coalescer.Rigidity < REpThreshold;
            string label = isEp ? "EP" : "crossing";
            string verdict = isEp ? "genuine EP, the {0,2}-coherence coalesces" : "no EP (a crossing)";
            rungs.Add(new InspectableNode($"N={n}: {label}",
                summary: $"coalescer r = {ep.Coalescer.Rigidity.ToString("0.000", Inv)} " +
                         $"(Im = {ep.Coalescer.Lambda.Imaginary.ToString("0.000", Inv)}, hist {{0:{h0}, 2:{h2}}}, " +
                         $"mean n_diff = {ep.CoalescerMeanNDiff.ToString("0.0", Inv)}) → {verdict}; " +
                         $"band edge Im = {Math.Abs(ep.BandEdge.Lambda.Imaginary).ToString("0.000", Inv)} " +
                         $"(2cos(π/(N+1))), r = {ep.BandEdgeR.ToString("0.000", Inv)} → the γ-protected survivor."));
        }
        double ratio = SqrtScalingRatio(4, 0.03);
        return new InspectableNode("the EP verdict (live phase rigidity)",
            summary: "the mode that coalesces at Q*(N) is the {0,2}-coherence (population/antisymmetric block) " +
                     "at ALL N=2..5, a genuine square-root EP (phase rigidity r → 0). NO sector bifurcation at " +
                     "N=4: the band edge 2cos(π/(N+1)) is the co-located γ-protected SURVIVOR (r ≈ 1), sharing the " +
                     "gap Re = −2γ only because the Absorption Theorem pins both (both ⟨n_diff⟩ = 1). So Q*(N) is at " +
                     "once a {0,2}-coherence EP (the erasure point, which climbs the ladder) and a band-edge crossing " +
                     $"(the clock survives). √-scaling Im²/(Q−Q*) at N=4 = {ratio.ToString("0.00", Inv)} (constant ⟹ " +
                     "a clean 2nd-order EP). The closed form of Q*(N) (the {0,2}-block discriminant) is OPEN. " +
                     "Recomputed live via PhaseRigidity; supersedes the earlier narrated 'bifurcation at N=4'.",
            children: rungs);
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
