using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

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

    public string DisplayName =>
        $"CoherenceHorizonWitness (Q*(N) live, J={J.ToString("0.#", Inv)}, ε={OmegaEps.ToString("0.###", Inv)})";

    public string Summary =>
        "the coherence horizon Q*(N) live (typed-home intent: docs/ANALYTICAL_FORMULAS.md F2b corollary " +
        "\"Coherence horizon Q*(N)\"): the Q below which the slowest non-zero Liouvillian mode stops " +
        "oscillating (the coherence hand freezes). Computed by bisecting γ on Symphony.Clock.Omega at J=1: " +
        $"Q*(2)={Horizon(2).ToString("0.####", Inv)}, Q*(3)={Horizon(3).ToString("0.####", Inv)}, " +
        $"Q*(4)={Horizon(4).ToString("0.####", Inv)}, Q*(5)={Horizon(5).ToString("0.####", Inv)}. " +
        "This equals the carbon Frost-Hückel coherent↔incoherent threshold (√2 / 1.879 / 2.372 at N=3/4/5, " +
        "FROST_CIRCLE_AS_THE_CLOCK_FACE.md) under the label swap J ↔ |β|.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheLadder();
            yield return TheEpBase();
            yield return TheBandEdgeCoincidence();
            yield return WhatIsOpen();
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

    /// <summary>what is open: the closed form of Q*(N) is OPEN; the freezing mode is half-lit, not dark.</summary>
    private InspectableNode WhatIsOpen()
    {
        return new InspectableNode("what is open",
            summary: "the closed form of Q*(N) is OPEN (the band edge 2cos(π/(N+1)) fits only N=2,3). The mode that " +
                     "freezes just below Q* is overdamped (real) AND half-lit: its light content ⟨n_XY⟩ ≈ ½, NOT the " +
                     "dark {I,Z} sector (⟨n_XY⟩ ≈ 0). The probe _carbon_quantum_same_mountain.py reads this directly. " +
                     "No mechanism is asserted here: the witness reports the live threshold and its carbon coincidence, " +
                     "not why the freezing mode sits where it does.");
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
