using System;
using System.Collections.Generic;
using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The topology generalization of <see cref="ClockHandLadderWitness"/> (typed home:
/// <c>TopologyBandEdgeClaim</c>). The single-excitation band edge = J × the hopping-graph adjacency
/// spectral radius ρ (chain 2cos(π/(N+1)), star √(N−1), ring 2). Whether the clock reads it as ω_mem (the
/// band edge is the strict gap mode) is topology-specific: a Q-horizon for chain / star N≤5 / odd rings, a
/// structural ceiling for star N≥6, and a co-occupied floor for ring N=4. The witness reads ω_mem at two
/// Q points (the operating Q and ×50) so it can label horizon (climbing) vs co-occupied honestly; N=6 (the
/// cleanest ceiling) is exercised in the test suite (4^6 dense, too slow for the live map). The Re=−2γ floor
/// is the Absorption Theorem; this witness reads the Im side (ω_mem) and where it is the gap mode.
///
/// <para><b>Clock engine</b>: the gap-mode reading is the same Takt/coherence clock <see cref="Symphony"/>
/// exposes — the slowest nonzero decay rate (the Takt gap) and the max |Im λ| among the modes AT that strict
/// gap (ω_mem). The witness builds the dense d²×d² Liouvillian directly via
/// <see cref="ChainSystem.BuildLiouvillian"/> and eigendecomposes it (rather than through
/// <see cref="Symphony"/>, whose <c>MaxN</c>=5 cap could not reach the N=6 ceiling case). For N≤5 this
/// produces the identical clock to <see cref="Symphony.Clock"/> (same L, same EVD, same tolerances).</para></summary>
public sealed class TopologyBandEdgeWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    public const double DefaultJ = 1.0;
    public const double DefaultGamma = 0.05;
    private const double Tol = 1e-6;

    private static readonly TopologyKind[] Topos = { TopologyKind.Chain, TopologyKind.Star, TopologyKind.Ring };
    private static readonly int[] Ns = { 3, 4, 5 };

    public double J { get; }
    public double Gamma { get; }
    public TopologyBandEdgeWitness(double j = DefaultJ, double gamma = DefaultGamma)
    {
        if (j <= 0) throw new ArgumentOutOfRangeException(nameof(j), j, "J must be > 0.");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be > 0.");
        J = j; Gamma = gamma;
    }

    private readonly Dictionary<(TopologyKind, int, long), (double Gap, double Omega)> _cache = new();

    /// <summary>The clock (Takt gap, ω_mem) at the given (topology, N, Q): the slowest nonzero decay rate
    /// and the max |Im λ| among the modes at that strict gap, read off the dense L spectrum. Same reading
    /// as <see cref="Symphony.Clock"/>; built directly so it reaches N=6 (above Symphony's MaxN=5). The
    /// coupling is J=Q·γ (the band edge scales with it), so the Q-sweep is a genuine sweep, not a rescale.</summary>
    private (double Gap, double Omega) ClockAt(TopologyKind topo, int n, double q)
    {
        var key = (topo, n, (long)Math.Round(q * 1e6));
        if (_cache.TryGetValue(key, out var c)) return c;

        var chain = new ChainSystem(N: n, J: q * Gamma, GammaZero: Gamma,
            HType: HamiltonianType.XY, Topology: topo);
        var L = chain.BuildLiouvillian();
        c = EigenvalueClockExtraction.ExtractClockFromSpectrum(L.Evd().EigenValues);
        _cache[key] = c;
        return c;
    }

    /// <summary>The band edge J·ρ for (topology, N) — the law (the bonds carry the topology).</summary>
    public double BandEdge(TopologyKind topo, int n)
    {
        var bonds = new ChainSystem(N: n, J: J, GammaZero: Gamma, HType: HamiltonianType.XY, Topology: topo).Bonds;
        return TopologyBandEdge.BandEdge(n, bonds, J);
    }

    /// <summary>ω_mem read at the strict gap for (topology, N) at Q (dimensionful, scales with J=q·γ).</summary>
    public double OmegaAtGap(TopologyKind topo, int n, double q) => ClockAt(topo, n, q).Omega;

    /// <summary>The dimensionless reading ω_mem / band-edge (= the verifier's w_strict/band): 1 when the
    /// band edge is the strict gap mode, 0 when a real mode holds the gap, √2 for the ring-N=4 (2,2) mode
    /// (Im=2√2·J ÷ band edge 2J). J=q·γ scales out.</summary>
    public double OmegaOverBand(TopologyKind topo, int n, double q)
    {
        double rho = BandEdge(topo, n) / J;        // dimensionless ρ (J-independent)
        double band = (q * Gamma) * rho;           // band edge at the actual coupling J=q·γ
        return band > Tol ? ClockAt(topo, n, q).Omega / band : 0.0;
    }

    /// <summary>Regime label from a two-Q-point read (operating Q and ×50).</summary>
    private string Regime(TopologyKind topo, int n, double qLo, double qHi)
    {
        double rho = BandEdge(topo, n) / J;
        double bandLo = (qLo * Gamma) * rho, bandHi = (qHi * Gamma) * rho;
        double oLo = OmegaAtGap(topo, n, qLo), oHi = OmegaAtGap(topo, n, qHi);
        if (Math.Abs(oLo - bandLo) < Tol) return "PROTECTED (band edge is the gap mode)";
        if (Math.Abs(oHi - bandHi) < Tol) return $"Q-horizon (ω_mem → band edge by Q={qHi.ToString("0", Inv)})";
        if (oHi > Tol) return $"co-occupied floor (a different mode at −2γ, Im/J = {(oHi / (qHi * Gamma)).ToString("0.###", Inv)} ≠ ρ)";
        return "horizon not reached at these Q (N≤5: rises slowly; N≥6 would be a structural ceiling — see the test suite)";
    }

    public string DisplayName =>
        $"TopologyBandEdgeWitness (band edge = J·ρ, J={J.ToString("0.###", Inv)}, γ={Gamma.ToString("0.###", Inv)})";

    public string Summary =>
        "the topology-general band edge (typed home: TopologyBandEdgeClaim): J × adjacency spectral radius " +
        "(chain 2cos(π/(N+1)), star √(N−1), ring 2). The Re=−2γ floor is the Absorption Theorem. Gap-dominance " +
        "(whether ω_mem reads the band edge) is topology-specific — see the law and the map below.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return LawNode();
            yield return MapNode();
        }
    }

    private InspectableNode LawNode()
    {
        var rows = new List<IInspectable>();
        foreach (var topo in Topos)
            foreach (var n in Ns)
            {
                double rho = BandEdge(topo, n) / J;
                double closed = topo switch
                {
                    TopologyKind.Chain => 2.0 * Math.Cos(Math.PI / (n + 1)),
                    TopologyKind.Star => Math.Sqrt(n - 1),
                    _ => 2.0,
                };
                rows.Add(new InspectableNode($"{topo} N={n}",
                    summary: $"ρ = {rho.ToString("0.#####", Inv)} (closed {closed.ToString("0.#####", Inv)}, " +
                             $"{(Math.Abs(rho - closed) < 1e-9 ? "match" : "MISMATCH")}); band edge = J·ρ = {(J * rho).ToString("0.#####", Inv)}"));
            }
        return new InspectableNode("the law (band edge = J × adjacency spectral radius)",
            summary: "the SE block of XY = J × adjacency, so the band edge = J·ρ; chain 2cos(π/(N+1)), star √(N−1), ring 2.",
            children: rows);
    }

    private InspectableNode MapNode()
    {
        const double qLo = 20.0, qHi = 1000.0;
        var rows = new List<IInspectable>();
        foreach (var topo in Topos)
            foreach (var n in Ns)
                rows.Add(new InspectableNode($"{topo} N={n}",
                    summary: $"at Q={qLo.ToString("0", Inv)}: ω_mem/J = {(OmegaAtGap(topo, n, qLo) / J).ToString("0.####", Inv)}, " +
                             $"band edge/J = {(BandEdge(topo, n) / J).ToString("0.####", Inv)} → {Regime(topo, n, qLo, qHi)}"));
        return new InspectableNode("the gap-dominance map (is the band edge the strict gap mode?)",
            summary: "is the band edge the strict gap mode? At Q=20: PROTECTED = chain (all N), star N≤4, odd rings. " +
                     "Q-horizon (broken at Q=20, climbs to the band edge by high Q) = star N=5. " +
                     "Co-occupied floor (a different mode at −2γ, never the band edge) = ring N=4 (a (2,2) mode, Im=2√2·J > band edge). " +
                     "Structural ceiling (gap saturates below 2γ) = star N≥6 (test suite only).",
            children: rows);
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
