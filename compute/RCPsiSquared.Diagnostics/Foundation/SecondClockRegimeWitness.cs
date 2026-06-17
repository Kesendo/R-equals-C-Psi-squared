using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The second clock's regime map computed live (typed home: <c>SecondClockRegimeClaim</c>).
/// The {0,2}/half-filling coherence (the SECOND clock) is ONE mode across the arcs; its regime is set by
/// the single-particle band through two knobs: DEGENERACY m (the high-Q fate: a ceiling g2 = 4/(m+1) below
/// the −2γ floor iff m ≥ 4) and DISPERSION (the low-Q character: a sharp √-EP coherence horizon on a
/// dispersive band, only asymptotic protection on a flat one).
///
/// <para>The live gate is the N=4 frame: it builds the full 4^N Liouvillian for all five topologies at
/// Q=2 and Q=8, reads the high-Q fate from the commutant rep structure (<see cref="StructuralCeilingWitness"/>,
/// exact, no Q=1000 eigenvalue-precision risk), classifies the regime, and asserts it matches the
/// gate-verified 2D map — all three regimes (EP-horizon: chain/disordered/star; gradual: ring; ceiling:
/// complete) are present at N=4. The C# twin of simulations/second_clock_regime_axis.py (which gates
/// 15/15 over N=4,5,6; the N=5,6 EP/gradual confirmation needs the 1024²/4096² Liouvillian and stays in
/// the Python verifier, while the high-Q ceiling ladder below is exact and cheap at every N).</para></summary>
public sealed class SecondClockRegimeWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Gamma = 0.05;   // the fixed dephasing the verifier reads Q = J/γ against

    // ---- topology + couplings (a faithful port of the verifier; supports the disordered chain) ----

    private static List<(int i, int j, double w)> Couplings(string topo, int n)
    {
        var bonds = new List<(int, int)>();
        switch (topo)
        {
            case "chain":
            case "chain_disordered": for (int i = 0; i < n - 1; i++) bonds.Add((i, i + 1)); break;
            case "ring": for (int i = 0; i < n; i++) bonds.Add((i, (i + 1) % n)); break;
            case "star": for (int i = 1; i < n; i++) bonds.Add((0, i)); break;
            case "complete": for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) bonds.Add((i, j)); break;
            default: throw new ArgumentException($"unknown topology '{topo}'", nameof(topo));
        }
        var cps = new List<(int, int, double)>(bonds.Count);
        for (int k = 0; k < bonds.Count; k++)
        {
            var (i, j) = bonds[k];
            // disordered chain: fixed deterministic positive couplings (no RNG), still m=1 dispersive
            double w = topo == "chain_disordered" ? 0.6 + 0.8 * (((7 * k + 3) % 11) / 11.0) : 1.0;
            cps.Add((i, j, w));
        }
        return cps;
    }

    /// <summary>The N×N single-particle adjacency band (the hopping graph spectrum), its max degeneracy m,
    /// distinct-level count, and bandwidth. The band is the selector both knobs read.</summary>
    public static (double[] Band, int M, int Distinct, double Bandwidth) AdjacencyBand(string topo, int n)
    {
        var A = Matrix<double>.Build.Dense(n, n);
        foreach (var (i, j, w) in Couplings(topo, n)) { A[i, j] = w; A[j, i] = w; }
        var ev = A.Evd(Symmetricity.Symmetric)
                  .EigenValues.Select(e => e.Real).OrderBy(x => x).ToArray();
        var rounded = ev.Select(x => Math.Round(x, 6)).ToArray();
        int distinct = rounded.Distinct().Count();
        int m = rounded.GroupBy(x => x).Max(g => g.Count());
        return (ev, m, distinct, ev[^1] - ev[0]);
    }

    /// <summary>g2 = strict_gap / 2γ of the FULL 4^N Liouvillian at coupling J = Q·γ (uniform rescale of
    /// all bonds). Built directly (no kron): the Hamiltonian part −i(H⊗I − I⊗H) for the real-symmetric XY
    /// hopping H, plus the Z-dephasing dissipator, which is the real diagonal −2γ·hamming(a,b) (the
    /// Absorption Theorem). The smallest positive decay rate over 2γ.</summary>
    public static double G2Full(string topo, int n, double q)
    {
        int d = 1 << n;
        double scale = q * Gamma;

        // H (real symmetric XY hopping, amplitude w per bond), scaled to J = Q·γ
        var H = new double[d, d];
        foreach (var (i, j, w) in Couplings(topo, n))
            for (int s = 0; s < d; s++)
            {
                int bi = (s >> i) & 1, bj = (s >> j) & 1;
                if (bi != bj) { int s2 = s ^ (1 << i) ^ (1 << j); H[s2, s] += w * scale; }
            }

        int dim = d * d;
        var L = Matrix<Complex>.Build.Dense(dim, dim);
        var negI = new Complex(0, -1);
        // −i(H⊗I): rows a*d+b ← c*d+b with −i·H[a,c]
        for (int a = 0; a < d; a++)
            for (int c = 0; c < d; c++)
            {
                double h = H[a, c];
                if (h != 0.0) for (int b = 0; b < d; b++) L[a * d + b, c * d + b] += negI * h;
            }
        // +i(I⊗H) (the −(I⊗Hᵀ) term, H symmetric): rows a*d+b ← a*d+e with +i·H[b,e]
        for (int b = 0; b < d; b++)
            for (int e = 0; e < d; e++)
            {
                double h = H[b, e];
                if (h != 0.0) for (int a = 0; a < d; a++) L[a * d + b, a * d + e] += -negI * h;
            }
        // Z-dephasing dissipator = real diagonal −2γ·hamming(a,b)
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
            {
                int hamm = BitOperations.PopCount((uint)(a ^ b));
                L[a * d + b, a * d + b] += -2.0 * Gamma * hamm;
            }

        var rates = L.Evd().EigenValues.Select(e => -e.Real).Where(r => r > 1e-9);
        return rates.Min() / (2.0 * Gamma);
    }

    /// <summary>The high-Q fate (knob 1) from the exact commutant rep structure: g2 saturates at
    /// min(1, darkest-(1,1), darkest-(2,2)) — the always-present band-edge floor 1, undercut by a darker
    /// degenerate-manifold coherence iff one exists. Reuses <see cref="StructuralCeilingWitness.CommutantDarkest"/>;
    /// avoids the Q=1000 full-Liouvillian eigenvalue-precision risk and IS the stitch to StructuralCeilingClaim.</summary>
    public static double HighQCeiling(string topo, int n)
    {
        // the high-Q fate is graph-structural (degenerate dark manifold or none); the disordered chain is a
        // path graph (m=1, all-distinct band) that never ceilings, identical to the uniform chain here.
        string repTopo = topo == "chain_disordered" ? "chain" : topo;
        double g = 1.0;   // the band-edge |vac⟩⟨ψ_1| floor, always present (rate exactly 2γ)
        double? c11 = StructuralCeilingWitness.CommutantDarkest(repTopo, n, 1, 1);
        double? c22 = StructuralCeilingWitness.CommutantDarkest(repTopo, n, 2, 2);
        if (c11 is { } v1 && v1 < g) g = v1;
        if (c22 is { } v2 && v2 < g) g = v2;
        return g;
    }

    /// <summary>Classify the second-clock regime: CEILING (high-Q ceiling below 1), else EP-HORIZON (g2
    /// reaches 1 sharply by Q=8 — a dispersive band), else GRADUAL (only asymptotic protection). The exact
    /// port of the verifier's classify(): high-Q from the commutant, g2(Q=8) from the full Liouvillian.</summary>
    public static string Classify(string topo, int n)
    {
        if (HighQCeiling(topo, n) < 1.0 - 1e-4) return "CEILING";
        if (G2Full(topo, n, 8.0) > 1.0 - 1e-3) return "EP";
        return "GRADUAL";
    }

    /// <summary>The gate-verified 2D regime map at N=4 (all three regimes present): the live full-Liouvillian
    /// classification must reproduce this.</summary>
    public static readonly IReadOnlyDictionary<string, string> ExpectedN4 = new Dictionary<string, string>
    {
        ["chain"] = "EP",
        ["chain_disordered"] = "EP",
        ["ring"] = "GRADUAL",
        ["star"] = "EP",
        ["complete"] = "CEILING",
    };

    private static readonly string[] N4Topologies =
        { "chain", "chain_disordered", "ring", "star", "complete" };

    // ---- IInspectable ----

    public string DisplayName => "SecondClockRegimeWitness (the {0,2}-coherence regime = map(degeneracy, dispersion))";

    public string Summary =>
        "the second clock's regime map computed live (typed home: SecondClockRegimeClaim): the {0,2}/half-filling " +
        "coherence is ONE mode whose regime is set by the band — knob 1 (degeneracy m) → high-Q ceiling 4/(m+1) " +
        "(below the floor iff m≥4), knob 2 (dispersion) → low-Q √-EP horizon (dispersive) vs gradual (flat). " +
        "The N=4 frame is a self-validating full-Liouvillian gate (all three regimes present). " +
        "Sector overview: inspect --root blockspectrum (this zooms the (2,2)/half-filling {0,2} sector).";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheN4Frame();
            yield return Knob1Ladder();
            yield return Knob2Dispersion();
            yield return TheStitch();
        }
    }

    /// <summary>The live gate: classify all five topologies at N=4 from the full Liouvillian and assert the
    /// 2D map (5/5). EP (chain/disordered/star), GRADUAL (ring), CEILING (complete) all appear at N=4.</summary>
    private InspectableNode TheN4Frame()
    {
        var rows = new List<IInspectable>();
        int pass = 0;
        foreach (var topo in N4Topologies)
        {
            double g2 = G2Full(topo, 4, 2.0);
            double g8 = G2Full(topo, 4, 8.0);
            double ghi = HighQCeiling(topo, 4);
            string got = Classify(topo, 4);
            string want = ExpectedN4[topo];
            bool ok = got == want;
            if (ok) pass++;
            var (_, m, distinct, width) = AdjacencyBand(topo, 4);
            rows.Add(new InspectableNode($"{topo} N=4",
                summary: $"m={m}, #distinct={distinct}, bandwidth={width.ToString("0.###", Inv)}; " +
                         $"g2(Q=2)={g2.ToString("0.0000", Inv)}, g2(Q=8)={g8.ToString("0.0000", Inv)}, " +
                         $"high-Q={ghi.ToString("0.0000", Inv)} → {got} " +
                         $"({(ok ? "matches map ✓" : $"MISMATCH, expected {want}")})"));
        }
        return new InspectableNode($"the N=4 frame (live gate, {pass}/{N4Topologies.Length})",
            summary: $"the full-Liouvillian regime of every topology at N=4, classified live and checked against the " +
                     $"gate-verified 2D map: {(pass == N4Topologies.Length ? "all match ✓" : "a row MISMATCHES — investigate")}. " +
                     "All three regimes are present at N=4: EP-horizon (chain/disordered/star), GRADUAL (ring), " +
                     "CEILING (complete). A firing gate is a finding (diagnose, don't loosen).",
            children: rows);
    }

    /// <summary>Knob 1: the high-Q ceiling 4/(m+1) across N=4,5,6 (the exact commutant), showing the star
    /// graduate from reaches-floor (N≤5) to a structural ceiling (N=6, 4/5 &lt; 1) as m grows past 3.</summary>
    private InspectableNode Knob1Ladder()
    {
        var rows = new List<IInspectable>();
        foreach (var topo in new[] { "complete", "star" })
        {
            var sub = new List<IInspectable>();
            for (int n = 4; n <= 6; n++)
            {
                double ghi = HighQCeiling(topo, n);
                var (_, m, _, _) = AdjacencyBand(topo, n);
                double mPlus1 = SecondClockRegimeClaim_SymmetricForm(m);
                string fate = ghi < 1.0 - 1e-9 ? "CEILING" : "reaches the floor";
                sub.Add(new InspectableNode($"{topo} N={n}",
                    summary: $"band degeneracy m={m}, high-Q g2 = {ghi.ToString("0.#####", Inv)} " +
                             $"(4/(m+1) = {mPlus1.ToString("0.#####", Inv)}); {fate}"));
            }
            rows.Add(new InspectableNode($"{topo}: the 4/(m+1) ladder",
                summary: topo == "complete"
                    ? "K_N: m=N−1, g2=4/N, a ceiling at every N≥5 (m≥4); deepens 4/5 → 2/3 → 4/7."
                    : "star: the leaf manifold m grows with N; g2=4/(N−1) reaches the floor at N≤5 (m≤3) and " +
                      "first ceilings at N=6 (m=4, 4/5 < 1) — the star graduating into the CEILING regime.",
                children: sub));
        }
        return new InspectableNode("knob 1: the high-Q ceiling 4/(m+1) (the CEILING fate)",
            summary: "the exact commutant high-Q g2 (the stitch to StructuralCeilingClaim): a degeneracy-m manifold " +
                     "ceilings at 4/(m+1), below the −2γ floor iff m ≥ 4. The star crosses into a ceiling at N=6.",
            children: rows);
    }

    private static double SecondClockRegimeClaim_SymmetricForm(int m) => 4.0 / (m + 1);

    /// <summary>Knob 2: the dispersion the regime reads once the second clock reaches the floor. A dispersive
    /// band (distinct cosine spectrum) gives a sharp √-EP coherence horizon; a flat/marginal band gives only
    /// asymptotic protection (GRADUAL).</summary>
    private InspectableNode Knob2Dispersion()
    {
        var rows = new List<IInspectable>();
        foreach (var topo in N4Topologies)
        {
            var (band, m, distinct, width) = AdjacencyBand(topo, 4);
            string bandStr = "[" + string.Join(", ", band.Select(x => x.ToString("0.##", Inv))) + "]";
            rows.Add(new InspectableNode($"{topo} N=4 band",
                summary: $"{bandStr}: m={m}, #distinct={distinct}, bandwidth={width.ToString("0.###", Inv)} " +
                         $"({(distinct >= 4 ? "dispersive → sharp √-EP horizon" : "low-dispersion → check the floor character")})"));
        }
        return new InspectableNode("knob 2: dispersion (EP vs GRADUAL)",
            summary: "once the second clock reaches the floor, the band's dispersion sets HOW: a dispersive band " +
                     "(real cosine spectrum: chain, disordered chain, ring, star N=4) → the two clocks coalesce at a " +
                     "sharp √-EP coherence horizon Q*(N) (CoherenceHorizonClaim); a flat/marginal band (star N=5, m=3) " +
                     "→ only asymptotic 1−c/Q² protection, no sharp horizon (the star-no-horizon reading).",
            children: rows);
    }

    private static InspectableNode TheStitch() =>
        new InspectableNode("the stitch (why this is one node)",
            summary: "CoherenceHorizonClaim (the EP regime), StructuralCeilingClaim (the CEILING regime, 4/(m+1) the " +
                     "bridge), and the star-no-horizon (the GRADUAL regime) are not three findings: they are ONE mode " +
                     "(the {0,2}/half-filling coherence) whose regime = map(band degeneracy, dispersion). The N=4 " +
                     "anomalies (ring-4 GRADUAL, complete-4 CEILING) are the (2,2) half-filling sector, the recurring " +
                     "n3/n4 specials. Everything is one object, seen at max zoom on a different facet.");

    public InspectablePayload Payload => InspectablePayload.Empty;
}
