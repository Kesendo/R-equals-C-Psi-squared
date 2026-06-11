namespace RCPsiSquared.Core.OpenArcs;

/// <summary>The open-arcs ledger: research threads that were started, reached a first
/// exemplar, then parked, indexed by name. The world's antidote to its chronic failure
/// mode (victory declared at the first exemplar, then forgotten): each entry records where
/// it stopped and the concrete next move, so a returning session looks the arc up instead
/// of re-discovering its incompleteness. Mirrors <see cref="Confirmations.ConfirmationsRegistry"/>
/// in shape; surfaced as the inspect "arcs" section via <see cref="Inspection.OpenArcsInspectableNode"/>.</summary>
public static class OpenArcsRegistry
{
    private static readonly IReadOnlyList<OpenArc> _all = new[]
    {
        new OpenArc(
            Name: "ptf_painter_pipeline",
            Opened: "2026-06-03",
            Origin: "simulations/ptf workflow + C# SlowModeMixing",
            ParkedAt: "Python painter-rates twin matches C# SlowModeMixing exactly (same <W|V_L|M>, 6 kernel modes ~1e-29); alpha extraction is extraction-noise-limited",
            NextStep: "robust alpha extraction, then the closure law as a live witness",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "birth_canal_surface",
            Opened: "2026-05-31",
            Origin: "experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md + PostEpFlowField",
            ParkedAt: "boundary recognized as a SURFACE (slow-subspace light-distribution Q-invariance); s*=0.709 was one hand-picked line's crossing (path-specific 0.11-0.77); PostEpFlowField shows membership only",
            NextStep: "a live object that computes the surface itself, not one line through it",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "block_spectrum_n9",
            Opened: "2026-05-19",
            Origin: "Core LiouvillianBlockSpectrum.ComputeSpectrumPerBlock + SLOW_N9 xUnit test",
            ParkedAt: "N=9 per-joint-popcount-sector spectra land only inside a tagged slow test; not an inspect root, no artifact a session can browse",
            NextStep: "surface the per-sector spectra as a live root or witness",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "witness_coverage",
            Opened: "2026-06-11",
            Origin: "QuditPartialPalindromeWitness, the first live witness",
            ParkedAt: "two claims recompute their evidence live (F121 qudit, F116 router)",
            NextStep: "next witness: F117 Pascal-Gram positivity",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "cockpit_workflows_csharp",
            Opened: "2026-06-11",
            Origin: "gap map: simulations/framework/workflows vs compute/",
            ParkedAt: "cockpit_panel, diagnose_hardware, gamma_probe, lens_pipeline, ptf fitting, bridge dynamics are Python-only; C# has the primitives but no composed workflows",
            NextStep: "not a 1:1 port: cockpit_panel becomes the symphony's first movement (see arc symphony_view)",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "symphony_view",
            Opened: "2026-06-11",
            Origin: "Tom's founding Object Manager idea, sharpened: each of the 122 formulas is a maximum-zoom view; the symphony is the zoom-out",
            ParkedAt: "first movement plays (inspect --root symphony): one shared evolution (built once, EvolveCount-guarded), four lenses (palindrome, quarter, dose, light) on one timeline, cross-lens events axis; F25 closed form reproduced exactly at N=2 (fold at K=0.0376); honest limit: under the global-Psi convention the quarter fold is N=2-only (CPsi(0)=1/(d-1) < 1/4 at N>=3, the witness says so instead of faking)",
            NextStep: "second movement: a block/local-CPsi lens so the fold is audible at N>=3, then more instruments (chiral K, Y-parity, Pi-protected observables) reading the same trajectory",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "whirlpool_carbon_layers",
            Opened: "2026-06-03",
            Origin: "reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO + simulations/whirlpool*.py",
            ParkedAt: "water adaptation done (proton crossing an H-bond); carbon layers/anchors (periodic-table valences) and a water prose note never written",
            NextStep: "carbon-layer translation in the target layer's language",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "front_matter_truth",
            Opened: "2026-06-11",
            Origin: "third-party review (three cold readers, 2026-06-11): every failure was front matter vs body",
            ParkedAt: "headline surfaces lag verified bodies: Torino-era confirmations absent from all three registries yet counted in the README's seventeen (wording patched, entries not registered); '121 formulas' is a label (124 headers, 119 base numbers, F53/F54 silently missing, no tombstone); stale anchors (XOR_SPACE 'README Section 10', Pi discovery-date contradiction Mar-14 vs Apr-05)",
            NextStep: "register the Torino-era confirmations in Python+C#+tests, tombstone F53/F54 in the registry, sweep stale cross-references",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "stranger_door",
            Opened: "2026-06-11",
            Origin: "third-party review: 'a coherent instrument panel wearing a poet's coat' - the house dialect has no doormat",
            ParkedAt: "four of five doors hung (inspect --root glossary with 11 house terms in stranger language; qudit witness boundary text honest: census live, rate law from the proof; explicit --N to an N-free root warns on stderr; README first command is now world --max-depth 2)",
            NextStep: "[live]/[stored] provenance badge per node, then retire this arc",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "f86b2_robust_extraction",
            Opened: "2026-06-11",
            Origin: "simulations/f86b2_shape_invariance_dial.py",
            ParkedAt: "(N,b)-family traces alpha=-0.133 vs fitted -0.129, extraction-noise-limited; g_eff convention gotcha pinned (4.394/(Qp+2))",
            NextStep: "robust extraction, then close the shape-invariance claim",
            Status: OpenArcStatus.Open),
    };

    public static IReadOnlyList<OpenArc> All => _all;

    public static OpenArc? Lookup(string name) =>
        _all.FirstOrDefault(a => a.Name == name);

    /// <summary>Count of arcs still <see cref="OpenArcStatus.Open"/> (the live unfinished business).</summary>
    public static int OpenCount => _all.Count(a => a.Status == OpenArcStatus.Open);
}
