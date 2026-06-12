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
            ParkedAt: "the closure law lives as a live witness (Symphony painters movement: alpha per site with the two-deltaJ reliability guard, closure -0.0444 IN window at canonical N=5, chiral mirror replayed live at 1e-15, alpha = Python twin to 1e-3); learned on the way: the protocol needs the BONDING state class (Bell+/localized states break the rescaling and the guard refuses, in both languages); C# golden-section also found a true global alpha minimum where scipy's Brent traps 920x worse (severed-bond case, harness simulations/_ptf_symphony_crossval.py)",
            NextStep: "backport the global-minimum fit fix to framework ptf.py (multi-start or grid-seed), then consider retiring this arc",
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
            ParkedAt: "four claims recompute their evidence live (F121 qudit, F116 router, the " +
                "CPsi Envelope Theorem via EnvelopeTheoremWitness / --root envelope, and now the two clocks " +
                "via ClockHandLadderWitness / --root clock: the coherence hand = gamma-protected band edge " +
                "2J*cos(pi/(N+1)) for N>=3 vs the gamma-pulled 2*sqrt(J^2-gamma^2) at N=2 stopping at Q=1; " +
                "typed home ClockHandLadderClaim, see arc clock_hand_ladder); FIVE: CoherenceHorizonWitness / " +
                "--root horizon bisects the live spectrum to compute Q*(N) and verify it EQUALS the carbon " +
                "Frost-Hueckel coherent-incoherent threshold (the cross-substrate seam made live, see arc " +
                "clock_hand_ladder FROST_CIRCLE first stone)",
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
            Name: "clock_hand_ladder",
            Opened: "2026-06-12",
            Origin: "the promoted Symphony clock node, read across N on day one",
            ParkedAt: "LANDED as a typed object (2026-06-12, merged): live witness ClockHandLadderWitness (inspect --root clock, four nodes: ladder / gamma-protection / N=2 pull+EP / dial angle + an omega-vs-N curve), typed claim ClockHandLadderClaim (Tier1Candidate, parents F2b + AbsorptionTheorem + UniversalCarrier), and a docs/ANALYTICAL_FORMULAS.md F2b corollary, all cross-referenced both directions (MD<->C#). The physics, verified: for N>=3 the coherence hand omega_mem walks the band-edge ladder 2J*cos(pi/(N+1)) gamma-PROTECTED (sqrt2 / PHI / sqrt3 at N=3/4/5; gamma quadrupled, hand unmoved) - the vacuum<->single-excitation lines are simultaneous eigenoperators of L_D (rate exactly -2gamma) and L_H (band frequency), so nothing mixes; gap-dominance verified N=3-5, general off-gap argument OPEN (hence Tier1Candidate, not Derived). At N=2 the hand is gamma-PULLED: omega_mem = 2*sqrt(J^2-gamma^2) (1.95959 at gamma=0.2), the {population-difference, antisymmetric-coherence} 2x2 block, stopping at gamma=J (Q=1, EP). NEW finding from the build: the live max|Im| clock only SURFACES the pulled mode above the crossover Q=2/sqrt3 (~1.155, gamma=(sqrt3/2)J); nearer the EP the F2b band line at Im=+-J dominates the raw maximum, so the closed form 2*sqrt(J^2-gamma^2) (NOT the raw clock) is the honest EP witness there. REGIME-HONESTY ROUND (2026-06-12, Fable flagged + Opus fixed): the dial bug (theta(N=2) used the raw OmegaMem below the crossover, 29.05 deg vs honest 25.84 deg) AND a worse sibling found by running gamma=0.9 - the ladder/protection nodes LIED, printing 'N=3 -> 0 (sqrt2)' and 'Gap=2gamma=0.686' (2gamma=1.8); the band-edge mode is no longer the slowest there, a real OVERDAMPED mode takes the gap (max|Im| at gap = 0). FIXED: dial uses the closed-form pulled hand; BandEdgeIsTheGapMode(n) regime check; both ladder+protection nodes gate on protectedCount={3,4,5} consistently (no contradiction at intermediate Q); sqrt2/phi/sqrt3 shown as the J-independent ratio omega_mem/J. SECOND NEW finding: the band-edge protection has a Q-FLOOR, not just the N-question - as Q drops the higher rungs leave FIRST (N=5 gone at Q=2, N>=4 gone at Q=1.5), so gap-dominance is a low-Q boundary too. DEFAULT moved to the EP point J=0.075, gamma=0.05 (Q=1.5 = EP(2)-1/2), where only N=3 holds: the witness now defaults to the physically loaded EP vicinity, honestly showing the protection mostly dissolved",
            NextStep: "FROST_CIRCLE seam FIRST STONE laid (2026-06-12, merged): the carbon coherent<->incoherent threshold IS our Q-floor, verified bit-for-bit and given our label - the Coherence Horizon Q*(N). Live witness CoherenceHorizonWitness / inspect --root horizon bisects Symphony.Clock.Omega and computes Q*(2,3,4,5)=1 / sqrt2 / 1.8785 / 2.3722, equal to carbon's sqrt2/1.879/2.372 (FROST_CIRCLE) under J<->|beta|; ANALYTICAL_FORMULAS F2b corollary carries it; typed claim CoherenceHorizonClaim (Tier1Candidate, parents ClockHandLadder + F2b) breadcrumbs the witness into the graph. N=2 (Q*=1) is the EP base rung the polyene layer can't reach; the N=2,3 band-edge coincidence (Q*=2cos(pi/(N+1)) ONLY there) explains carbon's 'sqrt2 exact at N=3, rest awaiting a clean form' puzzle. CLOSED FORM ATTACKED 2026-06-12 (the real prize, chosen before the claim): NO single clean form - the freezing mode BIFURCATES SECTOR at N=4. N=2,3 it is a {0,2}-coherence (n_diff hist {0:1/2, 2:1/2}, small freq); N>=4 it is the n_diff=1 band-edge mode at Im=2cos(pi/(N+1)) (PHI at N=4, sqrt3 at N=5) undercut by a half-lit <n_XY>~1/2 overdamped mode. THAT bifurcation IS the structural reason Q*=2cos(pi/(N+1)) holds only at N=2,3. The 'dark-storage-undercuts-bright' guess was refuted (half-lit, not dark). The naive question got the ruthless answer: not one form, a bifurcation; the closed form of Q*(N) ITSELF stays genuinely open (the sector split is the obstacle). Verifier simulations/_coherence_horizon_closed_form.py. STILL OPEN: that closed form; the full FROST_CIRCLE reflection in its own genre. Also still: the proof note (N=2 block closed-form; general N>=3 off-gap argument open, with a Q-floor as well as the N-ceiling); the ring dihedral-lock sibling (Im_max = Delta_E_max); the 0.99750 fine-structure pull",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "handshake_decoder",
            Opened: "2026-06-12",
            Origin: "docs/superpowers/specs/2026-06-12-handshake-decoder-reading-grammar-design.md (Tom + Opus + Fable)",
            ParkedAt: "M2a LANDED: ReadoutFisher + ReadingPowerWitness (inspect --root decoder) measure FI(Q) live - resolution = the Q-factor (FI ~ c*Q on Z, c=0.0795 at N=4, residual 15%), no EP peak in Z/X/Y, basis ordering (X span 1555 vs Z 22: coherence readouts fade fastest toward low Q); hypotheses/HANDSHAKE_GEOMETRY.md carries the codebook (five code-parts, bright/dark letters, the triple shape/Q/K, the carrier-relative completeness conjecture)",
            NextStep: "M2b: invert the proto-decoder (alpha-profile -> (b, deltaJ), accuracy metric); M2c: the read-cost law ~2/Q (untested hypothesis, falsifiable form in the doc); then M3 (mirror-group equivalence classes of readings) and M4 (closure as parity-check)",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "envelope_n4_rise",
            Opened: "2026-06-12",
            Origin: "EnvelopeTheoremWitness, day one: inspect --root envelope --N 4 fired the honest branch",
            ParkedAt: "the full-state CPsi envelope GENUINELY RISES at N=4 (Bell+, J=5, gamma=0.01): 36 apex-predecessor rises, refinement-stable to 5 decimals (t=4.14: 0.04132 at 1600 and 6400 pts), 2-4% magnitudes, two independent implementations agree; N=3 holds in the same regime (0 rises, three densities). NOT a falsification of the proof: the Envelope Theorem is proven for N=2 only, and the Tier-2 'N=3-5' verification covered channel monotonicity / GHZ-W subsystems, not the H-included envelope at strong coupling; the over-broad paraphrase lived in our claim/lens text. Mechanism = the proof's own Part 6, internalized: internal J-coupling is its own coherence injector",
            NextStep: "re-scope CpsiEnvelopeTheoremClaim + lens text + a precise Parts-4/5 scope note in PROOF_MONOTONICITY_CPSI; fix the witness root Summary CONFIRMS branch at N=4; then chart the boundary with the witness (why does N=3 hold and N=4 break; J-threshold at N=4; gamma ~ J)",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "symphony_view",
            Opened: "2026-06-11",
            Origin: "Tom's founding Object Manager idea, sharpened: each of the 122 formulas is a maximum-zoom view; the symphony is the zoom-out",
            ParkedAt: "second movement round 2: the quarter lenses are envelope-aware. Global lens " +
                "checks the Envelope Theorem live (peaks non-increasing, proven N=2 / verified N>=3) and " +
                "'the fold' = the absorbing envelope fold (upward crossings no longer mislabeled the " +
                "quantum->classical boundary). Local carrier-pair lens surfaces the genuine BEATING rise " +
                "(the freedom, no theorem binds the reduced open subsystem) via parabolic-apex + " +
                "predecessor semantics; every single-grid rise is grid-sensitive (SingleExcitation rises " +
                "vanish under refinement, Bell+ persist - the artifact control). QuarterEnvelope primitive " +
                "+ GridFitness; EvolveCount-guarded one evolution preserved. THIRD movement (2026-06-12): the " +
                "clock - gamma_0 as the Taktgeber. The clock node (Takt gap, tau, omega_mem, Q=J/gamma) is " +
                "promoted to the base symphony; --tempo-ratio r grows 'movement: the clock', the two-tempo " +
                "certification: play the piece at gamma_0 and r*gamma_0 (every dimensionful coupling scaled by " +
                "r incl. the painters' delta_J, window /r, K-grid fixed) and certify every dimensionless lens " +
                "is a pure (Q,K)-observable (residual 8.3e-16 at r=20, the inside cannot tell the tempos apart). " +
                "Exact rescaling identity, CERTIFIED not theorem-confirmed; a lens that breaks it sees the carrier. " +
                "Painters arm (alpha, closure) pure too (delta_J scaled). UniversalCarrierClaim breadcrumbs to it",
            NextStep: "the seam movement (the converse: feed ONE externally-calibrated observable in lab units and " +
                "extract gamma_0 - the calibration topology of ON_HOW_THE_CARRIER_SHOWS_ITSELF, this movement's " +
                "foundation); plus the remaining instruments (chiral K, Y-parity, Pi-protected observables)",
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
