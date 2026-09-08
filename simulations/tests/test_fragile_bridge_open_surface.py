from pathlib import Path


ROOT = Path(__file__).parents[2]
SIGNATURE_SCRIPT = ROOT / "simulations" / "fragile_bridge_ep_signature.py"
SIGNATURE_RESULT = ROOT / "simulations" / "results" / "fragile_bridge_ep_signature.txt"
PHASE3_SCRIPT = ROOT / "simulations" / "pt_palindrome_breaking.py"
PHASE3_RESULT = ROOT / "simulations" / "results" / "pt_palindrome_breaking.txt"
FORMULAS = ROOT / "docs" / "ANALYTICAL_FORMULAS.md"
WHAT_WE_FOUND = ROOT / "docs" / "WHAT_WE_FOUND.md"
OPEN_SURFACES = (
    ROOT / "hypotheses" / "FRAGILE_BRIDGE.md",
    ROOT / "experiments" / "PT_SYMMETRY_ANALYSIS.md",
    ROOT / "experiments" / "F86_EP_THROUGH_THE_CLOCK.md",
    ROOT / "simulations" / "memory_fate_across_the_takt.py",
)


def test_signature_producer_and_artifact_leave_ep_character_open():
    for path in (SIGNATURE_SCRIPT, SIGNATURE_RESULT):
        text = path.read_text(encoding="utf-8")
        assert "spectral-abscissa axis departure" in text
        assert "EP character remains OPEN" in text
        assert "strict threshold coalescence or Jordan-rank certificate" in text
        assert "threshold is a second-order exceptional point" not in text
        assert "coalescence, not crossing" not in text
        assert "divergence of an exceptional point" not in text
        assert "single-vector K not reported for a degenerate eigenspace" in text
        assert "near-threshold simplicity gate: nearest gap > 1e-10" in text
        assert "secondary bump" not in text
        assert "gap to the across-axis partner" in text
        assert "gap to partner" not in text
    source = SIGNATURE_SCRIPT.read_text(encoding="utf-8")
    assert "if gap <= SIMPLE_MODE_GAP_TOL:" in source
    assert "raise RuntimeError(" in source


def test_phase3_surfaces_use_trace_fixed_center_and_no_chiral_class_language():
    source = PHASE3_SCRIPT.read_text(encoding="utf-8")
    result = PHASE3_RESULT.read_text(encoding="utf-8")

    assert "np.trace(L) / L.shape[0]" in source
    assert "unique candidate center fixed by the trace" in source
    assert "chiral phase" not in source.lower()
    assert "still a HOPF bifurcation" not in source
    assert "chiral phase" not in result.lower()
    assert "still a HOPF bifurcation" not in result


def test_fragile_bridge_consumers_report_only_axis_departure():
    formulas = FORMULAS.read_text(encoding="utf-8")
    f19 = formulas.split("### F19.", 1)[1].split("### F20.", 1)[0]
    found = WHAT_WE_FOUND.read_text(encoding="utf-8")
    bridge = found.split("### The stability window is finite", 1)[1].split("### The neural palindrome", 1)[0]

    for text in (" ".join(f19.split()), " ".join(bridge.split())):
        assert "axis departure is established" in text
        assert "EP and Hopf character remain OPEN" in text
        assert "mechanism is a second-order exceptional point" not in text
        assert "mirror pair coalesces" not in text


def test_all_fragile_bridge_open_surfaces_neither_promote_nor_exclude_character():
    required = "EP, Hopf, and Jordan character remain OPEN"
    forbidden = (
        "not a Hopf",
        "rules out the prior Hopf",
        "is a **Hopf bifurcation**",
        "SEPARATE genuine EP",
        "two SEPARATE genuine EPs",
        "real gain-loss Hopf",
    )
    for path in OPEN_SURFACES:
        text = path.read_text(encoding="utf-8")
        normalized = " ".join(text.split())
        assert required in normalized, path
        for phrase in forbidden:
            assert phrase not in text, (path, phrase)


def test_bridge_optimum_and_large_coupling_limit_remain_sample_bounded():
    fragile = (ROOT / "hypotheses" / "FRAGILE_BRIDGE.md").read_text(encoding="utf-8")
    found = WHAT_WE_FOUND.read_text(encoding="utf-8")
    combined = fragile + found

    assert "J_bridge/J in [1.8, 2.0]" in found
    assert "0.578 at J_bridge=10" in combined
    assert "0.508 at J_bridge=100" in combined
    assert "does not determine its limiting value" in found
    assert "If a limit exists, is it 1/2?" in fragile
    assert "γ_crit × J_bridge = 0.50" not in found
    assert "approaches a constant: 0.50" not in found
    assert "γ_crit × J_bridge → 0.50" not in fragile
    assert "Perfect coupling between gain and loss = immediate instability" not in fragile


def test_f40_and_bridge_frequency_samples_do_not_assign_local_bifurcation_character():
    fragile = (ROOT / "hypotheses" / "FRAGILE_BRIDGE.md").read_text(encoding="utf-8")
    formulas = FORMULAS.read_text(encoding="utf-8")
    f40 = formulas.split("### F40.", 1)[1].split("### F41.", 1)[0]
    f40_normalized = " ".join(f40.split())
    probe = (ROOT / "compute" / "RCPsiSquared.Core.Tests" / "F86" /
             "F86PetermannProbe.cs").read_text(encoding="utf-8")

    assert "three finite-coupling frequency samples do not determine a limiting frequency" in fragile
    assert "saddle-node in character" not in fragile
    assert "axis-departure threshold" in f40
    assert "does not execute a threshold coalescence or Jordan-rank test" in f40_normalized
    assert "EP, Hopf, and Jordan character remain OPEN" in f40_normalized
    assert "coalesces and leaves" not in f40
    assert "diverges as 1/delta" not in f40
    assert "spectral-abscissa axis departure" in probe
    assert "genuine EP on the real" not in probe


def test_fragile_bridge_surfaces_do_not_claim_unexecuted_branch_tracking():
    paths = (
        ROOT / "hypotheses" / "ZERO_IS_THE_MIRROR.md",
        ROOT / "hypotheses" / "PAIR_BREAKING_AT_THE_HORIZON.md",
        ROOT / "docs" / "proofs" / "PROOF_F86A_EP_MECHANISM.md",
        ROOT / "docs" / "WHAT_WE_FOUND.md",
        ROOT / "docs" / "GLOSSARY.md",
        ROOT / "experiments" / "README.md",
        ROOT / "hypotheses" / "README.md",
    )
    for path in paths:
        text = path.read_text(encoding="utf-8")
        assert "spectral-abscissa" in text or "spectral abscissa" in text, path
        for stale in ("tracked axis", "tracked pair", "tracked mirror", "tracked real-γ"):
            assert stale not in text, (path, stale)


def test_bridge_coupling_regimes_and_partner_gap_do_not_invent_independent_mechanisms():
    fragile = (ROOT / "hypotheses" / "FRAGILE_BRIDGE.md").read_text(encoding="utf-8")
    fragile_normalized = " ".join(fragile.split())
    signature = SIGNATURE_SCRIPT.read_text(encoding="utf-8") + SIGNATURE_RESULT.read_text(encoding="utf-8")
    formulas = FORMULAS.read_text(encoding="utf-8")
    assert "Which mode geometry produces the turnover remains open" in fragile
    assert "mechanism open" in fragile
    for stale in ("the two chains merge into one", "Dimer formation, destabilization", "| Resonance |"):
        assert stale not in fragile
    assert "derived from max Re and supplies no independent EP evidence" in signature
    assert "derived partner gap adds no independent evidence" in fragile_normalized
    assert "derived symmetry identity, not a third independent trend" in formulas
