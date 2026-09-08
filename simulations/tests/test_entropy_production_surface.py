from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "entropy_production.py"
RESULT = ROOT / "results" / "entropy_production.txt"
DOCUMENT = Path(__file__).parents[2] / "experiments" / "ENTROPY_PRODUCTION.md"


def test_source_and_artifact_keep_thermodynamics_out_of_algebraic_diagnostics():
    combined = SCRIPT.read_text(encoding="utf-8") + RESULT.read_text(encoding="utf-8")

    for stale in (
        "maximum entropy production",
        "σ peak",
        "ρ → ρ_ss",
        "Jarzynski",
        "Crooks",
        "β_eff",
        "infinite-temperature bath",
        "T_eff",
        "Carnot",
        "η =",
        "No net entropy production",
        "heat engine",
    ):
        assert stale not in combined

    for required in (
        "chosen I/d reference",
        "largest sampled dS/dt",
        "algebraic decay-rate pair sum",
        "not a thermodynamic entropy-production scale",
        "formal gain-loss generator is not a physical Lindblad channel",
    ):
        assert required in combined


def test_document_does_not_promote_first_grid_point_to_continuum_peak():
    document = DOCUMENT.read_text(encoding="utf-8")
    assert "largest sampled dS/dt" in document
    assert "first point, t=0.01" in document
    assert "not a continuum-time peak" in document
    assert "grows logarithmically as t approaches zero from above" in document
    assert "reports a peak dS/dt" not in document
    assert "reported dS/dt peak" not in document


def test_zero_entropy_control_requires_hamiltonian_invariance_not_only_z_diagonality():
    document = DOCUMENT.read_text(encoding="utf-8")
    normalized = " ".join(document.split())
    assert "both a Z-dephasing pointer state and an eigenstate of this Hamiltonian" in normalized
    assert "Z-diagonality alone is not enough" in normalized
    assert "|01⟩⟨01| initially has no Z-basis coherence" in normalized
    assert "Only initial states with coherence" not in document


def test_bell_cpsi_crossing_is_refined_and_flat_variance_has_no_peak_time():
    result = RESULT.read_text(encoding="utf-8")
    document = DOCUMENT.read_text(encoding="utf-8")
    combined = result + document
    assert "refined CΨ=1/4 crossing: t = 0.747003" in result
    assert "sampled variance is constant; no peak time is defined" in result
    assert "t ≈ 0.747" in document
    assert "both Z-dephasing and the Hamiltonian preserve this Bell+ population" in document
    assert "CΨ crosses 1/4 at t ≈ 0.51" not in combined
    assert "Max variance: 0.250000 at t = 0.01" not in combined


def test_fragile_bridge_initial_state_is_four_qubit_ghz_not_a_bell_pair():
    source = SCRIPT.read_text(encoding="utf-8")
    result = RESULT.read_text(encoding="utf-8")
    document = DOCUMENT.read_text(encoding="utf-8")
    combined = source + result + document
    assert "four-qubit GHZ/cat coherence" in source
    assert "four-qubit GHZ/cat coherence" in result
    assert "four-qubit GHZ/cat coherence" in document
    assert "Bell pair across the bridge" not in combined


def test_phase2_all_pairs_are_distinct_from_phase3_positive_rate_population():
    source = SCRIPT.read_text(encoding="utf-8")
    result = RESULT.read_text(encoding="utf-8")
    document = DOCUMENT.read_text(encoding="utf-8")
    combined = source + result + document
    assert "all palindrome pairs" in source
    assert "all palindrome pairs" in result
    assert "All palindrome pairs" in document
    assert "stationary/max-rate endpoint pairs" in source
    assert "stationary/max-rate endpoint pairs" in document
    assert "nonzero pairs" not in combined
    assert "exclude zero-rate pairs" not in source
