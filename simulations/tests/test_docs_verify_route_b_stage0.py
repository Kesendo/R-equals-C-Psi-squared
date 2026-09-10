import copy
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).parents[2]
SIMULATIONS = ROOT / "simulations"
if str(SIMULATIONS) not in sys.path:
    sys.path.insert(0, str(SIMULATIONS))

import docs_verify


N6 = json.loads((SIMULATIONS / "results" / "route_b_a2_n6.json").read_text(encoding="utf-8"))
READOUT_DOCS = (
    ROOT / "experiments" / "ROUTE_B_N4_VIRTUAL_READOUT.md",
    ROOT / "experiments" / "ROUTE_B_N4_READOUT_SEARCH.md",
    ROOT / "experiments" / "ROUTE_B_N4_HISTOGRAM_FILTER.md",
)
EXPECTED_TEST = "ReadingPowerWitnessTests.FiIsMonotoneOnTheSevenPointSampledQGrid_InEveryBasis"


def test_n6_exact_box_contract_accepts_the_committed_schema3_inventory():
    assert docs_verify.route_b_n6_real_q_exclusion_errors(N6) == []


@pytest.mark.parametrize("field,value", [("schemaVersion", 2), ("n", 5)])
def test_n6_exact_box_contract_rejects_wrong_document_identity(field, value):
    changed = copy.deepcopy(N6)
    changed[field] = value
    assert docs_verify.route_b_n6_real_q_exclusion_errors(changed)


def test_n6_exact_box_contract_rejects_wrong_parameter_convention():
    changed = copy.deepcopy(N6)
    changed["conventions"]["parameter"] = "t=qCSharp"
    assert docs_verify.route_b_n6_real_q_exclusion_errors(changed)


def test_n6_exact_box_contract_rejects_boxes_detached_from_the_certified_carrier():
    changed = copy.deepcopy(N6)
    for locus in changed["loci"]:
        locus["tBox"]["real"] = {
            "lower": {"numerator": "1", "denominator": "1"},
            "upper": {"numerator": "2", "denominator": "1"},
        }
    assert docs_verify.route_b_n6_real_q_exclusion_errors(changed)


def test_n6_exact_box_contract_rejects_missing_or_duplicate_loci():
    missing = copy.deepcopy(N6)
    missing["loci"].pop()
    assert docs_verify.route_b_n6_real_q_exclusion_errors(missing)

    duplicate = copy.deepcopy(N6)
    duplicate["loci"][1]["id"] = duplicate["loci"][0]["id"]
    assert docs_verify.route_b_n6_real_q_exclusion_errors(duplicate)


@pytest.mark.parametrize("bad", [True, 1.0, "1.5"])
def test_n6_exact_box_contract_rejects_non_integer_rational_endpoints(bad):
    changed = copy.deepcopy(N6)
    changed["loci"][0]["tBox"]["real"]["lower"]["numerator"] = bad
    errors = docs_verify.route_b_n6_real_q_exclusion_errors(changed)
    assert any("malformed exact tBox.real interval" in error for error in errors)


def test_n6_exact_box_contract_rejects_float_inverted_interval():
    changed = copy.deepcopy(N6)
    real = changed["loci"][0]["tBox"]["real"]
    real["lower"]["numerator"] = 1.25
    real["upper"]["numerator"] = -1.25
    errors = docs_verify.route_b_n6_real_q_exclusion_errors(changed)
    assert any("malformed exact tBox.real interval" in error for error in errors)


def test_n6_exact_box_contract_rejects_wrong_parity_counts_and_multiplicity():
    changed = copy.deepcopy(N6)
    changed["loci"][0]["parity"] = "O"
    changed["loci"][1]["algebraicMultiplicity"] = 1
    assert docs_verify.route_b_n6_real_q_exclusion_errors(changed)


def test_readout_docs_and_verifier_require_the_actual_test_symbol():
    assert docs_verify.READING_POWER_MONOTONICITY_TEST == EXPECTED_TEST
    for path in READOUT_DOCS:
        assert EXPECTED_TEST in path.read_text(encoding="utf-8")


def test_f164_proof_scopes_the_only_claim_to_the_two_a2_inventories():
    proof = " ".join((ROOT / "docs" / "proofs" / "PROOF_N5_REAL_Q_DIABOLIC.md")
                     .read_text(encoding="utf-8").split())
    assert "only real-q points in the N=5/N=6 A2 inventories" in proof
    assert "q is the settable parameter" in proof
    assert "lambda is the resulting spectral eigenvalue" in proof
    assert "observed eigenvalue" not in proof
    assert "pinned by measurement" not in proof


def test_reading_power_current_truth_keeps_nonzero_coherence_readout_and_live_span():
    text = " ".join((ROOT / "hypotheses" / "HANDSHAKE_GEOMETRY.md")
                    .read_text(encoding="utf-8").split())
    assert "1555.3×" in text
    assert "X/Y remain nonzero" in text
    assert "population basis remains much stronger" in text
    assert "1670×" not in text
    assert "only the population basis still reads" not in text
