"""Task 6 source-local contracts for neutral calibration/hardware producers."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[2]


CONTRACTS = {
    "data/ibm_history/ibm_history_analysis.py": {
        "required": (
            "R_STAR = 0.21275477982200533",
            "normalized-purity proxy",
            "proxy_band",
            "normalized_purity_proxy_analysis.png",
            "normalized_purity_proxy_at_time",
            "find_proxy_quarter_crossing",
            "proxy = 1/4",
        ),
        "forbidden": (
            "Purity crosses C*Psi = 1/4",
            "NEVER reaches 1/4",
            "r* (kritisch)",
            "quarter_boundary_analysis.png",
            "C·Ψ = ¼",
            "below_rstar = r_vals < R_STAR",
            "pct = 100*(rv<R_STAR).mean()",
        ),
    },
    "simulations/ptf_clock_field.py": {
        "required": ("normalized-purity proxy", "both proxy bands"),
        "forbidden": ("Q3 (the quarter-boundary field)", "Q3 quarter boundary"),
    },
    "simulations/marrakesh_quarter_boundary_review.py": {
        "required": (
            "R_STAR = 0.21275477982200533",
            "is_below_rstar",
            "below-R*",
            "at-or-above-R*",
        ),
        "forbidden": (
            "R_THRESHOLD = 0.213",
            "quantum-side",
            "classical-side",
            "pure-quantum",
            "specific physical threshold",
            "flip_to_quantum",
            "stable_quantum",
        ),
    },
    "simulations/marrakesh_uniform_quantum_chain.py": {
        "required": ("below_rstar_stable", "below_fraction", "below %", "at-or-above-R*"),
        "forbidden": (
            "quantum_stable",
            "uniform-quantum",
            "uniform-classical",
            "quantum-side",
            "\"crossing\":",
            "cross %",
        ),
    },
    "simulations/marrakesh_kingston_fez_compare.py": {
        "required": ("below_rstar_rows", "below-R*", "at-or-above-R*"),
        "forbidden": (
            "stable_below_rstar",
            "r_threshold",
            "stable_quantum",
            "uniform-quantum",
            "uniform-classical",
        ),
    },
    "simulations/qubit_biography.py": {
        "required": (
            "R_STAR = 0.21275477982200533",
            "below_flags = rs < R_STAR",
            "stable-below",
            "stable-at-or-above",
            "r = t2 / (2.0 * t1)",
        ),
        "forbidden": (
            "np.sign(rs - R_STAR)",
            "pulse-stable",
            "silent-stable",
            "classic-stable",
            "quantum-side",
            "classical-side",
        ),
    },
    "simulations/marrakesh_may05_preflight.py": {
        "required": (
            "empirical composite heuristic",
            "volatility penalty",
            "variable-near penalty",
        ),
        "forbidden": (
            "volatility veto",
            "twitch penalty",
            "4. RECOMMENDATION",
            "uniform-quantum",
            "uniform-classical",
        ),
    },
    "simulations/marrakesh_path_biography.py": {
        "required": ("normalized-purity proxy", "R* bands", "below_fraction", "below %"),
        "forbidden": ("CΨ regime identity", "crossing =", "cross %"),
    },
    "simulations/f88b_lens_ibm_kingston_uniform_quantum.py": {
        "required": ("predominantly below-R*", "at-or-above-R*", "confounded association"),
        "forbidden": (
            "all below-R*",
            "all below R*",
            "all-below-R*",
            "uniform-quantum",
            "uniform-classical",
            "quantum-side dephasing dominates",
        ),
    },
    "simulations/kingston_f97_lens.py": {
        "required": ("radial diagnostic", "not the F97 parameter", "c=+1/4"),
        "forbidden": (
            "cross the cardioid boundary",
            "cardioid attractor",
            "F97-anchor crossings",
            "Quarter cusp",
            "Half cardioid",
        ),
    },
}


def source(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8-sig")


def findings(path: str, text: str) -> tuple[str, ...]:
    contract = CONTRACTS[path]
    lowered = text.lower()
    found = []
    for phrase in contract["required"]:
        if phrase.lower() not in lowered:
            found.append(f"{path}: missing positive control: {phrase}")
    for phrase in contract["forbidden"]:
        if phrase.lower() in lowered:
            found.append(f"{path}: restored unsupported label: {phrase}")
    return tuple(found)


HISTORY_LEGACY_KEYS = {
    "r_param",
    "cpsi_min",
    "distance_from_quarter",
    "crosses_quarter",
    "t_star_us",
    "t_star_over_T2",
    "t_coherence_crossing_us",
    "t_coh_over_T2",
}
HISTORY_NEUTRAL_KEYS = {
    "r",
    "proxy_band",
    "normalized_purity_proxy_min",
    "proxy_offset_from_quarter",
    "is_below_rstar",
    "proxy_quarter_time_us",
    "proxy_quarter_time_over_T2",
}

EXACT_RSTAR_HOSTS = {
    "data/ibm_history/ibm_history_analysis.py",
    "simulations/ptf_clock_field.py",
    "simulations/marrakesh_quarter_boundary_review.py",
    "simulations/qubit_biography.py",
}

ACTIVE_RSTAR_BOUNDARIES = {
    "data/ibm_history/ibm_history_analysis.py": (
        "proxy_band",
        'return "below-R*" if r < R_STAR else "at-or-above-R*"',
    ),
    "simulations/marrakesh_quarter_boundary_review.py": (
        "is_below_rstar",
        "return r_param(q) < R_STAR",
    ),
    "simulations/marrakesh_uniform_quantum_chain.py": (
        "main",
        "below_flags = rs < R_STAR",
    ),
    "simulations/marrakesh_kingston_fez_compare.py": (
        "below_rstar_rows",
        "(q.t2_us / (2 * q.t1_us)) < R_STAR",
    ),
    "simulations/qubit_biography.py": (
        "archetype_from_series",
        "below_flags = rs < R_STAR",
    ),
    "simulations/marrakesh_may05_preflight.py": (
        "band_count",
        "by_id[qid].t2_us / (2.0 * by_id[qid].t1_us) < R_STAR",
    ),
    "simulations/marrakesh_path_biography.py": (
        "main",
        "below_flags = rs < R_STAR",
    ),
}


def function(tree: ast.AST, name: str) -> ast.FunctionDef | ast.AsyncFunctionDef | None:
    return next(
        (node for node in ast.walk(tree)
         if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name),
        None,
    )


def literal_dict_keys(node: ast.AST) -> set[str]:
    return {
        key.value
        for candidate in ast.walk(node)
        if isinstance(candidate, ast.Dict)
        for key in candidate.keys
        if isinstance(key, ast.Constant) and isinstance(key.value, str)
    }


def history_schema_findings(text: str) -> tuple[str, ...]:
    tree = ast.parse(text)
    found = []
    compute = function(tree, "compute_qubit_record")
    if compute is None:
        found.append("missing compute_qubit_record")
    else:
        keys = literal_dict_keys(compute)
        missing = HISTORY_NEUTRAL_KEYS - keys
        if missing:
            found.append(f"neutral output keys missing: {sorted(missing)}")
        leaked = HISTORY_LEGACY_KEYS & keys
        if leaked:
            found.append(f"legacy output keys remain active: {sorted(leaked)}")

    names = {
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    retired_functions = {
        "cpsi_from_purity",
        "cpsi",
        "find_quarter_crossing",
        "normalized_purity_proxy",
        "find_proxy_quarter_time",
    }
    leaked_functions = retired_functions & names
    if leaked_functions:
        found.append(f"legacy active functions remain: {sorted(leaked_functions)}")

    adapter = function(tree, "_adapt_legacy_calibration_schema")
    if adapter is None:
        found.append("missing legacy calibration schema adapter")
    elif "LEGACY-CALIBRATION-SCHEMA" not in (ast.get_docstring(adapter) or ""):
        found.append("legacy adapter is not visibly marked")

    for node in tree.body:
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            call = node.value
            if isinstance(call.func, ast.Attribute) and call.func.attr == "mkdir":
                found.append("import-time directory creation remains")
    return tuple(found)


def ptf_writer_findings(text: str) -> tuple[str, ...]:
    tree = ast.parse(text)
    found = []
    main = function(tree, "main")
    if main is None or "output_path" not in {arg.arg for arg in main.args.args}:
        found.append("main(output_path=...) seam missing")
    if function(tree, "build_report") is None:
        found.append("deterministic build_report seam missing")
    adapter = function(tree, "_adapt_archived_proxy_band")
    if adapter is None:
        found.append("archived proxy-band adapter missing")
    elif "LEGACY-CALIBRATION-SCHEMA" not in (ast.get_docstring(adapter) or ""):
        found.append("archived proxy-band adapter is not visibly marked")
    if "--output" not in text:
        found.append("explicit --output CLI missing")
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if any(isinstance(target, ast.Name) and target.id == "RNG" for target in targets):
                found.append("module-global evolving RNG remains")
    return tuple(found)


def assigned_numeric_constant(text: str, name: str) -> float | None:
    tree = ast.parse(text)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == name for target in node.targets
        ):
            if isinstance(node.value, ast.Constant) and isinstance(node.value.value, (int, float)):
                return float(node.value.value)
    return None


def active_rstar_comparisons(text: str, function_name: str) -> tuple[type[ast.cmpop], ...]:
    tree = ast.parse(text)
    owner = function(tree, function_name)
    if owner is None:
        return ()
    comparisons = []
    for node in ast.walk(owner):
        if not isinstance(node, ast.Compare) or len(node.ops) != 1:
            continue
        operands = [node.left, *node.comparators]
        if any(isinstance(item, ast.Name) and item.id == "R_STAR" for item in operands):
            comparisons.extend(type(op) for op in node.ops)
    return tuple(comparisons)


@pytest.mark.parametrize("path", CONTRACTS)
def test_hardware_producer_labels_are_current(path):
    assert findings(path, source(path)) == ()


@pytest.mark.parametrize("path", CONTRACTS)
def test_hardware_producer_contracts_have_live_mutations(path):
    original = source(path)
    for phrase in CONTRACTS[path]["forbidden"]:
        mutant = original + "\n# " + phrase + "\n"
        assert any(phrase in item for item in findings(path, mutant)), (path, phrase)


@pytest.mark.parametrize("path", CONTRACTS)
def test_hardware_producers_remain_syntactically_valid(path):
    ast.parse(source(path), filename=path)


def test_history_writer_uses_neutral_output_schema_and_marked_legacy_adapter():
    text = source("data/ibm_history/ibm_history_analysis.py")
    assert history_schema_findings(text) == ()


def test_history_schema_gate_kills_active_legacy_output_and_unmarked_adapter():
    text = source("data/ibm_history/ibm_history_analysis.py")
    active_legacy = text.replace(
        "'proxy_band': band,",
        "'crosses_quarter': band,",
        1,
    )
    assert active_legacy != text
    assert any("legacy output keys" in item for item in history_schema_findings(active_legacy))
    unmarked = text.replace("LEGACY-CALIBRATION-SCHEMA", "legacy calibration schema", 1)
    assert "legacy adapter is not visibly marked" in history_schema_findings(unmarked)


def test_ptf_writer_has_explicit_deterministic_output_seam():
    text = source("simulations/ptf_clock_field.py")
    assert ptf_writer_findings(text) == ()


def test_ptf_writer_gate_kills_global_rng_and_missing_output_seam():
    text = source("simulations/ptf_clock_field.py")
    global_rng = text + "\nRNG = np.random.default_rng(135)\n"
    assert "module-global evolving RNG remains" in ptf_writer_findings(global_rng)
    missing_seam = text.replace("def main(output_path", "def main(destination", 1)
    assert "main(output_path=...) seam missing" in ptf_writer_findings(missing_seam)
    unmarked = text.replace("LEGACY-CALIBRATION-SCHEMA", "legacy calibration schema", 1)
    assert "archived proxy-band adapter is not visibly marked" in ptf_writer_findings(unmarked)


@pytest.mark.parametrize("path", sorted(EXACT_RSTAR_HOSTS))
def test_active_rstar_constants_are_exact_and_mutation_sensitive(path):
    text = source(path)
    assert assigned_numeric_constant(text, "R_STAR") == 0.21275477982200533
    mutant = text.replace("0.21275477982200533", "0.212754779822", 1)
    assert mutant != text
    assert assigned_numeric_constant(mutant, "R_STAR") != 0.21275477982200533


@pytest.mark.parametrize("path", sorted(ACTIVE_RSTAR_BOUNDARIES))
def test_active_rstar_boundaries_use_strict_less_than_and_reject_equality_mutation(path):
    text = source(path)
    function_name, fragment = ACTIVE_RSTAR_BOUNDARIES[path]
    assert fragment in text, path
    ops = active_rstar_comparisons(text, function_name)
    assert ast.Lt in ops, path
    assert ast.LtE not in ops, path
    mutated_fragment = fragment.replace("< R_STAR", "<= R_STAR", 1)
    mutant = text.replace(fragment, mutated_fragment, 1)
    assert mutant != text
    mutant_ops = active_rstar_comparisons(mutant, function_name)
    assert ast.LtE in mutant_ops, path


def test_history_runtime_schema_uses_exact_raw_rstar_classification():
    import importlib.util

    path = ROOT / "data/ibm_history/ibm_history_analysis.py"
    spec = importlib.util.spec_from_file_location("task6_ibm_history", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    record = module.compute_qubit_record("fixture", 1, 0.5, module.R_STAR)
    assert record["r"] == pytest.approx(module.R_STAR)
    assert record["proxy_band"] == "at-or-above-R*"
    assert record["is_below_rstar"] is False
    assert HISTORY_LEGACY_KEYS.isdisjoint(record)

    legacy = {
        "date": "fixture",
        "qubit": "1",
        "T1_us": "0.5",
        "T2_us": repr(module.R_STAR),
        "frequency_GHz": "5.0",
        "r_param": "0.212755",
        "cpsi_min": "0.25",
        "distance_from_quarter": "0.0",
        "crosses_quarter": "True",
        "t_star_us": "0.1",
        "t_star_over_T2": "0.2",
        "t_coherence_crossing_us": "0.3",
        "t_coh_over_T2": "1.386294",
    }
    adapted = module._adapt_legacy_calibration_schema(legacy)
    # The archive rounded T1/T2 after computing its flag.  A contradictory
    # equality-looking fixture must therefore carry the archived raw-derived
    # flag, not silently reclassify rounded display columns.
    assert adapted["proxy_band"] == "below-R*"
    assert adapted["is_below_rstar"] is True
    assert adapted["band_provenance"] == "legacy-raw-derived-flag"
    assert HISTORY_LEGACY_KEYS.isdisjoint(adapted)


def test_ptf_writer_repeats_in_fresh_processes_without_touching_tracked_output(tmp_path):
    script = ROOT / "simulations/ptf_clock_field.py"
    tracked = ROOT / "simulations/results/clock_field/ptf_clock_field_out.txt"
    tracked_before = hashlib.sha256(tracked.read_bytes()).hexdigest()
    outputs = []
    for name in ("a", "b"):
        destination = tmp_path / name / "ptf.txt"
        destination.parent.mkdir()
        subprocess.run(
            [sys.executable, str(script), "--output", str(destination)],
            cwd=tmp_path,
            check=True,
            capture_output=True,
            text=True,
        )
        assert tuple(destination.parent.iterdir()) == (destination,)
        outputs.append(destination.read_bytes())
    assert outputs[0] == outputs[1]
    assert outputs[0].endswith(b"\n") and not outputs[0].endswith(b"\n\n")
    assert b"\r" not in outputs[0]
    assert str(ROOT).encode() not in outputs[0]
    assert hashlib.sha256(tracked.read_bytes()).hexdigest() == tracked_before
