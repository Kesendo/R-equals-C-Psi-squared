"""Independent Bell+ crossing references and the Commit-A production-seam RED door."""
from __future__ import annotations

import ast
from collections import defaultdict
from contextlib import redirect_stdout
import inspect
import io
from pathlib import Path
import re
import runpy

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]


def _format_red_contract(marker, findings_by_category):
    normalized = {}
    for category, findings in findings_by_category.items():
        assert re.fullmatch(r"[A-Z][A-Z0-9_]*", category)
        details = tuple(sorted({str(finding) for finding in findings if str(finding)}))
        if details:
            normalized[category] = details
    categories = tuple(sorted(normalized))
    summary = f"{marker} categories={','.join(categories)} category_count={len(categories)}"
    details = tuple(f"{category}: {finding}" for category in categories for finding in normalized[category])
    return summary, details


def fail_red_contract(marker, findings_by_category):
    """Fail with one ASCII-safe printed line per finding (plain Windows console safe)."""
    summary, details = _format_red_contract(marker, findings_by_category)
    for detail in details:
        print(f"{marker}_DETAIL {detail}".encode("ascii", errors="backslashreplace").decode("ascii"))
    assert not details, summary


def test_red_contract_detail_printing_is_ascii_safe_and_identifiable(capsys):
    with pytest.raises(AssertionError, match="categories=WRITER"):
        fail_red_contract("RED_CONTRACT_TEST", {
            "WRITER": ["plot.py:VISIBLE:CΨ_com"],
        })
    output = capsys.readouterr().out
    assert output.isascii()
    assert "RED_CONTRACT_TEST_DETAIL WRITER: plot.py:VISIBLE:" in output
    assert r"C\u03a8_com" in output
    assert "CΨ_com" not in output


def test_red_contract_formatter_reassignment_changes_category_identity():
    before = _format_red_contract("PROBE", {"ALPHA": ["one"], "BETA": ["two"], "EMPTY": [""]})
    after = _format_red_contract("PROBE", {"WRONG": ["one"], "BETA": ["two"]})
    assert before == ("PROBE categories=ALPHA,BETA category_count=2",
                      ("ALPHA: one", "BETA: two"))
    assert after == ("PROBE categories=BETA,WRONG category_count=2",
                     ("BETA: two", "WRONG: one"))
    assert before != after


def test_red_contract_formatter_rejects_a_category_that_is_not_one_token():
    with pytest.raises(AssertionError):
        _format_red_contract("PROBE", {"shared book": ["one"]})


SOURCE = ROOT / "simulations/crossing_taxonomy_books.py"
EXPECTED = {
    "mutual_info": {
        "clean_K": 0.02965683339109007035,
        "feedback_K": 0.03264460387550044147,
    },
    "concurrence": {
        "clean_K": 0.03596025905647261593,
        "feedback_K": 0.03867513459481288225,
    },
    "correlation": {
        "clean_K": 0.07192051811294523186,
        "feedback_K": 0.07192051811294523186,
    },
}
NEVER = ("mutual_purity", "overlap")
GAMMA = 0.05


def copied_book(tmp_path, *, gamma=GAMMA, tighter=False, source=None):
    """Full source copy, run under a non-main name, with a disposable output cwd."""
    tree = ast.parse(SOURCE.read_text(encoding="utf-8") if source is None else source)
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "GAMMA" for t in node.targets):
            node.value = ast.Constant(gamma)
    # Before the parameter seam lands the convergence control changes only the
    # real solve_ivp call's tolerance literals in this temporary source copy.
    if tighter:
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "solve_ivp":
                for keyword in node.keywords:
                    if keyword.arg in ("rtol", "atol") and isinstance(keyword.value, ast.Constant):
                        keyword.value = ast.Constant(1e-13 if keyword.arg == "rtol" else 1e-15)
    target = tmp_path / "crossing_source_copy.py"
    target.write_text(ast.unparse(ast.fix_missing_locations(tree)), encoding="utf-8", newline="\n")
    return runpy.run_path(str(target))


def validate_bridge_definitions(mapping):
    for f in (0., .25, .5, .875, 1.):
        eigenvalues = np.array([(1+f)/2, (1-f)/2])
        nz = eigenvalues[eigenvalues > 0]
        entropy = -np.sum(nz*np.log2(nz))
        independent = {"mutual_info": (2-entropy)/2, "concurrence": f,
                       "correlation": min(1., .5 + f*f), "mutual_purity": .5, "overlap": .25}
        assert set(mapping) == set(independent)
        for name, expected in independent.items():
            assert abs(mapping[name](f)-expected) <= 2e-16, (name, f)


def actual_book(module, book, bridges=None):
    if "crossings_for_book" in module:
        return module["crossings_for_book"](book, bridges=module["BRIDGES"] if bridges is None else bridges)
    # Base-only numerical door: the missing shared door is separately RED.
    function = module[book+"_crossing"]
    old = function.__globals__["BRIDGES"]
    function.__globals__["BRIDGES"] = old if bridges is None else bridges
    try:
        return {name: function(name) for name in old}
    finally:
        function.__globals__["BRIDGES"] = old


def validate_crossing_book(actual, expected, book, *, gamma=GAMMA):
    assert set(actual) == set(expected) | set(NEVER)
    for name, references in expected.items():
        k, t = actual[name]
        budget = 1e-13 if book == "clean" or name == "correlation" else 1e-10
        assert abs(k-references[book+"_K"]) <= budget, (book, name, k)
        assert abs(t-k/gamma) <= 2e-9, (book, name, t)
    for name in NEVER:
        assert actual[name] == (None, None)


@pytest.mark.parametrize("book", ("clean", "feedback"))
def test_all_crossing_reference_values_and_never_bridges(book, tmp_path):
    module = copied_book(tmp_path)
    validate_bridge_definitions(module["BRIDGES"])
    validate_crossing_book(actual_book(module, book), EXPECTED, book)


@pytest.mark.parametrize("gamma", (.025, .05, .1))
@pytest.mark.parametrize("book", ("clean", "feedback"))
def test_gamma_scaling_stays_inside_fixed_bellplus_book(gamma, book, tmp_path):
    module = copied_book(tmp_path, gamma=gamma)
    validate_crossing_book(actual_book(module, book), EXPECTED, book, gamma=gamma)


def test_default_and_tighter_feedback_converge_empirically(tmp_path):
    module = copied_book(tmp_path)
    default = module["feedback_crossing"]
    if {"rtol", "atol"} <= set(inspect.signature(default).parameters):
        tightened = lambda name: default(name, rtol=1e-13, atol=1e-15)
    else:
        tight_module = copied_book(tmp_path, tighter=True)
        tightened = tight_module["feedback_crossing"]
    for name in EXPECTED:
        first, second = default(name)[0], tightened(name)[0]
        assert abs(first-second) <= 2e-11, name


@pytest.mark.parametrize("bridge", tuple(EXPECTED) + NEVER)
@pytest.mark.parametrize("book", ("clean", "feedback"))
def test_mutated_real_bridge_mapping_fails_fixed_reference_gate(bridge, book, tmp_path):
    module = copied_book(tmp_path)
    check_bridge_outcome_mutation(module, bridge, book)


def check_bridge_outcome_mutation(module, bridge, book):
    original = module["BRIDGES"]
    validate_bridge_definitions(original)
    mutated = dict(original)
    # Each legacy label's distinct object is sent through the real consumer door.
    alternatives = {"mutual_info": lambda f: 2-f,
                    "concurrence": lambda f: f*f,
                    "correlation": lambda f: f*f,
                    "mutual_purity": lambda f: 1.,
                    "overlap": lambda f: (1+f)/2}
    mutated[bridge] = alternatives[bridge]
    actual = actual_book(module, book, mutated)
    with pytest.raises(AssertionError):
        validate_crossing_book(actual, EXPECTED, book)


def test_outcome_gate_rejects_real_feedback_consumer_ignoring_injection(tmp_path):
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    cpsi = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "cpsi")
    frozen_cpsi = ast.parse(ast.unparse(cpsi)).body[0]
    frozen_cpsi.name = "_frozen_cpsi"

    class FreezeBridgeMapping(ast.NodeTransformer):
        def visit_Name(self, node):
            if node.id in {"BRIDGES", "bridges"}:
                return ast.copy_location(ast.Name(id="_frozen_bridges", ctx=node.ctx), node)
            if node.id == "cpsi":
                return ast.copy_location(ast.Name(id="_frozen_cpsi", ctx=node.ctx), node)
            return node

    feedback = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "feedback_crossing")
    FreezeBridgeMapping().visit(feedback)
    FreezeBridgeMapping().visit(frozen_cpsi)
    tree.body.extend(ast.parse("_frozen_bridges = dict(BRIDGES)").body + [frozen_cpsi])
    mutant = copied_book(tmp_path, source=ast.unparse(ast.fix_missing_locations(tree)))
    # Only the real feedback consumer is frozen; clean still consumes injection.
    # A contract checking clean alone silently accepts this broken consumer.
    with pytest.raises(pytest.fail.Exception, match="DID NOT RAISE"):
        check_bridge_outcome_mutation(mutant, "concurrence", "feedback")


@pytest.mark.parametrize("bridge", tuple(EXPECTED) + NEVER)
def test_mutated_real_bridge_definition_fails_definition_gate(bridge, tmp_path):
    module = copied_book(tmp_path)
    mutated = dict(module["BRIDGES"])
    alternatives = {"mutual_info": lambda f: 2-f,
                    "concurrence": lambda f: f*f,
                    "correlation": lambda f: f*f,
                    "mutual_purity": lambda f: .25,
                    "overlap": lambda f: (1+f)/2}
    mutated[bridge] = alternatives[bridge]
    with pytest.raises(AssertionError):
        validate_bridge_definitions(mutated)
    if bridge == "mutual_purity":
        # This wrong definition has the same never-crossing outcome. It belongs
        # exclusively to the definition door, not the outcome-mutation control.
        validate_crossing_book(actual_book(module, "clean", mutated), EXPECTED, "clean")


def generated_summary(module):
    stream = io.StringIO()
    with redirect_stdout(stream):
        module["main"]()
    return stream.getvalue()


def label_findings(source, *, summary=None):
    findings = []
    doc = ast.get_docstring(ast.parse(source)) or ""
    surfaces = [("module docstring", doc)]
    if summary is not None:
        surfaces.append(("generated summary", summary))
    for surface, text in surfaces:
        normalized = " ".join(text.split()).lower()
        if "model-independent" in normalized:
            findings.append(f"{surface}: model-independent Lindblad label: {normalized}")
        if "two books" not in normalized and "two named evolution-law families" not in normalized:
            findings.append(f"{surface}: missing two-book scope")
        for wrong in ("the whole feedback book is standard lindblad",
                      "the whole feedback book is linear lindblad",
                      "the whole feedback book is nonlinear"):
            if wrong in normalized:
                findings.append(f"{surface}: {wrong}")
    return findings


@pytest.mark.parametrize("wrong", (
    "The surviving Lindblad-scaling result is model-independent.",
    "The whole feedback book is standard Lindblad.",
    "The whole feedback book is linear Lindblad.",
    "The whole feedback book is nonlinear.",
))
def test_old_model_labels_fail_through_actual_module_docstring(wrong):
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    tree.body[0].value.value = "Two books: "+wrong
    mutated = ast.unparse(tree)
    assert label_findings(mutated)


@pytest.mark.parametrize("wrong", (
    "The surviving Lindblad-scaling result is model-independent.",
    "The whole feedback book is standard Lindblad.",
    "The whole feedback book is linear Lindblad.",
    "The whole feedback book is nonlinear.",
))
def test_old_model_labels_fail_through_actual_generated_summary(wrong, tmp_path):
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    tree.body[0].value.value = "Two books: finite Bell+ local Z-dephasing equations."
    main = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main")
    main.body.insert(0, ast.Expr(ast.Call(func=ast.Name(id="print", ctx=ast.Load()),
                                        args=[ast.Constant(wrong)], keywords=[])))
    source = ast.unparse(ast.fix_missing_locations(tree))
    module = copied_book(tmp_path, source=source)
    summary = generated_summary(module)
    assert wrong in summary.splitlines()
    assert any(item.startswith("generated summary:") and wrong.lower().rstrip(".") in item
               for item in label_findings(source, summary=summary))


def test_red_contract_requires_shared_book_and_solver_seams(tmp_path):
    module = copied_book(tmp_path)
    findings = defaultdict(list)
    shared = module.get("crossings_for_book")
    if shared is None:
        findings["SHARED_BOOK"].append("crossing_taxonomy_books.py: missing crossings_for_book")
    else:
        for book in ("clean", "feedback"):
            try:
                validate_crossing_book(actual_book(module, book), EXPECTED, book)
            except AssertionError as error:
                findings["SHARED_BOOK"].append(f"{book}: {error}")
    if not {"rtol", "atol"} <= set(inspect.signature(module["feedback_crossing"]).parameters):
        findings["FEEDBACK_SOLVER_TOLERANCE"].append("feedback_crossing lacks explicit rtol/atol parameters")
    findings["CURRENT_LABEL"].extend(label_findings(
        SOURCE.read_text(encoding="utf-8"), summary=generated_summary(module)))
    fail_red_contract("RED_CONTRACT_CROSSING_SEAMS", findings)
