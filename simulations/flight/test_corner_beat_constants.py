"""Tests for the corner-beat CONSTANTS SCHEMA (the freeze gate).

Honest provenance: unlike test_corner_beat_verdict.py, these are
tests written AFTER the code they cover. `_value_admissible` and its
value-form table grew under review pressure across rounds 13, 18 and
19, and the repo's rule for that case is to add tests for existing
code rather than pretend to a red-green cycle that did not happen.
They are characterization tests plus one parity test that is the real
point: the form table lives in TWO places, section 8a of the
pre-registration and VALUE_FORMS in the runner, and the arc's whole
pattern is a rule stated twice and repaired once.

Run: python -m pytest test_corner_beat_constants.py -q
"""
import hashlib
import importlib.util
import pathlib

import pytest

HERE = pathlib.Path(__file__).parent
DOC = pathlib.Path(
    r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
    r"\experiments\CORNER_BEAT_HARDWARE_PREDICTION.md")


def _load_runner():
    spec = importlib.util.spec_from_file_location(
        "rcb", HERE / "run_corner_beat.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rcb = _load_runner()


# ------------------------------------------------------- the forms
@pytest.mark.parametrize("key,value,ok", [
    # scalars, the default for every key not named in the table
    ("theta_D", 0.00253, True),
    ("theta_D", "0.00253", False),      # round 13: a string reached the fits
    ("theta_D", True, False),           # a bool is not a measurement
    ("theta_D", [0.001, 0.003], False),
    # round 19: every band and margin is a HALF-WIDTH, hence scalar
    ("dmag_band_C", 0.001, True),
    ("dmag_band_C", [0.001, 0.003], False),
    ("J_pass_band", 0.02, True),
    ("A_margin_C", 0.01, True),
    ("T1_clean_station_band", 0.05, True),
    # R_BOOT is a replicate COUNT
    ("R_BOOT", 2000, True),
    ("R_BOOT", 2000.5, False),
    ("R_BOOT", 0, False),
    # the profile is one number per deep cell
    # Amendment 1.8 moved the worst deep cell 87 -> 208 kept; the form
    # test does not read the value, but a stale number in test data is
    # how the next reader learns the wrong one
    ("deep_end_kept_profile", [208.0, 400.0, 900.0], True),
    ("deep_end_kept_profile", [], False),
    ("deep_end_kept_profile", 208.0, False),
    # named reals
    ("family_error_rates", {"P_any_void": 0.05, "P_false_confirmed": 0.01}, True),
    ("family_error_rates", 0.05, False),
    ("N0_clean_thresholds", {"spectrum": 1.0, "residual": 0.1, "batch": 0.2}, True),
    # the one dimensioned key: named axes, each an interval
    ("band_validity_window", {"T2echo_us": [150, 400], "p2": [0.0, 0.005]}, True),
    ("band_validity_window", {"T2echo_us": 150}, False),
    # the sentinel, on its two registered keys and no others
    ("dmag_band_Cp", "REPORTED", True),
    ("f_Cp_dressing", "REPORTED", True),
    ("A_margin_Cp", "REPORTED", False),   # round 19 blocker: never registered
    ("theta_D", "REPORTED", False),
    # the axis enumeration, closed since round 14
    ("time_axis", "nominal", True),
    ("time_axis", "realized_dose", True),
    ("time_axis", "realised_dose", False),
])
def test_value_admissibility(key, value, ok):
    assert rcb._value_admissible(key, value) is ok


def test_a_code_form_must_hash_to_the_committed_file():
    """Round 19: this accepted {"source":"todo","sha256":"todo"} at the
    last gate before a paid submission."""
    src = "run_corner_beat.py"
    good = hashlib.sha256((HERE / src).read_bytes()).hexdigest()
    assert rcb._value_admissible(
        "f_C_dressing", {"source": src, "symbol": "f_C", "sha256": good})
    assert not rcb._value_admissible(
        "f_C_dressing", {"source": src, "symbol": "f_C", "sha256": "0" * 64})
    assert not rcb._value_admissible(
        "f_C_dressing", {"source": src, "symbol": "f_C", "sha256": "todo"})
    assert not rcb._value_admissible(
        "f_C_dressing", {"source": "no_such_file.py", "symbol": "f",
                         "sha256": good})


def test_the_dressed_centres_cannot_be_frozen_as_numbers():
    """Section 7, verbatim: 'a frozen NUMBER cannot exist pre-flight'.
    Round 18's blocker was the guard demanding one anyway."""
    for key in ("f_C_dressing", "A_centres_dressed", "pi3_centre_dressed",
                "dayof_width_scaling"):
        assert not rcb._value_admissible(key, 0.5)


# ------------------------------------------------------- doc parity
def test_every_manifest_key_has_a_reader_or_a_ledger_entry():
    """The rule the ledger exists for: a frozen number without its
    evaluator must never release submission."""
    src = (HERE / "run_corner_beat.py").read_text(encoding="utf-8")
    # Cut the two declaration blocks out first. Without this every key
    # finds itself in its own declaration and the test is a fake green
    # (which it was, for one run, until the walk below caught it).
    def _strip_block(text, marker, closer):
        i = text.index(marker)
        j = text.index(closer, i)
        return text[:i] + text[j + len(closer):]

    body = _strip_block(src, "CONSTANTS_MANIFEST = [", "\n]")
    body = _strip_block(body, "DEFAULT_CONSTANTS = {", "\n}")
    ledger = " ".join(rcb.MACHINERY_PENDING)
    orphans = []
    for key in rcb.CONSTANTS_MANIFEST:
        if key in ledger or f'"{key}"' in body:
            continue
        orphans.append(key)
    assert orphans == [], f"manifest keys with no reader and no ledger entry: {orphans}"


def test_the_form_table_on_the_page_matches_the_form_table_in_the_code():
    """The arc's signature failure is a rule stated twice and repaired
    once. Section 8a names the non-scalar keys; VALUE_FORMS assigns
    them. They must agree, or the executor freezing by the document
    is stopped by the code."""
    if not DOC.exists():
        pytest.skip(
            "pre-registration not reachable; set CORNER_BEAT_DOC to run "
            "the parity check (a silent skip here is a green that means "
            "nothing, which is why the machine-independent half below "
            "runs everywhere)")
    text = DOC.read_text(encoding="utf-8")
    non_scalar = {k: v for k, v in rcb.VALUE_FORMS.items()
                  if v not in ("scalar",)}
    missing = [k for k in non_scalar if k not in text]
    assert missing == [], f"keys typed non-scalar in code, unnamed on the page: {missing}"
    # presence is not agreement: the FORM word must sit near the key,
    # or a key documented scalar and coded map would pass (review
    # finding: the check could detect absence only).
    disagreeing = []
    for key, form in non_scalar.items():
        # every occurrence, not the first: a key is named in the owed
        # lists long before the form paragraph, and taking text.index
        # measured distance to the wrong mention (my own defect, found
        # by this test's first run).
        starts = [i for i in range(len(text))
                  if text.startswith(key, i)]
        if not any(form in text[max(0, i - 1500):i + 1500] for i in starts):
            disagreeing.append((key, form))
    assert disagreeing == [], (
        f"keys whose code form is not stated beside them on the page: "
        f"{disagreeing}")


def test_every_typed_form_key_is_a_manifest_key():
    """Machine-independent half of the parity check: a form assigned to
    a key that is not in the manifest is a typo that would silently
    never apply."""
    stray = [k for k in rcb.VALUE_FORMS if k not in rcb.CONSTANTS_MANIFEST]
    assert stray == [], f"forms typed for non-manifest keys: {stray}"


def test_no_manifest_key_is_declared_twice_in_the_form_table():
    """The arc's signature failure, a rule stated twice and repaired
    once, once sat inside this very table: band_validity_window was
    declared 'map' and then 'map_of_intervals', and Python kept the
    later one."""
    import re
    src = (HERE / "run_corner_beat.py").read_text(encoding="utf-8")
    block = src[src.index("VALUE_FORMS = {"):]
    block = block[:block.index("\n}")]
    keys = re.findall(r'^\s*"([A-Za-z0-9_]+)":', block, re.M)
    dupes = sorted({k for k in keys if keys.count(k) > 1})
    assert dupes == [], f"declared twice in VALUE_FORMS: {dupes}"


def test_default_constants_and_the_shipped_file_agree_on_keys():
    import json
    shipped = json.load(open(HERE / "corner_beat_constants.json"))
    a = sorted(k for k in rcb.DEFAULT_CONSTANTS if k != "_source")
    b = sorted(k for k in shipped if k != "_source")
    assert a == b


def test_every_shipped_value_passes_its_own_form():
    """A shipped file that cannot pass the gate is a trap for the
    executor who freezes on flight day."""
    import json
    shipped = json.load(open(HERE / "corner_beat_constants.json"))
    bad = [k for k, e in shipped.items()
           if k != "_source" and e.get("value") is not None
           and not rcb._value_admissible(k, e["value"])]
    assert bad == [], f"shipped values failing their form: {bad}"


# ------------------------------------------- the p2 Class-1 guard bound
# Amendment 1.8 (Tom, 2026-08-17) tightened the bound 0.6% -> 0.5%. The
# inventory taken before that edit found the bound living as a bare
# module constant: absent from the manifest, absent from the sync-check
# loop, and stated a SECOND time as a hardcoded "0.6%" inside the
# hard-abort message. So the one-line change would have produced an
# abort that names the wrong rule, which is this arc's signature defect
# sitting on the very number being moved. These tests close both halves.
def test_the_p2_guard_bound_is_a_registered_manifest_key():
    """A threshold that hard-aborts a paid submission is a registered
    threshold, not a code constant, the same ruling round 12 made for
    kept_count_floor."""
    assert "p2_guard_bound" in rcb.CONSTANTS_MANIFEST
    assert "p2_guard_bound" in rcb.DEFAULT_CONSTANTS
    assert (rcb.DEFAULT_CONSTANTS["p2_guard_bound"]["value"]
            == rcb.RULE_MAX_P2_RZZ)


def test_the_p2_guard_bound_is_sync_checked_against_the_code():
    """Freezing the entry to a value the guard does not read is exactly
    the gap time_axis had: the manifest would record the tightening as
    landed while the guard still passed lines at the old bound."""
    src = (HERE / "run_corner_beat.py").read_text(encoding="utf-8")
    assert '("p2_guard_bound", RULE_MAX_P2_RZZ)' in src


def test_the_day_of_rule_message_states_the_bound_from_the_constant():
    """The second copy. The abort message must derive every threshold it
    names from the constant that enforces it, so no edit can leave the
    executor with an abort naming a rule the code no longer applies."""
    text = rcb._day_of_rule_text()
    assert f"{rcb.RULE_MAX_P2_RZZ * 100:g}%" in text
    assert f"{rcb.RULE_MIN_T2_US:g}" in text
    assert f"{rcb.RULE_MAX_READOUT * 100:g}%" in text
    assert f"{rcb.RULE_MAX_T2_RATIO:g}" in text


def test_no_percentage_literal_survives_beside_the_p2_guard():
    """Characterization of the fix: the rule text is built once. A
    literal percentage anywhere in the day-of abort is a second copy
    waiting to go stale."""
    import re
    src = (HERE / "run_corner_beat.py").read_text(encoding="utf-8")
    body = src[src.index("def _day_of_rule_text("):]
    body = body[:body.index("\ndef ", 10)]
    assert not re.search(r'"[^"]*\d\.\d\s*%', body), (
        "a hardcoded percentage inside the rule-text builder")


# ------------------------------- the weighted line budget (2026-08-17)
# The Class-1 guard is a MAXIMUM over the used edges; the quantity it
# protects is a WEIGHTED SUM, because a Strang step gives the odd bonds
# two XY blocks and the even bonds one. These tests pin the weights at
# the flown circuit, and pin the one structural fact that decides how
# far the budget may be used: under the 0.5% ceiling it can never gate.
def test_edge_weights_are_the_flown_two_qubit_gate_count():
    w = rcb.edge_gate_weights()
    assert [w[b] for b in rcb.CHAIN_BONDS] == [242, 122, 242, 122, 242]
    # 60 steps x (odd: 2 blocks, even: 1 block) x 2 rzz, plus 2 rzz of
    # preparation on every bond (5 Givens on overlapping pairs)
    assert sum(w.values()) == 970


def test_the_budget_of_a_uniform_line_is_the_documented_exponent():
    line = {b: 0.005 for b in rcb.CHAIN_BONDS}
    assert abs(rcb.line_error_budget(line) - 4.85) < 1e-12
    kept = rcb.projected_deep_end_kept(line)
    assert abs(kept - 208.28) < 0.05
    assert abs(kept / rcb.KEPT_COUNT_FLOOR - 4.17) < 0.01


def test_the_three_x_requirement_written_as_a_budget():
    """>= 3x is exactly an upper bound on the weighted budget, and on a
    uniform line it is the 0.538% the pre-registration states."""
    assert abs(rcb.P2_BUDGET_MAX - 5.2150) < 1e-3
    uniform_equivalent = rcb.P2_BUDGET_MAX / 970
    assert abs(uniform_equivalent - 0.005376) < 1e-6


def test_the_budget_can_never_gate_under_the_max_ceiling():
    """The reason no budget GUARD is built. A line passing the per-edge
    ceiling has budget <= 970*RULE_MAX_P2_RZZ, which is below the >= 3x
    bound with room to spare, so a budget guard would be machinery that
    cannot fire. If a future amendment raises the ceiling past this
    ratio, this test fails and the guard becomes real."""
    assert 970 * rcb.RULE_MAX_P2_RZZ < rcb.P2_BUDGET_MAX


def test_the_budget_orders_lines_the_unweighted_mean_gets_wrong():
    """Why the score changes. Two admissible lines where mean(p2) and
    the budget disagree: the mean prefers P, the deep end prefers Q."""
    P = dict(zip(rcb.CHAIN_BONDS, [0.004, 0.0005, 0.004, 0.0005, 0.004]))
    Q = dict(zip(rcb.CHAIN_BONDS, [0.0025, 0.0035, 0.0025, 0.0035, 0.0025]))
    import statistics
    assert statistics.mean(P.values()) < statistics.mean(Q.values())
    assert rcb.line_error_budget(P) > rcb.line_error_budget(Q)
    assert rcb.projected_deep_end_kept(Q) > rcb.projected_deep_end_kept(P)


def test_a_bad_odd_edge_costs_twice_a_bad_even_edge():
    D = dict(zip(rcb.CHAIN_BONDS, [0.005, 0.001, 0.005, 0.001, 0.005]))
    E = dict(zip(rcb.CHAIN_BONDS, [0.001, 0.005, 0.001, 0.005, 0.001]))
    assert max(D.values()) == max(E.values()) == rcb.RULE_MAX_P2_RZZ
    assert (rcb.projected_deep_end_kept(E)
            > 5 * rcb.projected_deep_end_kept(D))


# ------------------------------- the ledger's own promise (2026-08-17)
# Section 8a promises that "every machinery entry NAMES the manifest
# keys it discharges, in brackets, so that 'empty the ledger by its
# names' is an operation an executor can actually perform". When that
# promise was written the brackets were added to five entries of
# twenty-four, and twelve of the nineteen unbracketed ones plainly
# discharge a key. An executor emptying the ledger by name would have
# missed them, which is the same defect round 10 named and the
# composer's review round found again one level down.
def test_every_owed_manifest_key_is_reachable_from_the_ledger():
    """Either a ledger entry names the key in brackets, or the key's
    consumer already EXISTS and is declared as such. No third state:
    a key in neither set is a number nobody has to build anything for,
    which is how a freeze releases a flight with a dead arm."""
    import re
    named = set()
    for entry in rcb.MACHINERY_PENDING:
        for grp in re.findall(r"\[([^\]]+)\]", entry):
            named.update(k.strip() for k in grp.split(","))
    unreachable = [k for k in rcb.CONSTANTS_MANIFEST
                   if k not in named and k not in rcb.KEYS_WITH_LIVE_CONSUMERS]
    assert unreachable == [], (
        f"manifest keys named by no ledger entry and not declared as "
        f"having a live consumer: {unreachable}")


def test_the_live_consumer_declaration_holds_only_real_keys():
    stray = [k for k in rcb.KEYS_WITH_LIVE_CONSUMERS
             if k not in rcb.CONSTANTS_MANIFEST]
    assert stray == [], f"declared live consumers for non-manifest keys: {stray}"


def test_bracketed_ledger_names_are_real_keys():
    """A bracket naming a key that does not exist discharges nothing
    and reads as if it did (the composer round found `thetaF` for
    `theta_F` and `null_margin_12` spelled two ways)."""
    import re
    bad = []
    for entry in rcb.MACHINERY_PENDING:
        for grp in re.findall(r"\[([^\]]+)\]", entry):
            for k in (x.strip() for x in grp.split(",")):
                if k not in rcb.CONSTANTS_MANIFEST:
                    bad.append((entry.split(" [")[0], k))
    assert bad == [], f"ledger brackets naming non-manifest keys: {bad}"


# --------------------------- the fit time axis, gate side (2026-08-17)
# Machinery ledger item `realized_dose_time_axis_refit`, and the first
# step of the theta_D refreeze: section 5 registers the refit on the
# REALIZED-DOSE axis t_eff = max(n-1,0)*dt "in gate and runner
# together", the runner derives its axis from its TIME_AXIS constant,
# and the gate built the nominal axis inline in four places. --certify
# handed ONE array to both fit implementations, so it certified that
# the two FITS agree on a given axis and never that the gate
# CONSTRUCTS the same one. These tests are that missing check.
def _load_gate():
    spec = importlib.util.spec_from_file_location("cbgate", rcb.GATE_PATH)
    g = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g)
    return g


def test_the_gate_derives_its_fit_axis_from_a_switch():
    g = _load_gate()
    assert hasattr(g, "TIME_AXIS")
    assert hasattr(g, "fit_time_axis")


def test_gate_and_runner_build_the_same_axis_on_both_settings():
    """The parity that was missing. Exact equality, not a tolerance:
    both sides compute the same product of the same integers, so a
    residual would be a defect in the construction, not rounding."""
    import numpy as np
    g = _load_gate()
    saved = rcb.TIME_AXIS
    try:
        for axis in ("nominal", "realized_dose"):
            rcb.TIME_AXIS = axis
            mine = rcb.fit_time_axis()
            theirs = g.fit_time_axis(rcb.STEPS_MAX, rcb.K_GRID, rcb.JDT,
                                     axis=axis)
            assert np.array_equal(mine, theirs), (
                f"{axis}: gate {theirs} vs runner {mine}")
    finally:
        rcb.TIME_AXIS = saved


def test_the_gate_refactor_reproduces_the_axis_it_replaced():
    """The four inline copies read [i*k*dt for i in range(steps//k+1)].
    An exact route exists, so compare exactly (house rule)."""
    import numpy as np
    g = _load_gate()
    # every (steps, k, dt) the gate sweeps: 2 working points x 2
    # spacings, endpoint fixed at 60 in all sixteen configurations
    for steps, k, dt in ((60, 3, 0.15), (60, 5, 0.15),
                         (60, 3, 0.10), (60, 5, 0.10)):
        old = np.array([i * k * dt for i in range(steps // k + 1)])
        assert np.array_equal(g.fit_time_axis(steps, k, dt, axis="nominal"),
                              old)


def test_the_realized_dose_axis_moves_only_the_dosed_points():
    """Depth 0 carries no injection layer at all, so it stays at 0 and
    every dosed point loses exactly one layer (Amendment 1.4)."""
    import numpy as np
    g = _load_gate()
    nom = g.fit_time_axis(60, 3, 0.15, axis="nominal")
    real = g.fit_time_axis(60, 3, 0.15, axis="realized_dose")
    assert real[0] == nom[0] == 0.0
    assert np.allclose(nom[1:] - real[1:], 0.15)
