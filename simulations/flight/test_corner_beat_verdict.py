"""Tests for the corner-beat VERDICT-LINE COMPOSER.

Each test encodes one rule of the composition algorithm registered in
section 7 of experiments/CORNER_BEAT_HARDWARE_PREDICTION.md, and cites
the round that established it. The document is the specification; this
file is the executable form of it, built after eight review rounds in
which every round found its defects in the previous round's prose
repairs and none in the physics.

Run: python -m pytest test_corner_beat_verdict.py -q
"""
import pytest

from corner_beat_verdict import TRUE, FALSE, VOID, compose_verdict, classify_arm


def line(value, trigger=None):
    return {"value": value, "trigger": trigger}


def arms(dsign=TRUE, w=TRUE, dmag=None, a=None):
    """The default is the all-holding strong positive, one line per
    single-line arm and two per two-line arm."""
    return {
        "D-sign": {"": line(dsign)},
        "W": {"": line(w)},
        "D-mag": dmag if dmag is not None else {"C": line(TRUE), "C'": line(TRUE)},
        "A": a if a is not None else {"C": line(TRUE), "C'": line(TRUE)},
    }


# --------------------------------------------------------- step 2
def test_arm_with_all_lines_valued_and_one_failing_is_not_partial():
    """Round 17 blocker: PARTIAL is a statement about VOIDNESS only.
    An arm whose lines all valued, one of which fails, is not-arm."""
    assert classify_arm({"C": line(TRUE), "C'": line(FALSE)}) == "VALUED_FALSE"


def test_arm_with_some_lines_void_is_partial():
    assert classify_arm({"C": line(TRUE), "C'": line(VOID, "(2,3) null")}) == "PARTIAL"


def test_arm_with_all_lines_void_is_void():
    assert classify_arm({"C": line(VOID, "dose"), "C'": line(VOID, "dose")}) == "VOID"


def test_valued_arm_holds_only_if_every_line_holds():
    assert classify_arm({"C": line(TRUE), "C'": line(TRUE)}) == "VALUED_TRUE"


# --------------------------------------------------------- step 4
def test_all_four_arms_valued_and_holding_is_the_only_route_to_confirmed():
    assert compose_verdict(arms()) == "CONFIRMED"


def test_a_void_arm_forbids_the_registered_name_and_the_line_enumerates():
    """Round 18/19: no partial pattern is ever promoted. The
    enumeration IS the line, in the reporting order, and D-mag has one
    line under the reverted-band freeze this document expects."""
    state = arms(
        w=VOID,
        dmag={"C": line(TRUE)},
        a={"C": line(TRUE), "C'": line(VOID, "(2,3) null")},
    )
    state["W"][""] = line(VOID, "(2,3) null")
    assert compose_verdict(state) == (
        "D-sign detected; W VOID ((2,3) null); D-mag holds; "
        "A PARTIAL (C line holds; C′ line VOID ((2,3) null))"
    )


def test_partial_form_names_the_surviving_lines_truth_value():
    """Round 18: without it, a partial A with a holding C line and one
    with a failing C line printed the same string."""
    state = arms(a={"C": line(FALSE), "C'": line(VOID, "dose certificate")})
    assert "A PARTIAL (C line fails; C′ line VOID (dose certificate))" in compose_verdict(state)


def test_a_failing_two_line_arm_may_carry_its_line_verdicts():
    """The form is only reachable where NO registered name covers the
    pattern: with all four arms valued, D-sign holding and A failing,
    section 7 names the line 'Split observed, anchor failed'. A void
    elsewhere drops it to the enumeration, which is where the arm's
    own line verdicts get printed."""
    state = arms(a={"C": line(TRUE), "C'": line(FALSE)})
    state["W"] = {"": line(VOID, "(2,3) null")}
    assert "¬A (C′ line fails; C line holds)" in compose_verdict(state)


def test_a_named_line_wins_over_the_enumeration_when_all_arms_are_valued():
    state = arms(a={"C": line(TRUE), "C'": line(FALSE)})
    assert compose_verdict(state) == "Split observed, anchor failed"


# --------------------------------------------------------- step 1
def test_dmag_has_one_line_when_its_Cprime_band_reverted_to_reported():
    """Round 19 blocker 3: a (2,3) null cannot make D-mag PARTIAL in
    the regime where D-mag has no C' line at all."""
    state = arms(dmag={"C": line(TRUE)})
    state["W"][""] = line(VOID, "(2,3) null")
    state["A"]["C'"] = line(VOID, "(2,3) null")
    assert "D-mag holds" in compose_verdict(state)
    assert "D-mag PARTIAL" not in compose_verdict(state)


def test_a_single_line_arm_can_never_be_partial():
    assert classify_arm({"": line(VOID, "x")}) == "VOID"


# --------------------------------------------------------- step 3
def test_a_global_trigger_gives_a_flight_level_void_naming_every_trigger():
    """Round 18: step 3 appends ALL firing triggers, in precedence
    order, inside the parenthesis (round 19 bracketing)."""
    got = compose_verdict(arms(), global_triggers=["grid incomplete", "N0-CLEAN"])
    assert got == "VOID (N0-CLEAN; grid incomplete)"


def test_the_leading_trigger_is_chosen_by_the_registered_precedence():
    got = compose_verdict(arms(), global_triggers=["J pass band", "T1-CLEAN"])
    assert got.startswith("VOID (T1-CLEAN;")


def test_a_trigger_that_voids_D_sign_itself_is_flight_level():
    state = arms()
    state["D-sign"][""] = line(VOID, "(1,2) null")
    assert compose_verdict(state) == "VOID ((1,2) null)"


# --------------------------------------------------------- step 5
def test_the_partition_label_leads_and_the_enumeration_starts_at_W():
    """Round 19: the label carries D-sign, so W is first in the list."""
    state = arms(dsign=FALSE, dmag={"C": line(TRUE)})
    state["W"][""] = line(VOID, "(2,3) null")
    state["A"]["C'"] = line(VOID, "(2,3) null")
    got = compose_verdict(state, partition_label="FALSIFIED")
    assert got == (
        "FALSIFIED; W VOID ((2,3) null); D-mag holds; "
        "A PARTIAL (C line holds; C′ line VOID ((2,3) null))"
    )


def test_precondition_iii_failing_by_a_FAILING_A_C_line_is_anchor_failed():
    state = arms(dsign=FALSE, a={"C": line(FALSE), "C'": line(TRUE)})
    got = compose_verdict(state, partition_label="FALSIFIED")
    assert got.startswith("anchor failed, no split observed")


def test_precondition_iii_failing_by_a_VOID_A_C_line_is_partition_unavailable():
    """Round 18: a void line is neither true nor false and cannot
    satisfy not-A at all, so the two failure modes route differently."""
    state = arms(dsign=FALSE, a={"C": line(VOID, "dose certificate"), "C'": line(TRUE)})
    got = compose_verdict(state, partition_label="FALSIFIED")
    assert got == "NOT DETECTED; partition unavailable (A VOID (dose certificate))"


def test_the_partition_label_is_not_invented_by_the_composer():
    """kappa and the CI decide it; the composer only places it."""
    with pytest.raises(ValueError):
        compose_verdict(arms(dsign=FALSE))


# --------------------------------------------------------- interaction
def test_a_forcing_never_overwrites_a_void_line():
    """Round 18: the s2(C) floor forces not-W, but a void W line has no
    value to force, so it stays VOID."""
    state = arms(dsign=FALSE)
    state["W"][""] = line(VOID, "(2,3) null")
    got = compose_verdict(state, partition_label="FALSIFIED", forced_false=["W"])
    assert "W VOID ((2,3) null)" in got
    assert "¬W" not in got


def test_a_forcing_applies_to_a_non_void_line():
    """Forcing W false makes the pattern D-sign and not-W and A, which
    section 7 names; the forcing has taken effect precisely because
    the line is no longer CONFIRMED."""
    state = arms()
    got = compose_verdict(state, forced_false=["W"])
    assert got == "Split confirmed, class anomalous"


# ------------------------------------- findings of the code review
# Every test below was written from a defect the first reviewer of
# this module MEASURED by running it. They are the reason the
# composer is code: a reader cannot produce these outputs.

def test_step3_appends_triggers_that_fired_with_ARM_scope():
    """Blocker: the step collected only global triggers and D-sign's
    own, so a J-pass-band that voided D-mag and A vanished from the
    line. Section 7: appending, never replacing, keeps the RECORD's
    diagnosis complete."""
    state = {
        "D-sign": {"": line(VOID, "(1,2) null")},
        "W": {"": line(VOID, "(1,2) null")},
        "D-mag": {"C": line(VOID, "J pass band")},
        "A": {"C": line(VOID, "J pass band"), "C'": line(VOID, "J pass band")},
    }
    assert compose_verdict(state) == "VOID ((1,2) null; J pass band)"


def test_an_unregistered_trigger_spelling_is_refused_not_ranked_last():
    """Blocker: unknown strings sank to the lowest rank, so they could
    NAME the line, including the two triggers section 7 proves can
    never lead."""
    state = arms()
    with pytest.raises(ValueError):
        compose_verdict(state, global_triggers=["dose certificat"])


def test_the_registered_unicode_and_spaced_spellings_are_recognised():
    """The document writes s²(C) with a superscript and T1-CLEAN both
    hyphenated and spaced; all are the same registered triggers."""
    state = arms()
    got = compose_verdict(state, global_triggers=["s²(C) floor", "T1 CLEAN"])
    assert got == "VOID (T1 CLEAN; s²(C) floor)"


def test_the_floor_and_the_J_band_can_never_lead_a_void_line():
    state = arms()
    got = compose_verdict(state, global_triggers=["J pass band", "N0-CLEAN"])
    assert got.startswith("VOID (N0-CLEAN;")


@pytest.mark.parametrize("dmag,w,a,expected", [
    (FALSE, TRUE, TRUE, "Split confirmed, width off-prediction"),
    (TRUE, FALSE, TRUE, "Split confirmed, class anomalous"),
    (TRUE, TRUE, FALSE, "Split observed, anchor failed"),
])
def test_the_three_qualified_split_lines_are_emitted(dmag, w, a, expected):
    """Blocker: the code claimed a 'G1-frozen table' owed these names
    and emitted an enumeration instead. Section 7's verdict bullets
    ARE that table, and a pre-data deferral cannot be registered in a
    code comment."""
    state = arms(w=w, dmag={"C": line(dmag)}, a={"C": line(a)})
    assert compose_verdict(state) == expected


def test_the_partition_label_must_be_one_of_the_registered_four():
    with pytest.raises(ValueError):
        compose_verdict(arms(dsign=FALSE), partition_label="FALSIFYED")


def test_confirmed_is_unreachable_with_an_arm_missing():
    """Section 7: CONFIRMED = D-sign AND W AND D-mag AND A."""
    with pytest.raises(ValueError):
        compose_verdict({"D-sign": {"": line(TRUE)}, "W": {"": line(TRUE)},
                         "A": {"C": line(TRUE)}})


def test_a_two_line_D_sign_is_refused():
    """Section 7: D-sign, being single-line, can never be PARTIAL. A
    half-void detection arm had printed as a falsification."""
    state = arms()
    state["D-sign"] = {"C": line(TRUE), "C'": line(VOID, "(1,2) null")}
    with pytest.raises(ValueError):
        compose_verdict(state, partition_label="FALSIFIED")


def test_a_void_line_without_a_trigger_is_refused():
    """It printed 'W VOID (None)'; no registered form has a
    placeholder where the cause belongs."""
    state = arms()
    state["W"] = {"": line(VOID)}
    with pytest.raises(ValueError):
        compose_verdict(state)


def test_the_void_form_names_the_higher_precedence_cause():
    """It named whichever void line came first in the dict, so the
    caller's insertion order decided the printed cause."""
    state = arms(a={"C": line(VOID, "J pass band"),
                    "C'": line(VOID, "(2,3) null")})
    assert "A VOID ((2,3) null)" in compose_verdict(state)


def test_a_forcing_naming_an_unknown_arm_is_refused():
    with pytest.raises(ValueError):
        compose_verdict(arms(), forced_false=["Z"])


def test_the_printed_line_name_is_the_documents_prime_whatever_the_caller_typed():
    """The registered examples spell it C-prime with U+2032. The
    composer echoed the caller's key, so an ASCII apostrophe in the
    runner would have printed a string that is not the registered one
    while every test stayed green."""
    ascii_key = arms(a={"C": line(TRUE), "C'": line(VOID, "(2,3) null")})
    prime_key = arms(a={"C": line(TRUE), "C′": line(VOID, "(2,3) null")})
    expected = "A PARTIAL (C line holds; C′ line VOID ((2,3) null))"
    assert expected in compose_verdict(ascii_key)
    assert expected in compose_verdict(prime_key)


def test_an_unregistered_line_name_is_refused():
    state = arms(a={"C": line(TRUE), "D": line(TRUE)})
    with pytest.raises(ValueError):
        compose_verdict(state)


# ------------------------------------------------- exhaustive walk
# The reason this composer is code and not prose. Eight review rounds
# read the rules and each found the previous round's defect; a state
# machine can simply be walked.
import itertools

VALUES = [TRUE, FALSE, VOID]


def _all_states():
    """Every assignment of line values, in both D-mag regimes."""
    for dmag_lines in (["C"], ["C", "C'"]):
        for ds, w in itertools.product(VALUES, VALUES):
            for dm in itertools.product(VALUES, repeat=len(dmag_lines)):
                for av in itertools.product(VALUES, repeat=2):
                    yield {
                        "D-sign": {"": line(ds, "(1,2) null" if ds == VOID else None)},
                        "W": {"": line(w, "(2,3) null" if w == VOID else None)},
                        "D-mag": {n: line(v, "dose certificate" if v == VOID else None)
                                  for n, v in zip(dmag_lines, dm)},
                        "A": {n: line(v, "dose certificate" if v == VOID else None)
                              for n, v in zip(["C", "C'"], av)},
                    }


def test_every_reachable_state_produces_exactly_one_line():
    """No state may be unclassified, and none may raise except the one
    documented case: step 5 without the partition label the kappa/CI
    evaluator owns."""
    unclassified = []
    for state in _all_states():
        try:
            got = compose_verdict(state, partition_label="FALSIFIED")
        except Exception as e:  # noqa: BLE001 - the point is to catch any
            unclassified.append((state, repr(e)))
            continue
        if not isinstance(got, str) or not got.strip():
            unclassified.append((state, repr(got)))
    assert unclassified == [], f"{len(unclassified)} states without a line"


def test_no_state_prints_a_spelling_outside_the_registered_vocabulary():
    """Section 7: five forms are the whole vocabulary; any other
    spelling in the document is a defect. The same must hold here."""
    import re

    def split_top_level(s):
        """Split on '; ' at parenthesis depth 0. A PARTIAL clause
        carries a ';' INSIDE its parentheses, so a naive split would
        manufacture clauses that were never printed."""
        out, depth, cur = [], 0, ""
        i = 0
        while i < len(s):
            c = s[i]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
            if depth == 0 and s[i:i + 2] == "; ":
                out.append(cur)
                cur = ""
                i += 2
                continue
            cur += c
            i += 1
        out.append(cur)
        return out

    # Closed vocabularies. The earlier version used ".*" for every
    # variable segment, which accepted "W VOID (None)" and
    # "A PARTIAL ( holds;  VOID ())": it asserted the shape of string
    # constants and could not fail on a malformed line (review finding).
    TRIG = r"(?:band-validity window|N0-CLEAN|T1-CLEAN|grid incomplete"
    TRIG += r"|dose certificate|\(1,2\) null|\(2,3\) null|fit health"
    TRIG += r"|s2\(C\) floor|J pass band)"
    LN = r"(?:C|C′)"
    WHOLE_LINE = re.compile(
        r"^(CONFIRMED|Split confirmed, width off-prediction"
        r"|Split confirmed, class anomalous"
        r"|Split observed, anchor failed"
        rf"|VOID \({TRIG}(?:; {TRIG})*\))$")
    LABEL = re.compile(
        r"^(FALSIFIED|Anti-D|INCONCLUSIVE \((underpowered|indeterminate)\)"
        r"|anchor failed, no split observed"
        r"|D-sign detected|¬D-sign)$")
    CLAUSE = re.compile(
        rf"^((W|D-mag|A) holds"
        rf"|¬(W|D-mag|A)( \({LN} line (fails|holds)"
        rf"(; {LN} line (fails|holds))*\))?"
        rf"|(W|D-mag|A) PARTIAL \({LN} line (holds|fails); "
        rf"{LN} line VOID \({TRIG}\)\)"
        rf"|(W|D-mag|A) VOID \({TRIG}\))$")

    bad = []
    seen = set()
    for state in _all_states():
        try:
            got = compose_verdict(state, partition_label="FALSIFIED")
        except ValueError:
            continue
        if WHOLE_LINE.match(got):
            seen.add("whole-line")
            continue
        if re.fullmatch(
                rf"NOT DETECTED; partition unavailable \(A VOID \({TRIG}\)\)",
                got):
            seen.add("partition-unavailable")
            continue
        clauses = split_top_level(got)
        if not LABEL.match(clauses[0]):
            bad.append(got)
            continue
        for cl in clauses[1:]:
            if not CLAUSE.match(cl):
                bad.append(got)
                break
    assert bad == [], f"{len(bad)} lines outside the vocabulary, e.g. {bad[:3]}"
    # a form no state ever produces is itself a finding
    assert {"whole-line", "partition-unavailable"} <= seen, seen
