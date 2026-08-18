"""The corner-beat VERDICT-LINE COMPOSER.

Specification: the composition algorithm registered in section 7 of
`experiments/CORNER_BEAT_HARDWARE_PREDICTION.md` (R=CPsi-squared repo).
This module is that algorithm in executable form, and it exists for a
measured reason: eight consecutive review rounds found their defects
in the previous round's PROSE repairs of these rules and none in the
physics. A state machine that can be walked settles what a paragraph
that is read cannot. The first review of this module made the point
again in the other direction: it found four blockers by RUNNING the
code, none of which nineteen readings of the prose had produced.

The composer does NOT decide verdicts. It decides how a decided state
is WRITTEN. The partition label comes from the kappa and CI
evaluator, a separate machinery-ledger item; passing no label, or an
unregistered one, is an error and never a guess.

Line values are TRUE, FALSE or VOID. A line is a dict
{"value": <one of those>, "trigger": <str or None>}; the trigger is
what voided it and is REQUIRED on a void line. An arm is a dict of
named lines; single-line arms use the empty name. An arm's LINE SET
is data, not an assumption: under the reverted-band freeze the
pre-registration tells the executor to expect, `dmag_band_Cp` holds
the REPORTED sentinel and D-mag is passed with its C line alone.

This module refuses rather than guesses. Every refusal is a state the
pre-registration does not define, and on a paid one-shot a refusal
before submission is cheap where a plausible wrong line is not.
"""
TRUE = "TRUE"
FALSE = "FALSE"
VOID = "VOID"

# Section 7's registered trigger precedence, causal order: the line
# NAMES the first firing trigger and appends the rest in this order.
# Each entry lists every spelling the document uses for that trigger,
# because a spelling the table does not know must be REFUSED and not
# silently ranked last, where it could name the line (review finding:
# "s2(C) floor" never matched the document's own "s²(C) floor").
TRIGGER_ALIASES = [
    ("band-validity window", ("band-validity window", "band validity window")),
    ("N0-CLEAN", ("N0-CLEAN", "N0 CLEAN")),
    ("T1-CLEAN", ("T1-CLEAN", "T1 CLEAN", "T1 CLEAN conditions")),
    ("grid incomplete", ("grid incomplete",)),
    ("dose certificate", ("dose certificate", "dose certificates")),
    ("(1,2) null", ("(1,2) null", "(1,2) difference-null")),
    ("(2,3) null", ("(2,3) null", "(2,3) difference-null")),
    ("fit health", ("fit health",)),
    ("s2(C) floor", ("s2(C) floor", "s²(C) floor")),
    ("J pass band", ("J pass band", "G5 J pass band")),
]

# Amendment 1.6's brake rides at the grid-incomplete rank as a
# SUB-LABEL: it joins with a semicolon when it leads the line and with
# a comma when a higher-ranked trigger leads, which are the two forms
# the Amendment registers.
BRAKE_SUFFIX = "brake-truncated per Amendment 1.6"

ARM_ORDER = ["D-sign", "W", "D-mag", "A"]

# The registered line names. The document spells the second one with a
# PRIME (U+2032); a caller typing an ASCII apostrophe means the same
# line and must not produce a string that is not the registered one,
# so the name is canonicalised on the way out and anything else is
# refused rather than echoed.
LINE_NAME_CANON = {"": "", "C": "C", "C'": "C′", "C′": "C′"}

# Section 7's registered partition labels; the composer places one, it
# never invents or corrects one.
PARTITION_LABELS = (
    "FALSIFIED",
    "Anti-D",
    "INCONCLUSIVE (underpowered)",
    "INCONCLUSIVE (indeterminate)",
)

VALUED_TRUE = "VALUED_TRUE"
VALUED_FALSE = "VALUED_FALSE"
PARTIAL = "PARTIAL"


class VerdictStateError(ValueError):
    """A state the pre-registration does not define."""


def classify_arm(lines):
    """Step 2. VOID if all lines are void; PARTIAL if some but not
    all are; otherwise VALUED, holding iff every line holds.

    PARTIAL is a statement about VOIDNESS ONLY. An arm whose lines all
    valued and one of which FAILS is VALUED_FALSE, never PARTIAL: the
    distinction cost four review rounds, so it is one branch here.
    """
    if not lines:
        raise VerdictStateError("an arm with no lines has no state")
    values = [ln["value"] for ln in lines.values()]
    if all(v == VOID for v in values):
        return VOID
    if any(v == VOID for v in values):
        return PARTIAL
    return VALUED_TRUE if all(v == TRUE for v in values) else VALUED_FALSE


def _canonical_trigger(trigger):
    """The registered trigger this string names, or a refusal."""
    if trigger and trigger.startswith("grid incomplete") and BRAKE_SUFFIX in trigger:
        return "grid incomplete"
    for canon, spellings in TRIGGER_ALIASES:
        if trigger in spellings:
            return canon
    raise VerdictStateError(
        f"unregistered trigger spelling {trigger!r}: a trigger the "
        f"precedence table does not know could name the line")


def _rank(trigger):
    canon = _canonical_trigger(trigger)
    return [c for c, _ in TRIGGER_ALIASES].index(canon)


def _sorted_triggers(triggers):
    return sorted(dict.fromkeys(triggers), key=_rank)


def _void_lines(lines):
    return {n: ln for n, ln in lines.items() if ln["value"] == VOID}


def _leading_void_line(lines):
    """The void line whose trigger ranks highest, so the wider-reaching
    cause names the state rather than the caller's dict order."""
    voids = _void_lines(lines)
    name = min(voids, key=lambda n: _rank(voids[n]["trigger"]))
    return name, voids[name]


def _surviving_line(lines):
    for name, ln in lines.items():
        if ln["value"] != VOID:
            return name, ln
    return None, None


def _qualify(name):
    canon = LINE_NAME_CANON[name]
    return f"{canon} line" if canon else "line"


def _render_arm(arm, lines):
    """Step 4's five print forms, and no others."""
    state = classify_arm(lines)
    if state == VOID:
        _, ln = _leading_void_line(lines)
        return f"{arm} VOID ({ln['trigger']})"
    if state == PARTIAL:
        sname, sln = _surviving_line(lines)
        vname, vln = _leading_void_line(lines)
        verb = "holds" if sln["value"] == TRUE else "fails"
        return (f"{arm} PARTIAL ({_qualify(sname)} {verb}; "
                f"{_qualify(vname)} VOID ({vln['trigger']}))")
    if state == VALUED_TRUE:
        return f"{arm} detected" if arm == "D-sign" else f"{arm} holds"
    if len(lines) > 1:
        failing = [n for n, ln in lines.items() if ln["value"] == FALSE]
        holding = [n for n, ln in lines.items() if ln["value"] == TRUE]
        parts = ([f"{_qualify(n)} fails" for n in failing]
                 + [f"{_qualify(n)} holds" for n in holding])
        return f"¬{arm} ({'; '.join(parts)})"
    return f"¬{arm}"


def _validate(arms, forced_false):
    missing = [a for a in ARM_ORDER if a not in arms]
    if missing:
        raise VerdictStateError(
            f"every verdict arm must be present; missing {missing}. "
            f"CONFIRMED is the conjunction over all four.")
    if len(arms["D-sign"]) != 1:
        raise VerdictStateError(
            "D-sign is single-line by section 7 and can never be PARTIAL")
    for arm, lines in arms.items():
        for name, ln in lines.items():
            if name not in LINE_NAME_CANON:
                raise VerdictStateError(
                    f"{arm}: {name!r} is not a registered line name "
                    f"{tuple(sorted(set(LINE_NAME_CANON.values())))}")
            if ln["value"] not in (TRUE, FALSE, VOID):
                raise VerdictStateError(
                    f"{arm} {name!r}: {ln['value']!r} is not a line value")
            if ln["value"] == VOID:
                if not ln.get("trigger"):
                    raise VerdictStateError(
                        f"{arm} {name!r} is VOID with no trigger; no "
                        f"registered form has a placeholder where the "
                        f"cause belongs")
                _canonical_trigger(ln["trigger"])
    unknown = [a for a in forced_false if a not in arms]
    if unknown:
        raise VerdictStateError(f"forcing names no such arm: {unknown}")


def _apply_forcings(arms, forced_false):
    """A forcing never overwrites a VOID: a void line has no value to
    force. Registered for the s2(C) floor, which forces not-W and
    not-D-sign regardless of the computed exceedance."""
    out = {}
    for arm, lines in arms.items():
        if arm in forced_false:
            out[arm] = {n: (ln if ln["value"] == VOID
                            else {"value": FALSE, "trigger": None})
                        for n, ln in lines.items()}
        else:
            out[arm] = lines
    return out


def _named_line(states):
    """Section 7's four all-valued names. Returns None when no
    registered name covers the pattern, which is itself a finding."""
    holds = {a: states[a] == VALUED_TRUE for a in ARM_ORDER}
    if all(holds.values()):
        return "CONFIRMED"
    if holds["D-sign"] and holds["W"] and not holds["D-mag"] and holds["A"]:
        return "Split confirmed, width off-prediction"
    if holds["D-sign"] and not holds["W"] and holds["A"]:
        return "Split confirmed, class anomalous"
    if holds["D-sign"] and not holds["A"]:
        return "Split observed, anchor failed"
    return None


def _flight_level_void(triggers):
    ordered = _sorted_triggers(triggers)
    parts = []
    for t in ordered:
        # the brake is a sub-label at the grid-incomplete rank: it
        # leads with a semicolon and is appended with a comma
        parts.append(t)
    return f"VOID ({'; '.join(parts)})"


def compose_verdict(arms, global_triggers=(), partition_label=None,
                    forced_false=()):
    """Return the single verdict LINE for one decided state."""
    _validate(arms, forced_false)
    arms = _apply_forcings(arms, forced_false)

    # Step 3: a global trigger, or a trigger that voided D-sign itself,
    # is flight-level. This is the ONE place triggers are appended, and
    # it appends EVERY trigger that fired, including those whose scope
    # was a single arm (review blocker: those were dropped).
    dsign_state = classify_arm(arms["D-sign"])
    if global_triggers or dsign_state == VOID:
        fired = list(global_triggers)
        for lines in arms.values():
            for ln in _void_lines(lines).values():
                fired.append(ln["trigger"])
        return _flight_level_void(fired)

    states = {a: classify_arm(arms[a]) for a in ARM_ORDER}

    # Step 4: D-sign holds.
    if dsign_state == VALUED_TRUE:
        if all(s in (VALUED_TRUE, VALUED_FALSE) for s in states.values()):
            named = _named_line(states)
            if named:
                return named
        return "; ".join(_render_arm(a, arms[a]) for a in ARM_ORDER)

    # Step 5: D-sign is valued false.
    a_lines = arms["A"]
    if "C" not in a_lines:
        raise VerdictStateError(
            "the not-D-sign partition reads A's C line and A has none")
    a_c = a_lines["C"]
    rest = [a for a in ARM_ORDER if a != "D-sign"]
    if a_c["value"] == VOID:
        # (iii) fails by VOIDNESS: a void line cannot satisfy not-A.
        return (f"NOT DETECTED; partition unavailable "
                f"(A VOID ({a_c['trigger']}))")
    if a_c["value"] == FALSE:
        # (iii) fails by FAILURE: the anchor did not stand.
        body = "; ".join(_render_arm(a, arms[a]) for a in rest)
        return f"anchor failed, no split observed; {body}"
    if partition_label is None:
        raise VerdictStateError(
            "step 5 needs the partition label from the kappa/CI "
            "evaluator; the composer never invents one")
    if partition_label not in PARTITION_LABELS:
        raise VerdictStateError(
            f"{partition_label!r} is not one of the registered "
            f"partition labels {PARTITION_LABELS}")
    body = "; ".join(_render_arm(a, arms[a]) for a in rest)
    return f"{partition_label}; {body}" if body else partition_label
