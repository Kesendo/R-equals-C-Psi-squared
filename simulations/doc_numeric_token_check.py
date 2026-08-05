"""
Document numeric-token checker: does every number a document publishes still
appear in the run it was published from?

The question this answers. A script's output moves (a convention is repaired, a
parameter changes). The document quoting that output has to be re-synced. The
failure mode this tool exists to stop is sweeping only for the numbers you
already knew had moved, which can only confirm what you already knew: in one
session that left four whole tables stale, including a peak printed below two
trajectory values above it, and a follow-up check missed a row because its regex
matched only decimals and was structurally blind to the integers "119/60".

So: EVERY numeric token in the document, integers included, is checked against
the run. A token matches if the run prints a number that agrees with it at the
token's OWN printed precision, so "2.54" matches a run value of 2.5372 and
"0.0242" does not match 0.0121.

Usage:
    python simulations/doc_numeric_token_check.py <doc.md> <new_run.txt> \
        --baseline <old_run.txt> [--context]
    python simulations/doc_numeric_token_check.py <doc.md> <run.txt>

USE THE --baseline FORM DURING A REPAIR. It splits the unmatched tokens into
WENT STALE (agreed with the old run, does not agree with the new one: the
re-sync worklist) and NEVER FROM THIS RUN (unmatched under both, so imported
from another document, derived by the prose, or an address the tool mistook for
a measurement). Without a baseline, every derived number, commit hash, date and
section number in the document lands in one undifferentiated unmatched list and
the signal is buried. The baseline form is also what contains the collision
weakness below, since a spurious match has to survive both runs.

WHAT THIS TOOL CAN GET WRONG, which is the part to read before believing it.

1. It matches against ANY number the run printed, not the number at that
   position, so a token can collide with an unrelated run value. Measured
   against a real 4489-value run file, an invented one-decimal number in 0..10
   is accepted 63% of the time, an invented 2-digit integer 18%, an invented
   4-decimal number in 0..1 0.3%. Low-precision tokens are therefore close to
   uncheckable one at a time; precision is what makes this tool work.
2. It cannot see numbers written as words, inside a formula, as a ratio the
   prose computed itself from two printed values, or in notations the regex does
   not cover. The report ends with a scan for the notations known to be
   invisible (vulgar fractions, superscript exponents, leading-dot decimals) and
   states how many it found, so their absence is asserted rather than assumed.
3. A match is not a reading. A token can agree with the run and still be
   attached to the wrong claim in the sentence around it.

Deliberately NOT done: filtering out years, dates and small integers. An earlier
version skipped those to reduce noise and thereby stopped checking 73% of the
tokens in a real document, including "2048" (an N=11 dimension in this
repository, not a year) and every single-digit count the repo publishes
constantly ("zero exceptions", "3, 6, 11"). Noise in the unmatched list is the
price of checking everything, and --baseline mode is where that price is paid.

Exit code 0 if every token matched, 1 if any did not. In --baseline mode the
exit code is 1 whenever anything is unmatched, including a report whose entire
content is the expected NEVER FROM THIS RUN bucket, so read the WENT STALE
count rather than the exit code during a repair.
"""

import re
import sys
import unicodedata
from pathlib import Path

# The documents are Unicode throughout (gamma, mu, arrows), and the Windows
# console default is cp1252, which raises on the first one printed as context.
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# Thousands separators are only accepted in full 3-digit groups, so the comma in
# "March 30, 2026" is not swallowed into the number.
#
# The sign class carries U+2212 as well as ASCII "-". This repository's
# convention is Unicode throughout, so the minus a document actually contains is
# usually U+2212; a sign class of [+-] alone silently reads "−0.5" as +0.5 and
# calls a flipped sign clean. That is the unsafe spelling being the prescribed
# one, so it is the spelling that must work.
#
# The trailing lookahead is (?!\d) and NOT (?![\w]). This repository writes its
# numbers with units and suffixes attached: 2.86x, 22.4us, 5.0%. A lookahead
# that refuses a following letter skips every one of them, on the document side
# AND on the run side, and reports a confident clean. That is not hypothetical:
# an earlier version did exactly that, and silently passed a "2.54x" whose run
# had moved to 2.53x. The cost is that digit-led addresses (commit hash 22481a3,
# a "15:20" timestamp, "Section 2.3.1") also yield tokens; in --baseline mode
# those match neither run and land in NEVER FROM THIS RUN, out of the way.
SIGN = r'[+\-−]?'
NUM_RE = re.compile(
    rf'(?<![\w.]){SIGN}\d{{1,3}}(?:,\d{{3}})+(?:\.\d+)?(?!\d)'
    rf'|(?<![\w.]){SIGN}\d+(?:\.\d+)?(?:[eE][+\-−]?\d+)?(?!\d)'
)

# Notations the regex cannot represent. Counted and reported, never silently
# dropped: the point of the report is that "I saw nothing" and "there was
# nothing" are different statements.
INVISIBLE_NOTATIONS = [
    ('vulgar fraction (¼, ½, ...)', re.compile(r'[¼-¾⅐-⅞]')),
    ('superscript exponent (10⁻³)', re.compile(r'[⁰-₟]')),
    ('leading-dot decimal (.5)', re.compile(r'(?<![\w.\d])\.\d+')),
    ('dotted version/section (1.5.3)', re.compile(r'\d+\.\d+\.\d+')),
]


def normalize_sign(tok):
    """Map U+2212 to ASCII minus so float() will take the token."""
    return tok.replace('−', '-')


def parse_number(tok):
    """Return the float value of a token, or None if it will not parse."""
    try:
        return float(normalize_sign(tok).replace(',', ''))
    except ValueError:
        return None


def decimals_of(tok):
    """
    The token's own printed precision, in decimal places, EXPONENT INCLUDED.

    Reading only the mantissa is a false-clean route, not a rounding nicety:
    "1.5e-3" read as 1 decimal gets a half-band of 0.05 around a value of
    0.0015, so it matches almost anything the run prints, including the 0.0 that
    nearly every run prints somewhere. Every Ne-k token in the repository would
    be unfalsifiable.
    """
    t = normalize_sign(tok)
    exp = 0
    if 'e' in t or 'E' in t:
        body, _, e = re.split(r'([eE])', t, maxsplit=1)[0], None, None
        m = re.search(r'[eE]([+-]?\d+)', t)
        exp = int(m.group(1)) if m else 0
    else:
        body = t
    mantissa_decimals = len(body.split('.')[1]) if '.' in body else 0
    return mantissa_decimals - exp


def collect_run_values(paths):
    """Every number the run printed, as floats."""
    values = []
    for p in paths:
        text = Path(p).read_text(encoding='utf-8', errors='replace')
        for m in NUM_RE.finditer(text):
            v = parse_number(m.group())
            if v is not None:
                values.append(v)
    return values


def matches_run(value, decimals, run_values):
    """
    True if the run printed a number that agrees with this token at the token's
    own printed precision. Rounding both sides is what makes "2.54" accept a run
    value of 2.5372 without accepting 2.53.
    """
    target = round(value, decimals)
    for rv in run_values:
        if round(rv, decimals) == target:
            return True
    # The round() test above misses the half-way cases, where the document
    # rounded 16.95 up to "17.0" and Python's round() takes it down to 16.9. The
    # band test catches those, and its boundary must be inclusive for exactly
    # that reason: a strict "<" rejects the half-way case it exists to accept.
    # (The epsilon is for binary representation, not for slack: 16.95 is not
    # exactly representable, so the difference lands a few ulps either side.)
    half = 10.0 ** (-decimals) / 2
    for rv in run_values:
        if abs(rv - value) <= half + 1e-12:
            return True
    return False


def scan(doc_lines, run_values):
    """Return (checked_count, unmatched list, unparsable list) for one run."""
    checked = 0
    unmatched = []
    unparsable = []
    for lineno, line in enumerate(doc_lines, start=1):
        for m in NUM_RE.finditer(line):
            tok = m.group()
            value = parse_number(tok)
            if value is None:
                unparsable.append((lineno, tok))
                continue
            checked += 1
            if not matches_run(value, decimals_of(tok), run_values):
                unmatched.append((lineno, tok, line.strip()))
    return checked, unmatched, unparsable


def report_invisible(doc_lines):
    print("NOT SEEN BY THE REGEX (counted so their absence is a statement, not an")
    print("assumption; each is a number this tool did not check at all)")
    print("-" * 70)
    any_found = False
    for name, rx in INVISIBLE_NOTATIONS:
        hits = [(i, m.group()) for i, line in enumerate(doc_lines, start=1)
                for m in rx.finditer(line)]
        if hits:
            any_found = True
            sample = ', '.join(f"line {i}: {t}" for i, t in hits[:5])
            more = '' if len(hits) <= 5 else f", ... (+{len(hits) - 5})"
            print(f"  {len(hits):>5}  {name}: {sample}{more}")
    if not any_found:
        print("  none of the known-invisible notations occur in this document")
    print()
    print("Still outside reach and not counted: numbers stated in words, numbers")
    print("inside a formula, and ratios the prose computes from two printed values.")


def _print_block(title, items, show_context):
    print(title)
    print("-" * 70)
    for lineno, tok, context in items:
        print(f"  line {lineno:>5}  {tok}")
        if show_context:
            print(f"              | {context[:110]}")
    if not items:
        print("  none")
    print()


def check(doc_path, run_paths, show_context=False, baseline_paths=None):
    doc_lines = Path(doc_path).read_text(
        encoding='utf-8', errors='replace').splitlines()
    run_values = collect_run_values(run_paths)
    checked, unmatched, unparsable = scan(doc_lines, run_values)

    print(f"document : {doc_path}")
    for p in run_paths:
        print(f"run      : {p}")
    if baseline_paths:
        for p in baseline_paths:
            print(f"baseline : {p}")
    print(f"run values collected: {len(run_values)}")
    print(f"tokens checked      : {checked}")
    print(f"tokens unmatched    : {len(unmatched)}")
    print()

    if baseline_paths:
        base_values = collect_run_values(baseline_paths)
        _, base_unmatched, _ = scan(doc_lines, base_values)
        base_keys = {(ln, tok) for ln, tok, _ in base_unmatched}
        went_stale = [u for u in unmatched if (u[0], u[1]) not in base_keys]
        never = [u for u in unmatched if (u[0], u[1]) in base_keys]
        _print_block(
            f"WENT STALE ({len(went_stale)}): matched the baseline run, does not "
            f"match the new one", went_stale, show_context)
        _print_block(
            f"NEVER FROM THIS RUN ({len(never)}): unmatched under BOTH runs\n"
            "(imported from another document, derived in the prose, or an address\n"
            " the regex read as a measurement -- this tool cannot tell which)",
            never, show_context)
    elif unmatched:
        _print_block(
            "UNMATCHED (no number in the run agrees at the token's own precision)",
            unmatched, show_context)

    if unparsable:
        _print_block(f"UNPARSABLE ({len(unparsable)})",
                     [(ln, tok, '') for ln, tok in unparsable], False)

    report_invisible(doc_lines)
    return 0 if not unmatched else 1


def main(argv):
    show_context = '--context' in argv
    rest = [a for a in argv[1:] if a != '--context']

    baseline = []
    if '--baseline' in rest:
        i = rest.index('--baseline')
        baseline = rest[i + 1:]
        rest = rest[:i]
        if not baseline:
            print("--baseline was given with no file. Refusing to fall back to "
                  "plain mode silently: that is the mode you did not ask for.")
            return 2

    if len(rest) < 2:
        print(__doc__)
        return 2
    doc, runs = rest[0], rest[1:]
    for p in [doc] + runs + baseline:
        if not Path(p).exists():
            print(f"missing file: {p}")
            return 2
    return check(doc, runs, show_context, baseline or None)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
