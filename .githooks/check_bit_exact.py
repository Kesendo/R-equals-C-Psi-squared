#!/usr/bin/env python3
"""Warn when a staged (or, with --all, any tracked) file calls something
"bit-exact" on the same line as the tolerance that refutes it.

The repo's own rule is already written down, at docs/GLOSSARY.md:304: *if a
threshold accompanies it, the claim is not bit-exact*. The triage that decides
which routes may use the word at all is at
docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md:260: integer, rational,
BigInteger, GF(p), exact-permutation and symbolic routes may; anything out of an
eigensolver or a time evolution may not.

A hook cannot see which route produced a number, so this checks the narrower
thing it CAN see, and the one the 2026-08-06 census kept finding: the
self-refuting line, where the word sits beside a threshold, a residual, a
percentage or the phrase "machine precision". Every one of those is a sentence
that disproves itself without leaving the line:

    "bit-exact at machine precision"                      (7 sites, fixed 08-06)
    "holds bit-exact (0.000% relative error)"             (a printf width, fixed 08-06)
    "prints 'bit-exact ... max deviation 8.88e-16'"       (PopcountCoherenceClaim:48)
    "the bit-exact balance scale" over a 1e-10 gate       (the three F112 witnesses)

NEGATED USES ARE NOT FLAGGED, and that is deliberate: the repo's model wordings
say the word in order to deny it, quoting the residual right beside it, e.g.
BlockSpectrumWitness's "3.48e-13, i.e. ~174 eps·2σ and NOT bit-exact". Those are
the shape to copy, so the hook must stay quiet on them.

Scope = tracked .cs / .md / .py, which is the census's own denominator.

Used as part of the repo pre-commit hook (warns, never blocks) and runnable by
hand:
    python .githooks/check_bit_exact.py        # check staged changes
    python .githooks/check_bit_exact.py --all  # sweep the whole tree
"""
import re
import subprocess

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from _prose_scope import excluded
import sys

ALL = "--all" in sys.argv

SCOPE = re.compile(r"\.(cs|md|py)$")

# This file is a catalogue of the shape it hunts: its docstring quotes four real
# defects and SELF_TEST below carries five more as fixtures. Without this line the
# hook reports ten hits against itself on every commit that touches it, which is
# how a warn-only check gets ignored. It is the only self-exemption, and it is
# safe because nothing here is a claim about the physics.
SELF = "githooks/check_bit_exact.py"

# The word in its four spellings, plus the noun form.
WORD = re.compile(r"bit[-\s]?(exact(ly|ness)?|identical|for[-\s]bit)", re.I)

# What, on the same line, refutes it. Kept deliberately narrow: each of these is
# a THRESHOLD or a RESIDUAL, not merely a number, so "bit-exact at N=4, 5, 6"
# and "the 4^N basis" do not trip it.
TOLERANCE = re.compile(
    r"""(
        \d\s*(e|E)\s*[-+]\s*\d           # 1e-10, 8.88e-16, 3.48E-13
      | \d\.\d+\s*%                      # 0.000%
      | machine\s+precision
      | machine\s+zero
      | relative\s+(error|deviation)
      | max(imum)?\s+deviation
      | tolerance
      | \bwithin\s+\d
      | \beps\b
    )""",
    re.I | re.X,
)

# A denial in the 40 characters before the word. The window is short on purpose:
# "not" earlier in a long sentence usually negates something else. Quotes and
# emphasis markers may sit between the denial and the word, and DO: the model
# wording in StarImMaxBoundClaim is `says "machine precision" and not "bit-exact"`,
# which a bare \s*$ anchor misses and then reports as a defect.
# NOT re.X here, deliberately. Under re.X the literal spaces inside "no longer",
# "rather than" and "instead of" are stripped, so the engine sees "nolonger" and
# those three denials silently never match. The first version of this file had
# that bug, and it hit exactly the sentences the exemption exists for: the model
# wording in F113_T1_EXTRACTION_KINGSTON.md is "machine precision RATHER THAN
# bit-identical". The self-test at the bottom of this file covers all of them.
# An ARTICLE may stand between the denial and the word, and the plainest denial in
# English has one: "is not A bit-exact scale". The first version anchored the
# negation immediately against the phrase, so every "not a" / "not the" sentence was
# reported as a defect. Found 2026-08-06 by the hook firing on the commit that was
# rewriting exactly such a line. The determiner is optional and eats no more than one
# word, so "not the place; bit-exact holds" is still reported (the semicolon breaks it).
NEGATION = re.compile(
    r"""(not|n't|never|no longer|rather than|instead of|isn't|aren't)"""
    r"""[\s"'`*_(\[]*(a|an|the)?[\s"'`*_(\[]*$""",
    re.I,
)


def flags(line):
    """True if this line spends the word beside the threshold that refutes it."""
    m = WORD.search(line)
    if not m or not TOLERANCE.search(line):
        return False
    before = line[max(0, m.start() - 40):m.start()].rstrip()
    return not NEGATION.search(before)


# The cases the regexes must get right, runnable with --self-test. Every "quiet"
# row below is a real sentence from the repo that the hook must NOT report: they
# are the model wordings the campaign tells people to copy, so a hook that flags
# them trains the opposite of what it is for. Two of them were flagged by earlier
# drafts of this file, which is why they are pinned here rather than trusted.
SELF_TEST = [
    (False, 'pairing distance of 3.48e-13, i.e. ~174 eps*2sigma and NOT bit-exact,'),
    (False, 'why this claim says "machine precision" and not "bit-exact".'),
    (False, 'within 4e-15 of 1, i.e. machine precision rather than bit-identical'),
    (False, 'the identity is no longer bit-exact once the tolerance is 1e-9'),
    (False, 'reports the residual instead of bit-exactness (max deviation 8.9e-16)'),
    (False, 'the residual is *not* **bit-exact** (1e-12)'),
    (False, 'bit-exact at N=4, 5, 6 over the 4^N Pauli basis'),
    # The article gap, and its boundary: a denial may carry one determiner, no more.
    (False, 'Default 1e-10 is not a bit-exact scale (the break is at 7.7e-3)'),
    (False, 'this is not an bit-exact route, the gate is 1e-9'),
    (False, 'reports machine precision rather than the bit-exact claim (1e-12)'),
    (True, 'not the place to say it; bit-exact holds to 1e-9 here'),
    (True, "/// Default 1e-10 matches the bit-exact balance scale of F112-X's"),
    (True, "prints 'bit-exact ... max deviation 8.88e-16' in one summary"),
    (True, 'the Absorption Theorem holds bit-exact (0.000% relative error)'),
    (True, 'Verified bit-exactly (1e-9) against the Build zgeev path'),
    (True, 'bit-identical eigenvalue multiset (binned at 1e-7)'),
]

if "--self-test" in sys.argv:
    bad = [(w, s) for w, s in SELF_TEST if flags(s) != w]
    for want, s in bad:
        print(f"FAIL want={want} got={not want}: {s}")
    print(f"{len(SELF_TEST) - len(bad)}/{len(SELF_TEST)} pass")
    sys.exit(1 if bad else 0)


def git(*a):
    return subprocess.run(["git", *a], capture_output=True, text=True,
                          encoding="utf-8", errors="replace").stdout


if ALL:
    files = git("ls-files").splitlines()

    def content(f):
        try:
            return open(f, encoding="utf-8", errors="replace").read()
        except OSError:
            return ""
else:
    files = git("diff", "--cached", "--name-only", "--diff-filter=ACM").splitlines()

    def content(f):
        return git("show", f":{f}")  # the staged blob


hits = []  # (file, line_no, excerpt)
for f in files:
    if not SCOPE.search(f) or f.replace("\\", "/").endswith(SELF) or excluded(f):
        continue
    for i, line in enumerate(content(f).splitlines(), 1):
        if not flags(line):
            continue
        excerpt = line.strip()
        if len(excerpt) > 110:
            excerpt = excerpt[:107] + "..."
        hits.append((f, i, excerpt))

if hits:
    w = sys.stderr.write
    w("\n  \033[33mbit-exact warning\033[0m: these lines call something bit-exact "
      "and quote the threshold that refutes it.\n"
      "  The repo's rule is docs/GLOSSARY.md:304 -- if a threshold accompanies it, "
      "the claim is not bit-exact.\n")
    for f, i, excerpt in hits[:20]:
        w(f"    {f}:{i}\n        {excerpt}\n")
    if len(hits) > 20:
        w(f"    ... +{len(hits) - 20} more\n")
    w("  Say the measured number in its own units instead, and name why it is not "
      "bit-exact\n"
      "  (models: BlockSpectrumWitness, StarImMaxBoundClaim's WorstAnchorRelativeDeviation).\n"
      "  (warning only -- commit not blocked)\n\n")

sys.exit(0)  # never block
