#!/usr/bin/env python3
"""One owner for "which files do the prose hooks judge, and which do they not".

The prose hooks (em dashes, broken tables, spent words, bit-exact, WIP refs)
used to agree on their scope by each carrying its own copy of the same regex
and by relying on .gitignore to keep the review apparatus out. That second
half broke on 2026-08-19, when `docs/CAUGHT_ERRORS.md` was deliberately taken
out of .gitignore so that the store CLAUDE.md's Stage-0 rule names could
actually be read from a clone. The moment it became a tracked file it started
tripping four hooks at once, ~400 warnings, none of them actionable.

`check_em_dashes.py` had said in its own docstring that the review apparatus
"is gitignored and so never appears here". That was the exclusion: incidental,
undeclared, and invisible until it stopped holding. This module makes it
explicit, so that the next file that leaves .gitignore does not re-run the
same surprise.

A hook that always warns is a hook people learn to scroll past, which costs
more than the convention it enforces.
"""
import re

# The five documentation directories the review backlog covers.
DOC_DIRS = ("docs", "experiments", "hypotheses", "reflections", "recovered")
SCOPE = re.compile(r"^(" + "|".join(DOC_DIRS) + r")/.*\.md$")

# Tracked files the prose conventions deliberately do NOT judge, each with the
# reason, because an exclusion without a reason is indistinguishable from a
# file someone quietly got tired of fixing.
EXCLUDED = {
    "docs/CAUGHT_ERRORS.md":
        "an append-only ledger of dated minutes, not prose under revision. Its "
        "entries quote other documents verbatim (em dashes and all), its "
        "anchors are deliberately provenance rather than links (many point at "
        "gitignored scratch scripts, which is stated in the file's own "
        "preamble), and it is never given the full quality round whose absence "
        "the em-dash check is a fingerprint for. Judging it would rewrite the "
        "record it exists to keep. NEW entries still follow the conventions; "
        "that is on the writer, not on a hook that cannot tell an entry "
        "written today from one written in June.",
}


def excluded(path):
    """The exclusion reason for `path`, or None. Accepts either separator."""
    return EXCLUDED.get(path.replace("\\", "/"))


def in_prose_scope(path):
    """True if `path` is a repo doc the prose hooks judge."""
    p = path.replace("\\", "/")
    return bool(SCOPE.match(p)) and p not in EXCLUDED
