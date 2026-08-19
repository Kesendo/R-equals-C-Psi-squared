#!/usr/bin/env python3
"""Warn if a staged (or, with --all, any tracked) file references a gitignored
local-only WIP script `simulations/_<name>.py`.

Such a reference is dangling by construction: `simulations/_*.py` is local-only
(see .gitignore, "WIP exploration scripts"). Fix by promoting the script
(`git mv simulations/_foo.py simulations/foo.py`; commit it) or dropping the
link. Load-bearing verifier evidence should become a C# witness, not a script.

Used as the repo pre-commit hook (warns, never blocks) and runnable by hand:
    python .githooks/check_wip_refs.py        # check staged changes
    python .githooks/check_wip_refs.py --all  # sweep the whole tree
"""
import re
import subprocess

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from _prose_scope import excluded
import sys

ALL = "--all" in sys.argv
# A reference to a local-only WIP file, in the three shapes docs actually write.
#
# The original pattern demanded a CONCRETE `_<name>.py` and so had a blind spot
# that swallowed the worst cases, found 2026-08-06: a doc that cites a GLOB
# ("simulations/_polarity_probe_*.py", "simulations/_f89_edgeA_*") or a non-.py
# WIP file ("simulations/_atmosphere_cluster_notes.md") was invisible to it. The
# --all sweep therefore reported CLEAN while five tracked files, three of them in
# docs/proofs/ and one the formula registry, named glob-shaped WIP scripts as
# their stated empirical anchor. Those scripts are not in the repo and, checked
# that day, no longer on the machine either.
#
# Three alternatives, in order: a concrete _name.py; a glob _name_*(.py); a
# non-.py WIP file. The trailing (?![*\w]) on the first keeps a glob from being
# reported twice.
PATH = re.compile(
    r"simulations/_[A-Za-z0-9][A-Za-z0-9_]*\.py(?![*\w])"
    r"|simulations/_[A-Za-z0-9][A-Za-z0-9_]*\*(?:\.py)?"
    r"|simulations/_[A-Za-z0-9][A-Za-z0-9_]*\.(?:md|txt|json|ipynb|cmd)"
)
SKIP_EXT = (".png", ".jpg", ".jpeg", ".gif", ".pdf", ".csv", ".npy", ".npz",
            ".zip", ".dll", ".pyc", ".bin", ".ico", ".svg")
# exclude files that legitimately spell the pattern (the rule + this checker)
SKIP_PATHS = {".gitignore"}


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

# TWO CLASSES, and only one of them may ever be excluded (2026-08-19).
#
# A reference to `simulations/_foo.py` whose PROMOTED twin `simulations/foo.py`
# is tracked is not a local-only anchor at all: it is a citation left behind by
# a rename, pointing at nothing while the thing it names sits in the repo under
# its real name. That is a broken link for today's reader, so it is reported
# even in files the prose conventions exclude. The other class, an anchor with
# no promoted twin, is provenance: the script really is local-only, and a file
# whose preamble says so may be excluded from the nagging.
#
# The distinction was forced by a live case: docs/CAUGHT_ERRORS.md left
# .gitignore on 2026-08-19, was excluded from the prose hooks the same morning,
# and five stale-rename citations went silent with it. A reviewer then read the
# silence as a blind spot in this checker. The checker was fine; the exclusion
# was too blunt. Link integrity is not a style rule and does not become
# historical: a dead link is dead for the reader today.
tracked = set(git("ls-files").splitlines())

stale, anchors = [], []
for f in files:
    if f in SKIP_PATHS or f.startswith(".githooks/"):
        continue
    if f.lower().endswith(SKIP_EXT):
        continue
    for i, line in enumerate(content(f).splitlines(), 1):
        for m in set(PATH.findall(line)):
            twin = m.replace("simulations/_", "simulations/", 1)
            if twin in tracked:
                stale.append((f, i, m, twin))
            elif not excluded(f):
                anchors.append((f, i, m))

w = sys.stderr.write
if stale:
    w("\n  \033[33mstale-rename warning\033[0m: these cite a `_`-path whose "
      "PROMOTED twin is in the repo.\n  The script was renamed and the citation "
      "was not, so the link points at nothing.\n")
    for f, i, m, twin in stale:
        w(f"    {f}:{i}\n        cites {m}\n        but   {twin} is tracked\n")
    w("  Fix the citation. (warning only -- commit not blocked)\n\n")

if anchors:
    w("\n  \033[33mWIP-ref warning\033[0m: these reference gitignored, local-only "
      "scripts (simulations/_*.py):\n")
    for f, i, m in anchors:
        w(f"    {f}:{i}  ->  {m}\n")
    w("  Promote it (git mv simulations/_foo.py simulations/foo.py; commit it), "
      "drop the link,\n  or port the evidence to a C# witness. "
      "(warning only -- commit not blocked)\n\n")

sys.exit(0)  # never block
