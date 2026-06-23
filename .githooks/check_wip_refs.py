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
import sys

ALL = "--all" in sys.argv
# a concrete _<name>.py (won't match the glob `_*.py` in .gitignore itself)
PATH = re.compile(r"simulations/_[A-Za-z0-9][A-Za-z0-9_]*\.py")
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

hits = []
for f in files:
    if f in SKIP_PATHS or f.startswith(".githooks/"):
        continue
    if f.lower().endswith(SKIP_EXT):
        continue
    for i, line in enumerate(content(f).splitlines(), 1):
        for m in set(PATH.findall(line)):
            hits.append((f, i, m))

if hits:
    w = sys.stderr.write
    w("\n  \033[33mWIP-ref warning\033[0m: these reference gitignored, local-only "
      "scripts (simulations/_*.py):\n")
    for f, i, m in hits:
        w(f"    {f}:{i}  ->  {m}\n")
    w("  Promote it (git mv simulations/_foo.py simulations/foo.py; commit it), "
      "drop the link,\n  or port the evidence to a C# witness. "
      "(warning only -- commit not blocked)\n\n")

sys.exit(0)  # never block
