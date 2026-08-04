#!/usr/bin/env python3
"""Warn when a staged arc document is cited by tracked files it never cites back.

The missing BACK-EDGE (the 2026-08-04 m=0 hunt's structural find): an arc opener
gets edited long after its successors landed, still listing as open what they
closed, because nothing in the workflow ever walks back to the opener. The repo
held ~300 such one-directional pairs at the sweep. The successors are scrupulous;
the opener cannot feel them. So the trigger lives in the structure: when you
stage an edit to a doc in experiments/, hypotheses/, or reflections/, this warns
with the list of tracked files that cite it and are not cited back, and asks the
one question: do they close anything this file still calls open?

Non-blocking, like every hook here. Runnable by hand:
    python .githooks/check_back_edges.py            # staged docs
    python .githooks/check_back_edges.py --all      # sweep the whole tree
"""
import os
import re
import subprocess
import sys

ALL = "--all" in sys.argv
# arc documents owe back-edges; registries/indexes are MEANT to be cited
# one-way, so only these directories are checked as EDITED files:
ARC_DIRS = ("experiments/", "hypotheses/", "reflections/")
SKIP_BASENAMES = {"README.md"}
MAX_LIST = 10


def git(*a):
    return subprocess.run(["git", *a], capture_output=True, text=True,
                          encoding="utf-8", errors="replace").stdout


tracked_md = [f for f in git("ls-files", "*.md").splitlines() if f]

if ALL:
    targets = [f for f in tracked_md
               if f.startswith(ARC_DIRS) and os.path.basename(f) not in SKIP_BASENAMES]
else:
    staged = git("diff", "--cached", "--name-only", "--diff-filter=M").splitlines()
    targets = [f for f in staged
               if f.startswith(ARC_DIRS) and f.endswith(".md")
               and os.path.basename(f) not in SKIP_BASENAMES]

if not targets:
    sys.exit(0)

# one read of every tracked md (basename mentions are the citation currency here;
# markdown links, backticked paths, and bare filenames all carry the basename)
contents = {}
for f in tracked_md:
    try:
        contents[f] = open(f, encoding="utf-8", errors="replace").read()
    except OSError:
        contents[f] = ""

warned = False
for t in targets:
    base = os.path.basename(t)
    body = contents.get(t, "")
    citers = []
    for f, c in contents.items():
        if f == t or os.path.basename(f) in SKIP_BASENAMES:
            continue
        if base in c and os.path.basename(f) not in body:
            citers.append(f)
    if citers:
        if not warned:
            sys.stderr.write(
                "\n  \033[33mback-edge warning\033[0m: staged arc docs are cited by "
                "files they never cite back.\n  Before finishing, ask the one "
                "question: do any of these CLOSE something the staged doc still "
                "calls open?\n")
            warned = True
        sys.stderr.write(f"    {t}  <-  {len(citers)} un-returned citer(s):\n")
        for c in sorted(citers)[:MAX_LIST]:
            sys.stderr.write(f"        {c}\n")
        if len(citers) > MAX_LIST:
            sys.stderr.write(f"        ... and {len(citers) - MAX_LIST} more\n")

if warned:
    sys.stderr.write("  (Non-blocking. A back-edge is one line: cite the successor "
                     "and strike the stale open item.)\n\n")
sys.exit(0)
