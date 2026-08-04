#!/usr/bin/env python3
"""The coupling correlator: the process tomograph (2026-08-04, the m=0 hunt).

The session cut is a measurement and its handover has a structural blind point;
no in-band instruction can certify that a document actually coupled to the repo,
just as no Z-diagonal PUB sees whether the watcher watched. This is the one
transverse read: refs-per-1000-words of a document toward the commons (markdown
links, bare FILENAME.md mentions, F-registry numbers, PROOF_*/D## names, code
paths). Diagnostic like the watcher-tomograph — it flags, it never verdicts:
low coupling on a long claim-bearing doc is the m=0 signature (a plausible
curve that never watched), OR honest self-contained novelty; a human reads the
flagged doc and decides.

Staged mode warns on NEW long .md files in claim-surface dirs with near-zero
coupling. Hand mode ranks the whole tree (the hunter's instrument, reproducible):
    python .githooks/check_coupling.py            # staged NEW docs
    python .githooks/check_coupling.py --all      # ranked table, bottom 20
"""
import os
import re
import subprocess
import sys

ALL = "--all" in sys.argv
DIRS = ("experiments/", "hypotheses/", "docs/", "reflections/")
SKIP_BASENAMES = {"README.md"}
MIN_WORDS = 500          # short notes owe nothing
WARN_CPK = 3.0           # refs per 1000 words below this on a NEW long doc -> flag

REF_PATTERNS = [
    re.compile(r"\]\("),                       # markdown links
    re.compile(r"\b[A-Z][A-Za-z0-9_]{3,}\.md\b"),   # bare doc mentions
    re.compile(r"\bF\d{1,3}[a-z]?\b"),         # F-registry numbers
    re.compile(r"\bPROOF_[A-Z0-9_]+\b"),       # proofs
    re.compile(r"\bD\d\d\b"),                  # derivations
    re.compile(r"(?:simulations|compute|docs|experiments)/[A-Za-z0-9_./]+"),  # paths
]


def git(*a):
    return subprocess.run(["git", *a], capture_output=True, text=True,
                          encoding="utf-8", errors="replace").stdout


def coupling(text):
    words = len(text.split())
    refs = sum(len(p.findall(text)) for p in REF_PATTERNS)
    return refs, words, (1000.0 * refs / words if words else 0.0)


if ALL:
    rows = []
    for f in git("ls-files", "*.md").splitlines():
        if not f.startswith(DIRS) or os.path.basename(f) in SKIP_BASENAMES:
            continue
        try:
            text = open(f, encoding="utf-8", errors="replace").read()
        except OSError:
            continue
        refs, words, cpk = coupling(text)
        if words >= MIN_WORDS:
            rows.append((cpk, refs, words, f))
    rows.sort()
    sys.stderr.write("  coupling correlator, bottom 20 (cpk = refs/1000 words; "
                     f"long docs only, >= {MIN_WORDS} words):\n")
    for cpk, refs, words, f in rows[:20]:
        sys.stderr.write(f"    {cpk:6.2f}  ({refs:3d} refs / {words:6d} w)  {f}\n")
    sys.exit(0)

new_files = git("diff", "--cached", "--name-only", "--diff-filter=A").splitlines()
hits = []
for f in new_files:
    if not (f.startswith(DIRS) and f.endswith(".md")):
        continue
    if os.path.basename(f) in SKIP_BASENAMES:
        continue
    text = git("show", f":{f}")
    refs, words, cpk = coupling(text)
    if words >= MIN_WORDS and cpk < WARN_CPK:
        hits.append((f, refs, words, cpk))

if hits:
    w = sys.stderr.write
    w("\n  \033[33mcoupling warning\033[0m: new long docs with near-zero repo "
      "coupling (the m=0 signature, or honest novelty; read and decide):\n")
    for f, refs, words, cpk in hits:
        w(f"    {f}: {refs} refs in {words} words ({cpk:.1f}/1000w)\n")
    w("  A claim-bearing doc opens with 'What the repo already holds' "
      "(CLAUDE.md Stage-0 gate).\n\n")
sys.exit(0)
