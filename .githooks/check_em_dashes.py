#!/usr/bin/env python3
"""Warn if a staged (or, with --all, any tracked) repo doc contains an em dash
(U+2014, '—').

The repo convention (CLAUDE.md, "No em dashes") forbids the em dash in all
documentation: replace it context-dependently with a comma, colon, or semicolon
(en dashes, U+2013, are allowed for ranges, and the minus sign U+2212 is fine).

Beyond punctuation hygiene, an em dash is a reliable *fingerprint*: the
no-em-dash convention is applied during the full quality round (readability,
glosses, links, currency, punctuation), so an em dash surviving in a doc means
that round was never done -- even if a math/scope/honesty gate was. See the
2026-06-25 row in docs/superpowers/REVIEW_LOG.md (the breadcrumb Tom asked to
leave) and the em-dash re-read campaign it seeded.

Scope = the doc dirs named in that breadcrumb (docs/ experiments/ hypotheses/
reflections/ recovered/). The private review apparatus (docs/superpowers/,
docs/CAUGHT_ERRORS.md) is gitignored and so never appears here.

Used as the repo pre-commit hook (warns, never blocks) and runnable by hand:
    python .githooks/check_em_dashes.py        # check staged changes
    python .githooks/check_em_dashes.py --all  # sweep the whole tree
"""
import re
import subprocess
import sys

ALL = "--all" in sys.argv
EM_DASH = "—"
# .md under the five documentation dirs the convention/backlog covers
SCOPE = re.compile(r"^(docs|experiments|hypotheses|reflections|recovered)/.*\.md$")


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

hits = []  # (file, line_no, count)
for f in files:
    if not SCOPE.match(f):
        continue
    for i, line in enumerate(content(f).splitlines(), 1):
        n = line.count(EM_DASH)
        if n:
            hits.append((f, i, n))

if hits:
    by_file = {}
    for f, i, n in hits:
        by_file.setdefault(f, []).append((i, n))
    w = sys.stderr.write
    w(f"\n  \033[33mem-dash warning\033[0m: these docs contain the em dash "
      f"(U+2014), forbidden by the no-em-dash convention.\n"
      f"  An em dash also flags a doc whose full quality round was never done "
      f"(see REVIEW_LOG 2026-06-25).\n")
    for f in sorted(by_file):
        total = sum(n for _, n in by_file[f])
        lines = ", ".join(str(i) for i, _ in by_file[f][:12])
        more = "" if len(by_file[f]) <= 12 else f", +{len(by_file[f]) - 12} more"
        w(f"    {f}  ({total} on lines {lines}{more})\n")
    w("  Replace with a comma/colon/semicolon (en dash U+2013 ok for ranges). "
      "(warning only -- commit not blocked)\n\n")

sys.exit(0)  # never block
