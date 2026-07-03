#!/usr/bin/env python3
"""Warn if a staged (or, with --all, any tracked) repo doc contains a broken
GFM (Markdown) table.

Like the em dash (see check_em_dashes.py), a broken table is a reliable
*fingerprint*: a table that does not render was never actually *looked at*, so
the doc's full quality round was never done. The tell is orthogonal to the em
dash: a doc can be em-dash-clean and still carry a table that renders as a wall
of pipes. Both together cover "written fast, moved on, never looked back."

The check is a RAGGED-ROW check on separator-confirmed tables: a real GFM table
(one that has a |---|---| delimiter row) whose header or a data row has a cell
count different from the delimiter's. The usual cause is an UNESCAPED pipe
inside a cell or header (a ket like |01>, an absolute value like |rho|, a norm
||R||, a restriction |_w), which GFM reads as an extra column delimiter and which
breaks alignment or the table's recognition. The fix is to escape it as \\| (or
wrap the expression in a backtick code span).

Requiring a delimiter row keeps precision high: bra-ket math display lines and
ASCII art that happen to be pipe-bounded have no separator, so they are never
mistaken for a table.

Pipe-counting is escape- and code-span-aware: escaped pipes (\\|) and pipes
inside backtick code spans (`...`) are protected by GFM and are NOT counted, so
they never false-positive.

Scope = the doc dirs the review backlog covers (docs/ experiments/ hypotheses/
reflections/ recovered/), matching check_em_dashes.py.

Used as part of the repo pre-commit hook (warns, never blocks) and runnable by
hand:
    python .githooks/check_broken_tables.py        # check staged changes
    python .githooks/check_broken_tables.py --all  # sweep the whole tree
"""
import re
import subprocess
import sys

ALL = "--all" in sys.argv
SCOPE = re.compile(r"^(docs|experiments|hypotheses|reflections|recovered)/.*\.md$")
UPIPE = re.compile(r"(?<!\\)\|")     # an unescaped pipe
CODESPAN = re.compile(r"`[^`]*`")    # inline code span: GFM protects its pipes


def git(*a):
    return subprocess.run(["git", *a], capture_output=True, text=True,
                          encoding="utf-8", errors="replace").stdout


def mask(s):
    # blank out backtick code spans so their pipes are not counted as delimiters
    return CODESPAN.sub(lambda m: "x" * len(m.group(0)), s.strip())


def cells(line):
    parts = UPIPE.split(mask(line))
    if parts and parts[0].strip() == "":
        parts = parts[1:]
    if parts and parts[-1].strip() == "":
        parts = parts[:-1]
    return parts


def is_row(line):
    s = line.strip()
    if not (s.startswith("|") and s.endswith("|")):
        return False
    if "<--" in s or "-->" in s:      # ASCII-art arrows, not a table
        return False
    return len(UPIPE.findall(mask(s))) >= 2


def is_sep(line):
    cs = [c.strip() for c in cells(line)]
    return bool(cs) and all(c and re.fullmatch(r":?-{1,}:?", c) for c in cs)


def scan(text):
    out = []  # (line_no, reason)
    lines = text.split("\n")
    i, n = 0, len(lines)
    while i < n:
        if is_row(lines[i]):
            block = []
            while i < n and is_row(lines[i]):
                block.append((i, lines[i]))
                i += 1
            # only a separator-confirmed table is checked (keeps precision high)
            if len(block) >= 2 and is_sep(block[1][1]):
                ref = len(cells(block[1][1]))     # the delimiter = authoritative column count
                hln, hc = block[0]                # the header row
                if len(cells(hc)) != ref:
                    out.append((hln + 1, f"header has {len(cells(hc))} cells vs {ref} in the separator "
                                         f"(unescaped '|' in a header cell? write it '\\|')"))
                for ln, c in block[2:]:           # the data rows
                    k = len(cells(c))
                    if k != ref:
                        out.append((ln + 1, f"row has {k} cells vs {ref} in the separator "
                                            f"(unescaped '|' in a cell? write it '\\|')"))
        else:
            i += 1
    return out


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


hits = {}  # file -> [(line, reason), ...]
for f in files:
    if not SCOPE.match(f):
        continue
    found = scan(content(f))
    if found:
        hits[f] = found

if hits:
    w = sys.stderr.write
    w("\n  \033[33mbroken-table warning\033[0m: these docs contain a Markdown "
      "table that does not render.\n"
      "  A broken table also flags a doc whose full quality round was never done "
      "(the tell is orthogonal to the em dash).\n")
    for f in sorted(hits):
        rows = hits[f]
        lines = ", ".join(str(ln) for ln, _ in rows[:12])
        more = "" if len(rows) <= 12 else f", +{len(rows) - 12} more"
        w(f"    {f}  ({len(rows)} on lines {lines}{more})\n")
        # show the first reason per file as the actionable hint
        w(f"        e.g. L{rows[0][0]}: {rows[0][1]}\n")
    w("  Escape math pipes inside cells as '\\|', or add the missing '|---|' "
      "separator row. (warning only -- commit not blocked)\n\n")

sys.exit(0)  # never block
