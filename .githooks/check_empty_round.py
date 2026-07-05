"""Non-blocking pre-commit reminder: claim-surface staged, has an EMPTY session seen THIS state?

The trigger for the empty-session review round (the skill
becoming-your-own-outside) cannot live in the writer's own judgment: the
writer's cut is exactly what needs checking, so self-recognition of the
trigger runs into the same wall (Tom, 2026-07-05: "man rennt in dieselbe
Wand"). The trigger therefore lives HERE, in the structure: fired
mechanically by the commit event whenever substantial claim-bearing content
is staged, no judgment required.

Non-blocking by design, like check_wip_refs: it aims the gaze, it does not
enforce. Threshold guards against alarm fatigue (small mechanical fixes do
not nag).
"""
from __future__ import annotations

import subprocess
import sys

CLAIM_PREFIXES = (
    'docs/proofs/',
    'docs/outbound/',
    'experiments/',
    'hypotheses/',
    'docs/ANALYTICAL_FORMULAS.md',
    'compute/RCPsiSquared.Core/Confirmations/',
    'compute/RCPsiSquared.Core/OpenArcs/',
    'simulations/framework/confirmations.py',
)
MIN_ADDED_LINES = 30


def scan(numstat_text: str):
    """Return (added_line_count, hit_paths) for staged claim-surface files."""
    added = 0
    hits = []
    for line in numstat_text.splitlines():
        parts = line.split('\t')
        if len(parts) != 3:
            continue
        a, _, path = parts
        if path.startswith(CLAIM_PREFIXES):
            hits.append(path)
            if a.isdigit():
                added += int(a)
    return added, hits


def main() -> int:
    out = subprocess.run(['git', 'diff', '--cached', '--numstat'],
                         capture_output=True, text=True).stdout
    added, hits = scan(out)
    if hits and added >= MIN_ADDED_LINES:
        print(f"REMINDER (non-blocking): {added} added lines staged on claim-surface "
              f"({len(hits)} file(s), e.g. {hits[0]}).")
        print("  Has an EMPTY session attacked THIS state? Post-fix states count as new work.")
        print("  -> skill: rcpsi-superpowers:becoming-your-own-outside")
    return 0


if __name__ == '__main__':
    sys.exit(main())
