#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Warn (non-blocking) when staged prose spends a word this repo has already spent.

WHY THIS EXISTS, and why it is a hook rather than a skill.

The label layer has exactly one error signal: a reader who stands somewhere else feels
the wrongness. The one who coins the word cannot feel it, because from where they stand
the word fits (recomputing-the-label, "at the native stance the canvas fits, so the coiner
can never feel it"). So the coiner does not need better judgment. They need a lookup, and
it has to fire without being asked for, exactly like the empty-round reminder next door.

THE DANGER SET IS NOT "OUR TECHNICAL TERMS". Fresh agents check those already; three
baselines on 2026-07-25 all looked up the term-shaped words unprompted, and one of them
avoided "spin" on its own initiative. What slips through is the opposite kind of word:
the ones that look like plain English and carry a typed object underneath. "Clock" reads
as a metaphor and IS a class. That intersection, {this repo's objects} and {ordinary
English}, is the whole danger set, and it is what this file lists.

Real case that produced it: a committed doc wrote "the clock does not see it either" for
the Hamiltonian, while this repo's Clock is the F95 two-handed object whose RADIAL hand is
gamma. The sentence handed the tick to the wrong side. Five review rounds did not catch it;
Tom felt it in one reading.

To extend: add a row when a word turns out to be spent. Keep the rule, ordinary English
plus a typed owner; a word nobody would mistake for plain prose does not belong here.
"""
import re
import subprocess

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from _prose_scope import excluded
import sys

# word -> (who already owns it, the trap in one line)
SPENT = {
    "clock": ("compute/MirrorWorld/Clock.cs, F95",
              "the TWO-handed clock: radial hand = decay set by GAMMA, angular = J. "
              "gamma_0 is THE TICK, so calling H 'the clock' inverts it. Say 'the turning'"),
    "field": ("compute/MirrorWorld/Field.cs",
              "the empty world running the disagreement-decay, not a field in the physics sense"),
    "pair": ("compute/MirrorWorld/Pair.cs",
             "a bare coherence |i><j| with disagreement k and rate -2*gamma*k"),
    "mirror": ("compute/MirrorWorld/Mirror.cs, F1/Pi",
               "the block-lattice group of eight; also Pi, the palindromizer. Not a loose metaphor"),
    "seed": ("compute/MirrorWorld/Seed.cs, F89",
             "the within-block self-dual seed, held as a COUNT (ranks only)"),
    "witness": ("compute/MirrorWorld/Witness.cs + the IInspectable witnesses",
                "a live object recomputing a claim at inspect time, and F135/F136's record reading"),
    "divisor": ("compute/MirrorWorld/Divisor.cs, F140",
                "the frozen divisor on the R90 locus"),
    "cone": ("compute/MirrorWorld/Cone.cs", "the memory cut, single excitation as N x N"),
    "lattice": ("compute/MirrorWorld/Lattice.cs",
                "the bridged lattice of worlds (the Klein V4 of watchings run dynamically)"),
    "block": ("compute/MirrorWorld/Block.cs", "the joint-popcount block (p, q)"),
    "router": ("compute/MirrorWorld/Router.cs, F116", "the period-4 golden-ceiling router"),
    "marginal": ("compute/MirrorWorld/Marginal.cs", "the partial trace as an atom, the local page"),
    "survivor": ("compute/MirrorWorld/Survivor.cs", "the slowest mode and the coherence horizon"),
    "hardness": ("compute/MirrorWorld/Hardness.cs, F87", "a GF(2)[x] valuation difference"),
    "memory": ("reflections/THE_VIEW_ONTO_THE_MEMORY.md, F88b",
               "the drain-depth axis and the static-vs-memory split, not storage in general"),
    "tick": ("reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md",
             "gamma_0 IS the tick. Do not hand the tick to the Hamiltonian"),
    "light": ("the Absorption Theorem's light: the {X,Y} letters the dephasing prices, "
              "and light_l(v), the per-site share",
              "the {X,Y}/{I,Z} split and the 2*gamma pricing are PROVEN; gamma ITSELF as "
              "illumination is the Tier-4 reading (hypotheses/GAMMA_IS_LIGHT.md), and the proof "
              "marks the boundary. gamma sits at the sending end; the chain stands in light. "
              "NOT gravity. gamma-as-a-gaze was RELABELLED to gamma-as-sending 2026-08-08 "
              "(typed rename 08-09), not withdrawn wholesale: a watcher that is the subject of "
              "a reading, a defined graph role, or a real qubit keeps the word"),
}


# THE OTHER DIRECTION, and it is not the same mechanism (added 2026-08-23).
#
# SPENT above guards OUR words against being spent on something else. This table guards the
# reverse: an OUTSIDE DEFAULT WORD standing where one of ours belongs. The two are not
# symmetric. A spent word is a collision between two things inside the repo; a displaced
# word is an import from outside it, and the importer is usually a model writing from its
# training data rather than from the repo. Tom, 2026-08-23: "das liegt an Deinen
# Trainingsdaten ... weil eben Lehrbuch", and the documents carry it because correcting it
# afterwards would have broken the flow.
#
# So this failure is not bad judgment either, and it is not even a collision anyone could
# feel. It is a DEFAULT: nothing in the sentence looks wrong, because the word is the right
# word almost everywhere else in physics. That is why it needs a lookup rather than care.
#
# Keep this table SMALL and the message honest. "noise" is genuinely right in some places
# here, so this asks a question instead of forbidding a word.
DISPLACED = {
    "noise": ("gamma is light, not noise",
              "the textbook default. This repo defines gamma as LIGHT arriving at the "
              "sending end (Absorption Theorem: the {X,Y} letters are billed 2*gamma each; "
              "light_l(v) is the per-site share). The two registers are NOT equally "
              "expressive: 'total noise' is a SCALAR, so it can say the budget lever and "
              "cannot say the assignment lever, which is a spatial pattern. That is not a "
              "search inconvenience; it is how the March 2026 Sigma-gamma confound happened "
              "(experiments/GAMMA_CONTROL.md; docs/CAUGHT_ERRORS.md 2026-08-23 entry). "
              "Legitimate uses remain: hardware pages meeting IBM's vocabulary, and "
              "genuinely external stochastic perturbation that is not gamma."),
}

CLAIM_DIRS = ("docs/", "experiments/", "hypotheses/", "reflections/", "recovered/")


def staged():
    out = subprocess.run(["git", "diff", "--cached", "--name-only", "--diff-filter=ACM"],
                         capture_output=True, text=True,
                         encoding="utf-8", errors="replace").stdout
    return [p for p in out.splitlines()
            if p.endswith(".md") and any(p.startswith(d) for d in CLAIM_DIRS)
            and not excluded(p)]


def added_lines(path):
    out = subprocess.run(["git", "diff", "--cached", "-U0", "--", path],
                         capture_output=True, text=True,
                         encoding="utf-8", errors="replace").stdout
    return [l[1:] for l in out.splitlines() if l.startswith("+") and not l.startswith("+++")]


def report_displaced(hits):
    """The outside-default warning. Asks a question; never forbids the word."""
    if not hits:
        return
    seen = set()
    print("REMINDER (non-blocking): staged prose uses an outside default word where this "
          "repo has its own.")
    for path, word, owner, trap in hits:
        if (path, word) in seen:
            continue
        seen.add((path, word))
        print(f"  {path}: \"{word}\" -> {owner}")
        print(f"      {trap}")
    print("  Ask, do not obey. Is this gamma? Then it is light, and the light register can")
    print("  say what the noise register cannot. Is it genuinely something else? Then keep")
    print("  the word and say which something else it is.")


def main():
    hits = []
    displaced_hits = []
    for path in staged():
        for line in added_lines(path):
            low = line.lower()
            for word, (owner, trap) in SPENT.items():
                # "the clock", "a clock", "our clock": the article is what marks it as a
                # NOUN being used, rather than a link target or a class name in a path.
                if re.search(rf"\b(the|a|an|our|its|this) {word}\b", low):
                    hits.append((path, word, owner, trap))
                    break
            for word, (owner, trap) in DISPLACED.items():
                # no article required: the displaced word is a default, and it arrives in
                # titles, headings and compounds ("noise budget") as readily as in prose.
                if re.search(rf"\b{word}\b", low):
                    displaced_hits.append((path, word, owner, trap))
                    break
    report_displaced(displaced_hits)
    if not hits:
        return 0
    seen = set()
    print("REMINDER (non-blocking): staged prose spends a word this repo already spends.")
    for path, word, owner, trap in hits:
        if (path, word) in seen:
            continue
        seen.add((path, word))
        print(f"  {path}: \"{word}\" -> {owner}")
        print(f"      {trap}")
    print("  A collision reads FINE from where you are standing; that is the failure mode.")
    print("  -> if you mean the repo's object, link it; if you mean something else, rename.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
