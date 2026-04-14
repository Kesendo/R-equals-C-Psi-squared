# Merge Log: closed-form

**Batch:** closed-form (9 entries)
**Proposal file:** `OPEN_QUESTIONS_INDEX_PROPOSAL_closed-form.md`
**Status:** 2 / 9 entries resolved

---

## OQ-109 -- superseded

- **Cowork proposal:** superseded
- **Final decision:** superseded
- **Deviation from proposal:** no
- **Session:** 2026-04-14
- **Resolver:** existing repo state (both source and weakness document already carry the reformulation)

**Rationale.**

The entry text is a section heading, not a question: "Active Weakness
#4 (the natural variable u): REFORMULATED" in
`experiments/CRITICAL_SLOWING_AT_THE_CUSP.md` (Section 8). The section
is itself the resolution document, and it explicitly reformulates the
original weakness: u is a conjugation variable that linearises the
Mandelbrot iteration but does not provide a simpler dynamical clock
along real Bell+ trajectories (u ~ 0.61 * Psi^1.02).

The original weakness lives in `docs/WEAKNESSES_OPEN_QUESTIONS.md`
Section "Active weaknesses" item 4, and carries its own explicit
status field: "Partially reformulated (April 2026)" with the same
reformulation in place.

Both source and weakness document are self-documented. No intervention
needed in either file. Only the index needs to reflect the status.

**Category A pattern:** source document already carries the resolution
in-place. Index and merge log must be updated; source document stays
untouched. Distinguish from Category B where the source has an open
question that still needs an in-place status marker.

---

## OQ-209 -- resolved

- **Cowork proposal:** resolved
- **Final decision:** resolved
- **Deviation from proposal:** no
- **Session:** 2026-04-14
- **Resolver:** existing repo state (source document carries the fitted formula)

**Rationale.**

The entry text is a section heading: "5 Threshold formula: ANSWERED"
(numbered 8.5) in `experiments/STAR_TOPOLOGY_OBSERVERS.md`. The section
directly below it contains the resolution: the AB crossing threshold
J_th as a function of gamma was fitted from 10 verified data points
(gamma in {0.001, 0.010, 0.020, 0.050, 0.070, 0.100, 0.120, 0.150,
0.170, 0.200}) and yields:

    J_th(gamma) ~ 7.35 * gamma^1.08 + 1.18    (R^2 = 0.999)

A simple linear approximation J_th ~ 6.39*gamma + 1.16 (R^2 = 0.998)
is also reported. The relationship is smooth and monotonic, no
divergence near gamma = 0.2. The threshold genuinely exists and grows
smoothly with noise.

The entry explicitly says "ANSWERED" in the title, and the answer is
in-place. No intervention needed in the source document.

**Category A pattern** (same as OQ-109): source self-documents.

---

## OQ-040 -- pending
## OQ-044 -- pending
## OQ-048 -- pending
## OQ-100 -- pending
## OQ-122 -- pending
## OQ-139 -- pending
## OQ-264 -- pending
