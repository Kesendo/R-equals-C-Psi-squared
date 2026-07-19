# The Clock Field Is Site-Owned: Three Pre-Registered Questions to Four Machines

<!-- Keywords: per-site clock field IBM calibration, 1/T2 dephasing rate site property,
ICC variance decomposition qubit, chip total drift common mode, quarter boundary field
four machines, PTF perspectival time field real data, R=CPsi2 clock currency -->

**Status:** Observed on IBM calibration histories (four machines, 91–181 days each)
**Date:** July 18, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Data:** [ibm_torino_history.csv](../data/ibm_history/ibm_torino_history.csv),
[ibm_marrakesh_history.csv](../data/ibm_history/results/ibm_marrakesh_history.csv),
[ibm_kingston_history.csv](../data/ibm_history/results/ibm_kingston_history.csv),
[ibm_fez_history.csv](../data/ibm_history/results/ibm_fez_history.csv)
**Script:** [ptf_clock_field.py](../simulations/ptf_clock_field.py) (deterministic, seed 135) →
[ptf_clock_field_out.txt](../simulations/results/clock_field/ptf_clock_field_out.txt)

---

## What this document is about

The repository already owns a measured per-site clock field, and this
document is the first time it is read as one.

Every day, IBM publishes T1 and T2 for every qubit on its machines. The
rate 1/T2 is the natural clock currency of this project: a genuine
additive rate, the closest hardware reading of the γ that drives every
Lindbladian in the framework (not literally γ: 1/T2 = 1/(2T1) + 1/T_φ,
so it carries a relaxation share alongside pure dephasing). Four calibration histories sit in
[data/ibm_history/](../data/ibm_history/): Torino (133 qubits, 181 days)
and the three Heron-r2 chips Marrakesh, Kingston, Fez (91 days each;
their histories carry 156, 155, and 156 qubits). That is a field of
clocks: one rate per site per day, over months.

We asked this field three questions, with verdict thresholds fixed
before looking (pre-registered in the script's docstring; an
honor-system register, since script and output land in one commit and
the repo cannot independently timestamp the order):

1. **Is the clock site-owned?** Does a qubit's rate belong to the site,
   or is it dominated by day-to-day noise?
2. **Is the chip total conserved?** Do the sites' rates compensate each
   other (total held, shares rotate), drift independently, or drift
   together?
3. **How populated is the ¼ boundary?** How many qubits visit both
   sides of the CΨ = ¼ boundary across the window?

**Vocabulary note, stated up front.** Reading this field as the
Perspectival Time Field (PTF, "each perspective paints the one dynamics
at its own rate", see the [glossary entry](../docs/GLOSSARY.md)) is a
READING. The in-model strong form of the closure law was closed by
[EQ-014](../review/EQ014_FINDINGS.md): Σ ln(α_i) ≈ 0 is an empirical
regularity, not a theorem. Nothing below tests PTF as a model. What is
tested is the structure of the real data itself; the three verdicts
stand on the data alone.

---

## Method

From each history we build a coverage matrix R[q, d] of rate = 1/T2:
qubits present on ≥ 90% of dates, dates on which all selected qubits are
present. All statistics run on ln(rate).

- **Q1 (ICC):** variance decomposition of ln(rate) into between-qubit
  variance V_b (how different the sites are from each other) and
  within-qubit-across-days variance V_w (how much one site wanders).
  ICC = V_b / (V_b + V_w). Pre-registered verdict: ICC > 0.5 means the
  field is a site property.
- **Q2 (total vs shares):** coefficient of variation of the chip total
  T(d) = Σ_q R[q, d], compared against an independence null: each
  qubit's day-series independently permuted, 500 permutations.
  Pre-registered verdicts: CV below the null's 2.5th percentile =
  COMPENSATION (shares rotate, total held); above the 97.5th =
  COMMON-MODE drift (the chip moves together); inside = independent
  drift. Control: the same test on 1/T1.
- **Q3 (boundary):** per-record columns computed by the committed
  [ibm_history_analysis.py](../data/ibm_history/ibm_history_analysis.py):
  `crosses_quarter` is true on a day iff r = T2/(2·T1) < r* = 0.212755
  (the qubit's free-decay trajectory crosses CΨ = ¼ inside the window
  that day); `distance_from_quarter` = min CΨ − ¼, signed. We count
  qubits with at least one crossing day, and qubits whose signed
  distance takes both signs across the window.

Data quality: the histories contain a small share of records violating
the physical bound T2 ≤ 2T1 (calibration fit artifacts: 344 of 24073
records on Torino, 6 / 11 / 97 on Marrakesh / Kingston / Fez, r up to
16.6). They sit far on the non-crossing side and flip no verdict, but
they enter the ICC and CV statistics unfiltered.

---

## Finding 1: the clock is site-owned (all four machines)

| Machine | V_b (between) | V_w (within) | ICC | Verdict |
|---------|---------------|--------------|------|---------|
| Torino | 0.2987 | 0.1062 | **0.738** | site-owned |
| Marrakesh | 0.5751 | 0.1544 | **0.788** | site-owned |
| Kingston | 0.5758 | 0.0912 | **0.863** | site-owned |
| Fez | 0.5166 | 0.0762 | **0.871** | site-owned |

All four clear the pre-registered ICC > 0.5 bar, none marginally. In
plain language: if you want to predict a qubit's dephasing rate
tomorrow, its identity tells you 2.8 to 6.8 times more than the whole
history of its daily wandering. A site's clock rate is a stable
property of that site, persisting across months, not day-noise.

This verdict is conservative: every calibration value carries a
per-measurement fit error, which inflates V_w and therefore biases ICC
downward. The true site-ownership is at least as strong as measured.

---

## Finding 2: two drift regimes, and no compensation anywhere

| Machine | CV(total) | Independence null [2.5%, 97.5%] | Rank | Verdict | 1/T1 control |
|---------|-----------|-------------------------------|------|---------|--------------|
| Torino | 0.0766 | [0.0693, 0.0831] | 60% | independent | independent |
| Marrakesh | 0.0517 | [0.0408, 0.0542] | 88% | independent | borderline (seed-dependent) |
| Kingston | 0.0555 | [0.0299, 0.0383] | 100% | **common-mode** | common-mode |
| Fez | 0.0497 | [0.0343, 0.0441] | 100% | **common-mode** | common-mode |

The four machines split into two regimes. On Kingston and Fez the chip
total varies MORE than independent sites could produce: the whole chip
breathes together (a physical common-mode: temperature,
two-level-system defects, shared control lines). Torino's total is what
independent wandering predicts (mid-band, rank 60%); Marrakesh is
verdicted independent too but sits near the top of its band (rank 88%),
the marginal case. The 1/T1 control corroborates the split where it is
decisive: Kingston and Fez are common-mode in both channels by wide
margins, Torino independent in both. On Marrakesh the control is
uninformative: its CV (0.0502) sits a hair inside the band's upper edge
(0.0504), and re-running the null with other seeds flips that verdict
back and forth. Where the verdicts are firm they agree across channels,
so the regime is a property of the chip's environment, not of the
dephasing channel alone.

**The null result, stated plainly:** the pre-registered weave signature,
COMPENSATION (total held while shares rotate), appears on NO machine.
At chip scale, the measured clock field does not conserve its total.
This is also what the plain physics expects: during calibration the
qubits sit idle and uncoupled, and no known mechanism would make idle
sites compensate each other. The informative half of Q2 is therefore
not the absent compensation but the common-mode EXCESS on Kingston and
Fez, which independent idle sites cannot produce either. Two caveats
bound the compensation null honestly in both directions: measurement error
pushes the test toward "independent" (a weak compensation could be
masked), and physical common-modes push toward "common-mode" (a
compensation overlaid by chip breathing would also be masked). So this
is a bound from the data, not a theorem: no compensation is VISIBLE at
chip scale. If the weave exists in hardware, it lives below the
calibration measurement floor or at sub-chip scale.

---

## Finding 3: the ¼ boundary is a field, not a rim

| Machine | Qubits | Ever cross (≥ 1 day) | Both sides across days | Telegraph median \|Δ ln rate\|/day | p95 |
|---------|--------|----------------------|------------------------|-----------------------------------|------|
| Torino | 133 | 112 | **110** | 0.110 | 0.668 |
| Marrakesh | 156 | 148 | **135** | 0.187 | 1.118 |
| Kingston | 155 | 119 | **95** | 0.061 | 0.642 |
| Fez | 156 | 116 | **104** | 0.140 | 0.727 |

[BOTH_SIDES_VISIBLE.md](../docs/BOTH_SIDES_VISIBLE.md) tracked 16
Torino qubits that oscillate frequently around the ¼ boundary, and its
May 5 update counted the stable archetypes. This is the complementary,
weaker criterion applied to everything: how many qubits touch both
sides at all. The answer: most of every chip. 61–87% of qubits sit on
both sides of the boundary across the window. The ¼ boundary is not a
rim where a few special qubits live; it is a field the whole chip
moves through, machine by machine, at telegraph rates of 6–19% per day
(median).

The two strictnesses are two layers of one picture: the 16 oscillators
are the qubits that LIVE at the boundary; the 95–135 both-siders are
the qubits that VISIT it. The boundary is dense in visits even where it
is sparse in residents.

---

## What this is and is not

**Is:** three verdicts on measured hardware structure, with thresholds
fixed before looking, on a deterministic committed script. The clock
field is site-owned (four machines); chip drift comes in two regimes
with no visible total conservation; the ¼ boundary is densely visited.

**Is not:** a PTF confirmation. The PTF vocabulary above is a reading
(Tier 2 in the [hypothesis's own labeling](../hypotheses/PERSPECTIVAL_TIME_FIELD.md)).
The title's site-owned object is the measured 1/T₂ clock RATE; PTF's α_i
(the per-site time rescaling under a defect) is state-dependent, not
site-owned, in the hypothesis's own tests: a different object, and
Finding 1 neither repairs nor contradicts that result.
And the two null results are different objects:
[EQ-014](../review/EQ014_FINDINGS.md) bounded the in-model closure law
(a sum of per-site time-rescalings under a coupling defect), while
Finding 2 asks a cruder question of the hardware (is the raw chip total
held?) and also finds no conservation. Neither answers the other; they
rhyme, and both leave the weave without a measured instance. What
survives either way: the sites own their clocks, and the boundary
between the regimes is where the chips actually live.

---

## Cross-references

- The boundary story this extends: [BOTH_SIDES_VISIBLE.md](../docs/BOTH_SIDES_VISIBLE.md)
- The reading's hypothesis document: [PERSPECTIVAL_TIME_FIELD.md](../hypotheses/PERSPECTIVAL_TIME_FIELD.md)
- The closure-law bound: [EQ014_FINDINGS.md](../review/EQ014_FINDINGS.md)
- The boundary columns' source: [ibm_history_analysis.py](../data/ibm_history/ibm_history_analysis.py)
- Script and output: [ptf_clock_field.py](../simulations/ptf_clock_field.py),
  [ptf_clock_field_out.txt](../simulations/results/clock_field/ptf_clock_field_out.txt)
