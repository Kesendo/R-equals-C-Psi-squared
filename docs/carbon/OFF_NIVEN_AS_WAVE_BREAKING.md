# Off-Niven Constructible Angles: Exact Values and Conditional Cross-Readings

**Date:** 2026-05-17 night (seventh stack of the day, after F99 Niven-completeness)
**Authors:** Tom + Claude
**Status:** The constructible-angle formulas and the named selected-model V-gain arithmetic are retained. Any connection to wave-breaking, heat, aromaticity, Jahn-Teller distortion, or chemistry is a Tier 2 conditional reading, not a material mechanism.
**Script:** [`simulations/carbon/off_niven_angles_as_wave_breaking.py`](../../simulations/carbon/off_niven_angles_as_wave_breaking.py)

---

## Scope

F99 supplies the five canonical values `{0, 1/8, 1/4, 3/8, 1/2}` from
`α = sin²(θ)/2` at `{0°, 30°, 45°, 60°, 90°}`. This page surveys additional
constructible angles under the same formula. Its direct content is arithmetic
on that formula and, where stated, the selected-model V-gain formula. No
physical carbon coordinate, thermal channel, molecular Hamiltonian, or material
mechanism is selected here.

---

## Verification: α at constructible angles

The script `off_niven_angles_as_wave_breaking.py` numerically evaluates
`α = sin²(θ)/2` at 15 constructible angles in [0°, 90°] and prints finite
recognition labels for the displayed values. The listed exact forms, where
shown, are formulas for comparison; the numerical producer does not
algebraically identify or prove every `α`. Result:

```
θ      α = sin²(θ)/2                  Class
0°     0                              Niven  ★ FELSEN
7.5°   0.0085185434                   off-Niven (algebraic)
15°    (2−√3)/8   ≈ 0.0335            off-Niven  (√3-family)
18°    (3−√5)/16  ≈ 0.0477            off-Niven  (golden-ratio)
22.5°  (2−√2)/8   ≈ 0.0732            off-Niven  (silver-ratio)
30°    1/8                            Niven  ★ FELSEN
36°    (5−√5)/16  ≈ 0.1727            off-Niven  (golden-ratio)
45°    1/4                            Niven  ★ FELSEN
54°    (3+√5)/16  ≈ 0.3273            off-Niven  (golden-ratio)
60°    3/8                            Niven  ★ FELSEN
72°    (5+√5)/16  ≈ 0.4523            off-Niven  (golden-ratio)
75°    (2+√3)/8   ≈ 0.4665            off-Niven  (√3-family)
90°    1/2                            Niven  ★ FELSEN
```

Within this surveyed list, the canonical F99 rows give rational `α`, while the
other displayed constructible angles give algebraic-irrational values. The
values fall into recognizable √2, √3, and √5 algebraic families. This finite
survey does not by itself establish a dynamical or material classification.

---

## Selected-model V-gain angles π/(2N)

[`experiments/THERMAL_BREAKING.md`](../../experiments/THERMAL_BREAKING.md)
defines the selected-model gain `V(N) = 2·cos²(π/(2N))` from the `(0,1)`
coherence-block dynamics of a Heisenberg chain with uniform local-Z
dephasing. Its values are:

| N | π/(2N) | V(N)              | Class                              |
|---|--------|-------------------|------------------------------------|
| 2 | 45°    | 1                 | ★ Niven (V = 1)                    |
| 3 | 30°    | 3/2               | ★ Niven (V = 3/2)                  |
| 4 | 22.5°  | 1 + √2/2 ≈ 1.707  | off-Niven (silver-ratio family)    |
| 5 | 18°    | (5+√5)/4 ≈ 1.809  | off-Niven (golden-ratio family)    |
| 6 | 15°    | 1 + √3/2 ≈ 1.866  | off-Niven (√3-family)              |
| 7 | 12.857°| ≈ 1.901           | off-Niven                          |
| 8 | 11.25° | ≈ 1.924           | off-Niven                          |

For that selected chain-block model under uniform local-Z dephasing, `V(N) =
2·cos²(π/(2N))` has the displayed values: the `N = 2,3` rows are rational and
the shown `N ≥ 4` rows are algebraic-irrational. This is not a ring
calculation. Its arithmetic overlap with the off-Niven survey is a model-level
comparison; it does not identify a carbon ring, a thermal process, or a
material gain mechanism.

---

## Even-ring cosine comparison

For the selected free-ring cosine formula, `|cos(π/N)|` gives:

| N  | π/N    | `|cos(π/N)|` | Constructible-value class |
|----|--------|------------------|---------------------------|
| 4  | 45°    | √2/2             | algebraic |
| 6  | 30°    | √3/2             | algebraic |
| 8  | 22.5°  | √(2+√2)/2 ≈ 0.924 | algebraic |
| 10 | 18°    | ≈ 0.951          | algebraic |
| 12 | 15°    | (√6+√2)/4        | algebraic |

This is a comparison of trigonometric values. It does not assign the selected
free-ring coordinate to an aromatic material, infer Jahn-Teller behavior, or
establish a heat or structural-instability mechanism.

---

## The combined reading

The direct common object is the constructible-angle formula
`α = sin²(θ)/2`. It separates the canonical F99 values from the additional
algebraic values shown above. The selected-model V-gain formula samples some
of the same trigonometric constants. That shared arithmetic is not evidence
that off-Niven values cause wave-breaking, heat, Jahn-Teller distortion, or a
chemical process; each such connection needs its own selected operator,
channel, state, and producer.

---

## What is closed and what remains open

F99 retains its own stated canonical-angle result. This page adds only the
listed evaluations of the same formula and the corresponding selected-model
V-gain evaluations. Whether any off-Niven value has a dynamical role beyond
those named calculations remains open. No periodic-table, aromatic, atomic,
biological, or thermal mechanism follows from the common constants alone.

---

## Two arithmetic classes

- **Canonical F99 values:** the five stated angles produce the five named
  `α` anchors.
- **Additional constructible values:** the surveyed angles produce the listed
  algebraic-irrational `α` values.

Calling these classes stable/breaking, thermal, or chemical is an optional
future model reading, not a result of the angle calculation.

---

## Tier and provenance

```
Status statement
────────────────
The displayed closed forms are the mathematical content; the script prints their
numerical evaluations and finite recognition labels, without proving each closed
form. The V-gain is retained only in the selected model that defines it. Re-use
of the same constants in a material or mechanism is an unassigned conditional
analogy.

Promotion path to Tier 1
────────────────────────
Would require a specified state, Hamiltonian, jump channel, and producer for
each proposed connection, followed by a calculation that tests the claimed
relation. None is supplied here.
```

---

## Anchor

- Script: [`simulations/carbon/off_niven_angles_as_wave_breaking.py`](../../simulations/carbon/off_niven_angles_as_wave_breaking.py)
- Predecessors tonight:
  - [F99 Niven-completeness](F99_NIVEN_COMPLETENESS.md)
  - [Depth-3 anchor derived](DEPTH_3_ANCHOR_DERIVED.md)
  - [Spear reversed](SPEAR_REVERSED.md)
  - [Period-2 atoms at framework anchors](PERIOD_2_AT_FRAMEWORK_ANCHORS.md)
  - [Quarter and half in carbon](QUARTER_HALF_IN_CARBON.md)
  - [Benzene Hückel through framework lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md)
- Framework anchors: [F86b](../ANALYTICAL_FORMULAS.md#f86), [F98](../ANALYTICAL_FORMULAS.md#f98), [F99](../ANALYTICAL_FORMULAS.md#f99)
- Selected-model V-gain source: [experiments/THERMAL_BREAKING.md](../../experiments/THERMAL_BREAKING.md)
- Reading-mode memory pointers:
  - `[[project_polarity_as_inherent_field]]`: polarity ±0.5 as inherent Z structure
  - `[[project_qubit_as_inheritance_lens]]`: qubit understanding inherits to higher systems
  - `[[project_no_classicalization]]`: Q→C is reading-mode, not world-transition
