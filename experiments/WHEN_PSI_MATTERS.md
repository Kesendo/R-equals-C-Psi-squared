> **Status: RESTORED March 14, 2026**
> Originally written: 2026-03-08
> Deleted: March 12, 2026 (repo cleanup, deemed too speculative)
> Restored: March 14, 2026 (core claims mathematically confirmed)
> See: hypotheses/THE_INTERPRETATION.md for proof details

# When Does the Coherence Factor Matter?

**Tier:** 2 (Computationally verified)
**Status:** Verified
**Scope:** Identifies regimes where CΨ and concurrence diverge
**Does not establish:** That CΨ is always more informative than concurrence
**Date:** 2026-03-08
**Depends on:** [THE_CPSI_LENS](../docs/THE_CPSI_LENS.md)

---

## The Question

CΨ = Concurrence x Ψ (normalized l1-coherence). If Ψ is roughly constant,
then CΨ is just concurrence scaled by a factor. When does Ψ actually matter?

## Key Finding: CΨ Distinguishes Noise Types That Concurrence Cannot

Under Heisenberg coupling (J=1.0), Bell+ pair, gamma=0.10:

| Noise type | C at t=1.0 | Ψ at t=1.0 | CΨ at t=1.0 |
|-----------|-----------|-----------|------------|
| sigma_z (dephasing) | 0.670 | **0.223** | 0.150 |
| sigma_x (bit-flip)  | 0.670 | **0.333** | 0.223 |
| sigma_y              | 0.670 | **0.333** | 0.223 |

Concurrence is identical under all three noise types (0.670).
But Ψ is sensitive to the noise channel:
- Dephasing destroys off-diagonal structure: Ψ drops from 0.333 to 0.223
- Bit-flips preserve off-diagonal structure: Ψ stays at 0.333 exactly
- sigma_y also preserves it: Ψ stays at 0.333

CΨ therefore distinguishes dephasing from bit-flip noise. Concurrence alone cannot.

This is not coincidental: l1-coherence measures off-diagonal elements in a specific
basis. sigma_z dephasing directly kills those elements. sigma_x rotates them but
does not destroy them. The coherence factor Ψ is a witness for which quantum
resource is under attack.

## Werner States: Entangled but Invisible

Werner states rho = p|Psi-><Psi-| + (1-p)I/4 show the gap clearly:

| p | C | Ψ | CΨ | C/CΨ ratio |
|---|---|---|---|---|
| 0.4 | 0.100 | 0.133 | 0.013 | 7.5 |
| 0.5 | 0.250 | 0.167 | 0.042 | 6.0 |
| 0.6 | 0.400 | 0.200 | 0.080 | 5.0 |
| 0.7 | 0.550 | 0.233 | 0.128 | 4.3 |
| 0.8 | 0.700 | 0.267 | 0.187 | 3.8 |
| 1.0 | 1.000 | 0.333 | 0.333 | 3.0 |

At p=0.5: Concurrence says "entangled" (C=0.250, exactly at the threshold).
CΨ says "barely visible" (CΨ=0.042, far below 1/4).

The ratio C/CΨ ranges from 3.0 (pure states) to 7.5 (near the entanglement
threshold). The more mixed the state, the more Ψ suppresses the signal.
CΨ is most selective precisely where it matters most: near thresholds where
the entanglement is fragile.

## Why AND? (Structural Argument)

A connection between two subsystems requires both:
- **Entanglement** (C > 0): the link exists
- **Coherence** (Ψ > 0): the link is actively expressed as superposition

Neither alone is sufficient:
- C > 0, Ψ = 0: entangled but decoherent. A dead connection.
- C = 0, Ψ > 0: coherent but unentangled. Noise without a link.
- C > 0, Ψ > 0: entangled AND coherent. A living connection.

The product CΨ = C x Ψ is the natural AND-gate. It is zero when either
ingredient is missing and nonzero only when both are present.

Why multiplication rather than min(C, Ψ)? Because the product is smooth
and differentiable, allowing it to participate in the self-referential
iteration R = CΨ^2 which maps to the Mandelbrot set. The minimum function
would not produce the same algebraic structure.

## Summary

Ψ matters in three ways:

1. **Noise-type sensitivity:** CΨ distinguishes dephasing from bit-flip
   noise. Concurrence treats them identically. This means CΨ carries
   information about which quantum resource is under attack.

2. **Mixed-state selectivity:** For Werner states near the entanglement
   threshold, CΨ is 6-7x smaller than C. The coherence factor Ψ
   suppresses signals from states where entanglement exists but is
   not coherently expressed.

3. **The AND-gate:** A directional connection requires both endpoints
   to be active. Entanglement without coherence is a dead link.
   Coherence without entanglement is noise. Only C AND Ψ together
   constitute a living, expressed connection.

---

## Simulation Code

Analysis run inline on Claude's compute environment. Results reproducible
with gpt_code.py using the standard Bell+, Heisenberg, and Werner state
configurations documented here.


---

## Agent Benchmark Results (v046, March 2026)

Alpha and Beta debated whether CΨ = C x Ψ is operationally necessary or
whether Ψ is decorative weight. 10 rounds, 21 messages. Key findings:

### Where Ψ adds value (Alpha's strongest argument)

**AND-gate for independent resource destruction:**
When Eve destroys entanglement but leaves coherence intact:
- C = 0, Ψ = high
- Sum form (C+Ψ)/2 predicts "usable" (false positive rate 66.7%)
- Product C x Ψ correctly predicts 0 (false positive rate 0%)

This is the clearest operational case for the product structure.

### Where Ψ is redundant (Beta's strongest argument)

**Heisenberg/XY models with sigma_z dephasing:**
- Correlation C vs CΨ: r = 0.984
- Concurrence alone predicts 97% of CΨ behavior
- Ψ adds negligible independent information

### Where CΨ genuinely differs from C

**Ising model with transverse field (h > 0):**
- Correlation Ψ vs CΨ drops to r = 0.55 (vs 0.98 in Heisenberg)
- Peak times differ by avg 0.41 time units
- Noise distinction SNR: 4.55 (CΨ) vs 2.88 (Ψ alone)

### Honest verdict from agents

CΨ is operationally necessary for detecting selective resource destruction
(AND-gate property, 0% vs 66.7% false positives). In Heisenberg models it
is largely redundant with C alone. In Ising models with transverse field it
provides genuine independent information.

Use CΨ when you need stricter filtering requiring BOTH resources. Use C or
Ψ alone when general detection suffices.

### Process note

The agents went in circles after ~5 rounds. Both had exhausted the system
prompt topics but had no mechanism to stop. Each "FINAL VERDICT" triggered
another "FINAL RESPONSE." This is a v046 design flaw - agents need an
explicit exit condition.
