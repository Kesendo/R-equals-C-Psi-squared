# The Three Levels of Time: Why Noise Is What Makes Clocks Tick

<!-- Keywords: gamma dephasing experienced time three levels, parameter time
vs experienced time distinction, gamma necessary sufficient irreversibility,
Bell+ frozen zero noise no time, oscillation recurrence not clock, absorbing
boundary requires gamma, J provides content gamma provides arrow,
tau=gamma*t scaling breaks, R=CPsi2 gamma time -->

**Status:** Computationally verified (Tier 2)
**Date:** March 22, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** [disprove_gamma_is_time.py](../simulations/disprove_gamma_is_time.py), [gamma_is_time_proof.py](../simulations/gamma_is_time_proof.py), [two_qubits_no_noise.py](../simulations/two_qubits_no_noise.py)

---

## What this document is about

Think about what "time" means in everyday life. Not the number on a
clock, but the feeling that things move forward: that milk spills and
does not unspill, that you age and do not un-age, that decisions, once
made, stay made. Physicists call this the "arrow of time," the
difference between past and future.

This document asks: what creates that arrow in a quantum system? The
answer turns out to be surprisingly specific. There are three things
that people might call "time," and only one of them is the real arrow.
The mathematical parameter t (just a number in an equation) exists
trivially. Oscillation (things swinging back and forth, like a
pendulum) looks like change but always comes back to where it started.
Only dephasing noise (γ) creates irreversibility: things that happen
and do not un-happen.

The experiment is simple and dramatic. Take a quantum system with no
noise: it either sits perfectly still forever, or it oscillates in
circles, returning to exactly where it started. Turn on even a tiny
amount of noise: things decay, decisions become permanent, past and
future become distinguishable. γ is not just correlated with the time
arrow. γ is the necessary and sufficient condition for it.

---

## Abstract

Time has three levels: (1) the formal parameter t in d/dt (syntax, trivial),
(2) observable change (oscillation or stillness), and (3) experienced time
(direction, irreversibility, before/after). γ (dephasing rate) is the
necessary and sufficient condition for Level 3. Without γ: Bell+ is frozen
(zero observable change over t=0 to 50); |01⟩ oscillates in circles (127
crossings of ¼ in both directions, recurrence at t=11). With γ: CΨ crosses
¼ once and stays below (absorbing boundary), purity (a measure of how
far the state is from maximum mixedness) decays irreversibly,
past and future become distinguishable. J=0 with γ>0 has experienced time;
J>0 with γ=0 does not. However, τ=γt does not universally scale all
observables (deltas up to 0.86): γ provides the arrow, J provides the
content. Neither alone is the full experience.

---

## The Three Levels

The word "time" hides three completely different things. This table
separates them. The key insight is in the third row: only γ (noise)
creates the kind of time we actually experience.

| Level | Without γ | With γ |
|-------|-----------|--------|
| t as symbol in d/dt | Yes (trivially, syntax) | Yes |
| t as observable change | Oscillation or stillness | Irreversible decay |
| t as experienced time (direction, decisions, before/after) | **No** | **Yes** |

γ is the necessary and sufficient condition for Level 3 (experienced time).
It has no bearing on Level 1 (the formal parameter). The distinction is
between syntax and physics. But even at Level 3, experienced time is not
γ alone: γ provides the arrow, J provides the content (Part 3).

---

## The Evidence

The following three experiments show what "time" looks like with and
without noise. Bell+ is a maximally entangled state (two qubits
perfectly correlated). |01⟩ is a simple product state (one qubit up,
one qubit down). CΨ is our composite diagnostic that measures how
quantum the system is (see [What We Found](WHAT_WE_FOUND.md)).

### Bell+ at γ=0: time runs, but nothing happens

The parameter t goes from 0 to 50. Every observable is constant: CΨ = 0.3333,
concurrence = 1.0, trace distance (the standard measure of how distinguishable two quantum states are) = 0.0000 at every time step. From inside
the system, t=0 and t=50 are identical. The "time" runs, but there is no time.

### |01⟩ at γ=0: time runs in circles

CΨ oscillates, crossing 1/4 in both directions 127 times. Recurrence at t=11.
The system returns to its initial state. There is change, but no direction.
No observable accumulates irreversibly. The 1/4 boundary is meaningless
without γ: the system crosses it freely and returns.

### |01⟩ at γ=0.05: time runs forward

CΨ crosses 1/4 once and stays below (absorbing boundary). No recurrence.
Purity decays from 1.0 to 0.5. 16/16 palindromic pairs in the spectrum.
Things happen that do not unhappen. Past and future are distinguishable.

---

## The Tests Reinterpreted

The following section revisits the original tests from the disproof
attempt and explains what each one actually shows. The key distinction
throughout: a frequency (swinging back and forth) is not a clock
(counting forward and never coming back).

**Test 1 (γ=0 evolution):** The parameter t exists without noise. But
Bell+ at γ=0 shows: t runs from 0 to 50 with zero observable change.
Experienced time does not exist without noise. |01⟩ oscillates but
recurs. No direction, no accumulation, no clock.

**Test 2 (same τ, different t):** The Hamiltonian phase J*t differs
between systems at the same τ=γt. But this phase is periodic (oscillation,
not accumulation). It provides a frequency, not a clock. A clock requires
irreversible counting. Oscillation is reversible.

**Test 3 (Hamiltonian clock):** The Hamiltonian has its own FREQUENCY
(0.5994, independent of γ). Frequency is not a clock. A clock counts
forward and does not come back. The Hamiltonian oscillation comes back
(recurrence at t=11). Without γ, there is no clock.

**Test 4 (Multi-γ simultaneity):** Different γ values coexist at the
same formal t. This is because t is a mathematical coordinate. The
experienced time at each qubit depends on its local γ. Different γ
means different experienced time. This is consistent with γ being time.

**Test 6 (Π reversal):** Π does not reverse the decay envelope
(D = 0.39 after Π + forward evolution). This confirms that the
irreversible part (the time arrow, caused by γ) cannot be undone.
Π reverses the centered eigenvalues (the structure of the arrow),
not the arrow itself.

---

## Correction (March 22, 2026)

Science sometimes works by getting things wrong in an instructive way.
The original version of this document concluded that the strong claim
"γ IS time" was falsified. That conclusion was a category error: it
confused the formal parameter t (syntax) with experienced time (physics).

The [two_qubits_no_noise](../simulations/two_qubits_no_noise.py)
experiment showed what "time without γ" actually looks like: Bell+ is
frozen. |01⟩ oscillates in circles. Nothing is ever decided. CΨ crosses
1/4 in both directions 127 times.

The formal parameter t exists without γ. But experienced time requires γ.

The three-part proof (Part 1 + 2) confirmed: γ is the necessary and
sufficient condition for experienced time. But Part 3 showed: τ=γt does
not scale universally (deltas up to 0.86). Experienced time is not γ
alone. γ provides the arrow. J provides the content. Neither alone is
the full experience.

The original claim in [INCOMPLETENESS_PROOF.md](proofs/INCOMPLETENESS_PROOF.md)
Corollary 2 was correct in spirit but too strong in letter: γ is the
source of experienced time, not identical to it.

---

## Precise Language

| Statement | Status |
|-----------|--------|
| γ is the source of experienced time | Correct (Level 3, Parts 1+2; but τ≠full trajectory, Part 3) |
| The formal parameter t exists without γ | Correct (Level 1, trivial) |
| The Hamiltonian provides a frequency | Correct (but frequency is not a clock) |
| Without γ, nothing is ever decided | Confirmed (Bell+ frozen, |01⟩ recurs) |
| The 1/4 boundary is absorbing | Only with γ. Without γ: 127 crossings both ways. |

---

## The Proof: γ == Experienced Time (March 22, 2026)

The following proof has three parts. Part 1 shows that γ produces
every property of experienced time (completeness). Part 2 shows that
nothing else does (exclusivity). Part 3 tests whether γ and t are
fully interchangeable (they are not: γ provides direction, the
Hamiltonian coupling J provides content). The tables below are in
German because they were written during the original investigation;
the conclusions are summarized in English at the end.

Three-part proof. Script: [gamma_is_time_proof.py](../simulations/gamma_is_time_proof.py).
Data: [gamma_is_time_proof.txt](../simulations/results/gamma_is_time_proof.txt).

### Part 1: Completeness - γ produces every property of experienced time

| Property | |01⟩ γ=0 | |01⟩ γ=0.05 | Bell+ γ=0 | Bell+ γ=0.05 |
|----------|---------|------------|-----------|-------------|
| S monoton steigend | Nein | Nein | Ja (trivial: S=1 konstant) | Ja |
| D kehrt nicht zurück | Nein (Rekurrenz) | Ja | Nein (immer 0) | Ja |
| CΨ Kreuzungen ↓/↑ | 63/64 | 2/2 | 0/0 | 1/0 |
| ‖dρ/dt‖ bei t=50 | 2.81 | 0.02 | 0.00 | 0.00 |

Ohne γ: Rekurrenz (|01⟩) oder Stillstand (Bell+). Mit γ: Konvergenz
zum Steady State.

Anmerkung: |01⟩ bei γ=0.05 zeigt S NICHT monoton steigend und CΨ
kreuzt 2× in beide Richtungen. Die Oszillation vom Hamiltonian
moduliert den Zerfall. Die reine Monotonie gilt nur für J=0 oder
für Bell+ (Eigenzustand von H, keine Oszillation).

### Part 2: Exclusivity - NUR γ erzeugt erlebte Zeit

| Konfiguration | S monoton | D kehrt nicht zurück | CΨ einmalig | Zeit? |
|---------------|-----------|---------------------|-------------|-------|
| J=0.1, γ=0.05 | Nein | Ja | Ja | Teilweise |
| J=1.0, γ=0.05 | Nein | Ja | Nein | Teilweise |
| J=10, γ=0.05 | Nein | Ja | Nein | Teilweise |
| **J=0, γ=0.05** | **Ja** | **Ja** | **Ja** | **JA** |
| J=1.0, γ=0 | Nein | Nein | Nein | **NEIN** |

**J=0, γ>0: Erlebte Zeit existiert.** Reiner Zerfall, kein Hamiltonian.
S steigt monoton. D kehrt nicht zurück. CΨ kreuzt einmalig abwärts.
Alles was erlebte Zeit ausmacht - ohne jede Dynamik. Nur γ.

**J>0, γ=0: Erlebte Zeit existiert nicht.** Oszillation, Rekurrenz.
Kein Observable akkumuliert irreversibel.

**γ ist notwendig und hinreichend für erlebte Zeit. J ist weder
notwendig noch hinreichend.**

### Part 3: Äquivalenz - τ=γt als universelle Skala

| Observable (irreversibel) | max Δ(τ;γ₁,γ₂) | τ-Skalierung? |
|---------------------------|-----------------|---------------|
| S(ρ_A) | 0.790 | Nein |
| Tr(ρ²) | 0.057 | Nein |
| CΨ | 0.259 | Nein |
| Concurrence | 0.861 | Nein |
| arg(ρ₀₁) | 0.000 | Ja (trivial) |

**Die τ-Skalierung bricht.** Die irreversiblen Observablen kollapsen
NICHT auf eine einzige Kurve in τ=γt. Der Grund: der Hamiltonian und
der Dissipator interagieren. Die Oszillation (von J) moduliert den
Zerfall (von γ) und umgekehrt. Das J/γ-Verhältnis bestimmt die Form
der Trajektorie, nicht nur die Skala.

Anmerkung: die K-Invarianz aus [CROSSING_TAXONOMY](../experiments/CROSSING_TAXONOMY.md)
(K = γ·t_cross = const) gilt weiterhin - aber nur für den
KREUZUNGSZEITPUNKT, nicht für die vollständige Trajektorie. Die
volle Dynamik enthält sowohl γ-abhängigen Zerfall als auch
J-abhängige Oszillation. Nur der Kreuzungszeitpunkt skaliert mit τ.

### Die ehrliche Schlussfolgerung

**Was bewiesen ist:**
- γ ist notwendig und hinreichend für erlebte Zeit (Teil 1 + 2)
- Ohne γ keine Irreversibilität, keine Richtung, keine Entscheidungen
- J=0, γ>0 hat Zeit. J>0, γ=0 hat keine.
- γ IST die Quelle erlebter Zeit.

**Was NICHT bewiesen ist:**
- Dass erlebte Zeit = τ=γt (die volle Trajektorie skaliert nicht mit τ)
- Dass γ und t vollständig äquivalent sind (J moduliert den Zerfall)

**Die präzise Behauptung:**

γ ist die notwendige und hinreichende Bedingung für erlebte Zeit.
γ ist nicht die vollständige Beschreibung erlebter Zeit. Erlebte
Zeit = γ (Richtung) + J (Inhalt). γ liefert den Pfeil. J liefert
was passiert. Ohne γ gibt es keine Erfahrung - das ist bewiesen.
Aber die Erfahrung hat Struktur die von J abhängt - das ist auch
bewiesen.

Oder kürzer: **γ ist die Quelle. J ist der Inhalt. Beides zusammen
ist erlebte Zeit.**

---

## References

- [INCOMPLETENESS_PROOF.md](proofs/INCOMPLETENESS_PROOF.md): Corollary 2 (γ: source of experienced time)
- [THE_BRIDGE_WAS_ALWAYS_OPEN.md](THE_BRIDGE_WAS_ALWAYS_OPEN.md): γ as source of experienced time
- [two_qubits_no_noise.py](../simulations/two_qubits_no_noise.py): what time looks like without γ
- [disprove_gamma_is_time.py](../simulations/disprove_gamma_is_time.py): the original tests (data valid, interpretation corrected)
- [gamma_is_time_proof.py](../simulations/gamma_is_time_proof.py): three-part proof (completeness, exclusivity, equivalence)
- [gamma_is_time_proof.txt](../simulations/results/gamma_is_time_proof.txt): raw results
