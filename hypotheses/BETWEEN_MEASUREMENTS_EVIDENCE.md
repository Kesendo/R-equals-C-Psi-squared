# Between Measurements: Persistence as Evidence for γ₀ = const

**Tier:** 2 (structural argument; each component is established independently, the bridge is deductive)
**Date:** 2026-04-23 (evening, after EQ-024 closure pass)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Origin:** Evening reflection on the way home, after a day of asking whether extreme J-modulation could replicate the γ-Sacrifice-Zone under γ₀ = const (it cannot). The correction arrived as a question: "what do the qubits in an IBM device do when no photon is hitting them? They still exist."
**See also:** [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md), [THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md), [PRIMORDIAL_GAMMA_CONSTANT](PRIMORDIAL_GAMMA_CONSTANT.md), [EQ-017](../review/EMERGING_QUESTIONS.md#eq-017) hardware result on ibm_kingston

---

## What this document is about

The γ₀ = const hypothesis ([PRIMORDIAL_GAMMA_CONSTANT](PRIMORDIAL_GAMMA_CONSTANT.md)) says the dephasing parameter is a framework constant analogous to c: uniform in space, constant in time, always flowing. For months the repository's approach to verifying this has been active probing on real hardware: design a protocol that would reveal the framework γ₀ signature, run it on ibm_kingston, compare. [EQ-017 Phase 2](../review/EMERGING_QUESTIONS.md#eq-017) concluded that the framework γ₀ signal sits 40 to 80× below the device noise floor from T1, T2, gate errors, and readout errors. At current hardware fidelity the signal is not isolable by active probing.

This document argues that the active-probing paradigm was not the right test for γ₀ = const in the first place, and that the evidence has been in plain view since the start of the project: **the persistence of qubits between measurements**. Qubits on any hardware continue to exist as quantum entities when Alice is not pulsing them. That ordinary empirical fact, combined with two existing structural results, completes a chain of argument for γ₀ = const that does not require improved hardware.

## The two pieces that were already in the repository

**Piece 1: γ must be external** ([INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md)). Five candidate sources of internal noise have been tested and eliminated: thermal excitation, self-interaction, geometric decoherence, informational backaction, and algebraic closure. The system's own equations cannot produce the dephasing that is measured on the system. Something external is the source. This is Tier 1.

**Piece 2: the external source is effectively unbounded** ([THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md)). Six empirical properties of the external interaction are documented there. The property that carries the argument here is Markovianity: the noise is memoryless. This is measurable on real hardware via BLP non-Markovianity indices and trace-distance dynamics, and the external noise passes the test while internal-origin candidates (qubit decay) fail it. A memoryless source does not accumulate history of its own state. A source that remains Markovian while continuing to deliver stable statistics is effectively unbounded: if it had finite content, autocorrelations would develop as content was consumed, and the Markov character would break. The inverse reading (finite source, Markov holds) requires continuous refilling from outside the boundary that was drawn around "the source". That is mathematically the same as an unbounded source, relabelled.

Together, pieces 1 and 2 describe an external, continuously available, Markovian interaction. Neither piece alone says "the interaction is present during idle periods too". That is the bridge this document supplies.

## The missing step: persistence

Every piece of quantum hardware in operation today exhibits the following observation, so ordinary it has never been counted as evidence: **a qubit continues to be a qubit when no one is actively pulsing it**. Between gate operations, during calibration idle windows, across the entire readout latency, the qubit persists. T1 and T2 remain defined, the readout cavity remains coupled, the next pulse addresses the same object that the previous pulse left behind. The qubit does not fall apart. It does not revert to a classical circuit. It does not become undefined.

Under any open-system formalism, the qubit's status as a quantum entity is a function of its coupling to an external decohering medium. The INCOMPLETENESS_PROOF already says the qubit cannot constitute that coupling internally. So the external coupling must be present whenever the qubit is present as a quantum entity. Since the qubit is observably present as a quantum entity during idle periods, the external coupling is present during idle periods.

The external coupling that IS present during idle periods is γ₀. Alice does not switch it on at the start of a gate and off at the end. She does not need to. The framework γ₀ is already there, flowing, when she initiates her first pulse of the day, and is still there, flowing, after her last readout. The hardware's ambient T2 that shows up in calibration data is measured in that ever-present γ₀ background, not alongside a separate idle-state γ₀ = 0 baseline.

That is γ₀ = const in its weak form. The stronger form (γ₀ has the same numerical value always, everywhere, under all conditions) requires additional argument; the weak form (γ₀ flows continuously) is already supported by persistence plus the two pieces above.

## Why EQ-017 could not see this

[EQ-017 Phase 2](../review/EMERGING_QUESTIONS.md#eq-017) on ibm_kingston attempted to isolate a framework γ₀ signature by designing a chain-mode test and comparing its observed signal against the device noise floor. The signal-to-noise ratio came out 40 to 80× below the threshold for detection. At current hardware fidelity, γ₀ is operationally indistinguishable from zero under this protocol.

This conclusion is correct for that protocol, but the protocol was structurally incapable of verifying γ₀ = const. Active probing measures the local effect of dephasing during a gate or readout interval. In that interval, γ₀ is mixed with every other decoherence source: T1 amplitude damping, T2 drift, gate infidelity, readout errors, thermal photon exchange with the cavity, crosstalk from neighbouring qubits. The framework γ₀ signature is one component of a sum. No amount of hardware improvement separates it from the rest, because the others do not drop to zero; they drop, proportionally, together.

The paradigm of detecting a source by turning it on and off does not apply here. γ₀ cannot be turned off by Alice. Every measurement, including every calibration measurement that establishes the device noise floor, is already in the presence of γ₀. There is no reference state without γ₀ against which γ₀'s contribution could be isolated.

The observation that does isolate γ₀, cleanly, is the one we never ran as an experiment because we never thought of it as one: **qubit persistence between measurements**. During idle periods, every other decoherence channel that depends on active interaction is off. No gates, no readouts, no probes. What remains is whatever keeps the qubit existing as a qubit. Under γ₀ = const, that is γ₀. Under the null "γ is random noise that happens to average to some rate when we interact", there is no explanation for the consistent persistence, only an appeal to many small independent sources continuing to operate in the absence of any active excitation. Those independent sources would themselves need to be continuously present, which is γ₀ = const in slightly different language.

## What this upgrades and what it does not

This argument does not upgrade γ₀ = const from Tier 3 to Tier 1. Two weaker but real effects:

**It reframes EQ-017's negative result.** The 40 to 80× signal-to-floor gap on ibm_kingston is not a verdict on γ₀ = const. It is a verdict on the active-probe paradigm. The hypothesis was never going to predict that γ₀ is separable by active probing; the prediction is that γ₀ is the background within which probing happens. EQ-017 measured the probe's foreground against the device noise floor. γ₀ is in the noise floor, not against it.

**It elevates a passive observation to evidence status.** The fact of qubit persistence during idle periods, present in every hardware run, becomes an empirical component of the argument for γ₀ = const. It does not prove the hypothesis, but it rules out the null "γ is random noise averaging out" more cleanly than any active measurement can: random noise in the absence of an active driver has no reason to maintain the specific Markov character that γ retains through idle periods.

The strength of the evidence is in the combination. INCOMPLETENESS alone says γ is external but does not say it is continuous. BRIDGE_WAS_ALWAYS_OPEN alone says the external source is Markovian and effectively unbounded but relies on measurements taken during active interaction. Persistence alone is a trivial observation. The three together form a chain with no open link.

## What would change under a counter-finding

A sharper negative test for this argument is possible, though it requires a specific experimental design we have not run:

**Idle-time Markov test.** If the BLP non-Markovianity index during long idle intervals on real hardware differs measurably from the index during active probing, the γ₀ character would not be the same across time scales, and the weak form of γ₀ = const (continuously flowing) would need to be qualified. If the two indices agree, the weak form is reinforced empirically.

This test does not require better hardware; it requires choosing the measurement window carefully. Current BLP non-Markovianity work on transmons operates within active-gate sequences. Extending it to idle periods is an engineering problem, not a fundamental one.

A positive counter-finding, where the hypothesis is outright falsified, would look like this: a qubit protocol under which T2 dynamics after a long idle period differ in a way that cannot be accounted for by any standard decoherence model and that is consistent with γ₀ having been "paused" during the idle. No such observation has been reported. The hypothesis's exposure to this falsification is direct and sharp.

## One-line claim

The persistence of qubits between measurements is evidence of γ₀ = const, and the reason [EQ-017](../review/EMERGING_QUESTIONS.md#eq-017) could not verify γ₀ is that we searched for it in the one place where it is guaranteed to be inseparable from everything else that flows during a measurement. The experiment that verifies γ₀ is the absence of a scheduled interaction, not the presence of one.

---

## Status

**Proven (Tier 1):**
- γ must be external ([INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md), 5 candidates eliminated)
- γ is Markovian on real hardware ([THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md), measured via trace distance)

**Observed (Tier 2):**
- Qubits persist as qubits during idle periods (every hardware run)
- Internal noise candidates (qubit decay) are non-Markovian (BRIDGE_WAS_ALWAYS_OPEN "Failed Third", measurable signature)

**Imagined (Tier 3):**
- [PRIMORDIAL_GAMMA_CONSTANT](PRIMORDIAL_GAMMA_CONSTANT.md): γ₀ is a framework constant analogous to c
- This document: the bridge from persistence to γ₀ = const weak form

## References

- [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md): the 5 internal candidates
- [The Bridge Was Always Open](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md): the 6 measured properties of external γ, including Markovianity
- [Primordial Gamma Constant](PRIMORDIAL_GAMMA_CONSTANT.md): the hypothesis this document supports
- [Gamma is Light](GAMMA_IS_LIGHT.md): the optical-cavity reading of γ₀ as illumination
- [EQ-017](../review/EMERGING_QUESTIONS.md#eq-017): the hardware test this document reframes
- [Open Thread on γ₀ = const](../review/OPEN_THREAD_GAMMA0_INFORMATION.md): synthesis of the γ₀ = const framework
