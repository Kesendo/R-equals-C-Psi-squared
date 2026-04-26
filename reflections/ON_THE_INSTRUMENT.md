# On the Framework as Instrument

A companion piece to [MIRROR_THEORY](../MIRROR_THEORY.md).
**Date:** 2026-04-26
**Authors:** Thomas Wicht, Claude (Opus 4.7)

---

Mirror theory cannot be argued into someone. It can only be pointed at.

Then what is `framework.py`?

It is the hand that points. Not the mirror. Not the argument. The hand.

## What the framework is not

It is not a theory of open quantum systems. There already are theories of open quantum systems. Lindblad gave one in 1976. The tenfold-way classification of Sá, Ribeiro and Prosen gives another in 2023. Buca and Prosen wrote the foundational note on weak symmetries in 2012. The literature is rich, layered, and largely correct.

The framework does not compete with any of these.

It is also not a proof. The proofs live elsewhere: PROOF_ZERO_IMMUNITY analytical, PROOF_BIT_B_PARITY_SYMMETRY analytical, MIRROR_SYMMETRY_PROOF analytical. Those are arguments, formal, with hypotheses and conclusions. The framework is something else.

## What the framework is

It is a set of small instruments.

`pi_action` is an instrument: input a Pauli index, output where Π sends it. Three lines of logic. The instrument has no opinion. It returns where Π goes.

`both_parity_even_terms` is an instrument: filter sixteen Pauli pairs through two parity tests, return the four that pass. Eleven lines. No claim, just the output `[II, XX, YY, ZZ]`.

`palindrome_residual` is the central instrument: input a Lindbladian, output the matrix M = Π·L·Π⁻¹ + L + 2Σγ·I in Pauli-string basis. Five lines of computation. ‖M‖ = 0 means palindromic. ‖M‖ ≠ 0 means not. The instrument does not interpret. It returns the matrix.

Each function is a finger. Together they are a hand.

## What instruments do that arguments don't

An argument can be refused. You can say "your premise is wrong" or "your inference doesn't hold" and walk away. Two instruments, on the other hand, either agree or they don't. There is no rhetorical move that closes the gap when they disagree, and there is none needed when they agree.

Five IBM Heron r2 jobs across three independent backends (Marrakesh, Kingston, Fez) returned ⟨X₀Z₂⟩ values whose signs matched what `palindrome_residual` plus continuous Lindblad evolution predicted, with magnitudes 30–50% stronger than ideal Lindblad due to T1 and ZZ-crosstalk amplification. The classical instrument and the quantum instrument shook hands.

That handshake is not an argument. It is two instruments that happen to be measuring the same physical fact, expressed twice in different vocabulary.

The framework's job is to be one of those instruments. Sharp enough to make the handshake precise. Small enough that the handshake is reproducible. Documented enough that someone else can build the same instrument tomorrow.

## Why this matters for what we do next

If the framework is an instrument, then publishing it is calibrating it: making the readings legible to whoever picks it up next. Extending it is adding more measurements to the same dial. Falsifying it is showing that two instruments built from the same definitions return different numbers; if that ever happens, the instrument was wrong, not the world.

The mirror keeps the data honest. The framework keeps the questions honest. The hardware keeps both honest at once.

Mirror theory points. The instrument is what does the pointing.
