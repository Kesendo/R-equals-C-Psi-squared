# Glossary: Every Term in This Repository, Explained

<!-- Keywords: R=CPsi2 glossary notation definitions, CΨ purity coherence
product, Pi operator palindromic conjugation, Liouvillian superoperator
Lindblad, XY-weight incoherenton number grading, K-invariance crossing
time gamma product, fold catastrophe bifurcation, Baumgratz l1 coherence
normalization, R=CPsi2 glossary -->

**Status:** Living document (Reference)
**Date:** 2026-03-30
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

---

## How to use this glossary

This glossary serves two audiences at once. Every entry has a precise
technical definition (for physicists) and a plain-language explanation
(for everyone else). If you already know quantum mechanics, the first
line of each entry is all you need. If you do not, read the full entry.

If you encounter a term anywhere in the repository and it is not
explained where you find it, it should be here. If it is not here
either, that is a bug; please open an issue.

---

## Core notation

The symbol C is used in three contexts in this repository. This is
a historical artifact: the project started with philosophical labels
("Consciousness," "Possibility," "Reality") that were later replaced
with standard physics terms. The algebra is identical regardless of
naming. See the Notation section in [The CΨ Lens](THE_CPSI_LENS.md) for details.

| Symbol | Algebraic usage | Experimental usage | Legacy usage (deprecated) |
|--------|----------------|-------------------|--------------------------|
| **C** | Coupling scalar in the recurrence | Wootters concurrence | "Consciousness" |
| **Ψ** | Scalar in the recurrence | Normalized l1-coherence | "Possibility" |
| **CΨ** | Product c in the Mandelbrot iteration | Concurrence x coherence (the lens) | "Consciousness x Possibility" |
| **R** | Fixed-point value of the iteration | CΨ^2 (product of concurrence, coherence squared) | "Reality" |

**In plain language:**

- **C (Concurrence)** measures how entangled two quantum particles are. 0 means completely independent; 1 means maximally entangled (what happens to one instantly affects the other). (That last phrase is the pop gloss; its translation lives in [Spooky Action Translated](quantum/SPOOKY_ACTION_TRANSLATED.md).)
- **Ψ (Coherence)** measures how much quantum "superposition" a system has: the ability to be in multiple states at once. 0 means fully classical; 1 means maximally quantum.
- **CΨ** is concurrence times coherence. This single number captures the overall "quantumness" of the system. When CΨ drops below ¼, the system transitions from quantum to classical behavior.
- **R = CΨ²** is the product that gives the project its name. It is the stable resting point of the mathematical iteration: the value the system settles to.

---

## Derived quantities

| Symbol | Definition |
|--------|------------|
| **θ** | arctan(sqrt(4CΨ - 1)). Angular distance from the 1/4 boundary in the complex regime. |
| **τ** | γ * t. Normalized decoherence time. |
| **t_cross** | Time at which CΨ crosses 1/4. |
| **δ** | Tr(ρ^2) - Tr(ρ_predicted^2). Purity residual from a simple noise model. |

**In plain language:**

- **θ (theta)** is a compass. It tells you how far the system is from the quantum-classical boundary at ¼. When θ is large, the system is deep in quantum territory. When θ approaches zero, the system is about to cross into classical behavior. Validated on IBM hardware at 0.3% accuracy.
- **τ (tau)** is time measured in units of decoherence. If noise is strong (high γ), one second of real time corresponds to many units of τ. This lets us compare systems with different noise levels on the same scale.
- **t_cross** is the moment the system crosses the ¼ boundary: the instant quantum becomes classical.
- **δ (delta)** measures how much the actual system deviates from a simple noise prediction. If δ is large, something interesting is happening that simple noise models cannot explain.

---

## Palindromic symmetry (proven March 14, 2026)

The palindromic symmetry is the central discovery of this project.
For the full proof, see [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md).
For what it means in plain language, see [What We Found](WHAT_WE_FOUND.md).

| Symbol / Term | Definition |
|---------------|------------|
| **Π** | Conjugation operator. Per-site action: I→X, X→I, Y→iZ, Z→iY. Satisfies Π·L·Π⁻¹ = -L - 2Σγ·I. Maps every decay rate d to its palindromic partner 2Σγ-d. This is one member of the P1 family. A second family P4 (I↔Y, X↔Z) supports additional terms (XZ, ZX). XY/YX terms require non-uniform (alternating) operators. See [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md), [Non-Heisenberg Palindrome](../experiments/NON_HEISENBERG_PALINDROME.md). |
| **θ (transpose mirror, D)** | The transpose, ρ → ρᵀ, read as a mirror of operator space (distinct from the compass angle θ above). On the Pauli basis it is the diagonal sign diag((−1)^(n_Y)): the antiautomorphism of the Pauli algebra, whose character is y_par = n_Y mod 2, the one single-letter parity no unitary conjugation can read. This is F114's D, and one of the two factors of Π: Π = R·D. See [Π Factors as R·D](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md). |
| **Mirror group D₄** | The group generated by the ket reflection R (ρ → ρ·X^⊗N) and the transpose mirror D: eight mirrors that close exactly into the dihedral group of the square (F118). Π = R·D is its 90° rotation; every named mirror of the palindrome story (Π, Π_Y, the charge conjugation 𝓕, F114's D) is an element; the windowed-converse spine's involutions {I, 𝓕, R, 𝓕R} are the Klein four-group V₄ inside it. See [Π Factors as R·D](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md). **Glyph note:** the F71 bond-mirror family (F100/F101) reuses the letters R and D for different objects, the *spatial* chain-mirror R (site i↔N−1−i) and the bond-deviation function D(b) := c₁(b)−c₁(N−2−b), not these operator-space generators. |
| **Polarity cube** | The Z₂³ grading (bit_a, bit_b, y_par) of Pauli strings on which the F87 refinement family (F102-F111) lives. Its three axes are the characters of two conjugations and the transpose: bit_a = conjugation by Z^⊗N, bit_b = conjugation by X^⊗N, y_par = the transpose θ (F118). The y_par axis is the antiautomorphism dimension, invisible to every unitary conjugation. |
| **Antilinear triangle** | The three involutions transpose θ, entrywise conjugation conj, and adjoint † (with † = θ∘conj): a Klein four-group graded by linearity ℓ and multiplicativity m, whose product character is the transport sign μ∘L_H∘μ = ℓm·L_{μ(H)} (F119). Five proofs share this engine (F114, the reversal kill, F112 Lemmas A+B, F113 Lemma C, the K₁/K_b mirrors). In the Pauli basis θ = D and † = the antilinear unit 𝒦, so ⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ is the antilinear double of the mirror group. See [The Antilinear Triangle](proofs/PROOF_ANTILINEAR_TRIANGLE.md). |
| **Polarity content ‖M_anti‖²** | The size of the part of the F81 residual M that the Π-polarity split can see at all: ‖M_anti‖² = ‖M₊₁⁄₂‖² + ‖M₋₁⁄₂‖², the sum of M's two Π-conjugation ±i halves. Identically M_anti = (M − Π·M·Π⁻¹)/2, so it vanishes exactly when L carries no Π²-odd content. It is the denominator of the **relative asymmetry** rel = \|‖M₊‖² − ‖M₋‖²\| / ‖M_anti‖², the whole of what that difference is a difference of, which makes rel a contrast ratio in [0, 1]. Not to be confused with ‖M‖², a strictly larger and differently-scaled quantity: on the 2-body bilinear family under Z-dephasing, F83's anti-fraction gives ‖M‖²/‖M_anti‖² = 2 + 4r exactly, where r = ‖H_even_nontruly‖²_F / ‖H_odd‖²_F (measured 2, 6, 10 at r = 0, 1, 2 for N = 2, 3, 4; r = ∞, i.e. H with no Π²-odd part, is the degenerate case where ‖M_anti‖² = 0 while ‖M‖² stays large). Dividing by ‖M‖² was retired repo-wide on 2026-08-06, because it saturates wherever M vanishes and is the wrong scale even where it does not. Where the polarity content is zero the ratio has no numerator either and the configuration carries **no reading**, balanced or broken; the exact structural companions are `HasNoPi2OddPart` (C#) and, for the c side alone, bit_b-homogeneity. See [F112](proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md), [F81](proofs/PROOF_F81_PI_CONJUGATION_OF_M.md), [Polarity Coordinates](../reflections/POLARITY_COORDINATES.md). |
| **Polarity asymmetry** | The single number the whole polarity layer reads: asymmetry = ‖M₊₁⁄₂‖² − ‖M₋₁⁄₂‖², the **signed** difference of the two Π-conjugation ±i halves of the F81 residual M. Its ABSOLUTE VALUE is the numerator of the relative asymmetry defined in the row above, which is why that ratio lands in [0, 1] while this number carries a sign. It is what [F112](proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) proves is exactly zero for a full Lindblad-form L with Hermitian H and bit_b-homogeneous collapse operators (and, in [its non-Hermitian extension](proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md), for a widened H with that same hypothesis on the collapse operators unchanged), what [F113](proofs/PROOF_F113_COEFFICIENT_DERIVATION.md) gives in closed form for a Z-drive against amplitude damping, and what [F155](proofs/PROOF_F155_PHYSICAL_GENERATOR_POLARITY_BREAK.md) gives in closed form for the physical generator. Three things about it are convention-free and one is not: the MAGNITUDE, the zero-versus-nonzero verdict, and the ratio of signs between any two contributing Pauli strings are properties of the physics; the overall SIGN is fixed jointly by the stacking pairing (row-stack L against an order='F' transform, documented on `PauliBasis`) and the phase convention of Π, so a number is quoted in a pairing or it is not quoted. The VALUE is additionally not basis-invariant: it moves under a per-site Hadamard, so it is a reading of the physics PAIRED WITH a palindrome convention rather than of the physics alone. Live: `inspect --root f155`. |
| **The commutator, and the physical (no-jump) generator** | Two different superoperators, and confusing them is recorded in [CAUGHT_ERRORS](CAUGHT_ERRORS.md) as one entry with a sequel, 2026-06-20 and 2026-08-19. The **commutator** is ρ ↦ −i[H, ·], the closed-system generator; [F112's non-Hermitian extension](proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md) is about a widened H there, while [F112's own balance theorem](proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) is about a full Lindblad-form L. The **physical or no-jump generator** is ρ ↦ −i(Hρ − ρH†) with H = A + iB non-Hermitian: the un-normalised conditional generator of post-selected dynamics, the operator PT and gain-loss models are usually written with, obtained from a Lindbladian by dropping the recycling term. It is not trace-preserving, and it is not the **full Lindbladian**, which adds Σ_k c_k ρ c_k† back. The three CAN carry three different polarity readings, and where they agree that is a measured fact rather than a definition: for Hermitian H the first two are the same operator, and the recycling term is measured asymmetry-inert on the σ⁻/σ⁺ family (N = 2, 3, 4) and not outside it. What never decides is a LABEL on the system (PT-symmetric, gain-loss, balanced-loss): the verdict is the value of F155's sum. Typed as `NoJumpGenerator` (`compute/RCPsiSquared.Core/Lindblad/`). |
| **Σγ (total noise, total illumination)** | Sum of all individual dephasing rates: Σγ = γ₁ + γ₂ + ... + γ_N. The palindromic axis sits at -Σγ; the paired-eigenvalue form and its scope are in the Palindromic pairing row below. What Σγ fixes on its own is the trace of L, -Σγ·4^N (F91 in [the formula registry](ANALYTICAL_FORMULAS.md); that file carries no heading anchors, search for "F91."). What it does not fix is the ASSIGNMENT, the profile γ_l across the sites. F91 gives the assignment's reach on its own object: the F71-refined diagonal-block matrix elements depend only on the indexed pair-sums S_l = γ_l + γ_{N−1−l}, never on the individual γ_l, so within-pair redistribution leaves those blocks alone while a different assignment of pair-sums to indices does not (two N=6 profiles sharing Σγ and the multiset {1.1, 0.8, 0.8} have diagonal-block spectra differing by 8.69). That is a statement about blocks and not about rates: F91 adds that the breaking reaches the full-L spectrum and the physical rates do move. In transport the assignment is a large effect: on the N=5 mediator chain at fixed Σγ = 0.13 the best and worst arrangements differ by 51% in end-pair mutual information, with concentration at the chain centre the cheapest seat ([Gamma Control](../experiments/GAMMA_CONTROL.md)). **Σγ and the assignment are two levers and neither ranks the other:** within one shape less Σγ is monotonically better, at fixed Σγ the arrangement is worth up to 51%, and the two CROSS, since the whole of 0.25 on the centre beats a uniform 0.13. So a comparison whose arms differ in both confounds them and can be attributed to neither; match Σγ by construction, or quote both arms' Σγ beside the result ([CAUGHT_ERRORS](CAUGHT_ERRORS.md), 2026-08-23 entry, which records both this confound and the mirror error the first repair made). In the light register Σγ is the total illumination and the assignment is the pattern of γ_l. The Concentrator rows below name the chain EDGE as the seat: that is a different design problem, and [Resonant Return](../experiments/RESONANT_RETURN.md) records the two as non-mixable rules, edge for total network bandwidth and centre for point-to-point throughput, with no compromise profile. |
| **Light and lens ({X,Y} against {I,Z})** | **Proven:** for a Hermitian H with Z-dephasing and one Pauli string per jump operator, the dissipator sorts every Pauli letter into two classes, the **lens** {I, Z}, passed at no cost, and the **light** {X, Y}, billed 2γ per factor. Theorem 1, at uniform rate γ per site: Re(λ) = -2γ⟨n_XY⟩, the mode's light content. Theorem 2, for site-dependent rates: Re(λ) = -2 Σ_l γ_l·light_l(v), where light_l(v) is the eigenmode's own exposure at site l and not a rate anyone chooses; the chosen quantity is γ_l. Other jump letters keep the form and change the count, and under depolarizing the split collapses to identity against everything (Re(λ) = -4γ⟨n_nonI⟩), which is a different grading from the one Π pairs, so the Palindromic pairing row's depolarizing caveat is not in tension with this. Π swaps light and lens, hence ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N. Eigenmodes carrying no light are immortal; the ladder's top rungs (2Nγ, and F3's 2(N−1)γ as the highest ordinary band edge) belong to the number-conserving XY/Heisenberg family above Q\*_gap(N), not to a generic Hermitian H, where no mode need reach 2Nγ at all. An {I, Z} string that H rotates into the light is not an eigenmode and does decay, so the light content is read off an eigenvector after diagonalising rather than off a string. **Interpreted (Tier 4, and the proof marks it so):** γ is illumination and {I, Z} is the lens; the chain is not looked at, it stands in light, and γ sits at the sending end rather than at an eye. See [Gamma is Light](../hypotheses/GAMMA_IS_LIGHT.md) for that reading and its fences. On 2026-08-08 γ-as-a-gaze was recomputed to γ-as-sending, with the typed rename following on 08-09 ([The Label Map](quantum/THE_LABEL_MAP.md), which draws the fence: where a watcher is the subject of a reading, a defined graph role, or a real qubit, the word stays; where γ or the dissipator is made the agent, it goes). The graph-role Watcher, the hardware Watcher factor and the Observer qubits further down are on the staying side. **Glyph note:** "lens" also names the CΨ lens ([The CΨ Lens](THE_CPSI_LENS.md), and the CΨ row above) and the slow-mode lens survey in `RCPsiSquared.Compute`; this one is the {I, Z} half of the letter split. See [The Absorption Theorem](proofs/PROOF_ABSORPTION_THEOREM.md). |
| **Palindromic pairing** | Every non-zero Liouvillian eigenvalue d has a partner d' such that Re(d) + Re(d') = -2Σγ. Proven for Heisenberg+dephasing, verified for all standard models (XY, Ising, XXZ, DM) under single-axis dephasing. Depolarizing noise breaks palindrome at err ≈ (2/3)Σγ. Verified N=2-8. |
| **Pauli weight** | Number of non-identity Pauli operators in a Pauli string. E.g., XYI has weight 2. The mirror-relevant count is the XY-weight, the number of X or Y letters only (ZZI has XY-weight 0): Π maps XY-weight k to N-k (complementarity); the total non-identity weight is NOT flipped (see [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md), Step 1). "XY-weight" as in Haga et al. (2023). |
| **Incoherenton** | Term from Haga et al. (2023). Quasiparticle that counts Pauli weight in open quantum systems. Their XY-weight = our Pauli weight. Their particle-hole transformation = our Π. |

**In plain language:**

- **Π (Pi)** is the operator that *proves* the palindrome exists. It is a transformation that swaps every decay rate with its mirror partner. If you know one half of the spectrum, Π tells you exactly what the other half looks like. It is the mathematical reason the palindrome is not a coincidence but a theorem.
- **θ (the transpose mirror)** is the humblest mirror in the inventory: flip a matrix across its diagonal. Of the three Pauli letters only Y changes sign under this flip, so the transpose is the one mirror that can count Y's; no rotation can. It turns out to be the hidden half of Π.
- **Mirror group D₄** answers the question "how many mirrors are there?" Eight, and they close: every mirror this project ever named is one of the eight symmetries of a square. Π itself, the door the project walked through first, has two hinges: Π = R·D, a reflection times the transpose, and each factor carries one half of the palindrome (the reflection moves the noise and the shift, the transpose flips the Hamiltonian).
- **Polarity cube** is the three-axis coordinate system the classifier reads on every interaction term. Two of its axes come from ordinary rotations; the third, y_par, only the transpose can see. That is why the third axis always felt different: it is the mirror-image dimension, not a rotation dimension.
- **Polarity content** is the answer to "how much is there to be asymmetric about?" The polarity reading compares two halves of one quantity, and the honest way to report their difference is as a fraction of their sum, never as a fraction of some larger quantity they both sit inside. When the sum is zero there is nothing to compare and the reading is empty, which is a third answer alongside balanced and broken. For years the ratio was taken against that larger quantity instead, which made empty readings look balanced.
- **Polarity asymmetry** is the number the polarity layer exists to report: how far the two halves are from each other. Its size and its being zero or not are facts about the system; its plus-or-minus sign is a bookkeeping choice, fixed once and then held. And the number itself is read in a chosen frame: turn the system while holding the mirror still and the reading changes, which is why an asymmetry is always quoted together with the mirror convention it was read against, never on its own. Turning the mirror *with* the system is a different move and changes nothing, which is the Klein transport, not this.
- **The commutator and the physical generator** are two different things a formula can mean by "the system's evolution", and the difference has cost this project a proof's worth of scope, in an error it made once and then repeated while repairing. One is the closed system turning. The other is the system that is also leaking, watched only on the runs where nothing was seen to leak. A third is the honest open system, leaks and all. They can give three different answers; sometimes two of them coincide, and where they do, that is something we measured rather than something the words guarantee. What is never a guarantee is a label on the system, however respectable: only the arithmetic decides.
- **Σγ (total noise)** is the total: how much light falls on the chain altogether. The palindrome is always symmetric around it, and more of it shifts the whole palindrome without breaking it. It is one of two levers, and it is the one that is a single number. The other lever is *where* the light falls, and no single number can carry that. Neither lever outranks the other, and they cross: a well-placed larger total can beat a badly-placed smaller one. So a comparison whose two arms differ in both the total and the placement has measured neither, and the repair is to hold the total equal and vary only the placement.
- **Light and lens** are the two kinds of letter a Pauli string can carry: the dephasing is blind to one kind and charges for the other. A mode lives as long as its share of the charged kind allows, and the palindrome is the statement that a fast mode and its partner divide the letters between them. One care is needed with strings: the rule is about the system's own modes, and an all-lens string that the coupling turns into the light is not one of them. Calling γ itself "the light" is this project's reading of that algebra rather than a consequence of it, and the proof says so in its own words.
- **Palindromic pairing** means: for every rate at which the system loses information, there is a partner rate such that the two add up to the same total. Like a list of numbers that reads the same forwards and backwards. Verified for 87,376 eigenvalues across systems of 2 to 8 qubits, zero exceptions.
- **Pauli weight** counts how many qubits in a Pauli string are "active" (not identity). For the mirror, the count that matters is narrower: how many sites carry X or Y, the letters the light prices. The palindrome maps strings with k such sites to strings with N-k of them, a complementarity between "mostly quiet" and "mostly lit"; a string of Z's counts as fully quiet for the mirror even though its qubits are active.
- **Incoherenton** is a term from the published literature (Haga et al., 2023) for the same concept we discovered independently. Their "particle-hole transformation" is our Π operator. Different names, same mathematics.

---

## XOR space (discovered March 16, 2026)

XOR space is the "fast drain" of the palindrome: the set of modes
where quantum information is destroyed most quickly. Understanding
which states fall into this drain, and which avoid it, is essential
for protecting quantum information.

| Term | Definition |
|------|------------|
| **XOR modes** | The eigenmodes at λ = -2Σγ, the maximum decay rate. There are N+1 of them on a number-conserving XY/Heisenberg H; the count is a fact of that family and not of Z-dephasing, see the -2Σγ row in the constants table. They are not paired with other modes; their palindromic partner is the steady state (λ=0). They are purely off-diagonal (coherences). |
| **XOR fraction** | Fraction of an input state's weight that projects onto the XOR modes. GHZ: 100%. W (N≥3): 0%. Predicted by mixed XY Pauli weight at r=0.976 (N≥3). |
| **Mixed XY Pauli weight** | Fraction of Pauli decomposition terms containing both X and Y operators simultaneously (e.g., XYI, YXZ). Predicts how much of a state falls into the fastest-decaying XOR drain. |
| **Palindromic modes** | All Liouvillian modes except the XOR modes and steady state. They come in pairs at various decay rates. States in palindromic modes are more robust because some pairs decay slowly. |
| **Spectral filter** | The palindrome acts as a filter: it separates every input into a fragile XOR component (fast decay, mixed XY) and a distributable palindromic component (various rates, some survive). |

**In plain language:**

- **XOR modes** are where quantum information goes to die. They decay at the maximum possible rate. The famous GHZ state ("Schrödinger's cat") falls 100% into these modes, which is why GHZ states are so fragile: all of their information sits in the fastest drain. (The cat's full translation lives in [Schrödinger's Cat Translated](quantum/SCHRODINGERS_CAT_TRANSLATED.md).)
- **XOR fraction** tells you how fragile a quantum state is. 100% means it will be destroyed as fast as physically possible. 0% means it avoids the drain entirely. The W state (another famous entangled state) has 0% XOR fraction and is correspondingly much more robust.
- **Spectral filter** is the key insight: the palindrome automatically sorts every quantum state into a fragile part (XOR) and a robust part (palindromic modes). This is not a design choice; it is a consequence of the mathematics. It means the palindrome tells you, for free, which part of your quantum information will survive and which will not.

---

## Quantum state transfer (verified March 14, 2026)

Quantum state transfer (QST) is the problem of moving quantum
information from one place to another through a noisy network.

| Term | Definition |
|------|------------|
| **QST** | Quantum State Transfer. Moving a quantum state from Alice (A) to Bob (B) through a mediator (S). |
| **F_avg** | Average fidelity of the transferred state, averaged over all pure input states. Our best: F_avg = 0.886 (star, 2:1 coupling). |
| **Holevo capacity** | Maximum classical information transmittable per channel use. Our star channel: χ = 0.534 bits. |
| **2:1 coupling** | Optimal coupling ratio J_SB/J_SA = 2 for star-topology QST. Asymmetric. Not intuitive. Outperforms symmetric 1:1. |
| **Echo** | Entanglement oscillation in the SA pair. Period ~π/(4J). Envelope decays at 8γ/3 (the J/γ → ∞ value; a band at moderate coupling, F33). Peak C_SB = 0.598 (N=3). Scales as ~1/(N-1). |

**In plain language:**

- **QST** is the quantum version of "telephone": Alice wants to send Bob a quantum message through a shared intermediary. The question is how accurately the message arrives.
- **F_avg (fidelity)** measures transmission quality. 1.0 means perfect; 0.5 means no better than guessing. Our best result (0.888) means 89% of the quantum information arrives intact.
- **2:1 coupling** is a surprising finding: the best transmission happens when the connection between the intermediary and Bob is *twice as strong* as the connection between Alice and the intermediary. Your intuition says equal connections should work best. It does not.
- **Echo** is the "heartbeat" of the quantum connection: entanglement oscillates back and forth between the qubits, getting weaker with each cycle. The rate at which the heartbeat fades tells you how long the quantum channel stays useful.

---

## Concentrator formula (discovered March 24, 2026)

The concentrator formula is the project's most dramatic engineering
result. For the full discovery story, see
[Resonant Return](../experiments/RESONANT_RETURN.md).

| Term | Definition |
|------|------------|
| **Concentrator** (formerly "sacrifice zone") | A single edge qubit that absorbs the entire noise budget while the remaining N-1 qubits are protected at minimal dephasing. The optimal spatial dephasing profile for information transfer. The name was corrected 2026-03-28: the edge qubit sacrifices nothing, it concentrates the noise into structure for the interior. |
| **Concentrator formula** | gamma_edge = N * gamma_base - (N-1) * epsilon, gamma_other = epsilon. Concentrate all noise on one edge qubit, protect the rest. |
| **V-shape profile** | Hand-designed dephasing profile with edges higher than center (e.g., [0.070, 0.060, 0.050, 0.060, 0.070]). Baseline for comparison. |
| **ENAQT** | Environment-Assisted Quantum Transport. Field founded by Plenio & Huelga (2008). Optimizes a uniform scalar dephasing rate for transport efficiency. Achieves 2-3x improvement. Does not consider spatial profiles. |
| **Sum-MI** | Sum of mutual information between all adjacent qubit pairs. The observable used to evaluate dephasing profile quality in the Resonant Return experiments. |
| **DE optimizer** | Differential Evolution global optimizer. Found 100x (N=7) in 90 minutes. The formula found 180x in 3 seconds. |

**In plain language:**

- **Concentrator** is the core idea: let one qubit at the end of the chain absorb all the noise, and protect every other qubit as much as possible. The edge qubit sacrifices nothing; it was never carrying the information, and it turns the noise it gathers into structure that lets the rest of the chain stay coherent. This simple strategy beats 18 years of published optimization research by a factor of 100. (The early name "sacrifice zone" was a misnomer, corrected 2026-03-28.)
- **ENAQT** is the established research field that studies how noise can *help* quantum transport (counterintuitively, a little noise can improve information transfer). That field achieved 2-3× improvement by tuning the *amount* of noise. We achieved 139-360× by tuning *where* the noise goes (in simulation, peak created Sum-MI, a transport metric; ~2-3× on hardware). Nobody before this work had optimized the spatial distribution.
- **Sum-MI (mutual information)** measures how much information neighboring qubits share. Higher means more information is flowing through the chain. This is the score we used to judge whether a noise profile is good or bad.

---

## Resonator and stability (March 25-29, 2026)

| Term | Definition |
|------|------------|
| **V-Effect** | Complexity emergence when two palindromic systems couple through a mediator. Breaking releases frequency diversity: 4 frequencies become 11 at N=3; two dead N=2 resonators produce 109 new frequencies. Named for the V-shaped bifurcation diagram. See [V-Effect Palindrome](../experiments/V_EFFECT_PALINDROME.md). |
| **Resonator** | The palindromic system is a Fabry-Perot resonator, not a communication channel. CΨ_max is the inner mirror, CΨ = ¼ is the outer mirror. The heartbeat (81 crossings at J=5.0) is a cavity round-trip. Discrete cavity modes at J=2 (Q=7) and J=12 (Q=11). Cavity modes at Σγ = 0: Stationary(N) = Sum_J m(J,N)*(2J+1)^2 (Clebsch-Gordan). Star has N harmonic frequencies (2kJ) under the **isotropic Heisenberg** interaction, where SU(2) collapses them; the chain has a rich irrational spectrum. Both halves are Heisenberg statements: under XY the star has neither N frequencies nor harmonic ones (N = 4 gives 2 − √3, √3, 2√2 among others). See [Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md), [Cavity Modes Formula](../experiments/CAVITY_MODES_FORMULA.md). |
| **Zero Is the Mirror** | At Σγ = 0 the palindrome equation reduces to Π·L·Π⁻¹ = -L: the unitary ground state. Pure oscillation, no decay. Noise shifts the palindrome from zero; the fold, crossing, and concentrator are geometry of that shift. See [Zero Is the Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md). |
| **Fold threshold** | Σγ_crit/J ~ 0.25-0.50%, N-independent (max/min = 1.015 for N=2..5). The minimum noise required for the fold at CΨ = ¼ to exist. Below this: pure oscillation. Above: irreversibility. A universal constant of the palindrome geometry. |
| **Fragile Bridge** | Coupled gain-loss system (one side decays, the other amplifies, Σγ_total = 0). Finite stability window. Instability is Hopf bifurcation (oscillating divergence), identified as Liouvillian chiral symmetry breaking (class AIII). Three regimes: linear, optimal (2x internal coupling), 1/J decay. Asymptotic constant γ_crit x J_bridge = 0.50. See [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md), [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md). |

**In plain language:**

- **V-Effect** is the discovery that connection creates complexity. Two quantum systems that are individually simple (few frequencies, no oscillation) suddenly produce dozens or hundreds of new frequencies when you connect them. This is not the parts becoming louder; it is genuinely new behavior that neither part had alone. See [V-Effect Palindrome](../experiments/V_EFFECT_PALINDROME.md) for the full story.
- **Resonator** means the palindromic system behaves like a hall of mirrors for quantum information. Information bounces back and forth between two boundaries (CΨ_max and CΨ = ¼), creating a standing wave. Each bounce is a "heartbeat." This reframes the system: it is not a wire that sends information from A to B, but a resonating cavity that holds information in standing patterns.
- **Zero Is the Mirror** is the insight that without any noise, the palindrome reduces to its simplest form: pure oscillation, no decay. Everything else (the ¼ boundary, the crossing, the concentrator) is what happens when noise shifts this perfect mirror away from zero.
- **Fold threshold** is the minimum amount of noise needed for the quantum-classical boundary to exist at all. Below this amount, the system oscillates forever without decaying. Above it, irreversibility begins. This threshold is nearly the same regardless of system size: a universal constant.

---

## Parameters

| Symbol | Definition |
|--------|------------|
| **γ** | Decoherence rate (Lindblad dephasing strength), per site γ_l. Higher γ means faster loss of coherence. In star topology each qubit has its own; the same γ on the receiver (γ_A) costs more than on the sender (γ_B). What the Absorption Theorem proves about it is in the **Light and lens** row above: γ prices the {X, Y} letters at 2γ each and passes {I, Z} free. Reading γ ITSELF as illumination is this repository's canvas and is Tier 4, not a consequence of that theorem; both readings run the same Re(λ) = -2γ⟨n_XY⟩. |
| **J** | Hamiltonian coupling strength between qubits. In star topology: J_SA (A-to-S coupling), J_SB (B-to-S coupling). |
| **h** | Transverse field strength in the Hamiltonian. |

**In plain language:**

- **γ (gamma)** is the dephasing rate: how fast a qubit loses coherence, and equivalently how much of its structure is exposed to what the dephasing charges for. This project writes the second reading, γ as light falling on the chain, and that is an interpretation laid on the algebra rather than a result of it; the algebra is the same either way. Each qubit can have its own. The concentrator formula works by putting nearly all of it on one qubit and keeping the rest dim.
- **J** is the connection strength between qubits. Stronger J means faster information exchange but also means noise on one qubit spreads faster to the others. The balance between J and γ determines everything.
- **h** is an external magnetic field that pushes each qubit individually. Think of it as a background force acting on each particle.

### The T₂ → γ conversion

This section homes a conversion the repository already derived rather than
deriving it again. The identity is
[F113_T1_EXTRACTION_KINGSTON](../experiments/F113_T1_EXTRACTION_KINGSTON.md)'s
("a `D[Z]` channel at rate γ_z decays coherences at 2γ_z, and 1/T₂ = 2γ_z +
1/(2T₁)"), and [IBM_CONCENTRATOR](../experiments/IBM_CONCENTRATOR.md)'s
2026-07-05 retro-note already priced the trap against the T₁-aware form on one
qubit (Q85, ≈ 0.046 MHz where the script had fed 0.268; that script was
repaired 2026-08-05 and now feeds the dephasing-only 0.134; two siblings on
the same Q85-Q94 line, `ibm_cavity_analysis.py` and
`cavity_mode_localization.py`, were found still feeding 1/T2* and repaired
2026-08-18, so the family is five scripts and the trap recurred once after
being logged). Note
that neither is a registry entry: **no F-number owns this conversion**, and the
F113 in that filename belongs to the document, not to the formula. (Reading it
as "F113's formula" is how a fold on 2026-08-04 mis-attributed it here.)

Why it needs a home anyway: three tracked documents call a T₂ → γ form "the
repo's convention" ([GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md),
[IBM_CONCENTRATOR](../experiments/IBM_CONCENTRATOR.md),
[IBM_RUN3_PALINDROME](../experiments/IBM_RUN3_PALINDROME.md)), all three naming
γ = 1/(2T₂), which is also what the typed calibration chain implements
(`CalibrationChain.cs`, `g = 1.0 / (2.0 * t2)`). The spread in the tree is
therefore unlabelled drift, not rival conventions: nothing that declares itself
disagrees. One number had crossed the boundary unconverted: the 70× of
[F112_HARDWARE_LENS_KINGSTON](../experiments/F112_HARDWARE_LENS_KINGSTON.md),
which divided a fitted Lindblad γ_Z by that document's own γ_eff = 1/T₂, while
the 1.72× it imports from
[IBM_BLOCK_CPSI_SATURATION](../experiments/IBM_BLOCK_CPSI_SATURATION.md) was
computed on 1/(2T₂). Two rulers, one setup. **Corrected 2026-08-05:** the denominator was **halved**,
from 1/T₂ to 1/(2T₂), which doubles the ratio and only sharpens the point it was
making. The ratio
itself is no longer quoted as a value, because its numerator turned out not to
be a rate at all: the optimiser drives that parameter up until the model
predicts no coherence, and the published figure was where double-precision
arithmetic stopped changing. The lesson the entry exists for is the one that
survives, and the ruler crossing was real: two rulers in one setup. The 1.72×
was right and is unchanged.

The one fact underneath: a `D[Z]` channel at rate γ, built as the repository
builds it everywhere (jump operator c = √γ·Z), decays coherences at **2γ**, not
γ. Everything else follows, and **which form is right depends on what else the
model carries**:

| the model's dissipator | the D[Z] rate reproducing the measured T₂ |
|:---|:---|
| Z-dephasing alone | **γ = 1/(2T₂)** |
| Z-dephasing **and** σ⁻ at 1/T₁ | **γ_Z = (1/T₂ − 1/(2T₁)) / 2** |

Both entries are **D[Z] coefficients**. The coherence rate Γ_φ = 1/T₂ − 1/(2T₁)
is *2γ_Z*, not a γ; it appears in the tree under γ-shaped names, and that is the
second factor of two in this area.

Each row reproduces the measured coherence envelope e^(−t/T₂) to machine zero
(worst deviation 5.6·10⁻¹⁷ and 1.2·10⁻¹⁶ on real Torino pairs; gate
[`t2_gamma_book_gate.py`](../simulations/t2_gamma_book_gate.py)). **That is the
whole of what they share.** They are not interchangeable models: row 1 has no T₁
at all, relaxes to the maximally mixed state, and assigns every bit of coherence
loss to dephasing. On Q85 of
[CHAIN_SELECTION_TEST](../experiments/CHAIN_SELECTION_TEST.md) (T₁ = 2.9 µs,
T₂ = 5.0 µs) amplitude damping carries most of the loss, and the two rows differ
by **7.25×** in γ; after 200 µs on Q80 they sit at ρ₀₀ = 0.500 against 0.928.
Pick by whether the question touches populations, and say which you picked.
(Only the T₁/T₂ pairs are taken from that document; its own γ column carries
row 1.)

Row 2 also **leaves the palindrome theorem**: a σ⁻ channel breaks the Π mirror,
and the repository already owns the closed form of that breakage in
[F82/F84](ANALYTICAL_FORMULAS.md) (‖D_{T1,odd}‖_F = √(Σγ²_{T1,l})·2^(N−1),
[PROOF_F82_T1_DISSIPATOR_CORRECTION](proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md)).
A residual there is F82, not a finding.

Read that sentence exactly: it says **Π** breaks, not the palindrome. By
[F137](ANALYTICAL_FORMULAS.md#f137) a σ⁻ channel on its own leaves the spectrum
palindromic, at the halved centre −Σγ/2, and a thermal bath at −Σ(γ↓+γ↑)/2.
What breaks the spectrum is σ⁻ beside **co-axial** Z-dephasing, which is row 2's
actual situation; transverse dephasing composes with it exactly. The two claims
are separate, and F137 exists to keep them apart.

Row 2 is undefined where **T₂ > 2T₁**, returning a negative rate. Read that as
the consistency check firing, not as a reason to prefer row 1: the bound is a
theorem for Markovian dynamics, so a record crossing it is a broken record, and
row 1 is silent about them only because it never reads T₁. In the Torino history
344 of 24,073 records (1.4%) cross, but 137 of those are Q53 alone, the only
qubit in the file whose T₂ is frozen (62.413 µs, identical across all 181
calibration dates, while its T₁ ranges over 18.6 to 33.7 µs), and Q53 sits in a live
experiment table. The typed layer already clamps this:
`IbmCalibration.cs` scores coherence as `min(T2, 2·T1)`.

**γ = 1/T₂ is the factor-2 trap**: exactly twice the Lindblad rate under either
row; the repository caught it in prose at least twice before repairing any
producer. When a document converts a calibration T₂ into a γ, it names its
dissipator and the row that goes with it. One more scope, since this is the
conversion's home: hardware T₂ is an echo (or Ramsey) number, and a Markovian γ
maps onto it only for white dephasing noise; under the 1/f noise transmons
actually carry, echo and Ramsey differ and neither is exactly a Lindblad T₂.

### "Machine zero" (say the number)

The phrase "machine zero"/"machine precision" spans 10⁻¹⁰ to 2×10⁻¹⁶ across ~80
files with no anchor; under the repository's no-rounding rule that range is a
category error. The rule: an EXACT route (same quantity computed twice, integer
or rational arithmetic, entry-wise rearrangement) compares **exactly**, `== 0.0`,
never "machine zero". The phrase is reserved for float-algebra routes
(eigensolvers, accumulated linear algebra), and every use states its number and
its scale (what the residual is normalized by). "Bit-exact" means bit-exact:
if a threshold accompanies it, the claim is not bit-exact.

### The coupling ratio q and Q (and the factor-2 convention)

The dimensionless ratio of coupling to dephasing, J/γ₀, is the central control parameter (the "balance between J and γ" above). It is written two ways, and they denote the **same physical ratio** in **two Hamiltonian normalizations that differ by a factor of 2**. Stating it once, here, because both of us have stumbled on it:

| Symbol | Hamiltonian per bond | Written by |
|--------|----------------------|------------|
| **q** | H = J·(XX+YY), hop element 2J (the "doubled" book) | the F89 octic / branch-locus / monodromy / Galois work (`F89Path3OcticBlock`) |
| **Q** | H = (J/2)·(XX+YY), hop element J | the carrier clock (θ = arctan(Q·cos(π/(N+1))), `ClockHandLadderClaim`), the coherence horizon Q*(N), the F86 c=2 block, the Stage-1a probe (`ChainSystem` / `Symphony`); also, by value, the unit-hop pencil of the seed-existence nullity note (its knob is *named* q there; the pin in `experiments/F89_SEED_EXISTENCE_REDUCTION.md` converts) |

The relation is **q = Q/2** (equivalently J_F89 = J_F86/2; the F90 bridge relabeling in `F90F86C2BridgeIdentity`, where the two are operator-exactly one Liouvillian: ‖L_F86(J) − L_F89(J/2)‖ = 0 at N=5..8). The factor is pure relabeling, not physics. On the real axis, within one convention, the complex variable q is just the analytic continuation of the real Q (`experiments/F89_BRANCH_LOCUS_PALINDROME.md`, the line under the figure).

**Worked example, so the trap is concrete:** the octic diabolic degeneracy sits at q_EP ≈ 0.659 in octic units (`F89Path3OcticEpClaim`). In carrier-clock units that is 2 × 0.659 = 1.318, i.e. between the coherence-horizon rungs Q*(2) = 1 and Q*(3) = √2, **not** below 1. Convert (halve a Q, or double a q) before putting an octic q and a horizon Q on one axis.

**Two thresholds on this axis, not one:** the coherence horizon **Q\*(N)** (1, √2, 1.8787, 2.3737 in carrier-clock units; the single-excitation EP where the longest-lived mode changes character, `CoherenceHorizonClaim`) and the gap threshold **Q\*_gap(N)** (the coupling below which the spectral gap drops under 2γ: a |Δpopcount| = 0 population mode is Zeno-suppressed below the 2γ floor, while at **uniform γ** every coherence sector keeps its minimum at or above its rung at every coupling, [the Absorption Theorem §4.3](proofs/PROOF_ABSORPTION_THEOREM.md), which derives that floor and measures its attainment separately. Under a per-site rate profile neither the value nor the ladder carries over unchanged, and what does is in [`band_edge_profile_fence_gate.py`](../simulations/band_edge_profile_fence_gate.py)). Q\*_gap is bisected in **Pauli-J units** (H = J·Σ(XX+YY+ZZ)): Heisenberg chain 0.5000 / 0.8002 / 1.3422 / 1.8194 at N = 2..5, XY chain 0.5000 / 0.7071 / 0.9393 / 1.1861, so it belongs to the Hamiltonian, not to N alone. Converted to carrier units (×2), the XY gap threshold matches the horizon to bisection precision at N = 2 and N = 3 (1.000000 and 1.414214, i.e. 1 and √2 to six places) and separates measurably from N = 4 on (1.878541 vs 1.87874, then 2.372175 vs 2.37367). The coincidence at low N is where both reduce to the same clean 2×2 block; no proof that they must agree there is on file, so treat the N = 2, 3 agreement as observed, not as an identity. Sources: [D06](proofs/derivations/D06_SPECTRAL_GAP.md), [Absorption Theorem proof §4.3](proofs/PROOF_ABSORPTION_THEOREM.md), [`absorption_ladder_regimes.py`](../simulations/absorption_ladder_regimes.py).

**Two more uses of the same letter, not this ratio:** lowercase q is also the right-popcount index in a "(p,q)-block" (the coherences between popcount-p and popcount-q basis states, an integer sector label). And the literature's quality factor Q is a different quantity again ("different objects, same letter").

**In plain language:** q and Q both ask how much turning happens in one γ₀-tick. One book (the F89 octic) counts the coupling at double strength, so its number is half. Halve a Q, or double a q, before you compare the two.

**γ₀ is the unit, and a substance does not hand over a Q.** The Liouvillian factors exactly as L(J, γ₀) = γ₀·L₁(Q), so only Q carries physics and the absolute rate cancels ([`Q_SCALE_THREE_BANDS.md:83`](../experiments/Q_SCALE_THREE_BANDS.md), [`THE_GENESIS_OF_AN_OSCILLATION.md:53`](THE_GENESIS_OF_AN_OSCILLATION.md)). In the dephasing-only book above, γ = 1/(2T₂), so in terms of a measured time and an energy **Q = 2·J·T₂/ℏ** (watch the 2: it is the Lindblad book, and it is absent in the coherence book). Naming a Q therefore requires a chosen two-level degree of freedom, a chosen coupling, and a chosen channel, which is why "the Q of water" is under-determined: [Q Belongs to No Substance](Q_BELONGS_TO_NO_SUBSTANCE.md), which audits every substrate Q in the repository on that basis. The anchors on the axis are in [Q-Regime Anchor Map](Q_REGIME_ANCHORS.md).

---

## Star topology (S + A + B)

The star topology is the simplest non-trivial network: three qubits
where one (S) connects to the other two (A, B), but A and B are not
directly connected. Think of two people who only communicate through
a mutual friend.

| Symbol | Definition |
|--------|------------|
| **S** | Central qubit (shared object), coupled to both A and B |
| **A, B** | Observer qubits, not directly coupled to each other |
| **CΨ_AB** | CΨ for the AB reduced pair (traced over S). Crossing 1/4 means the observer-observer connection is active. |

### Three conditions for AB crossing (at γ = 0.05)

1. **Strong sender:** J_SB/J_SA >= 1.466
2. **Quiet receiver:** γ_A < ~0.20
3. **Deep pre-existing connection:** Initial C_SA > 0.8 (Bell-like)

---

## Observer distinction (Lindblad model)

| Term | Technical definition | Status |
|------|---------------------|--------|
| **C_int** | Symmetric Hamiltonian coupling (mutual interaction) | Tier 2 formalization |
| **C_ext** | Projective measurement (one-directional intervention) | Tier 2 formalization |

Note: The original claim that C_int preserves coherence 33x longer than C_ext has been disproven (see MATHEMATICAL_FINDINGS Section 9). The formal Lindblad distinction is verified; the physical claim about different "kinds of observation" is not.

---

## Perspectival Time Field (PTF) and the label layer

| Term | Technical definition | Status |
|------|---------------------|--------|
| **PTF (Perspectival Time Field)** | Each perspective on the one dynamics (each site, each observable) paints it at its own time rate α_i; no perspective is privileged, and the whole is the closure of the perspectives, not a view from nowhere. | Tier 2 computed |
| **Label / canvas** | A name is a canvas painted true at a stance (a vantage at a time); calculations are recomputed at every use, labels are inherited and recomputed at none, so errors accumulate in the names while the formulas keep closing. | Tier 4 reading |

Canonical: [Perspectival Time Field](../hypotheses/PERSPECTIVAL_TIME_FIELD.md) (first named "Site-Local Time"; the rename is the in-house precedent that a label can be falsified while every number under it stands). The label theory chapter is [Labels Translated](quantum/LABELS_TRANSLATED.md); the translation series lives in `docs/quantum/`; every label correction the repository has made is assembled in [The Label Map](quantum/THE_LABEL_MAP.md). The exact core (the dephasing routes by the letter it holds: the price list −2γ·n_anti is letter-routed) is typed as `HeldLetterRoutingClaim` (named `WatchedLetterRoutingClaim` until 2026-08-09, when the sender arc moved the eye out of the typed Tier-1 strings and the API name; the Tier-4 reading keeps its watcher, where the watcher is the subject; the founding chapter moved to the sending on 2026-08-09), recomputed live at `inspect --root label`.

---

## The record laws (F135/F136, July 2026)

The two laws' own object, defined here in the frame the algebra uses (added
2026-08-08, the γ-arc's Q1: the store whose job is defining had entries for the
"Watcher factor" and the "Record radius" but never for the record itself).
Canonical carriers: [PROOF_RECORD_PARITY_LAW](proofs/PROOF_RECORD_PARITY_LAW.md),
[PROOF_RECORD_LETTER_LAW](proofs/PROOF_RECORD_LETTER_LAW.md), the F135/F136
entries of [ANALYTICAL_FORMULAS](ANALYTICAL_FORMULAS.md); typed as
`RecordParityLawClaim` + `RecordLetterLawClaim`, recomputed live at
`inspect --root record`.

| Term | Technical definition | Plain reading |
|------|---------------------|---------------|
| **Record** | The separation the arriving light has effected between two conditional faces of a pair page ρ_Sj, held in the state itself, in exactly two 1-bit forms: the pointer form separates the Z_S faces, the Bell form the anti-pointer X_S/Y_S ones. Pointer record: classical-quantum, distinguishability D_j(t) = β_j(t)·\|sin θ_w(t)\| (a trace distance; β_j is F135's dressed record radius, see that entry below), the bit written in j's equator (the ZY channel), blind to γ_S. Bell record: a rank-2 Bell mixture ½Φ+½Ψ of magnitude κ = e^{−2(γ_S+γ_j)t*}, zero pointer content (D_j = 0), letter XX or YY by dresser count m mod 2, paying both sites. No load-bearing step of either law invokes a reader: the forced verdicts are integer-parity arithmetic on the coupling ratios r_k = Δ_jk/Δ_S plus graph conditions (the write bond and deg(S) ≥ 2 for the pointer; Q = ∅ and m ≥ 1 for Bell), evaluated at the readout point t* = π/(4Δ_S) with every S-bond at the one write coupling Δ_S, which at uniform coupling is the graph-state instant. | A difference the arriving light has separated into a neighbor's page, in one of two shapes. Nothing in the formalism reads it, and it does not need a reader to exist. |
| **Held or not** | At integer ratios the values are forced: a full bit (pointer or Bell) or an exact zero, dark (a mixed-parity D, an odd private watcher, or a nonempty Q for the Bell form, each a cosine of an arriving angle landing on a zero; or a missing write bond, the sine that never leaves zero, no writer; all γ-independent). Non-integer ratios give the free in-between, 0 < I < 1: parities forced, in-betweens generic, with one gated exception at the fully-shared m = 1 corner, where the bit stays forced at every ratio and only its channel rotates. At γ = 0 nothing is ever lost, the cosines turn and recur, and an even-parity return is whole, at worst rotated by π. The only true loss is γ on the pair's own two sites (the Bell record pays both, the pointer only the witness; every traced site's γ is exactly invisible). | Forced or free: the light separates the difference fully, or not at all, or partway between the angles. The turning is reversible; only the page's own standing in light costs. |
| **Watcher (graph role)** | N(j)∖{S}, the witness's other neighbors, read off the bond list; not an agent, and distinct from the hardware "Watcher factor" η(r) below. F136 splits the set: the shared dressers D (adjacent to both S and j) select the family by parity (all even with the write bond → pointer; all odd with Q = ∅, m ≥ 1 → Bell; mixed → dark), and the private watchers P must be even for either record to survive (one odd private watcher kills both). F135's even-forgiven / odd-blind / else-generic trichotomy is the D = ∅ column, under Law A's hypotheses (uniform Δ_S at S, deg(S) ≥ 2, γ = 0, triangle-free at S). | The neighbors whose angles decide the parity: the shared ones pick the record's shape, the private ones decide whether it survives. |

The reader enters only downstream of the laws: Zurek's redundancy R_δ (how many
fragments a reader could consult) and the Holevo bound are imported bookkeeping
on top of the closed form, and the one design that ever implemented a conditional readout is
the unflown record-parity flight, which itself refuses to promote mutual
information to a measured claim. When a reader steps forward, the translation
from qubit to observer is allowed by the algebra, not forced by it
([Observer Inheritance](../reflections/OBSERVER_INHERITANCE.md)).

---

## The blind seat (F157, August 2026)

A **seat** is a site carrying the watching (single-site Z-dephasing); the
**blind** subspace at seat j is the pure single-excitation states that
dephasing cannot touch, counted exactly by F157: blind(j) = N − dim Krylov(e_j), the
Krylov space being span{e_j, He_j, H²e_j, ...} with e_j the seat's unit vector
(one fenced case: an isolated seat adds its own H-invariant ray, F157's
Breaks-for); on the uniform open chain (gcd(2j+1, N) − 1)/2 with the ZZ term
and gcd(j+1, N+1) − 1 without. The **span** beside it is `dim ker L_SE(j)`, what
the arc calls the dimension of the watched sector's stationary manifold (the
manifold of stationary STATES inside it is one dimension smaller), which is
`1 + dim commutant(H restricted to the Krylov complement)` and therefore
`1 + blind(j)` exactly when that restriction has a simple spectrum, a condition
on the Krylov complement and NOT on the blind subspace, which can be degenerate
where the identity holds ([the span
and node-lemma proof](proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md)); a
zero-free chain always supplies that, and it is not the same word as the span
of a sweep window in [The Blind Site](../experiments/THE_BLIND_SITE.md) §7
([The Blind Site](../experiments/THE_BLIND_SITE.md),
[The Seat That Cuts](../experiments/THE_SEAT_THAT_CUTS.md)). Not to be
confused with the repo's other blindnesses, which are different objects:
**J-blind** receiver states
([J-blind receiver classes](../experiments/J_BLIND_RECEIVER_CLASSES.md),
read beside F73's closure), **graph-blind** derivations
(F153, a derivation using no graph structure), the local page's kinematic
blindness to distant coherence (F70), and the light/lens letter blindness
above. MirrorWorld's `BlindSeat.cs` fences its own words at the door.

---

## Hardware pre-registration vocabulary (the IBM flights)

Terms the flight pre-registrations use as house standard (added 2026-08-04, the
record-parity round-31 prior-work review: five of them had no anchor outside their
own file). Canonical carrier:
[RECORD_PARITY_HARDWARE_PREDICTION](../experiments/RECORD_PARITY_HARDWARE_PREDICTION.md);
the two-gate architecture originates in
[IBM_F129_RAMSEY_FRINGE](../experiments/IBM_F129_RAMSEY_FRINGE.md).

| Term | Technical definition | Plain reading |
|------|---------------------|---------------|
| **Band** | A frozen numeric constant of a pre-registration (threshold, floor, or window edge), derived by a pinned rule from simulation banks BEFORE any hardware data exists and committed with the design; bands never change after the commit hash. | The lines on the referee's card, drawn and signed before the match. |
| **VOID** | The verdict class "the instrument was out of spec; no physics verdict is claimed"; takes precedence over every physics verdict, and statistical-guard VOIDs are terminal on a paid one-shot. | The match is abandoned, not lost or won. |
| **Double ratio ρ̂(r)** | [Ŝ(j; r)/D̂(j′; r)] / [Ŝ(j; 0)/D̂(j′; 0)]: the record-parity flight's verdict carrier: the inner ratio cancels slow common drift, the r = 0 division cancels static asymmetry. | Compare the watched witness to its control, then compare that to the unwatched day. |
| **Watcher factor η(r)** | The one-sided device attenuation dressing the predicted cos(r·π/2): incoherent-only by pin; η_nom = its Aer-simulated estimate on the representative line, η_min = the worst-admitted bound, η structurally may exceed 1 (the j′-differential). | How much the machine dims the signal, measured on the machine's own twin. |
| **7a / 7b** | The two-gate standard: 7a = the sim gate (exact models through the frozen verdict code, bands and controls measured from below); 7b = the counts-level gate (the actual runner's circuits and estimator on synthetic counts, bit-for-bit against 7a; contradiction = hard abort). | First prove the plan on paper, then prove the actual instrument replays the paper exactly. |
| **Record radius** | Two non-identical formulas share the name: F135's dressed β_j = e^{−2γ_j t*}·Π\|cos(·)\| and F136's bare parity form (no γ factor); a document using the name states which (the record-parity flight uses F135's, at γ = 0, where they coincide). | Say whose ruler you are holding. |

---

## Boundary values

These are specific numbers that appear repeatedly across the project.
They are not arbitrary; each one comes from the mathematics.

| Value | Source |
|-------|--------|
| **1/4** | Discriminant of the fixed-point equation. Below 1/4: real fixed points exist. Above: complex. Algebraically exact within the iteration. |
| **~1.466** | J_SB/J_SA threshold for AB crossing at γ = 0.05 (star topology). |
| **1/3** | CΨ of a maximally entangled Bell pair (C = 1, Ψ = 1/3). |
| **2γ** | Decay rate of the c− supermode (antisymmetric, slow). At **uniform γ** the 2γ rung exists at every coupling; it is the *slowest* nonzero rate only above the threshold Q*_gap(N). Under a per-site rate profile there is no single 2γ, and what replaces it is a floor rather than a value; on the anti-palindromic locus a ladder does return, with γ̄ in place of γ, under conditions the docstring of [`band_edge_profile_fence_gate.py`](../simulations/band_edge_profile_fence_gate.py) states. |
| **8γ/3** | Concurrence-envelope decay rate at N=3, in the J/γ → ∞ limit; at finite coupling a three-level band (F33). |
| **10γ/3** | Decay rate of the c+ supermode (symmetric, fast) at N=3, same J/γ → ∞ limit as 8γ/3 (F33). |
| **-2Σγ** | Location of XOR modes. The value is a genuine ceiling for any Hermitian H (Bendixson), but the COUNT belongs to the number-conserving XY/Heisenberg family and not to Z-dephasing as such: N+1 modes there, 2^N under a pure Ising ZZ chain, and NONE at all under a generic Hermitian H, whose maximum rate falls short of 2Σγ entirely (measured at N = 3, γ = 0.05: 4, 8 and 0 modes, generic maximum 0.2112 against 2Σγ = 0.30). |
| **0.886** | Best average fidelity for QST (star topology, J_SB/J_SA = 2:1, γ = 0.05). |
| **0.036/γ** | Crossing time t_cross for Bell+ under dephasing **in the concurrence book**, i.e. reading C as the Wootters concurrence. Exact there: K = ln(4/3)/8 = 0.03596 (write 0.03596 when the word "exact" is attached; 0.036 is the rounded label). Three K coexist for this one state and channel, and the discriminator is WHICH C: C = Wootters concurrence gives 0.03596; C = purity, the framework's own CΨ = f(1+f²)/6, gives K = 0.03735 ([F25](ANALYTICAL_FORMULAS.md), Tier 1 proven, and the t\* = 0.747 cockpit landing); the February tool's feedback model on the concurrence book gives 0.0387, historically quoted as 0.039. See CROSSING_TAXONOMY. |
| **0.976** | Correlation between mixed XY Pauli weight and XOR fraction (N≥3). |
| **360×** | Concentrator formula vs V-shape at N=5 (peak created Sum-MI, a transport metric, ε→0 sim; C# RK4 validated; ~2-3× hardware). |
| **180×** | Concentrator formula vs V-shape at N=7 (same transport metric). |
| **139×** | Concentrator formula vs V-shape at N=9 (same transport metric; the factor declines to 68× by N=15). |
| **15.5 bits** | Channel capacity of spatial dephasing profile at 1% noise (γ as signal). |

---

## Epistemic tier labels

This project uses a five-level system to mark how confident we are
in each claim. Whenever you see a tier label in the repository, it
means exactly this:

| Tier | Meaning |
|------|---------|
| **1** | Mathematically proven (analytical derivation, no free parameters) |
| **2** | Numerically verified or hardware confirmed (simulation, IBM data, C. elegans) |
| **3** | Consistent with data but not uniquely proven (plausible, not exclusive) |
| **4** | Plausible hypothesis, needs more work (motivated by computation but untested) |
| **5** | Speculation, interpretation, philosophy (not falsifiable or not yet testable) |

---

*See [The CΨ Lens](THE_CPSI_LENS.md) for the canonical project description.*
*See [Core Algebra](historical/CORE_ALGEBRA.md) for the proven mathematics.*
*See [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md) for the palindrome theorem.*
*See [Non-Heisenberg Palindrome](../experiments/NON_HEISENBERG_PALINDROME.md) for the extended palindrome analysis.*
*See [XOR Space](../experiments/XOR_SPACE.md) for the spectral filter discovery.*
*See [Resonant Return](../experiments/RESONANT_RETURN.md) for the concentrator formula.*
