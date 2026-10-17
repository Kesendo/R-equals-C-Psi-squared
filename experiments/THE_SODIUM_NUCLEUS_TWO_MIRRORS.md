# The sodium nucleus: two mirrors in the quadrupolar relaxation of a spin 3/2

**Date:** 2026-10-17
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5.5)
**Verifier:** [`sodium_two_mirrors.py`](../simulations/sodium_two_mirrors.py) (`python simulations/sodium_two_mirrors.py`; exact, sympy 1.14, every parameter a symbol where it can be; 36 checks)
**Rates, not ours:** A. M. Souza, A. Gavini-Viana, I. S. Oliveira, R. S. Sarthour, R. Auccaise, J. Teles, E. R. deAzevedo, T. J. Bonagamba, *Nuclear spin 3/2 electric quadrupole relaxation as a quantum computation*, [arXiv:1105.2539](https://arxiv.org/abs/1105.2539), Eqs. (15) to (27), which cite R. Auccaise et al., J. Magn. Reson. 192, 17 (2008); the spectral densities follow J. R. C. van der Maarel, J. Chem. Phys. 94, 4765 (1991), and G. Jaccard, S. Wimperis, G. Bodenhausen, J. Chem. Phys. 85, 6282 (1986).
**A word fence:** *order* p = m − m′ is the coherence order of |m⟩⟨m′|; J₀, J₁, J₂ are the spectral densities of the electric field gradient at 0, ω₀ and 2ω₀ as they enter the secular Redfield generator; *dark* here means untouched by the J₀ dissipator, an operator |m⟩⟨m′| that commutes with T₂₀: the operator sense of [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)'s DARK (commuting with every jump), restricted to the one jump T₂₀, and not a state annihilated by a jump (the singlet page's sense); *bright* the rest; a *mirror about −c* is the map λ ↦ −λ̄ − 2c on the eigenvalues of one order block, and c is called its centre. Inside a block every eigenvalue carries the block's frequency with one sign, so the repository's F1 form λ ↦ −λ − 2σ would send it to the block of order −p; since the order −p block has the complex-conjugate characteristic polynomial of the order p block (checked), the map used here is F1's map followed by that conjugation, read inside one block.

**Conventions.** Spin operators I_z, I_± for I = 3/2; T₂₀ = (3I_z² − I(I + 1))/√6, T₂,±₁ = ∓(I_zI_± + I_±I_z)/2, T₂,±₂ = I_±²/2; Lindblad dissipator Σ_m J_|m|·(T₂ₘρT₂ₘ† − ½{T₂ₘ†T₂ₘ, ρ}); Hamiltonian ω_Q·T₂₀ + ω·I_z, so ω_Q is the coefficient of T₂₀ in this normalization, not the conventional quadrupolar splitting. With these conventions Souza et al.'s rates are reproduced with their constant set to C = 3; that is a normalization of the J here, not a value of their physical constant. Dynamic frequency shifts (the imaginary parts of the spectral densities) and the finite polarization of the equilibrium state are not in the model.

## What this is about

Sodium in tissue is a spin 3/2: four levels, not two. It relaxes through the fluctuating electric field gradient around it, not through a magnetic field, and sodium MRI uses that relaxation to tell tissue compartments apart, most directly through the triple-quantum coherence between m = +3/2 and m = −3/2. The repository's mirror results live at d = 2, and for d > 2 they hold only partially ([F121](../docs/ANALYTICAL_FORMULAS.md), under equidistant dephasing); [Qubit Necessity](../docs/QUBIT_NECESSITY.md) §8a leaves every other dissipator at d > 2 open. This page reads one physical one, the quadrupolar dissipator of a single spin 3/2.

When all three spectral densities are present it has no single mirror in the cases checked. It has two. The slow part of the gradient, J₀, sees one property of a level: outer (|m| = 3/2) or inner (|m| = 1/2). Everything that agrees on that property, the populations, the central transition and the triple quantum, is dark to it, and this dark half of the spectrum is palindromic about −3(J₁ + J₂), with the central transition and the triple quantum sitting exactly at its centre. The bright half, the satellite and double-quantum coherences, is palindromic about a second centre 3J₀ further out. Both mirrors are exact, for every static quadrupolar splitting and every offset. The dark centre is the mean of the two longitudinal rates 6J₁ and 6J₂, and it is also the triple-quantum rate; that identity was already a line away from the published rates, and the two mirrors give it a place.

## Abstract

One spin 3/2 under the secular Redfield generator of quadrupolar relaxation at infinite temperature (jumps T₂,ₘ, m = −2 … 2, at rates J₀, J₁, J₁, J₂, J₂), with a static splitting ω_Q·T₂₀ and an offset ω·I_z. The generator is block diagonal in the coherence order. Exactly, with every parameter a symbol: the populations decay at 0, 6J₁, 6J₂, 6(J₁ + J₂); the central transition and the triple quantum at 3(J₁ + J₂), shifted only by their frequencies; at ω_Q = 0 the satellites at 3(J₀ + J₁) and 3(J₀ + J₁) + 6J₂, the double quantum at 3(J₀ + J₂) and 3(J₀ + J₂) + 6J₁ (Souza et al.'s rates with their constant set to 3); a splitting draws each pair together along the real axis, and past 6ω_Q² = 9J₂² (satellites) or 6ω_Q² = 9J₁² (double quantum) turns it into a conjugate pair on the same real part. The dark set pairs under λ ↦ −λ̄ − 2c₁ with c₁ = 3(J₁ + J₂), the bright set under the same map with c₂ = 3(J₀ + J₁ + J₂) = c₁ + 3J₀; with ω_Q and ω present, the satellite and double-quantum polynomials keep real coefficients after the frequency shift p·ω and the root sum −2c₂. The two centres coincide exactly when J₀ = 0. At ω_Q = ω = 0, in every case checked, the whole spectrum has one centre exactly when one of J₀, J₁, J₂ vanishes; that rule is not proved.

## The sweep: what the repository held

Swept 2026-10-17 by one survey agent, stores by name.

- `docs/ANALYTICAL_FORMULAS.md`: [F121](../docs/ANALYTICAL_FORMULAS.md), the qudit partial palindrome under full-Cartan (equidistant) dephasing; at d = 4, N = 1 it counts 8 of 16 modes paired about one centre. No entry on quadrupolar relaxation, spectral densities J₀/J₁/J₂, Redfield, ²³Na, triple quantum or central transition.
- `docs/proofs/`: [PROOF_QUDIT_PARTIAL_PALINDROME](../docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md), the equidistant ladder, the ceiling, the product cap and its clock-shift mirror, checked at d = 4, N = 1, 2. Nothing on spin 3/2 as an angular momentum or on non-equidistant levels.
- `docs/`: [Qubit Necessity](../docs/QUBIT_NECESSITY.md), the trunk d² − 2d = 0, and §8a: 0 of 236 sampled qutrit dissipators paired fully, non-dephasing dissipators at d > 2 left open.
- The typed layer: `QuditPartialPalindromeCeiling`, `QuditProductMirrorCap`, `QuditMirrorProtectionScalingClaim` (Tier 1) and the live `QuditPartialPalindromeWitness` (`inspect --root qudit`), all on the equidistant model. MirrorWorld: `Formulas.F121_*`, d = 2 and 3 only.
- `experiments/`, the OpenArcs registry, `docs/GLOSSARY.md`, `docs/CAUGHT_ERRORS.md`, `fw.Confirmations`: nothing on spin > 1/2 physics or quadrupolar relaxation; CAUGHT_ERRORS holds the qudit counting corrections (the product cap holds for d ≤ 5, so at d = 4).

Checked for adjacency: F121 and this page are about the same Hilbert space and different dissipators; F121's single-centre count (8 of 16) and this page's two-centre pairing (16 of 16) do not contradict each other. [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md) found the same shape for qubit chains under correlated Z noise: the global mirror breaks, mirrors survive block by block about centres that differ by what the watching charges. [The Shared Field in a Crystal](THE_SHARED_FIELD_IN_A_CRYSTAL.md) is the companion page on multiple-quantum NMR.

**Prior art.** The rates are Souza et al.'s Eqs. (15) to (27), and the reading of the four levels as two qubits, |3/2⟩, |1/2⟩, |−1/2⟩, |−3/2⟩ as |00⟩, |01⟩, |10⟩, |11⟩, is theirs too, with a quantum-circuit model of the same relaxation; in their labelling the outer/inner property below is the parity Z⊗Z. They also state that the triple quantum and the central transition are the only coherences free of J₀, a decoherence-free subspace of the J₀ channel spanned by Bell states: this page's dark coherences are theirs. That the two share the slow rate J₁ + J₂, that the longitudinal relaxation is biexponential, and that the satellites decay at the bright centre's rate, is the standard picture behind triple-quantum-filtered sodium NMR; both centres are known rates. What this page adds is reading them as the centres of two exact palindromes, and the global-centre rule below; the literature was not searched for either, and the mean-rate identity follows in one line from the cited rates.

## The prediction and how it fared

**Written before computing:** with J₀ alone the spectrum is an exact palindrome, since T₂₀ is a Z on one bit (held: two levels of rate, 0 and 3J₀, eight modes each); with J₁, J₂ present the palindrome breaks in general, and a guess that extreme narrowing (J₀ = J₁ = J₂) restores it. The guess failed: at extreme narrowing the two centres sit at 6J and 9J and the global mirror breaks. What survives is the block structure below.

## The one bit the slow field reads

T₂₀ = (3I_z² − 15/4)/√6 is diagonal, (3, −3, −3, 3)/√6 on m = 3/2, 1/2, −1/2, −3/2. It sees one thing about a level, whether it is outer or inner. Label each level by two bits, a = outer/inner and b = the sign of m; then T₂₀ = (3/√6)·Z_a, one Z on bit a (in Souza et al.'s labelling the same operator is Z⊗Z). Its dissipator charges |m⟩⟨m′| at J₀·½·(t_m − t_m′)², with t the diagonal entries: zero when m and m′ agree on bit a, ½·(6/√6)²·J₀ = 3J₀ when they disagree. That is the repository's Hamming price, on one bit of two.

| coherence | order p | outer/inner agree? | J₀ cost |
|---|---|---|---|
| populations | 0 | yes | 0 |
| ±1/2 central transition | 1 | yes (both inner) | 0 |
| ±3/2 triple quantum | 3 | yes (both outer) | 0 |
| satellites, ±3/2 ↔ ±1/2 | 1 | no | 3J₀ |
| double quantum | 2 | no | 3J₀ |

## The two mirrors

The verifier computes the characteristic polynomial of every order block with J₀, J₁, J₂, ω_Q, ω as symbols; the blocks of order −p have the conjugate polynomials of those of order p, so the statements below for p ≥ 0 carry over.

- **Dark, centre c₁ = 3(J₁ + J₂).** The populations have roots 0, −6J₁, −6J₂, −6(J₁ + J₂) at every ω_Q and ω, two pairs about −c₁. The central transition is the factor x + c₁ + iω of the order-1 polynomial and the triple quantum the whole order-3 polynomial x + c₁ + 3iω: both sit at Re λ = −c₁, each its own mirror image.
- **Bright, centre c₂ = 3(J₀ + J₁ + J₂).** The remaining order-1 quadratic and the order-2 quadratic, written in y = x + p·iω, have real coefficients and the linear coefficient 6(J₀ + J₁ + J₂). Their roots are either a real pair about −c₂ or a complex pair on Re = −c₂ (the static splitting enters only the constant term, as 6ω_Q²), and in both cases λ ↦ −λ̄ − 2c₂ maps the pair to itself.
- **The distance between the centres** is c₂ − c₁ = 3J₀, the J₀ cost of a bright cell. The two centres coincide exactly when J₀ = 0, so then the whole spectrum has one centre, for every ω_Q and ω.
- **One centre for the whole spectrum otherwise.** At ω_Q = ω = 0 the verifier decides it exactly at rational rates: a global centre exists for J₀ alone, J₁ alone, J₁ with J₂, J₀ with J₁ and J₀ with J₂, pairing across blocks, and does not exist at extreme narrowing or at three other points with all three densities positive; 24 further random rational points (outside the verifier) follow the same rule, a global centre exactly when one density vanishes. That rule is not proved, and special ratios with all three positive are not excluded. A static splitting breaks the cross-block pairing even for J₀ alone.

## What a sodium measurement can read

Under this generator the triple-quantum rate, equal to the slow transverse rate, is the mean of the two biexponential longitudinal rates: c₁ = (6J₁ + 6J₂)/2. In Souza et al.'s units (rates C·J) the longitudinal rates are 2CJ₁ and 2CJ₂ (the third population rate, 2C(J₁ + J₂), is their sum and not a separate longitudinal rate), and the slow transverse and triple-quantum rate is C(J₁ + J₂). It is a property of the model, not a new measurement: a relaxation study that reports these rates of one sample can check it, and a deviation beyond the measurement error and the fit bias says that something outside the model is present (more than one sodium environment in the voxel, dynamic frequency shifts, motion that the secular Redfield form does not capture).

## What this does not show

No sodium data were compared. The model is a single isolated spin 3/2 at infinite temperature in the secular Redfield picture; multiple compartments, exchange, finite polarization and non-secular terms are outside it. Nothing here improves a sodium image. The mean-rate identity follows directly from rates published before this page; only its place as a mirror centre is ours.

## Open

- The identity against measured sodium relaxation rates (gels, cartilage, brain), whether it is already in use as a consistency check, and which deviation a two-compartment voxel produces.
- The global-centre rule (one centre exactly when one spectral density vanishes): a proof, or a ratio with all three positive that breaks it.
- Whether the two-centre structure extends to spin 5/2 and 7/2 (more bits, more centres), and to an ensemble of spins with a shared gradient.
- The general d > 2 question of Qubit Necessity §8a: which other physical dissipators split a palindrome into centres instead of breaking it.

## Links

- [F121](../docs/ANALYTICAL_FORMULAS.md) and [PROOF_QUDIT_PARTIAL_PALINDROME](../docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md), the equidistant d > 2 case.
- [Qubit Necessity](../docs/QUBIT_NECESSITY.md), the d² − 2d = 0 trunk and the open d > 2 dissipators.
- [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md), blockwise centres for qubit chains.
- [The Shared Field in a Crystal](THE_SHARED_FIELD_IN_A_CRYSTAL.md), the multiple-quantum companion page.
