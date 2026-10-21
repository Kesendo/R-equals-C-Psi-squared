# The sodium chain: the two mirrors under a coupling between nuclei

**Date:** 2026-10-21
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 5.5)
**Verifier:** [`sodium_chain_mirrors.py`](../simulations/sodium_chain_mirrors.py) (`python simulations/sodium_chain_mirrors.py`, 56 checks at two and three nuclei, about a minute and a half; `--n3` adds the sharp subspaces at three nuclei, 60 checks; integer and Gaussian-integer arithmetic where a statement is exact, two primes where it is read modulo p)
**Parent page:** [The Sodium Nucleus: Two Mirrors](THE_SODIUM_NUCLEUS_TWO_MIRRORS.md), whose conventions, labelling and word fence (*dark*, *bright*) are used here unchanged. One difference: the parent reads a mirror inside one coherence order as λ ↦ −λ̄ − 2c; here a mirror is the map λ ↦ −λ − 2s on a whole invariant part, the form of the parent's paragraph "One map for both mirrors" in its section on the Dirac matrices.

## What this is about

One sodium nucleus, a spin 3/2, relaxes with two mirrors: the half of its modes that the slow component J₀ of the field-gradient fluctuations cannot see is palindromic about one centre, the other half about a second centre 3J₀ further out, and both are one map. This page asks what a static coupling between nuclei does to the two mirrors, for the forms such a coupling takes: the Ising part I_z·I_z, the secular dipolar coupling, and its truncations when the quadrupolar splitting dominates it.

The answer has three layers. Whenever the coupling keeps each nucleus's outer/inner bit, the chain carries a ladder of exact mirrors, one centre for each number of bright nuclei. When two neighbours have the same quadrupolar splitting, their satellites exchange polarization resonantly (satellite-transition spin diffusion), brightness moves from one nucleus to the next, and the ladder is exact only on the part of the spectrum where brightness stays sharp. The full secular dipolar coupling, untruncated, breaks the map itself, at exactly the terms that couple a satellite of one nucleus to the central transition of the other.

The physical regime is an assumption, not established here: a static coupling between specific nuclei together with the motionally narrowed relaxation of the parent page. Where molecular motion is fast enough for that relaxation, it usually averages the coupling between nuclei; ordered media, where residual couplings survive the motion, are the nearest candidate. No sample is known to this page in which both hold.

## Abstract

Spins 3/2, each with its own field gradient (independent quadrupolar relaxation, rates J₀,k, J₁,k, J₂,k per nucleus), a static splitting and an offset, coupled by a real two-body Hamiltonian on any graph. In the parent page's two-qubit reading each nucleus is a qubit pair, its relaxation is dephasing by the strings ZZ, ZX, ZY, XI, YI, and S(ρ) = (⊗_k ZZ_k)·ρᵀ. **Lemma** (any number of nuclei, any graph): if H is real and commutes with ⊗_k ZZ_k, then S·L(J₀)·S⁻¹ = −L(−J₀) − 6·Σ_k(J₁,k + J₂,k): the mirror image of the chain is the chain with the J₀ rates reversed in sign. **Ladder**: if H commutes with every ZZ_k, the part whose bright nuclei form the set B is palindromic about −s_B, s_B = 3·Σ_k(J₁,k + J₂,k) + 3·Σ_{k∈B} J₀,k, up to 2ⁿ centres; for identical rates n + 1, s_m = 3n(J₁ + J₂) + 3m·J₀ with m bright nuclei, the parent's two at n = 1. The Ising coupling and the secular dipolar coupling between nuclei of unequal, non-zero splittings both qualify. **Equal splittings:** the lemma holds and the ladder is exact on K_m, the largest invariant subspace with m bright nuclei, and not forced off it; for identical nuclei on an open chain with one coupling, K_m has dimensions at most 24 + 128 + 24 = 176 of 256 at two nuclei and 64 + 140 + 84 + 16 = 304 of 4096 at three (upper bounds over two primes at the verifier's parameters, and so at generic ones; at two nuclei the same at two rate points). At J₀ = 0 the whole spectrum is palindromic. **The full secular dipolar coupling** breaks S at its eight √3 terms, and at J₀ = 0 the global mirror (about −σ, the only candidate centre) breaks exactly, by F158 at every length and again by its characteristic polynomial at two nuclei. At two nuclei its K₁ and K₂ (12 and 1 directions) read palindromic over two primes, and K₀ is not (proved); nothing found here forces those palindromic readings.

## The sweep: what the repository held

Swept 2026-10-20 and 2026-10-21 by three survey agents, stores by name.

- The parent page: one nucleus only; its Open list keeps the shared gradient (collective jumps), and the C# `EndCount` witness at N = 2 it names is not built. At the time of the sweep nothing in the repository cited the parent page except the experiments index and its verifier.
- `docs/ANALYTICAL_FORMULAS.md`, `docs/proofs/`, the OpenArcs registry, `docs/CAUGHT_ERRORS.md`, `fw.Confirmations`, the glossary, the typed layer in both halves (the Claim graph in Core, the witnesses in Diagnostics): no coupling between quadrupolar nuclei, no dipolar or J coupling of spins 3/2, no sodium entry. The spin operators I_x, I_y as two-qubit strings are written nowhere; the dictionary below is gated here.
- [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md) and MirrorWorld's `Atmosphere`: qubit chains under correlated Z noise, blocks that pair about their own centres when the excess price is one constant on the block and break when it is not. It holds Z jumps only and places collective Z jumps outside F158; the parent page draws the sodium consequence.
- [The Shared Field in a Crystal](THE_SHARED_FIELD_IN_A_CRYSTAL.md): multiple-quantum NMR of spin-1/2 clusters; not this object. `reflections/`: the portfolio page below; `hypotheses/` was not searched.
- The instruments: MirrorWorld's `EndCount` takes multi-site string jumps; the Sextant, the LetterCube, the Theodolite and the Python cube layer take one-site jumps only.
- The shape of S: [F118](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) builds the palindromizers as a transpose and a multiplication, Π_Y(ρ) = X^⊗N·ρᵀ among them; S is that shape with ⊗ZZ_k in place of X^⊗N, as the parent page notes for one nucleus. [The One Square](../docs/THE_ONE_SQUARE.md#9-the-moves-and-the-two-questions) §9's corner-shift row: a Pauli string commuting with H gives a copy where it is dark to every jump and a reflection where it is lit to every jump; ⊗ZZ_k is lit to four jumps per nucleus and dark to ZZ_k, and the lemma's kept part is that dark remainder. At J₀ = 0 it is lit to every jump, and for H commuting with it a colouring in the sense of [the palindrome as a colouring](THE_PALINDROME_AS_A_COLOURING.md). The string itself, ⊗ZZ_k = Z^⊗2n, is the global qubit parity: F137's S_IZ multiplies by it, §9's dark shift by Z^N is it, and it is the strong symmetry of [the Z_N partnership](Z_N_PARTNERSHIP.md); S is S_IZ's string on the left, composed with the transpose.
- The figure "a mirror maps the system onto a neighbouring one": [F137](../docs/ANALYTICAL_FORMULAS.md#f137) (an element carrying L onto the reflection of the Liouvillian with every jump A replaced by A†; a jump even under the carrier is kept, not reflected); MirrorWorld's `GammaFold` (L_anti(γ) = L(−γ) − 2σ for the whole profile, composed of the per-site turns s_l, each moving exactly the cells where site l differs by 4γ_l); [PROOF_R90_FROZEN_DIVISOR](../docs/proofs/PROOF_R90_FROZEN_DIVISOR.md) (right multiplication by X^N sends L(γ) to L(−γ) − 2σ, and the theorem is applied at −γ); [The Pairing Condition](THE_PAIRING_CONDITION.md), whose criterion asks the carrier to anticommute with every jump, so a jump it commutes with is outside it. None of them states the single-channel form below.
- Kept channels elsewhere: F82 and F84 in the registry and [PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM](../docs/proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md) give the Π residual of dissipator channels the carrier does not reflect; F80/F81 give the Hamiltonian-side residual. [The mirror proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)'s uniform longitudinal field is a different case: it pairs each block with its Π image, the real parts about the one centre, and moves the complex centre off the axis, a lab-frame shift, not a real shift of the centre.
- The construction of K_m below, the largest invariant subspace inside a coordinate subspace by Krylov stabilisation, is the one F157's blind subspace uses (the largest H-invariant subspace inside a kernel; MirrorWorld's `BlindSeat`); same construction, different operator, not joined here.
- The per-channel portfolio: [The View onto the Memory](../reflections/THE_VIEW_ONTO_THE_MEMORY.md) and Theorem 2 of [the Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (a mode's rate is twice its portfolio paired with the channel rates); its §4.2 swaps the portfolio for the whole mirror, summed over channels; `CarrierVectorPortfolio` (Diagnostics) names the sharp bit Δ ∈ {0, 1}. No store reverses the portfolio channel by channel.

Checked for adjacency: the lemma's kept part is the dark remainder of One Square §9's corner row, F137's kept jump, and GammaFold's gain turn read on the J₀ channel alone, the composite of its per-site turns there; none of them cites the parent page, nor it them. The ladder and the sharp/blurred split are the shared-field page's blockwise centres in another dissipator. These pages are linked from here; whether one statement holds them all is not decided on this page.

## Setting

Nucleus k is a qubit pair, levels |3/2⟩, |1/2⟩, |−1/2⟩, |−3/2⟩ = |00⟩, |01⟩, |10⟩, |11⟩. The spin operators, gated exactly against the spin-3/2 matrices:

  I_z = ZI + ½·IZ,  I_x = (√3/2)·IX + ½(XX + YY),  I_y = (√3/2)·IY + ½(YX − XY),  T₂₀ = (3/√6)·ZZ.

The √3 part of I_x and I_y is the satellite (3/2 ↔ 1/2 and −1/2 ↔ −3/2) and anticommutes with ZZ; the rest is the central transition and commutes with it. Per nucleus: H_k = ω_Q,k·T₂₀ + ω_k·I_z, jumps ZZ at 3J₀/2, ZX and ZY at 3J₁/2, XI and YI at 3J₂/2. Between two nuclei:

- **Ising**, d·I_z⊗I_z: Z strings only. A J or dipolar coupling reduces to it between nuclei all of whose transition frequencies differ by much more than d (another isotope, or sodium nuclei with different shifts and splittings); such nuclei have their own rates as well, which the lemma and the ladder allow.
- **Secular dipolar** between nuclei of one resonance frequency, d·(3I_zI_z − **I**·**I**): 22 strings, of which eight carry √3. The √3 strings are exactly those that anticommute with ZZ₁ZZ₂, the satellite of one nucleus against the central transition of the other. Of the rational strings, IXIX and IYIY anticommute with each ZZ_k (satellite against satellite); the rest commute with both, among them the central-transition flip-flop.
- **Its truncations** when the first-order quadrupolar splittings dominate the coupling, with equal offsets ω₁ = ω₂: the part commuting with ω_Q,1·T₂₀,₁ + ω_Q,2·T₂₀,₂. The √3 strings change one ZZ_k and drop (both splittings non-zero). IXIX and IYIY change both, by −√6(ω_Q,1·z₁ + ω_Q,2·z₂) with z_k = ±1. For **unequal** splittings, |ω_Q,1| ≠ |ω_Q,2|, they drop too, and what is left, the central-transition spin diffusion among them, commutes with every ZZ_k. For **equal** splittings, ω_Q,1 = ω_Q,2, they survive on z₁ = −z₂, as (IXIX + IYIY)·(1 − ZZ₁ZZ₂)/2, which adds ZXZX and ZYZY (for opposite splittings the resonance sits on z₁ = z₂ instead, a case not computed here). A difference of splittings comparable to d is neither truncation. Both forms are gated against the projection onto the energy shells.

All three are real (Hᵀ = H), and so are second-order quadrupolar shifts, which are diagonal and commute with every ZZ_k.

## The lemma: the mirror image is the J₀-reversed chain

Let Q = ⊗_k ZZ_k and S(ρ) = Q·ρᵀ. If H is real and commutes with Q, then

  S·L(J₀)·S⁻¹ = −L(−J₀) − 6·Σ_k(J₁,k + J₂,k),

with L(−J₀) the same chain with every J₀,k reversed in sign, a negative spectral density: not a physical, completely positive generator, the place the mirror image lives. *Proof.* The transpose reverses −i[H, ·] because H is real, and multiplying by Q keeps it because Q commutes with H. A jump P is a Pauli string, so Pᵀ = ±P, and the sign enters twice: (PρP)ᵀ = PρᵀP. With Q·P = ±P·Q, S carries the dissipator γ(PρP − ρ) to itself when P commutes with Q and to −γ(PρP − ρ) − 2γ when it anticommutes. Only ZZ_k commutes with Q; the others, at total rate 3(J₁,k + J₂,k) on nucleus k, are reflected. Reversing the sign of the kept J₀ terms is the same as reversing the J₀,k. The graph of the coupling is never used. ∎

Written as the parent page writes it, S·L·S⁻¹ = −L + 2·L_{J₀} − 6·Σ_k(J₁,k + J₂,k) with L_{J₀} = −3·Σ_k J₀,k·P_k, P_k the projector on the cells bright for nucleus k. The verifier gates the identity with residual exactly 0 at two and three nuclei for the Ising coupling and both truncations, for the Ising coupling also with rates per nucleus, and its control: the √3 part alone breaks it.

## The ladder

Every Pauli jump keeps each cell's brightness. If H commutes with every ZZ_k as well, each P_k commutes with L, the cells split by the set B of bright nuclei, and on the cells of B, L_{J₀} = −3·Σ_{k∈B} J₀,k. There the lemma reads S·L·S⁻¹ = −L − 2s_B:

  **every part with bright set B is palindromic about −s_B, s_B = 3·Σ_k(J₁,k + J₂,k) + 3·Σ_{k∈B} J₀,k.**

For identical rates s_B depends only on m = |B|, s_m = 3n(J₁ + J₂) + 3m·J₀, m = 0, …, n, and the rungs are 3J₀ apart, the J₀ cost of one bright nucleus; with rates per nucleus there are up to 2ⁿ centres. For n = 1 they are the parent page's dark and bright centres, c₁ and c₂. The Ising coupling commutes with every ZZ_k, and so does the secular dipolar coupling between nuclei of unequal, non-zero splittings, because its surviving flip-flop is the central one, which moves no brightness. The verifier gates [L, P_k] = 0 exactly at two and three nuclei, and at two nuclei the characteristic polynomial of every pattern block against its own reflection over two primes, for identical rates and, with the Ising coupling, for rates per nucleus, with the exact controls that a block is not palindromic about a neighbouring centre. The rungs with m ≥ 2 belong to multi-spin modes, which a plain spectrum does not show; reading them would need multiple-quantum or filtered experiments.

## Equal splittings: the ladder where brightness stays sharp

With equal splittings the satellite flip-flop between neighbours is resonant: one nucleus goes outer to inner while the other goes inner to outer. The coupling still commutes with Q, so the lemma holds, but it moves brightness: it changes the number m of bright nuclei by 0 or ±2, P_k no longer commutes with L, and only the parity of m is kept. The rungs are then exact on a part of the spectrum.

Take identical rates here. Let K_m be the largest subspace that L leaves invariant inside the cells with m bright nuclei, W = Σ_k P_k. *If the lemma holds, S maps K_m onto itself and L restricted to K_m is palindromic about −s_m.* Proof: S keeps every P_k, so it keeps the cells with m bright nuclei. From the lemma, L = −S·L·S⁻¹ − 6n(J₁ + J₂) − 6J₀·W, so for x in K_m, L·Sx = S(−Lx) − 2s_m·Sx lies in S·K_m; S·K_m is invariant and inside the m cells, hence contained in K_m, and equal to it by dimension. On it S·L·S⁻¹ = −L − 2s_m. ∎ Off the K_m the rungs are not forced; whether any of those modes pair is not tested here.

For identical nuclei (one splitting, one offset, identical rates) on an open nearest-neighbour chain with one coupling strength, at the verifier's parameters:

| nuclei | K₀ | K₁ | K₂ | K₃ | at most, in all |
|---|---|---|---|---|---|
| 2 | 24 of 64 | 128 of 128 | 24 of 64 | | 176 of 256 |
| 3 | 64 of 512 | 140 of 1536 | 84 of 1536 | 16 of 512 | 304 of 4096 |

The dimensions are computed over two primes and are upper bounds on the complex ones at these parameters, and so at generic ones (a dimension can only rise at special points); at two nuclei they are the same at the rate points (6, 4, 2) and (8, 6, 2). K₁ at two nuclei is exactly the whole block with one bright nucleus, which the paired flip-flop cannot leave. The verifier also reads each K_m palindromic about −s_m over both primes; that follows from the lemma over GF(p) and checks the code, not the proposition. The sharp share is at most 69 % at two nuclei and at most 7.4 % at three. At J₀ = 0 all rungs coincide, and the whole spectrum is palindromic at every length, by the lemma.

## The full secular dipolar coupling

Untruncated, the coupling keeps its eight √3 strings, which anticommute with Q, and S is no symmetry of the form above. This is the Hamiltonian where the splitting does not dominate the coupling, at a cubic site in a solid; it is computed at zero splitting. At J₀ = 0, with J₁, J₂ > 0, the global mirror, about −σ, the only candidate centre, breaks exactly. Every jump is Hermitian and squares to 1, so [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) applies: the far end of its count lies in the span of the strings that anticommute with every jump, which on each nucleus must read ZZ, so at every length the span is Q alone; Q does not commute with H, the far end is empty, the near end holds the identity, and that holds for any non-zero coefficient of the √3 strings, √3 included. The verifier enumerates that span at two and three nuclei, with the control that it is empty once the J₀ jumps are present, and gives the break a second route at two nuclei: the characteristic polynomial is not palindromic about −σ modulo a prime, which is exact.

What remains: the sharp subspaces have dimensions at most 10, 12 and 1 at two nuclei, and 12, 16, 0, 0 at three. K₀ is not palindromic about −s₀, at every length, and that is proved: the identity lies in K₀ with eigenvalue 0, so −2s₀ would have to be an eigenvalue of L on the cells with no bright nucleus. The real part of an eigenvalue is ⟨x, Dx⟩/⟨x, x⟩ for its eigenvector x, the Hamiltonian part being anti-Hermitian, and on those cells the dissipator D is diagonal in the strings with −2s₀ as its most negative entry, reached on Q alone (gated at two and three nuclei). So the eigenvector would be Q, and Q is no eigenvector, since [H, Q] ≠ 0. Over both primes K₁ and K₂ read palindromic about their rungs at two nuclei, and K₁ at three. Without S nothing found here forces those readings, and a reading modulo p is not a proof of them.

## A reading

The Absorption Theorem prices a mode by its portfolio, the share it holds in each channel. S reverses the share in every channel it reflects and keeps the share in the channel it cannot see, J₀, whose jump ZZ_k is S's own carrier. A mode and its image therefore have rates summing to the reflected rates plus the kept channel's share counted twice. Where that share is sharp, a whole number of bright nuclei, the kept channel adds a constant and the mirror lives on with a shifted centre: the ladder. Where the coupling blurs it, the share differs from mode to mode, and the mirror keeps its element but loses its centre; its image is then a mode of the J₀-reversed chain, not of this one. Where the element itself breaks, at the √3 strings, that is gone as well. This is a reading of the lemma, not a further result.

## What this does not show

No sodium data were compared. The model is the parent page's secular Redfield picture at infinite temperature, with an independent gradient on every nucleus, plus a static coupling; that the two coexist in a real sample is assumed, and no such sample is named. The independent jumps per nucleus are the local form of the master equation, valid for d·τ_c ≪ 1 with τ_c the correlation time of the gradient, and the relaxation is secular with respect to each nucleus's own levels, not in the coupled eigenbasis, which would give jumps spread over the nuclei. Outside the model and outside F158's class, each of these would in general break the ladder or the lemma, and none is computed: dipolar relaxation from the motion that modulates the coupling, and its cross-correlation with the quadrupolar relaxation, whose jumps (I_z⊗I_z-like, flip-flop-like) are sums of strings; relaxation through I_z or I_± from other spins, chemical-shift anisotropy or paramagnetic centres, again sums of strings; field gradients that are partly correlated between nearby nuclei, the generic case between independent and shared. The lemma and the ladder are proved for any number of nuclei and any graph; the dimensions of K_m and the readings under the full coupling are computed at two and three nuclei, for identical nuclei on an open chain, at the rates and couplings of the verifier. The literature on dipolar couplings between quadrupolar nuclei was not searched. For isotropic motion J₀ is the largest density and the statements at J₀ = 0 are empty; they hold for motion that only turns the gradient about the field, where J₀ = 0 (the parent page).

## Open

- What carries the palindromic readings under the full secular dipolar coupling (K₁ and K₂ at two nuclei, K₁ at three), and whether they hold over the complex numbers.
- The sharp subspaces K_m for equal splittings beyond three nuclei: a count, or a formula; opposite splittings.
- A shared or partly correlated gradient, and the relaxation channels named above.
- A sample in which a static coupling between sodium nuclei and their motionally narrowed relaxation coexist; for the full coupling, a cubic site in a solid with a relaxation of this form.
- Whether one statement holds this ladder, F137's kept jump, GammaFold's gain turn, the corner-shift row of One Square §9 and the shared-field page's blockwise centres.
- A typed home with a live C# witness, if the ladder becomes a claim.

## Links

- [The Sodium Nucleus: Two Mirrors](THE_SODIUM_NUCLEUS_TWO_MIRRORS.md), the single nucleus, its theorem and its conventions.
- [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md), the two-end count, [the palindrome as a colouring](THE_PALINDROME_AS_A_COLOURING.md) and [F118](../docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md), the palindromizer as a transpose and a multiplication.
- [The One Square](../docs/THE_ONE_SQUARE.md#9-the-moves-and-the-two-questions) §9, the corner shift.
- [F137](../docs/ANALYTICAL_FORMULAS.md#f137) and MirrorWorld's `GammaFold` (`compute/MirrorWorld/GammaFold.cs`), mirrors onto a neighbouring Liouvillian.
- [The Mirror under a Shared Field](THE_MIRROR_UNDER_A_SHARED_FIELD.md), blockwise centres for qubit chains.
- [The View onto the Memory](../reflections/THE_VIEW_ONTO_THE_MEMORY.md) and [the Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), the portfolio of a mode.
