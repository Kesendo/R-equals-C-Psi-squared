# Chain Dissipation-Gap Sector Diagnostic: the slow mode is a near-stationary magnon-admixture

**Status:** Tier 1 candidate (3 structural findings, verified at N=4, 5, 6; empirical prefactor 0.55·Q²/N² at Q = 2, matching to ~1%, its Q → 0 value N²(1 − cos(π/N))/8 closed). Closes the "which weight sector hosts the slow mode" question from F1_DISSIPATION_GAP_PATTERN.md Q5.
**Date:** 2026-05-19
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:** [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md), [the weight-1 degeneracy proof](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md) (F50's 2γ floor on the distance-1 coherences, and its Zeno end for the Q → 0 values), [the F1 dissipation-gap pattern](../hypotheses/F1_DISSIPATION_GAP_PATTERN.md)

**Verification:** [`simulations/chain_gap_sector_diagnostic.py`](../simulations/chain_gap_sector_diagnostic.py) (the Heisenberg chain, H = (J/4)·Σ_b(X_iX_j + Y_iY_j + Z_iZ_j) with uniform Z-dephasing γ, N=4, 5, 6 at γ=0.5, J=1, Q=2)

---

## What this document is

The 2026-05-19 Q-sweep + Absorption Theorem reading proposed that the chain dissipation gap lives in a mixed sector with fractional `⟨n_XY⟩ ≪ 1`, not in the pure-w=1 sector that F50 pins at 2γ. This experiment opens the black box: it block-diagonalises L by joint-popcount sectors `(p_col, p_row)`, runs eigendecomposition per block, identifies the slow mode by smallest `|Re(λ)|`, and reads off both its sector and its Pauli-basis light content. Three structural findings emerge.

## Three structural findings

### Finding 1. Slow mode lives in the **central diagonal popcount block**

At chain N=4, 5, 6 with γ=0.5, J=1 (Q=2), the slow mode sits in the `(⌈N/2⌉, ⌈N/2⌉)` joint-popcount sector:

| N | slow-mode sector | block size | gap |
|---|---|---:|---:|
| 4 | (2, 2) | 36  | 0.13616 |
| 5 | (3, 3) | 100 | 0.08837 |
| 6 | (3, 3) | 400 | 0.06069 |

This is the **largest** joint-popcount block, of dimension `C(N, ⌈N/2⌉)²` (the binomial-squared central peak). The slow mode is NOT in any off-diagonal-popcount sector `(k, k±1)`, NOT in any boundary `(0, 0)` or `(N, N)` sector.

Every off-diagonal-popcount sector `(k, k±1)` has its **slowest** mode at exactly `2γ`, and that is a floor rather than a description of the sector. The floor follows: an operator in `(k, k±1)` is built from coherences `|α⟩⟨β|` whose popcounts differ, so their Hamming distance is at least 1, so `⟨n_XY⟩ ≥ 1`, so `Re ≤ −2γ` by the Absorption Theorem, with equality exactly on the distance-1 coherences F50 pins. The sectors spread well above their floor: at N=5, γ=0.05 the `(1,2)` sector's 50 eigenvalues run `Re ∈ [−0.300, −0.100]`. Only the two END sectors, `(0,1)` and `(N−1,N)`, are distance-1 throughout and therefore sit at `−2γ` as a whole. Checked at N=3, 4, 5 over four (J, γ) pairs: every off-diagonal sector is exactly closed, every one attains `2γ` as its smallest `|Re|`, and the flat ones are the two end sectors and no others. Either way the diagonal sectors carry the actual gap structure, because the gap is below `2γ` and no off-diagonal sector can go there.

The X⊗N pairing `(k, k) ↔ (N−k, N−k)` (Π², F1²) holds at the sector level, exactly, since H commutes with X⊗N, and to the four decimals quoted in the readings: e.g. for N=6, sectors `(2, 2)` and `(4, 4)` both report slow eigenvalue −0.0626; `(1, 1)` and `(5, 5)` both report −0.0681. The two per-sector eigensolver runs were never compared below that precision, so this is agreement at the quoted width and not a measured residual. The slow mode of the full L is the smallest of these X⊗N-paired per-sector minima, which is always the central `(⌈N/2⌉, ⌈N/2⌉)` block.

### Finding 2. Absorption Theorem holds on the slow mode, to the floor of the two float routes

The Absorption Theorem prediction `Re(λ_slow) = −2γ·⟨n_XY⟩_slow` is exact as a theorem. What this table compares is two independent floating-point routes to it: the eigensolver's `Re λ` on the block, and the Pauli-basis projection's `⟨n_XY⟩`. They agree to

| N | gap | 2γ·⟨n_XY⟩_slow | relative deviation |
|---|---:|---:|---:|
| 4 | 0.13616 | 0.13616 | 1.8e-15 |
| 5 | 0.08837 | 0.08837 | 6.8e-14 |
| 6 | 0.06069 | 0.06069 | 1.0e-13 |

which is a verification that the Absorption Theorem is the right reading of F3 for the gap question: the slow-mode decay rate is `2γ` times its Pauli-basis light content.

Read the third column as an error floor between two float routes, not as a quality of the theorem. It is not bit-exactness, and it was reported as such here until 2026-08-06: the producing script formatted the deviation as `{:.3f}%`, so everything under 5e-6 printed `0.000%`, and this document read that printf width as exactness. The script now prints the deviation in scientific notation for exactly that reason.

Do not read a trend into the third column, and in particular not the "grows with N" one that a first draft of this repair put here. It is a RELATIVE deviation and the gap it divides by falls like 1/N², so most of the apparent growth is the denominator. In absolute terms the deviation is 2.5e-16 / 6.0e-15 / 6.3e-15 at N = 4 / 5 / 6: one step up between N=4 and N=5, then flat to 4% between N=5 and N=6. Three points at a single (γ, J) do not carry a scaling law, and no error model has been stated here, so the honest reading is a floor of order 1e-15 in absolute size with its N-dependence unmeasured.

The empirical form `⟨n_XY⟩_slow ≈ 0.55·Q²/N²` at Q = 2 agrees with the measurements to ~1%:

| N | ⟨n_XY⟩_slow observed | 0.55·Q²/N² predicted | error |
|---|---:|---:|---:|
| 4 | 0.13616 | 0.13750 | 1.0% |
| 5 | 0.08837 | 0.08800 | 0.4% |
| 6 | 0.06069 | 0.06111 | 0.7% |

The 0.55 coefficient is approximately N-independent; the small finite-N drift suggests a sub-leading `1/N²` correction that vanishes as N → ∞. Its Q → 0 value is closed: at strong dephasing the populations relax under the exclusion graph's Laplacian, whose gap in every block of the chain is the path's ([PROOF_WEIGHT1_DEGENERACY](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-zeno-end), the Zeno end), so D06's asymptote gives, in this document's book H = (J/4)·Σ σ_i·σ_j, `⟨n_XY⟩_slow·N²/Q² → N²(1 − cos(π/N))/8` = 0.586, 0.597, 0.603 at N = 4, 5, 6 and π²/16 = 0.617 as N → ∞. The 0.55 at Q = 2 is that value moved by the finite coupling, and the form of that drift is open.

### Finding 3. The slow mode is a **near-stationary magnon-admixture**

The Pauli-basis weight distribution of the slow mode is sharply peaked at `n_XY = 0`:

| N | weight n_XY=0 | weight n_XY=2 | weight n_XY=4 |
|---|---:|---:|---:|
| 4 | **93.20%** | 6.80% | 0.00% |
| 5 | **95.59%** | 4.41% | 0.00% |
| 6 | **96.97%** | 3.03% | 0.00% |

The slow mode is **93-97% pure I/Z Pauli strings** (operators diagonal in the computational basis, the "dark" / n_XY=0 sector that F50 pins at the kernel = Re=0). The remaining 3-7% is **n_XY=2** (two X/Y operators, typically an XX or YY pair: a single magnon excitation). The n_XY=4 weight is below 10⁻⁴ at all N.

This explains everything: `⟨n_XY⟩_slow = 0·w_0 + 2·w_2 + 4·w_4 + ... ≈ 2·w_2` since w_4 vanishes. The empirical amplitude is `w_2 ≈ 0.275·Q²/N²` (so that `⟨n_XY⟩ ≈ 2·w_2 = 0.55·Q²/N²` matches the table above). The `Q²/N²` scaling is the perturbation-theoretic order of magnitude (mixing amplitude `~ J·k_min/2γ ~ Q/N`, squared); at Q → 0 the 0.275 becomes N²(1 − cos(π/N))/16, half of the closed value above, and its drift with Q is what a further calculation would need to derive in closed form.

The physical picture: in the limit Q → 0 (no Hamiltonian, only dephasing), the slow mode is **exactly stationary** (n_XY=0, in the kernel of L). Turning on H mixes a small magnon excitation (n_XY=2, since one nearest-neighbour XY-bond flip introduces two X or Y letters), carried by `k_min = π/N`, the lowest Neumann cosine of the population diffusion, into the otherwise-stationary mode. The dephasing dissipator only "sees" this small magnon admixture; the decay rate of the mode is therefore `2γ⟨n_XY⟩ = 4γ × w_2 ≈ γ·Q²/N²`. The slow mode is a **nearly-conserved operator dressed with a small magnon component**, and the magnon-mixing amplitude is the gap prefactor.

## Structural role of the admixture

The 3-7% magnon admixture is small in amplitude but plays two structurally large roles that the bare gap-prefactor analysis does not surface.

### Role 1: the loophole channel for an otherwise-conserved quantity

In the absence of H (Q = 0 limit), the I/Z-only Pauli content of the slow mode would be exactly stationary: every diagonal-popcount operator commutes with the dephasing dissipator (`L_D` is diagonal in the computational basis) and Z-magnetisation is the conserved charge of the system. The kernel of `L_D` restricted to a diagonal popcount block is the full set of populations on that block; nothing decays.

Turning on H opens precisely one decay channel for that conserved subspace: the XX+YY pair-flip term mixes a small n_XY=2 magnon coherence into the otherwise-stationary mode. The dephasing dissipator then acts on this small admixture, not on the (still-protected) population content. The slow-mode decay rate is therefore

    gap  =  2γ · w_admixture · (n_XY of the admixture)  ≈  2γ · w_2 · 2  =  4γ · w_2.

The empirical amplitude is `w_2 ≈ 0.275·Q²/N²` (chain plateau N ≥ 4 from the Q-sweep; the plateau itself has ~10% sub-Q² drift across Q ∈ [0.5, 2.5], so `c = 1.10` should be read as the plateau-mean of `c(Q)`, not a perfect constant), giving `gap ≈ 4γ · 0.275 · Q²/N² = 1.10·γ·Q²/N²` consistent with the observed chain plateau. The `Q²/N²` scaling is the perturbation-theoretic order of magnitude (mixing amplitude `~ J·k_min/2γ ~ Q/N`, squared); the prefactor 0.275 is the finite-Q value, its Q → 0 limit N²(1 − cos(π/N))/16 is closed as above, and its drift with Q awaits a closed form.

The admixture is therefore **the unique decay channel** for the otherwise-conserved population content: without it, gap = 0 exactly, and the slow mode would be a true zero-mode of L_H + L_D. The magnon admixture is the structural "loophole" that lets dissipation reach an operator that is otherwise protected by conservation. The small size of the admixture (3-7%) is what makes the slow mode slow: a 50% admixture would land at the 2γ floor F50 pins (no slow mode at all), and zero admixture would make the mode a kernel addition (infinite lifetime). The empirical `Q²/N²` scaling is the scaling of the loophole's opening.

### Role 2: the synthesis point for F1², F2, F3, F50, and the Absorption Theorem

Each of the May 2026 typed claims plays an explicit role in the admixture's structure, and together they form a self-consistent decomposition of the slow mode:

- **F50** (`PROOF_WEIGHT1_DEGENERACY`): if the admixture lived alone in an off-diagonal popcount sector, it would be pinned at Re = −2γ exactly. The relation is an inclusion and not an equivalence: a weight-1 Pauli string carries one X or Y letter, so on a computational-basis state |α⟩⟨β| it flips exactly one bit, and the weight-1 operators span the coherences at Hamming distance 1. Those sit inside the off-diagonal popcount sectors `(k, k±1)` but do not fill them, because a popcount step of ±1 also admits Hamming distance 3, 5, and so on: at N=5 the weight-1 span is 160-dimensional against the sectors' 420. Only the two end sectors, `(0,1)` and `(N−1,N)`, lie at distance 1 throughout, and only they sit at Re = −2γ as a whole; the interior ones do not. Measured at N=5, γ=0.05, the `(1,2)` sector's 50 eigenvalues run Re ∈ [−0.300, −0.100]. So F50 pins the distance-1 coherences, not whole popcount sectors, and D10 Step 6 carries the scope. The admixture inherits this 2γ scale; the slow-mode rate is `2γ·⟨n_XY⟩ = 4γ × (admixture weight)`.
- **F2** (`F2W1DispersionPi2Inheritance`): the magnon component carries the open-chain dispersion `ω_k = 4J·(1 − cos(πk/N))`. The slowest magnon mode is at k = 1 with ω_1 ≈ 2π²·J/N², which is the "kinetic" frequency that drives the mixing. F2 explains why the admixture-amplitude scales as `Q/N`: the mixing is set by the ratio of H's hopping rate (k_min × J) to the dissipator's decay rate (2γ). F2 does not describe the slow mode itself, in two separate ways. Its object is the `(0,1)` coherence block, not any off-diagonal sector as a whole, and the slow mode is not in an off-diagonal sector at all; and the slow mode has `Im(λ) ≈ 0` to machine precision, because it lives in a diagonal-popcount sector where nothing oscillates. What F2 supplies here is the magnon's intrinsic frequency scale, borrowed.
- **F3 / Absorption Theorem**: the operator-level identity `Re(λ) = −2γ·⟨n_XY⟩` reads the decay rate of any Lindblad eigenmode directly from its Pauli-basis light content. Applied to the slow mode (`⟨n_XY⟩ ≈ 2·w_2`), it gives the gap.
- **X⊗N pairing**: the slow mode at sector `(k, k)` is partnered by X⊗N (Π², F1²) with a mode at sector `(N−k, N−k)` with identical decay rate. The per-block analysis confirms this to the precision the readings are quoted at: e.g. for N=6, sectors `(2, 2)` and `(4, 4)` both give slow eigenvalue −0.0626; `(1, 1)` and `(5, 5)` both give −0.0681. The admixture obeys the X⊗N symmetry by inheriting it from its host population.

Read together: **X⊗N, the square of F1, pairs the diagonal sectors around the center, F50 puts a 2γ floor under the off-diagonal popcount sectors, F2 governs the magnon's intrinsic frequency, F3 / Absorption Theorem reads decay from light content, and the admixture is where all four meet.** The slow mode at the central diagonal popcount block (⌈N/2⌉, ⌈N/2⌉) is where, at finite Q, each formula contributes one structural ingredient and all four compose; that it sits at the centre is none of the four's doing but the ZZ term's, the XY chain, which holds the same four with F2b's spectrum in F2's place, tying every filling (at the Zeno end, Theorem E (d) of [the weight-1 proof](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-zeno-end)).

The 0.55 is non-trivial in its drift, not in its limit: at Q → 0 two of the four inputs suffice and the X⊗N pairing and the F50 floor drop out, the Zeno generator being the exclusion graph's Laplacian, blind to the anisotropy, whose gap read through the Absorption Theorem gives the N²(1 − cos(π/N))/8 above, that gap being F2's slowest frequency times J/γ in the Pauli book, one Laplacian running the (0,1) block and the populations (F152). The drift from there to 0.55 at Q = 2 is what remains, the ZZ term turning it downward (item 1 below).

---

## Extensions (resolved 2026-05-19) and remaining open work

Items 2-5 from the original open list were closed in a sector-diagnostic sweep on 2026-05-19 (`simulations/slow_mode_sector_sweep.py`). Item 1, the closed form of `c ≈ 0.55`, is closed at Q → 0 (above), the ZZ term's share of its drift at the first step by Theorem E (d), and open in the XY chain's own drift with Q.

### Item 2 resolved: Ring N=4..6 sector

Ring slow mode lives in the **central diagonal popcount sector**, same as chain, with a magnon admixture the Absorption Theorem reads the same way. The admixture amplitude is 3-5× larger than chain at each N (Q=2 anchor):

| N | sector | gap | ⟨n_XY⟩ | w_0 | w_2 | w_4 |
|---|---|---|---|---|---|---|
| 4 | (2, 2) | 0.379 | 0.379 | 0.812 | 0.187 | 0.001 |
| 5 | (3, 3) | 0.317 | 0.317 | 0.842 | 0.157 | 0.001 |
| 6 | (3, 3) | 0.230 | 0.230 | 0.885 | 0.114 | 0.000 |

Predicted `⟨n_XY⟩ = 2·Q²/N²` (4× chain coefficient) gives 0.500 / 0.320 / 0.222 for N=4/5/6; observed 0.379 / 0.317 / 0.230 (N=5 within 1%; N=4 has finite-size deviation, N=6 within 4%). The "4× ring/chain prefactor matches cyclic-vs-open k_min² ratio" reading from the F1_DISSIPATION_GAP_PATTERN doc is structurally confirmed.

### Item 3 resolved: Star N=3..6 sector (surprise: NOT central)

Star slow mode lives at **boundary popcount sectors** `(1, 1)` or `(N−1, N−1)` (X⊗N partners), NOT central:

| N | sector | gap | ⟨n_XY⟩ | w_0 | w_2 |
|---|---|---|---|---|---|
| 3 | (1, 1) | 0.270 | 0.270 | 0.865 | 0.135 |
| 4 | (3, 3) | 0.210 | 0.210 | 0.895 | 0.105 |
| 5 | (1, 1) | 0.164 | 0.164 | 0.918 | 0.082 |
| 6 | (5, 5) | 0.130 | 0.130 | 0.935 | 0.065 |

This is the structural signature of the star's separate scaling family at Q = 2 (`gap ~ 1/N` rather than `1/N²` there; toward Q → 0 the star is a tree, its Zeno gap N-independent and shared by every block, and the ZZ term chooses at the next order, detuning a hop between the hub and an arm by the imbalance of the other arms, most at the edges of filling, by the detuning form of [the Zeno end](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-zeno-end), the hops' own correction read shared on the stars N = 4 to 7 and not proved there): the slow mode sits at the popcount-boundary sector `(1, 1)` (or its X⊗N partner `(N−1, N−1)`), i.e. the sector of single-excitation operators on either bra or ket. The admixture-as-channel picture still holds (gap = 2γ·⟨n_XY⟩), but the channel content sits at the popcount boundary rather than the centre.

Promotion implication: the chain reading "slow mode in central diagonal popcount sector" was N-universal for chain and ring, but NOT for star. Future "slow mode lives at the central popcount block" statements need a topology qualifier (the chain and the ring central, the star at the boundary, each where the ZZ term detunes the hops most at the Zeno end; without ZZ the chain ties every filling, and the star keeps the boundary at the order after, (J/γ)⁶).

### Item 4 resolved: N=7, 8, 9 chain sector confirmed

| N | source | sector | gap | ⟨n_XY⟩ | 0.55·Q²/N² (Q = 2 form) | match |
|---|---|---|---|---|---|---|
| 4 | dense N=4 | (2, 2) | 0.1362 | 0.1362 | 0.1375 | 1.0% |
| 5 | dense N=5 | (3, 3) | 0.0884 | 0.0884 | 0.0880 | 0.4% |
| 6 | dense N=6 | (3, 3) | 0.0607 | 0.0607 | 0.0611 | 0.7% |
| 7 | dense N=7 | **(4, 4)** | 0.0450 | 0.0450 | 0.0449 | 0.2% |
| 8 | SLOW_N8 sweep + AT | (4, 4) | 0.0344 | 0.0344 | 0.0344 | 0.06% |
| 9 | MklDirect bridge + AT | (4, 4) ≡ (5, 5) X⊗N-paired | 0.0273 | 0.0273 | 0.0272 | 0.4% |

Central-popcount-block reading holds at N=4..9 chain. The reading is the identification of a sector INDEX, an integer, so there is no tolerance in it either way. The N=8 and N=9 numbers are read via Absorption Theorem `⟨n_XY⟩ = gap/(2γ)` from existing JSON metric files (no new compute required) plus MaxBlockSectorPCol/PRow which both report the central popcount block.

### Item 5 resolved: c=0.55 drifts ~10% with Q at fixed N

Q-sweep at chain N=5, γ₀=0.05, across the six canonical Q-anchors gives `⟨n_XY⟩(Q) / (0.55·Q²/N²)` from 1.08 at Q=0.5 down to 0.97 at Q=2.5:

| Q | ⟨n_XY⟩ | 0.55·Q²/N² | ratio |
|---|---|---|---|
| 0.5  | 0.00593 | 0.00550 | **1.079** |
| 1.0  | 0.02334 | 0.02200 | **1.061** |
| 1.5  | 0.05123 | 0.04950 | **1.035** |
| √3   | 0.06739 | 0.06600 | **1.021** |
| 2.0  | 0.08837 | 0.08800 | **1.004** |
| 2.5  | 0.13352 | 0.13750 | **0.971** |

The "0.55" coefficient is therefore Q-specific: ~0.59 at Q=0.5, 0.55 at Q=2, ~0.53 at Q=2.5. The drift matches the ~10% sub-Q² drift in the chain plateau f(Q)/Q² documented separately in `F1_DISSIPATION_GAP_PATTERN.md`. Closed-form derivation needs to produce a c(Q) function, not just a single number.

The Q-sweep also surfaces that the slow-mode sector at N=5 alternates between `(2, 2)` and `(3, 3)` across Q values. Both are X⊗N-paired (N=5 central is `⌈5/2⌉ = 3` so `(2, 2)` and `(3, 3)` are X⊗N partners carrying the same spectrum up to eigensolver noise, the two blocks being diagonalised separately), so the "winner" is numerical chance from the eigensolver. The sector identity is "(2,2)+(3,3) X⊗N pair", not a single block.

### Framework-convention cross-check at Q=1.5 (γ₀=0.05, J=0.075)

Re-running the sector diagnostic at the F86 Q_peak c=2 canonical anchor (`Q=1.5` from `docs/Q_REGIME_ANCHORS.md`) gives a consistent reading across all three topologies:

| topology | N | sector | ⟨n_XY⟩ | predicted | ratio |
|---|---|---|---|---|---|
| chain | 3 | (1, 1) X⊗N=(2,2) | 0.1462 | 0.55·Q²/N² = 0.1375 | 1.063 |
| chain | 4 | (2, 2) | 0.0788 | 0.0773 | 1.019 |
| chain | 5 | (3, 3) | 0.0512 | 0.0495 | 1.035 |
| chain | 6 | (3, 3) | 0.0355 | 0.0344 | 1.032 |
| ring | 3 | (1, 1) | 0.4665 | 2·Q²/N² = 0.5000 | 0.933 |
| ring | 4 | (2, 2) | 0.2413 | 0.2813 | 0.858 (the Heisenberg ring-N=4 Im-max bound interferes) |
| ring | 5 | (2, 2) X⊗N=(3,3) | 0.1858 | 0.1800 | 1.033 |
| ring | 6 | (3, 3) | 0.1337 | 0.1250 | 1.070 |
| star | 3 | (2, 2) X⊗N=(1,1) | 0.1462 | (boundary, no Q²/N² form) | – |
| star | 4 | (1, 1) | 0.1269 | – | – |
| star | 5 | (4, 4) X⊗N=(1,1) | 0.1074 | – | – |
| star | 6 | (5, 5) X⊗N=(1,1) | 0.0902 | – | – |

Chain ratio stays around 1.03 at Q=1.5 across N (matches the Q-sweep prediction). Ring N≥5 sits between 1.03 and 1.07 (consistent with cyclic-vs-open k_min² factor 4); ring N=3,4 sit below 1.0 due to dihedral-lock finite-size interference (a separate Im-max bound for N=4 that the gap doesn't follow cleanly). Star at every N=4..6 reports a non-central popcount sector (one of the X⊗N-paired boundary sectors), confirming the topology-distinct scaling family at the framework's canonical Q anchor.

The framework's `lebensader.py + cockpit_panel` workflow defaults run at this convention; the chain N=5, ring N=5, star N=5 numbers above are therefore directly comparable to any hardware data taken under the same convention.

## Item 1: closed at Q → 0, open in its drift

At N = 5 the coefficient runs 0.593, 0.584, 0.569, 0.552, 0.534 at Q = 0.5, 1, 1.5, 2, 2.5, down from the Zeno value N²(1 − cos(π/N))/8 = 0.597 above and even in Q, as the Zeno expansion requires (c(Q) = c(0) + O(Q²), no linear term); the form of that drift is open. MEP 2016 (`arXiv:1606.09122`) gives `2π² · Q² · γ/N²` for periodic XX; in this document's book that is `2π²/4 = π²/2 = 4.93`, the ring's Zeno-end coefficient `gap·N²/(γQ²) → N²(1 − cos(2π/N))/4 → π²/2`, which the Zeno generator gives for XX and XXX alike (the ring's gap equality holds by the outside route the Zeno end recognizes, and gate row Z7 of [`f50_zeno_end_ferromagnet.py`](../simulations/f50_zeno_end_ferromagnet.py) reads it in every block of the rings N = 4 to 9). The factor of about 9 between it and the chain's `c(Q=2) ≈ 0.55` is 4 (ring against chain) times 2 (a gap against ⟨n_XY⟩), neither of which sees the ZZ term, times 1.12, which is finite N and finite Q together (at N = 5, 0.617/0.597 = 1.034 times 0.597/0.552 = 1.081). The ZZ term turns that drift downward: at Q = 2 the XX chain sits at 0.632, 0.624, 0.620 for N = 4, 5, 6, above its Zeno values 0.586, 0.597, 0.603, the Heisenberg chain at 0.545, 0.552, 0.546, below them (gate row Z8 prints these). That downward turn is Theorem E (d)'s at its first step: the hops' own correction is the same with and without the ZZ term, so c_Heis(Q)/c_XY(Q) = 1 − δ_p/|λ₀| + O(Q⁴) = 1 − b·(Q/4)² + O(Q⁴) in this book, b = [(N − 4 sin²(π/N))·φ(p) + sin²(π/N)]/N with φ(p) = 2(p − 1)(N − p − 1)/((N − 2)(N − 3)), 0.5515 in the central blocks at N = 5 ([the Zeno end](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-zeno-end), δ_p itself read on the Liouvillian in gate row Z12). What stays open is the XY chain's own drift, upward from its Zeno value. The 4× chain-to-ring ratio in our data matches the open vs cyclic k_min² ratio.

## Open extensions (not closed today)

- **Closed-form derivation of the XY chain's drift of c(Q)** away from its Zeno value (item 1 above; the ZZ term's share is Theorem E (d)'s at the first step): the dominant remaining analytical item.
- **Star sector beyond N=6**: at N=7 the dense N=7 method costs ~50s; star N=7,8 sector via the block-spectrum bridge could confirm the boundary-popcount-sector reading at scale.
- **Ring N≥7 sector**: similar; would verify that the 4× chain-to-ring prefactor and central-popcount-sector picture persist.
- **The drift of c(Q) at other N**: the Q-sweep of item 5 is at N=5 only; an N=4, 6 sweep would show how the Q² drift away from the Zeno value depends on N.

## Reproduction

```
python simulations/chain_gap_sector_diagnostic.py
```

Runs N=4, 5, 6 chain at γ=0.5, J=1; outputs the slow-mode sector + light content + per-block slow eigenvalues + Pauli-basis weight distribution to stdout. Total wall time: ~5 seconds on a standard desktop.

## Cross-references

- Parent doc (the open question that surfaced this finding): [`hypotheses/F1_DISSIPATION_GAP_PATTERN.md`](../hypotheses/F1_DISSIPATION_GAP_PATTERN.md) Q5.
- Companion Tier-1-derived Im_max bounds from the same sprint: [`RingN4DihedralLockClaim`](../compute/RCPsiSquared.Core/Symmetry/RingN4DihedralLockClaim.cs), [`StarImMaxBoundClaim`](../compute/RCPsiSquared.Core/Symmetry/StarImMaxBoundClaim.cs).
- Absorption Theorem (the engine of Finding 2): [`docs/proofs/PROOF_ABSORPTION_THEOREM.md`](../docs/proofs/PROOF_ABSORPTION_THEOREM.md).
- F50's 2γ floor (the engine of the per-sector hierarchy in Finding 1): [`docs/proofs/PROOF_WEIGHT1_DEGENERACY.md`](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md).
- F2 dispersion (the Im(λ) structure of the (0,1) coherence block): [`compute/RCPsiSquared.Core/Symmetry/F2W1DispersionPi2Inheritance.cs`](../compute/RCPsiSquared.Core/Symmetry/F2W1DispersionPi2Inheritance.cs).
- Q-anchor canonical table: [`docs/Q_REGIME_ANCHORS.md`](../docs/Q_REGIME_ANCHORS.md).
- Earlier slow-mode work (different Q regime, different sector concept): [`experiments/SLOW_MODE_R_PARITY.md`](SLOW_MODE_R_PARITY.md), [`simulations/slow_modes_r_parity.py`](../simulations/slow_modes_r_parity.py), [`compute/RCPsiSquared.Compute/LensAnalysis.cs`](../compute/RCPsiSquared.Compute/LensAnalysis.cs).
- External literature anchors: MEP 2016 (`arXiv:1606.09122`, exact Bethe-ansatz for periodic XX + dephasing); Žnidarič 2024 (`arXiv:2311.07375`, local-dephasing → diffusive 1/N²); Bortz-Stolze 2008 (`arXiv:cond-mat/0612382`, central-spin model for the star scaling family).
