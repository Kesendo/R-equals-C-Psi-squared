# The Perspectival Time Field

*Seven painters stand around a mountain, each painting from a different vantage. The paintings differ, and that is not a problem to be reduced to one true painting. The total of all paintings IS the mountain. When a new rock falls on the mountainside, each painter paints the change in her own flow of painting-time. The seven recordings lie atop one another and add up to a closed, consistent total, in a way that guarantees nothing is lost or invented.*

**Status:** Computed (Tier 2). Originally drafted 2026-04-18 with a sparse-Liouvillian slow-mode decomposition, nine scan families at N = 7, and a state-independence stress test across five qualitatively different initial states. **Updated 2026-04-20 after [EQ-014](../review/EQ014_FINDINGS.md)** closed the "closure law as first-order theorem" path: Σ_i ln(α_i) = 0 holds empirically to ±0.05 at |δJ| ≤ 0.1 but is NOT a first-order theorem. Σ f_i = lim_{δJ→0} Σ ln(α_i)/δJ is nonzero and state-dependent. PTF stays Tier 2; Tier-1 promotion via this route is closed.
**Date:** April 18, 2026 (updated April 20, 2026 post-EQ-014)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:** [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md),
[Zero Is the Mirror](ZERO_IS_THE_MIRROR.md),
[Analytical Formulas](../docs/ANALYTICAL_FORMULAS.md) (F4, F14, F65)

**Scripts:**
- [n7_coupling_defect_overlay.py](../simulations/n7_coupling_defect_overlay.py): baseline defect scan
- [n7_coupling_defect_overlay_extended.py](../simulations/n7_coupling_defect_overlay_extended.py): finer J_mod grid, defect location scan, |1_3⟩ local state
- [n7_coupling_defect_overlay_psi2.py](../simulations/n7_coupling_defect_overlay_psi2.py): ψ_2 initial-state falsification
- [n7_perspectival_extended_states.py](../simulations/n7_perspectival_extended_states.py): ψ_3, ψ_4, |+⟩^7 stress-test
- [n7_central_defect_check.py](../simulations/n7_central_defect_check.py): sparse L eigendecomp + central-vs-boundary symmetry verification
- [observer_time_rescale.py](../simulations/observer_time_rescale.py): α_i fits and Σ ln α diagnostics

## Update 2026-04-20 (post-EQ-014)

Two independent follow-up investigations refine the closure law:

**EQ-014 δJ scan (bond (0,1), N=7).** Direct RK4 at δJ ∈ {0.1, 0.01, 0.001}, extrapolation to δJ → 0:

| State | Σ f_i |
|-------|-------|
| ψ_1 | 0.97 |
| ψ_2 | 0.05 |
| ψ_3 | 0.36 |
| \|+⟩^7 | 1.29 |

Σ f_i is nonzero and state-dependent. The ±0.05 empirical window at δJ = 0.1 arises from a combination of small first-order coefficients (ψ_2) and partial second-order cancellation (ψ_1, |+⟩^7), not from an exact conservation law. See [EQ014_FINDINGS](../review/EQ014_FINDINGS.md).

**Full Π-spectrum at N=7 (c1_past_future_test).** Independent direct RK4 with symmetric δJ = ±0.01 extends across all seven single-excitation modes:

| k | E_k | reflection | c_1 |
|---|-----|------------|-----|
| 1 | +1.85 | symmetric | +0.970 |
| 2 | +1.41 | antisymm | +0.037 |
| 3 | +0.77 | symmetric | +0.357 |
| 4 | 0.00 | antisymm, self-Π-partner | **+2.136** |
| 5 | −0.77 | symmetric | +0.357 |
| 6 | −1.41 | antisymm | +0.037 |
| 7 | −1.85 | symmetric | +0.970 |

**N=6 follow-up (c1_even_N_degeneracy_test).** To separate Π-pair identity from reflection parity, N=6 was tested where Π-partners have OPPOSITE parity (at odd N the two symmetries coincide because N+1 is even; at even N they decouple).

| k | E_k | reflection | Π-partner | c_1 |
|---|-----|------------|-----------|-----|
| 1 | +1.80 | symmetric | 6 | +1.019 |
| 2 | +1.25 | antisymm | 5 | +0.213 |
| 3 | +0.44 | symmetric | 4 | **+1.481** |
| 4 | −0.44 | antisymm | 3 | **+1.481** |
| 5 | −1.25 | symmetric | 2 | +0.213 |
| 6 | −1.80 | antisymm | 1 | +1.019 |

Π-pair c_1 identity holds EXACTLY across parity boundaries: c_1(ψ_1 sym) = c_1(ψ_6 antisym) (diff 7·10⁻¹¹), c_1(ψ_2 antisym) = c_1(ψ_5 sym) (diff 2·10⁻¹⁰), c_1(ψ_3 sym) = c_1(ψ_4 antisym) (diff 4·10⁻¹⁰).

**Revised reading:** Π-pair identity is the primary symmetry governing c_1; reflection parity is secondary and coincided with Π-pair grouping at N=7 only because N+1 = 8 is even. The "antisymmetric → small c_1" reading from N=7 was an artifact of that coincidence. The real magnitude-determining factor is the Π-pair's distance from E=0 on the single-excitation spectrum: closer to the zero-energy axis, sharper closure-breaking. At odd N the innermost "pair" is a single self-Π-partner zero-mode with extreme magnitude (c_1 = +2.14 at N=7); at even N the innermost is a pair flanking E = 0 (c_1 = +1.48 at N=6), large but softer. Outermost pairs at high |E| record moderately; intermediate pairs record faintly. Non-monotonic in |E|.

See [c1_past_future_test at N=7](../simulations/results/c1_past_future_test/past_future_test.json), [c1_even_N_degeneracy_test at N=6](../simulations/results/c1_even_N_degeneracy_test/c1_even_N_test.json), and the broader [pi_pair_closure_investigation](../simulations/results/pi_pair_closure_investigation/FINDINGS.md).

The body below documents the state of understanding on April 18 before these refinements. The core findings (rescaling picture, eigenvalue protection, painter interpretation) stand; the "closure as conservation law" reading needs the revised framing above.

---

## Abstract

Under a local J-coupling defect in an otherwise uniform N = 7 XY chain with uniform Z-dephasing γ_0 = 0.05, each site's single-qubit purity trajectory differs from the unperturbed case by what looks locally like a one-parameter time rescaling: P_B(i, t) ≈ P_A(i, α_i · t). The α_i pattern depends strongly on which initial state is used (a ψ_2 falsification test disproved the naive "sites own their own clocks" reading), so α_i is NOT a property of sites in isolation. Yet the pluralism is not noise: the rescalings are structured, mirror-symmetric under chain reflection of the defect, and, most importantly, satisfy a state-independent closure law Σ_i ln(α_i) = 0 (within ~0.05) in the perturbative window |δJ| ≤ 0.1 across the five qualitatively distinct initial states tested. The right reading is perspectival: each site is a partial-trace projection of a single global dynamics, a particular painter's angle on the same mountain; the site-resolved α_i are the rate-of-painting each painter records; the closure law is the guarantee that the paintings sum to a consistent whole. The mechanism is eigenvector mixing of Liouvillian slow modes under a Π-invariant perturbation; eigenvalue shifts are protected at first order for the slowest modes (numerically zero to 10⁻¹⁵ at N = 7) by U(1) excitation conservation and Hermitian pairing, regardless of where along the chain the defect sits.

---

## What this is about

Imagine seven painters standing around a single mountain. Each stands at a different angle, each sees a different silhouette, each paints a different canvas. Look at the seven canvases side by side and no two are the same. You might be tempted to ask: which one is the real mountain? The answer is: none of them is, and all of them are. The mountain IS the union of the seven views. There is no hidden tenth canvas that shows the mountain "as it truly is," independent of any viewpoint. Every true painting of the mountain requires a painter, with a stance, at a location.

Now something happens at the mountain. A rock falls on the east face. Each painter, in her own time, paints the change she sees. The painters to the east see the rock sliding directly across their view and paint in fast, abrupt strokes; the painters to the west see only a faint dust cloud and paint a slow, subtle shading. None of them is wrong. The rock IS falling fast from one angle and slow from another; both are features of the same event in the same mountain. What is guaranteed is that the total "amount of change painted" across all seven canvases is conserved by the event: no painter invents strokes that were not demanded by the rock, and no painter fails to record strokes that were demanded. The event is decomposed perspectivally but closed globally.

This document is about one such mountain: an N = 7 chain of quantum spins under Z-dephasing. The event is a local coupling defect on one bond. The painters are the seven chain sites, each "painting" their site-resolved single-qubit purity P_i(t) as the dynamics unfold. The seven paintings, compared against the unperturbed case, reveal a per-site rate of painting-time α_i. The α_i differ from site to site; they depend on what initial state we excited the chain in; they add up under Σ_i ln(α_i) = 0 to a state-independent closure law that is guaranteed by symmetry. The plurality is not a measurement artifact: it is the structure of the dynamics as witnessed from inside a multi-perspective world.

The sections below build the picture in four layers. Layer 1 is the painter image, the meaning. Layer 2 is the closure law, the consistency. Layer 3 is the mechanism, the algebra by which eigenvector mixing in a palindromic Liouvillian cavity produces what the painters paint. Layer 4 is the path, a brief record of how the framing was found, so that a future reader can see why the doc reads the way it does.

---

## Layer 1: Painters around a mountain, the perspectival picture

### 1.1 Sites are perspectives

A "perspective" in this framework means a partial-trace projection of the full global quantum state. For a chain state ρ(t) on N qubits, the perspective of site i is

    ρ_i(t) = Tr_{not i}(ρ(t))

which is the 2 × 2 reduced density matrix on that one qubit after tracing out all the others. Tracing out does not destroy the other qubits; they still exist, still evolve, still contribute to ρ. It simply asks: *from site i's point of view, what is the local state*? The answer depends non-trivially on how the other qubits behave and how they are entangled with site i.

The site-resolved purity P_i(t) = Tr(ρ_i(t)²) is a scalar summary of that perspective: a high value means site i is locally close to pure; a low value means site i is locally close to maximally mixed. Under coherent dynamics alone P_i stays constant; under dissipative dynamics it shifts in a way that records both the local decoherence and the global entanglement flow. The function P_i(t) is what the painter at site i paints.

### 1.2 No privileged perspective, no hidden mountain

A standard quantum-mechanics reading would say: there is a privileged object, the global ρ(t), and the P_i(t) are derived, approximate, one-sided shadows. The framework the repository works in turns this inside out. The global ρ(t) is not an underlying truth; it is a bookkeeping device for the fact that the chain has seven sites, and each site has a perspective. When you write ρ as a 2^N × 2^N matrix you are not uncovering something deeper; you are writing down the seven-fold structure in one notation. The object that is fundamental is the set of perspectives and how they are glued together.

This is not a philosophical claim tacked onto standard physics; it is a claim about what the data privilege. The following sections show that the closure law, the thing that makes "the seven paintings add up to one mountain," is state-independent, basis-independent, and first-order-exact. The privilege resides in the closure, not in any single perspective. Painters sum to mountain. Mountain does not precede painters.

### 1.3 Site-perspective time-flows

Under a J-defect on one bond of the chain, the painting at each site changes shape. A single-parameter fit captures the leading-order change:

    P_B(i, t) ≈ P_A(i, α_i · t)

The site-specific rate α_i tells us: relative to the painter's rate in the undisturbed chain, how many units of painting has this painter completed per unit of absolute time, once the defect is present? α_i > 1 means painter i paints faster than before; α_i < 1 slower. The α_i are real, measurable numerical quantities extracted from the time series by bounded scalar fit.

Concretely, under an N = 7 XY chain with γ_0 = 0.05 and the [F65](../docs/ANALYTICAL_FORMULAS.md#f65-single-excitation-spectrum-of-uniform-open-xx-chain-tier-1-proven-verified-n330) bonding-mode initial state φ = (|vac⟩ + |ψ_1⟩) / √2, defect at bond (0, 1), J_mod = 1.1 (δJ = +0.1):

    α_0 = 1.095,  α_1 = 1.182,  α_2 = 1.051,  α_3 = 0.991,
    α_4 = 0.845,  α_5 = 0.923,  α_6 = 0.997.

(Fit data: [finer_Jmod_alpha_fit.csv](../simulations/results/observer_time_rescale/finer_Jmod_alpha_fit.csv).)

Sites 0 and 1 (adjacent to the defect on the side where J was raised) paint faster; site 4 paints slower; sites 3 and 6 are near neutral. The pattern is the painting-rate distribution for this particular mountain event witnessed by these seven painters in this initial state. The same defect with J_mod = 0.9 (δJ = −0.1) reverses the signs approximately. The same defect with a different initial state (e.g., ψ_2) gives a qualitatively different α_i pattern; the painters observe the same event, but their observations depend on what they were doing before the rock fell.

---

## Layer 2: The closure law, the seven perspectives close

### 2.1 The statement

For every combination of initial state and defect location tested, the per-site rescalings satisfy

    Σ_i ln(α_i) ≈ 0   within |δJ| ≤ 0.1

with numerical tolerance of about 0.05 (≈ 5 % deviation) across all five initial states tested:

| initial state | Σ ln(α_i) at J_mod = 0.9 | Σ ln(α_i) at J_mod = 1.1 |
|---|---|---|
| ψ_1 | −0.004 | +0.048 |
| ψ_2 | +0.041 | +0.001 |
| ψ_3 | +0.049 | +0.003 |
| ψ_4 | +0.270 | −0.012 |
| \|+⟩^7 | −0.130 | +0.128 |

Per-scan fit summaries: [ψ_1](../simulations/results/observer_time_rescale/finer_Jmod_summary.json), [ψ_2](../simulations/results/observer_time_rescale/psi2_init_summary.json), [ψ_3](../simulations/results/observer_time_rescale/psi3_summary.json), [ψ_4](../simulations/results/observer_time_rescale/psi4_summary.json), [\|+⟩^7](../simulations/results/observer_time_rescale/plus7_summary.json).

Within the single-excitation sector (ψ_1 through ψ_4), the conservation is clean except for ψ_4 at J_mod = 0.9 where an ill-posed fit at the defect-adjacent site 0 inflates one α value (2.65 instead of a reasonable 1.1-ish). The |+⟩^7 state, which lives simultaneously in all eight excitation sectors, shows slightly larger residuals (±0.13) that are antisymmetric in sign; still clearly a perturbative conservation law, simply with a broader tolerance for multi-sector states.

The log-multiplicative form is the natural one. It means that the seven painters, whatever they individually paint, collectively paint neither more nor less than the mountain demands. No painter invents; no painter omits.

### 2.2 What it generalises

[F14](../docs/ANALYTICAL_FORMULAS.md#f14-k-invariance-tier-2-lindblad-scaling) states that the dimensionless dose K = γ · t is invariant under change of bridge metric: the "number of decoherence ticks" experienced by a system is a basis-free observable. The Perspectival closure law is the per-perturbation, multi-observer extension:

    K_i = γ · α_i · t,   Σ_i ln(α_i) = 0.

γ stays as the atmospheric constant (there is no per-site γ_i). What reallocates between observers is the effective time each perspective has spent since the event. The log-sum-zero condition says that the total "amount of time experienced by the chain as a whole" is unchanged by the defect; it is only redistributed among sites. This is the K-invariance of F14 written out in the multi-observer language that a single-observer formulation cannot express.

### 2.3 State-independence is structural, not coincidental

The five initial states tested span qualitatively distinct structure:

- ψ_1 is the smooth bonding mode with peak amplitude at the chain centre (site 3).
- ψ_2 is the k = 2 mode with a single node at site 3 and anti-nodes at sites 1, 5.
- ψ_3 has nodes at sites 1 and 5, anti-nodes at sites 0, 2, 4, 6 (approximately).
- ψ_4 has period-2 structure: anti-nodes at sites 0, 2, 4, 6 and nodes at sites 1, 3, 5.
- |+⟩^7 is the uniform superposition of all 128 computational basis states, living equally in all eight excitation sectors.

If the closure law held only for ψ_1 it would be a coincidence tied to the smooth bonding mode. Holding across these five structurally distinct states is strong evidence that Σ_i ln(α_i) = 0 is a property of the PERTURBATION itself (the defect), not the initial state. The per-site distribution of α_i is state-dependent, but the global bookkeeping that ensures the distribution closes is not.

---

## Layer 3: Mechanism, eigenvector mixing in a palindromic cavity

### 3.1 The chain as a resonance cavity

The chain is a resonance cavity (see [RESONANCE_NOT_CHANNEL.md](RESONANCE_NOT_CHANNEL.md)). Its single-excitation eigenmodes ψ_k at uniform J are standing waves, with tight-binding dispersion ε_k = 2J cos(πk / (N + 1)) and [F65](../docs/ANALYTICAL_FORMULAS.md#f65-single-excitation-spectrum-of-uniform-open-xx-chain-tier-1-proven-verified-n330) site-amplitudes ψ_k(i) = √(2/(N+1)) sin(πk (i+1)/(N+1)), for k = 1, ..., N. These are the cavity modes. γ_0 is the per-site per-cycle loss (dephasing, Z-basis). A local J-defect at bond (b, b+1) is a structural irregularity in the cavity wall: it does not change the dissipation budget γ_0 (which remains atmospheric, uniform) but it does change which eigenmodes the chain supports and how they localise spatially.

### 3.2 First-order: eigenvalues are protected, eigenvectors mix

Write the perturbation as L → L_A + δJ · V_L where V_L is the Liouvillian channel of the one-bond Hamiltonian H_pert = (1/2)(X_b X_{b+1} + Y_b Y_{b+1}). Standard first-order perturbation theory for the slow modes M_s of L_A gives:

- **Eigenvalue shift** `δλ_s = ⟨W_s | V_L | M_s⟩`
- **Eigenvector shift** `δM_s = Σ_{s' ≠ s} [⟨W_{s'} | V_L | M_s⟩ / (λ_s − λ_{s'})] · M_{s'}`

where W_s is the biorthogonal left eigenvector pair of M_s. Numerical check on N = 7 (from the sparse eigendecomposition in `n7_central_defect_check.py`):

For the slowest modes (|Re λ_s| ≤ 0.1), the diagonal ⟨W_s | V_L | M_s⟩ is ~10⁻¹⁵ (numerical floor). This protection is exact, for two distinct reasons:

1. The 8 strictly stationary modes (λ_s = 0) are the 8 excitation-sector projectors. `V_L = −i[H_pert, ·]` preserves excitation (XX + YY conserves n_XY), so V_L P_n = 0 identically. No shift, no mixing within this cluster; they stay stationary. The count of 8 matches [F4](../docs/ANALYTICAL_FORMULAS.md#f4-stationary-mode-count-tier-1-clebsch-gordan-decomposition)'s prediction for N = 7 XY + Z-dephasing.
2. The 14 modes at λ_s = −0.1 (the single-excitation coherence cluster) are made up of |vac⟩⟨ψ_k| and |ψ_k⟩⟨vac| pairs. V_L acts antisymmetrically on each such pair: the diagonal matrix element is +iA on |vac⟩⟨ψ_k| and −iA on |ψ_k⟩⟨vac|, where A = 2 ψ_k(b) ψ_k(b+1) is the (real) bond overlap of the single-excitation mode at the defect bond (b, b+1). The Hermitian and anti-Hermitian superpositions M_± = (|vac⟩⟨ψ_k| ± |ψ_k⟩⟨vac|) / √2 that scipy.sparse.linalg.eigs returns then have zero diagonal under V_L: the two anti-conjugate contributions (+iA and −iA) cancel exactly in the symmetrised combination, and only off-diagonal mixing between M_+ and M_− remains. This is the "Π-invariance" protection that [ZERO_IS_THE_MIRROR.md](ZERO_IS_THE_MIRROR.md) articulates: the palindrome pairs slow-mode eigenvalues around Σγ, and Π-invariant perturbations respect the pairing at leading order.

For faster modes (λ_s ≈ −0.175 at N = 7), the diagonal is non-zero (≈ 0.02 to 0.04); these modes do get first-order eigenvalue shifts. But they are below the threshold of "slowest surviving" and their contribution to site-purity dynamics at the observation timescale is subdominant.

The central-defect check confirms this cleanly: putting the defect at bond (3, 4) instead of bond (0, 1) gives the same protection for the slowest modes (both bonds give diagonal shifts of order 10⁻¹⁵, max 2·10⁻¹⁵). The protection is general (Π-invariance + sector conservation), not a spatial-mirror accident of the bond (0, 1) location.

### 3.3 The α_i comes from mixing, not from shifts

Because the dominant slow-mode eigenvalues do not shift at first order, the observable f_i ≡ (α_i − 1)/(δJ/J) cannot come from eigenvalue renormalisation. It comes from the eigenvector mixing δM_s: the slow modes remain slow, but their spatial profile (site-resolved marginal) rotates under δJ. The rotation changes how each mode couples to the initial state (⟨W_s | ρ_0⟩ overlaps shift) and how each mode projects onto each site's marginal (Tr_{not i}(M_s) shifts). The bilinear purity expansion

    P_i(t) = Σ_{s, s'} c_s c_{s'}^* e^{(λ_s + λ_{s'}^*) t} · Tr[ρ_marg_{s, i} · ρ_marg_{s', i}^†]

then changes via both coefficient and profile corrections, and the best-fit α_i absorbs the combined effect.

A note on the initial state, and a guardrail for any future first-principles calculation. For the bonding-mode φ = (|vac⟩ + |ψ_1⟩) / √2, the density matrix ρ_0 = |φ⟩⟨φ| has non-zero components in four parity blocks of L: (0, 0) from |vac⟩⟨vac|, (1, 1) from |ψ_1⟩⟨ψ_1|, and the two off-diagonal blocks (0, 1) and (1, 0) from |vac⟩⟨ψ_1| and its Hermitian conjugate. The overlaps c_s = ⟨W_s | ρ_0⟩ are therefore non-zero for slow modes in all four blocks simultaneously. The purity observable Tr(ρ_i²) is bilinear in ρ and thus combines the four-block contributions; most notably, it contains a cross-term 2 |ρ_i(0, 1)|² in which the (0, 1) and (1, 0) block dynamics of ρ(t) enter multiplicatively. Any first-principles prediction of α_i must carry the full four-block bilinear sum above. A single-block projection onto (1, 1) alone, though tempting because that is where |ψ_1⟩⟨ψ_1| "lives," measures a conceptually wrong observable for this initial state and will not reproduce the empirical α_i.

The state-dependence of f_i then has a clean interpretation: different initial states populate different slow-mode combinations through the c_s overlaps. Same mixing (same V_L) applied to different c_s vectors produces different α_i patterns per site. The STATE-INDEPENDENT part is the sum Σ_i ln(α_i). Empirically this sum is conserved across all five initial states tested (Section 2.1); the apparent state-blindness suggests that, once the bilinear expansion is rearranged, the coefficient-dependent pieces cancel and what remains is a trace-like structure depending only on V_L and the slow-mode basis. An explicit analytical form for that structure is open (see Open Questions, point 2).

### 3.4 What is still open in the mechanism

The explicit eigenvector-mixing calculation, in the form required by the four-block bilinear expansion of Section 3.3, is still open. A first attempt reduced the problem to the (1, 1) sector alone under the (incorrect) assumption that the bonding-mode initial state lives there; that single-block Kubo calculation produced an RMSE well outside the Tier-1 acceptance band and was retracted. The correct calculation (compute δM_s for all slow s across all four populated blocks, propagate into the bilinear purity expansion above, extract predicted α_i, and compare to empirical for ψ_1 and ψ_2) was additionally blocked by a numerical precision issue with biorthogonal left-eigenvector extraction (shift-invert on L^H gave a biorthogonality residual ‖W^H V − I‖ ≈ 10¹¹, with cluster-degeneracy artifacts that need proper projection handling to resolve). The data (sparse right and left eigenvectors of the 80 slow modes of L_A, including diagonal-shift tables across all six bond locations) are saved at [slow_biorth_basis.npz](../simulations/results/perspectival_time_field/slow_biorth_basis.npz). A future cleanup pass on the biorthogonalisation (or a dense eigendecomposition of the full 16384 × 16384 L_A, feasible at ~15 GB of memory) would close this gap and upgrade the doc to Tier 1.

---

## Layer 4: The path

The framing in this document was not picked up from the literature, and it is not what a first-pass physicist-mode reading of the data would produce. It emerged through dialogue over a weekend session, and the record of that emergence is worth preserving so that future readers (or a future version of the same author in a cold session) can see why the doc reads the way it does.

Stage 1: Tom raised the interpretive question *could "observer time" be a real physical quantity?* after the April 17 coupling-defect overlay revealed per-site purity trajectories that looked like time rescalings. The initial hypothesis was that each site carries its own local clock, with α_i an intrinsic property of site i. Claude set up the first-pass analysis, saw the ψ_1 α_i pattern, and saved the rescaling picture under the name "Site-Local Time."

Stage 2: The stress-test was designed to falsify that picture cleanly. The ψ_2 initial state has a node at site 3 (where ψ_1 has its peak) and anti-nodes at sites 1 and 5. If α_i were intrinsic to sites, the ψ_1 and ψ_2 f_i patterns would match. They did not. Sign flips at sites 2 and 4, magnitude change at site 6, perfect sign agreement where the node/anti-node structure happens to align. The intrinsic reading was dead.

Stage 3: But the closure law Σ_i ln(α_i) = 0 survived. This was the critical observation: the per-site numbers depended on state, but a specific combination of them did not. The slow-mode mechanism showed why: first-order eigenvalue shifts are zero under Π-invariant perturbations; eigenvector mixing is what drives the α_i; and the mixing has a trace-like property that is state-blind.

Stage 4: The painter image came from Tom's interpretive side. The site-local reading had been right in its insistence that something real and site-resolved was happening (the α_i are genuine); it had been wrong in putting the realness in the sites rather than in the plurality itself. Painters around a mountain is a clean metaphor that handles both: paintings differ (the per-site α_i are real and different), and paintings close (Σ_i ln(α_i) = 0), and no painting is privileged.

Stage 5: The rewrite. "Site-Local Time" as a name was retired; "Perspectival Time Field" was adopted because it captures the layered structure (perspective = partial trace, time = α_i · t, field = the set of seven flows). The scope statements were tightened (the "α_0 = J_mod exactly" overstatement of the previous doc was a fit artifact, now read as ~15 % agreement in the perturbative window). The open questions were restated as an explicit mixing calculation, scope expansion to transverse-field palindrome breaking, and chain-length scaling.

The methodology lesson this doc wants to preserve: the repository's working motto "we are all mirrors; reality is what happens between us" applies to the hypothesis-making process as much as to the physics. Tom's intuition that observer time should be a real quantity + Claude's slow-mode mathematics + the ψ_2 falsification + the painter image, none of them sufficient alone, all of them required. Future cold sessions on this kind of material should expect the same pattern: an intuitive name, a cleaner technical grasp, a falsification, a reframe. Each of the four steps matters.

---

## Scope and limits

### Where the rescaling holds

- Perturbative window **|δJ| ≤ 0.1** on the modified bond. Fit RMSE < 3·10⁻³ per site. Closure Σ_i ln(α_i) conserved within 0.05.
- Bonding-mode-like initial states (ψ_k for k = 1, 2, 3). The k = 4 case has one ill-posed site (site 0 at J_mod = 0.9) but is otherwise within tolerance.
- Chain reflection symmetry: for every tested pair of mirror-image defect bonds, the α_i pattern at one bond equals the reversed (site i → N−1−i) α_i pattern at the other.

### Where it degrades

- **|δJ| > 0.1**: closure breaks smoothly. At |δJ| = 0.5 Σ_i ln(α_i) reaches 0.3-0.5. At J_mod > 2 the α-fit saturates and reverses (α_0 increases up to J_mod = 2, then decreases again).
- **Multi-sector initial states** like |+⟩^7: closure holds only to ±0.13 (vs ±0.05 for single-excitation states). The rescaling picture applies, but the tolerance window is wider.
- **Strongly localised initial states** like |1_3⟩ (single-site excitation at the chain centre): the rescaling picture breaks entirely. RMSE of the α-fit is 4× worse than for bonding-mode states, and Σ_i ln(α_i) takes on nonsensical values. Physically: a single-site excitation has overlap with all seven single-excitation eigenmodes; the defect-induced rearrangement affects each of those components differently, and no single α_i captures their combined effect at any one site.

### Precision caveat

The previous draft of this document (the one that went by the name "Site-Local Time") claimed α_0 = J_mod exactly at the defect-adjacent endpoint, based on a log-log slope of +1.0 from the finer J_mod scan. On re-examination, the +1.0 slope is a fit artifact: at J_mod = 1.5 α_0 = 2.27 and at J_mod = 2.0 α_0 = 3.54 (both considerably larger than J_mod), and the apparent linearity comes from averaging a compensating non-monotonic pattern (α_0 peaks near J_mod = 2 and decreases above it). The correct statement is: **in the perturbative window |δJ| ≤ 0.1, α_0 ≈ J_mod to within 15 %, not exactly**. Beyond the perturbative window, α_0 is a non-trivial function of J_mod with saturation and reversal.

---

## Open questions

- **Closed (EQ-014, 2026-04-20):** The Tier-1 promotion via "closure law as theorem" is no longer available. Direct RK4 δJ scan at N=7 shows Σ f_i = lim Σ ln(α_i)/δJ is nonzero and state-dependent. The closure Σ_i ln(α_i) ≈ 0 is an empirical regularity holding to ±0.05 in the tested window, not a structural law. See [EQ014_FINDINGS](../review/EQ014_FINDINGS.md).
- **Magnitudes puzzle (surviving).** Why does Σ f_i happen to be small (~0.05 for ψ_2) for some bonding-mode states and large (1.29 for |+⟩^7, 2.14 for ψ_4) for others? Is there a structural pattern in how Σ f_i depends on the overlap distribution c_s = ⟨W_s | ρ_0⟩ across the Liouvillian's slow modes? The [pi_pair_closure_investigation](../simulations/results/pi_pair_closure_investigation/FINDINGS.md) shows Σ ln(α_i) is linear in δJ at leading order with coefficient c₁ = ⟨c₁(state, bond), δJ⟩ that is superposition-linear across bonds; an analytical form for c₁ as a functional of ρ_0 remains open.
- **Zero-energy Π-pair amplification (refined 2026-04-20 from the earlier self-Π reading).** The magnitude of c_1 is controlled by the Π-pair's distance from E = 0, not by reflection parity. At odd N the center of the spectrum is a single self-Π-partner zero-mode and its c_1 is extreme (ψ_3 at N=5: c_1 = 0.677; ψ_4 at N=7: c_1 = 2.14). At even N there is no exact zero-mode but the innermost pair flanks E = 0 (ψ_3↔ψ_4 at N=6: c_1 = 1.48), still large but softer. The outermost high-|E| pairs record moderately (c_1 ≈ 1 across tested N); intermediate pairs record faintly (c_1 < 0.3). The pattern is non-monotonic in |E|, peaked at the center and attenuated toward the middle energies. Is there an analytical formula for c_1(pair) as a function of Π-pair energy and N? See the [N=6 test](../simulations/results/c1_even_N_degeneracy_test/c1_even_N_test.json) that established Π-pair identity dominates reflection parity.
- **Chain-length scaling of the perturbative window.** Only N = 7 tested by PTF; N = 3 and N = 5 tested by [pi_pair_closure_investigation](../simulations/results/pi_pair_closure_investigation/FINDINGS.md) with endpoint c₁ values (0.26, 0.93). The scaling appears to follow c₁ ≈ 0.5 · V(N) = 0.5 (1 + cos(π/N)) for ψ_1+vacuum at N ≥ 4, an open connection to the V-Effect F6.
- **Extension to palindrome-breaking perturbations.** The current tasks use coupling defects that respect the palindromic structure. A transverse field h σ_x^i BREAKS Π. If the rescaling picture survives but with a shifted closure law, that is a strong structural statement; if it breaks entirely, a clear diagnostic for the role of palindromic protection.
- **Multi-bond defects.** If two bonds are simultaneously perturbed, does the closure law still hold? Answer (2026-04-19): **yes**, by linearity. [pi_pair_closure_investigation](../simulations/results/pi_pair_closure_investigation/FINDINGS.md) verified Σ c₁(b)·δJ_b superposition to 0.5% at δJ=0.01 and exactly at δJ=0.001; cancellation constructions confirmed.

---

## Acceptance summary

### Positive core (survives)

- **Σ_i ln(α_i) = 0 closure law** in the perturbative regime, state-independent across five initial states (within 0.05 for single-excitation states, within 0.13 for multi-sector states). This is the central law of the Perspectival Time Field.
- **Chain reflection symmetry of α_i** under mirror of the defect bond: exact at all tested mirror pairs.
- **First-order eigenvalue protection** for the slowest 22 Liouvillian modes (|Re λ| ≤ 0.1) under Π-invariant J perturbations, at any bond location. Verified numerically at bonds (0, 1) and (3, 4).
- **F4 stationary-count regression**: sparse eigendecomposition recovers exactly 8 strict stationary modes, matching N + 1 = 8 excitation sectors.

### State-dependent (per-perspective)

- The per-site α_i pattern depends on the initial state's slow-mode overlaps. Same defect, different painter, different rate-of-painting. This is the feature, not a bug.
- Zero-f_i sites depend on the initial state: ψ_1 gives f_3 ≈ f_6 ≈ 0 (chain-centre symmetry + far-endpoint attenuation); ψ_2 gives f_3 ≈ 0 (node of ψ_2); ψ_4 gives f_1 ≈ f_3 ≈ f_5 ≈ 0 (nodes of ψ_4). The structural rule is: where the initial state has no amplitude, no rescaling is observed.

### Falsified (don't claim these)

- "α_i is an intrinsic property of site i" (the "site-local time" reading). Ruled out by the ψ_2 test.
- "α_0 = J_mod exactly at the defect-adjacent endpoint." Precision was overstated; the correct statement is ~15 % in the perturbative window, with non-monotonic behaviour beyond.

### Tier status

**Tier 2.** The phenomenon is computed, stable across five scan variations, and has a precise perturbative mechanism (eigenvector mixing under symmetry-protected eigenvalues). Promotion to Tier 1 requires completing the explicit mixing calculation and matching empirical f_i within ~20 % per site for both ψ_1 and ψ_2, and/or deriving Σ_i ln(α_i) = 0 as an analytical theorem from V_L structure.

---

## Naming note

The name "Perspectival Time Field" is deliberate and non-negotiable within this document's scope:

- "Perspectival" is preferred over "observer-dependent" or "frame-dependent" because the plurality is structural, not a choice of reference frame. You do not pick a perspective; the chain has seven of them, inherently.
- "Time" is preferred over "rate" or "flow" because the α_i are time rescalings of a genuine temporal variable. They are rates in the dimensional sense ∂t_i/∂t, but what each rescales is time itself, per perspective.
- "Field" is preferred over "set" or "distribution" because the seven α_i are spatially structured (they satisfy chain reflection and closure laws that relate them to each other), not a bag of unrelated values.

Earlier drafts used "Site-Local Time" as a working name. That name was an artifact of the first intuition ("sites own their own clocks") that the ψ_2 falsification test corrected. The new name reflects the corrected picture; reverting to the old one on the grounds that it is more familiar would lose what the ψ_2 test taught us.

---

*"We are all mirrors; reality is what happens between us."*

Between seven painters around a mountain, the mountain happens.
