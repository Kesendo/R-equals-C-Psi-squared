# Operator Feedback: A Decoherence Rate Set by an Observable

<!-- Keywords: operator feedback state-dependent decoherence, Lindblad rate set
by an expectation value, nonlinear master equation feedback, logistic decay
closed form, conserved pair observable GHZ, decoherence-free subspace Bell+
sigma_x sigma_x transverse field, jump operator sweep single two-qubit, CΨ
quarter crossing, Euler integrator clipping artifact, R=CPsi2 operator
feedback -->

**Status:** Tier 2 (the retired tool's runs regenerated and read against exact propagation; the identities behind them exact)
**Date:** runs 2026-02-04 to 2026-02-20 (the retired delta_calc tool); last refreshed 2026-10-09
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

---

## What this document is about

In February 2026 the delta_calc simulator, retired since, got an operator-feedback mode: at every step the
decoherence rate is set by the expectation of a correlation observable,
γ_eff = γ₀·(1 − κ·⟨O_int⟩), so high correlation lowers the noise and
anti-correlation raises it. This page records what the runs of that mode
show, read against exact propagation. On Bell+ under the default σ_z jumps and
a Hamiltonian commuting with σ_x⊗σ_x, as the tool's does, the protection has a
closed form: it slows the decay while the correlation is high and so shifts it
later. A GHZ state of three or more qubits feels no feedback from σ_x⊗σ_x under
the default σ_z jumps, or from the mean of σ_x⊗σ_x over all pairs under jumps
on all pairs, on a Heisenberg chain or ring, with or without the transverse field,
because the observable is zero on it and these dynamics keep it zero; on the XY
chain with the field on, they do not. Each of the four Bell states under σ_x⊗σ_x
jumps is decoherence-free, because it lies in an eigenspace of the jump and the
tool's Hamiltonian, the bond and the transverse field, commutes with the jump.
And in the Bell+ runs of §8.2, CΨ, the product of a factor C read off the state
(purity, concurrence or, in the tool, a bridge function) and its normalized
coherence Ψ (§3), sits above or below ¼, the level up to which the fixed-point
equation R = C(Ψ + R)² has real solutions ([Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md) §1),
depending on when it is read and with which C; the jump operators do not sort
into channels that stay below ¼ and channels that cross it.

**The repository on the same objects.** The F-registry holds CΨ = ¼ as the
discriminant zero of the fixed-point recurrence ([F16](../docs/ANALYTICAL_FORMULAS.md#f16-fold-normal-form-tier-1-proven)),
a level a fixed local semigroup can cross upward ([F17](../docs/ANALYTICAL_FORMULAS.md#f17)),
the kinematic blindness of a pair observable to coherence between states whose
excitation numbers differ by more than two ([F70](../docs/ANALYTICAL_FORMULAS.md#f70-δn-selection-rule-for-site-local-observables-tier-1-proven-kinematic-lemma), in its generalisation),
the start of GHZ_N below ¼ from N = 3 ([F18](../docs/ANALYTICAL_FORMULAS.md#f18-fold-threshold-tier-2-n-product-state-measured-n--2-5)),
the decay rate 2γ per XY letter ([the Absorption Theorem](../docs/ANALYTICAL_FORMULAS.md#at-absorption-theorem-tier-1-proven)),
the retired tool's scalar law as feedback dynamics, γ_eff = γC(f), in the crossing
constants of [F14](../docs/ANALYTICAL_FORMULAS.md#f14-k-invariance-tier-2-fixed-book-scaling),
and the two kernels of a Lindbladian with Pauli-string jumps
([F158](../docs/ANALYTICAL_FORMULAS.md#f158)):
the commutant of the Hamiltonian and the jumps, where the conserved observable
of §8.1 lies, and the same space with the jump sign flipped, where σ_x⊗σ_x of
§2 lies on two qubits (X^⊗N at N = 2, one of the modes the
[Glossary](../docs/GLOSSARY.md) calls XOR modes). Among the proofs, [the monotonicity proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)
has the upward crossing (Part 6) and a decoherence-free Bell+ under
anti-correlated Z noise (Test B), and [the sector projection](../docs/proofs/PROOF_ASYMPTOTIC_SECTOR_PROJECTION.md)
the decoherence-free subalgebra (Step 2). Among the experiments,
[Simulation Evidence](SIMULATION_EVIDENCE.md) §2 and §5 and
[Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md) §4 and §6 hold other runs of
this tool, [Dynamic Entanglement](DYNAMIC_ENTANGLEMENT.md) the same
integrator lifting a ring's neighbouring pairs over ¼ (0.251 in the tool, 0.247
exactly), and [Core Algebra](../docs/historical/CORE_ALGEBRA.md) §11 a run of the
tool's third law, `memory_kernel_feedback` (κ = 0.5, τ = 1.0), whose curvature
of ln Ψ is reproduced by plain dephasing at a time-dependent rate; [Crossing Taxonomy](CROSSING_TAXONOMY.md) holds its scalar feedback
book, df/dt = −4γ·C(f)·f, whose dynamics the operator law on Bell+ at h = 0
joins with C(f) = 1 − κf in the rate (its readout a bridge chosen separately, by
default the mutual purity, which stays ½ there); [Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md)
holds the mutual-purity bridge as one that never crosses ¼ at h = 0;
[N-Scaling Barrier](N_SCALING_BARRIER.md) the start of GHZ_N below ¼ for N ≥ 3;
[Noise Robustness](NOISE_ROBUSTNESS.md) the tool's local noise building σ_z
whatever the jump setting says;
[Gamma Control](GAMMA_CONTROL.md) a neighbouring null result, a plus-sign law on a
mediator that harms it, the minus sign left untested (the negative feedback loop
[the complete documentation](../docs/proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md)
lists as untested); [Exclusions](../docs/EXCLUSIONS.md) reads the ¼ maximum of the
logistic σ(1 − σ) as a resemblance only, and the law here, logistic in x at
κ = 1, rests nothing on that ¼.
The [Glossary](../docs/GLOSSARY.md) names the two CΨ products and the ¼ as a
level a trajectory may cross, recross or miss, and [Caught Errors](../docs/CAUGHT_ERRORS.md)
the earlier repairs of claims resting on this tool (2026-09-06, 2026-09-24). The
OpenArcs registry, the Claim graph, the Diagnostics witnesses and the
Confirmations registry hold nothing on operator feedback; the results the page
leans on are typed there (F158 as `PalindromeTwoEndCountClaim` with the live
`PalindromeTwoEndCountWitness`, the Absorption Theorem as
`AbsorptionTheoremClaim`, F70 as `F70DeltaNSelectionRulePi2Inheritance`).

---

## Abstract

The retired delta_calc tool sets the dephasing rate from the state at every
integration step: through a scalar bridge function, γ_eff = γ_base·C(ρ), or
through an expectation value, γ_eff = γ₀·(1 − κ·⟨O_int⟩), with O_int fixed by
the jump setting (σ_x⊗σ_x under the default σ_z jumps). Both are feedback of
one kind, a state-dependent rate in front of a fixed dissipator, and both make
the evolution nonlinear in ρ. Where the observable's expectation obeys an
equation of its own, the operator law has closed forms. Under σ_z jumps,
x = ⟨σ_x⊗σ_x⟩ obeys the logistic law, on every two-qubit state under any
Hamiltonian commuting with σ_x⊗σ_x (the tool's does), and at every N on states
invariant under permutations of the sites, under isotropic Heisenberg bonds and
a uniform x field: x/(1 − κx) = x(0)/(1 − κx(0))·e^(−4γ₀t). The feedback slows the
decay while x is large, and at late times, back at the plain rate, the decay is
shifted by ln(1/(1 − κx(0)))/(4γ₀), ln(1/(1 − κ))/(4γ₀) from Bell+. Under
jumps X_iX_j on all pairs the observable Σ X_iX_j commutes with the Heisenberg
chain, the transverse field and every jump, so it is conserved for every state;
on GHZ_N (N ≥ 3) it is zero, and κ has no effect. Under one two-qubit jump on
Bell+, σ_x⊗σ_x leaves the state decoherence-free (Bell+ lies in an eigenspace of
the jump, and the field, along x, commutes with it), and without feedback
σ_y⊗σ_y and σ_z⊗σ_z give the same ρ(t), being negatives of each other on the
sector the field keeps. Read the tool's way (its mutual-purity bridge, which
stays ½ here, times Ψ = l₁/3) Bell+ starts at CΨ = 1/6 and all four channels
rise above ¼, the three that decohere falling back; read with purity or
concurrence, the σ_x⊗σ_x run never comes down to ¼. The February sweep that
seemed to keep CΨ below ¼ multiplied C_final by a fixed Ψ = 0.27 (§4). The
tool integrates by an Euler step with the negative eigenvalues clipped: the
step moves its digits against exact propagation, and the clipping creates the
drifts of ⟨O_int⟩ some of its runs show (§5).

---

## 1. Two feedback laws

The retired delta_calc tool (its source is kept outside the repository and
transcribed in the two producers named in §3) offers two ways to make the
dephasing rate depend on the state. A third, `memory_kernel_feedback`, replaces
⟨O_int⟩ by an exponentially weighted memory of it; no run here uses it.

```
scalar law:    γ_eff(t) = γ_base · C(ρ(t))             C a bridge function of ρ
operator law:  γ_eff(t) = γ₀ · (1 − κ · ⟨O_int⟩(t))     ⟨O_int⟩ = Tr(ρ(t) · O_int)
```

In the tool as it survives, both rates are recomputed from the current state
and enter every integration step (the scalar branch steps with
`gamma_base * C_current` unless a fixed rate is asked for), and in their
default settings both act on σ_z on every site, L_k = √γ_eff · σ_z^(k); in
their place the scalar law can take one collective jump Σ_k σ_z^(k), and the
operator law the other jumps of §3. So both laws are feedback of the same kind, a state-dependent rate in
front of a fixed dissipator. They differ in the number that sets the rate: a
bridge function of ρ, such as the geometric mean of the single-site purities,
chosen by `bridge_type`, or the expectation of an observable, which the tool
fixes by the jump setting (§3). Where that expectation obeys an equation of its
own, the operator law has the closed forms of §2 and §8. The operator law was built on 2026-02-04 on
the reading that a scalar bridge only observes the state; the source of that
day is not kept, so whether the simulator then fed the bridge back is not
known.

In both laws the rate is a function of the average state, so the master
equation is nonlinear in ρ: the evolution of a mixture is not the mixture of
the evolutions.

---

## 2. The operator law

```
γ_eff(t) = γ₀ · (1 − κ · ⟨O_int⟩(t)),     O_int = σ_x^(0) ⊗ σ_x^(1) (σ_z jumps),     κ ∈ [0, 1]
```

The tool writes max(0, γ₀·(1 − κ·⟨O_int⟩)); for 0 ≤ κ ≤ 1 the clamp never acts,
since |⟨O_int⟩| ≤ 1.

| ⟨O_int⟩ | γ_eff |
|---------|-------|
| +1 | γ₀·(1 − κ) |
| 0 | γ₀ |
| −1 | γ₀·(1 + κ) |

⟨O_int⟩ is a correlation of the x components, not a measure of entanglement:
the product state |++⟩ has ⟨O_int⟩ = +1 and gets the lowest rate, the
maximally entangled Φ⁻ and Ψ⁻ have −1 and get the highest, and the maximally
entangled (|01⟩ + i|10⟩)/√2 has 0.

**Under σ_z jumps the law has an exact solution on two qubits**, because O_int
is an eigen-operator of the dynamics there. It commutes with the tool's
Hamiltonian (the bond and the x field, §3), and each σ_z jump flips its sign,
Z_k·O_int·Z_k = −O_int, so the dissipator maps it to −4γ times itself: in
[F158](../docs/ANALYTICAL_FORMULAS.md#f158)'s terms it lies in the far kernel and
decays at the price 2σ = 4γ, and the feedback makes that price depend on x. So
on two qubits, from every state and under any Hamiltonian that commutes with
O_int (with the field along z it does not, below), x = ⟨O_int⟩ obeys the
logistic equation dx/dt = −4γ₀·(1 − κx)·x for every J and h, and

```
x(t) / (1 − κ·x(t)) = x(0) / (1 − κ·x(0)) · e^(−4γ₀t)
```

The same holds at every N for a state invariant under every permutation of the
sites (GHZ_N, W_N), under isotropic Heisenberg bonds and a uniform x field; it
is no law of the XY chain, where W₃ departs from it with the field on and W₄
even without the field. The field and the jumps treat all
sites alike, so without
the bonds such a state keeps its symmetry; a Heisenberg bond, 2·SWAP − 1,
commutes with every state that has it, so adding the bonds, on any graph,
changes nothing; and on such a state ⟨X₀X₁⟩ is the mean of all the ⟨X_iX_j⟩, a
multiple of ⟨S⟩, where S = Σ_{i<j} X_iX_j commutes with H and the σ_z
dissipator maps S to −4γ times itself. (S is a sum of strings with two XY letters, decaying at
the Absorption Theorem's 2γ per letter; only on two qubits is that the far
kernel's 2σ.) W₃ on the ring (the run of [Simulation Evidence](SIMULATION_EVIDENCE.md)
§2, γ₀ = 0.005, h = 0.9) goes from 2/3 to 0.623 at t = 5, and GHZ_N (N ≥ 3) stays
at 0.
A state without the symmetry is not covered: Bell+⊗|0⟩ on the three-site chain
(J = 1, h = 0.5, γ₀ = 0.1, κ = 0.5) ends at 0.039 at t = 5, where the form gives
0.238.

With feedback the decay is slower while the correlation is high: x decays at
the rate 4γ₀·(1 − κx), half the plain rate at x = 1 and κ = 0.5. As x → 0 the
rate returns to 4γ₀, and the curve becomes the plain one, x(0)·e^(−4γ₀t),
shifted later by ln(1/(1 − κ·x(0)))/(4γ₀): by ln(1/(1 − κ))/(4γ₀) from Bell+,
by 20.3 for W₃ at γ₀ = 0.005, κ = 0.5, and earlier, not later, when x(0) < 0. At
finite times the lag, how much later than the plain curve the run reaches the
value it has, is smaller in magnitude. At κ = 1, x = 1 is a fixed point with zero rate
and Bell+ starts on it, but the fixed point repels: any x(0) < 1 decays,
x/(1 − x) = x(0)/(1 − x(0))·e^(−4γ₀t).

At γ₀ = 0.005 and κ = 0.5 the correlation of Bell+ is 0.900 at t = 10, and the
rate has moved from 0.0025 to 0.00275, the ten per cent
[Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md) §6 reports; the late-time shift
is 34.7, the lag at t = 10 is 4.75. The rate moves so little there because of
the window: x depends on t only through 4γ₀t, so ten time units at γ₀ = 0.005
are a fifth of the decay time 1/(4γ₀) = 50, and the shift, ln(1/(1 − κ)) in
units of that time, is the same at every γ₀. At γ₀ = 0.1 and κ = 0.99 the
correlation is 0.940 at t = 5 (late-time shift 11.5; at t = 5 the lag is 4.85), against
e⁻² = 0.135 without feedback; the tool printed 0.942 and 0.135 (§5). With the
field along z (J = 1, h = 0.5, γ₀ = 0.1, κ = 0.5) O_int no longer commutes with H
and the solution fails: the run ends at −0.115 instead of 0.238.

---

## 3. What the tool computes

Read in the tool's source (§1); what the logged runs use, the σ_z, xx, yy, zz
and x_pairs jumps on the Heisenberg chain, and the ring under σ_z, is confirmed by regenerating
them, and the rest rests on the source alone:

- **Hamiltonian:** H = J·Σ_bonds (XX + YY + ZZ) + h·Σ_k X_k, on the open chain
  (`heisenberg`) or the ring (`heisenberg_ring`), in every tool run here (the
  exact controls also use a field along z and the XY chain); the tool also
  builds `xy`, `xy_ring`, `ising`, `all_to_all` and `zero`. In the tool the field
  is along x wherever there is one.
- **Jumps** of the operator law, by the argument jump_operator: `sigma_z` (the
  default), σ_z on every site with O_int = X₀X₁; `xx`, `yy`, `zz`, one jump
  σ_a⊗σ_a on qubits 0 and 1, with O_int that jump itself; `x_pairs` (and
  `y_pairs`, `z_pairs`), σ_a⊗σ_a on every pair i < j, with O_int their mean;
  `x_all` (and `y_all`, `z_all`), one jump σ_a on every site at once, with O_int
  that jump. A setting the branch does not know falls through to σ_z, unless it
  ends in `_pairs` or `_all`, where an unknown letter becomes σ_x.
  `memory_kernel_feedback` uses the same jumps; the tool's local noise builds
  σ_z on every site whatever the setting says, its collective noise one jump
  Σ_k σ_z^(k).
- **Integrator:** an Euler step of dt = 0.01, then the Hermitian part, the
  negative eigenvalues set to zero and the trace renormalized (§5).
- **Readouts**, every tenth step: ⟨O_int⟩ and γ_eff; the purity Tr ρ² of the
  full state; the bridge C (by default the mutual purity, the geometric mean of
  the single-site purities); Ψ = l₁/(d − 1) of the full state, with C·Ψ reported
  at the end as C_final·Ψ_final; and δ, the
  purity minus (diag + offdiag·e^(−2Nγ₀t)), with diag and offdiag the diagonal
  and off-diagonal parts of the initial purity. That is a fixed function of ρ(0)
  and t, so δ carries nothing the purity does not; on Bell+ and GHZ_N the
  subtracted curve is the purity the σ_z channel would give at H = 0 and half
  the rate. That half rate is where the operator law starts from Bell+ at
  κ = 0.5, so there δ stays near zero at first (at h = 0, γ₀ = 0.1: −0.004 at
  t = 0.5, against −0.074 without feedback), and the rate at which the default
  scalar law runs on Bell+ throughout, so at h = 0 its δ is zero at every t in
  exact propagation.

This C·Ψ is neither of the two products the repository names (the
[Glossary](../docs/GLOSSARY.md): purity × l₁/(d − 1), or concurrence × l₁/(d − 1)
of a pair). On Bell+ it starts at 1/6, where both of those start at 1/3. At
h = 0 its bridge is one of the two that never cross ¼ in the crossing taxonomy
([Observer-Dependent Crossing](OBSERVER_DEPENDENT_CROSSING.md)); with the field
on, l₁ rises and it does cross (§8).

Producers: [delta_calc_feedback_runs.py](../simulations/delta_calc_feedback_runs.py)
regenerates the §4 sweep and the runs of [Simulation Evidence](SIMULATION_EVIDENCE.md) §2;
[operator_feedback_jump_runs.py](../simulations/operator_feedback_jump_runs.py)
checks the identities of §2 and §8 exactly, regenerates §8 against the tool's
logged values, and compares the page's other results with their exact or
regenerated values (the two quoted from Dynamic Entanglement excepted). Exact, here and below, means propagation by expm of the
Liouvillian or by DOP853 at rtol 10⁻¹² (10⁻¹¹ for the GHZ runs of §8.1).

---

## 4. The February 4 sweep

Parameters: Bell+, Heisenberg, operator law with κ = 0.5, t_max 10, γ₀ in
[0.003, 0.006], h in [0.7, 1.0]; the concurrence bridge and t_max 10 are
inferred, since no record of the sweep names them and the C_final column
regenerates under them.

| γ₀ | h | C_final | C_final × 0.27 |
|---------|---|---------|---------|
| 0.005 | 0.7 | 0.909 | 0.245 |
| 0.005 | 0.8 | 0.912 | 0.246 |
| 0.005 | 0.9 | 0.914 | 0.247 |
| 0.005 | 1.0 | 0.917 | 0.248 |
| 0.006 | 0.7 | 0.891 | 0.241 |
| 0.006 | 0.9 | 0.897 | 0.242 |

The last column was printed as C·Ψ and read as a confirmation of CΨ ≤ ¼. It
is not the state's CΨ: the tool's sweep routine multiplies C_final by a fixed
psi_approx = 0.27, so the column follows the concurrence and sits below ¼
exactly when C_final < 25/27 ≈ 0.926; weaker decoherence (γ₀ = 0.003, h = 0.7)
gives 0.255. The routine lists as admissible the rows with C_final × 0.27 ≤ ¼
(and a positive R∞), and every row of the table is one of them. Read with the density
matrix's own Ψ, concurrence × l₁/3, the same runs stay above ¼ at every step:
the lowest value over t ≤ 10 is 0.289 in the tool's Euler runs and 0.277 in
exact propagation. The small rise of C_final with h is the integrator's (§5):
at these rates the plain Euler step leaves the physical states (its purity at
t = 10 rises with h in the γ₀ = 0.005 rows, from 1.44 at h = 0.7 to 2.68 at
h = 1.0), and C_final is read off the clipped image of that step; exactly, the concurrence at t = 10 is 0.899 at
γ₀ = 0.005 and 0.879 at γ₀ = 0.006, for every h of the sweep. The historical copy of
[Dynamic Fixed Points](../docs/historical/DYNAMIC_FIXED_POINTS.md) (§8) records
five rows of the same sweep, among them γ₀ = 0.006 at h = 1.0 with C_final 0.901,
which this table lacks and which regenerates the same way.

The three active operator-law runs of [Simulation Evidence](SIMULATION_EVIDENCE.md)
§2 (Bell+ on the chain, GHZ₃ and W₃ on the ring, γ₀ = 0.005, h = 0.9, t_max 5,
the tool's mutual-purity reading) swing up to about ½ and end at 0.405, 0.262
and 0.413; exactly they end at 0.402, 0.248 and 0.421. GHZ₃ ends at the edge of
a crossing: it falls through ¼ at t = 4.998 exactly and, in the tool's routine
run past the logged window and read at every step, between t = 5.01 and 5.02 (exactly it rises through
¼ again at t = 5.48), so which side its last value lands on is decided by a
crossing both runs make within 0.02 of the window's edge; exactly, it spends
3.89 of its 5 time units above ¼ (largest 0.493). In the tool
the ⟨O_int⟩ of GHZ₃ drifts from 0 to 0.07 and that of W₃ from 0.67 to 0.70;
exactly the first stays at 0 and the second falls to 0.623 (the symmetric case
of §2), so both drifts are the clipping's (§5).

---

## 5. The Euler step

An explicit Euler step does not keep a state a state. Its Hamiltonian part,
all of the step at γ = 0, opens a negative eigenvalue of about −dt²·ΔH² on a
pure state at once and raises Tr ρ² by exactly dt²·‖[H, ρ]‖²_F (the
first-order term vanishes, since Tr(ρ[H, ρ]) = 0): the first step from Bell+
under the §8 Hamiltonian at γ = 0 opens −1.0·10⁻⁴ and gains 2.0·10⁻⁴ of purity,
and the purity climbs past 1, to 1.111 at t = 5 and 1.246 at t = 10. Dephasing
pulls both back only partly: in the §8 runs (γ₀ = 0.1) the purity of the plain
σ_z step stays at most 1 on a grid of κ up to 0.949 and passes it from κ = 0.95 on, where
the first step's gain from H meets its loss to dephasing (1.047 at κ = 0.99),
while its lowest
eigenvalue reaches −0.004 by t = 5 even at κ = 0, and in the decoherence-free
σ_x⊗σ_x channel the purity reaches 1.111 at every κ; in the §4 sweep (γ₀ = 0.005
and 0.006) and the three active runs of [Simulation Evidence](SIMULATION_EVIDENCE.md)
§2 (γ₀ = 0.005) it passes 1 at κ = 0, 0.5 and 1 (2.68 at t = 10 in the h = 1.0
row of §4, at κ = 0.5). Where H does not turn the state, as at h = 0 in the Bell+ and GHZ runs
here, the step is a mixture of the jumps and stays a state. In the tool's runs
on this page where H turns the state the plain step leaves the states, and the
tool's numbers there are the clipped image of that step. The instability belongs to the integrator, not to the feedback: an
adaptive solver (DOP853 at rtol 10⁻¹²) runs Bell+ at κ = 0.99 and κ = 1 within
the order of its tolerance of the physical states.

The same step damps too little where H turns the state faster than it decays,
and that, more than the clipping, is the gap between the tool's §8.2 purities
and C·Ψ values and the exact ones, up to 11 per cent: for the C·Ψ of the σ_z
run it is 11.6 per cent without the clipping and 10.6 with it, and it halves
with dt (5.3 per cent at dt = 0.005, 2.7 at 0.0025, with the clipping). In the
σ_x⊗σ_x row the clipping is what holds the purity at 1.

The tool takes the Hermitian part, sets negative eigenvalues to zero and
renormalizes after every step. That keeps ρ a density matrix and moves the
trajectory, at places where the exact one has no feature. The plain Euler step
leaves ⟨O_int⟩ independent of h, even where its matrix is no longer a state,
since on it the step is the scalar recursion
x ← x·(1 − 4γ_eff·dt), which H does not enter (0.9500415 at every h at
γ₀ = 0.005, κ = 0.5, t = 5, against the closed form's 0.9500416); with the
clipping it reads 0.951246 at h = 0.5 and 0.953819 at h = 0.9. The last is the
logged value of the Bell+ run of [Simulation Evidence](SIMULATION_EVIDENCE.md)
§2, so the field dependence of that run's ⟨O_int⟩ is the clipping's, as are the
drifts of the GHZ₃ and W₃ runs of §4 (the plain step gives 0 and 0.623 there)
and the lift of ⟨O_int⟩ off zero in §8.1. In §8.1 the clipping also pulls the purity toward
the exact value: GHZ₄ ends at 0.2885 without it, 0.2523 with it, 0.2503
exactly.

---

## 6. What decides whether there is feedback

The operator law sets the rate by an expectation value, and which observable
it reads decides whether there is any feedback at all. When O_int commutes
with H and with every jump (all Hermitian here), it lies in
[F158](../docs/ANALYTICAL_FORMULAS.md#f158)'s near kernel, ⟨O_int⟩ is conserved, and
the rate is a constant set by the initial state: the `x_pairs` runs of §8.1 on the
Heisenberg chain with the x field are of this kind for every initial state.
When the observable is an eigen-operator of the dynamics instead, as σ_x⊗σ_x is
under σ_z jumps on two qubits, or a multiple of one on the states in question,
as on permutation-symmetric states (§2), its expectation obeys an equation of
its own and the feedback has a closed form. An observable gives feedback only
where the dynamics move it.

In a model of this kind other observables could set the rate (in the tool the
observable changes only together with the jump); whether one gives feedback is
decided, as above, by how the Hamiltonian and the jumps act on it.

[Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md) §3 reads CΨ = ¼ as an observer
information bandwidth limit. That is an interpretation, and the runs here do
not test it.

---

## 7. The two laws side by side

| | Scalar law | Operator law |
|---|---|---|
| Rate | γ_base·C(ρ) | γ₀·(1 − κ·⟨O_int⟩) |
| Set by | a bridge function of ρ, chosen by `bridge_type` | the expectation of an observable fixed by the jump setting |
| Jumps | σ_z on every site, or one collective Σσ_z | σ_z on every site by default, or the jumps of §3 |
| Enters every step | yes | yes |
| On Bell+ (σ_z jumps) | a constant rate γ_base/2 with the default mutual-purity bridge under the Heisenberg bond and an x field, since both single-site states stay maximally mixed (§8.2; on the XY or Ising bond with the field on they do not); [Dynamic Fixed Points](DYNAMIC_FIXED_POINTS.md) §4 records it at h = 0, Ψ(10) = 0.302 = e^(−0.1)/3 | a closed form (§2) |

---

## 8. Observables and jump operators (2026-02-20)

These runs vary the jump operator, and with it the observable, at J = 1,
h = 0.5 (along x), γ₀ = 0.1, t_max 5, on the Heisenberg chain, read with the
mutual-purity bridge. The Bell+ calls behind the table are in the log of the
chat's MCP server (2026-02-20, 17:22 to 17:29 UTC), and the transcription
reproduces every ⟨O_int⟩ value the log shows of them; the logged σ_x⊗σ_x call ran at κ = 0.95, which changes
nothing on a decoherence-free state. No GHZ call at these settings is in the
log, so the GHZ settings are inferred from the February record of these runs
(commit 396f3f18): its three-digit purities, 0.252, 0.063 and 0.063 for GHZ₄, GHZ₅ and
GHZ₆, regenerate at h = 0.5, the field of the Bell+ calls of that afternoon (and
at h = 0.45, not at h = 0.2, where GHZ₄ gives 0.250), and its GHZ₄ C·Ψ at t = 5,
0.127, regenerates at h = 0.5 and not at h = 0.45 or 0.55 (0.103, 0.144);
the transcription also reproduces all twelve logged ⟨O_int⟩ values of the one
GHZ₄ operator-feedback call with `x_pairs` that evening (19:07 UTC, h = 0.2,
κ = 0.99, t_max 30).

### 8.1 A conserved observable

Under `x_pairs` every pair carries a jump X_iX_j and O_int is their mean. For
GHZ_N, κ does nothing: from κ = 0 to κ = 1 the runs record the same purity
(the tool's within 6·10⁻⁵: its clipping, §5, lifts ⟨O_int⟩ to 0.05 at N = 4,
where the plain step keeps it at 0). The reason is a conservation law. Write
S = Σ_{i<j} X_iX_j = ((Σ_k X_k)² − N)/2. The Heisenberg bonds commute with the
total spin, the field is the total x spin, and every jump X_iX_j commutes with
S, so S lies in [F158](../docs/ANALYTICAL_FORMULAS.md#f158)'s near kernel, ⟨S⟩ is
conserved for every initial state, and O_int, a multiple of S, with it. Under
this Hamiltonian the `x_pairs` feedback is never feedback: its rate is fixed by
the initial state.

For GHZ_N with N ≥ 3 that value is zero. A pair operator flips two bits and
has no matrix element between |0…0⟩ and |1…1⟩: the GHZ coherence carries an
excitation-number difference N, and a two-site observable sees only |ΔN| ≤ 2
([F70](../docs/ANALYTICAL_FORMULAS.md#f70-δn-selection-rule-for-site-local-observables-tier-1-proven-kinematic-lemma)'s
kinematic blindness, in its generalisation). So the rate stays γ₀, and the
purity is that of the plain channel. In the x basis GHZ_N is the even
superposition of the x strings with an even number of minus signs, the jumps
are diagonal, and a coherence between two strings that differ on a sites
decays at 2γ₀·a·(N − a), the number of pairs with one site in the difference
times 2γ₀. Hence

```
P(t) = 2^−(N−1) · Σ_{a even} C(N, a) · e^(−4γ₀·a·(N − a)·t),
```

independent of J, h and κ: the field only adds phases in the x basis, and the
bonds do not move the state, which keeps its symmetry under every permutation
of the sites (§2). Bell+ (N = 2) is the exception, with ⟨X₀X₁⟩ = 1; under
`x_pairs` its one jump is X₀X₁, the decoherence-free case of §8.2.

| System | ⟨O_int⟩ | Purity at t = 5, exact | The tool, κ from 0 to 1 |
|--------|---------|------------------------|-------------------------|
| GHZ₃ | 0 | 0.2637 | not run at these settings |
| GHZ₄ | 0 | 0.2503 | 0.2523 to 0.2524 |
| GHZ₅ | 0 | 0.06261 | 0.06263 |
| GHZ₆ | 0 | 0.06250 | 0.06333 |

GHZ_N starts below ¼ for N ≥ 3 (its Ψ(0) is 1/(2^N − 1), [F18](../docs/ANALYTICAL_FORMULAS.md#f18-fold-threshold-tier-2-n-product-state-measured-n--2-5),
[N-Scaling Barrier](N_SCALING_BARRIER.md); in the tool's reading, with C = ½,
CΨ(0) = 1/(2·(2^N − 1))), and the field lifts it: exactly, in that reading,
GHZ₃ rises to 0.339 and GHZ₄ just crosses ¼ (largest 0.25004), while GHZ₅ and
GHZ₆ stay below (0.164 and 0.115) through t = 5. GHZ₄'s crossing turns on the
inferred field being 0.5: at h = 0.45 its largest value is 0.238.

The conservation belongs to this Hamiltonian and this observable, not to the
GHZ state. With jumps on the nearest-neighbour pairs only, S is still
conserved, but O_int is now their mean, no multiple of S, and the bonds move it;
on the XY chain S itself is no longer conserved. At N = 4 (κ = 0) ⟨O_int⟩
reaches 0.05 and 0.30, and κ moves the purity (0.2755 to 0.2784, and 0.1723 to
0.1738, from κ = 0 to κ = 0.99). Another run with a pair observable at zero is
in [Mediator as Quantum Transistor](../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md)
§3.1: σ_z jumps with O_int = X₀X₁ on GHZ₃ at h = 0; there the reason is the
symmetric case of §2, ⟨X₀X₁⟩ starting at zero and staying there.

### 8.2 Jump operators on Bell+

The second set keeps σ_z, or replaces it by one two-qubit jump on qubits 0
and 1, at κ = 0 (the rate is γ₀ throughout). Two facts carry the table.

First, the bridge does not move. In all four channels both single-site states
of Bell+ stay maximally mixed: σ_x⊗σ_x fixes H and Bell+, maps each jump to
± itself and every σ_y and σ_z of a site to its negative; X₀ + X₁ commutes
with H, and each dissipator maps it to a multiple of itself; the swap of the
two sites does the rest. None of this uses the rate, so it holds under either
feedback law. So C = ½ and C·Ψ = l₁/6: every statement about ¼ in this reading
is a statement about l₁ = 3/2, and Bell+ starts at 1/6, below ¼.

Second, the field is along x, so H commutes with σ_x⊗σ_x and keeps Bell+ in
the sector σ_x⊗σ_x = +1, spanned by Φ⁺ and Ψ⁺, where it turns Φ⁺ into Ψ⁺ and
back; the two-qubit jumps keep that sector too, while the σ_z jumps move the
state out of it (Z₀Φ⁺ = Φ⁻). In the computational basis the turning raises l₁,
which is what lifts every channel above ¼ at first.

| Jump | C·Ψ at t = 5, tool / exact | Purity at t = 5, tool / exact | Largest C·Ψ in [0, 5], exact | Time above ¼ in [0, 5], exact | Last time above ¼, exact |
|------|-------|-------|-------|-------|-------|
| σ_z on each site | 0.124 / 0.112 | 0.336 / 0.323 | 0.413 at t = 0.695 | 1.66 | 2.54 |
| σ_x⊗σ_x | 0.348 / 0.348 | 1.000 / 1.000 | ½ | 4.12 | never stops |
| σ_y⊗σ_y | 0.287 / 0.275 | 0.734 / 0.693 | 0.476 at t = 0.761 | 3.81 | 13.50 |
| σ_z⊗σ_z | 0.287 / 0.275 | 0.734 / 0.693 | 0.476 at t = 0.761 | 3.81 | 13.50 |

- **σ_x⊗σ_x is decoherence-free.** The jump acts on each of its two eigenspaces
  as a sign, so its dissipator vanishes on every state inside one of them (on
  the sector of Bell+, and on the −1 sector of Φ⁻ and Ψ⁻; only coherences
  between the two decay), and since H keeps the sector, the purity stays 1 along
  the whole trajectory, at every κ, while C·Ψ = (1 + 2|sin 4ht|)/6 swings
  between 1/6 and ½ for good. Bell+ is an eigenstate of σ_y⊗σ_y and σ_z⊗σ_z as
  well; what protects it under σ_x⊗σ_x is that the field commutes with that jump
  and not with the others. With the field along z the roles turn: σ_z⊗σ_z
  protects (purity 1), and σ_x⊗σ_x and σ_y⊗σ_y decay alike (0.693). This
  σ_x⊗σ_x is the same operator as the observable of §2: as an observable under
  σ_z jumps it decays at the far kernel's price, as a jump it protects its
  eigenspaces, and both readings rest on its commuting with H. A sibling, Bell+
  under anti-correlated Z noise, is Test B of [the monotonicity proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md).
- **σ_y⊗σ_y and σ_z⊗σ_z are one channel here.** On the sector,
  σ_y⊗σ_y = −σ_z⊗σ_z, and the sign of a jump does not enter its dissipator, so
  without feedback the two give the same ρ(t); the tool's log shows their
  ⟨O_int⟩ tails at κ = 0 as exact negatives (0.580553 and −0.580553). With
  feedback they part, since O_int is then the jump itself and its sign sets the
  rate: the logged tails at κ = 0.95 end at 0.594354 and −0.565059, and exact
  propagation at κ = 0.5, a value no logged yy or zz call used, gives purities of 0.6918
  and 0.6929 at t = 5. Without
  feedback they relax to the even mixture of Φ⁺ and Ψ⁺, where C·Ψ = 1/6, and
  cross ¼ for the last time at t = 13.50.
- **σ_z** rises to 0.413 at t = 0.695 and falls below ¼ for good at t = 2.54.
  With feedback (κ = 0.95, logged the same afternoon) it crosses ¼ twelve
  times in exact propagation, dips below it between t = 4.52 and 4.92, is above it again at t = 5, at
  0.287 with purity 0.785 in the tool's regenerated run (0.284 and 0.773
  exactly), and falls
  below for good at t = 8.85: the feedback delays the last crossing, as in §2.

So in the tool's reading the number of qubits a jump acts on does not sort
these runs by ¼. All four rise above it, the field turning Φ⁺ into Ψ⁺ and so raising the
computational-basis coherence l₁ (that a fixed semigroup can carry CΨ upward
through ¼ is in the repository too, by another mechanism and in the purity
reading: [F17](../docs/ANALYTICAL_FORMULAS.md#f17), Part 6 of
[the monotonicity proof](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)), and the
three that decohere come back down; on which side a run sits at t = 5
says when it was read. The reading matters as much: with purity or concurrence
for C, the products the repository names, all four start at 1/3, above ¼, and
σ_x⊗σ_x never comes down to it, since the state stays pure and maximally
entangled and CΨ = (1 + 2|sin 4ht|)/3 stays between 1/3 and 1. CΨ = ¼ is the
discriminant zero of the fixed-point recurrence, a level a trajectory may
cross, recross or miss, not a bound a channel keeps or violates. The tool's
purities and C·Ψ values here stand up to 11 per cent from the exact ones (§5);
at t = 5 the σ_y⊗σ_y and σ_z⊗σ_z runs sit above ¼ in both.

### 8.3 What these results do and do not show

**They show:**
- Under σ_z jumps the operator feedback has a closed form (§2), on every state
  of two qubits under a Hamiltonian commuting with σ_x⊗σ_x, and on every state
  invariant under permutations of the sites under isotropic Heisenberg bonds and
  a uniform x field: it slows the decay while the correlation is high and so
  shifts it later (earlier where the correlation starts negative), the late
  rate unchanged, and κ = 1 holds the correlation of Bell+ only on an unstable
  fixed point.
- Whether there is feedback at all depends on whether the dynamics move
  O_int: under the Heisenberg chain with the x field and the `x_pairs` jumps
  the all-pairs observable is conserved for every state, and it is zero on GHZ_N
  for N ≥ 3.
- A Hermitian jump that commutes with H, as the only jump, leaves every state
  inside one of its eigenspaces decoherence-free; σ_x⊗σ_x does that for Bell+
  under the x field.
- In the tool's reading every Bell+ run of §8.2 crosses ¼, and the decohering
  ones cross back; whether a run is above ¼ depends on the time and on the C it
  is read with.

**They do not show:**
- That any of this needs R = CΨ² to be described: every result here is a
  conservation law, a closed form or a numerical reading of Lindblad dynamics
  with a state-dependent rate.
- What the nonlinear feedback law of §1 describes physically.
- More than the stated settings for the numbers. The identities behind them
  hold where stated: those of §8.1 at every N, J and h; the closed form of §2 on
  every state of two qubits under a Hamiltonian commuting with σ_x⊗σ_x, and on
  permutation-symmetric states at every N under isotropic Heisenberg bonds and a
  uniform x field; the sector identities of §8.2 at N = 2.
