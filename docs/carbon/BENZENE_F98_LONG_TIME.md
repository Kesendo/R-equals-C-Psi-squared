# Selected C₄/C₆ XX+YY Rings and the F98 Long-Time Bridge

**Date:** 2026-05-22
**Authors:** Tom + Claude
**Status:** F98 is a Tier-1 derived theorem for its specified KIntermediate
initial state, a W-conserving Hamiltonian, and all-site Z-dephasing. The C₄/C₆
calculations are exact selected-model instances; they do not assign a material
carbon degree of freedom, β-to-J convention, bath, γ, T₂, or Q (see
[README](README.md#conditional-c4-and-c6-working-model) and
[Q audit](../Q_BELONGS_TO_NO_SUBSTANCE.md)).
**Script:** [`simulations/carbon/benzene_f98_long_time.py`](../../simulations/carbon/benzene_f98_long_time.py)
**Tested:** selected C₄ and C₆ XX+YY spin rings.
**Answers:** [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md) open
question 2 / [README](README.md) open question 5.

---

## The question

F98 is the framework's long-time bridge. The KIntermediate Dicke superposition
ψ = (|D_{N/2−1}⟩ + |D_{N/2}⟩)/√2, evolved under a magnetization-conserving
Hamiltonian (`[H, Ŵ] = 0`, with `Ŵ = Σ_l (I − Z_l)/2`) plus Z-dephasing on every
site, ends at a Π²-odd Frobenius² fraction

```
α(∞) = (N + 2) / [4(N + 1)]      →  1/4  as N → ∞
```

It is the dynamical partner of the static F86b 3/8 anchor: the same KIntermediate
state begins at α(0) = 3/8 and decays along an N-dependent curve to
(N+2)/[4(N+1)]. F98 is valid for any bond graph under its stated premise; the
source records bit-exact N = 4..16 checks on the Heisenberg chain
([`proton_chain_dicke_anchor.py`](../../simulations/water/proton_chain_dicke_anchor.py)).

The question here is narrower: do selected C₄/C₆ XX+YY rings, prepared in that
specified KIntermediate state and given all-site Z-dephasing, instantiate the
F98 premise?

## The selected C₄/C₆ model instances

Choose `n_l = (I − Z_l)/2` as a site-occupation coordinate and XX+YY as the
free-hopping spin model. XX+YY commutes with `Ŵ`, so these selected ring
Hamiltonians meet F98's Hamiltonian condition. The Hückel `β → J` relation is a
Tier-2 structural translation within that model; it does not select a
material-carbon Hamiltonian or convention.

If the selected local jump is the density `n_l`, then
`D[n_l] = ¼·D[Z_l]`. The script evaluates the corresponding all-site
Z-dephasing model. This conditional algebra does not identify a physical carbon
degree of freedom or a material bath, rate, γ, T₂, or Q.

The KIntermediate initial state is load-bearing. F98 does not state an asymptote
for arbitrary initial states, and W conservation is likewise a premise rather
than a consequence of arbitrary XX/YY terms.

## Result

The specified KIntermediate state was evolved under the selected C₄/C₆ XX+YY
plus all-site-Z Liouvillian. α(∞) is the exact t → ∞ limit obtained by projecting
ρ onto `ker L`.

| Ring | α(t = 0) | α(∞) measured | F98 (N+2)/[4(N+1)] | Match |
|------|----------|---------------|--------------------|-------|
| C₄ selected ring | 3/8 | 0.30000000 | 3/10 | numerical agreement (\|diff\| = 4.44 × 10⁻¹⁶) |
| C₆ selected ring | 3/8 | 0.28571429 | 2/7 | numerical agreement (\|diff\| = 2.05 × 10⁻¹⁵) |

The direct selected-model runs use `J = 1` and `γ = 0.5`. Their float-eigensolver
outputs agree with F98 at the stated numerical precision; F98's closed form remains
the exact theorem. Two further selected-model checks pass: `ker L` has dimension N + 1 (the F4 prediction,
`ker L = span(P_0, …, P_N)`); and the long-time state ρ_∞ equals the
F98-predicted `½·[P_m/C(N,m) + P_{m+1}/C(N,m+1)]`, `m = N/2 − 1`, to 10⁻¹⁵.

The selected C₆ trace reads 0.375 at t = 0, 0.327 at t = 1, and settles onto
0.2857 by t ≈ 3 at `γ = 0.5`: the F86b static anchor 3/8 relaxing onto F98's
long-time value 2/7.

## The answer to Question 5

**Yes.** With F98's specified KIntermediate state, the selected C₄/C₆ XX+YY
rings and all-site-Z channel yield the F98 bridge: 3/10 for C₄ and 2/7 for C₆.
This is an explicit selected-model instance of the theorem, not a claim about
arbitrary benzene dynamics, a material bath, or a carbon material mapping.

## What the selected instances check

- **Ring topology within F98's premise.** F98 is valid for any bond graph with a
  W-conserving Hamiltonian and all-site Z-dephasing; the C₄/C₆ rings are direct
  selected instances of those conditions.
- **XX+YY within the W-conserving class.** The selected ring Hamiltonian commutes
  with W. This condition, rather than a Hückel or carbon label, supplies the
  F98 connection; a Hamiltonian that does not conserve W is outside F98's scope.
- **The specified state and long-time operator.** The selected propagation and
  `ker L` projection recover both the closed-form α(∞) and its stated ρ_∞ for
  the KIntermediate initial state.

## Framework-vocabulary translation

| Selected model item | Framework object | Status |
|---------------------|------------------|--------|
| free-hopping Hückel/framework comparison | XX+YY spin ring, N=6 | Tier 2 structural translation |
| selected site occupation and density jump | local Z-dephasing, `D[n_l] = ¼·D[Z_l]` | exact after `n_l` and the jump are selected |
| KIntermediate Dicke (\|D₂⟩+\|D₃⟩)/√2 | F98 / F86b KIntermediate anchor | Tier 1 specified initial state |
| C₆ long-time Π²-odd fraction 2/7 | F98 α(∞) = (N+2)/[4(N+1)] at N=6 | Tier 1 selected-model result |
| static Π²-odd fraction 3/8 | F86b Dicke anchor α(0) | Tier 1 |

## Open follow-ups

- F98 fixes the asymptote only for its KIntermediate state. Other selected
  initial states require their own state-class analysis.
- A selected mixed `D[Z] + D[B]` jump model is outside F98's all-site-Z premise.
  Its result would be a model-channel question, not a material-bath
  classification.
- [The Frost circle clock comparison](FROST_CIRCLE_AS_THE_CLOCK_FACE.md) supplies
  Q* values only for its selected XX+YY/Z-dephasing model; it supplies no carbon
  `γ`, `T₂`, or Q.

## Anchor

- Script: [`simulations/carbon/benzene_f98_long_time.py`](../../simulations/carbon/benzene_f98_long_time.py)
- Companion doc: [Selected C₄/C₆ ring Liouvillians](BENZENE_LIOUVILLIAN_PALINDROME.md)
  (the F1 selected-channel comparison), [README.md](README.md)
- Framework anchors: F98 long-time bridge, F86b Dicke anchor, F4 kernel
  decomposition, all in [`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md);
  [`simulations/water/proton_chain_dicke_anchor.py`](../../simulations/water/proton_chain_dicke_anchor.py)
  (the Heisenberg-chain verification)
