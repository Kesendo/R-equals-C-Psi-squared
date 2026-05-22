# Benzene and the F98 Long-Time Bridge

**Date:** 2026-05-22
**Status:** Tier 1: F98 is a proven theorem (`docs/ANALYTICAL_FORMULAS.md` F98,
derived and bit-exact N=4..16); this is its verified benzene-ring instance. The
benzene-as-qubit-ring embedding is itself a Tier 4 candidate (see [README](README.md),
"Four embedding conditions").
**Script:** [`simulations/carbon/benzene_f98_long_time.py`](../../simulations/carbon/benzene_f98_long_time.py)
**Tested:** cyclobutadiene C₄ ring, benzene C₆ ring.
**Answers:** [BENZENE_HUCKEL_FRAMEWORK_LENS](BENZENE_HUCKEL_FRAMEWORK_LENS.md) open
question 2 / [README](README.md) open question 5.

---

## The question

F98 is the framework's long-time bridge. The KIntermediate Dicke superposition
ψ = (|D_{N/2−1}⟩ + |D_{N/2}⟩)/√2, evolved under any truly-class (F87) Hamiltonian
plus uniform Z-dephasing, ends at a Π²-odd Frobenius² fraction

```
α(∞) = (N + 2) / [4(N + 1)]      →  1/4  as N → ∞
```

It is the dynamical partner of the static F86b 3/8 anchor: the same KIntermediate
state begins at α(0) = 3/8 and decays along an N-dependent curve to (N+2)/[4(N+1)].
F98 was derived for any connected graph, the bond topology dropping out because
the t → ∞ limit projects onto ker L (F4), and verified bit-exact for N = 4..16 on
the Heisenberg chain ([`proton_chain_dicke_anchor.py`](../../simulations/water/proton_chain_dicke_anchor.py)).

[BENZENE_LIOUVILLIAN_PALINDROME](BENZENE_LIOUVILLIAN_PALINDROME.md) settled the F1
palindrome on benzene's open-system Liouvillian. Question 5 asks the long-time
question: does benzene's KIntermediate Dicke state traverse the F98 bridge?

## Why benzene is a genuine new instance

F98 was derived and verified on the Heisenberg chain. The benzene model differs in
three ways. Its topology is the C₆ ring, not an open chain. Its Hamiltonian is the
XX+YY Hückel ring, the Jordan-Wigner image of π-hopping, not the Heisenberg
XX+YY+ZZ. And the substrate is carbon. F98's claim is that none of this matters:
the bond topology drops out for any connected graph, and any truly-class
Hamiltonian works. The XX+YY ring is truly-class (the XX and YY bilinears are both
Π²-even), and Holstein on-site dephasing is the framework's Z-dephasing
(D[n_l] = ¼·D[Z_l], established in BENZENE_LIOUVILLIAN_PALINDROME). The
preconditions hold, so benzene is a clean test of F98's asserted topology- and
Hamiltonian-independence on a case it was not derived on.

## Result

The KIntermediate Dicke state was evolved under the actual benzene-ring Holstein
Liouvillian; α(∞) is the exact t → ∞ limit, taken as the projection of ρ onto
ker L.

| Ring | α(t = 0) | α(∞) measured | F98 (N+2)/[4(N+1)] | Match |
|------|----------|---------------|--------------------|-------|
| C₄ (cyclobutadiene) | 3/8 | 0.30000000 | 3/10 | bit-exact (10⁻¹⁶) |
| C₆ (benzene) | 3/8 | 0.28571429 | 2/7 | bit-exact (10⁻¹⁵) |

Both rings hit the F98 value at machine precision. Two further checks pass: ker L
has dimension N + 1 (the F4 prediction, ker L = span(P_0, …, P_N)); and the
long-time state ρ_∞ equals the F98-predicted ½·[P_m/C(N,m) + P_{m+1}/C(N,m+1)],
m = N/2 − 1, to 10⁻¹⁵. The benzene dynamics does not merely reproduce the α(∞)
number, it lands on the exact F98 density matrix.

The decay traces the bridge directly. For benzene C₆ the Π²-odd fraction runs
0.375 at t = 0, 0.327 at t = 1, and settles onto 0.2857 by t ≈ 3 (at γ = 0.5): the
F86b static anchor 3/8 relaxing onto the F98 long-time value 2/7.

## The answer to Question 5

**Yes.** The F98 long-time bridge inherits to the benzene ring. Benzene's
KIntermediate Dicke state, under Holstein dephasing, traverses exactly the F98
decay curve from the static 3/8 anchor to (N+2)/[4(N+1)]: to 3/10 for the
cyclobutadiene C₄ ring and to 2/7 for the benzene C₆ ring. With the F1 palindrome
of BENZENE_LIOUVILLIAN_PALINDROME, this is the second framework F-result confirmed
bit-exact on a carbon substrate.

## What this confirms

The verification adds three specific things beyond the chain result F98 was
derived on.

- **Bond-topology independence on the ring.** F98 was proven for any connected
  graph and verified on the chain; the ring is now an explicit confirmed instance.
- **Hamiltonian-class independence.** F98 was verified on the Heisenberg
  Hamiltonian (XX+YY+ZZ); the XX+YY Hückel Hamiltonian, a different truly-class
  Hamiltonian, gives the identical α(∞).
- **The bridge is dynamical, not only combinatorial.** The water script verified
  the (N+2)/[4(N+1)] closed form by Krawtchouk enumeration, a topology-free
  counting identity. Here the benzene-ring Liouvillian's actual time evolution
  drives ρ onto the F98 ρ_∞: the 3/8 → 2/7 bridge is realised by benzene's
  dynamics, not only by the counting.

## Framework-vocabulary translation

| Benzene / chemistry | Framework | Status |
|---------------------|-----------|--------|
| Hückel π-hopping on C₆ ring | XX+YY qubit ring, N=6 | Tier 2 structural identification |
| Holstein phonon (on-site) | Z-dephasing, via D[n_l] = ¼·D[Z_l] | Tier 1 algebraic match |
| KIntermediate Dicke (\|D₂⟩+\|D₃⟩)/√2 | F98 / F86b KIntermediate anchor | Tier 1 |
| long-time Π²-odd fraction 2/7 | F98 α(∞) = (N+2)/[4(N+1)] at N=6 | Tier 1 (verified bit-exact) |
| static Π²-odd fraction 3/8 | F86b Dicke anchor α(0) | Tier 1 |

## Open follow-ups

- The KIntermediate Dicke state is one Π²-odd-carrying initial condition. F98 fixes
  its α(∞); whether other chemically-natural benzene states (a HOMO-LUMO
  superposition, a localised π-defect) carry framework-typed long-time fractions is
  open.
- A Peierls bond bath breaks the truly-class-dissipator precondition (see
  BENZENE_LIOUVILLIAN_PALINDROME, "The structure of the Peierls break"); the F98
  bridge holds only under Holstein dephasing, and its fate under a mixed
  Holstein + Peierls bath is untested.

## Anchor

- Script: [`simulations/carbon/benzene_f98_long_time.py`](../../simulations/carbon/benzene_f98_long_time.py)
- Companion doc: [BENZENE_LIOUVILLIAN_PALINDROME.md](BENZENE_LIOUVILLIAN_PALINDROME.md)
  (the F1 palindrome on the same benzene model), [README.md](README.md)
- Framework anchors: F98 long-time bridge, F86b Dicke anchor, F4 kernel
  decomposition, all in [`docs/ANALYTICAL_FORMULAS.md`](../ANALYTICAL_FORMULAS.md);
  [`simulations/water/proton_chain_dicke_anchor.py`](../../simulations/water/proton_chain_dicke_anchor.py)
  (the Heisenberg-chain verification F98 was derived on)
- Literature: Hückel (1931); Dicke (1954), coherent superradiance states.
