# The band-edge transition matrix: a Dirichlet-edge coupling fixes both the Frobenius deficit and the spectral floor

**Status:** Tier 1 derived; gate-exact N = 3..20. Open chain, band-edge carrier. 2026-06-20.
**Verifier:** `simulations/_handshake_M_checksum.py`, `simulations/_handshake_M_topology.py`, `simulations/_handshake_F124_adversarial.py`.
**Registry:** [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md) F124.
**Typed:** `BandEdgeTransitionInvariantClaim` (Tier1Derived, dual parents `KPartnerSelectionRuleClaim` + `ClockHandLadderClaim`) + live witness `inspect --root transition` (`compute/RCPsiSquared.Diagnostics/Foundation/BandEdgeTransitionInvariantWitness.cs`), which recomputes the identity, the genuine-minimum carrier selection, the frame identities, and the object/topology breakages across N. 2026-06-20.
**Framing corrected** after a five-lens panel (2026-06-20): `physics-first-review` caught the object-misnaming and the oversold "conservation law"; `mathematical-review` confirmed the proof (no error); `grounding-in-the-quantum` + `borrowing-a-discipline` supplied the frame-theory reading. This file replaces the earlier "location-matrix invariant" draft.

> **Object guard.** `M` here is the FULL single-excitation bond-transition matrix `M[b,k] = ⟨ψ_k|V_b|ψ_1⟩` over **all** N modes `k = 1..N` (including `k = 1`, the carrier / "strength" channel). It is **not** the handshake decoder's *location dictionary* (`hypotheses/HANDSHAKE_GEOMETRY.md`, which runs `k = 2..N`, excluding the strength channel); on that `k = 2..N` object the identity FAILS — the sum is not 2 and `λ_min = 0` (the K-partner ψ_N is a null column). The clean "2" **requires** the strength column. `M` is also distinct from the F1/F49 palindrome-residual `M`.

## What this is about

A band-edge standing wave `ψ_1` (the smoothest single-excitation mode, the F67 receiver) sits on an open chain. Perturb one bond `b`: `V_b` drives a unit of current across it, scattering `ψ_1` into the other modes, and `M[b,k] = ⟨ψ_k|V_b|ψ_1⟩` is that scattering amplitude. Two numbers come from the full scattering matrix — its Frobenius energy and the smallest eigenvalue of its Gram `M Mᵀ` — and they sum to exactly 2, the chain's coordination number. One boundary fact is the reason: the part of the carrier the bulk does not double-count is its weight on the two free ends, and that *same* end-weight is what the most contrarian (staggered, zone-boundary) bond-modulation couples to, the entire bulk cancelling.

## The real content (lead with this, not the "2")

The non-trivial fact is **`λ_min = E`**. For the BAND-EDGE carrier, the staggered (zone-boundary, `q = π`) bond modulation `Σ_b (−1)^b V_b` scatters `ψ_1` with a strength `λ_min = E` set **entirely by the Dirichlet ends**: a conserved standing-wave envelope `Q` of the carrier's recurrence makes the bulk telescope away *in the eigenvalue*, leaving only the boundary value `2Q_0 = E` (the scattered state itself is a bulk staggered-cosine — it is the coupling *strength*, not the wavefunction, that the boundary fixes). This **reads as** an **SSH / Peierls** edge effect: a dimerization's grip on a band-edge state is boundary-dominated. And the *same* `E` is the deficit of the carrier's degree-weighted norm from the coordination number. One boundary quantity, `c₀²`, governs both the spectral floor and the trace deficit.

## Statement

Open chain of `N` sites; carrier `ψ_1` = the band-edge mode, `c_i := ψ_1(i) = √(2/(N+1))·sin(π(i+1)/(N+1))`; bonds `b = (a, a+1)`, `V_b = |a⟩⟨a+1| + |a+1⟩⟨a|`; full transition matrix `M[b,k] = ⟨ψ_k|V_b|ψ_1⟩` (`N−1` bonds × `N` modes). For all `N ≥ 3`:
```
   ‖M‖_F² = z − E ,    λ_min(M Mᵀ) = E ,    ‖M‖_F² + λ_min = z = 2 ,
```
where `z = 2` is the chain's coordination number and `E = c₀² + c_{N-1}² = (4/(N+1))·sin²(π/(N+1))` is the carrier's weight on the two free ends. The `λ_min` eigenvector is the staggered bond wave `(−1)^b`.

## Proof

**Part 1 — `‖M‖_F² = z − E` (degree counting; basis-independent).** By completeness of the `N` modes, `‖M‖_F² = Σ_b Σ_k ⟨ψ_k|V_bψ_1⟩² = Σ_b ‖V_bψ_1‖²` — **independent of the mode basis** (the whole mode/transition column structure washes out). With `‖V_bψ_1‖² = c_a² + c_{a+1}²`, summing over bonds weights each site by its degree:
```
   ‖M‖_F² = Σ_site deg(site)·c_site² = ψ_1ᵀ D ψ_1 = 2·Σc² − c₀² − c_{N-1}² = 2 − E ,
```
using `Σc² = 1` and `deg = 2` interior, `1` at the two free ends. The "2" is the bulk coordination number `z` — it equals the bulk degree only because `ψ_1` is normalized, and it rides on `‖V_b‖² = 2` (scaling as `2α²` under `V_b → αV_b`). This half carries **no** transition content; it is the carrier's degree-weighted norm. ∎(Part 1)

**Part 2 — `λ_min(M Mᵀ) = E` (the band-edge edge-coupling).** `M Mᵀ` is the Gram of the bond-scattered carriers `{V_bψ_1}`, tridiagonal with `diag_a = c_a² + c_{a+1}²`, `offdiag(a,a+1) = c_a c_{a+2}`, and — for the band-edge (**nodeless**) carrier — **strictly positive** off-diagonals (`c_a c_{a+2} > 0`). The staggered vector `w_a = (−1)^a` gives
```
   (M Mᵀ w)_a = (−1)^a [(c_a² + c_{a+1}²) − c_{a-1}c_{a+1} − c_a c_{a+2}] = (−1)^a · 2Q_a ,
```
with `Q_a := c_a² + c_{a+1}² − E₁ c_a c_{a+1}`, `E₁ = 2cos(π/(N+1))` the band edge. `Q_a` is the **conserved discrete-energy envelope** of the carrier's recurrence `c_{a-1} + c_{a+1} = E₁ c_a`: `Q_a − Q_{a-1} = (c_{a+1}−c_{a-1})[(c_{a+1}+c_{a-1}) − E₁ c_a] = 0`. The bulk telescopes away; only the boundary survives. At `a = 0` (Dirichlet `c_{-1} = 0 ⟹ c_1 = E₁ c_0`), `Q_0 = c₀²`, and the boundary row gives `(c₀² + c_1²) − c_0 c_2 = 2c₀²` directly. So the staggered mode is an eigenvector with eigenvalue `2c₀² = c₀² + c_{N-1}² = E`. Because the off-diagonals are **strictly positive** (the band-edge condition), `M Mᵀ` is an *irreducible* Jacobi matrix: Perron-Frobenius gives the **largest** eigenvalue the all-positive eigenvector, hence the eigenvector of **maximal** sign changes carries the **least** eigenvalue. The staggered vector has the maximal `N−2` sign changes, so `E = λ_min`. ∎(Part 2)

Adding: `‖M‖_F² + λ_min = (2 − E) + E = 2`. ∎

## What singles out the band edge (the carrier-selecting step)

Part 2's conserved envelope `Q` holds for **any** eigenmode-carrier `ψ_k`: the staggered vector is always an eigenvector, with eigenvalue `2Q_k`. But only for the band-edge (nodeless) carrier are the Gram off-diagonals `c_a c_{a+2}` positive everywhere, which is what makes the staggered (maximal sign-change) mode the **minimum**. For an interior carrier the off-diagonals change sign, the staggered is no longer the least eigenvalue, and the sum drops strictly below 2 (verified: `N=7`, `ψ_2 → 1.75`, `ψ_3 → 1.58`). **The band-edge / F67-receiver carrier is the unique one for which "= 2" holds**; the conserved envelope is necessary but not the selector — the **positivity** is.

## The frame reading (grounding + borrowing-a-discipline, converged independently)

The bond-scattered carriers `{V_bψ_1}` form a **frame** for their span; `M Mᵀ` is the frame Gram, `MᵀM = Σ_b |V_bψ_1⟩⟨V_bψ_1|` the frame operator. They are a **deficient, non-tight Riesz basis**: rank `N−1`, with one exact kernel direction = the K-partner `ψ_N = (−1)^i ψ_1` (`⟨ψ_N|V_b|ψ_1⟩ ≡ 0`, the typed `KPartnerSelectionRuleClaim`). Then, term-for-term:
- `λ_min` = the optimal **lower frame bound** `A = σ_min²(M)` = the **Eckart-Young squared distance** from `M` to the nearest rank-deficient matrix = the **worst-case reconstruction floor** `1/‖S⁻¹‖`. The end-leakage `E` is not a loss; it is the **conditioner** that keeps the bond→mode map invertible.
- F124 is then a **completeness ⊕ conditioning** identity: the Parseval deficit `(z − trace)` (a trace quantity) equals the lower bound `A` (a single extremal eigenvalue) — both `= E`.
- The **condition number** `λ_max/E` grows with `N` (`3.4, 5.3, 7.6, …` for `N = 4,5,6`): shorter chains are better-conditioned; the long-chain / ring limit (`E → 0`) is ill-conditioned, the staggered K-direction going singular.

**K-partner caveat (two spaces, not one).** The kernel `ψ_N = (−1)^i ψ_1` is a **site-space** null mode of `MᵀM`; the staggered `λ_min` eigenvector `(−1)^b` is a **bond-space** vector. They share the bipartite `Z₂` grading but are **distinct** objects. On the location dictionary `k = 2..N` the strength column is absent and `λ_min` *is* the K-partner at 0; the strength column is exactly what lifts the floor to `E`.

## Scope and breakages

- **Carrier:** band-edge (nodeless) only. Other carriers → sum `< 2`.
- **Topology:** the "2" is the coordination number `z` (chain/ring degree-2 → 2; a degree-3 backbone → larger, weighted by `ψ_1²`; star hub → `‖M‖_F² = N/2`). The `λ_min = E` half needs the staggering to close: the chain (Dirichlet) and the **even** ring hold — but the even ring is *degenerate* (`E = 0`, `λ_min = 0`: no boundary, so it does **not** exhibit the deficit↔floor cancellation, not a second instance of the mechanism); the **odd** ring frustrates the staggering (`λ_min > 0`, sum `> 2`); the **star** breaks the trace half (the hub).
- **Object:** the full transition matrix `k = 1..N`; the decoder's location dictionary `k = 2..N` does *not* satisfy it.

## Honest calibration

The "2" is the coordination number — a *contract*, not a conservation law, and it scales with `‖V_b‖²`. The finding is the band-edge edge-coupling `λ_min = E` and its coincidence with the degree-weighted-norm deficit: **one boundary quantity `c₀²` fixing both the spectral floor and the trace deficit** — an SSH/Peierls edge effect with a frame-theoretic face (the lower frame bound = the distance to singularity = the chain's margin from rank-collapse).
