# The star's frozen survivor is the structural ceiling, read dynamically

**What this is about.** Watch the longest-lived coherence of a dephased spin network as you weaken the
dephasing (raise `Q = J/γ`, the coupling-to-observation ratio). Does that survivor stay *frozen* (a pure
decay, no oscillation, `|Im λ| = 0`) or does it *un-freeze* (start ringing)? The answer separates the
topologies, and it gives the star its own place in a trichotomy the typed layer had named only for the
chain and the ring. (`|Im λ|` is the imaginary part of the survivor's Liouvillian eigenvalue: its
oscillation frequency. `g2 = 4/(N−1)` is the star's structural ceiling, the darkest reachable decay-gap.)

## The finding

A single `Q` does not tell the topologies apart: *below* its coherence horizon every topology's slowest
mode is overdamped (real, `|Im| = 0`). At `Q = 1.5` the chain `(2,2)`, ring `(2,2)`, and star `(1,1)`
survivors are all frozen. The signature is the `|Im|(Q)` *curve* (gate-verified, `simulations/star_frozen_seam.py`):

| topology | survivor as `Q` grows | un-freezes? |
|---|---|---|
| **chain** | `(p,p)` interior, frozen, then the `(0,1)` band edge takes over above `Q*(N)` | **yes** (the coherence horizon) |
| **ring** | `(2,2)` frozen seam, then the oscillating band edge above the handover | **yes** (the handover) |
| **star** | `(1,1)` boundary survivor, frozen at **every** `Q` (for `N ≥ 5`) | **no** |

So the star has its own frozen survivor, and it **never un-freezes**: the third case of the trichotomy.

## Why, and the threshold

The star's survivor is the darkest `[H,A] = 0` **commutant** coherence (the structural-ceiling mode).
A coherence that commutes with `H` has no coherent evolution, `−i[H,ρ] = 0`, so it cannot oscillate: it
is frozen by construction. But it is only the *survivor* (the slowest mode) when it is darker than the
`−2γ` Absorption floor, i.e. exactly when the ceiling `g2 = 4/(N−1) ≤ 1`:

- `N ≥ 5` (`g2 ≤ 1`): the commutant coherence undercuts the floor, is the survivor, and is frozen at all
  `Q`. Gate-verified frozen at `N = 5, 6, 7, 8` (max `|Im|` over the `Q` sweep `< 1e-15`).
- `N = 4` (`g2 = 4/3 > 1`): the commutant mode is not dark enough; an oscillating band-edge mode wins and
  the star **un-freezes** (max `|Im| = 1.73`). This is the known star outlier (the `(2,2)`/`K₄` half-filling
  special case).

**So the star's frozen seam IS the structural ceiling `g2 = 4/(N−1) ≤ 1`, read dynamically.** The high-Q
ceiling and the all-Q frozenness of the survivor are the same fact: the commutant coherence sits below the
floor. This is the survivor-level reading of `PROOF_STRUCTURAL_CEILING.md` §7's "the star has no coherence
horizon" — the global slowest mode never acquires a frequency.

## What this closes

The typed layer characterized *frozen vs oscillating* for the chain (oscillates, the SE-EP horizon) and
the ring (frozen `(2,2)` level crossing), but the star was absent from those statements (only "flat band /
no horizon / boundary" was typed). This names the star's case and ties it to an already-typed quantity
(the ceiling), so it is a consolidation, not a new mechanism. The star's freeze is a *different route* to
`|Im| = 0` than the ring's: the ring freezes by a **level crossing** (two real eigenvalues coincide), the
star by a **commutant** (the mode commutes with `H`). Both are frozen; the chain is neither.

## See also

- Verifier (self-validating): `simulations/star_frozen_seam.py` (the pin at `Q=1.5`; the `|Im|(Q)` sweep;
  the `g2 = 4/(N−1) ≤ 1` threshold across `N = 4..8`). Reuses `simulations/carbon/incompleteness_survivor.py`
  (the sector-projected survivor, validated bit-for-bit vs the full `4^N` `L` at `N = 4`).
- The ceiling: `docs/proofs/PROOF_STRUCTURAL_CEILING.md` (`g2 = 4/(N−1)`, §7 no-horizon); `StructuralCeilingClaim`,
  `SecondClockRegimeClaim` (the regimes), `CoherenceHorizonClaim` (the chain's un-freezing).
- The seam picture this extends: `docs/STERILE_BIRTHCANAL_AND_THE_JUNCTION.md` (chain/ring frozen vs
  oscillating); the sibling negative result `experiments/THE_HUB_KILLS_THE_HORIZON.md` (a hub removes the
  horizon however dispersive the leaves).
