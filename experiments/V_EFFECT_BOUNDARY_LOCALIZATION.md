# V-Effect Boundary Localization: Where the Palindrome Can Break

**Status:** Computational result (Tier 1-2), N=3 and N=4 to machine precision (residual 10⁻¹⁵); the extreme-block immunity is proven in [the zero-immunity proof](../docs/proofs/PROOF_ZERO_IMMUNITY.md). Analytical interpretation Tier 2.
**Date:** 2026-04-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Pipeline:** `simulations/veffect_weight_immunity.py`
**See also:** [the V-Effect palindrome](V_EFFECT_PALINDROME.md), [the zero-immunity proof](../docs/proofs/PROOF_ZERO_IMMUNITY.md), [the bit-b parity symmetry proof](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md), [the primordial qubit](../hypotheses/PRIMORDIAL_QUBIT.md)

---

## What this finding establishes

For an N-qubit chain with two-body bond Hamiltonians H and Z-dephasing rate γ_l per site, the palindromic relation

```
Π · L · Π⁻¹ + L + 2Σγ · I = 0                     (*)
```

can be read block by block in XY-weight (w = total number of X or Y operators across all sites in a Pauli string): the residual M = Π·L·Π⁻¹ + L + 2Σγ·I splits into blocks (w_row, w_col). Three structural facts emerge:

1. **Extreme-block immunity (Π-fixed).** The diagonal blocks (w=0, w=0) and (w=N, w=N) of M vanish **for every two-body Hamiltonian H**, regardless of whether H respects or violates either Z₂ parity. This is the theorem of [the zero-immunity proof](../docs/proofs/PROOF_ZERO_IMMUNITY.md): the Hamiltonian has no part inside either extreme sector, the dissipator is zero on w=0 and the uniform −2Σγ on w=N, and Π carries that −2Σγ back as the shift (*) needs. "Two-body" means two-site bilinears; a single-site Z field breaks both blocks.

2. **Boundary diagonal blocks (H-fixed).** The diagonal blocks with 0 < w < N vanish when every term of H is truly (an even number of Y and an even number of Z, F85/F87). They break when H holds a non-truly term with two X/Y letters (XY or YX, the case of Tests 3 and 5). A non-truly term with a single X/Y letter (YZ, ZY, XZ, ZX) also breaks (*), but never on the diagonal: a letter X or Y flips its site's weight bit, so such a term moves w by exactly one, and its breaking sits entirely in the blocks (w, w ± 1). bit_b parity alone is therefore not the criterion: YZ is bit_b-even and breaks, off the diagonal.

3. **Off the diagonal the break reaches the extremes.** The blocks joining w=0 or w=N to other weights are not protected. In Tests 3 and 5 the blocks (0,2), (2,0), (1,3) and (3,1) carry norm 4√2 at N=3; at N=4 the blocks (0,2) and (2,4) and their transposes carry 8.

So the operator-level break behind the V-Effect "14 of 36 mode-pairs break, 22 remain palindromic" reading has two parts: on the diagonal, **all** of it lives in 0 < w < N, and the diagonal blocks of the extreme sectors w=0 (all I, Z) and w=N (all X, Y) are untouched by any bilinear; off the diagonal the boundary sectors and the extremes are coupled, and the coupling carries the break too, as [the V-Effect palindrome](V_EFFECT_PALINDROME.md) reads in its cross blocks w = 0↔2 and 1↔3. Which modes the spectrum then orphans is a separate, spectral question.

## Numerical results

`simulations/veffect_weight_immunity.py`, N ∈ {3, 4}, γ_l = 0.1 per site, Σγ = N · 0.1.

Five test Hamiltonians on each chain (bonds (0,1), (1,2), and at N=4 also (2,3)):

| Test | Hamiltonian per bond | Every term truly? |
|------|---------------------|-------------------|
| 1 | XX + YY + ZZ (Heisenberg) | Yes |
| 2 | XX + YY (XY model) | Yes |
| 3 | XX on (0,1), XY on (1,2) | **No** (XY has one Y) |
| 4 | XX + ZZ | Yes |
| 5 | ZZ on (0,1), XY on (1,2) | **No** |

**N=3 results (diagonal-block residual norms):**

| Sector | Block size | Test 1 | Test 2 | Test 3 | Test 4 | Test 5 |
|--------|------------|--------|--------|--------|--------|--------|
| w=0 | 8×8 | 0 | 0 | **0** | 0 | **0** |
| w=1 | 24×24 | ~10⁻¹⁶ | ~10⁻¹⁶ | **8.0** | ~10⁻¹⁶ | **8.0** |
| w=2 | 24×24 | ~10⁻¹⁶ | ~10⁻¹⁶ | **8.0** | ~10⁻¹⁶ | **8.0** |
| w=3 | 8×8 | 0 | 0 | **0** | 0 | **0** |

**N=4 results:**

| Sector | Block size | Test 1 | Test 2 | Test 3 | Test 4 | Test 5 |
|--------|------------|--------|--------|--------|--------|--------|
| w=0 | 16×16 | ~10⁻¹⁶ | ~10⁻¹⁶ | **~10⁻¹⁶** | ~10⁻¹⁶ | **~10⁻¹⁶** |
| w=1 | 64×64 | ~10⁻¹⁷ | ~10⁻¹⁷ | **11.314** | ~10⁻¹⁷ | **11.314** |
| w=2 | 96×96 | ~10⁻¹⁵ | ~10⁻¹⁵ | **16.000** | ~10⁻¹⁵ | **16.000** |
| w=3 | 64×64 | ~10⁻¹⁷ | ~10⁻¹⁷ | **11.314** | ~10⁻¹⁷ | **11.314** |
| w=4 | 16×16 | ~10⁻¹⁶ | ~10⁻¹⁶ | **~10⁻¹⁶** | ~10⁻¹⁶ | **~10⁻¹⁶** |

In both cases more than ten orders of magnitude separate "exact zero" (≤ 10⁻¹⁵) from "broken" (≥ 8). The break is structural, not perturbative.

The off-diagonal blocks, the field cases, YZ and the census of the 36 two-term Hamiltonians below come from [`veffect_residual_all_blocks.py`](../simulations/veffect_residual_all_blocks.py), which reads every block exactly.

**The break norms are a string count.** One unit XY bond in an N-site chain gives

```
‖M|(w,w)‖² = 32 · C(N−2, w−1) · 2^(N−2)
```

(8, 8√2, 16 and 8√2 in the tables; computed exactly at N = 3, 4, 5). The residual acts on the bond's two sites and as the identity on the N−2 spectators, so a diagonal block is the bond's own two-site residual times the spectator strings of the matching weight. On the bond the extreme diagonal blocks vanish (immunity at N = 2) and its blocks (0,2), (2,0) shift w by two and never land on the diagonal, which leaves its weight-1 block, of squared norm 32; the spectators then supply w−1 transverse sites out of N−2, with two letters to choose at each of the N−2 sites, C(N−2, w−1)·2^(N−2) strings.

## Π-pairing structure

Π acts per site as I↔X (sign 1), Y↔Z (sign i). On total XY-weight, this maps w → N − w. The Π-paired sectors:

- **N=3:** (w=0 ↔ w=3), (w=1 ↔ w=2)
- **N=4:** (w=0 ↔ w=4), (w=1 ↔ w=3), (w=2 ↔ w=2 self-pair)

The script also reads the Π-paired off-diagonal blocks (w_row + w_col = N). The extreme pair (w=0, w=N) is zero for every tested H, since no bilinear moves w by N ≥ 3. The boundary pair follows the parity of the move: XY moves w by 0 or ±2, so in Tests 3 and 5 the pair (w=1, w=2) at N=3 stays zero, while (w=1, w=3) at N=4 breaks at 11.314, as its diagonal neighbours do.

## The structural reading

The palindrome-relation (*) splits into two independent facts:

**(I) Extreme-block immunity.** The diagonal blocks of the pair (w=0, w=N) are fixed by Π alone. For every two-body bilinear H, Π·L·Π⁻¹ on these blocks returns exactly the negative of L, shifted by 2Σγ ([the proof](../docs/proofs/PROOF_ZERO_IMMUNITY.md)). The reasons:

- w=0 strings (all I, Z): commute with the Z-dephasing operators; the dissipator is zero on this sector.
- w=N strings (all X, Y): each anticommutes with every Z-dephasing operator; the dissipator gives the uniform decay rate −2Σγ.
- Π pairs these two extreme behaviors exactly.
- The Hamiltonian has no part inside either sector. A bond letter X or Y flips its site's weight bit, so only a bond built of I and Z keeps a string in w=0 or in w=N; among the bilinears that is ZZ, which commutes with every I/Z string and anticommutes twice, so commutes, with every all-transverse string.

**(II) Boundary-block dependence.** Inside 0 < w < N, whether every term is truly (even Y count, even Z count) controls whether Π·L·Π⁻¹ matches −L − 2Σγ·I. If every term of H is truly, the palindrome holds in every block. A non-truly term breaks it in the blocks it reaches: two X/Y letters (XY) reach the diagonal and the blocks two weights apart, one X/Y letter (YZ) only the blocks one weight apart. Even total Y+Z alone (bit_b parity) is not enough: YZ breaks (*).

The **break magnitude** (8.0 at N=3, 11.314 ≈ 8√2 and 16.0 = 8·2 at N=4) is the string count above, scaled by the strength of the non-truly term.

## Why this matters

**For V-Effect:** On the residual's diagonal the operator-level break behind the "14 of 36 mode-pairs break" finding is not distributed across all weights but **strictly confined to boundary sectors**, while off the diagonal it couples the boundary sectors to the extremes. The boundary sectors are "where chemistry happens" in HIERARCHY_OF_INCOMPLETENESS terms (orphaned valence modes, the carbon-like incomplete shells). The extreme sectors are the inert cores.

**For Heisenberg-form selection:** The C²⊗C² parity argument from PRIMORDIAL_QUBIT says Heisenberg is the unique both-parity-even 2-body bilinear. This finding does not single Heisenberg out by itself: **every truly coupling (even Y and even Z count per term) preserves the palindrome in every block through V-Effect transitions**, Heisenberg, XX + YY and XX + ZZ alike (Tests 1, 2, 4 in the tables above). Coupling forms with non-truly terms (XY, ZX, and also the bit_b-even YZ) couple sites but break the palindrome, XY on the boundary diagonal, ZX and YZ off it.

**For inheritance up the level stack:** When V-Effect bridges Level 0 to Level 1, the algebra that gets transferred is the algebra of the boundary sectors. If the bridge is truly (Heisenberg, XX + YY, XX + ZZ), the inheritance is clean. If it is not, the bridge breaks the boundary palindrome and that level becomes algebraically inconsistent with the next one above. This is the structural reason why truly coupling is the form that propagates: it does not damage what it inherits.

## What this does and does not establish

**Establishes:**

- The diagonal blocks (w=0, w=0) and (w=N, w=N) of M vanish for every two-body bilinear H; numerically at N=3 and N=4 to machine precision (10⁻¹⁵), analytically in [the zero-immunity proof](../docs/proofs/PROOF_ZERO_IMMUNITY.md). A single-site Z field breaks them.
- On the diagonal, the break is confined exactly to 0 < w < N and is discrete in magnitude, a string count.
- Off the diagonal it is not confined: the blocks joining w=0 and w=N to the boundary sectors carry it (4√2 at N=3, 8 at N=4, Tests 3 and 5).
- Among the Hamiltonians run here (the five tests, the 36 two-term census, and XY − YX and XX + YY + XY at N=3), the diagonal blocks break exactly when H contains a non-truly term with two X/Y letters (Tests 3 and 5 break; Tests 1, 2, 4 do not); a non-truly term with one X/Y letter, such as YZ, never reaches them and breaks only the blocks one weight apart.
- The Π-paired off-diagonal blocks (w_row + w_col = N) break only where the move's parity lets the term reach them: zero at (1,2) for N=3, broken at (1,3) for N=4.

**Does not establish:**

- Scaling of the numerical tables beyond N=4. The immunity proof and the parity argument hold at every N; the tables and the self-pair (w=N/2, w=N/2), which breaks at 16.0 at N=4 when the boundary breaks, were not run at N ≥ 5. The string-count law was checked at N = 5.
- Connection to specific atomic-level observables. The "boundary modes" inheriting to Level 1 valence chemistry is an asserted structural correspondence; making it operationally testable is open work.

## References

- [the V-Effect palindrome](V_EFFECT_PALINDROME.md): the original numerical finding (14 of 36 mode-pairs break, boundary modes orphan) that this localizes structurally.
- [the zero-immunity proof](../docs/proofs/PROOF_ZERO_IMMUNITY.md): the analytical proof of the extreme-block immunity, with where it stops (Z fields, the off-diagonal blocks, terms with one X/Y letter).
- [the bit-b parity symmetry proof](../docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md): the analytical proof [L, Π²] = 0 for any N.
- [the primordial qubit](../hypotheses/PRIMORDIAL_QUBIT.md): the C²⊗C² parity structure (bit_a, bit_b) that determines which Hamiltonian terms preserve which Z₂.
- [the hierarchy of incompleteness](../docs/HIERARCHY_OF_INCOMPLETENESS.md): the level-stack picture into which this finding embeds.
- Simulation: `simulations/veffect_weight_immunity.py`.
