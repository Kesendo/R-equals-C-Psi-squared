# The ζ² anti-protection law: the chiral symmetry lifted to the Floquet step

**Date:** 2026-07-16
**Status:** derived law under the F129 surface (no own F number, by decision; precedent: the F129 family inventory). Gate: [`simulations/zeta2_anti_protection.py`](../../simulations/zeta2_anti_protection.py), exit 0.
**Owner of the static symmetry:** `ChiralKClaim` (class BDI) + [PROOF_K_PARTNERSHIP](PROOF_K_PARTNERSHIP.md). This document adds only the Floquet-step lift and its consequence for the flown F129 Ramsey fringe.
**General home (added 2026-07-16):** this law is Theorem B, the antiunitary column, of the mirror's order-sorting law, [PROOF_MIRROR_ORDER_SORTING](PROOF_MIRROR_ORDER_SORTING.md): the response-order sorting by q·σ_eff, of which the pair difference here is the even cell and the pair sum the odd cell.

Before the algebra, feel the shape of it: the same mirror that makes the two
notes equal is the one that forbids the noise from telling them apart at first
order, and then, one order deeper, the very same mirror turns around and pushes
them apart, twice one branch's shift, exactly where the trained intuition
expects the protected pair to feel nothing. Protection and anti-protection
are not two symmetries. They are one symmetry read at odd and at even order.

## 1. Setting and vocabulary

The instrument is the flown F129 standing-fringe experiment
([IBM_F129_RAMSEY_FRINGE](../../experiments/IBM_F129_RAMSEY_FRINGE.md),
Confirmation 24). N = 8 sites, open chain; the comb index is **n = N + 1 = 9**
(continuum mode energies cos(kπ/9), k = 1..8; the step's Floquet phases sit
near 2θ·cos(kπ/9)). One Trotter step is the brickwork XX
propagator U(θ) at θ = 0.5 (even bonds then odd bonds; one-particle blocks
B(θ) with cos θ on the diagonal and i·sin θ off it). Floquet modes are the
eigenvectors of U(θ); quasi-energies are the eigenphases. The A0 arm prepares
the superposition of the two 3-magnon Slater branches

τ = modes (1, 5, 7), ν = modes (2, 4, 8) = (9−7, 9−5, 9−1),

a **Θ-mirror pair**: ν is the image of τ under the mode mirror k ↔ N+1−k.
The perturbation is the always-on ZZ crosstalk, injected per step as
rzz(2πζ_i·τ_step) on each bond: an **occupation-diagonal** phase operator A
(diagonal in the site-occupation basis), with site-dependent rates ζ_i.
The observable is the fringe phase difference of the pair; per step it is the
quasi-energy difference Δ(ζ) = θ_τ(ζ) − θ_ν(ζ).

## 2. What was already owned

- **The static symmetry:** K = diag((−1)^ℓ) with K·H·K = −H for any NN
  tight-binding H, spectrum inversion E_{N+1−k} = −E_k; K is linear (chiral,
  NOT time reversal); with T = complex conjugation the class is BDI and the
  particle-hole line of the table is the antiunitary **Θ = T·K with Θ² = +I**.
  All of this is [PROOF_K_PARTNERSHIP](PROOF_K_PARTNERSHIP.md) and the typed
  `ChiralKClaim`; the mode-comb pairing lives in `XyJordanWignerModes`
  (ChiralPairingResidual ≤ 1e−12) and the triple mirror in
  `JwDispersionStructure`.
- **The first-order certificates (committed 2026-07-15):**
  `chiral_certificate()` in
  [`simulations/f129_ramsey_fringe_design.py`](../../simulations/f129_ramsey_fringe_design.py)
  already certifies, at machine zero, (i) equality of every Z-diagonal
  observable between the two branches, per site and per pair at all ranges,
  with site-dependent couplings (zero first-order drift from ZZ and detuning),
  and (ii) the Slater-Condon branch-mixing zero (the branches differ in three
  orbitals; the operators are at most 2-body).
- **The numerical law (committed 2026-07-15, not derived then):** the 7a gate's
  `zz2_scan` found the A0 slope bias 0.00257·(ζ/3.8 kHz)²·(τ/1.2 µs)² under
  the flown estimator and funded the one-sided budget b_zz2.

New in this document: the lift of Θ from H to the Floquet step, the all-orders
identity, the evenness, the exact factor 2, and the clean-limit coefficient
that the flown constant dresses.

## 3. The lift

**Lemma (Θ commutes with the brickwork step).** Write C = K = diag((−1)^i),
the same single-site diagonal as §2's K (the lemma keeps the letter C for the
one-particle matrix). Each bond block satisfies
B(θ)* = B(−θ), and conjugation by C flips the sign of every
NN off-diagonal element, so C·B(−θ)·C⁻¹ = B(θ). C is a single on-site
diagonal, so it distributes across both brickwork layers and every disjoint
bond independently; open ends and even N are irrelevant. Hence

C·U(θ)*·C⁻¹ = U(θ), i.e. Θ = T·K commutes with U(θ).

Being antiunitary, Θ negates every quasi-energy; on the sorted comb it
realizes k ↔ N+1−k and maps the Slater branch τ = (1,5,7) to ν = (2,4,8).
(Gate G0: residual exactly 0.0; mode antisymmetry 2e−16; branch overlap
1.000000000000.)

**Theorem (the lift with the perturbation).** Let A be ANY occupation-diagonal
real phase operator (site-dependent ZZ of any range, single-site detunings),
and W(ζ) = e^{−iA(ζ)}·U(θ) the perturbed step with A linear in the global
scale ζ. C is diagonal in occupation, so Θ·A·Θ⁻¹ = A, and

Θ·W(ζ)·Θ⁻¹ = e^{+iA(ζ)}·U(θ) = W(−ζ), exactly.

Antiunitary conjugation NEGATES eigenphases: if W(ζ)·ψ = e^{iφ}·ψ, then
W(−ζ)·(Θψ) = e^{−iφ}·(Θψ). That negation is the whole mechanism. So Θ maps
the eigenvector of W(ζ) on the τ branch to an eigenvector of W(−ζ) with
negated eigenphase. For the **isolated branch** (τ is a simple
eigenphase whose nearest neighbour is ν itself, at gap 0.0213, with coupling
exactly zero by Slater-Condon, so the analytic branch is well defined and no
crossing relabels it at the flown scales):

**θ_τ(ζ) = −θ_ν(−ζ), exactly, to all orders in ζ.**

## 4. Corollaries: one symmetry, opposite action at even and odd order

Write θ_τ(ζ) = Σ_m c_m ζ^m and θ_ν(ζ) = Σ_m d_m ζ^m. The identity gives
d_m = −(−1)^m c_m: **odd orders transfer to both branches equally** (they
cancel in the difference: the first-order protection the certificates pin),
and **even orders are exactly opposite** (they double in the difference).
Hence:

- **Evenness:** Δ(ζ) − Δ(0) is an even function of the global sign of ζ
  (gate G4: bias(+ζ) − bias(−ζ) = 0.0 exactly).
- **The exact factor 2:** the leading observable bias is 2·ε⁽²⁾_τ, twice one
  branch's second-order Floquet shift, instead of the zero an untrained
  intuition expects from "the protected pair". That is the anti-protection.
  The factor 2 is exact; the O(ζ⁴) tail is the truncation of θ_τ's even part,
  not of the factor.
- **The mechanism, completed:** gap negation alone is not sufficient. The
  second-order shift is ε⁽²⁾_s = ½·Σ_{m≠s} |A_ms|²·cot((ε_s − ε_m)/2), and the
  branches' shifts are opposite because Θ BOTH negates every gap (cot is odd)
  AND equalizes every coupling modulus |A_{Θm,ν}| = |A_{m,τ}| (Θ-invariance of
  the occupation-diagonal A). The committed one-liner "mirrored Floquet gaps
  make the shifts opposite" named only the first half.
- **The clean-limit coefficient:** at the 7a seed profile and the τ = 1.2 µs
  reference, the ζ → 0 limit of the exact bias is **C_phys = 0.002522**, and
  the PT2 cot formula reproduces it to 3·10⁻⁵ relative (gate G2). This is the
  estimator-free physical value of the flown configuration (this profile,
  θ = 0.5, N = 8; not a universal constant). The flown C2LAW = 0.00257 is C_phys
  dressed by the predicted-V WLS estimator weights (and the earlier 0.00242 by
  the equal-weight padding-era estimator); both dressed values live in the
  instrument doc, the physics value lives here.
- **The sign is NOT fixed by the symmetry.** Θ fixes evenness and the factor
  2; the sign of ε⁽²⁾_τ is a property of the profile and regime (same-sign
  always-on ZZ gives the positive direction; the 40-draw check stays in the
  instrument doc with b_zz2).
- **Consistency bonus:** the unperturbed Trotter split of the pair is
  ε_τ − ε_ν = 2ε_τ (ε_ν = −ε_τ exactly), numerically +0.02127 = the
  PREDICTED A0 drift center of §5 of the fringe doc (+0.0213; the RECORD's
  measured slope is +0.0326, the center plus the budgeted excess) (gate G5).

## 5. The fences (load-bearing)

1. **Θ-mirror pairs ONLY.** Every statement above requires ν to be the mode
   mirror of τ. Non-mirror collision pairs carry FIRST-order ZZ shifts: the
   flown control arms A1 (2,5,6) and A2 (1,6,7) do (Wick values −0.0024 and
   +0.0010 in the fringe doc), and so does the **n = 12 pair
   (1,2,10) ~ (3,5,6)**, which is not a mirror pair. **An n = 12 flight must
   budget ζ¹, not ζ².** This fence is the practical content of the law for any
   future collision lab.
2. **Open chain, NN hopping.** The exact fermionic-compound identification
   (and with it the machine-zero gate statements) relies on the Jordan-Wigner
   map being string-free: open boundary, nearest-neighbour hopping only. A
   ring or beyond-NN hop reintroduces parity strings and the statements need
   revisiting.
3. **Quasi-energy level.** "Exact" throughout means the quasi-energies of the
   ideal step-plus-diagonal-phase model. At the flown-observable level (GHZ
   prep, Givens network, decoherence, readout correction, WLS slope) what
   survives exactly is the first-order protection and the evenness; the
   factor 2 and the coefficient arrive estimator-dressed. The always-on (not
   stroboscopic) character of real ZZ does not break the symmetry (the
   continuous generator H_hop + H_ZZ commutes with Θ the same way), but the
   coefficient value is specific to the rzz-per-step model and timing.
4. **Isolated branch.** The all-orders identity uses the simplicity and
   isolation of the τ eigenphase (verified: nearest neighbour is ν, coupling
   zero). The PT2 helper in the gate is fenced to this Slater-Condon-protected
   pair; the comb globally has exactly degenerate pairs with real coupling
   where the naive formula would silently drop an O(ζ) term.

## 6. Hardware context, honestly sized

The flight (ibm_kingston, job d9br4vmg26ic73dgbgk0, Confirmation 24) measured
the A0 deviation +0.0114 = +3.36σ_a, positive, inside the pre-registered
budgeted window. That is **sign-consistency within a blended budget**, not a
measurement of the ζ² sign: b_zz2 = +0.0043 is about 38 % of the observed
excess, the remainder sits in the two-sided quasi-static draw b_qs = ±0.008
whose sign on one device is a single frozen unknown. The flight banks the
standing fringe (clauses a ∧ b ∧ c), and the excess direction is consistent
with this law in the same-sign regime; no more is claimed.

## 7. Where the pieces live

- This proof + the gate [`simulations/zeta2_anti_protection.py`](../../simulations/zeta2_anti_protection.py) (G0 lift, G1 cited certificate, G2 clean limit + PT2, G3 exact antisymmetry, G4 evenness, G5 Trotter-split consistency; all PASS, exit 0).
- The static owner: [PROOF_K_PARTNERSHIP](PROOF_K_PARTNERSHIP.md) + `ChiralKClaim` (class BDI; K linear, T = conjugation, Θ = T·K the PHS line).
- The instrument constants (b_zz2, the 1.22 estimator uplift, the 2× transferred-ζ policy, one-sidedness, the 40-draw same-sign check): [IBM_F129_RAMSEY_FRINGE](../../experiments/IBM_F129_RAMSEY_FRINGE.md) §5, unchanged and normative for the flight record.
- Typed surface: deferred; when it lands it is a child of `ChiralKClaim` (the Floquet lift), not of the crosstriple family and not a MirrorWorld adoption.
