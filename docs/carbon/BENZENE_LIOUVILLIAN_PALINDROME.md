# Selected C₄/C₆ Ring Liouvillians: F1 Palindrome and Bond-Jump Comparison

**Date:** 2026-05-22
**Authors:** Tom + Claude
**Status:** F1 proves the selected XX+YY + all-site-Z model palindrome; the
selected bond-jump comparison is numerically evaluated at C₄ and C₆. The
material-carbon degree of freedom, β-to-J convention, bath channel and rate, and
therefore material γ, T₂, and Q remain unassigned (see [README](README.md#conditional-c4-and-c6-working-model)
and [Q audit](../Q_BELONGS_TO_NO_SUBSTANCE.md)).
**Script:** [`simulations/carbon/benzene_liouvillian_palindrome.py`](../../simulations/carbon/benzene_liouvillian_palindrome.py)
**Tested:** selected C₄ and C₆ spin rings.
**Answers:** [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md) open
question 1 / [README](README.md) open question 2.

---

## The question

[Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md) showed that the
**closed-system** Hückel MO spectrum of an alternant hydrocarbon is palindromic
(Coulson-Rushbrooke, 1940) and that this is a Z₂ palindrome of the same shape as
the framework's F1, with a different trigger. But Coulson-Rushbrooke is the closed-system, Hamiltonian-level
statement. F1 is the **open-system** claim: the Liouvillian spectrum of an XY-class
Hamiltonian under Z-dephasing is closed under λ → −λ − 2Σγ, palindromic about the
centre −Σγ.

Does the selected C₄/C₆ open-system model satisfy F1, and what changes when its
local jump is replaced by the selected bond jump?

## The selected model and jump comparison

Choose a site occupation `n_l = (I − Z_l)/2` and the XX+YY ring as a free-hopping
mathematical model. The Hückel `β → J` relation is a Tier-2 structural translation
within that selected model, not a material-carbon assignment. F1 applies to the
selected XX+YY/XY spin graph with all-site local Z-dephasing; it does not identify
the physical carbon degree of freedom or a molecular bath.

Two deliberately selected jump channels make the scope comparison:

- **Local density / Z.** If the selected jump is the local density `n_l`, then
  `D[n_l] = ¼·D[Z_l]`; equivalently the script evaluates local Z-dephasing. F1
  applies to this selected channel.

- **Bond B.** The selected jump is `B_b = X_aX_b + Y_aY_b`, a two-body operator.
  This is outside F1's Z-dephasing hypothesis.

These two selected channels are not an exhaustive classification of physical bath
couplings and do not establish a material Holstein or Peierls bath.

## Result

The selected Liouvillian `L = −i[H,·] + Σ D[√γ·jump]` was built for C₄ and C₆.
For local-Z rows the spectrum is tested against the strict F1 involution
λ → −λ − 2Σγ about the predicted centre −Σγ; for bond-B rows, which have no F1
prediction, it is tested against the most generous reflection (about the spectrum
mean). The residual is the largest distance from a reflected eigenvalue to the
nearest actual one.

| Ring | Coupling | Reflection centre | Residual | Palindrome |
|------|----------|-------------------|----------|------------|
| C₄ selected ring | local Z, D[Z_l] | −Σγ = −4 | 3.5 × 10⁻⁸ | holds |
| C₄ selected ring | bond B, D[B_b] | n/a (no palindrome) | ≈ 11 | broken |
| C₆ selected ring | local Z, D[Z_l] | −Σγ = −6 | 1.2 × 10⁻⁷ | holds |
| C₆ selected ring | bond B, D[B_b] | n/a (no palindrome) | ≈ 14 | broken |

For the selected local-Z channel, the strict-F1 residual is at the numerical floor
of the non-Hermitian eigendecomposition, with centres −4 and −6 at `γ = 1` per
site; F1 proves this model palindrome exact. For the selected bond-B channel,
even the most generous reflection leaves a residual of order the spectral width.

## Answer within the selected comparison

For the selected XX+YY plus all-site-local-Z model, F1 supplies the exact
palindrome. For the selected bond-B jump, the listed finite C₄/C₆ simulations do
not show that F1 symmetry. This compares two chosen jump operators; it neither
classifies a molecular environment nor says that either channel is exhaustive.

Coulson-Rushbrooke/K and F1 remain structural siblings on distinct objects:
Coulson-Rushbrooke and K pair a single-particle Hückel hopping problem, whereas
F1 pairs a Liouvillian operator space. K is the formal framework partner of the
Hückel bipartite condition; it fails on odd rings, while F1 remains
topology-blind within its selected Hamiltonian-and-Z-channel scope. No material
embedding follows from their shared pairing shape (see [the framework lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md)).

## Structure of the selected bond-B break

The selected bond-B jump gives a clear numerical break in its C₄/C₆ γ-scan.

**It is linear in γ.** For small γ the palindrome residual, measured against the
most generous reflection (the spectrum mean), grows as 4(N+2)·γ: ratio 24.00 at
C₄, 32.00 at C₆, log-log slope 1.00. Above γ ≈ 0.3 the growth bends sub-linear. As
γ → 0 the Liouvillian returns to −i[H,·], palindromic about 0, and the residual
vanishes.

**It is total.** At γ = 1 on C₆, even under the best-fit reflection centre not one
of the 4096 eigenvalues has a reflected partner within 10⁻⁶. The break is not
carried by a few outliers: 95% of eigenvalues sit more than 1% of the spectral
width from their nearest partner, and the per-eigenvalue gap rises smoothly across
the spectrum. No real centre rescues it; the only surviving involution is complex
conjugation, the conjugate-pairing every Lindbladian spectrum has, not the F1
mirror. The Lindblad zero mode (the steady state) is intact.

**It is outside F1's jump premise.** The selected bond-B Liouvillian does not
close under λ → −λ − 2Σγ (F1-strict residual 16.5 at C₄, 22.5 at C₆, against
the 10⁻⁷ floor of the selected local-Z case). A like-for-like selected
depolarising-noise scan on the same ring is also γ-linear (slope 1.00, residual
(2/3)·Σγ). These are model-channel comparisons, not a statement about a
material bath.

This says something scoped about the selected bond operator
`B = X_aX_b + Y_aY_b`. It is a truly-class Hamiltonian bilinear in the selected
`H = Σ_b B_b`, yet `D[B]` is outside F1's Z-dephasing hypothesis. F87 classifies
Hamiltonian terms; it does not make the same operator an F1-compatible dissipator
jump. The simulation establishes that model distinction only.

## Framework-vocabulary translation

| Hückel/framework model item | Selected framework object | Status |
|------------------------------|---------------------------|--------|
| free Hückel hopping on C₆ graph | XX+YY spin ring, N=6 | Tier 2 structural translation |
| selected site occupation and density jump | local Z-dephasing, `D[n_l] = ¼·D[Z_l]` | exact once `n_l` and the jump are selected |
| selected bond operator | bond-jump dephasing `D[B_b]` | finite comparison channel |
| local-Z selected model | F1, proven for XY-on-graph + Z-dephasing | Tier 1 model result |
| bond-B selected model | outside F1; numerical C₄/C₆ comparison | selected-model simulation |

## Open follow-ups

- A mixed selected `D[Z] + D[B]` jump model: what interpolation does its finite
  spectrum show? This would remain a model-channel scan until a material channel
  is independently specified.
- F98 ((N+2)/[4(N+1)] → 1/4 long-time bridge): F98 needs its specified
  KIntermediate initial state, a magnetization-conserving H, and all-site
  Z-dephasing. [The C₄/C₆ selected-ring instances](BENZENE_F98_LONG_TIME.md)
  give 3/10 and 2/7; they are not a material-bath assignment.
- The [selected-model clock comparison](FROST_CIRCLE_AS_THE_CLOCK_FACE.md) reports
  frequencies and Q* values only after choosing its XX+YY Hamiltonian and
  Z-dephasing rate; it supplies no carbon `γ`, `T₂`, or Q.

## Anchor

- Scripts: [`simulations/carbon/benzene_liouvillian_palindrome.py`](../../simulations/carbon/benzene_liouvillian_palindrome.py)
  (the selected channel comparison), [`simulations/carbon/peierls_break_structure.py`](../../simulations/carbon/peierls_break_structure.py)
  (the selected bond-jump γ-scan)
- Companion doc: [Benzene Hückel through the Framework Lens](BENZENE_HUCKEL_FRAMEWORK_LENS.md)
  (the closed-system half), [README.md](README.md)
- Framework anchors: [F1 palindrome](../ANALYTICAL_FORMULAS.md#f1-palindrome-equation-tier-1-proven),
  [`compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs`](../../compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs),
  [`docs/proofs/MIRROR_SYMMETRY_PROOF.md`](../proofs/MIRROR_SYMMETRY_PROOF.md)
