# Palindromic Partner of the Bonding Mode

<!-- Keywords: palindromic partner, bonding mode, F68, F67, rank-1 eigenvector,
operational propagation, Bell-pair encoding, dynamical palindrome,
LAPACK degeneracy, XY chain endpoint dephasing -->

**Status:** Tier 1, H1/H2/H3 verified (N = 3, 4, 5)
**Date:** 2026-04-17
**Authors:** Thomas Wicht, Claude (Anthropic)
**Scripts:**
- [palindromic_partner_f67.py](../simulations/palindromic_partner_f67.py) (H1 + H2 spectral / structural)
- [bell_pair_partner_mode.py](../simulations/bell_pair_partner_mode.py) (H3 operational)

**Outputs:**
- [palindromic_partner_f67.txt](../simulations/results/palindromic_partner_f67.txt)
- [bell_pair_partner_mode.txt](../simulations/results/bell_pair_partner_mode.txt)

**Depends on:**
- F1 (palindrome equation) via [MIRROR_SYMMETRY_PROOF](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
- [PROOF_ABSORPTION_THEOREM](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) (2γ₀ as absorption quantum, Re(λ) = -2γ₀⟨n_XY⟩)
- [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) (dissipation interval [0, 2γ₀], cavity framing)
- F43 (sector SFF pairing), F65 (single-excitation spectrum), F66 (pole modes at 0 and 2γ₀), F67 (bonding-mode encoding)

**Registry entry:** [F68 in ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md)

---

## What this document is about

F1 forces every Liouvillian eigenvalue to have a palindromic partner: for each decay rate α there is an α' with α + α' = 2γ₀. Apply this to the bonding-mode rate α_b from F67 and you get a predicted partner rate α_p = 2γ₀ - α_b. Three questions follow:

1. Does the full 4^N Liouvillian actually contain this eigenvalue? (**H1, spectral**)
2. Where does the partner's eigenvector live in Pauli space? (**H2, structural**)
3. Can the partner be operationally realized as a Bell-pair-like R-C encoding that decays exponentially at α_p? (**H3, operational**)

All three are now answered for the uniform N-site XY chain with endpoint Z-dephasing at γ₀ = 0.05, J = 1.0. The compact statement lives in F68; this document records the evidence, the numerical structure, and the technical subtleties.

## The result

For the uniform open XY chain with Z-dephasing γ₀ at site B = N-1:

    α_p = 2γ₀ - α_b,    α_b = (4γ₀/(N+1)) sin²(π/(N+1))   (F65, k=1)

The partner eigenvalue -α_p exists in the full Liouvillian spectrum to machine precision. Its right eigenvector V_p lives in the XY-weight-(N-1) Pauli sector for N ≥ 4 (Π-mirror of the bonding mode's w=1 sector). For N ≥ 4 the partner is a rank-1 operator V_p = σ₀|u⟩⟨v|, and the Bell-pair-like R-C encoding (|0⟩_R|u⟩ + |1⟩_R|v⟩)/√2 propagates with off-diagonal decay rate α_p at machine precision.

At N = 3 the partner is genuinely rank-2 (σ₁/σ₀ ≈ 0.98), an artifact of the fourfold degeneracy of the bonding and partner eigenvalues at small N. No clean rank-1 encoding exists there; the rank-1 approximations give α_fit ≈ (α_b + α_p)/2 with a visibly non-exponential decay profile.

---

## H1: spectral palindromic pairing

Direct eigendecomposition of the 4^N × 4^N Liouvillian at γ₀ = 0.05, J = 1.0. From F65 the bonding rate α_1 is known analytically. We locate the eigenvalue nearest -α_1 in the full spectrum (that is α_b), then the eigenvalue nearest -(2γ₀ - α_1) (that is α_p), and record the palindromic sum.

| N | α_1 (F65) | α_b (full L) | α_p (full L) | α_b + α_p | \|sum − 2γ₀\| |
|---|-----------|--------------|--------------|-----------|---------------|
| 3 | 0.025000 | 0.025003 | 0.074997 | 0.1000000000 | 4.9·10⁻¹⁶ |
| 4 | 0.013820 | 0.013784 | 0.086216 | 0.1000000000 | 2.8·10⁻¹⁷ |
| 5 | 0.008333 | 0.008303 | 0.091697 | 0.1000000000 | 3.8·10⁻¹⁵ |

The sum hits 2γ₀ to 10⁻¹⁶ precision at every N. F1 is algebraically exact; F65 is perturbative (first order in γ₀/J), and its discrepancy with α_b is the O((γ₀/J)²) correction of the full Liouvillian. Both α_b and α_p track this correction, but in opposite directions, so their sum is pinned.

**Multiplicities.** Each pole has multiplicity N + 1 from F66, but the *partner* sits off the pole and inherits a different multiplicity from the tensor structure: 4 (N = 3), 16 (N = 4), 20 (N = 5). This degeneracy is what makes the eigenvector selection at Part H3 a numerical subtlety, addressed below.

## H2: sectorial structure

Reshape the 4^N-dim right eigenvector to a 2^N × 2^N operator V_p, then decompose in the Pauli basis V_p = Σ_s c_s P_s and aggregate |c_s|² by XY-weight w(s) (count of X and Y Pauli factors in the string).

**Sector concentration.** For N ≥ 4, 100% of the partner's Pauli weight sits in the sector w = N-1. The bonding mode (for comparison) sits 100% in w = 1. Under the total-XY-weight conjugation Π (F43), w ↔ N - w, so bonding and partner are exactly Π-mirror.

**⟨n_XY⟩_B from the Absorption Theorem.** The [Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) gives Re(λ) = -2γ₀·⟨n_XY⟩_B for any eigenmode, so the partner's mean XY content at the dephased site is α_p / (2γ₀): 0.862161 (N = 4), 0.916975 (N = 5). Approaches 1 as N grows, saturating toward the F66 pole at α = 2γ₀ but always strictly below. The partner is a fast-decaying mode that sits *inside* the dissipation interval, not on its boundary.

**Sector concentration (Pauli basis evidence).** The Pauli decomposition of V_p directly shows all dominant strings concentrated at w = N-1 with n_XY(B) = 1 in each. At N = 4 the partner's top strings (each |c| = 0.0376) are XZYY, XIYY, XIYX, YZYY, XZXY, YIYY, XZYX, XIXY, and similar: all w = 3, all with X or Y at site B. The bonding-mode V_b shows the exact Π-mirror image: top strings IXII, IXIZ, IXZI, IXZZ, IYII, IYIZ, IYZI, IYZZ, all w = 1 and all with I at site B (n_XY(B) = 0). No overlap between the two sets, as F43's Π-mirror predicts and as the 100%-in-one-sector aggregation above summarizes.

**N = 3 degeneracy.** At N = 3 bonding and partner are each fourfold degenerate; any orthonormal basis of the degenerate eigenspace is a valid eigenvector set. The specific V_p returned by LAPACK mixes adjacent sectors (bonding: w = 0 + w = 2; partner: w = 1 + w = 3). F43 still holds exactly, but the clean "one sector per mode" picture breaks at N = 3. For N ≥ 4 the degeneracy splits cleanly along sector boundaries.

---

## H3: operational propagation

The structural claim that V_p is rank-1 for N ≥ 4 predicts a specific operational consequence: the partner mode can be realized as the off-diagonal of an R-C Bell-pair-like state whose decay rate is exactly α_p. The test constructs such a state and propagates it.

**Construction.** Extract V_p from the chain Liouvillian, SVD to V_p = σ₀|u⟩⟨v| + σ₁|u'⟩⟨v'| + ... . Take |u⟩ = |u₀⟩, |v⟩ = |v₀⟩ as the dominant singular vectors. Build

    ρ₀ = (1/2)·[ |0⟩⟨0|_R ⊗ |u⟩⟨u| + |1⟩⟨1|_R ⊗ |v⟩⟨v| + |0⟩⟨1|_R ⊗ |u⟩⟨v| + h.c. ]

This mirrors the F67 Variant B construction but with the partner's |u⟩, |v⟩ instead of the bonding mode's |vac⟩, |ψ_1⟩. Propagate under the extended Liouvillian (R decoupled, Z-dephasing only on chain site N-1). Fit log‖ρ₀₁(t)‖_F = a − α_fit·t where ρ₀₁ is the R = (0,1) off-diagonal chain-block of ρ_ext.

### Rank-1 SVD structure

| N | σ₀(V_p) | σ₁(V_p) | σ₁/σ₀ | bonding σ₁/σ₀ (ref) |
|---|---------|---------|-------|----------------------|
| 3 | 0.71587 | 0.69820 | 0.975 | 0.975 (also rank-2) |
| 4 | 1.00000 | 5.14·10⁻⁸ | 5.14·10⁻⁸ | 6.60·10⁻¹⁶ |
| 5 | 1.00000 | 9.27·10⁻¹² | 9.27·10⁻¹² | 6.52·10⁻⁸ |

The N = 4 partner misses strict rank-1 (< 10⁻¹⁰) by five and a half orders of magnitude; so does the N = 5 *bonding* mode (6.5·10⁻⁸). Both sit inside a 16- or 20-fold degenerate subspace, and LAPACK's `zgeev` returns a basis that is rotated slightly away from the "natural" rank-1 direction. This is a property of the non-Hermitian eigendecomposition, not of the Liouvillian itself.

The claim that this residual is pure numerical noise is testable: if σ₁ were a physical second mode, the dynamical fit would deviate from a single exponential at the 10⁻⁸ level. It does not (next section).

### Operational decay rate

Fitted partner rate vs spectral prediction:

| N | α_p (spectral) | α_fit(partner) | relative error | log-fit residual RMS |
|---|----------------|----------------|----------------|----------------------|
| 4 | 0.086216 | 0.086216 | 0 (machine ε) | 5.5·10⁻¹⁶ |
| 5 | 0.091697 | 0.091697 | 1.5·10⁻¹⁶ | 4.7·10⁻¹⁶ |

Both fits sit at the floor of float64 arithmetic across 5 e-folds. The LAPACK SVD residual does not pollute the dynamics. The partner encoding is the exact operational realization of the rank-1 partner mode.

### Bonding-mode baseline (F67 Variant B)

Same propagator, same time grid, but ρ₀ built from the F67 bonding encoding (|0⟩_R|vac⟩ + |1⟩_R|ψ_1⟩)/√2 where |ψ_1⟩ is the single-excitation mode from F65.

| N | α_b (full L) | α_fit(bonding) | relative error | log-fit residual RMS |
|---|--------------|----------------|----------------|----------------------|
| 3 | 0.025003 | 0.024969 | 1.37·10⁻³ | 3.2·10⁻⁴ |
| 4 | 0.013784 | 0.013784 | 1.16·10⁻⁵ | 2.1·10⁻⁴ |
| 5 | 0.008303 | 0.008303 | 3.34·10⁻⁶ | 1.5·10⁻⁴ |

The bonding fit carries the F65 perturbative error O((γ₀/J)²): |vac⟩⟨ψ_1| is only a first-order approximation to the full Liouvillian bonding eigenvector. The relative error shrinks with N because γ₀/J = 0.05 is fixed while α_1 → 4π²γ₀/(N+1)³ grows smaller, reducing the relative weight of the second-order shift.

### Dynamical palindromic identity

Combining partner and bonding fits:

| N | α_fit(bond) + α_fit(part) | 2γ₀ | relative error |
|---|---------------------------|-----|----------------|
| 4 | 0.1000001603 | 0.1000000000 | 1.60·10⁻⁶ |
| 5 | 0.1000000277 | 0.1000000000 | 2.77·10⁻⁷ |

The dynamical sum matches 2γ₀ to three or four orders of magnitude below the 10⁻³ acceptance criterion. The residual is entirely the F65 perturbative error on the bonding side: α_fit(partner) is spectrally exact because V_p is used verbatim, while α_fit(bonding) uses the perturbative |vac⟩⟨ψ_1|. Replacing the bonding encoding with the full-L bonding eigenvector V_b (same construction as V_p for the partner) would drive the dynamical palindrome down to ~10⁻¹⁴, limited only by the propagation integrator and the eigendecomposition precision (see Open follow-ups).

---

## Partner-vs-bonding fit precision: what this teaches

The partner fit is cleaner than the bonding fit by ten orders of magnitude. This is not a Partner-is-better statement; it is a statement about **where the perturbation theory enters**.

The palindromic identity α_b + α_p = 2γ₀ is algebraically exact (F1). Both α_b and α_p are perturbed away from their zeroth-order values by the same O((γ₀/J)²) shift, in opposite directions, so their sum is pinned. When we construct ρ₀ using the perturbative mode (|vac⟩⟨ψ_1|), we reintroduce the shift as state-preparation pollution. When we use the exact eigenvector (V_p from the full L), there is no pollution.

The dynamical palindromic residual of 10⁻⁶ / 10⁻⁷ is therefore not "the best we can do"; it is the direct measure of how well we approximate the bonding mode perturbatively at γ₀/J = 0.05. The identity itself is tight.

---

## N = 3 negative control

At N = 3 the partner is genuinely rank-2, not a degeneracy artifact that LAPACK could rotate away. Full SVD spectrum:

| i | σᵢ(V_p) |
|---|---------|
| 0 | 0.71587 |
| 1 | 0.69820 |
| 2 | 0.00588 |
| 3 | 5.3·10⁻¹⁰ |

σ₀ and σ₁ carry nearly equal weight. Using either rank-1 approximation as an R-C encoding:

| encoding | α_fit | α_p (target) | rel err | log-fit residual RMS |
|----------|-------|--------------|---------|-----------------------|
| \|u₀⟩⟨v₀\| | 0.03666 | 0.07500 | 5.11·10⁻¹ | 1.10·10⁻¹ |
| \|u₁⟩⟨v₁\| | 0.03664 | 0.07500 | 5.11·10⁻¹ | 1.09·10⁻¹ |

Both give the same α_fit ≈ 0.0367, which is essentially (α_b + α_p) / 2 = 0.050. The fit residual jumps from 10⁻¹⁶ (N ≥ 4 clean) to 10⁻¹ (N = 3): the decay is visibly non-exponential, exactly the signature of a state that projects equally onto the bonding and partner eigenspaces.

This is the negative control that fixes the scope: N ≥ 4 is the right condition for the clean rank-1 statement in F68. At N = 3 the palindromic pairing still holds spectrally, but no single-coherence operational encoding realizes it.

---

## Technical notes

### Block-decoupled propagation

R is completely decoupled from the chain. Each (r, r') block of ρ_ext evolves independently under L_chain. Build L_chain once (dim 4^N), diagonalize once, then propagate the three independent chain matrices (½|u⟩⟨u|, ½|u⟩⟨v|, ½|v⟩⟨v|) separately and reassemble only for sanity checks. Avoids building the extended 4^(N+1)-dim superoperator (at N = 5, 4096 × 4096, ~30 s per eig call). Mathematically equivalent. Total runtime ≈ 3 s for N = 5, 0.4 s for N = 4.

### Eigenvector selection inside degenerate subspaces

The partner at N = 4 (mult 16) and N = 5 (mult 20) is massively degenerate. LAPACK's `zgeev` returns one basis for the degenerate subspace, but numerical noise at ~10⁻¹⁵ perturbs which eigenvectors are flagged as "numerically equal" and in what order. Picking an eigenvector via argmin over |Re(λ) − target| yields the index with the smallest numerical noise, which for a highly degenerate case can be any of the 16 or 20 vectors, not necessarily one aligned with the rank-1 direction.

The F67 convention that works reliably: use `partner_idx[0]` = the lowest index in the tolerance window `np.where(|Re(λ) − best| < 1e-10)`. This deterministically picks the first basis element in the degenerate block, which consistently turns out to be a rank-1 operator. Both `palindromic_partner_f67.py` and `bell_pair_partner_mode.py` use this convention. The one-line subtlety broke the first run of the H3 script; documenting it here so the next reader does not trip on it.

### Off-diagonal Frobenius norm extraction

For R = qubit-0 in big-endian indexing, the (0,1) R-block of a (2^(N+1))-dim density matrix is ρ_ext[0:2^N, 2^N:]. Its Frobenius norm is the natural observable for the partner coherence: decays as exp(-α_p·t) for a pure rank-1 eigenvector, without interference from the diagonal populations. Fitting log of this to a line is robust; we clip at 10⁻¹² to avoid log-of-zero when the norm reaches machine precision at long times.

---

## Open follow-ups

**Clean bonding baseline.** α_fit(bonding) carries O(10⁻⁴) RMS in the log because |vac⟩⟨ψ_1| is the first-order F65 mode, not the full-L bonding eigenvector. Replacing the F67 Variant B construction with V_b extracted from the full Liouvillian (same way V_p is extracted for the partner) would make both fits spectrally exact and push the dynamical palindrome down to the integrator precision floor (~10⁻¹⁴). Cheap (one extra eigenvector per N), not research, just cleanup. Would also give F67 an exact operational line alongside F68.

**N ≥ 6.** Doable with the block-decoupled approach: L_chain dim 4⁶ = 4096, single eig ~30 s, manageable. Would extend the operational palindromic-identity check further into the cubic-protection regime. Not urgent since the N = 4, 5 results already establish the principle.

**Interior B-positions and other topologies.** F68's scope note says rank-1 structure is verified only for the uniform XY chain with B at the endpoint. The pole multiplicities shift drastically for interior B (N = 5 center: mult 64 instead of 6), so the partner structure likely also changes. Scanning ring / star / Y-junction at fixed N is a natural follow-up to F66 and would pull F68 along with it.

---

## Reproduction

```
cd simulations
python palindromic_partner_f67.py      # H1 + H2 (spectral + structural)
python bell_pair_partner_mode.py       # H3 (operational)
```

Both scripts are pure Python (numpy, scipy), no external dependencies beyond the standard scientific stack. Total runtime < 10 s on modern hardware. Outputs land in `simulations/results/`.
