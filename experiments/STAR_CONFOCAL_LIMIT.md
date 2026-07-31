# The Star Spread: Im_max = J·N/2, and What That Does and Does Not Distinguish

<!-- Keywords: star topology Liouvillian spectrum, optical cavity point focus,
SU(2) Schur-Weyl Heisenberg, imaginary spectral spread, hub-spoke topology,
minimal energy gap connected graph, 24 anchors Q-sweep, R=CPsi2 star -->

**Status:** Tier 1 derived. Formal proof [`PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md`](../docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md), typed as [`StarImMaxBoundClaim`](../compute/RCPsiSquared.Core/Symmetry/StarImMaxBoundClaim.cs).
**Date:** 2026-05-19 (the closed form and the anchors), 2026-07-31 (the scope fences and the minimiser search)
**Authors:** Thomas Wicht, Claude
**Depends on:** [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md),
[Proof: Weight-1 Degeneracy](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md) (the Heisenberg + Schur-Weyl substrate),
[F50 typed claim](../compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs),
[Q-Regime Anchor Map](../docs/Q_REGIME_ANCHORS.md) (the 10-anchor canonical Q-table)

**This page is the reading; the proof is the primary source.** The derivation, the scope fences and the counterexamples live in [the proof](../docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md); this page states what they amount to and places the result in the optical-cavity picture. Everything numeric on both pages is produced by one runnable gate, [`simulations/star_saturation_gate.py`](../simulations/star_saturation_gate.py).

---

## Conventions

Fix these once, because two of them are easy to trip over.

- **The Hamiltonian.** `H = J · Σ_bonds S_i · S_j = (J/4) · Σ_bonds (X_i X_j + Y_i Y_j + Z_i Z_j)`, with **J > 0**. In this sign the fully polarised computational state |0…0⟩ is the **maximum**-energy eigenstate, not the ground state: each bond term has maximum eigenvalue +J/4 on the triplet, and |0…0⟩ is triplet on every bond at once, so it attains all of them simultaneously. Below, "fully polarised extreme" always means this state (or |1…1⟩), and it sits at the top of the spectrum.
- **The dissipation.** Uniform Z-dephasing at rate γ on every site: `D[ρ] = γ Σ_l (Z_l ρ Z_l − ρ)`. `σ := N·γ` is the palindrome center, and `ΔE_max(H)` denotes the spectral spread `E_max − E_min` of H.
- **The ratio.** `Q := J/γ`, the dimensionless coupling-to-dephasing ratio. Note the canonical anchor table [`docs/Q_REGIME_ANCHORS.md`](../docs/Q_REGIME_ANCHORS.md) defines `Q = J/γ₀` against a **fixed** code-convention substrate `γ₀ = 0.05`; the two agree exactly when γ = γ₀, which is the setting of the Q-sweep below. The "Marrakesh convention" is the other setting used here, `J = 1, γ = 0.5` (so `Q = 2`), named for the hardware run the repo calibrated against.

---

## What this document is

A sharpening of [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md). In April 2026 we read the qubit chain under Heisenberg + Z-dephasing as a Fabry-Perot optical cavity: weight sectors are transverse planes, the Hamiltonian is free-space propagation, the degeneracy profile is the beam intensity. That framework was developed for the chain.

Turning it on the star gives a closed form:

    Im_max(star, N, J)  =  J·N/2      for every N ≥ 3 and every (J, γ), at uniform γ

Three statements are tangled together in that line, and keeping them apart is the whole content of this page.

1. **The Casimir gap** (star-specific, exact, every N). `ΔE_max(H_star) = J·N/2`. A fact about H alone: no Liouvillian, no dephasing, no γ.
2. **The saturation** (universal, not star-specific). `max |Im(λ_L)| = ΔE_max(H)` holds for **every** isotropic Heisenberg graph, so meeting one's own spread carries no star information. Chain, ring, complete, an asymmetric six-vertex graph and a disconnected one all do it.
3. **The minimum** (star-specific, searched at N ≤ 6). Among all connected graphs on N sites, the star is the **unique** minimiser of `ΔE_max`, and its value is exactly J·N/2.

The star's distinction is (1) and (3): the smallest Hamiltonian spread, with a clean closed form. It is not (2).

---

## Empirical anchors

Heisenberg J=1, uniform Z-dephasing γ=0.5, σ = N·γ. Dense numpy eigvals at small N, `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock` at N=8.

| N | σ = N·γ | max \|Im(λ)\| | Im/σ | Other topologies (Im/σ) |
|---|---:|---:|---:|---|
| 3 | 1.5 | 1.5000 | **1.0000** | chain=1.000 (= star, isomorphic at N=3); ring=1.000 (triangle) |
| 4 | 2.0 | 2.0000 | **1.0000** | chain=1.183, ring=1.500 |
| 5 | 2.5 | 2.5000 | **1.0000** | chain=1.171, ring=1.247 |
| 6 | 3.0 | 3.0000 | **1.0000** | chain=1.248, ring=1.434 |
| 8 | 4.0 | 4.0000 | **1.0000** | chain=1.281, ring=1.413, K₄ ⊔ P₄=1.342 |

At N=3 the star equals the chain by graph isomorphism (path on 3 sites = star with 2 leaves), and the N=3 ring is K₃; all three coincide. For N ≥ 4 the star is the only one of these with Im/σ = 1 at J = 2γ.

25 distinct (N, Q) anchors over 29 runs: a 24-point Q-sweep, the N=8 point, and four Python re-runs of the Q=2 column the sweep already covers. They agree with J·N/2 to a worst relative deviation of **1.98e-14**. None of them is exactly J·N/2 in floating point, so these are machine-precision agreements, not bit-exact ones.

---

## Q-sweep: the value does not depend on γ (24 anchors)

The 2026-05-19 Q-sweep (`f1_q_sweep_anchor.py`) tests `Im/σ = Q/2` at γ₀ = 0.05 across six Q values. All 24 (N, Q) combinations match to machine precision:

| N \ Q | 0.5 | 1.0 | 1.5 | √3 ≈ 1.732 | 2.0 | 2.5 |
|---|---:|---:|---:|---:|---:|---:|
| pred Q/2 | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **3** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **4** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **5** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |
| **6** | 0.2500 | 0.5000 | 0.7500 | 0.8660 | 1.0000 | 1.2500 |

Five of the six Q values (1.0, 1.5, √3, 2.0, 2.5) are canonical anchors from the Q-regime table; Q = 0.5 is not in that table, sitting between its onset band (0.2 to 0.35) and its peak band (from 1.2).

Two cautions about what this sweep is evidence for. First, `Im/σ = Q/2` and `Im_max = J·N/2` are the same statement: σ = N·γ makes N cancel identically, so the *ratio* Q/2 carries no N-dependence to confirm. The four N rows are still four genuine values of N for `Im_max = J·N/2`; what they are not is four independent confirmations of the dimensionless law, which is one identity. Second, the sweep varies **J** at fixed γ₀ = 0.05 (`f1_q_sweep_anchor.py` sets `J = Q · GAMMA_SUBSTRATE`); γ itself is never varied there. The γ-independence is checked separately in the gate, which holds J = 1 fixed and moves γ across γ ∈ {0.005, 0.05, 0.5, 5.0, 50.0}, four decades, with `Im_max` unchanged at J·N/2 throughout.

The Marrakesh convention J = 2γ is the column Q = 2, where Im/σ = 1. That reading is the Q=2 specialization, not the law.

---

## The result in outline

The full derivation is in [the proof](../docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md); this is the shape of it.

**The gap.** Every star bond touches the hub, so `H_star = J · S_0 · S_L` with `S_L` the total leaf spin, and the Casimir identity gives two levels per `S_L ≥ 1/2` sector split by `J·(S_L + 1/2)`. The maximum at `S_L = (N−1)/2` is `ΔE_max = J·N/2`, between `E_max = +J(N−1)/4` and `E_min = −J(N+1)/4`.

**The bound.** Write `L = K + D` with `K = −i[H, ·]`. Because the jump operators `Z_l` are Hermitian, `D` is self-adjoint in the Hilbert-Schmidt inner product; `K` is anti-Hermitian and normal. For an eigenvector v, `λ = ⟨v, Lv⟩/‖v‖²` with `⟨v, Dv⟩` real, so `Im λ` lies in the numerical range of `K/i`, which for a normal operator is the convex hull of its spectrum, `[−ΔE_max, +ΔE_max]`. Hermiticity of the jumps is what carries this, not the informal "the dissipator adds only real decay": see the scope section for a dissipator that adds only real decay in that informal sense and violates the bound outright.

**The attainment.** Equality forces v into the extremal `K`-eigenspace *and* makes v an eigenvector of `D`. Both hold for `|β_k⟩⟨ferro|`, where |ferro⟩ is a fully polarised extreme and `|β_k⟩` the `E_min` eigenstate in the Hamming-weight-k sector. Uniform dephasing acts on that operator as the **scalar** `−2γk`, giving `λ = −2γk + i·J·N/2` exactly (the sign follows from |β_k⟩ sitting at E_min and |ferro⟩ at E_max, so [H, |β_k⟩⟨ferro|] = −ΔE_max·|β_k⟩⟨ferro|; the conjugate mode |ferro⟩⟨β_k| carries −i). For the star, `E_min` occupies every rung k = 1..N−1 (its multiplet has `S_tot = (N−2)/2`, whose `S_z` values span exactly those rungs), which is why there are `4(N−1)` such eigenvalues counted with multiplicity: 8, 12, 16, 20 at N = 3, 4, 5, 6. The distinct values number 2(N−1), each doubly degenerate.

The choice of the polarised extreme is not cosmetic. Products between the two extremal multiplets generally fail to be eigenoperators at all: at N = 3, 4, 5 the gate finds most of them failing, with worst residuals 0.94, 1.15, 1.20. What makes the family above work is that a fully polarised state is a Z-eigenstate, so uniform dephasing cannot mix it.

**On SU(2).** [DEGENERACY_HUNT.md](DEGENERACY_HUNT.md) records that Z-dephasing breaks SU(2): the jump operators Z_l do not commute with S², so L is not SU(2)-symmetric. Nothing above needs it to be. The Casimir argument uses the SU(2) symmetry of **H alone**, which is intact on any graph, and the bound and attainment arguments use only Hermiticity of the jumps and the Z-eigenbasis. Broken SU(2) for L would matter if the modes were being classified by total spin; they are not.

---

## Scope: what the law actually needs

Three hypotheses are load-bearing. Each one, dropped, breaks something, and the gate holds all three fences.

**Uniform γ is required.** With site-dependent γ_l the dissipator becomes diag(−2γ_l) on the weight-1 block, stops commuting with the one-magnon Hamiltonian, and the saturation fails:

| N | γ profile | max\|Im\| | J·N/2 |
|---|---|---:|---:|
| 3 | (2.0, 0.5, 0.5) | 1.032100016 | 1.5 |
| 4 | (2.0, 0.5, 0.5, 0.5) | 1.592689809 | 2.0 |
| 4 | (0.1, 0.2, 0.3, 0.4) | 1.981222484 | 2.0 |
| 5 | (0.11, 0.29, 0.47, 0.65, 0.83) | 2.440907393 | 2.5 |

Every sampled profile falls strictly below. The hub is the most expensive site to detune, and that comparison has to be made at equal perturbation to mean anything: raising one γ by δ = 0.05 from a uniform 0.5 gives 1.999063 (hub) against 1.999527 (leaf) at N=4, and at δ = 0.2, 1.985041 against 1.992691; the ordering repeats at N=5. What is measured is that these sampled profiles fall below, not a proof that no profile can attain the bound.

The *bound* survives non-uniform γ, since the numerical-range argument needs only D self-adjoint. The *saturation* does not. Note the scope of this fence: it constrains statement (2), the Liouvillian equality. Statements (1) and (3), the Casimir gap and the minimality, are facts about H with no γ in them and are untouched.

**The polarised extreme is required.** The realising modes need the fully polarised state to be an extreme of H. Under the same uniform dephasing, models without that property do not saturate: the XY chain gives max|Im| = 1.4398736 against a spread of 2.2360680 (N=4) and 1.8812144 against 2.7320508 (N=5); the transverse-field Ising chain gives 3.7344319 against 4.1883993 (N=4). So the exact saturation belongs to isotropic Heisenberg, not to "dephasing plus any H".

**Hermitian jump operators are required, or the bound itself fails.** Take a single qubit, H = 0 (so max H gap = 0), and the jump operator c = I + iY. The generator is trace preserving, unital, has all Re λ ≤ 0, and has no Hamiltonian part at all, so it "adds only real decay" in every informal sense. Its spectrum is {0, 0, −2+2i, −2−2i}: max|Im| = 2 against a bound of 0. What is missing is self-adjointness, ‖D − D†‖ = 5.6569.

---

## Cross-topology comparison

Maximum H eigenvalue gap in the spin normalisation H = J·Σ S_i·S_j, with Im/σ read at J = 2γ, so σ = N·J/2. Both the Im/σ columns and the closed forms in the second column are checked in the gate.

| Topology | Max H gap | Im/σ at N = 4, 5, 6 | Equals J·N/2? |
|---|---|---|---|
| Star (hub + N-1 leaves) | J·N/2 (Casimir) | 1.000, 1.000, 1.000 | ✓ exactly, every N |
| Chain (open) | above J·N/2 from N = 4 | 1.183, 1.171, 1.248 | No, larger |
| Ring (closed) | c_N·J·N, c_4 = 3/4, c_6 = (5+√13)/12 | 1.500, 1.247, 1.434 | No, larger |
| Complete K_N | J·N(N+2)/8 (even N), J·(N−1)(N+3)/8 (odd N) | 1.500, 1.600, 2.000 | No, larger |
| Disconnected | sum over components | varies | varies |

At N = 3 all rows coincide at Im/σ = 1: there are only two distinct **connected** graphs on three vertices, the star is the path, and the ring is K₃.

The ring coefficient c_N approaches ln 2 ≈ 0.693147 in two branches, even N from above and decreasing (0.750000, 0.717129, 0.706387, 0.701545, 0.698949 at N = 4..12) and odd N from below and increasing (0.500000, 0.623607, 0.657883, 0.671922, 0.678994 at N = 3..11), so the c_5 = 0.6236 behind the table's own 1.247 sits below the limit while c_4 and c_6 sit above it.

The disconnected row composes by **addition**, since H = H₁ ⊕ H₂ makes the spread additive. The K₄ ⊔ P₄ instance at N=8 has spread 3.0000000000 + 2.3660254038 = 5.3660254038, which divided by σ = 4 gives the Im/σ = 1.3415063509 recorded in the anchor table above (shown there rounded to 1.342).

### The star is the minimum, and uniquely so

Exhaustive search over every connected labelled graph on N vertices, comparing ΔE_max:

| N | connected graphs searched | minimum ΔE_max | J·N/2 | minimisers |
|---|---:|---:|---:|---|
| 4 | 38 | 2.0000000000 | 2.0 | the 4 stars, nothing else |
| 5 | 728 | 2.5000000000 | 2.5 | the 5 stars, nothing else |
| 6 | 26704 | 3.0000000000 | 3.0 | the 6 stars, nothing else |

At N ≤ 6 the star is the unique connected graph attaining the minimal Hamiltonian spread, and that minimum is exactly J·N/2. This is a searched result, not a proved one, and the search grows too fast to continue by brute force past N = 6.

---

## The cavity reading

An interpretive layer, and worth being explicit about how much of it carries content.

In the April cavity picture the star is the point-focus member of the family: all N−1 bonds converge on the hub, as rays converge on a focal point, and the resulting Hamiltonian spread is the smallest available. That much is a fair description of the graph, and the minimiser search is what stands behind it.

Two things the cavity picture does not deliver here. The first is any conversion story: Im_max is exactly independent of γ while σ = N·γ is exactly independent of J, so the two coincide only on the line J = 2γ, and at γ = 100, J = 1, N = 4 the "illumination" σ = 400 stands against an oscillation of 2. Nothing is being converted; the two numbers share units and a convention. The second is the numerical aperture, which in the parent analysis is computed from the chain eigenvalue exports and is defined for the chain only; it does not carry over to the star, and no star NA has been computed.

What the dictionary does carry for the star is the weight-sector-as-transverse-plane correspondence, which is topology-independent, and the point-focus reading of the hub.

---

## Open questions

1. **Does the star remain the unique minimiser at every N?** Exhaustive at N = 4, 5, 6. A proof would presumably run through the Casimir structure: the star is the only connected graph whose Heisenberg Hamiltonian factors as J·S_A·S_B with |A| = 1.

2. **Is the fully polarised extreme also necessary for the universal saturation?** It is sufficient, and that direction is derived. XY and transverse-field Ising fail it and do not saturate, but that is two data points, not a converse. The clean question is whether some H without a polarised extreme can still attain `ΔE_max`.

3. **How far past isotropic Heisenberg does the saturation extend?** The attainment argument needs an H that commutes with Σ Z_l and has a fully polarised extreme. XXZ satisfies both and should inherit; the boundary is unmapped.

4. **Non-uniform γ.** The saturation breaks under the site-dependence sampled here, with the hub the most expensive site. The size of the deviation as a function of the γ profile is unstudied, and it is the quantity a hardware realisation would actually see. Whether *any* non-uniform profile can still attain J·N/2 is open.

5. **Connection to the K₃ weight-1 anomaly** (PROOF_WEIGHT1_DEGENERACY appendix 2026-05-17). K₃ = ring = star at N=3 also shows a +2 SWAP-invariant excess at weight-1 (S₃ standard 2-dim irrep). Whether that shares a source with the N=3 coincidence here is unexamined.

---

## Reproduction

- The gate: `python simulations/star_saturation_gate.py` (runs in a few minutes, the N=6 exhaustive search dominating).
- Python anchors N=3..6 chain/ring/star: `python simulations/f1_topology_heisenberg_small_n_anchor.py`; outputs `star_N{3,4,5,6}_python.json`.
- C# N=8 star: `dotnet test --filter "FullyQualifiedName~F1GeneralTopologyN8BlockSpectrumTests.Star"` (SLOW_N8 trait, opt-in). That test asserts the F1 palindromic-pairing identity; the `MaxImag` value used here is one of the statistics it logs, not something it asserts.
- Data: `simulations/results/f1_n8_n9_metrics/star_N{3..6}_python.json` and `star_N8.json`; `simulations/results/q_sweep_anchor/star_N{3..6}_Q*.json`. In every star JSON the field `MaxImag` equals `N · GammaValue` at the Marrakesh convention. (There is no `σ` field; `SigmaShift` is the palindrome center −2Nγ, a different quantity.)
