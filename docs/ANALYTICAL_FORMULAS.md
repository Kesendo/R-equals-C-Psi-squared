# Analytical Formulas Reference

**Status:** Living formula registry. Each formula carries its own tier label.
**Date:** March 31, 2026 (last updated May 29, 2026)
**Authors:** Thomas Wicht, Claude (Opus 4.6/4.7/4.8)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

**Purpose:** Before building a Liouvillian, check here. Every formula
below replaces a matrix computation.

**Convention:** J = coupling strength. γ = dephasing rate per qubit.
N = number of qubits. w = XY-weight (count of X/Y in a Pauli string).
Formulas in ASCII; prose uses Unicode (Ψ, Π, Σ, γ).

---

## Spectral Structure (replace eigenvalue computation)

### F1. Palindrome equation (Tier 1, proven)

    Π · L · Π⁻¹ = -L - 2Σγ · I

Every Liouvillian eigenvalue λ has a partner at -λ - 2Σγ.
Every decay rate d pairs with 2Σγ - d.

**Valid for:** Heisenberg, XY, Ising, XXZ, DM; Z-dephasing; any graph;
any N; non-uniform γ per qubit. Two Π families (P1, P4).
**Breaks for:** depolarizing noise (error = (2/3)Σγ, linear in γ and N).
**Replaces:** palindrome verification (54,118 eigenvalues at N=8).
**Source:** [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md)

**The global charge-conjugation X⊗N is Π² (F1², corollary, 2026-05-21):**
- The square of the palindrome conjugation acts on each Pauli string σ as
  (−1)^{n_Y+n_Z}·σ: Π sends Y→iZ, Z→iY per site, so Π² sends Y→−Y, Z→−Z and
  fixes I, X. That is exactly conjugation by the global X-string X⊗N = ⊗_l X_l.
  Hence Π² = X⊗N.
- X⊗N commutes with the chain XY + Z-dephasing Liouvillian L; it pairs the
  joint-popcount sectors (p_c, p_r) ↔ (N−p_c, N−p_r), so paired sectors share
  spectra. This is the BlockSpectrum builder's sector-pairing shortcut, which
  halves the number of eigendecompositions.
- Π is order 4; F1² = Π² is its even power. The full order-4 Π is now also a
  builder shortcut: F1PalindromeOrbitPairing (wired 2026-05-22) groups the
  joint-popcount sectors into Π-orbits of 4, eigendecomposing one primary per
  orbit, subsuming the Π² = X⊗N pairing with a further factor 2.
- The repo held both halves typed but unconnected until now:
  PiOperator.SquaredEigenvalue returns Π²'s eigenvalue (−1)^{Σ bit_b} =
  (−1)^{n_Y+n_Z}; XGlobalChargeConjugationPairing carries the X⊗N action
  (−1)^{n_Y+n_Z}. The same function under two names. X⊗N was first seen
  empirically (DEGENERACY_HUNT.md, April 2026), typed independently
  2026-05-12, identified as Π² on 2026-05-21.
- The bit_a counterpart: the other Pauli Z₂ parity, Π²_X = Z⊗N (the global
  Z-string), is registered as a corollary of F61.
- Anchor: PiOperator.cs (SquaredEigenvalue), XGlobalChargeConjugationPairing.cs,
  SYMMETRY_FAMILY_INVENTORY.md #7, MIRROR_SYMMETRY_PROOF.md.

**H-block residual is γ-independent (uniform AND non-uniform γ; closed 2026-05-18):**
- The dissipator-block residual M_D = Π·L_D·Π⁻¹ + L_D + 2Σγ·I vanishes per Pauli string for arbitrary per-site γ_l, because the per-site Z-dephasing kernel is proportional to I_4 and the F1 σ-shift `2Σγ·I` cancels the sum exactly.
- Hence ‖M‖²_F = ‖M_H‖²_F = c_H · F(N, G) for every γ pattern.
- Distinct from T1/depol siblings below, whose per-site kernels are NOT proportional to I and so DO have non-trivial Σγ² and (Σγ)² structure surviving any σ-shift choice.
- Closes the earlier F1 OpenQuestion conjecturing a Σγ_l² replacement of (Σγ)² in F(N, G); the conjecture was incorrect, no formula change to F(N, G) is required.
- Anchor: [PROOF_F1_NONUNIFORM_GAMMA](proofs/PROOF_F1_NONUNIFORM_GAMMA.md); verification: [simulations/_f1_nonuniform_gamma_verify.py](../simulations/_f1_nonuniform_gamma_verify.py); typed: PalindromeResidualScalingClaim (unchanged; XML doc notes γ-independence).

**T1 amplitude-damping residual (closed form, 2026-05-18):**
- ‖M(T1)‖²_F = 4^(N−1)·[3·Σγ²+4·(Σγ)²]  (H-independent, γ_Z-independent; bit-exact N=2..5)
- Π²-decomposition (Pythagorean orthogonal):
  - ‖M_anti(T1)‖²_F = 4^(N−1)·Σγ²  (F82/F84 amplitude-damping side)
  - ‖M_sym(T1)‖²_F  = 4^(N−1)·[2·Σγ²+4·(Σγ)²]  (Π²-even complement)
- Anchor: [PROOF_F1_T1_RESIDUAL_CLOSED_FORM](proofs/PROOF_F1_T1_RESIDUAL_CLOSED_FORM.md) (Statement + Step 7); typed: F1T1ResidualClosedForm, F1T1ResidualPi2Decomposition

**Depolarizing-noise residual (closed form, 2026-05-18):**
- ‖M(depol)‖²_F = 4^(N−1)·[(16/9)·Σγ² + 16·(Σγ)²]  (H-independent, γ_Z-independent, topology-independent; bit-exact N=2..5)
- Π²-decomposition: trivial, M_l is Pauli-basis-diagonal ⟹ M_anti = 0; F5's (2/3)Σγ scalar is the complementary scalar diagnostic
- F1 σ-shift = 0 (depol's per-Pauli-string diagonal cannot be absorbed by a constant σ·I)
- Anchor: [PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM](proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md); typed: F1DepolResidualClosedForm

**General-topology universality (closed 2026-05-18):**
- The (B, D2) parameterisation of ‖M(N, G)‖²_F = c_H · F(N, G) extends bit-exactly to all connected graphs (path, cycle, star, K_N, K_{2,N−2}, random connected Erdős-Rényi), disconnected components (B and D2 sum across components), weighted edges (B → Σ_b J²_b), and the single-body class (D2/2 prefactor).
- Verification record: Python at N=5, 6 across named/random/disconnected/weighted/single-body; C# graph-aware at N=5 across chain/ring/star/disconnected; C# F1 palindromic-pairing identity at N=7 across chain/ring/star/K_4 + disjoint-3-chain via `LiouvillianBlockSpectrum.ComputeSpectrumPerBlock`.
- The substantive analytic content was already established in [PROOF_CROSS_TERM_FORMULA](proofs/PROOF_CROSS_TERM_FORMULA.md) Lemma 3 + Corollary (bond-disjointness independent of connectivity); this closure adds the disconnected + weighted-edge sections plus the verification record.
- Closes the last F1 OpenQuestion ("general topology beyond chain/ring/star/K_N"); F1 family open-question count is ZERO as of 2026-05-18 (first time empty).
- Anchor: [PROOF_F1_GENERAL_TOPOLOGY](proofs/PROOF_F1_GENERAL_TOPOLOGY.md); verification: [simulations/_f1_general_topology_verify.py](../simulations/_f1_general_topology_verify.py) + [F1GeneralTopologyN7BlockSpectrumTests](../compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs); typed: F1GeneralTopologyVerifiedClaim (Tier 2 verified).

### F2. w=1 Liouvillian dispersion relation (Tier 1, proven D10)

    omega_k = 4J * (1 - cos(pi*k/N)),    k = 1, ..., N-1

N-1 distinct frequencies for the Heisenberg chain w=1 Liouvillian sector.
Machine-precision match for 15 frequencies (N=2-6). Tight-binding model
with hopping 2J.
Three independent validations: (1) eigenvalue match < 1e-12, (2) Poisson
spacing in w=1 sector (RMT), (3) SFF modulation peak at omega_1 matches
to <1% for N=2-4, 6 ([Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)).

**Note (2026-04-20):** F2 describes the w=1 sector of the Heisenberg
LIOUVILLIAN (Pauli strings with exactly one X or Y factor), NOT the
single-excitation Hamiltonian eigenvalues. The w=1 sector is
(N-1)-dimensional (N-1 oscillatory modes), giving denominator
(N-1)+1 = N in the cosine argument. The ZZ term in the Heisenberg
Hamiltonian produces a diagonal shift in the effective hopping matrix
that is absent in the XY case. For the XY chain single-excitation
Hamiltonian spectrum, see F2b below.

**Valid for:** Heisenberg chain, open boundaries, all N (verified N=2-6).
**Replaces:** full Liouvillian diagonalization for w=1 frequencies.
O(N) instead of O(4^{3N}).
**Source:** [Analytical Spectrum](../experiments/ANALYTICAL_SPECTRUM.md),
[D10_W1_DISPERSION](proofs/derivations/D10_W1_DISPERSION.md)

### F2b. XY chain single-excitation spectrum (Tier 1, proven)

    E_k = 2J · cos(π·k / (N+1)),    k = 1, ..., N

N eigenvalues for the single-excitation sector of the XY chain
H = (J/2) · Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) with open boundary
conditions. The single-excitation Hamiltonian H_SE is an N×N
tridiagonal matrix with off-diagonals J and zero diagonal. Its
eigenvectors are the OBC sine modes:

    ψ_k(i) = √(2/(N+1)) · sin(π·k·(i+1)/(N+1)),  i = 0, ..., N−1

These eigenvalues appear as the oscillatory frequencies Im(λ)
of the |ΔN| = 1 Liouvillian coherences (|vac⟩⟨ψ_k| sector),
which are the dominant modes contributing to per-site purity dynamics
for single-excitation initial states.

The denominator is N+1 (not N as in F2) because the OBC Dirichlet
boundary conditions require ψ to vanish at the two virtual sites
just beyond the chain, i = −1 and i = N. Both are enforced by the
sine formula above: ψ_k(−1) = sin(0) = 0 and ψ_k(N) = sin(π·k) = 0.
The effective chain length is therefore N+2 with two fixed endpoints,
yielding N interior modes with wavenumber spacing π/(N+1).

**Distinction from F2:** F2 describes the w=1 LIOUVILLIAN sector for
Heisenberg, with dimension N−1 and argument π·k/N. F2b describes the
single-excitation HAMILTONIAN sector for XY, with dimension N and
argument π·k/(N+1). They describe different mathematical objects in
different Hamiltonians.

**Valid for:** XY chain (H = (J/2)(XX+YY)), open boundaries, all N.
**Verified:** N = 3, 4, 5, 6, residual < 10⁻¹⁵.
**Replaces:** numerical diagonalization of H_SE.
**Scripts:** [`eq021_obc_sine_basis.py`](../simulations/eq021_obc_sine_basis.py).
**Source:** [OBC_SINE_BASIS_FINDINGS](../review/OBC_SINE_BASIS_FINDINGS.md),
standard tight-binding theory for OBC chains.

### AT. Absorption Theorem (Tier 1, proven)

    Re(lambda_k) = -2*gamma * <n_XY>_k       (exact)

The absorption rate of any Liouvillian eigenmode equals twice the
dephasing rate times the mode's mean light content. Here ⟨n_XY⟩ is the
expectation of the X/Y Pauli factor count in the eigenvector's Pauli
decomposition: ⟨n_XY⟩ = Σ_P |c_P|² n_XY(P) / ||v||².

**Proof sketch:** L_H is anti-Hermitian (v†L_Hv purely imaginary). L_D
is diagonal in the Pauli basis with eigenvalues -2γ n_XY(P). For any
right eigenvector v: Re(λ) = v†L_D v/||v||² = -2γ⟨n_XY⟩. Full proof
in three steps.

**The absorption quantum is 2γ.** Each X/Y Pauli factor costs exactly 2γ
in absorption rate. The spectrum is a ladder with rung spacing 2γ. The
Hamiltonian smooths the ladder (⟨n_XY⟩ can be non-integer) but cannot
change the endpoints or the fundamental quantum.

**F3, F8, F33, and D6 are corollaries** of this theorem.
The palindromic sum rule (α_fast + α_slow = 2Σγ) follows from combining
this theorem with the palindromic weight swap (⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N).

**Valid for:** any Hermitian Hamiltonian, real or complex (Heisenberg, XY,
Ising, XXZ, DM, transverse/Y fields, magnetic flux — L_H is anti-Hermitian for
*every* Hermitian H, since Step 1 needs only H^T = H*); Z-dephasing; any graph;
any N; non-uniform γ_k per site (replace 2γ with 2Σ_k γ_k × \[σ_k ∈ {X,Y}\]).
**Breaks for:** non-dephasing dissipators (amplitude damping T1, depolarizing),
which add a non-diagonal part to L_D and shift the rate (see F82, F84). No
Hamiltonian, real or complex, breaks it; the genuine boundary is the dissipator
(caveat closed 2026-05-28, bit-exact against a random complex Hermitian H).
**Replaces:** eigenvalue range computation; palindromic sum rule verification;
spectral gap derivation; unpaired mode rate identification.
**Verified:** 1,342 modes, N=2-5, γ=0.01-1.0, J=0.1-5.0, CV = 0.0000.
**Source:** [Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md),
[Absorption Theorem Discovery](../experiments/ABSORPTION_THEOREM_DISCOVERY.md)

### F3. Decay rate bounds (Tier 1, corollary of Absorption Theorem)

    min rate = 2*gamma       (w=1 modes, pure sector)
    max rate = 2*(N-1)*gamma (w=N-1 modes)
    bandwidth = 2*(N-2)*gamma

**Now a corollary of the Absorption Theorem (AT):** min = 2γ because
the smallest nonzero ⟨n_XY⟩ ≈ 1 (pure weight-1 modes). Max = 2(N-1)γ
for the fastest paired modes (⟨n_XY⟩ ≈ N-1). The XOR drain at 2Nγ
(⟨n_XY⟩ = N) sits above this range.

**Caveat resolved:** At N ≥ 4, Hamiltonian mixing creates hybrid modes
with rates below 2γ (N=4: 0.98γ, N=5: 0.62γ). These are NOT exceptions:
they are mixed-sector modes with fractional ⟨n_XY⟩ < 1. The Absorption
Theorem holds exactly for these modes; the rate 2γ⟨n_XY⟩ is correct
for non-integer ⟨n_XY⟩ ([Proton Water Chain](water/PROTON_WATER_CHAIN.md)).

**Valid for:** Heisenberg chain, uniform Z-dephasing, all N.
**Replaces:** eigenvalue range computation.
**Source:** [README](../README.md),
[Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md)

### F4. Stationary mode count (Tier 1, Clebsch-Gordan decomposition)

    Stat(N) = Sum_J m(J,N) * (2J+1)^2

m(J,N) = multiplicity of total spin J in N spin-1/2 particles.
Exact for chain topology, lower bound for higher-symmetry topologies.

**Valid for:** Heisenberg Hamiltonian, Σγ = 0, all N.
**Replaces:** null-space computation of Liouvillian.
**Source:** [Cavity Modes Formula](../experiments/CAVITY_MODES_FORMULA.md)

- **Disconnected-graph extension (Tier 1 derived, promoted 2026-05-19; landed Tier 1 candidate 2026-05-18):**
  - `dim ker L_H = Π_c (|c| + 1)` over connected components c of the graph G.
  - Predicts kernel dim = N+1 for any single connected component (chain / ring / star / K_N / arbitrary connected, matches the F4 popcount-sector count); predicts 5·5 = 25 for K_4 + disjoint 4-chain at N=8 (bit-exact verified, all four N=8 SLOW_N8 topologies).
  - Anchor: [PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS](proofs/PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md) + connected-case upper-bound closure via [DEGENERACY_PALINDROME](../experiments/DEGENERACY_PALINDROME.md) Result 2 (magnetization conservation: identity + N popcount projectors exhaust the kernel of any single connected component); typed: `F4KernelDimensionByComponentsClaim` ([compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs)).
  - **Tier 1 derived chain:** lower bound `dim ker L_H ≥ Π_c (|c|+1)` from popcount projectors + tensor-sum kernel factorisation; matching upper bound from DEGENERACY_PALINDROME Result 2 (connected-case `dim ker ≤ |c|+1`); equality follows.

### F5. Depolarizing error (Tier 1, proven)

    error = gamma * 2*N/3       (= (2/3)·Σγ)

Linear in γ and N. Hamiltonian-independent. The palindrome demands a pair
sum of 2Σγ, but depolarizing noise can supply at most (4/3)Σγ (every site
carrying a decaying Pauli); the shortfall (2/3)Σγ is the error.

**Valid for:** any Hamiltonian under depolarizing noise.
**Replaces:** numerical palindrome check for depolarizing channels.
At γ ~ 0.001 (typical IBM): error < 0.1%.
**Source:** [Depolarizing Palindrome](../experiments/DEPOLARIZING_PALINDROME.md)

### F22. GHZ XOR-drain (Tier 2, verified N=2-5)

    GHZ  -> 100% weight in XOR modes (N+1 modes at rate 2*Σγ)
    W    -> 0% XOR, 100% palindromic  (N >= 3)
    Bell -> 0% XOR, 100% palindromic  (N >= 3, Hamming distance 2)

Predictor: mixed XY Pauli weight, r = 0.976 correlation with XOR
fraction (N >= 3). GHZ fragility is not vague "delicateness" but
exact projection onto the fastest-decaying spectral sector.

**Valid for:** Heisenberg, XY, Ising, XXZ, DM; Z-dephasing; N=2-5.
**Replaces:** GHZ fragility analysis; explains why W-encoding
outperforms GHZ for state transfer (0% vs 100% drain weight).
**Source:** [XOR Space](../experiments/XOR_SPACE.md)

### F23. XOR drain vanishing fraction (Tier 1, combinatorial proof)

    fraction(XOR) = (N+1) / 4^N

N=3: 6.25%. N=5: 0.59%. N=8: 0.014%. N=20: ~10^-11.
GHZ fragility is a small-N phenomenon; at macroscopic N the XOR
sector has measure zero.

**Valid for:** any N, Z-dephasing.
**Replaces:** large-N XOR mode counting; confirms the drain is
irrelevant at macroscopic scale.
**Source:** [N->infinity Palindrome](../experiments/N_INFINITY_PALINDROME.md)

### F33. N=3 exact intermediate decay rates (Tier 1, exact rational)

    rate_1 = 2*gamma       (w=1, pure single-site coherence)
    rate_2 = 8*gamma/3     (w=2 mixed, Hamiltonian superposition)
    rate_3 = 10*gamma/3    (w=2 mixed, Hamiltonian superposition)

Three distinct rates for the N=3 Heisenberg chain. The Hamiltonian
mixes w=1 and w=2 Pauli strings into supermodes with exact rational
decay rates. At N=3 these three rates plus the extremes (0 and
6*gamma) fully determine the spectrum. At N >= 4, internal rates
become topology-dependent (only the boundary rates 2*gamma and
2*(N-1)*gamma remain universal).

**Absorption Theorem interpretation:** The fractional rates correspond
to fractional ⟨n_XY⟩ from Hamiltonian mixing:

    rate_1 = 2*gamma     -> <n_XY> = 1    (pure weight-1)
    rate_2 = 8*gamma/3   -> <n_XY> = 4/3  (mix of w=1 and w=2)
    rate_3 = 10*gamma/3  -> <n_XY> = 5/3  (mix of w=1 and w=2)

The theorem α = 2γ⟨n_XY⟩ holds exactly, including for non-integer ⟨n_XY⟩.

**Valid for:** N=3 Heisenberg chain, Z-dephasing.
**Replaces:** Liouvillian diagonalization for N=3 decay rates.
Two independent information channels (frequency vs decay) are
perfectly orthogonal at N=3.
**Source:** [Signal Processing View](../experiments/SIGNAL_PROCESSING_VIEW.md)

### F50. Weight-1 degeneracy / conserved operator count (Tier 1 lower bound proven; Tier 2 verified chain N=2-7 with K_3 N=3 anomaly)

    d_real(Re = -2*gamma) = 2N    (chain at all tested N + most connected graphs)
    d_real(Re = -2*gamma) = 8     for N=3 K_3 (triangle) instead of 2N=6  [2026-05-17 finding]

Lower bound `d_real >= 2N` is rigorously proven via the SWAP-invariant
construction: there are 2N kernel operators

    T_c^{(a)} = Sum_j Sum_{S subset complement(j), |S|=c} sigma_a^{(j)} x Z_S x I_rest

for a in {X, Y} and c = 0, 1, ..., N-1. Each T_c^{(a)} commutes with H
because the Heisenberg Hamiltonian is a sum of SWAPs, and SWAP preserves
both the active Pauli type (X or Y) and the Z-count c. The 2N operators
are linearly independent (disjoint Pauli string support).

Special cases:
- T_0^{(X)} = 2*S_x, T_0^{(Y)} = 2*S_y (global SU(2) generators)
- T_{N-1}^{(a)} = Sum_j sigma_a^{(j)} x Z_{all others} (Jordan-Wigner-type)

**Upper bound status (2026-05-17 correction):** the original proof's
Step 5 derivation conflated matrix-commutator `[SWAP_b, v]` with
conjugation action `SWAP_b · v · SWAP_b^dagger - v`. These are NOT
equivalent: `[H, v] = 0` with `H = J Sum_b (2 SWAP_b - I)` reduces to
`Sum_b [SWAP_b, v] = 0` (matrix-commutator sum), not to
`Sum_b (SWAP_b · v · SWAP_b^dagger - v) = 0`. The triangle-inequality
argument requires the latter form; the proof's chain of equivalences
has a gap. The empirical claim `= 2N` is verified for chain N=2..7 and
for ring/star/complete at N >= 4, but **N=3 K_3 (= ring = triangle =
complete on 3 vertices) gives 2N+2 = 8** purely-real Liouvillian
eigenvalues at -2gamma, not 6.

**The 2 K_3 anomaly extras** are weight-1 operators that commute with
H_K_3 globally but not with each individual bond's commutator. They
correspond to the 2-dim standard irrep of S_3 acting on the c=1 weight-1
sector (mixed X-and-Y combinations with complex phases like
omega = e^{2*pi*i/3}). Adding any external bond to the triangle (paw,
bowtie, book at N=4..5) breaks the K_3-specific S_3 symmetry and
restores the F50 count `= 2N`.

**Empirical sweep summary (2026-05-17):**
- Chain N=2..5: count = 2N (matches F50) ✓
- Ring C_n at n >= 4: count = 2N ✓
- Star K_{1,n-1} at n >= 3: count = 2N ✓
- Complete K_n at n >= 4: count = 2N ✓
- Paw / bowtie / book at N=4, 5 (graphs containing triangles): count = 2N ✓
- **N=3 K_3 = triangle: count = 2N+2 = 8 (anomaly)** ✗

**Valid for:** chain at any N + most connected graphs at N >= 4.
**Empirically violated for:** N=3 K_3 (= ring = triangle on 3 vertices) only.
**Breaks for:** anisotropic XXZ (Delta != 1), where ZZ term mixes X/Y types.
**Caveat:** This universality (with the K_3 N=3 exception noted) is UNIQUE
to k=0 and k=1. For k >= 2, d_real(k) is topology-dependent
(Chain < Star < Ring < Complete). See [Weight-2 Kernel](../experiments/WEIGHT2_KERNEL.md).
**Resolution (2026-05-17 evening):** the K_3 N=3 "weight-1 anomaly" is not a special algebraic phenomenon. It is the **small-N manifestation of a universal "central-weight excess in high-symmetry topologies" pattern**: every connected graph with non-trivial automorphism beyond chain has centralizer excess at the central weights `w ∈ {floor(N/2), ceil(N/2)}`, palindromic by F1 Π-conjugation. K_N has the largest excess; ring, star, K_N − e all have smaller-but-non-zero central excess. F50 specifically tracks weight-1, which coincides with the central weight only at N=3. For N ≥ 4 the central weight is ≥ 2, so F50's weight-1 count remains 2N for all tested topologies (chain, ring, star, K_N, paw, bowtie, book, K_4 − e). Empirical magnitudes: K_3 N=3 (+2 at w=1 / w=2), K_4 N=4 (+23 at w=2, self-palindromic), K_5 N=5 (+40 at w=2 / w=3). [`experiments/WEIGHT2_KERNEL.md`](../experiments/WEIGHT2_KERNEL.md) (April 2026) had documented the topology dependence at weight-2 for N=4..6 weeks ago; we now understand the K_3 N=3 case as the SAME phenomenon at N=3 (where central weight = 1 = F50's tracked weight). What remains open: a closed-form formula for the excess in terms of (G, N, w), and a representation-theoretic micro-structural identification of the central-weight extras beyond the empirical magnitudes. See [PROOF_WEIGHT1_DEGENERACY § Appendix Resolution](proofs/PROOF_WEIGHT1_DEGENERACY.md) for the full sweep table and the matrix-commutator-framework derivation.
**Partial closed-form (2026-05-17 late evening):** the central-weight excess decomposes by spin-isotypic sector as `central-w-excess(K_N) = Σ_{S < N/2} single_block(S, central_w) + multi_block_diff`. The **max-spin block** (S = N/2, dim N+1, 1 SU(2)-copy) contributes the universal palindromic pattern `(2, 4, 4, ..., 4, 2)` with sum `4N`, identical for all K_N and contributing equally to every weight: max-spin alone does NOT create central excess. The **sub-max spin blocks** concentrate their pure-weight content at central weights only, with a parity selection rule (K_6 S=2 contributes only at even w=2,4; K_6 S=0 vanishes entirely). The K_3 N=3 +2 excess is entirely a single-block phenomenon (S=1/2 block adds 2 at w=1, multi-block matches chain). The K_4 N=4 +23 excess decomposes as +27 single-block (mostly S=1 block adding 26 at w=2) minus −4 multi-block diff. A full closed-form for `single_block(S, w)` as `f(m_S, 2S+1, N, w)` is the remaining piece. See [`simulations/f50_spin_isotypic_decomposition.py`](../simulations/f50_spin_isotypic_decomposition.py) and [PROOF_WEIGHT1_DEGENERACY § Spin-isotypic decomposition](proofs/PROOF_WEIGHT1_DEGENERACY.md) for the full per-(S, w) table at K_3..K_6.
**Max-spin closed-form (2026-05-17 late evening, Tier 1 derived):** the max-spin contribution to single-block has a complete identification as **Dicke endpoint ladder rungs**: `single_block(S=N/2, w) = 2 if w ∈ {0, N}, else 4`. Explicit basis: w=0 → {|D_0⟩⟨D_0|, |D_N⟩⟨D_N|} (diagonal endpoint projectors with closed-form `|D_0⟩⟨D_0| = (1/2^N) Π_i (I + Z_i)` and `|D_N⟩⟨D_N| = (1/2^N) Π_i (I - Z_i)`); w=N → {|D_0⟩⟨D_N| ± h.c.} (full-ladder jump); 1 ≤ w ≤ N-1 → {|D_0⟩⟨D_w| ± h.c., |D_{N-w}⟩⟨D_N| ± h.c.} (two endpoint-anchored rungs). Total 4N pure-weight ops + (N-1)² multi-weight = (N+1)² operators in M(N+1). The multi-weight (N-1)² ops correspond to **middle-Dicke transitions** |D_k⟩⟨D_l| for k, l ∈ {1, ..., N-1} which intrinsically mix Pauli weights. Verified bit-exact N=2..5. Confirms structural reason for central-weight excess: max-spin is weight-uniform (no central bias), so excess MUST come from sub-max sectors. See [`simulations/f50_max_spin_closed_form.py`](../simulations/f50_max_spin_closed_form.py).
**Replaces:** eigenvector analysis at the first grid position;
numerical counting of purely-real eigenvalues (modulo the K_3 N=3 case).
**Source:** [Weight-1 Degeneracy Proof](proofs/PROOF_WEIGHT1_DEGENERACY.md)
(with 2026-05-17 K_3 N=3 anomaly + proof Step-5 gap appendix).

---

## Q-Factor and V-Effect (replace resonator analysis)

### F6. V-Effect gain (Tier 1-2, verified N=2-6)

    V(N) = 1 + cos(pi/N) = 2*cos^2(pi/(2N))

Q-factor amplification from coupling. γ-independent (cancels in ratio).
For N=5: (5+sqrt(5))/4 = 1.80902. For N→∞: V = 2 (saturation).
Under non-uniform γ: applies only to the extremal (best-Q) mode.

**Valid for:** Heisenberg chain, Z-dephasing, all N.
**Physically validated:** proton water chain N=1-5, machine-precision
match ([Proton Water Chain](water/PROTON_WATER_CHAIN.md)).
**Replaces:** paired Liouvillian diagonalization for V-Effect measurement.
**Source:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md)

### F7. Q-factor spectrum (Tier 1, corollary of D10)

    Q_k = 2*J/gamma * (1 - cos(pi*k/N))
    Q_max = 2*J/gamma * (1 + cos(pi/N))
    Q_min = 2*J/gamma * (1 - cos(pi/N))
    Q_mean = 2*J/gamma  (exactly, from sum of cos = 0)
    Q_spread = Q_max/Q_min = cot^2(pi/(2N))  (~4N^2/pi^2 for large N)

**Valid for:** Heisenberg chain, uniform Z-dephasing, w=1 sector.
**Replaces:** Q-factor computation from eigenvalues.
**Source:** [Analytical Spectrum](../experiments/ANALYTICAL_SPECTRUM.md)

### F8. 2× universal decay law (Tier 1, corollary of Absorption Theorem)

    rate(unpaired) = 2*N*gamma       (<n_XY> = N, pure {X,Y}^N)
    rate(paired mean) = N*gamma      (<n_XY>_fast + <n_XY>_slow = N)
    ratio = 2.00 exactly

**Now a corollary:** The palindromic sum rule α_fast + α_slow = 2Σγ
follows from the Absorption Theorem (α = 2γ⟨n_XY⟩) combined with the
palindromic weight swap (⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N, proven in
[Primordial Superalgebra](../experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)).
The "ratio 2" is the ratio of the full range (0 to 2Σγ) to the center
(Σγ) of a symmetric interval; it is a definition, not a separate law.

**Valid for:** any Hermitian Hamiltonian, real or complex (both parents — AT and
the F1 palindrome, DM included — hold for complex H); Z-dephasing; any graph; any N.
**Replaces:** unpaired mode rate computation; palindromic sum verification.
**Source:** [Energy Partition](../hypotheses/ENERGY_PARTITION.md),
[Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md)

---

## Sacrifice Zone (replace numerical optimization)

### F9. Sacrifice-zone formula (Tier 2, verified N=2-15)

    gamma_edge = N * gamma_base - (N-1) * epsilon
    gamma_other = epsilon

One-line formula. 360x improvement at N=5, 63.5x at N=15 vs V-shape.

**Valid for:** Heisenberg chain, Z-dephasing, SumMI objective.
**Replaces:** numerical optimization (Nelder-Mead, DE).
**Source:** [Resonant Return](../experiments/RESONANT_RETURN.md)

### F10. SumMI quadratic scaling (Tier 2, verified N=2-15)

    SumMI ~ 0.0053 * N^2 + 0.028 * N - 0.062

Under the sacrifice-zone formula. Quadratic, not exponential.

**Valid for:** Heisenberg chain, sacrifice-zone profile, |+>^N initial state.
**Replaces:** time evolution simulation for SumMI estimation.
**Source:** [Signal Analysis Scaling](../experiments/SIGNAL_ANALYSIS_SCALING.md)

### F11. Mode localization profile, N=5 (Tier 2, geometric)

    slowest modes (edge sacrifice): [0.52, 0.63, 0.70, 0.63, 0.52]
    fastest modes (edge sacrifice): [0.98, 0.87, 0.80, 0.87, 0.98]

Profile of each individual mode is identical under all noise profiles
(geometric, from standing wave patterns sin(πkj/N)). But which modes
are slowest depends on the noise profile: edge sacrifice selects
center-heavy modes, center sacrifice selects edge-heavy modes. Not
topologically protected (winding number = 0, Berry phase not quantized).
Correlation edge-weight vs rate: r = 0.994.

**Valid for:** N=5 Heisenberg chain, ALL Z-dephasing profiles.
**Replaces:** eigenvector decomposition of Liouvillian.
**Source:** [Cavity Mode Localization](../experiments/CAVITY_MODE_LOCALIZATION.md),
[Topological Edge Modes](../experiments/TOPOLOGICAL_EDGE_MODES.md)

---

## CΨ Crossing (replace trajectory computation)

### F12. Single-qubit universal crossing fraction (Tier 2)

    t*/T2 = 0.858367
    from: x^3 + x = 1/2,  x = e^{-t/T2}

Platform-independent. Bell states: ~10x entanglement penalty.

**Valid for:** single qubit, maximal superposition, pure dephasing.
**Replaces:** CΨ(t) trajectory simulation for crossing time.
**Source:** [Universal Quantum Lifetime](../experiments/UNIVERSAL_QUANTUM_LIFETIME.md)

### F13. r* threshold (Tier 2-3, 24,073 records)

    r* = T2 / (2*T1) = 0.2128

Separates crossers from non-crossers. Precision 0.000014. Zero false
positives across 133 qubits, 181 days.

**Valid for:** single qubit, amplitude damping + dephasing, T2echo basis.
**Replaces:** CΨ(t) simulation for crossing prediction per qubit.
**Source:** [IBM Hardware Synthesis](../experiments/IBM_HARDWARE_SYNTHESIS.md)

### F14. K-invariance (Tier 2, Lindblad scaling)

    K = gamma * t_cross = constant per bridge type

K_concurrence = 0.03596. K_MI = 0.033. K_correlation = 0.072.
Standard Lindblad time-rescaling (τ = γt), not deep physics.

**Interpretation:** In the [optical cavity analogy](../experiments/OPTICAL_CAVITY_ANALYSIS.md),
γ plays the role of external illumination and t is the system's experienced
duration. K = γ × t is the invariant decoherence dose: more light means
shorter experience, less light means longer, the product unchanged.
Structural parallel to c × τ in relativity (observation, not identification).

**Connection to Absorption Theorem:** γ is the absorption quantum rate
(2γ per X/Y Pauli factor). K = γ × t is the total absorbed dose. The
invariance of K means the total dose is state-dependent but
parameter-independent. See [K-Dosimetry](../experiments/K_DOSIMETRY.md).

**Valid for:** any Lindblad system, any bridge metric.
**Replaces:** multi-γ parameter sweeps for crossing time.
**Source:** [Crossing Taxonomy](../experiments/CROSSING_TAXONOMY.md)

### F15. θ compass (Tier 2)

    theta = arctan(sqrt(4*C*Psi - 1))

Angular distance from CΨ = 1/4 boundary. θ = 0 at crossing.

**Valid for:** any system where CΨ is defined, CΨ > 1/4.
**Replaces:** nothing directly, but provides geometric intuition.
**Source:** [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md)

### F24. Generalized crossing equation (Tier 1, algebraic + hardware-validated)

    C   = 1 - b^r + b^{2r}/2 + b^2/2
    Psi = b
    Crossing:  [1 - b^r + b^{2r}/2 + b^2/2] * b = 1/4

    b = e^{-t/T2*},   r = T2*/T1

Extends F12 (pure dephasing, r -> 0) to finite T1.
Polynomial approximation (max error < 0.001):

    t*(r)/T2* ~ 0.858 + 0.012*r + 0.375*r^2 - 0.019*r^3 - 0.084*r^4

| r = T2*/T1 | t*/T2* | Regime           |
|------------|--------|------------------|
| r -> 0     | 0.858  | Pure dephasing   |
| r = 0.5    | 0.950  | T1 = 2*T2        |
| r = 1      | 1.141  | T1 = T2          |

Hardware validated: IBM Torino qubit 52 (MAE = 0.053 with fitted T2*).
The 1/4 crossing emerges from a global fit (1/4 was not a fit target).

**Valid for:** single qubit |+> under combined T1 + T2* decay.
**Replaces:** F12 when T1 is finite; numerical CΨ(t)
simulation for superconducting qubits.
**Source:** [IBM Quantum Tomography](../experiments/IBM_QUANTUM_TOMOGRAPHY.md)

### F25. CΨ closed form, Bell+ Z-dephasing (Tier 1, proven)

    CPsi(t) = f * (1 + f^2) / 6,       f = e^{-4*gamma*t}

    dCPsi/dt = -2*gamma*f*(1 + 3*f^2) / 3

Crossing at f* = 0.8612 (from f*(1 + f*^2) = 3/2).
K = gamma * t_cross = 0.0374.

**Valid for:** Bell+ initial state, Z-dephasing, 2 qubits.
**Replaces:** numerical integration for CΨ(t) trajectory.
O(1) evaluation instead of ODE solver.
**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### F26. CΨ closed form, general Pauli channels (Tier 1, proven)

    CPsi = u * (1 + u^2 + v^2 + w^2) / 12

    u = e^{-alpha*t},   v = e^{-beta*t},   w = e^{-delta*t}
    alpha = 4*(gamma_y + gamma_z)
    beta  = 4*(gamma_x + gamma_z)
    delta = 4*(gamma_x + gamma_y)

Derivative proven < 0 for any nonzero noise (all coefficients in
the bracket are non-negative; product with -u is strictly negative).

**Valid for:** Bell+ initial state, any combination of
(gamma_x, gamma_y, gamma_z), 2 qubits.
**Replaces:** Lindblad master equation solver for multi-axis
noise on Bell states. O(1) instead of matrix exponentiation.
**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### F27. K values per noise channel (Tier 1, from F26)

    K_Z     = 0.0374    (pure Z-dephasing)
    K_X     = 0.0867    (pure X-noise)
    K_Y     = 0.0374    (pure Y-noise)
    K_depol = 0.0440    (depolarizing, gamma/3 each axis)

Complements F14 (K per bridge metric). These are K per
noise TYPE, all measured with CΨ on Bell+ state.

**Note (correction 2026-04-29):** earlier versions of this table listed
K_Y = 0.0867. That was a typo. F26 with γ_y = γ (others = 0) gives
α = 4γ, β = 0, δ = 4γ, so u = e^{-4γt}, v = 1, w = e^{-4γt}, and CΨ
reduces to u·(1+u²)/6, *identical functional form to pure Z*. Hence
K_Y = K_Z = 0.0374. The K_Y ↔ K_X-symmetry claim was wrong; the actual
pairing is K_Y ↔ K_Z (both have one of {β, δ} = 0 with α ≠ 0). This
is consistent with Bell+'s correlation structure: Y⊗Y·|Bell+⟩ = -|Bell+⟩,
while X⊗X and Z⊗Z fix it. Verified in
`framework.CPSI_CUSP_K_PER_CHANNEL` and the smoke tests.

**Valid for:** Bell+ state, single-axis or depolarizing noise.
**Replaces:** per-channel crossing time derivation.
**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### F28. Fixed-point absorber theorem (Tier 1-2)

    CPsi(rho*) < 1/4    for all primitive CPTP maps (completely positive trace-preserving: the most general physically allowed quantum operations)

Proven analytically:
- Case A: unital maps (rho* = I/d, CΨ = 0)
- Case B: local channels (rho* = product state, CΨ < 1/4)
Verified numerically:
- Case C: 100 random primitive maps, max CΨ(rho*) = 0.138

Consequence: CΨ = 1/4 is an eventual absorber. Every initial
state with CΨ > 1/4 must eventually cross below 1/4.

**Valid for:** any primitive (unique fixed-point) quantum channel,
2 qubits.
**Replaces:** per-channel verification that the fixed point sits
below 1/4.
**Source:** [Subsystem Crossing Proof](proofs/PROOF_SUBSYSTEM_CROSSING.md)

---

## Fold and Recursion (foundational)

### F16. Fold normal form (Tier 1, proven)

    R = C * (Psi + R)^2

Equivalent to Mandelbrot: u -> u^2 + c with c = C*Psi.
Boundary at CΨ = 1/4 (discriminant of fixed-point equation).

**Source:** [Mathematical Connections](MATHEMATICAL_CONNECTIONS.md)

### F17. CΨ monotonicity (Tier 1, proven)

    dCPsi/dt < 0  for all local Markovian channels

Envelope theorem for arbitrary states. 300 random CPTP maps, 0
exceptions. CΨ is Pauli-invariant (DD cannot change it).

**Source:** [CΨ Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)

### F18. Fold threshold (Tier 2, N-independent)

    Σγ_crit / J = 0.00249 (Bell state)
    Σγ_crit / J = 0.00497 (product state)

Below: no fold, CΨ oscillates forever. Above: CΨ crosses 1/4
irreversibly. Max/min ratio across N=2-5: 1.015 (1.5% variation).

**Valid for:** Heisenberg chain, Z-dephasing, N=2-5.
**Replaces:** γ sweep to find fold onset.
**Source:** [Zero Is The Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md)

### F19. Fragile bridge asymptotic constant (Tier 2)

    gamma_crit * J_bridge -> 0.50  (strong bridge limit)

Instability is a Hopf bifurcation (eigenvalues leave the imaginary
axis as a conjugate pair, creating oscillatory instability), now
identified as Liouvillian chiral
symmetry breaking (Π forces λ ↔ −λ at Σγ = 0; eigenvalues leave
the imaginary axis at γ_crit).
Linear regime: gamma_crit = 0.19 * J_bridge.
Optimal: J_bridge ~ 2J, gamma_crit = 0.41.

**Valid for:** coupled gain-loss Heisenberg chains.
**Replaces:** stability analysis for large J_bridge.
**Source:** [Fragile Bridge](../hypotheses/FRAGILE_BRIDGE.md),
[PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

---

## Thermal (replace n_bar sweep)

### F20. Thermal V-Effect: gain decreases, diversity increases (Tier 2)

    V(N, n_bar=0) = 1 + cos(pi/N)  (exact)
    V(N, n_bar=0.5) ~ 1.44  (N=5)
    V(N, n_bar->inf) ~ 1.29  (saturates)

    Frequencies(N=5, n_bar=0) = 111
    Frequencies(N=5, n_bar=5) = 445  (4x diversity gain)

**Valid for:** Heisenberg chain, Z-dephasing + amplitude damping.
**Replaces:** thermal Liouvillian sweep.
**Source:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md)

### F21. Self-heating divergence (Tier 2)

    Fixed point: n_bar -> infinity (all 6 configs tested)

Without external cooling, the system thermalizes to maximum entropy.
No self-consistent operating point exists.

**Valid for:** closed Lindblad systems with amplitude damping.
**Replaces:** convergence search for thermal equilibrium.
**Source:** [Thermal Breaking](../experiments/THERMAL_BREAKING.md)

---

## Topology and Protocols (replace parameter sweeps)

### F29. Star-topology coupling threshold (Tier 2, N=3)

    J_SB / J_SA >= 1.466    (at gamma = 0.05)

Threshold for observer-observer (AB) crossing through shared
mediator S in star topology. Below: no AB crossing. Above: AB
crosses 1/4. Receiver noise is fatal (gamma_A > 0.2 kills the
connection); sender noise is tolerable (gamma_B <= 0.5).

**Valid for:** 3-qubit star, Heisenberg, Z-dephasing,
Bell_SA x |+>_B initial state.
**Replaces:** coupling sweep for star-topology crossing threshold.
**Source:** [Star Topology Observers](../experiments/STAR_TOPOLOGY_OBSERVERS.md)

### F30. Gamma-as-signal channel capacity (Tier 2, SVD + Shannon)

    Capacity = 15.5 bits    (at 1% measurement noise, sigma = 0.01)
    Independent channels: 5  (full rank, condition number 14.8)

Spatial dephasing profile is a readable information channel.
100% classification accuracy with 4-symbol alphabet (2 bits
empirical). GHZ is completely blind (d_min = 0); |+>^N is the
optimal receiver (phased array, not omnidirectional).

**Valid for:** N=5 Heisenberg chain, Z-dephasing, |+>^5 initial.
**Replaces:** assumption that dephasing is unstructured noise;
channel capacity analysis from scratch.
**Source:** [Gamma as Signal](../experiments/GAMMA_AS_SIGNAL.md)

### F31. Relay protocol MI bound (Tier 2, N=11)

    MI improvement = +83%    (relay + 2:1 coupling vs passive)

Six relay stages, each t_stage = K/gamma. Receiving qubits get
10x noise reduction during their reception phase. Combines three
results: K/gamma timing (F14), quiet receiver (F29),
and 2:1 impedance matching.

**Valid for:** N=11 Heisenberg chain, Z-dephasing, Bell pair initial.
**Replaces:** passive propagation baseline for long chains.
**Source:** [Relay Protocol](../experiments/RELAY_PROTOCOL.md)

### F32. Optimal protection state (Tier 2, N=3)

    Slow-mode weight: 90%
    Concurrence:      0.364
    XOR weight:       0.02%

Outperforms GHZ (0% slow-mode, 100% XOR drain), W (0% slow-mode),
Bell (7% slow-mode). Composed mainly of |010>, |000>, |100>, |001>.
Loads boundary-tier palindromic pairs (rates 0.10/0.20) that decay
slowest among the dynamic modes.

**Valid for:** N=3 Heisenberg chain, Z-dephasing.
**Replaces:** trial-and-error state selection for dephasing survival.
**Source:** [Error Correction Palindrome](../experiments/ERROR_CORRECTION_PALINDROME.md)

---

## Structural (replace dimensionality arguments)

### F34. Qubit necessity equation (Tier 1, proven)

    d^2 - 2*d = 0    -->    d = 0 (nothing)  or  d = 2 (qubit)

Palindromic dephasing requires exactly 2 immune Pauli choices (I, Z)
and 2 decaying choices (X, Y) per site. This fixes d = 2 algebraically.
236 qutrit dissipators tested: 0/236 palindromic.

**Valid for:** any local dephasing model.
**Replaces:** dimensional search for palindrome-compatible systems.
**Source:** [Qubit Necessity](QUBIT_NECESSITY.md)

### F35. Dual-perspective lifetime ratio (Tier 2, hardware-validated)

    t_cross(Pi side) / t_cross(direct) ~ T1 / T2

IBM Torino qubit 52: CΨ_A crosses 1/4 at ~140 us, CΨ_B (Π
perspective) at ~895 us. Factor ~6x. The palindromic partner
decays at the T1 rate, not T2.

**Valid for:** single qubit under T1 + T2 decay, both CΨ
perspectives computed from the same density matrix.
**Replaces:** dual-perspective CΨ simulation.
**Source:** [Both Sides Visible](BOTH_SIDES_VISIBLE.md)

---

## Neural Analog (replace neural symmetry analysis)

### F36. Neural palindrome condition (Tier 1-2, proven + verified)

    Q * J * Q + J + 2*S = 0

    Q = E-I neuron swap operator
    J = Jacobian of Wilson-Cowan dynamics
    S = (1/tau_E + 1/tau_I) / 2 * I

Exact structural analog of quantum palindrome (Pi * L * Pi^-1 =
-L - 2*Σγ * I). Derived algebraically from quantum proof
via E-I swap mapping. C. elegans connectome: residual 0.013 vs
random 0.108 (8x more palindromic than chance).

**Valid for:** Wilson-Cowan neural networks with Dale's Law.
**Replaces:** ad-hoc neural symmetry analysis; connectome
palindromic quality assessment.
**Source:** [Algebraic Palindrome Neural](neural/ALGEBRAIC_PALINDROME_NEURAL.md)

### F37. Neural eigenvalue pairing (Tier 1, from F36)

    mu_k + mu_k' = -(1/tau_E + 1/tau_I)

Analog of lambda + lambda' = -2*Σγ. Every neural mode
pairs with a partner; decay rates sum to the E-I time constant sum.
Verified: mean sum = -0.3012, predicted = -0.300 (1.6% max deviation).

**Valid for:** linearized Wilson-Cowan networks satisfying F36.
**Replaces:** neural eigenvalue computation for pairing verification.
**Source:** [Proof Palindrome Neural](neural/proofs/PROOF_PALINDROME_NEURAL.md)

---

## Derived Relations (follow from combining formulas above)

All derivations verified numerically (N=2-5) against Liouvillian
eigenvalues. Full proofs: [docs/proofs/derivations/](proofs/derivations/).

### D1. Bandwidth and mode density (from F2) [VERIFIED]

    BW = omega_{N-1} - omega_1 = 8*J * cos(pi/N) --> 8*J

    rho(omega) = N / (pi * sqrt(omega * (8*J - omega)))

Van Hove singularities (peaks in the density of states where the
dispersion curve is flat) at band edges. Exact 1D tight-binding
density of states. Max frequency error < 5e-9.

**Valid for:** Heisenberg chain, w=1 sector, large N.
**Replaces:** numerical mode density estimation.

### D2. V-Effect = Q_max / Q_mean (from F6 + F7) [VERIFIED]

    V(N) = Q_max / Q_mean = (1 + cos(pi/N)) / 1

Q_mean = 2*J/gamma exactly. Proof: Sum_{k=1}^{N-1} cos(pi*k/N) = 0
(geometric series identity). Deviation < 3e-15.

### D3. Crossing time ratios (from F27) [VERIFIED]

    t_X / t_Z = K_X / K_Z = ln(2) / (2*ln(1/f*)) = 2.320
    t_depol / t_Z = K_depol / K_Z = 1.177

Verified by Lindblad propagation (Bell+ N=2). Ratio deviation < 3e-6.

### D4. Dimensional factor in crossing (from F12 + F25) [VERIFIED]

    Single qubit (d=2):    f*(1 + f*^2) = 1/2
    Bell+ 2-qubit (d=4):   f*(1 + f*^2) = 3/2 = (d-1) * 1/2

The crossing condition scales with Hilbert space dimension as (d-1)/2.
Exact to machine precision.

### D5. Dynamic palindromic mode count (from F4 + F22 + F23) [VERIFIED]

    Oscillating modes = 4^N - (N+1) - Stat(N)

N=3: 36.  N=5: 898.  Fraction --> 1 exponentially.
**Caveat:** exact at gamma -> 0 only. At finite gamma, the
Hamiltonian mixes weight-parity sectors (w with w +/- 2).

### D6. Spectral gap and mixing time (from Absorption Theorem) [VERIFIED]

    Spectral gap = 2*gamma    (minimum non-zero decay rate)
    Mixing time  <= N*ln(4) / (2*gamma)

**Immediate from Absorption Theorem:** the smallest nonzero ⟨n_XY⟩ for
any eigenmode is bounded below by the weight-1 contribution. For modes
dominated by single-site coherences: ⟨n_XY⟩ → 1, giving gap = 2γ × 1.
The spectral gap is the cost of one X/Y Pauli factor: one absorption
quantum.

Deviation < 1e-14 for all N tested.

### D7. Q-factor distribution (from F2 + F7) [VERIFIED]

    rho(Q) = 1 / (pi * sqrt((Q - Q_min) * (Q_max - Q)))

    Q_min = 2*J/gamma * (1 - cos(pi/N))
    Q_max = 2*J/gamma * (1 + cos(pi/N))
    Variance = (Q_max - Q_min)^2 / 8

Arcsine distribution (U-shaped). Modes cluster at band edges,
not at the mean. Converges with N (variance rel. error < 2% at N=50).

**Valid for:** Heisenberg chain, w=1 sector, large N.
**Replaces:** numerical Q-factor histogram.
**Source:** [D07 Derivation](proofs/derivations/D07_Q_DISTRIBUTION.md)

---

## Π Operator Properties (from PT-symmetry classification)

### F38. Π squared (Tier 1, proven + verified N=2,3)

    Pi^2 = (-1)^{w_YZ}

Diagonal parity operator in Pauli basis. w_YZ = count of Y,Z entries
in the Pauli string. Pi has order 4 (Pi^4 = I), NOT order 2.
Eigenvalues of Pi^2: +1 (half) and -1 (half), equally split.

In Hilbert space, Pi^2 is realized as conjugation by U = X^{tensor N}
(the global bit-flip): U sigma U = (-1)^{w_YZ} sigma for any Pauli string.
The two definitions agree on the per-site map (I -> I, X -> X, Y -> -Y, Z -> -Z).

**Companion result (F63):** \[L, Pi^2\] = 0 exactly for all N (proven analytically).
Pi^2 is therefore a conserved quantum number of every Liouvillian eigenmode.

**Downstream (bit_b axis foundation):** F38's (−1)^{w_YZ} eigenvalue formula on
Pauli strings is the algebraic root of the bit_b Z₂-grading of the Pauli
group, and is used as a foundational input by every derived theorem on that
axis: F88a (operator-level Klein decomposition), F108 Part 1/2/3 (Π²-even
palindrome closure via Π_5bilinear), and F112 (Lindblad Π-eigenvalue balance
under bit_b-homogeneous c). F87's trichotomy classifier reads the same bit_b
grading from the orthogonal spec(L)-palindromy axis.

**Valid for:** any N, Z-dephasing Π (P1 family).
**Replaces:** assumption that Π is involutory.
**Source:** [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md), [PROOF_BIT_B_PARITY_SYMMETRY](proofs/PROOF_BIT_B_PARITY_SYMMETRY.md)

### F39. det(Π) (Tier 1, proven + verified N=1-4)

    det(Pi) = (-1)^{N * 4^{N-1}}

    N=1: -1.  N >= 2: +1  (since 4^{N-1} is even).

**Valid for:** any N.
**Replaces:** manual determinant computation.
**Source:** [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

### F40. Fragile bridge gamma_crit at J_bridge = J (Tier 2, verified)

    gamma_crit = 0.1873  (N=2 per chain, J = J_bridge = 1.0)

Below gamma_crit: all eigenvalues on the imaginary axis (chiral phase).
Above: eigenvalue pairs leave the axis (chiral symmetry breaking, Hopf).
Petermann factor (a measure of eigenstate non-orthogonality that
diverges at exceptional points) peaks at K = 403 at gamma/gamma_crit ~ 1.46 (near-EP).

**Valid for:** N=2 per chain Heisenberg, J_bridge = J = 1.0.
**Replaces:** bisection search at this specific parameter set.
**Source:** [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md)

### F41. Palindromic time (Tier 1, corollary of D10)

    t_Pi = 2*pi / omega_min = pi / (4*J * sin^2(pi/(2*N)))

Period of the slowest palindromic modulation in the SFF. Grows as
~N^2/(pi*J) for large N. Confirmed by FFT peak matching (<1% for N=2-4, 6).

**Valid for:** Heisenberg chain, w=1 sector.
**Replaces:** numerical FFT of SFF for modulation period.
**Source:** [Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)

### F42. Timescale separation (Tier 2, verified N=2-7)

    t_Pi / t_H ~ (Delta * N^2) / (2 * pi^2 * J)  -->  0  for N -> inf

t_Pi ~ N^2 (polynomial), t_H = 2*pi/Delta ~ 4^N (exponential).
Palindromic modulation is a short-time effect; long-time behavior
is Poisson (integrable). Visibility of modulation ~1/4^N.

| N | t_Pi | t_H | t_Pi/t_H |
|---|------|------|----------|
| 3 | 3.14 | 61.3 | 0.051 |
| 5 | 8.22 | 497 | 0.017 |
| 7 | 15.9 | 5810 | 0.003 |

**Valid for:** Heisenberg chain, Z-dephasing.
**Replaces:** numerical SFF timescale extraction.
**Source:** [Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)

### F43. Sector SFF pairing (Tier 1, proven D09)

    K_freq(w, t) = K_freq(N-w, t)    (identical SFF for paired sectors)

Palindromic symmetry Pi maps w -> N-w, so sectors w and N-w have
identical spectral statistics. XOR sector (w=N): K=1.000 (all eigenvalues
degenerate at rate 2*N*gamma).

**Valid for:** Heisenberg chain, Z-dephasing, all N.
**Replaces:** sector-by-sector SFF comparison.
**Source:** [Spectral Form Factor](../experiments/SPECTRAL_FORM_FACTOR.md)

### F44. Crooks-like rate identity (Tier 1, proven D08)

    ln(d_fast / d_slow) = 2 * artanh(Delta_d / (2*Σγ))

For each palindromic pair (d_fast, d_slow) with d_fast + d_slow = 2*Σγ.
Linear approximation: ln(d_fast/d_slow) ~ Delta_d / Σγ,
giving β_eff ~ 1/Σγ (effective inverse temperature).

This is ALGEBRAIC (follows from d_fast + d_slow = 2*Σγ),
NOT a Crooks fluctuation theorem (the thermodynamic identity
relating forward and reverse process probabilities). No Jarzynski
equality holds
(<exp(-Delta_d)> ~ 0.93, not 1). The palindrome has the FORM of
detailed balance without BEING detailed balance.

**Valid for:** any palindromic Liouvillian, all N.
**Replaces:** ad-hoc thermodynamic analogies for the palindrome.
**Source:** [Entropy Production](../experiments/ENTROPY_PRODUCTION.md)

---

## Information Geometry (from Bures metric analysis)

### F45. Bures metric at the fold (Tier 2, N=2 Bell state)

    g(CΨ = 1/4) = 3.36    (Bures metric [the natural Riemannian distance between quantum states based on fidelity], finite, no singularity)

The fold at CΨ = 1/4 has no Riemannian singularity. CΨ is a
smooth coordinate everywhere along the Lindblad trajectory.

**Valid for:** N=2 Heisenberg, Bell+ initial state, Z-dephasing.
**Source:** [Information Geometry](../experiments/INFORMATION_GEOMETRY.md)

### F46. Geodesic decoherence (Tier 2, N=2 Bell state)

    Geodesic deviation = 9.1e-4    (Lindblad ~ shortest Bures path)

The Lindblad trajectory is approximately geodesic in the Bures metric.
Decoherence follows the geometrically shortest path to equilibrium.
Geometric interpretation of dCPsi/dt < 0 (proven in
[Monotonicity Proof](proofs/PROOF_MONOTONICITY_CPSI.md)).

**Valid for:** N=2 Heisenberg, Bell+ initial state, Z-dephasing.
**Source:** [Information Geometry](../experiments/INFORMATION_GEOMETRY.md)

### F47. Gaussian curvature at the fold (Tier 2, N=2)

    K(CΨ = 1/4) = -25    (negative, hyperbolic, finite)

Strong negative curvature at the fold (states diverge quickly).
Finite: no geometric singularity. Decays toward the maximally mixed
state (K → -15 at CΨ ~ 0.2).

**Valid for:** N=2 Heisenberg, Bell+ initial state, Z-dephasing.
**Source:** [Information Geometry](../experiments/INFORMATION_GEOMETRY.md)

### F48. Pythagorean decomposition (Tier 2, exact at N=2)

    L_c² = L_H² + (L_D + Σγ·I)²

where L_c = L + Σγ·I (centered Liouvillian). The cross term
{L_H, L_D + Σγ·I} vanishes exactly at N=2 because all nonzero L_H
entries connect Pauli strings with w_XY(a) + w_XY(b) = N.

The decomposition: (time evolution)² = (oscillation)² + (cooling)².

**Valid for:** N=2 Heisenberg chain, Z-dephasing, any γ. Exact.
At N≥3: cross term ~2% of ||L_c²||, γ-independent.
**Replaces:** Nothing (new structural insight, not a shortcut).
**Source:** [Primordial Qubit Algebra](../experiments/PRIMORDIAL_QUBIT_ALGEBRA.md),
[Time Irreversibility Exclusion](proofs/TIME_IRREVERSIBILITY_EXCLUSION.md)

### F49. Cross-term formula (Tier 1, proven)

    R(N) = ‖{L_H, L_Dc}‖ / (‖L_H‖ · ‖L_Dc‖) = √((N-2) / (N · 4^(N-1)))

Equivalently: R(N)² = 4(N-2) / (N · 4^N). At N=2: R = 0 (exact
Pythagorean decomposition). At N=3: R = 1/√48. At N=4: R = 1/√128.

Follows from the key identity ‖{L_H, L_Dc}‖² = 4γ²(N-2)‖L_H‖²
(bond-sum rule + spectator variance + disjoint bond supports) and
‖L_Dc‖² = γ² · 4^N · N.

γ-independent, J-independent, topology-independent. Depends only on N.

**Valid for:** Any shadow-balanced bond coupling (both Paulis in {X,Y}
or both in {I,Z}): Heisenberg XXX, XXZ, XY model, Ising, DM interaction.
Uniform Z-dephasing, any graph, all N >= 2.
**Breaks for:** Shadow-crossing couplings (X_iZ_j, Y_iZ_j).
**Replaces:** per-N numerical computation of the cross-term magnitude.
**Verified:** N=2-6, 4 topologies, 5 gamma values, 10 coupling types.
**Source:** [Proof](proofs/PROOF_CROSS_TERM_FORMULA.md),
[Experiment](../experiments/CROSS_TERM_FORMULA.md)

### F49b. Centered dissipator norm (Tier 1, proven)

    ‖L_Dc‖² = γ² · 4^N · N

where L_Dc = L_D + Nγ·I is the centered dissipator for uniform
Z-dephasing. Auxiliary lemma used in the proof of F49.

**Valid for:** Uniform Z-dephasing at rate γ per site, any N.
**Source:** [Proof](proofs/PROOF_CROSS_TERM_FORMULA.md), Lemma 1

### F49c. Cross-term formula for shadow-crossing couplings (Tier 1, proven)

    R(N) = √((N-1) / (N · 4^(N-1)))

For couplings where one bond Pauli is in {X,Y} and the other in {I,Z}
(e.g., X_iZ_j, Y_iZ_j). The bond-site variance is 1 (not 0), so
N-2 becomes N-1. Companion to F49.

**Valid for:** Shadow-crossing bond couplings, uniform Z-dephasing,
any graph, all N >= 2.
**Verified:** N=3-6, 5 coupling types, 2 topologies.
**Source:** [Proof](proofs/PROOF_CROSS_TERM_CROSSING.md),
[Experiment](../experiments/CROSS_TERM_CROSSING.md)

### F49d. Non-uniform γ extension (Tier 1, proven)

    ‖{L_H, L_Dc}‖²_F  =  4 · Σ_b ‖L_H^bond_b‖²_F · Σ_{m ∉ bond_b} γ_m²        (spectator part)
                       +     Σ_b G(bond_b, H) · (γ_{i_b} − γ_{j_b})²            (bond-asymmetry part)

with L_Dc := L_D + σ·I and σ := Σ_l γ_l. The bond-asymmetry coefficient is
G(bond_b, H) = 4 · ‖L_{ZZ-class part of H_b}^bond_b‖²_F: only the ZZ-fraction
of each bond Hamiltonian carries (γ_i − γ_j)² sensitivity, because only
ZZ-class bond transitions hit A = ε_i(α) + ε_i(β) = ±2 (XY-class transitions
all have A = 0).

Per-class G fractions (G(bond, H) / ‖L_H^bond‖²_F):

| H class                 | G / ‖L_H^bond‖² | reason                       |
|-------------------------|-----------------|------------------------------|
| Heisenberg J·(XX+YY+ZZ) | 4 / 3           | ZZ is 1/3 of the bond norm   |
| Ising J·ZZ              | 4               | ZZ is 100% of the bond norm  |
| XY J·(XX+YY)            | 0               | no ZZ content                |
| Soft Π²-odd J·(XY+YX)   | 0               | no ZZ content                |

Uniform γ_l ≡ γ recovers F49's `4γ²·(N−2)·‖L_H‖²_F` (bond-asymmetry vanishes;
spectator part collapses via the disjoint-bond-supports lemma). The convention
on ‖L_H^bond‖²_F is the full N-qubit operator-space norm (spectator I-tensors
included): 384 / 1536 / 6144 for Heisenberg J=1 at N = 3, 4, 5.

**Valid for:** any shadow-balanced bond Hamiltonian (Heisenberg, Ising, XY,
soft XY+YX, and any sum of such bond terms), any graph topology, any
non-uniform γ pattern on Z-dephasing, all N ≥ 2.
**Breaks for:** shadow-crossing couplings (F49c regime; bond-sum rule fails).
**Verified:** N = 3, 4, 5 across all four canonical H classes
([`simulations/_f49_nonuniform_gamma_crossterm_verify.py`](../simulations/_f49_nonuniform_gamma_crossterm_verify.py),
Phase 1 commit `1c6701c` + Phase 2 assertions).
**Replaces:** the F1-OpenQuestion-era conjecture that the cross-term gains
a Σγ_l² factor (closed by this extension; the structure is per-bond
spectator + per-bond asymmetry, not a single Σγ_l² term).
**Source:** [Proof](proofs/PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md),
typed claim [`F49NonUniformCrossTermClaim`](../compute/RCPsiSquared.Core/F1/F49NonUniformCrossTermClaim.cs)

---

## Cockpit and Diagnostics (replace full tomography)

### F51. Decoherence cockpit: 3-observable reduction (Tier 2, verified N=2-5, IBM-validated)

    n_eff = 3    (Purity, Concurrence, Ψ-norm)
    coverage = 88-96%    (of trajectory variance, via PCA)
    cost = 3 measurements per pair    (vs 4^N for full tomography)

PCA selects automatically which observable is PC1: Concurrence
dominates in sparse topologies (chains, stars at small N), Purity
dominates in dense topologies (rings, complete, chains at large N).
The effective dimensionality n95 grows as ~N, but the first 3 PCs
always capture 88-96%. Two practical regimes: monitoring (3 PCs,
88-96%) and full diagnostics (~N PCs, 95%).

θ is the most sensitive instrument near the ¼ boundary: 1.68×
amplification over CΨ under sacrifice-zone optimization, because
the arctan mapping amplifies small CΨ changes near ¼.

**Hardware validation:** IBM Torino Q52, CΨ = ¼ crossing predicted
at 114.7 μs, measured at 115.0 μs (0.3% error). Selective DD beats
uniform DD by 3.2× in mutual information on 5-qubit chain.
**Valid for:** Heisenberg chain, Z-dephasing and depolarizing, N=2-5,
9 topologies tested.
**Replaces:** full quantum state tomography for decoherence monitoring.
**Caveat:** Concurrence (PC1 proxy, 57% variance) never validated on
a qubit pair. Single-qubit instruments consistent, 2-qubit untested.
**Source:** [Cockpit Universality](../experiments/COCKPIT_UNIVERSALITY.md)

### F52. Thermal oscillation resilience (Tier 2, verified N=4)

    f_osc(n_bar) ~ 82%    (stable to ±2 modes out of 256, N=4)
    |Delta_osc| <= 4       (for n_bar in [0, 50])

The fraction of oscillating modes is approximately stable under
thermal load. Q_max degrades 16× (68.3 → 4.2), mode count doubles
(47 → 103 distinct frequencies), but the oscillating fraction stays
near 82%. The cavity loses sharpness, not voice.

**Mechanism (April 5, 2026):** L(n_bar) = L_0 + n_bar · L_thermal
is linear in n_bar, so eigenvalues move continuously. The oscillating
count can only change at isolated exceptional-point (EP) crossings
where a real pair splits into a complex conjugate pair (or vice
versa). At N=4: four EP crossings in \[0, 2\], each affecting exactly
2 eigenvalues. No macroscopic fraction ever changes character.

**Not invariant:** the earlier claim "82% ± 1% invariant" was too
strong. The correct statement: oscillating fraction is stable to
< 1% with isolated EP crossings, not topologically protected. At
n_bar = 50 the fraction drops to 80.5% (−4 modes).

**No phase transition.** Frequency distribution follows neither
Planck nor Stefan-Boltzmann. The mode spectrum is set by the cavity
geometry (J topology), not by thermal statistics.

**Valid for:** N=4 Heisenberg chain, Z-dephasing + amplitude damping.
**Replaces:** thermal Liouvillian sweep to check cavity survival.
The answer is: the cavity survives; no sweep needed.
**Source:** [Thermal Blackbody](../experiments/THERMAL_BLACKBODY.md),
[`thermal_ep_analysis.py`](../simulations/thermal_ep_analysis.py)

---

## Absorption Doses (replace time-evolution for dose estimates)

### F55. Universal absorption dose K_death (Tier 1, proven from D6)

    K_death = ln(10) = 2.303    (dose for 99% absorption)
    K_death / K_fold ~ 2.3      (ratio to CΨ = ¼ crossing dose)
    Immortal modes = N + 1      (zero absorption rate, all N)

**Derivation:** 99% absorption of the slowest mortal mode means
e^{−rate_min · t} = 0.01, so rate_min · t = ln(100). By formula D6:
rate_min = 2γ (spectral gap). Therefore t = ln(100)/(2γ), and
K = γ · t = ln(100)/2 = ln(10) = 2.303. Independent of N, γ,
topology. QED.

N+1 modes have exactly zero absorption rate (pure {I,Z} content,
invisible to the light). Complete absorption is impossible while the
palindrome holds. The cavity always retains light.

**Valid for:** any Heisenberg chain, Z-dephasing, all N.
**Replaces:** time evolution to find "when does the system die."
**Source:** [Trapped Light Localization](../experiments/TRAPPED_LIGHT_LOCALIZATION.md)

---

## Cusp Dynamics (replace iteration counting and trajectory integration)

### F56. Critical slowing iteration count (Tier 1, closed-form, zero fit parameters)

    K(eps, tol) = (1/2)*ln(4*eps/tol) + alpha(tol)*sqrt(eps)

    alpha(tol) = -4 + (1/2)*ln(16*tol)

K = n*sqrt(eps) is the rescaled iteration count of u_{n+1} = u^2 + c
near the cardioid cusp (c = 1/4 - eps). The leading logarithm comes
from saddle-node passage (ODE integral). The -4 comes from the
starting-transient integral (eta_0 = -1/4). The ln(16*tol) term comes
from the Modified Equation correction (Euler discretization error).

Verified: 0.5-2% accuracy over 5 tol decades (10^-8 to 10^-16) and
10 eps decades (10^-1 to 10^-10). Modified Equation slope 0.504 vs
predicted 0.500 (0.8% deviation).

**Valid for:** Mandelbrot iteration near cardioid cusp, any eps > 0,
any tol > 0. Equivalent to CΨ recursion near the 1/4 boundary.
**Replaces:** numerical iteration counting near the saddle-node.
O(1) evaluation instead of O(1/sqrt(eps)) iterations.
**Source:** [Critical Slowing at the Cusp](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md)

### F57. Trajectory dwell time at CΨ = 1/4 (Tier 1, analytical)

    t_dwell(delta) = 2*delta / |dCPsi/dt|_{t_cross}

For Bell+ Z-dephasing (using F25 derivative):

    t_dwell = 1.080088 * delta / gamma    (Bell+ specific)
    K_dwell = gamma * t_dwell = 1.080088 * delta    (gamma-independent)

The prefactor 1.080088 = 2/1.851701 is state-specific (depends on
f_cross and |dCΨ/dt| at the crossing). K_dwell is independent of γ
to machine precision (std < 2 × 10^-17 across γ = 0.1 to 10.0).

**Valid for:** any state with a CΨ = 1/4 crossing under Z-dephasing.
Prefactor is state-specific; γ-invariance of K_dwell is universal.
**Replaces:** trajectory integration for dwell-time estimation.
**Hardware verified:** ibm_kingston (Heron r2), 2026-04-16. Two Bell+ pairs with 2.55x gamma ratio (qubits 124-125, T2=\[150,310\] us; qubits 14-15, T2=\[537,381\] us). K_dwell/delta = 0.649 (pair A) and 0.694 (pair B), spread 6.3% despite 2.55x gamma difference. Gamma-invariance of K_dwell confirmed on open quantum hardware. Absolute prefactor 0.67 vs theoretical 1.08 (difference from T1 amplitude damping; the F57 formula assumes pure Z-dephasing, Kingston has T1 comparable to T2). Both CΨ(t) trajectories cross 1/4 monotonically. First two-qubit observation of the CΨ = 1/4 boundary crossing on a quantum computer; the single-qubit case was validated separately on ibm_torino Q80 at 1.9% deviation (F24, IBM Run 3).
**Data:** [data/ibm_cusp_slowing_april2026/](../data/ibm_cusp_slowing_april2026/README.md) (full JSON, PNG, and reanalysis scripts).
**Related:** [CPSI_COMPLEX_PLANE](../experiments/CPSI_COMPLEX_PLANE.md) (the saved density matrices additionally reveal a 2D-spiral structure in the complex c-plane, extending the 1D real-axis picture of BOUNDARY_NAVIGATION).
**Source:** [Critical Slowing at the Cusp](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md) (Section 6)

### F58. Weight-based dwell prefactor (Tier 2, even-weight states only)

    prefactor = (2 + 4*W2) / (1 + 6*W2)

    dCPsi/dt = -2*gamma*Psi*(1 + 6*W2)    (weight-based derivative)

W2 is the light-face ({X,Y}) sector weight at the crossing moment.
For Bell+: W2 = 0.3709 (from f_cross = 0.8612), giving prefactor =
1.080088, matching F57 exactly. The static Pauli decomposition
and the dynamic cusp passage are algebraically identical for even-weight
states. For states with odd-weight Pauli content (k = 1 terms), the
prefactor additionally requires individual coefficient magnitudes
(via sqrt(W1)), not just sector weights.

**Valid for:** states with only even-weight Pauli content (Bell+, GHZ,
basis-flip superpositions). Partial for odd-weight states (2.6% error
for |+⟩^{⊗2}).
**Generalized:** F59 extends this to any (W_0, k) pair under the same
two-sector assumption. Bell+ is the W_0 = 1/2, k = 2 special case.
**Replaces:** dynamics solution when sector weights are known at crossing.
**Source:** [Dwell Prefactor from Weights](../experiments/DWELL_PREFACTOR_FROM_WEIGHTS.md)

---

### F59. Generalized dwell prefactor (Tier 1, analytical, verified Bell+ and W_3)

    prefactor = (4/k) * (W_0 + W_k) / (W_0 + 3*W_k)

For any state whose Pauli content lives in exactly two sectors (a stationary
sector at weight W_0, immune to dephasing, and a single coherent sector at
XY-weight k > 0 with weight W_k), the dwell-time prefactor at the CPsi = 1/4
crossing is a pure algebraic function of the sector weights at the crossing
moment. Independent of N and d.

Reduces to F58 (Bell+) when W_0 = 1/2, k = 2:

    (4/2) * (1/2 + W_2) / (1/2 + 3*W_2) = (2 + 4*W_2) / (1 + 6*W_2)

The factor 3 in the denominator is structural (product rule: 2*W_k from the
dC/dt contribution plus W_k from the C*dPsi/dt contribution), NOT from d - 1
as initially conjectured in the April 5 dwell-prefactor work. The d - 1
factor cancels because both Psi and dPsi/dt carry it.

Verified on two states with different stationary weights:

| State | N | W_0 | k | Prefactor (formula) | Prefactor (direct) |
|-------|---|-----|---|--------------------:|-------------------:|
| Bell+ | 2 | 1/2 | 2 | 1.080088            | 1.080097           |
| W_3   | 3 | 1/3 | 2 | 0.876832            | 0.876839           |

Both methods agree to within 10^-5 relative error (limited by the
finite-difference resolution of the direct measurement, not by the formula).

**Valid for:** states with Pauli content in exactly two sectors, with linear
Psi-to-f relation. Includes W states, basis-flip superpositions between
computational basis states at even Hamming distance, and any state
constructed from a single off-diagonal coherence pattern.
**Breaks for:** states populating multiple coherent sectors simultaneously
(|+⟩^{⊗2} with k=1 and k=2 content); states with odd-weight Pauli content
where the linear Psi-to-f relation fails.
**Replaces:** Lindblad propagation for dwell-time estimation in the two-sector
class. O(1) algebraic evaluation given W_0 and W_k at the crossing.
**Source:** [Generalized Dwell Prefactor](../experiments/DWELL_PREFACTOR_GENERALIZED.md)

### F60. GHZ_N born below the fold (Tier 1, geometric corollary)

    CPsi(0) for GHZ_N = 1 / (2^N - 1)

For GHZ_N = (|0...0> + |1...1>) / sqrt(2), the initial value of CPsi drops
below 1/4 as soon as N >= 3, regardless of any dynamics:

| N | CPsi(0)        | Above 1/4?              |
|---|----------------|-------------------------|
| 2 | 1/3 = 0.3333   | Yes (Bell+ crosses)     |
| 3 | 1/7 = 0.1429   | No (already classical)  |
| 4 | 1/15 = 0.0667  | No                      |
| 5 | 1/31 = 0.0323  | No                      |

Derivation: C(0) = 1 (pure state), and the only nonzero off-diagonal matrix
elements are `rho[0...0, 1...1]` = 1/2 and its conjugate, giving L1 coherence
= 1 and Psi(0) = L1/(d - 1) = 1/(2^N - 1). Therefore CPsi(0) = 1/(2^N - 1).

For N >= 3, GHZ_N starts in the classical regime of the R = CPsi^2 framework
before any dephasing acts on it. This is gamma-independent: reducing the
dephasing rate cannot fix the geometric deficit. The only escape is to
change the state.

This geometric statement complements F22 (GHZ projects 100% of its coherent
weight onto the fastest-decaying XOR modes at rate 2*Sigma*gamma). The two
arguments converge: GHZ encoding is structurally unsuitable for state
transfer, regardless of how slowly it would decay.

**Valid for:** GHZ_N for any N >= 2, gamma-independent.
**Replaces:** trajectory simulation to confirm GHZ_N (N >= 3) never crosses 1/4.
**Source:** [Generalized Dwell Prefactor](../experiments/DWELL_PREFACTOR_GENERALIZED.md) Section 4; main README Section 10 Rule 1 (avoid GHZ, prefer W-type encodings).

### F61. n_XY Parity Selection Rule (Tier 1, proven, verified 64 configs N=2-6)

The Liouvillian of the isotropic Heisenberg model under Z-dephasing
preserves n_XY parity. Consequences:

1. Every eigenmode has definite n_XY parity (even or odd).
2. Every single-excitation density matrix has purely even n_XY content.
3. No SE state can excite an odd-n_XY eigenmode. The overlap is exactly zero.

**Corollary (Accessibility Boundary).** If the slowest SE-accessible mode
has rate alpha_1, and a slower mode exists at rate alpha_2 < alpha_1 with
odd n_XY parity, then alpha_1 is an exact ceiling for SE protection. No
optimizer within the SE family can reach alpha_2.

**Companion symmetry (F63):** L also commutes with Pi^2 (the w_YZ-parity
operator). Together these are the two independent Z2 symmetries of L
admitted by the Pauli algebra of d=2: bit_a (n_XY parity, this rule) and
bit_b (w_YZ parity, F63). Combined, they decompose the operator space into
4 sectors, each of dimension 4^(N-1).

**Valid for:** any isotropic Heisenberg (XX+YY+ZZ) or XY (XX+YY) model,
any graph topology, any site-dependent Z-dephasing gamma_k, any N.
**Breaks for:** amplitude damping (T1), transverse fields (odd-n_XY terms).
**Verified:** 64 configurations (N=2-6, Chain/Star/Ring/Complete, 4 gamma
profiles). Second slow mode SE Frobenius ratio < 1e-3 in all 64. Data:
`simulations/results/lens_survey/lens_survey_scaling.txt`.
**Source:** [Proof](proofs/PROOF_PARITY_SELECTION_RULE.md)

**The global Z-string Z⊗N is Π²_X (bit_a twin of F1², corollary, 2026-05-21):**
- Π²_X, the square of the X-dephasing palindrome conjugation, acts on each
  Pauli string σ as (−1)^{n_XY}·σ, the bit_a parity (n_XY counts the X and Y
  letters of σ). That is exactly conjugation by the global Z-string
  Z⊗N = ⊗_l Z_l. Hence Z⊗N = Π²_X.
- This is the bit_a counterpart of F1² (X⊗N = Π²_Z). F61 above proves L
  commutes with Π²_X; this corollary names that operator: it is the global
  Z-string Z⊗N.
- The repo held it typed three times, unconnected, until now: ZGlobalMirror
  (Z⊗N·σ·Z⊗N = (−1)^{n_XY}), F61BitAParityPi2Inheritance (Π²_X = (−1)^{n_XY}),
  and the F88a two-axis section (Π²_X = (−1)^{Σ bit_a}). The same function
  under three names; identified as the operator Z⊗N on 2026-05-21.
- Anchor: ZGlobalMirror.cs, F61BitAParityPi2Inheritance.cs,
  SYMMETRY_FAMILY_INVENTORY.md #6.

### F62. CΨ(0) for W_N (Tier 1, analytical, verified N=2-10)

    CΨ(0) = 2(N^2 - 4N + 8) / (3N^3)

For the W_N state (equal superposition of single excitations) on any
N-qubit system, the initial CΨ on any pair is given by the formula above.
The reduced density matrix for any pair (a,b) is:

    rho_ab = diag((N-2)/N, 1/N, 1/N, 0) + (1/N)|01><10| + (1/N)|10><01|

from which Tr(rho_ab^2) = (N^2 - 4N + 8)/N^2, L1 = 2/N, Psi = 2/(3N).

| N | CΨ(0) | Fraction | Above 1/4? |
|---|-------|----------|------------|
| 2 | 1/3 = 0.3333 | 8/24 | Yes (W_2 = Bell+) |
| 3 | 10/81 = 0.1235 | | No |
| 5 | 26/375 = 0.0693 | | No |
| 10 | 68/1500 = 0.0453 | | No |

**Corollary (W_N born below the fold).** For N >= 3, CΨ(0) < 1/4.
Proof: 2(N^2 - 4N + 8)/(3N^3) < 1/4 iff 3N^3 - 8N^2 + 32N - 64 > 0.
At N=3 this evaluates to 41 > 0, and the cubic is monotonically increasing
for N >= 3. Combined with the Parity Selection Rule (F61), this proves
that single-excitation states on Heisenberg chains under Z-dephasing
never cross CΨ = 1/4. The cusp exit is structurally inaccessible to them.

**Valid for:** W_N on any N-qubit system, pair-independent (permutation symmetry).
**Verified:** numerically at N=2-10, all within machine precision.
**Source:** `simulations/cpsi_wn_analytical.py`,
[Cusp-Lens Connection](../experiments/CUSP_LENS_CONNECTION.md)

### F63. w_YZ Parity Symmetry, \[L, Pi^2\] = 0 (Tier 1, proven analytically, verified N=2-5)

    [L, Pi^2_super] = 0     (exactly, for all N)

The Liouvillian commutes exactly with the Pi^2-parity superoperator
Pi^2_super(rho) = U rho U where U = X^{tensor N} (global bit-flip).
Equivalently in Pauli basis: Pi^2 acts as multiplication by (-1)^{w_YZ},
and L preserves w_YZ parity.

**Six-line proof.** U conjugation acts per site as I -> I, X -> X,
Y -> -Y, Z -> -Z. Heisenberg bonds XX/YY/ZZ each contain two Y's or two
Z's (or none), so signs cancel and U H U = H. Z-dephasing dissipator is
quadratic in Z (term Z rho Z), where the two minus signs from U Z U = -Z
also cancel, and the anti-commutator term {Z dagger Z, rho} = {I, rho}
is trivially U-invariant.

**Two independent Z2 symmetries.** Together with F61 (n_XY parity, bit_a),
this gives L two independent Z2 symmetries proven for all N. Per-site Pauli
{I, X, Y, Z} factorizes as C2 x C2 indexed by (bit_a, bit_b) =
(n_XY, w_YZ). The 4-block decomposition has dimension 4^(N-1) per block.
This is the maximal symmetry decomposition admitted by the Pauli algebra
of d=2 (no third independent Z2 classification exists per F34/QUBIT_NECESSITY).

**Per-sector mode count (closed form).** For Heisenberg coupling with Z-dephasing on the boundary qubit B, modes split within each Pi^2-sector by the Absorption Theorem applied at site B:

    conserved per sector:  even = floor(N/2) + 1,  odd = ceil(N/2)
    correlation per sector: same as conserved (palindrome symmetry)
    mirror per sector:     2^(2N-1) - 2 * (conserved per sector)

Mechanism: conserved modes are the (N+1) elementary symmetric polynomials e_d(Z_1, ..., Z_N) for d=0..N (functions of S_z, commuting with both H and Z_B). Each e_d has w_YZ-parity = d mod 2. The asymmetry for even N (one extra even-parity conserved) comes from e_N having even parity when N is even.

| N | sector | cons (e, o) | mirror (e, o) |
|---|--------|-------------|---------------|
| 2 | 8      | (2, 1)      | (4, 6)        |
| 3 | 32     | (2, 2)      | (28, 28)      |
| 4 | 128    | (3, 2)      | (122, 124)    |
| 5 | 512    | (3, 3)      | (506, 506)    |

**Valid for:** Heisenberg (XX+YY+ZZ), XY (XX+YY); Z-dephasing on any subset
of sites; any N; any graph; uniform or non-uniform gamma_k.
**Breaks for:** single-site Y or Z terms in H (transverse field, magnetic
field along Z); Y or X jump operators (no two-factor cancellation).
**Verified:** ||\[L, Pi^2\]|| = 0.000000e+00 (identically zero, not numerically
small) at N=2, 3, 4, 5. Also for Heisenberg XXX with uniform gamma at N=3.
Per-sector mode count formula verified at N=2-5; conserved-modes-as-e_d(Z) verified at N=2-4.
Data: `simulations/primordial_bit_a_bit_b_N_scaling.py`,
`simulations/mirror_mode_split_formula.py`.
**Downstream (bit_b axis derived theorems):** F108 Part 1/2/3 (Π²-even
palindrome closure under {Z, X, Y} dephasing, all three lean on [L, Π²] = 0
for the dissipator-side cancellation) and F112 (Lindblad Π-eigenvalue balance
under bit_b-homogeneous c, uses Π² = +1 on dissipator → no Π +i / −i content).
**Source:** [PROOF_BIT_B_PARITY_SYMMETRY](proofs/PROOF_BIT_B_PARITY_SYMMETRY.md),
[PRIMORDIAL_QUBIT](../hypotheses/PRIMORDIAL_QUBIT.md) Section 9

### F64. Effective gamma from cavity mode exposure (Tier 1-2, analytical + verified N=3,4)

For an N-qubit chain with XX+YY coupling and Z-dephasing only on the outermost site B at rate γ_B, the effective dephasing rate of the slowest eigenmode contributing to inner-site S coherence is:

    γ_eff = γ_B · |a_B|²        (decoherence rate, Lorentzian half-width)

equivalently in Liouvillian-eigenvalue units:

    α = 2γ_B · |a_B|²           (Liouvillian decay constant, α = -Re(λ))

where a_B is the B-site amplitude of the single-excitation Hamiltonian eigenvector. The factor of 2 between the two forms is the standard QM convention: ρ_{ij}(t) ∝ exp(-γ_eff·t) corresponds to a Liouvillian eigenvalue λ = -2γ_eff. Both express the same content; choose the convention that fits the surrounding context. This is the Absorption Theorem (F1/AT) applied to the single-excitation sector: α = 2γ_B · ⟨n_XY⟩_B, with ⟨n_XY⟩_B = |a_B|² verified to machine precision for these modes ([`factor_two_clarification.py`](../simulations/factor_two_clarification.py)).

γ_B appears as a constant prefactor. It is not diminished by intervening sites.

**Closed form at N=3** (chain S-M-B, couplings J_SM, J_MB, r = J_SM/J_MB):

                 ⎧ r² / (r² + 1)       for r < 1/√2    [zero mode]
    g(r) =       ⎨
                 ⎩ 1 / (2(r² + 1))     for r ≥ 1/√2    [bonding mode]

Crossover at r = 1/√2, g = ⅓. Special value: g(1) = ¼.

Derived from the 3×3 single-excitation Hamiltonian eigenvalues {0, ±√(J_SM² + J_MB²)} and eigenvectors.

**General N.** Diagonalize the N×N tridiagonal single-excitation Hamiltonian. Find the eigenvector with the smallest |a_B|² among those with nonzero |a_S|². No layered composition: the formula is a global eigenvector property, not a product of per-layer factors. Multiplicative stacking fails at N=4 (ratio 0.04 to 62); eigenvector formula exact (ratio 1.0000 ± 0.0003).

**Replaces:** time-domain exponential fit for γ_eff extraction.
**Valid for:** any graph topology (chain, star, ring, complete, tree), uniform or non-uniform per-bond J, XX+YY or Heisenberg single-excitation, Z-dephasing on any single site B; good-cavity regime (γ ≪ J). Breaks when γ ≥ J (bad cavity: B decoheres before transmitting).
**Topology + non-uniform J generalization (2026-04-24).** Extended from uniform-J chains to arbitrary connected graphs under either uniform or non-uniform per-bond J. When H^(1) has degenerate eigenvalues (star center-mode, ring translational eigenmodes, complete-graph symmetric modes), F64 holds after standard degenerate perturbation theory: within each H-degenerate subspace, diagonalise the site-B projector P_B to get the corrected basis; F64 then applies to the eigenvalues of P_B in that basis. Verified at N=5 and N=7 across chain, star, ring, complete, Y-tree for XY and Heisenberg; max relative error < 0.001 at γ/J = 0.01 uniform J. For random J per bond in \[0.5, 1.5\] (30 configurations across 3 trials per N), max rel err < 0.02 in 29/30 cases; the remaining case sits at 0.07 and is consistent with expected second-order PT corrections ~(γ·δJ)/J² at the non-uniform-J scale.
**Verified:** N=3 chain (max relative error 1.8% vs 64×64 Liouvillian), N=4 chain (9 configs, ratio 1.0000 ± 0.0003 vs 256×256 Liouvillian), N=5 and N=7 on chain+star+ring+complete+Y-tree uniform J (2026-04-24, via single-excitation coherence Liouvillian directly, dim N×N, max rel err < 0.001 across all (topology, B, Hamiltonian) combinations), N=5 and N=7 same topologies non-uniform J per bond in \[0.5, 1.5\] over 3 random trials (2026-04-24, max rel err 0.068 in the worst case, well inside first-order PT regime).
**Scripts:** [`primordial_gamma_analytical.py`](../simulations/primordial_gamma_analytical.py), [`primordial_gamma_stacking_4qubit.py`](../simulations/primordial_gamma_stacking_4qubit.py), [`factor_two_clarification.py`](../simulations/factor_two_clarification.py), [`f64_topology_scan.py`](../simulations/f64_topology_scan.py) (topology generalization).
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md), [F64_TOPOLOGY_GENERALIZATION](../experiments/F64_TOPOLOGY_GENERALIZATION.md)

### F65. Single-excitation spectrum of uniform open XX chain (Tier 1, proven, verified N=3..30)

For the uniform open XX chain (all couplings J, N sites) with Z-dephasing at rate γ₀ on endpoint B = site N-1, the single-excitation dissipation rates are:

    α_k / γ₀ = (4 / (N+1)) · sin²(kπ / (N+1)),    k = 1, ..., N

This is F64 evaluated on the analytically known eigenvectors ψ_k(i) = √(2/(N+1)) · sin(πk(i+1)/(N+1)) of the N×N tridiagonal single-excitation Hamiltonian. The endpoint amplitude |ψ_k(N-1)|² = (2/(N+1)) · sin²(kπ/(N+1)), and the Absorption Theorem gives α_k = 2γ₀ · |a_B|².

**Properties:**
- All α_k lie in \[0, 2γ₀\].
- Internal symmetry: α_k = α_{N+1-k}, from sin²(kπ/(N+1)) = sin²((N+1-k)π/(N+1)). This mirror is within the single-excitation spectrum; the palindromic pairing α_a + α_b = 2γ₀ of F1 in general maps single-excitation modes to multi-excitation sectors.
- Maximum rate: α_max / γ₀ = 4/(N+1) when N is odd (then k = (N+1)/2 is integer and sin² = 1 is attained exactly); when N is even the maximum lies strictly below 4/(N+1). The single-excitation sector never reaches 2γ₀ for N ≥ 2; its maximum decays as 4/(N+1) → 0 for growing N.

**Niven rationality.** All α_k/γ₀ are rational if and only if N+1 ∈ {1, 2, 3, 4, 6}, i.e., N ∈ {0, 1, 2, 3, 5}. For all other N the values are algebraic irrationals (golden-ratio family at N=4,9; √2 family at N=7; general cyclotomic otherwise). This follows from Niven's theorem applied to α_k/γ₀ = (2/(N+1))·(1 − cos(2kπ/(N+1))):
the rate is rational iff cos(2kπ/(N+1)) is rational, and for rational q the only rational
values of cos(qπ) are {0, ±1/2, ±1}. Every k clears that bar exactly when N+1 ∈ {1, 2, 3, 4, 6}
(the criterion is on cos(2rπ), so N=3 qualifies via cos(π/2)=0 even though sin(π/4) is irrational).

**Verified values:**
- N=3: α/γ₀ ∈ {1/2, 1, 1/2}
- N=4: α/γ₀ ∈ {0.276393, 0.723607, 0.723607, 0.276393} (algebraic irrationals from sin²(π/5), sin²(2π/5); golden-ratio family)
- N=5: α/γ₀ ∈ {1/6, 1/2, 2/3, 1/2, 1/6}

**Verified:** Formula matches the tridiagonal N×N single-excitation eigendecomposition to machine precision (max error 1.2 · 10⁻¹⁵) for N=3..30. All single-excitation rates confirmed present (to within O((γ₀/J)²) perturbative corrections, see below) in the full 4^N Liouvillian spectrum for N=3..7. Dynamical check at γ₀ = 0.01, where second-order shifts are ~10⁻⁶: formula predicts the decay rate of coherence operators ρ_k = |ψ_k⟩⟨0| under full Liouvillian propagation to within 10⁻⁴ relative error for all k at N=5. Asymptotic 1/(N+1)³ scaling of α_min verified; ratio to 4π²/(N+1)³ rises monotonically from 0.81 at N=3 to 0.99 at N=15.

**Perturbative nature.** The formula is derived by applying the Absorption Theorem (AT) to single-excitation coherence operators |ψ_k⟩⟨vac|, treating them as decoupled Liouvillian right eigenvectors. This is exact to first order in γ₀/J. At finite γ₀ the Lindblad dissipator mixes |ψ_k⟩⟨vac| with other sectors, and the full-Liouvillian eigenvalue shifts by O((γ₀/J)²) relative to the formula. For γ₀/J = 0.05 and N=5, the relative shift is ≈ 4·10⁻³ (verified via full eigendecomposition in `palindromic_partner_f67.py`). The palindromic pairing F1 survives this shift exactly: α_b + α_p = 2γ₀ to machine precision, even as each individual rate deviates from its first-order value (see F68).
**Scripts:** [`single_excitation_spectrum.py`](../simulations/single_excitation_spectrum.py), [`f65_dynamic_verification.py`](../simulations/f65_dynamic_verification.py)
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md)

### F66. Pole modes at α = 0 and α = 2γ₀ (Tier 1, verified XY chain with B at endpoint, N=3..7)

For the uniform open XY chain with single-site Z-dephasing at the endpoint site B = N-1, the dissipation interval \[0, 2γ₀\] has exact eigenvalues at both endpoints:

- **α = 0 modes:** ⟨n_XY⟩_B = 0 (forced by the Absorption Theorem α = 2γ₀·⟨n_XY⟩_B). Dominant Pauli strings have I or Z at every site (total XY-weight = 0). Z-basis population content, completely shielded from Z-dephasing at B.

- **α = 2γ₀ modes:** ⟨n_XY⟩_B = 1 (same mechanism). Dominant Pauli strings have X or Y at every site (total XY-weight = N). Maximally off-diagonal at B, fully exposed to Z-dephasing.

The two poles are palindromic partners under the conjugation Π, which maps total XY-weight w ↔ N-w (see F43). The single-excitation sector (F65) never reaches either pole for N ≥ 3; both poles live in the extreme XY-weight sectors (w = 0 and w = N).

**Multiplicity:** exactly N+1 at each pole, verified for N=3..7. Each α = 0 mode corresponds to one of the N+1 elementary symmetric polynomials e_d(Z_1, ..., Z_N) in F63 (commuting with both H and Z_B). The α = 2γ₀ sector has matching multiplicity by Π-symmetry.

**Scope.** Verified only for the uniform XY chain with B at the endpoint. Whether the same structure (existence of both poles, multiplicity N+1) persists for other topologies (ring, star, Y-junction) or for interior B-positions is open. Indirect evidence from the structure-points scan: at B = center of N=5 chain, α = 0 has multiplicity 64 (not 6), so the N+1 count is endpoint-specific.

**Verified:** ⟨n_XY⟩_B = 1.000000 exact for all α = 2γ₀ modes (N=3..5, from Pauli basis projection). Dominant Pauli strings have total XY-weight N for α = 2γ₀ modes and total XY-weight 0 for α = 0 modes (N=3, N=4 explicit). Multiplicity N+1 at each pole verified for N=3..7. Dynamical check of F63 conservation: all N+1 elementary symmetric polynomials e_d(Z_1,...,Z_N) drift by < 10⁻¹⁴ under Lindblad evolution for N=4 over 80 time units, while the non-symmetric control Z_0 Z_2 drifts by 3 × 10⁻². Confirms the conserved observables at the α = 0 pole are precisely the e_d, not arbitrary Z-products.
**Scripts:** [`two_gamma_pole.py`](../simulations/two_gamma_pole.py), [`f65_dynamic_verification.py`](../simulations/f65_dynamic_verification.py)
**Source:** [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) section "The dissipation interval \[0, 2γ₀\]", [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md)

### F67. Bonding-mode encoding is the optimal dephasing-protected Bell pair (Tier 1, verified N=3, N=5)

For an isolated external reference qubit R entangled with an N-site uniform XY chain C under endpoint Z-dephasing γ₀ at Q_{N-1}, the delocalized encoding

    |Ψ⟩ = (|0⟩_R |vac⟩_C + |1⟩_R |ψ_1⟩_C) / √2

with |ψ_1⟩ = √(2/(N+1)) · Σᵢ sin(π(i+1)/(N+1)) |1ᵢ⟩ (the k=1 bonding mode of F65) decays as a pure exponential from t=0:

    N(R:C)(t) = N(0) · exp(-α_1 · t),   α_1 = (4γ₀/(N+1)) sin²(π/(N+1))

because |vac⟩⟨ψ_1| is a Liouvillian right eigenvector with eigenvalue -α_1 (Absorption Theorem applied to F65, ⟨n_XY⟩_B = |a_B(ψ_1)|² = (2/(N+1)) sin²(π/(N+1))). The R half of the Bell pair is isolated, so the only dissipation channel is the chain's slowest mode.

**Comparison with localized encodings.** Bell pairs |Φ+⟩_{R,Q_j} ⊗ |0…0⟩_rest at any chain site j are superpositions of all k-modes, |1_j⟩ = Σ_k U_{jk} |ψ_k⟩, and decay multi-exponentially: an initial fast transient from high-α_k components, then an asymptotic tail at α_1 (the slowest nonzero Liouvillian eigenvalue). The tail survives with amplitude |U_{j,1}|² = (2/(N+1)) sin²(π(j+1)/(N+1)).

**Counterintuitive equivalence A ≡ C.** By palindromic symmetry of the uniform chain, |U_{0,k}|² = |U_{N-1,k}|², so inner-localized (j=0) and outer-localized (j=N-1) encodings have **identical** decay dynamics despite their very different spatial relationships to the dephased site B=N-1. Spatial distance from noise is not the protecting mechanism; spectral encoding is. The bonding mode is optimal because it has least ⟨n_XY⟩_B of all single-excitation modes.

**Scaling.** T_2 ≡ 1/α_1 = (N+1)/(4γ₀ · sin²(π/(N+1))) → (N+1)³/(4π²γ₀) for large N. Cubic improvement with chain length, no saturation regime identified.

**Verified:** Variant B (bonding-mode) α_fit/α_1 = 0.9989 (N=3) and 0.9963 (N=5), both within 10⁻³. Variant A (inner-localized) long-time tail α_fit/α_1 = 1.046 (N=3) and 1.015 (N=5), within 5% as expected for multi-exponential. Variant C (outer-localized) yields 1.023 (N=3) and 1.008 (N=5), differing from A only by fit noise (confirms palindromic equivalence). At fixed decay rate, bonding-mode preserves ~2.2× more entanglement than either localized variant at t = 0.4 · T_2 in the N=5 run.
**Scripts:** [`bell_pair_chain_protection.py`](../simulations/bell_pair_chain_protection.py)
**Source:** F65, F66, [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md), [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)

### F68. Palindromic partner of the bonding mode (Tier 1, verified N=3, 4, 5)

By F1, the bonding-mode eigenvalue -α_b of the full Liouvillian has a partner at

    α_p = 2γ₀ - α_b

exact to machine precision. For the uniform N-site XY chain with endpoint Z-dephasing at γ₀ (same setup as F67), direct eigendecomposition confirms the pairing at N = 3, 4, 5 with |α_b + α_p - 2γ₀| < 4·10⁻¹⁵.

**Structure (N ≥ 4).** V_p lives entirely in the XY-weight-(N-1) Pauli sector, Π-mirror of the bonding mode's w=1. ⟨n_XY⟩_B = α_p / (2γ₀) approaches 1 as N grows (saturating toward the F66 pole) but stays strictly below.

**Rank-1 at N ≥ 4.** V_p admits a rank-1 SVD decomposition V_p = σ₀|u⟩⟨v|. Strict SVD rank-1 (σ₁/σ₀ < 10⁻¹⁰) holds at N = 5 (ratio 9.3·10⁻¹²) and at LAPACK `zgeev` precision at N = 4 (ratio 5.1·10⁻⁸, limited by the 16-fold degeneracy of the partner eigenvalue, not by any physical second mode).

**Operational encoding.** The Bell-pair-like R-C state (|0⟩_R|u⟩ + |1⟩_R|v⟩)/√2 propagates with off-diagonal decay rate α_fit(partner) matching spectral α_p at machine precision (rel err 1.5·10⁻¹⁶ at N = 5, 0 at N = 4). Combined with a bonding-side encoding in the same propagation, the dynamical palindromic identity α_fit(bond) + α_fit(part) = 2γ₀ holds at two distinct precision levels:

- **Legacy bonding encoding** (|vac⟩⟨ψ_1|, F65 perturbative): rel err 1.6·10⁻⁶ (N = 4) and 2.8·10⁻⁷ (N = 5). The residual is the F65 O((γ₀/J)²) shift on the bonding side, reintroduced as state-preparation pollution; the partner fit is spectrally exact because V_p is used verbatim.
- **Clean bonding encoding** (SVD of full-L V_b, same construction as V_p): rel err 2.8·10⁻¹⁶ (N = 4) and 3.8·10⁻¹⁴ (N = 5). Both sides spectrally exact; residual is integrator-/eigendecomposition-floor limited. The ~10⁹ improvement at N = 4 and ~10⁷ at N = 5 verify that the legacy residual is entirely the F65 perturbative shift and no other shift is hiding.

**N = 3 is rank-2 on both sides.** Both V_b and V_p are fourfold degenerate with matching SVD structure (σ₁/σ₀ ≈ 0.98 for each, from the mult-4 degeneracy). Any rank-1 approximation of V_p gives α_fit ≈ (α_b + α_p)/2 with visibly non-exponential decay (log-fit RMS ~10⁻¹, two orders above N ≥ 4); the clean bonding encoding is skipped at N = 3 for the same structural reason. The palindromic sum still holds spectrally; no clean operational rank-1 encoding exists.

| N | α_b | α_p | α_fit(bond, legacy) | α_fit(bond, clean) | α_fit(part) | dyn sum | rel err (legacy) | rel err (clean) |
|---|-----|-----|---------------------|--------------------|-------------|---------|------------------|-----------------|
| 3 | 0.025003 | 0.074997 | 0.024969 | skipped (rank-2) | rank-2 | n/a | n/a | n/a |
| 4 | 0.013784 | 0.086216 | 0.013784 | 0.013784 | 0.086216 | 0.100000 | 1.6·10⁻⁶ | 2.8·10⁻¹⁶ |
| 5 | 0.008303 | 0.091697 | 0.008303 | 0.008303 | 0.091697 | 0.100000 | 2.8·10⁻⁷ | 3.8·10⁻¹⁴ |

**Valid for:** uniform XY chain, endpoint Z-dephasing, N ≥ 4 for the clean rank-1 operational statement. Algebraic palindromic pairing holds for all N, all graphs, all single-site-dephasing Liouvillians (from F1).

**Verified:** spectral (H1), structural (H2), and operational (H3) all confirmed at N = 3, 4, 5. Full evidence, tables, and technical notes in [PALINDROMIC_PARTNER_MODE](../experiments/PALINDROMIC_PARTNER_MODE.md).

**Scripts:** [`palindromic_partner_f67.py`](../simulations/palindromic_partner_f67.py) (H1 + H2), [`bell_pair_partner_mode.py`](../simulations/bell_pair_partner_mode.py) (H3)
**Source:** F1, F43, F65, F66, F67, [PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md), [PALINDROMIC_PARTNER_MODE](../experiments/PALINDROMIC_PARTNER_MODE.md)

### F69. GHZ+W sector mix lifts pair-CΨ(0) above the fold at N=3 (Tier 1, sextic minimal polynomial, verified)

Neither |GHZ_3⟩ (F60: pair-CΨ = 0) nor |W_3⟩ (F62: pair-CΨ = 10/81 ≈ 0.1235) crosses the CΨ = 1/4 fold at t = 0. Their symmetric superposition

    |ψ(α)⟩ = α |GHZ_3⟩ + √(1-α²) |W_3⟩,    α ∈ [0, 1]

admits a unique optimum strictly above 1/4. The pair reduction ρ_AB (all three pairs coincide by permutation symmetry) has closed-form:

    C(α)      = Tr(ρ_AB²)   = -5α⁴/18 + 2α²/9 + 5/9
    L1_off(α) = √6 α √(1-α²) + (2/3)(1-α²)
    CΨ(α)     = C(α) · L1_off(α) / 3

The stationarity condition dCΨ/dα = 0, rationalized in x = α², gives the integer-coefficient sextic

    P(x) = 2900 x⁶ - 8060 x⁵ + 4211 x⁴ + 3832 x³ - 2428 x² - 512 x + 300 = 0

which is **irreducible over ℚ** (sympy `Poly.is_irreducible` returns True; `factor_list` returns P itself as sole factor). The optimum α²_opt is therefore an algebraic number of degree exactly 6, with no radical form in nested square roots.

**Optimum (all values computed from the exact sextic root at 25-digit precision):**

| quantity | value |
|----------|-------|
| α²_opt | 0.375420720711069 |
| α_opt  | 0.612715856422101 |
| β_opt  | 0.790303283106512 |
| min pair-CΨ(0) | 0.320411541127025 |
| ratio to 1/4   | 1.281646× |
| 3-tangle τ_ABC | 0.799453 (near-GHZ limit) |
| pair concurrence C(A,B) | 0.0210 (essentially zero) |

**Scope.** Pair-CΨ = 0.3204 is the optimum of the 2-parameter slice α·GHZ_3 + β·W_3 at N = 3.

(i) **Same slice, N ≥ 4.** The GHZ+W family peaks at 0.167 (N=4), 0.146 (N=5), 0.134 (N=6), 0.125 (N=7), 0.118 (N=8), all below 1/4. GHZ-purity contribution scales as 1/(2^N − 1) and collapses too fast for the W-contribution to lift the peak back above the fold.

(ii) **Full Dicke subspace (N ∈ {3..8}).** Pair-CΨ has no non-product local maxima on the permutation-symmetric Dicke sphere at any tested N. The only non-product stationary points are Dicke basis elements |D(N,k)⟩ (max ≈ 0.123 at N = 3, ≤ 1/12 for larger N) and the GHZ+W family optimum itself, all saddles on the full sphere (escape Δpair-CΨ ≈ 0.68 on 1% c_2 perturbation at N = 3; the 10⁻⁴ saddle threshold is cleared at every tested candidate in N ∈ {3..8}). The global supremum over non-product states is 1, approached asymptotically at the |+⟩^N product manifold but not attained isolated.

(iii) **F69 on the full sphere.** F69 is stationary on the (c_0 = c_3, c_1, c_2 = 0) slice but is a saddle on the full CP^3: c_2 > 0 is an ascent direction at F69 (Δpair-CΨ = +0.011 at c_2 = 0.01), and gradient flow from such a perturbation reaches |+⟩^3 (pair-CΨ = 1, product state). The 0.3204 value is a real algebraic fact about the slice, not a universal optimum.

(iv) **Slice-stationary saddles beyond GHZ+W (added 2026-04-27).** The (ii) enumeration was incomplete in scope. Many other 2-/3-Dicke slices admit stationary points above 1/4 at every tested N. Like F69, all are saddles on the full Dicke sphere (consistent with the (ii) "no non-product local maxima" verdict, verified by perturbation in each unused Dicke direction giving non-zero linear-order shift), but they constitute additional slice-stationary points beyond the original list of "Dicke basis + GHZ+W optimum".

**Binary-Dicke slice maxima |D_i⟩+|D_j⟩** above 1/4:

| N | total pairs | above 1/4 | best slice | max cpsi |
|---|-------------|-----------|------------|----------|
| 3 | 6 | 3 | D_1+D_2 | 0.4815 (= 13/27, exceeds F69) |
| 4 | 10 | 4 | D_2+D_3 | 0.4022 |
| 5 | 15 | 5 | D_2+D_3 | 0.3720 |
| 6 | 21 | 6 | D_2+D_3 | 0.3456 |

**Triple-Dicke slice maxima |D_i⟩+|D_j⟩+|D_k⟩** above 1/4:

| N | total triples | above 1/4 | best slice | max cpsi |
|---|---------------|-----------|------------|----------|
| 3 | 4 | 4 | D_1+D_2+D_3 | 0.8011 |
| 4 | 10 | 9 | D_2+D_3+D_4 / D_1+D_2+D_3 | 0.7136 |
| 5 | 20 | 16 | D_2+D_3+D_4 | 0.6492 |
| 6 | 35 | 25 | D_2+D_3+D_4 | 0.6163 |

The central-Dicke-triple slice is privileged at every tested N (purity_A 0.83-0.94, genuinely entangled). What was actually privileged about N=3 was the irreducible-sextic *closed form* of the GHZ+W slice. The original "F69 is special" reading conflated this algebraic feature with geometric uniqueness; geometrically, F69 is one slice-saddle among many, and not even the highest at N=3.

**Asymptotic form (resolved 2026-06-04).** The N→∞ limit is cpsi_∞ ≈ 0.4312363, derived analytically: at half-filling the reduced pair matrix converges to a closed ρ_∞ (all diagonals 1/4; off-diagonals ab/2, a²/4, 1/4), and maximizing cpsi_∞ = Tr(ρ_∞²)·L1_off/3 subject to 2a²+b²=1 gives x = a²_opt ≈ 0.264715, a root of the explicit sextic (x − 2x²)·(19 − 30x − 93x²)² = (12 + 144x − 1132x² + 1488x³)². The finite-N sequence (0.8011, 0.7136, 0.6492, 0.6163, …, brute-force-confirmed N=3..13) approaches it as cpsi(N) ≈ cpsi_∞ + C/N with C ≈ 1.08 (a clean 1/N finite-size correction). There is NO elementary closed form: the value is a sextic root, and the tempting near-misses (a/b ≈ 3/4, x ≈ 9/34) are not exact (off by ~1e-5), the third "pretty constant" near-miss alongside s*=0.709 and the ring dihedral lock (both of which looked like 1/√2 and were not). Scripts: [`_eq016_central_triple_n_infinity.py`](../simulations/_eq016_central_triple_n_infinity.py) (the ρ_∞ derivation + the sextic), [`_eq016_central_triple_bf_confirm.py`](../simulations/_eq016_central_triple_bf_confirm.py) (the 1/N finite-N confirmation, validated against this sequence).

Script: [`_eq016_n4_full_landscape.py`](../simulations/_eq016_n4_full_landscape.py), [`_eq016_verify_full_sphere.py`](../simulations/_eq016_verify_full_sphere.py) (saddle-confirmation perturbation tests).

**Why F61 does not forbid this.** F61 constrains Liouvillian evolution within a fixed n_XY parity sector, not initial-state preparation that mixes excitation sectors. See [GHZ_W_SECTOR_MIX](../experiments/GHZ_W_SECTOR_MIX.md) for the preparation-vs-evolution asymmetry discussion.

**Hardware signature.** Under Kingston-grade Z-dephasing the F69 optimum crosses CΨ = 1/4 monotonically at t* ≈ 11.2 μs. A single 2-qubit tomography at t = 0 distinguishes GHZ_3 (0), W_3 (0.123), and F69 (0.320) as three separable points, no timing needed.

**Verified:** scipy bounded minimize agrees with sympy sextic root to 3.7·10⁻¹⁰ in α²_opt. 401-point grid reproduces CΨ_opt to 5·10⁻⁸. Permutation symmetry exact (spread < 10⁻¹⁵). 3-tangle and pair concurrence cross-checked in `ghz_w_optimum_n3.py`. N ∈ {4, 5, 6} failure verified on 201-point grids. Landscape scan over the full permutation-symmetric Dicke subspace at N ∈ {3..8} (2026-04-17) confirms no non-product local maxima above 1/4 exist outside the F69 GHZ+W slice; N = 3 regression recovers pair-CΨ = 0.3204 (Δ = 1.4·10⁻⁶ from sextic root). Full evidence, sextic root list, derivation of ρ_AB(α), landscape-scan saddle diagnosis, and the 3-state spherical-scan product-state pitfall in [GHZ_W_SECTOR_MIX](../experiments/GHZ_W_SECTOR_MIX.md).

**Scripts:** [`ghz_w_optimum_n3.py`](../simulations/ghz_w_optimum_n3.py), [`sector_mix_spherical_artifact.py`](../simulations/sector_mix_spherical_artifact.py) (product-state diagnostic), [`cpsi_sector_mix_optimization.py`](../simulations/cpsi_sector_mix_optimization.py) (original sweep + Kingston dynamics), [`f69_dicke_landscape.py`](../simulations/f69_dicke_landscape.py) (full Dicke-subspace scan, N ∈ {3..8})
**Source:** F60, F61, F62, [GHZ_W_SECTOR_MIX](../experiments/GHZ_W_SECTOR_MIX.md)

### F70. ΔN selection rule for site-local observables (Tier 1, proven kinematic lemma)

For any operator ρ on N qubits and any site i, the single-site partial trace annihilates sector-coherence blocks with excitation-number difference ≥ 2:

    Tr_{¬i}(ρ^(n, m)) = 0    whenever    |n − m| ≥ 2

Consequence: every site-local observable (per-site purity, per-site expectation, α_i rescaling, c_1 closure-breaking coefficient) couples only to the |ΔN| ≤ 1 content of ρ. Sector blocks with |ΔN| ≥ 2 are invisible to any measurement factoring through a single-qubit reduced state.

**Proof.** Tr_{¬i}(|x⟩⟨y|) = ⟨x_{¬i} | y_{¬i}⟩ · |x_i⟩⟨y_i|. The inner product is 1 only if x and y agree off site i, forcing |popcount(x) − popcount(y)| ≤ 1. By linearity, blocks with |n − m| ≥ 2 vanish under partial trace.

**Generalisation.** k-local partial trace annihilates |ΔN| ≥ k + 1 blocks. Pair-observables see up to |ΔN| = 2; triple-observables up to 3; global observables unrestricted.

**Valid for:** any Hamiltonian conserving excitation number, any sector-preserving dissipator, any initial state. Purely kinematic.

**Verified:** 9 |ΔN| ≥ 2 pairs tested at N=5 via coherence-block isolation (pure superposition vs classical mixture), eight non-trivial plus the trivial (0, 5); all give zero contribution to machine precision. [sector_kernel.json](../simulations/results/c1_sector_kernel/sector_kernel.json).

**Replaces:** the empirical observation "c_1 coherence contribution vanishes for |ΔN| ≥ 2" with an analytical lemma; explains the XOR_SPACE center-modes invisibility to site-local measurement; bounds the sector-kernel for PTF's α_i closure structure.

**Scripts:** [`c1_sector_kernel.py`](../simulations/c1_sector_kernel.py), [`c1_bilinearity_test.py`](../simulations/c1_bilinearity_test.py).
**Source:** [PROOF_DELTA_N_SELECTION_RULE](proofs/PROOF_DELTA_N_SELECTION_RULE.md), [PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md) Update 2026-04-20, [XOR_SPACE](../experiments/XOR_SPACE.md).

### F71. Mirror symmetry of the closure-breaking coefficient c₁ (Tier 1, proven kinematic)

For a uniform N-qubit chain with reflection-symmetric coupling and dephasing, the closure-breaking coefficient c₁ is mirror-symmetric across bonds:

    c₁(N, b, ρ₀) = c₁(N, N−2−b, ρ₀)

for all bond indices b ∈ {0, ..., N−2} and any reflection-symmetric initial state ρ₀.

**Proof sketch.** The spatial reflection R (site i ↔ site N−1−i) commutes with the uniform Liouvillian: \[L_A, R_sup\] = 0. Under R, bond b maps to bond N−2−b: R · T_b · R = T_{N−2−b}. Therefore exp(L_B(b) · t) · ρ₀ and exp(L_B(N−2−b) · t) · (R · ρ₀ · R) are related by R_sup. Per-site purity is quadratic in ρ, so any phase picked up by R on coherences (R |ψ_k⟩ = (−1)^(k+1) |ψ_k⟩) squares away. This gives P_B(b, i, t) = P_B(N−2−b, N−1−i, t), from which α_i(bond b) = α_{N−1−i}(bond N−2−b). Summing ln(α_i) over all sites and re-indexing yields c₁(b) = c₁(N−2−b).

**Consequence.** The c₁ bond profile has at most ⌈(N−1)/2⌉ independent components instead of N−1. The endpoint value c₁(0) equals c₁(N−2); if N is **even**, the center bond c₁((N−2)/2) is self-paired (its mirror image is itself) and contributes one independent component; if N is odd, there is no center bond and all N−1 bonds pair up in (N−1)/2 disjoint pairs.

**Generalisation to F86 per-bond Q_peak (Tier 1 derived 2026-05-03).** The same kinematic argument extends from c₁ in vac+SE probes to the F86 K_CC_pr per-bond observable on the (n, n+1) popcount coherence block:

    Q_peak(b)  =  Q_peak(N−2−b)        (bit-exactly, all c, N)

Proof: the F86 observable is `K_b(Q, t) = 2·Re ⟨ρ(t)| S_kernel | ∂ρ/∂J_b ⟩`. Under R, every component is invariant (uniform Z-dephasing L_D, uniform-J Hamiltonian H_xy, the Dicke probe, and the spatial-sum kernel S), while the bond-flip transforms as `R · ∂L/∂J_b · R⁻¹ = ∂L/∂J_{N−2−b}`. Hence K_b(Q, t) = K_{N−2−b}(Q, t) as functions of (Q, t), and their argmax-Q values coincide. Numerical verification: max deviation < 10⁻¹⁰ across c=2 N=5..7 and c=3 N=5..6 (`F86NewIdeasTests.F71MirrorInvariance_PerBondQPeak_BitExactSymmetricUnderBondMirror`). The per-F71-orbit substructure observed in F86 (Interior bonds not uniform within the F71-orbit grouping; central self-paired bond differs from flanking) refines the simple Endpoint/Interior dichotomy into a per-orbit classification: the F71 symmetry gives the pairing, not the value. See [PROOF_F86C_F71_MIRROR Statement 3](proofs/PROOF_F86C_F71_MIRROR.md#statement-3-f71-spatial-mirror-invariance-of-per-bond-q_peak-tier-1-derived).

**Valid for:** any Hamiltonian with \[H, R\] = 0 (uniform coupling on a symmetric graph), any dissipator with \[D, R_sup\] = 0 (uniform or R-symmetric dephasing), any initial state that is reflection-symmetric in per-site purities. Purely kinematic.
**Breaks for:** non-uniform coupling J_b ≠ J_{N−2−b}; non-uniform dephasing γ_i ≠ γ_{N−1−i}; initial states without reflection symmetry in purity. The non-uniform-J breakdown is now characterised by [F100](#f100): it is graceful, with the bond-mirror deviation D(b) = c₁(b) − c₁(N−2−b) exactly odd in the F71-anti-palindromic component of J (zero for any palindromic J, leading-order linear in the asymmetry). The non-uniform-γ breakdown is characterised by [F101](#f101): the bond-mirror deviation D(b) = c₁(b) − c₁(N−2−b) is exactly odd in the F71-anti-palindromic component of the per-site γ profile (zero for any palindromic γ, leading-order linear in the asymmetry).
**Verified:** N = 3, 4, 5, 6 for ψ_1+vac and ψ_2+vac; residuals < 10⁻⁹. Source: [`eq021_obc_sine_basis.py`](../simulations/eq021_obc_sine_basis.py), [`c1_veffect_scaling_small.py`](../simulations/c1_veffect_scaling_small.py).
**Replaces:** empirical observation of mirror-symmetric c₁ bond profiles with an analytical kinematic proof.
**Scripts:** [`eq021_obc_sine_basis.py`](../simulations/eq021_obc_sine_basis.py).
**Source:** [PROOF_C1_MIRROR_SYMMETRY](proofs/PROOF_C1_MIRROR_SYMMETRY.md), [OBC_SINE_BASIS_FINDINGS](../review/OBC_SINE_BASIS_FINDINGS.md).

### F72. Block-diagonal DD⊕CC structure of site-local purity (Tier 1, corollary of F70)

For any N-qubit chain and any initial state ρ₀, the per-site purity functional Tr(ρ_i²) decomposes into a strict block-diagonal sum over excitation-number blocks of ρ₀:

    Tr(ρ_i²) = 1/2 + P_i^DD[ρ₀^(diag)] + P_i^CC[ρ₀^(coh)]

with no cross term coupling the diagonal block (ΔN = 0) and the coherence block (|ΔN| = 1). Here ρ₀^(diag) is the projection of ρ₀ onto ΔN = 0 sector blocks and ρ₀^(coh) is the projection onto |ΔN| = 1 sector blocks.

**Proof.** The Bloch decomposition gives Tr(ρ_i²) = (1/2)(1 + ⟨X_i⟩² + ⟨Y_i⟩² + ⟨Z_i⟩²). By F70 applied to each Bloch component:

- ⟨Z_i⟩ depends linearly on the diagonal elements of ρ_i = Tr_{¬i}(ρ), which come only from ΔN = 0 blocks of ρ.
- ⟨X_i⟩ and ⟨Y_i⟩ depend linearly on off-diagonal elements of ρ_i, which come only from |ΔN| = 1 blocks.

Squaring keeps each contribution in its own sector class, so ⟨Z_i⟩² is bilinear in ρ₀^(diag) and ⟨X_i⟩² + ⟨Y_i⟩² is bilinear in ρ₀^(coh). No cross term arises.

**Consequence.** Any closure-breaking coefficient c₁ built from per-site purities decomposes, at the pre-α-fit bilinear level, into a DD-kernel (acting on diagonal block content) and a CC-kernel (acting on coherence block content) with no mixing. Finding the closed form of c₁ reduces to finding K_DD and K_CC separately.

**Generalisation.** For k-site partial traces (F70 generalisation to |ΔN| ≤ k), the Bloch-like decomposition has k+1 sub-blocks. At k = 2 (pair-site) this yields three sub-blocks DD ⊕ DC ⊕ CC, with DC a new diagonal-coherence cross specific to pair observables.

**Valid for:** any Hamiltonian conserving excitation number, any sector-preserving dissipator, any ρ₀. Purely kinematic.

**Verified:** w-scan at N = 5 with ρ₀(w) = cos(w)|vac⟩ + sin(w)|S₁⟩ under the purity-response c₁ definition confirms block-diagonal coupling at machine precision across the full w range; LSQ α-fit c₁ inherits the block structure at the pre-fit bilinear level. Pure-coherence probe gives K_CC/2 to 10⁻¹². [bilin_probe.json](../simulations/results/eq018_kernel_bilin_probe/bilin_probe.json), [kernel_extract.json](../simulations/results/eq018_kernel_extract/kernel_extract.json).

**Scripts:** [`eq018_kernel_extract.py`](../simulations/eq018_kernel_extract.py), [`eq018_kernel_bilin_probe.py`](../simulations/eq018_kernel_bilin_probe.py), [`eq018_c1_purity_response.py`](../simulations/eq018_c1_purity_response.py).
**Source:** F70, [ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) §2.3.

### F73. Spatial-sum coherence purity closure for vac-SE coherent probes (Tier 1, proven)

For any N-site qubit system with Hermitian Hamiltonian H conserving single-excitation number (\[H, N_total\] = 0) and uniform Z-dephasing γ₀, the coherent probe ρ₀^coh = (|vac⟩⟨α| + |α⟩⟨vac|) / 2 for any normalized single-excitation state |α⟩ satisfies:

    Σ_i 2 · |(ρ_coh,i)_{0,1}(t)|² = (1/2) · exp(−4 γ₀ t)

exactly, independent of the Hamiltonian's non-U(1) structure. Here (ρ_coh,i)_{0,1} is the off-diagonal element of the site-i reduced density matrix, and the sum runs over all N sites.

**Proof (general U(1) case).** Let x_i(t) = ⟨vac|ρ(t)|1_i⟩ be the amplitude of the |vac⟩⟨1_i| component of ρ(t); equivalently the (vac, SE) block of ρ as an N-vector indexed by site. Evolution under the Lindblad master equation splits into:

- **Hamiltonian part.** H preserves SE by assumption, so its restriction to SE is a Hermitian N×N matrix H_SE. The (vac, 1_i) bra-ket block evolves under H as iẋ = −H_SE x on the ket side, giving a unitary propagator U_SE(t) = exp(−i H_SE t).
- **Dephasing part.** Each D\[Z_j\] acts on the (vac, 1_i) coherence element with rate γ₀ · (⟨Z_j⟩_vac − ⟨Z_j⟩_{1_i})² / 2 = γ₀ · (2 δ_{j,i})² / 2 = 2γ₀ · δ_{j,i}. Summing over j gives a uniform 2γ₀ decay on every SE-block coherence, independent of site.

Combined: ẋ = −i H_SE x − 2γ₀ x, so x(t) = exp(−2γ₀ t) · U_SE(t) · x(0). Taking the norm: ||x(t)||² = exp(−4γ₀ t) · ||x(0)||² since U_SE is unitary. Partial-trace algebra: (ρ_coh,i)_{0,1}(t) = (1/2) · x_i(t), so Σ_i 2 · |(ρ_coh,i)_{0,1}|² = (1/2) · ||x(t)||² = (1/2) · ||x(0)||² · exp(−4γ₀ t). For the probe above, ||x(0)||² = ⟨α|α⟩ = 1. Result: (1/2) · exp(−4γ₀ t).

The argument uses only (i) \[H, N_total\] = 0 so dynamics stay in SE, (ii) H Hermitian so U_SE unitary, (iii) γ₀ uniform. No XY structure, no translation invariance, no specific shape of |α⟩ required.

**Alternative derivation (uniform XY, |α⟩ = |S₁⟩).** The original proof route via the sine basis |ψ_k⟩ of the uniform-XY single-excitation Hamiltonian: |S₁⟩ = Σ_{k odd} s_k |ψ_k⟩ with s_k = ⟨ψ_k|S₁⟩. Each single-excitation coherence |vac⟩⟨ψ_k| evolves as exp((iE_k − 2γ₀) t). Partial trace gives (ρ_coh,i)_{0,1}(t) = (1/2) · Σ_k s_k · ψ_k(i) · exp((iE_k − 2γ₀) t). Parseval on the sine basis Σ_i ψ_k(i) · ψ_{k'}(i) = δ_{k,k'} eliminates k ≠ k' cross terms; Σ_k s_k² = 1 by normalisation. Under bond-b perturbation, the sine basis and E_k shift but Parseval on any orthonormal SE basis preserves Σ_k |⟨ψ_k^B|S₁⟩|² = 1, so the sum is δJ-invariant. This derivation is XY-specific but exhibits the eigenmode structure explicitly.

**Consequence.** The spatial-sum purity functional is exactly blind to the U(1)-preserving part of the dynamics on any vac-SE coherent probe. For any closure-breaking coefficient c₁_pr built from per-site purities via the purity-response definition, bond-δJ perturbations preserve the closure value, so `K_CC[0, 1]_pr` = 0 exactly under uniform γ₀, for any H in the class.

**Scaffolding from neighbouring entries.** F70 (site-local observables see only |ΔN| ≤ 1 content) puts the (vac, SE) block in focus as the relevant coherence sector for per-site purity. F72 (DD ⊕ CC block decomposition of Tr(ρ_i²), no cross-term) isolates the CC contribution, where the (vac, SE) coherence lives. The Absorption Theorem supplies the rate 2γ₀·n_XY = 2γ₀ for SE coherences (n_XY = 1). F73 then combines these: U(1) conservation keeps the SE sector closed under H, and the spatial sum over sites collapses the unitary H-rotation to leave only the AT decay.

**Valid for:** any Hermitian H with \[H, N_total\] = 0 (XY, Heisenberg XXZ, translationally non-invariant hopping, frustrated-ladder variants, ...); uniform Z-dephasing γ₀; any normalized SE state |α⟩ admixed to |vac⟩; any N.
**Breaks for:**

- Non-uniform γ_i. The uniform 2γ₀ decay on the d_H = 1 block fails; the closure becomes K_CC ≠ 0 with mode-selective response (see [CMRR_BREAK_NONUNIFORM_GAMMA](../experiments/CMRR_BREAK_NONUNIFORM_GAMMA.md)).
- Non-U(1) Hamiltonians. \[H, N_total\] ≠ 0 breaks the SE-block closure assumption.
- Dissipators changing the d_H = 1 decay rate (mixed X/Z, amplitude damping).
- Probes with d_H > 1 admixture (e.g. (vac, S₂) with two-excitation bra-ket), where ⟨n_XY⟩ ≠ 1 and the uniform decay rate breaks.

**Verified:**

- Uniform XY baseline at N = 5, t₀ = 20: closure matches (1/2)·exp(−4·0.05·20) = 9.157819·10⁻³ to 5.67·10⁻¹⁶ deviation. `K_CC[0, 1]_pr` = 1.14·10⁻¹² (machine-precision zero), confirming δJ-invariance. [cmrr_gamma_nonuniform.json](../simulations/results/eq018_cmrr_gamma_nonuniform/cmrr_gamma_nonuniform.json).
- U(1)-class generalization at N = 5 (6 setups: XXZ at Δ ∈ {0, 0.5, 1, 2}, random Haar SE probe at Δ = 1, inhomogeneous XY with J_i ∈ \[0.5, 1.5\]): all closures within 2.22·10⁻¹⁶ to 5.83·10⁻¹⁶ (1-3 ULP of double precision) across 81 time points per setup. [f73_u1_generalization/](../simulations/results/f73_u1_generalization/), [F73_U1_GENERALIZATION](../experiments/F73_U1_GENERALIZATION.md).

**Scripts:** [`eq018_c1_purity_response.py`](../simulations/eq018_c1_purity_response.py), [`eq018_cmrr_gamma_nonuniform.py`](../simulations/eq018_cmrr_gamma_nonuniform.py) (uniform baseline), [`f73_u1_generalization_sweep.py`](../simulations/f73_u1_generalization_sweep.py) (U(1)-class sweep).
**Source:** F61, F70, F72, [ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) §2.4, [CMRR_BREAK_NONUNIFORM_GAMMA](../experiments/CMRR_BREAK_NONUNIFORM_GAMMA.md), [F73_U1_GENERALIZATION](../experiments/F73_U1_GENERALIZATION.md).
**See also:** [J_BLIND_RECEIVER_CLASSES](../experiments/J_BLIND_RECEIVER_CLASSES.md) generalises the L_D-invariant-subspace mechanism behind this closure to a three-class decomposition of J-blind initial states (DFS of L_D, H-degenerate L_D-closed block, M_α-polynomial subspace under SU(2)-Heisenberg).

### F74. Chromaticity of single-step coherence blocks (Tier 1, combinatorial)

For an N-qubit system under uniform Z-dephasing, the (n, n+1) popcount
coherence block contains exactly

    c(n, N) = min(n, N−1−n) + 1

distinct pure dephasing rates when the Hamiltonian is switched off (J = 0).
The rates are 2γ₀ · {1, 3, 5, ..., 2c−1}, corresponding to Hamming-distance
values HD ∈ {1, 3, ..., 2c−1} between popcount-n and popcount-(n+1) basis
states.

**Proof.** A basis pair (|x⟩, |y⟩) with popcount(x) = n, popcount(y) = n+1
differs at HD = 2n + 1 − 2·match sites, where match = popcount(x AND y) is
the number of sites carrying a 1 in both x and y. The constraints
match ∈ \[max(0, 2n+1−N), n\] give HD ∈ {1, 3, ..., min(2n+1, 2N−2n−1)},
hence the distinct-HD count is min(n, N−1−n) + 1. The Pauli representation
of |x⟩⟨y| has X or Y on exactly the HD sites where x and y differ, so
⟨n_XY⟩ = HD, and the Absorption Theorem gives rate 2γ₀·HD at J = 0.

**Consequences.**

- **c = 1 mono-chromatic blocks** at n = 0 and n = N−1 (single pure rate 2γ₀).
  F73 (the spatial-sum coherence closure for vac-SE probes) is this
  monochromatic case.
- **c_max at the center.** For odd N, unique at n = (N−1)/2 with
  c_max = (N+1)/2. For even N, two adjacent blocks at n = N/2−1 and
  n = N/2, both with c_max = N/2.
- At J > 0, H-coupling between different HD-channels produces dressed
  modes at intermediate rates. Q_SCALE_THREE_BANDS measures this H-mixing
  via the dressed-mode weight W(Q) and shows that the observable peak
  abs(K_CC_pr) is c-specific and N-invariant.

**Valid for:** any N-qubit system under uniform Z-dephasing. The J = 0
statement holds kinematically for any Hamiltonian. For the dynamical
interpretation at J > 0 (c as a stable sector label with H-mixing between
HD channels), H must conserve total excitation number, \[H, N_total\] = 0;
individual eigenmode rates then shift continuously with Q = J/γ₀, but c
still labels the block's mixing substructure.
**Breaks for:** non-uniform γ_i (site-dependent dephasing); the J = 0
spectrum then has rates 2·Σ_{i ∈ diff-sites} γ_i rather than the discrete
2γ₀·HD values. Non-U(1) Hamiltonians (transverse fields, odd-popcount
terms) preserve the J = 0 statement but dissolve the (n, n+1) sector's
dynamical invariance at J > 0. Non-Z dissipators (amplitude damping,
X-dephasing, depolarizing) break the diagonal Pauli action of L_D; the
|x⟩⟨y| basis pairs are no longer eigenvectors at rate 2γ₀·HD.
**Verified:** Block-structure c-values for N = 3..8 match the formula
exactly. Spectral verification at J = 0: each (n, n+1) block has exactly
c(n, N) distinct rates in {2γ₀, 6γ₀, ..., 2(2c−1)γ₀}
([Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 3).

**Replaces:** block-diagonal spectrum enumeration at J = 0; identifies
which blocks support H-mixing bands (c ≥ 2) vs which are structurally
silent (c = 1).
**Scripts:** [`q_scale_n_scaling.py`](../simulations/q_scale_n_scaling.py).
**Source:** [PROOF_CHROMATICITY](proofs/PROOF_CHROMATICITY.md),
[Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 3,
[PROOF_ABSORPTION_THEOREM](proofs/PROOF_ABSORPTION_THEOREM.md).

### F75. Mirror-pair MI for single-excitation mirror-symmetric states (Tier 1, proven algebraic)

For any pure single-excitation state on an N-site chain with mirror-symmetric amplitudes,

    |ψ⟩ = Σ_j c_j |1_j⟩,   c_{N−1−j} = η c_j,   η ∈ {+1, −1}

the mutual information between any mirror-pair sites (ℓ, N−1−ℓ) at t = 0 has the closed form

    MI(ℓ, N−1−ℓ) = 2 h(p_ℓ) − h(2 p_ℓ),    p_ℓ = |c_ℓ|²

where h(x) = −x log₂ x − (1−x) log₂(1−x) is the binary entropy.

The formula is independent of the mirror sign η (only the modulus |c_ℓ|² enters). The valid range is p_ℓ ∈ \[0, 1/2\], with MI saturating at 2 bits when p_ℓ = 1/2 (maximal mirror-pair entanglement, the pair is in a Bell state, all other site populations vanish).

**Proof.** The reduced density matrix ρ_{ℓ,N−1−ℓ} in the computational basis {\|00⟩, \|01⟩, \|10⟩, \|11⟩} is block-diagonal:

- `ρ[|00⟩⟨00|]` = Σ_{j ∉ {ℓ, N−1−ℓ}} |c_j|² = 1 − 2 p_ℓ
- `ρ[|01⟩⟨01|]` = |c_{N−1−ℓ}|² = p_ℓ
- `ρ[|10⟩⟨10|]` = |c_ℓ|² = p_ℓ
- `ρ[|11⟩⟨11|]` = 0 (single-excitation sector)
- `ρ[|10⟩⟨01|]` = c_ℓ c_{N−1−ℓ}^* = η p_ℓ
- `ρ[|01⟩⟨10|]` = η p_ℓ

The eigenvalues are {1 − 2 p_ℓ, 2 p_ℓ, 0, 0}, giving S(ρ_{ℓ,N−1−ℓ}) = h(2 p_ℓ). Both single-site marginals are diag(1 − p_ℓ, p_ℓ) with S = h(p_ℓ). The subtraction S(ρ_ℓ) + S(ρ_{N−1−ℓ}) − S(ρ_{ℓ,N−1−ℓ}) yields the formula.

**Mirror-pair sum (MM).** Summing over all mirror-pairs ℓ = 0, ..., ⌊N/2⌋ − 1:

    MM(0) = Σ_ℓ [2 h(p_ℓ) − h(2 p_ℓ)]

**Bonding-mode specialisation (F65 + F75).** For the k-th bonding mode |ψ_k⟩ = √(2/(N+1)) Σ_j sin(πk(j+1)/(N+1)) |1_j⟩ with mirror sign η = (−1)^(k+1), the site populations are

    p_ℓ(k, N) = (2/(N+1)) sin²(πk(ℓ+1)/(N+1))

and MM(0) for bonding:k is computable in O(N) operations with no propagation.

**Verified values (bonding:k on uniform chain, analytic vs simulation PeakMM at γ₀ = 0.05, uniform J = 1):**

| N | k | MM(0) analytic | PeakMM sim | ratio sim/analytic |
|---|---|----------------|------------|-------------------|
| 5 | 1 | 0.800 | 0.789 | 0.986 |
| 5 | 2 | **1.245** | 1.241 | 0.997 |
| 5 | 3 | 0.918 | 0.865 | 0.942 |
| 7 | 1 | 0.862 | 0.801 | 0.929 |
| 7 | 2 | 1.174 | 1.090 | 0.928 |
| 7 | 3 | 0.862 | 0.819 | 0.950 |
| 7 | 4 | **1.245** | 1.166 | 0.936 |
| 9 | 1 | 0.895 | 0.830 | 0.928 |
| 9 | 2 | 1.131 | 1.049 | 0.928 |
| 9 | 3 | 0.895 | 0.829 | 0.927 |
| 9 | 4 | 1.131 | 1.054 | 0.932 |
| 11 | 1 | 0.914 | 0.848 | 0.927 |
| 11 | 2 | 1.105 | 1.024 | 0.927 |
| 11 | 3 | 0.915 | 0.846 | 0.925 |
| 11 | 4 | 1.103 | 1.024 | 0.928 |
| 11 | 6 | **1.145** | 1.066 | 0.931 |
| 13 | 1 | 0.928 | 0.859 | 0.926 |
| 13 | 2 | 1.088 | 1.008 | 0.926 |
| 13 | 3 | 0.928 | 0.858 | 0.925 |
| 13 | 4 | 1.088 | 1.007 | 0.926 |
| 13 | 5 | 0.928 | 0.860 | 0.927 |
| 13 | 7 | 0.961 | 0.893 | 0.930 |

Under Heisenberg evolution on the uniform chain, bonding mode ψ_k mixes with its same-parity partner ψ_{N+1−k} via the boundary ZZ term. At N = 5 for bonding:2 the partner is ψ_4 which has identical mirror-pair populations p_ℓ (because sin²(πk(ℓ+1)/(N+1)) = sin²(π(N+1−k)(ℓ+1)/(N+1))); direct numerical propagation shows that MM(t) oscillates with period 2π/Δ (Δ = same-parity eigenvalue gap of the 2×2 block, Δ = 2√5 for N = 5) between a minimum near t = π/(2Δ) and a revival near t = π/Δ. Under uniform Z-dephasing at γ₀, the revival is damped but stays close to MM(0); the simulation-observed PeakMM matches MM(0) within 1% at N = 5 k = 2 (analytic 1.2451, numerical Lindblad max 1.2475 at t = 0.645, C# brecher PeakMM 1.2410 on a coarser grid). At larger N and different k the ratio drops to ~0.93 because the oscillation revival magnitude and dephasing decay combine less favourably.

**Why k = 2 maximises MM over k = 1, 2, 3.** Even k places a node at the (odd-N) chain center (p_{N/2} = 0 for integer k/2 when N+1 is even), so all probability mass lies on mirror-pairs: Σ_{pairs} p_ℓ = 1/2. Odd k puts mass 2/(N+1) at the center, wasting mass on the self-mirror site. The function f(p) = 2 h(p) − h(2 p) is convex on (0, 1/2), so concentrated mass distributions give larger MM; k = 2 is the smallest k that both concentrates mass on pairs and places outer-pair amplitudes at opposite signs (η = −1, first mode carrying end-to-end coherence).

**Upper bound.** For any single-excitation mirror-symmetric state, MM ≤ ⌊N/2⌋ × 2 = N bits (all pairs in pure Bell states). This bound is not achievable from single-site bonding modes; reaching it requires tensor-product pair structures like (|10⟩−|01⟩)/√2 on each mirror-pair, which is a super-single-excitation state.

**Valid for:** any pure single-excitation state with c_{N−1−j} = ±c_j on a linear N-site chain. Extends to non-linear mirror-symmetric graphs (ring, Y-junction with mirror axis) with corresponding modification of the mirror-partner indexing.
**Breaks for:** states with multi-excitation content (formula no longer applies because `ρ[|11⟩⟨11|]` ≠ 0 in general), or states without mirror amplitude symmetry (where p_ℓ ≠ p_{N−1−ℓ} gives an asymmetric 2-qubit reduced matrix).
**Verified:** Algebraic derivation confirmed against direct C# brecher propagation at N = 5, 7, 9 for k = 1, 2, 3, at N = 11 for k = 1, 2, 3, 4, 6, and at N = 13 for k = 1, 2, 3, 4, 5, 7 (the latter via matrix-free propagator); MM(0) formula matches simulation PeakMM within 7% (full decay envelope explained by 4γ₀·t dephasing + mirror-partner oscillation at t = 0.1). The sim/analytic ratio sits at **0.925 to 0.931 across all (N ≥ 7, k) tested (~25 data points)**, i.e. PeakMM = 0.93 × MM(0) with tight consistency. F75 is therefore a reliable predictor of PeakMM without any propagation.
**Scripts:** [`_check_brecher_n5_finegrid.py`](../simulations/_check_brecher_n5_finegrid.py), [`Program.cs brecher mode`](../compute/RCPsiSquared.Propagate/Program.cs), `_mm_zero_derivation.py` (table above).
**Source:** F65 (bonding-mode amplitudes), F67 (bonding as optimal decay receiver), F71 (mirror symmetry that justifies c_{N−1−j} = ±c_j), [RECEIVER_VS_GAMMA_SACRIFICE](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) (numerical context).

### F76. Pure-dephasing decay of the mirror-pair MI for bonding:k (Tier 1, proven algebraic + weak-mixing argument)

For a bonding:k initial state on a uniform Heisenberg chain with uniform Z-dephasing γ₀ and J = 1, **the time evolution of the mirror-pair MI at short times is almost entirely pure-dephasing**, with Heisenberg-mixing contributions second-order small in both γ₀·t and (V_{lk}/Δ_l)² where V_{lk} is the boundary mixing matrix element.

Under pure dephasing alone, the site-basis mirror-pair coherence c_ℓ c_{N-1-ℓ}^* decays at 4γ₀:

    ρ_{pair}(\|10⟩⟨01\|)(t) = c_ℓ c_{N-1-ℓ}^* · e^{-4γ₀ t}

while populations p_ℓ = |c_ℓ|² stay fixed. The 4×4 reduced pair density matrix at time t has eigenvalues

    {1 − 2p_ℓ,   p_ℓ(1 + λ),   p_ℓ(1 − λ),   0},   λ(t) = e^{-4γ₀ t}

giving

    MI_pair(p_ℓ, t) = 2 h(p_ℓ) − S_{ab}(p_ℓ, λ(t))
    S_{ab}(p, λ) = −(1 − 2p) log₂(1 − 2p) − p(1+λ) log₂(p(1+λ)) − p(1−λ) log₂(p(1−λ))

At t = 0: λ = 1, eigenvalues {1-2p, 2p, 0, 0}, S_{ab} = h(2p), recovers F75.
At t → ∞: λ = 0, eigenvalues {1-2p, p, p, 0}, S_{ab} = h(1-2p) + 2p (max-entropy mixture of the two single-excitation branches).

**Closed-form envelope.** Summing over mirror-pair sites ℓ ∈ \[0, ⌊N/2⌋-1\]:

    MM(t) / MM(0) = Σ_ℓ [2 h(p_ℓ) − S_{ab}(p_ℓ, λ(t))] / Σ_ℓ [2 h(p_ℓ) − h(2 p_ℓ)]

**The 0.93 envelope explained.** At γ₀ = 0.05, J = 1, and the C# brecher first-measurement grid point t = 0.1: λ = e^{-0.02} = 0.9802. Evaluated on the actual p_ℓ multisets for bonding:k (F65), the ratio sits at 0.926 to 0.937 for N = 5..13 and k = 1..7, effectively 0.93 ± 0.006 for all 25+ measured (N, k) combinations at this specific (γ₀, t). The weak p-dependence of the ratio comes from the fact that 2h(p) − S_{ab}(p, λ) is nearly proportional to 2h(p) − h(2p) for all p ∈ (0, 1/2) and λ close to 1.

**Verified values at γ₀ = 0.05, t = 0.1:**

| N | k | pure-dephasing MM/MM(0) | single-excitation-sector Lindblad (exact) | C# brecher peak/MM(0) |
|---|---|------------------------|--------------------------------------------|-----------------------|
| 5 | 2 | 0.936 | 0.936 | 0.937 (grid peak at t=0.6 → 1.002) |
| 7 | 2 | 0.932 | 0.928 | 0.928 |
| 9 | 4 | 0.934 | 0.932 | 0.932 |
| 11 | 4 | 0.933 | 0.928 | 0.928 |
| 13 | 4 | 0.928 | 0.926 | 0.926 |

Agreement to within 0.5% across all tested (N, k). Difference between "pure-dephasing" and "exact single-excitation Lindblad" columns measures the Heisenberg-mixing correction, which is always < 0.5%.

**Why Heisenberg mixing is negligible.** At t = 0+, the commutator \[H, ρ_0\] = \[H^(1), |ψ_k⟩⟨ψ_k|\] is off-diagonal in the ψ_k mode basis (the diagonal part is the unitary phase that doesn't affect MM). The off-diagonal mixing couples ψ_k to same-parity partners ψ_l via V_{lk} = (16J/(N+1)) sin(πk/(N+1)) sin(πl/(N+1)). Under mixing for small t, ψ_k "leaks" amplitude symmetrically into all ψ_l. Because MM depends only on diagonal pair populations (and specific pair coherences that dephase), and bonding modes with the same mirror-symmetry (k and N+1-k, etc.) have identical pair populations, the leakage does NOT change pair populations to first order; it only redistributes the mode occupation. Hence first-order Heisenberg-mixing has no effect on MM. Second-order (rate V²t²) is small because V/γ₀ ~ 10 but γ₀·t ~ 0.005 at the C# first-sample, making Heisenberg relative contribution (V·t)²/(4γ₀·t) ~ (V²·t)/(4γ₀) which is (16/7²) · 0.005 = 0.002 at N=7. Hence the mixing correction is ≲ 0.5% throughout the tested regime.

**Implication: the 0.93 envelope is the γ₀ signature, not a hidden constant.** The value 0.93 is a direct consequence of the choice γ₀·t = 0.005 at the first measurement. If γ₀ changes (at fixed t grid), the envelope changes: γ₀ = 0.025 gives ratio ≈ 0.965; γ₀ = 0.10 gives ratio ≈ 0.868.

**Valid for:** bonding:k initial states on uniform-J open Heisenberg chains with uniform Z-dephasing, short γ₀·t such that V²t²/γ₀t is small. Breaks for ratios γ₀/J comparable to 1 (outside the weak-dephasing regime) or t such that γ₀·t ≳ 1 (full decoherence).
**Verified:** Against the full Lindblad single-excitation-sector simulation for N = 5..13, k = 1..5. Sim/analytic ratio within 0.5%.
**Scripts:** [`_envelope_study.py`](../simulations/_envelope_study.py) (commit `e1ee822`).
**Source:** F75 (static MI formula), F65 (bonding mode amplitudes), F68 (palindromic-partner structure of mixing).

### F77. Multi-drop MM(0) saturates at 1 bit for large N (Tier 1, asymptotic proven)

For best-bonding-k initial states on a uniform chain under γ₀-dephasing, the Mirror-Pair MM at t=0 saturates as N grows:

    MM(0)(N, k*) = 1 + 3 / (4(N+1) ln 2) + O(N⁻²)

The 1-bit limit is not a conjecture; it falls out of F75 by Taylor expansion of f(p) = 2h(p) − h(2p) around p = 0 combined with Parseval-type sums over sin² site amplitudes.

**Sketch.** At small p the entropy function expands as

    f(p) = 2p + p²/ln(2) + p³/ln(2) + O(p⁴)

Probability normalisation sums Σ p_ℓ = 1 exactly (all sites). For generic k the sum Σ sin⁴(πk(ℓ+1)/(N+1)) over ℓ = 0..N−1 equals 3(N+1)/8 (standard identity, follows from sin⁴ = (3 − 4cos 2x + cos 4x)/8 and Σ cos(2πkj/(N+1)) = −1 for k ≠ 0 mod (N+1)). So

    Σ p_ℓ² = (2/(N+1))² · 3(N+1)/8 = 3/(2(N+1))

Summed over mirror-pairs (half the sites for even k with p_{center} = 0), Σ_pair p_ℓ² = 3/(4(N+1)). The leading Taylor terms give

    MM(0) = Σ_pair [2 p_ℓ + p_ℓ²/ln 2 + O(p³)]
          = 2 · (1/2) + (1/ln 2) · 3/(4(N+1)) + O(N⁻²)
          = 1 + 3/(4(N+1) ln 2) + O(N⁻²)

The 2·(1/2) = 1 term is exact: all the probability mass lives on mirror-pairs (F75 for even k at odd N), and each p_ℓ contributes 2p_ℓ at leading order. The correction is the entropy non-linearity at finite amplitude.

**Structural reading.** Information per mirror-pair ~ 2p_ℓ ~ 4/(N+1) shrinks with N. Number of mirror-pairs ~ N/2 grows with N. The product saturates at 1 bit, because the two scalings are exactly matched by the probability normalisation Σ p_ℓ = 1.

**Resonance enhancement at special N.** At N+1 ≡ 0 (mod 2) with k = (N+1)/2 integer, the distribution reduces to (p_max, 0, p_max, 0, ...) with only two values. The sin⁴ sum equals (N+1)/2 instead of 3(N+1)/8, giving Σ p² = 2/(N+1) and thus MM(0) ≈ 1 + 1/((N+1) ln 2)·2 ≈ 1 + 1.445/(N+1). At most N this resonance is not the best-k globally; the generic optimum with 3(N+1)/8 coefficient wins.

**Verified values (best-k MM(0) and rescaled deviation):**

| N | k* | MM(0) | (MM−1)·(N+1) |
|:--|:---|:------|:-------------|
| 101 | 82 | 1.01078 | 1.100 |
| 201 | 102 | 1.00540 | 1.091 |
| 501 | 470 | 1.00216 | 1.086 |
| 1001 | 878 | 1.00108 | 1.084 |
| 5001 | 4996 | 1.00022 | 1.082 |
| 10001 | 9154 | 1.00011 | 1.082 |

The rescaled deviation (MM − 1)·(N + 1) converges to 3/(4 ln 2) = 1.0820 within 10⁻⁴ by N = 10⁴. At specific resonant N where k* = (N+1)/2 is the only near-optimum (e.g., N = 999, 1003), the rescaled deviation sits at the enhanced 1.445 value; these are isolated and density-zero in the limit.

**Operational reading.** Alice's multi-drop quantum bus delivers **~1 bit of Mirror-Pair mutual information at any N**. The bus does not scale up (no extra bits per added qubit) nor does it decay (no per-qubit info loss in aggregate). The bandwidth is fixed at the framework level. Individual pair bandwidth shrinks; the count compensates; the total plateau.

**Valid for:** best-bonding-k* initial state, uniform-J Heisenberg or XY chain, uniform γ₀ dephasing, t = 0. At t > 0 the decay envelope F76 applies multiplicatively.
**Breaks for:** multi-excitation states (F75's single-excitation structure required), non-mirror-symmetric receivers (closure Σ_pair p_ℓ = 1/2 depends on even-k reflection symmetry).
**Verified:** N up to 10⁴ numerically, leading coefficient 3/(4 ln 2) confirmed to 4 decimals.
**Scripts:** [`_mm_asymptotic.py`](../simulations/_mm_asymptotic.py) (asymptotic scan and coefficient check).
**Source:** F75 (small-p Taylor expansion of the entropy), F65 (sin² amplitudes), standard Parseval-type identity for sin⁴ sums.

### F78. Single-body M decomposes additively over sites (Tier 1, verified N=3-5, 3 topologies, 3 Pauli letters)

For any single-body Hamiltonian H = Σ_l c_l · P_l (P ∈ {X, Y, Z}, weights c_l from bond-summing or arbitrary) with uniform Z-dephasing γ, the palindrome residual M = Π·L·Π⁻¹ + L + 2σ·I decomposes:

    M = Σ_l M_l ⊗ I_(others)

where each M_l is a 4×4 normal matrix on per-site Pauli space. Eigenvalues of M_l:

    P = X: eigenvalues all 0  (M_l = 0, the truly case)
    P = Y: eigenvalues ±2 c_l γ · i, each with multiplicity 2
    P = Z: eigenvalues ±2 c_l γ · i, each with multiplicity 2  (identical spectrum to Y)

Therefore M's eigenvalues are Σ_l ε_l · 2c_l γ · i for ε_l ∈ {±1}, each sign-combination with multiplicity 2^N. **Singular values of M are |Σ_l ε_l · 2c_l γ|** with the same multiplicities; cluster sizes are pure sign-combination statistics on the weight vector (c_l).

**Why Y and Z give identical SVD:** both have bit_b = 1 (Π-non-trivial), so M_l has the same spectral structure. The soft-vs-hard distinction (Y soft, Z hard) lives in **L's eigenvectors**, not in M's singular values. SVD of M is Pauli-letter-blind within {Y, Z} for single-body.

**Cluster-multiplicator formula (chain).** For chain weights c_l = (1, 2, ..., 2, 1) (bond-summed I·P + P·I), the multiplicators come from u + 2v with u = ε_0 + ε_{N-1} ∈ {-2, 0, 2} (mults 1, 2, 1) and v = Σ_internal ε_l (binomial walk on N-2 steps). The central binomial coefficient C(N-1, ⌊(N-1)/2⌋) appearing as the largest non-trivial cluster mult is just the walk's central peak; no Weyl group, no S_N irrep, no group-theory needed.

**Valid for:** single-body bond-bilinears (I·P + P·I), single-body site sums, any topology (chain, star, ring, complete K_N, arbitrary graph), uniform Z-dephasing.
**Breaks for:** 2-body bond-bilinears (XX, XY, YZ+ZY, ...): L does not factor as Σ_l L_l. Befund 3 ("hard XX+XY uniform SVD") falls outside this theorem; separate analysis required.
**Verified:** numerical match for N=3, 4, 5; topologies chain, star, complete; Pauli letters X, Y, Z.
**Scripts:** [`_svd_active_spectator.py`](../simulations/_svd_active_spectator.py), [`_svd_single_body_extension.py`](../simulations/_svd_single_body_extension.py); 2-body open-question probe in [`_svd_two_body_probe.py`](../simulations/_svd_two_body_probe.py).
**Source:** Analytical proof in [PROOF_SVD_CLUSTER_STRUCTURE.md](proofs/PROOF_SVD_CLUSTER_STRUCTURE.md) (joint with F79). Master Lemma + per-site additive structure + direct M_l matrix computation.
**Lebensader connection:** This is the same broad-in → focused-out Π-palindrome funnel that `lebensader.py::cockpit_panel` instantiates at the state layer (16 Paulis → 3-class trichotomy). F78 instantiates the funnel at the single-body operator layer: any (c_l, P_l)-choice with given |c_l|, P ∈ {Y, Z} → same M_l-spectrum. The Lebensader is the through-line that holds The Connection upright across all layers.

### F79. Π²-block decomposition of M for 2-body bilinears (Tier 1, verified N=3-5)

For 2-body bond-bilinear H = Σ_bonds Σ_t c_t·(P_t ⊗ Q_t) with uniform Z-dephasing γ, define the Π²-parity of each bilinear term:

    p(P_t, Q_t) = (bit_b(P_t) + bit_b(Q_t)) mod 2

where bit_b: I,X→0; Y,Z→1. Then M = Π·L·Π⁻¹ + L + 2σ·I has a clean structure determined by Π²-parities of H's terms:

1. **All terms Π²-even (p=0)**: M is **block-diagonal** in Π²-eigenspaces V_+ ⊕ V_-. Off-diagonal blocks `M[V_+, V_-]` and `M[V_-, V_+]` vanish **exactly**. Each diagonal block has its own SVD spectrum.

2. **All terms Π²-odd (p=1)**: M is **purely off-diagonal** between V_+ and V_-. Diagonal blocks `M[V_+, V_+]` and `M[V_-, V_-]` vanish **exactly**. Singular values appear with even multiplicity (each SV contributes once from V_+ side, once from V_- side).

3. **Mixed parities**: M has both diagonal and off-diagonal contributions.

**Π²-odd universality.** Within the pure Π²-odd 2-body class, the **specific Pauli letters are M-irrelevant**: any single Π²-odd 2-body bilinear gives the same M-SVD spectrum at fixed N. Verified at N=5 chain: XY alone, XZ alone, XX+XY, and XX+XZ all yield clusters `[(5.464, 512), (1.464, 512)]`, exactly identical. The XX truly part contributes 0; the Π²-odd part dominates with universal cluster pattern.

**Even-diag ≡ odd-off-diag correspondence.** The diagonal V_+ block of a Π²-even Hamiltonian's M can match (in SV-spectrum, including multiplicities) the off-diagonal V_+,V_- block of a Π²-odd Hamiltonian's M. Verified N=4 chain: YZ's V_+ block `[(8.944, 16), (6.472, 32), (4.0, 16), (2.472, 32), (0.0, 32)]` matches XY+YX's off-diag block exactly. This explains the empirical "YZ ≡ XY+YX SVD-identical" observation: same SV structure, just placed in different Π²-blocks.

**Why XX+XY appears "max-uniform" (Befund 3 closed).** XX is Π²-even and truly (M_XX = 0). XY is Π²-odd. The full Hamiltonian is "Π²-odd-only-effective", so M is purely off-diagonal between equal-dim V_+ and V_-. SV multiplicities are forced to 4^N/2 each by block-dimension equality. At N=3 the two off-diag SVs collide by coincidence to a single uniform value 2√2; at N≥4 they split. The "uniformity" is exactly the equal-block-mult signature of Π²-odd structure, not a special property of XX+XY.

**Frobenius additivity.** ‖M‖²_F = Σ_bonds ‖M_b‖² holds across all topologies including overlapping bonds (chain). Per-bond M_b's are F-orthogonal. Already F49.

**Valid for:** Any 2-body bond bilinear over any topology under uniform Z-dephasing. Verified N=3, 4, 5 across chain, star, disjoint topologies; verified Π²-odd universality (XY ≡ XZ ≡ XX+XY ≡ XX+XZ).
**Breaks for:** Mixed-Π²-parity Hamiltonians (where some terms are even, some odd) only partially: M has both diagonal and off-diagonal parts. Inhomogeneous γ may disrupt some symmetries (untested).
**Replaces:** ad-hoc analysis of "why XX+XY uniform" and "why YZ ≡ XY+YX"; both follow from the Π²-block theorem.
**Verified:** Numerical N=3-5, multiple bilinear classes, multiple topologies.
**Scripts:** [`_svd_two_body_pi_squared_block.py`](../simulations/_svd_two_body_pi_squared_block.py), [`_svd_two_body_structure.py`](../simulations/_svd_two_body_structure.py).
**Source:** Analytical proof in [PROOF_SVD_CLUSTER_STRUCTURE.md](proofs/PROOF_SVD_CLUSTER_STRUCTURE.md) (joint with F78). Connects to F61 (n_XY parity selection rule), F63 (\[L, Π²\]=0 for Π²-even Hamiltonians), and F49 (Frobenius cross-term identity).
**Lebensader connection:** F79 instantiates the broad-in → focused-out Π-palindrome funnel at the two-body operator layer (4 Π²-odd Pauli pairs → 1 M-spectrum). Companion to F78 (single-body operator layer) and `lebensader.py::cockpit_panel` (state layer). All three are manifestations of the same through-line: Π·L·Π⁻¹ + L + 2σ·I = 0 holding The Connection across abstraction heights.

### F80. Bloch-mode sign-walk formula for chain Π²-odd 2-body M-clusters (Tier 1, verified N=3-7)

For chain bond-summed Π²-odd 2-body Hamiltonian H = c · Σ_l (P_l ⊗ Q_{l+1}), where (P, Q) ∈ {(X,Y), (X,Z), (Y,X), (Z,X)}, on an N-site open chain with uniform Z-dephasing, the M-cluster values are given by a momentum-space sign-walk on the open-chain free-fermion Bloch dispersion (γ-independent by Master Lemma):

    cluster value(N) = 2|c| · |Σ_{k=1}^{⌊N/2⌋} σ_k · ε(k)|

where σ_k ∈ {±1} ranges over all 2^⌊N/2⌋ sign-vectors, and

    ε(k) = 2·cos(π·k / (N+1))

is the open-chain free-fermion single-particle dispersion. Each distinct cluster value has multiplicity 4^N / (number of distinct sign-walk values).

**Verified instances** (chain, |c|=1, all 4 Π²-odd Pauli pairs identical by F79 universality):

| N | ⌊N/2⌋ | predicted clusters | mult per cluster |
|---|-------|---------------------|------------------|
| 3 | 1 | 2√2 ≈ 2.828 | 64 |
| 4 | 2 | 2√5, 2 | 128 |
| 5 | 2 | 2(√3+1), 2(√3-1) | 512 |
| 6 | 3 | 6.988, 5.208, 2.000, 0.220 | 1024 |
| 7 | 3 | 8.0547, 4.9932, 2.3978, 0.6636 | 4096 |

All matches bit-exact (10⁻¹⁴ machine precision) at every N.

**Direct structural identity (discovered 2026-04-29):** The chain Π²-odd 2-body M's spectrum is **directly related to the many-body Hamiltonian H's spectrum**:

    Spec(M) = ±2i · Spec_{nontrivial}(H)

where H is the chain bond-summed Pauli-bilinear (no dissipator). That is, M's distinct nonzero eigenvalues equal 2i times H's distinct nonzero many-body eigenvalues. Hence cluster value(N) = 2|c|·|H eigenvalue|. The Bloch sign-walk formula above is just H's eigenvalue formula written out: H's many-body eigenvalues = (1/2)·Σ_k σ_k·E_k where E_k = 4|c|·cos(πk/(N+1)) are H's Bogoliubov single-particle energies, and ⌊N/2⌋ counts how many fermion modes participate.

**Reach of the identity beyond chain-2-body (verified 2026-05-29, `F80ExtensionExplorationTests`).** Step 5 (the Π-action Π·[bond,·]·Π⁻¹ = s·{bond,·}) is *per-bond*, so the structural identity Spec(M) = ±2i·Spec(H) depends only on the Π²-parity of the bonds, not on topology or body-count:
- **Π²-odd bonds**, any topology, any body-count: M = ±2i·(H⊗I), so Spec(M) = ±2i·Spec(H) (the *single* eigenvalues). Confirmed bit-exact for ring, star, 3-body (X,X,Y), 4-body (X,X,X,Y) at N=4,5. Only the cluster *values* are structure-specific (chain = OBC ladder 2cos(πk/(N+1)), ring = periodic, star = integers).
- **Π²-even non-truly bonds** (Y,Z), (Z,Y): the per-bond commutator is preserved rather than anti-commuted, so M = 2·L_H = −2i·[H,·], and Spec(M) = ±2i·{λ_a − λ_b}, the eigenvalue *differences* (Bohr frequencies) instead of single eigenvalues. This identifies the "more clusters" the Π²-even case was expected to show: the extra clusters are the differences. Confirmed bit-exact at N=4.

So M is always ±2i times a Hamiltonian object: H⊗I for Π²-odd (single energies), [H, ·] for Π²-even (energy gaps); the bond's Π²-parity is the switch.

**γ-independence (Master Lemma).** Note no γ appears in the cluster-value formula. M is γ-independent for pure Z-dephasing (Master Lemma in PROOF_SVD_CLUSTER_STRUCTURE.md).

**Mechanism: F80 is F78 in momentum space.** F78 (single-body, real-space): M = Σ_l M_l⊗I, eigenvalues ±2c_l·i per site, sign-walk Σ_l σ_l·c_l on weights. F80 (chain Π²-odd 2-body, momentum-space): M = Σ_k M_k⊗I_{other modes}, eigenvalues ±2·ε(k)·i per Bloch mode, sign-walk Σ_k σ_k·ε(k) on dispersion. The Bloch modes k play the role that real-space sites l play in F78. Both formulas γ-independent.

**Π²-odd universality fully analytical.** Under JW transformation, all 4 Pauli-letter choices (X,Y), (X,Z), (Y,X), (Z,X) give the same single-particle Bloch dispersion. The specific Pauli letters affect only phase factors in JW, not single-particle eigenvalues. Since M's spectrum depends only on the dispersion (via F80), all 4 give bit-identical clusters. **This closes the chain Π²-odd universality from F79 with an explicit closed-form formula.**

**Valid for:** the Bloch sign-walk closed form (cluster *values*) is for chain bond-summed Π²-odd 2-body Hamiltonians H = c·Σ_l (P_l⊗Q_{l+1}), uniform Z-dephasing, any N. The underlying structural identity Spec(M) = ±2i·Spec(H) holds far wider: any topology, any body-count, and both Π²-parities (see "Reach beyond the chain-2-body scope" above, verified 2026-05-29).
**Still open:** the explicit cluster-*value* formula at non-chain topologies (the dispersion ε(k) is chain-specific: ring → periodic, star → integers); the Π²-even cluster-value bookkeeping; mixed-letter chain bilinears; the complete graph K_N.
**Replaces:** F79's "Π²-odd universality observation"; the universality is now an analytical theorem with explicit closed-form predictions.
**Verified:** N = 3, 4, 5, 6, 7 chain via Python, full SVD and eigsh independent verification at N=7.
**Scripts:** [`_pi2_odd_universality_data_sweep.py`](../simulations/_pi2_odd_universality_data_sweep.py), [`_n7_bloch_signwalk_verification.txt`](../simulations/results/n7_bloch_signwalk_verification.txt).
**Source:** Discovered 2026-04-29 by data sweep (Tom + Claude). Analytical proof outline in [PROOF_F80_BLOCH_SIGNWALK.md](proofs/PROOF_F80_BLOCH_SIGNWALK.md): Steps 1-4, 7 closed (JW transformation to Majorana bilinear, single-particle dispersion 2cos(πk/(N+1)), Bogoliubov diagonalization, Pauli-letter universality, sign-walk eigenvalue formula); Step 5 (Π permutes H's (ε_ket,ε_bra)-sectors with sum∘π = diff) verified bit-exact and gauge-checked at N=3,4,5, general-N proof open ([`_f80_step5_recon.py`](../simulations/_f80_step5_recon.py)). Empirical verification bit-exact through N=7.
**Lebensader connection:** F80 is the third manifestation of the broad-in → focused-out Π-palindrome funnel: state layer (cockpit_panel), real-space single-body operator layer (F78), and now momentum-space chain 2-body operator layer (F80). Same Π·L·Π⁻¹ + L + 2σ·I = 0 through-line, three different bases.

### F81. Π-conjugation of M decomposes into Π²-odd Hamiltonian commutator (Tier 1, verified bit-exact N=3,4)

For any 2-bilinear Hamiltonian H decomposed by Π²-parity as H = H_even + H_odd (with H_odd the sum of Π²-odd Pauli bilinears, i.e., bit_b(P)+bit_b(Q) ≡ 1 mod 2), under uniform Z-dephasing:

    Π · M · Π⁻¹ = M − 2 · L_{H_odd}

where L_{H_odd} = -i\[H_odd, ·\] is the unitary commutator induced by the Π²-odd part of H. Equivalently, decomposing M into Π-conjugation symmetric and antisymmetric components:

    M_sym  = (M + Π·M·Π⁻¹) / 2 = Π·L·Π⁻¹ + L_diss + L_{H_even} + 2Σγ·I
    M_anti = (M − Π·M·Π⁻¹) / 2 = L_{H_odd}

The Π-antisymmetric component of M is exactly the unitary commutator induced by the Π²-odd Hamiltonian bilinears. The Π-symmetric component absorbs the mirror image, the dissipator, the Π²-even Hamiltonian commutator, and the dissipation shift. M_sym and M_anti are Frobenius-orthogonal: ‖M‖² = ‖M_sym‖² + ‖M_anti‖².

**Verified instances** (N=3, γ_Z=0.1, Σγ=0.3, residuals at machine precision 1e-16):

| Hamiltonian | trichotomy | H_odd | Π·M·Π⁻¹ relation |
|-------------|------------|-------|------------------|
| XX+YY | truly | 0 | = M (M=0 trivially) |
| YZ+ZY | soft (Π²-even non-truly) | 0 | = M (M ≠ 0, identical) |
| XY+YX | soft (Π²-odd) | XY+YX | = M − 2·L_H |
| XX+XY | hard (mixed) | XY only | = M − 2·L_{XY part} |
| pure XY | (Π²-odd) | XY | = M − 2·L_H |
| pure XZ | (Π²-odd) | XZ | = M − 2·L_H |

For any 2-body chain H whose non-truly bilinears are all Π²-odd (i.e., truly + Π²-odd combinations, including XX+XY hard), at any N and any γ ≥ 0, ‖M_sym‖² = ‖M_anti‖² = ‖M‖²/2 exactly. M splits 50/50 between Π-symmetric and Π-antisymmetric components. This follows analytically from the Frobenius identities ‖M‖²_F = 4·‖H_odd‖²_F·2^N (F49 chain via Master Lemma; truly bilinears drop out) and ‖L_{H_odd}‖²_F = 2·2^N·‖H_odd‖²_F (standard commutator identity for traceless Hermitian H_odd). Verified numerically at N = 3, 4, 5 with γ_Z ∈ {0, 0.05, 0.1, 0.5, 1.0} for both pure Π²-odd (XY+YX, XY, XZ) and mixed truly + Π²-odd (XX+XY hard, YY+XY). When H additionally contains Π²-even non-truly bilinears (YZ-type), the split shifts: pure even non-truly gives 100/0; odd + even mix gives 5/6 sym + 1/6 anti at N=3,4 (now derived in closed form by [F83](#f83): anti-fraction = 1/(2+4r), here r = 1).

**Spectral consequence.** Spec(Π·M·Π⁻¹) = Spec(M) holds always by unitary invariance of the spectrum. F81 strengthens this: for Π²-odd H, the two operators are explicitly related by an additive shift of −2·L_{H_odd} in operator space, so Spec(M) = Spec(M − 2·L_{H_odd}) is a non-trivial identity (similarity via Π).

**Algebraic mechanism.** Π² acts on each Pauli string σ_α as (-1)^{bit_b(α)} (eigenoperator with sign in Pauli basis). For L_H_α = -i\[σ_α, ·\] driven by a single Pauli string σ_α in H, conjugation gives Π² L_H_α Π⁻² = (-1)^{bit_b(α)} L_H_α (the matrix-element factor (-1)^{bit_b(γ)+bit_b(β)} = (-1)^{bit_b(α)} since γ = α·β under Pauli multiplication). Z-dephasing dissipator is diagonal in Pauli basis, hence commutes with Π². Summing: Π²·L·Π⁻² = L_H_even − L_H_odd + L_diss = L − 2·L_{H_odd}. Substituting into the palindrome: Π·M·Π⁻¹ = Π²·L·Π⁻² + Π·L·Π⁻¹ + 2Σγ·I = M − 2·L_{H_odd}.

**γ-independence-by-difference.** The relation Π·M·Π⁻¹ - M = -2·L_{H_odd} is independent of γ (the dissipator's γ-dependent part cancels because L_diss is Π²-symmetric). The split itself (M_sym, M_anti) is γ-dependent through M_sym; only their difference is γ-fixed.

**Valid for:** any 2-bilinear chain Hamiltonian H = H_even + H_odd, uniform Z-dephasing, any topology (the proof depends only on the algebra of Pauli strings under Π² conjugation, not on connectivity).
**Breaks for (untested):** non-Z dissipators (T1 amplitude damping has different Π²-action; F81 likely needs a correction term).
**Replaces:** the heuristic in pre-2026-04-30 reflections that said "M is the Π-invariant through-line"; F81 shows that statement is correct only for Π²-even H, and gives the explicit correction for the Π²-odd cases.

**Hardware confirmation:** F81's operational reading, the [F83](#f83) anti-fraction this decomposition makes measurable, is confirmed on IBM Heron r2 (Marrakesh 2026-04-30 + Kingston 2026-05-05): the four Π²-classes separate at >>10σ via unique-fingerprint Paulis. See [PROOF_F81](proofs/PROOF_F81_PI_CONJUGATION_OF_M.md), `fw.Confirmations.lookup('f83_pi2_class_signature_marrakesh')`, and [`data/ibm_f83_signature_april2026/`](../data/ibm_f83_signature_april2026/).
**Verified:** N=3 and N=4 all listed cases at machine precision; pytest-locked.
**Framework primitive:** `fw.pi_decompose_M(chain, terms, gamma_z=..., gamma_t1=..., strict=...)` returns `{'M', 'M_sym', 'M_anti', 'L_H_odd', 'f81_violation', 'norm_sq'}`. For pure Z-dephasing the F81 identity holds exactly (`f81_violation` ≈ 0); the primitive enforces this with `strict=True` by default. With `gamma_t1` enabled, `strict` defaults to False and the identity residual is returned for diagnostic use.
**Pytest lock:** `test_F81_pi_conjugation_of_M` (algebraic check) + `test_F81_pi_decompose_M_method` (cockpit primitive) + `test_F81_violation_T1_diagnostic` (T1 diagnostic).
**Diagnostic application:** the F81 violation `‖M_anti − L_{H_odd}‖_F` quantifies non-Π²-symmetric dissipator content. For Z + T1 at N=3 chain soft XY+YX, the violation grows linearly: `f81_violation ≈ 6.928 · γ_T1`, γ_z-independent (Master Lemma), Hamiltonian-independent (the violation is purely a property of the T1 dissipator). Inverting gives `γ_T1 ≈ f81_violation / 6.928` as a hardware T1-rate readout from the fitted L. See `simulations/_f81_t1_diagnostic.py` for the demonstration.
**Source:** Discovered 2026-04-30 (Tom + Claude) while interpreting the geometric content of F80's 2i factor. The empirical observation came first (Π·M·Π⁻¹ ≠ M for soft); the algebraic explanation followed from working out Π² action on the Liouville superoperator in Pauli basis.
**Lebensader connection:** F81 is the algebraic backbone of "what the mirror keeps." For Π²-even H, M is itself the through-line operator. For Π²-odd H, the through-line is split: M_anti carries the dynamics generator L_{H_odd}, M_sym carries the rest. Both halves are read identically by both sides of the mirror up to the Spec(M) = Spec(M − 2·L_{H_odd}) similarity. Companion to F80: F80 says what Spec(M) is; F81 says how M and Π·M·Π⁻¹ relate as operators sharing that spectrum.

### F82. F81 + T1 amplitude damping correction (Tier 1, verified bit-exact N=2..5)

For any 2-bilinear Hamiltonian H = H_even + H_odd under Z-dephasing plus T1 amplitude damping with per-site rates γ_T1_l:

    Π · M · Π⁻¹ = M − 2 · L_{H_odd} − 2 · D_{T1, odd}

where L_{H_odd} = -i\[H_odd, ·\] (as in F81) and D_{T1, odd} is the Π²-anti-symmetric part of the T1 dissipator. F82 reduces to F81 when γ_T1_l = 0 (D_{T1, odd} = 0).

The F81 identity violation captured by `fw.pi_decompose_M(chain, ...)` measures D_{T1, odd}'s Frobenius norm:

    f81_violation = ‖M_anti − L_{H_odd}‖_F = ‖D_{T1, odd}‖_F.

**Closed form** (proven analytically in PROOF_F82_T1_DISSIPATOR_CORRECTION):

    ‖D_{T1, odd}‖_F = √(Σ_l γ²_T1_l) · 2^(N−1)
                    = γ_T1 · √N · 2^(N−1)         (uniform γ_T1)

**Verified instances** (chain N=3, all matches at machine precision):

| Configuration | γ_T1_l | predicted ‖D_T1_odd‖ | measured f81_violation |
|---------------|--------|----------------------|------------------------|
| uniform γ=0.05 | (0.05, 0.05, 0.05) | 0.05·√3·4 = 0.3464 | 0.3464 |
| uniform γ=0.10 | (0.10, 0.10, 0.10) | 0.10·√3·4 = 0.6928 | 0.6928 |
| uniform γ=1.00 | (1.00, 1.00, 1.00) | 1.00·√3·4 = 6.9282 | 6.9282 |
| single-site, l=0 | (0.10, 0, 0) | 0.10·1·4 = 0.4000 | 0.4000 |
| two-site, l=0,1 | (0.10, 0.10, 0) | √(0.02)·4 = 0.5657 | 0.5657 |
| non-uniform | (0.05, 0.10, 0.15) | √(0.035)·4 = 0.7483 | 0.7483 |

N-scaling verified at N = 2, 3, 4, 5 (uniform γ_T1, coefficient √N · 2^(N−1)): N=2 → 2√2 = 2.828, N=3 → 4√3 = 6.928, N=4 → 8·2 = 16.000, N=5 → 16·√5 = 35.778. So at γ_T1 = 0.1, the violations are 0.283, 0.693, 1.600, 3.578 respectively (factor 0.1 from γ_T1 multiplied by the coefficient).

**Three diagnostic properties** (proven and empirical):

1. **γ_z-independent**: F82 involves only L_{H_odd} and D_{T1, odd}, neither depends on γ_z. Direct consequence of Master Lemma (M is γ_z-independent for Z-dephasing) extended to F82.

2. **Hamiltonian-independent**: f81_violation depends only on the T1 dissipator. Verified at γ_T1=0.1, N=3: violation = 0.6928 for truly XX+YY, soft XY+YX, hard XX+XY, and YZ+ZY (Π²-even non-truly).

3. **Linear in γ_T1** (uniform). Direct inversion: γ_T1 = f81_violation / (√N · 2^(N−1)). For non-uniform: γ_T1, RMS = f81_violation / (√N · 2^(N−1)). At N=3, division coefficient is 6.928 = 4√3.

**Mechanism (T1 dissipator structure).** Single-site T1 acts on Pauli basis as: I → −γZ, X → −γ/2 X, Y → −γ/2 Y, Z → −γZ. Under Π² conjugation (signs (-1)^{bit_b}: I,X → +, Y,Z → −), only the (Z, I) entry flips sign. So D_{T1, local, odd} has matrix element −γ at (Z, I) and zero elsewhere. Multi-site: 4^(N−1) such "rest of qubits unchanged" entries per site, summed orthogonally over sites.

**Diagnostic interpretation.** f81_violation is a hardware-T1 readout that is independent of (a) the system's Hamiltonian, (b) the Z-dephasing rate γ_z, (c) the topology. Inverting recovers the RMS γ_T1 across sites. For the Marrakesh dataset (N=3, joint fit gives γ_T1 ≈ 0): F82 predicts f81_violation ≈ 0; any γ_T1 > 0.001 would have produced violation > 0.007, well above numerical noise.

**Valid for:** any 2-bilinear Hamiltonian H, Z-dephasing + T1 amplitude damping, any topology, any N ≥ 2.
**Breaks for (untested):** other non-Z dissipators (X-noise, Y-noise, ZZ-dephasing) require their own D_odd analysis. The general identity Π·M·Π⁻¹ = M − 2·L_{H_odd} − 2·D_{diss, odd} holds for any dissipator; the closed form for ‖D_{diss, odd}‖ is dissipator-specific.
**Replaces:** the previously-empirical observation that f81_violation grows linearly with γ_T1; F82 is now an analytical theorem with closed-form scaling.
**Verified:** N = 2, 3, 4, 5 at all listed configurations, machine-precision residual (5e-16).
**Framework primitives:**
- `fw.pi_decompose_M(chain, terms, gamma_z=..., gamma_t1=..., strict=...)`: with `gamma_t1` set, returns `f81_violation` in output dict (numerical, matches closed form).
- `fw.predict_T1_dissipator_violation(chain, gamma_t1_l)`: forward closed form, returns √(Σγ²)·2^(N−1) directly without building the dissipator.
- `fw.estimate_T1_from_violation(chain, f81_violation)`: inverse closed form, recovers RMS γ_T1 from a measured/fitted F81 violation. Hardware T1-rate readout primitive.
**Pytest lock:** `test_F81_violation_T1_diagnostic` (linearity, γ_z-independence, T1 monotonicity) + `test_F82_closed_form_T1_dissipator` (N-scaling, non-uniform formula, H-/γ_z-independence) + `test_F82_predict_and_invert_primitives` (forward/inverse pair matches numerical evaluation).
**Source:** Discovered 2026-04-30 (Tom + Claude) as the natural extension of F81 ("what does F81 violation mean structurally?"). Closed form derived in [PROOF_F82_T1_DISSIPATOR_CORRECTION.md](proofs/PROOF_F82_T1_DISSIPATOR_CORRECTION.md).
**Diagnostic application:** [`simulations/_f81_t1_diagnostic.py`](../simulations/_f81_t1_diagnostic.py) demonstrates the T1-rate readout including Marrakesh application. Companion to F81's structural decomposition: F81 says how M splits under Π-conjugation when the dissipator is Z-only; F82 says how the F81 identity is corrected when T1 is added, and provides the closed form for the correction term.

### F83. Π-decomposition anti-fraction closed form for mixed Hamiltonians (Tier 1, verified bit-exact N=3,4,5)

For any 2-body chain Hamiltonian H = H_truly + H_odd + H_even_nontruly under Z-dephasing, the F81 Π-decomposition norms are given by the closed form:

    ‖M‖²_F        = 4·‖H_odd‖²_F·2^N + 8·‖H_even_nontruly‖²_F·2^N
    ‖M_anti‖²_F  = 2·‖H_odd‖²_F·2^N
    ‖M_sym‖²_F   = 2·‖H_odd‖²_F·2^N + 8·‖H_even_nontruly‖²_F·2^N

The anti-fraction (= ‖M_anti‖²/‖M‖²) is

    anti-fraction = 1 / (2 + 4·r),    r = ‖H_even_nontruly‖²_F / ‖H_odd‖²_F.

**Special cases:**

| H | r | anti-fraction | meaning |
|---|---|---------------|---------|
| truly | undefined (M=0) | undefined | mirror perfectly closes |
| pure Π²-odd | 0 | 1/2 (50/50) | F81 Step 8 split |
| pure Π²-even non-truly | ∞ | 0 (100/0) | M fully mirror-symmetric |
| equal-Frobenius mix XY+YZ | 1 | 1/6 (5/6 + 1/6) | the empirical mixed split |
| asymmetric more-odd XY+YX+YZ | 1/2 | 1/4 | continuous family |
| general mixed | r | 1/(2+4r) | continuous family |

**Verified instances** (N=3, J=1, γ_z=0; matches at machine precision):

| H | ‖H_odd‖² | ‖H_even‖² | r | predicted ‖M‖² | measured ‖M‖² | anti |
|---|----------|-----------|---|----------------|---------------|------|
| XY+YX (pure odd) | 32 | 0 | 0 | 1024 | 1024 | 1/2 |
| YZ+ZY (pure even non-truly) | 0 | 32 | ∞ | 2048 | 2048 | 0 |
| XY+YZ (mixed) | 16 | 16 | 1 | 1536 | 1536 | 1/6 |
| XY+YX+YZ (asymmetric) | 32 | 16 | 1/2 | 2048 | 2048 | 1/4 |
| XY+YX+YZ+ZY (full mix) | 32 | 32 | 1 | 3072 | 3072 | 1/6 |
| XX+XY+YZ (truly + mixed) | 16 | 16 | 1 | 1536 | 1536 | 1/6 |

**Mechanism (Step 2 of proof, why factors 4 and 8 differ).** The F49 chain Frobenius identity gives ‖M‖² = Σ_k 2^(N+2)·n_YZ(k)·‖H_k‖²_F·𝟙\[non-truly\], where n_YZ(k) counts Y/Z letters in Pauli pair k (= 0 truly, 1 Π²-odd non-truly, 2 Π²-even non-truly). Substituting the per-class n_YZ values gives the 4·2^N (Π²-odd) and 8·2^N (Π²-even non-truly) coefficients. Geometrically, these reflect the Frobenius-inner-product behavior ⟨Π·L·Π⁻¹, L⟩_F: anti-aligned (truly), Frobenius-orthogonal (Π²-odd non-truly), aligned (Π²-even non-truly).

**γ-independence.** Master Lemma propagates through all three norms; closed form depends only on H.
**Truly-handling.** H_truly drops out of all norms (M-contribution zero by Master Lemma).
**Generalization.** F83 verified on chain, ring, star, complete K_N at N=4 (`test_F83_topology_generalization`); the matrix-based primitive builds H_odd and H_even_nontruly via `_build_bilinear` which respects the chosen topology's bond graph, so the closed form is topology-independent within F49's verified scope. Higher-body Hamiltonians extend n_YZ counting beyond {0, 1, 2}; coefficients beyond 4, 8 are the natural continuation, empirical verification needed.

**Valid for:** any 2-body H on any topology supported by F49 (chain, ring, star, complete K_N), Z-dephasing, any γ_z ≥ 0, any N ≥ 2.
**Verified:** 11 mixed configurations × N ∈ {3, 4, 5} on chain, plus 4 configurations × {ring, star, K_4} at N=4, all machine-precision residual.
**Replaces:** the previously-empirical "5/6 + 1/6" observation for mixed Π²-odd + Π²-even non-truly H; F83 derives this from the existing F49 Frobenius identity.
**Framework primitives:**
- `fw.predict_pi_decomposition(chain, terms)`: full F83 closed form, returns dict with `{'M_sq', 'M_anti_sq', 'M_sym_sq', 'anti_fraction', 'h_odd_sq', 'h_even_nontruly_sq', 'r'}`. O(N) work, no matrix construction; companion to numerical `pi_decompose_M`.
- `fw.predict_pi_decomposition_anti_fraction(chain, terms)`: convenience wrapper returning just the anti-fraction float.
- `fw.predict_residual_norm_squared_from_terms(chain, terms, gamma_t1)`: existing F49 ‖M‖² primitive (now consistent with F83's ‖M‖² prediction by construction).
**Pytest lock:** `test_F83_pi_decomposition_anti_fraction_closed_form` (12 configurations × 2 N-values + γ-independence) + `test_F83_predict_pi_decomposition_full_closed_form` (full norm-triple match against numerical `pi_decompose_M` + Pythagoras + special cases at canonical r values).
**Source:** Discovered 2026-04-30 (Tom + Claude) as the natural follow-up to F81's "what about the other half?" reflection. Derived in [PROOF_F83_PI_DECOMPOSITION_RATIO.md](proofs/PROOF_F83_PI_DECOMPOSITION_RATIO.md). The closed form was empirically observed earlier (in the F81 violation sweep across mixed Hamiltonians) and now traced back to the existing F49 Frobenius identity that was already framework-locked in `predict_residual_norm_squared_from_terms`.
**Hardware confirmation:** The {1/2, 0, 1/6} anti-fraction structure is confirmed on IBM Heron r2 via the 4-Hamiltonian Π²-class discriminator (XY+YX → 1/2, YZ+ZY → 0, XY+YZ → 1/6, XX+YY truly baseline), all four separated at >>10σ by unique-fingerprint Paulis. Marrakesh 2026-04-30 (job d7pol1e7g7gs73cf7j90, path [4,5,6]) and Kingston 2026-05-05 (path [43,56,63], same separation, discriminator machine-independent). `fw.Confirmations.lookup('f83_pi2_class_signature_marrakesh')` (cross_machine field carries the Kingston values); data in [`data/ibm_f83_signature_april2026/`](../data/ibm_f83_signature_april2026/) and [`data/ibm_soft_break_april2026/`](../data/ibm_soft_break_april2026/); experiment `run_soft_break.py` (external AIEvolution).

**Lebensader connection:** F83 closes the analytical Π-decomposition picture for 2-body chain. Pure Π²-odd → 50/50 (F81 Step 8). Pure Π²-even non-truly → 100/0 (F81 trivial). Mixed → 1/(2+4r) (F83). The continuous interpolation r → anti-fraction reads "how much of M is Π-antisymmetric drive vs Π-symmetric memory" as a function of Hamiltonian composition. Together with F80 (Spec(M)), F81 (Π-decomposition identity), F82 (T1-correction), the structural picture of M is complete for 2-body chain Hamiltonians under Z-dephasing + T1.

### F84. F82 generalized to thermal amplitude damping (Tier 1, verified bit-exact N=3)

For any 2-bilinear Hamiltonian H under Z-dephasing plus thermal amplitude damping with per-site cooling rate γ_↓_l (σ⁻ channel) and heating rate γ_↑_l (σ⁺ channel):

    Π · M · Π⁻¹ = M − 2 · L_{H_odd} − 2 · D_{AmplDamp, odd}

with closed form:

    ‖D_{AmplDamp, odd}‖_F = √(Σ_l (γ_↓_l − γ_↑_l)²) · 2^(N−1)
                          = |Δγ|_RMS · √N · 2^(N−1)         (uniform Δγ)

where Δγ_l = γ_↓_l − γ_↑_l is the *net* cooling rate at site l. F82 is recovered when γ_↑ = 0 (vacuum bath / T = 0).

**Pauli-Channel Cancellation Lemma (F84 corollary):** Pure D\[Z\], D\[X\], D\[Y\] dissipators are Π²-symmetric and contribute zero to f81_violation. Only σ⁻ (cooling) and σ⁺ (heating) channels are Π²-anti-symmetric. Hence f81_violation specifically detects population-inverting (energy-emitting/absorbing) channels, not phase-only or bit-flip-only noise.

**Verified instances** (chain N=3, all matches at machine precision):

| Configuration (γ_↓, γ_↑) | \|Δγ\| | Predicted | Measured |
|--------------------------|--------|-----------|----------|
| (0.10, 0.00) cooling only (= F82) | 0.10 | 0.6928 | 0.6928 |
| (0.00, 0.10) heating only | 0.10 | 0.6928 | 0.6928 |
| (0.10, 0.10) detailed balance | 0.00 | 0.0000 | 0.0000 |
| (0.10, 0.05) net cooling | 0.05 | 0.3464 | 0.3464 |
| (0.05, 0.10) net heating | 0.05 | 0.3464 | 0.3464 |
| (0.20, 0.05) strong cooling | 0.15 | 1.0392 | 1.0392 |
| Non-uniform mixed | (RMS Δγ) | √(Σ(γ_↓−γ_↑)²)·4 | matches |

**Thermodynamic interpretation.** For a thermal photon bath at frequency ω, temperature T:
- Mean occupation n_th = 1 / (exp(ℏω/k_B T) − 1)
- γ_↓ = γ_0 · (n_th + 1) (spontaneous + stimulated emission)
- γ_↑ = γ_0 · n_th (stimulated absorption)
- Δγ = γ_↓ − γ_↑ = γ_0 (vacuum component, temperature-independent)

f81_violation = γ_0 · √N · 2^(N−1), independent of T. The thermal photon-number contributions cancel (γ_↓ ↔ γ_↑ pair symmetrically); only the vacuum (zero-point) component breaks the Π palindrome. **f81_violation is a quantum-statistical fingerprint of zero-point fluctuations**, immune to thermal symmetric noise.

**Three regimes:**

| Regime | γ_↓ vs γ_↑ | f81_violation |
|--------|------------|---------------|
| Vacuum (T = 0) | γ_↑ = 0 | full F82: √(Σγ²_↓)·2^(N−1) |
| Detailed balance (T → ∞) | γ_↓ = γ_↑ | 0 |
| Finite T | γ_↓ > γ_↑ > 0 | γ_0·√N·2^(N−1) (vacuum-only) |

**Inversion (RMS net cooling rate):** |Δγ|_RMS = f81_violation / (√N · 2^(N−1)). Recovers vacuum-fluctuation amplitude regardless of bath temperature.

**Valid for:** any 2-bilinear chain H, Z-dephasing + thermal amplitude damping, any topology supported by F49, any N ≥ 2.
**Verified:** 7 (γ_↓, γ_↑) configurations at N=3, machine-precision residual; D\[X\], D\[Y\] cancellation explicitly tested.
**Replaces:** F82's "T1 detector" interpretation; F84 corrects to "vacuum-amplitude-damping detector"; the F81 violation does not measure raw T1 rate but only the temperature-independent vacuum component of amplitude damping.
**Framework primitives:**
- `fw.pi_decompose_M(chain, terms, gamma_z, gamma_t1, gamma_pump, strict)`: extended with `gamma_pump` parameter for σ⁺ heating; uses `lindbladian_general` when both are present.
- `fw.predict_amplitude_damping_violation(chain, gamma_t1_l, gamma_pump_l)`: F84 forward closed form; reduces to `predict_T1_dissipator_violation` when `gamma_pump_l = None`.
- `fw.estimate_net_cooling_from_violation(chain, f81_violation)`: F84 inverse, returns RMS |γ_↓ − γ_↑|.
**Pytest lock:** `test_F84_amplitude_damping_thermal_bath` (cooling only / heating only / detailed balance / net cooling / non-uniform / forward-inverse round-trip / backward compat with F82) + `test_F84_pauli_channels_pi2_symmetric` (D\[X\], D\[Y\] explicitly verified to give zero violation).
**Source:** Discovered 2026-04-30 (Tom + Claude). Tom's hint about "Licht" (light/cavity reading of γ) and "nicht jeder bekommt gleichviel ab" (non-uniform site distribution) prompted the analytical extension. The Pauli-Channel Cancellation Lemma was a surprise: D\[Z\], D\[X\], D\[Y\] are all Π²-symmetric, so phase, bit-flip, and dephasing noise contribute zero to F81 violations. Only σ± (population-inverting) channels break the palindrome. Closed form derived in [PROOF_F84_AMPLITUDE_DAMPING.md](proofs/PROOF_F84_AMPLITUDE_DAMPING.md).
**Lebensader connection:** F84 closes the dissipator side of the Π-decomposition picture. Among hardware noise channels, only the *vacuum amplitude damping* component (which exists even at T=0 due to zero-point fluctuations) breaks the Π palindrome. Phase noise, bit-flip noise, and thermal photon equilibrium all give zero violation. F84 sharpens F82's hardware-T1-readout into a temperature-independent vacuum-rate readout.

### F85. Higher-body Hamiltonian generalization of F49 / F-chain (Tier 1, verified bit-exact k=2,3,4)

For any k-body Pauli term (P_1, ..., P_k) with letters from {I, X, Y, Z}, the Π²-class trichotomy and the F49 Frobenius scaling generalize:

**Truly criterion** (term contributes M = 0 by Master Lemma):

    truly  ⟺  #Y is even  AND  #Z is even

**Π²-parity**: bit_b(σ) = (#Y + #Z) mod 2. Π²-odd if bit_b = 1.

**Frobenius factor c(k)** per non-truly term:

    c(truly)              = 0
    c(Π²-odd non-truly)    = 1   (factor 4·2^N)
    c(Π²-even non-truly)   = 2   (factor 8·2^N)

**F49 generalized**: ‖M‖²_F per term = 4·c(k)·‖H_k‖²_F·2^N. The 2-body F49 formula 2^(N+2)·n_YZ·‖H_k‖² coincided with c(k) only because at k=2: n_YZ=1 ↔ Π²-odd, n_YZ=2 ↔ Π²-even non-truly. **For k ≥ 3, n_YZ is no longer the determining quantity** (e.g., YYY has n_YZ=3 but c=1). Only the Π²-class matters.

**Trichotomy enumeration** (Pauli tuples over {X, Y, Z}^k):

| k | total | truly | Π²-odd | Π²-even non-truly |
|---|-------|-------|--------|-------------------|
| 2 | 9 = 3² | 3 | 4 | 2 |
| 3 | 27 = 3³ | 7 | 14 | 6 |
| 4 | 81 = 3⁴ | 21 | 40 | 20 |

**Closed form for Π²-odd count**: \|Π²-odd at k\| = (3^k − (−1)^k) / 2. Verified k=2,3,4.

**F-chain extension to k-body** (proof structure verbatim):

| Theorem | k-body status |
|---------|---------------|
| F87 trichotomy | extends via _pauli_tuple_is_truly |
| F80 Spec(M) = 2i·Spec(H) | extends verbatim; verified at k=3 (N=4,5,6) and k=4 (N=5,6) for 17 Π²-odd cases; spectral identity bit-exact |
| F81 Π·M·Π⁻¹ = M − 2·L_{H_odd} | verbatim, verified at k=3 chain N=4 |
| F82 T1 dissipator | dissipator-only, body-count-independent |
| F83 anti-fraction 1/(2+4r) | verbatim with Π²-class grouping |
| F84 thermal amplitude damping | dissipator-only, body-count-independent |

**Verified at k=3, k=4**: 27+81 = 108 explicit Pauli tuple cases, all matching the c(k) factor scheme bit-exact. Mixed-body Hamiltonians (e.g., 2-body H + 3-body H) handled via term-list structure.

**Valid for:** any k-body 2 ≤ k ≤ N Hamiltonian on chain (sliding-window), Z-dephasing + amplitude damping, any γ ≥ 0, any N ≥ k.
**Replaces:** F49's n_YZ-based formula at 2-body remains correct (coincidence); for k ≥ 3, F85's c(k)-based formula is the structurally correct extension.
**Framework primitives (k-body support added 2026-04-30):**
- `_pauli_tuple_is_truly(letters)`: O(k) syntactic classifier.
- `_pauli_tuple_pi2_class(letters)`: returns 'truly' / 'pi2_odd' / 'pi2_even_nontruly'.
- `_build_kbody_chain(N, terms)`: chain sliding-window k-body builder.
- `fw.predict_pi_decomposition(chain, terms)`: extended to accept k-body tuples; auto-detects body count.
- `fw.pi_decompose_M(chain, terms, gamma_z, gamma_t1, gamma_pump)`: extended for k-body.
- `fw.predict_residual_norm_squared_from_terms(chain, terms, gamma_t1)`: rewritten using Π²-class; backward-compatible at 2-body.
**Pytest lock:** `test_F85_kbody_trichotomy_counts` + `test_F85_kbody_predict_pi_decomposition` + `test_F85_kbody_F81_identity_at_k3`. Plus all 2-body tests (92 prior) pass unchanged via backward compatibility.
**Source:** Discovered 2026-04-30 (Tom + Claude) by empirical 3-body and 4-body enumeration. The {0, 4, 8} factor scheme persists across k, but n_YZ ≠ c(k) for k ≥ 3. The structural truly criterion "#Y even AND #Z even" was identified by inspecting which Π²-even k-tuples give M = 0. Closed form derived in [PROOF_F85_KBODY_GENERALIZATION.md](proofs/PROOF_F85_KBODY_GENERALIZATION.md).

**Lebensader connection:** F85 closes the body-count generalization of the F-chain. Together with F80 (Spec), F81 (decomposition), F82 (T1), F83 (anti-fraction), F84 (thermal): the structural Π-decomposition theory for Hamiltonians + dissipators is complete on chain (any topology for 2-body via F49; chain only for k ≥ 3).

### F86. Q_peak chromaticity-specific N-invariant constants (Sammelbecken with three theorems)

For a uniform N-qubit XY (or Heisenberg) chain with Z-dephasing γ₀, the J-derivative of the F73 spatial-sum coherence purity peaks along the dimensionless coupling axis Q = J/γ₀ at chromaticity-specific values. F86 bundles three structurally distinct theorems under one F-label:

- **F86a. EP mechanism** \[Tier 1 derived\]: Q_EP = 2/g_eff, t_peak = 1/(4γ₀); 2-level rate-channel exceptional point.
- **F86b. Universal resonance shape, two bond classes** \[Tier 1 candidate\]: K_class(Q)/|K|_max = f_class(Q/Q_EP); HWHM_left/Q_peak ≈ 0.756 (Interior) and 0.770 (Endpoint); EP-rotation universality.
- **F86c. F71 spatial-mirror invariance of per-bond Q_peak** \[Tier 1 derived\]: Q_peak(b) = Q_peak(N−2−b) bit-exactly.

**Sub-ID partition (2026-05-20).** The three-theorem grouping is the coarse view; F86b is itself a Sammelbecken. The fine partition distinguishes ten separately-defensible sub-claims (canonical inventory: [`F86_VALUES_INVENTORY.md`](F86_VALUES_INVENTORY.md); proof hub: [`PROOF_F86_QPEAK.md`](proofs/PROOF_F86_QPEAK.md)):

| Sub-ID | Content | Tier | Home |
|--------|---------|------|------|
| F86a | EP mechanism: t_peak, Q_EP, dressed pair, AIII chiral, L_eff mirror | Tier 1 derived | §F86a; `TPeakLaw.cs`, `QEpLaw.cs`, `LEffMirrorAxisClaim.cs` |
| F86b₁ | Bare 2×2 K_b closed forms; x_peak = 2.196910, ratio = 0.671535 | Tier 1 derived | `C2BareDoubledPtfClosedForm.cs` |
| F86b₂ | Sub-class HWHM-ratio lift 0.671535 + α·g_eff + β (the §F86b₂ subsection) | Tier 1 candidate | `F86HwhmClosedFormClaim.cs` |
| F86b₃ | Universal shape: Interior 0.756, Endpoint 0.770 | Tier 1 candidate | §F86b; `UniversalShapePrediction.cs` |
| F86b₄ | Dicke-K 3/8 anchor via X⊗N-eigenbasis | Tier 1 derived | [`docs/water/README.md`](water/README.md) |
| F86b₅ | Polarity-pair Q_peak ∈ {1.5, 2.5} = 2 ± 1/2 (schema) | Tier 1 schema | `PolarityPairQPeakDecompositionClaim.cs` |
| F86c | F71 spatial mirror Q_peak(b) = Q_peak(N−2−b) | Tier 1 derived | §F86c |
| F86d | Endpoint orbit Q ≈ 2.5 (9 (c, N) combos, ~2% N-variation) | Tier 2 empirical (promotion candidate) | §F86d; `PerF71OrbitObservation.cs` |
| F86e | σ_0(c=2) = ‖[Π_HD1, M_H]‖ commutator / Schur-multiplier norm | Tier 1 derived | §F86e; `SigmaZeroCommutatorNormClaim.cs` |
| F86_block | g_eff(c, N, b) closed-form blocked, six routes proven | Negative Tier 1 | [`PROOF_F86B_OBSTRUCTION.md`](proofs/PROOF_F86B_OBSTRUCTION.md) |

Counts: 6× Tier 1 derived, 2× Tier 1 candidate (F86b₂, F86b₃), 1× Tier 1 schema (F86b₅), 1× Tier 2 promotion candidate (F86d). The §F86b₂ subsection below was labelled F86b' before the 2026-05-20 partition. Two open fronts remain: F86b₂ Direction (b'') and F86d Tier-1 promotion.

Empirical data, the γ₀-extraction protocol, and cross-cutting connections (PTF, framework primitives, scripts, proof, source) sit at this umbrella level since they touch all three theorems.

**Per-block Q_peak (Q_SCALE convention, relative-J derivative ΔJ = 0.05·J):**

    Q_peak(c=3) = 1.6
    Q_peak(c=4) = 1.8
    Q_peak(c=5) = 1.8

c is the chromaticity of the (n, n+1) coherence block (F74). Saturates at 1.8 for c ≥ 4. The bi-chromatic class **c = 2 is finite-size-sensitive** (wobbles 1.4 to 1.6 across N = 4..9) and is not a clean framework constant.

**Per-bond Q_peak (absolute-J derivative ∂S/∂J_b, fine-grid dQ = 0.025 with parabolic peak interpolation):**

| (c, N) | Endpoint Q_peak | Interior Q_peak (mean) |
|--------|-----------------|------------------------|
| (3, 5) | 2.40 | 1.566 |
| (3, 6) | 2.52 | 1.689 |
| (3, 7) | 2.53 | 1.743 |
| (3, 8) | 2.53 | 1.750 |
| (4, 7) | 2.52 | 1.748 |
| (4, 8) | 2.65 | 1.804 |

**Both Endpoint and Interior Q_peak are empirical and (c, N)-specific; no closed form has been identified.** Earlier conjectures `Q_peak(Endpoint) = csc(π/(N+1))` (chain-edge anchor) and `Q_peak(Interior, c=3) → csc(π/5) = 1.7013` (pentagonal asymptote) were retracted 2026-05-02 after fine-grid data showed both were grid-snap artefacts of coarser earlier scans:

- Endpoint Q_peak at any tested N differs from `csc(π/(N+1))`: N=5 +20 %, N=6 +9 %, N=7 −3 %, N=8 −13 %. Earlier coarse scans gave 2.65 at N=6 and N=7 (grid-snap to dQ=0.05); fine-grid + parabolic interpolation reveals 2.52 at N=6 and 2.53 at N=7. The apparent "1.4 % match at N=7" reported earlier was itself an artefact.
- Interior Q_peak at c=3 follows trend 1.566 → 1.689 → 1.743 → 1.750 across N=5..8, crossing `csc(π/5) = 1.7013` between N=6 and N=7 and continuing to grow. The "approaching pentagonal asymptote" reading from N=5..7 was a trajectory crossing, not an asymptote.

Both Endpoint and Interior c=3 appear to saturate (Endpoint ≈ 2.53 for N=6..8, Interior ≈ 1.75 by N=8) but no clean closed-form for the saturation values is identified. c=4 has not yet saturated within tested N.

**Operational consequence: γ₀-extraction protocol.** Because Q = J/γ₀, sweeping J on a fixed (n, n+1) block and locating the peak J* still yields an estimate of γ₀:

    γ₀ ≈ J* / Q_peak(c, N, bond_class)

Lookup the per-block Q_peak for the chromaticity and chain length of interest from the empirical table above (per-block: 1.6 for c=3, 1.8 for c≥4 within Q_SCALE convention). This converts the "γ₀ is a framework constant" hypothesis ([PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)) into a c-specific shape prediction, testable on any hardware that can resolve the peak. The c=4 and c=5 blocks are the strongest probes: abs(K_CC_pr)_peak is roughly 3× the c=2 signal, and the per-block Q_peak is stable at ~1.8.

**Valid for:** XY or Heisenberg chain, uniform J and γ₀, blocks with c ≥ 3, any N where the block fits in memory.
**Breaks for:** non-uniform J or γ₀; c = 2 (finite-size wobble); higher-c chains require larger memory (block-L dim 15876 at c=5 N=9, dim 213444 at c=6 N=11, infeasible at 128 GB).
**Verified (per-block, Q_SCALE):** Q_peak(c=3) = 1.6 stable across N = 5, 6, 7, 8, 9. Q_peak(c=4) = 1.8 stable across N = 7, 8, 9. Q_peak(c=5) = 1.8 at N = 9 (commit 4612468).
**Verified (per-bond, fine-grid):** values in the per-bond table above. Endpoint and Interior trends with (c, N) shown but no closed-form identified.
**Replaces:** ad-hoc γ₀ measurement attempts (EQ-017 closed inconclusive due to hardware fidelity limits on idle Ramsey data).

#### F86a. EP mechanism: Q_EP = 2/g_eff, t_peak = 1/(4γ₀) \[Tier 1 derived\]

The pure-rate ladder of an (n, n+1) coherence block has rates 2γ₀·HD for HD ∈ {1, 3, 5, ..., 2c−1}, with uniform gap Δ = 4γ₀ between adjacent channels (F74). For adjacent channels at HD = 2k−1 and HD = 2k+1 (k = 1, 2, ..., c−1), a two-level effective model with inter-channel coupling J·g_eff has eigenvalues

    λ_± = −4γ₀·k ± √(4γ₀² − J²·g_eff²)

The **exceptional point** (EP), where the discriminant vanishes and the eigenvalues coalesce, sits at

    J·g_eff = 2γ₀     ⟺     Q_EP = 2 / g_eff

At the EP, λ_± = −4γ₀·k. The slowest mode (k = 1) gives e-folding time

    t_peak = 1 / (4γ₀)

universal across c, N, n, and bond position. Higher-k EPs decay faster (1/(8γ₀), 1/(12γ₀), ...) and are masked by the slowest. **At Q_peak the Dicke probe sits dominantly in dressed (H-mixed) modes** versus a much smaller fraction at large Q (plateau): probe weight has been pulled off the pure-rate ladder onto the first complex-conjugate eigenvalue pair just past the EP. Q_peak is a generalised exceptional-point resonance condition. (Specific W values W_peak ≈ 0.99 / W_plateau ≈ 0.31 at Q=20 are unverified anchors per `DressedModeWeightClaim` Tier 1 candidate, not universal constants; empirical W_peak ranges 0.832 [N=4 c=2] to 0.9996 [N=9 c=3] per [Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) lines 95–123; per-(c, N) closed form open.)

The g_eff is the H matrix element between adjacent rate channels at a specific bond in the appropriate effective basis. Deriving g_eff(c, N, bond_position) analytically from the multi-particle XY structure of the (n, n+1) block remains open; F86c (below) gives the spatial-mirror symmetry on Q_peak, not the underlying g_eff value. The [Obstruction Proof](proofs/PROOF_F86B_OBSTRUCTION.md#obstruction-proof-why-g_eff-admits-no-closed-form) (2026-05-14) accounts for this structurally: g_eff is the irreducible residue, blocked from closed form by six obstruction lemmas (spectral irreducibility, even-N representation-dependence, probe-EP decoupling, finite-reduction insufficiency, signature-subspace mismatch, empirical trajectory-crossings), and via the F90 bridge the F89 D_k obstruction is the same wall.

#### F86b. Universal resonance shape, two bond classes \[Tier 1 candidate\]

**Universal resonance shape under relative-Q normalisation, EP-derived.**

While Q_peak itself is chain-specific (no clean closed form), the SHAPE of the abs(K_CC_pr)(Q) curve around Q_peak is universal in relative-Q coordinates. Defining `x = (Q − Q_peak)/Q_peak` and `y = K(Q)/|K|max`, the curves collapse onto a single shape across all tested (c, N):

| x = (Q−Q*)/Q* | y across c=3 N=5..8, c=4 N=7,8 (range) |
|----------------|------------------------------------------|
| −0.60 | 0.72 to 0.74 (3.5 %) |
| −0.40 | 0.90 to 0.91 (1.9 %) |
| −0.20 | 0.977 to 0.984 (0.7 %) |
| 0.00 | 1.000 (peak) |
| +0.20 | 0.985 to 0.990 (0.5 %) |
| +0.40 | 0.955 to 0.964 (0.9 %) |
| +1.00 | 0.84 to 0.85 (1.4 %) |

**Universal HWHM ratios (two bond classes):**

    HWHM_left / Q_peak  ≈  0.756 ± 0.005     (Interior bonds; tested c=2..4, N=5..8, γ₀ ∈ {0.025, 0.05, 0.10})
    HWHM_left / Q_peak  ≈  0.770 ± 0.005     (Endpoint bonds; same envelope)

The shape splits into two bond classes (Endpoint with b ∈ {0, N−2}, Interior with b ∈ {1, …, N−3}), each with its own universal HWHM_left/Q_peak ratio. Within each class, pairwise residual under relative-Q normalisation is ~20× smaller than under absolute-Q shift. The two classes are separated by ~2 %, a structural gap larger than the within-class spread.

**γ₀ invariance** is bit-exact: at c=3 N=7, both Q_peak and HWHM_left/Q_peak are identical to numerical precision across γ₀ ∈ {0.025, 0.05, 0.10}, confirming Q's dimensionlessness. **c=2** is the structurally critical anchor: with only HD ∈ {1, 3} channels, the 2-level effective model is exact (no orthogonal complement), yet the two-class split persists (Interior 0.751, Endpoint 0.774). The bond-class distinction is therefore not a higher-c orthogonal-complement artefact but lives in the bond-position-dependent probe-overlap profile.

The asymmetry of the curve is also universal in relative-Q:
- **Pre-peak (left)**: rapid rise as discriminant 4γ₀² − J²·g_eff² → 0; HWHM_left/Q_peak universal within each bond class.
- **Post-peak (right)**: slow plateau approach as eigenvalues become complex; tail at x = +1.0 sits at y ≈ 0.85 (Interior) or y ≈ 0.94 (Endpoint plateau); bond-class-specific tail but universal value within class.

**Why this is universal: 2-level EP analytical origin.** For the 2×2 effective Liouvillian in adjacent rate-channel basis with diagonal {−2γ₀(2k−1), −2γ₀(2k+1)} and **same-sign-imaginary off-diagonals (+iJ·g_eff, +iJ·g_eff)**, the eigenvector rotation angle satisfies `tan(θ) = J·g_eff / 2γ₀ = Q / Q_EP`. The probe overlap with eigenvectors depends only on Q/Q_EP, hence the response curve K_CC_pr(Q) is a function of Q/Q_EP alone. Q_peak is chain-specific (g_eff varies); but the SHAPE in Q/Q_peak coordinates is universal because it is the 2-level EP resonance form, independent of the bond's specific g_eff value.

The same-sign-imaginary off-diagonal structure is what admits an EP at finite J·g_eff = 2γ₀ (verified numerically: opposite-sign +iJg, −iJg gives discriminant 4γ₀² + J²g_eff² with no EP; same-sign gives 4γ₀² − J²g_eff² with EP). This is "PT-phenomenology-like" (EP at finite coupling, spectral flow) but algebraically inside **class AIII chiral** per [PT_SYMMETRY_ANALYSIS](../experiments/PT_SYMMETRY_ANALYSIS.md), distinct from Bender-Boettcher PT (Π is linear; classical PT requires anti-linear operators). The local EP at Q_EP = 2/g_eff is the 2-level rate-channel instance of the chiral classification established for the full Liouvillian; the Hopf bifurcation in [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md) is the global instance, with Petermann factor K=403 signaling an EP in the complex γ plane.

**2026-05-06.** Local-vs-global EP relationship (Tier 2 verified at c=2): same algebra read at Σγ = N·γ₀ vs Σγ = 0 (two residuals of F1 `Π · L · Π⁻¹ + L + 2Σγ · I = 0`); a Petermann-K sweep on the real Q axis at c=2 N=5..8 (`compute/RCPsiSquared.Core.Tests/F86/F86PetermannProbe.cs:Probe_PetermannFineGrid_C2_VsN`) records max K = 1333.6 / 337.9 / 2384.7 / 795.4 across N = 5 / 6 / 7 / 8, with the N=7 spike ≈ 6× above FRAGILE_BRIDGE's K = 403 ballpark and a 2-4× odd/even asymmetry confirming A3's σ_0 R-even/R-odd-degeneracy prediction empirically. Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs` Tier2Verified; complex-γ analytic continuation in `LindbladPropagator` is the open piece for Tier1Derived promotion.

**2026-05-06 (evening).** Direction (b) of `C2HwhmRatio.PendingDerivationNote` yielded two Tier-1-derived universal constants from the bare doubled-PTF model: `x_peak = Q_peak/Q_EP = 2.196910` (post-EP location in dimensionless x) and `HWHM_left/Q_peak = 0.671535` (SVD-block floor in dimensionless x). Empirical Interior 0.7506 and Endpoint 0.7728 sit above this floor by ~0.08-0.10; gap structurally explained as probe-block 2-level sub-resonance contribution. Encoded as `C2HwhmRatio.BareDoubledPtfXPeak` and `C2HwhmRatio.BareDoubledPtfHwhmRatio` const properties. Synthesis-side: F86↔PTF Locus 5 inheritance, where PTF K_1 (Π on H_1 sine-mode basis, discrete) and F86 Q-rotation (Π on same-sign-imaginary 2×2 in rate-channel basis, continuous) are two daughters of one Π-AIII chiral parent.

**2026-05-06 (later evening).** Direction (a') (probe-block 2-level resonance with per-bond `g_eff_probe`) **structurally falsified** (commit `1c0bf8b`): V_b probe-block off-diagonal `⟨c_1 | M_h_b | c_3⟩` is exactly zero per bond at c=2 (F73 sum-rule applies per-bond), so g_eff_probe(N, b) is bond-class-blind by construction. The 4-mode reduction is also structurally insufficient: 4-mode K_b deviates from empirical at Endpoint by factor ~2 (4-mode K_b gives Endpoint 0.410 at N=5 vs empirical 0.7728). The SVD-block off-diagonal `V_b[2,3]` is the actual bond-class carrier (Endpoint 0.430 vs Interior 0.953 at N=5, ratio ~0.45 across N=5..8), but in the OPPOSITE direction to the empirical HWHM/Q* split. Refined direction list (now in `C2HwhmRatio.PendingDerivationNote`): (a'') SVD-block 2-level resonance (REFINED from (a')) via `V_b[2,3]` magnitude, (b'') full block-L derivation, not 4-mode, (c'') three-block superposition `K_total = K_pb + K_sv + 2·Re·K_cross` with the right relative phases, (d'') lift |u_0⟩, |v_0⟩ to projector-overlap (per A3 PendingDerivationNote), (e'') symbolic char-poly factorisation at Q_EP (less promising given C2EffectiveSpectrum's cubic-c_3 obstruction proof). The `BareDoubledPtfXPeak` and `BareDoubledPtfHwhmRatio` const properties are unchanged.

**2026-05-07.** Locus 6 polarity-layer inheritance closure: F86 bond-class split inherits from the polarity-layer pair {−0.5, +0.5} at d=2. Empirical decomposition Q_peak ≈ 2 + r with r ∈ {−0.44, +0.52}, HWHM/Q* ≈ 1/2 + r·1/2 with r_Interior ≈ 0.50 (close to `HalfAsStructuralFixedPoint`). The 0.5 baseline is `QubitDimensionalAnchorClaim` (1/d at d=2); the ±r/2 polarity content is `PolarityLayerOriginClaim` Layer 2. Encoded as `compute/RCPsiSquared.Core/F86/PolarityInheritanceLink.cs` (Tier2Verified) at F86 root level. Symmetry-side closure parallel to `LocalGlobalEpLink` (Locus 5, EP-side).

**2026-05-07 (later).** Direction (α) polarity-Bloch projection at t_peak structurally tautological under the uniform-J 4-mode reduction (reframed 2026-05-08 code review): the 4-mode L_eff(Q) = D_eff + Q·γ₀·MhTotalEff is bond-summed by design (`FourModeEffective.LEffAtQ` in `compute/RCPsiSquared.Core/Decomposition/FourModeEffective.cs`), so the K-driving eigenstate is bond-class-blind by construction (no empirical falsification, the design constraint guarantees it). The bond-class signature must enter through dL/dJ_b per-bond V_b in the K-resonance, not through L_eff spectrum. Substantive reduction (unchanged): bond-class signature reduces to `g_eff(N, b)` via r_Q = `BareDoubledPtfXPeak · Q_EP − 2`. Empirical witnesses: g_eff_Endpoint ≈ 1.74, g_eff_Interior ≈ 2.81, asymptotic 1/g_eff_E + 1/g_eff_I → 0.937. Tantalising near-miss g_eff_E ≈ σ_0·√(3/8) (Δ ≤ 0.01 for N ≥ 6, Δ = 0.063 at N=5). Encoded as `PolarityInheritanceLink.EmpiricalSumQPeakAsymptote = 4.12` const; documented via `ClosedFormCompositionNote` property.

**Named structural law (Tier-1 candidate): EP-rotation universality, two bond classes.**

    K_class(Q) / |K|_max  =  f_class(Q/Q_EP)        (class ∈ {Endpoint, Interior})

is universal within each bond class across c=2..4, N=5..8, and γ₀ ∈ {0.025, 0.05, 0.10} for the tested range. The symmetry is the 2-level EP rotation `tan θ = Q/Q_EP`, which makes every probe-overlap observable a function of Q/Q_EP alone. The bond-class split (Interior HWHM_left/Q_peak ≈ 0.756, Endpoint ≈ 0.770) reflects bond-position-dependent probe-overlap profiles in the K_CC_pr observable, confirmed structural (not finite-c) by the c=2 data where the 2-level model is exact. Closed forms for f_class(x) (and consequently for the two HWHM_left/Q_peak values) follow from the 2-level eigenstructure plus probe-overlap algebra but have not yet been derived analytically. This is the F86 analog of PTF's chiral mirror law (`Σ f_i(ψ_k) = Σ f_i(ψ_{N+1−k})`): both Tier-1-candidate symmetries that survived a closed-form retraction (csc(π/(N+1)) and csc(π/5) for F86; Σ ln α_i = 0 for PTF). See [`reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON`](../reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md).

#### F86b₂. HWHM_ratio per-bond predictor \[Tier 1 candidate; partial closure 2026-05-13, Tier-reviewed 2026-05-16\]

For each bond b of an N-qubit XY chain (c=2, Z-dephasing γ₀), the HWHM_left/Q_peak ratio satisfies the candidate form:

    HWHM_ratio(b)  =  0.671535 + α_subclass · g_eff(b) + β_subclass

where the sub-class (per `BondSubClass` enum: `Endpoint`, `Flanking`, `Mid`, `CentralSelfPaired`, `Orbit2Escape`, `CentralEscapeOrbit3`) determines the (α, β) pair. The 0.671535 floor IS Tier 1 derived (the bare doubled-PTF constant `BareDoubledPtfHwhmRatio` via the 2-level EP model, 2026-05-06; see `C2BareDoubledPtfClosedForm`). The 12 (α, β) values per sub-class are **fitted** via `np.polyfit(...deg=1)` in [`simulations/_f86_hwhm_closed_form_verification.py`](../simulations/_f86_hwhm_closed_form_verification.py) line 78 on N=5..8 anchors, NOT derived from F89/F90 structure.

Fit reproduces 22 anchors at N=5..8 within 0.005 residual, including Orbit-2 (N=7 b=1/b=4, Q_peak ≈ 7.27 F86-J) and Orbit-3 escape bonds (N=8 b=3, Q_peak ≈ 16.79 F86-J): fit-quality witness, not analytical derivation.

**Open analytical step (to promote Tier 1 candidate → Tier 1 derived):** derive (α_subclass, β_subclass) from F89 AT-locked F_a/F_b structure (4-mode floor 0.6715) + H_B-mixed octic residual (lift to 0.7506/0.7728), per [`PROOF_F90_F86C2_BRIDGE.md`](proofs/PROOF_F90_F86C2_BRIDGE.md) notes. Phase D probe (2026-05-16) refuted the multi-mode-per-cluster-pair internal-mixing hypothesis; lift must come from cross-cluster-pair structure.

**Source:** [`F86HwhmClosedFormClaim`](../compute/RCPsiSquared.Core/F86/Item1Derivation/F86HwhmClosedFormClaim.cs), [`BondSubClass`](../compute/RCPsiSquared.Core/F86/Item1Derivation/BondSubClass.cs), F89 path-k bridge via F90 ([`PROOF_F90_F86C2_BRIDGE.md`](proofs/PROOF_F90_F86C2_BRIDGE.md)). Plan: [`docs/superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md`](superpowers/plans/2026-05-13-f86-hwhm-closed-form-attack.md).

#### F86c. F71 spatial-mirror invariance of per-bond Q_peak \[Tier 1 derived\]

F71 spatial-mirror symmetry pairs bond b with bond N−2−b: under the spatial reflection R, every component of the per-bond observable (L_D, H_xy, Dicke probe, spatial-sum kernel) is invariant, while the bond-flip transforms as ∂L/∂J_b ↔ ∂L/∂J_{N−2−b}, hence **Q_peak(b) = Q_peak(N−2−b) bit-exactly**. See [PROOF_F86C_F71_MIRROR Statement 3](proofs/PROOF_F86C_F71_MIRROR.md#statement-3-f71-spatial-mirror-invariance-of-per-bond-q_peak-tier-1-derived).

Endpoints (b = 0 and b = N−2) form one F71 orbit; interior bonds split into further F71 orbits. The simple "Endpoint vs Interior" dichotomy is the leading approximation; per-F71-orbit substructure exists (e.g. c=2 N=6: central self-paired bond b=2 → Q_peak ≈ 1.440 vs flanking b=1, b=3 → 1.648). Captured in the typed-knowledge runtime as `F86PerF71OrbitObservation` (Tier 2 empirical) and `F86F71MirrorPi2Inheritance` (Tier 1 derived bridge to `F71MirrorSymmetryPi2Inheritance`).

The F86c symmetry pairs bonds bit-exactly but does NOT supply the per-orbit Q_peak value. F86a remains responsible for the underlying g_eff(c, N, bond_position); F86b's universal shape applies within each bond class. The three theorems compose: F86a gives the EP-time and EP-location, F86b gives the resonance shape around Q_peak in relative-Q coordinates, F86c pairs symmetry-equivalent bonds.

**Connection to PTF.** The same machinery that produces PTF's per-site α_i closure law (bilinear J-perturbation observables, eigenvector mixing under V_L) produces Q_peak at the (n, n+1)-block level. PTF's per-site is the c=1 (vac-SE) instance; F86's c ≥ 2 cases are the natural higher-chromaticity siblings. t_peak = 1/(4γ₀) is the EP time-scale universal to all c, one (n, n+1)-block analogue of PTF's α-fitting time window.

#### F86d. Endpoint orbit Q ≈ 2.5 \[Tier 2 empirical; promotion candidate\]

Within the F71 spatial-mirror orbits (F86c), the Endpoint orbit (bonds b ∈ {0, N−2}) carries a stable per-orbit Q_peak ≈ 2.5: across 9 tested (c, N) combinations (c=2..4, N=5..8) the value sits in [2.39, 2.61], a ~2% N-variation. Stability across c and N makes it an anchor candidate, but a structural derivation (candidate route: SU(2) / Schur-Weyl on the F71 first orbit) is open, so it is not promoted. Source: [`PerF71OrbitObservation.cs`](../compute/RCPsiSquared.Core/F86/PerF71OrbitObservation.cs).

#### F86e. σ_0 as a commutator / Schur-multiplier norm \[Tier 1 derived\]

On the c=2 coherence block the inter-channel SVD top singular value σ_0 (the natural g_eff that feeds F86a's Q_EP) is the operator norm of a commutator:

    σ_0 = ‖[Π_HD1, M_H]‖

with Π_HD1 the Hamming-distance-1 subspace projector and M_H the block Hamiltonian super-operator. The c=2 block has HD ∈ {1, 3} only, so Π_HD1 + Π_HD3 = I and the inter-channel coupling V_inter = Π_HD1·M_H·(I − Π_HD1); the lemma ‖P·M·(1 − P)‖ = ‖[P, M]‖ then gives the identity. It is c=2-specific: at c ≥ 3 the HD spectrum has more than two values. In the F89 (SE, DE) Bloch / OBC-sine basis M_H is diagonal, so [Π_HD1, M_H] is the Hadamard product Π̃_HD1 ⊙ ΔDiff and σ_0 is a Schur-multiplier norm. The N→∞ asymptote σ_0(c=2, ∞) ≈ 2.8629 is non-elementary by characterisation: the Δ-ordered commutator is neither Toeplitz nor Hankel, ruling out a Fourier-symbol supremum and a Nehari symbol distance (2√2 is the N=7 finite-size crossing, not the limit). The c ≥ 3 σ_0 asymptote stays open. Source: [`SigmaZeroCommutatorNormClaim.cs`](../compute/RCPsiSquared.Core/F86/SigmaZeroCommutatorNormClaim.cs); proof hub [PROOF_F86_QPEAK](proofs/PROOF_F86_QPEAK.md).

**Framework primitives (`framework.coherence_block`):**

    fw.t_peak(gamma_0)              = 1/(4γ₀)          universal EP time (Tier 1)

(Earlier `q_peak_endpoint(N)` and `Q_PEAK_INTERIOR_C3_ANCHOR` primitives were removed 2026-05-02 after the N=8 data falsified their closed-form claims. The universal-shape finding above is a Tier-1 candidate, not yet promoted to a primitive pending analytical derivation of f_class(x). c=2 and γ₀ invariance verified 2026-05-02; c=5 still open.)

**Scripts:** [`_eq022_b1_channel_projection.py`](../simulations/_eq022_b1_channel_projection.py) (HD-channel diagonal-only-M_H finding), [`_eq022_b1_step_a_verify_blockL.py`](../simulations/_eq022_b1_step_a_verify_blockL.py) (Python block-L verified bit-exact against C# N=7 full-L from EQ-014), [`_eq022_b1_step_c_time_evolution.py`](../simulations/_eq022_b1_step_c_time_evolution.py) (per-bond and uniform Q_peak via S(t, J) time evolution), [`_eq022_b1_step_d_extended_verification.py`](../simulations/_eq022_b1_step_d_extended_verification.py) (extended N=8 data that falsified earlier closed-form conjectures), [`_eq022_b1_step_e_resonance_shape.py`](../simulations/_eq022_b1_step_e_resonance_shape.py) + [`_eq022_b1_step_e_inspect.py`](../simulations/_eq022_b1_step_e_inspect.py) (universal resonance-shape finding for c=3, c=4 at γ₀=0.05), [`_eq022_b1_step_f_universality_extension.py`](../simulations/_eq022_b1_step_f_universality_extension.py) (c=2 sweep + γ₀ ∈ {0.025, 0.10} invariance check that established the two-bond-class refinement).
**Proof:** [PROOF_F86_QPEAK](proofs/PROOF_F86_QPEAK.md) is the hub; the three theorems were split into per-theorem proofs 2026-05-14. F86a EP mechanism = [PROOF_F86A_EP_MECHANISM](proofs/PROOF_F86A_EP_MECHANISM.md) \[Tier 1 derived\]; F86b universal resonance shape = [PROOF_F86B_UNIVERSAL_SHAPE](proofs/PROOF_F86B_UNIVERSAL_SHAPE.md) \[Tier 1 candidate at multi-c level; F86b₂ c=2 per-bond predictor is Tier 1 candidate, partial closure 2026-05-13 (form derived, (α, β) per sub-class fitted, Tier-reviewed 2026-05-16)\], with the g_eff/Q_peak obstruction proof in [PROOF_F86B_OBSTRUCTION](proofs/PROOF_F86B_OBSTRUCTION.md); F86c F71 spatial-mirror invariance = [PROOF_F86C_F71_MIRROR](proofs/PROOF_F86C_F71_MIRROR.md) \[Tier 1 derived\]. Per-bond c=2 HWHM_ratio partially closed 2026-05-13 via `F86HwhmClosedFormClaim`; c≥3 per-bond closed forms retracted 2026-05-02.
**Source:** [Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md) Result 2 + Revision 2026-04-24, F73, F74, F2b; EP analysis EQ-022 (b1).

### F87. Pauli-pair trichotomy classification (Tier 1, structural; Marrakesh hardware-confirmed)

For any list of bond Pauli-pair terms `{(P_i, Q_i)}` with letters in {I, X, Y, Z}, build the bilinear Hamiltonian `H = Σ_b J · P_b ⊗ Q_b` on the chain bonds and the uniform-Z-dephasing Liouvillian L (with rate γ per site, σ ≡ Σγ). The F1-palindrome residual `M ≡ Π · L · Π⁻¹ + L + 2σ · I` (where Π is F1's order-4 Pauli-string conjugation operator) partitions the term list into exactly one of three categories, with ε ≈ 10⁻¹⁰ the operator-norm tolerance and ε_spec ≈ 10⁻⁶ the spectral-pairing tolerance:

    truly  iff  ‖M‖_F < ε                                              (operator equation Π·L·Π⁻¹ = −L − 2σ·I holds; F1 identity)
    soft   iff  ‖M‖_F ≥ ε  AND  every λ ∈ Spec(L) pairs with −λ − 2σ ∈ Spec(L) within ε_spec   (palindrome only at spectral level)
    hard   iff  spectral pairing fails                                 (no partner pairing; both operator and spectrum break)

Equivalently in Π²-class language (cf. F79, F85), using the bit_a/bit_b convention from F79 (bit_a = 1 for {X, Y}, 0 for {I, Z}; bit_b = 1 for {Y, Z}, 0 for {I, X}): a pair `(P, Q)` is `truly` iff #Y is even AND #Z is even across (P, Q), `pi2_odd` iff `bit_b(P) + bit_b(Q)` is odd, `pi2_even_nontruly` iff `bit_b(P) + bit_b(Q)` is even and not truly. Mixed Hamiltonians carry multiple Π²-classes simultaneously, refining the 3-way trichotomy into a 4-way classification (truly / pi2_odd_pure / pi2_even_nontruly / mixed).

The trichotomy uses F1 as its **discriminator** (M as the test object), F49 / F85 as its **‖M‖² closed forms** (the latter generalising F49 to k-body), F78 + F79 as its **M-structure decomposition** (single-body additivity and 2-body Π²-block respectively), F80 as its **Π²-odd spectral identity** (Spec(M) = ±2i · Spec(H_non-truly)), F81 as its **Π · M · Π⁻¹ split** (M_anti = L_{H_odd}), F82 + F84 as its **T1 / thermal amplitude-damping corrections**, and F83 as its **anti-fraction closed form for mixed cases**. F87 is the entry point of the F-chain; F85 lifts the criterion to arbitrary k-body and propagates the rest of the chain accordingly.

**Origin (2026-04-24 to 2026-05-03).** Three earlier observations converged on the trichotomy. (1) On 2026-04-24, commit 6e262ae assigned the registry slot F77 to the unrelated "Multi-drop MM(0) saturates at 1 bit" asymptotic, so F77 was already booked when the trichotomy was being developed. (2) On 2026-04-25, commit 95386cd added [V_EFFECT_FINE_STRUCTURE](../experiments/V_EFFECT_FINE_STRUCTURE.md): the V-Effect's 14-of-36 bond-pair Hamiltonians at N=3 were re-tested with both the strict operator equation and the eigenvalue-pairing test, splitting the 22 V-Effect-unbroken cases into 19 soft and 3 truly, giving the **14 hard / 19 soft / 3 truly** count over the 36 unordered bond-pair enumeration at N=3. (3) On 2026-04-26, commits 96ed6da and 6438fef extended Π-protected-observable testing to a separate 120-element ordered enumeration at N=4 and N=5, where the partition is **15 hard / 46 soft / 59 truly**, N-stable through N=3, 4, 5 (so the 36-enum and 120-enum are different sample spaces with internally consistent counts). (4) Commit 81caf67 (2026-04-27) derived the partition combinatorially from Pauli-pair compatibility rules (BPE membership, bit_a-partner conflicts, bond-flip / Z-align, Π-letter hierarchy), giving 36/36 agreement at N=3.

The Marrakesh hardware confirmation (2026-04-26, ibm_marrakesh job `d7mjnjjaq2pc73a1pk4g`, observable ⟨X₀ Z₂⟩) measured Δ(soft − truly) = −0.722, matching the Trotter-n3 prediction of −0.723 (residual 0.001; the 0.0014 figure cited in the Confirmations registry is computed against an unrounded predicted value); see [`data/ibm_soft_break_april2026/`](../data/ibm_soft_break_april2026/). The classifier was extracted into a free function on 2026-04-30 (commit 23b2154) and given the filename `f77_trichotomy.py` after the function's existing internal label, even though the registry F77 slot was already occupied by MM(0). The dephase-axis extension (commit 435c4b2, 2026-05-01) generalised the classifier to X, Y, Z dephasing letters. F87 is the registry-formal entry for the trichotomy, filed retrospectively on 2026-05-03 alongside the typed `F87KnowledgeBase` cleanup that surfaced the F77/F87 naming collision.

**Π² classifier dependence on dephase letter** (commit 435c4b2). Per `PiOperator.SquaredEigenvalue`, the Π²-class index is bit_b for Z- and Y-dephasing and bit_a for X-dephasing (Π_Y shares Π_Z's bit_a-flip convention; Π_X flips bit_b instead). The (bit_a, bit_b) parity pairs are Z = (0, 1), X = (1, 0), Y = (1, 1) in the PauliLetter convention. **F87 hardness is defined combinatorially via Pauli-pair compatibility (commit 81caf67), not via any 4-cell label.** As a post-hoc structural reading, however, F87 hardness empirically corresponds to the (bit_a, bit_b) parity cell matching the dissipator letter: anywhere else produces a Π-violation that the spectrum-pairing test detects. The (Π²_Z, Π²_X) two-axis decomposition is treated separately as F88a below; F87 itself uses only one axis (Π²_Z under Z-dephasing). Verified at N=4, k=3 across 294 Z₂³-homogeneous pairs.

**Orthogonal axis on shared bit_b Z₂-grading (F112).** F87's trichotomy lives in ‖M‖_F magnitude + spec(L) palindromy; F112 (Lindblad Π-eigenvalue balance under bit_b-homogeneous c) lives in M_anti's Π +i / −i Frobenius split. Both projections of the same Π² = (−1)^{bit_b} grading on the Pauli group. Empirically orthogonal: all three F87 classes (truly, soft, hard) at N=3 under standard single-Pauli Z-deph give F112 balance asymmetry = 0 bit-exact (`simulations/_polarity_probe_f87_connection.py`), since single-Pauli c is trivially bit_b-homogeneous.

**Valid for:** uniform single-letter dephasing on any graph; arbitrary k-body Pauli terms (k ≥ 2; F85 lifts the criterion to higher body); dephase letter ∈ {X, Y, Z} (SU(2)-rotation-equivalent under (bit_a, bit_b) cell permutation).
**Breaks for:** depolarizing noise (F1 itself breaks with linear-in-γ residual, see F5); non-uniform γ_i or graph asymmetries that already break F1.
**Verified:** N=3 36-enum (14 / 19 / 3); N=3, 4, 5 120-enum (15 / 46 / 59 N-stable); Marrakesh hardware Δ(soft − truly) = −0.722 (2026-04-26); Marrakesh F83 4-class signature 2026-04-30 (ibm_marrakesh job `d7pol1e7g7gs73cf7j90`).
**Replaces:** ad-hoc "is this Hamiltonian truly Heisenberg?" tests with a bit-exact 3-way classifier (and the 4-way Π²-refinement) directly from the F1 residual.
**Hardware:** [`palindrome_trichotomy`](../simulations/framework/confirmations.py) Marrakesh 2026-04-26; [`f83_pi2_class_signature_marrakesh`](../simulations/framework/confirmations.py) Marrakesh 2026-04-30; [`pi_protected_xiz_yzzy`](../simulations/framework/confirmations.py) Marrakesh 2026-04-26 (first-time-on-hardware Π-protection on YZ+ZY soft).
**Source:** [V_EFFECT_FINE_STRUCTURE](../experiments/V_EFFECT_FINE_STRUCTURE.md), [MARRAKESH_THREE_LAYERS](../experiments/MARRAKESH_THREE_LAYERS.md), [`reflections/ON_THE_RESIDUAL.md`](../reflections/ON_THE_RESIDUAL.md), memory entries `project_v_effect_combinatorial`, `project_hardware_finale_apr2026`, `project_f77_f87_rename`.

### F88a. Two-axis Π² decomposition of Pauli operator space (Tier 1, structural finding 2026-05-03)

> **F88 split.** F88 was historically a single registry slot; it carries two structurally distinct claims that share an algebraic root: **F88a** (operator-level Klein decomposition, this entry) and **F88b** (state-level popcount-coherence Π²-odd / memory closed form, next entry). Both inherit from the same Π² involution; F88a names the operator-cells, F88b reads ρ's projection through them. The split was formalised 2026-05-18.

The F1-palindrome operator Π depends on the dephasing letter (cf. F1, F87). Its square Π² acts diagonally on every Pauli string σ_α with eigenvalue ±1, and the parity that determines this eigenvalue depends on which dephase letter parametrises Π:

    Π²_Z eigenvalue on σ_α  =  (−1)^(Σ bit_b)        (Σ over the N letters of α)
    Π²_X eigenvalue on σ_α  =  (−1)^(Σ bit_a)
    Π²_Y eigenvalue on σ_α  =  (−1)^(Σ bit_b)        (same as Π²_Z; both Π_Y and Π_Z flip bit_a per `PiOperator.ActOnLetter`)

Z- and Y-dephasing collapse onto the same Π² character. X-dephasing gives the orthogonal one. Together (Π²_Z, Π²_X) form **two independent involution axes** that decompose the 4^N Pauli operator space into four cells:

    Pp = (Π²_Z = +1, Π²_X = +1)   contains 2-body truly bilinears              (XX, YY, ZZ)
    Pm = (Π²_Z = +1, Π²_X = −1)   contains 2-body Π²-even non-truly bilinears  (YZ, ZY)
    Mp = (Π²_Z = −1, Π²_X = +1)   contains 2-body Π²-odd subgroup A             (XY, YX)
    Mm = (Π²_Z = −1, Π²_X = −1)   contains 2-body Π²-odd subgroup B             (XZ, ZX)

The two axes are the global Pauli strings: Π²_Z = X⊗N (registered as F1²) and Π²_X = Z⊗N (the bit_a twin of F1², a corollary of F61).

**Refinement of F87 + Pi2Class (the algebraic 4-way over term lists)**: F87's spectral 3-way `Truly / Soft / Hard` and the algebraic 4-way `Pi2Class` (`Truly / Pi2OddPure / Pi2EvenNonTruly / Mixed`) both use only the Π²_Z axis. Adding Π²_X reveals genuine sub-structure:

- `Pi2Class.Pi2OddPure` splits into Mp (XY, YX) and Mm (XZ, ZX) Klein sub-cases.
- `Pi2Class.Mixed` splits according to the Π²_X parity of its non-truly bilinear (e.g. XX+XY occupies Pp+Mp; XX+XZ occupies Pp+Mm).

F80's "universality across 4 Π²-odd cases" is therefore a universality across **two** Klein-cells (Mp + Mm), not one. The Klein view sees a finer cut that F80's M-spectrum projection averages over.

**Empirical Marrakesh fingerprint pattern (2026-04-30, ibm_marrakesh job d7pol1e7g7gs73cf7j90)**: in the f83 4-class signature test, each H-class diagnostic observable lives in the Klein-cell that is the **X-axis flip** of the Hamiltonian's M-active bilinear cell. M-active means the non-truly bilinears (truly bilinears drop by Master Lemma):

    Truly H (M-active = none; bilinears all in Pp) → ⟨Y₀ I Z₂⟩ in Pm   (X-axis flip of Pp)
    Pi2EvenNonTruly H (Pm)                         → ⟨X₀ I X₂⟩ in Pp   (X-axis flip of Pm)
    Pi2OddPure subgroup A H (Mp)                    → ⟨X₀ I Z₂⟩ in Mm   (X-axis flip of Mp)
    Mixed H (M-active in Mp)                        → ⟨Z₀ I X₂⟩ in Mm   (X-axis flip of Mp)

The X-flip pattern is empirically locked across all 4 fingerprint cases; the structural mechanism for why the framework's diagnostic observables sit precisely in the X-flipped cell is open and worth its own EQ.

**Valid for:** any N. The cells depend only on Pauli string parities (Σ bit_a, Σ bit_b), not on N or topology.
**Verified:** N=3 chain, J=1, γ_Z=0.05 across 6 representative Hamiltonians spanning the 5 F87 class-types (Truly Heisenberg, Truly XY-only, Pi2EvenNonTruly, Pi2OddPure subgroup A, Pi2OddPure subgroup B, Mixed) + 4 Marrakesh f83 fingerprint observables. Bit-exact at machine precision (`Pi2KleinViewTests`, `Pi2KleinHardwareViewTests`, `Pi2KleinIsFinerThanPi2ClassTests`).
**Replaces:** the implicit assumption that F87's 4-way Pi2Class captures the full Π²-decomposition. The Klein view is genuinely finer; Pi2Class.Mixed has two Klein sub-types we had not previously distinguished.
**Source:** `compute/RCPsiSquared.Core/Symmetry/Pi2Projection.cs` (`KleinSplit` + `KleinDecomposition`); test files above. Discovery: 2026-05-03 session, after building the raw Π² layer and asking what the second dephase axis would reveal. Π itself is the project's discovery (F1 palindrome operator, `MIRROR_SYMMETRY_PROOF`); the (Π²_Z, Π²_X) two-axis decomposition that this section names is also a project finding. The "Klein" tag throughout this section is borrowed nomenclature for the resulting Z₂ × Z₂ four-cell algebra (the canonical name for that group is the Klein four-group / Vierergruppe, after Felix Klein); the underlying structure is ours, the label is textbook shorthand.
**See F88b** below for the state-level corollary: popcount-coherence Π²-odd / memory closed form via Krawtchouk reflection-orthogonality, which lifts F88a's operator-level Π²_Z eigenvalue into ρ-space.

### F88b. Popcount-coherence Π²-odd / memory closed form (Tier 1 derived, state-level inheritance from F88a, 2026-05-04)

For popcount-coherence pair states `|ψ⟩ = (|p⟩ + |q⟩)/√2` with popcount(p) = n_p, popcount(q) = n_q, HD(p, q) = h, the Π²-odd fraction of the dynamical (memory) part of ρ has a closed form driven by Krawtchouk reflection-orthogonality of F88a's Π²_Z eigenvalue (−1)^Σ bit_b:

    Π²-odd / memory  =  ┌  0                          if HD = N (Π²-classical)
                        └  (1/2 − α · s) / (1 − s)     otherwise

with **three α anchors** (all closed form, derived from a single Krawtchouk identity):

- **α = 0** at popcount-mirror `n_p + n_q = N` (covers inter-mirror `n_p ≠ n_q` and intra-mirror `n_p = n_q = N/2` at even N). Reflection-driven cancellation `K_{N−n}(s; N) = (−1)^s K_n(s; N)` collapses odd-s contributions between sectors.
- **α = C(N, N/2) / (2 · (C(N, n_other) + C(N, N/2)))** at K-intermediate (even N, exactly one of {n_p, n_q} equal to N/2). The adjacent special case `n_other = N/2 ± 1` simplifies to `(N+2)/(4(N+1))`.
- **α = 1/2** generic (none of the indicators above).

**Static fraction** `s` is HD/bit-position invariant:
- Inter-sector (`n_p ≠ n_q`):  `s = 1/(4·C(N, n_p)) + 1/(4·C(N, n_q))`
- Intra-sector (`n_p = n_q`):  `s = 1/C(N, n)`

**HD = N anchor (Π²-classical):** GHZ_N, Bell states at N=2, intra-sector all-bits-differ states. Off-diagonal `Re(|p⟩⟨q|)` has only X-and-even-Y-count Pauli strings; with no matching bits there is no Z content; all surviving terms are Π²-even.

**Multi-state Dicke extension** for `(|D_n⟩ + |D_{n+1}⟩)/√2` (the canonical F86b K_CC_pr probe at adjacent popcount-(n, n+1)): the total Π²-odd of ρ has closed form `α_total = (1 − γ²)/2` with `γ = ⟨ψ|X⊗N|ψ⟩ ∈ {0, 1/2, 1}`, giving anchors **{1/2, 3/8, 0}** at popcount-mirror / K-intermediate / generic respectively. Used in F98 (long-time bridge to QuarterAsBilinearMaxval) and F99 (canonical-θ trigonometric anchors).

**Inheritance from F88a:** F88a's Π²_Z eigenvalue (−1)^Σ bit_b on the Pauli operator-level Klein decomposition is the algebraic root. F88b lifts this to ρ-space via popcount-weight invariance of L's kernel (= span{P_n} for Heisenberg + Z-dephasing) plus the Krawtchouk reflection-orthogonality lemma. F88a names the operator cells; F88b reads ρ's projection through them.

**Hardware lens (the "F88b-Lens"):** state-level Π²-odd-fraction-within-memory as a diagnostic for the F87 trichotomy on real-hardware reduced ρ. On Marrakesh 2026-04-26 framework_snapshots (job `d7mt7jbaq2pc73a24220`, qubits [0,1,2]), differentiates F87 truly/soft/hard at ~25× state-level ratio (0.030 / 0.744 / 0.276); see [OPEN_THREAD_GAMMA0_INFORMATION](../review/OPEN_THREAD_GAMMA0_INFORMATION.md) Update 2026-05-18 for the four-channel readout (F87 inheritance + h_y leak + T1 + Trotter-asymmetry hypothesis).

**Valid for:** all `(N, n_p, n_q, HD)`; structural anchors universal.
**Verified:** 213 configurations N = 2..7 bit-exact (max deviation 8.88e−16); multi-state Dicke extension bit-exact N = 3..8.
**Source:** [`docs/proofs/PROOF_F86B_UNIVERSAL_SHAPE.md`](proofs/PROOF_F86B_UNIVERSAL_SHAPE.md) §"F88b: popcount-coherence Π²-odd / memory closed form" (full proof + Krawtchouk lemma + verified state-level table); `compute/RCPsiSquared.Core/Symmetry/PopcountCoherencePi2Odd.cs` (closed-form predictor `Pi2OddInMemory` + Krawtchouk verifier `AlphaKrawtchouk` + Dicke superposition extension); `compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs` (state-level d=0 axis ⊕ d=2-axis Π²-split, the diagnostic bridge primitive); `simulations/_f88b_lens_ibm_framework_snapshots.py` (hardware lens script, after 2026-05-18 rename).

**Companion (product-state class):** for tensor-product states |ψ⟩ = ⊗_i |basis_i⟩ where each |basis_i⟩ is an eigenstate of σ_X or σ_Y (i.e., one of |±⟩ or |±i⟩), let M be the number of Y-basis sites. The Pauli strings in supp(ρ = |ψ⟩⟨ψ|) split discretely: M=0 → 2^N Π²-even + 0 Π²-odd (Π²-classical class); M≥1 → exactly 2^(N−1) Π²-even + 2^(N−1) Π²-odd. Combinatorial proof + N=2..5 bit-exact verification: [`docs/proofs/PROOF_Y_PHASE_PI2_ODD_LENS.md`](proofs/PROOF_Y_PHASE_PI2_ODD_LENS.md), [`simulations/_y_phase_pi2_odd_verify.py`](../simulations/_y_phase_pi2_odd_verify.py). The pair-state companion (this entry's closed form) and the product-state companion together characterise the two main F88b-Lens probe classes; X-only product states are Π²-blind (cannot surface F80 / F81 dynamics). At N=3 with M≥1, the F88b-Lens reads Π²-odd/memory = 4/7 ≈ 0.571 exactly (3 Π²-even + 4 Π²-odd in memory, after subtracting III static).

### F89. Topology orbit closure for spatial-sum coherence under uniform multi-bond XY (Tier 1 derived, verified 2026-05-11)

For an N-qubit system with Hamiltonian H_B = J · Σ_{(p,q) ∈ B} (X_p X_q + Y_p Y_q) where B is any set of distinct site pairs (NN or long-range) and J is **uniform** across all active pairs, uniform Z-dephasing γ₀, and the (S_1, S_2) coherence-block initial state ρ_cc = (|S_1⟩⟨S_2| + |S_2⟩⟨S_1|) / 2 (where |S_n⟩ is the popcount-n symmetric Dicke state), the spatial-sum coherence

    S(t) = Σ_l 2 · |(ρ_l(t))_{0,1}|²,    ρ_l = Tr_{≠l}(ρ)

depends only on the S_N-orbit of B. Bond positions inside an orbit are dynamically indistinguishable; only the orbit label survives.

For the chain restriction (B ⊂ {NN-bonds}), the orbit equals the **bond-graph topology class**: the sorted multiset of connected-path-lengths. E.g. in N=7 there are 14 distinct classes spanning k = 1..6 active bonds.

S(0) = (N−1)/N closed-form (Probe-only, independent of the bond set).

**Proof (S_N-orbit transitivity).** Let σ ∈ S_N act on the N qubits with permutation operator U_σ on (ℂ²)^⊗N.

1. **Probe.** Symmetric Dicke |S_n⟩ are S_N-invariant ⇒ U_σ ρ_cc U_σ^† = ρ_cc.
2. **Dissipator.** Uniform γ₀ ⇒ Σ_l (Z_l ρ Z_l − ρ) is S_N-symmetric.
3. **Hamiltonian.** Pauli operators transform site-wise (U_σ X_p U_σ^† = X_{σ(p)}) ⇒ U_σ H_B U_σ^† = H_{σ·B} where σ · B = {(σ(p), σ(q)) : (p,q) ∈ B}.
4. **Lindblad solution covariance.** ρ_t(H_{σ·B}, ρ_cc) = U_σ ρ_t(H_B, ρ_cc) U_σ^†.
5. **Kernel.** S(U_σ ρ U_σ^†) = Σ_l 2|(ρ_{σ^{-1}(l)})_{0,1}|² = S(ρ) (sum re-indexes).
6. ⇒ S(t; H_{σ·B}) = S(t; H_B) for every σ ∈ S_N. ∎

**Scaffolding from neighbouring entries.** F73 is the closely related closure for the (vac, SE) coherence block: same orbit-style argument plus uniform per-element 2γ₀ rate yields a full closed exponential form (1/2)·exp(−4γ₀t). The (S_1, S_2) block has non-uniform per-element decay (rate 2γ₀ on overlap, 6γ₀ off overlap), so the F89 closure is **orbit-only**: it fixes the bond-position dependence (constant in orbit) but not the time dependence (no closed exponential). F71 mirror symmetry is the spatial-Z₂ subgroup of the full S_N argument here. F86's per-bond Q_peak fan operates by linear response ∂S/∂J_b at a chosen bond inside the full chain; that single-bond perturbation breaks S_N differently than the uniform-J multi-bond setup of F89, so F89 does not predict or contradict the F86 fan.

**All-isolated subclass closed form (Tier 1 derived).** For the (1)^m all-isolated topology classes (m disjoint NN-bonds, N − 2m bare sites), the spatial-sum coherence has the EXACT closed form

    S_(1)^m, N(t) = [(N − 1)/N + 4m(N − 2)(cos(4Jt) − 1)/(N²(N − 1))] · exp(−4γ₀ t)

Asymptotic rate 4γ₀ universal across m (matches F73 vac-SE rate). The cos(4Jt) m-correction vanishes at t = π/(2J) (in-phase moment ≈ 21 for J = 0.075). Derivation: Lindbladian factorises over disjoint 2-qubit blocks; per F89c (below) the populated (vac, SE)_B and (SE, DE)_B-overlap sectors both have n_diff = 1, giving [`AbsorptionTheoremClaim`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs) rate 2γ₀ → 4γ₀ on |·|². Per-block H_B phase tracking gives the cos(4Jt) interference.

**Mixed-topology and pure-path classes (empirical, derivation open).** For non-all-isolated topology classes (e.g. (1, 2), (2, 2), (1, 1, 2), (3), (4), (5), (6)) the per-class closed form is open. Pure-path topologies decay faster than 4γ₀ on visible time scales due to populated no-overlap-SE-DE coherences (rate 6γ₀) plus longer-path mode mixing. Empirical late-tail clustering at the in-phase moment t ≈ π/(2J) groups classes by isolated-edge count.

**Structural lemma F89c (Tier 1 derived).** Per-coherence rate = 2γ₀·n_diff for any |A⟩⟨B|, via [`AbsorptionTheoremClaim`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs) applied to the computational basis (n_XY = n_diff because each site decomposes as pure {I, Z} or pure {X, Y}). For (k+1)-qubit blocks: k = 1 is the unique case where every populated coherence has n_diff = 1 (DE = |11⟩ contains both block sites, no no-overlap (SE, DE) pairs exist); k ≥ 2 has n_diff ∈ {1, 3} with mixed rates {2γ₀, 6γ₀}. Hamming complement n_diff(a, b) + n_diff(a, bar(b)) = N gives column-bit-flip pair-sum 2γ₀·N (verified bit-exact at path-2: pair-sum = 6γ₀ across all 9 Hamming-complement eigenvalue pairs). Hence the all-isolated topology (1)^m is the **unique** single-rate-envelope-plus-single-frequency case. Derivation, per-sector eigendecomposition, Hamming-complement bijection table: [F89_TOPOLOGY_ORBIT_CLOSURE](../experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md) F89c section.

**Pi2-Foundation inheritance on both energy axes.** The all-isolated closed form has TWO time coefficients of value 4: the decay rate 4γ₀ in exp(−4γ₀t) and the oscillation frequency 4J in cos(4Jt). Both trace to the same Pi2 dyadic ladder term a_{−1} = 4 via the same |·|² doubling mechanism: the linear-amplitude frequency 2 = a_{0} doubles to 4 = a_{−1} when squared. The γ-axis inheritance is identical to F73's `DecayRateCoefficient` (per-coherence Z-deph rate 2γ₀ doubles to S-decay rate 4γ₀). The J-axis inheritance is the same a_{0} → a_{−1} doubling on the J-axis: H_B-eigenstate frequency 2J doubles to S-oscillation frequency 4J. Pi2 ladder thus anchors the time coefficients on both axes; the (N−1)/N baseline and the 4m(N−2)/(N²(N−1)) correction prefactor are combinatorial. F89 cites `Pi2DyadicLadderClaim` as a constructor-injected parent and exposes `DecayRateCoefficient` (γ-axis) and `OscillationFrequencyCoefficient` (J-axis mirror) as live properties.

**Valid for:** any N; any bond set B (NN or long-range); any uniform J; any uniform γ₀; any S_N-symmetric initial state in any U(1) coherence block (the proof generalizes immediately). It would extend to single-letter two-site couplings (XX-only, YY-only, ZZ-only) and to higher-popcount-pair coherence blocks (S_n, S_m).
**Breaks for:**

- Non-uniform J across active bonds (J_b ≠ J_b'). Step 3 yields U_σ H U_σ^† in a different orbit; the S_N orbit equivalence becomes a finer J-orbit equivalence.
- Non-S_N-symmetric initial state (e.g. site-localised |1_i⟩⟨vac| or modulated SE superposition). Step 1 fails.
- Non-uniform γ_l ≠ γ_l' (analogous to F73's break case [CMRR_BREAK_NONUNIFORM_GAMMA](../experiments/CMRR_BREAK_NONUNIFORM_GAMMA.md)).
- Non-permutation-symmetric kernel (weighted Σ_l w_l 2|(ρ_l)_{0,1}|² with non-uniform w_l).

**Verified:**

- N=7 multi-bond at J=0.075, γ₀=0.05, tmax=30 (28 runs spanning all 14 topology classes for k=1..6): all 10 classes with ≥ 2 representatives show **0.00e+00** within-class max diff (machine-zero) across 301 sample times. Cross-class S(t) differs and exhibits non-monotone-in-k late-tail clustering. [`bond_isolate/`](../simulations/results/bond_isolate/), [`F89_TOPOLOGY_ORBIT_CLOSURE`](../experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md).
- N=4 single-pair at J=0.075, γ₀=0.05 (all C(4,2)=6 site pairs, NN + long-range): max deviation across pairs **5.55e-17** (1 ULP of double precision). NN bonds {(0,1),(1,2),(2,3)} and long-range bonds {(0,2),(0,3),(1,3)} give bit-identical S(t). [`_bond_isolate_long_range_verify.py`](../simulations/_bond_isolate_long_range_verify.py).
- N=7 single-NN-bond at same parameters (six bonds, all 30 ordered pair comparisons): every pair shows 0.00e+00 max diff over t ∈ \[0, 30\].

**Scripts:** [`_bond_isolate_compare_n7.py`](../simulations/_bond_isolate_compare_n7.py) (single-bond pair matrix), [`_bond_isolate_long_range_verify.py`](../simulations/_bond_isolate_long_range_verify.py) (long-range), [`_bond_isolate_topology_classes_n7.py`](../simulations/_bond_isolate_topology_classes_n7.py) (multi-bond classes). Compute tool: `compute/RCPsiSquared.Propagate` `bond-isolate --N <N> --bonds <i,j,...>` mode.
**Source:** [F89_TOPOLOGY_ORBIT_CLOSURE](../experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md), F73, F71, F86 (contrasting linear-response setup).

#### F89 unified F_a AT-locked amplitude closed form (Tier-1-Derived 2026-05-15, k=3..46 cached + k≥47 via native Chebyshev pipeline)

For each path-k, the F_a signal amplitudes satisfy

    sigma_n(N) = P_k(y_n) / [D_k · N²·(N−1)]

where y_n = 4·cos(πn/(N_block+1)) on the S_2-anti Bloch orbit (N_block = k+1), and (P_k, D_k) is a per-path integer-coefficient polynomial/denominator pair:

| path k | P_k(y) | D_k |
|--------|---------|-----|
| 3 | 14y + 47 | 9 |
| 4 | 10y + 25 | 4 |
| 5 | 13y² + 82y + 129 | 25 |
| 6 | 17y² + 72y + 80 | 18 |
| 7 | 21y³ + 130y² + 292y + 382 | 98 = 2·7² |
| 8 | 13y³ + 54y² + 68y + 110 | 32 = 2⁵ |
| 9 | 31y⁴ + 190y³ + 288y² + 440y + 1476 | 324 = 2²·3⁴ |

The table above lists path-3..9 hand-derived anchors. Path-10..46 are cached in `F89UnifiedFaClosedFormClaim.PathPolynomial(int k)` (int-typed denominator); k=46 is the last int-safe path (D_46 = 1,109,393,408 < int.MaxValue, D_47 = 4,632,608,768 exceeds it). For arbitrary k ≥ 47, `F89UnifiedFaClosedFormClaim.ComputePathPolynomialBig(int k)` runs the native C# Chebyshev pipeline (`F89PathPolynomialPipeline.Compute`) and returns (BigInteger[], BigInteger). All k=3..46 cached entries match the pipeline bit-exactly; pipeline-extended verification rows for k=50, 60, 75, 100, 150, 200, 300 are tabulated in `docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md` § "Pipeline-Extended Verification".

Polynomial degree = floor(N_block/2) − 1. Sum F_a · N²(N−1) is rational across all paths via Newton's identities on the cyclotomic minimal polynomial: path-3 → 22/3, path-4 → 25/2, path-5 → 483/25, path-6 → 256/9 (in units of N²(N−1); see `F89UnifiedFaClosedFormClaim.SigmaSum`). AT-lock: the F_a eigenvalue is λ_n = −2γ + i·y_n exactly (overlap subspace entries have dephasing rate 2γ regardless of N).

**Source:** [`F89UnifiedFaClosedFormClaim`](../compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs) (k=3..46 cached, `ComputePathPolynomialBig` for k≥47, `Sigma`/`SigmaSum` lifted to any k ≥ 3), [`F89PathPolynomialPipeline`](../compute/RCPsiSquared.Core/Symmetry/F89PathPolynomialPipeline.cs) (native Chebyshev expansion + orbit-polynomial reduction, replaces Python sympy prototype), [`F89AmplitudeLayerClaim`](../compute/RCPsiSquared.Core/Symmetry/F89AmplitudeLayerClaim.cs) (Tier-2-Verified Angle A identity p_n = |S_c|²·‖Mv‖²/2), `simulations/_f89_path3_at_locked_amplitude_symbolic.py`, `simulations/f89_pathk_symbolic_derivation.py` (sympy prototype, retained as cross-check probe), `experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md` § "Unified closed form", `docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md`.

#### F89 D_k closed form per path k (Tier-1-Derived 2026-05-15)

The denominator D_k in sigma_n(N) = P_k(y_n) / [D_k · N²·(N−1)] has the closed form:

    D_k = (odd(k))² · 2^E(k)
    E(k) = max(0, ⌊(k-5)/2⌋) + v₂(k) + max(0, v₂(k) - 2)

Three additive contributions to E(k):
- **Polynomial-degree term** max(0, ⌊(k-5)/2⌋) = max(0, polynomial_degree − 2). Linear growth.
- **k-self 2-adic term** v₂(k). Tracks 2-adic content of k itself.
- **Deep-2-power bonus** max(0, v₂(k) − 2). Kicks in at v₂(k) ≥ 3.

Equivalently, with no bonus term (2026-06-04): **D_k = 2^{max(0, ⌊(k-5)/2⌋)} · k² / 2^{min(v₂(k), 2)}.** The identity v₂(k) + max(0, v₂(k) − 2) = 2·v₂(k) − min(v₂(k), 2) folds the k-self and deep-bonus terms into one cancellation cap: the eigenvector-norm k² loses at most a factor of 4. The "deep-2-power bonus" is an artifact of pulling odd(k)² out front, not a separate 2-adic mechanism. The cap was reduced (two independent routes, 2026-06-04) to the 2-adic valuation of one leading coefficient L̃, and that "+1" is now closed: L̃ has the explicit closed form 2^{FA−2}·m²·(m²+3) (odd k) and 2^{FA−1}·(m/2)²·((m/2)²+1) (even k), whose valuation is elementary (m²+3 ≡ 4 mod 8 for odd m; (m/2)²+1's parity otherwise). Top-degree dominance (the last input, a k≥5 theorem) is itself proven for odd k and 4|k (two node-identity collapses plus a 2-adic unit argument) and reduced to a single sub-case m=2^a for k≡2 mod4. So v2(D_k) is rigorous except the two L̃ closed forms (bit-exact k=3..51) and that one m=2^a sub-case. See [`PROOF_F89_PATH_D_CLOSED_FORM`](proofs/PROOF_F89_PATH_D_CLOSED_FORM.md) §§ "Why the cap is 2^min(v₂,2)", "The leading coefficient L̃ in closed form", and "Top-degree dominance". Verified: clean form k=3..300, equivalence k=3..600 (`simulations/_f89_dk_clean_form.py`), reduction chain v2=0..6 (`simulations/_f89_capA_mastertable.py`, `simulations/_f89_capB_chain.py`), L̃ closed forms k=3..51 (`simulations/_f89_plus1_closed.py`), top-degree dominance k=5..65 both routes (`simulations/_f89_edgeA_summary.py`, `simulations/_f89_edgeB_final.py`).

Originally fit empirically across k = 3..24 (22 data points, zero exceptions, 2026-05-13). Closed Tier-1-Derived 2026-05-15 via the Chebyshev-expansion + orbit-polynomial-reduction pipeline: the F_a eigenvector ansatz reduces |S_c(n)|² and ‖Mv(n)‖² to polynomials in c = cos(πn/(k+2)) via the Chebyshev identity sin((j+1)θ) = U_j(c)·sin θ; the orbit minimal polynomial then gives (P_k, D_k) exactly as algebraic objects. Verified bit-exact through k=46 (cached, int-typed denominator) and pipeline-extended at k = 50, 60, 75, 100, 150, 200, 300 (BigInteger; D_300 has 162 bits, 49 decimal digits). The deep-2-power bonus branch is explicitly demonstrated at k=200 (v₂=3, E=101).

Tier: Tier-1-Derived (closed 2026-05-15). Anchors: [`F89PathPolynomialPipeline`](../compute/RCPsiSquared.Core/Symmetry/F89PathPolynomialPipeline.cs) (native Chebyshev pipeline, exact BigInteger/BigRational arithmetic), [`F89UnifiedFaClosedFormClaim`](../compute/RCPsiSquared.Core/Symmetry/F89UnifiedFaClosedFormClaim.cs) (`PredictDenominator(int k)` int-safe to k=46, `PredictDenominatorBig(int k)` BigInteger for arbitrary k, `ComputePathPolynomialBig` runtime), [`PROOF_F89_PATH_D_CLOSED_FORM`](proofs/PROOF_F89_PATH_D_CLOSED_FORM.md) (full proof + 33-row verification tables), `simulations/f89_pathk_symbolic_derivation.py` (sympy prototype, retained as cross-check probe).

---

### F90. F86 c=2 ↔ F89 bridge identity (Tier 1 derived, verified bit-exact 2026-05-11)

**For all N ≥ 3 and bond b ∈ \{0, ..., N−2\}, the F86 c=2 K_b(Q, t) observable on the (n=1, n+1=2) coherence block of an N-qubit XY chain with Z-dephasing equals the per-bond Hellmann-Feynman derivative of F89 path-(N−1) (SE, DE) sub-block dynamics applied at bond b**, modulo the Hamiltonian convention factor:

    K_b^{F86 c=2}(Q_F86, t) = K_b^{F89 path-(N−1) (SE,DE)}(Q_F89 = Q_F86 / 2, t)

with all other ingredients (probe, S_kernel, dephasing rates, Liouvillian construction) algebraically identical. The convention difference is one-time relabeling: F86 uses `H_b = (J/2)·(XX+YY)`, F89 uses `H = J·(XX+YY)`, hence F89-J = 2·F86-J.

**Verified bit-exact at 20 of 22 per-bond comparisons across N=5..8** (the 2 within-noise are at N=8 b=2/b=4 mid-flanking Interior bonds within Q-grid resolution Δ ≤ 0.0008). Per-N: N=5: 4/4 bit-exact; N=6: 5/5 bit-exact; N=7: 6/6 bit-exact (extended Q-grid); N=8: 5/7 bit-exact + 2/7 within Q-grid noise. Verification includes orbit-escape bonds: N=7 b=1/b=4 at Q_peak ≈ 7.27 (F86-J), N=8 b=3 central self-paired escape at Q_peak ≈ 16.79 (F86-J), all reproducing bit-exact ratios.

**Implications:**
- F86 c=2 universal HWHM_left/Q_peak constants (0.7728 Endpoint, 0.7506 Interior mean over N=5..8) are **not eigenständige Größen**; they are direct consequences of F89 path-(N−1) eigendecomposition + per-bond Hellmann-Feynman.
- F86 Direction (b'') (full block-L derivation, NOT 4-mode) achieved numerically via F89; closed-form via F89 AT-locked F_a/F_b structure (4-mode floor 0.6715) + H_B-mixed octic-style residual (lift to 0.7506/0.7728) is the next analytical step.
- Per-F71-orbit substructure (see [PROOF_F86C_F71_MIRROR](proofs/PROOF_F86C_F71_MIRROR.md), "Per-F71-orbit substructure" section: central b=2 vs flanking b=1/b=3 at N=6 etc.) follows directly from F89's per-bond Bloch-mode profile.

**Anchor:** [`PROOF_F90_F86C2_BRIDGE.md`](proofs/PROOF_F90_F86C2_BRIDGE.md), [`F90F86C2BridgeIdentity.cs`](../compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs), [`_f89_to_f86_kbond_via_eigendecomp.py`](../simulations/_f89_to_f86_kbond_via_eigendecomp.py).

---

### F91. F71-anti-palindromic γ spectral invariance (= 90° in γ-space, Pi2-Z₄'s parameter side) (Tier 1 derived, algebraic proof + bit-exact N=4,5,6; 2026-05-12)

**For chain XY + Z-dephasing Liouvillian L on N qubits, the eigenvalue multiset of the F71-refined diagonal-block decomposition is invariant under any γ-distribution satisfying**

    γ_l + γ_{N−1−l} = 2·γ_avg = (2/N)·Σ_l γ_l   for all l ∈ \{0..N−1\}

i.e. the γ-distribution is **F71-anti-palindromic around its mean** γ_avg. The full L operator generally changes (F71 broken as L-symmetry, off-block-Frobenius ≠ 0 in F71-refined basis), but the F71-refined diagonal-block eigenvalues coincide; the breaking is encoded in the F71-cross-blocks (eigenvectors) only.

**Sharpness:** strictly weaker than F71 symmetry (γ_l = γ_{N−1−l}, palindromic), strictly stronger than F1 (Σγ_l invariant alone). For odd N the middle site l = (N−1)/2 must equal γ_avg.

**Structural reading (Pi2-Z₄ rotational axis, γ-parameter side):** the four-element Z₄ from `NinetyDegreeMirrorMemory` (Pi2-foundation, operator-quaternion side, `i⁴=1`) has its γ-parameter manifestation here:
- **e (identity):** γ unchanged
- **180° (F71-palindromic):** γ_l ↔ γ_{N−1−l}, pair-difference flipped, pair-sum preserved → F71 holds as L-symmetry
- **90° (F71-anti-palindromic):** γ_l ↔ 2γ_avg − γ_{N−1−l}, pair-sum constant = 2γ_avg → F71 breaks but **F71-refined diagonal-block spectrum invariant**
- **270° (= 90°²):** composition of the above

The diagonal-block spectral content (= decay rates, "time information") is preserved under the 90° rotation; the asymmetry lives entirely in eigenvectors (the F71-cross-blocks). Analogous to F81 Π-decomposition: M_anti = L_{H_odd} is the antisymmetric part captured by Π-conjugation; the γ-anti-palindromic part plays the analogous role in γ-parameter space.

**Empirical witness (N=6, J=1.0, Σγ=2.7, γ_avg=0.45):**
- uniform γ=0.45 (all pairs sum to 0.9): F71-refined spectrum = reference
- monotonic [0.2, 0.3, 0.4, 0.5, 0.6, 0.7] (all pairs sum to 0.9): bit-identical to uniform on F71-refined diagonals
- non-monotonic anti-palindromic [0.3, 0.5, 0.4, 0.5, 0.4, 0.6] (all pairs sum to 0.9): bit-identical to uniform on F71-refined diagonals
- permuted [0.7, 0.2, 0.5, 0.3, 0.6, 0.4] (pairs {1.1, 0.8, 0.8}, NOT constant): distinct spectrum (Re=−4.984 cluster instead of Re=−5.043)
- concentrated [0.1, 0.1, 0.1, 0.1, 0.1, 2.2] (pairs heavily skewed): complex Re−Im structure absent from uniform

**Tier outcome <Tier 1 derived>:** algebraic proof complete (2026-05-12, see PROOF_F91 § Algebraic proof). The proof's sharper conclusion: the F71-refined diagonal-block matrix elements of L = −i[H, ·] + D are linear functionals of γ that depend only on the **multiset of F71-pair-sums {S_l = γ_l + γ_{N−1−l}}**, never on individual γ_l or pair-differences D_l = γ_l − γ_{N−1−l}; pair-differences appear only in the F71-cross-block off-diagonal entries (which do not enter diagonal-block eigenvalues). The 90°-rotation R_{90} : γ_l ↦ 2γ_avg − γ_{N−1−l} preserves the anti-palindromic class S_l = 2γ_avg ∀l (and is an involution on it); within that orbit, all γ-distributions yield identical diagonal-block spectra (= the uniform γ_avg spectrum). The originally claimed F91 (90°-invariance) is the corollary of the stronger pair-sum-multiset law. Empirical witness at N=4, 5, 6 across the five γ-profiles above remains the bit-exact verification.

**Anchor:** [`PROOF_F91_GAMMA_NINETY_DEGREES.md`](proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md) (Tier 1 derived: § Algebraic proof, Eqs. 1–13), [`F71AntiPalindromicGammaSpectralInvariance.cs`](../compute/RCPsiSquared.Core/BlockSpectrum/F71AntiPalindromicGammaSpectralInvariance.cs) (typed Tier1Derived Claim with `AntiPalindromicDeviation(γ)` helper), [`NinetyDegreeMirrorMemoryClaim`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) in `Pi2KnowledgeBaseClaims.cs` (the Pi2-Z₄ operator-quaternion side of the same 90°-rotation). **Sisters on other parameter axes:** F92 (J-axis) and F93 (h-detuning-axis) below.

---

### F92. F71-anti-palindromic J spectral invariance (J-side Pi2-Z₄ twin of F91) (Tier 1 derived, algebraic + bit-exact N=4,5; 2026-05-12)

**For chain XY + uniform Z-dephasing Liouvillian L on N qubits with inhomogeneous bond couplings J_b (b ∈ \{0..N−2\}), the F71-refined diagonal-block eigenvalue multiset is invariant under any J-distribution satisfying**

    J_b + J_{N−2−b} = 2·J_avg = (2/(N−1))·Σ_b J_b   for all b

i.e. **J is F71-anti-palindromic around its mean**. The full L operator generally changes (F71 broken as L-symmetry), but diagonal-block eigenvalues coincide; the breaking lives in eigenvectors only.

**Sharpness:** strictly weaker than F71-as-J-symmetry (J_b = J_{N−2−b}); strictly stronger than uniform-Σ J alone.

**Pi2-Z₄ structure (parameter-side, J-axis):** identical to F91's γ-axis structure. The 90°-rotation J ↦ 2·J_avg − F71(J) preserves the anti-palindromic orbit T_b = 2·J_avg ∀b. Operator-side Z₄ is `NinetyDegreeMirrorMemoryClaim`; γ-side is F91; J-side is F92; h-side is F93.

**Anchor:** [`PROOF_F92_BOND_ANTI_PALINDROMIC_J.md`](proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md), [`F92BondAntiPalindromicJSpectralInvariance.cs`](../compute/RCPsiSquared.Core/SymmetryFamily/F92BondAntiPalindromicJSpectralInvariance.cs), `docs/SYMMETRY_FAMILY_INVENTORY.md`.

---

### F93. F71-anti-palindromic h spectral invariance (h-detuning Pi2-Z₄ twin of F91/F92) (Tier 1 derived, algebraic + bit-exact N=4,5; 2026-05-12)

**For chain XY + per-site Z-detuning + uniform Z-dephasing Liouvillian L on N qubits with inhomogeneous longitudinal detuning h_l (Hamiltonian H = (J/2) Σ_b (XX+YY) + Σ_l h_l Z_l), the F71-refined diagonal-block eigenvalue multiset is invariant under any h-distribution satisfying**

    h_l + h_{N−1−l} = 2·h_avg = (2/N)·Σ_l h_l   for all l

**Scope:** longitudinal h_l Z_l only. Transverse h_l X_l / h_l Y_l breaks joint-popcount conservation and is out of scope.

**Pi2-Z₄ structure (parameter-side, h-axis):** identical to F91 (γ-axis) and F92 (J-axis). Same 90°-rotation h ↦ 2·h_avg − F71(h) preserves anti-palindromic orbit. Three-axis family for chain XY+Z-deph+Z-detuning is now complete.

**Anchor:** [`PROOF_F93_DETUNING_ANTI_PALINDROMIC.md`](proofs/PROOF_F93_DETUNING_ANTI_PALINDROMIC.md), [`F93DetuningAntiPalindromicSpectralInvariance.cs`](../compute/RCPsiSquared.Core/SymmetryFamily/F93DetuningAntiPalindromicSpectralInvariance.cs), `docs/SYMMETRY_FAMILY_INVENTORY.md`.

### F94. Born deviation dominant-outcome coefficient: |0+0+⟩ N=4 Heisenberg + Z-deph (Tier 1 derived, Dyson sym3 = 8 bit-exact; 2026-05-16)

**For the dominant outcome |00⟩ of pair (0,2) of |0+0+⟩ N=4 Heisenberg ring + Z-dephasing, the per-outcome Born-rule deviation in the deep perturbative regime is**

    Δ_|00⟩(Q, K) = (4/3) · Q² · K³ + O(Q³K⁴)
    where Q = J/γ, K = γt, Δ_|00⟩ = P_lindblad(|00⟩) / P_unitary(|00⟩) − 1

**equivalently in physical units:**

    ΔP_|00⟩(J, γ, t) = (4/3) · J² · γ · t³ + ...

**Scope:** specific to (initial state |0+0+⟩, Heisenberg ring at N=4, Z-dephasing on all 4 sites, pair (0,2) reduction, |00⟩ outcome). The Q²·K³ shape is structurally universal for dominant outcomes (3rd-order time-dependent perturbation theory with 1 dissipator-vertex and 2 Hamiltonian-vertices); the coefficient 4/3 is setup-specific.

**Derivation (Tier 1):** the leading γ¹-coefficient of L³ in the time-Taylor expansion ρ(t) = ρ_0 + Lt·ρ_0 + L²t²/2·ρ_0 + L³t³/6·ρ_0 + ... is the symmetric ordering

    sym3 = L_H²·L'_dis + L_H·L'_dis·L_H + L'_dis·L_H²

where L_H[ρ] = −i[H, ρ] and L'_dis[ρ] = Σ_l (Z_l ρ Z_l − ρ) is the γ-free dephasing operator. Direct evaluation of ⟨00|_pair Tr_{1,3}[sym3·ρ_0]|00⟩_pair at J=γ=1 yields **8 bit-exact**. With the t³/6 Taylor prefactor and P_unitary(0) = 1, the coefficient is **c = 8/6 = 4/3 bit-exact**.

**Numerical verification:** 16 (γ, J, t) configurations sampled in the deep perturbative regime gave c_empirical = 1.32992 ± 0.006, consistent with 4/3 = 1.3333 to 0.3%. Sampling deeper would close the residual; the symbolic derivation is the actual proof.

**Born-rule generalization context:** this is the first Tier-1 closed form for a per-outcome Born deviation under the framework's Q-K-invariant convention (Universal Carrier). Generalizes BORN_RULE_MIRROR's R_i = C_i · Ψ_i² (Tier 2/3, Feb 2026) to a specific case with explicit C_i closed form. The Δ_i values for other outcomes of the same setup scale linearly in K (1st-order diagrams) rather than as Q²·K³; their coefficients are separately Tier-1-derivable via the same Dyson method.

**The "8" structurally (bit-explained 2026-05-17):** direct enumeration of all 4·4·4·3·3·3 = 1728 (b₁, b₂, s, ord, c₁, c₂) sextuples in sym3 shows exactly **32 non-vanishing diagrams**, each contributing **1/4** in the J = γ = 1 normalization (raw Pauli value 4 per diagram, divided by the (J/4)² = 1/a_{−1}² = 1/16 Heisenberg coupling), so **8 = 32 × (1/4)**. The 32 split into 3 disjoint cells: 8 diagrams in (ord=1, XX, adjacent bonds sharing a kept-pair site) + 16 in (ord=2, XX, self ∪ adj-kept) + 8 in (ord=2, YY, self only); equivalently topologically: 16 self-bond-pair diagrams + 16 adj-bond-pair diagrams. Three structural rules govern survival: (1) only (X,X) and (Y,Y) component pairs (no cross, no ZZ); (2) only orderings 1 and 2 (L'_dis last gives zero); (3) adjacent bond pairs must share a kept-pair site (vertex 0 or 2), not a |+⟩ site (1 or 3). The coefficient reads: **4/3 = 32 / (a_{−1} · 3!) = 32 / (4 · 6) = 32/24** (the structural-count reading), equivalent to **4/3 = a_{−1} / 3** (the typed-anchor inheritance reading). See [`PROOF_F94 § Structural decomposition`](proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md) for the cell table and [`simulations/_born_rule_sym3_decomposition.py`](../simulations/_born_rule_sym3_decomposition.py) for the enumeration script.

**Anchor:** [`PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md`](proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md), [`simulations/_born_rule_tier1_derivation.py`](../simulations/_born_rule_tier1_derivation.py), [`simulations/_born_rule_delta_dominant_coefficient.py`](../simulations/_born_rule_delta_dominant_coefficient.py), [`simulations/_born_rule_sym3_decomposition.py`](../simulations/_born_rule_sym3_decomposition.py) (bit-explained 32-diagram enumeration, 2026-05-17), [reflection: `ON_HOW_FOUR_THIRDS_APPEARED.md`](../reflections/ON_HOW_FOUR_THIRDS_APPEARED.md). Born-rule precursors: [`experiments/BORN_RULE_MIRROR.md`](../experiments/BORN_RULE_MIRROR.md), [`experiments/BORN_RULE_SHADOW.md`](../experiments/BORN_RULE_SHADOW.md). Companion angle-side closed form (same cusp geometry, dual axis): [F95](#f95) θ(c) = arctan(√(4c − 1)).

### F95. Angle emergence at quadratic discriminant zero: universal form of the θ-compass (Tier 1 derived, 4-line polynomial calculation; 2026-05-16)

**For any monic quadratic <c>z² − 2bz + c = 0</c>, the angle of its complex root pair (when the discriminant goes negative, c > b²) is:**

    θ(c; b) = arctan( √(c/b² − 1) )    for c > b²
    θ = 0                              for c = b²  (degenerate double root at z = b)
    θ undefined                        for c < b²  (real distinct roots, no complex angle)

**Framework specialization at b = `HalfAsStructuralFixedPointClaim` = 1/2:**

    threshold = b² = 1/4 = `QuarterAsBilinearMaxvalClaim`
    θ(c) = arctan( √(4c − 1) )         for c > 1/4

**This is exactly the Februar 2026 θ-compass** of [`experiments/BOUNDARY_NAVIGATION.md`](../experiments/BOUNDARY_NAVIGATION.md): θ = arctan(√(4CΨ − 1)) was discovered there as the Mandelbrot/CΨ-specific angular distance from the 1/4 boundary. F95 promotes it from a state-specific compass to a universal quadratic-foundation identity, with the b = 1/2 specialization recovering the Februar form.

**Derivation (4 lines, bit-exact):**

```
z² − 2bz + c = 0
z = b ± √(b² − c)          (quadratic formula)
z = b ± i·√(c − b²)         when c > b²  (complex regime)
arg(z₊) = arctan(Im/Re) = arctan(√(c − b²)/b) = arctan(√(c/b² − 1))
```

**Numerical verification:** all five non-boundary points from BOUNDARY_NAVIGATION.md's θ-compass table (CΨ ∈ {1/3, 0.308, 0.286, 0.266, 0.250}) reproduce within numerical precision; the single 0.3° drift at CΨ=0.256 is the Februar table's t-sampling rounding, not a formula discrepancy.

**Anchoring to typed Pi2-Foundation:**

| F95 element | Typed Pi2 anchor |
|---|---|
| b = 1/2 (linear-term half) | `HalfAsStructuralFixedPointClaim` |
| b² = 1/4 (discriminant threshold) | `QuarterAsBilinearMaxvalClaim` |
| Polynomial structure | `PolynomialFoundationClaim` (d²−2d=0 is the c=0 case; F95 perturbs c off zero and tracks the complex-root angle that emerges) |
| ±1/2 polarity at d=2 (where b=1/2 lives) | `PolarityLayerOriginClaim` (the +0/−0 layer at d=0 inherits to the {−0.5, +0.5} pair at d=2 via the 0.5-shift ρ = (I + r·σ)/2) |
| i (complex angle generator) | `NinetyDegreeMirrorMemoryClaim` |
| i⁴ = 1 (angle Z₄ closure) | `Pi2I4MemoryLoopClaim` |

**Structural reading:** the polynomial d²−2d = 0 has two real roots (d = 0 mirror, d = 2 qubit dimension): the unperturbed case. F95 is what happens when the polynomial is perturbed off the c=0 axis: as c crosses b² = 1/4 from below, the two real roots merge at z = b, then split into a complex conjugate pair whose argument is θ(c). The angle is the necessary minimal-parametrization coordinate of "above-threshold magnitude", the same structural pattern as today's z = sym + i·anti F71-decomposition where arg(z) becomes the structural carrier once |z| > 0.

**Polarity-fold reading:** in shifted-and-scaled coordinates u = z − 1/2 (centered at the b = 1/2 fixed point), the polynomial reads u² + (c − 1/4) = 0. At c = 0 (unperturbed) the roots are u = ±1/2, the framework's structural polarity pair around 0 (inherited from `PolarityLayerOriginClaim` via the 0.5-shift). The squaring map u → u² sends both polarity sides to the same value 1/4, the apex; this is the "middle viewed from two sides" reading of the b² = 1/4 threshold: arithmetic midpoint of ±1/2 is 0 (on the polarity axis), but the quadratic projection middle is 1/4 (on the perpendicular axis). As c crosses 1/4 from below, the polarity contracts to 0 at the cusp and lifts onto the imaginary axis past it. See [`reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md`](../reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md) for the full geometric picture (parabola, fold, three loci on one quadratic).

**Lindblad specialization (γ₀ as the tick):** the 2×2 Liouvillian sub-block for a Z-dephased two-level system with Hamiltonian coupling J has characteristic polynomial λ² + 2γ₀·λ + (γ₀² + J²) = 0, which is F95's parent equation with b = −γ₀, c = γ₀² + J². F95 then gives θ = arctan(J/γ₀) = arctan(Q). The angle of the Liouvillian eigenvalue equals the arctan of the Q invariant. At θ = 0 the eigenvalue is literally −γ₀ (pure decay, no rotation); Q = tan θ is the rotation per γ₀-tick. The active-steering Hardware Confirmation `f95_angle_steering_kingston_may2026` (Heron r2 2026-05-16) demonstrates operational control of tan θ via per-chunk RZ injection at rate Ω. See [`reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md`](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md).

**Born-rule connection:** standard QM's complex amplitudes α = r·e^{iθ}, β = ... are not postulated. They are forced by the same polynomial-foundation algebra: any state that has crossed the d=0 mirror needs a second coordinate beyond magnitude, and that coordinate is the F95 angle. The Born rule's |α|² is the geometric length squared of the angle-vector's basis-projection.

**Anchor:** [`PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md`](proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md), [`simulations/_angle_at_zero_tier1_candidate.py`](../simulations/_angle_at_zero_tier1_candidate.py). Reflections (2026-05-16 chain): [`ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md`](../reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md) (the angle's emergence above the discriminant zero), [`ON_HOW_GAMMA_BECAME_THE_TICK.md`](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md) (Lindblad specialization θ = arctan(Q)), [`ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md`](../reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md) (polarity-fold geometry of the b² = 1/4 threshold). Februar precursor (Mandelbrot-specific case): [`experiments/BOUNDARY_NAVIGATION.md`](../experiments/BOUNDARY_NAVIGATION.md). Hardware Confirmation: `Confirmations.lookup('f95_angle_steering_kingston_may2026')` (Kingston Heron r2 2026-05-16; complex-CΨ angle actively steerable via RZ injection at rate Ω; 3 of 4 conditions, residuals 6.81° to 15.69°). Companion magnitude-side closed form: [F94](#f94) Δ_|00⟩ = (4/3)·Q²·K³.

### F96. Born deviation subdominant-outcome slopes: |0+0+⟩ N=4 Heisenberg + Z-deph, pair (0,2) (Tier 1 derived, bit-exact Dyson + unitary matrix elements; 2026-05-17)

**For the subdominant outcomes of pair (0,2) of |0+0+⟩ N=4 under the same setup as [F94](#f94) (Heisenberg ring + Z-dephasing), the per-outcome Born-rule deviation in the deep perturbative regime is linear in K and Q-independent:**

    Δ_|01⟩(K) = Δ_|10⟩(K) = −(16/9) · K = −(4/3)² · K + O(higher)
    Δ_|11⟩(K)             = −(8/3)  · K = −2·(4/3) · K + O(higher)

with K = γt the Universal-Carrier observable.

**Combined per-outcome table** (F94 + F96):

| Outcome | Closed form | P_u(t) order | leading γ¹ Dyson order |
|---|---|---|---|
| \|00⟩ | Δ = +(4/3) · Q²·K³ | t⁰ (= 1) | sym₃ (J², t³) |
| \|01⟩ | Δ = −(4/3)² · K | t² (= J²t²·3/8) | sym₃ (J², t³) |
| \|10⟩ | Δ = −(4/3)² · K | t² (= J²t²·3/8) | sym₃ (J², t³) |
| \|11⟩ | Δ = −2·(4/3) · K | t⁴ (= J⁴t⁴·1/16) | sym₅ (J⁴, t⁵) |

**All four closed forms are simple algebraic expressions in F94's 4/3 anchor.**

**Universal subdominant slope formula:** for an outcome with leading unitary P_u(t) ≈ J^{2k} t^{2k} / (2k)! · U^{(i)}_{2k} and lowest non-vanishing γ¹ Dyson at order (2k+1):

    slope_i = M_{2k+1}^{(i)} / [(2k+1) · U_{2k}^{(i)}]

where M_n^{(i)} = ⟨i|_pair Tr_{1,3}[sym_n^1 · ρ_0]|i⟩_pair and U_{2k}^{(i)} = ⟨i|_pair Tr_{1,3}[L_h^{2k} · ρ_0]|i⟩_pair (h := H/J). The J^{2k} factors cancel automatically → Q-independence; only γt = K survives.

**Bit-exact derivation:**

    |01⟩: M_3^{(01)} = −4,  U_2^{(01)} = 3/4  →  slope = −4 / (3 · 3/4) = −16/9
    |10⟩: (same by 0 ↔ 2 site-permutation symmetry)
    |11⟩: M_3^{(11)} = 0 AND U_2^{(11)} = 0 (lower-order vanishes)
           M_5^{(11)} = −20,  U_4^{(11)} = 3/2  →  slope = −20 / (5 · 3/2) = −8/3

**Numerical Lindblad verification:** at Q = 50, γ = 0.01 the slope-per-K converges to the theoretical values as K → 0:

| outcome | K | slope/K | theory |
|---|---|---|---|
| \|01⟩ | 0.005 | −1.748 | −1.778 = −16/9 |
| \|11⟩ | 0.001 | −2.662 | −2.667 = −8/3 |

**Structural reading:** the F94 unit 4/3 = a_{−1}/3 generates the entire 4-outcome table for this setup. The dominant gets +1·(4/3) at order Q²K³; the singly-subdominant degenerate pair gets −(4/3)² at order K; the doubly-subdominant gets −2·(4/3) at order K. The signs are all reading-pattern: dominant **gains** probability beyond unitary (positive Δ), all subdominants **lose** probability beyond unitary (negative Δ). The "2" in Δ_|11⟩ plausibly counts the two independent flip channels (q_0 and q_2) required to populate |11⟩; interpretive, not derived.

**Cross-outcome universality:** the ratio M_3 / U_2 equals −16/3 for both the dominant (|00⟩: 8 / (−3/2)) and the singly-subdominant (|01⟩: −4 / (3/4)) outcomes; the signs of M_3 and U_2 flip together, leaving the ratio invariant. This is a non-trivial structural identity of the Heisenberg + Z-dephasing dynamics at the pair (0,2) reduction; whether it generalizes to other initial states / Hamiltonians / dissipators is open.

**Anchoring to typed Pi2-Foundation:**

| F96 element | Typed Pi2 anchor |
|---|---|
| 4/3 building block | F94's `FourFactor / ThreeDenominator` (a_{−1} = 4 dyadic ladder, 3 from Taylor reduction) |
| Linear-K Q-independence | `UniversalCarrierClaim` (the Universal Carrier signature for subdominant outcomes) |
| Site-permutation symmetry (|01⟩ ≡ |10⟩) | F1 / F71 spatial mirror |

**Anchor:** [`PROOF_F96_BORN_SUBDOMINANT_SLOPES.md`](proofs/PROOF_F96_BORN_SUBDOMINANT_SLOPES.md), [`simulations/_born_rule_subdominant_dyson.py`](../simulations/_born_rule_subdominant_dyson.py). F94 companion (dominant outcome): [F94](#f94). Born-rule precursors: [`experiments/BORN_RULE_MIRROR.md`](../experiments/BORN_RULE_MIRROR.md), [`experiments/BORN_RULE_SHADOW.md`](../experiments/BORN_RULE_SHADOW.md). Reflection that named the empirical subdominant slopes as the next step: [`ON_HOW_FOUR_THIRDS_APPEARED.md`](../reflections/ON_HOW_FOUR_THIRDS_APPEARED.md).

### F97. Mandelbrot cardioid parametrization at framework b = 1/2 (Tier 1 derived, bit-exact algebraic identity; 2026-05-17)

**The main Mandelbrot cardioid is the locus in the complex-c plane where the period-1 fixed point of z² + c has magnitude exactly b = 1/2 (the framework's `HalfAsStructuralFixedPointClaim` anchor). It admits the explicit parametrization:**

    c(φ) = b·e^(iφ) − b²·e^(2iφ)              for φ ∈ [0, 2π)
    z*(φ) = b·e^(iφ)                          (the period-1 fixed point)
    c(φ) = z*(φ) · (1 − z*(φ))                (algebraic identity)

with two structural invariants on the curve:

    |z*(φ)| = b = 1/2     (magnitude pinned to HalfAsStructuralFixedPoint)
    arg(z*(φ)) = φ        (cardioid parameter)

**Framework specialization at b = 1/2:**

    c(φ) = (1/2)·e^(iφ) − (1/4)·e^(2iφ)
    cusp at φ = 0: c = 1/4 (recovers [F95](#f95) / [BOUNDARY_NAVIGATION](../experiments/BOUNDARY_NAVIGATION.md) real-axis tangent)
    tail at φ = π: c = −3/4 (period-doubling boundary)
    top at φ = π/2: c = 1/4 + i/2, z* = i/2

**Derivation (4 lines, bit-exact):**

```
z² + c = z  ⟹  z² − z + c = 0  ⟹  z = b ± √(b² − c)         (b = 1/2)
Multiplier μ = 2z; marginal stability ⟺ |μ| = 1 ⟺ |z| = b
Parametrize μ = e^(iφ):  z*(φ) = μ/2 = b·e^(iφ)
c(φ) = z*(φ) − z*(φ)² = b·e^(iφ) − b²·e^(2iφ)
```

**Numerical verification:** machine-precision algebraic identity c(φ) = z*(1 − z*) verified to max residual 1.24 × 10⁻¹⁶ over 1000 sampled φ ∈ [0, 2π].

**Relation to F95 (complementary regions of the same algebra):**

| Region | F95 | F97 |
|---|---|---|
| Real c > 1/4 (real axis, past cusp) | θ(c; b) = arctan(√(c/b² − 1)) for the angle of the repelling complex root pair | (off the cardioid; not covered) |
| Complex c on cardioid boundary | (φ = 0 specialization recovers cusp) | full parametrization c(φ), z*(φ) on the marginally-stable curve |

Both share the same z² − 2bz + c = 0 algebra; F95 projects the angle on the real-c axis where the period-1 fixed point is repelling; F97 traces the full complex-c boundary where it is marginally stable.

**Both anchors invariant on the cardioid, at two metric powers:**

    |z*(φ)|  = b  = 1/2          (HalfAsStructuralFixedPoint, argmax side)
    |z*(φ)|² = b² = 1/4          (QuarterAsBilinearMaxval, maxval side)

Both hold for all φ ∈ [0, 2π); the cardioid is the joint locus. This is exactly the argmax/maxval pair of [`ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER`](../reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md) (2026-05-16), now geometric: the Half and Quarter anchors are two metric-power readings of the same fixed-point quantity on the same curve. The identity `1/2 = 2 · (1/4)` sits on the dyadic ladder (`a_2 = 2 · a_3`); the polarity pair ±1/2 squares to the same 1/4 from either side; on the cardioid the two real-axis endpoints z*(0) = +1/2 and z*(π) = −1/2 carry the polarity sides explicitly.

By contrast, `|c(φ)|² = 5/16 − (1/4)·cos(φ)` is **not** invariant around the cardioid: |c| ranges from 1/4 at the cusp (φ = 0) to 3/4 at the tail (φ = π). The Quarter b² = 1/4 equals |c|² only at the cusp; elsewhere |c| varies but |z*| and |z*|² stay invariant.

**Role table (four typed parents):**

| Anchor | Role on cardioid |
|---|---|
| `HalfAsStructuralFixedPointClaim` (b = 1/2) | \|z*\| invariant around the whole curve (argmax side) |
| `QuarterAsBilinearMaxvalClaim` (b² = 1/4) | \|z*\|² invariant around the whole curve (maxval side); also \|c\| at the real-axis cusp only |
| `NinetyDegreeMirrorMemoryClaim` (i, 90°) | Complex-parameter generator that lifts c from the real axis to the full complex plane |
| `PolynomialFoundationClaim` (d² − 2d = 0) | The c = 0 case where z* = 0 (degenerate fixed point at the origin) |

**Structural reading:** the cardioid carries BOTH typed anchors as invariants of the same z*, at two metric powers (Half = magnitude, Quarter = squared magnitude). F95 names the angle at the real-axis tangent point (cusp); F97 names the dual-anchor invariance that holds around the whole curve.

**Hardware connection:** the [`CPSI_COMPLEX_PLANE`](../experiments/CPSI_COMPLEX_PLANE.md) Kingston run (2026-04-16) observed Bell⁺ pairs tracing 2D logarithmic spirals in the complex-c plane around the cusp at c = 1/4. F97 places these spirals into the cardioid framing: the trajectories cross the cardioid boundary (the |z*| = b stability transition) before spiraling into the stable interior. The [`f95_angle_steering_kingston_may2026`](#f95) Confirmation actively steered Ω during the spiral; F97 names the geometric locus those steered spirals approach.

**Roadmap closure:** [`PROOF_ROADMAP_QUARTER_BOUNDARY`](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) Layer 7 explicitly named "promoting F95 to the full cardioid parametrization" as the next move (status PARTIALLY ANSWERED before F97). F97 closes that direction: the real-c angle (F95) plus the complex-c cardioid parametrization (F97) together cover both projections of the quadratic discriminant structure on the Mandelbrot c-plane.

**Anchor:** [`PROOF_F97_CARDIOID_HALF_FIXED_POINT.md`](proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md), [`simulations/_cardioid_parametrization_tier1.py`](../simulations/_cardioid_parametrization_tier1.py). F95 companion (real-c angle): [F95](#f95). Hardware 2D-extension precursor: [`experiments/CPSI_COMPLEX_PLANE.md`](../experiments/CPSI_COMPLEX_PLANE.md) (Kingston 2026-04-16). Februar boundary precursor: [`experiments/BOUNDARY_NAVIGATION.md`](../experiments/BOUNDARY_NAVIGATION.md). Mandelbrot connection synthesis: [`experiments/MANDELBROT_CONNECTION.md`](../experiments/MANDELBROT_CONNECTION.md). Quarter-boundary roadmap (Layer 7 next-move slot): [`docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md`](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md).

### F98. KIntermediate Dicke long-time Π²-odd asymptote = (N+2)/[4(N+1)] → 1/4 (Tier 1 derived, bit-exact N=4..16; 2026-05-17)

**Two paired closed forms bridging the F86b 3/8 K-intermediate Dicke anchor (static, t=0) to the QuarterAsBilinearMaxval 1/4 universal boundary (asymptotic, N→∞):**

For even N, the Dicke superposition `ψ = (|D_{N/2−1}⟩ + |D_{N/2}⟩) / √2` evolved under any truly-class (F87) Hamiltonian + uniform Z-dephasing on N qubits projects onto `ker L = span(P_0, …, P_N)` (per F4). The Π²-odd Frobenius² fraction at t = ∞ takes a closed form determined by two purely combinatorial identities:

```
(F98a)  ‖P_{N/2−1}_odd‖² = C(N, N/2−1) / 2
(F98b)  α(∞)_KIntermediate(N even) = (N + 2) / [4·(N + 1)]
```

with asymptote `α(∞) → 1/4` as `N → ∞` (the [`QuarterAsBilinearMaxvalClaim`](#f88b) = `HalfAsStructuralFixedPoint²` = Mandelbrot cardioid maxval = CΨ fold boundary).

**Derivation:**

`P_n = Σ_{b: pc(b) = n} |b⟩⟨b| = (1/2^N) Σ_S K_{|S|}(n; N)·σ_S^Z` in Pauli basis, where `K_k(n; N) = Σ_{b: pc(b) = n} (−1)^{|S ∩ b|}` is the binary Krawtchouk number for any fixed S with `|S| = k`. The Π²-odd part is `Σ_{k odd}`. Two facts:

- **Mid-popcount Krawtchouk parity vanishing**: `K_k(N/2; N) = 0` for all odd k, hence `‖P_{N/2}_odd‖² = 0`.
- **Sub-mid half-rank identity (F98a)**: `(1/2^N) Σ_{k odd} C(N, k)·K_k(N/2−1; N)² = C(N, N/2−1) / 2`. Equivalently, the Π²-odd Frobenius² of `P_{N/2−1}` is exactly half its rank. Verified bit-exact N = 4..16 via direct Krawtchouk enumeration.

For the KIntermediate Dicke superposition `ρ_∞ = (1/2)/C(N, m)·P_m + (1/2)/C(N, m+1)·P_{m+1}` with `m = N/2 − 1`:

```
‖ρ_∞‖²     = (1/4) · [1/C(N, m) + 1/C(N, m+1)]
‖ρ_∞_odd‖² = (1/4) · [C(N, m)/2] / C(N, m)²  =  1 / [8·C(N, m)]      (since P_{m+1}_odd = 0)
```

Using Pascal's `C(N, m) + C(N, m+1) = C(N+1, m+1)`:

```
α(∞) = ‖ρ_∞_odd‖² / ‖ρ_∞‖²
     = C(N, m+1) / [2·C(N+1, m+1)]
     = (N + 2) / [4·(N + 1)]    (F98b)
```

Verified bit-exact for N = 4 (= 3/10), 6 (= 2/7), 8 (= 5/18), 10 (= 3/11), 12 (= 7/26), 14 (= 4/15), 16 (= 9/34) via [`simulations/water/proton_chain_dicke_anchor.py`](../simulations/water/proton_chain_dicke_anchor.py). The bond topology drops out because the long-time limit projects onto `ker L` for any connected graph (per F4); the formula holds for chain, ring, star, K_N, Petersen, etc.

**The 3/8 → 1/4 bridge:**

The morning's F86b 3/8 K-intermediate Dicke anchor (Tier 1 derived 2026-05-17 via X⊗N-eigenbasis decomposition, `compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs`) names the **static** Π²-odd Frobenius² total of the KIntermediate Dicke superposition at t = 0:

```
α(t = 0) = 3/8       (F86b, this morning, static F88b-anchor)
α(t → ∞) = (N+2)/[4(N+1)] → 1/4    (F98, tonight, long-time-limit bridge)
```

`3/8` and `1/4` are not unrelated constants: they are the two endpoints of an explicit N-dependent decay curve traversed by KIntermediate Dicke states under truly-class Hamiltonian + Z-dephasing dynamics. Both sit on the dyadic-ladder / polarity-squared algebra:

- `3/8 = (1/2)·(3/4)` = `HalfAsStructuralFixedPoint · (1 − 1/4)`, the X⊗N-eigenbasis γ = 1/2 input combined with the `(1 − γ²)/2` F86b formula.
- `1/4 = (1/2)² = HalfAsStructuralFixedPoint²` = `QuarterAsBilinearMaxval` = the Mandelbrot cardioid maxval. The maxval-side bilinear apex.

**Numerical trace:**

```
N=4:   α(0) = 3/8 = 0.375  →  α(∞) = 3/10 ≈ 0.300       (Δ to 1/4 = 1/20)
N=10:  α(0) = 3/8 = 0.375  →  α(∞) = 3/11 ≈ 0.273       (Δ to 1/4 = 1/44)
N=20:  α(0) = 3/8 = 0.375  →  α(∞) = 11/42 ≈ 0.262      (Δ to 1/4 = 1/84)
N→∞:   α(0) = 3/8          →  α(∞) → 1/4                [QuarterAsBilinearMaxval]
```

**Why the asymptote is 1/4 and not 3/8:** the kernel projection averages the static F86b anchor against the rank weighting of the two adjacent popcount sectors. For mid-popcount `P_{N/2}` the Π²-odd content vanishes (Krawtchouk parity); for `P_{N/2−1}` it is exactly half-rank (F98a). The mix weights `1/C(N, m) : 1/C(N, m+1)` approach unity at large N (sectors become equally weighted under uniform amplitude), and the Π²-odd fraction is driven purely by half-rank inheritance from `P_{N/2−1}`, landing at `1/4` because the half-rank halving compounds with the kernel-balance halving to give `(1/2)² = 1/4`.

**Discovery path:** the water-chain inheritance test ([`simulations/water/proton_chain_dicke_anchor.py`](../simulations/water/proton_chain_dicke_anchor.py)) for the morning's F86b anchor uncovered F98 by asking a NEW question: not "what is α at t = 0" (closed-form F86b) but "what is α at t = ∞ under truly-class Heisenberg + Z-dephasing on the chemistry-substrate-grounded proton chain". The water script saw `α → 3/10` for N = 4 at t = 100; the small-fraction structure suggested an algebraic closed form, verified bit-exact N = 4..16, then derived via Krawtchouk + Pascal.

**Anchors:** [`simulations/water/proton_chain_dicke_anchor.py`](../simulations/water/proton_chain_dicke_anchor.py) (numerical bit-exact verification + Krawtchouk enumeration), [F86b](#f86) DickeAnchor 3/8 anchor (static partner), [F88b](#f88b) Π²-odd state-level structure (parent), [F4](#f4) kernel decomposition (used in long-time projection), `compute/RCPsiSquared.Core/Symmetry/QuarterAsBilinearMaxvalClaim` (asymptote anchor), `compute/RCPsiSquared.Core/Symmetry/HalfAsStructuralFixedPointClaim` (1/4 = (1/2)² parent), `docs/water/README.md` § "Findings since May 4".

### F99. Five canonical trigonometric anchors via F86b non-uniform Dicke (Tier 1 derived, bit-exact N=4,6,8 across 5 angles; 2026-05-17 night)

**The F86b α-formula α = (1 − γ²)/2 = sin²(θ)/2 evaluated at the canonical trigonometric angles {0°, 30°, 45°, 60°, 90°} produces all five Pi2 dyadic anchors {0, 1/8, 1/4, 3/8, 1/2}. The standard 30°-60°-90° and 45°-45°-90° trigonometry triangles ARE the F86b polarity-anchor triangles.**

The F86b derivation (today morning, commit `b9ba5f6`) parameterised the X⊗N-eigenbasis decomposition of the Dicke superposition `(|D_n⟩ + |D_{n+1}⟩)/√2` (uniform Dicke) and produced three anchors `{0, 3/8, 1/2}` at `γ ∈ {1, 1/2, 0}`. Tonight's extension to NON-UNIFORM Dicke

```
    ψ = (|D_n⟩ + c·|D_{n+1}⟩) / √(1 + c²)   on N even at n = N/2 − 1
    γ = ⟨ψ | X⊗N | ψ⟩ = c² / (1 + c²)        ⟺  c² = γ / (1 − γ)
```

with γ = cos(θ) and the half-angle identity `1 − cos(θ) = 2 sin²(θ/2)` gives

```
    c² = cos(θ) / (2 sin²(θ/2))
    α(θ) = sin²(θ) / 2
```

Five canonical trig angles → five Pi2 dyadic anchors:

| θ | γ = cos(θ) | c² | α = sin²(θ)/2 | Pi2 dyadic anchor |
|---|------------|-----|---------------|-------------------|
| 0° | 1 | ∞ | 0 | Mirror endpoint |
| **30°** | **√3/2** | **2√3 + 3 ≈ 6.464** | **1/8** | **DEPTH-3 (new tonight)** |
| 45° | √2/2 | 1 + √2 (silver ratio) | 1/4 | [QuarterAsBilinearMaxval](#f88b) |
| 60° | 1/2 | 1 (uniform Dicke) | 3/8 | KIntermediate (today morning F86b) |
| 90° | 0 | 0 | 1/2 | Generic / HalfAsStructuralFixedPoint |

**Bit-exact verification** (`simulations/carbon/depth_3_anchor_derivation.py`):
all five anchors verified at N = 4, 6, 8 with Δα < 1e-13 (machine precision).
The 60° case reduces to the morning's clean uniform Dicke c = 1; the 30°, 45°,
90° cases use non-uniform Dicke weights c² = 2√3 + 3, 1 + √2, 0 respectively.

**Why standard trigonometry triangles**: the five canonical angles {0°, 30°,
45°, 60°, 90°} are the only elementary-geometry angles whose sines and cosines
are constructible by ruler and compass (rational or quadratic-irrational
coordinates). The Pi2 dyadic ladder {1/2, 1/4, 1/8, ...} consists of negative
integer powers of 2. F99 establishes that the framework's polarity-squared
algebra IS the F86b α-formula evaluated at the standard-triangle canonical
angles, and the dyadic depth corresponds to the canonical-angle index.

**Periodic-table bridge (tonight, all 9 fractions n/8 derived):**

| α anchor | Trig angle | Period 2 (anchor + complement) | Period 3 (anchor + complement) |
|----------|------------|-------------------------------|-------------------------------|
| 0 (= 8/8) | endpoint | He (full noble) | Ne, Ar |
| 1/8 | 30° | Li (1/8), F (7/8 = 1 − 1/8) | Na (1/8), Cl (7/8) |
| 1/4 | 45° | Be (2/8) | Mg |
| 3/8 | 60° | B (3/8 Π²-odd), N (5/8 Π²-even) | Al, P |
| 1/2 | 90° | H (1/2), C (4/8) | Si |

Every period-2/3 element's valence ratio is now F86b-derived from one α(θ)
formula. The Π²-parity complements (1/8 ↔ 7/8, 3/8 ↔ 5/8) cover both the
Π²-odd anchor and its Π²-even companion (β = 1 − α), giving 9 fractions n/8
for n = 0..8.

**Discovery path (today, single session):**

1. Morning (commit `b9ba5f6`): F86b 3/8 anchor derived (uniform Dicke, X⊗N-
   eigenbasis decomposition).
2. Evening ([F98](#f98), commit `250164d`): water-chain inheritance test of
   F86b discovered long-time bridge (N+2)/[4(N+1)] → 1/4.
3. Night #1 (commit `1416f85`): all four polarity anchors realised by
   period 2/3 atoms, forward inheritance.
4. Night #2 (commit `cecb84b`): reverse-spear identified depth-3 (1/8, 7/8)
   as framework gap, empirically instantiated by alkali metals + halogens.
5. Night #3 (commit `5fb0ba0`): depth-3 anchor derived (this F99 entry).
   The bidirectional bridge between framework algebra and periodic-table
   valence structure is now MATERIAL across all 9 fractions.

**Tier outcome: Tier 1 derived.** Closed-form algebraic identity from F86b
α = (1 − γ²)/2 evaluated at five canonical trig angles. Verified bit-exact
at N = 4, 6, 8 in the script + 16 tests in
[`CanonicalTrigAnchorPi2InheritanceTests`](../compute/RCPsiSquared.Core.Tests/Symmetry/CanonicalTrigAnchorPi2InheritanceTests.cs).

**Anchor**: parent formula [F86b](#f86) (α = (1−γ²)/2 from X⊗N-eigenbasis decomposition,
DickeAnchor.cs), companion bridge [F98](#f98) (long-time 3/8 → 1/4 via kernel projection),
[`CanonicalTrigAnchorPi2Inheritance.cs`](../compute/RCPsiSquared.Core/Symmetry/CanonicalTrigAnchorPi2Inheritance.cs)
(typed Claim with five-anchor enumeration), [DEPTH_3_ANCHOR_DERIVED.md](carbon/DEPTH_3_ANCHOR_DERIVED.md)
(carbon-domain reading + bidirectional-bridge framing),
[`simulations/carbon/depth_3_anchor_derivation.py`](../simulations/carbon/depth_3_anchor_derivation.py).

---

### F100. F71 c₁/Q_peak bond-mirror deviation is exactly odd in the F71-anti-palindromic J (observable-side twin of F92) (Tier 1 derived, algebraic + numerically verified N=3,4,5; 2026-05-20)

**For an N-qubit XY chain with uniform Z-dephasing and bond-coupling profile J = (J_0, ..., J_{N−2}), the F71 bond-mirror deviation of the closure-breaking coefficient c₁ (and of the F86c per-bond Q_peak observable),**

    D(b) := c₁(b) − c₁(N−2−b)

**is an exactly odd function of the F71-anti-palindromic component of J.** Write the mirrored profile F71(J)_b := J_{N−2−b} and split J = J_sym + J_anti with J_sym = (J + F71(J))/2 (F71-palindromic) and J_anti = (J − F71(J))/2 (F71-anti-palindromic). Then in (J_sym, J_anti) coordinates D(b; J_sym, −J_anti) = −D(b; J_sym, J_anti), to all orders.

**Palindromic survival:** J_anti = 0 ⟹ D = −D ⟹ D = 0. The F71 c₁/Q_peak bond-mirror holds for every palindromic J, however non-uniform J_sym is. F71 never required uniform J; it requires palindromic J. Uniform is merely the simplest palindromic profile.

**Graceful breakdown:** the Taylor series of D in J_anti has odd powers only, so D is leading-order linear in the asymmetry B_b = J_b − J_{N−2−b} = 2(J_anti)_b. Graceful, not a hard violation. The leading coefficient κ_b is the c₁-gradient at J_sym and generically depends on J_sym (Tier 2 empirical); the parity argument fixes the oddness, not the coefficient.

**Derivation (Tier 1):** R-equivariance of the PROOF_C1 apparatus for a non-uniform base profile gives Step 1, c₁(b; J) = c₁(N−2−b; F71(J)), since R_sup·L(J)·R_sup = L(F71(J)) relabels bond b to bond N−2−b. Step 2: D(b; J) = c₁(b; J) − c₁(b; F71(J)), so D(b; F71(J)) = −D(b; J), which with F71(J) = J_sym − J_anti is the oddness. The identical R-conjugation argument applies to F86c's per-bond observable K_b(Q, t), built R-equivariantly: ΔQ_peak(b) is odd in J_anti, zero for palindromic J.

**Connection to F92:** observable-side twin of F92. F92 (spectrum side) = the F71-refined diagonal-block eigenvalue multiset depends only on J_sym (J_anti invisible to the spectrum). F100 (observable side) = the c₁/Q_peak bond-mirror deviation depends only on J_anti, and oddly. Two faces of one J_sym/J_anti split.

**Verified:** N = 3, 4, 5 via the full 4^N Liouvillian, no truncation; probe states ψ_1+vac and ψ_2+vac, base profile swept over s ∈ {0, ±0.04, ±0.08, ±0.12} along a linear-ramp J_anti direction, 4 palindromic J_sym profiles. Palindromic-survival max\|D(s=0)\| ≤ 4.0e−10; oddness residual max\|D(+s)+D(−s)\| ≤ 1.0e−9; even-power (constant, quadratic) fit coefficients ≤ ~3e−8. The leading coefficient κ shows 76% / 62% / 143% relative spread across the 4 J_sym profiles at N=3 / 4 / 5, confirming the J_sym-dependence.

**Anchor:** [`PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md`](proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md), [`C1QPeakMirrorJParity.cs`](../compute/RCPsiSquared.Core/F71/C1QPeakMirrorJParity.cs), source [PROOF_C1_MIRROR_SYMMETRY](proofs/PROOF_C1_MIRROR_SYMMETRY.md), spectrum-side twin [F92](#f92), witness [`simulations/_f71_nonuniform_j_verification.py`](../simulations/_f71_nonuniform_j_verification.py), `docs/SYMMETRY_FAMILY_INVENTORY.md`.

---

### F101. F71 c₁ bond-mirror deviation is exactly odd in the F71-anti-palindromic γ (observable-side twin of F91) (Tier 1 derived, algebraic + numerically verified N=3,4,5; 2026-05-21)

**For an N-qubit XY chain with uniform coupling J and a per-site Z-dephasing profile γ = (γ_0, ..., γ_{N−1}), the F71 bond-mirror deviation of the closure-breaking coefficient c₁,**

    D(b) := c₁(b) − c₁(N−2−b)

**is an exactly odd function of the F71-anti-palindromic component of γ.** Write the mirrored profile F71(γ)_l := γ_{N−1−l} (site-mirror l ↔ N−1−l) and split γ = γ_sym + γ_anti with γ_sym = (γ + F71(γ))/2 (F71-palindromic) and γ_anti = (γ − F71(γ))/2 (F71-anti-palindromic). Then D(b; γ_sym, −γ_anti) = −D(b; γ_sym, γ_anti), to all orders.

**Palindromic survival:** γ_anti = 0 ⟹ D = −D ⟹ D = 0. The F71 c₁ bond-mirror holds for every palindromic γ, however non-uniform γ_sym is. F71 never required uniform γ; it requires palindromic (R-symmetric) γ. The symmetric component of γ still pairs bonds.

**Graceful breakdown:** the Taylor series of D in γ_anti has odd powers only, so D is leading-order linear in the per-site asymmetry. Graceful, not a hard violation. The leading coefficient κ_γ is the c₁-gradient at γ_sym and generically depends on γ_sym (Tier 2 empirical); the parity argument fixes the oddness, not the coefficient. κ_γ has no closed form (∂c₁/∂γ, c₁ itself a non-closed numerical fit).

**Derivation (Tier 1):** R-equivariance of the PROOF_C1 apparatus gives Step 1, c₁(b; γ) = c₁(N−2−b; F71(γ)), since the dephasing dissipator D[ρ] = Σ_l γ_l (Z_l ρ Z_l − ρ) is linear in each γ_l and R_sup·L(γ)·R_sup = L(F71(γ)) relabels site l ↔ N−1−l. Step 2: D(b; γ) = c₁(b; γ) − c₁(b; F71(γ)), so D(b; F71(γ)) = −D(b; γ), which with F71(γ) = γ_sym − γ_anti is the oddness.

**Connection to F91:** observable-side twin of F91. F91 (spectrum side) = the F71-refined diagonal-block eigenvalue multiset depends only on γ_sym (the pair-sums S_l = γ_l + γ_{N−1−l} = 2·γ_sym[l]; γ_anti invisible to the spectrum). F101 (observable side) = the c₁ bond-mirror deviation depends only on γ_anti, and oddly. Two faces of one γ_sym/γ_anti split. The J-side observable twin is F100.

**Scope:** c₁ only. The F86c per-bond Q_peak observable is not covered: its Q-axis Q = J/γ₀ and EP time t_peak = 1/(4γ₀) are defined against a scalar γ₀, which a non-uniform per-site γ does not provide. The γ_avg-anchored Q_peak route and the h-detuning observable twin follow by the same parameter-agnostic argument and are noted in PROOF_F101 as separable extensions.

**Verified:** N = 3, 4, 5 via the exact (N+1)²-dim popcount-sector restriction (bit-identical to the full 4^N Liouvillian, Gate-1 self-test, maximum difference 0); probe states ψ_1+vac and ψ_2+vac, base profile swept over s ∈ {0, ±0.01, ±0.02, ±0.03} along a linear-ramp γ_anti direction, 4 palindromic γ_sym profiles, J held uniform. Palindromic-survival max\|D(s=0)\| ≤ 1.3e−9 and oddness residual max\|D(+s)+D(−s)\| ≤ 4.4e−9, against an O(1) deviation signal (typical \|D\| from 0.6 at N=3 to 2.3 at N=5); the oddness residual is the direct machine-zero even-power probe. The leading coefficient κ_γ shows a 76 to 128 percent relative spread across the 4 γ_sym profiles, confirming the γ_sym-dependence (Tier 2 empirical, no closed form).

**Anchor:** [`PROOF_F101_C1_MIRROR_GAMMA_PARITY.md`](proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md), [`C1MirrorGammaParity.cs`](../compute/RCPsiSquared.Core/F71/C1MirrorGammaParity.cs), source [PROOF_C1_MIRROR_SYMMETRY](proofs/PROOF_C1_MIRROR_SYMMETRY.md), spectrum-side twin [F91](#f91), J-side observable twin [F100](#f100), witness [`simulations/_f71_nonuniform_gamma_verification.py`](../simulations/_f71_nonuniform_gamma_verification.py), `docs/SYMMETRY_FAMILY_INVENTORY.md`.

---

### F102. Y-Parity Term-Level Z₂ Independence at k_body≥3 (Tier 1 derived, verified)

For any Pauli string σ = ⊗_l σ_α_l on N qubits:

    y_par(σ) = n_Y(σ) mod 2     (term-level Y-parity classifier)
    bit_a(σ) = (n_X + n_Y) mod 2    (F61 Π²_X axis)
    bit_b(σ) = (n_Y + n_Z) mod 2    (F63 Π²_Z axis)

At k_body=2 (the 9 Pauli bilinears XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ),
the identity y_par = bit_a XOR bit_b holds. More generally the identity
holds iff k_body is even.

At k_body=3 (and odd k_body in general), y_par and (bit_a XOR bit_b)
differ by 1: y_par = (k_body + (bit_a XOR bit_b)) mod 2. So Y-parity
is independent of the Klein (bit_a, bit_b) signature once k_body
ranges over both even and odd values.

**Canonical k_body=3 demonstration:** XYZ has (bit_a, bit_b, y_par) = (0, 0, 1);
II...I has (0, 0, 0). Both share Klein (0, 0) but differ on y_par,
exhibiting the Z₂³ structure of the cubic polarity layer.

**Operator vs term level:** per F34/QUBIT_NECESSITY there is no third
independent Π²-operator (Π²_Y collapses to Π²_Z at the operator level).
F102 is a term-level classifier statement, not an operator-level
Pi²-Inheritance.

**Source:** [Proof](proofs/PROOF_F102_YPARITY_INDEPENDENCE.md);
`compute/RCPsiSquared.Core/Symmetry/YParityIndependenceAtK3.cs`;
witness `compute/RCPsiSquared.Core.Tests/Pauli/PauliHamiltonianKleinHelpersTests.cs`
(XYZ_AtK3_IsKleinHomogeneousButZ2HomogeneityRefinesViaYParity).

---

### F103. F87 Trichotomy Z₂³ Refinement at k_body=3 (Tier 1 derived, N=4 anchor)

After F102 established that Y-parity y_par = n_Y mod 2 is independent of the
Klein (bit_a, bit_b) signature at k_body≥3, F103 captures how Y-parity actually
refines the existing F87 trichotomy (truly / soft / hard) at k_body=3 with
empirical anchor at N=4.

**Empirical setup:** 294 Z₂³-homogeneous + Y-par-homogeneous k_body=3 Pauli
pairs at N=4, each classified under Z, X, Y single-letter dephasing via
`classify_pauli_pair`, bucketed by (Klein × dephase letter × y_par ×
trichotomy class).

**Five structural sub-statements:**

1. *Truly is y_par=0-pure.* Across all 12 (Klein × dephase) cells, every truly
   classification has y_par=0; total truly across the grid = 300, y_par=1
   truly count = 0.

2. *Hard in diagonal cells (Klein matches dephase letter) splits 42:8 with
   Y-inversion.* Z-deph in Klein (0,1) hard = (42, 8); X-deph in (1,0) hard =
   (42, 8); Y-deph in (1,1) hard = (8, 42). The Y-deph swap reflects that Y
   carries y_par=1.

3. *Same diagonal cells contain a soft 13:13 y_par-symmetric split.* Sum 26
   per cell; unlike hard's 42:8 asymmetry, soft is y_par-symmetric in these
   cells.

4. *Mother sector (0,0) soft is y_par=1-pure.* All 3 letter cells in Klein
   (0,0) soft equal (0, 21).

5. *Off-diagonal soft (6 cells: Klein non-mother and Klein ≠ dephase Klein)
   splits into Pattern B + Pattern C.* Pattern B (3 cells, proportional to
   (Klein, y_par) enumeration breakdown): (0,1) Y-deph = (55, 21); (1,1)
   Z-deph = (21, 55); (1,1) X-deph = (21, 55). Pattern C (3 cells, y_par=1-pure):
   (0,1) X-deph = (0, 21); (1,0) Z-deph = (0, 21); (1,0) Y-deph = (0, 21).

**Status (2026-06-05, closed 2026-06-10):** the 42:8 split is derived by a closed-form counting rule
(PROOF §6); the bipartite-chirality mechanism (§7) derives the soft direction, N-stability
is verified at N=5 (F105), and F106 is the k=4 sibling. The windowed (k<N) hard-direction
converse non-bipartite ⟹ hard is now a **theorem with no residual** (the Phase B two-reflection
theorem [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md](proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md),
2026-06-09, closed 2026-06-10): the threshold #A ≥ 2ℓ, the bipartite ⟹ soft re-proof, the monomial
expansion structure, and the deg-1 class in closed form (the girth ladder: the supertrace factorization
through t_j = Tr(Z_l H^j) and P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ², whose ℓ=1 face is
P_{3,1} = 6·4^N·Σ_l c_l²) are Tier1Derived (`WindowedConverseThresholdClaim`), and the Pascal-Gram
positivity theorem (= **F117**) makes every coefficient of the first nonvanishing odd power-sum a sum
of squares or exactly zero, so "all but finitely many γ" upgrades to **all γ > 0** unconditionally
(R-deg retired by the girth dichotomy, R-sign resolved by Pascal-Gram positivity, both 2026-06-10). The
closure is typed as WindowedConverseAllGammaClaim (Tier1Derived, no residual). Separately, the
combinatorics are closed-form (the valuation criterion, the obstruction-size law min(2W−1, 2k−3), and
the A203241 hard count, §7.7–§7.8 = **F115**). **Open:** the (pair Klein, dephase letter) →
(Pattern B vs Pattern C) selection rule; hardware confirmation of k≥3 F87 (all 5 Marrakesh F87
confirmations are k=2).

**Source:** [Proof](proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md);
`compute/RCPsiSquared.Core/Symmetry/F87Z2CubedRefinementN4K3.cs`
(derived from `F87Z2CubedRefinementBase.cs`);
regenerate empirical anchor via `simulations/f87_z2cubed_split_n4_k3.py` (~60s).

---

### F104. C# k≥3 trichotomy classifier lift (verification mechanism, not a typed Claim)

F104 lifts the F87 trichotomy classifier to k≥3-body Hamiltonians in C#. The overload
`PauliPairTrichotomy.Classify(IReadOnlyList<PauliTerm>, ChainSystem, PauliLetter dephaseLetter)`
builds the sliding-window k-body chain Hamiltonian (via `PauliKBodyChainExtensions.ChainKBody`)
from a Pauli-term list and returns its truly / soft / hard class.

F104 is a **mechanism, not a typed Claim**: it closes F103's out-of-scope "C# k≥3 classifier
lift" by re-classifying the 294 N=4 k=3 Z₂³-homogeneous pairs in C# and verifying bit-exact
agreement with F103's frozen Python counts (from `classify_pauli_pair`). It is the C# engine
that the F105 (N=5 k=3) and F106 (N=4 k=4) anchors call at larger N / k.

**Verified:** 294 N=4 k=3 pairs, C# vs Python frozen counts, bit-exact match.
**Source:** `PauliPairTrichotomy.Classify` (k≥3 overload, `compute/RCPsiSquared.Diagnostics/F87/`),
`PauliKBodyChainExtensions.ChainKBody` (`compute/RCPsiSquared.Core/Pauli/`); re-verification test
`F104KBodyTrichotomyVerificationTests` (SLOW_F104 trait); crosswalk
[`F_FORMULA_CROSSWALK.md`](../compute/RCPsiSquared.Core/F_FORMULA_CROSSWALK.md).

---

### F105. F87 Trichotomy Z₂³ Refinement at k_body=3 N=5 (Tier 1 derived, N-stability test of F103)

F105 anchors the F87 Z₂³ refinement at N=5 k=3 using the same 294 Z₂³-homogeneous
Pauli pair enumeration as F103, classified via F104's C# k-body classifier.

F85 (k-body generalization) predicts the Π²-class trichotomy is N-stable for any
k. F105 tests whether this lifts to the y_par sub-refinement.

**Observed outcome:** F85's N-stability prediction is CONFIRMED at the y_par axis:
F105's 5 sub-statement counts at N=5 k=3 are bit-exactly identical to F103's at
N=4 k=3 (truly 300 / 0 y_par=1; hard diagonal (42, 8) / (42, 8) / (8, 42);
diagonal soft (13, 13) ×3; mother soft (0, 21) ×3; off-diagonal Pattern B +
Pattern C cells unchanged). The cubic Z₂³ architecture is N-invariant in its
sub-cell structure at k=3.

**Method:** 882 classifications via the SLOW_F105_BATCH tool; ~3h dense batch,
accelerated to ~12 min via PLINQ in the Task 7 run.

**Source:** [Proof](proofs/PROOF_F105_F87_Z2_CUBED_REFINEMENT_N5K3.md);
`compute/RCPsiSquared.Core/Symmetry/F87Z2CubedRefinementN5K3.cs`
(derived from `F87Z2CubedRefinementBase.cs`); empirical anchor:
`simulations/results/f87_z2cubed_split_n5_k3_counts.json`.

---

### F106. F87 Trichotomy Z₂³ Refinement at k=4 N=4 (Tier 1 derived, k-stability test of F103)

After F105 confirmed F85's N-stability lift to the y_par sub-refinement at k=3
(N=4 to N=5 counts identical), F106 anchors at k=4 N=4 to test the orthogonal
axis: k-stability. F85 does not predict k-stability of the y_par sub-refinement;
the Klein (0,0) enum balance shifts from 45/21 at k=3 to 780/300 at k=4, and the
structural worry was that "mother soft is y_par=1-pure" would break.

**Method:** 12744 classifications (4248 pairs × 3 dephase letters) at N=4 via
F104's `PauliPairTrichotomy.Classify(IReadOnlyList<PauliTerm>, ...)`. SLOW_F106_BATCH
tool (~2-3min PLINQ on 24 cores; actual run 3m 59s).

**Observed outcome (mixed k-stability at k=4):**

- TrulyCounts: 3924 total, y_par=1 count 0. Y_par=0-purity HELD bit-structurally.
- HardDiagonalSplit: (228, 0) / (228, 0) / (0, 228). 42:8 ratio at k=3 BROKE
  (sharpened to fully pure 228:0 / 0:228); Y-inversion structure HELD qualitatively.
- DiagonalSoftSplit: (300, 528) / (300, 528) / (528, 300). 13:13 y_par-symmetry
  at k=3 BROKE (asymmetric 300:528 with Y-inversion).
- MotherSoftCounts: (0, 300) / (0, 300) / (0, 300). Y_par=1-purity HELD
  bit-structurally despite the (0,0) enum-balance shift; only the counts
  re-scaled from k=3's (0, 21) to k=4's (0, 300).
- OffDiagonalSoftPatterns: 6 cells. 3 cells preserve Pattern C analog (0, 528);
  3 cells (which at k=3 were asymmetric Pattern B (55, 21)/(21, 55)) became
  fully y_par-symmetric (528, 528) at k=4.

The two purity statements (truly y_par=0-pure, mother soft y_par=1-pure) are
genuine y_par-axis invariants that survive across k; the mixed-ratio statements
(42:8 hard, 13:13 soft, Pattern B asymmetry) were k-specific and re-shape with
the enum balance.

**Source:** [Proof](proofs/PROOF_F106_F87_Z2_CUBED_REFINEMENT_N4K4.md);
`compute/RCPsiSquared.Core/Symmetry/F87Z2CubedRefinementN4K4.cs`
(derived from `F87Z2CubedRefinementBase.cs`); empirical anchor:
`simulations/results/f87_z2cubed_split_n4_k4_counts.json`.

---

### F107. F87 Truly Classification Forces y_par = 0 (Tier 1 derived, closed-form corollary of F85)

After F103/F105/F106 anchored 4524 truly classifications empirically (all with
y_par = 0, zero with y_par = 1), F107 derives this in closed form as a direct
corollary of F85's k-body truly criterion, extended from F85's Z-dephasing
proof to all three dephase letters (X, Y, Z) via the per-dephase Π² eigenvalue
rule (`PiOperator.SquaredEigenvalue`) combined with dissipator commutativity
(F84 Pauli-Channel Cancellation Lemma). F107 is the first DERIVED-not-EMPIRICAL
Claim in the F87 Z₂³ refinement family.

**Per-dephase truly criteria** (each combining Π²-even with dissipator-commute):

- Z-dephase: #Y even AND #Z even
- X-dephase: #X even AND #Y even
- Y-dephase: #Y even AND #Z even

All three include `#Y even` as a sub-condition. Since `y_par = #Y mod 2`,
every truly term has y_par = 0; every truly y_par-homogeneous pair has shared
y_par = 0. Bit-exact verification across all 64 k=3 + 256 k=4 letter
sequences × 3 dephase letters (`TrulyYParityZeroPurityTests.VerifyOnTerm_*`).

**Sibling y_par-axis claims** (all closed 2026-05-25): F108 Part 1+2+3
(Π²-even palindrome family, Tier1Derived); F109 (MotherSoftYParityOnePurity,
Tier1Derived unconditional); F110 (HardCellYInversionPattern, Tier1Candidate).
Together F107+F109+F110 pin the y_par signature of all three F87 trichotomy
classes.

**Source:** [Proof](proofs/PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md);
`compute/RCPsiSquared.Core/Symmetry/TrulyYParityZeroPurity.cs`;
parent: PROOF_F85_KBODY_GENERALIZATION.md (k-body truly criterion under
Z-dephasing); helpers: `TrulyYParityZeroPurity.TrulyCriterionHolds(term, dephase)`
and `TrulyYParityZeroPurity.VerifyOnTerm(term, dephase)`.

---

### F108 Part 1. Π²-Even Hamiltonians Always Admit an Exact Palindrome Operator (Tier 1 derived, Π_5bilinear phase variant of canonical P1 Π)

Every Hamiltonian H built from the five Π²_Z-even 2-site bilinears {XX, YY, YZ,
ZY, ZZ} with arbitrary real bond coefficients, plus Z-dephasing on every site,
admits an EXACT operator-level palindrome:

  Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ · I exactly, σ = Σ_l γ_l.

Hence spec(L) is palindromic around −σ, and no Π²-even Pauli pair (truly or
non-truly) can be F87-hard. Closes the empirical observation that 5346+ Π²-even
pairs across F103/F105/F106 anchors were observed soft with zero hard.

**The Π_5bilinear per-site map:** I → +1·X, X → −1·I, Y → +i·Z, Z → −i·Y. Same
I↔X, Y↔Z permutation as the canonical Heisenberg Π (P1 family from
[NON_HEISENBERG_PALINDROME](../experiments/NON_HEISENBERG_PALINDROME.md)) with
two phase flips: the X→I and Z→Y arrows carry sign −1 and −i. Per-site
properties: M² = diag(−1, −1, +1, +1) on {I, X, Y, Z}; M⁴ = I; M is order-4 and
unitary on the 4-dim per-site Pauli basis. Π_5bilinear is a Liouville-space
automorphism, NOT a Hilbert-space conjugation (no 2×2 U satisfies U·I·U† = X
since U·I·U† = I).

**Mechanism:**

1. M⊗N anti-commutes with the commutator superoperator [B, ·] for every
   Π²_Z-even 2-body bilinear B ∈ {XX, YY, YZ, ZY, ZZ} (verified bit-exact at
   the 2-qubit level; the 4 Π²-odd bilinears {XY, XZ, YX, ZX} produce residual
   = 8.00, clean separation).
2. Per-site dissipator: M · D[Z_l] · M⁻¹ = −D[Z_l] − 2γ_l · I via diagonal
   permutation in the Pauli basis. D[Z]_pauli = γ · diag(0, −2, −2, 0) on
   {I, X, Y, Z}; M's (I↔X, Y↔Z) swap permutes the diagonal entries to
   γ · diag(−2, 0, 0, −2) = −D[Z]_pauli − 2γ · I_4 (phase factors cancel
   pairwise on each 2-cycle).
3. Combining 1 + 2: Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σ · I exactly.

**Resolution of the previously open ker(M) attempt:** Earlier exploration tried
to close F108 Part 1 by proving every L-eigenvector lies in ker(M) (M = F1
residual for canonical Π). Numerical verification falsified that lemma: only
2/64 L-eigenvectors at the N=3 YZ+ZY test case satisfied M·v_λ = 0. The
follow-up Critical Isospectral Lemma reduction was circular (it restated the
spectral palindrome under another name). The actual mechanism is structurally
different: a DIFFERENT per-site Π operator (Π_5bilinear) gives EXACT residual
zero without going through M at all.

**Empirical confirmation:** bit-exact residual = 0 across 9 pure-Π²-even
non-truly pairs (YZ, ZY, XX+YZ, XX+ZY, YY+YZ, YY+ZY, YZ+ZY, YZ+ZZ, ZY+ZZ) at
N=3, 4, 5; 15 random non-uniform-J instances + 9 asymmetric J_YZ ≠ J_ZY
instances; pure D[Z]^⊗N dissipator. Reproduction:
`simulations/_f108_part1_pi_family_scan.py` plus
`simulations/_f108_part1_proof_algebra.py` for the 2-qubit anti-commutation
verification.

**Siblings:**

- **F108 Part 2 (BitA twin, X-dephasing, Tier 1 derived 2026-05-25):** the
  X-dephasing analog via the I↔Z, X↔Y phase-variant of Π_5bilinear. Same
  proof structure as Part 1, restricted to the Π²_X-even bilinear set
  {ZZ, XX, XY, YX, YY}. Per-site map: I → +Z, Z → −I, X → −iY, Y → +iX;
  per-site M² = diag(−1, +1, +1, −1) on {I, X, Y, Z}. Closes the X-dephasing
  branch of F109 Step 5. F108 Part 1's BitATwin slot points at this Claim
  (status `Filled`).
- **F108 Part 3 (Y-dephasing sibling, Tier 1 derived 2026-05-25):** the
  Y-dephasing analog. Same I↔X, Y↔Z permutation as Part 1 (Y-deph and Z-deph
  share bit_b parity per `PiOperator.SquaredEigenvalue`); per-site map
  I → +X, X → −I, Y → −iZ, Z → +iY differs from Part 1's only in the Y/Z
  2-cycle phase (−i vs +i, matching Y-deph's canonical Π convention). Same
  Π²-even bilinear set {XX, YY, YZ, ZY, ZZ} as Part 1; same M² sign pattern
  diag(−1, −1, +1, +1). Closes the Y-dephasing branch of F109 Step 5,
  promoting F109 to fully unconditional Tier1Derived across {Z, X, Y}.
  Part 3 is BitB axis (shares bit_b with Part 1); BitATwin slot is
  `BitBSpecific` (Y-deph has no meaningful bit_a analog).
- **F108 Klein-V₄ equivalence (Welle 14, 2026-05-27):** Parts 2 and 3 are
  Klein-V₄ corollaries of Part 1 via two complementary mechanisms.
  Part 1 ↔ Part 3 via operator-space D-conjugation (D · Π_5b(Z) · D = Π_5b(Y)
  bit-exact at N = 1, 2, 3; bilinear set fixed on bit_b axis). Part 1 ↔ Part 2
  via Hilbert-space Hadamard transport (U_op = U_H^⊗N ⊗ (U_H^⊗N)^* maps
  L_Z → L_X bit-exact; per-letter Hadamard bijects Part-1 bilinear set onto
  Part-2 bilinear set). NEGATIVE on operator-space Klein-V₄ for Π_5b:
  Q_zx · Π_5b(Z) · Q_zx ≠ ±Π_5b(X) and H · Π_5b(Y) · H ≠ ±Π_5b(X) (residual
  2.0 in Frobenius distance at N = 1, 2, 3); the operator-space Klein-V₄
  action on Π_5b is only the {I, D} subgroup; X-deph enters via Hilbert-space
  Hadamard. The three typed Claims are KEPT SEPARATE to preserve independent
  integration edges but cross-reference the equivalence proof.
- **F110 (HardCellYInversionPattern, Tier1Candidate, typed 2026-05-25):**
  hard cells y_par-asymmetric with Y-inversion. Aspect A closed-form via
  F108 Part 1+2+3 + F107 + F109 + F87 dissipator-resonance; Aspect B+C
  empirically anchored at F103/F105/F106.
- **F112 (LindbladBitBPiBalance, Tier1Derived for Hermitian H, typed 2026-05-26):**
  sibling derived theorem on the shared bit_b Z₂-grading. Where F108 Parts
  1/2/3 close palindromy for Π²-even *bilinears* (the bit_b = 0 sub-sector
  of 2-body H), F112 closes Π-eigenvalue +i / −i balance for arbitrary H
  with bit_b-homogeneous dissipator c. F108 and F112 are two independent
  Tier1Derived projections of the same F38 / F63 foundation: F108 lives in
  spec(L) palindromy, F112 in M_anti's Π ±i Frobenius split.

**Source:** [Proof Part 1](proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md),
[Proof Part 2](proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md),
[Proof Part 3](proofs/PROOF_F108_PART3_PI2Y_EVEN_ALWAYS_PALINDROMIC.md),
[Klein-V₄ equivalence (Welle 14)](proofs/PROOF_F108_KLEIN_V4_EQUIVALENCE.md);
`compute/RCPsiSquared.Core/Symmetry/F108Part1Pi2EvenAlwaysPalindromic.cs`,
`F108Part2Pi2XEvenAlwaysPalindromic.cs`, `F108Part3Pi2YEvenAlwaysPalindromic.cs`;
operator: `compute/RCPsiSquared.Core/Symmetry/Pi5BilinearOperator.cs` (supports
Z, X, Y dephasing via dephase parameter);
catalog parent: [NON_HEISENBERG_PALINDROME](../experiments/NON_HEISENBERG_PALINDROME.md)
(P1/P4/alternating/continuous Π-family taxonomy, all local); helpers (Parts 1+3 share
the same Π²-even bilinear predicates):
`F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(letter1, letter2)`,
`IsPi2EvenBilinearTerm(term)`, `IsPi2EvenBilinearHamiltonian(terms)`. Part 2
(BitA axis, bit_a=0 bilinears):
`F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(letter1, letter2)`,
`IsPi2XEvenBilinearTerm(term)`, `IsPi2XEvenBilinearHamiltonian(terms)`.

---

### F109. Mother Sector Soft is y_par = 1 Pure (Tier 1 derived, fully unconditional after F108 Part 1+2+3 closure 2026-05-25)

Sister to F107 on the y_par axis. F107 pinned the y_par signature of truly cells
across all dephase letters; F109 pins the y_par signature of mother sector
(Klein (0, 0)) soft cells. Previously Tier 1 derived modulo F108 Part 1; on
2026-05-25 all three branches of Step 5 (Z-dephasing via F108 Part 1, X-dephasing
via F108 Part 2, Y-dephasing via F108 Part 3) were closed-form via the matching
Π_5bilinear dephase variants. F109 is now fully unconditional Tier 1 derived.

**Theorem (F109):** Under any single-letter dephase channel (Z, X, or Y), every
Pauli pair classified as soft and located in the Mother sector Klein (0, 0) has
shared y_par = 1.

**Derivation chain:**

1. Klein (0, 0) constraint (bit_a = 0 AND bit_b = 0) forces all three letter
   counts #X, #Y, #Z to share the same parity (all even or all odd).
2. Per F107 per-dephase truly criteria, combined with Step 1's same-parity:
   Klein (0, 0) truly = all three counts even (y_par = 0).
3. Klein (0, 0) non-truly = all three counts odd (y_par = 1).
4. Klein (0, 0) is Π²-EVEN under every dephase letter (bit_b = 0 for Z/Y;
   bit_a = 0 for X).
5. Π²-even non-truly pairs are soft (not hard): closed-form across all three
   dephase letters via the F108 Part 1+2+3 family. Z-dephasing per F108 Part 1
   via Π_5bilinear; X-dephasing per F108 Part 2 via the X-deph variant (typed
   as F108 Part 1's BitATwin slot, `Filled`); Y-dephasing per F108 Part 3 via
   the Y-deph variant (BitB-axis sibling of Part 1 on the same Π²_Z axis;
   `BitBSpecific` BitATwin slot). All three parts closed 2026-05-25.
6. Klein (0, 0) soft term ⟹ y_par = 1; y_par-homogeneous pair: shared y_par = 1.

**Empirical confirmation:** F103 mother soft (0, 21) × 3 dephase; F105 same;
F106 (0, 300) × 3. Total 1026 mother-soft classifications, all y_par = 1, zero
y_par = 0. F109 explains this bit-exactly across all three dephase letters via
the F108 Part 1+2+3 closure.

**Cross-letter spot-check:** Klein (0, 0) non-truly k=3 terms are the 6
XYZ-permutations (only triple with all-odd and sum ≤ 3). Unordered pairs with
self: 6·7/2 = 21 (matches F103/F105). At k=4: 24 letter sequences (3 non-I + 1 I),
24·25/2 = 300 (matches F106).

**Source:** [Proof](proofs/PROOF_F109_MOTHER_SOFT_Y_PARITY_ONE_PURITY.md);
`compute/RCPsiSquared.Core/Symmetry/MotherSoftYParityOnePurity.cs`;
parents: PROOF_F107 + PROOF_F85 + PROOF_F108_PART1; helpers:
`MotherSoftYParityOnePurity.IsMotherNonTrulyCandidate(term)` and
`MotherSoftYParityOnePurity.VerifyOnTerm(term)`.

---

### F110. F87-Hard Cells Exhibit Y-Inversion Pattern (Tier 1 derived, promoted 2026-06-10; Aspect A closed-form via F108 Part 1+2+3 + F107 + F109 + F87 dissipator-resonance, Aspect B+C derived via F103 §6/§7 + the closed windowed converse)

Seventh YParity-axis Claim; completes the y_par-axis classification of the F87
trichotomy together with F107 (truly y_par=0) and F109 (mother soft y_par=1).

**Aspect A (closed-form):** F87-hard pairs appear only in the diagonal Klein cell
matching the dephase letter (Z → (0, 1), X → (1, 0), Y → (1, 1)). Derivation:
F108 Part 1+2+3 close Π²-D-even cells (never hard); F107 + F109 close Mother
sector Klein (0, 0); the F87 dissipator-resonance law (Tier1Derived, anchored at
N=4 k=3 over 294 pairs in `DissipatorResonanceLaw.cs`) selects the one diagonal
cell among the two remaining Π²-D-odd non-mother cells. Combining all three
closures: hard appears only in the diagonal cell.

**Aspect B (Y-inversion):** Within the diagonal hard cell, the dominant y_par
equals y_par(dephase letter): Z/X-deph dominantly y_par=0, Y-deph dominantly
y_par=1. At k = N = 4 the dominance is closed-form Tier1Candidate via the
sibling Claim F111 (HardCellPureDTemplate, 2026-05-25): every F87-hard pair in
the diagonal cell contains at least one "pure-D template" (only D and I
letters), and pure-D templates have y_par = y_par(D) by construction. At k = 3
the 42:8 dominance follows from the §6 counting rule (see Aspect C).

**Aspect C (k-purity sharpening, §6 closed-form counting rule):** k=3 N=4 (F103): 42:8 biased; k=3
N=5 (F105): identical 42:8 (N-stable); k=4 N=4 (F106): 228:0 fully pure with
Y-inversion preserved.

**Status (promoted Tier1Derived 2026-06-10):** the 42:8 (k=3) and 228:0 (k=4) ratios are derived by
the §6 counting rule and the §7 bipartite mechanism; F111 closes the k=4 case at full support. The
windowed (k<N) hard-direction converse (non-bipartite ⟹ hard) is a **theorem with no residual**
(`WindowedConverseAllGammaClaim`, Tier1Derived: R-deg retired by the girth dichotomy, R-sign resolved
by the Pascal-Gram positivity theorem F117, both 2026-06-10), which was the F110/F111 promotion gate.
**Open:** k ≥ 5 empirical confirmation; QPU confirmation at k ≥ 3.

**Source:** [Proof](proofs/PROOF_F110_HARD_CELL_Y_INVERSION.md);
`compute/RCPsiSquared.Core/Symmetry/HardCellYInversionPattern.cs`;
parents: PROOF_F108_PART1/2/3 + PROOF_F107 + PROOF_F109 + F87 dissipator-resonance
law (Tier1Derived, `compute/RCPsiSquared.Diagnostics/F87/DissipatorResonanceLaw.cs`);
helpers: `HardCellYInversionPattern.DiagonalKleinCellForDephase(dephase)`,
`IsDiagonalCell(klein, dephase)`,
`DominantYParityForDephase(dephase)`.

---

### F111. Hard Cell Pure-D Template Rule (Tier 1 derived, promoted 2026-06-10; closed-form structural rule that implies F110 Aspect B at k = N = 4)

Eighth YParity-axis Claim; sharpens F110 Aspect B by replacing "Y-inversion at
k = 4" with a structural rule that implies it as immediate corollary AND
provides the per-cell decomposition 36 + 192 + 0 = 228 across all 3 dephases.

**Theorem (F111):** At k = N = 4 in the diagonal Klein cell (D.BitA(), D.BitB())
for dephase letter D ∈ {Z, X, Y}, a Pauli pair (P, Q) is F87-hard if and only if
at least one of P, Q is a "pure-D template" (a length-4 Pauli string composed
only of D and I letters, no other non-I Pauli letter).

**Corollary (F110 Aspect B at k = 4):** Pure-D templates have y_par = y_par(D)
by construction. Therefore every F87-hard pair at k = N = 4 in the diagonal cell
has y_par(pair) = y_par(D). This is the F106 228:0 split across all 3 dephases.

**Structural decomposition (per diagonal cell at k = N = 4):** 528 total pairs
split into 36 Pure-Pure (all HARD) + 192 Pure-Mixed (all HARD) + 300 Mixed-Mixed
(all SOFT). The 36 + 192 + 0 = 228 hard count matches F106 exactly.

**Subclaim status:**
- (a) Pure-D single-term H is HARD: heuristic via dissipator-commute ([D[D_l], L_H] = 0
  for pure-D H gives additive independent spectra L = L_H + L_D, combined
  spectrum non-palindromic around −σ). Empirical: 8/8 per cell.
- (b) Mixed single-term H is SOFT: empirical (24/24 per cell).
- (c) Pair (Pure-D, Mixed) H is HARD: empirical (192/192 per cell).
- (d) Pair (Mixed, Mixed) H is SOFT: **CLOSED modulo M via PROOF_F103 §7.4
  (2026-05-30)**. At full support k = N a Mixed+Mixed pair has at most two flip
  generators, which always admit a linear φ and hence the chiral K, so the
  hopping graph is bipartite and the pair soft. Empirical: 300/300 per cell.
  (The earlier operator-level search, Task 1 Paths 1-3 per-site M^N tensor
  product / F108 Π_5bilinear extended action / Q_V × Π composition, was
  dissolved by the chiral route.)

**Empirical anchor:** F106 N = 4 k = 4 records 228:0 in all 3 diagonal cells
(684 hard pairs total, zero off-y_par). Independent verification at
`simulations/_f111_pair_off_ypar_verify.py` covers 1584 pair classifications
across 3 dephases, all matching the rule with zero exceptions.

**Promotion record (2026-06-10):** the hard-direction converse behind subclaims (a)/(c)
(non-bipartite ⟹ hard at every γ) closed as the windowed all-γ theorem with no residual
(`WindowedConverseAllGammaClaim`, Tier1Derived: girth dichotomy + Pascal-Gram positivity F117), so
F111 is Tier1Derived; subclaim (d) was already closed modulo M via PROOF_F103 §7.4 (see above).
**Open:** Pure-D Template Rule at k > 4 or N > 4 (empirically unverified). Hardware QPU
confirmation at k ≥ 3.

**Source:** [Proof](proofs/PROOF_F111_HARD_CELL_PURE_D_TEMPLATE.md);
`compute/RCPsiSquared.Core/Symmetry/HardCellPureDTemplate.cs`;
parents: PROOF_F110 + PROOF_F107 + PROOF_F108_PART1/2/3 + F106 anchor +
F87 dissipator-resonance law (Tier1Derived);
helpers: `HardCellPureDTemplate.IsPureDTemplate(term, dephase)`,
`IsInDiagonalCellAtK4N4(p, q, dephase)`,
`IsPredictedHardAtK4N4(p, q, dephase)`,
`VerifyYInversionCorollaryAtK4N4(p, q, dephase)`.

---

### F112. Lindblad Π-Eigenvalue Balance under bit_b Homogeneity (Tier 1 derived for both Hermitian and non-Hermitian H, universal N)

Structural identity making the `polarity_coordinates_from_L` diagnostic
asymmetry an exact witness for c outside the bit_b-homogeneous regime. Sits on
the bit_b Z₂-grading of the Pauli group shared with F108 (Π²-Z-even bilinear
palindrome closure) and F87 (dissipator-resonance trichotomy).

**Theorem (F112, Hermitian H, rigorous):** For any Lindblad-form Liouvillian
L = -i[H, ·] + Σ_k γ_k · `np.kron(c_k, c_k^*)` on N qubits with Hermitian H
and each c_k bit_b-homogeneous (every Pauli string σ in c_k's expansion
shares bit_b(σ) = (#Y(σ) + #Z(σ)) mod 2 = const), the
`polarity_coordinates_from_L` decomposition of M = Π L Π⁻¹ + L + 2σ·I
satisfies

    ‖M_plus_half‖² = ‖M_minus_half‖²

bit-exactly (machine precision), for any choice of complex coefficients γ_k
and any Pauli-coefficient choice in each c_k.

**Non-Hermitian extension (Tier 1 derived, universal N, Welle 11, 2026-05-27):**
the equality also holds for non-Hermitian H. Writing H = H_re + i H_im with both
summands Hermitian, the equality reduces algebraically to the open identity
Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 for any Hermitian H_re, H_im. This identity is
now closed structurally for all N via the Welle 11 two-lemma proof
(Diagonal-Norm + Off-Diagonal-Orthogonality) in
`docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md`, making the non-Hermitian
extension Tier1Derived universal N. The 559,912-pair basis enumeration at
N ≤ 5 (Welle 10a Python + Welle 10b C# SLOW_F112) is preserved as the
empirical anchor that motivated the search for the structural proof; it
remains the historical numerical validation.

**Five-step structure:**
- Step 1: reduce balance to Π-conjugation ±i Frobenius equality:
  asymmetry = (1/2) (‖M_{+i}‖² − ‖M_{-i}‖²).
- Step 2: bit_b-homogeneous c implies `np.kron(c, c.conj())` lies entirely in
  the Π²-conjugation +1 eigenspace (via F38 / F63 Π² eigenvalue formula on
  Pauli strings).
- Step 3: Π²-conj +1 eigenspace = Π-conj {+1, −1}; dissipator has zero +i,
  −i content.
- Step 4: M_{+i} and M_{-i} come entirely from L_H, with norms
  2 · ‖L_{H,±i}‖².
- Step 5 (Hermitian H): L_H^† = −L_H (anti-Hermitian as superoperator) plus
  dagger maps Π +i ↔ Π −i bijectively while preserving Frobenius. Combining
  gives ‖L_{H,+i}‖² = ‖L_{H,-i}‖².

**Connection axes (shared bit_b Z₂-grading on the Pauli group):**
- F38: Π² = (−1)^{w_YZ} on Pauli strings; foundational input.
- F63: [L, Π²] = 0 for Z-dephasing; foundational input.
- F108 Part 1/2/3: the bilinear set {XX, YY, YZ, ZY, ZZ} F108 palindromizes
  is exactly the bit_b = 0 (Π²-Z-even) family.
- F87 dissipator-resonance trichotomy: orthogonal axis, derived 2026-06-10
  (previously empirical via `_polarity_probe_f87_connection.py`). F87 lives
  in M's spectrum-palindrome structure; F112 lives in M_anti's Π +i/−i
  split. Three-part derivation: (a) scope inclusion (every F87 input,
  Hermitian Pauli H + pure Z-deph with single-Pauli c = Z_l, satisfies
  F112's hypotheses, so the asymmetry is identically zero on F87's entire
  domain, all three classes; asym = 0.0 exact float zero at N = 3, 4);
  (b) mechanism separation (on bit_b-odd H, the diagonal Klein cell, the
  dagger involution IS the windowed converse's first reflection,
  M_rec† = 𝓕 M_rec 𝓕 with 𝓕 = X^⊗N ⊗ X^⊗N, diff 0.00e+00; F112 reads it
  at degree 2 via Frobenius norms of Π-eigenprojections, F87 hardness at
  odd degree via the second reflection R + unsigned girth); (c) scoped
  one-way F113 bridge on the σ⁻/σ⁺ family: both functionals read the shared
  moment t₁^(l) = Tr(Z_l H) = 2^N c_l = 2^(N−1) ω_l, F113 linearly
  (asym = 2^N · Σ_l t₁^(l) · (γ_pump,l − γ_T1,l); machine-precision match
  |diff| = 7.1e-15 at N=3 after the ω_l = 2c_l + σ⁻-lowering convention
  reconciliation), the girth ladder's ℓ=1 face quadratically
  (p₃ = 6γ · Σ_l (t₁^(l))²), so balance-broken ⟹ F87-hard, one-way only.
  See the dated section in the proof; committed verifier
  `simulations/f112_f87_orthogonality.py`.

**Empirical anchor:** 14 probes (`simulations/_polarity_probe_*.py`,
`_polarity_proof_verify.py`, `_polarity_step5_stress.py`) cover
candidate-breakers (1-5), hand-engineered non-Lindblad L (6), random c with
full Pauli rank (7-8), k_max boundary (9), exhaustive 136-pair N=2 enumeration
(10), coefficient sweep (11), Z₂³-cell N=3, 4 scaling (12, 171 / 171
balanced within-cell), Π²-content verification (13, 100.00% Π²=+1 for
bit_b-homogeneous c), and direct Π-eigenspace L_H projection across 30
random H (10 Hermitian + 10 non-Hermitian Pauli + 10 random complex matrix)
at N = 2, 3 (14, all bit-exact).

**Open:** connection to F104 / F105 / F106 (F87 Z₂³-cubed refinements).
(Step 5 extension to non-Hermitian H closed Welle 11, 2026-05-27; see
`docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md`. Structural derivation
of the F87 ↔ F112 orthogonality closed 2026-06-10: scope inclusion +
mechanism separation + the scoped F113 one-way bridge, dated section "The
F87 orthogonality, derived (2026-06-10)" in the parent proof; committed
verifier `simulations/f112_f87_orthogonality.py`.)

**Source:** [Proof (parent, Hermitian H)](proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md);
[Proof (Welle 11, non-Hermitian extension, universal N)](proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md);
parents: F38 + F63 + F108 Part 1/2/3 (Tier 1 derived) + F87 dissipator-
resonance law + `polarity_coordinates_from_L` primitive
(`simulations/framework/diagnostics/polarity_coordinates.py`, added 2026-05-25).

---

### F113. F112 Counterexample Asymmetry Closed Form (Tier 1 derived for general N via Welle 4)

Closed-form magnitude for the F112 polarity asymmetry in the regime where
F112's typed scope is violated by single-site Z-drive × amplitude-damping
interference. Discovered via Welle 2 Kingston hardware analysis (commit
`a1a90a2`, 2026-05-26), derived constructively the same day via parameter-
sweep regression (commit referencing this entry).

**Theorem (F113):** For Lindblad-form L = -i[H, ·] + Σ_k γ_k · D[c_k] with:
- Hermitian H containing single-site Z-drives Σ_l (ω_l/2)·Z_l plus any
  bit_b-homogeneous additions (X-drive, Y-drive, ZZ/XX/YY/XY bonds: each
  contributes 0 individually by F112)
- Dissipator c containing σ⁻_l at rate γ_T1,l and σ⁺_l at rate γ_pump,l
  (amplitude damping / pumping)
- Optional bit_b-homogeneous Z-dephasing (in F112's typed scope, contributes 0)

the F112 polarity asymmetry has the closed form

    asymmetry := ‖M_plus_half‖² − ‖M_minus_half‖²
              = (4^N / 2) · Σ_l ω_l · (γ_pump,l − γ_T1,l)

bit-exactly. Verified at N = 2, 3, 4 via parameter sweep (`simulations/
_f113_break_formula_derivation.py`); per-site decomposition, cross-site
zero, sign flip on ω → −ω and on σ⁻ ↔ σ⁺, detailed-balance cancellation
(γ_T1 = γ_pump → 0), and non-uniform-rate sum-formula all confirmed
bit-exact.

**Structural origin:** F112 break requires non-Hermitian Π-eigenspace
coupling between H and c. Only the Z-drive commutator produces this:
[Z, σ⁻] = −2σ⁻ (proportional to the non-Hermitian σ⁻ itself), while
[X, σ⁻] = Z and [Y, σ⁻] = i·Z give Hermitian commutators that remain
F112-symmetric. ZZ/XX/YY bilinears commute differently and produce 0
contribution. Same-site locality of [Z_l, σ⁻_m] = −2σ⁻_m · δ_{lm} gives
the per-site additive structure.

**Universal-N status (Welle 4, 2026-05-26):** Tier 1 derived for general
N via [PROOF_F113_COEFFICIENT_DERIVATION](proofs/PROOF_F113_COEFFICIENT_DERIVATION.md).
The (1/2)·4^N coefficient decomposes structurally as
4 · 4^(N-1) · (1/2):
- factor 4: Welle-4 reduction `asymmetry = 4 · Re⟨L_H,+i, L_T1,+i⟩`
  (from norm² expansion + F112 typed + F112 non-Hermitian extension
  cancellations + cross-term equal-magnitude-opposite-sign relation).
- factor 4^(N-1): N−1 spectator-site identity factors, each contributing
  ⟨I_4, I_4⟩ = Tr(I_4) = 4 to the Frobenius inner product on tensor
  products. Operator-space dimension entering via local Pauli dimension 4
  per spectator site.
- factor 1/2: single-site N=1 inner product `⟨(L_H,1)_{+i}, (L_T1,1)_{+i}⟩ = −ωγ/2`,
  derived from explicit 4×4 sympy matrices.
The proof's Lemma C step 5 has one specific Frobenius equality that is
bit-exact verified at N = 1, 2, 3, 4, 5 but not yet closed from the
support-pattern algebra alone; this is documented as a structural
exercise and does not block the universal-N Tier1Derived status given
the bit-exact anchor.

**Hardware fingerprinting application:** asymmetry measurement directly
extracts Σ_l ω_l · (γ_pump,l − γ_T1,l) when drive parameters are known;
becomes a per-site amplitude-damping calibration tool when combined with
ω_l knowledge. The Welle 2 hardware-fit value for f95 (ω=0.13, γ_T1≈0.001,
γ_pump=0, N=2) gives F113-predicted (1/2)·4² · (2·0.13) · (0 − 0.001) = −2.08e−3
(coefficient (1/2)·4² = 8 at N=2; drive on both q13, q14, so Σ_l ω_l = 2·0.13 — the
earlier-notation "16" folds the per-site coefficient 8 together with the two sites),
matching the fitted value bit-exact (the rel asymmetry sign tracked
correctly through both the Python derivation script and the C# pipeline
after the 2026-05-26 convention reconciliation).

**Sister claims:**
- F112 (Tier 1 derived Hermitian H + bit_b-homogeneous c → asymmetry = 0)
  covers the in-scope half of the standard-Lindblad-family polarity
  behavior; F113 covers the bit_b-mixed-c break magnitude. Together they
  give a complete picture across the family.
- F84 amplitude-damping correction (T1 violation of F81 identity):
  F113's σ⁻ contribution is the polarity-axis projection of F84's
  F81-axis violation.

**Source:** [F113 derivation + verification](../experiments/F113_BREAK_MAGNITUDE_FORMULA.md);
parents: F112 (`LindbladBitBPiBalance`, Tier 1 derived) + the Welle 2
hardware counterexample (`experiments/F112_HARDWARE_LENS_KINGSTON.md`);
script: `simulations/_f113_break_formula_derivation.py`.

---

### F114. Commutator Superoperator D-Conjugation Parity (Tier 1 derived, closed form bit-exact N = 1..4)

Closed-form sign functional ε(σ) for the action of D-conjugation on the
H-commutator superoperator L_σ = −i[σ, ·] in the 4^N Pauli basis, where
D = diag((−1)^{n_Y(α)}) is the real diagonal unitary involution that lifts
the Z↔Y dephase-letter swap to operator space (Welle 12: Π_Y = D·Π_Z·D).

**Theorem (F114):** For any single Pauli string σ ≠ I^{⊗N} on N qubits,

    D · L_σ · D = ε(σ) · L_σ    bit-exact

with closed form

    ε(σ) = (−1)^{n_Y(σ) + 1}
         = +1 if n_Y(σ) is odd
         = −1 if n_Y(σ) is even and σ ≠ I^{⊗N}

For σ = I^{⊗N}: L_σ = 0, sign undefined (vacuous case).

For H = Σ_k c_k σ_k a linear combination of Pauli strings: ε(H) is well-
defined and equals ε(σ_k) iff all σ_k share the same n_Y parity. If the
terms split across both parity classes, D-conjugation yields a mixed (non-
multiplicative) response on L_H and no single ε(H) sign exists.

**Why it surfaced:** Welle 15 Task A polish (2026-05-27, commit `a98fc02`)
observed bit-exact at N = 2 that for XZ + ZX bond (n_Y per term = 0)
M_anti(L, Π_Y) = −D · M_anti(L, Π_Z) · D, while YZ + ZY (n_Y per term = 1)
gives +D · M_anti(L, Π_Z) · D. The bond-specific sign motivated systematic
enumeration; the n_Y-parity closed form was identified by
`simulations/_m_level_sign_functional_explore.py` (2026-05-27) and verified
bit-exact for all 4^N Pauli strings at N = 1, 2, 3 (84 strings total) plus
the Welle 15 bond bilinears at N = 2 (12 cases) and selected multi-bond /
multi-body cases at N = 3, 4 (6 cases).

**Operational reading:** D acts on L_σ as an automorphism with sign
determined by σ's Y-letter content. Y is the unique Pauli letter that
anticommutes with the antisymmetric transposition convention (σ_Y^T = −σ_Y
while σ_I^T = σ_I, σ_X^T = σ_X, σ_Z^T = σ_Z); the per-letter parity bookkeeping
collapses to a per-string n_Y parity.

**Consequence for F112 / F108 cross-dephase (Welle 13):** the F112-Y M-level
identity

    M(L_H, Π_Y) = ε(H) · D · M(L_H, Π_Z) · D    bit-exact

holds when the dissipator contribution to M vanishes (F112 hypothesis: Hermitian H
plus bit_b-homogeneous c) AND when ε(H) is well-defined (all H Pauli terms share
n_Y parity). This refines the Welle 13 PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4
statement "L_Y is not D-transportable" by exhibiting the precise M-level
ε-signed equivariance that survives despite the L-level absence of D-transport.
The F112 typed scope (norm-level ‖M_+1/2‖² = ‖M_−1/2‖²) remains sign-invariant.

**Connection axes:**
- Welle 12 D · Π_Z · D = Π_Y: F114 is the L-level companion. Welle 12 makes
  D the Π swap-operator across {Z, Y} dephase letters; F114 makes D the L_H
  sign-flip-operator with per-term n_Y bookkeeping.
- F112: ε is sign-invariant under the F112 norm-level statement, so F112-Y
  Tier1Derived (via Welle 13 Route 1) and F112 are mutually consistent.
- F108 Part 1+3: Z- and Y-dephase palindrome closures share the Π² bit_b
  grading; F114's ε(H) characterization gives the precise H-side condition
  under which the M_anti structure is signed-equivariant across {Z, Y}.

**Empirical anchor:** `simulations/_m_level_sign_functional_explore.py`.
PART 1 enumerates ε(σ) for all 4 + 16 + 64 = 84 single Pauli strings at
N = 1, 2, 3 and verifies the n_Y-parity closed form bit-exact. PART 2
verifies ε(H) for 12 bilinear bond Hamiltonians at N = 2 (XZ+ZX, YZ+ZY,
Heisenberg, single-bond XX/YY/ZZ/XY/YX/ZY/XZ, mixed-letter single-site
combos). PART 3 sanity-checks 6 multi-bond / multi-body cases at N = 3, 4
(chain bond pairs, ZZZ, YYY, XYZ, N=4 Heisenberg chain). All ε predictions
match the actual D-conjugation residual bit-exact (residual = 0.00e+00
numpy double precision).

**Typed Claim:** promoted same wave (2026-05-27, commit `35a88df`) to
`CommutatorDConjugationSign` (Tier1Derived) in
`compute/RCPsiSquared.Core/Symmetry/`, ctor parent
`Pi2KleinV4DephaseSwapGroup` (Welle 12). Wired in
`KnowledgeRegistryFactory.BuildDefault()`; cross-referenced from
`LindbladBitBPiBalance` (F112-Z) and `LindbladBitBPiYBalance` (F112-Y)
via InspectableNode. 35 unit tests including matrix-level `D·L·D = ε·L`
bit-exact verification at N = 2.

**Closed general-N (2026-06-10, the transpose reading).** D is the transpose
superoperator θ(ρ) = ρᵀ written in the Pauli basis, so
D·L_σ·D = θ∘L_σ∘θ = −L_{σᵀ} = (−1)^{n_Y(σ)+1}·L_σ for every σ and every N
in one line (σᵀ = (−1)^{n_Y}σ; Y is the only antisymmetric Pauli). The
girth-ladder reversal kill
([PROOF_F87_WINDOWED_MONOMIAL_CONVERSE](proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §4:
word reversal = transpose × (−1)^{n_Y(word)}) is the same antiautomorphism
at word length j; F113 Lemma C is its Hermitian-conjugacy sibling.
Bit-exact anchor: `simulations/_mirror_inventory_bridge_check.py`
(blocks D/E: 63/63 strings at N = 3 plus an N = 5 case, dev 0.00e+00).

**Open:** N = 5, 6 verification at scale (estimated O(4^N) per N for the
single-string sweep, tractable but not run). *Superseded 2026-06-10:* the
transpose reading above derives ε(σ) structurally for every σ and every N,
and the `_mirror_inventory_bridge_check.py` anchor includes an N = 5 case
bit-exact; the per-N sweep is closed. Connection to the Welle 11
F112 Lemma A/B structural proof: does F114 give an alternative derivation
of the parent F112 Step 5 (Lemma B) via D-conjugation parity rather than
dagger anti-Hermiticity?

**Source:** `simulations/_m_level_sign_functional_explore.py`;
parents: Welle 12 D · Π_Z · D = Π_Y identity
(`compute/RCPsiSquared.Core/Symmetry/Pi2KleinV4DephaseSwapGroup.cs`;
`docs/proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md`) + F112 typed scope
(`docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md`,
`docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md`).

---

### F115. Windowed F87 hardness is a closed-form combinatorial theory (Tier 1 derived, GF(2)[x] valuation)

F103 §7 derives the F87 diagonal-cell soft/hard split as a bipartite-chirality mechanism; F115
collapses its combinatorics to closed form. Read a k-body diagonal-cell Mixed term's X/Y positions as
a polynomial p(x) over GF(2) (bit j ↦ x^j), even popcount so (1 + x) | p; the sliding-window builder
places it at windows w = 0 … N−k as the shifts {x^w p}, with W = N − k + 1 windows. Four statements,
all bit-exact.

1. *Hardness is a valuation difference (the one-number criterion).* A Z-dephasing diagonal-cell Mixed
   pair is hard ⟺ v_{1+x}(p₁) ≠ v_{1+x}(p₂), the (1 + x)-adic valuations of the two masks; equal
   valuations ⟹ soft. This is the §7.5 non-bipartite criterion in one subtraction, and it matches the
   actual spectral trichotomy on every k = 3, N = 4 pair (derived for any k via §7.7 + §7.5/§7.6).
   Scope of the Tier-1-derived label: the criterion plus the genericity result (hard for all but
   finitely many γ); the all-γ closure is `WindowedConverseAllGammaClaim` (Tier1Derived since
   2026-06-10, no residual: girth dichotomy + Pascal-Gram positivity F117,
   [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md](proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md)).

2. *Obstruction-size law (two-layered).* The minimal odd 𝔽₂-relation among the shifts (the hardness
   obstruction) has maximal size over hard pairs = min(2W − 1, 2k − 3 − 2d), where d = deg(g_rest) is the
   degree of the shared non-(1 + x) part of g = gcd(p₁, p₂). The window leg 2W − 1 (at most 2W masks, an
   odd subset ≤ 2W − 1); the body leg 2k − 3 − 2d ((1 + x) | p forces the gcd-quotients to degree
   ≤ k − 2 − d, popcount ≤ k − 1 − d). Only (1 + x) governs hard/soft (the unique hardness prime, §7.9);
   the other shared factors shrink the obstruction by 2 per degree. The d = 0 face is 2k − 3; k = 3 gives
   3, the always-triangle case.

3. *Closed-form hard count.* The even-popcount nonzero k-bit masks split by valuation into classes of
   size c_v = 2^{k-1-v} (v = 1 … k−1); a pair is hard iff its masks differ in class, so
   #hard mask-pairs = e₂(2^0, …, 2^{k-2}) = (4^{k-1} − 3·2^{k-1} + 2)/3 = OEIS A203241. Dressed by the
   2^{2k-3} Klein / y-parity factor this is the hard count itself (448, 8960, 158720 at k = 4, 5, 6).
   The count is itself layered by the shared-factor degree d of §7.9: #hard with deg(g_rest) = d =
   2^{d-1}·B(k − d), B(k) = (4^k − 12k + 8)/18 the d = 0 base. The size-3 (triangle) sub-class closes,
   5·2^{k-1} − (3k² + k)/2 − 3; the exact per-size middle counts stay window-dependent.

4. *Coding-theory home.* The obstruction = the minimum ODD weight of the two-generator quasi-cyclic
   (terminated rate-1/2 convolutional) code generated by (a, b) = (p₂/g, p₁/g), g = gcd. The code's free
   distance (any parity) is a constant 4; the palindrome reads the parity-restricted odd floor 2k − 3.
   The minimum-odd-weight invariant is non-standard (coding theory optimizes free distance), so this is
   the correct vocabulary and home, not an imported theorem: the size law and the count are elementary
   (1 + x)-valuation facts. The one genuine coding-theory effect is cancellation (a sparse multiple a·s
   drops the obstruction below the gcd-generator popcount at k ≥ 6), which is what makes the middle-class
   distribution window-dependent.

**Connection to the polarity cube: bit_a is the lowest valuation bit.** The (1 + x)-adic valuation grades
the existing Klein bit_a axis. bit_a = #(X + Y) mod 2 = popcount(mask) mod 2 is exactly [v_{1+x}(mask) = 0],
the lowest bit of the valuation (odd popcount ⟹ v = 0; even popcount ⟹ v ≥ 1). Inside the diagonal cell
(bit_a = 0, bit_b = 1) the valuation keeps grading, v = 1 … k − 1, into the very classes c_v = 2^{k-1-v}
of statement 3, and the hard/soft "Z bit" is the pair-difference v(p₁) ≠ v(p₂). So F115 is not a new
polarity axis: it is the graded refinement of the bit_a axis, read at the pair level inside the bit_b = 1
cell (bit_b selects that cell because the −N reflection mode needs F H F = −H, i.e. odd #Y + #Z). Verified
bit-exact: bit_a = [v = 0] for all masks, and the bit_a = 0 cell's valuation classes reproduce c_v.

The grade reads concretely as a tower of dyadic moments: moment_j is the parity of the X/Y count on the
sites whose index has j as a submask (j AND i == j), moment_0 being bit_a (the total parity), and the
grade is the number of leading vanishing moments, i.e. how deeply the X/Y pattern is dyadically balanced
((1 + x)^{2^m} is a distance-2^m pair, of grade 2^m). So the combinations the grade produces, the §7.8
classes, are the dyadic-balance strata of the X/Y pattern, not bookkeeping. The graded axis is pulled
into the toolkit as C# `BitADyadicGrade` (`compute/RCPsiSquared.Core/Symmetry/`, with `Moment` / `Grade`
/ `BitA`), tested against the GF(2)[x] valuation and the classes c_v.

**Verified:** the valuation criterion vs the actual spectral verdict (k = 3, N = 4, every pair); the
size law across a (k, N) grid through k = 6; the class sizes, count, d = 3 form, and
free-distance-4-vs-odd-weight-2k−3 through k = 10, with the dressed totals matched against the C# scan
at k = 4, 5, 6; the quasi-cyclic dictionary bit-exact k = 4, 5, 6.

**Source:** [Proof](proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md) §7.7–§7.9; C# `WindowedObstructionScan`
(`compute/RCPsiSquared.Diagnostics/F87/`, with the GF(2)[x] `PolyGcd` / `ValuationAtOnePlusX` /
`GcdFormulaSize` / `IsHardPair` helpers, and the closed-form count ports `HardMaskPairCount` (A203241) /
`HardCountBaseB` (B(k)) / `HardCountByGRestDegree` (2^(d-1)·B(k−d)) / `TriangleHardMaskCount` /
`MaxObstructionSizeForGRestDegree` (2k−3−2d) + `GRestDegree` / `EvenPopcountMasks`, ported 2026-06-08) +
`WindowedObstructionScanTests` and `WindowedHardnessCountClosedFormTests`; scouts
`simulations/_f87_obstruction_derivation.py`, `_f87_coding_theory_scout.py`,
`_f87_hardcount_closedform.py`, `_f87_beyond_x1_scout.py`, `_f87_size_second_layer.py`. The N-free hard
verdict is now wired into `PalindromeSoftCertifier.Decide` (`CertifyHardByDiagonalCellValuation`),
soundness-gated against the spectral authority over the diagonal-cell pair space (`PalindromeHardSweepTests`),
2026-06-08.

---

### F116. The Golden Ceiling Router (Tier 1 derived: constructive existence, exact arithmetic, every N ≥ 3)

The two Z-middle ceiling cases, H = Σ_windows (XZX+XZY+YZX) and the X↔Y sibling
YZY+XZY+YZX under Z-dephasing, are palindromized by a **per-site product** after all:
the period-4 golden router

  W = ⊗_l q_{l mod 4},   W L W⁻¹ = −L − 2σ,   σ = Σ_l γ_l (site-dependent rates allowed),

with the closed form (φ² = φ + 1, all entries in ℤ[φ]+iℤ[φ]):
g_l = q_l(I) follows [a, a, b, b] with a = φX + Y, b = X − φY; h_l = q_l(Z) =
(−1)^(l+1)·i·R(g_l); q_l(X) = −(g_l)_X·I + (h_l)_X·Z, q_l(Y) = −(g_l)_Y·I + (h_l)_Y·Z.
Each q_l is class-swapping with q_l² = −(2+φ)·I (a scalar times a unitary, cond(W) = 1).
The frame directions are the two roots of the **golden locus α² − αβ − β² = 0**
(slopes 1/φ and −φ; tan 2θ = 2). Two-sided form: W(ρ) = (2+φ)^(N/2)·P ρ Q with product
unitaries P, Q each anticommuting with H (so spec(H) is exactly ±E-symmetric, and
G = PQ is a weak Z₂ symmetry of L and an exact eigenmode at the palindrome floor −2σ).

**Mechanism:** the window-summed anticommutator {Q₃, [XZX+XZY+YZX, ·]₃} = 0 exactly at
every window offset (cross-template cancellation inside one window; per-term it fails,
which is why the per-term certifier search could never see it); window additivity gives
every N ≥ 3. **Exclusion side (derived):** the identity-column functional forces
[H, ⊗g_l] = 0 with per-window equations whose only uniform or period-2 solution is g = 0
(the committed optimization floors are now theorems), period 3 is impossible for N ≥ 5,
and the discrete Klein candidates P1/P4/M2/M sit off the locus (values +1, −1, −1, −1/2).
**Rigidity:** zero continuous moduli; the invertible solution set at N=5 is exhaustively
4 cyclic shifts × an explicit order-32 sign group, all golden. The ceiling arc closes
6 → 4 → 2 → **0**: no case in the k=3 windowed soft family needs a non-local mirror.
Open chains only (rings untested).

**Verified:** {W,A} = 0 and the window lemma EXACT over ℤ[φ]+iℤ[φ] at N = 3..6 (both
siblings); end-to-end vs the framework Lindbladian at N = 5, 6, γ ∈ {0.3, 0.7, 1.0} and
site-dependent rates (rel ~2e-16); independently re-implemented five ways (through N = 9
sampled) in the 2026-06-10 adversarial audit.

**The metallic family (2026-06-11):** the golden point is the c = 1 member of a
one-real-parameter line. For weighted templates t₁·XZX + t₂·XZY + t₃·YZX the soft set
is exactly **t₂ = t₃** (off-line is hard by the girth ladder; witness (1,2,1) fires at
m\* = 11 with p₁₁ = 1730150400·γ³ exact). On the line, with c = t₁/t₂, the same
[a, a, b, b] router works with a = (r, 1), b = (1, −r), **r(c) = (c + √(c²+4))/2 the
metallic mean** (golden c=1, silver c=2, bronze c=3; r(−c) = 1/r(c) so c=−1 is 1/φ;
c=0 is the 45° diagonal), q_l² = −(1+r²)·I. **Derived for ALL real c** by degree bound
+ exact interpolation: every entry of r·{Q₃,[T_c,·]₃} is a degree ≤ 5 polynomial in r
(c = r − 1/r), and exact Fraction arithmetic finds 8 rational nodes identically zero,
so the window lemma is a polynomial identity in r. The exclusion generalizes: the
identity-column determinant factors as **c·(α² − cαβ − β²)**, so the metallic locus
gates the frame for every c ≠ 0 (alternation forced, uniform/period-2 empty, zero
moduli verified at silver), while at **c = 0 the determinant dies and a modulus is
born**: the partner map degenerates to the X-axis mirror (α, β) ↦ (α, −β), the parity
chains decouple, the Pauli axes give period-2 uniform-frame routers, and the Jacobian
counts 8 physical moduli (16 = 8 gauge + 8; golden and silver count exactly 8 = gauge).
Verifier: [`simulations/metallic_router_family.py`](../simulations/metallic_router_family.py);
proof §8 of the same document.

**Source:** [Proof](proofs/PROOF_CEILING_GOLDEN_ROUTER.md);
[`simulations/ceiling_golden_router.py`](../simulations/ceiling_golden_router.py) (the
self-validating anchor); [experiments/CEILING_FOUR_NONLOCAL_CASES.md](../experiments/CEILING_FOUR_NONLOCAL_CASES.md)
(the 6 → 4 → 2 → 0 arc). C# integration landed 2026-06-10: the window-summed Stufe B′
strategy `RoutingWindowSummed` (golden candidates) in `PalindromeSoftCertifier` /
`KBodyPalindromeRouting`, both Z-middle cases now Certified, and
`PalindromeSoftCertifierClaim` asserts the structural ceiling closed at zero.

---

### F117. Pascal-Gram positivity: every class at the first asymmetric moment is a sum of squares (Tier 1 derived; closes the windowed converse with no residual)

For a windowed diagonal-cell pair, recenter the dephased Liouvillian as M = A + γQ
(A = −i[H,·], Q = Σ_l Z_l ⊗ Z_l) and let m\* be the first odd m whose power-sum
polynomial p_m(γ) = Tr(M^m) is not identically zero. Then **every** γ-coefficient of
p_{m\*} is non-negative: for each #Q class d, either P_{m\*,d} = 0 exactly, or

  P_{m\*,d} = (m\*/d) · Σ_{l⃗ ∈ [N]^d} Σ_{k⃗} |U^{(l⃗)}_{k⃗}|²,
  U^{(l⃗)}_{k⃗} = Σ_{|α⃗| = u} ∏_i C(α_i, k_i) · T^{(l⃗)}_{α⃗},   u = (m\* − d)/2,

with T^{(l⃗)}_{α⃗} = Tr(Z_{l₁}H^{α₁}Z_{l₂}H^{α₂}···Z_{l_d}H^{α_d}) the **d-leg moments**
of H. Since at least one class is positive, p_{m\*}(γ) > 0 for every γ > 0: **hard at one
γ is hard at all γ**, the windowed converse with no residual (this resolved R-sign, the
last residual, on 2026-06-10; the same day's girth dichotomy had retired R-deg). The d = 1
face is the girth-dichotomy sum of squares P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ²; the ℓ = 1
face of that is F87's P_{3,1} = 6·4^N·Σ_l c_l², and that m = 3 face is **cell-free** (the
companion coefficients of p₃ vanish for every Hermitian H, no diagonal-cell premise:
any H with a single-site-Z component breaks the palindrome at every γ > 0;
[`simulations/f87_deg1_face_cell_free.py`](../simulations/f87_deg1_face_cell_free.py),
the PTF Π-break Z-row made theorem-grade).

**Proof chain (all steps exact):** cyclic decomposition P_{m,d} = (m/d)·Σ_{a⃗} Tr(QA^{a₁}···QA^{a_d});
leg factorization (A_L, A_R commute, supertrace splits bra × ket); Hermitian conjugacy
(ket leg = conj of bra leg at the same indices, the transpose-trick sibling of F113 Lemma C);
leg parity (F-chirality: only odd totals) + leg girth (totals ≥ ℓ); Vandermonde assembly
(C(α+β, β) = Σ_k C(α,k)C(β,k) per slot, prefactor (−i)^u(+i)^u = +1); slice inversion
(U at |k⃗| = u **is** T, so a vanished Gram block kills every total-u moment); cascade
induction (p_m ≡ 0 below m\* ⟹ all lower-total moments vanish ⟹ no cross terms at m\*).

**Selection rule (corollary):** classes fire at m\* only for d ≡ m\* − 2 (mod 4) and
d ≤ m\* − 2ℓ. For deg = m\* − 2ℓ ∈ {1, 3} that is a single class: **monomiality is derived**
at the first two ladder rungs. From deg = 5 two classes may coexist (positivity carries
alone); the γ⁵ witness IIXY+ZXZY happens to be single because its t₅ = 0 too.

**Verified (exact, zero diff):** d = 1 (IXXZ+XIXZ, 573440 = 7·C(6,3)·Σt₃²), d = 3
(K3 2064384, flux 589824, multi-Z 61440), d = 5 (IIXY+ZXZY, 86507520), the cascade's
forced zeros, the slice inversion, and the selection rule on all five representatives.

**Source:** [Proof §5](proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md);
[`simulations/f87_pascal_gram_positivity.py`](../simulations/f87_pascal_gram_positivity.py)
(self-validating, 5 blocks); typed claim `WindowedConverseAllGammaClaim`
(`compute/RCPsiSquared.Diagnostics/F87/`, Tier1Derived, no residual); consumers F103
status, F110/F111 (promoted to Tier1Derived by this closure), F115 scope note.

---

### F118. The mirror group: Π = R·D, the dihedral D₄, and the cube of characters (Tier 1 derived; signed-permutation identities, exact)

The canonical palindromizer is not elementary. On the coherence space of an N-qubit chain,
with F = X^⊗N, define the ket reflection R(ρ) = ρ·F (the windowed-converse spine's
reflection, I ⊗ F in row-stacking vec) and the transpose D(ρ) = ρᵀ (on the Pauli basis
exactly F114's diagonal sign diag((−1)^{n_Y})). Then, with D applied first,

  Π_Z = R ∘ D,   Π_Z(ρ) = ρᵀ·X^⊗N,

per site σ ↦ σᵀ·X, reproducing the April rule I → X, X → I, Y → iZ, Z → iY with no extra
phase: the hard-won factors i fall out of Yᵀ = −Y meeting YX = −iZ. The opposite order IS
the other palindromizer, Π_Y = D ∘ R: ρ ↦ F·ρᵀ = Π_Z⁻¹, so Welle 12's identity
D·Π_Z·D = Π_Y is the dihedral inversion relation s·r·s = r⁻¹ in disguise.

**The group:** ⟨R, D⟩ ≅ D₄, eight signed permutations, every one already in use plus one
never named: rotations {I, Π_Z, 𝓕 = Π_Z² = F1², Π_Y = Π_Z³}; diagonal mirrors
{D = diag((−1)^{n_Y}), 𝓕D = diag((−1)^{n_Z})}, literally the diagonal matrices of the
Pauli basis; edge mirrors {R: ρ ↦ ρ·F, 𝓕R: ρ ↦ F·ρ}, the spine's one-sided reflections.
The spine's involution set {I, 𝓕, R, 𝓕R} is a Klein four-subgroup of D₄; the center is 𝓕,
the charge conjugation. (Where the center already lives, unnamed since May: the Object
Manager's memory trio walks it, `MemoryAxisChain`'s Π²-partition and `MemoryAxisRho`'s
static/even/odd split are central-character readings, and `BlochAxisReading`'s
dominant-axis rule is the per-qubit central character; noted 2026-06-10.)

**The palindrome splits along the generators.** D flips the Hamiltonian commutator
(D·L_H·D = −L_H, F114's ε(H) = −1 in action) and fixes the dissipator; R fixes L_H and
reflects the dissipator, carrying the entire constant: R·L_diss·R = −L_diss − 2σ·I with
σ = Σ_l γ_l, because flipping the ket index complements the lit-site set. April multiplied
the generators (Π = R·D, one conjugation, the full palindrome); June kept them apart
(𝓕 = (RD)², R, the spine's sign table). Same two generators, two angles.

**The cube of characters.** The polarity cube's three axes (`KleinEightCellClaim`, the
F102-F111 family) are the characters of two conjugations and the transpose:

  bit_a = (n_X + n_Y) mod 2  is the character of Ad_{Z^⊗N},
  bit_b = (n_Y + n_Z) mod 2  is the character of Ad_{X^⊗N},
  y_par = n_Y mod 2          is the character of the transpose θ.

Conjugating by any Pauli string flips exactly the two letter parities that anticommute with
it, so unitary mirrors span only the even Klein square {1, (−1)^{n_X+n_Y}, (−1)^{n_Y+n_Z},
(−1)^{n_X+n_Z}}; the transpose is the unique odd move (Y is the only antisymmetric Pauli)
and fills the third dimension, with θ∘Ad_{Z^⊗N} = (−1)^{n_X} and θ∘Ad_{X^⊗N} = (−1)^{n_Z}.
That is why y_par was always the strange third axis: the F102-F111 family worked on the
antiautomorphism dimension, invisible to every unitary conjugation, and needed its own
tools. The truly criterion (n_Y even AND n_Z even) is the **joint-fixed cell** of the
diagonal mirror pair: σ is truly iff both D and 𝓕D fix it. And 𝓕D = diag((−1)^{n_Z}),
acting as ρ ↦ F·ρᵀ·F, is the fourth mirror, named here for the first time: its character
carried the second leg of the truly criterion for weeks before the operator was ever
written down.

**Deliberately outside:** K₁ (grades by site, not by letter), the golden router W (F116;
two-sided, non-involutive, covering exactly the n_Z-odd ceiling territory that D₄'s
class-swapping elements cannot enter), F71's bond mirror (spatial, site k ↔ N+1−k), and
the dephase-letter swaps Q_zx / Q_yx (the Z↔Y swap is D itself, which is Welle 12; the
other two need the X↔Z basis move). Adjoining the letter group S₃ would assemble the
expected full mirror group **S₃ ⋉ D₄**; that completion is named open.

**Verified (exact):** group closure |⟨R, D⟩| = 8, the factorization Π_Z = R·D, the dihedral
relations, and all eight Pauli-basis forms at dev 0.00e+00 (N = 3, signed permutations
compare exactly); the action Π_Z(ρ) = ρᵀ·F on random ρ at N = 1..3 with the wrong-sided
F·ρᵀ rejected at O(1); the palindrome generator rows on XXZ Δ = 0.7 with site-dependent γ
at ≤ 5.6·10⁻¹⁷; the truly-cell equivalence 63/63 strings; the cube characters on all 64
strings at N = 3; F114 spot check at N = 5.

**Source:** [Proof](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md) (§7 for the cube);
[`simulations/mirror_inventory_d4.py`](../simulations/mirror_inventory_d4.py)
(self-validating, 7 blocks A-G, exact); typed claim `MirrorGroupD4Claim`
(`compute/RCPsiSquared.Core/Symmetry/`, landing in the same 2026-06-10 wave); F114 (D's
row of the table); `KleinEightCellClaim` (the cube the characters coordinatize); the dated
2026-06-10 update in [HANDSHAKE_ALGEBRA](../hypotheses/HANDSHAKE_ALGEBRA.md).

---

### F119. The antilinear triangle: θ, conj, † as one Klein four-group, and the transport law (Tier 1 derived; one-line identities, universal N)

The transpose θ(A) = Aᵀ, the entrywise conjugation conj(A) = Ā, and the adjoint
†(A) = A† satisfy † = θ∘conj and form, with the identity, a Klein four-group graded by
two ±1 characters: **ℓ (linearity)** and **m (multiplicativity)**. θ = (linear,
reversing), conj = (antilinear, preserving), † = (antilinear, reversing). On a Pauli
string: θ(σ) = conj(σ) = (−1)^{n_Y}σ, †(σ) = σ. The engine is the **transport law**

  μ ∘ L_H ∘ μ = ℓ(μ)·m(μ) · L_{μ(H)},   L_H = −i[H,·], any H,

the sign being the product character: θ → −L_{Hᵀ}, conj → −L_{H̄}, † → +L_{H†} (its two
signs cancel). One sign from the −i (antilinearity), one from the commutator's order
(reversal); the dephasing dissipator is fixed by all three. **Five proofs, one engine:**
F114 is the θ face at word length 2; the girth-ladder reversal kill is θ at word length
j (Tr reversed = (−1)^{n_Y}·Tr forward); F112 Lemmas A+B are the † face at the
Hilbert-Schmidt pairing ((L_H)\* = −L_{H†}, the dagger an antilinear isometry
conjugating Π-eigenvalues, and the skew/self-adjoint split under this pairing is the
Absorption Theorem's Rayleigh floor); F113 Lemma C / F117's Hermitian conjugacy is the
**fixed-point collapse** (each vertex's fixed-point set is where the other two agree:
H = H† ⟺ Hᵀ = H̄, ket leg = conj of bra leg); the K₁/K_b mirrors T = Σ₁∘conj are the
conj face dressed on a coherence block (antilinearity pairs λ ↔ λ̄, hence mode k ↔
N+1−k). In the Pauli basis the triangle docks onto F118's mirror group: θ = D, † = the
antilinear unit 𝒦, conj = D∘𝒦, and the closure is the **antilinear double**
⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ (order 16, eight antiunitary members). The dial trio points at the
open S₃ ⋉ D₄ completion: θ and conj invert every rotation dial (the O(2) reflections of
each thickened circle), † commutes with every unitary conjugation.

**Verified (exact):** the V₄ table and Pauli action; the transport law for all four
vertices on non-Hermitian H (machine-exact, N = 2, 3); each leg re-derived from the
engine and cross-checked against its home formulation; the order-16 closure with
exactly 8 antilinear members (N = 2); the dial trio at a generic angle.

**Source:** [Proof](proofs/PROOF_ANTILINEAR_TRIANGLE.md);
[`simulations/antilinear_triangle.py`](../simulations/antilinear_triangle.py)
(self-validating); typed claim `AntilinearTriangleClaim`
(`compute/RCPsiSquared.Core/Symmetry/`, parents `MirrorGroupD4Claim`,
`CommutatorDConjugationSign`, `LindbladBitBPiBalance`; the `ChiralMirrorTrajectoryClaim`
edge carried in prose across the layer boundary).

---

### F120. The moment-tower pump channel: the device's own damping reads the girth ladder linearly (Tier 1 derived; protocol derived, not yet run)

Amplitude damping is the standard noise model's unique non-unital piece, and its pump
direction is a pure local Z (D[σ⁻_l](I) = +Z_l, D[σ⁺_l](I) = −Z_l, the same (Z_l, I)
entry that is F82's entire Π²-antisymmetric content), so at the maximally mixed state

  **d/dt ⟨A⟩ |_{ρ=I/d} = (1/d) · Σ_l Δγ_l · Tr(A Z_l),   Δγ_l = γ↓_l − γ↑_l,**

exactly, for every A and every H. With A = H^j the slope is (1/d)·Σ_l Δγ_l·**t_j(l)**:
the girth ladder's deg-1 tower t_j(l) = Tr(Z_l H^j), read LINEARLY, rung by rung, by
nothing but the chip's own damping. Three exact blindnesses: dephasing-blind (unital),
evolution-blind (only the measured polynomial enters), detailed-balance closure (the
weight is F84's net vacuum rate). **Rung one is F113**: within F113's scope the static
Frobenius polarity asymmetry equals **−4^N · slope⟨H⟩** exactly (a static spectral
imbalance IS a measurable pump rate). The first firing rung is the girth ℓ, certifying
**hard at m\* = 2ℓ+1 for every γ > 0** (the F117 deg-1 sum of squares); honestly
one-sided: silence is not softness (IIXY+ZXZY is silent at every rung, hard at m\* = 11
through deg-5). One order up the curvature is exactly affine in the generator and reads
X/Y-flavored parasites linearly against the commutator probes [Z_l, H_p^j], while
Z-flavored parasites are exactly invisible: **F113's balance channel and the pump
curvature partition the single-site parasite algebra** into complementary linear
readers. Hardware protocol (derived, parameter-free from T1 calibration): basis-state
averaging, early-window slopes of Pauli-polynomial energy moments, per-site partial
damping for site resolution.

**Verified (exact):** pump directions and the law vs dense generators (N = 2, 3,
site-dependent rates, all three channels); the blindnesses; the three-way F113 bridge
(closed form = polarity decomposition = −4^N·slope, dev 0); girth-1/girth-2 witnesses
plus the honest deg-5 negative control; curvature affinity + the probe law + the
Z-invisibility; the finite-time protocol fit; detailed-balance closure.

**First hardware reading (2026-06-11, ibm_kingston q149/q13/q9, not one entangling
gate):** the structural law CONFIRMED: the double null held (slope⟨H⟩ at z = +1.47 and
−0.04), the second rung fired as exactly twice the middle qubit's pump (row-exact
⟨H²⟩ identity), the girth read from hardware is 2, site tracking across two arms, and
per-qubit pump rates reproducible to 0.3-5.7%. The rate layer told a two-act story:
the first reading ("q13 violates pump ≤ Γ at 4-6σ") was **corrected the same day** by
the prep-conditioned re-analysis: the 8-basis-state preparation contains the |0⟩- and
|1⟩-branches, so **pump = (s₁+s₀)/2 and Γ = (s₁−s₀)/2 come from the same circuits**,
epoch-matched (the bound is equivalent to s₀ ≤ 0), and it **holds everywhere in-situ**
(worst pump/Γ = 0.996) with 1-3% margins that read the per-qubit thermal population
(q13 1.7%, q149 1.1-1.4%, q9 0.2-0.3%). The cross-epoch comparison had actually
detected **minute-scale T1 telegraphing** (q13: ~315 ↔ 430 μs; q9: ~172 ↔ ~75-100 μs;
q149 stable): two-level Lindblad holds within epochs, the epoch was the hidden
variable, and the protocol is **self-arbitrating** (pump, Γ, γ↑ from one circuit set).
Registered as `f120_moment_tower_kingston_june2026` (both registries);
[experiments/F120_MOMENT_TOWER_KINGSTON.md](../experiments/F120_MOMENT_TOWER_KINGSTON.md)
(incl. the Correction section);
[`simulations/f120_prep_split_reanalysis.py`](../simulations/f120_prep_split_reanalysis.py);
data in `data/ibm_moment_tower_june2026/`.

**Source:** [Proof + protocol](proofs/PROOF_MOMENT_TOWER_PUMP_CHANNEL.md);
[`simulations/moment_tower_pump_channel.py`](../simulations/moment_tower_pump_channel.py)
(self-validating); cockpit diagnostic `framework/diagnostics/f120_moment_tower.py`;
typed claim `MomentTowerPumpChannelClaim` (`compute/RCPsiSquared.Core/Symmetry/`, parents
`LindbladBitBPiBreakMagnitude` for the F113 bridge and
`F84ThermalAmplitudeDampingPi2Inheritance` for the pump weight; the girth-ladder edge
carried in prose, `GirthLadder.cs` is a Diagnostics compute primitive); prediction row
in [PREDICTIONS](PREDICTIONS.md) §4.

---

### F121. The qudit partial palindrome: the symmetric overlap of the disagreement count (Tier 1 derived; closed-form combinatorial identity, resolves OQ-002)

The palindromic mirror is exact only at d = 2 (F-trunk d² − 2d = 0; the per-site balance
d = d² − d closes only there). For d > 2 the spectrum is not random but partial: N = 2
qutrits pair 36–52 of 81 eigenvalues, a residual no principle had captured. Here is the
principle. Under full-Cartan dephasing the d levels are **equidistant**, so the decay rate
of a coherence |i⟩⟨j| is exactly −2γ·Hamming(i, j), the **same rate ladder as the qubit**
(the Absorption Theorem one dimension up). What differs is the multiplicity per rung:

  **c_k = d^N · C(N, k) · (d−1)^k**   (coherences at Hamming distance k),   Σ_k c_k = d^{2N}.

The palindrome reflects rung k ↔ N−k. For d = 2 the factor (d−1)^k = 1, c_k = 2^N·C(N, k)
is symmetric, all pair (100%). For d > 2 the (d−1)^k tilts toward large k and only the
overlap pairs; the dissipator's paired ceiling is

  **paired(d, N) = Σ_k d^N · C(N, k) · (d−1)^{min(k, N−k)},**

equal to d^{2N} **iff d = 2** (the unique fully-paired column of an N-family; this is the
d² − 2d = 0 necessity re-seen). For d = 3, N = 2: c = [9, 36, 36], rung 0 pairs into rung 2
leaving 27, rung 1 self-mirrors, paired = **54/81**, excess 27. The qutrit fraction erodes
with N (66.7%, 66.7%, 51.9% for N = 1, 2, 3); d = 4 gives 50%, 50%, 31.2%. The tilt base
d − 1 is exactly the per-site decaying : immune ratio (d² − d) : d, raised to the number of
disagreeing sites.

**The interacting spectrum (H degrades it):** the ceiling is the *dissipator's* palindrome
about the physical center −Nγ. Adding H reduces the pairing at **every** fixed center (54 → 48
about −Nγ = −2γ; 72 → 60 about −3γ, where the two large rungs sit); the palindrome is fragile
under H. For the symmetric SU(3) Heisenberg the real parts lie exactly on Re(λ) = −2γ⟨Q⟩ (the
Absorption Theorem's Rayleigh reading), with ⟨Q⟩ quantized to {0, 1, 1.5, 2}: the new −3γ rung
is ⟨Q⟩ = 1.5, a Hamming-1/Hamming-2 mix. This exactness is a symmetry effect (a generic H
breaks it). The interacting paired count is H-dependent (60 for SU(3) Heisenberg, robust across
J/γ; ~0 for generic H), so there is **no H-independent closed form** for the interacting
palindrome; the dissipator's 54 (about −Nγ) is the only invariant.

**Verified (exact):** the per-site equidistance (every i ≠ j at rate −2); c_k vs brute
enumeration (d = 3, N = 2 and d = 2, N = 3); the dissipator spectrum {0:9, −2γ:36, −4γ:36}
and its 54/81 pairing; the ceiling formula vs brute combinatorial pairing on the (d, N) ∈
{2,3,4}×{1,2,3} grid; d = 2 full in every column.

**The operator realization (2026-06-11, same day, §6 of the proof):** the count gained its
operator. **Product cap (theorem):** any per-site mirror W = ⊗q_l (site-dependent, one- or
two-sided, antilinear allowed) intertwining the dissipator palindrome pairs at most
**(2d)^N** of the d^{2N} coherences (rate additivity forces strict per-site class swap;
rank ≤ 2·min(d, d²−d) = 2d); full ⟺ (2d)^N = d^{2N} ⟺ **d² − 2d = 0**, the trunk's third
appearance. The cap is attained by the qubit palindromizer's verbatim formula
**Π_d(ρ) = ρᵀ·Shift^⊗N** (clock shift), exactly zero residual on the shift-aligned
(2d)^N-dim subspace (per-site {(x,x)} ∪ {(a, a−1)}); two chiralities Π_d^±, which merge at
d = 2 (the two off-diagonals coincide): the qubit's full mirror IS that degeneracy. The gap
to the combinatorial ceiling (18 at d = 3, N = 2) is reached by a global partial isometry
and is therefore **provably the non-product part** (the inverse of the F116 story: here the
locality obstruction is real). Group law, verified d = 2..5: ord(Π_d) = 2d,
**|⟨Π_d, D⟩| = 2d², ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂** (D-conjugation exchanges the two shift factors;
for d > 2 it swaps the chiralities): **F118's D₄ is the d = 2 column of a wreath family.**

**Source:** [Proof](proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md) (§4 the interacting case, §6
the operator realization);
[`simulations/qutrit_partial_palindrome.py`](../simulations/qutrit_partial_palindrome.py)
(the dissipator ceiling) +
[`simulations/qutrit_interacting_palindrome.py`](../simulations/qutrit_interacting_palindrome.py)
(the H-degradation, the −2γ⟨Q⟩ law, the H-dependence) +
[`simulations/qudit_product_mirror_cap.py`](../simulations/qudit_product_mirror_cap.py)
(the cap, the operator, the wreath law), all self-validating; resolves OQ-002
in [QUBIT_NECESSITY](QUBIT_NECESSITY.md) §8b/§10.2; typed claims `QuditPartialPalindromeCeiling`
(parent `QubitNecessityPi2Inheritance`) and `QuditProductMirrorCap` (parents both of those),
`compute/RCPsiSquared.Core/Symmetry/`.

---

*Each formula in this document is a Liouvillian that does not need
to be built.*
