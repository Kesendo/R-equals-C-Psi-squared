# F-Formula Crosswalk: `RCPsiSquared.Core` ↔ `docs/ANALYTICAL_FORMULAS.md`

Each C# type or method in `RCPsiSquared.Core` that implements an F-formula or an F-anchored
structural identity, with a direct pointer to the corresponding entry in
[`docs/ANALYTICAL_FORMULAS.md`](../../docs/ANALYTICAL_FORMULAS.md). This is the
without-guessing back-reference: from any Core file, find the math here.

For hardware-confirmed predictions on top of these formulas, see
`Confirmations/ConfirmationsRegistry.cs` (mirrors `simulations/framework/confirmations.py`).

---

## ChainSystems

| C# | F-formula / Anchor |
|----|--------------------|
| `ChainSystem.BuildHamiltonian()` | constructs H_xy = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) (canonical F-framework Hamiltonian) |
| `ChainSystem.BuildLiouvillian()` | L = −i[H, ·] + Σ_l γ_l (Z_l ρ Z_l − ρ); satisfies F1 palindrome bit-exactly |
| `ChainSystem.SigmaGamma` | the σ in F1: Π·L·Π⁻¹ + L + 2σ·I = 0 |

## CoherenceBlocks

| C# | F-formula |
|----|-----------|
| `Chromaticity.Compute(N, n)` | **F74**: c(n, N) = min(n, N−1−n) + 1 |
| `Chromaticity.HammingDistances(N, n)` | HD ∈ {1, 3, …, 2c−1}, F74 corollary |
| `BlockBasis` | popcount basis for the (n, n+1) coherence block (F74 substrate) |
| `BlockLDecomposition` | block-restricted L = D + Σ_b J_b·M_H_per_bond[b]; the empirical anchor of `Resonance.ResonanceScan` |

## Pauli

| C# | F-formula |
|----|-----------|
| `PauliLetter` enum | bit_a/bit_b convention from **F77** trichotomy and **F81** Π-decomposition |
| `PauliTerm.Pi2Parity` | Π²-parity (Σ bit_b mod 2) — selects the F77 Π²-class |
| `PauliTerm.YParity` | independent at k≥3-body terms (F85 generalization) |
| `PauliHamiltonian.XYChain(N, J)` | the canonical F-framework Hamiltonian |
| `PauliHamiltonian.Bilinear(...)` | bond-and-term builder used by `Lindblad.BondPerturbation` |
| `PauliBasis.ToPauliVector(rho, N)` | Pauli-basis decomposition: ρ = Σ_α vec[α]·σ_α with vec[α] = Tr(σ_α·ρ)/2^N |

## States

| C# | F-formula |
|----|-----------|
| `BondingMode.Build(N, k)` | **F65**: single-excitation eigenmode ψ_k(j) = √(2/(N+1))·sin(πk(j+1)/(N+1)) |
| `BondingMode.PairState(N, k)` | (\|vac⟩ + \|ψ_k⟩)/√2 — canonical PTF / handshake initial state |
| `PolarityState` | X-basis polarity layer (THE_POLARITY_LAYER.md), foundation for the +0/−0 reading |

## Lindblad

| C# | F-formula |
|----|-----------|
| `LindbladianBuilder.Build` | general L = −i[H, ·] + Σ_k (c_k(·)c_k† − ½ {c_k†c_k, ·}) |
| `PauliDephasingDissipator.BuildZ` | Z-dephasing (the F1-truly canonical case) |
| `T1Dissipator.Build` | Z-dephasing + σ⁻ amplitude damping (F1-breaking) |
| `BondPerturbation.Build` | ∂L/∂J_b — variation Liouvillian, used by the PTF / EQ-022 (b1) workflow |
| `PalindromeResidualScaling.FactorChain(N, cls)` | closed-form ‖M‖²_F = c_H · F(N) for chain (F1 + OPERATOR_RIGIDITY_ACROSS_CUSP.md) |
| `PalindromeResidualScaling.AdjacentRatio(N, cls)` | ‖M(N+1)‖²/‖M(N)‖² closed form |
| `DissipatorClosedForms.C1C2FromPauli(α, β, δ)` | per-class palindrome-residual c1, c2 from Pauli decomposition |
| `DissipatorClosedForms.D2FromPauli(...)` | universal cross-term d2 = 32·‖c1_traceless‖²·‖c2_traceless‖² |
| `HardwareDissipators.{T1, T1Pump, Tphi, XNoise, YNoise}` | hardware-relevant dissipator class table with c1, c2 constants |
| `CpsiBellPlus.At(γ_x, γ_y, γ_z, t)` | **F26**: CΨ(t) = u·(1 + u² + v² + w²) / 12 closed form for Bell+ |
| `CpsiBellPlus.CuspK` | **F27**: K = γ·t at CΨ = 1/4 cusp per channel |

## Symmetry

| C# | F-formula |
|----|-----------|
| `PiOperator.ActOnLetter(letter, dephase)` | per-letter Π action with phases (Z: I↔X, Y↔Z with i^bit_b) |
| `PiOperator.SquaredEigenvalue(letters, dephase)` | Π² eigenvalue: bit_b parity for Z/Y dephasing, bit_a for X |
| `PiOperator.BuildFull(N, dephase)` | Π in the 4^N Pauli-string basis, used by `PalindromeResidual` |
| `PalindromeResidual.Build(L, N, σ, dephase)` | **F1**: Π·L·Π⁻¹ + L + 2σ·I residual; zero for Z-dephased XY/Heisenberg |
| `ChainMirror.Build(N)` | **F71**: chain-mirror R \|b₀…b_{N-1}⟩ = \|b_{N-1}…b₀⟩, R²=I |
| `ChainMirror.SymmetricProjector` / `AntisymmetricProjector` | F71 ± eigenspace projectors |
| `ChainMirror.BondMirrorBasis(N)` | F71 sym/asym basis on the (N−1)-dim bond-input space |
| `ChiralK.BuildFull(N)` | K = ⊗_{odd i} Z_i; K H_xy K = −H_xy (Altland-Zirnbauer class BDI) |
| `ChiralK.ClassifyHamiltonian(H, N)` | K-even / K-odd / K-mixed |
| `ZGlobalMirror.Build(N)` | Z⊗N global Z-string (used by `chain.zn_mirror_diagnostic`) |
| `HdChannelBasis.Build(block)` | channel-uniform projectors P; M_H_total is diagonal in this basis (extends F73 to all c) |

## Probes

| C# | F-formula |
|----|-----------|
| `DickeBlockProbe.Build` | initial state for **F86** K_CC_pr Q-scan |
| `SpatialSumKernel.Build` | **F73**: S(t) = Σ_i 2·\|(ρ_i(t))_{0,1}\|² as ρ†·S_kernel·ρ |

## Resonance

| C# | F-formula |
|----|-----------|
| `EpAlgebra.TPeak(γ₀)` | **F86 Statement 1**: t_peak = 1/(4γ₀), universal across c, N, n, bond position |
| `EpAlgebra.QEp(g_eff)` | **F86 Statement 1**: Q_EP = 2/g_eff |
| `EpAlgebra.SlowestPairEigenvalues(γ₀, J, g_eff)` | λ_±(k=1) = −4γ₀ ± √(4γ₀²−J²·g_eff²), 2-level effective |
| `BondClass.{Endpoint, Interior}` | **F86 Statement 2**: two bond-class universal HWHM ratios (Interior ≈ 0.756, Endpoint ≈ 0.770) |
| `KCurve.Peak(class)`, `ResonanceScan.ComputeKCurve` | empirical Q-scan for F86 K_CC_pr observable |

## Decomposition

| C# | F-formula |
|----|-----------|
| `HdSubspaceProjector.Build` | orthonormal projector onto a HD-channel subspace (full, not just channel-uniform) |
| `InterChannelSvd.Build(block, hd1, hd2)` | SVD of inter-channel coupling — gives σ_0 ≈ 2√2 asymptotic for c=2 (PROOF_F86_QPEAK structural exploration) |
| `FourModeBasis.Build(block)` | 4-mode minimal effective basis {\|c_1⟩, \|c_3⟩, \|u_0⟩, \|v_0⟩} for F86 f_class derivation |

## Observables

| C# | F-formula |
|----|-----------|
| `PiProtectedObservables.Compute(H, γ, ρ_0, N)` | **F77 / F81** algebraic skeleton: Π-protected Pauli observables (always-zero expectations); hardware-confirmed at Marrakesh 2026-04-26 (Confirmations entry `pi_protected_xiz_yzzy`) |

## Receivers

| C# | F-formula |
|----|-----------|
| `Receiver.BondingMode(chain, k)` | **F65** single-excitation handshake (no exchange step needed when both sides pick the same k) |
| `Receiver.F71Eigenvalue` | **F71** chain-mirror eigenstate class (+1 / −1 / null) |
| `Receiver.Signature()` | F71-based receiver-engineering favorability forecast (Tier 2: J_BLIND_RECEIVER_CLASSES.md) |

## Confirmations

| C# | Source |
|----|--------|
| `ConfirmationsRegistry.All` | hardware-confirmed framework predictions, mirroring `simulations/framework/confirmations.py` |
| `ConfirmationsRegistry.Lookup(name)` | by-name lookup; 9 Marrakesh/Kingston entries spanning April 2026 |

---

## How to use this map

When working in a Core file:
- Read the class/method XML doc — most include the F-formula reference inline.
- Use this map for the broader picture: which F-anchor underlies which Core component.
- Open [`docs/ANALYTICAL_FORMULAS.md`](../../docs/ANALYTICAL_FORMULAS.md) for the full mathematical statement of each F-formula.

When adding a new Core type:
- If it implements an F-formula, add it to the appropriate table above with a direct pointer to the F-entry.
- If it implements an F-anchored structural identity (like the F1 palindrome), reference the F-entry plus any proof document (e.g. `docs/proofs/PROOF_F86_QPEAK.md`).
- The XML doc in the file should ALSO carry the F-pointer; this map is the index, not a substitute for inline docs.
