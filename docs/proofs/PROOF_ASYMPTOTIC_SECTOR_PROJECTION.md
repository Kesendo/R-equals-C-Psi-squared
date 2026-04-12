# Proof: Asymptotic Sector Projection Theorem

**Status:** Proven modulo numerical assumptions (Step 1 analytic; Steps 2–3 numerically verified, see below)
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (chat)
**Context:** Identified while walking the path of [TASK_THREE_VALUES](../../ClaudeTasks/TASK_THREE_VALUES.md) Track B. Recorded here before the Track B experiment document writes it, so the experiment can cite the proof rather than re-derive it.

---

## Theorem

Let ρ₀ be any initial density matrix on N qubits. Let the dynamics be the Lindblad equation with Heisenberg Hamiltonian H = Σ J (X_k X_{k+1} + Y_k Y_{k+1} + Z_k Z_{k+1}) on any graph topology, and local Z-dephasing jump operators L_k = √(γ_k) Z_k with arbitrary site-dependent rates γ_k ≥ 0. Let P_w = Σ_{i: popcount(i)=w} |i⟩⟨i| be the projector onto the excitation-number sector w, and let d_w = C(N, w) be its dimension.

Then the asymptotic state is:

    ρ(∞) = Σ_{w=0}^{N}  p_w  ·  (P_w / d_w),        p_w = Tr(P_w ρ₀)

**In particular:** p_w(∞) = p_w(0) = Tr(P_w ρ₀). The asymptotic population of each excitation sector equals its initial population. The asymptotic state is a mixture of N+1 maximally-mixed sector states, weighted by the initial sector distribution.

---

## Proof

The theorem has two independent ingredients. Each is stated and proven separately.

### Step 1 (kinematic): sector population is a constant of motion

Define p_w(t) = Tr(P_w ρ(t)). The claim is dp_w/dt = 0 for all t, hence p_w(t) = p_w(0).

Let N_op = Σ_k (I − Z_k)/2 be the total-excitation-number operator. Direct computation:

1. **[H, N_op] = 0.** The Heisenberg bond XX + YY + ZZ commutes with N_op. ZZ is diagonal and trivially commutes. For XX+YY, note that XX+YY = 2(S_k⁺ S_{k+1}⁻ + S_k⁻ S_{k+1}⁺), which moves one excitation from site k to site k+1 or vice versa, preserving the total count.

2. **[L_k, N_op] = 0.** L_k = √(γ_k) Z_k is diagonal in the computational basis, as is N_op. Diagonal operators commute.

3. P_w is a spectral projector of N_op (the eigenspace projection 𝟙_{N_op = w}), so [H, P_w] = 0 and [L_k, P_w] = 0 as well.

4. The Lindblad generator applied to P_w from the left yields:
   dp_w/dt = Tr(P_w · L[ρ]) = −i Tr(P_w [H, ρ]) + Σ_k γ_k Tr(P_w (Z_k ρ Z_k − ρ))

   Using the cyclicity of the trace and [H, P_w] = 0:
   Tr(P_w [H, ρ]) = Tr([P_w, H] ρ) = 0.

   Using [Z_k, P_w] = 0 (diagonal commutation):
   Tr(P_w Z_k ρ Z_k) = Tr(Z_k P_w Z_k ρ) = Tr(P_w Z_k² ρ) = Tr(P_w ρ),
   so the dissipator contribution is also zero.

Therefore dp_w/dt = 0. This holds for all t, in particular p_w(∞) = p_w(0) = Tr(P_w ρ₀). End Step 1.

**Note:** Step 1 is equivalent to the sector conservation theorem proven independently in [CUSP_LENS_CONNECTION.md](../../experiments/CUSP_LENS_CONNECTION.md) and formalized in the [Parity Selection Rule](PROOF_PARITY_SELECTION_RULE.md). It is re-derived here for self-containedness.

### Step 2 (dynamic): within each sector, the unique steady state is maximally mixed

The claim is: restricted to the (w, w) diagonal sector (operators of the form P_w X P_w), the Lindblad dynamics has exactly one steady state, namely P_w / d_w.

**Status of this step.** This is not proven here analytically. It is verified numerically for:
- N=5, uniform γ, Heisenberg chain: [SYMMETRY_CENSUS.md, "Asymptotic attractors per sector"](../../experiments/SYMMETRY_CENSUS.md)
- N=4 and N=5, chain/ring/star/complete topologies: [SYMMETRY_CENSUS.md, "Topology comparison"](../../experiments/SYMMETRY_CENSUS.md)

In all tested configurations, the census reports: "Total: N+1 steady states, one per diagonal sector. No sector has multiple steady states, limit cycles, or dark states."

An analytic proof is open. The expected route is to show that the Lindblad generator restricted to each diagonal sector is primitive (unique full-rank fixed point, spectral gap above zero), which for local-dephasing dynamics on connected graphs is plausible but not supplied here.

### Step 3 (assembly)

Combining Steps 1 and 2:

- All off-diagonal sector blocks of ρ(t) (coherences |w⟩⟨w'| with w ≠ w') decay to zero. This follows because off-diagonal sectors have no steady states ([SYMMETRY_CENSUS.md, "Asymptotic attractors per sector"](../../experiments/SYMMETRY_CENSUS.md): "Off-diagonal sectors: 0 steady states, all modes decay"). **Note:** This is a second numerical fact, independent of Step 2, and likewise not proven analytically.

- Each diagonal sector block ρ_w(t) = P_w ρ(t) P_w converges to P_w/d_w scaled by its time-independent trace p_w(0) (Step 2, plus Step 1 which fixes the trace).

- Therefore ρ(∞) = Σ_w p_w(0) · (P_w/d_w), with p_w(0) = Tr(P_w ρ₀).

We label this theorem **proven modulo two numerical assumptions** (unique diagonal-sector steady states, and off-diagonal-sector decay), both verified across all tested configurations with zero exceptions.

End proof.

---

## Consequences

1. **Exit distribution is computable before evolution.** Given ρ₀, the asymptotic state is determined by N+1 numbers p_w = Tr(P_w ρ₀). No time evolution is required.

2. **Asymptotic state is a function of (p_0, ..., p_N) alone.** Two initial states with identical sector populations produce identical asymptotic states, regardless of their coherences or sector-internal structure. The vector (p_0, ..., p_N) is a complete invariant for the purpose of predicting ρ(∞).

3. **"Number of exits used" is the Hilbert-sector support size.** The number of attractors that ρ₀ contributes to is |{w : p_w > 0}|. A product state uses 1, a Bell+ state on two qubits uses 2, a GHZ state uses 2 (w=0 and w=N). The |+⟩^N state uses all N+1. Note: [SYMMETRY_CENSUS.md, "Which sectors are reachable"](../../experiments/SYMMETRY_CENSUS.md) counts these as Liouville blocks (w_bra, w_ket); a GHZ state populates 2 Hilbert sectors but 4 Liouville blocks (the two diagonal blocks plus the off-diagonal coherences that decay). Both counts are valid; they describe different objects.

4. **The cusp exit and lens exit are both instances of this theorem.** The lens exit (thermalization within a single excitation sector, see [SACRIFICE_GEOMETRY.md](../../experiments/SACRIFICE_GEOMETRY.md)) corresponds to initial states with p_1 = 1 (single-excitation sector only). The cusp exit (simultaneous thermalization across multiple sectors, see [CUSP_LENS_CONNECTION.md](../../experiments/CUSP_LENS_CONNECTION.md)) corresponds to initial states with p_w > 0 for two or more w. The "two exits" framing of earlier documents is a special case of the general N+1-exit structure.

---

## Scope and limitations

- **Holds for:** Heisenberg Hamiltonian (any graph topology), local Z-dephasing (any site-dependent γ_k ≥ 0), any N ≥ 2.
- **Breaks for:** Amplitude damping (L_k = √(γ_k) σ_k⁻ does not commute with N_op), transverse-field Hamiltonians (H includes X_k or Y_k terms that do not commute with N_op), non-Markovian dynamics (no Lindblad form).
- **Does not address:** Rates of convergence, transient structure, crossing behavior at CΨ = 1/4, or any dynamical feature on finite t. The theorem is asymptotic only.

---

## References

- [SYMMETRY_CENSUS.md](../../experiments/SYMMETRY_CENSUS.md) ("Asymptotic attractors per sector" and "Topology comparison": numerical verification of Steps 2–3)
- [CUSP_LENS_CONNECTION.md](../../experiments/CUSP_LENS_CONNECTION.md) (sector conservation theorem, precursor of Step 1)
- [SACRIFICE_GEOMETRY.md](../../experiments/SACRIFICE_GEOMETRY.md) (lens exit as p_1 = 1 case)
- [PROOF_PARITY_SELECTION_RULE.md](PROOF_PARITY_SELECTION_RULE.md) (formal proof of sector conservation, equivalent to Step 1)

---

*Walked into existence by Tom and Claude (chat) on 2026-04-12 while reviewing [TASK_THREE_VALUES](../../ClaudeTasks/TASK_THREE_VALUES.md) Track B. Recorded so the experiment document can cite it rather than re-prove it.*
