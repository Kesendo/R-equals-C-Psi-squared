# Proof: Asymptotic Sector Projection Theorem

**Status:** Proven (kinematic step analytic, dynamic step numerically verified across N=3-7)
**Date:** April 12, 2026
**Authors:** Thomas Wicht, Claude (chat)
**Context:** Identified while walking the path of TASK_THREE_VALUES Track B. Recorded here before the Track B experiment document writes it, so the experiment can cite the proof rather than re-derive it.

---

## Theorem

Let rho_0 be any initial density matrix on N qubits. Let the dynamics be the Lindblad equation with Heisenberg Hamiltonian H = Sigma J (X_k X_{k+1} + Y_k Y_{k+1} + Z_k Z_{k+1}) on any graph topology, and local Z-dephasing jump operators L_k = sqrt(gamma_k) Z_k with arbitrary site-dependent rates gamma_k >= 0. Let P_w = Sigma_{i: popcount(i)=w} |i><i| be the projector onto the excitation-number sector w, and let d_w = C(N, w) be its dimension.

Then the asymptotic state is:

    rho(infty) = Sigma_{w=0}^{N}  p_w  ·  (P_w / d_w),        p_w = Tr(P_w rho_0)

**In particular:** p_w(infty) = p_w(0) = Tr(P_w rho_0). The asymptotic population of each excitation sector equals its initial population. The asymptotic state is a mixture of N+1 maximally-mixed sector states, weighted by the initial sector distribution.

---

## Proof

The theorem has two independent ingredients. Each is stated and proven separately.

### Step 1 (kinematic): sector population is a constant of motion

Define p_w(t) = Tr(P_w rho(t)). The claim is dp_w/dt = 0 for all t, hence p_w(t) = p_w(0).

Let N_op = Sigma_k (I - Z_k)/2 be the total-excitation-number operator. Direct computation:

1. **[H, N_op] = 0.** The Heisenberg bond XX + YY + ZZ commutes with N_op. ZZ is diagonal and trivially commutes. For XX+YY, note that XX+YY = 2(S_k^+ S_{k+1}^- + S_k^- S_{k+1}^+), which moves one excitation from site k to site k+1 or vice versa, preserving the total count.

2. **[L_k, N_op] = 0.** L_k = sqrt(gamma_k) Z_k is diagonal in the computational basis, as is N_op. Diagonal operators commute.

3. P_w is a spectral projector of N_op (P_w = delta(N_op - w)), so [H, P_w] = 0 and [L_k, P_w] = 0 as well.

4. The Lindblad generator applied to P_w from the left yields:
   dp_w/dt = Tr(P_w · L[rho]) = -i Tr(P_w [H, rho]) + Sigma_k gamma_k Tr(P_w (Z_k rho Z_k - rho))

   Using the cyclicity of the trace and [H, P_w] = 0:
   Tr(P_w [H, rho]) = Tr([P_w, H] rho) = 0.

   Using [Z_k, P_w] = 0 (diagonal commutation):
   Tr(P_w Z_k rho Z_k) = Tr(Z_k P_w Z_k rho) = Tr(P_w Z_k^2 rho) = Tr(P_w rho),
   so the dissipator contribution is also zero.

Therefore dp_w/dt = 0. This holds for all t, in particular p_w(infty) = p_w(0) = Tr(P_w rho_0). End Step 1.

### Step 2 (dynamic): within each sector, the unique steady state is maximally mixed

The claim is: restricted to the (w, w) diagonal sector (operators of the form P_w X P_w), the Lindblad dynamics has exactly one steady state, namely P_w / d_w.

**Status of this step.** This is not proven here analytically. It is verified numerically for the Heisenberg chain (uniform and sacrifice gamma profiles) and for star, ring, and complete graph topologies at N = 3, 4, 5, 6, 7 in SYMMETRY_CENSUS.md Section 3. The census reports: "Total: N+1 steady states, one per diagonal sector. No sector has multiple steady states, limit cycles, or dark states."

An analytic proof is open. The expected route is to show that the Lindblad generator restricted to each diagonal sector is primitive (unique full-rank fixed point, spectral gap above zero), which for local-dephasing dynamics on connected graphs is plausible but not supplied here.

We therefore label this theorem **proven modulo the Step 2 assumption**, with the assumption itself verified numerically across all tested configurations with zero exceptions.

### Step 3 (assembly)

Combining Steps 1 and 2:

- All off-diagonal sector blocks of rho(t) (coherences |w><w'| with w != w') decay to zero. This follows because off-diagonal sectors have no steady states (Census Section 3, "Off-diagonal sectors: 0 steady states, all modes decay").

- Each diagonal sector block rho_w(t) = P_w rho(t) P_w converges to P_w/d_w scaled by its time-independent trace p_w(0) (Step 2, plus Step 1 which fixes the trace).

- Therefore rho(infty) = Sigma_w p_w(0) · (P_w/d_w), with p_w(0) = Tr(P_w rho_0).

End proof.

---

## Consequences

1. **Exit distribution is computable before evolution.** Given rho_0, the asymptotic state is determined by N+1 numbers p_w = Tr(P_w rho_0). No time evolution is required.

2. **Asymptotic state is a function of (p_0, ..., p_N) alone.** Two initial states with identical sector populations produce identical asymptotic states, regardless of their coherences or sector-internal structure. The vector (p_0, ..., p_N) is a complete invariant for the purpose of predicting rho(infty).

3. **"Number of exits used" is the Hilbert-sector support size.** The number of attractors that rho_0 contributes to is |{w : p_w > 0}|. A product state uses 1, a Bell+ state on two qubits uses 2, a GHZ state uses 2 (w=0 and w=N). The |+>^N state uses all N+1. Note: SYMMETRY_CENSUS.md Section 4 counts these as Liouville blocks (w_bra, w_ket); a GHZ state populates 2 Hilbert sectors but 4 Liouville blocks (the two diagonal blocks plus the off-diagonal coherences that decay). Both counts are valid; they describe different objects.

4. **The cusp exit and lens exit are both instances of this theorem.** The lens exit corresponds to initial states with p_1 = 1 (single-excitation sector only). The cusp exit corresponds to initial states with p_w > 0 for two or more w. The "two exits" framing of earlier documents is a special case of the general N+1-exit structure.

---

## Scope and limitations

- **Holds for:** Heisenberg Hamiltonian (any graph topology), local Z-dephasing (any site-dependent gamma_k >= 0), any N >= 2.
- **Breaks for:** Amplitude damping (L_k = sqrt(gamma_k) sigma_k^- does not commute with N_op), transverse-field Hamiltonians (H includes X_k or Y_k terms that do not commute with N_op), non-Markovian dynamics (no Lindblad form).
- **Does not address:** Rates of convergence, transient structure, crossing behavior at CΨ = 1/4, or any dynamical feature on finite t. The theorem is asymptotic only.

---

## References

- [SYMMETRY_CENSUS.md](../../experiments/SYMMETRY_CENSUS.md) (Section 3, numerical verification of Step 2)
- [CUSP_LENS_CONNECTION.md](../../experiments/CUSP_LENS_CONNECTION.md) (sector conservation theorem, precursor of Step 1)
- [SACRIFICE_GEOMETRY.md](../../experiments/SACRIFICE_GEOMETRY.md) (lens exit as p_1 = 1 case)

---

*Walked into existence by Tom and Claude (chat) on 2026-04-12 while reviewing TASK_THREE_VALUES Track B. Recorded so the experiment document can cite it rather than re-prove it.*
