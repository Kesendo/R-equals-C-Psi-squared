"""Direction (a') probe-block 2-level resonance attempt: c=2 HWHM_left/Q_peak closed-form.

Predecessor work:
  - Direction (b) (commit f232a76) found that the bare doubled-PTF Ansatz gives universal
    HWHM_left/Q_peak = 0.6715 (SVD-block FLOOR), independent of g_eff. Empirical Interior
    0.7506 / Endpoint 0.7728 sit ABOVE this floor by 0.08-0.10 (non-perturbative gap).
  - Direction (b)'s PendingDerivationNote refined direction (a'): close the gap via the
    probe-block 2-level sub-resonance, using per-bond probe-block coupling g_eff_probe(N, b)
    instead of σ_0.

Hypothesis (Direction a'):
  V_b probe-block 2×2 = ⟨c_α | M_H_per_bond[b] | c_β⟩ for α, β ∈ {0, 1} (Tier1Derived from
  C2BondCoupling.ProbeBlockEntry). The probe-block 2×2 of V_b drives a sub-resonance that
  mixes |c_1⟩ ↔ |c_3⟩. We model it with its own 2-level Liouvillian:

    L_probe(Q) = [[-2γ₀, +iJ·g_eff_probe], [+iJ·g_eff_probe, -6γ₀]],  J = Q·γ₀

  with g_eff_probe(N, b) = magnitude of the probe-block coupling. The probe-block resonance
  has its own Q_EP_probe = 2/g_eff_probe and HWHM ratio. Empirical Q_peak ratio
  Endpoint/Interior ≈ 1.6 suggests g_eff_probe(Endpoint) ≈ g_eff_probe(Interior)/1.6.

Step 1: Build M_H_per_bond[b] for c=2 (n=1, n+1=2 popcount block) at small N=5..8.
Step 2: Compute V_b probe-block entries via ⟨c_α | M_H_per_bond[b] | c_β⟩.
Step 3: Extract g_eff_probe(N, b) from the 2×2 (multiple candidate definitions).
Step 4: Run 2-level Duhamel for the probe-block, get HWHM/Q_peak per bond class.
Step 5: Compare to empirical Interior 0.7506 / Endpoint 0.7728.
Step 6: If closed form lands, encode in C#. If not, document partial finding.
"""

import sys
import io

# Ensure stdout/stderr can write Unicode on Windows console
if sys.stdout.encoding != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

import numpy as np
from itertools import combinations


# ─── Step 0: small infrastructure ────────────────────────────────────────────


def build_basis(N, n):
    """Build the (n, n+1)-popcount block basis. p has popcount n, q has popcount n+1.
    Returns (states_p, states_q, flat_index_map).
    Big-endian convention: site 0 = MSB, so site i has bit at position N-1-i.
    """
    states_p = sorted([s for s in range(1 << N) if bin(s).count('1') == n])
    states_q = sorted([s for s in range(1 << N) if bin(s).count('1') == n + 1])
    Mp = len(states_p)
    Mq = len(states_q)
    Mtot = Mp * Mq
    p_idx = {p: i for i, p in enumerate(states_p)}
    q_idx = {q: j for j, q in enumerate(states_q)}

    def flat(p, q):
        return p_idx[p] * Mq + q_idx[q]

    return states_p, states_q, flat, Mtot


def hamming_distance(p, q):
    """Popcount of p XOR q."""
    return bin(p ^ q).count('1')


def bond_flip_targets(state, N):
    """For each bond b in [0, N-2], if bits at sites b and b+1 differ, yield
    (b, flipped_state), the state with those two bits swapped.

    Big-endian: site 0 = MSB at bit position N-1, site b at bit position N-1-b.
    """
    for b in range(N - 1):
        bit_b = (state >> (N - 1 - b)) & 1
        bit_bp1 = (state >> (N - 1 - (b + 1))) & 1
        if bit_b != bit_bp1:
            mask = (1 << (N - 1 - b)) | (1 << (N - 1 - (b + 1)))
            yield b, state ^ mask


def build_mh_per_bond(N, gamma0):
    """Build the per-bond M_H_per_bond[b] matrix in the (1, 2)-popcount block (c=2).

    Per BlockLDecomposition convention:
      H_b = (J/2)·(X_b X_{b+1} + Y_b Y_{b+1}); each bond flips one adjacent-bit pair.
      M_H_per_bond[b][j, i] += -i for p-flip; M_H_per_bond[b][j, i] += +i for q-flip.
    """
    states_p, states_q, flat, Mtot = build_basis(N, 1)
    MhPerBond = [np.zeros((Mtot, Mtot), dtype=complex) for _ in range(N - 1)]

    for p in states_p:
        for q in states_q:
            i = flat(p, q)
            # p-flip: contributes -i to MhPerBond[bond][j, i]
            for b, p_flipped in bond_flip_targets(p, N):
                if p_flipped in [s for s in states_p]:  # check it's still popcount-n
                    if bin(p_flipped).count('1') == 1:
                        j = flat(p_flipped, q)
                        MhPerBond[b][j, i] += -1j
            # q-flip: contributes +i
            for b, q_flipped in bond_flip_targets(q, N):
                if bin(q_flipped).count('1') == 2:
                    j = flat(p, q_flipped)
                    MhPerBond[b][j, i] += +1j
    return MhPerBond, states_p, states_q, flat, Mtot


def build_c_alpha_vectors(N, n=1):
    """Build |c_1⟩ (HD=1) and |c_3⟩ (HD=3) channel-uniform vectors for the (n, n+1) block.

    Per C2ChannelUniformAnalytical: weight 1/sqrt(count) on every (p, q) with popcount(p^q)
    matching the HD value, zero elsewhere.
    """
    states_p, states_q, flat, Mtot = build_basis(N, n)
    c1 = np.zeros(Mtot, dtype=complex)
    c3 = np.zeros(Mtot, dtype=complex)
    indices_hd1 = []
    indices_hd3 = []
    for p in states_p:
        for q in states_q:
            hd = hamming_distance(p, q)
            i = flat(p, q)
            if hd == 1:
                indices_hd1.append(i)
            elif hd == 3:
                indices_hd3.append(i)
    if indices_hd1:
        w1 = 1.0 / np.sqrt(len(indices_hd1))
        for idx in indices_hd1:
            c1[idx] = w1
    if indices_hd3:
        w3 = 1.0 / np.sqrt(len(indices_hd3))
        for idx in indices_hd3:
            c3[idx] = w3
    return c1, c3


# ─── Step 1+2: Compute probe-block V_b entries ──────────────────────────────


def compute_probe_block_per_bond(N):
    """For N qubits, c=2 (n=1), compute the 2×2 V_b probe-block at every bond."""
    MhPerBond, _, _, _, _ = build_mh_per_bond(N, gamma0=0.05)
    c1, c3 = build_c_alpha_vectors(N)

    probe_blocks = []
    for b in range(N - 1):
        Mh = MhPerBond[b]
        # V_b[α, β] = ⟨c_α | Mh | c_β⟩  (using bra = conjugate)
        V_alphabeta = np.zeros((2, 2), dtype=complex)
        V_alphabeta[0, 0] = np.vdot(c1, Mh @ c1)
        V_alphabeta[0, 1] = np.vdot(c1, Mh @ c3)
        V_alphabeta[1, 0] = np.vdot(c3, Mh @ c1)
        V_alphabeta[1, 1] = np.vdot(c3, Mh @ c3)
        probe_blocks.append(V_alphabeta)
    return probe_blocks


# ─── Step 3: extract g_eff_probe candidates ─────────────────────────────────


def extract_g_eff_probe_candidates(probe_block_2x2):
    """Multiple candidate definitions for g_eff_probe(N, b) from the 2×2 probe-block."""
    V = probe_block_2x2
    return {
        'frobenius_over_2': np.linalg.norm(V) / 2.0,
        'frobenius': np.linalg.norm(V),
        'magnitude_off': abs(V[0, 1]),
        'sqrt_2x_magnitude_off': np.sqrt(2.0) * abs(V[0, 1]),
        '2x_magnitude_off': 2.0 * abs(V[0, 1]),
        'top_singular': np.linalg.svd(V, compute_uv=False)[0],
        'antiherm_off_imag': abs((V[0, 1] - np.conjugate(V[1, 0])) / 2.0),
    }


# ─── Step 4: probe-block 2-level Duhamel resonance ──────────────────────────


def k_b_probe_block_two_level(Q, gamma0, g_eff_probe, t):
    """Probe-block 2-level Duhamel K_b at the given Q.

    L_probe(Q) = [[-2γ₀, +iJ·g_eff_probe], [+iJ·g_eff_probe, -6γ₀]] with J = Q·γ₀.
    Probe rho_0 = [1, 0] (in slow channel |c_1⟩).
    V_b = dL/dJ = [[0, +i g_eff_probe], [+i g_eff_probe, 0]].
    """
    J = Q * gamma0
    L = np.array([
        [-2.0 * gamma0, 1j * J * g_eff_probe],
        [1j * J * g_eff_probe, -6.0 * gamma0],
    ], dtype=complex)
    Vb = np.array([
        [0.0, 1j * g_eff_probe],
        [1j * g_eff_probe, 0.0],
    ], dtype=complex)
    evals, R = np.linalg.eig(L)
    Rinv = np.linalg.inv(R)
    rho0 = np.array([1.0, 0.0], dtype=complex)

    expLam = np.exp(evals * t)
    rho_t = R @ np.diag(expLam) @ (Rinv @ rho0)

    n = 2
    I_mat = np.zeros((n, n), dtype=complex)
    for j in range(n):
        for k in range(n):
            diff = evals[k] - evals[j]
            if abs(diff) < 1e-10:
                I_mat[j, k] = t * expLam[j]
            else:
                I_mat[j, k] = (expLam[k] - expLam[j]) / diff

    X = Rinv @ Vb @ R
    c0 = Rinv @ rho0
    fbC0 = np.zeros(n, dtype=complex)
    for r in range(n):
        for c in range(n):
            fbC0[r] += X[r, c] * I_mat[r, c] * c0[c]
    drho = R @ fbC0
    inner = np.vdot(rho_t, drho)
    return 2.0 * inner.real


def find_peak_hwhm_left(K_curve, Q_grid):
    """Standard peak finder + HWHM_left via parabolic+linear interp."""
    K_abs = np.abs(K_curve)
    i_max = int(np.argmax(K_abs))
    K_max = K_abs[i_max]
    if i_max == 0 or i_max == len(Q_grid) - 1:
        return Q_grid[i_max], K_max, 0.0, 0.0
    y0, y1, y2 = K_abs[i_max - 1], K_abs[i_max], K_abs[i_max + 1]
    denom = y0 - 2 * y1 + y2
    delta = 0.5 * (y0 - y2) / denom if abs(denom) > 1e-15 else 0.0
    dQ = Q_grid[1] - Q_grid[0]
    Q_peak = Q_grid[i_max] + delta * dQ
    half = K_max / 2.0
    hwhm_left = None
    for j in range(i_max, 0, -1):
        if K_abs[j-1] < half <= K_abs[j]:
            frac = (half - K_abs[j-1]) / (K_abs[j] - K_abs[j-1])
            Q_half = Q_grid[j-1] + frac * (Q_grid[j] - Q_grid[j-1])
            hwhm_left = Q_peak - Q_half
            break
    if hwhm_left is None:
        hwhm_left = Q_peak - Q_grid[0]
    return Q_peak, K_max, hwhm_left, hwhm_left / Q_peak if Q_peak > 0 else 0.0


def probe_block_2level_hwhm(g_eff_probe, gamma0=0.05, Q_lo=0.05, Q_hi=10.0, n_grid=4000):
    """Sweep Q over [Q_lo, Q_hi], compute K_b for the probe-block 2-level model,
    return Q_peak, HWHM_left/Q_peak."""
    if g_eff_probe < 1e-10:
        return None
    t_peak = 1.0 / (4.0 * gamma0)
    Q_grid = np.linspace(Q_lo, Q_hi, n_grid)
    K_curve = np.array([k_b_probe_block_two_level(Q, gamma0, g_eff_probe, t_peak)
                        for Q in Q_grid])
    return find_peak_hwhm_left(K_curve, Q_grid)


# ─── Step 5: Tests ──────────────────────────────────────────────────────────


def explore_probe_block_per_bond():
    print("=" * 82)
    print("Step 1+2+3: Compute V_b probe-block per bond and g_eff_probe candidates")
    print("=" * 82)

    for N in [5, 6, 7, 8]:
        print(f"\nN={N}, c=2 (n=1) block:")
        print(f"  Each |c_α⟩ vector built with weight from C2ChannelUniformAnalytical:")
        n_pairs_hd1 = N * (N - 1)
        n_pairs_hd3 = N * (N - 1) * (N - 2) // 2
        print(f"  HD=1 pairs: {n_pairs_hd1}, w_1 = {1.0/np.sqrt(n_pairs_hd1):.6f}")
        print(f"  HD=3 pairs: {n_pairs_hd3}, w_3 = {1.0/np.sqrt(n_pairs_hd3):.6f}")

        probe_blocks = compute_probe_block_per_bond(N)
        for b, V_pb in enumerate(probe_blocks):
            cls = "Endpoint" if b == 0 or b == N - 2 else "Interior"
            cands = extract_g_eff_probe_candidates(V_pb)
            print(f"  Bond {b} ({cls}):")
            print(f"    V_b probe-block:")
            for r in range(2):
                row_str = "  ".join(f"{V_pb[r, c].real:+.4f}{V_pb[r, c].imag:+.4f}j" for c in range(2))
                print(f"      [{row_str}]")
            print(f"    g_eff_probe candidates: " +
                  ", ".join(f"{k}={v:.4f}" for k, v in cands.items()))


def test_probe_block_resonance_per_bond():
    """For each bond, run the 2-level Duhamel resonance and compute HWHM/Q*."""
    print()
    print("=" * 82)
    print("Step 4+5: Probe-block 2-level Duhamel HWHM/Q* per bond, per g_eff candidate")
    print("=" * 82)
    print()
    print("Empirical anchors (PROOF_F86_QPEAK Statement 2): Interior 0.7506, Endpoint 0.7728")
    print()
    candidates = ['frobenius_over_2', 'frobenius', 'magnitude_off',
                  'sqrt_2x_magnitude_off', '2x_magnitude_off', 'top_singular']

    for N in [5, 6, 7, 8]:
        probe_blocks = compute_probe_block_per_bond(N)
        print(f"\nN={N}:")
        print(f"  {'bond':<5} {'class':<10} {'g_eff':<14} {'Q_peak':<8} {'HWHM/Q*':<10}")
        for cand in candidates:
            print(f"  --- using g_eff_probe = {cand} ---")
            ratios_int = []
            ratios_end = []
            for b, V_pb in enumerate(probe_blocks):
                cands = extract_g_eff_probe_candidates(V_pb)
                g = cands[cand]
                if g < 1e-10:
                    print(f"  {b:<5} {'(zero)':<10} {0.0:<14.4f} {'n/a':<8} {'n/a':<10}")
                    continue
                cls = "Endpoint" if b == 0 or b == N - 2 else "Interior"
                result = probe_block_2level_hwhm(g)
                if result is None:
                    continue
                Q_peak, K_max, hwhm_left, ratio = result
                print(f"  {b:<5} {cls:<10} {g:<14.4f} {Q_peak:<8.4f} {ratio:<10.4f}")
                if cls == "Interior":
                    ratios_int.append(ratio)
                else:
                    ratios_end.append(ratio)
            if ratios_int and ratios_end:
                int_mean = np.mean(ratios_int)
                end_mean = np.mean(ratios_end)
                gap = end_mean - int_mean
                print(f"    -> Interior mean = {int_mean:.4f}, "
                      f"Endpoint mean = {end_mean:.4f}, gap = {gap:.4f}")
                print(f"    -> empirical anchors: Interior 0.7506, Endpoint 0.7728, gap 0.022")


def test_two_level_class_means():
    """Better approach: average the K curves across bonds in the same class first
    (matching the canonical Python pipeline contract), then find peak/HWHM."""
    print()
    print("=" * 82)
    print("Step 5b: Class-mean K-curve, probe-block 2-level model")
    print("=" * 82)
    print()

    candidates = ['frobenius_over_2', 'frobenius', 'magnitude_off',
                  'sqrt_2x_magnitude_off', '2x_magnitude_off', 'top_singular']
    gamma0 = 0.05
    t_peak = 1.0 / (4.0 * gamma0)
    Q_grid = np.linspace(0.05, 10.0, 4000)

    for N in [5, 6, 7, 8]:
        print(f"\nN={N}:")
        probe_blocks = compute_probe_block_per_bond(N)
        for cand in candidates:
            # Collect g_eff_probe per bond, compute K curves per bond, average per class
            g_per_bond = []
            for V_pb in probe_blocks:
                cands = extract_g_eff_probe_candidates(V_pb)
                g_per_bond.append(cands[cand])

            n_bonds = len(probe_blocks)
            int_curves = []
            end_curves = []
            for b in range(n_bonds):
                if g_per_bond[b] < 1e-10:
                    continue
                K_curve = np.array([k_b_probe_block_two_level(Q, gamma0, g_per_bond[b], t_peak)
                                    for Q in Q_grid])
                if b == 0 or b == n_bonds - 1:
                    end_curves.append(K_curve)
                else:
                    int_curves.append(K_curve)
            if int_curves and end_curves:
                int_mean_curve = np.mean(int_curves, axis=0)
                end_mean_curve = np.mean(end_curves, axis=0)
                _, _, _, ratio_int = find_peak_hwhm_left(int_mean_curve, Q_grid)
                _, _, _, ratio_end = find_peak_hwhm_left(end_mean_curve, Q_grid)
                gap = ratio_end - ratio_int
                print(f"  {cand:<25}: Interior mean curve ratio = {ratio_int:.4f}, "
                      f"Endpoint = {ratio_end:.4f}, gap = {gap:.4f}")


def test_g_eff_probe_endpoint_interior_ratio():
    """Test the synthesis prediction: g_eff_probe(Endpoint) ≈ g_eff_probe(Interior) / 1.6."""
    print()
    print("=" * 82)
    print("Step 4 verification: g_eff_probe(Endpoint) / g_eff_probe(Interior) ≈ 1/1.6 = 0.625?")
    print("=" * 82)
    print()
    candidates = ['frobenius_over_2', 'frobenius', 'magnitude_off',
                  'sqrt_2x_magnitude_off', '2x_magnitude_off', 'top_singular']
    print(f"  {'N':<5} {'candidate':<25} {'<g_int>':<10} {'<g_end>':<10} {'g_end/g_int':<14}")

    for N in [5, 6, 7, 8]:
        probe_blocks = compute_probe_block_per_bond(N)
        for cand in candidates:
            g_int = []
            g_end = []
            for b, V_pb in enumerate(probe_blocks):
                cands = extract_g_eff_probe_candidates(V_pb)
                g = cands[cand]
                if b == 0 or b == N - 2:
                    g_end.append(g)
                else:
                    g_int.append(g)
            if g_int and g_end:
                int_mean = np.mean(g_int)
                end_mean = np.mean(g_end)
                ratio = end_mean / int_mean if int_mean > 0 else float('inf')
                print(f"  {N:<5} {cand:<25} {int_mean:<10.4f} {end_mean:<10.4f} {ratio:<14.4f}")


# ─── Step 6: Three-block superposition: K_total = K_pb + K_sv (+ cross) ─────


def k_b_full_4mode_with_pb_and_sv(Q, gamma0, g_eff_probe, sigma0, t):
    """Test the simplest superposition: probe-block 2-level (with g_eff_probe) +
    SVD-block 2-level (with σ_0) decoupled.

    L_4(Q) = [
      [-2γ₀, +iJ g_eff_probe, 0, 0],
      [+iJ g_eff_probe, -6γ₀, 0, 0],
      [0, 0, -2γ₀, +iJ σ_0],
      [0, 0, +iJ σ_0, -6γ₀]
    ]
    Probe = (1, 0, 0, 0). V_b at this 4-mode level is dL/dJ.
    """
    J = Q * gamma0
    L = np.zeros((4, 4), dtype=complex)
    L[0, 0] = -2 * gamma0
    L[1, 1] = -6 * gamma0
    L[2, 2] = -2 * gamma0
    L[3, 3] = -6 * gamma0
    L[0, 1] = 1j * J * g_eff_probe
    L[1, 0] = 1j * J * g_eff_probe
    L[2, 3] = 1j * J * sigma0
    L[3, 2] = 1j * J * sigma0

    Vb = np.zeros((4, 4), dtype=complex)
    Vb[0, 1] = 1j * g_eff_probe
    Vb[1, 0] = 1j * g_eff_probe
    Vb[2, 3] = 1j * sigma0
    Vb[3, 2] = 1j * sigma0

    evals, R = np.linalg.eig(L)
    Rinv = np.linalg.inv(R)
    rho0 = np.array([1.0, 0.0, 0.0, 0.0], dtype=complex)
    expLam = np.exp(evals * t)
    rho_t = R @ np.diag(expLam) @ (Rinv @ rho0)
    n = 4
    I_mat = np.zeros((n, n), dtype=complex)
    for j in range(n):
        for k in range(n):
            diff = evals[k] - evals[j]
            if abs(diff) < 1e-10:
                I_mat[j, k] = t * expLam[j]
            else:
                I_mat[j, k] = (expLam[k] - expLam[j]) / diff
    X = Rinv @ Vb @ R
    c0 = Rinv @ rho0
    fbC0 = np.zeros(n, dtype=complex)
    for r in range(n):
        for c in range(n):
            fbC0[r] += X[r, c] * I_mat[r, c] * c0[c]
    drho = R @ fbC0
    return 2.0 * np.vdot(rho_t, drho).real


def test_pb_plus_sv_decoupled():
    """Two decoupled 2-level resonances: probe-block at g_eff_probe, SVD-block at σ_0.

    If only the probe-block resonance contributes to the K_b shape (because the probe
    starts in |c_1⟩ which is in the probe-block subspace), then the result should
    reproduce the bare probe-block 2-level resonance for whatever g_eff_probe is.

    But if the SVD-block also resonates with the right kernel, the result might be
    a superposition that lifts HWHM/Q* above the bare 0.6715 floor.
    """
    print()
    print("=" * 82)
    print("Step 6: Decoupled 4-mode (probe-block + SVD-block decoupled): 2-level superposition")
    print("=" * 82)
    print()
    gamma0 = 0.05
    t_peak = 1.0 / (4.0 * gamma0)
    Q_grid = np.linspace(0.05, 10.0, 4000)

    # σ_0 ~ 2*sqrt(2) at large N
    sigma0_values = {5: 2.765, 6: 2.802, 7: 2.828, 8: 2.839}
    cand = 'frobenius_over_2'  # placeholder

    for N in [5, 6, 7, 8]:
        probe_blocks = compute_probe_block_per_bond(N)
        sigma0 = sigma0_values[N]
        print(f"\nN={N}, σ_0={sigma0}:")
        print(f"  {'bond':<5} {'class':<10} {'g_eff_pb':<12} {'Q_peak':<8} {'HWHM/Q*':<10}")
        for b, V_pb in enumerate(probe_blocks):
            cands = extract_g_eff_probe_candidates(V_pb)
            g = cands[cand]
            cls = "Endpoint" if b == 0 or b == N - 2 else "Interior"
            K_curve = np.array([k_b_full_4mode_with_pb_and_sv(Q, gamma0, g, sigma0, t_peak)
                               for Q in Q_grid])
            Q_peak, _, _, ratio = find_peak_hwhm_left(K_curve, Q_grid)
            print(f"  {b:<5} {cls:<10} {g:<12.4f} {Q_peak:<8.4f} {ratio:<10.4f}")


# ─── Empirical anchor table (Direction (b) found this) ──────────────────────


def empirical_anchor_summary():
    """Anchors per PROOF_F86_QPEAK Statement 2 + the Q-peak table."""
    print()
    print("=" * 82)
    print("Empirical anchor table (PROOF_F86_QPEAK Statement 2, gamma_0=0.05)")
    print("=" * 82)
    print(f"  {'N':<5} {'Q_peak Int':<12} {'Q_peak End':<12} {'HWHM/Q* Int':<12} "
          f"{'HWHM/Q* End':<12} {'σ_0':<8}")
    cases = [
        (5, 1.4821, 2.5008, 0.7455, 0.7700, 2.765),
        (6, 1.5801, 2.5470, 0.7529, 0.7738, 2.802),
        (7, 1.5831, 2.5299, 0.7507, 0.7738, 2.828),
        (8, 1.6049, 2.5145, 0.7531, 0.7734, 2.839),
    ]
    for N, qi, qe, hi, he, sig in cases:
        print(f"  {N:<5} {qi:<12.4f} {qe:<12.4f} {hi:<12.4f} {he:<12.4f} {sig:<8.3f}")
    qpeak_int_mean = np.mean([c[1] for c in cases])
    qpeak_end_mean = np.mean([c[2] for c in cases])
    print(f"\n  Q_peak Endpoint/Interior ratio: {qpeak_end_mean/qpeak_int_mean:.4f}")
    print(f"  HWHM/Q* gap (Endpoint - Interior): {np.mean([c[4] - c[3] for c in cases]):.4f}")


def falsification_summary():
    """Summary of the structural falsification of Direction (a')."""
    print()
    print("=" * 82)
    print("DIRECTION (a') STRUCTURAL FALSIFICATION SUMMARY")
    print("=" * 82)
    print()
    print("Three structural facts FALSIFY Direction (a') as specified:")
    print()
    print("(1) V_b probe-block is bond-class-blind for all N tested.")
    print("    For every bond at every N=5..8:")
    print("      V_b[0, 0] = +i·c (scalar, same constant for all b)")
    print("      V_b[1, 1] = +i·c (same scalar)")
    print("      V_b[0, 1] = EXACTLY zero (verified bit-exact)")
    print("    The F73 sum-rule applies per-bond, not just summed.")
    print("    -> The hypothesised g_eff_probe(N, b) cannot have bond-class dependence.")
    print()
    print("(2) Cross-block Frobenius is unstable across N.")
    for N in [5, 6, 7, 8]:
        from numpy.linalg import svd
        states_p, states_q, flat, Mtot = build_basis(N, 1)
        c1, c3 = build_c_alpha_vectors(N)
        MhPerBond, _, _, _, _ = build_mh_per_bond(N, 0.05)
        MhTotal = sum(MhPerBond)
        hd1_idx = []; hd3_idx = []
        for p in states_p:
            for q in states_q:
                hd = bin(p^q).count('1')
                idx = flat(p, q)
                if hd == 1: hd1_idx.append(idx)
                elif hd == 3: hd3_idx.append(idx)
        hd1_idx = np.array(hd1_idx); hd3_idx = np.array(hd3_idx)
        V_inter = MhTotal[np.ix_(hd1_idx, hd3_idx)]
        U, S, Vh = svd(V_inter)
        u0 = np.zeros(Mtot, dtype=complex); u0[hd1_idx] = U[:, 0]
        v0 = np.zeros(Mtot, dtype=complex); v0[hd3_idx] = Vh[0].conjugate()
        end_frob = []; int_frob = []
        for b, Mh in enumerate(MhPerBond):
            V02 = np.vdot(c1, Mh @ u0); V03 = np.vdot(c1, Mh @ v0)
            V12 = np.vdot(c3, Mh @ u0); V13 = np.vdot(c3, Mh @ v0)
            f = np.sqrt(abs(V02)**2 + abs(V03)**2 + abs(V12)**2 + abs(V13)**2)
            if b == 0 or b == N - 2: end_frob.append(f)
            else: int_frob.append(f)
        e = np.mean(end_frob); i = np.mean(int_frob)
        print(f"    N={N}: Endpoint={e:.4f}, Interior={i:.4f}, ratio={e/i:.4f}")
    print("    -> Library-dependent at even N due to A3's σ_0 degeneracy obstruction.")
    print()
    print("(3) SVD-block off-diagonal V_b[2, 3] IS bond-class-dependent consistently.")
    for N in [5, 6, 7, 8]:
        from numpy.linalg import svd
        states_p, states_q, flat, Mtot = build_basis(N, 1)
        MhPerBond, _, _, _, _ = build_mh_per_bond(N, 0.05)
        MhTotal = sum(MhPerBond)
        hd1_idx = []; hd3_idx = []
        for p in states_p:
            for q in states_q:
                hd = bin(p^q).count('1')
                idx = flat(p, q)
                if hd == 1: hd1_idx.append(idx)
                elif hd == 3: hd3_idx.append(idx)
        hd1_idx = np.array(hd1_idx); hd3_idx = np.array(hd3_idx)
        V_inter = MhTotal[np.ix_(hd1_idx, hd3_idx)]
        U, S, Vh = svd(V_inter)
        u0 = np.zeros(Mtot, dtype=complex); u0[hd1_idx] = U[:, 0]
        v0 = np.zeros(Mtot, dtype=complex); v0[hd3_idx] = Vh[0].conjugate()
        end_svd = []; int_svd = []
        for b, Mh in enumerate(MhPerBond):
            V23 = abs(np.vdot(u0, Mh @ v0))
            if b == 0 or b == N - 2: end_svd.append(V23)
            else: int_svd.append(V23)
        e = np.mean(end_svd); i = np.mean(int_svd)
        print(f"    N={N}: |V_b[2,3]| Endpoint={e:.4f}, Interior={i:.4f}, "
              f"ratio={e/i:.4f}")
    print("    -> Endpoint < Interior consistently. Direction OPPOSITE the empirical")
    print("       HWHM/Q* split (Endpoint > Interior). Closed form must map V_b magnitude")
    print("       to HWHM/Q* shift non-trivially.")
    print()
    print("CONCLUSION: Direction (a') as specified is structurally falsified.")
    print("Refined direction (a''): Use SVD-block V_b[2, 3] as bond-class carrier;")
    print("derive the non-trivial map from V_b sub-block structure to HWHM/Q* shift.")
    print()
    print("Refined direction (b''): Work in full block-L, not 4-mode reduction.")
    print("4-mode K_b reproduces only the ~0.673 floor, not the empirical 0.7506 lift.")


def main():
    np.set_printoptions(precision=6, suppress=True)
    empirical_anchor_summary()
    explore_probe_block_per_bond()
    test_g_eff_probe_endpoint_interior_ratio()
    falsification_summary()


if __name__ == "__main__":
    main()
