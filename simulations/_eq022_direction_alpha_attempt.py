"""Direction (alpha) polarity-Bloch projection at t_peak: c=2 per-bond r(N, b) closed-form attempt.

Predecessor work (PolarityInheritanceLink Tier2Verified, commit 7886946 / 12f3181):
  - F86 c=2 bond-class split (Endpoint vs Interior in Q_peak and HWHM/Q*) inherits structurally
    from the polarity-layer pair {-0.5, +0.5} at d=2 named in PolarityLayerOriginClaim
    (Pi2KnowledgeBase): rho = (I + r sigma)/2 maps Pi-spectrum {+1, -1} onto Bloch-diagonal
    1/2 +- r/2.
  - Empirical witnesses (gamma_0 = 0.05, c=2 N=5..8):
      N | r_Q Endpoint | r_Q Interior | r_H Endpoint | r_H Interior
      5 | +0.5008      | -0.5179      | 0.5400       | 0.4910
      6 | +0.5470      | -0.4199      | 0.5476       | 0.5058
      7 | +0.5299      | -0.4169      | 0.5476       | 0.5014
      8 | +0.5145      | -0.3951      | 0.5468       | 0.5062
    where r_Q = Q_peak - 2 (from Q_peak axis), r_H = 2*(HWHM/Q* - 1/2) (from HWHM/Q* axis).

Goal of Direction (alpha):
  Derive r(N, b) analytically by:
  1) Computing the K-driving 4-mode eigenstate rho_K(Q_peak, t_peak, b) per bond.
  2) Defining a "polarity Bloch axis" in the 4-mode basis and projecting rho_K onto it.
  3) Comparing r_polarity to empirical r_Q and r_H.
  4) If a clean structural form lands, derive r(N, b) closed form.

Approach:
  - Use the 4-mode L_eff(Q) = D_eff + Q gamma_0 * sum_b V_b in basis B = [c_1, c_3, u_0, v_0]
  - Probe lives in span{c_1, c_3} (Tier 1 structural fact; SVD components are 0).
  - At c=2 the spectrum is bond-INDEPENDENT; the bond enters via per-bond V_b in dL/dJ_b.
  - Time evolution: rho(t) = exp(L_eff(Q) * t) * rho_probe.
  - Restriction: at t_peak = 1/(4 gamma_0), peak Q = Q_peak(b), evaluate rho_K(t_peak).

Polarity axis candidates explored (3 categories):
  (A) Channel-uniform polarity: project rho_K onto (|c_1> +- |c_3>)/sqrt(2)
      ("X-like" split between HD=1 and HD=3 channels).
  (B) Bond-coupling-influenced polarity: project rho_K onto bond-coupling-sensitive direction
      (e.g. M_h_total_eff axis of L_eff or per-bond V_b axis).
  (C) Mixed-block polarity: project onto u_0 +- v_0 (SVD-block X-like split, parallel to A).

Step 0: stdout fix + imports.
Step 1: build the c=2 4-mode L_eff(Q) numerically (using same Pauli/popcount basis).
Step 2: per-bond compute rho_K(Q_peak, t_peak), with Q_peak = empirical anchor.
Step 3: project onto the polarity axes (A, B, C), tabulate r_polarity per bond.
Step 4: compare r_polarity to empirical r_Q and r_H; check tolerance.
Step 5: if structural form lands, write r(N, b) closed form; else document partial finding.
"""

import sys
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import os

import numpy as np
from numpy.linalg import eig, inv, svd, norm

# Use the framework primitive for the F73 spatial-sum coherence kernel rather
# than a placeholder identity matrix. (Code-review fix, 2026-05-08.)
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from framework.coherence_block import spatial_sum_coherence_kernel


# Empirical anchor table from PolarityInheritanceLink witnesses (gamma_0 = 0.05, c=2 N=5..8)
# Values in (Q_peak_Interior, Q_peak_Endpoint, HwhmRatio_Interior, HwhmRatio_Endpoint).
EMPIRICAL = {
    5: dict(Q_int=1.4821, Q_end=2.5008, H_int=0.7455, H_end=0.7700),
    6: dict(Q_int=1.5801, Q_end=2.5470, H_int=0.7529, H_end=0.7738),
    7: dict(Q_int=1.5831, Q_end=2.5299, H_int=0.7507, H_end=0.7738),
    8: dict(Q_int=1.6049, Q_end=2.5145, H_int=0.7531, H_end=0.7734),
}


def empirical_r(N):
    """Return (r_Q_end, r_Q_int, r_H_end, r_H_int) for the given N."""
    e = EMPIRICAL[N]
    return (
        e["Q_end"] - 2.0,
        e["Q_int"] - 2.0,
        2.0 * (e["H_end"] - 0.5),
        2.0 * (e["H_int"] - 0.5),
    )


# ---------- Block construction ----------------------------------------------

def build_basis(N):
    """Build the (n=1, n+1=2) popcount block basis for c=2."""
    states_p = sorted([s for s in range(1 << N) if bin(s).count("1") == 1])
    states_q = sorted([s for s in range(1 << N) if bin(s).count("1") == 2])
    Mp, Mq = len(states_p), len(states_q)
    p_idx = {p: i for i, p in enumerate(states_p)}
    q_idx = {q: j for j, q in enumerate(states_q)}

    def flat(p, q):
        return p_idx[p] * Mq + q_idx[q]

    return states_p, states_q, flat, Mp * Mq


def hamming_distance(p, q):
    return bin(p ^ q).count("1")


def bond_flip_targets(state, N):
    """Big-endian bond convention: site 0 = MSB; bond b spans sites (b, b+1)."""
    for b in range(N - 1):
        bit_b = (state >> (N - 1 - b)) & 1
        bit_bp1 = (state >> (N - 1 - (b + 1))) & 1
        if bit_b != bit_bp1:
            mask = (1 << (N - 1 - b)) | (1 << (N - 1 - (b + 1)))
            yield b, state ^ mask


def build_full_block_operators(N, gamma0):
    """Construct M_h_per_bond[b], M_h_total, D, S_kernel, probe in the c=2 (1,2) block.

    Conventions:
      H_b = (J/2)(X_b X_{b+1} + Y_b Y_{b+1}); each bond flips exactly one differing pair.
      M_h_per_bond[b][j, i] += -i for p-flip; +i for q-flip.
      D = -2 gamma_0 on HD=1 indices, -6 gamma_0 on HD=3 indices, etc.; in the (1,2) block
        D[i, i] = -gamma_0 * (n_p^2 + n_q^2 - n_p - n_q)... actually for the c=2 block the
        relevant rates by HD are: HD=1: 2 gamma_0, HD=3: 6 gamma_0 (from F86 Statement 1).
    """
    states_p, states_q, flat, Mtot = build_basis(N)

    M_h_per_bond = [np.zeros((Mtot, Mtot), dtype=complex) for _ in range(N - 1)]
    for p in states_p:
        for q in states_q:
            i = flat(p, q)
            for b, p_flipped in bond_flip_targets(p, N):
                if bin(p_flipped).count("1") == 1:
                    j = flat(p_flipped, q)
                    M_h_per_bond[b][j, i] += -1j
            for b, q_flipped in bond_flip_targets(q, N):
                if bin(q_flipped).count("1") == 2:
                    j = flat(p, q_flipped)
                    M_h_per_bond[b][j, i] += +1j

    M_h_total = sum(M_h_per_bond)

    # D: the full Z-dephasing dissipator restricted to the (1, 2) block.
    # Convention from BlockLDecomposition.cs / framework.coherence_block.block_L_split_xy:
    #   D[i, i] = -2·gamma_0·HD(p, q)
    # NOT -gamma_0·HD². (Code-review bug fix, 2026-05-08; was -gamma_0·HD² before.)
    D = np.zeros((Mtot, Mtot), dtype=complex)
    for p in states_p:
        for q in states_q:
            i = flat(p, q)
            hd = bin(p ^ q).count("1")
            D[i, i] = -2.0 * gamma0 * hd

    # Probe: Dicke (n=1) state to (n=2) state coupling. The Dicke probe in the (1,2) block
    # is the (proportional) symmetric-superposition factored as |Dicke_1> tensor |Dicke_2>:
    # uniform weight on every (p, q) basis vector with p of popcount 1 and q of popcount 2.
    probe = np.ones(Mtot, dtype=complex)
    probe = probe / norm(probe)

    # S_kernel: the F73 spatial-sum coherence kernel computed via the framework primitive.
    # (Code-review bug fix, 2026-05-08; was identity before, which silently changed the
    # observable being maximised. With the proper kernel, K_b(Q, t_peak) = ρ(t)† S ∂_J ρ.)
    # The (1, 2) block corresponds to n = 1 in the framework convention.
    S_kernel = spatial_sum_coherence_kernel(N, n=1)

    return states_p, states_q, flat, Mtot, M_h_per_bond, M_h_total, D, probe, S_kernel


# ---------- 4-mode basis and effective L --------------------------------------


def build_4mode_basis(states_p, states_q, flat, Mtot, M_h_total):
    """Build the 4-mode basis B = [c_1, c_3, u_0, v_0]."""
    indices_hd1 = []
    indices_hd3 = []
    for p in states_p:
        for q in states_q:
            i = flat(p, q)
            hd = bin(p ^ q).count("1")
            if hd == 1:
                indices_hd1.append(i)
            elif hd == 3:
                indices_hd3.append(i)

    c1 = np.zeros(Mtot, dtype=complex)
    c3 = np.zeros(Mtot, dtype=complex)
    if indices_hd1:
        w1 = 1.0 / np.sqrt(len(indices_hd1))
        for idx in indices_hd1:
            c1[idx] = w1
    if indices_hd3:
        w3 = 1.0 / np.sqrt(len(indices_hd3))
        for idx in indices_hd3:
            c3[idx] = w3

    # u_0, v_0 from inter-channel SVD of M_h_total restricted to HD=1 x HD=3 block.
    hd1_arr = np.array(indices_hd1)
    hd3_arr = np.array(indices_hd3)
    V_inter = M_h_total[np.ix_(hd1_arr, hd3_arr)]
    U, S, Vh = svd(V_inter)
    u0 = np.zeros(Mtot, dtype=complex)
    v0 = np.zeros(Mtot, dtype=complex)
    u0[hd1_arr] = U[:, 0]
    v0[hd3_arr] = Vh[0].conjugate()

    B = np.column_stack([c1, c3, u0, v0])
    # Verify orthonormality
    G = B.conj().T @ B
    off_max = max(abs(G[i, j]) for i in range(4) for j in range(4) if i != j)
    if off_max > 1e-8:
        print(f"[warning] basis off-orthonormality residual = {off_max:.3e}")

    return B, c1, c3, u0, v0


def project_4mode(B, M):
    """Project a Mtot x Mtot matrix onto the 4-mode basis."""
    return B.conj().T @ M @ B


def project_vec_4mode(B, v):
    """Project a Mtot vector onto the 4-mode basis."""
    return B.conj().T @ v


# ---------- K-driving eigenstate ---------------------------------------------


def k_driving_subspace_at_q(L_eff, probe_eff):
    """Return the 4 eigenvalues of L_eff sorted Re desc / Im asc, the eigenvector matrix R
    in matching column order, and the 2 K-driving indices ranked by |<probe | w_i>|^2."""
    evals_raw, R_raw = eig(L_eff)
    # Sort by Re desc, Im asc
    idx = sorted(range(4), key=lambda k: (-evals_raw[k].real, evals_raw[k].imag))
    evals = evals_raw[idx]
    R = R_raw[:, idx]

    overlaps = np.zeros(4, dtype=float)
    for i in range(4):
        inner = np.vdot(probe_eff, R[:, i])
        overlaps[i] = abs(inner) ** 2

    # Top 2 by overlap (descending)
    pair_indices = sorted(range(4), key=lambda k: -overlaps[k])[:2]
    # Order so that LamPlus has larger Re
    if evals[pair_indices[0]].real < evals[pair_indices[1]].real:
        pair_indices = pair_indices[::-1]

    return evals, R, pair_indices, overlaps


def evolve_probe_to_t(L_eff, probe_eff, t):
    """rho(t) = exp(L_eff t) * probe_eff (4-mode)."""
    evals, R = eig(L_eff)
    Rinv = inv(R)
    expLam = np.exp(evals * t)
    return R @ (expLam * (Rinv @ probe_eff))


def project_to_k_driving_pair(rho_t_4mode, R, pair_indices):
    """Restrict rho_t to the 2D K-driving eigenvector subspace by keeping only the
    pair_indices components (in the eigenbasis), zeroing the rest, then transforming back."""
    Rinv = inv(R)
    coords = Rinv @ rho_t_4mode
    # Restrict
    restricted = np.zeros_like(coords)
    restricted[pair_indices[0]] = coords[pair_indices[0]]
    restricted[pair_indices[1]] = coords[pair_indices[1]]
    return R @ restricted


# ---------- Polarity axes ----------------------------------------------------


def polarity_axis_channel_uniform_plus(B, c1, c3):
    """(c_1 + c_3)/sqrt(2);channel-uniform 'positive polarity'."""
    return (c1 + c3) / np.sqrt(2.0)


def polarity_axis_channel_uniform_minus(B, c1, c3):
    """(c_1 - c_3)/sqrt(2);channel-uniform 'negative polarity'."""
    return (c1 - c3) / np.sqrt(2.0)


def polarity_axis_4mode_pp(B):
    """In 4-mode basis, the (e_0 + e_1)/sqrt(2) direction.  Projected back: (c_1 + c_3)/sqrt(2)."""
    p = np.array([1.0, 1.0, 0.0, 0.0]) / np.sqrt(2.0)
    return p


def polarity_axis_4mode_pm(B):
    """In 4-mode basis, the (e_0 - e_1)/sqrt(2) direction.  Projected back: (c_1 - c_3)/sqrt(2)."""
    p = np.array([1.0, -1.0, 0.0, 0.0]) / np.sqrt(2.0)
    return p


# ---------- Per-bond evaluation ----------------------------------------------


def per_bond_polarity_evaluation(N, gamma0, Q_peak_per_class, t_peak):
    """For the given N at c=2, build the 4-mode picture, evolve probe to t_peak at Q = Q_peak,
    project onto the polarity Bloch axes per bond, and report r values.

    Q_peak_per_class: dict with 'Endpoint' -> Q_peak_end, 'Interior' -> Q_peak_int.

    Returns one row per bond with (bond, class, r_polarity_pp, r_polarity_pm).
    """
    states_p, states_q, flat, Mtot, M_h_per_bond, M_h_total, D, probe_full, S_kernel \
        = build_full_block_operators(N, gamma0)

    B, c1, c3, u0, v0 = build_4mode_basis(states_p, states_q, flat, Mtot, M_h_total)

    D_eff = project_4mode(B, D)
    M_h_total_eff = project_4mode(B, M_h_total)
    probe_eff = project_vec_4mode(B, probe_full)

    rows = []
    for b in range(N - 1):
        is_end = (b == 0 or b == N - 2)
        cls = "Endpoint" if is_end else "Interior"
        Q_peak = Q_peak_per_class["Endpoint"] if is_end else Q_peak_per_class["Interior"]
        L_eff_at_Q = D_eff + Q_peak * gamma0 * M_h_total_eff

        evals, R, pair_indices, overlaps = k_driving_subspace_at_q(L_eff_at_Q, probe_eff)
        rho_t = evolve_probe_to_t(L_eff_at_Q, probe_eff, t_peak)
        rho_K = project_to_k_driving_pair(rho_t, R, pair_indices)

        # Polarity axes in the 4-mode basis directly
        p_pp = polarity_axis_4mode_pp(B)
        p_pm = polarity_axis_4mode_pm(B)

        # Project rho_K onto the polarity axis: r = <p | rho_K> / norm(rho_K)
        # The 'state-Bloch' interpretation: rho is a state vector (in operator space here),
        # the polarity axis is one of the d^2 'directions'. r is the (real-part of the
        # complex projection) / norm;the signed Bloch coordinate along the axis.
        nrm = norm(rho_K)
        if nrm > 1e-12:
            r_pp = np.vdot(p_pp, rho_K).real / nrm
            r_pm = np.vdot(p_pm, rho_K).real / nrm
            # Also unsigned absolute components:
            r_pp_signed = np.vdot(p_pp, rho_K).real
            r_pm_signed = np.vdot(p_pm, rho_K).real
        else:
            r_pp = r_pm = r_pp_signed = r_pm_signed = 0.0

        # Alternate: Bloch coordinate of rho(t)/rho(0) along (e_0 - e_1)/sqrt(2)
        # in the 4-mode basis (raw, not subspace-restricted).
        rho_t_unrestricted = rho_t
        # decomposition of probe_eff
        p_e0 = np.array([1.0, 0.0, 0.0, 0.0])
        p_e1 = np.array([0.0, 1.0, 0.0, 0.0])
        amp_c1 = np.vdot(p_e0, rho_t_unrestricted).real
        amp_c3 = np.vdot(p_e1, rho_t_unrestricted).real
        amp_c1_init = np.vdot(p_e0, probe_eff).real
        amp_c3_init = np.vdot(p_e1, probe_eff).real

        r_diag_signed = (amp_c1 - amp_c3)
        r_diag_norm = (amp_c1 - amp_c3) / max(abs(amp_c1) + abs(amp_c3), 1e-12)

        rows.append({
            "N": N,
            "bond": b,
            "class": cls,
            "Q_peak": Q_peak,
            "evals": evals,
            "pair_indices": pair_indices,
            "overlaps": overlaps,
            "r_pp_normalized": r_pp,
            "r_pm_normalized": r_pm,
            "r_pp_signed": r_pp_signed,
            "r_pm_signed": r_pm_signed,
            "r_diag_signed": r_diag_signed,
            "r_diag_norm": r_diag_norm,
            "amp_c1": amp_c1,
            "amp_c3": amp_c3,
            "amp_c1_init": amp_c1_init,
            "amp_c3_init": amp_c3_init,
        })

    return rows


# ---------- Main exploration -------------------------------------------------


def main():
    print("=" * 90)
    print("Direction (alpha) attempt: per-bond polarity-Bloch projection at t_peak.")
    print("=" * 90)
    print()

    gamma0 = 0.05
    t_peak = 1.0 / (4.0 * gamma0)  # = 5.0
    print(f"gamma_0 = {gamma0}, t_peak = 1/(4 gamma_0) = {t_peak}")
    print()
    print("Empirical anchor table (PolarityInheritanceLink witnesses, c=2 N=5..8):")
    print(f"  {'N':<3} {'r_Q End':<10} {'r_Q Int':<10} {'r_H End':<10} {'r_H Int':<10}")
    for N in (5, 6, 7, 8):
        rQE, rQI = EMPIRICAL[N]["Q_end"] - 2, EMPIRICAL[N]["Q_int"] - 2
        rHE = 2 * (EMPIRICAL[N]["H_end"] - 0.5)
        rHI = 2 * (EMPIRICAL[N]["H_int"] - 0.5)
        print(f"  {N:<3} {rQE:<+10.4f} {rQI:<+10.4f} {rHE:<+10.4f} {rHI:<+10.4f}")
    print()

    print("=" * 90)
    print("Step 1: per-bond polarity projection of rho_K(Q_peak, t_peak)")
    print("=" * 90)
    print()
    print("Polarity axes tried in 4-mode basis (e_0 = c_1, e_1 = c_3):")
    print("  (A_pp) (e_0 + e_1)/sqrt(2) ;channel-uniform '+' polarity")
    print("  (A_pm) (e_0 - e_1)/sqrt(2) ;channel-uniform '-' polarity")
    print("  (D)    diag-amplitude difference (amp_c1 - amp_c3)")
    print()

    for N in (5, 6, 7, 8):
        rQE, rQI, rHE, rHI = empirical_r(N)
        Q_peak_per_class = {
            "Endpoint": EMPIRICAL[N]["Q_end"],
            "Interior": EMPIRICAL[N]["Q_int"],
        }
        rows = per_bond_polarity_evaluation(N, gamma0, Q_peak_per_class, t_peak)

        print(f"--- N = {N} ---")
        print(f"  empirical: r_Q_end = {rQE:+.4f}, r_Q_int = {rQI:+.4f}, "
              f"r_H_end = {rHE:+.4f}, r_H_int = {rHI:+.4f}")
        print(f"  {'b':<3} {'class':<10} {'Q_peak':<8} {'amp_c1(t)':<12} {'amp_c3(t)':<12} "
              f"{'r_diag':<10} {'r_pp':<10} {'r_pm':<10}")
        for r in rows:
            print(f"  {r['bond']:<3} {r['class']:<10} {r['Q_peak']:<8.4f} "
                  f"{r['amp_c1']:<+12.4f} {r['amp_c3']:<+12.4f} "
                  f"{r['r_diag_signed']:<+10.4f} {r['r_pp_signed']:<+10.4f} "
                  f"{r['r_pm_signed']:<+10.4f}")
        print()

    print("=" * 90)
    print("Step 2: aggregated per-class polarity projection (mean across bonds in class)")
    print("=" * 90)
    print()
    print(f"  {'N':<3} {'r_Q_E':<10} {'r_Q_I':<10} {'r_H_E':<10} {'r_H_I':<10} "
          f"{'<r_diag>_E':<12} {'<r_diag>_I':<12} {'<r_pm>_E':<12} {'<r_pm>_I':<12}")
    for N in (5, 6, 7, 8):
        rQE, rQI, rHE, rHI = empirical_r(N)
        Q_peak_per_class = {
            "Endpoint": EMPIRICAL[N]["Q_end"],
            "Interior": EMPIRICAL[N]["Q_int"],
        }
        rows = per_bond_polarity_evaluation(N, gamma0, Q_peak_per_class, t_peak)

        end_diag = [r["r_diag_signed"] for r in rows if r["class"] == "Endpoint"]
        int_diag = [r["r_diag_signed"] for r in rows if r["class"] == "Interior"]
        end_pm = [r["r_pm_signed"] for r in rows if r["class"] == "Endpoint"]
        int_pm = [r["r_pm_signed"] for r in rows if r["class"] == "Interior"]
        eDe = float(np.mean(end_diag)) if end_diag else float("nan")
        iDe = float(np.mean(int_diag)) if int_diag else float("nan")
        ePm = float(np.mean(end_pm)) if end_pm else float("nan")
        iPm = float(np.mean(int_pm)) if int_pm else float("nan")
        print(f"  {N:<3} {rQE:<+10.4f} {rQI:<+10.4f} {rHE:<+10.4f} {rHI:<+10.4f} "
              f"{eDe:<+12.4f} {iDe:<+12.4f} {ePm:<+12.4f} {iPm:<+12.4f}")
    print()

    print("=" * 90)
    print("Step 3: try OBC sine-mode candidate r_polarity(N, b) closed form")
    print("=" * 90)
    print()
    print("OBC sine modes: psi_k(j) = sqrt(2/(N+1)) * sin(pi k (j+1)/(N+1)).")
    print("Bond-position dispersion: 2 cos(pi k / (N+1)).")
    print("F80 'sign-walk' shows up as the OBC dispersion on chain.")
    print()
    print("Candidate (1): r(N, b) = +/- sin(pi (b+1)/(N+1)) - sin(pi (N-1-b)/(N+1)) ?")
    print("Candidate (2): r(N, b) = (-1)^b * 2 * cos(pi/(N+1)) ?")
    print("Candidate (3): r(N, b) = +/- f(N) - g(b)?")
    print()
    print("Compute the OBC sine first-mode amplitude per bond and tabulate vs empirical r_Q.")
    print()
    print(f"  {'N':<3} {'b':<3} {'class':<10} {'2cos(pi/(N+1))':<16} {'sin(pi(b+1)/(N+1))':<22} "
          f"{'OBC mode at b':<20}")
    for N in (5, 6, 7, 8):
        for b in range(N - 1):
            cls = "Endpoint" if (b == 0 or b == N - 2) else "Interior"
            two_cos = 2.0 * np.cos(np.pi / (N + 1))
            sin_b = np.sin(np.pi * (b + 1) / (N + 1))
            obc = np.sqrt(2.0 / (N + 1)) * sin_b
            print(f"  {N:<3} {b:<3} {cls:<10} {two_cos:<16.4f} {sin_b:<22.4f} {obc:<20.4f}")
        print()

    print()
    print("=" * 90)
    print("Step 3b: per-bond g_eff(b) via M_h_total_eff vs per-bond V_b_eff")
    print("=" * 90)
    print()
    # In the 4-mode picture, Q_EP = 2 / g_eff where g_eff = M_h_total_eff[u_0, v_0]
    # (the SVD-block off-diagonal entry). For per-bond analysis the relevant quantity is
    # V_b_eff[u_0, v_0]. Let's compute these and look at Q_peak / Q_EP per bond.
    print(f"  {'N':<3} {'sigma_0(global)':<18} {'Q_EP_global':<14}")
    for N in (5, 6, 7, 8):
        states_p, states_q, flat, Mtot, M_h_per_bond, M_h_total, D, probe_full, S_kernel \
            = build_full_block_operators(N, gamma0)
        B, c1, c3, u0, v0 = build_4mode_basis(states_p, states_q, flat, Mtot, M_h_total)
        M_h_total_eff = project_4mode(B, M_h_total)
        sigma0 = abs(M_h_total_eff[2, 3])
        Q_EP_global = 2.0 / sigma0
        print(f"  {N:<3} {sigma0:<18.6f} {Q_EP_global:<14.6f}")
    print()
    print(f"  {'N':<3} {'b':<3} {'class':<10} {'|V_b[u0,v0]|':<14} {'Q_EP(b)':<10} "
          f"{'Q_peak/Q_EP(b)':<18}")
    for N in (5, 6, 7, 8):
        states_p, states_q, flat, Mtot, M_h_per_bond, M_h_total, D, probe_full, S_kernel \
            = build_full_block_operators(N, gamma0)
        B, c1, c3, u0, v0 = build_4mode_basis(states_p, states_q, flat, Mtot, M_h_total)
        for b in range(N - 1):
            V_b_eff = project_4mode(B, M_h_per_bond[b])
            mag = abs(V_b_eff[2, 3])
            Q_EP_b = 2.0 / mag if mag > 1e-10 else float("inf")
            cls = "Endpoint" if (b == 0 or b == N - 2) else "Interior"
            Q_peak = EMPIRICAL[N]["Q_end"] if cls == "Endpoint" else EMPIRICAL[N]["Q_int"]
            ratio = Q_peak / Q_EP_b
            print(f"  {N:<3} {b:<3} {cls:<10} {mag:<14.4f} {Q_EP_b:<10.4f} {ratio:<18.4f}")
        print()
    print()
    print("=" * 90)
    print("Step 3c: K-driving eigenvalue and its polarity at Q_peak (per bond class)")
    print("=" * 90)
    print()
    print(f"  {'N':<3} {'class':<10} {'Q_peak':<8} {'Re(lam_K+)':<14} {'Im(lam_K+)':<14} "
          f"{'Re(lam_K-)':<14} {'Im(lam_K-)':<14}")
    for N in (5, 6, 7, 8):
        states_p, states_q, flat, Mtot, M_h_per_bond, M_h_total, D, probe_full, S_kernel \
            = build_full_block_operators(N, gamma0)
        B, c1, c3, u0, v0 = build_4mode_basis(states_p, states_q, flat, Mtot, M_h_total)
        D_eff = project_4mode(B, D)
        M_h_total_eff = project_4mode(B, M_h_total)
        probe_eff = project_vec_4mode(B, probe_full)
        for cls, Q_peak_key in (("Endpoint", "Q_end"), ("Interior", "Q_int")):
            Q_peak = EMPIRICAL[N][Q_peak_key]
            L_eff_at_Q = D_eff + Q_peak * gamma0 * M_h_total_eff
            evals, R, pair_indices, overlaps = k_driving_subspace_at_q(L_eff_at_Q, probe_eff)
            lam_plus = evals[pair_indices[0]]
            lam_minus = evals[pair_indices[1]]
            print(f"  {N:<3} {cls:<10} {Q_peak:<8.4f} {lam_plus.real:<+14.6f} {lam_plus.imag:<+14.6f} "
                  f"{lam_minus.real:<+14.6f} {lam_minus.imag:<+14.6f}")
        print()

    print()
    print("=" * 90)
    print("Step 3d: K_b shape per bond at t_peak; Q_peak from LIVE Duhamel scan vs empirical")
    print("=" * 90)
    print()
    # Per-bond K_b(Q, t_peak) Duhamel, with the actual M_h_total_eff dynamics + V_b in dL/dJ_b.
    Q_grid = np.linspace(0.20, 4.00, 153)
    print(f"  {'N':<3} {'b':<3} {'class':<10} {'Q_peak (sim)':<15} {'Q_peak (emp)':<15} "
          f"{'r=Q_peak(sim)-2':<18}")
    for N in (5, 6, 7, 8):
        states_p, states_q, flat, Mtot, M_h_per_bond, M_h_total, D, probe_full, S_kernel \
            = build_full_block_operators(N, gamma0)
        B, c1, c3, u0, v0 = build_4mode_basis(states_p, states_q, flat, Mtot, M_h_total)
        D_eff = project_4mode(B, D)
        M_h_total_eff = project_4mode(B, M_h_total)
        probe_eff = project_vec_4mode(B, probe_full)
        S_kernel_eff = project_4mode(B, S_kernel)
        for b in range(N - 1):
            cls = "Endpoint" if (b == 0 or b == N - 2) else "Interior"
            V_b_eff = project_4mode(B, M_h_per_bond[b])
            # K_b(Q, t_peak) via Duhamel
            K_curve = []
            for Q in Q_grid:
                L_eff_at_Q = D_eff + Q * gamma0 * M_h_total_eff
                evals, R = eig(L_eff_at_Q)
                Rinv = inv(R)
                expLam = np.exp(evals * t_peak)
                rho_t = R @ (expLam * (Rinv @ probe_eff))
                I_mat = np.zeros((4, 4), dtype=complex)
                for j in range(4):
                    for k in range(4):
                        diff = evals[k] - evals[j]
                        if abs(diff) < 1e-10:
                            I_mat[j, k] = t_peak * expLam[j]
                        else:
                            I_mat[j, k] = (expLam[k] - expLam[j]) / diff
                X = Rinv @ V_b_eff @ R
                c0 = Rinv @ probe_eff
                fbC0 = np.zeros(4, dtype=complex)
                for r in range(4):
                    for c in range(4):
                        fbC0[r] += X[r, c] * I_mat[r, c] * c0[c]
                drho = R @ fbC0
                sDrho = S_kernel_eff @ drho
                K = 2.0 * np.vdot(rho_t, sDrho).real
                K_curve.append(K)
            K_arr = np.array(K_curve)
            i_max = int(np.argmax(np.abs(K_arr)))
            if 0 < i_max < len(Q_grid) - 1:
                y0, y1, y2 = abs(K_arr[i_max-1]), abs(K_arr[i_max]), abs(K_arr[i_max+1])
                denom = y0 - 2*y1 + y2
                delta = 0.5 * (y0 - y2) / denom if abs(denom) > 1e-15 else 0.0
                dQ = Q_grid[1] - Q_grid[0]
                Q_peak_sim = Q_grid[i_max] + delta * dQ
            else:
                Q_peak_sim = Q_grid[i_max]
            Q_peak_emp = EMPIRICAL[N]["Q_end"] if cls == "Endpoint" else EMPIRICAL[N]["Q_int"]
            r_sim = Q_peak_sim - 2.0
            print(f"  {N:<3} {b:<3} {cls:<10} {Q_peak_sim:<15.4f} {Q_peak_emp:<15.4f} {r_sim:<+18.4f}")
        print()

    print()
    print("=" * 90)
    print("Step 3e: Tier1 closed-form composition r_Q(N,b) = 2*(BareDoubledPtfXPeak/g_eff(b) - 1)")
    print("=" * 90)
    print()
    print("Two proven Tier1 ingredients:")
    print("  - BareDoubledPtfXPeak = 2.196910 (Tier1Derived universal constant in C2HwhmRatio)")
    print("  - Q_EP(N,b) = 2 / g_eff(N,b) (F86 Statement 1, Tier1Derived)")
    print()
    print("Composition: Q_peak(N,b) = BareDoubledPtfXPeak * Q_EP(N,b)")
    print("           = 2.196910 / g_eff(N,b) * 2 = 4.393820 / g_eff(N,b)")
    print("Therefore: r_Q(N,b) = Q_peak(N,b) - 2 = 4.393820 / g_eff(N,b) - 2")
    print()
    print("Solve for g_eff(N,b) from empirical Q_peak:")
    BARE_X_PEAK = 2.196910
    print(f"  {'N':<3} {'class':<10} {'Q_peak (emp)':<14} {'g_eff(b) inferred':<22} "
          f"{'r_Q (formula)':<16} {'r_Q (emp)':<12} {'|delta|':<10}")
    for N in (5, 6, 7, 8):
        for cls, Q_peak_key in (("Endpoint", "Q_end"), ("Interior", "Q_int")):
            Q_peak_emp = EMPIRICAL[N][Q_peak_key]
            g_eff = 2.0 * BARE_X_PEAK / Q_peak_emp
            r_Q_formula = BARE_X_PEAK * (2.0 / g_eff) - 2.0
            r_Q_emp = Q_peak_emp - 2.0
            diff = abs(r_Q_formula - r_Q_emp)
            print(f"  {N:<3} {cls:<10} {Q_peak_emp:<14.4f} {g_eff:<22.6f} "
                  f"{r_Q_formula:<+16.4f} {r_Q_emp:<+12.4f} {diff:<10.6f}")
        print()
    print("Note: this composition is mathematically identity (we solved for g_eff(N,b) from")
    print("Q_peak, so Q_peak(N,b) = BareDoubledPtfXPeak * Q_EP(N,b) reproduces Q_peak exactly).")
    print("The Tier1Candidate piece is that the COMPOSITION is structural; what is OPEN is the")
    print("closed-form g_eff(N,b) per bond class.")
    print()
    print("Empirical g_eff(N,b) values:")
    print(f"  {'N':<3} {'g_eff_E':<12} {'g_eff_I':<12} {'ratio E/I':<12}")
    for N in (5, 6, 7, 8):
        g_E = 2.0 * BARE_X_PEAK / EMPIRICAL[N]["Q_end"]
        g_I = 2.0 * BARE_X_PEAK / EMPIRICAL[N]["Q_int"]
        print(f"  {N:<3} {g_E:<12.6f} {g_I:<12.6f} {g_E/g_I:<12.6f}")
    print()
    print("Asymptotic behavior of g_eff:")
    print("  g_eff_E -> 4.39382 / 2.5 = 1.757")
    print("  g_eff_I -> 4.39382 / 1.6 = 2.746")
    print("  Both classes have N-stable g_eff;bond-class signature is in g_eff(class), not")
    print("  in N-scaling. This matches the FRAGILE_BRIDGE / EP-rotation hypothesis.")
    print()

    print()
    print("=" * 90)
    print("Step 3f: harmonic-mean closure check on g_eff(N,b)")
    print("=" * 90)
    print()
    print("If 1/g_E + 1/g_I = const, then Q_E + Q_I = 2·BareDoubledPtfXPeak·(1/g_E + 1/g_I)")
    print("            = 4.394·(1/g_E + 1/g_I)")
    print()
    print(f"  {'N':<3} {'1/g_E':<10} {'1/g_I':<10} {'1/g_E+1/g_I':<14} "
          f"{'Q_E+Q_I':<10} {'(Q_E+Q_I)/4.394':<18}")
    for N in (5, 6, 7, 8):
        Q_E = EMPIRICAL[N]["Q_end"]
        Q_I = EMPIRICAL[N]["Q_int"]
        g_E = 2.0 * BARE_X_PEAK / Q_E
        g_I = 2.0 * BARE_X_PEAK / Q_I
        recip_sum = 1.0/g_E + 1.0/g_I
        sum_Q = Q_E + Q_I
        print(f"  {N:<3} {1.0/g_E:<10.4f} {1.0/g_I:<10.4f} {recip_sum:<14.6f} "
              f"{sum_Q:<10.4f} {sum_Q/(2*BARE_X_PEAK):<18.6f}")
    print()
    print("Observation: 1/g_E + 1/g_I appears to converge to ~0.937 across N=6..8 (a")
    print("possible asymptotic structural constant; at N=5 the deviation is ~0.03 which")
    print("may be a finite-size correction).")
    print()
    print("Test against pi-related candidates:")
    candidates = {
        "1 - 1/sqrt(N+1)": [1 - 1/np.sqrt(N+1) for N in (5,6,7,8)],
        "2/pi": [2/np.pi]*4,
        "1 - 1/N": [1 - 1/N for N in (5,6,7,8)],
        "(N+1)/(N+3)": [(N+1)/(N+3) for N in (5,6,7,8)],
    }
    print(f"  {'N':<3} {'empirical':<12} " + "".join(f"{name:<22}" for name in candidates))
    for N in (5, 6, 7, 8):
        Q_E = EMPIRICAL[N]["Q_end"]
        Q_I = EMPIRICAL[N]["Q_int"]
        g_E = 2.0 * BARE_X_PEAK / Q_E
        g_I = 2.0 * BARE_X_PEAK / Q_I
        recip_sum = 1.0/g_E + 1.0/g_I
        i = N - 5
        cand_strs = [f"{candidates[name][i]:<22.6f}" for name in candidates]
        print(f"  {N:<3} {recip_sum:<12.6f} " + "".join(cand_strs))
    print()

    print()
    print("=" * 90)
    print("Step 4: hybrid candidate;1/2 inheritance + perturbation")
    print("=" * 90)
    print()
    print("Note from PolarityInheritanceLink:")
    print("  r_H ~ 0.5 across all bonds (close to HalfAsStructuralFixedPoint)")
    print("  r_Q has bond-class-distinct values around +-0.5")
    print("So r_H reads the FIXED POINT (the 'unsigned 1/2 baseline of the polarity layer')")
    print("and r_Q reads the SIGNED component (the +-0.5 polarity content).")
    print()
    print("Candidate: at the ratio level, if Q_peak / Q_EP ~ 1.x, then the bond-class")
    print("structure is in g_eff(b) -> Q_EP(b) = 2/g_eff(b). Looking at the empirical")
    print("r_Q split:")
    for N in (5, 6, 7, 8):
        rQE = EMPIRICAL[N]["Q_end"] - 2.0
        rQI = EMPIRICAL[N]["Q_int"] - 2.0
        sumR = rQE + rQI
        gap = rQE - rQI
        print(f"  N = {N}: r_Q_E = {rQE:+.4f}, r_Q_I = {rQI:+.4f}, "
              f"sum = {sumR:+.4f}, gap = {gap:+.4f}")
    print()


if __name__ == "__main__":
    main()
