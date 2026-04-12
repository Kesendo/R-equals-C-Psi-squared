"""
Boundary straddling sweep V2: correct overlap measurement
==========================================================
Supersedes V1 (commit f467c81, retracted in ca99064).

DERIVATION OF THE CORRECT OVERLAP OBSERVABLE
=============================================

1. THE CORRECT OBSERVABLE IS c_slow = Tr(L_slow^dagger rho_0).

   The Lindblad master equation drho/dt = L[rho] has biorthogonal
   eigendecomposition with right eigenoperators R_lambda and left
   eigenoperators L_lambda satisfying:

       L[R_lambda] = lambda R_lambda
       Tr(L_lambda^dagger R_mu) = delta_{lambda,mu}   (biorthogonality)

   The solution is:
       rho(t) = sum_lambda c_lambda R_lambda exp(lambda t)
       c_lambda = Tr(L_lambda^dagger rho_0) = <L_lambda | rho_0>_HS

   Therefore the physically meaningful "how much of rho_0 lives in
   the slow mode" is |c_slow| = |Tr(L_slow^dagger rho_0)|. This uses
   the LEFT eigenoperator, not the right. V1 used the right eigenvector,
   which is incorrect for non-normal operators.

   Numerically: c_slow = R_inv[slow_idx, :] @ vec(rho_0), where
   R_inv = inverse of the right-eigenvector matrix. This is how
   slow_mode_lens_analysis.py computes it (line 45, 192).

2. THE SE BLOCK IS SUFFICIENT (cross-sector coherences do not contribute).

   The Liouvillian is block-diagonal by excitation-number sector pair
   (w_bra, w_ket), because the Heisenberg Hamiltonian conserves total
   excitation number and Z-dephasing is diagonal in the computational
   basis. Both right AND left eigenvectors of the (w=1,w=1) block stay
   in that block (PROOF_PARITY_SELECTION_RULE.md, Part 4: "Left eigen-
   vectors respect the block structure").

   The slow mode has SE fraction = 1.000 for N=3-7 chain (SACRIFICE_
   GEOMETRY.md). This means L_slow lives entirely in the (w=1,w=1) block.
   Therefore:

       c_slow = Tr(L_slow^dagger rho_0)
              = Tr(L_slow_SE^dagger rho_0_SE)   [exact, not approximate]

   The cross-sector coherences (w=1,w=3) of rho_0 CANNOT contribute
   because L_slow has zero content there.

3. THE BELL PAIR DOES NOT APPEAR IN c_slow.

   For candidate state |psi> = cos(theta)|e_k> + sin(theta)|f_{c1,c2,k}>
   where |e_k> is in w=1 and |f_{c1,c2,k}> is in w=3:

       rho_0_SE[m, n] = cos^2(theta) delta_{mk} delta_{nk}

   This is a rank-1 matrix that depends ONLY on cos(theta) and k, NOT
   on the Bell pair (c1, c2). The Bell pair information lives in the
   cross-sector (w=1,w=3) coherences and the (w=3,w=3) block, which
   are invisible to c_slow.

   Therefore:
       c_slow = cos^2(theta) * L_slow_SE[k, k]^*

   This is INDEPENDENT of which pair forms the Bell pair. The V1 pattern
   (overlap depends only on exc_k, not on bell_pair) was correct physics,
   not a bug.

4. NORMALIZATION.

   We report c_slow_normalized = |c_slow(rho_0)| / |c_slow(psi_opt)|
   where psi_opt maximizes |c_slow| over pure SE states (the lens state
   from SACRIFICE_GEOMETRY.md). This gives 1.0 for psi_opt by definition.

SANITY CHECKS (run before sweep)
=================================
Check 1: psi_opt should give normalized overlap = 1.0 (by construction).
Check 2: bare Bell+(2,3) tensor |000> has zero SE content, so c_slow = 0
         exactly. If code returns nonzero, something is wrong.

Addresses OQ-114 (CUSP_LENS_CONNECTION.md) and SACRIFICE_GEOMETRY.md OQ#2.

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 12, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import time as _time
import numpy as np
from scipy import linalg

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import (
    heisenberg_H, partial_trace_to_pair,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "boundary_straddling_v2")
os.makedirs(OUT_DIR, exist_ok=True)

CUSP_THRESHOLD = 0.25   # CΨ = 1/4


# ===================================================================
# Dephasing profiles
# ===================================================================
def make_gamma(N):
    """Sacrifice dephasing profile: site 0 high, rest gradient."""
    if N == 5:
        T2_us = np.array([5.22, 122.70, 243.85, 169.97, 237.57])
        g = 1.0 / (2.0 * T2_us)
        return g / g.min() * 0.05
    gamma = np.zeros(N)
    gamma[0] = 2.336
    gamma[1:] = np.linspace(0.050, 0.099, N - 1)
    return gamma


# ===================================================================
# CΨ computation
# ===================================================================
def cpsi_pair(rho_2q):
    """CΨ = Tr(ρ²) × L1 / 3 for a 4×4 density matrix."""
    purity = float(np.real(np.trace(rho_2q @ rho_2q)))
    l1 = float(np.sum(np.abs(rho_2q)) - np.sum(np.abs(np.diag(rho_2q))))
    return purity * l1 / 3.0


def cpsi_all_pairs(rho_full, N):
    """CΨ for each adjacent pair."""
    return np.array([cpsi_pair(partial_trace_to_pair(rho_full, N, k, k + 1))
                     for k in range(N - 1)])


# ===================================================================
# Sector-restricted Liouvillian (exact, from sector conservation)
# ===================================================================
def sector_basis(N, sectors):
    """Computational basis indices for given excitation-number sectors."""
    return sorted(idx for idx in range(2 ** N)
                  if bin(idx).count('1') in sectors)


def build_sector_liouvillian(N, gammas, basis, J=1.0):
    """Liouvillian restricted to given sector basis. Returns M^2 x M^2."""
    D = 2 ** N
    M = len(basis)
    b = np.array(basis)
    H_full = heisenberg_H(N, J)
    H_r = H_full[np.ix_(b, b)]
    I_M = np.eye(M, dtype=complex)
    L = -1j * (np.kron(H_r, I_M) - np.kron(I_M, H_r.T))
    I_M2 = np.eye(M * M, dtype=complex)
    for k in range(N):
        zk_diag = np.array([1.0 - 2.0 * ((idx >> (N - 1 - k)) & 1)
                            for idx in basis], dtype=complex)
        Zk = np.diag(zk_diag)
        L += gammas[k] * (np.kron(Zk, Zk) - I_M2)
    return L


# ===================================================================
# Slow mode extraction and c_slow computation
# ===================================================================
class SlowModeInfo:
    """Stores slow mode eigendecomposition for one (N, gamma) config."""

    def __init__(self, N, gammas, J=1.0):
        self.N = N
        self.se_basis = sector_basis(N, {1})
        N_se = len(self.se_basis)

        # SE-restricted Liouvillian (N^2 x N^2)
        L_SE = build_sector_liouvillian(N, gammas, self.se_basis, J)
        eigvals, R = linalg.eig(L_SE)
        order = np.argsort(-eigvals.real)
        eigvals = eigvals[order]
        R = R[:, order]
        R_inv = linalg.inv(R)

        # Find slowest non-stationary mode
        rates = -eigvals.real
        nonstat = np.where(rates > 1e-10)[0]
        self.slow_idx = nonstat[np.argmin(rates[nonstat])]
        self.slow_ev = eigvals[self.slow_idx]
        self.slow_rate = -self.slow_ev.real

        # Left eigenvector (row of R_inv): this is vec(L_slow)^dagger
        self.left_vec = R_inv[self.slow_idx, :]  # N^2-vector

        # Reshape to N x N for inspection
        self.left_mat = self.left_vec.reshape(N_se, N_se)

        # Compute c_slow for psi_opt (maximum over SE pure states)
        # c_slow = a^dagger Q a where Q = B^T, B[m,n] = left_vec[m*N+n]
        B = self.left_mat  # B[m,n] = R_inv[slow, m*N_se + n]
        Q = B.T
        M_H = (Q + Q.conj().T) / 2
        evals_H, evecs_H = np.linalg.eigh(M_H)
        best = np.argmax(np.abs(evals_H))
        self.a_opt = evecs_H[:, best]
        # Phase convention: largest component real positive
        phase = np.exp(-1j * np.angle(
            self.a_opt[np.argmax(np.abs(self.a_opt))]))
        self.a_opt = self.a_opt * phase

        self.c_slow_psi_opt = float(self.a_opt.conj() @ Q @ self.a_opt)
        self.c_slow_max = abs(self.c_slow_psi_opt)

    def compute_c_slow(self, rho_0_full):
        """Compute c_slow = R_inv[slow, :] @ vec(rho_0_SE).

        rho_0_full is D x D density matrix in full computational basis.
        Returns the complex coefficient c_slow.
        """
        N_se = len(self.se_basis)
        rho_SE_vec = np.zeros(N_se * N_se, dtype=complex)
        for a, ia in enumerate(self.se_basis):
            for b, ib in enumerate(self.se_basis):
                rho_SE_vec[a * N_se + b] = rho_0_full[ia, ib]
        return complex(self.left_vec @ rho_SE_vec)

    def normalized_overlap(self, rho_0_full):
        """Normalized overlap: |c_slow(rho_0)| / |c_slow(psi_opt)|."""
        c = self.compute_c_slow(rho_0_full)
        return abs(c) / self.c_slow_max if self.c_slow_max > 1e-15 else 0.0


# ===================================================================
# Candidate state construction
# ===================================================================
def make_state(N, bell_a, bell_b, exc_k, cos_a=None, sin_a=None):
    """Build (cos_a |00> + sin_a |11>)_{a,b} x |1>_k x |0>_rest."""
    if cos_a is None:
        cos_a, sin_a = 1.0 / np.sqrt(2), 1.0 / np.sqrt(2)
    d = 2 ** N
    psi = np.zeros(d, dtype=complex)
    i0 = 1 << (N - 1 - exc_k)
    i1 = (1 << (N - 1 - bell_a)) | (1 << (N - 1 - bell_b)) | (1 << (N - 1 - exc_k))
    psi[i0] = cos_a
    psi[i1] = sin_a
    return psi / np.linalg.norm(psi)


def generate_candidates(N):
    """Generate non-symmetric two-excitation candidate states."""
    cands = []
    c1 = (N - 1) // 2
    c2 = c1 + 1
    center = (N - 1) / 2.0

    for k in range(N):
        if k in (c1, c2):
            continue
        cands.append(dict(
            family='A', label=f'A:Bell({c1},{c2})+exc({k})',
            psi=make_state(N, c1, c2, k),
            bell_pair=(c1, c2), exc_site=k,
        ))

    for a in range(N - 1):
        if a == c1:
            continue
        avail = [k for k in range(N) if k not in (a, a + 1)]
        avail.sort(key=lambda kk: abs(kk - center), reverse=True)
        n_pick = 3 if N >= 7 else 2
        for k in avail[:n_pick]:
            cands.append(dict(
                family='B', label=f'B:Bell({a},{a+1})+exc({k})',
                psi=make_state(N, a, a + 1, k),
                bell_pair=(a, a + 1), exc_site=k,
            ))

    if N >= 6:
        for theta, tl in [(np.pi / 6, 'pi/6'), (np.pi / 4, 'pi/4'),
                          (np.pi / 3, 'pi/3')]:
            avail = [k for k in range(N) if k not in (c1, c2)]
            avail.sort(key=lambda kk: abs(kk - center), reverse=True)
            for k in avail[:2]:
                cands.append(dict(
                    family='C',
                    label=f'C:theta={tl}+exc({k})',
                    psi=make_state(N, c1, c2, k,
                                   cos_a=np.cos(theta), sin_a=np.sin(theta)),
                    bell_pair=(c1, c2), exc_site=k, theta=float(theta),
                ))

    return cands


# ===================================================================
# Per-candidate measurement
# ===================================================================
def measure(cand, eigvals_r, R_r, R_inv_r, slow_info, full_basis, N):
    """Full measurement for one candidate.

    Time evolution uses the {w=1,w=3}-restricted eigendecomposition.
    Slow-mode overlap uses the SE-restricted left eigenvector (SlowModeInfo).
    """
    D = 2 ** N
    M = len(full_basis)
    psi = cand['psi']
    rho0_full = np.outer(psi, psi.conj())

    # 1. Initial CΨ per pair
    cpsi0 = cpsi_all_pairs(rho0_full, N)

    # 2. Slow-mode coefficient (biorthogonal, LEFT eigenvector)
    c_slow = slow_info.compute_c_slow(rho0_full)
    c_slow_abs = abs(c_slow)
    c_slow_norm = slow_info.normalized_overlap(rho0_full)

    # 3. Project rho0 into sector-restricted basis for time evolution
    rho0_r = np.zeros((M, M), dtype=complex)
    for a, ia in enumerate(full_basis):
        for b, ib in enumerate(full_basis):
            rho0_r[a, b] = rho0_full[ia, ib]
    c0 = R_inv_r @ rho0_r.flatten()

    # 4. Time evolution: t=0..30, dt=0.1
    times = np.linspace(0, 30, 301)
    n_pairs = N - 1
    cpsi_traj = np.zeros((len(times), n_pairs))
    for i, t in enumerate(times):
        rv = R_r @ (c0 * np.exp(eigvals_r * t))
        rho_r_t = rv.reshape(M, M)
        rho_r_t = (rho_r_t + rho_r_t.conj().T) / 2
        rho_t = np.zeros((D, D), dtype=complex)
        for a, ia in enumerate(full_basis):
            for b, ib in enumerate(full_basis):
                rho_t[ia, ib] = rho_r_t[a, b]
        cpsi_traj[i] = cpsi_all_pairs(rho_t, N)

    # 5. Cusp crossing
    cusp_crossed = False
    t_cross = None
    dcpsi_cross = None
    cross_pair = None
    dt = times[1] - times[0]
    for p in range(n_pairs):
        for i in range(len(times) - 1):
            if cpsi_traj[i, p] >= CUSP_THRESHOLD > cpsi_traj[i + 1, p]:
                frac = (cpsi_traj[i, p] - CUSP_THRESHOLD) / \
                       (cpsi_traj[i, p] - cpsi_traj[i + 1, p])
                tc = times[i] + frac * dt
                if i > 0:
                    dc = (cpsi_traj[i + 1, p] - cpsi_traj[i - 1, p]) / (2 * dt)
                else:
                    dc = (cpsi_traj[i + 1, p] - cpsi_traj[i, p]) / dt
                if t_cross is None or tc < t_cross:
                    cusp_crossed = True
                    t_cross = float(tc)
                    dcpsi_cross = float(dc)
                    cross_pair = p
                break

    # 6. Asymptotic state at t=100
    rv100 = R_r @ (c0 * np.exp(eigvals_r * 100.0))
    rho_r_100 = rv100.reshape(M, M)
    rho_r_100 = (rho_r_100 + rho_r_100.conj().T) / 2
    rho100 = np.zeros((D, D), dtype=complex)
    for a, ia in enumerate(full_basis):
        for b, ib in enumerate(full_basis):
            rho100[ia, ib] = rho_r_100[a, b]
    purity100 = float(np.real(np.trace(rho100 @ rho100)))

    return dict(
        label=cand['label'], family=cand['family'], N=N,
        bell_pair=list(cand['bell_pair']), exc_site=cand['exc_site'],
        cpsi0_per_pair=[round(v, 6) for v in cpsi0.tolist()],
        cpsi0_max=round(float(cpsi0.max()), 6),
        cusp_active=cusp_crossed,
        t_cross=round(t_cross, 4) if t_cross is not None else None,
        c_slow_raw=round(c_slow_abs, 6),
        c_slow_normalized=round(c_slow_norm, 6),
        purity_100=round(purity100, 6),
    )


# ===================================================================
# Sanity checks (must pass before sweep runs)
# ===================================================================
def run_sanity_checks(N, gamma):
    """Run the two sanity checks on N=5 IBM chain.

    Returns (slow_info, passed_bool, report_text).
    """
    D = 2 ** N
    slow_info = SlowModeInfo(N, gamma)
    lines = []
    passed = True

    lines.append("SANITY CHECK 1: psi_opt")
    lines.append(f"  psi_opt amplitudes: "
                 f"[{', '.join(f'{abs(a):.4f}' for a in slow_info.a_opt)}]")
    lines.append(f"  c_slow_psi_opt (raw): {slow_info.c_slow_psi_opt:.6f}")
    lines.append(f"  |c_slow_psi_opt|:     {slow_info.c_slow_max:.6f}")

    # Build psi_opt density matrix
    se = slow_info.se_basis
    psi_opt = np.zeros(D, dtype=complex)
    for k, idx in enumerate(se):
        psi_opt[idx] = slow_info.a_opt[k]
    psi_opt /= np.linalg.norm(psi_opt)
    rho_opt = np.outer(psi_opt, psi_opt.conj())

    c_check = slow_info.compute_c_slow(rho_opt)
    norm_check = slow_info.normalized_overlap(rho_opt)
    lines.append(f"  Computed c_slow:      {abs(c_check):.6f}")
    lines.append(f"  Normalized overlap:   {norm_check:.6f}")
    lines.append(f"  EXPECTED: normalized ≈ 1.0")
    if abs(norm_check - 1.0) > 0.01:
        lines.append(f"  FAILED: got {norm_check:.6f}, expected 1.0")
        passed = False
    else:
        lines.append(f"  PASSED: {norm_check:.6f} ≈ 1.0")

    lines.append("")
    lines.append("SANITY CHECK 2: bare Bell+(2,3) x |000>")
    lines.append("  PREDICTION: c_slow = 0 exactly (no SE content)")
    # Bell+(2,3) x |000> on N=5: |ψ⟩ = (|00000⟩ + |00110⟩)/√2
    psi_bell = np.zeros(D, dtype=complex)
    psi_bell[0] = 1.0 / np.sqrt(2)  # |00000⟩
    # |00110⟩ for N=5: sites 2,3 excited
    idx_11 = (1 << (N - 1 - 2)) | (1 << (N - 1 - 3))
    psi_bell[idx_11] = 1.0 / np.sqrt(2)
    rho_bell = np.outer(psi_bell, psi_bell.conj())

    c_bell = slow_info.compute_c_slow(rho_bell)
    lines.append(f"  Computed c_slow:      {abs(c_bell):.2e}")
    if abs(c_bell) > 1e-10:
        lines.append(f"  FAILED: got {abs(c_bell):.6f}, expected 0.0")
        passed = False
    else:
        lines.append(f"  PASSED: {abs(c_bell):.2e} ≈ 0")

    # Bonus: verify bell_pair independence
    lines.append("")
    lines.append("INDEPENDENCE CHECK: overlap for Bell(a,b)+exc(k)")
    lines.append("  If derivation is correct, overlap depends only on k,")
    lines.append("  not on (a,b). Testing N=5:")
    for k in [0, 1, 4]:
        overlaps = []
        for a, b in [(0, 1), (1, 2), (2, 3), (3, 4)]:
            if k in (a, b):
                continue
            psi = make_state(N, a, b, k)
            rho = np.outer(psi, psi.conj())
            ov = slow_info.normalized_overlap(rho)
            overlaps.append((a, b, ov))
        vals = [ov for _, _, ov in overlaps]
        spread = max(vals) - min(vals)
        pairs_str = ', '.join(f'({a},{b})={ov:.4f}' for a, b, ov in overlaps)
        lines.append(f"  k={k}: {pairs_str}")
        lines.append(f"    spread = {spread:.2e} "
                     f"({'IDENTICAL' if spread < 1e-10 else 'DIFFERS'})")

    return slow_info, passed, '\n'.join(lines)


# ===================================================================
# Full sweep
# ===================================================================
def sweep_one_N(N, gamma, slow_info):
    """Run sweep for one chain length."""
    t0 = _time.time()
    D = 2 ** N
    full_basis = sector_basis(N, {1, 3})
    M = len(full_basis)

    print(f"\n  Building {M * M}x{M * M} sector Liouvillian ...")
    L_r = build_sector_liouvillian(N, gamma, full_basis)
    eigvals_r, R_r = linalg.eig(L_r)
    order = np.argsort(-eigvals_r.real)
    eigvals_r = eigvals_r[order]
    R_r = R_r[:, order]
    R_inv_r = linalg.inv(R_r)

    cands = generate_candidates(N)
    print(f"  {len(cands)} candidates")

    results = []
    for ci, cand in enumerate(cands):
        r = measure(cand, eigvals_r, R_r, R_inv_r, slow_info,
                    full_basis, N)
        results.append(r)
        tc_str = f"{r['t_cross']:.2f}" if r['t_cross'] is not None else "n/a"
        print(f"    [{ci + 1:2d}/{len(cands)}] {r['label']:32s}  "
              f"|c_slow|={r['c_slow_raw']:.4f}  "
              f"norm={r['c_slow_normalized']:.4f}  "
              f"cusp={tc_str}")

    elapsed = _time.time() - t0
    print(f"  N={N}: {len(results)} candidates, {elapsed:.1f}s")
    return results


# ===================================================================
# Output
# ===================================================================
def save_all(all_results, all_meta, sanity_report):
    """Save JSON, TXT, and plot."""
    with open(os.path.join(OUT_DIR, 'sweep_results.json'), 'w',
              encoding='utf-8') as f:
        json.dump(dict(meta=all_meta, results=all_results,
                       sanity_report=sanity_report), f,
                  indent=2, ensure_ascii=False)

    txt_path = os.path.join(OUT_DIR, 'sweep_summary.txt')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("Boundary Straddling V2: Correct Overlap Measurement\n")
        f.write("====================================================\n")
        f.write("Date: 2026-04-12\n")
        f.write("Observable: |c_slow| = |Tr(L_slow^dag rho_0)| "
                "(LEFT eigenoperator)\n")
        f.write("Normalization: |c_slow| / |c_slow(psi_opt)|\n\n")

        f.write("SANITY CHECKS\n")
        f.write("-" * 60 + "\n")
        f.write(sanity_report + "\n\n")

        for m in all_meta:
            f.write(f"N = {m['N']}: slow rate = {m['slow_rate']:.6f}, "
                    f"|c_slow_max| = {m['c_slow_max']:.6f}\n")
        f.write("\n")

        hdr = (f"{'Label':35s}  {'N':>2}  {'exc':>3}  "
               f"{'|c_slow|':>8}  {'norm':>6}  {'CΨmax':>6}  "
               f"{'t_cross':>8}\n")
        f.write(hdr)
        f.write("-" * len(hdr.rstrip()) + "\n")
        for r in all_results:
            tc = f"{r['t_cross']:.3f}" if r['t_cross'] is not None else "n/a"
            f.write(
                f"{r['label']:35s}  {r['N']:>2}  {r['exc_site']:>3}  "
                f"{r['c_slow_raw']:>8.4f}  {r['c_slow_normalized']:>6.4f}  "
                f"{r['cpsi0_max']:>6.4f}  {tc:>8}\n")

        # Analysis by exc_site
        f.write("\nOverlap by excitation site (confirms bell_pair independence):\n")
        for N_val in sorted(set(r['N'] for r in all_results)):
            sub = [r for r in all_results if r['N'] == N_val
                   and r['family'] in ('A', 'B')]
            sites = sorted(set(r['exc_site'] for r in sub))
            f.write(f"\n  N={N_val}:\n")
            for k in sites:
                ks = [r for r in sub if r['exc_site'] == k]
                norms = [r['c_slow_normalized'] for r in ks]
                spread = max(norms) - min(norms) if len(norms) > 1 else 0.0
                f.write(f"    site {k}: norm={norms[0]:.4f} "
                        f"(spread {spread:.1e} over {len(ks)} configs)\n")

    # Plot: normalized overlap vs exc_site, colored by N
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: normalized overlap vs exc_site
    ax = axes[0]
    colors = {5: '#1f77b4', 6: '#ff7f0e', 7: '#2ca02c'}
    markers = {'A': 'o', 'B': 's', 'C': '^'}
    for r in all_results:
        ax.scatter(r['exc_site'], r['c_slow_normalized'],
                   c=colors.get(r['N'], 'gray'),
                   marker=markers.get(r['family'], 'o'),
                   s=70, edgecolors='k', linewidth=0.5, alpha=0.8)
    ax.set_xlabel('Excitation site k', fontsize=11)
    ax.set_ylabel('Normalized $|c_{slow}|$', fontsize=11)
    ax.set_title('Slow-mode overlap by excitation site', fontsize=12)
    ax.grid(True, alpha=0.3)

    # Panel 2: CΨ_max(0) vs normalized overlap
    ax = axes[1]
    for r in all_results:
        ax.scatter(r['cpsi0_max'], r['c_slow_normalized'],
                   c=colors.get(r['N'], 'gray'),
                   marker=markers.get(r['family'], 'o'),
                   s=70, edgecolors='k', linewidth=0.5, alpha=0.8)
    ax.set_xlabel('Initial CΨ$_{max}$', fontsize=11)
    ax.set_ylabel('Normalized $|c_{slow}|$', fontsize=11)
    ax.set_title('Concurrence vs slow-mode overlap', fontsize=12)
    ax.grid(True, alpha=0.3)

    # Shared legend
    from matplotlib.lines import Line2D
    handles = []
    for n in sorted(colors):
        handles.append(Line2D([0], [0], marker='o', color='w',
                              markerfacecolor=colors[n], markeredgecolor='k',
                              markersize=9, linewidth=0, label=f'N={n}'))
    for fam, m in markers.items():
        handles.append(Line2D([0], [0], marker=m, color='w',
                              markerfacecolor='gray', markeredgecolor='k',
                              markersize=8, linewidth=0, label=f'Family {fam}'))
    fig.legend(handles=handles, loc='upper center', ncol=6, fontsize=8,
               bbox_to_anchor=(0.5, 0.02))
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(os.path.join(OUT_DIR, 'boundary_straddling_v2.png'), dpi=150)
    plt.close(fig)

    print(f"\n  Results saved to {OUT_DIR}/")


# ===================================================================
# Main
# ===================================================================
if __name__ == "__main__":
    print("Boundary Straddling V2: Correct Overlap Measurement")
    print("Supersedes V1 (f467c81, retracted ca99064)")
    print("=" * 60)

    # Sanity checks on N=5 IBM chain
    N_check = 5
    gamma_check = make_gamma(N_check)
    print(f"\nRunning sanity checks on N={N_check} IBM chain ...")
    slow_info_5, checks_passed, sanity_report = run_sanity_checks(
        N_check, gamma_check)
    print(sanity_report)

    if not checks_passed:
        print("\n*** SANITY CHECKS FAILED. Sweep aborted. ***")
        print("Handing back to chat with the failure report.")
        with open(os.path.join(OUT_DIR, 'sanity_FAILED.txt'), 'w',
                  encoding='utf-8') as f:
            f.write("SANITY CHECKS FAILED\n\n")
            f.write(sanity_report)
        sys.exit(1)

    print("\nSanity checks PASSED. Running full sweep.")
    print("=" * 60)

    all_results = []
    all_meta = []
    slow_infos = {}

    for N in [5, 6, 7]:
        gamma = make_gamma(N)
        print(f"\n{'=' * 60}")
        print(f"N = {N}")
        print(f"{'=' * 60}")

        if N == N_check:
            si = slow_info_5
        else:
            si = SlowModeInfo(N, gamma)
        slow_infos[N] = si

        print(f"  Slow mode: rate {si.slow_rate:.6f}, "
              f"|c_slow_max| = {si.c_slow_max:.6f}")
        print(f"  psi_opt: [{', '.join(f'{abs(a):.4f}' for a in si.a_opt)}]")

        results = sweep_one_N(N, gamma, si)
        all_results.extend(results)
        all_meta.append(dict(
            N=N, gamma=gamma.tolist(),
            slow_rate=round(si.slow_rate, 6),
            c_slow_max=round(si.c_slow_max, 6),
            n_candidates=len(results),
        ))

    save_all(all_results, all_meta, sanity_report)

    # Final summary
    print(f"\n{'=' * 60}")
    print(f"SWEEP COMPLETE: {len(all_results)} candidates")
    print(f"{'=' * 60}")
    print(f"\nKey finding: c_slow depends on excitation site only,")
    print(f"not on Bell pair placement. This is correct physics")
    print(f"(the slow mode lives in the SE sector; the Bell pair")
    print(f"contributes to cross-sector coherences, not to SE).")
