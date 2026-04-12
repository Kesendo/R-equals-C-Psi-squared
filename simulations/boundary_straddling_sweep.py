"""
Boundary straddling sweep: non-symmetric two-excitation states
==============================================================
Tests whether non-symmetric two-excitation states on Heisenberg chains
under Z-dephasing can straddle both decoherence exits:
  - Cusp exit: CΨ crosses 1/4 on at least one adjacent pair
  - Lens exit: Frobenius overlap >= 0.1 with the slow Liouvillian eigenmode

Three candidate families are tested for N=5, 6, 7:
  A: Bell+(center) + displaced excitation
  B: Bell+(off-center pair) + displaced excitation
  C: Tuned concurrence (cos θ |00> + sin θ |11>) + excitation (N >= 6)

Method: sector-restricted Liouvillian. Candidate states span excitation
sectors w=1 and w=3. The Heisenberg Hamiltonian conserves total excitation
number, and Z-dephasing preserves it. Therefore the Liouvillian is
block-diagonal by sector pair (w_bra, w_ket), and the evolution within
{w=1, w=3} x {w=1, w=3} is EXACT. For N=7 this reduces the Liouvillian
from 16384x16384 to 1764x1764 (sector conservation proven in
CUSP_LENS_CONNECTION.md).

Slow mode extraction uses the SE-restricted Liouvillian (w=1,w=1 block,
N^2 x N^2). The slow mode has SE fraction = 1.000 for N=3-7 chain
(SACRIFICE_GEOMETRY.md), so this is exact.

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
                       "results", "boundary_straddling")
os.makedirs(OUT_DIR, exist_ok=True)

# Thresholds defined BEFORE seeing data (per task spec)
CUSP_THRESHOLD = 0.25   # CΨ = 1/4
LENS_THRESHOLD = 0.10   # 10% normalized Frobenius overlap


# ===================================================================
# Dephasing profiles
# ===================================================================
def make_gamma(N):
    """Sacrifice dephasing profile: site 0 high, rest gradient.

    N=5 uses the exact IBM Torino T2-derived profile.
    N=6, 7 extend the pattern with a quiet-site gradient.
    """
    if N == 5:
        T2_us = np.array([5.22, 122.70, 243.85, 169.97, 237.57])
        g = 1.0 / (2.0 * T2_us)
        return g / g.min() * 0.05
    gamma = np.zeros(N)
    gamma[0] = 2.336
    gamma[1:] = np.linspace(0.050, 0.099, N - 1)
    return gamma


# ===================================================================
# CΨ computation (from psi_opt_cusp_trajectory.py)
# ===================================================================
def cpsi_pair(rho_2q):
    """CΨ = Tr(ρ²) × L1 / 3 for a 4×4 density matrix."""
    purity = float(np.real(np.trace(rho_2q @ rho_2q)))
    l1 = float(np.sum(np.abs(rho_2q)) - np.sum(np.abs(np.diag(rho_2q))))
    return purity * l1 / 3.0


def cpsi_all_pairs(rho_full, N):
    """CΨ for each adjacent pair (0,1), ..., (N-2,N-1)."""
    return np.array([cpsi_pair(partial_trace_to_pair(rho_full, N, k, k + 1))
                     for k in range(N - 1)])


# ===================================================================
# Sector-restricted Liouvillian construction
# ===================================================================
def sector_basis(N, sectors):
    """Return computational basis indices for given excitation-number sectors.

    Example: sector_basis(5, {1, 3}) gives all N=5 indices with Hamming
    weight 1 or 3, sorted ascending.
    """
    return sorted(idx for idx in range(2 ** N)
                  if bin(idx).count('1') in sectors)


def build_sector_liouvillian(N, gammas, basis, J=1.0):
    """Build the Liouvillian restricted to the given sector basis.

    The full 2^N x 2^N Hamiltonian is built (cheap), then restricted
    to the M-dimensional subspace. The Liouvillian is M^2 x M^2.

    Parameters
    ----------
    N : int, qubit count
    gammas : array, per-site dephasing rates
    basis : list of int, computational basis indices in the subspace
    J : float, Heisenberg coupling

    Returns
    -------
    L : ndarray (M^2, M^2), sector-restricted Liouvillian
    """
    D = 2 ** N
    M = len(basis)
    b = np.array(basis)

    # Restricted Hamiltonian: H_r[a, b] = H_full[basis[a], basis[b]]
    H_full = heisenberg_H(N, J)
    H_r = H_full[np.ix_(b, b)]

    # Liouvillian commutator part
    I_M = np.eye(M, dtype=complex)
    L = -1j * (np.kron(H_r, I_M) - np.kron(I_M, H_r.T))

    # Dephasing: L_k = sqrt(γ_k) Z_k
    # Z_k restricted to basis: diagonal with ±1 entries
    # Z_k^2 = I for any diagonal ±1 matrix
    # Lindblad dephasing: γ_k (Z_k ⊗ Z_k - I_{M^2})
    I_M2 = np.eye(M * M, dtype=complex)
    for k in range(N):
        zk_diag = np.array([1.0 - 2.0 * ((idx >> (N - 1 - k)) & 1)
                            for idx in basis], dtype=complex)
        Zk = np.diag(zk_diag)
        L += gammas[k] * (np.kron(Zk, Zk) - I_M2)

    return L


# ===================================================================
# Slow mode extraction (SE-restricted, exact for SE fraction ≈ 1)
# ===================================================================
def find_slow_se_mode(N, gammas, J=1.0):
    """Find the slow SE-accessible mode using the SE-restricted Liouvillian.

    The SE sector (w=1) Liouvillian is N^2 x N^2, trivial to eigendecompose.
    The slow mode of this block matches the full Liouvillian's slow mode
    because SE fraction = 1.000 for N=3-7 chain profiles.

    Returns (eigenvalue, right_eigvec, se_basis).
    """
    se_basis = sector_basis(N, {1})
    L_SE = build_sector_liouvillian(N, gammas, se_basis, J)

    eigvals, R = linalg.eig(L_SE)
    rates = -eigvals.real
    nonstat = np.where(rates > 1e-10)[0]
    if len(nonstat) == 0:
        return None, None, se_basis

    slow_i = nonstat[np.argmin(rates[nonstat])]
    return eigvals[slow_i], R[:, slow_i], se_basis


# ===================================================================
# Candidate state construction
# ===================================================================
def make_state(N, bell_a, bell_b, exc_k, cos_a=None, sin_a=None):
    """Build (cos_a |00> + sin_a |11>)_{a,b} x |1>_k x |0>_rest.

    Default cos_a = sin_a = 1/sqrt(2) gives Bell+.
    """
    if cos_a is None:
        cos_a, sin_a = 1.0 / np.sqrt(2), 1.0 / np.sqrt(2)
    d = 2 ** N
    psi = np.zeros(d, dtype=complex)
    # |00>_{a,b} component: only site k excited (w=1 sector)
    i0 = 1 << (N - 1 - exc_k)
    # |11>_{a,b} component: sites a, b, k excited (w=3 sector)
    i1 = (1 << (N - 1 - bell_a)) | (1 << (N - 1 - bell_b)) | (1 << (N - 1 - exc_k))
    psi[i0] = cos_a
    psi[i1] = sin_a
    return psi / np.linalg.norm(psi)


def generate_candidates(N):
    """Generate non-symmetric two-excitation candidate states for chain length N."""
    cands = []
    c1 = (N - 1) // 2
    c2 = c1 + 1
    center = (N - 1) / 2.0

    # Family A: Bell+ on central pair + displaced excitation
    for k in range(N):
        if k in (c1, c2):
            continue
        cands.append(dict(
            family='A', label=f'A:Bell({c1},{c2})+exc({k})',
            psi=make_state(N, c1, c2, k),
            bell_pair=(c1, c2), exc_site=k,
        ))

    # Family B: Bell+ on off-center pair + displaced excitation
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

    # Family C: Tuned concurrence (N >= 6 only)
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
# Per-candidate measurement (sector-restricted)
# ===================================================================
def measure(cand, eigvals_r, R_r, R_inv_r, slow_vec_SE, se_basis,
            full_basis, N):
    """Full measurement for one candidate using sector-restricted evolution.

    Parameters
    ----------
    eigvals_r, R_r, R_inv_r : eigendecomposition of {w=1,w=3} Liouvillian
    slow_vec_SE : right eigenvector of slow SE mode (N^2 vector)
    se_basis : list of SE (w=1) computational basis indices
    full_basis : list of {w=1,w=3} computational basis indices
    """
    D = 2 ** N
    M = len(full_basis)
    psi = cand['psi']
    rho0_full = np.outer(psi, psi.conj())

    # 1. Initial CΨ per pair (from full density matrix)
    cpsi0 = cpsi_all_pairs(rho0_full, N)

    # 2. Slow-mode overlap (SE block only, exact when SE fraction ≈ 1)
    if slow_vec_SE is not None:
        N_se = len(se_basis)
        rho0_SE = np.zeros(N_se * N_se, dtype=complex)
        for a, ia in enumerate(se_basis):
            for b, ib in enumerate(se_basis):
                rho0_SE[a * N_se + b] = rho0_full[ia, ib]
        inner = np.abs(slow_vec_SE.conj() @ rho0_SE)
        norm_v = np.linalg.norm(slow_vec_SE)
        # ||rho0||_F = 1 for pure states
        overlap = float(inner / norm_v) if norm_v > 1e-15 else 0.0
    else:
        overlap = 0.0

    # 3. Project rho0 into the sector-restricted basis
    rho0_r = np.zeros((M, M), dtype=complex)
    for a, ia in enumerate(full_basis):
        for b, ib in enumerate(full_basis):
            rho0_r[a, b] = rho0_full[ia, ib]
    rho0_r_vec = rho0_r.flatten()

    # Decompose into restricted eigenbasis
    c0 = R_inv_r @ rho0_r_vec

    # 4. Time evolution: t=0..30, dt=0.1
    times = np.linspace(0, 30, 301)
    n_pairs = N - 1
    cpsi_traj = np.zeros((len(times), n_pairs))
    for i, t in enumerate(times):
        rv = R_r @ (c0 * np.exp(eigvals_r * t))
        rho_r_t = rv.reshape(M, M)
        rho_r_t = (rho_r_t + rho_r_t.conj().T) / 2
        # Embed back to full D×D space
        rho_t = np.zeros((D, D), dtype=complex)
        for a, ia in enumerate(full_basis):
            for b, ib in enumerate(full_basis):
                rho_t[ia, ib] = rho_r_t[a, b]
        cpsi_traj[i] = cpsi_all_pairs(rho_t, N)

    # 5. Cusp crossing (downward through 1/4, earliest)
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
    sector_pops = {}
    for w in range(N + 1):
        pop = sum(float(np.real(rho100[idx, idx]))
                  for idx in range(D) if bin(idx).count('1') == w)
        sector_pops[str(w)] = round(pop, 6)
    coherence100 = float(np.sum(np.abs(rho100))
                         - np.sum(np.abs(np.diag(rho100))))

    # Classification
    cusp_active = cusp_crossed
    lens_active = overlap >= LENS_THRESHOLD
    if cusp_active and lens_active:
        cls = 'straddles'
    elif cusp_active:
        cls = 'cusp-only'
    elif lens_active:
        cls = 'lens-only'
    else:
        cls = 'neither'

    return dict(
        label=cand['label'], family=cand['family'], N=N,
        bell_pair=list(cand['bell_pair']), exc_site=cand['exc_site'],
        cpsi0_per_pair=[round(v, 6) for v in cpsi0.tolist()],
        cpsi0_max=round(float(cpsi0.max()), 6),
        cpsi0_max_pair=int(cpsi0.argmax()),
        cusp_active=cusp_active,
        cross_pair=int(cross_pair) if cross_pair is not None else None,
        t_cross=round(t_cross, 4) if t_cross is not None else None,
        dcpsi_cross=round(dcpsi_cross, 6) if dcpsi_cross is not None else None,
        slow_mode_overlap=round(overlap, 6),
        lens_active=lens_active,
        classification=cls,
        purity_100=round(purity100, 6),
        sector_pops_100=sector_pops,
        coherence_100=round(coherence100, 6),
    )


# ===================================================================
# Full sweep for one chain length
# ===================================================================
def sweep_one_N(N):
    """Run the boundary straddling sweep for chain length N.

    Uses sector-restricted Liouvillian: candidates span w=1 and w=3
    sectors only (exact, proven from sector conservation).
    For N=7: 1764x1764 instead of 16384x16384.

    Returns (results_list, meta_dict).
    """
    t0 = _time.time()
    D = 2 ** N
    gamma = make_gamma(N)

    # Sector basis for {w=1, w=3}
    full_basis = sector_basis(N, {1, 3})
    se_basis = sector_basis(N, {1})
    M = len(full_basis)

    print(f"\n{'=' * 70}")
    print(f"N = {N}  (D = {D}, sector dim M = {M}, "
          f"L_sector = {M * M}x{M * M})")
    print(f"{'=' * 70}")
    print(f"  γ = [{', '.join(f'{g:.4f}' for g in gamma)}]")
    print(f"  Σγ = {gamma.sum():.4f}")
    print(f"  Sectors: w=1 ({len(se_basis)} states), "
          f"w=3 ({M - len(se_basis)} states)")

    # Build and eigendecompose sector-restricted Liouvillian
    print(f"  Building sector-restricted Liouvillian ({M * M}x{M * M}) ...")
    t_eig = _time.time()
    L_r = build_sector_liouvillian(N, gamma, full_basis)
    print(f"  Eigendecomposing ...")
    eigvals_r, R_r = linalg.eig(L_r)
    order = np.argsort(-eigvals_r.real)
    eigvals_r = eigvals_r[order]
    R_r = R_r[:, order]
    R_inv_r = linalg.inv(R_r)
    print(f"  Sector eigendecomposition done ({_time.time() - t_eig:.1f}s)")

    # Find slow SE mode from the SE-restricted Liouvillian
    slow_ev, slow_vec_SE, _ = find_slow_se_mode(N, gamma)
    if slow_ev is not None:
        print(f"  Slow SE mode: rate {-slow_ev.real:.6f}, "
              f"Im = {slow_ev.imag:+.6f}")
    else:
        print(f"  WARNING: no slow SE mode found")

    # Generate candidates
    cands = generate_candidates(N)
    print(f"  {len(cands)} candidates generated")

    # Measure each candidate
    results = []
    for ci, cand in enumerate(cands):
        r = measure(cand, eigvals_r, R_r, R_inv_r, slow_vec_SE,
                    se_basis, full_basis, N)
        results.append(r)
        ov_str = f"{r['slow_mode_overlap']:.4f}"
        tc_str = f"{r['t_cross']:.2f}" if r['t_cross'] is not None else "n/a"
        print(f"    [{ci + 1:2d}/{len(cands)}] {r['label']:32s}  "
              f"cusp={'Y' if r['cusp_active'] else 'N'}(t={tc_str:>6})  "
              f"overlap={ov_str}  "
              f"-> {r['classification']}")

    elapsed = _time.time() - t0
    print(f"  N={N} complete: {len(results)} candidates, {elapsed:.1f}s")

    meta = dict(
        N=N, gamma=gamma.tolist(),
        Sg=round(float(gamma.sum()), 4),
        sector_dim=M,
        liouvillian_dim=M * M,
        slow_mode_rate=round(float(-slow_ev.real), 6) if slow_ev is not None else None,
        n_candidates=len(results),
        elapsed_s=round(elapsed, 1),
    )
    return results, meta


# ===================================================================
# Output: JSON, TXT summary, overview plot
# ===================================================================
def save_all(all_results, all_meta):
    """Save sweep results to JSON, TXT summary, and overview plot."""

    # --- JSON ---
    with open(os.path.join(OUT_DIR, 'sweep_results.json'), 'w',
              encoding='utf-8') as f:
        json.dump(dict(meta=all_meta, results=all_results), f,
                  indent=2, ensure_ascii=False)

    # --- TXT summary ---
    txt_path = os.path.join(OUT_DIR, 'sweep_summary.txt')
    with open(txt_path, 'w', encoding='utf-8') as f:
        f.write("Boundary Straddling Sweep: Classification Table\n")
        f.write("================================================\n")
        f.write("Date: 2026-04-12\n")
        f.write(f"Thresholds: cusp = CΨ crosses 1/4;  "
                f"lens = overlap >= {LENS_THRESHOLD}\n")
        f.write("Method: sector-restricted Liouvillian "
                "(exact, sector conservation)\n\n")

        for m in all_meta:
            f.write(f"N = {m['N']}: "
                    f"γ = [{', '.join(f'{g:.4f}' for g in m['gamma'])}]\n")
            f.write(f"  Sector dim M = {m['sector_dim']}, "
                    f"L_sector = {m['liouvillian_dim']}x{m['liouvillian_dim']}\n")
            if m['slow_mode_rate'] is not None:
                f.write(f"  Slow SE mode rate = {m['slow_mode_rate']:.6f}\n")
            else:
                f.write(f"  No slow SE mode found\n")
        f.write("\n")

        hdr = (f"{'Label':35s}  {'N':>2}  {'Fam':>3}  "
               f"{'CΨmax(0)':>8}  {'Cusp':>4}  {'t_cross':>8}  "
               f"{'Overlap':>8}  {'Lens':>4}  {'Classification':>14}\n")
        f.write(hdr)
        f.write("-" * len(hdr.rstrip()) + "\n")
        for r in all_results:
            tc = f"{r['t_cross']:.3f}" if r['t_cross'] is not None else "n/a"
            f.write(
                f"{r['label']:35s}  {r['N']:>2}  {r['family']:>3}  "
                f"{r['cpsi0_max']:>8.4f}  "
                f"{'Y' if r['cusp_active'] else 'N':>4}  {tc:>8}  "
                f"{r['slow_mode_overlap']:>8.4f}  "
                f"{'Y' if r['lens_active'] else 'N':>4}  "
                f"{r['classification']:>14}\n")

        # Summary statistics
        classes = [r['classification'] for r in all_results]
        f.write(f"\nClassification counts:\n")
        for cls in ['straddles', 'cusp-only', 'lens-only', 'neither']:
            f.write(f"  {cls:>12}: {classes.count(cls)}\n")

        max_ov = max(r['slow_mode_overlap'] for r in all_results)
        best = max(all_results, key=lambda r: r['slow_mode_overlap'])
        f.write(f"\nMax slow-mode overlap: {max_ov:.6f} ({best['label']})\n")

        if max_ov < LENS_THRESHOLD:
            f.write(
                f"\nNo non-symmetric two-excitation state in the tested "
                f"families reaches\n"
                f"lens-accessibility above {LENS_THRESHOLD * 100:.0f}% "
                f"overlap. Combined with the symmetric-case\n"
                f"exclusion (SACRIFICE_GEOMETRY Level 2 line 112), this "
                f"substantially\n"
                f"resolves OQ-114 with negative evidence for the "
                f"two-excitation candidate\n"
                f"class. The boundary, if it exists, requires ansatze "
                f"beyond two-excitation\n"
                f"states.\n")
        else:
            straddlers = [r for r in all_results
                          if r['classification'] == 'straddles']
            if straddlers:
                f.write(f"\n{len(straddlers)} candidate(s) straddle both "
                        f"exits:\n")
                for s in straddlers:
                    f.write(f"  {s['label']} (N={s['N']}, "
                            f"overlap={s['slow_mode_overlap']:.4f}, "
                            f"t_cross={s['t_cross']:.3f})\n")

    # --- Overview plot ---
    fig, ax = plt.subplots(figsize=(10, 7))
    colors_N = {5: '#1f77b4', 6: '#ff7f0e', 7: '#2ca02c'}
    markers_fam = {'A': 'o', 'B': 's', 'C': '^'}

    seen = set()
    for r in all_results:
        x = r['t_cross'] if r['t_cross'] is not None else -0.5
        y = r['slow_mode_overlap']
        key = (r['N'], r['family'])
        ax.scatter(x, y,
                   c=colors_N.get(r['N'], 'gray'),
                   marker=markers_fam.get(r['family'], 'o'),
                   s=80, edgecolors='k', linewidth=0.5, zorder=3)
        seen.add(key)

    ax.axhline(LENS_THRESHOLD, color='red', linestyle=':', alpha=0.7,
               linewidth=1.5)

    from matplotlib.lines import Line2D
    handles, labels = [], []
    for key in sorted(seen):
        n, fam = key
        h = Line2D([0], [0], marker=markers_fam[fam], color='w',
                   markerfacecolor=colors_N[n], markeredgecolor='k',
                   markersize=9, linewidth=0)
        handles.append(h)
        labels.append(f'N={n} Family {fam}')
    handles.append(Line2D([0], [0], color='red', linestyle=':', linewidth=1.5))
    labels.append(f'Lens threshold ({LENS_THRESHOLD})')
    ax.legend(handles, labels, fontsize=8, loc='upper right')

    ax.set_xlabel('Cusp crossing time $t_{cross}$', fontsize=11)
    ax.set_ylabel('Slow-mode Frobenius overlap', fontsize=11)
    ax.set_title('Boundary straddling test: cusp crossing vs slow-mode overlap',
                 fontsize=12)
    ax.set_ylim(bottom=-0.005)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'boundary_straddling_overview.png'),
                dpi=150)
    plt.close(fig)

    print(f"\n  Results saved to {OUT_DIR}/")
    print(f"    sweep_results.json")
    print(f"    sweep_summary.txt")
    print(f"    boundary_straddling_overview.png")


# ===================================================================
# Main
# ===================================================================
if __name__ == "__main__":
    print("Boundary Straddling Sweep")
    print("Addresses OQ-114 + SACRIFICE_GEOMETRY OQ#2")
    print("Method: sector-restricted Liouvillian (exact)")
    print("=" * 70)

    all_results = []
    all_meta = []

    for N in [5, 6, 7]:
        try:
            results, meta = sweep_one_N(N)
            all_results.extend(results)
            all_meta.append(meta)
            # Incremental save after each N
            save_all(all_results, all_meta)
        except (MemoryError, np.linalg.LinAlgError) as e:
            print(f"\n  N={N}: failed ({e})")

    # Final summary
    print(f"\n{'=' * 70}")
    print(f"SWEEP COMPLETE: {len(all_results)} candidates across "
          f"N = {[m['N'] for m in all_meta]}")
    print(f"{'=' * 70}")

    classes = [r['classification'] for r in all_results]
    for cls in ['straddles', 'cusp-only', 'lens-only', 'neither']:
        n = classes.count(cls)
        if n > 0:
            print(f"  {cls:>12}: {n}")

    max_ov = max(r['slow_mode_overlap'] for r in all_results)
    best = max(all_results, key=lambda r: r['slow_mode_overlap'])
    print(f"\n  Max slow-mode overlap: {max_ov:.6f} ({best['label']}, N={best['N']})")

    if max_ov < LENS_THRESHOLD:
        print(f"\n  NEGATIVE RESULT: No candidate reaches the "
              f"{LENS_THRESHOLD * 100:.0f}% lens threshold.")
    else:
        straddlers = [r for r in all_results
                      if r['classification'] == 'straddles']
        print(f"\n  POSITIVE RESULT: {len(straddlers)} candidate(s) straddle "
              f"both exits.")
