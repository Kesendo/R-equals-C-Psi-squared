#!/usr/bin/env python3
"""PTF chiral mirror: the K1 site-wise trajectory identity (Edge 3 of the PTF chain, 2026-06-10).

THE LAW (EQ-014, review/EQ014_FINDINGS.md; "PTF's surviving Tier-1 law" in
reflections/ON_THE_Q_AXIS_AND_THE_PTF_LESSON.md): for the uniform XY chain (N sites, J=1)
under uniform Z-dephasing gamma_0 with a single J-bond defect dJ, the per-site first-order
rates f_i = (alpha_i - 1)/(dJ/J) from the PTF time-rescaling fits satisfy

    Sum_i f_i(psi_k) = Sum_i f_i(psi_{N+1-k})    (machine-exact at N = 5, 7, 8).

THE DERIVATION (verified here): the Sigma-law is a corollary of a strictly stronger
SITE-WISE TRAJECTORY IDENTITY. Let K1 = Prod_{l odd} Z_l (odd 0-indexed sublattice).
  Step 1 (algebra):     K1 H K1 = -H, K1 V K1 = -V (every XY bond has exactly one odd
                        endpoint), [K1, Z_l] = 0 (K1 is a Z string).
  Step 2 (conjugation): rho -> K1 rho K1 maps the Lindblad dynamics of H+V to that of
                        -(H+V) with the same dissipator; site purities are K1-invariant.
  Step 3 (reality):     H, V real in the computational basis and the Z-dephasing action
                        is real, so complex conjugation maps -(H+V) dynamics back to
                        +(H+V) dynamics; for real initial states purities are unchanged.
  Step 4 (modes):       K1 psi_k = psi_{N+1-k} EXACTLY (sine identity, 0-indexed odd
                        sublattice; the complementary even-sublattice product picks up
                        a minus sign, absorbed by the U(1) phase exp(i*pi*Nhat) =
                        Prod_l Z_l which commutes with H, V, and the dissipator).
HENCE  P_i(t; phi_k) = P_i(t; phi_{N+1-k})  for EVERY site i and ALL t, exactly, where
phi_k = (|vac> + |psi_k>)/sqrt(2) is the PTF pair state. The alpha_i fits are functionals
of the per-site purity trajectories, so alpha_i and f_i are SITE-WISE equal between k and
N+1-k, and the published Sigma f_i mirror law follows by summing over i. (No alpha re-fit
is needed here: the trajectory identity is strictly stronger than the fitted Sigma-law.)

Blocks (self-validating, every block raises on failure):
  1  algebra at N = 5: K1 H K1 = -H, K1 V K1 = -V, [K1, Z_l] = 0, H and V real,
     K1 psi_k = +psi_{N+1-k} for all k, even-sublattice product gives -psi_{N+1-k}.
  2  trajectory identity at N = 5 (gamma_0 = 0.05, dJ = 0.1 on bond (0,1)): site purities
     of phi_k vs phi_{N+1-k} for k = 1, 2 agree < 1e-10 at t in {0, 0.5, 1, 2, 4, 8},
     via dense expm of the 1024x1024 Liouvillian (vec convention self-checked).
  3  the U(1) minus-sign branch: chi = (|vac> - |psi_4>)/sqrt(2) (= K1c phi_2 with K1c
     the even-sublattice Z product, which maps psi_2 -> -psi_4) has the same site
     purities as BOTH phi_4 (the exp(i*pi*Nhat) absorption) and phi_2 (the full chain).
  4  N = 7 spot check (sparse expm_multiply, t in {0, 2, 4, 6}): phi_1 vs phi_7 and
     phi_2 vs phi_6 agree < 1e-9 (sparse builder self-checked against the framework
     dense Lindbladian at N = 4).

Run: python simulations/ptf_chiral_mirror_trajectory.py   (~30 s). 2026-06-10.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import (  # noqa: E402
    site_op,
    bonding_mode_state,
    bonding_mode_pair_state,
)
from framework.lindblad import lindbladian_z_dephasing  # noqa: E402

GAMMA = 0.05
DJ = 0.1


def xy_chain(N, J=1.0):
    """Uniform XY chain H = (J/2) Sum_l (X_l X_{l+1} + Y_l Y_{l+1}) (ChainSystem convention)."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for l in range(N - 1):
        H += (J / 2.0) * (site_op(N, l, 'X') @ site_op(N, l + 1, 'X')
                          + site_op(N, l, 'Y') @ site_op(N, l + 1, 'Y'))
    return H


def bond_defect(N, b, dJ):
    """V = dJ * (1/2)(X_b X_{b+1} + Y_b Y_{b+1}), the single J-bond defect."""
    return dJ * 0.5 * (site_op(N, b, 'X') @ site_op(N, b + 1, 'X')
                       + site_op(N, b, 'Y') @ site_op(N, b + 1, 'Y'))


def z_product(N, sites):
    """Prod_{l in sites} Z_l."""
    out = np.eye(2 ** N, dtype=complex)
    for l in sites:
        out = out @ site_op(N, l, 'Z')
    return out


def site_purities(rho, N):
    """[Tr(rho_i^2) for i in range(N)] via per-site partial trace."""
    t = rho.reshape((2,) * (2 * N))
    out = []
    for i in range(N):
        others = [a for a in range(N) if a != i]
        perm = [i] + others + [N + i] + [N + o for o in others]
        t2 = t.transpose(perm).reshape(2, 2 ** (N - 1), 2, 2 ** (N - 1))
        red = np.einsum('akbk->ab', t2)
        out.append(float(np.trace(red @ red).real))
    return np.array(out)


def block1():
    print('BLOCK 1 -- algebra at N = 5', flush=True)
    N = 5
    H = xy_chain(N)
    V = bond_defect(N, 0, DJ)
    K1 = z_product(N, [l for l in range(N) if l % 2 == 1])      # odd sublattice
    K1c = z_product(N, [l for l in range(N) if l % 2 == 0])     # even sublattice

    assert np.max(np.abs(K1 @ H @ K1 + H)) < 1e-14, 'K1 H K1 != -H'
    assert np.max(np.abs(K1 @ V @ K1 + V)) < 1e-14, 'K1 V K1 != -V'
    assert np.max(np.abs(K1c @ H @ K1c + H)) < 1e-14, 'K1c H K1c != -H'
    for l in range(N):
        Zl = site_op(N, l, 'Z')
        assert np.max(np.abs(K1 @ Zl - Zl @ K1)) < 1e-14, f'[K1, Z_{l}] != 0'
    assert np.max(np.abs(H.imag)) < 1e-14 and np.max(np.abs(V.imag)) < 1e-14, \
        'H or V not real in the computational basis'
    print('  K1 H K1 = -H, K1 V K1 = -V, [K1, Z_l] = 0, H and V real  OK', flush=True)

    for k in range(1, N + 1):
        psi_k = bonding_mode_state(N, k)
        psi_p = bonding_mode_state(N, N + 1 - k)
        dev_plus = np.max(np.abs(K1 @ psi_k - psi_p))
        dev_minus = np.max(np.abs(K1c @ psi_k + psi_p))
        assert dev_plus < 1e-14, f'K1 psi_{k} != +psi_{N + 1 - k} (dev {dev_plus:.2e})'
        assert dev_minus < 1e-14, f'K1c psi_{k} != -psi_{N + 1 - k} (dev {dev_minus:.2e})'
        print(f'  K1 psi_{k} = +psi_{N + 1 - k} ({dev_plus:.1e}), '
              f'K1c psi_{k} = -psi_{N + 1 - k} ({dev_minus:.1e})  OK', flush=True)
    print('BLOCK 1 PASS', flush=True)


def _dense_setup(N):
    """Dense Liouvillian for H + V (bond (0,1) defect) + vec-convention self-check."""
    H = xy_chain(N) + bond_defect(N, 0, DJ)
    L = lindbladian_z_dephasing(H, [GAMMA] * N)
    # Self-check the row-major vec convention: L @ vec_C(rho) == vec_C(RHS(rho)).
    rng = np.random.default_rng(14)
    d = 2 ** N
    R = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    rho = R + R.conj().T
    rhs = -1j * (H @ rho - rho @ H)
    for l in range(N):
        Zl = site_op(N, l, 'Z')
        rhs += GAMMA * (Zl @ rho @ Zl - rho)
    dev = np.max(np.abs((L @ rho.flatten()).reshape(d, d) - rhs))
    assert dev < 1e-10, f'vec convention self-check failed (dev {dev:.2e})'
    return L


def _purity_trajectories(U_half, states, N, n_steps=16, record={1, 2, 4, 8, 16}):
    """Propagate vec(rho) by repeated U_half = expm(L/2); return {label: {t: purities}}."""
    d = 2 ** N
    out = {}
    for label, psi in states.items():
        v = np.outer(psi, psi.conj()).flatten()
        traj = {0.0: site_purities(v.reshape(d, d), N)}
        for step in range(1, n_steps + 1):
            v = U_half @ v
            if step in record:
                traj[0.5 * step] = site_purities(v.reshape(d, d), N)
        out[label] = traj
    return out


def block2():
    print('BLOCK 2 -- trajectory identity at N = 5 (gamma=0.05, dJ=0.1 on bond (0,1))',
          flush=True)
    N = 5
    L = _dense_setup(N)
    print('  vec convention self-check OK', flush=True)
    U_half = expm(L * 0.5)
    states = {k: bonding_mode_pair_state(N, k) for k in (1, 2, 4, 5)}
    trajs = _purity_trajectories(U_half, states, N)
    worst = 0.0
    for k, kp in [(1, 5), (2, 4)]:
        dev = max(np.max(np.abs(trajs[k][t] - trajs[kp][t])) for t in trajs[k])
        worst = max(worst, dev)
        assert dev < 1e-10, f'phi_{k} vs phi_{kp}: max site-purity deviation {dev:.3e}'
        print(f'  phi_{k} vs phi_{kp}: max_i,t |P_i(t) deviation| = {dev:.3e}  OK', flush=True)
    print(f'BLOCK 2 PASS (worst {worst:.3e})', flush=True)
    return U_half, trajs


def block3(U_half, trajs):
    print('BLOCK 3 -- U(1) minus-sign branch: chi = (|vac> - |psi_4>)/sqrt(2)', flush=True)
    N = 5
    # chi = K1c phi_2 exactly (K1c = even-sublattice Z product maps psi_2 -> -psi_4,
    # vac -> +vac). Its purities must equal phi_4's (exp(i*pi*Nhat) absorption) and,
    # through the full chain, phi_2's.
    psi4 = bonding_mode_state(N, 4)
    chi = np.zeros(2 ** N, dtype=complex)
    chi[0] = 1.0
    chi -= psi4
    chi /= np.linalg.norm(chi)
    K1c = z_product(N, [0, 2, 4])
    dev_state = np.max(np.abs(K1c @ bonding_mode_pair_state(N, 2) - chi))
    assert dev_state < 1e-14, f'chi != K1c phi_2 (dev {dev_state:.2e})'
    traj_chi = _purity_trajectories(U_half, {'chi': chi}, N)['chi']
    dev4 = max(np.max(np.abs(traj_chi[t] - trajs[4][t])) for t in traj_chi)
    dev2 = max(np.max(np.abs(traj_chi[t] - trajs[2][t])) for t in traj_chi)
    assert dev4 < 1e-10, f'chi vs phi_4 deviation {dev4:.3e}'
    assert dev2 < 1e-10, f'chi vs phi_2 deviation {dev2:.3e}'
    print(f'  chi = K1c phi_2 exactly ({dev_state:.1e})', flush=True)
    print(f'  chi vs phi_4 (U(1) absorption): {dev4:.3e}  OK', flush=True)
    print(f'  chi vs phi_2 (full chain):      {dev2:.3e}  OK', flush=True)
    print('BLOCK 3 PASS', flush=True)


def _sparse_liouvillian(N):
    """Sparse row-major-vec Liouvillian for H + V with uniform Z-dephasing GAMMA."""
    d = 2 ** N
    H = sp.csr_matrix(xy_chain(N) + bond_defect(N, 0, DJ))
    I = sp.identity(d, format='csr', dtype=complex)
    L = -1j * (sp.kron(H, I, format='csr') - sp.kron(I, H.T, format='csr'))
    popc = np.array([bin(x).count('1') for x in range(d)])
    A, B = np.meshgrid(np.arange(d), np.arange(d), indexing='ij')
    mask = -2.0 * GAMMA * popc[(A ^ B)].ravel()          # row-major vec index a*d + b
    return (L + sp.diags(mask)).tocsc()


def block4():
    print('BLOCK 4 -- N = 7 spot check (sparse expm_multiply)', flush=True)
    # Self-check the sparse builder against the framework dense Lindbladian at N = 4.
    N_chk = 4
    L_dense = _dense_setup(N_chk)
    dev_b = np.max(np.abs(_sparse_liouvillian(N_chk).toarray() - L_dense))
    assert dev_b < 1e-12, f'sparse builder mismatch at N=4 (dev {dev_b:.2e})'
    print(f'  sparse builder vs framework dense at N=4: {dev_b:.1e}  OK', flush=True)

    N = 7
    d = 2 ** N
    L = _sparse_liouvillian(N)
    worst = 0.0
    for k, kp in [(1, 7), (2, 6)]:
        purities = {}
        for label in (k, kp):
            psi = bonding_mode_pair_state(N, label)
            v0 = np.outer(psi, psi.conj()).flatten()
            traj = expm_multiply(L, v0, start=0.0, stop=6.0, num=4)   # t = 0, 2, 4, 6
            purities[label] = [site_purities(traj[i].reshape(d, d), N) for i in range(4)]
        dev = max(np.max(np.abs(pa - pb))
                  for pa, pb in zip(purities[k], purities[kp]))
        worst = max(worst, dev)
        assert dev < 1e-9, f'phi_{k} vs phi_{kp} at N=7: deviation {dev:.3e}'
        print(f'  phi_{k} vs phi_{kp}: max_i,t |P_i(t) deviation| = {dev:.3e}  OK', flush=True)
    print(f'BLOCK 4 PASS (worst {worst:.3e})', flush=True)


if __name__ == '__main__':
    block1()
    U_half, trajs = block2()
    block3(U_half, trajs)
    block4()
    print('\nALL BLOCKS PASS -- the K1 chiral mirror is a site-wise trajectory identity:')
    print('P_i(t; phi_k) = P_i(t; phi_{N+1-k}) exactly for every site and time, so the')
    print('EQ-014 Sigma f_i mirror law is a corollary (site-wise f_i equality, summed).')
