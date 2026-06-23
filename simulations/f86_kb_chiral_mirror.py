#!/usr/bin/env python3
"""F86 Edge 3 (WIP): mode-resolved K_b under the K1 chiral mirror.

THE QUESTION: the ChiralMirrorTrajectoryClaim (PROOF_PTF_CHIRAL_MIRROR_RATE_LAW,
2026-06-10) gives P_i(t; phi_k) = P_i(t; phi_{N+1-k}) exactly for the pair states
phi_k = (|vac> + |psi_k>)/sqrt(2). Does the F86 susceptibility
K_b = 2*Re<rho(t)| S_kernel |d rho/dJ_b> (ResonanceScan convention, peak |.| over the
t-grid [0.6, 1.6]*t_peak with t_peak = 1/(4*gamma_0)) inherit the same k <-> N+1-k
mirror when computed from mode-resolved initial data?

DERIVATION (on a coherence block; the repo computes K_b on the (n, n+1) block):
L(J) = D + Sum_b J_b*M_b with D the REAL dephasing diagonal (-2*gamma*HD) and M_b
PURELY IMAGINARY (entries +-i; real hopping). K1 = Prod_{l odd} Z_l acts on the block
as the diagonal sign Sigma1[(p,q)] = (-1)^(n_odd(p)+n_odd(q)). Define the antiunitary

    T = Sigma1 o complex conjugation.

  (i)   generator even:     T L T = L         (Sigma1 M_b Sigma1 = -M_b because the
        bond is K1-odd; conj(M_b) = -M_b because the hopping is real; the two minus
        signs cancel -- the block shadow of "K1 flips H, conjugation flips it back").
  (ii)  perturbation even:  T M_b T = M_b     (same cancellation; the K_b direction
        d/dJ_b is itself K1-odd AND real, so it is T-even). Duhamel then gives
        d rho/dJ_b (t; T v0) = T d rho/dJ_b (t; v0).
  (iii) observable even:    Sigma1 S Sigma1 = S and S real. The K1 sign on the
        site-i kernel pair (p, p|bit_i) is the constant (-1)^[i odd], which squares
        away inside A_i^dag A_i; hence <Tv|S|Tw> = conj(<v|S|w>) and 2*Re is invariant.
  (iv)  initial data:       Sigma1 psi_k = psi_{N+1-k} (the ChiralKClaim gauge); the
        U(1)/parity signs drop because K_b is bilinear in v0 (phase-invariant).

PREDICTION: K_b(b, Q, t; mode-resolved probe) = K_b(b, Q, t; T-mirrored probe)
EXACTLY, at the SAME bond b, pointwise in Q and t. No bond reflection is induced:
K1 is site-diagonal. The F71 spatial mirror b <-> N-2-b is an INDEPENDENT exact
symmetry (F86c); sine modes are parity eigenstates, so it acts trivially on the mode
index, and the composed mirror (k -> N+1-k, b -> N-2-b) holds as the product. On a
non-palindromic J background only the K1 mode mirror survives.

THE c=1 SURPRISE (derived first, verified in block 2): the literal phi_k pair state
puts its K_b-relevant coherence in the (0,1) block, where L = i*h(J) - 2*gamma*I
with h Hermitian, so S(t) = 2*||v(t)||^2 = 2*e^(-4*gamma*t)*||v0||^2 is J-INDEPENDENT:
K_b(b; phi_k) == 0 identically (F73's c=1 spatial-sum closure). The trajectory
identity holds there vacuously. The non-trivial mode-resolved K_b lives on F86's
canonical c=2 block (1,2), probed with |psi_k><Psi_{k1,k2}|/2 (chi_k states), where
Psi_{k1,k2} is the 2-excitation sine-mode pair and T maps
(k; k1,k2) -> (N+1-k; N+1-k2, N+1-k1).

Blocks (self-validating):
  1  algebra: Sigma1/Mb/S/D T-evenness on the (1,2) block (N=5,6), R-evenness,
     Sigma1-mapping of the mode-resolved probes, single-particle K1 relations.
  2  conventions at N=5: block builders vs the framework full Liouvillian
     (row restriction), Duhamel vs central finite difference dJ=1e-6 (block space
     AND full 1024x1024 Lindblad propagation of the chi state), and the c=1
     K_b == 0 theorem for the literal phi_k.
  3  MAIN TABLES: (1,2) block, N = 5, 6, 7, Q in {1.5, 2.0}, all bonds, probes
     (k, (1,2)) for k = 1..N plus mixed spot probes; mode mirror / F71 bond
     mirror / composed mirror deviations (peak-K and pointwise-in-t).
  4  controls at N=6: J-defect background (mode mirror survives, bond mirror
     breaks), site-0 Z field (mode mirror breaks), uniform Z field (invisible
     to K_b entirely: constant E_p - E_q = global phase).
  5  Q-grid sweep (default 153-point grid) at N=6: max mode-mirror deviation.

Honest caveats: t_peak = 1/(4*gamma_0) is the F86a carrier law for the Dicke probe;
mode-resolved probes may peak elsewhere in (or at the edge of) the t-window, so the
peak-K numbers carry the window convention. The mirror itself is pointwise in t and
does not depend on it. Run: python simulations/f86_kb_chiral_mirror.py (~1 min).
2026-06-10.
"""
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import site_op, bonding_mode_pair_state  # noqa: E402
from framework.lindblad import lindbladian_z_dephasing  # noqa: E402

GAMMA = 0.05
T_PEAK = 1.0 / (4.0 * GAMMA)                       # F86a carrier law, = 5.0
T_GRID = np.linspace(0.6 * T_PEAK, 1.6 * T_PEAK, 21)   # ResonanceScan convention


def psi_sine(N, k):
    """F2b sine mode psi_k(j) = sqrt(2/(N+1))*sin(pi*k*(j+1)/(N+1)), j = 0..N-1."""
    j = np.arange(N)
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (j + 1) / (N + 1))


def pair_psi(N, k1, k2, qpairs):
    """Normalized 2-excitation sine-mode pair Psi_{k1,k2}(a,b) over the qpairs list."""
    p1, p2 = psi_sine(N, k1), psi_sine(N, k2)
    v = np.array([p1[a] * p2[b] - p1[b] * p2[a] for (a, b) in qpairs])
    return v / np.linalg.norm(v)


# ----------------------------------------------------------------------
# (1,2) coherence block: python twin of BlockLDecomposition + SpatialSumKernel
# ----------------------------------------------------------------------

def build_block12(N, gamma):
    """Basis (p=1exc site j, q=2exc pair (a,b)), D diag, per-bond M_b, S kernel,
    Sigma1 sign vector, R (site-reversal) permutation, bond reflection map."""
    qpairs = list(combinations(range(N), 2))
    basis = [(j, pq) for j in range(N) for pq in qpairs]
    flat = {bk: i for i, bk in enumerate(basis)}
    n = len(basis)
    nb = N - 1

    d_diag = np.array([-2.0 * gamma * (1 if j in pq else 3) for (j, pq) in basis],
                      dtype=complex)

    M = [np.zeros((n, n), dtype=complex) for _ in range(nb)]
    for col, (j, (a, b)) in enumerate(basis):
        for bnd in range(nb):
            # p-side hop (sign -i, matching BlockLDecomposition's hopP)
            if j == bnd:
                M[bnd][flat[(bnd + 1, (a, b))], col] += -1j
            elif j == bnd + 1:
                M[bnd][flat[(bnd, (a, b))], col] += -1j
            # q-side hops (sign +i, hopQ)
            for x, o in ((a, b), (b, a)):
                if x == bnd and bnd + 1 != o:
                    M[bnd][flat[(j, tuple(sorted((bnd + 1, o))))], col] += 1j
                elif x == bnd + 1 and bnd != o:
                    M[bnd][flat[(j, tuple(sorted((bnd, o))))], col] += 1j

    S = np.zeros((n, n), dtype=complex)
    site_groups = []
    for i in range(N):
        g = [flat[(j, tuple(sorted((j, i))))] for j in range(N) if j != i]
        site_groups.append(g)
        for x in g:
            for y in g:
                S[x, y] += 2.0

    sigma1 = np.array([(-1) ** ((j % 2) + (a % 2) + (b % 2)) for (j, (a, b)) in basis],
                      dtype=float)
    rperm = np.array([flat[(N - 1 - j, tuple(sorted((N - 1 - a, N - 1 - b))))]
                      for (j, (a, b)) in basis])

    def assemble(jvec, hz=None):
        L = np.diag(d_diag).astype(complex)
        for bnd in range(nb):
            L += jvec[bnd] * M[bnd]
        if hz is not None:
            # -i(E_p - E_q) = -2i*(sum_{l in q} h_l - h_j)
            f = np.array([-2j * (hz[a] + hz[b] - hz[j]) for (j, (a, b)) in basis])
            L += np.diag(f)
        return L

    def probe(k, kpair):
        psi = psi_sine(N, k)
        Psi = pair_psi(N, kpair[0], kpair[1], qpairs)
        v = np.zeros(n, dtype=complex)
        for idx, (j, pq) in enumerate(basis):
            v[idx] = psi[j] * Psi[qpairs.index(pq)] / 2.0
        return v

    return dict(N=N, basis=basis, flat=flat, qpairs=qpairs, n=n, nb=nb,
                d_diag=d_diag, M=M, S=S, sigma1=sigma1, rperm=rperm,
                assemble=assemble, probe=probe)


def mirror_label(N, k, kpair):
    return (N + 1 - k, tuple(sorted((N + 1 - kpair[0], N + 1 - kpair[1]))))


# ----------------------------------------------------------------------
# Duhamel K_b(t) engine (python twin of ResonanceScan.ScanAtQ)
# ----------------------------------------------------------------------

def k_curves(L, M_list, S, v0_list, t_grid=T_GRID):
    """Signed K[probe, bond, t] = 2*Re(<rho(t)| S |d rho/dJ_b>), Duhamel via eig."""
    evals, R = np.linalg.eig(L)
    Rinv = np.linalg.inv(R)
    c0s = [Rinv @ v for v in v0_list]
    Xs = [Rinv @ Mb @ R for Mb in M_list]
    K = np.zeros((len(v0_list), len(M_list), len(t_grid)))
    diff = evals[None, :] - evals[:, None]              # lambda_c - lambda_r
    deg = np.abs(diff) <= 1e-10
    denom = np.where(deg, 1.0, diff)
    for it, t in enumerate(t_grid):
        e = np.exp(evals * t)
        iM = np.where(deg, t * e[:, None], (e[None, :] - e[:, None]) / denom)
        for ip, c0 in enumerate(c0s):
            rho_t = R @ (e * c0)
            for ib, X in enumerate(Xs):
                drho = R @ ((X * iM) @ c0)
                K[ip, ib, it] = 2.0 * np.vdot(rho_t, S @ drho).real
    return K


def peak(K):
    """ResonanceScan peak convention: kAbs = max_t |K|, plus t at peak."""
    idx = np.argmax(np.abs(K), axis=-1)
    kpk = np.take_along_axis(np.abs(K), idx[..., None], axis=-1)[..., 0]
    tpk = T_GRID[idx]
    return kpk, tpk


# ----------------------------------------------------------------------
# Full-space helpers (validation only, N=5)
# ----------------------------------------------------------------------

def h_full(N, jvec, hz=None):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        H += (jvec[b] / 2.0) * (site_op(N, b, 'X') @ site_op(N, b + 1, 'X')
                                + site_op(N, b, 'Y') @ site_op(N, b + 1, 'Y'))
    if hz is not None:
        for l in range(N):
            H += hz[l] * site_op(N, l, 'Z')
    return H


def spatial_sum_full(rho, N):
    """S = Sum_i 2*|(rho_i)_{01}|^2 via per-site partial trace."""
    t = rho.reshape((2,) * (2 * N))
    s = 0.0
    for i in range(N):
        others = [a for a in range(N) if a != i]
        perm = [i] + others + [N + i] + [N + o for o in others]
        t2 = t.transpose(perm).reshape(2, 2 ** (N - 1), 2, 2 ** (N - 1))
        red = np.einsum('akbk->ab', t2)
        s += 2.0 * abs(red[0, 1]) ** 2
    return s


def bit(N, j):
    return 1 << (N - 1 - j)        # big-endian, framework/C# convention


def chi_state(N, k, kpair, qpairs):
    """chi = (|1exc psi_k> + |2exc Psi_{k1,k2}>)/sqrt(2); its (1,2) coherence block
    equals the block probe psi_k Psi^T / 2."""
    psi, Psi = psi_sine(N, k), pair_psi(N, kpair[0], kpair[1], qpairs)
    v = np.zeros(2 ** N, dtype=complex)
    for j in range(N):
        v[bit(N, j)] = psi[j] / np.sqrt(2)
    for (a, b), w in zip(qpairs, Psi):
        v[bit(N, a) | bit(N, b)] = w / np.sqrt(2)
    return v


# ----------------------------------------------------------------------
# Block 1: algebra
# ----------------------------------------------------------------------

def block1():
    print('BLOCK 1 -- algebra: T-evenness on the (1,2) block, probe mapping', flush=True)
    for N in (5, 6):
        blk = build_block12(N, GAMMA)
        s1 = blk['sigma1']
        worst_modd = max(np.max(np.abs(s1[:, None] * Mb * s1[None, :] + Mb))
                         for Mb in blk['M'])
        worst_mreal = max(np.max(np.abs(Mb.real)) for Mb in blk['M'])
        worst_teven = max(np.max(np.abs(s1[:, None] * np.conj(Mb) * s1[None, :] - Mb))
                          for Mb in blk['M'])
        S = blk['S']
        dev_s = np.max(np.abs(s1[:, None] * S * s1[None, :] - S)) \
            + np.max(np.abs(S.imag))
        assert worst_modd < 1e-14, 'Sigma1 M_b Sigma1 != -M_b'
        assert worst_mreal < 1e-14, 'M_b not purely imaginary'
        assert worst_teven < 1e-14, 'T M_b T != M_b'
        assert dev_s < 1e-14, 'S not Sigma1-even/real'
        # R-evenness: P_R M_b P_R = M_{N-2-b}, P_R S P_R = S, D R-invariant
        rp = blk['rperm']
        dev_r = max(np.max(np.abs(blk['M'][b][np.ix_(rp, rp)] - blk['M'][N - 2 - b]))
                    for b in range(N - 1))
        dev_rs = np.max(np.abs(S[np.ix_(rp, rp)] - S))
        dev_rd = np.max(np.abs(blk['d_diag'][rp] - blk['d_diag']))
        assert dev_r < 1e-14 and dev_rs < 1e-14 and dev_rd < 1e-14, 'R-evenness fails'
        # Probe mapping: Sigma1 v0(k, pair) = +- v0(mirror)
        worst_map = 0.0
        for (k, kp) in [(1, (1, 2)), (2, (1, 3)), (3, (2, 3))]:
            mk, mkp = mirror_label(N, k, kp)
            v, vm = blk['probe'](k, kp), blk['probe'](mk, mkp)
            dev = min(np.max(np.abs(s1 * v - vm)), np.max(np.abs(s1 * v + vm)))
            worst_map = max(worst_map, dev)
        assert worst_map < 1e-14, 'Sigma1 probe mapping fails'
        print(f'  N={N}: Sigma1 M_b Sigma1 = -M_b ({worst_modd:.1e}), '
              f'T M_b T = M_b ({worst_teven:.1e}), Sigma1 S Sigma1 = S real ({dev_s:.1e}),')
        print(f'        R: M_b -> M_(N-2-b) ({dev_r:.1e}), S, D invariant '
              f'({max(dev_rs, dev_rd):.1e}); Sigma1 probe(k,pair) = ±probe(mirror) '
              f'({worst_map:.1e})  OK', flush=True)
    # single-particle layer: Sigma1 h Sigma1 = -h, Sigma1 psi_k = psi_{N+1-k}
    N = 6
    sgn = np.array([(-1) ** l for l in range(N)], dtype=float)
    h = np.diag(np.ones(N - 1), 1) + np.diag(np.ones(N - 1), -1)
    assert np.max(np.abs(sgn[:, None] * h * sgn[None, :] + h)) < 1e-14
    for k in range(1, N + 1):
        assert np.max(np.abs(sgn * psi_sine(N, k) - psi_sine(N, N + 1 - k))) < 1e-13
    print(f'  single-particle: Sigma1 h Sigma1 = -h, Sigma1 psi_k = psi_(N+1-k) '
          f'for all k (N={N})  OK', flush=True)
    print('BLOCK 1 PASS\n', flush=True)


# ----------------------------------------------------------------------
# Block 2: conventions at N=5 (+ the c=1 zero theorem)
# ----------------------------------------------------------------------

def block2():
    print('BLOCK 2 -- conventions at N = 5 (Q = 2.0, J = Q*gamma = 0.1)', flush=True)
    N, Q = 5, 2.0
    J = Q * GAMMA
    jvec = [J] * (N - 1)
    d = 2 ** N
    blk = build_block12(N, GAMMA)

    # (a) block builder vs framework full Liouvillian, restricted to (1,2) rows
    L_full = lindbladian_z_dephasing(h_full(N, jvec), [GAMMA] * N)
    rows = [bit(N, j) * d + (bit(N, a) | bit(N, b)) for (j, (a, b)) in blk['basis']]
    dev_restrict = np.max(np.abs(L_full[np.ix_(rows, rows)] - blk['assemble'](jvec)))
    other = [i for i in range(d * d) if i not in set(rows)]
    dev_leak = np.max(np.abs(L_full[np.ix_(other, rows)]))
    assert dev_restrict < 1e-12, f'(1,2) restriction mismatch {dev_restrict:.2e}'
    assert dev_leak < 1e-12, f'(1,2) block leaks {dev_leak:.2e}'
    print(f'  (1,2) block = full-L restriction ({dev_restrict:.1e}), '
          f'no out-of-block coupling ({dev_leak:.1e})  OK', flush=True)

    # (b) Duhamel vs block-space central FD (dJ = 1e-6), chi probe (2,(1,2))
    k, kp = 2, (1, 2)
    v0 = blk['probe'](k, kp)
    L = blk['assemble'](jvec)
    Kd = k_curves(L, blk['M'], blk['S'], [v0])
    delta = 1e-6
    it5 = int(np.argmin(np.abs(T_GRID - T_PEAK)))     # t = 5.0 exactly on grid
    worst_fd = 0.0
    for b in (0, 2):
        sv = []
        for sgn in (+1, -1):
            jv = list(jvec)
            jv[b] += sgn * delta
            vt = expm(blk['assemble'](jv) * T_PEAK) @ v0
            sv.append(np.vdot(vt, blk['S'] @ vt).real)
        kfd = (sv[0] - sv[1]) / (2 * delta)
        dev = abs(kfd - Kd[0, b, it5])
        worst_fd = max(worst_fd, dev)
        assert dev < 1e-7 + 1e-5 * abs(kfd), f'Duhamel vs FD bond {b}: {dev:.2e}'
        print(f'  bond {b}: Duhamel K = {Kd[0, b, it5]:+.8e}, block FD = {kfd:+.8e} '
              f'(dev {dev:.1e})  OK', flush=True)

    # (c) full-space FD with the physical chi state (validates probe + observable)
    chi = chi_state(N, k, kp, blk['qpairs'])
    b = 0
    sv = []
    for sgn in (+1, -1):
        jv = list(jvec)
        jv[b] += sgn * delta
        Lf = lindbladian_z_dephasing(h_full(N, jv), [GAMMA] * N)
        rho_t = (expm(Lf * T_PEAK) @ np.outer(chi, chi.conj()).flatten()).reshape(d, d)
        sv.append(spatial_sum_full(rho_t, N))
    kfd_full = (sv[0] - sv[1]) / (2 * delta)
    dev = abs(kfd_full - Kd[0, b, it5])
    assert dev < 1e-7 + 1e-5 * abs(kfd_full), f'full-space FD mismatch {dev:.2e}'
    print(f'  full 1024x1024 Lindblad FD on chi state, bond 0: {kfd_full:+.8e} '
          f'(dev vs block Duhamel {dev:.1e})  OK', flush=True)

    # (d) the c=1 zero theorem: literal phi_k has K_b == 0 (F73 closure)
    phi = bonding_mode_pair_state(N, 2)
    sv = []
    for sgn in (+1, -1):
        jv = list(jvec)
        jv[0] += sgn * delta
        Lf = lindbladian_z_dephasing(h_full(N, jv), [GAMMA] * N)
        rho_t = (expm(Lf * T_PEAK) @ np.outer(phi, phi.conj()).flatten()).reshape(d, d)
        sv.append(spatial_sum_full(rho_t, N))
    kfd_phi = (sv[0] - sv[1]) / (2 * delta)
    # block-level: (0,1) generator is i*h - 2*gamma*I, unitary x scalar decay
    h1 = np.zeros((N, N))
    for bb in range(N - 1):
        h1[bb, bb + 1] = h1[bb + 1, bb] = J
    v0_01 = psi_sine(N, 2) / 2.0
    vt = (expm((1j * h1 - 2 * GAMMA * np.eye(N)) * T_PEAK) @ v0_01)
    s_pred = 2.0 * np.exp(-4 * GAMMA * T_PEAK) * np.linalg.norm(v0_01) ** 2
    dev_s = abs(2.0 * np.linalg.norm(vt) ** 2 - s_pred)
    assert abs(kfd_phi) < 1e-8, f'phi_k K_b not zero: {kfd_phi:.2e}'
    assert dev_s < 1e-14, f'c=1 closure violated: {dev_s:.2e}'
    print(f'  literal phi_2: full-space FD K_b = {kfd_phi:+.2e} (== 0; F73 c=1 '
          f'closure S(t) = 2 e^(-4 gamma t) ||v0||^2, dev {dev_s:.1e})  OK', flush=True)
    print('BLOCK 2 PASS\n', flush=True)


# ----------------------------------------------------------------------
# Block 3: main mirror tables on the (1,2) block
# ----------------------------------------------------------------------

def probe_set(N):
    """Base probes: (k, (1,2)) for all k, plus mixed spot probes."""
    base = [(k, (1, 2)) for k in range(1, N + 1)]
    base += [(2, (1, 3)), (1, (2, 3))]
    full, seen = [], set()
    for lab in base:
        for cand in (lab, mirror_label(N, *lab)):
            if cand not in seen:
                seen.add(cand)
                full.append(cand)
    return base, full


def run_case(N, jvec, hz=None):
    blk = build_block12(N, GAMMA)
    base, full = probe_set(N)
    pidx = {lab: i for i, lab in enumerate(full)}
    v0s = [blk['probe'](k, kp) for (k, kp) in full]
    L = blk['assemble'](jvec, hz)
    K = k_curves(L, blk['M'], blk['S'], v0s)
    kpk, tpk = peak(K)
    return blk, base, pidx, K, kpk, tpk


def mirror_devs(N, base, pidx, K, kpk):
    """Max deviations (relative to global peak K) for the three mirrors."""
    kmax = np.max(kpk)
    dev_mode = dev_bond = dev_comp = 0.0
    dev_mode_pt = 0.0
    for lab in base:
        i, im = pidx[lab], pidx[mirror_label(N, *lab)]
        for b in range(N - 1):
            dev_mode = max(dev_mode, abs(kpk[i, b] - kpk[im, b]))
            dev_mode_pt = max(dev_mode_pt, np.max(np.abs(K[i, b] - K[im, b])))
            dev_bond = max(dev_bond, abs(kpk[i, b] - kpk[i, N - 2 - b]))
            dev_comp = max(dev_comp, abs(kpk[i, b] - kpk[im, N - 2 - b]))
    return (dev_mode / kmax, dev_mode_pt / kmax, dev_bond / kmax,
            dev_comp / kmax, kmax)


def block3():
    print('BLOCK 3 -- mode-resolved K_b mirrors on the (1,2) block '
          '(gamma_0 = 0.05, t-grid [0.6, 1.6]*t_peak)', flush=True)
    summary = []
    for N in (5, 6, 7):
        for Q in (1.5, 2.0):
            jvec = [Q * GAMMA] * (N - 1)
            blk, base, pidx, K, kpk, tpk = run_case(N, jvec)
            dm, dmpt, db, dc, kmax = mirror_devs(N, base, pidx, K, kpk)
            summary.append((N, Q, dm, dmpt, db, dc, kmax))
            if N == 6 and Q == 2.0:
                print(f'\n  Full table at N=6, Q=2.0 (J=0.1); columns: peak K_b per '
                      f'bond [t_at_peak]; bonds 0,4 = Endpoint, 1-3 = Interior')
                hdr = ' | '.join(f'b={b}' for b in range(N - 1))
                print(f'    probe (k; k1,k2)         | {hdr}')
                for lab in base:
                    i = pidx[lab]
                    cells = ' | '.join(f'{kpk[i, b]:.6e}[{tpk[i, b]:.2f}]'
                                       for b in range(N - 1))
                    print(f'    ({lab[0]}; {lab[1][0]},{lab[1][1]})'.ljust(29)
                          + f'| {cells}')
                    mlab = mirror_label(N, *lab)
                    im = pidx[mlab]
                    cells = ' | '.join(f'{kpk[im, b]:.6e}[{tpk[im, b]:.2f}]'
                                       for b in range(N - 1))
                    print(f'    ({mlab[0]}; {mlab[1][0]},{mlab[1][1]}) [mirror]'
                          .ljust(29) + f'| {cells}')
                    dev = max(abs(kpk[i, b] - kpk[im, b]) for b in range(N - 1))
                    print(f'      -> mode-mirror dev (same b): {dev:.2e}')
                print()
    print('  Summary (deviations relative to the global peak K of each case):')
    print('    N | Q   | mode k<->N+1-k (peak) | mode (pointwise t) | '
          'F71 bond b<->N-2-b | composed | max K')
    for (N, Q, dm, dmpt, db, dc, kmax) in summary:
        print(f'    {N} | {Q:.1f} | {dm:.2e}             | {dmpt:.2e}'
              f'           | {db:.2e}           | {dc:.2e} | {kmax:.3e}')
        assert dm < 1e-10, f'mode mirror broken at N={N}, Q={Q}: {dm:.2e}'
        assert dmpt < 1e-10, f'pointwise mode mirror broken: {dmpt:.2e}'
        assert db < 1e-10, f'F71 bond mirror broken: {db:.2e}'
        assert dc < 1e-10, f'composed mirror broken: {dc:.2e}'
    # non-triviality: K depends on the mode index
    N, Q = 6, 2.0
    blk, base, pidx, K, kpk, tpk = run_case(N, [Q * GAMMA] * (N - 1))
    col = [kpk[pidx[(k, (1, 2))], 0] for k in range(1, N + 1)]
    spread = (max(col) - min(col)) / max(col)
    assert spread > 1e-3, f'K does not resolve the mode index (spread {spread:.2e})'
    print(f'  Non-triviality: K(b=0; k,(1,2)) spread over k at N=6 is '
          f'{spread:.1%} of max -- the mirror is not vacuous.')
    print('BLOCK 3 PASS\n', flush=True)


# ----------------------------------------------------------------------
# Block 4: controls at N=6
# ----------------------------------------------------------------------

def block4():
    print('BLOCK 4 -- controls at N = 6, Q = 2.0', flush=True)
    N, Q = 6, 2.0
    J = Q * GAMMA

    # (a) J-defect background: K1 mirror survives, F71 bond mirror breaks
    jvec = [J] * (N - 1)
    jvec[0] *= 1.3
    blk, base, pidx, K, kpk, tpk = run_case(N, jvec)
    dm, dmpt, db, dc, kmax = mirror_devs(N, base, pidx, K, kpk)
    assert dm < 1e-10 and dmpt < 1e-10, f'mode mirror broke on defect: {dm:.2e}'
    assert db > 1e-4, f'bond mirror unexpectedly survived the defect: {db:.2e}'
    print(f'  (a) J-defect (+30% on bond 0): mode mirror {dm:.2e} / pointwise '
          f'{dmpt:.2e} (EXACT), bond mirror {db:.2e} (BROKEN), composed {dc:.2e}')

    # (b) single-site Z field on site 0: K1-even term, mode mirror breaks
    hz = [0.0] * N
    hz[0] = 0.2
    blk, base, pidx, K, kpk, tpk = run_case(N, [J] * (N - 1), hz)
    dm, dmpt, db, dc, kmax = mirror_devs(N, base, pidx, K, kpk)
    assert dm > 1e-6, f'mode mirror unexpectedly survived the Z field: {dm:.2e}'
    print(f'  (b) Z field h=0.2 on site 0: mode mirror {dm:.2e} / pointwise '
          f'{dmpt:.2e} (BROKEN), bond mirror {db:.2e} (broken too: field breaks R)')

    # (c) uniform Z field: E_p - E_q = 2h constant on any (n, n+1) block ->
    #     global phase -> K_b identical to the no-field case
    blk = build_block12(N, GAMMA)
    base, full = probe_set(N)
    v0s = [blk['probe'](k, kp) for (k, kp) in full]
    K0 = k_curves(blk['assemble']([J] * (N - 1)), blk['M'], blk['S'], v0s)
    Kh = k_curves(blk['assemble']([J] * (N - 1), [0.2] * N), blk['M'], blk['S'], v0s)
    dev_u = np.max(np.abs(K0 - Kh)) / np.max(np.abs(K0))
    assert dev_u < 1e-10, f'uniform field not invisible: {dev_u:.2e}'
    print(f'  (c) uniform Z field h=0.2: K_b(h) vs K_b(0) deviation {dev_u:.2e} '
          f'(INVISIBLE -- constant E_p - E_q is a global block phase)')
    print('BLOCK 4 PASS\n', flush=True)


# ----------------------------------------------------------------------
# Block 5: Q-grid sweep
# ----------------------------------------------------------------------

def block5():
    print('BLOCK 5 -- mode-mirror over the default 153-point Q grid, N = 6', flush=True)
    N = 6
    blk = build_block12(N, GAMMA)
    lab, mlab = (2, (1, 2)), mirror_label(N, 2, (1, 2))
    v0s = [blk['probe'](*lab), blk['probe'](*mlab)]
    qgrid = np.linspace(0.20, 4.00, 153)
    worst = 0.0
    kmax = 0.0
    for Q in qgrid:
        K = k_curves(blk['assemble']([Q * GAMMA] * (N - 1)), blk['M'], blk['S'], v0s)
        worst = max(worst, np.max(np.abs(K[0] - K[1])))
        kmax = max(kmax, np.max(np.abs(K)))
    rel = worst / kmax
    assert rel < 1e-10, f'Q-sweep mode mirror broken: {rel:.2e}'
    print(f'  probes {lab} vs {mlab}: max pointwise |Delta K(b, Q, t)| = {worst:.2e}'
          f' ({rel:.2e} of max K = {kmax:.3e}) across Q in [0.20, 4.00]  OK')
    print('BLOCK 5 PASS\n', flush=True)


if __name__ == '__main__':
    block1()
    block2()
    block3()
    block4()
    block5()
    print('ALL BLOCKS PASS -- verdict:')
    print('  1. The literal phi_k pair states give K_b == 0 identically (F73 c=1')
    print('     closure: the (0,1) block evolves as e^(-2 gamma t) x unitary, so the')
    print('     spatial sum has no J-susceptibility). The trajectory identity holds')
    print('     there vacuously; F86\'s K_b needs c >= 2 precisely because dephasing')
    print('     must act non-uniformly (HD 1 vs 3) inside the block.')
    print('  2. On the canonical c=2 block, mode-resolved K_b obeys the K1 chiral')
    print('     mirror EXACTLY at the SAME bond b, pointwise in Q and t:')
    print('     K_b(b; k, (k1,k2)) = K_b(b; N+1-k, (N+1-k2, N+1-k1)).')
    print('     The F71 bond mirror b <-> N-2-b and the composed mirror are')
    print('     independent exact symmetries; a J-defect background kills the bond')
    print('     mirror but not the K1 mirror; a non-uniform Z field kills the K1')
    print('     mirror; a uniform Z field is invisible to K_b entirely.')
