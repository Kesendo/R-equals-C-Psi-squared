"""
Pi Time Reversal Verification -- March 19, 2026
=================================================
Verifies that Pi maps eigenvectors of L to their palindromic partners.

Steps:
  1. Build Pi (64x64) and verify fundamental identity
  2. Apply Pi to each eigenvector, verify it lands in partner eigenspace
  3. Verify XY-weight swap: w(r_k) + w(Pi|r_k>) = N
  4. Demonstrate on Bell(0,1)'s oscillating modes (omega ~ 2, 4, 6)
  5. Demonstrate Pi(ZZZ) = -i*YYY

Script: simulations/pi_time_reversal_verify.py
Output: simulations/results/pi_time_reversal_verify.txt
"""
import numpy as np
from itertools import product as iproduct
from datetime import datetime

OUT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results\pi_time_reversal_verify.txt"
f = open(OUT, "w", buffering=1)

def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n"); f.flush()

# ============================================================
# PAULI BASICS
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
NAMES = ['I', 'X', 'Y', 'Z']

PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}   # I<->X, Y<->Z
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}  # phase factors


def xy_weight(indices):
    return sum(1 for i in indices if i in (1, 2))


def plabel(indices):
    return ''.join(NAMES[i] for i in indices)


def avg_xy_wt(v, all_idx):
    """Weighted average XY-weight of a Pauli-basis vector."""
    w = np.array([xy_weight(idx) for idx in all_idx], dtype=float)
    v2 = np.abs(v)**2
    s = np.sum(v2)
    return np.dot(w, v2) / s if s > 1e-30 else 0.0


def wt_dist(v, all_idx, N):
    """Fraction of |v|^2 in each XY-weight sector (0..N)."""
    v2 = np.abs(v)**2
    s = np.sum(v2)
    if s < 1e-30:
        return np.zeros(N + 1)
    d = np.zeros(N + 1)
    for a, idx in enumerate(all_idx):
        d[xy_weight(idx)] += v2[a]
    return d / s


def fmt(z, p=6):
    r, i = z.real, z.imag
    if abs(i) < 10**(-p):
        return f"{r:+.{p}f}"
    if abs(r) < 10**(-p):
        return f"{i:+.{p}f}i"
    return f"{r:+.{p}f}{i:+.{p}f}i"


def fmtp(z):
    if abs(z - 1) < 1e-10: return "+1"
    if abs(z + 1) < 1e-10: return "-1"
    if abs(z - 1j) < 1e-10: return "+i"
    if abs(z + 1j) < 1e-10: return "-i"
    return fmt(z, 4)


# ============================================================
# BUILD SYSTEM IN PAULI BASIS
# ============================================================
def build_H_xxz(N, bonds, J=1.0, delta=1.0):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        for pidx, pauli in enumerate([sx, sy, sz]):
            term = np.eye(1, dtype=complex)
            for k in range(N):
                term = np.kron(term, pauli if k in (i, j) else I2)
            H += J * (delta if pidx == 2 else 1.0) * term
    return H


def build_system(N, H, gammas):
    d = 2**N
    num = 4**N
    all_idx = list(iproduct(range(4), repeat=N))

    # Build Pauli operator matrices
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for i in idx[1:]:
            m = np.kron(m, PAULIS[i])
        pmats.append(m)

    # Hamiltonian part of Liouvillian in Pauli basis
    L_H = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)
        for a in range(num):
            L_H[a, b] = np.trace(pmats[a] @ comm) / d

    # Dephasing diagonal
    L_D_diag = np.zeros(num)
    for a, idx in enumerate(all_idx):
        L_D_diag[a] = -sum(2 * gammas[s] for s in range(N) if idx[s] in (1, 2))

    L = L_H.copy()
    for a in range(num):
        L[a, a] += L_D_diag[a]

    # Pi operator
    Pi = np.zeros((num, num), dtype=complex)
    for b, idx_b in enumerate(all_idx):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b:
            sign *= PI_SIGN[i]
        Pi[all_idx.index(mapped), b] = sign

    return L, L_D_diag, Pi, all_idx, pmats


# ============================================================
# SETUP
# ============================================================
N = 3
gamma = 0.05
gammas = [gamma] * N
Sg = sum(gammas)   # 0.15
bonds = [(i, i + 1) for i in range(N - 1)]
num = 4**N         # 64
dim = 2**N         # 8

H = build_H_xxz(N, bonds)
L, L_D_diag, Pi, all_idx, pmats = build_system(N, H, gammas)

log("=" * 90)
log("PI TIME REVERSAL VERIFICATION")
log(f"Date: {datetime.now()}")
log(f"N={N}, gamma={gamma}, Heisenberg chain, bonds={bonds}")
log(f"Sg = sum(gamma) = {Sg}, 2*Sg = {2*Sg}")
log("=" * 90)


# ############################################################
# SECTION 1: Pi construction
# ############################################################
log()
log("=" * 90)
log("SECTION 1: Pi construction")
log("=" * 90)

log(f"\n  {num}x{num} matrix, Pauli-vectorized basis")

# Unitary check
uerr = np.max(np.abs(Pi @ Pi.conj().T - np.eye(num)))
log(f"  Unitary check: max|Pi*Pi^H - I| = {uerr:.2e}")
log(f"    -> {'UNITARY' if uerr < 1e-10 else 'NOT UNITARY'}")

# Pi^2 phases by XY-weight sector
Pi_sq = Pi @ Pi
log(f"\n  Pi^2 phases by XY-weight sector:")
for w in range(N + 1):
    indices = [a for a, idx in enumerate(all_idx) if xy_weight(idx) == w]
    phase = Pi_sq[indices[0], indices[0]]
    log(f"    w_xy = {w}: Pi^2 = {fmtp(phase):>3}  ({len(indices)} Pauli strings)")

# Fundamental identity
Pi_inv = np.linalg.inv(Pi)
ferr = np.max(np.abs(Pi @ L @ Pi_inv + L + 2 * Sg * np.eye(num)))
log(f"\n  Fundamental identity: Pi*L*Pi^(-1) = -L - 2*Sg*I")
log(f"    Max error: {ferr:.2e}")
log(f"    -> {'VERIFIED' if ferr < 1e-10 else 'FAILED'}")

# Per-site action
log(f"\n  Per-site map:")
for i in range(4):
    log(f"    {NAMES[i]} -> {fmtp(PI_SIGN[i]):>3} * {NAMES[PI_PERM[i]]}")


# ############################################################
# SECTION 2: Eigenvector mapping
# ############################################################
log()
log("=" * 90)
log("SECTION 2: Eigenvector mapping")
log("  For each right eigenvector |r_k>, verify Pi|r_k> lies in")
log("  the eigenspace of the palindromic partner lambda_k' = -(lambda_k + 2*Sg)")
log("=" * 90)

evals, R = np.linalg.eig(L)
mu = evals + Sg

# Check eigendecomposition quality
eig_errs = [np.linalg.norm(L @ R[:, k] - evals[k] * R[:, k])
            for k in range(num)]
log(f"\n  Eigendecomposition: max|L*r_k - lambda_k*r_k| = {max(eig_errs):.2e}")

# Palindromic pairing (greedy)
paired = np.zeros(num, dtype=bool)
pair_map = -np.ones(num, dtype=int)

for k in range(num):
    if paired[k]:
        continue
    target = -(evals[k] + 2 * Sg)
    best_j, best_d = -1, np.inf
    for j in range(k + 1, num):
        if paired[j]:
            continue
        d = abs(evals[j] - target)
        if d < best_d:
            best_d, best_j = d, j
    if best_j >= 0 and best_d < 1e-6:
        paired[k] = paired[best_j] = True
        pair_map[k] = best_j
        pair_map[best_j] = k

n_paired = int(np.sum(paired))
log(f"  Palindromic pairing: {n_paired}/{num} eigenvalues = {n_paired // 2} pairs")
if n_paired < num:
    log(f"  WARNING: {num - n_paired} eigenvalues unpaired!")

# Verify eigenvector mapping for each pair
hdr_lp = "lambda_k'"
log(f"\n  {'#':>3}  {'lambda_k':>24}  {hdr_lp:>24}  "
    f"{'eval_err':>10}  {'vec_resid':>10}")
log(f"  {'-'*78}")

n_p = 0
max_res = 0
seen = set()

for k in range(num):
    j = pair_map[k]
    if j < 0 or k in seen:
        continue
    seen.add(k)
    seen.add(j)
    n_p += 1

    ev_err = abs(evals[j] - (-(evals[k] + 2 * Sg)))

    # Apply Pi to eigenvector k
    Pi_rk = Pi @ R[:, k]

    # Partner eigenspace = all eigenvectors with eigenvalue close to evals[j]
    # (handles degenerate eigenvalues correctly)
    pidx = [m for m in range(num) if abs(evals[m] - evals[j]) < 1e-6]
    V = R[:, pidx]

    # Project Pi|r_k> onto partner eigenspace, compute residual
    proj = V @ (np.linalg.pinv(V) @ Pi_rk)
    res = np.linalg.norm(Pi_rk - proj) / np.linalg.norm(Pi_rk)
    max_res = max(max_res, res)

    log(f"  {n_p:>3}  {fmt(evals[k]):>24}  {fmt(evals[j]):>24}  "
        f"{ev_err:>10.2e}  {res:>10.2e}")

log(f"\n  All {n_p} pairs verified: {'YES' if max_res < 1e-6 else 'NO'}")
log(f"  Max eigenvector residual: {max_res:.2e}")


# ############################################################
# SECTION 3: XY-weight swap
# ############################################################
log()
log("=" * 90)
log("SECTION 3: XY-weight swap")
log("  For each pair: w(r_k) + w(Pi|r_k>) = N (algebraic identity)")
log("  Also show XY-weight distribution of each eigenvector")
log("=" * 90)

seen = set()
sums = []

hdr_wdp = "weight dist k'"
log(f"\n  {'#':>3}  {'w_k':>7}  {'w_Pi':>7}  {'sum':>8}  "
    f"{'weight dist k':>24}  {hdr_wdp:>24}")
log(f"  {'-'*82}")

for k in range(num):
    j = pair_map[k]
    if j < 0 or k in seen:
        continue
    seen.add(k)
    seen.add(j)

    wk = avg_xy_wt(R[:, k], all_idx)
    w_pi = avg_xy_wt(Pi @ R[:, k], all_idx)
    sums.append(wk + w_pi)

    dk = wt_dist(R[:, k], all_idx, N)
    dj = wt_dist(R[:, j], all_idx, N)
    dks = ' '.join(f'{v:.2f}' for v in dk)
    djs = ' '.join(f'{v:.2f}' for v in dj)

    log(f"  {len(sums):>3}  {wk:>7.3f}  {w_pi:>7.3f}  {wk + w_pi:>8.5f}  "
        f"[{dks}]  [{djs}]")

max_dev = max(abs(s - N) for s in sums)
log(f"\n  Average w_k + w(Pi|r_k>): {np.mean(sums):.10f} (should be {N})")
log(f"  Max |sum - N|: {max_dev:.2e}")
log(f"  XY-weight swap: {'VERIFIED' if max_dev < 1e-8 else 'FAILED'}")

# Summary by dominant weight sector
log(f"\n  Dominant weight sector pairing:")
dom_pairs = {}
seen = set()
for k in range(num):
    j = pair_map[k]
    if j < 0 or k in seen:
        continue
    seen.add(k)
    seen.add(j)
    dk = wt_dist(R[:, k], all_idx, N)
    dj = wt_dist(R[:, j], all_idx, N)
    dom_k = int(np.argmax(dk))
    dom_j = int(np.argmax(dj))
    key = (min(dom_k, dom_j), max(dom_k, dom_j))
    dom_pairs[key] = dom_pairs.get(key, 0) + 1

for (w1, w2), count in sorted(dom_pairs.items()):
    log(f"    w={w1} <-> w={w2}: {count} pairs  "
        f"(sum={w1 + w2}, {'= N' if w1 + w2 == N else '!= N'})")


# ############################################################
# SECTION 4: Bell oscillating modes
# ############################################################
log()
log("=" * 90)
log("SECTION 4: Bell(0,1) oscillating modes")
log("  Pi maps forward (Re(mu) > 0) to backward (Re(mu) < 0)")
log("  Forward = slow decay = {I,Z}-heavy; Backward = fast decay = {X,Y}-heavy")
log("=" * 90)

# Bell(0,1) = (|000> + |110>)/sqrt(2) -> Pauli vector
psi = np.zeros(dim, dtype=complex)
psi[0] = 1 / np.sqrt(2)    # |000>
psi[6] = 1 / np.sqrt(2)    # |110>
rho = np.outer(psi, psi.conj())
v_bell = np.array([np.trace(P @ rho) for P in pmats])

# Eigenbasis decomposition
c = np.linalg.solve(R, v_bell)
recon_err = np.linalg.norm(R @ c - v_bell) / np.linalg.norm(v_bell)
log(f"\n  Bell(0,1) = (|000> + |110>)/sqrt(2)")
log(f"  Eigenbasis reconstruction error: {recon_err:.2e}")

# Show oscillating pairs by frequency band
for bname, lo, hi in [("omega ~ 2", 1.5, 2.5),
                       ("omega ~ 4", 3.5, 4.5),
                       ("omega ~ 6", 5.5, 6.5)]:
    # Collect palindromic pairs in this band excited by Bell
    band_pairs = []
    seen2 = set()
    for k in range(num):
        j = pair_map[k]
        if j < 0 or k in seen2:
            continue
        omega = abs(np.imag(mu[k]))
        if not (lo < omega < hi):
            continue
        if max(abs(c[k]), abs(c[j])) < 0.001:
            continue
        seen2.add(k)
        seen2.add(j)
        band_pairs.append((k, j))

    if not band_pairs:
        log(f"\n  {bname}: no excited pairs")
        continue

    log(f"\n  {bname}: {len(band_pairs)} excited palindromic pairs")

    # Sort by total weight
    band_pairs.sort(key=lambda kj: -(abs(c[kj[0]])**2 + abs(c[kj[1]])**2))

    for pi_idx, (k, j) in enumerate(band_pairs[:3]):
        # Identify forward (Re(mu) > 0) and backward (Re(mu) < 0)
        if np.real(mu[k]) > np.real(mu[j]):
            fwd, bwd = k, j
        else:
            fwd, bwd = j, k

        log(f"\n    Pair {pi_idx + 1}:")
        log(f"      Forward  (k={fwd:>2}): mu = {fmt(mu[fwd])}, "
            f"Re(mu) = {np.real(mu[fwd]):+.6f}, |c| = {abs(c[fwd]):.6f}")
        log(f"      Backward (k={bwd:>2}): mu = {fmt(mu[bwd])}, "
            f"Re(mu) = {np.real(mu[bwd]):+.6f}, |c| = {abs(c[bwd]):.6f}")

        # XY-weight analysis
        wf = avg_xy_wt(R[:, fwd], all_idx)
        wb = avg_xy_wt(R[:, bwd], all_idx)
        log(f"      XY-weight: forward = {wf:.3f}, backward = {wb:.3f}, "
            f"sum = {wf + wb:.4f} (expect {N})")

        # Weight distribution
        df = wt_dist(R[:, fwd], all_idx, N)
        db = wt_dist(R[:, bwd], all_idx, N)
        df_str = '  '.join(f'w{w}={v:.3f}' for w, v in enumerate(df) if v > 0.005)
        db_str = '  '.join(f'w{w}={v:.3f}' for w, v in enumerate(db) if v > 0.005)
        log(f"      Forward  weight dist: {df_str}")
        log(f"      Backward weight dist: {db_str}")

        # Dominant Pauli strings for each
        for label, ki in [("Forward ", fwd), ("Backward", bwd)]:
            v_mode = R[:, ki]
            va = np.abs(v_mode)
            order = np.argsort(-va)
            tops = []
            for idx in order[:6]:
                if va[idx] < 0.05 * va[order[0]]:
                    break
                w = xy_weight(all_idx[idx])
                iz_type = all(all_idx[idx][s] in (0, 3) for s in range(N))
                xy_type = all(all_idx[idx][s] in (1, 2) for s in range(N))
                tag = ""
                if iz_type:
                    tag = " {IZ}"
                elif xy_type:
                    tag = " {XY}"
                tops.append(f"{plabel(all_idx[idx])}(w={w}{tag})")
            log(f"      {label} dominant: {', '.join(tops)}")

    # Summary for this frequency band
    all_fwd = []
    all_bwd = []
    for k, j in band_pairs:
        if np.real(mu[k]) > np.real(mu[j]):
            all_fwd.append(k)
            all_bwd.append(j)
        else:
            all_fwd.append(j)
            all_bwd.append(k)

    if all_fwd and all_bwd:
        avg_wf = np.mean([avg_xy_wt(R[:, k], all_idx) for k in all_fwd])
        avg_wb = np.mean([avg_xy_wt(R[:, k], all_idx) for k in all_bwd])
        avg_re_f = np.mean([np.real(mu[k]) for k in all_fwd])
        avg_re_b = np.mean([np.real(mu[k]) for k in all_bwd])
        log(f"\n    Band summary ({len(band_pairs)} pairs):")
        log(f"      Forward:  avg Re(mu) = {avg_re_f:+.6f}, avg w_xy = {avg_wf:.3f}")
        log(f"      Backward: avg Re(mu) = {avg_re_b:+.6f}, avg w_xy = {avg_wb:.3f}")
        log(f"      Forward is {'low-w (classical)' if avg_wf < avg_wb else 'high-w (quantum)'}")
        log(f"      Backward is {'high-w (quantum)' if avg_wf < avg_wb else 'low-w (classical)'}")


# ############################################################
# SECTION 5: ZZZ node demonstration
# ############################################################
log()
log("=" * 90)
log("SECTION 5: ZZZ -> YYY")
log("=" * 90)

# Direct computation: Pi applied to ZZZ basis vector
zzz = (3, 3, 3)
yyy = (2, 2, 2)
zzz_i = all_idx.index(zzz)
yyy_i = all_idx.index(yyy)

e_zzz = np.zeros(num, dtype=complex)
e_zzz[zzz_i] = 1.0
Pi_zzz = Pi @ e_zzz

log(f"\n  Pi applied to ZZZ (Pauli basis index {zzz_i}):")
for a in range(num):
    if abs(Pi_zzz[a]) > 1e-10:
        log(f"    Pi(ZZZ) = {fmtp(Pi_zzz[a])} * {plabel(all_idx[a])}  "
            f"(index {a})")

# Verify
expected = -1j
actual = Pi_zzz[yyy_i]
log(f"\n  Expected: Pi(ZZZ) = -i * YYY")
log(f"  Actual:   Pi(ZZZ) = {fmtp(actual)} * YYY")
log(f"  Match: {'EXACT' if abs(actual - expected) < 1e-10 else 'MISMATCH'}")

# All other components should be zero
other_max = max(abs(Pi_zzz[a]) for a in range(num) if a != yyy_i)
log(f"  Other components: max = {other_max:.2e} "
    f"({'zero' if other_max < 1e-10 else 'NONZERO'})")

# Decay rates
rate_zzz = 2 * gamma * xy_weight(zzz)
rate_yyy = 2 * gamma * xy_weight(yyy)
log(f"\n  Decay rates:")
log(f"    ZZZ: XY-weight = {xy_weight(zzz)}, "
    f"rate = {rate_zzz:.6f}")
log(f"    YYY: XY-weight = {xy_weight(yyy)}, "
    f"rate = {rate_yyy:.6f} = 2*Sg = {2 * Sg:.6f}")

if rate_zzz == 0:
    log(f"    Ratio: infinity (ZZZ is immortal, YYY decays at max rate)")
else:
    log(f"    Ratio: {rate_yyy / rate_zzz:.1f}")

log(f"\n  Physical interpretation:")
log(f"    ZZZ: zero dephasing rate -> most persistent -> purely classical")
log(f"    YYY: maximum dephasing (2*Sg={2*Sg:.2f}) -> most fragile -> purely quantum")
log(f"    Pi(past) = future: the most decided maps to the most undecided")

# Extended Pauli string demonstrations
log(f"\n  Extended Pi mappings:")
log(f"  {'Input':>6}  {'w':>2}  {'rate':>6}  {'-->':>5}  "
    f"{'Output':>8}  {'w':>2}  {'rate':>6}  {'sum_w':>5}")
log(f"  {'-'*50}")

demo_strings = [
    (3, 3, 3), (0, 0, 0), (3, 3, 0), (3, 0, 3), (0, 3, 3),
    (1, 1, 0), (2, 2, 0), (1, 1, 3), (2, 2, 3),
    (1, 0, 1), (2, 0, 2), (1, 2, 0), (0, 1, 2),
]

for idx in demo_strings:
    li = plabel(idx)
    wi = xy_weight(idx)
    ri = 2 * gamma * wi
    mp = tuple(PI_PERM[i] for i in idx)
    ph = 1
    for i in idx:
        ph *= PI_SIGN[i]
    lo2 = plabel(mp)
    wo = xy_weight(mp)
    ro = 2 * gamma * wo
    log(f"  {li:>6}  {wi:>2}  {ri:>6.3f}  {fmtp(ph):>5}  "
        f"{lo2:>8}  {wo:>2}  {ro:>6.3f}  {wi + wo:>5}")


# ############################################################
# DONE
# ############################################################
log()
log("=" * 90)
log("VERIFICATION COMPLETE")
log(f"Date: {datetime.now()}")
log("=" * 90)
f.close()
print(f"\n>>> Results written to {OUT}")
