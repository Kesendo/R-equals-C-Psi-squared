# Zero-QPU audit of the A-vs-B concentrator test's created-MI observable.
# Evidence for Downgrade 1 of experiments/CONCENTRATOR_AB_MECHANISM_TEST.md
# (the reckoning section).
#
# Question: is the "created interior MI" the A-vs-B concentrator test measures a
# genuine concentrator effect, or a classical-mixing artifact of the frozen K=16
# RZ-phase sink construction?
#
# Method (noiseless, exact density-matrix sim of the flown circuit, N=5):
#   (a) struct-K16 : the ACTUAL flown 16 frozen phase paths, pooled  -> Sum-MI
#   (b) true channel: per-step edge phase-damping at the same per-step retention
#                     r (the K->inf / Markovian dephasing the sink is meant to be)
#   (c) iid-K      : independent-per-step +/-a phases (cos a = r_s), pooled,
#                    K = 16..1024 -> must converge to (b) if the channel ref is right
#
# Result (run 2026-07-05, seed 20260705), all three flown doses; the artifact
# fraction is DOSE-DEPENDENT:
#   phi=0 sanity -> Sum-MI = 0 (fixed-point theorem, convention correct)
#   run-3 dose S=0.378105 (injected gamma_edge = N*gamma_0):
#     struct-K16 = 0.319 / 0.398 / 0.491 / 0.592 / 0.505   (t=1..5)
#     true channel = 0.122 / 0.073 / 0.024 / 0.030 / 0.027
#     => classical-mixing floor = 62% / 82% / 95% / 95% / 95% of the signal
#   run-2 dose S=0.270756 (half dose): floor = 56% / 81% / 95% / 96% / 95%
#   run-1 ceiling S=1.0 (r=0): floor = 3.5% / 27% / 20% / 24% / -9%, i.e. the
#     ceiling construction is CHANNEL-DOMINATED (only the per-step marginal is
#     exact scramble at r=0; the frozen multi-step correlation survives, hence
#     the 20-27% mid-t floor and the unexplained -9% at t=5); the artifact
#     owns the partial doses only. "floor = struct - channel" is an additive
#     split of a nonlinear functional (MI), so the fractions are heuristic.
#   iid-K converges to the channel (within noise of it at t=1, ~3x above at
#   t=5 for K=16, 21-29% above at K=64, tight by K>=256), so the partial-dose
#   artifact is the frozen intra-path correlation of the STRUCTURE, NOT
#   finite-K (more shots/K does not fix it). At the partial doses the flown
#   "created MI" is dominated by classical mixing; note also that the
#   channel-only signal clears the binding band (+-0.0486) at t=1-2 only, so
#   the late-t persistence that fired the run-2/3 rule is not reachable by the
#   intended channel at those doses.
#
# Run:  python simulations/ab_classical_mixing_check.py

import sys, numpy as np, json
sys.path.insert(0, r"D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography")
import run_ab_test as R
from scipy.linalg import expm

N = 5
theta = 2 * R.J_COUPLING * R.DT          # = 1.0
SEED = 20260705
DOSE = 0.378105                          # run-3 corrected dose (r*=e^-0.25); iid + S2 sections
DOSES = [("run-1 ceiling S=1.0", 1.0),
         ("run-2 half dose S=0.270756", 0.270756),
         ("run-3 gamma_0-anchored dose S=0.378105", 0.378105)]
HW = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\data\ibm_ab_test_july2026\ab_test_hardware_20260705_203200.json"

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)

def op(s, q):
    m = 1
    for i in range(N):
        m = np.kron(m, s if i == q else I2)
    return m

def op2(a, b, qa, qb):
    m = 1
    for i in range(N):
        m = np.kron(m, a if i == qa else (b if i == qb else I2))
    return m

U_H = np.eye(2 ** N, dtype=complex)      # one Trotter step of the Heisenberg bonds, b=0..3
for b in range(N - 1):
    U_H = expm(-1j * (theta / 2) * (op2(X, X, b, b + 1) + op2(Y, Y, b, b + 1) + op2(Z, Z, b, b + 1))) @ U_H
Z0 = op(Z, 0)

def RZ0(phi):
    return op(np.diag([np.exp(-1j * phi / 2), np.exp(1j * phi / 2)]).astype(complex), 0)

plus = np.array([1, 1], complex) / np.sqrt(2)
psi0 = 1
for _ in range(N):
    psi0 = np.kron(psi0, plus)
psi0 = psi0.astype(complex)
bonds = [(i, i + 1) for i in range(N - 1)]

def r1(rho, q):
    t = rho.reshape([2] * N + [2] * N); L = 'abcdefghij'; br = list(L[:N]); ke = list(L[N:2 * N])
    for o in range(N):
        if o != q: ke[o] = br[o]
    return np.einsum(''.join(br) + ''.join(ke) + '->' + br[q] + ke[q], t).reshape(2, 2)

def r2(rho, i, j):
    t = rho.reshape([2] * N + [2] * N); L = 'abcdefghij'; br = list(L[:N]); ke = list(L[N:2 * N])
    for o in range(N):
        if o not in (i, j): ke[o] = br[o]
    return np.einsum(''.join(br) + ''.join(ke) + '->' + br[i] + br[j] + ke[i] + ke[j], t).reshape(4, 4)

def vn(r):
    w = np.linalg.eigvalsh(r); w = w[w > 1e-13]
    return float(-np.sum(w * np.log2(w)))

def sumMI(rho):
    return sum(vn(r1(rho, i)) + vn(r1(rho, j)) - vn(r2(rho, i, j)) for (i, j) in bonds)

def mix(phases):                          # pooled mixture of unitary trajectories
    K, nt = phases.shape; rho = np.zeros((2 ** N, 2 ** N), complex)
    for k in range(K):
        psi = psi0.copy()
        for s in range(nt):
            psi = RZ0(phases[k, s]) @ (U_H @ psi)
        rho += np.outer(psi, psi.conj())
    return sumMI(rho / K)

def chan(nt, ret):                        # true per-step edge phase-damping channel
    rho = np.outer(psi0, psi0.conj())
    for s in range(nt):
        rho = U_H @ rho @ U_H.conj().T
        r = ret[s]; rho = 0.5 * (1 + r) * rho + 0.5 * (1 - r) * (Z0 @ rho @ Z0)
    return sumMI(rho)

_, meta = R.freeze_sink_phases(SEED, DOSE)
paths_by_t, _ = R.freeze_sink_phases(SEED, DOSE)

print("phi=0 sanity (fixed-point theorem, must be ~0):")
for ti, nt in enumerate(R.TROTTER_STEPS):
    psi = psi0.copy()
    for _ in range(nt):
        psi = U_H @ psi
    print(f"  t={ti+1}: {sumMI(np.outer(psi, psi.conj())):.4f}")

for dose_label, dose_s in DOSES:
    d_paths, d_meta = R.freeze_sink_phases(SEED, dose_s)
    print(f"\nS1 [{dose_label}] -- struct-K16 (flown) vs true channel; classical-mixing floor:")
    print(f"{'t':>3}{'struct-K16':>12}{'channel':>10}{'floor':>9}{'floor/K16':>11}")
    for ti, nt in enumerate(R.TROTTER_STEPS):
        ret = [row['r'] for row in d_meta[ti]['retention']]
        mk = mix(d_paths[ti]); ch = chan(nt, ret); fl = mk - ch
        print(f"{ti+1:>3}{mk:>12.4f}{ch:>10.4f}{fl:>9.4f}{fl/mk if abs(mk)>1e-9 else 0:>11.3f}")

print("\niid-K convergence (binary +/-a, cos a = r_s) -> validates the channel ref:")
rng = np.random.default_rng(7)
for ti in (0, 2, 4):
    nt = R.TROTTER_STEPS[ti]; ret = np.array([row['r'] for row in meta[ti]['retention']])
    a = np.arccos(np.clip(ret, -1, 1)); struct = mix(paths_by_t[ti]); ch = chan(nt, ret)
    row = []
    for K in (16, 64, 256, 1024):
        signs = rng.choice([-1, 1], size=(K, nt)); row.append(mix((signs * a[None, :]) % (2 * np.pi)))
    print(f"  t={ti+1}: struct-K16={struct:.4f} | iid K=16/64/256/1024=" +
          "/".join(f'{v:.4f}' for v in row) + f" | channel={ch:.4f}")

print("\nS2 -- hardware-native null (2.33*bootstrap-SE, shot-noise-only) vs sim band 0.0486, both legs:")
d = json.load(open(HW))
for leg in ("Delta", "Delta_u"):
    pd = d["primary_diffs"][leg]
    for ti in range(5):
        lo, hi = pd["delta_ci95"][ti]; se = (hi - lo) / (2 * 1.96); hw = 2.33 * se
        m = abs(pd['delta'][ti]) / hw
        print(f"  {leg:>7} t={ti+1}: {pd['delta'][ti]:+.4f}  hw-null~{hw:.4f}  margin={m:.2f}x  clears={'yes' if m>1 else 'NO'}")
