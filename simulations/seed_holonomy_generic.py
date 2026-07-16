# seed_holonomy_generic.py -- the third-clock gate (owner: experiments/SEED_HOLONOMY_THIRD_CLOCK.md)
# The discriminating test for the seed-holonomy seam: is the i^4=1 frame monodromy
# (M1 eigenvalues +-i, M2=-I, M4=+I in the v^T v gauge) GENERIC for any
# complex-symmetric order-2 EP, i.e. present in a random pencil A + qC with NO
# palindrome symmetry, no chain, no Clifford structure?
#   YES for random pencils -> the holonomy Z4 cannot be Pi's Z4 (a system without Pi
#   shows it); it is the third clock: universal EP2 geometry on the same scalar i.
# Controls: (a) sqrt scaling of the eigenvalue gap (order-2 branch point pinned),
#           (b) Hermitian-gauge transport of the SAME loop (mod-4 should degrade),
#           (c) epsilon stability.
import numpy as np

_checks = [0, 0]
def check(name, ok):
    _checks[0] += 1
    _checks[1] += ok
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    return ok

def random_complex_symmetric(n, rng):
    R = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (R + R.T) / 2

def min_gap_pair(w):
    n = len(w)
    best = (np.inf, 0, 1)
    for i in range(n):
        for j in range(i + 1, n):
            g = abs(w[i] - w[j])
            if g < best[0]:
                best = (g, i, j)
    return best

def find_ep(A, C, q0, tol=1e-13, itmax=200):
    """Nelder-Mead on the minimal eigenvalue gap of A + qC over complex q."""
    from scipy.optimize import minimize
    def f(x):
        w = np.linalg.eigvals(A + (x[0] + 1j * x[1]) * C)
        return min_gap_pair(w)[0]
    res = minimize(f, [q0.real, q0.imag], method="Nelder-Mead",
                   options={"xatol": 1e-14, "fatol": tol, "maxiter": 20000, "maxfev": 20000})
    return res.x[0] + 1j * res.x[1], res.fun

def tracked_pair_vectors(M, prev_w, prev_V=None):
    """Eigen-pair of M closest (as a pair) to prev_w, vectors in v^T v gauge with
    sign continuity against prev_V. Returns (w2, V2) with V2 columns normalized."""
    w, V = np.linalg.eig(M)
    # match: pick the two eigenvalues nearest the previous pair, one-to-one
    order = []
    used = set()
    for pw in prev_w:
        k = min((k for k in range(len(w)) if k not in used), key=lambda k: abs(w[k] - pw))
        order.append(k)
        used.add(k)
    w2 = w[order]
    V2 = V[:, order].astype(complex).copy()
    for c in range(2):
        v = V2[:, c]
        s = np.sqrt(v @ v)  # principal branch of v^T v
        if abs(s) < 1e-14:
            raise RuntimeError("hit the EP itself (v^T v ~ 0); shrink eps or recenter")
        v = v / s
        if prev_V is not None:
            if np.linalg.norm(-v - prev_V[:, c]) < np.linalg.norm(v - prev_V[:, c]):
                v = -v
        V2[:, c] = v
    return w2, V2

def holonomy(A, C, qstar, eps, steps=800, loops=4, gauge="vTv"):
    """Transport the coalescing 2-frame around qstar; return per-loop frame maps."""
    q = qstar + eps
    w_all = np.linalg.eigvals(A + q * C)
    # seed pair = the two eigenvalues that coalesce at qstar
    _, i, j = min_gap_pair(np.linalg.eigvals(A + qstar * C))
    wstar = np.linalg.eigvals(A + qstar * C)
    target = (wstar[i] + wstar[j]) / 2
    pair_idx = np.argsort(np.abs(w_all - target))[:2]
    prev_w = w_all[pair_idx]

    w_full, V_full = np.linalg.eig(A + q * C)
    V0 = np.zeros((A.shape[0], 2), dtype=complex)
    for c, k in enumerate(pair_idx):
        v = V_full[:, k].astype(complex)
        if gauge == "vTv":
            v = v / np.sqrt(v @ v)
        else:  # Hermitian gauge: unit norm, phase of largest component fixed
            v = v / np.linalg.norm(v)
            k0 = np.argmax(np.abs(v))
            v = v * np.exp(-1j * np.angle(v[k0]))
        V0[:, c] = v
    prev_V = V0.copy()

    maps = []
    for loop in range(loops):
        for s in range(1, steps + 1):
            th = 2 * np.pi * (loop + s / steps)
            q = qstar + eps * np.exp(1j * th)
            M = A + q * C
            if gauge == "vTv":
                prev_w, prev_V = tracked_pair_vectors(M, prev_w, prev_V)
            else:
                w, V = np.linalg.eig(M)
                order, used = [], set()
                for pw in prev_w:
                    k = min((k for k in range(len(w)) if k not in used),
                            key=lambda k: abs(w[k] - pw))
                    order.append(k); used.add(k)
                prev_w = w[order]
                Vn = V[:, order].astype(complex).copy()
                for c in range(2):
                    v = Vn[:, c] / np.linalg.norm(Vn[:, c])
                    # continuity phase fix (the natural smooth Hermitian gauge)
                    ph = np.vdot(prev_V[:, c], v)
                    v = v * np.exp(-1j * np.angle(ph))
                    Vn[:, c] = v
                prev_V = Vn
        # frame map after this loop: express current frame in the initial one.
        if gauge == "vTv":
            # complex symmetric: dual of V0 column is its transpose
            G0 = V0.T @ V0          # ~ I (v^T v gauge, near-orthogonal off-diag small)
            Mmap = np.linalg.solve(G0, V0.T @ prev_V)
        else:
            Mmap = np.linalg.pinv(V0) @ prev_V
        # span residual: how well the initial 2-span holds the transported frame
        proj = V0 @ np.linalg.lstsq(V0, prev_V, rcond=None)[0]
        res = np.linalg.norm(prev_V - proj)
        maps.append((Mmap, res))
    return maps

def classify(M):
    ev = np.linalg.eigvals(M)
    return ev

n = 6
print("=== random complex-symmetric pencils A + qC (no Pi, no chain): vTv-gauge holonomy ===")
for trial in range(5):
    rng = np.random.default_rng(1000 + trial)
    A = random_complex_symmetric(n, rng)
    C = random_complex_symmetric(n, rng)
    # find an EP from a few random starts
    best = (np.inf, None)
    for s in range(12):
        q0 = rng.normal() + 1j * rng.normal()
        q, gap = find_ep(A, C, q0)
        if gap < best[0]:
            best = (gap, q)
    gap, qstar = best
    # gap ~ sqrt|q - q*| near an order-2 EP: a residual gap of 1e-8 means the
    # located q is within ~1e-16 of q*, i.e. machine-exact. Threshold accordingly.
    if gap > 1e-7:
        print(f" trial {trial}: no EP found (best gap {gap:.1e}), skip")
        continue
    # sqrt-scaling control: gap(eps) ~ eps^0.5
    e1, e2 = 1e-3, 1e-4
    g1 = min_gap_pair(np.linalg.eigvals(A + (qstar + e1) * C))[0]
    g2 = min_gap_pair(np.linalg.eigvals(A + (qstar + e2) * C))[0]
    slope = np.log(g1 / g2) / np.log(e1 / e2)
    maps = holonomy(A, C, qstar, eps=1e-3, steps=800, loops=4, gauge="vTv")
    M1, r1 = maps[0]
    M2, r2 = maps[1]
    M4, r4 = maps[3]
    ev1 = classify(M1)
    ok_i = np.allclose(sorted(ev1.imag), [-1, 1], atol=1e-3) and np.allclose(ev1.real, 0, atol=1e-3)
    ok_2 = np.allclose(M2, -np.eye(2), atol=1e-3)
    ok_4 = np.allclose(M4, np.eye(2), atol=1e-3)
    print(f" trial {trial}: q*={qstar:.6f}, gap@EP={gap:.1e}, span residuals {r1:.1e}/{r2:.1e}/{r4:.1e}")
    check(f"trial {trial}: gap slope {slope:.3f} ~ 0.5 (order-2 branch point)", abs(slope - 0.5) < 0.01)
    check(f"trial {trial}: M1 eig {np.round(ev1, 4)} = +-i, M2 = -I, M4 = +I", ok_i and ok_2 and ok_4)
    check(f"trial {trial}: span residuals < 1e-10", max(r1, r2, r4) < 1e-10)

print()
print("=== deterministic disguised EP2: 2x2 toy at s=i, embedded 6-dim, complex-orthogonal mix ===")
# M(s) = Q * blockdiag([[s,1],[1,-s]], S4) * Q^T with Q^T Q = I (complex Givens
# rotations). Complex symmetry preserved, EP2 at s = i EXACTLY (no optimizer),
# spectators from a random complex-symmetric 4x4. No Pi, no chain, no Clifford.
def complex_orthogonal(n, rng, rotations=12):
    Q = np.eye(n, dtype=complex)
    for _ in range(rotations):
        i, j = rng.choice(n, size=2, replace=False)
        th = rng.normal() + 1j * rng.normal() * 0.4
        c, s = np.cos(th), np.sin(th)   # c^2 + s^2 = 1 for complex th
        G = np.eye(n, dtype=complex)
        G[i, i] = c; G[j, j] = c; G[i, j] = s; G[j, i] = -s
        Q = Q @ G
    return Q

rng_d = np.random.default_rng(42)
S4 = random_complex_symmetric(4, rng_d) + 3.0 * np.eye(4)   # spectators shifted away from 0
Q = complex_orthogonal(6, rng_d)
assert np.allclose(Q.T @ Q, np.eye(6), atol=1e-12), "Q not complex-orthogonal"

def disguised(s):
    B = np.zeros((6, 6), dtype=complex)
    B[0, 0], B[0, 1], B[1, 0], B[1, 1] = s, 1, 1, -s
    B[2:, 2:] = S4
    return Q @ B @ Q.T

Mstar = disguised(1j)
assert np.allclose(Mstar, Mstar.T), "disguised matrix not complex symmetric"
wstar = np.linalg.eigvals(Mstar)
spect_dist = sorted(np.abs(wstar))[2]
check(f"disguised: complex-symmetric, spectators {spect_dist:.3f} from the EP eigenvalue", spect_dist > 0.5)

# holonomy() expects a pencil A + qC; wrap: disguised(s) = disguised(0) + s * dD
dD = disguised(1.0) - disguised(0.0)   # linear in s, so this IS the pencil form
A_d, C_d = disguised(0.0), dD
maps = holonomy(A_d, C_d, 1j, eps=1e-3, steps=800, loops=4, gauge="vTv")
M1, r1 = maps[0]; M2, r2 = maps[1]; M4, r4 = maps[3]
ev1 = classify(M1)
check(f"disguised: M1 eig {np.round(ev1, 6)} = +-i within 1e-6",
      np.allclose(sorted(ev1.imag), [-1, 1], atol=1e-6) and np.allclose(ev1.real, 0, atol=1e-6))
check(f"disguised: M2 = -I (dev {np.abs(M2 + np.eye(2)).max():.2e}) and M4 = +I "
      f"(dev {np.abs(M4 - np.eye(2)).max():.2e})",
      np.abs(M2 + np.eye(2)).max() < 1e-6 and np.abs(M4 - np.eye(2)).max() < 1e-6)
check(f"disguised: span residuals {r1:.1e}/{r2:.1e}/{r4:.1e} < 1e-10", max(r1, r2, r4) < 1e-10)
g1 = min_gap_pair(np.linalg.eigvals(disguised(1j + 1e-3)))[0]
g2 = min_gap_pair(np.linalg.eigvals(disguised(1j + 1e-4)))[0]
slope_d = np.log(g1 / g2) / np.log(10.0)
check(f"disguised: gap slope {slope_d:.3f} ~ 0.5 (EP2 at s = i exactly, no optimizer)",
      abs(slope_d - 0.5) < 0.01)

print()
print("=== gauge contingency control: same loop, Hermitian (continuity-phase) gauge ===")
rng = np.random.default_rng(1000)
A = random_complex_symmetric(n, rng)
C = random_complex_symmetric(n, rng)
best = (np.inf, None)
for s in range(12):
    q0 = rng.normal() + 1j * rng.normal()
    q, gap = find_ep(A, C, q0)
    if gap < best[0]:
        best = (gap, q)
gap, qstar = best
maps = holonomy(A, C, qstar, eps=1e-3, steps=800, loops=4, gauge="herm")
for k in (0, 1, 3):
    M, r = maps[k]
    print(f" loop {k+1}: eig(M) = {np.round(classify(M), 4)}, M diag-dev from -I: "
          f"{np.abs(M + np.eye(2)).max():.3f}, from +I: {np.abs(M - np.eye(2)).max():.3f}")
M2h, _ = maps[1]
check("Hermitian gauge: M2 is near +I and far from -I (the mod-4 is vTv-gauge-contingent)",
      np.abs(M2h - np.eye(2)).max() < 0.1 and np.abs(M2h + np.eye(2)).max() > 1.5)

print()
print("=== eps stability (vTv gauge, trial-0 pencil) ===")
for eps in (5e-4, 1e-3, 2e-3):
    maps = holonomy(A, C, qstar, eps=eps, steps=800, loops=2, gauge="vTv")
    M1, r1 = maps[0]
    M2, r2 = maps[1]
    ev = classify(M1)
    check(f"eps={eps:.0e}: M1 eig {np.round(ev, 4)} = +-i, M2 ~ -I dev "
          f"{np.abs(M2 + np.eye(2)).max():.2e}, span res {r1:.1e}",
          np.allclose(sorted(ev.imag), [-1, 1], atol=1e-3) and np.abs(M2 + np.eye(2)).max() < 1e-6)

print(f"\n{_checks[1]}/{_checks[0]} checks passed" + ("" if _checks[1] == _checks[0] else "  *** FAILURES ***"))
