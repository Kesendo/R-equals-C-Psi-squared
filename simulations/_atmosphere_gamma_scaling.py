"""OQ#5: is the cluster |Im|_min a structural number or a gamma0-artefact?

The atmosphere notes flag this: at N=6 the cluster bottoms at ~7e-5 relative to
gamma0=0.05, ratio ~1.4e-3. Is the ratio characteristic (a real F-formula-style
number) or does it scale with gamma0 itself?

Test: vary gamma0 (J fixed) so Q = J/gamma0 changes; for each gamma0 scan eps
around the N=5 dip (~-0.34) at fine grid; track dip depth |Im|_min and dip
location eps*. Two scaling hypotheses:
  (a) |Im|_min / gamma0 = const          -> dimensional artefact
  (b) |Im|_min depends on Q = J/gamma0   -> structural, encoded in Q
  (c) eps* shifts with gamma0/Q          -> the dip is a function of Q, not eps

N=5 is sufficient: same F1-mirror-pair structure, faster than N=6.

Investigation only.
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
N = 5
J = 0.075


def site_op(op, k, N):
    m = np.array([[1.0 + 0j]])
    for i in range(N):
        m = np.kron(m, op if i == k else I2)
    return m


def chain_H(N):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        for P in (X, Y):
            H += site_op(P, b, N) @ site_op(P, b + 1, N)
    return H


def build_L(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(Z, k, N)
        L += gamma[k] * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def sym_shape(N):
    i = np.arange(N, dtype=float)
    u = (i - (N - 1) / 2.0) ** 2
    u = u - u.mean()
    return u / np.max(np.abs(u))


u = sym_shape(N)
H = J * chain_H(N)

# Scan over gamma0 (J fixed): varies Q = J/gamma0
gamma0_list = [0.025, 0.05, 0.075, 0.1, 0.15, 0.2]
EPS = np.linspace(-0.40, -0.25, 16)  # around N=5 dip ~ -0.34

print(f"N={N}, J={J} (fixed)")
print(f"scanning gamma0 in {gamma0_list}, eps grid Δ={(EPS[-1]-EPS[0])/(len(EPS)-1):.4f}")
print(f"\n{'gamma0':>8} {'Q=J/g0':>8} {'eps*':>10} {'|Im|_min':>14} "
      f"{'|Im|_min/g0':>14} {'|Im|_min/J':>14}")
print("-" * 80)

results = []
for g0 in gamma0_list:
    ims = []
    for eps in EPS:
        gamma = g0 * (1.0 + eps * u)
        ev = np.linalg.eigvals(build_L(H, gamma, N))
        aim = np.sort(np.abs(ev.imag))
        nonreal = aim[aim >= 1e-13]
        smin = float(nonreal[0]) if nonreal.size else float('nan')
        ims.append(smin)
    imin = int(np.argmin(ims))
    eps_star = EPS[imin]
    im_star = ims[imin]
    Q = J / g0
    results.append((g0, Q, eps_star, im_star))
    print(f"{g0:>8.4f} {Q:>8.4f} {eps_star:>+10.4f} {im_star:>14.4e} "
          f"{im_star/g0:>14.4e} {im_star/J:>14.4e}")
    sys.stdout.flush()

print("\n--- analysis ---")
print("if |Im|_min/gamma0 is roughly constant -> structural number scales with gamma0")
print("if |Im|_min/J is roughly constant -> structural number scales with J")
print("if eps* shifts -> dip location is Q-dependent")
print()
g0s = np.array([r[0] for r in results])
ims = np.array([r[3] for r in results])
# log-log fit |Im|_min ~ g0^p
slope_g0, intercept = np.polyfit(np.log(g0s), np.log(ims), 1)
print(f"|Im|_min ~ gamma0^{slope_g0:.3f} (fit)")
slope_J = 0  # by construction J is fixed
print(f"|Im|_min at gamma0=0.05 is {ims[1]:.3e}")
print(f"ratio at gamma0=0.05: |Im|_min/gamma0 = {ims[1]/0.05:.3e}")
