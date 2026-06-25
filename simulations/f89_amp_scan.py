"""F89 path-3 octic amplitude ('loudness') across q (zeros_connecting_structure arc, piece 3).

Result (death = transition): the physical biorthogonal amplitude a_k = (O.R_k)(L_k.rho) DIVERGES at the
defective EPs (eigenvectors coalesce, R singular, Petermann amplification) = LOUD transitions, but stays
smooth at the diabolic (semisimple, eigenvectors stay independent) = the SILENT transition. Confirmed: a
real-q scan spikes near the genuine EPs (q~0.857 ~10x, ~1.74) and passes smoothly through the diabolic
(q=0.659, like generic). Both are degeneracies (eigenvalue gap -> 0) but only the defective one is loud:
defective vs diabolic seen in the loudness/time layer.

Usage: python simulations/f89_amp_scan.py [out.png]   (default out: f89_amp_scan.png in the cwd).
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "f89_amp_scan.png"


def build(J, gamma):
    de = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    basis = [(i, jk) for i in range(4) for jk in de]
    M_SE = np.zeros((4, 4))
    for a in range(4):
        for b in range(4):
            if abs(a - b) == 1:
                M_SE[a, b] = 2 * J
    M_DE = np.zeros((6, 6))
    for idx, (j, k) in enumerate(de):
        for nj in (j - 1, j + 1):
            if 0 <= nj <= 3 and nj != k:
                p = tuple(sorted([nj, k]))
                if abs(nj - j) == 1 and p in de:
                    M_DE[de.index(p), idx] += 2 * J
        for nk in (k - 1, k + 1):
            if 0 <= nk <= 3 and nk != j:
                p = tuple(sorted([j, nk]))
                if abs(nk - k) == 1 and p in de:
                    M_DE[de.index(p), idx] += 2 * J
    L = np.zeros((len(basis), len(basis)), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(4):
            if M_SE[i2, i] != 0:
                L[basis.index((i2, jk)), idx] += -1j * M_SE[i2, i]
        ji = de.index(jk)
        for j2 in range(6):
            if M_DE[ji, j2] != 0:
                L[basis.index((i, de[j2])), idx] += 1j * M_DE[ji, j2]
        L[idx, idx] += -2 * gamma if i in jk else -6 * gamma
    return L, basis, de


def s2(basis):
    perm = {0: 3, 1: 2, 2: 1, 3: 0}
    pde = lambda jk: tuple(sorted([perm[jk[0]], perm[jk[1]]]))
    cols, hd = [], set()
    for idx, (i, jk) in enumerate(basis):
        if idx in hd:
            continue
        idx2 = basis.index((perm[i], pde(jk)))
        v = np.zeros(len(basis), dtype=complex)
        if idx == idx2:
            v[idx] = 1.0
        else:
            v[idx] = v[idx2] = 1 / np.sqrt(2)
            hd.add(idx2)
        cols.append(v)
        hd.add(idx)
    return np.column_stack(cols)


def wsite(basis, nb=4):
    w = np.zeros((nb, len(basis)))
    for idx, (i, jk) in enumerate(basis):
        for l in jk:
            if l == i:
                continue
            other = jk[1] if jk[0] == l else jk[0]
            if other == i:
                w[l, idx] = 1.0
    return w


gamma = 1.0
L0, basis, de = build(1.0, gamma)        # structure only; rebuild per q
P = s2(basis)
rho_s = P.conj().T @ np.ones(len(basis), dtype=complex)
o_s = wsite(basis).sum(axis=0) @ P

qs = np.linspace(0.40, 3.0, 2600)
maxa, suma, gapmin = [], [], []
for q in qs:
    L, _, _ = build(q * gamma, gamma)
    Ls = P.conj().T @ L @ P
    lam, R = np.linalg.eig(Ls)
    Rinv = np.linalg.inv(R)
    a = (o_s @ R) * (Rinv @ rho_s)
    rate = -lam.real
    octic = ~((np.abs(rate - 2) < 1e-6) | (np.abs(rate - 6) < 1e-6))
    aa = np.abs(a[octic])
    maxa.append(aa.max())
    suma.append(aa.sum())
    # min eigenvalue gap among octic (cusp -> 0 near an EP)
    lo = lam[octic]
    g = min(abs(lo[i] - lo[j]) for i in range(len(lo)) for j in range(i + 1, len(lo)))
    gapmin.append(g)

maxa, suma, gapmin = np.array(maxa), np.array(suma), np.array(gapmin)

fig, ax = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
fig.patch.set_facecolor("#05060a")
for a_ in ax:
    a_.set_facecolor("#05060a"); a_.tick_params(colors="#39ff77")
    for s in a_.spines.values():
        s.set_color("#39ff77")
    a_.grid(color="#0c2a16")

ax[0].semilogy(qs, maxa, color="#39ff14", lw=1.8, label="max |a_k|  (loudest octic mode)")
ax[0].semilogy(qs, suma, color="#08f7fe", lw=1.2, label="sum |a_k|")
ax[0].axvline(0.659, color="#ffae00", ls="--", lw=1.5, label="diabolic q_EP=0.659 (silent transition)")
for qe, lab in [(0.857, "EPs (loud transitions)"), (1.74, None)]:
    ax[0].axvline(qe, color="#ff2bd6", ls=":", lw=1.5, label=lab)
ax[0].set_ylabel("amplitude (loudness)", color="#39ff77")
ax[0].set_title("F89 path-3 octic 'loudness' a_k across q: do the amplitudes spike at the transitions?", color="#7ef9ff")
ax[0].legend(facecolor="#0c1018", labelcolor="#eaffea", edgecolor="#39ff77", fontsize=9)

ax[1].semilogy(qs, gapmin, color="#fe53bb", lw=1.5, label="min octic eigenvalue gap (->0 at an EP)")
ax[1].axvline(0.659, color="#ffae00", ls="--", lw=1.5)
ax[1].axvline(0.857, color="#ff2bd6", ls=":", lw=1.5)
ax[1].axvline(1.74, color="#ff2bd6", ls=":", lw=1.5)
ax[1].set_xlabel("q = J/gamma", color="#39ff77"); ax[1].set_ylabel("min gap", color="#39ff77")
ax[1].legend(facecolor="#0c1018", labelcolor="#eaffea", edgecolor="#39ff77", fontsize=9)

fig.tight_layout(); fig.savefig(OUT, dpi=110, facecolor=fig.get_facecolor())

# report the peaks
def near(qval):
    i = np.argmin(np.abs(qs - qval)); return maxa[i]
print(f"max|a| at q=0.659 (diabolic): {near(0.659):.3e}")
print(f"max|a| at q=0.857 (EP):       {near(0.857):.3e}")
print(f"max|a| at q=1.74  (EP):       {near(1.74):.3e}")
print(f"max|a| at q=2.0   (generic):  {near(2.0):.3e}")
top = np.argsort(maxa)[-6:][::-1]
print("top-6 loudness peaks at q =", [f"{qs[i]:.3f}({maxa[i]:.1e})" for i in top])
print("saved", OUT)
