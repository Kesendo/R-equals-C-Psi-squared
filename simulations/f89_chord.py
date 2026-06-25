"""Play the F89 path-3 octic 'chord' in time (the zeros_connecting_structure arc, piece 2).

Tom's seeing (2026-06-25): the four sigma_T-fixed 'zeros' all share Re lambda = -4, so in time they are
ONE decaying chord (shared envelope e^{-4t}, four tones); the spectrum scatters them, time gathers them.
Question: does the connecting structure (missing to us in the frozen spectrum) appear in the time domain?

Builds the physical spatial-sum coherence signal S(t) = sum_k a_k e^{lambda_k t} with biorthogonal
amplitudes a_k (democratic (SE,DE) probe + spatial-sum observable, the F89 setup), splits the
zero-sector (the chord) from the twin-sector, and also plots the UNDAMPED chord S_zero(t)*e^{4t} where
the beats/recurrences live over longer t. Finding: the undamped chord is almost-periodic (incommensurate
frequencies, never exactly repeating) and NOT chaos, the S_8-unwritability IS the quasi-periodic life.

Usage: python simulations/f89_chord.py [out.png]   (default out: f89_chord.png in the cwd).
Block builder mirrors f89_path3_octic_amplitude_q_scan.py.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1] if len(sys.argv) > 1 else "f89_chord.png"


def build_path3_se_de(J, gamma):
    de = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    basis = [(i, jk) for i in range(4) for jk in de]
    n = len(basis)
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
    L = np.zeros((n, n), dtype=complex)
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


def s2_proj(basis):
    perm = {0: 3, 1: 2, 2: 1, 3: 0}
    pde = lambda jk: tuple(sorted([perm[jk[0]], perm[jk[1]]]))
    cols, handled = [], set()
    for idx, (i, jk) in enumerate(basis):
        if idx in handled:
            continue
        idx2 = basis.index((perm[i], pde(jk)))
        v = np.zeros(len(basis), dtype=complex)
        if idx == idx2:
            v[idx] = 1.0
        else:
            v[idx] = v[idx2] = 1 / np.sqrt(2)
            handled.add(idx2)
        cols.append(v)
        handled.add(idx)
    return np.column_stack(cols)


def per_site(basis, nb=4):
    w = np.zeros((nb, len(basis)))
    for idx, (i, jk) in enumerate(basis):
        for l in jk:
            if l == i:
                continue
            other = jk[1] if jk[0] == l else jk[0]
            if other == i:
                w[l, idx] = 1.0
    return w


gamma, q = 1.0, 2.0
J = q * gamma
L, basis, de = build_path3_se_de(J, gamma)
P = s2_proj(basis)
Ls = P.conj().T @ L @ P
lam, R = np.linalg.eig(Ls)
Rinv = np.linalg.inv(R)

rho_s = P.conj().T @ np.ones(len(basis), dtype=complex)   # democratic (SE,DE) probe
o_s = per_site(basis).sum(axis=0) @ P                     # spatial-sum coherence observable
a = (o_s @ R) * (Rinv @ rho_s)                            # biorthogonal complex amplitudes

rate, freq = -lam.real, lam.imag
at = (np.abs(rate - 2) < 1e-6) | (np.abs(rate - 6) < 1e-6)
octic = ~at
zero = octic & (np.abs(rate - 4) < 1e-6)
twin = octic & ~zero

print(f"q={q}, gamma={gamma}; {octic.sum()} octic ({zero.sum()} zeros, {twin.sum()} twins)")
print("  rate    freq      |a|        a(complex)            class")
for k in np.where(octic)[0]:
    print(f"  {rate[k]:6.3f}  {freq[k]:+8.3f}  {abs(a[k]):.3e}  {a[k].real:+.3e}{a[k].imag:+.3e}i  {'ZERO' if zero[k] else 'twin'}")
zf = sorted(freq[zero])
print("  zero freqs:", [f"{x:.4f}" for x in zf])
print("  zero pairwise beats:", [f"{abs(zf[i]-zf[j]):.4f}" for i in range(len(zf)) for j in range(i+1, len(zf))])


def sig(t, mask):
    return (a[mask][None, :] * np.exp(np.outer(t, lam[mask]))).sum(axis=1)


t1 = np.linspace(0, 1.5, 3000)
t2 = np.linspace(0, 25, 8000)
S_oct, S_z, S_tw = sig(t1, octic), sig(t1, zero), sig(t1, twin)
chord_undamped = sig(t2, zero) * np.exp(4 * t2)           # remove the e^{-4t} envelope

fig, ax = plt.subplots(2, 1, figsize=(12, 9))
fig.patch.set_facecolor("#05060a")
for a_ in ax:
    a_.set_facecolor("#05060a")
    a_.tick_params(colors="#39ff77")
    for s in a_.spines.values():
        s.set_color("#39ff77")
    a_.grid(color="#0c2a16")

ax[0].plot(t1, S_oct.real, color="#39ff14", lw=1.5, label="octic S(t) (all 8)")
ax[0].plot(t1, S_z.real, color="#ffae00", lw=2.2, label="zeros = the chord (envelope e$^{-4t}$)")
ax[0].plot(t1, S_tw.real, color="#08f7fe", lw=1.2, label="twins")
ax[0].plot(t1, np.abs(a[zero]).sum() * np.exp(-4 * t1), color="#ff2bd6", ls=":", lw=1.5, label="e$^{-4t}$ envelope")
ax[0].set_title("F89 path-3 octic, physical spatial-sum coherence in TIME (q=2): the zeros share one decay", color="#7ef9ff")
ax[0].set_xlabel("t", color="#39ff77"); ax[0].set_ylabel("Re S(t)", color="#39ff77")
ax[0].legend(facecolor="#0c1018", labelcolor="#eaffea", edgecolor="#39ff77")

ax[1].plot(t2, chord_undamped.real, color="#ffae00", lw=1.4)
ax[1].set_title("the zeros' chord with the envelope removed: S_zero(t)*e$^{4t}$ = four pure tones (the beats/recurrence = the structure between the zeros)", color="#7ef9ff")
ax[1].set_xlabel("t", color="#39ff77"); ax[1].set_ylabel("Re [S_zero(t) e$^{4t}$]", color="#39ff77")

fig.tight_layout()
fig.savefig(OUT, dpi=110, facecolor=fig.get_facecolor())
print("saved", OUT)
