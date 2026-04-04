"""
Correction: surviving mode energy per site (not populations)
=============================================================
Computes the spatial distribution of trapped light (surviving {X,Y}
coherence energy), the correct observable for "mass = trapped light."

Output: printed to console (results go into TRAPPED_LIGHT_LOCALIZATION.md)
"""

import numpy as np
from pathlib import Path
import sys, os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

J = 1.0; GAMMA = 0.05; EPS = 0.001
I2 = np.eye(2, dtype=complex)
Xm = np.array([[0,1],[1,0]], dtype=complex)
Ym = np.array([[0,-1j],[1j,0]], dtype=complex)
Zm = np.array([[1,0],[0,-1]], dtype=complex)

def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def build_liouvillian(N, gammas):
    d = 2**N; Id = np.eye(d, dtype=complex)
    H = np.zeros((d, d), dtype=complex)
    for i in range(N-1):
        for P in [Xm, Ym, Zm]:
            ops=[I2]*N; ops[i]=P; ops[i+1]=P; H+=J*kron_chain(ops)
    L = -1j*(np.kron(H,Id)-np.kron(Id,H.T))
    for k in range(N):
        ops=[I2]*N; ops[k]=Zm; Lk=np.sqrt(gammas[k])*kron_chain(ops)
        LdL=Lk.conj().T@Lk
        L+=np.kron(Lk,Lk.conj())-0.5*(np.kron(LdL,Id)+np.kron(Id,LdL.T))
    return L

def xy_weight_per_site(eigvec, N):
    """For an eigenvector in |i><j| basis, compute the {X,Y} coherence
    weight localized at each site.

    XY-weight at site q: the eigenvector component contributes to site q
    if the density matrix element (i,j) has i_q != j_q (i.e., the q-th
    bit differs between row and column, meaning an off-diagonal in that qubit).
    """
    d = 2**N
    site_weight = np.zeros(N)
    for idx in range(d*d):
        i = idx // d  # row
        j = idx % d   # col
        amp2 = np.abs(eigvec[idx])**2
        if amp2 < 1e-30:
            continue
        for q in range(N):
            iq = (i >> q) & 1
            jq = (j >> q) & 1
            if iq != jq:  # off-diagonal at site q = X or Y coherence
                site_weight[q] += amp2
    return site_weight

print("=" * 70)
print("CORRECTION: SURVIVING MODE ENERGY PER SITE")
print("=" * 70)
print()

for N in [4, 5]:
    d = 2**N; d2 = d*d

    for profile_name, gammas in [
        ("uniform", [GAMMA]*N),
        ("sacrifice", [N*GAMMA-(N-1)*EPS]+[EPS]*(N-1))
    ]:
        L = build_liouvillian(N, gammas)
        eigvals, R = np.linalg.eig(L)
        _, Lf = np.linalg.eig(L.T)
        ov = Lf.conj().T @ R
        for j in range(d2): Lf[:,j] /= ov[j,j]

        # Initial state |+>^N
        plus = np.array([1,1],dtype=complex)/np.sqrt(2)
        psi = plus
        for _ in range(N-1): psi = np.kron(psi, plus)
        rho0 = np.outer(psi, psi.conj())
        coeffs = Lf.conj().T @ rho0.ravel()

        # Surviving mode energy at t = K_fold / gamma ~ 20
        t_eval = 20.0

        # For each mode: survival factor and XY weight per site
        total_xy_per_site = np.zeros(N)
        total_weight = 0

        for mi in range(d2):
            # Mode properties
            re_lam = eigvals[mi].real
            im_lam = abs(eigvals[mi].imag)

            if im_lam < 1e-6:  # non-oscillating mode, skip
                continue

            # Survival at t_eval
            survival = np.abs(coeffs[mi] * np.exp(eigvals[mi] * t_eval))**2
            if survival < 1e-20:
                continue

            # Energy weight: frequency × survival × initial amplitude
            energy_weight = im_lam * survival

            # XY weight per site from eigenvector
            eigvec = R[:, mi]
            site_xy = xy_weight_per_site(eigvec, N)

            # Normalize site_xy
            total_xy = np.sum(site_xy)
            if total_xy > 1e-20:
                site_xy /= total_xy

            total_xy_per_site += energy_weight * site_xy
            total_weight += energy_weight

        # Normalize
        if total_weight > 0:
            total_xy_per_site /= total_weight

        print(f"N={N}, {profile_name}: surviving mode energy per site (t={t_eval})")
        print(f"  gamma = [{', '.join(f'{g:.4f}' for g in gammas)}]")
        print(f"  XY energy: [{', '.join(f'{w:.4f}' for w in total_xy_per_site)}]")

        # Is it center-localized?
        if N >= 4:
            center = total_xy_per_site[N//2-1:N//2+1].mean()
            edge = (total_xy_per_site[0] + total_xy_per_site[-1]) / 2
            ratio = center / edge if edge > 0 else 0
            print(f"  Center/edge ratio: {ratio:.3f} ({'center-localized' if ratio > 1.1 else 'uniform'})")
        print()

print("Reference from CAVITY_MODE_LOCALIZATION (March 30):")
print("  Sacrifice zone Q0 weight: [0.519, 0.631, 0.700, 0.631, 0.519]")
print("  r = 0.994 correlation between edge weight and decay rate")
