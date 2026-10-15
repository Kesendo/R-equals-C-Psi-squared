"""
LIOUVILLIAN DIAGONALIZATION - Generator-side spectral analysis
Original: GPT (March 2026), paths adjusted for local repo.

Diagonalizes the full 3-qubit Liouvillian superoperator,
identifies dominant oscillatory modes, and projects observables
onto eigenmodes. This is the spectral verification test.
"""
import math
import numpy as np
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as st

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_all(*ops):
    out = np.array([[1]], dtype=complex)
    for op in ops:
        out = np.kron(out, op)
    return out


def liouvillian_super(H, L_ops):
    d = H.shape[0]
    I = np.eye(d, dtype=complex)
    L = -1j * (np.kron(I, H) - np.kron(H.T, I))
    for Lop in L_ops:
        LdL = Lop.conj().T @ Lop
        L += np.kron(Lop.conj(), Lop)
        L += -0.5 * np.kron(I, LdL)
        L += -0.5 * np.kron(LdL.T, I)
    return L


def build_biorthogonal_eigensystem(L):
    """Right eigenvectors and their dual (left) basis from R^-1, biorthonormal also
    inside degenerate eigenspaces, so cluster sums in observable_amplitudes are exact."""
    vals_r, VR = np.linalg.eig(L)
    VL = np.linalg.inv(VR).conj().T
    return vals_r, VR, VL


def observable_amplitudes(vals_r, VR, VL, rho0_vec, obs_vec, tol=1e-6):
    """Project an observable onto Liouvillian eigenmodes."""
    clusters = []
    for i, v in enumerate(vals_r):
        placed = False
        for cl in clusters:
            if abs(v - cl['rep']) < tol:
                cl['idx'].append(i)
                placed = True
                break
        if not placed:
            clusters.append({'rep': v, 'idx': [i]})
    rows = []
    for cl in clusters:
        amp = 0.0 + 0.0j
        for i in cl['idx']:
            amp += (obs_vec.conj().T @ VR[:, i]) * (VL[:, i].conj().T @ rho0_vec)
        v = cl['rep']
        rows.append({
            'lambda_real': float(v.real),
            'lambda_imag': float(v.imag),
            'freq': float(v.imag / (2 * math.pi)),
            'decay': float(-v.real),
            'amp': float(abs(amp)),
            'mult': len(cl['idx']),
        })
    rows.sort(key=lambda r: -r['amp'])
    return rows
