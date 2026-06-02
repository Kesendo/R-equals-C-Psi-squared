"""The closed-form local mirror for the two continuous-crossover pairs (XZ+YZ, ZX+ZY).

[Non-Heisenberg Palindrome](../../../experiments/NON_HEISENBERG_PALINDROME.md) once
classified XZ+YZ and ZX+ZY as the only two of 36 two-term bilinear Hamiltonians whose
palindrome mirror Π is "genuinely non-local" (entangled across sites, operator-Schmidt
rank 8-9). That verdict was a lens artifact: it came from restricting the per-site map
to the discrete signed-permutation crossovers {P1, P4, M2}, and from the eigenvector-
pairing construction, which under the heavy spectral degeneracy of these Liouvillians
returns an entangled representative even when a product solution exists.

Widen the per-site map from a discrete permutation to a continuous block rotation and a
product mirror appears. A single uniform per-site unitary M, the SAME on every site,
gives Π = M^⊗N that conjugates the full Liouvillian to machine precision at every N
tested (2..6). The two cases are LOCAL.

The closed form, in the framework Pauli order [I, X, Z, Y] (rows = output, cols = input):

    I ↦ −(X + Y)/√2        X ↦ (I + iZ)/√2
    Z ↦  i(X − Y)/√2       Y ↦ (I − iZ)/√2

M is unitary with M² = −I (an order-4 element, eigenvalues ±i). It swaps the dark bus
{I, Z} and the light bus {X, Y} with a 45° rotation inside each: this is the diplexer
turned to the symmetric point between the X-router (P1) and the Y-router (P4), the one
setting that serves both bands at once. The light combination (X ± Y) echoes the bond
operator itself, since XZ + YZ = (X + Y)·Z.

Public API:
  crossover_map()                              -> the 4×4 closed-form M
  is_crossover_pair(term1, term2)              -> bool
  CROSSOVER_PAIRS                              -> the two frozenset pairs
  product_pi_residual(combo, maps, N, gamma)   -> ‖Π L Π⁻¹ + L + 2Σγ·I‖ for Π = ⊗maps
  verify_crossover_local(combo, N, gamma)      -> residual using M^⊗N (≈ 0 ⇒ local)
"""
from __future__ import annotations

import numpy as np

from ..pauli import _vec_to_pauli_basis_transform, site_op
from ..lindblad import lindbladian_z_dephasing

# Framework single-qubit Pauli order is [I, X, Z, Y] (index a + 2b).
_INV_SQRT2 = 1.0 / np.sqrt(2.0)

# The two pairs whose local mirror is continuous (not a discrete permutation):
# same-site X and Y collide over a shared dark Z, so no single P1/P4 crossover routes
# both bands, but the 45°-rotated crossover does.
CROSSOVER_PAIRS = frozenset({
    frozenset({'XZ', 'YZ'}),
    frozenset({'ZX', 'ZY'}),
})


def crossover_map():
    """The closed-form continuous per-site map M in [I, X, Z, Y] order (M² = −I, unitary).

    Returns the same M for both XZ+YZ and ZX+ZY (the uniform product Π = M^⊗N is
    site-symmetric, so the two site-swapped cases share it).
    """
    s = _INV_SQRT2
    M = np.zeros((4, 4), dtype=complex)
    M[0, 1] = s;  M[0, 3] = s        # output I  ←  X, Y      (i.e. X,Y ↦ +I/√2)
    M[1, 0] = -s; M[1, 2] = 1j * s   # output X  ←  I, Z      (I ↦ −X/√2, Z ↦ iX/√2)
    M[2, 1] = 1j * s; M[2, 3] = -1j * s  # output Z  ←  X, Y
    M[3, 0] = -s; M[3, 2] = -1j * s  # output Y  ←  I, Z
    return M


def is_crossover_pair(term1, term2):
    """True iff {term1, term2} is one of the two continuous-crossover pairs."""
    t1 = term1 if isinstance(term1, str) else ''.join(term1)
    t2 = term2 if isinstance(term2, str) else ''.join(term2)
    return frozenset({t1, t2}) in CROSSOVER_PAIRS


def _liouvillian_pauli(combo, N, gamma):
    """Full Z-dephased Liouvillian for H = Σ_bonds(t1 + t2), in the Pauli-string basis."""
    (a1, b1), (a2, b2) = combo
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        H = H + site_op(N, i, a1) @ site_op(N, i + 1, b1)
        H = H + site_op(N, i, a2) @ site_op(N, i + 1, b2)
    L = lindbladian_z_dephasing(H, [gamma] * N)
    V = _vec_to_pauli_basis_transform(N)
    return (V.conj().T @ L @ V) / d, N * gamma


def _product_pi(maps):
    """Π = maps[0] ⊗ maps[1] ⊗ … in the Pauli-string basis (site 0 most significant)."""
    out = maps[0]
    for M in maps[1:]:
        out = np.kron(out, M)
    return out


def product_pi_residual(combo, maps, N, gamma=0.5):
    """‖Π·L·Π⁻¹ + L + 2Σγ·I‖ for the custom product mirror Π = ⊗ maps.

    Unlike the canonical `palindrome_residual` (which builds the fixed signed-permutation
    Π), this verifies an arbitrary per-site product, the capability needed to show the
    continuous M closes these two cases.

    Args:
        combo: ((a1, b1), (a2, b2)) the two bond bilinears, e.g. (('X','Z'), ('Y','Z')).
        maps:  list of N single-site 4×4 maps (the per-site Π factors).
        N:     chain length.
        gamma: per-site Z-dephasing rate (default 0.5).
    """
    L_pauli, sigma = _liouvillian_pauli(combo, N, gamma)
    Pi = _product_pi(maps)
    resid = Pi @ L_pauli @ np.linalg.inv(Pi) + L_pauli + 2 * sigma * np.eye(L_pauli.shape[0])
    return float(np.linalg.norm(resid))


def verify_crossover_local(combo, N, gamma=0.5):
    """Residual of the uniform closed-form mirror Π = M^⊗N on a crossover combo.

    ≈ 0 (machine precision) confirms the case is LOCAL: a single continuous per-site
    unitary mirrors the full Liouvillian. Raises if `combo` is not a crossover pair.
    """
    (a1, b1), (a2, b2) = combo
    if not is_crossover_pair(a1 + b1, a2 + b2):
        raise ValueError(f"{combo} is not one of the continuous-crossover pairs {CROSSOVER_PAIRS}")
    M = crossover_map()
    return product_pi_residual(combo, [M] * N, N, gamma)
