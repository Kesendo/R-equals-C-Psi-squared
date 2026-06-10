#!/usr/bin/env python3
"""Constructive verification of the ceiling 4->2 correction (committed evidence in C# as the
SingleSiteField strategy; see compute/RCPsiSquared.Diagnostics/F87/PalindromeSoftCertifier.cs and
experiments/CEILING_FOUR_NONLOCAL_CASES.md).

The two I-heavy cases IXI+IIY+YII and IYI+IIX+XII, once counted in the non-local ceiling, are in fact
LOCAL. Summed over the windows of an N-chain each is a sum of weight-1 TRANSVERSE single-site fields:

    H = Σ_i (a_i X_i + b_i Y_i),   a transverse field on each site.

Single-site Paulis on different sites commute, so L = Σ_i L_i over commuting single-site Liouvillians, and
the per-site product Q = ⊗_i M_i palindromizes the whole chain. Each M_i is the per-site crossover map
Ad_{R_z(θ_i)} of the single-site X-field router, with θ_i = atan2(b_i, a_i): R_z rotates the lit field to
the X axis, where the single-site X-router R_X palindromizes it. This is CONSTRUCTIVE (Q is exhibited in
closed form), N-INDEPENDENT (additivity over commuting sites), and verified to machine precision here at
N = 4, 5, 6.

Soundness is specific to TRANSVERSE fields. A single-site X/Y field is soft: its Liouvillian spectrum
{0, −2γ, −γ ± 2i} is palindromic about the centre −γ. A single-site Z (LONGITUDINAL) field is HARD: its
spectrum {0, 0, −2γ ± 2i} has no partner for the 0 eigenvalue about −γ. The certifier therefore certifies
only transverse fields and excludes Z; the script confirms the Z exclusion as a contrast.

Self-validating: asserts the two I-heavy cases palindromize below 1e-12 at N=4,5,6 (LOCAL), the Z-middle
pair lies outside this single-site-field family (weight-3; it resists per-term routing and routes instead
via the period-4 golden router, docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md), and the single-site Z field
is hard.
Run: python simulations/ceiling_4to2_iheavy_local.py
"""
from __future__ import annotations
import sys
import numpy as np
if sys.platform == "win32":
    try: sys.stdout.reconfigure(encoding="utf-8")
    except Exception: pass

# Pauli single-letter product table P[a][b] = (letter_index, phase) for sigma_a . sigma_b.
# Letter indices: 0=I, 1=X, 2=Y, 3=Z. Matches simulations/ceiling_6to4_verification.py.
_P = [[(0, 1), (1, 1), (2, 1), (3, 1)], [(1, 1), (0, 1), (3, 1j), (2, -1j)],
      [(2, 1), (3, -1j), (0, 1), (1, 1j)], [(3, 1), (2, 1j), (1, -1j), (0, 1)]]
GAMMA = 0.05


def sprod(t, s):
    res, ph = [], 1 + 0j
    for a, b in zip(t, s):
        c, p = _P[a][b]; res.append(c); ph *= p
    return tuple(res), ph


def comm(t):
    """[T, .] as a 4^k x 4^k matrix in the Pauli basis (commutator superoperator of one term)."""
    k = len(t); dim = 4 ** k
    C = np.zeros((dim, dim), dtype=complex)
    for si in range(dim):
        s = tuple((si >> (2 * (k - 1 - j))) & 3 for j in range(k))
        ts, pts = sprod(t, s); _, pst = sprod(s, t)
        co = pts - pst
        if abs(co) > 1e-15:
            r = sum(ts[j] << (2 * (k - 1 - j)) for j in range(k)); C[r, si] = co
    return C


def embed(term, w, N):
    f = [0] * N
    for j, a in enumerate(term): f[w + j] = a
    return tuple(f)


def n_xy(idx, N):
    return sum(1 for j in range(N) if ((idx >> (2 * (N - 1 - j))) & 3) in (1, 2))


def build_L(terms, N):
    """Full Z-dephased Liouvillian L = -i[H, .] + D on the 4^N Pauli basis; centre sigma = N.gamma."""
    k = len(terms[0]); dim = 4 ** N
    C = np.zeros((dim, dim), dtype=complex)
    for w in range(N - k + 1):
        for t in terms: C += comm(embed(t, w, N))
    D = np.diag([-2 * GAMMA * n_xy(i, N) for i in range(dim)]).astype(complex)
    return -1j * C + D, N * GAMMA


def single_site_field_terms(name, N):
    """Window-sum of an I-heavy k=3 term-set over an N-chain, returned as the equivalent per-site
    transverse field {site: (a, b)} (coefficient of X_i, Y_i). IXI+IIY+YII places an X, a Y, and a Y as
    the window slides, so every interior site collects a transverse field; the chain Hamiltonian is
    H = Σ_i (a_i X_i + b_i Y_i)."""
    # (letter-per-slot) templates, letters 1=X 2=Y 3=Z, 0=I.
    templates = {
        "IXI+IIY+YII": [(0, 1, 0), (0, 0, 2), (2, 0, 0)],
        "IYI+IIX+XII": [(0, 2, 0), (0, 0, 1), (1, 0, 0)],
    }[name]
    field = {i: [0.0, 0.0] for i in range(N)}  # site -> [a (X coeff), b (Y coeff)]
    for w in range(N - 3 + 1):
        for tmpl in templates:
            for slot, letter in enumerate(tmpl):
                site = w + slot
                if letter == 1: field[site][0] += 1.0   # X
                elif letter == 2: field[site][1] += 1.0  # Y
    return {i: (ab[0], ab[1]) for i, ab in field.items() if ab[0] != 0.0 or ab[1] != 0.0}


def x_router():
    """The single-site X-field router R_X: the 4x4 per-site map (Pauli basis [I,X,Y,Z]) that palindromizes
    a single-site X-field Liouvillian L_X about its centre -gamma. It is the crossover map's representative
    for a lit X site, found in closed form by solving R L_X R^-1 = -L_X - 2.gamma over the class-swapping
    per-site maps ({I,Z} <-> {X,Y}, the -2.gamma dissipator shift). The signed permutation I<->X, Y<->Z is
    one such representative (verified by the main() assertion to residual 0)."""
    R = np.zeros((4, 4), dtype=complex)
    # Columns act on [I, X, Y, Z]. I->X, X->I, Y->Z, Z->Y. This sends L_X to -L_X-2.gamma exactly.
    R[1, 0] = 1.0    # I -> X
    R[0, 1] = 1.0    # X -> I
    R[3, 2] = 1.0    # Y -> Z
    R[2, 3] = 1.0    # Z -> Y
    return R


def rz_pauli(theta):
    """Ad_{R_z(theta)} on the single-site Pauli basis [I,X,Y,Z]: rotates X,Y by 2.theta, fixes I,Z."""
    c, s = np.cos(2 * theta), np.sin(2 * theta)
    return np.array([[1, 0, 0, 0],
                     [0, c, -s, 0],
                     [0, s, c, 0],
                     [0, 0, 0, 1]], dtype=complex)


def per_site_map(a, b):
    """The per-site crossover map M_i for a transverse field a.X + b.Y: rotate the field to the X axis by
    Ad_{R_z(-phi)}, apply the X-router R_X, rotate back. The adjoint Ad_{R_z(phi)} rotates the (X, Y) plane
    by 2.phi, so the field a.X + b.Y sits at plane-angle atan2(b, a) = 2.phi, giving phi = atan2(b, a) / 2.
    M_i = R_z(phi) R_X R_z(-phi). Because Ad_{R_z} commutes with the Z-dephasing dissipator (R_z commutes
    with Z) and fixes the identity, conjugating R_X's palindrome of L_X back through R_z(phi) yields the
    palindrome of L_{a.X + b.Y} exactly: M_i L_field M_i^-1 = -L_field - 2.gamma."""
    phi = np.arctan2(b, a) / 2.0
    Rz = rz_pauli(phi)
    return Rz @ x_router() @ np.linalg.inv(Rz)


def kron_sites(maps_by_site, N):
    """Q = ⊗_i M_i over the chain (site 0 the high-order factor, matching the basis bit order)."""
    Q = np.array([[1.0 + 0j]])
    for i in range(N):
        Q = np.kron(Q, maps_by_site[i])
    return Q


def field_resid(name, N):
    """Constructive palindrome residual ||Q L Q^-1 -(-L-2.sigma)|| / ||target|| for the single-site-field
    router Q = ⊗_i M_i of the I-heavy case `name` on an N-chain."""
    field = single_site_field_terms(name, N)
    # Build the Liouvillian directly from the window-sum of the k=3 templates (identical physics).
    templates = {"IXI+IIY+YII": [(0, 1, 0), (0, 0, 2), (2, 0, 0)],
                 "IYI+IIX+XII": [(0, 2, 0), (0, 0, 1), (1, 0, 0)]}[name]
    L, sg = build_L(templates, N)
    maps = {}
    for i in range(N):
        a, b = field.get(i, (0.0, 0.0))
        # A bare site (no field) still needs a class-swapping map for the dissipator; use the X-router
        # (theta = 0). Any class-swapping per-site map works on a site the Hamiltonian does not touch.
        maps[i] = per_site_map(a, b) if (a, b) != (0.0, 0.0) else x_router()
    Q = kron_sites(maps, N)
    Qi = np.linalg.inv(Q)
    RHS = -L - 2 * sg * np.eye(4 ** N)
    return float(np.linalg.norm(Q @ L @ Qi - RHS) / np.linalg.norm(RHS))


def single_site_spectrum(letter):
    """Spectrum of a single-site field Liouvillian (N=1) for letter in {1:X, 2:Y, 3:Z}, sorted by Re."""
    L, _ = build_L([(letter,)], 1)
    ev = np.linalg.eigvals(L)
    return ev[np.argsort(ev.real)]


def is_palindromic(ev, sigma, tol=1e-9):
    """True iff the spectrum is symmetric about the centre -sigma (every lambda has a partner -2.sigma-lambda)."""
    target = -2 * sigma - ev
    used = [False] * len(ev)
    for z in target:
        k = min((i for i in range(len(ev)) if not used[i]), key=lambda i: abs(ev[i] - z), default=None)
        if k is None or abs(ev[k] - z) > tol: return False
        used[k] = True
    return True


# The Z-middle pair: resists the UNIFORM family (ceiling_6to4_verification.py) and per-term routing, but
# is LOCAL via the period-4 golden router (docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md). Kept here for
# contrast: the single-site-field router does not apply (these are weight-3, not weight-1).
ZMIDDLE = {"XZX+XZY+YZX": [(1, 3, 1), (1, 3, 2), (2, 3, 1)],
           "YZY+XZY+YZX": [(2, 3, 2), (1, 3, 2), (2, 3, 1)]}


def main():
    # 1. Soundness ground truth: the single-site router actually palindromizes a lit X field at N=1.
    L1, s1 = build_L([(1,)], 1)
    R = x_router()
    r1 = np.linalg.norm(R @ L1 @ np.linalg.inv(R) - (-L1 - 2 * s1 * np.eye(4))) / np.linalg.norm(-L1 - 2 * s1 * np.eye(4))
    assert r1 < 1e-12, f"X-router must palindromize a single-site X field, got {r1:.1e}"

    print("SINGLE-SITE FIELD SPECTRA (soundness: transverse soft, longitudinal hard):", flush=True)
    for name, letter in (("X", 1), ("Y", 2), ("Z", 3)):
        ev = single_site_spectrum(letter)
        pal = is_palindromic(ev, GAMMA)
        evs = ", ".join(f"{z.real:+.2f}{z.imag:+.2f}i" for z in ev)
        print(f"  {name}-field: spectrum [{evs}]  palindromic about -gamma: {pal}", flush=True)
    assert is_palindromic(single_site_spectrum(1), GAMMA), "X-field must be soft"
    assert is_palindromic(single_site_spectrum(2), GAMMA), "Y-field must be soft"
    assert not is_palindromic(single_site_spectrum(3), GAMMA), "Z-field must be hard (no partner for 0)"

    print("\nTHE 2 I-HEAVY CASES, NOW LOCAL (constructive per-site product Q; verify N=4,5,6):", flush=True)
    for name in ("IXI+IIY+YII", "IYI+IIX+XII"):
        r4, r5, r6 = (field_resid(name, N) for N in (4, 5, 6))
        print(f"  {name}: single-site-field Q residual  N=4/5/6  {r4:.1e}/{r5:.1e}/{r6:.1e}", flush=True)
        assert r4 < 1e-12 and r5 < 1e-12 and r6 < 1e-12, f"{name} must palindromize (LOCAL via single-site fields)"

    print("\nTHE 2 Z-MIDDLE CASES, OUTSIDE THIS STRATEGY (weight-3, not single-site; local via the golden"
          " router):", flush=True)
    for name, terms in ZMIDDLE.items():
        # The single-site-field router does not apply (terms are weight-3). Confirm the gate rejects them:
        # there is no per-site transverse-field decomposition, so no Q = ⊗ M_i from this construction.
        applies = all(sum(1 for x in t if x != 0) == 1 and any(x in (1, 2) for x in t) for t in terms)
        print(f"  {name}: single-site-field strategy applies: {applies}  => declined here (routes via the"
              " golden router, PROOF_CEILING_GOLDEN_ROUTER)", flush=True)
        assert not applies, f"{name} is not a single-site-field set"

    print("\nOK: the 2 I-heavy are LOCAL (single-site fields); the 2 Z-middle are outside this strategy"
          " (golden-router local). 4->2 confirmed.", flush=True)


if __name__ == "__main__":
    main()
