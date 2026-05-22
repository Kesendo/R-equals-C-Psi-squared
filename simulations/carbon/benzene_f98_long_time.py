"""docs/carbon Question 5: a carbon analog of the F98 long-time bridge.

F98 (docs/ANALYTICAL_FORMULAS.md): for even N, the KIntermediate Dicke
superposition psi = (|D_{N/2-1}> + |D_{N/2}>)/sqrt(2), evolved under any
truly-class (F87) Hamiltonian + uniform Z-dephasing, has a long-time Pi^2-odd
Frobenius^2 fraction

    alpha(inf) = (N + 2) / [4 (N + 1)]   ->  1/4  as N -> inf

F98 was derived for any connected graph (the bond topology drops out: the
t -> inf limit projects onto ker L per F4) and verified bit-exact N = 4..16 on
the Heisenberg CHAIN, in simulations/water/proton_chain_dicke_anchor.py.

Question 5 asks for the BENZENE instance. Benzene's Hueckel pi-system is the
XX+YY RING (not the Heisenberg chain), and Holstein on-site dephasing is the
framework's Z-dephasing (BENZENE_LIOUVILLIAN_PALINDROME.md: D[n_l] = 1/4 D[Z_l]).
Both differ from the chain F98 was verified on. The benzene ring is still
truly-class (XX, YY bilinears are Pi^2-even) and the Holstein bath is still
Z-dephasing, so F98's preconditions hold and the long-time fraction should be
(N+2)/[4(N+1)]: 3/10 for the C4 ring, 2/7 for the C6 benzene ring.

This script evolves the KIntermediate Dicke state under the actual benzene-ring
Holstein Liouvillian (the exact t -> inf limit, via projection onto ker L) and
confirms alpha(inf). The Pi^2 machinery (the Lindbladian builder, pi2_split, the
propagator) is reused from the water F98 script; the only new piece is the XX+YY
ring Hamiltonian in place of the Heisenberg chain.
"""
import sys
import numpy as np
from math import comb
from fractions import Fraction

sys.stdout.reconfigure(encoding="utf-8")

# ---- primitives, reused from simulations/water/proton_chain_dicke_anchor.py ----
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, SX, SY, SZ]
BIT_B = [0, 0, 1, 1]                  # Pi^2 parity: Y, Z carry bit_b = 1


def kron_n(mats):
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def site_op(P, k, N):
    return kron_n([P if i == k else I2 for i in range(N)])


def z_dephasing_lindblad(H, gamma, N):
    """L for d(rho)/dt = -i[H,rho] + gamma Sum_l (Z_l rho Z_l - rho)."""
    d = 2 ** N
    eye = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, eye) - np.kron(eye, H.T))
    for k in range(N):
        Zk = site_op(SZ, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.kron(eye, eye))
    return L


def dicke_state(N, n):
    """|D_n^N>: normalised equal-amplitude superposition of all popcount-n states."""
    d = 2 ** N
    psi = np.zeros(d, dtype=complex)
    norm = 1.0 / np.sqrt(comb(N, n))
    for b in range(d):
        if bin(b).count("1") == n:
            psi[b] = norm
    return psi


def pi2_odd_total(rho, N):
    """The Pi^2-odd half of the water F98 script's pi2_split: returns
    (||rho||_F^2, ||rho_odd||_F^2) where rho_odd sums the Pauli strings whose
    bit_b letter-sum is odd."""
    d = 2 ** N
    rho_odd = np.zeros_like(rho)
    inv = 1.0 / d
    for k in range(4 ** N):
        kk, idxs = k, []
        for _ in range(N):
            idxs.append(kk & 3)
            kk >>= 2
        if sum(BIT_B[i] for i in idxs) & 1:
            sigma = kron_n([PAULIS[i] for i in idxs])
            rho_odd = rho_odd + (np.trace(sigma @ rho) * inv) * sigma
    return (np.linalg.norm(rho, "fro") ** 2,
            np.linalg.norm(rho_odd, "fro") ** 2)


# ---- the only new piece: the benzene Hueckel pi-ring Hamiltonian ---------------
def benzene_ring_hamiltonian(N):
    """XX+YY ring on N sites: H = Sum_{ring bonds} (X_a X_b + Y_a Y_b), the
    Jordan-Wigner image of Hueckel pi-hopping. Truly-class (F87): XX, YY are
    both Pi^2-even. This replaces the water script's Heisenberg chain."""
    bonds = [(l, (l + 1) % N) for l in range(N)]          # ring, incl. wrap-around
    return sum(site_op(SX, a, N) @ site_op(SX, b, N)
               + site_op(SY, a, N) @ site_op(SY, b, N)
               for a, b in bonds)


def popcount_projector(n, N):
    """P_n = Sum_{popcount(b) = n} |b><b|. F4: ker L = span(P_0, ..., P_N) for a
    truly-class Hamiltonian + Z-dephasing."""
    d = 2 ** N
    diag = [1.0 if bin(b).count("1") == n else 0.0 for b in range(d)]
    return np.diag(diag).astype(complex)


# ---- verification --------------------------------------------------------------
def verify(N, gamma=0.5):
    d = 2 ** N
    m = N // 2 - 1                                # KIntermediate: n in {m, m+1}
    alpha_f98 = Fraction(N + 2, 4 * (N + 1))

    H = benzene_ring_hamiltonian(N)
    L = z_dephasing_lindblad(H, gamma, N)

    psi = (dicke_state(N, m) + dicke_state(N, m + 1)) / np.sqrt(2.0)
    rho0 = np.outer(psi, psi.conj())
    tot0, odd0 = pi2_odd_total(rho0, N)

    print(f"=== C{N} benzene ring: F98 long-time bridge "
          f"(d={d}, L is {d * d}x{d * d}) ===")
    print(f"  KIntermediate Dicke (|D_{m}> + |D_{m + 1}>)/sqrt(2)")
    print(f"  alpha(t=0)   = {odd0 / tot0:.8f}   (F86b static anchor: 3/8 = 0.375)")

    # exact t -> inf limit: project rho0 onto ker L (the lambda = 0 eigenspace)
    lam, R = np.linalg.eig(L)
    R_inv = np.linalg.inv(R)
    coeff = R_inv @ rho0.flatten(order="F")
    ker = np.abs(lam) < 1e-9
    rho_inf = (R @ (ker * coeff)).reshape((d, d), order="F")
    rho_inf = (rho_inf + rho_inf.conj().T) / 2.0
    tot_inf, odd_inf = pi2_odd_total(rho_inf, N)
    alpha_inf = odd_inf / tot_inf

    print(f"  ker L dim    = {int(ker.sum())}   (F4 predicts N+1 = {N + 1})")
    print(f"  tr(rho_inf)  = {np.trace(rho_inf).real:.8f}   (trace is conserved: 1)")
    print(f"  alpha(inf)   = {alpha_inf:.8f}")
    print(f"  F98 (N+2)/[4(N+1)] = {alpha_f98} = {float(alpha_f98):.8f}   "
          f"|diff| = {abs(alpha_inf - float(alpha_f98)):.2e}")

    # cross-check: does the benzene dynamics land on the exact F98 rho_inf?
    rho_formula = 0.5 * (popcount_projector(m, N) / comb(N, m)
                         + popcount_projector(m + 1, N) / comb(N, m + 1))
    print(f"  ||rho_inf - F98 formula rho_inf||_F = "
          f"{np.linalg.norm(rho_inf - rho_formula, 'fro'):.2e}")

    verdict = ("F98 HOLDS on the benzene ring"
               if abs(alpha_inf - float(alpha_f98)) < 1e-9 else "MISMATCH")
    print(f"  -> {verdict}")

    # the 3/8 -> alpha(inf) decay curve: the F98 bridge, traced explicitly
    print(f"  3/8 -> {alpha_f98} decay trace:")
    for t in (0.0, 1.0, 3.0, 8.0, 20.0, 50.0):
        vec_t = R @ (np.exp(lam * t) * coeff)
        rho_t = vec_t.reshape((d, d), order="F")
        rho_t = (rho_t + rho_t.conj().T) / 2.0
        tt, oo = pi2_odd_total(rho_t, N)
        print(f"    t={t:>6.1f}   alpha(t) = {oo / tt:.8f}")
    print()


if __name__ == "__main__":
    print("Question 5: does the F98 (N+2)/[4(N+1)] long-time bridge inherit to")
    print("the benzene ring? F98 was verified on the Heisenberg chain; benzene")
    print("is the XX+YY ring under Holstein (= Z) dephasing.\n")
    verify(4)        # cyclobutadiene C4 ring: alpha(inf) = 6/20 = 3/10
    verify(6)        # benzene C6 ring:        alpha(inf) = 8/28 = 2/7
