"""
The XXZ axis as a relay: the slowest Liouvillian mode hands the baton from the bright
band-edge coherence to the Lebensader.

Walk H = J*(XX+YY) + Delta*ZZ under uniform Z-dephasing and track the SLOWEST
nonstationary mode's decay rate, oscillation frequency |Im lambda|, and n_XY composition
(weight by number of X/Y Pauli letters) as Delta goes from 0 (XY) through 1 (Heisenberg)
toward the Ising/Neel end.

What it shows:
  - XY end (Delta=0): the survivor is a pure single-magnon coherence (n_XY=1), fast-beating
    (|Im| large), decay 2*gamma. The bright band-edge coherence the clock's Rotation hand
    reads (FROST_CIRCLE_AS_THE_CLOCK_FACE).
  - Ising end (Delta large): the ZZ term slows the near-conserved I/Z population mode until,
    past a handover Delta*, it becomes the slowest mode: I/Z-dominated (n_XY=0) with a small
    magnon admixture (n_XY=2), NON-rotating (Im=0), sub-2*gamma. That is the Lebensader
    (reflections/ON_THE_ADMIXTURE_AS_LEBENSADER, experiments/CHAIN_GAP_SECTOR_DIAGNOSTIC):
    the near-conserved survivor kept alive by its magnon channel.

So the slowest mode RELAYS the baton from band-edge coherence (charge / bright) to
Lebensader (spin-Ising / still) at a handover Delta*, here located bit-by-bit.

Carbon map (docs/carbon): each carbon pi-site = a qubit; XX+YY = Huckel hopping; Z-dephasing
= Holstein phonon; ZZ = density-density (the Hubbard/PPP correlation Huckel lacks). Walking
Delta dials from the free-fermion (Huckel) charge clock to the correlated Ising Lebensader.
"""
import sys
import itertools
import warnings
sys.path.insert(0, 'simulations')
import numpy as np
import framework as fw

warnings.filterwarnings('ignore')  # ChainSystem(N=2) structural-degeneracy notice
GAMMA = 0.05
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {'I': I, 'X': X, 'Y': Y, 'Z': Z}


def site_op(M, i, N):
    mats = [I] * N
    mats[i] = M
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def H_xxz(N, Delta, J=1.0):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for i in range(N - 1):
        H += J * (site_op(X, i, N) @ site_op(X, i + 1, N) + site_op(Y, i, N) @ site_op(Y, i + 1, N))
        H += Delta * (site_op(Z, i, N) @ site_op(Z, i + 1, N))
    return H


def pauli_basis(N):
    """All 4^N Pauli strings as (4^N, d, d), with their n_XY (X/Y-letter count)."""
    d = 2 ** N
    mats = np.empty((4 ** N, d, d), dtype=complex)
    nxy = np.empty(4 ** N, dtype=int)
    for idx, letters in enumerate(itertools.product('IXYZ', repeat=N)):
        S = PAULI[letters[0]]
        for l in letters[1:]:
            S = np.kron(S, PAULI[l])
        mats[idx] = S
        nxy[idx] = sum(1 for c in letters if c in 'XY')
    return mats, nxy


def slowest_mode(N, Delta, basis_mats, basis_nxy, gamma=GAMMA):
    """(rate, |omega|, {n_XY: weight}) of the slowest nonstationary Liouvillian mode."""
    d = 2 ** N
    L = fw.lindbladian_z_dephasing(H_xxz(N, Delta), [gamma] * N)
    ev, evec = np.linalg.eig(L)
    rate = -ev.real
    nz = np.where(rate > 1e-9)[0]
    s = nz[np.argmin(rate[nz])]
    Op = evec[:, s].reshape(d, d)
    coeffs = np.einsum('sij,ij->s', basis_mats.conj(), Op) / d   # Tr(S^dag Op)/d
    w = np.abs(coeffs) ** 2
    by = {k: float(w[basis_nxy == k].sum()) for k in range(N + 1)}
    tot = sum(by.values())
    return float(rate[s]), float(abs(ev[s].imag)), {k: v / tot for k, v in by.items() if v / tot > 1e-4}


def is_lebensader(by):
    """The robust regime test: the Lebensader is I/Z-dominated (n_XY=0 outweighs the
    single-magnon n_XY=1). Im=0 alone is unreliable: at the SU(2) point the degenerate
    band-edge +/-omega pair yields a real (Im=0) combination that is still n_XY=1."""
    return by.get(0, 0.0) > by.get(1, 0.0)


def find_handover(N, basis_mats, basis_nxy, lo=0.5, hi=2.5):
    """Smallest Delta at which the slowest mode becomes the (I/Z-dominated) Lebensader."""
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        _, _, by = slowest_mode(N, mid, basis_mats, basis_nxy)
        if is_lebensader(by):
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


if __name__ == '__main__':
    for N in (4, 5):
        basis_mats, basis_nxy = pauli_basis(N)
        print(f"=== N={N}  (J=1, gamma={GAMMA}, 2*gamma={2 * GAMMA}) ===")
        print("  Delta   rate     |omega|   regime              n_XY composition")
        for Delta in [0.0, 0.5, 1.0, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0, 2.5]:
            rate, om, by = slowest_mode(N, Delta, basis_mats, basis_nxy)
            comp = '  '.join(f'{k}:{100 * v:2.0f}%' for k, v in sorted(by.items()))
            regime = 'LEBENSADER (still)' if is_lebensader(by) else 'band-edge (bright)'
            print(f"  {Delta:4.1f}   {rate:.4f}   {om:6.3f}   {regime:18s}  {comp}")
        print(f"  handover Delta* (slowest mode becomes the Lebensader) = "
              f"{find_handover(N, basis_mats, basis_nxy):.3f}")
        print()
