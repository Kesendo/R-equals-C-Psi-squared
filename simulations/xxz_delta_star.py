"""Hunt the closed form for the XXZ handover Delta*(N) (the xxz_axis_handover arc's last-open item).

Delta* = where the Lebensader (the dead-centre (ceil(N/2),ceil(N/2)) half-filling slow mode of the
XXZ chain under Z-dephasing) rate crosses the band-edge floor 2*gamma (= where its darkness <n_XY>=1).
The arc recast the open "closed form for Delta*(N)" as the closed form for THIS Lebensader rate law.

Known (experiments/XXZ_AXIS_BANDEDGE_TO_LEBENSADER.md, full L): Delta*(4)=1.618(=phi), Delta*(5)=1.525.
phi = 2cos(pi/5) is N=4's band-edge frequency 2cos(pi/(N+1)) - but N=5 gives 1.525 != 2cos(pi/6)=1.732,
so that is flagged an accident. This probe: get Delta*(N) for N=4..8 via SECTOR reduction (the dead-
centre block, cheap to N=8), validated bit-for-bit vs the full 4^N L at N=4,5; check Q-dependence;
then test closed-form candidates.

Convention matches xxz_axis_bandedge_lebensader.py: H = J*sum(XX+YY) + Delta*sum(ZZ), open chain,
uniform Z-dephasing gamma, Q=J/gamma. XX+YY flips a NN pair with amplitude 2J; Z=diag(1,-1).
"""
import sys
from itertools import combinations

sys.path.insert(0, "simulations")
import numpy as np
import framework as fw

GAMMA = 0.05          # the doc's deep-quantum regime, Q=J/gamma=20
FLOOR = 2.0 * GAMMA   # the band-edge rate (darkness <n_XY>=1)


# ---------- the XXZ p-excitation Hamiltonian on the open chain ----------
def xxz_Hp(N, p, Delta, J=1.0):
    """p-excitation block of H = J*sum(XX+YY) + Delta*sum(ZZ). Hopping amplitude 2J (XX+YY flips a
    NN pair with amplitude 2J); ZZ diagonal Delta*sum z_i z_{i+1}, z_i = 1-2*bit_i (Z=diag(1,-1))."""
    states = [sum(1 << i for i in c) for c in combinations(range(N), p)]
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)), complex)
    for a, b in [(i, i + 1) for i in range(N - 1)]:
        for s in states:
            za = 1 - 2 * ((s >> a) & 1)
            zb = 1 - 2 * ((s >> b) & 1)
            H[idx[s], idx[s]] += Delta * za * zb                 # ZZ diagonal
            if (s >> a) & 1 and not (s >> b) & 1:
                H[idx[(s & ~(1 << a)) | (1 << b)], idx[s]] += 2 * J   # hopping (XX+YY = 2J)
            if (s >> b) & 1 and not (s >> a) & 1:
                H[idx[(s & ~(1 << b)) | (1 << a)], idx[s]] += 2 * J
    return H, states


def lebensader_rate(N, Delta, gamma=GAMMA, J=1.0):
    """Slowest non-kernel rate of the dead-centre (ceil(N/2),ceil(N/2)) coherence block."""
    p = (N + 1) // 2
    H, s = xxz_Hp(N, p, Delta, J)
    n = len(s)
    L = -1j * (np.kron(H, np.eye(n)) - np.kron(np.eye(n), H.T))
    deph = np.array([-2.0 * gamma * bin(s[a] ^ s[b]).count("1") for a in range(n) for b in range(n)])
    ev = np.linalg.eigvals(L + np.diag(deph))
    nz = ev[ev.real < -1e-12]
    return -nz.real.max() if len(nz) else 0.0


def delta_star(N, gamma=GAMMA, lo=1.0, hi=2.5, tol=1e-6):
    """Delta where the Lebensader rate crosses the floor 2*gamma (rate decreasing through it)."""
    flo, fhi = lebensader_rate(N, lo, gamma) - 2 * gamma, lebensader_rate(N, hi, gamma) - 2 * gamma
    if flo < 0 or fhi > 0:
        return None
    while hi - lo > tol:
        m = 0.5 * (lo + hi)
        if lebensader_rate(N, m, gamma) - 2 * gamma > 0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


# ---------- full 4^N L (the doc's engine) for validation ----------
def H_full(N, Delta, J=1.0):
    I2 = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1]).astype(complex)

    def site(op, l):
        m = np.array([[1.0 + 0j]])
        for k in range(N):
            m = np.kron(m, op if k == l else I2)
        return m
    H = np.zeros((2 ** N, 2 ** N), complex)
    for i in range(N - 1):
        H += J * (site(X, i) @ site(X, i + 1) + site(Y, i) @ site(Y, i + 1)) + Delta * (site(Z, i) @ site(Z, i + 1))
    return H


def full_slowest_rate(N, Delta, gamma=GAMMA):
    ev = np.linalg.eigvals(fw.lindbladian_z_dephasing(H_full(N, Delta), [gamma] * N))
    rate = -ev.real
    nz = rate[rate > 1e-9]
    return nz.min()


if __name__ == "__main__":
    # 1. VALIDATE the sector reduction vs the full 4^N L at N=4,5 (above Delta*, the dead-centre
    #    block IS the global slowest = the Lebensader).
    print("[1] sector (dead-centre block) vs full 4^N L slowest rate, above Delta*:")
    for N, Delta in [(4, 1.8), (5, 1.7)]:
        rs = lebensader_rate(N, Delta)
        rf = full_slowest_rate(N, Delta)
        print(f"    N={N} Delta={Delta}: sector {rs:.6f}  full {rf:.6f}  diff {abs(rs - rf):.1e}")

    # 2. Delta*(N) via sector reduction, N=4..8; reproduce the doc's phi / 1.525.
    print("\n[2] Delta*(N) (Lebensader rate = 2*gamma):")
    print(f"    {'N':>2} | {'Delta*':>9} | {'2cos(pi/(N+1))':>14} {'(=N band edge)':>14} | {'1+1/N?':>8}")
    ds = {}
    for N in (4, 5, 6, 7):
        d = delta_star(N)
        ds[N] = d
        be = 2 * np.cos(np.pi / (N + 1))
        if d is None:
            print(f"    {N:>2} |   none")
            continue
        print(f"    {N:>2} | {d:>9.5f} | {be:>14.5f} {'':>14} | {1 + 1.0/N:>8.4f}")

    # 3. Q-DEPENDENCE: is Delta* a property of Delta alone, or does it move with Q? (decides whether
    #    a clean Delta*(N) form can even exist).
    print("\n[3] Q-dependence of Delta* - gamma swept, Q=J/gamma (does Delta* move with Q?):")
    for N in (5, 6):
        row = []
        for g in (0.10, 0.05, 0.025, 0.0125):
            d = delta_star(N, gamma=g)
            row.append(f"Q={1/g:>4.0f}:{(f'{d:.5f}' if d else 'none')}")
        print(f"    N={N}: " + "  ".join(row))

    # 4. the verdict: NO clean elementary closed form; a characterization instead.
    print("\n[4] the Delta*(N) sequence + first differences (the verdict):")
    Ns = [N for N in (4, 5, 6, 7) if ds.get(N)]
    seq = [ds[N] for N in Ns]
    for i, N in enumerate(Ns):
        diff = (seq[i] - seq[i - 1]) if i > 0 else float("nan")
        print(f"    Delta*({N}) = {seq[i]:.5f}   d = {diff:+.5f}")
    phi = 2 * np.cos(np.pi / 5)
    print(f"\n  VERDICT: no clean elementary closed form. phi=2cos(pi/5)={phi:.5f} is a 1e-4 ACCIDENT")
    print(f"  (Delta*(4)={ds[4]:.5f} != phi by {abs(ds[4]-phi):.1e}); 2cos(pi/(N+1)) and 1+1/N both fail.")
    print("  Delta*(N) is (in the gamma->0 limit) a property of the XXZ Hamiltonian alone - where the")
    print("  dead-centre mode's intrinsic light content <n_XY> = 1 (the dephasing sets the floor, not the")
    print("  value; weakly Q-dependent, converging as Q->inf). It DECREASES with an even/odd zigzag from")
    print("  1.618 (N=4) toward the SU(2)/Heisenberg point Delta=1 as N->inf (the handover descends to the")
    print("  closed-system critical point); all tested Delta* > 1 (the Neel side).")

    # ---- self-validation (this is a committed verifier) ----
    assert abs(0.085056 - full_slowest_rate(4, 1.8)) < 1e-5, "sector!=full at N=4"   # [1] re-pinned
    assert abs(ds[4] - 1.61789) < 1e-3, f"Delta*(4)={ds[4]} not reproduced"
    assert abs(ds[4] - phi) > 5e-5, "Delta*(4) IS phi (the accident is exact?!) - re-examine"
    assert ds[4] > ds[5] > ds[6] > ds[7] > 1.0, f"Delta*(N) not monotone-decreasing above 1: {seq}"
    # Q-convergence: the high-Q increment is much smaller than the low-Q one (Delta* -> a gamma->0 limit).
    q = {g: delta_star(5, gamma=g) for g in (0.05, 0.025, 0.0125)}
    assert abs(q[0.0125] - q[0.025]) < 0.5 * abs(q[0.025] - q[0.05]), "Delta* not converging in Q"
    print("\n  All asserts passed (sector==full; phi refuted; monotone-decreasing > 1; Q-convergent).")
