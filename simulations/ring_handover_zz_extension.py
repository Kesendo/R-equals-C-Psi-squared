"""Does the sqrt3 handover derivation EXTEND to the XXZ / Delta axis? Gate-first check.

My derivation Q_h = N*sqrt3/(2pi) rests on the FREE-FERMION (XY) coherence-ladder dispersion
lambda^2+8g*lambda+4J^2 q^2, whose darkness=1 point is Qq=sqrt3. The Delta* handover (ANALYTICAL_FORMULAS
"the handover Delta") is XXZ (ZZ term -> interacting fermions) and DESCENDS to Delta=1, not linear.

Structural prediction: the ZZ term breaks the free-fermion dispersion, so the (2,2) handover's Qq=sqrt3
should MOVE once Delta>0 -- i.e. sqrt3 is XY-specific, the derivation does NOT extend to the XXZ Delta axis.
GATE: if Qq at the handover stays sqrt3 for Delta>0, the extension holds (surprise); if it moves, sqrt3 is
XY-specific (the expected, honest negative). A firing gate (it moves) is the find: bounds sqrt3 to XY.
"""
import sys
import numpy as np
from itertools import combinations
from math import sqrt, pi

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
J = 1.0


def two_states(N):
    return [frozenset(c) for c in combinations(range(N), 2)]


def H_xxz_2exc(N, Delta):
    """(2,2)-sector Hamiltonian of the XXZ ring: hopping J (off-diag) + ZZ coupling Delta*J (diagonal).
    ZZ per bond = +1 if the two sites have equal occupation, -1 if they differ (Z=+/-1)."""
    states = two_states(N)
    idx = {s: a for a, s in enumerate(states)}
    D = len(states)
    H = np.zeros((D, D))
    bonds = [(i, (i + 1) % N) for i in range(N)]
    for a, s in enumerate(states):
        # diagonal ZZ
        zz = 0
        for (i, j) in bonds:
            same = (i in s) == (j in s)
            zz += 1 if same else -1
        H[a, a] += Delta * J * zz
        # off-diagonal hopping
        for (i, j) in bonds:
            if (i in s) ^ (j in s):
                moved = (s - {i}) | {j} if i in s else (s - {j}) | {i}
                H[idx[moved], a] += J
    return H, states


def hamming(s1, s2):
    return len(s1 ^ s2)


def darkness_22_xxz(N, Delta, Q):
    """darkness <n_XY> = -Re/(2g) of the slowest interior (2,2) mode of the XXZ ring at J=1, gamma=1/Q."""
    g = 1.0 / Q
    H, states = H_xxz_2exc(N, Delta)
    D = len(states)
    Id = np.eye(D)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    deph = np.array([hamming(states[a], states[b]) for a in range(D) for b in range(D)], float)
    L = L - 2 * g * np.diag(deph)
    ev = np.linalg.eigvals(L)
    interior = ev[ev.real < -1e-9]
    return -interior.real.max() / (2 * g) if interior.size else 0.0


def qh(N, Delta, lo=0.5, hi=8.0):
    flo = darkness_22_xxz(N, Delta, lo) - 1
    fhi = darkness_22_xxz(N, Delta, hi) - 1
    if flo * fhi > 0:
        return float("nan")
    for _ in range(30):
        if hi - lo < 0.01:
            break
        mid = 0.5 * (lo + hi)
        fm = darkness_22_xxz(N, Delta, mid) - 1
        if flo * fm < 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def main():
    GATE = {"fired": []}

    def gate(name, ok, detail=""):
        flag = "ok " if ok else "GATE-FIRE"
        if not ok:
            GATE["fired"].append(name)
        print(f"   [{flag}] {name}" + (f"   {detail}" if detail else ""))

    print("=" * 96)
    print("Does the (2,2) handover stay at Q*q_min = sqrt3 when the ZZ term (Delta) is turned on?")
    print("  Delta=0 = XY (free fermions, the derivation); Delta>0 = XXZ (interacting). q_min=2pi/N.")
    print("=" * 96)
    print(f"{'N':>3} {'Delta':>6} {'Q_h':>9} {'Q_h*q_min':>10} {'sqrt3=1.732':>12}")
    for N in (6, 8):
        qm = 2 * pi / N
        base = None
        for Delta in (0.0, 0.5, 1.0):
            q = qh(N, Delta)
            prod = q * qm
            if Delta == 0.0:
                base = prod
            print(f"{N:>3} {Delta:>6.2f} {q:>9.5f} {prod:>10.5f} {'':>12}")
        # gate: Delta=0 is the XY (2,2) handover (Q_h*q_min approaching sqrt3=1.732 from above with finite-N,
        # so in [sqrt3, sqrt3+0.5] at N=6,8); the ZZ term then MOVES it (a clear shift, or out of range = nan).
        q0 = qh(N, 0.0) * qm
        qh1 = qh(N, 1.0)
        q1 = qh1 * qm
        gate(f"N={N}: Delta=0 (XY) handover is the (2,2) result approaching sqrt3 from above (finite-N)",
             sqrt(3) - 0.05 < q0 < sqrt(3) + 0.5, f"Q_h*q_min(Delta=0)={q0:.4f} (-> sqrt3={sqrt(3):.4f})")
        moved = np.isnan(qh1) or abs(q1 - q0) > 0.05      # ZZ moved the handover (shift, or crossing left the range)
        gate(f"N={N}: the ZZ term MOVES the handover off the XY value (so sqrt3 is XY-specific, no Delta extension)",
             moved, f"Delta=0:{q0:.4f} vs Delta=1:{('out of range (nan)' if np.isnan(qh1) else format(q1,'.4f'))}")

    print("\n" + "=" * 96)
    print("VERDICT (honest, whatever the gates say):")
    if not GATE["fired"]:
        print("  Delta=0 sits at sqrt3 (the XY derivation) AND the ZZ term moves it -> sqrt3 is XY/free-fermion")
        print("  -specific; the derivation does NOT extend to the XXZ Delta axis. The two handovers share the")
        print("  darkness=1 floor but are distinct mechanisms (XY free-fermion dispersion vs XXZ interacting;")
        print("  Q_h grows ~N, Delta* descends to 1). The Delta* closed form stays the open Bethe-ansatz problem.")
    else:
        print(f"  GATES FIRED: {GATE['fired']} -- diagnose (did Delta=0 not give sqrt3? did ZZ NOT move it?).")
    print("=" * 96)


if __name__ == "__main__":
    main()
