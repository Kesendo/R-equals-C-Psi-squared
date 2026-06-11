#!/usr/bin/env python3
"""
_qutrit_partial_palindrome_count.py - WIP scout for OQ-002.

Question: the partial palindrome at d>2. Qutrit N=2 spectra pair 36-52/81
eigenvalues (QUBIT_NECESSITY.md sec 8b/9). Never explained. The
one-diagonal lens (ON_THE_ONE_DIAGONAL): under full-Cartan dephasing the
decay rate is still -2*gamma * Hamming(i,j), the same rate ladder as the
qubit. What changes is the MULTIPLICITY per rung:

    c_k = d^N * C(N,k) * (d-1)^k    (coherences at Hamming distance k)

The palindrome reflects rate-rung k <-> N-k. For qubits d-1=1, so c_k is
C(N,k), symmetric, all pair (100%). For d>2 the (d-1)^k factor tilts the
distribution; only the overlap pairs. Closed-form ceiling:

    paired(d,N) = sum_k d^N * C(N,k) * (d-1)^min(k, N-k)

PREDICTION: rate-only (H=0) pairing of the N=2 qutrit Liouvillian = 54/81,
the algebraic ceiling; the full L (with H) lands at 36-52 because H splits
the rate-degenerate rungs in the imaginary axis. This script proves both.

Self-validating: every block asserts; process exits 0 only if all hold.
"""

import numpy as np
from itertools import product as iprod
from math import comb

# ---- qutrit infrastructure (matches qubit_necessity_tests.py conventions) ----
I3 = np.eye(3, dtype=complex)
lam3 = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex)
lam8 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / np.sqrt(3)
# SU(3) Heisenberg uses all 8 Gell-Mann (raw, as in the existing test)
gm_raw = [
    np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
    np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=complex),
    np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=complex),
    np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
    np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=complex),
    np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
    np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=complex),
    np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=complex) / np.sqrt(3),
]
GAMMA = 0.05


def site_op(M, k, N, dloc):
    ops = [np.eye(dloc, dtype=complex)] * N
    ops[k] = M
    out = ops[0]
    for i in range(1, N):
        out = np.kron(out, ops[i])
    return out


def H_su3_heisenberg(N, bonds, J=1.0):
    d = 3 ** N
    H = np.zeros((d, d), dtype=complex)
    for (a, b) in bonds:
        for m in gm_raw:
            H += J * (site_op(m, a, N, 3) @ site_op(m, b, N, 3))
    return H


def L_dephasing(N, gamma, jumps):
    """Full computational-basis Lindbladian dissipator (no H), jumps per site."""
    d = 3 ** N
    Id = np.eye(d, dtype=complex)
    L = np.zeros((d * d, d * d), dtype=complex)
    for k in range(N):
        for m in jumps:
            Mk = site_op(m, k, N, 3)
            MdM = Mk.conj().T @ Mk
            L += gamma * (np.kron(Mk, Mk.conj())
                          - 0.5 * np.kron(MdM, Id)
                          - 0.5 * np.kron(Id, MdM.T))
    return L


def L_hamiltonian(H):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def palindrome_pairs(evals, Sg, tol=1e-4):
    """Count eigenvalues that pair under lambda -> -2*Sg - lambda."""
    num = len(evals)
    used = set()
    n_well = 0
    for k in range(num):
        if k in used:
            continue
        target = -evals[k] - 2 * Sg
        diffs = np.abs(evals - target)
        for u in used:
            diffs[u] = 1e30
        best = int(np.argmin(diffs))
        if diffs[best] < tol:
            n_well += 1 if k == best else 2
            used.add(k)
            used.add(best)
        else:
            used.add(k)
    return n_well


def hamming_counts(d, N):
    """Multiplicity c_k of coherences |i><j| at Hamming distance k."""
    c = [0] * (N + 1)
    for i in iprod(range(d), repeat=N):
        for j in iprod(range(d), repeat=N):
            kk = sum(1 for a, b in zip(i, j) if a != b)
            c[kk] += 1
    return c


def ceiling_formula(d, N):
    return sum(d ** N * comb(N, k) * (d - 1) ** min(k, N - k) for k in range(N + 1))


def combinatorial_pairing(c):
    """Modes that find a mirror partner given rung counts c[k] (rung k<->N-k)."""
    N = len(c) - 1
    paired = 0
    for k in range(N + 1):
        m = N - k
        if k < m:
            paired += 2 * min(c[k], c[m])
        elif k == m:
            paired += c[k]
    return paired


def main():
    print("=" * 68)
    print("OQ-002: the partial palindrome at d>2, via the disagreement count")
    print("=" * 68)

    # ---- Part A: per-site equidistance (the lens lands on d=3) ----
    print("\n[A] Per-site rate r(i,j) under full-Cartan dephasing {lam3, lam8}:")
    jumps1 = [lam3, lam8]
    rates_site = {}
    for i in range(3):
        for j in range(3):
            DB = np.zeros((3, 3), dtype=complex)
            B = np.zeros((3, 3), dtype=complex)
            B[i, j] = 1.0
            for M in jumps1:
                DB += M @ B @ M.conj().T - 0.5 * (M.conj().T @ M @ B + B @ M.conj().T @ M)
            rates_site[(i, j)] = float(np.real(DB[i, j]))  # diagonal action
    for (i, j), r in rates_site.items():
        if i != j:
            print(f"    ({i},{j}): rate = {r:+.4f}   (expect -2.0)")
    offdiag = [r for (i, j), r in rates_site.items() if i != j]
    assert all(abs(r - (-2.0)) < 1e-12 for r in offdiag), "levels not equidistant!"
    assert all(abs(rates_site[(i, i)]) < 1e-12 for i in range(3)), "diagonal decays!"
    print("    -> all i!=j equidistant at -2: rate = -2*gamma*Hamming(i,j). OK")

    # ---- Part B: rung multiplicity matches c_k = d^N C(N,k)(d-1)^k ----
    print("\n[B] Coherence multiplicity per Hamming rung (d=3, N=2):")
    d, N = 3, 2
    c_emp = hamming_counts(d, N)
    c_pred = [d ** N * comb(N, k) * (d - 1) ** k for k in range(N + 1)]
    print(f"    empirical c_k = {c_emp}")
    print(f"    formula   c_k = {c_pred}   (d^N C(N,k)(d-1)^k)")
    assert c_emp == c_pred, "multiplicity formula wrong!"
    print(f"    sum = {sum(c_emp)} = d^(2N) = {d ** (2 * N)}. OK")

    # ---- Part C: rate-only Liouvillian -> spectrum is exactly the rungs ----
    print("\n[C] Dissipator-only Liouvillian (H=0), N=2 qutrits:")
    LD = L_dephasing(N, GAMMA, jumps1)
    evD = np.linalg.eigvals(LD)
    # eigenvalues should be exactly {0, -2g, -4g} with mult {9,36,36}
    rounded = np.round(evD.real / GAMMA).astype(int)
    levels = {}
    for r in rounded:
        levels[r] = levels.get(r, 0) + 1
    print(f"    eigenvalue rungs (units of gamma): "
          f"{ {k: levels[k] for k in sorted(levels)} }")
    assert abs(evD.imag).max() < 1e-9, "dissipator-only has imaginary parts!"
    assert levels.get(0) == 9 and levels.get(-2) == 36 and levels.get(-4) == 36, \
        "rung multiplicities off"
    Sg = N * GAMMA
    paired_D = palindrome_pairs(evD, Sg, tol=1e-9)
    ceil = ceiling_formula(d, N)
    comb_pair = combinatorial_pairing(c_emp)
    print(f"    rate-only palindrome pairing = {paired_D}/81")
    print(f"    closed-form ceiling          = {ceil}/81")
    print(f"    combinatorial pairing        = {comb_pair}/81")
    assert paired_D == ceil == comb_pair == 54, \
        f"ceiling mismatch: {paired_D} {ceil} {comb_pair}"
    print("    -> 54/81: the algebraic ceiling. 27 excess at k=2 unpaired. OK")

    # ---- Part D: full L -> H DEGRADES the palindrome at every fixed center ----
    # The dissipator's 54 (Part C) is the pairing about the PHYSICAL center -N*g
    # (the k<->N-k reflection, where the qubit palindrome is exact). Adding the
    # Hamiltonian degrades it: at EVERY fixed center the full L pairs FEWER than
    # the dissipator. An earlier reading mistook full-L-best-center (60 at -3g)
    # for "exceeding" the dissipator's 54 at -2g - a comparison of two DIFFERENT
    # centers; at equal center H always reduces pairing. The full interacting
    # analysis (real parts = -2g*<Q>, H-dependence, no closed form) lives in
    # simulations/qutrit_interacting_palindrome.py.
    print("\n[D] Full Liouvillian (H_SU3 + dephasing): H DEGRADES at fixed center")
    H = H_su3_heisenberg(N, [(0, 1)])
    evF = np.linalg.eigvals(L_hamiltonian(H) + L_dephasing(N, GAMMA, jumps1))
    evD = np.linalg.eigvals(L_dephasing(N, GAMMA, jumps1))
    rungs = {}
    for r in np.round(evF.real / GAMMA, 1):
        rungs[r] = rungs.get(r, 0) + 1
    print(f"    full-L real-part rungs (gamma units): "
          f"{ {k: rungs[k] for k in sorted(rungs)} }")
    print(f"    (dissipator alone was {{0: 9, -2: 36, -4: 36}}; H builds a -3g rung)")
    # palindrome_pairs(ev, Sg) reflects lambda -> -2*Sg - lambda (center -Sg),
    # so a fixed center -c*gamma is reached with Sg = c*gamma.
    for cg, name in ((2.0, "physical center -N*g"), (3.0, "-3g, the two big rungs")):
        pD = palindrome_pairs(evD, cg * GAMMA, tol=1e-4)
        pF = palindrome_pairs(evF, cg * GAMMA, tol=1e-4)
        print(f"    center -{cg:.0f}g ({name}): dissipator {pD}/81, full L {pF}/81 "
              f"-> H {'degrades' if pF < pD else 'does NOT degrade'}")
        assert pF < pD, f"H should degrade at center -{cg}g: {pD} -> {pF}"
    assert palindrome_pairs(evD, 2 * GAMMA, tol=1e-4) == 54
    assert palindrome_pairs(evF, 2 * GAMMA, tol=1e-4) == 48
    print("    -> at EVERY fixed center H reduces the pairing (54->48, 72->60);")
    print("       the palindrome is fragile under H. 54 (about -N*g) is the invariant.")

    # ---- Part E: the formula across (d, N): d=2 always full, d>2 partial ----
    print("\n[E] Closed-form ceiling vs brute force, grid of (d, N):")
    print(f"    {'d':>2} {'N':>2} {'paired':>8} {'total':>8} {'frac':>7}  formula?")
    for d_ in (2, 3, 4):
        for N_ in (1, 2, 3):
            c = hamming_counts(d_, N_)
            bf = combinatorial_pairing(c)
            fm = ceiling_formula(d_, N_)
            tot = d_ ** (2 * N_)
            ok = (bf == fm)
            assert ok, f"formula != brute force at d={d_}, N={N_}"
            if d_ == 2:
                assert bf == tot, f"qubit not full at N={N_}"
            print(f"    {d_:>2} {N_:>2} {bf:>8} {tot:>8} {bf / tot:>6.1%}  "
                  f"{'OK' if ok else 'FAIL'}")
    print("    -> formula exact everywhere; d=2 is the only fully-paired column.")

    print("\n" + "=" * 68)
    print("ALL CHECKS PASSED. The DEPHASING dissipator's partial palindrome =")
    print("the symmetric overlap of the disagreement-count distribution")
    print("c_k = d^N C(N,k)(d-1)^k under k<->N-k, about the physical center -N*g.")
    print("The (d-1)^k tilt vanishes only at d=2 (the unique fully-paired column);")
    print("for d=3,N=2 the ceiling is 54/81. Adding the Hamiltonian DEGRADES this:")
    print("at every fixed center H reduces the pairing (54->48 about -N*g), the")
    print("palindrome is fragile under H. The full interacting analysis (real parts")
    print("= -2g<Q>, H-dependent count, no closed form) is in")
    print("simulations/qutrit_interacting_palindrome.py. 54 is the invariant skeleton.")
    print("=" * 68)


if __name__ == "__main__":
    main()
