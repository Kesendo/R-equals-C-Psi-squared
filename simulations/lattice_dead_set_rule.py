"""Gate for the dead-set rule (experiments/LATTICE_DEAD_SET_RULE.md): which Pauli
readouts are identically zero at every h on the h-thread setup, as a closed rule.

Setup (experiments/LATTICE_H_THREAD.md): open XY chain + longitudinal field h_l Z_l,
local Z-dephasing (any watching profile), prep = the operator-X-even population pair
rho(0) = (P_s + P_~s)/2, optionally + a cross coherence c|s><~s| + h.c. (real or
imaginary c). Scan h along a fixed generic direction; a readout O is ALIVE if
<O>(t) is not identically zero over the scan, DEAD otherwise.

THE RULE (three conserved structures; each necessity is exact, joint sufficiency is
the gated observation):

  K  (popcount blocks)  H conserves total popcount at every h and the dephasing mask
     moves nothing, so rho(t) is supported exactly on the (p,q) blocks populated at
     t = 0. O with XY-mask weight w must connect a supported block:
     exists a = (p + w - q)/2 with 0 <= a <= w, a <= p, p - a <= N - w.

  V  (the doubly-mirrored kill)  For each sublattice gauge U_g (Z on the even or the
     odd sites), W_g = U_g X^N sends H(h) -> -H(h) as an operator; entrywise
     conjugation leaves the real H untouched but is ANTIUNITARY (it reverses the i,
     so conj(rho) solves the flow of -H); the dephasing mask, real and a function of
     i XOR j only, is invariant under all three ingredients. So the antiunitary
     V_g(rho) = W_g conj(rho) W_g^dag is an exact symmetry of the flow AT EVERY
     FIXED h (two h-flipping mirrors compose to an h-preserving one). If V_g rho(0) = rho(0), then for every Pauli readout
     <O> = eps_g <O> with eps_g = (-1)^(n_Z + xy_g), xy_g = #{X/Y letters on
     sublattice g}; eps_g = -1 kills O identically. Stabilizer bookkeeping (checked
     directly here): the population part always; the cross-coherence part iff the
     sublattice size |g| is even (real AND imaginary c alike).

  F  (fermion degree)  Under the left Jordan-Wigner map every Pauli string is one
     Majorana monomial with a definite degree d(O); the XY + field Hamiltonian is
     quadratic and the adjoint dissipator is DIAGONAL on Pauli strings, so the flow
     preserves d exactly. rho(0) touches only even degrees (populations expand in
     pure-Z strings) plus degree N (the cross coherence expands in all-XY strings).
     So O with odd d is dead unless the coherence is present and d(O) = N.
     Visible only where K and V do not already kill: the first sighting is N = 5
     with a popcount-2 seed, where the 2N Jordan-Wigner-linear strings (d = 1) and
     their Gamma-duals (d = 2N-1) slip through K and V and are dead by F alone.

  THE MOD-4 SHADOW (T3/T4, the sharpened form): V's kill sign is a pure function of
  the Majorana degree, N-free: eps_odd-gauge = (-1)^(d(d-1)/2) and eps_even-gauge =
  (-1)^(d(d+1)/2) (see eps_degree for the one-paragraph derivation; checked against
  the letter formula AND direct matrix conjugation for all strings, N = 2..6). In
  this language V_g is ALWAYS a symmetry; the prep splits into V-eigen-sectors
  (population +1 under both gauges, coherence (-1)^|g|), each sector's contribution
  to <O> dies unless eps_g matches the sector sign for both g, and at d = N the
  coherence signs match AUTOMATICALLY (both N parities). The whole rule collapses to

      alive  iff  (K_pop and d = 0 mod 4)  or  (coherence and K_coh and d = N)

  (necessity derived; sufficiency, as everywhere in this gate, is the gated
  observation). The global-stabilizer form above is the special case where the readable sectors
  share their signs; the two forms diverge at N = 2 mod 4 with coherence (no gauge
  stabilizes globally: at N = 6 the global form over-predicts 2047 alive vs the
  actual 1055, the collapsed form stays exact). Alive-by-degree counts come out as
  binomial coefficients cut by kinematics (asserted via math.comb).

  Boundary (zz control): with a ZZ coupling the chiral gauge no longer flips H
  (kills V) and H is quartic (kills F); layer K persists, and the one fermionic
  survivor is the conserved parity Z^N, whose expectation starts (and stays) 0 for
  the odd-N population pair. Gated: at N = 4 the dead set collapses to exactly
  layer K; at N = 3 (pop) it is layer K plus {ZZZ}.

Exactness = predicted set == alive set over the full 4^N - 1 census, per config.
Runtime ~5 min (the N = 6 censuses dominate). Scout history: seven local scouts (the
candidate rule of the 2026-07-16 handover falsified as stated, then rebuilt as layer
V; the mod-4 shadow conjectured from the Majorana bookkeeping, then pinned exactly).
"""
import sys
from itertools import product

import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
s0 = np.eye(2, dtype=complex)
PAULI = {"I": s0, "X": sx, "Y": sy, "Z": sz}

FAILURES = []


def check(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILURES.append(name)


def eps_degree(d, gauge):
    """The mod-4 shadow, exact: eps_odd = (-1)^(d(d-1)/2), eps_even = (-1)^(d(d+1)/2).
    N-free. Derivation: X^N sends a_2l -> (-1)^l a_2l, a_2l+1 -> (-1)^(l+1) a_2l+1;
    conj fixes a_2l and flips a_2l+1; U_g flips both Majoranas of the sites in g. The
    combined per-site sign is (-1)^(l + [l in g]): uniformly -1 for the even-sites
    gauge, uniformly +1 for the odd-sites gauge. A Hermitian degree-d string is
    (+-) i^(d(d-1)/2) times an ordered Majorana monomial, and the antilinear map
    conjugates that phase to a sign. Hence eps_g = (-1)^(d(d-1)/2) * sigma_g^d."""
    return (-1) ** ((d * (d - 1)) // 2) if gauge == "odd" else (-1) ** ((d * (d + 1)) // 2)


def majorana_degree(name):
    """Left-JW Majorana degree of a Pauli string (a Majorana at site m puts Z on
    every site left of m), computed right-to-left via the tail parity tau."""
    tau = 0
    d = 0
    for l in range(len(name) - 1, -1, -1):
        ch = name[l]
        if ch in "XY":
            n = 1
        elif ch == "Z":
            n = 2 if tau % 2 == 0 else 0
        else:
            n = 0 if tau % 2 == 0 else 2
        d += n
        tau += n
    return d


def block_readable(N, w, support):
    for p, q in support:
        if (p + w - q) % 2:
            continue
        a = (p + w - q) // 2
        if 0 <= a <= w and a <= p and p - a <= N - w:
            return True
    return False


def census(N, seed, coh_kind, gammas, zz):
    """Run the h-scan and classify every readout; return per-string records."""
    D = 2 ** N

    def op(single, site):
        ops = [s0] * N
        ops[site] = single
        out = ops[0]
        for o in ops[1:]:
            out = np.kron(out, o)
        return out

    def bit(i, l):
        return (i >> (N - 1 - l)) & 1

    def hamiltonian(hs):
        H = np.zeros((D, D), dtype=complex)
        for b in range(N - 1):
            H += op(sx, b) @ op(sx, b + 1) + op(sy, b) @ op(sy, b + 1)
            if zz:
                H += zz * op(sz, b) @ op(sz, b + 1)
        for l in range(N):
            H += hs[l] * np.diag([1 - 2 * bit(i, l) for i in range(D)]).astype(complex)
        return H

    XN = op(sx, 0)
    for site in range(1, N):
        XN = XN @ op(sx, site)
    Ug = {}
    for gname, par in (("even", 0), ("odd", 1)):
        U = np.eye(D, dtype=complex)
        for l in range(N):
            if l % 2 == par:
                U = U @ op(sz, l)
        Ug[gname] = U

    mn = np.zeros((D, D))
    for i in range(D):
        for j in range(D):
            mn[i, j] = -2 * sum(g * (bit(i, l) ^ bit(j, l)) for l, g in enumerate(gammas))

    def rk4(rho, H, dt, ticks):
        def rhs(r):
            return -1j * (H @ r - r @ H) + mn * r
        out = [rho.copy()]
        for _ in range(ticks):
            k1 = rhs(rho)
            k2 = rhs(rho + dt / 2 * k1)
            k3 = rhs(rho + dt / 2 * k2)
            k4 = rhs(rho + dt * k3)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            out.append(rho.copy())
        return out

    rb = np.zeros((D, D), dtype=complex)
    rb[seed, seed] = 0.5
    rb[D - 1 - seed, D - 1 - seed] = 0.5
    ps, pt = bin(seed).count("1"), N - bin(seed).count("1")
    support = {(ps, ps), (pt, pt)}
    coh = coh_kind is not None
    if coh_kind == "real":
        rb[seed, D - 1 - seed] = 0.3
        rb[D - 1 - seed, seed] = 0.3
    elif coh_kind == "imag":
        rb[seed, D - 1 - seed] = 0.3j
        rb[D - 1 - seed, seed] = -0.3j
    if coh:
        support |= {(ps, pt), (pt, ps)}

    stab = []
    for gname, U in Ug.items():
        W = U @ XN
        if np.max(np.abs(W @ rb.conj() @ W.conj().T - rb)) < 1e-14:
            stab.append(gname)

    t_scan = 0.4
    hdir = ([1.0, -0.6, 0.3] + [0.5] * N)[:N]
    dt, ticks = 0.05, 40
    rp = rk4(rb, hamiltonian([t_scan * v for v in hdir]), dt, ticks)
    rm = rk4(rb, hamiltonian([-t_scan * v for v in hdir]), dt, ticks)

    # the V symmetry itself, pinned dynamically for every stabilizing gauge
    v_worst = 0.0
    for gname in stab:
        W = Ug[gname] @ XN
        v_worst = max(v_worst, max(float(np.max(np.abs(W @ r.conj() @ W.conj().T - r)))
                                   for r in rp))

    rows = []
    for letters in product("IXYZ", repeat=N):
        name = "".join(letters)
        if name == "I" * N:
            continue
        O = op(PAULI[letters[0]], 0)
        for site in range(1, N):
            O = O @ op(PAULI[letters[site]], site)
        fp = [float(np.real(np.trace(r @ O))) for r in rp]
        fm = [float(np.real(np.trace(r @ O))) for r in rm]
        mag = max(max(abs(a) for a in fp), max(abs(b) for b in fm))
        alive = mag > 1e-10

        nz = sum(1 for ch in name if ch == "Z")
        w = sum(1 for ch in name if ch in "XY")
        k_ok = block_readable(N, w, support)
        v_ok = True
        for gname, par in (("even", 0), ("odd", 1)):
            if gname in stab:
                xy_g = sum(1 for l, ch in enumerate(name) if ch in "XY" and l % 2 == par)
                if (nz + xy_g) % 2 == 1:
                    v_ok = False
        d = majorana_degree(name)
        f_ok = (d % 2 == 0) or (coh and d == N)
        rows.append((name, alive, mag, k_ok, v_ok, f_ok, d, w))
    return rows, stab, v_worst


def run_config(N, seed, coh_kind, gammas, zz, label, use_f=True):
    rows, stab, v_worst = census(N, seed, coh_kind, gammas, zz)
    pred_alive = {r[0] for r in rows if r[3] and r[4] and (r[5] or not use_f)}
    alive = {r[0] for r in rows if r[1]}
    total = len(rows)
    check(label, pred_alive == alive and v_worst < 1e-12,
          f"alive {len(alive)}/{total}, predicted {len(pred_alive)}, stab={stab}, "
          f"V-invariance worst {v_worst:.1e}")
    return rows


def t4_mod4_identity(maxN=6):
    """Check eps_letters == eps_degree == the direct matrix conjugation, ALL strings."""
    worst_bad = 0
    total = 0
    for N in range(2, maxN + 1):
        D = 2 ** N

        def op(single, site):
            ops = [s0] * N
            ops[site] = single
            out = ops[0]
            for o in ops[1:]:
                out = np.kron(out, o)
            return out

        XN = op(sx, 0)
        for site in range(1, N):
            XN = XN @ op(sx, site)
        W = {}
        for gname, par in (("even", 0), ("odd", 1)):
            U = np.eye(D, dtype=complex)
            for l in range(N):
                if l % 2 == par:
                    U = U @ op(sz, l)
            W[gname] = U @ XN
        for letters in product("IXYZ", repeat=N):
            name = "".join(letters)
            O = op(PAULI[letters[0]], 0)
            for site in range(1, N):
                O = O @ op(PAULI[letters[site]], site)
            d = majorana_degree(name)
            nz = sum(1 for ch in name if ch == "Z")
            for gname, par in (("even", 0), ("odd", 1)):
                xy_g = sum(1 for l, ch in enumerate(name) if ch in "XY" and l % 2 == par)
                e_let = (-1) ** (nz + xy_g)
                e_deg = eps_degree(d, gname)
                Wm = W[gname]
                res = float(np.max(np.abs((Wm @ O @ Wm.conj().T).conj() - e_deg * O)))
                total += 1
                if e_let != e_deg or res > 1e-12:
                    worst_bad += 1
    return worst_bad, total


def run_collapsed(N, seed, coh_kind, gammas, label, expect_hist):
    """The per-sector (collapsed fermionic) rule:
    alive iff (K_pop and d = 0 mod 4) or (coherence and K_coh and d = N)."""
    rows, stab, v_worst = census(N, seed, coh_kind, gammas, 0.0)
    ps, pt = bin(seed).count("1"), N - bin(seed).count("1")
    sup_pop = {(ps, ps), (pt, pt)}
    sup_coh = {(ps, pt), (pt, ps)}
    coh = coh_kind is not None
    pred_alive = set()
    hist = {}
    for name, alive, mag, k_ok, v_ok, f_ok, d, w in rows:
        pred = (block_readable(N, w, sup_pop) and d % 4 == 0) or \
               (coh and block_readable(N, w, sup_coh) and d == N)
        if pred:
            pred_alive.add(name)
        if alive:
            hist[d] = hist.get(d, 0) + 1
    alive = {r[0] for r in rows if r[1]}
    check(label, pred_alive == alive and hist == expect_hist,
          f"alive {len(alive)}/{len(rows)}, collapsed-rule predicted {len(pred_alive)}, "
          f"alive-by-degree {dict(sorted(hist.items()))}")
    return rows


def main():
    print("T1  the rule is exact, free world (XY + field), full censuses")
    one3, one4, one5 = [0.0, 0.0, 0.5], [0.0] * 3 + [0.5], [0.0] * 4 + [0.5]
    r3pop = run_config(3, 1, None, one3, 0.0, "N=3 s=1 population")
    run_config(3, 1, "real", one3, 0.0, "N=3 s=1 real coherence")
    run_config(3, 1, "imag", one3, 0.0, "N=3 s=1 imaginary coherence")
    run_config(3, 2, "real", one3, 0.0, "N=3 s=2 real coherence")
    r4pop = run_config(4, 1, None, one4, 0.0, "N=4 s=1 population")
    run_config(4, 1, "real", one4, 0.0, "N=4 s=1 real coherence")
    run_config(5, 1, "real", one5, 0.0, "N=5 s=1 real coherence")
    run_config(5, 3, None, one5, 0.0, "N=5 s=3 (popcount 2) population")
    r5coh = run_config(5, 3, "real", one5, 0.0, "N=5 s=3 (popcount 2) real coherence")
    run_config(5, 1, "real", [0.3, 0.1, 0.5, 0.2, 0.4], 0.0,
               "N=5 s=1 real coherence, FULL non-uniform watching")

    print("T2  layer attribution: each layer's own kill, witnessed")
    k_kill = {r[0] for r in r4pop if not r[3] and r[4] and r[5]}
    pure_xy4 = {"".join(p) for p in product("XY", repeat=4)}
    check("K alone kills the 16 pure-XY strings at N=4 population",
          k_kill == pure_xy4 and all(not r[1] for r in r4pop if r[0] in pure_xy4),
          f"{len(k_kill)} strings, all dead, all w=4 (no popcount-2 block populated)")
    f_kill = {r[0] for r in r5coh if r[3] and r[4] and not r[5]}
    f_degrees = {r[6] for r in r5coh if r[0] in f_kill}
    check("F alone kills the 20 fermion-linear strings at N=5 s=3 coherence",
          len(f_kill) == 20 and f_degrees == {1, 9}
          and all(not r[1] for r in r5coh if r[0] in f_kill),
          f"{len(f_kill)} strings, degrees {sorted(f_degrees)} = the a_k and their "
          "Gamma-duals, all dead")
    zii = next(r for r in r3pop if r[0] == "ZII")
    check("V alone kills the single-site Z at N=3 population (K and F both allow)",
          (not zii[1]) and zii[3] and (not zii[4]) and zii[5] and zii[6] == 2,
          "ZII: k_ok, d=2 (F allows), eps=-1 (V kills), dead: the doubly-mirrored "
          "zero of the h thread, now one cell of the closed rule")

    print("T3  the mod-4 shadow: eps is a pure function of Majorana degree, N-free")
    bad, total = t4_mod4_identity(6)
    check("eps_letters == eps_degree == matrix conjugation, all strings N=2..6",
          bad == 0, f"{total} (string, gauge) pairs, {bad} mismatches; "
          "eps_odd=(-1)^(d(d-1)/2), eps_even=(-1)^(d(d+1)/2)")

    print("T4  the collapsed per-sector rule (V_g always a symmetry; the prep splits "
          "into V-eigen-sectors; at d=N the coherence channel's signs are automatic)")
    one6 = [0.0] * 5 + [0.5]
    run_collapsed(5, 3, "real", one5, "N=5 s=3 coherence via the collapsed rule",
                  {4: 210, 5: 252, 8: 45})
    run_collapsed(6, 1, None, one6, "N=6 population via the collapsed rule",
                  {4: 255, 8: 255, 12: 1})
    r6 = run_collapsed(6, 1, "real", one6,
                       "N=6 coherence via the collapsed rule (N=2 mod 4: no gauge "
                       "stabilizes globally, the per-sector form still exact)",
                       {4: 255, 6: 544, 8: 255, 12: 1})
    # the honest divergence pin: the GLOBAL-stabilizer form of the rule (T1's) is
    # blind here (stab=[] with coherence at N=2 mod 4) and over-predicts
    global_pred = {r[0] for r in r6 if r[3] and r[4] and r[5]}
    alive6 = {r[0] for r in r6 if r[1]}
    check("the global-stabilizer form genuinely diverges at N=6 coherence",
          alive6 < global_pred,
          f"global form predicts {len(global_pred)} alive, actual {len(alive6)}: "
          "the per-sector refinement is load-bearing, not cosmetic")
    # binomial anatomy of the histograms (comb ties, not literals)
    from math import comb
    check("alive-by-degree counts are binomials cut by kinematics",
          210 == comb(10, 4) and 252 == comb(10, 5) and 45 == comb(10, 8)
          and 255 == comb(12, 4) - comb(6, 4) * 2 ** 4
          and 544 == comb(12, 6) - comb(6, 3) - comb(6, 2) * 4 * comb(4, 2)
          and 1 == comb(12, 12),
          "N=5 s=3: C(10,4)+C(10,5)+C(10,8); N=6: C(12,4)-240 (w=4 caps unreadable), "
          "C(12,6)-20-360 (w<4 cannot reach the coherence blocks), C(12,12)=Z^N")

    print("T5  the zz boundary: V and F are properties of the free world")
    rows3, _, _ = census(3, 1, None, one3, 0.7)
    dead3 = {r[0] for r in rows3 if not r[1]}
    k_dead3 = {r[0] for r in rows3 if not r[3]}
    check("N=3 population, zz=0.7: dead = layer K + the conserved Z^N",
          dead3 == k_dead3 | {"ZZZ"},
          f"{len(dead3)} dead; V's and F's kills revive, only ZZZ stays "
          "(parity conserved, starts at 0 for the odd-N pair)")
    rows4, _, _ = census(4, 1, "real", one4, 0.7)
    dead4 = {r[0] for r in rows4 if not r[1]}
    k_dead4 = {r[0] for r in rows4 if not r[3]}
    check("N=4 real coherence, zz=0.7: dead = exactly layer K",
          dead4 == k_dead4, f"{len(dead4)} dead, all odd-w (Z^N starts nonzero "
          "at even N, so no parity zero)")

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: the dead set is closed by the three conserved structures "
          "(popcount blocks, the doubly-mirrored kill, fermion degree), and the "
          "sharpened per-sector form collapses them to one fermionic line: alive iff "
          "(K_pop and d = 0 mod 4) or (coherence and K_coh and d = N); with zz the "
          "V and F layers collapse and only the blocks plus conserved parity remain.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
