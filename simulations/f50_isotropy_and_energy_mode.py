"""F50 Door 1, its shape: the SU(2) isotropy floor 4/3 in the full space, and what the sector projector costs.

The open statement of docs/proofs/PROOF_WEIGHT1_DEGENERACY.md ("The count as a plane crossing", the Hamiltonian end):
on the Heisenberg chain, inside block (p, p), every traceless function X = f(H_p) of the sector Hamiltonian has height
h(X) = sum_l ||[n_l, X]||^2 / ||X||^2 > 1, i.e. mu_2(W_p) < p - 1/2. This script pins what IS exact around it.

E1  Isotropy (exact, every graph, every N). For an operator A that commutes with the global SU(2) (every function of
    H and S^2 of a Heisenberg Hamiltonian on any graph), the three dephasing forms agree, h_X(A) = h_Y(A) = h_Z(A)
    with h_alpha(A) = sum_l ||[sigma^alpha_l, A]||^2 / (4 ||A||^2), since the global rotation carrying Z to X fixes A;
    and for EVERY operator their sum is 2 <A, w A> / ||A||^2 with w the Pauli weight (a non-identity letter at a site
    anticommutes with two of the three Paulis there, each contributing 4/4, and distinct strings stay distinct so the
    cross terms vanish). So h_Z(A) = (2/3) <w>_A, and for traceless invariant A, whose weight components are each
    invariant (local rotations keep the weight grading) and whose weight-1 component is therefore zero (a vector
    a_l fixed by every rotation is zero), h_Z(A) >= 4/3 with equality exactly when A is purely weight 2. EXACT route
    (the repo's case 1): A = d * B - tr(B) * I for B in {H, S^2, H^2, H S^2, H^3} is an integer matrix (Y_l = i X_l Z_l,
    so ||[Y_l, A]|| = ||[X_l Z_l, A]|| over the integers), every norm an integer, h_alpha and <w> Fractions, compared
    with ==, N = 3..6. Beside it the float reading on random traceless invariants (deviations printed, not gated) and
    the control: a traceless function of (H, S_z) that is NOT SU(2)-invariant breaks h_X = h_Z (the sum identity holds
    for every operator and needs no control).
E2  The minimum of h_Z over the invariant commutant (the span of the joint (E, S) eigenprojectors, which is the
    invariant commutant when no multiplet is degenerate, true on the chain at N <= 7 where the projector count equals
    the multiplet count) is 4/3 at N = 3..7, twice on the CHAIN: H - <H> and S^2 - <S^2> are the two weight-2 invariants
    among the functions of (H, S^2) there (on the complete graph H = 2 S^2 - 3N/2 and the multiplicity is one). Read
    against the generalized eigensolver's resolution.
E3  The sector's energy mode, exact in rationals: for X = H_p - <H_p> in block (p, p) the height is the Rayleigh quotient
    2 ||H_XY||_p^2 / (||H_XY||_p^2 + Var_p(H_ZZ)), both terms integers / rationals over the sector configurations
    (||H_XY||_p^2 = 8 (N-1) C(N-2, p-1), Var the SUM of squared deviations), and it equals
    4 N (N - 1) / (3 N^2 - 4 N - (N - 2p)^2) at every p, compared with == in Fraction arithmetic for N = 3..12 and
    every 1 <= p <= N - 1 (read, not derived); at fixed filling x = p/N the limit is 4 / (3 - (1 - 2x)^2), the isotropy
    floor 4/3 exactly at half filling.
    Its limit 4/3 is the full-space floor of E1; the sector projector P_p costs the difference.
E4  The energy mode is a Ritz vector, so h(H_p - <H_p>) is an UPPER bound on the minimum height 2 (p - mu_2); the
    reading: the gap times N and the overlap of the minimizer with the energy mode, N = 4..12 at half filling.
E5  The same matrix on the Delta axis: the golden-rule rate matrix R of simulations/xxz_delta_star_descent.py
    (experiments/XXZ_AXIS_BANDEDGE_TO_LEBENSADER.md, Delta* <=> gap(R) = 2) equals 4 W_p - 4 p I entry for entry, two
    routes to one matrix (rounding, ratio gated), so gap(R) = 2 h_min and the handover Delta* is mu_2 = p - 1/2 moved
    along Delta. READING beside it, ungated: the two lowest heights and the lowest mode's overlap with the energy mode
    on Delta = 1.00..1.32 in steps of 0.02 for N = 6..12, printed at five of the points with the minimum gap over the
    grid. What it reads to N = 12: the lowest height rises with Delta at first (N >= 8) and then falls; the gap to the
    second height does not close but shrinks with N (0.283, 0.244, 0.214, 0.195), an avoided crossing with a closing
    trend; the overlap rises to 0.99 and then swaps within about 0.06 in Delta; the fall with N shows from Delta = 1.2.
Prints ALL GATES PASS when E1, E2, E3, E4 and E5 hold. About one minute.
"""
import sys
from fractions import Fraction
from math import comb
import numpy as np
from scipy.linalg import eigh as geigh

EPS = np.finfo(float).eps
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
id2 = np.eye(2, dtype=complex)

ok_all = True


def gate(cond, msg):
    global ok_all
    ok_all &= bool(cond)
    print(("  PASS  " if cond else "  FAIL  ") + msg)


def kron_all(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def site_op(N, l, op):
    return kron_all([op if m == l else id2 for m in range(N)])


def full_heisenberg(N, bonds):
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for (l, m) in bonds:
        for s in (sx, sy, sz):
            H += site_op(N, l, s) @ site_op(N, m, s)
    return H


def h_alpha(A, ops):
    nrm = np.real(np.vdot(A, A))
    return sum(np.real(np.vdot(O @ A - A @ O, O @ A - A @ O)) for O in ops) / 4 / nrm


def weight_mean(A, N):
    """<A, w A> / ||A||^2 by partial traces: ||A||^2 - ||(1/2) tr_l A (x) I||^2 = ||A||^2 - (1/2) ||tr_l A||^2 per site."""
    nrm = np.real(np.vdot(A, A))
    T = A.reshape([2] * N + [2] * N)
    tot = 0.0
    for l in range(N):
        Tl = np.trace(T, axis1=l, axis2=N + l)
        tot += nrm - 0.5 * np.real(np.vdot(Tl, Tl))
    return tot / nrm


def joint_projectors(N, H, S2, Sz=None):
    """Projectors onto the joint eigenspaces of H and S^2 (and S_z if given)."""
    mix = H + np.pi * S2 + (np.e * Sz if Sz is not None else 0)
    E, V = np.linalg.eigh(mix)
    labels = []
    for a in range(len(E)):
        v = V[:, a]
        lab = [round(float(np.real(v.conj() @ H @ v)), 7), round(float(np.real(v.conj() @ S2 @ v)), 7)]
        if Sz is not None:
            lab.append(round(float(np.real(v.conj() @ Sz @ v)), 7))
        labels.append(tuple(lab))
    groups = {}
    for a, lab in enumerate(labels):
        groups.setdefault(lab, []).append(a)
    return [sum(np.outer(V[:, a], V[:, a].conj()) for a in idxs) for idxs in groups.values()]


def int_site_op(N, l, op):
    out = np.array([[1]], dtype=object)
    for m in range(N):
        out = np.kron(out, op if m == l else np.eye(2, dtype=object))
    return out


def int_norm2(A):
    return int(np.sum(A * A))


def e1_exact(N):
    """h_X = h_Y = h_Z = (2/3)<w> as Fractions on integer invariant matrices, compared with ==."""
    X = np.array([[0, 1], [1, 0]], dtype=object)
    Z = np.array([[1, 0], [0, -1]], dtype=object)
    XZ = X @ Z
    I2 = np.eye(2, dtype=object)
    d = 2 ** N
    Id = np.eye(d, dtype=object)
    H = sum(int_site_op(N, l, s) @ int_site_op(N, l + 1, s) for l in range(N - 1) for s in (X, Z)) + \
        sum(-(int_site_op(N, l, XZ) @ int_site_op(N, l + 1, XZ)) for l in range(N - 1))   # Y Y = (i XZ)(i XZ) = -XZ XZ
    S2x4 = sum((sum(int_site_op(N, l, s) for l in range(N))) @ (sum(int_site_op(N, l, s) for l in range(N))) for s in (X, Z)) + \
        (-(sum(int_site_op(N, l, XZ) for l in range(N))) @ (sum(int_site_op(N, l, XZ) for l in range(N))))   # 4 S^2
    ops = {"X": [int_site_op(N, l, X) for l in range(N)], "Y": [int_site_op(N, l, XZ) for l in range(N)],
           "Z": [int_site_op(N, l, Z) for l in range(N)]}
    worst_float = 0.0
    for name, B in (("H", H), ("S2", S2x4), ("H^2", H @ H), ("H S2", H @ S2x4), ("H^3", H @ H @ H)):
        A = d * B - int(np.trace(B)) * Id
        n2 = int_norm2(A)
        h = {}
        for a, lst in ops.items():
            h[a] = Fraction(sum(int_norm2(O @ A - A @ O) for O in lst), 4 * n2)
        # <w> by partial traces, exact: sum_l (||A||^2 - (1/2)||tr_l A||^2)
        T = A.reshape([2] * N + [2] * N)
        wsum = Fraction(0)
        for l in range(N):
            Tl = np.trace(T, axis1=l, axis2=N + l)
            wsum += Fraction(n2) - Fraction(int(np.sum(Tl * Tl)), 2)
        wmean = wsum / n2
        gate(h["X"] == h["Y"] == h["Z"] and h["X"] + h["Y"] + h["Z"] == 2 * wmean and h["Z"] >= Fraction(4, 3),
             f"E1 N={N} exact: B={name:4s} h_X = h_Y = h_Z = {h['Z']} = (2/3)<w> with <w> = {wmean}, >= 4/3")
        if name in ("H", "S2"):
            gate(h["Z"] == Fraction(4, 3), f"E1 N={N} exact: B={name:4s} h_Z == 4/3 (weight 2)")


def e1_e2(N):
    bonds = [(l, l + 1) for l in range(N - 1)]
    H = full_heisenberg(N, bonds)
    Sa = [sum(site_op(N, l, s) for l in range(N)) / 2 for s in (sx, sy, sz)]
    S2 = sum(s @ s for s in Sa)
    Sz = Sa[2] * 1.0
    Xs = [site_op(N, l, sx) for l in range(N)]
    Ys = [site_op(N, l, sy) for l in range(N)]
    Zs = [site_op(N, l, sz) for l in range(N)]
    d = 2 ** N
    I = np.eye(d)
    projs = joint_projectors(N, H, S2)
    rng = np.random.default_rng(N)
    worst_iso, worst_w, hmin = 0.0, 0.0, 9.0
    for _ in range(12):
        c = rng.normal(size=len(projs))
        A = sum(ci * Pi for ci, Pi in zip(c, projs))
        A -= np.trace(A) / d * I
        hx, hy, hz = h_alpha(A, Xs), h_alpha(A, Ys), h_alpha(A, Zs)
        worst_iso = max(worst_iso, abs(hx - hz), abs(hy - hz))
        worst_w = max(worst_w, abs(hx + hy + hz - 2 * weight_mean(A, N)))
        hmin = min(hmin, hz)
    scale = EPS * d
    print(f"  read  E1 N={N}: float, {len(projs)} invariant projectors, 12 random traceless: max |h_X - h_Z| {worst_iso:.1e}, "
          f"max |sum - 2<w>| {worst_w:.1e}, min h_Z {hmin:.4f} (the floor is E2's, not this sample's)")
    # negative control: a traceless function of (H, S_z) that is not invariant
    projs_z = joint_projectors(N, H, S2, Sz)
    c = rng.normal(size=len(projs_z))
    A = sum(ci * Pi for ci, Pi in zip(c, projs_z))
    A -= np.trace(A) / d * I
    dev = abs(h_alpha(A, Xs) - h_alpha(A, Zs))
    gate(dev > 1e-3, f"E1 N={N}: control, a traceless function of (H, S_z) breaks h_X = h_Z by {dev:.4f}")
    # E2: the minimum over the invariant commutant, generalized eigenproblem on the projector span
    G = np.array([[np.real(np.trace(Pi @ Pj)) for Pj in projs] for Pi in projs])
    Kz = np.array([[sum(np.real(np.vdot(Zl @ Pi - Pi @ Zl, Zl @ Pj - Pj @ Zl)) for Zl in Zs) / 4 for Pj in projs] for Pi in projs])
    vals = np.sort(geigh(Kz, G, eigvals_only=True))
    gate(abs(vals[0]) < 64 * scale and abs(vals[1] - 4 / 3) < 64 * scale and abs(vals[2] - 4 / 3) < 64 * scale
         and (len(vals) < 4 or vals[3] > 4 / 3 + 1e-6),
         f"E2 N={N}: spectrum of K_Z on the invariant commutant starts 0, 4/3, 4/3, then {vals[3] if len(vals) > 3 else float('nan'):.6f}: "
         f"the floor 4/3 attained exactly twice (H and S^2), deviations {abs(vals[1]-4/3):.1e}, {abs(vals[2]-4/3):.1e}")
    hH = h_alpha(H - np.trace(H) / d * I, Zs)
    hS = h_alpha(S2 - np.trace(S2) / d * I, Zs)
    gate(abs(hH - 4 / 3) < 64 * scale and abs(hS - 4 / 3) < 64 * scale, f"E2 N={N}: h_Z(H - <H>) = {hH:.12f}, h_Z(S^2 - <S^2>) = {hS:.12f}")


# ---- the sector (functions copied from f50_hamiltonian_end_heights.py; importing it would run the whole gate) ----
def sector(N, p):
    confs = [c for c in range(2 ** N) if bin(c).count("1") == p]
    return confs, {c: i for i, c in enumerate(confs)}


def h_sector(N, p, bonds, delta=1.0):
    confs, idx = sector(N, p)
    d = len(confs)
    H = np.zeros((d, d))
    for i, c in enumerate(confs):
        for (l, m) in bonds:
            b1, b2 = (c >> l) & 1, (c >> m) & 1
            H[i, i] += delta * (1.0 if b1 == b2 else -1.0)
            if b1 != b2:
                H[idx[c ^ (1 << l) ^ (1 << m)], i] += 2.0
    return confs, H


def w_closed(N, V, confs):
    d = len(confs)
    W = np.zeros((d, d))
    for l in range(N):
        bit = np.array([(c >> l) & 1 for c in confs], dtype=float)
        n_l = V.T @ (bit[:, None] * V)
        W += n_l * n_l
    return W


def energy_mode_exact(N, p):
    """2 ||H_XY||^2 / (||H_XY||^2 + Var(H_ZZ)) over the sector, exact in Fraction (Pauli book, J = 1)."""
    confs, _ = sector(N, p)
    d = len(confs)
    bonds = [(l, l + 1) for l in range(N - 1)]
    nxy = 0
    diag = []
    for c in confs:
        s = 0
        for (l, m) in bonds:
            if ((c >> l) & 1) != ((c >> m) & 1):
                nxy += 4      # the hop has matrix element 2, and each unequal bond gives one entry of 2 in the column
                s -= 1
            else:
                s += 1
        diag.append(s)
    mean = Fraction(sum(diag), d)
    var = sum((Fraction(x) - mean) ** 2 for x in diag)
    nxy = Fraction(nxy)
    return 2 * nxy / (nxy + var)


def e3_e4(nmax):
    bad = 0
    for N in range(3, nmax + 1):
        for p in range(1, N):
            if energy_mode_exact(N, p) != Fraction(4 * N * (N - 1), 3 * N * N - 4 * N - (N - 2 * p) ** 2):
                bad += 1
    gate(bad == 0, f"E3: h(H_p - <H_p>) == 4N(N-1)/(3N^2-4N-(N-2p)^2) in Fractions, every N = 3..{nmax}, every p = 1..N-1 ({bad} failures)")
    print("E3/E4  N  p   dim   h_min     h(E) exact        h(E)      gap*N   |<u2|E>|")
    for N in range(4, nmax + 1):
        p = N // 2
        bonds = [(l, l + 1) for l in range(N - 1)]
        confs, H = h_sector(N, p, bonds)
        d = len(confs)
        E, V = np.linalg.eigh(H)
        W = w_closed(N, V, confs)
        mu, U = np.linalg.eigh(W)
        order = np.argsort(mu)[::-1]
        mu, U = mu[order], U[:, order]
        hmin = 2 * (p - mu[1])
        c = E - E.mean()
        c /= np.linalg.norm(c)
        hE = 2 * (p - c @ W @ c)
        ov = abs(U[:, 1] @ c)
        ex = energy_mode_exact(N, p)
        print(f"      {N:2d} {p:2d} {d:5d}  {hmin:.5f}  {str(ex):>14}  {hE:.6f}  {(hE-hmin)*N:.4f}  {ov:.4f}")
        gate(abs(float(ex) - hE) < 64 * EPS * d, f"E3 N={N}: exact Rayleigh quotient vs the compression W_p, deviation {abs(float(ex)-hE):.1e}")
        closed = Fraction(4 * N * (N - 1), 3 * N * N - 4 * N - (N - 2 * p) ** 2)
        gate(ex == closed, f"E3 N={N} p={p}: h(E) == {closed} exactly")
        gate(hE > hmin - 1e-12, f"E4 N={N}: the energy mode is a Ritz value, h(E) >= h_min")


def rate_matrix_R(N, p, Delta):
    """Copied in spirit from simulations/xxz_delta_star_descent.py (gain sum_k |<a|Z_k|b>|^2, loss -4 sum_k Var_a(n_k))."""
    bonds = [(l, l + 1) for l in range(N - 1)]
    confs, H = h_sector(N, p, bonds, Delta)
    E, V = np.linalg.eigh(H)
    ns = len(confs)
    nk = np.array([[(s >> k) & 1 for k in range(N)] for s in confs], dtype=float)
    R = np.zeros((ns, ns))
    for k in range(N):
        zk = 1.0 - 2.0 * nk[:, k]
        Mk = V.T @ (zk[:, None] * V)
        R += Mk ** 2
    w = V ** 2
    mean_n = w.T @ nk
    var = (mean_n - mean_n ** 2).sum(axis=1)
    R[np.diag_indices(ns)] = -4.0 * var
    return R, V, confs


def e5():
    for N, Delta in ((6, 1.0), (7, 1.3), (8, 1.0)):
        p = (N + 1) // 2
        R, V, confs = rate_matrix_R(N, p, Delta)
        W = w_closed(N, V, confs)
        M = 4 * W - 4 * p * np.eye(len(confs))
        dev = np.max(np.abs(R - M))
        scale = EPS * len(confs) * 4 * p
        gate(dev < 64 * scale, f"E5 N={N} Delta={Delta}: R (xxz_delta_star_descent) = 4 W_p - 4p I entry-wise, max deviation {dev:.1e} = {dev/scale:.1f} eps*d*4p")
    print("E5 reading (ungated): half filling, Delta grid 1.00..1.32 step 0.02; h1 = gap(R)/2, h2 the next height,")
    print("   ov = |<u1|E - Ebar>|; five points shown, then the minimum of h2 - h1 over the grid and where")
    shown = [1.0, 1.1, 1.14, 1.2, 1.3]
    grid = [round(1.0 + 0.02 * k, 2) for k in range(17)]
    print("      N  " + "   ".join(f"D={D:<4} h1/h2/ov     " for D in shown) + "   min gap @ D")
    for N in (6, 8, 10, 12):
        p = N // 2
        bonds = [(l, l + 1) for l in range(N - 1)]
        rows = {}
        for D in grid:
            confs, H = h_sector(N, p, bonds, D)
            E, V = np.linalg.eigh(H)
            mu, U = np.linalg.eigh(w_closed(N, V, confs))
            o = np.argsort(mu)[::-1]
            mu, U = mu[o], U[:, o]
            c = E - E.mean()
            c /= np.linalg.norm(c)
            rows[D] = (2 * (p - mu[1]), 2 * (p - mu[2]), abs(U[:, 1] @ c))
        gmin = min(grid, key=lambda D: rows[D][1] - rows[D][0])
        cells = [f"{rows[D][0]:.4f}/{rows[D][1]:.4f}/{rows[D][2]:.2f}" for D in shown]
        print(f"     {N:2d}  " + "   ".join(f"{cl:22s}" for cl in cells) + f"   {rows[gmin][1]-rows[gmin][0]:.3f} @ {gmin:.2f}")


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    print("E1/E2  the full space, Heisenberg chain, N = 3..7")
    for N in range(3, 7):
        e1_exact(N)
    for N in range(3, 8):
        e1_e2(N)
    e3_e4(nmax)
    e5()
    print("\nALL GATES PASS" if ok_all else "\nSOME GATE FAILED")
