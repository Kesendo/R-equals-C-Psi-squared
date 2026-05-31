#!/usr/bin/env python3
"""Two excitations launched from opposite ends of a chain: do they meet, and is the meeting
a pass-through or a collision?

Setup: N=6 chain (two N=3 halves), one excitation on site 0 and one on site N-1, the
popcount-2 sector. Closed system first (gamma=0) so the collision is not blurred by
forgetting. We compare two Hamiltonians:

  XY        H = sum_b (X_b X_{b+1} + Y_b Y_{b+1})            -> free fermions (Jordan-Wigner),
                                                                excitations are transparent
  Heisenberg H = XY + Delta * sum_b Z_b Z_{b+1}              -> the ZZ term is a contact
                                                                interaction, excitations scatter

Read per-site occupation <n_l>(t), plus two discriminators:
  - end revival: max <n_0>(t) for t>0  (did the left excitation bounce back to its own end?)
  - center pile-up: max of <n_{c} n_{c+1}>(t), the joint two-site density at the middle bond
    (do the two excitations preferentially co-locate -> a bound/scattered pair?)

By the chain's mirror symmetry the state stays mirror-symmetric, so any meeting happens at
the center. The question is what KIND of meeting. The numbers tell it.
"""
import sys

import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
BLOCKS = "▁▂▃▄▅▆▇█"


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def bond(N, b, P, Qop):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == b else (Qop if i == b + 1 else I2))
    return o


def H_xy(N):
    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        H += bond(N, b, X, X) + bond(N, b, Y, Y)
    return H


def H_heisenberg(N, delta=1.0):
    H = H_xy(N)
    for b in range(N - 1):
        H += delta * bond(N, b, Z, Z)
    return H


def sparkline(ys):
    ys = np.asarray(ys, float)
    lo, hi = ys.min(), ys.max()
    if hi - lo < 1e-12:
        return BLOCKS[0] * len(ys)
    idx = np.clip(((ys - lo) / (hi - lo) * 7).round().astype(int), 0, 7)
    return "".join(BLOCKS[i] for i in idx)


def run(N, H, label, ts):
    d = 2 ** N
    # |1 0 ... 0 1>: excitation on site 0 and site N-1 (popcount 2)
    psi0 = np.array([1], complex)
    for i in range(N):
        bit = 1 if i in (0, N - 1) else 0
        psi0 = np.kron(psi0, np.array([0, 1], complex) if bit else np.array([1, 0], complex))
    # closed-system unitary evolution via eigendecomposition of H
    w, V = np.linalg.eigh(H)
    c = V.conj().T @ psi0
    nl = [(np.eye(d) - op_at(N, l, Z)) / 2.0 for l in range(N)]
    cmid = N // 2 - 1  # center bond (cmid, cmid+1)
    pair_center = nl[cmid] @ nl[cmid + 1]

    occ = {l: [] for l in range(N)}
    n0_series, paircenter_series = [], []
    for t in ts:
        psit = V @ (np.exp(-1j * w * t) * c)
        rho_diag = np.abs(psit) ** 2  # not needed directly; use expectations
        for l in range(N):
            occ[l].append(float((psit.conj() @ (nl[l] @ psit)).real))
        n0_series.append(occ[0][-1])
        paircenter_series.append(float((psit.conj() @ (pair_center @ psit)).real))

    end_revival = max(n0_series[1:])
    center_pileup = max(paircenter_series)
    print(f"\n  {label}")
    for l in range(N):
        ys = occ[l]
        print(f"    site {l}  {sparkline(ys)}  {ys[0]:.2f}->{ys[-1]:.2f}  peak {max(ys):.2f}")
    print(f"    end revival  max<n_0>(t>0) = {end_revival:.3f}")
    print(f"    center pile-up  max<n_{cmid} n_{cmid+1}>(t) = {center_pileup:.3f}"
          f"   (independent sites would give ~<n_{cmid}>*<n_{cmid+1}>)")
    return end_revival, center_pileup


def sweep_delta(N, deltas, ts):
    """For each Delta = ZZ coefficient, track the center bond and report the peak center
    pile-up and the peak connected correlation  C = <n_c n_{c+1}> - <n_c><n_{c+1}>
    (~0 for genuinely free particles; > 0 attraction/binding; < 0 repulsion)."""
    d = 2 ** N
    cmid = N // 2 - 1
    nc, nc1 = (np.eye(d) - op_at(N, cmid, Z)) / 2.0, (np.eye(d) - op_at(N, cmid + 1, Z)) / 2.0
    pair = nc @ nc1
    psi0 = np.array([1], complex)
    for i in range(N):
        bit = 1 if i in (0, N - 1) else 0
        psi0 = np.kron(psi0, np.array([0, 1], complex) if bit else np.array([1, 0], complex))

    print(f"Delta-sweep, two excitations from opposite ends, closed system, N={N}, center bond ({cmid},{cmid+1})")
    print(f"  Delta=0 is the free XY reference; connected C isolates the genuine two-body effect")
    print(f"\n  {'Delta':>6}   {'peak pile-up':>12}   {'peak C(+)':>10}   {'peak C(-)':>10}   reading")
    for delta in deltas:
        H = H_xy(N) + (delta * sum(bond(N, b, Z, Z) for b in range(N - 1)) if delta != 0 else 0)
        w, V = np.linalg.eigh(H)
        c = V.conj().T @ psi0
        pile, cmax, cmin = [], [], []
        for t in ts:
            p = V @ (np.exp(-1j * w * t) * c)
            e_pair = float((p.conj() @ (pair @ p)).real)
            e_c = float((p.conj() @ (nc @ p)).real)
            e_c1 = float((p.conj() @ (nc1 @ p)).real)
            conn = e_pair - e_c * e_c1
            pile.append(e_pair)
            cmax.append(conn)
            cmin.append(conn)
        peak_pile, peak_cp, peak_cm = max(pile), max(cmax), min(cmin)
        if peak_cp > 0.02 and peak_cp >= -peak_cm:
            tag = "attract (bind)"
        elif peak_cm < -0.02:
            tag = "repel"
        else:
            tag = "~free"
        print(f"  {delta:>6.2f}   {peak_pile:>12.3f}   {peak_cp:>10.3f}   {peak_cm:>10.3f}   {tag}")


def liouvillian(N, H, gamma):
    """Open-system Liouvillian L (F-order vec): -i[H, .] + gamma * sum_l (Z_l . Z_l - .)."""
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += gamma * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def open_meeting(N, gammas, ts):
    """Two excitations from opposite ends under dephasing: does the interaction signature
    (XY vs Heisenberg) survive forgetting, while both flow into 2/N? Tracks the peak center
    connected correlation |C| over the transient and the final per-site occupation."""
    d = 2 ** N
    cmid = N // 2 - 1
    nl = [(np.eye(d) - op_at(N, l, Z)) / 2.0 for l in range(N)]
    psi0 = np.array([1], complex)
    for i in range(N):
        bit = 1 if i in (0, N - 1) else 0
        psi0 = np.kron(psi0, np.array([0, 1], complex) if bit else np.array([1, 0], complex))
    rho0 = np.outer(psi0, psi0.conj())
    vec0 = rho0.flatten(order="F")

    def trajectory(H, gamma):
        L = liouvillian(N, H, gamma)
        w, V = np.linalg.eig(L)
        c = np.linalg.solve(V, vec0)
        peak_absC, occ_last = 0.0, None
        for t in ts:
            rho = (V @ (np.exp(w * t) * c)).reshape(d, d, order="F")
            es = [float(np.trace(nl[l] @ rho).real) for l in range(N)]
            epair = float(np.trace(nl[cmid] @ nl[cmid + 1] @ rho).real)
            conn = epair - es[cmid] * es[cmid + 1]
            peak_absC = max(peak_absC, abs(conn))
            occ_last = es
        return peak_absC, occ_last

    target = 2.0 / N
    print(f"Open-system meeting, N={N}, two excitations on sites 0 and {N-1}, center bond ({cmid},{cmid+1})")
    print(f"  flow target (open): 2/N = {target:.4f} per site; does the XY-vs-Heisenberg meeting survive forgetting?")
    print(f"\n  {'gamma':>6}   {'peak|C| XY':>11}   {'peak|C| Heis':>13}   {'ratio H/XY':>10}   {'sites->':>8} (XY final, Heis final)")
    for g in gammas:
        cxy, oxy = trajectory(H_xy(N), g)
        chs, ohs = trajectory(H_heisenberg(N, 1.0), g)
        ratio = chs / cxy if cxy > 1e-9 else float("nan")
        oxy_s = " ".join(f"{x:.2f}" for x in oxy)
        ohs_s = " ".join(f"{x:.2f}" for x in ohs)
        print(f"  {g:>6.2f}   {cxy:>11.4f}   {chs:>13.4f}   {ratio:>10.2f}   XY[{oxy_s}]  Heis[{ohs_s}]")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "open":
        N = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        ts = np.linspace(0.0, 40.0, 200)
        open_meeting(N, [0.0, 0.05, 0.2, 0.5, 1.0], ts)
        return
    if len(sys.argv) > 1 and sys.argv[1] == "sweep":
        N = int(sys.argv[2]) if len(sys.argv) > 2 else 6
        ts = np.linspace(0.0, 6.0, 120)
        sweep_delta(N, [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0], ts)
        return
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    ts = np.linspace(0.0, 4.0, 80)
    print(f"Two excitations from opposite ends, closed system (gamma=0), N={N}")
    print(f"  initial |1 0...0 1>: site 0 and site {N-1}; watch them meet at the center")
    print(f"  target if it were open + 1/N-flow: 2/N = {2.0/N:.4f} per site (here gamma=0, no relaxation)")
    rv_xy, cp_xy = run(N, H_xy(N), "XY  (free fermions, expect transparent pass-through)", ts)
    rv_h, cp_h = run(N, H_heisenberg(N, 1.0), "Heisenberg (XY + ZZ, expect interaction)", ts)
    print("\n  DISCRIMINATOR (XY vs Heisenberg):")
    print(f"    end revival : XY {rv_xy:.3f}  vs  Heis {rv_h:.3f}   (higher Heis = reflection back to own end)")
    print(f"    center pile : XY {cp_xy:.3f}  vs  Heis {cp_h:.3f}   (higher Heis = bound/scattered pair in the middle)")


if __name__ == "__main__":
    main()
