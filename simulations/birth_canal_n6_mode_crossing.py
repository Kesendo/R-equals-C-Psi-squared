#!/usr/bin/env python3
"""Why does the (1,0) block miss the full 4^6 global slowest for the N=6 deep-edge profile?

Part C of birth_canal_vacuum_block_verifier.py flagged a 3.5e-01 full-vs-block gap at N=6 deep-edge.
This pins which block carries the full global-slowest mode by decomposing its eigenvector by
(popcount_ket, popcount_bra). FINDING (refining the initial "a higher ODD (2,1) block" guess): at
low Q the global slowest is the EVEN {0,2}-coherence in the 2-excitation DENSITY block (n_diff hist
{0:0.78, 2:0.22}, <n_XY>=0.44 < 1) - the coherence-horizon mode - NOT a number-changing block.
"""
import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)


def op_at(N, s, P):
    o = np.array([[1]], complex)
    for i in range(N):
        o = np.kron(o, P if i == s else I2)
    return o


def H_xy_unit(N):
    H = np.zeros((2 ** N, 2 ** N), complex)
    for b in range(N - 1):
        for P in (X, Y):
            t = np.array([[1]], complex)
            for i in range(N):
                t = np.kron(t, P if i in (b, b + 1) else I2)
            H += t
    return H


def full_L(N, Q, profile):
    d = 2 ** N
    Id = np.eye(d)
    H1 = H_xy_unit(N)
    L = -1j * Q * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    return L


def block_slowest(N, Q, profile):
    h = np.zeros((N, N), complex)
    for i in range(N - 1):
        h[i, i + 1] = 2.0
        h[i + 1, i] = 2.0
    w = np.linalg.eigvals(-1j * Q * h - 2.0 * np.diag(profile))
    nz = w[np.abs(w) > 1e-7]
    return -float(np.max(nz.real))


def analyze(N, Q, profile, label):
    d = 2 ** N
    L = full_L(N, Q, profile)
    w, V = np.linalg.eig(L)
    # global slowest non-kernel
    mask = np.abs(w) > 1e-7
    idx = np.where(mask)[0]
    slow = idx[np.argmax(w[idx].real)]
    lam = w[slow]
    v = V[:, slow]
    rho = v.reshape(d, d)                      # row-major vec: rho[a,b] = vec[a*d+b]
    wts = np.abs(rho) ** 2
    wts /= wts.sum()
    pc = np.array([bin(x).count("1") for x in range(d)])
    # (popcount_ket, popcount_bra) block weights and n_diff histogram
    blk = {}
    ndiff = {}
    for a in range(d):
        for b in range(d):
            wab = wts[a, b]
            if wab < 1e-9:
                continue
            key = (pc[a], pc[b])
            blk[key] = blk.get(key, 0.0) + wab
            nd = bin(a ^ b).count("1")
            ndiff[nd] = ndiff.get(nd, 0.0) + wab
    top = sorted(blk.items(), key=lambda kv: -kv[1])[:4]
    rate_full = -lam.real
    rate_block = block_slowest(N, Q, profile)
    print(f"  {label}  Q={Q:>6.1f}: full slowest -Re={rate_full:.6f}  block(1,0) -Re={rate_block:.6f}"
          f"  gap={abs(rate_full-rate_block):.3f}")
    print(f"      (ket#,bra#) blocks: " + ", ".join(f"{k}:{w:.2f}" for k, w in top))
    print(f"      n_diff hist: " + ", ".join(f"{k}:{v:.2f}" for k, v in sorted(ndiff.items())))
    return rate_full, rate_block


def main():
    N = 6
    deep = np.array([0.25, 1.375, 1.375, 1.375, 1.375, 0.25])
    print(f"N={N} deep-edge {list(deep)}  (two protected gamma=0.25 edges)\n")
    for Q in (1.5, 1000.0):
        analyze(N, Q, deep, "deep-edge")
    print()
    # control: the N=5 flat-bulk-edge anchor (where part A said block==full)
    N5 = 5
    fbe = np.array([0.25, 1.5, 1.5, 1.5, 0.25])
    print(f"control N={N5} flat-bulk-edge {list(fbe)} (part A: block==full):\n")
    for Q in (1.5, 1000.0):
        analyze(N5, Q, fbe, "fbe-N5   ")


if __name__ == "__main__":
    main()
