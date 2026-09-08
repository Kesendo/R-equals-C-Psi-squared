#!/usr/bin/env python3
"""Read the N=5 slow tolerance-subspace light distribution at canonical Q=3 and Q=2000.

The slow real-rate tolerance cluster can contain several modes and need not be exactly degenerate.
This producer therefore uses the orthogonal projector onto the selected right-invariant subspace,
not an arbitrary eigenvector. The Absorption Theorem reproduces the CLUSTER-MEAN decay rate,
2*sum_k gamma_k*light_k; the spectral-edge rate is reported separately. Distribution drift is a
separate diagnostic: zero spectral-edge drift does not by itself prove a frozen distribution.
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


def slow_subspace(N, Q, profile, H1, cluster_tolerance=1e-6):
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (Q / 2.0) * (np.kron(Id, H1) - np.kron(H1.T, Id))
    for l in range(N):
        Zl = op_at(N, l, Z)
        L += profile[l] * (np.kron(Zl, Zl) - np.kron(Id, Id))
    w, V = np.linalg.eig(L)
    nonkernel = np.abs(w) > 1e-7
    edge = float(np.max(np.where(nonkernel, w.real, -np.inf)))
    cluster = np.where(nonkernel & (np.abs(w.real - edge) <= cluster_tolerance))[0]
    selected = V[:, cluster]
    U, singular, _ = np.linalg.svd(selected, full_matrices=False)
    rank_tol = max(selected.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.sum(singular > rank_tol))
    basis = U[:, :rank]
    return -edge, -float(np.mean(w[cluster].real)), basis, len(cluster)


def projector_light(basis, N):
    """Basis-free Tr(Pi_V Delta_k)/dim(V) for the selected right-invariant subspace V."""
    d = 2 ** N
    xy_site = np.zeros(N)
    indices = np.arange(d * d)
    row = indices % d                    # column-major vec: index = row + d*column
    col = indices // d
    weights = np.sum(np.abs(basis) ** 2, axis=1) / basis.shape[1]
    for site in range(N):
        bit = N - 1 - site
        differs = ((row >> bit) & 1) != ((col >> bit) & 1)
        xy_site[site] = float(np.sum(weights[differs]))
    return xy_site, float(np.sum(xy_site))


def main():
    N = 5
    H1 = H_xy_unit(N)
    profiles = {
        "peaked-V  (sterile)":   [0.25, 0.75, 3.0, 0.75, 0.25],
        "flat-bulk (birth canal)": [0.25, 1.5, 1.5, 1.5, 0.25],
    }
    print(f"N={N}. Basis-free slow tolerance-subspace light at canonical Q=3 and Q=2000.\n")
    for name, p in profiles.items():
        p = list(np.array(p, float) * N / np.sum(p))
        print(f"  {name}   profile {np.round(p,3).tolist()}")
        xy_lo = xy_hi = None
        for Q in [3.0, 2000.0]:
            edge_rate, mean_rate, basis, cluster_size = slow_subspace(N, Q, p, H1)
            xy_site, nxy = projector_light(basis, N)
            pred = 2.0 * sum(p[k] * xy_site[k] for k in range(N))
            print(f"    Q={Q:>7.1f}  edge={edge_rate:.5f}  cluster-mean={mean_rate:.5f}  "
                  f"2*Sum(g*XY)={pred:.5f}  (mean theorem err {abs(mean_rate-pred):.1e}, g={cluster_size})"
                  f"   <n_XY>={nxy:.5f}   per-site XY={np.round(xy_site,4).tolist()}")
            if Q == 3.0:
                xy_lo = xy_site
            else:
                xy_hi = xy_site
        drift = float(np.max(np.abs(xy_hi - xy_lo)))
        print(f"    max per-site distribution drift (Q: 3 -> 2000) = {drift:.5f}   "
              f"{'FROZEN DISTRIBUTION' if drift < 1e-4 else 'DISTRIBUTION DRIFTS'}\n")
    print("  reading: cluster-mean absorption and spectral-edge drift are distinct observables;")
    print("  the per-site projector distribution is required to diagnose distribution freeze.")


if __name__ == "__main__":
    main()
