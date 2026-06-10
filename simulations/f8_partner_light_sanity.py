# Sanity check for F8PartnerLightComplementarityTests (C# Diagnostics.Tests/Ptf).
# Verifies index conventions for: (a) absorption identity Re λ = −2γ·light(v),
# (b) F8 pairing λ_s + λ_f = −2Nγ with light complementarity light_s + light_f = N,
# (c) the April standing-wave null max |⟨W_f|V_L|M_s⟩| over partner pairs.
# Conventions match the C# layer: row-major vec (x = a·d + b), site 0 = most
# significant bit, H = Σ_bonds (J/2)(XX+YY), dissipator diag −2γ·popcount(a⊕b).
import sys

import numpy as np

if sys.stdout.encoding and sys.stdout.encoding.lower() != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def site_pair(N, i, j, P, Q):
    ops = [I2] * N
    ops[i] = P
    ops[j] = Q
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def run(N, gamma=0.05, J=1.0, tol_pair=1e-8):
    d = 1 << N
    H = sum(J / 2.0 * (site_pair(N, b, b + 1, X, X) + site_pair(N, b, b + 1, Y, Y))
            for b in range(N - 1))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    w = np.array([bin((x // d) ^ (x % d)).count("1") for x in range(d * d)], dtype=float)
    L += np.diag(-2.0 * gamma * w)

    vals, R = np.linalg.eig(L)
    W = np.linalg.inv(R)  # rows globally biorthogonal to R's columns, ⟨W_k|M_k⟩ = 1

    # (a) absorption identity
    prob = np.abs(R) ** 2
    light = (w @ prob) / prob.sum(axis=0)
    abs_resid = np.max(np.abs(vals.real + 2.0 * gamma * light))

    # (b) greedy multiset pairing λ_i + λ_j ≈ −2Nγ
    target = -2.0 * N * gamma
    unpaired = list(range(d * d))
    pairs = []
    worst_pair_resid = 0.0
    while unpaired:
        i = unpaired.pop(0)
        resid = [abs(vals[i] + vals[j] - target) for j in unpaired]
        k = int(np.argmin(resid))
        worst_pair_resid = max(worst_pair_resid, resid[k])
        j = unpaired.pop(k)
        pairs.append((i, j))
    complete = worst_pair_resid < tol_pair
    comp_resid = max(abs(light[i] + light[j] - N) for i, j in pairs)

    # (c) standing-wave null: V = ½(X0X1 + Y0Y1), V_L = −i[V,·], cross-partner elements
    V = 0.5 * (site_pair(N, 0, 1, X, X) + site_pair(N, 0, 1, Y, Y))
    VL = -1j * (np.kron(V, Id) - np.kron(Id, V.T))
    K = W @ VL @ R
    null_max = max(max(abs(K[j, i]), abs(K[i, j])) for i, j in pairs)

    print(f"N={N}: dim={d * d}, pairs={len(pairs)}")
    print(f"  (a) max |Re λ + 2γ·light|        = {abs_resid:.3e}")
    print(f"  (b) pairing complete: {complete} (worst |λ_s+λ_f+2Nγ| = {worst_pair_resid:.3e})")
    print(f"      max |light_s + light_f − N|  = {comp_resid:.3e}")
    print(f"  (c) max cross-partner |⟨W_f|V_L|M_s⟩| = {null_max:.3e}")
    print(f"      biorth check ‖W·R − I‖_max   = {np.max(np.abs(W @ R - np.eye(d * d))):.3e}")


if __name__ == "__main__":
    run(3)
    run(4)
