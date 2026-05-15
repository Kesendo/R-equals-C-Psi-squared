"""Step 6 attempt with JW infrastructure: V_inter R-parity split on c=2 stratum.

Restricts V_inter to the c=2 coherence sub-block (popcount-1 ↔ popcount-2,
which is exactly the stratum covered by RCPsiSquared.Core/F86/JordanWigner/
JwBlockBasis). Checks whether the parity-of-N effect (σ_0⁺ = σ_0⁻ at even N,
≠ at odd N) holds on this restricted sub-block too.

If yes, the JwClusterDEigenstructure 'same-size W_c matrices are unitarily
equivalent via F71-mirror invariance' (Tier1Derived) argument structurally
matches the V⁺⁺ ≅ V⁻⁻ structure on c=2; Step 6's open formal step reduces
to extending the cluster-equivalence argument from W_c (D-matrix clusters)
to V_inter (HD-block cross-section).

Run: python simulations/_step6_c2_v_inter.py
"""
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

PAULI = [
    np.array([[1, 0], [0, 1]], dtype=complex),
    np.array([[0, 1], [1, 0]], dtype=complex),
    np.array([[0, -1j], [1j, 0]], dtype=complex),
    np.array([[1, 0], [0, -1]], dtype=complex),
]


def pauli_string(alpha):
    result = PAULI[alpha[0]]
    for a in alpha[1:]:
        result = np.kron(result, PAULI[a])
    return result


def reverse_bits(b, n):
    out = 0
    for i in range(n):
        if (b >> i) & 1:
            out |= 1 << (n - 1 - i)
    return out


def analyse(N, J_ref=1.0):
    dim = 2 ** N
    pc1_states = sorted([s for s in range(dim) if bin(s).count('1') == 1])
    pc2_states = sorted([s for s in range(dim) if bin(s).count('1') == 2])

    hd1_c2_ops = []
    for p1 in pc1_states:
        for p2 in pc2_states:
            if bin(p1 ^ p2).count('1') == 1:
                hd1_c2_ops.append((p1, p2))
    hd3_c2_ops = []
    for p1 in pc1_states:
        for p2 in pc2_states:
            if bin(p1 ^ p2).count('1') == 3:
                hd3_c2_ops.append((p1, p2))

    H = np.zeros((dim, dim), dtype=complex)
    for l in range(N - 1):
        for axis in (1, 2):
            op = [0] * N
            op[l] = axis
            op[l + 1] = axis
            H += (J_ref / 2) * pauli_string(op)

    V_inter = np.zeros((len(hd1_c2_ops), len(hd3_c2_ops)), dtype=complex)
    for i, (a, b) in enumerate(hd1_c2_ops):
        for j, (p1, p2) in enumerate(hd3_c2_ops):
            val = complex(0)
            if b == p2:
                val += H[a, p1]
            if a == p1:
                val -= H[p2, b]
            V_inter[i, j] = -1j * val

    hd1_idx = {op: i for i, op in enumerate(hd1_c2_ops)}
    hd3_idx = {op: i for i, op in enumerate(hd3_c2_ops)}

    R_hd1 = np.zeros((len(hd1_c2_ops), len(hd1_c2_ops)), dtype=float)
    for i, (p1, p2) in enumerate(hd1_c2_ops):
        rp1, rp2 = reverse_bits(p1, N), reverse_bits(p2, N)
        j = hd1_idx[(rp1, rp2)]
        R_hd1[j, i] = 1.0

    R_hd3 = np.zeros((len(hd3_c2_ops), len(hd3_c2_ops)), dtype=float)
    for i, (p1, p2) in enumerate(hd3_c2_ops):
        rp1, rp2 = reverse_bits(p1, N), reverse_bits(p2, N)
        j = hd3_idx[(rp1, rp2)]
        R_hd3[j, i] = 1.0

    r_hd1_vals, r_hd1_vecs = np.linalg.eigh(R_hd1)
    r_hd3_vals, r_hd3_vecs = np.linalg.eigh(R_hd3)
    even_hd1 = r_hd1_vecs[:, np.abs(r_hd1_vals - 1) < 1e-6]
    odd_hd1 = r_hd1_vecs[:, np.abs(r_hd1_vals + 1) < 1e-6]
    even_hd3 = r_hd3_vecs[:, np.abs(r_hd3_vals - 1) < 1e-6]
    odd_hd3 = r_hd3_vecs[:, np.abs(r_hd3_vals + 1) < 1e-6]

    V_pp = even_hd1.T @ V_inter @ even_hd3
    V_mm = odd_hd1.T @ V_inter @ odd_hd3
    V_pm = even_hd1.T @ V_inter @ odd_hd3
    V_mp = odd_hd1.T @ V_inter @ even_hd3

    s_pp = np.linalg.svd(V_pp, compute_uv=False)
    s_mm = np.linalg.svd(V_mm, compute_uv=False)

    return {
        "N": N,
        "hd1_c2_dim": len(hd1_c2_ops),
        "hd3_c2_dim": len(hd3_c2_ops),
        "even_hd1": even_hd1.shape[1],
        "odd_hd1": odd_hd1.shape[1],
        "even_hd3": even_hd3.shape[1],
        "odd_hd3": odd_hd3.shape[1],
        "V_pp_norm": np.linalg.norm(V_pp),
        "V_mm_norm": np.linalg.norm(V_mm),
        "V_cross_norm": np.linalg.norm(V_pm) + np.linalg.norm(V_mp),
        "s_pp": s_pp,
        "s_mm": s_mm,
    }


for N in [3, 4, 5, 6]:
    r = analyse(N)
    print(f"\n# N = {N}")
    print(f"#   c=2 HD=1 dim: {r['hd1_c2_dim']}, HD=3 dim: {r['hd3_c2_dim']}")
    print(f"#   R-even HD=1: {r['even_hd1']}, R-odd HD=1: {r['odd_hd1']}")
    print(f"#   R-even HD=3: {r['even_hd3']}, R-odd HD=3: {r['odd_hd3']}")
    print(f"#   ‖V_++‖_F = {r['V_pp_norm']:.4f}, ‖V_--‖_F = {r['V_mm_norm']:.4f}")
    print(f"#   ‖V_cross‖_F = {r['V_cross_norm']:.2e} (should ~ 0)")
    n_show = min(6, len(r['s_pp']))
    print(f"#   Top σ V_++: {[f'{s:.4f}' for s in r['s_pp'][:n_show]]}")
    n_show_m = min(6, len(r['s_mm']))
    print(f"#   Top σ V_--: {[f'{s:.4f}' for s in r['s_mm'][:n_show_m]]}")
    if r['s_pp'][0] > 1e-9 and r['s_mm'][0] > 1e-9:
        print(f"#   |σ_0⁺ − σ_0⁻| = {abs(r['s_pp'][0] - r['s_mm'][0]):.6e}")
