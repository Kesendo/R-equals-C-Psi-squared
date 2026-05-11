"""F89 path-6 (single 7-qubit chain at N=7, no bare sites): numerical multi-exp.

Topology (6) at N=7: bonds {0,1,2,3,4,5} all active = 7 connected sites = full
chain, zero bare sites. Single block, m=1.

L_super dim = 4^7 = 16384. np.linalg.eig at this dim: ~12 GB RAM, ~30-60 min.

This completes the path-k inventory (k=1..6) for the F89 closed-form program.
Trivially satisfies the additive identity since m=1 → no subtraction term:
S_(6)(t) = block contribution = full result (no bare component).

Verification target: bond-isolate `N7_b0-1-2-3-4-5_J0.0750_gamma0.0500_probe-coherence.csv`
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO = Path(__file__).resolve().parents[1]
CSV_DIR = REPO / "simulations" / "results" / "bond_isolate"

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(P: np.ndarray, site: int, n: int) -> np.ndarray:
    op = np.array([[1.0]], dtype=complex)
    for q in range(n):
        op = np.kron(op, P if q == site else I2)
    return op


def main() -> None:
    J, gamma = 0.075, 0.05
    N = 7
    n_block = 7
    D = 2**n_block

    print(f"# F89 path-6 (full chain at N={N}): J={J}, γ={gamma}\n")
    print(f"# L_super dim: {D*D} ({D}×{D} block)")

    # Build path-6 H_B
    t0 = time.time()
    print("# Building H_B...", flush=True)
    H = np.zeros((D, D), dtype=complex)
    for b in range(n_block - 1):
        H += J * (
            kron_at(X, b, n_block) @ kron_at(X, b + 1, n_block)
            + kron_at(Y, b, n_block) @ kron_at(Y, b + 1, n_block)
        )
    print(f"# H_B built in {time.time() - t0:.1f}s", flush=True)

    t1 = time.time()
    print("# Building L_super (16384×16384)...", flush=True)
    Id = np.eye(D, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(n_block):
        Zl = kron_at(Z, l, n_block)
        L += gamma * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    print(f"# L_super built in {time.time() - t1:.1f}s; size {L.nbytes / 1024**3:.2f} GB", flush=True)

    # Build ρ_block(0) for topology (6) at N=7 (no bare sites: N_E=0)
    bit_pos = [2 ** (n_block - 1 - i) for i in range(n_block)]

    def state_idx(bits):
        return sum(bit_pos[i] * bits[i] for i in range(n_block))

    N_E = 0  # no bare sites
    pre = 1.0 / np.sqrt(N * N * (N - 1) / 2)
    rho = np.zeros((D, D), dtype=complex)
    print("# Computing ρ_block(0)...", flush=True)
    for i in range(n_block):
        bits = [0] * n_block
        bits[i] = 1
        idx_se = state_idx(bits)
        for j in range(n_block):
            for k in range(j + 1, n_block):
                bits_de = [0] * n_block
                bits_de[j] = 1
                bits_de[k] = 1
                rho[idx_se, state_idx(bits_de)] += pre
    # N_E = 0, so no term-2 contribution
    rho = (rho + rho.conj().T) / 2.0

    vec = rho.flatten(order="F")

    t2 = time.time()
    print("# Eigendecomposing L_super (this is the slow part, ~30-60 min)...", flush=True)
    eigvals, R = np.linalg.eig(L)
    print(f"# Eigendecomposition done in {time.time() - t2:.1f}s", flush=True)

    t3 = time.time()
    print("# Solving R c = vec(ρ_0)...", flush=True)
    c = np.linalg.solve(R, vec)
    print(f"# Solve done in {time.time() - t3:.1f}s", flush=True)

    # Per-site reduction setup
    def reduce_to_site(rho_block, l):
        other = [s for s in range(n_block) if s != l]
        val = 0.0 + 0.0j
        for cc in range(2 ** (n_block - 1)):
            bits_other = [(cc >> (n_block - 2 - i)) & 1 for i in range(n_block - 1)]
            idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(n_block - 1))
            idx_1 = idx_0 + bit_pos[l]
            val += rho_block[idx_0, idx_1]
        return val

    # Verification
    csv_path = CSV_DIR / f"N{N}_b0-1-2-3-4-5_J0.0750_gamma0.0500_probe-coherence.csv"
    if not csv_path.exists():
        print(f"\n# CSV not found: {csv_path}")
        return

    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    t_csv, S_csv = data[:, 0], data[:, -1]

    t4 = time.time()
    print("# Computing S(t) at all CSV time points...", flush=True)
    S_pred = np.zeros_like(t_csv)
    for ti, t in enumerate(t_csv):
        vec_t = R @ (np.exp(eigvals * t) * c)
        rho_t = vec_t.reshape((D, D), order="F")
        S_pred[ti] = sum(2.0 * abs(reduce_to_site(rho_t, l)) ** 2 for l in range(n_block))
    print(f"# S(t) computed in {time.time() - t4:.1f}s", flush=True)

    diff = S_pred - S_csv
    print()
    print(f"## N=7 path-6 verification vs bond-isolate CSV")
    print(f"# CSV: {csv_path.name}")
    print(f"# max |diff|: {np.max(np.abs(diff)):.3e} (CSV write precision is ~5e-7)")
    print(f"# mean |diff|: {np.mean(np.abs(diff)):.3e}")
    print(f"# S(0): pred={S_pred[0]:.6f}, csv={S_csv[0]:.6f}, expect (N-1)/N={(N-1)/N:.6f}")
    print()
    print("| t  | S_csv (RK4) | S_pred (closed form) | diff |")
    print("|---|---|---|---|")
    for i in [0, 30, 50, 100, 150, 200, 250, 300]:
        if i < len(t_csv):
            print(
                f"| {t_csv[i]:5.2f} | {S_csv[i]:.6f} | {S_pred[i]:.6f} | "
                f"{diff[i]:+.3e} |"
            )

    # Mode count
    sig = np.zeros(D * D)
    # Per-site reduction matrix w[l, idx]
    w = np.zeros((n_block, D * D), dtype=complex)
    for l in range(n_block):
        other = [s for s in range(n_block) if s != l]
        for cc in range(2 ** (n_block - 1)):
            bits_other = [(cc >> (n_block - 2 - i)) & 1 for i in range(n_block - 1)]
            idx_0 = sum(bit_pos[other[i]] * bits_other[i] for i in range(n_block - 1))
            idx_1 = idx_0 + bit_pos[l]
            w[l, idx_1 * D + idx_0] = 1.0
    M = w @ R
    a = M * c[None, :]
    sig = np.sum(np.abs(a) ** 2, axis=0)
    contributing = np.where(sig > 1e-12)[0]
    rates = -eigvals.real / gamma
    freqs = np.abs(eigvals.imag) / J
    groups = {}
    for k in contributing:
        key = (round(rates[k], 4), round(freqs[k], 4))
        groups[key] = groups.get(key, 0.0) + sig[k]

    print(f"\n## Path-6 mode-group survey")
    print(f"# 16384 total L_super modes, {len(contributing)} populated, {len(groups)} distinct (rate, |freq|) groups")
    print(f"# (Compare path-2..5: 4, 10, 12, 35 mode-groups respectively)")


if __name__ == "__main__":
    main()
