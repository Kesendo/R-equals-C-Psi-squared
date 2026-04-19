#!/usr/bin/env python3
"""eq014_step23_biorth.py

EQ-014 Step 2+3: classify modes + build biorthogonal basis.

Input: simulations/results/eq014_*.bin (+ eq014_metadata.json)
Produced by C# `dotnet run -c Release -- ptf` in RCPsiSquared.Compute.

Pipeline:
  1. Load eigenvalues + right eigenvectors R + left eigenvectors L_left.
     LAPACK pairs them in the same order (column j of both corresponds
     to eigenvalue j).
  2. Diagonal pairing check: c_j = <u_j | v_j> = Σ_k conj(L[k,j]) * R[k,j].
     Verify |c_j| is not near zero (would indicate defective eigenpair).
  3. Normalize: R_j /= c_j so that <u_j | v_j> = 1.
  4. Biorthogonality residual: ||W^H V - I||_F. Target < 1e-8.
     If the eigenvalues have degenerate clusters, additional in-cluster
     biorthogonalization (SVD-based) is applied.
  5. Classify slow modes: |Re λ| <= SLOW_CUTOFF (0.15 for N=7, γ=0.05).
     For each slow mode, compute Pauli projection and dominant label.

Outputs:
  simulations/results/eq014_biorthogonal_check.txt
  simulations/results/eq014_slow_modes.json
  simulations/results/eq014_biorth_right.bin    (rescaled R, same shape as input)
  simulations/results/eq014_biorth_left.bin     (rescaled L_left, same shape as input)
  simulations/results/eq014_slow_indices.npz    (indices + eigenvalues of slow modes)

Date: 2026-04-19
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding='utf-8')

RESULTS_DIR = Path(__file__).parent / "results"
SLOW_CUTOFF = 0.15            # |Re(λ)| <= SLOW_CUTOFF -> slow mode
DEGEN_TOL = 1e-6              # cluster-degeneracy threshold (in Δ eigenvalue)
RESIDUAL_TARGET = 1e-8        # acceptable ||W^H V - I||_F / d2

# Pauli single-site operators, big-endian convention (site 0 = MSB)
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = [I2, X, Y, Z]
PAULI_NAME = "IXYZ"


def load_metadata(path):
    with open(path) as f:
        return json.load(f)


def load_bin_matrix(path, shape):
    """Load binary file as complex128 column-major (Fortran-order) matrix."""
    t0 = time.time()
    data = np.fromfile(path, dtype=np.complex128)
    assert data.size == shape[0] * shape[1], \
        f"expected {shape[0] * shape[1]} elements, got {data.size}"
    mat = data.reshape(shape, order='F')
    print(f"  Loaded {path.name} {shape} in {time.time() - t0:.1f} s "
          f"({data.nbytes / 1e9:.2f} GB)")
    return mat


def classify_mode_sector(eigvec, d, N, popcount):
    """Classify an eigenvector M_s by its dominant (n_a, n_b) U(1) sector.
    Convention: vec(rho)[a*d + b] = rho[a, b] with MSB site-indexing
    (site 0 = MSB, matching PauliOps.At and BuildDirectRaw in the C# engine).

    Returns dict with dominant sector, its weight fraction, and total norm.
    """
    M = eigvec.reshape((d, d))  # row-major: M[a, b] = rho[a, b]
    sector_weights = np.zeros((N + 1, N + 1), dtype=float)
    # Group indices by popcount and sum |M|^2 per block.
    buckets = [np.where(popcount == n)[0] for n in range(N + 1)]
    for n_a in range(N + 1):
        rows = buckets[n_a]
        if len(rows) == 0:
            continue
        for n_b in range(N + 1):
            cols = buckets[n_b]
            if len(cols) == 0:
                continue
            sub = M[np.ix_(rows, cols)]
            sector_weights[n_a, n_b] = float(np.sum(np.abs(sub) ** 2))
    total = float(np.sum(sector_weights))
    if total < 1e-30:
        return {"n_a": -1, "n_b": -1, "weight_frac": 0.0, "total_norm": 0.0,
                "sector_weights_flat": sector_weights.flatten().tolist()}
    idx = np.unravel_index(int(np.argmax(sector_weights)), sector_weights.shape)
    frac = float(sector_weights[idx] / total)
    # XY-weight = |n_a - n_b| (Hamming distance for single-particle blocks)
    xy_weight = abs(int(idx[0]) - int(idx[1]))
    return {
        "n_a": int(idx[0]),
        "n_b": int(idx[1]),
        "xy_weight": xy_weight,
        "weight_frac": frac,
        "total_norm": float(np.sqrt(total)),
    }


def classify_modes_efficient(slow_indices, R, d, N):
    """Sector classification for slow modes. Returns list of dicts."""
    popcount = np.array([bin(i).count('1') for i in range(d)])
    results = []
    for s in slow_indices:
        ev = R[:, s]
        cls = classify_mode_sector(ev, d, N, popcount)
        results.append({
            "index": int(s),
            "sector": cls,
        })
    return results


def biorth_residual_sample(W, V, n_samples=200, rng=None):
    """Sample-based biorth residual check. Avoids full d^2 x d^2 multiply.
    Computes ||W[:, idx]^H V[:, idx] - I_sub||_F on a random subset."""
    if rng is None:
        rng = np.random.default_rng(42)
    d2 = W.shape[1]
    idx = rng.choice(d2, size=min(n_samples, d2), replace=False)
    W_sub = W[:, idx]  # (d2, n_samples)
    V_sub = V[:, idx]
    prod = W_sub.conj().T @ V_sub  # (n_samples, n_samples)
    I_sub = np.eye(len(idx), dtype=complex)
    res = np.linalg.norm(prod - I_sub, ord='fro') / len(idx)
    return float(res), idx


def find_degenerate_clusters(vals, tol=DEGEN_TOL):
    """Group eigenvalue indices into clusters of nearly-degenerate values.
    Returns list of arrays of indices (each cluster)."""
    d2 = len(vals)
    # Sort by complex eigenvalue (lex: real then imag)
    order = np.argsort(vals.real + 1j * vals.imag * 1e-8)  # break ties by imag
    clusters = []
    current = [order[0]]
    for i in range(1, d2):
        prev = order[i - 1]
        curr = order[i]
        if abs(vals[curr] - vals[prev]) < tol:
            current.append(curr)
        else:
            if len(current) > 1:
                clusters.append(np.array(current))
            current = [curr]
    if len(current) > 1:
        clusters.append(np.array(current))
    return clusters


def biorthogonalize_cluster(W_sub, V_sub):
    """Re-biorthogonalize a degenerate cluster.
    Given W_sub, V_sub with shape (d2, k) where k is cluster size,
    returns new W_sub, V_sub such that (W_sub)^H V_sub = I_k.

    Uses SVD: Let M = W_sub^H V_sub (k x k). Compute M = UΣV^H.
    Then define:
      W_new = W_sub @ U @ Σ^{-1/2}   (shape (d2, k))
      V_new = V_sub @ V @ Σ^{-1/2}   (shape (d2, k))
    check:  W_new^H V_new = Σ^{-1/2} U^H M V Σ^{-1/2} = Σ^{-1/2} Σ Σ^{-1/2} = I.
    """
    M = W_sub.conj().T @ V_sub
    # SVD of M
    U, s, Vh = np.linalg.svd(M)
    sqrt_s_inv = np.diag(1.0 / np.sqrt(s))
    W_new = W_sub @ (U @ sqrt_s_inv)
    V_new = V_sub @ (Vh.conj().T @ sqrt_s_inv)
    return W_new, V_new


def main():
    t_total = time.time()
    meta = load_metadata(RESULTS_DIR / "eq014_metadata.json")
    N = meta["N"]
    d = meta["d"]
    d2 = meta["d2"]
    gamma_0 = meta["gamma"]
    print(f"=== EQ-014 Step 2+3: biorthogonal basis ===")
    print(f"N={N}, d={d}, d2={d2}, γ={gamma_0}, topology={meta['topology']}")
    print()

    # Load eigenvalues
    vals = np.fromfile(RESULTS_DIR / "eq014_eigvals_n7.bin", dtype=np.complex128)
    assert len(vals) == d2
    print(f"Eigenvalues: {len(vals)}")
    print(f"  Re range: [{vals.real.min():.4f}, {vals.real.max():.4f}]")
    print(f"  |Im| max: {np.abs(vals.imag).max():.4f}")
    stationary = int(np.sum(np.abs(vals) < 1e-10))
    print(f"  Stationary: {stationary} (expected {N + 1})")
    print()

    # Load eigenvectors (column-major). For N=7, each is 4 GB.
    print("Loading right eigenvectors...")
    R = load_bin_matrix(RESULTS_DIR / "eq014_right_eigvecs_n7.bin", (d2, d2))
    print("Loading left eigenvectors...")
    W = load_bin_matrix(RESULTS_DIR / "eq014_left_eigvecs_n7.bin", (d2, d2))
    print()

    # Diagonal pairing check: c_j = <u_j | v_j>
    print("Computing diagonal c_j = <u_j | v_j>...")
    t0 = time.time()
    c_diag = np.einsum('kj,kj->j', W.conj(), R)
    print(f"  Done in {time.time() - t0:.1f} s")
    print(f"  |c_j| range: [{np.abs(c_diag).min():.2e}, {np.abs(c_diag).max():.2e}]")
    small_mask = np.abs(c_diag) < 1e-6
    print(f"  Modes with |c_j| < 1e-6: {int(np.sum(small_mask))}")
    if np.any(small_mask):
        small_idx = np.where(small_mask)[0][:5]
        for j in small_idx:
            print(f"    j={j}, λ={vals[j]:.4e}, c_j={c_diag[j]:.3e}")
    print()

    # Normalize: R[:, j] /= c_diag[j]  so that <u_j | v_j> = 1
    # To avoid division by zero for degenerate pairs, use safe normalization
    safe_c = np.where(np.abs(c_diag) < 1e-14, 1.0, c_diag)
    print("Normalizing right eigenvectors (R /= c_diag)...")
    t0 = time.time()
    R /= safe_c[np.newaxis, :]
    print(f"  Done in {time.time() - t0:.1f} s")
    print()

    # Sample-based biorth residual
    print("Sample-based biorthogonality residual (200 random columns)...")
    t0 = time.time()
    res_sample, _ = biorth_residual_sample(W, R, n_samples=200)
    print(f"  ||W_s^H R_s - I||_F / n = {res_sample:.3e} (done in {time.time() - t0:.1f} s)")
    print()

    # Degenerate cluster detection and in-cluster biorthogonalization
    print(f"Finding degenerate clusters (tol={DEGEN_TOL})...")
    clusters = find_degenerate_clusters(vals, tol=DEGEN_TOL)
    total_in_clusters = sum(len(c) for c in clusters)
    print(f"  Found {len(clusters)} clusters, {total_in_clusters} modes total")
    for k, cl in enumerate(clusters[:10]):
        print(f"    cluster {k}: size={len(cl)}, λ≈{vals[cl[0]]:.4e}")
    if len(clusters) > 10:
        print(f"    ... and {len(clusters) - 10} more")
    print()

    print("In-cluster SVD biorthogonalization...")
    t0 = time.time()
    clusters_fixed = 0
    for cl in clusters:
        if len(cl) < 2:
            continue
        W_sub = W[:, cl].copy()
        R_sub = R[:, cl].copy()
        W_new, R_new = biorthogonalize_cluster(W_sub, R_sub)
        W[:, cl] = W_new
        R[:, cl] = R_new
        clusters_fixed += 1
    print(f"  Fixed {clusters_fixed} clusters in {time.time() - t0:.1f} s")
    print()

    # Re-check biorthogonality residual after fix
    print("Post-fix sample-based biorth residual (200 random)...")
    t0 = time.time()
    res_post, _ = biorth_residual_sample(W, R, n_samples=200)
    print(f"  ||W_s^H R_s - I||_F / n = {res_post:.3e} (done in {time.time() - t0:.1f} s)")
    print()

    # Classify slow modes
    print(f"Identifying slow modes (|Re λ| <= {SLOW_CUTOFF})...")
    slow_mask = np.abs(vals.real) <= SLOW_CUTOFF
    slow_indices = np.where(slow_mask)[0]
    print(f"  Found {len(slow_indices)} slow modes")
    print()

    print("Classifying slow modes by Pauli decomposition...")
    t0 = time.time()
    classifications = classify_modes_efficient(slow_indices, R, d, N)
    print(f"  Done in {time.time() - t0:.1f} s")
    print()

    # Save outputs
    print("Saving outputs...")

    # Slow indices + eigenvalues
    np.savez(
        RESULTS_DIR / "eq014_slow_indices.npz",
        indices=slow_indices,
        eigenvalues=vals[slow_indices],
        all_eigenvalues=vals,
    )

    # Slow mode classification JSON
    slow_modes_json = {
        "N": N,
        "d": d,
        "d2": d2,
        "gamma": gamma_0,
        "slow_cutoff": SLOW_CUTOFF,
        "count": len(slow_indices),
        "modes": [
            {
                "rank": i,
                "index": int(s),
                "eigenvalue_re": float(vals[s].real),
                "eigenvalue_im": float(vals[s].imag),
                "pairing_c_abs": float(abs(c_diag[s])),
                "pairing_c_phase": float(np.angle(c_diag[s])),
                "sector": cls["sector"],
            }
            for i, (s, cls) in enumerate(zip(slow_indices, classifications))
        ],
    }
    with open(RESULTS_DIR / "eq014_slow_modes.json", "w", encoding="utf-8") as f:
        json.dump(slow_modes_json, f, indent=2)
    print(f"  Wrote eq014_slow_modes.json ({len(slow_indices)} modes)")

    # Biorthogonal check log
    with open(RESULTS_DIR / "eq014_biorthogonal_check.txt", "w", encoding="utf-8") as f:
        f.write(f"EQ-014 Biorthogonal Check (Step 2+3)\n")
        f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"N={N}, d2={d2}, γ={gamma_0}\n\n")
        f.write(f"Eigenvalue count: {len(vals)}\n")
        f.write(f"Stationary (|λ|<1e-10): {stationary} (expected {N + 1})\n")
        f.write(f"Diag pairing |c_j| range: [{np.abs(c_diag).min():.2e}, "
                f"{np.abs(c_diag).max():.2e}]\n")
        f.write(f"Modes with |c_j| < 1e-6: {int(np.sum(small_mask))}\n\n")
        f.write(f"Degenerate clusters (tol={DEGEN_TOL}):\n")
        f.write(f"  Cluster count: {len(clusters)}\n")
        f.write(f"  Modes in clusters: {total_in_clusters}\n")
        f.write(f"  Cluster sizes: {[len(c) for c in clusters[:20]]}\n\n")
        f.write(f"Biorthogonality residual (sampled 200 random columns):\n")
        f.write(f"  Before SVD fix: {res_sample:.3e}\n")
        f.write(f"  After SVD fix:  {res_post:.3e}\n")
        f.write(f"  Target:         < {RESIDUAL_TARGET:.0e}\n\n")
        verdict = "PASS" if res_post < RESIDUAL_TARGET else "INVESTIGATE"
        f.write(f"Verdict: {verdict}\n")
        f.write(f"Slow modes found: {len(slow_indices)}\n")
    print(f"  Wrote eq014_biorthogonal_check.txt")

    # Save cluster indices so Step 4-7 can reproduce the in-cluster fix
    # without storing two 4 GB biorth binaries on disk.
    # Flat concat + offsets (CSR-style).
    cluster_flat = np.concatenate([cl for cl in clusters]) if clusters else np.array([], dtype=int)
    cluster_offsets = np.array([0] + list(np.cumsum([len(c) for c in clusters])))
    np.savez(
        RESULTS_DIR / "eq014_biorth_fix.npz",
        c_diag=c_diag,
        cluster_flat=cluster_flat,
        cluster_offsets=cluster_offsets,
        slow_indices=slow_indices,
    )
    print(f"  Wrote eq014_biorth_fix.npz (c_diag + {len(clusters)} clusters)")

    # Copy metadata + augment
    meta["biorth_residual_pre"] = res_sample
    meta["biorth_residual_post"] = res_post
    meta["degenerate_cluster_count"] = len(clusters)
    meta["slow_mode_count"] = int(len(slow_indices))
    meta["slow_cutoff"] = SLOW_CUTOFF
    with open(RESULTS_DIR / "eq014_biorth_metadata.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    print()
    print(f"=== Step 2+3 complete in {time.time() - t_total:.1f} s ===")
    print(f"  Biorthogonality residual: {res_post:.3e} ({verdict})")
    print(f"  Slow modes: {len(slow_indices)}")


if __name__ == "__main__":
    main()
