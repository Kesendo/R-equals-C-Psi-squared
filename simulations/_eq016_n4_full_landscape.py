#!/usr/bin/env python3
"""EQ-016: binary-Dicke saddles above 1/4 across N.

Reframing of the original EQ-016 question. The naive reading was "N=3 is
structurally privileged for sector mixing" because F69 (GHZ+W binary at
N=3) is the only known SADDLE above 1/4. The framework doc says:

  "Pair-CPsi has no non-product local maxima on the permutation-symmetric
   Dicke sphere at any tested N."

But "no local maxima" leaves room for SADDLES. Today's finding: many
binary-Dicke pairs lift above 1/4 at every tested N, both at and beyond
N=3. F69's GHZ+W = (D_0+D_N)/sqrt(2) + D_1 is just ONE such saddle, not
the only one.

Strategy: enumerate all C(N+1, 2) Dicke pairs |D_i> + |D_j>, find their
binary-mix maxima. Report which lift above 1/4 and at what alpha. The
question "is N=3 privileged" reduces to: does the SET of binary-Dicke
saddles above 1/4 shrink, grow, or change character with N?

Open question (still requires non-symmetric scan or analytical work):
do any of these saddles correspond to LOCAL maxima of pair-CPsi within
some restricted sub-manifold (e.g., entanglement-class-fixed slice)?
Today's investigation maps the binary-Dicke landscape; full Hilbert
space saddle search is a separate harder problem.

Test A (clean): for each N in {3, 4, 5, 6}, enumerate all binary-Dicke
pair maxima. Report which are above 1/4.

Output: print + JSON report.
"""
from __future__ import annotations

import json
import math
import sys
import time
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "results" / "eq016_n4_full_landscape"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

FOLD = 0.25
PURITY_PRODUCT_THRESHOLD = 1.0 - 1e-3   # purity > this -> product flagged


def dicke_state_vector(N, k):
    v = np.zeros(2 ** N, dtype=complex)
    if k < 0 or k > N:
        return v
    norm = 1.0 / math.sqrt(comb(N, k))
    for idx in range(2 ** N):
        if bin(idx).count("1") == k:
            v[idx] = norm
    return v


def reduced_rho_AB(psi, N):
    M = psi.reshape(2 ** (N - 2), 4).T
    return M @ M.conj().T


def reduced_rho_A(rho_AB):
    t = rho_AB.reshape(2, 2, 2, 2)
    return np.trace(t, axis1=1, axis2=3)


def pair_cpsi(rho_AB):
    C = float(np.real(np.trace(rho_AB @ rho_AB)))
    diag = np.diag(np.diag(rho_AB))
    L1 = float(np.sum(np.abs(rho_AB - diag)))
    return C * L1 / 3.0


def single_qubit_purity(rho_AB):
    rho_A = reduced_rho_A(rho_AB)
    return float(np.real(np.trace(rho_A @ rho_A)))


def cpsi_from_psi(psi, N):
    rho_AB = reduced_rho_AB(psi, N)
    return pair_cpsi(rho_AB), single_qubit_purity(rho_AB), rho_AB


# -----------------------------------------------------------------------
# Test A: all binary-Dicke pair maxima at N=4
# -----------------------------------------------------------------------
def binary_dicke_optimum(N, i, j, n_grid=401):
    """For |psi(alpha)> = alpha |D_i> + sqrt(1-alpha^2) |D_j>, find max
    pair-CPsi over alpha in [0, 1]."""
    Di = dicke_state_vector(N, i)
    Dj = dicke_state_vector(N, j)

    def cpsi_alpha(alpha):
        psi = alpha * Di + math.sqrt(max(0, 1 - alpha ** 2)) * Dj
        return pair_cpsi(reduced_rho_AB(psi, N))

    # Grid search
    alphas = np.linspace(0, 1, n_grid)
    cpsis = np.array([cpsi_alpha(a) for a in alphas])
    best_idx = int(np.argmax(cpsis))
    best_alpha = alphas[best_idx]
    # Refine via bounded scalar minimization on -cpsi
    lb = max(0.0, best_alpha - 2.0 / (n_grid - 1))
    ub = min(1.0, best_alpha + 2.0 / (n_grid - 1))
    res = minimize_scalar(lambda a: -cpsi_alpha(a),
                            bounds=(lb, ub),
                            method="bounded",
                            options={"xatol": 1e-12})
    if res.success and lb <= res.x <= ub:
        best_alpha = res.x
    best_cpsi = cpsi_alpha(best_alpha)
    # Diagnostics on optimum
    psi_opt = best_alpha * Di + math.sqrt(max(0, 1 - best_alpha ** 2)) * Dj
    rho_AB = reduced_rho_AB(psi_opt, N)
    return {
        "i": i, "j": j,
        "alpha_opt": float(best_alpha),
        "cpsi_opt": float(best_cpsi),
        "purity_A_opt": float(single_qubit_purity(rho_AB)),
    }


def test_A_binary_dicke(N):
    print("=" * 80)
    print(f"Test A: binary-Dicke pair maxima at N={N}")
    print("=" * 80)
    print(f"  All C({N+1}, 2) = {comb(N+1, 2)} unordered Dicke pairs.")
    print(f"  Family: |psi(alpha)> = alpha|D_i> + sqrt(1-alpha^2)|D_j>")
    print()

    results = []
    for i, j in combinations(range(N + 1), 2):
        r = binary_dicke_optimum(N, i, j)
        results.append(r)

    results.sort(key=lambda r: -r["cpsi_opt"])
    print(f"  {'pair':<14} {'alpha_opt':>10} {'cpsi_opt':>12} "
          f"{'purity(A)':>11} {'above 1/4?':>11}")
    print("  " + "-" * 65)
    for r in results:
        flag = "YES" if r["cpsi_opt"] > FOLD else "no"
        # Mark known: GHZ = D(N,0)+D(N,N) at alpha=1/sqrt(2)
        label = ""
        if (r["i"] == 0 and r["j"] == N):
            label = " (GHZ family)"
        print(f"  D{r['i']}+D{r['j']:<10} {r['alpha_opt']:>10.4f} "
              f"{r['cpsi_opt']:>12.6f} {r['purity_A_opt']:>11.6f} "
              f"{flag:>11}{label}")
    print()

    above_fold = [r for r in results if r["cpsi_opt"] > FOLD]
    if above_fold:
        print(f"  Found {len(above_fold)} binary-Dicke pairs lifting above 1/4.")
    else:
        print(f"  NO binary-Dicke pair lifts above 1/4 at N={N}.")
    print()
    return results


# -----------------------------------------------------------------------
# Test C: 3-Dicke family maxima (includes F69's GHZ+W as D_0+D_N+D_1 at α²=0.375)
# -----------------------------------------------------------------------
def triple_dicke_optimum(N, i, j, k, n_grid=51):
    """Max pair-CPsi on |psi(c)> = c_i|D_i> + c_j|D_j> + c_k|D_k>,
    c real with c_i² + c_j² + c_k² = 1.

    Parameterise by spherical coords (theta, phi):
      c_i = cos(theta), c_j = sin(theta) cos(phi), c_k = sin(theta) sin(phi)
    """
    Di = dicke_state_vector(N, i)
    Dj = dicke_state_vector(N, j)
    Dk = dicke_state_vector(N, k)

    def cpsi_tp(theta_phi):
        t, p = theta_phi
        ci = math.cos(t)
        cj = math.sin(t) * math.cos(p)
        ck = math.sin(t) * math.sin(p)
        psi = ci * Di + cj * Dj + ck * Dk
        return pair_cpsi(reduced_rho_AB(psi, N))

    # Grid search
    thetas = np.linspace(0, math.pi, n_grid)
    phis = np.linspace(0, 2 * math.pi, n_grid)
    best = (None, -1.0)
    for t in thetas:
        for p in phis:
            c = cpsi_tp((t, p))
            if c > best[1]:
                best = ((t, p), c)
    # Refine via Nelder-Mead
    res = minimize(lambda x: -cpsi_tp(x), best[0], method="Nelder-Mead",
                    options={"xatol": 1e-10, "fatol": 1e-12})
    if res.success and -res.fun > best[1]:
        best = (tuple(res.x), -res.fun)
    t_opt, p_opt = best[0]
    cpsi_opt = best[1]
    ci = math.cos(t_opt)
    cj = math.sin(t_opt) * math.cos(p_opt)
    ck = math.sin(t_opt) * math.sin(p_opt)
    psi_opt = ci * Di + cj * Dj + ck * Dk
    rho_AB = reduced_rho_AB(psi_opt, N)
    return {
        "i": i, "j": j, "k": k,
        "c_i": float(ci), "c_j": float(cj), "c_k": float(ck),
        "cpsi_opt": float(cpsi_opt),
        "purity_A_opt": float(single_qubit_purity(rho_AB)),
    }


def test_C_triple_dicke(N):
    print("=" * 80)
    print(f"Test C: triple-Dicke family maxima at N={N}")
    print("=" * 80)
    print(f"  All C({N+1}, 3) = {comb(N+1, 3)} triples.")
    print()

    results = []
    for i, j, k in combinations(range(N + 1), 3):
        r = triple_dicke_optimum(N, i, j, k)
        results.append(r)
    results.sort(key=lambda r: -r["cpsi_opt"])

    print(f"  Top 10 triples (sorted by cpsi_opt):")
    print(f"  {'triple':<14} {'(c_i,c_j,c_k)':<28} {'cpsi':>10} "
          f"{'purity':>9} {'>1/4?':>6}")
    print("  " + "-" * 75)
    for r in results[:10]:
        flag = "YES" if r["cpsi_opt"] > FOLD else "no"
        c_str = f"({r['c_i']:+.3f},{r['c_j']:+.3f},{r['c_k']:+.3f})"
        # Mark F69 if applicable: GHZ_N + W_N = (D_0+D_N)/sqrt(2) + D_1
        # Triple (0, 1, N) at c_0=c_N (GHZ structure) hits this
        label = ""
        if N == 3 and r["i"] == 0 and r["j"] == 1 and r["k"] == 3:
            label = " (F69 GHZ+W)"
        print(f"  D{r['i']}+D{r['j']}+D{r['k']:<6} {c_str:<28} {r['cpsi_opt']:>10.6f} "
              f"{r['purity_A_opt']:>9.4f} {flag:>6}{label}")
    print()
    above_fold = [r for r in results if r["cpsi_opt"] > FOLD]
    print(f"  {len(above_fold)}/{len(results)} triples lift above 1/4.")
    print()
    return results


# -----------------------------------------------------------------------
# Test B: full Hilbert space random scan at N=4
# -----------------------------------------------------------------------
def random_complex_unit(dim, n_samples, rng):
    """Uniform on the (2*dim - 1)-sphere of C^dim states."""
    re = rng.normal(size=(n_samples, dim))
    im = rng.normal(size=(n_samples, dim))
    psi = re + 1j * im
    norms = np.linalg.norm(psi, axis=1, keepdims=True)
    return psi / norms


def neg_cpsi_real_im(x, N, dim):
    """Negate pair-CPsi (for minimization), with x = [Re; Im] flat array."""
    n2 = float(np.dot(x, x))
    if n2 < 1e-24:
        return 0.0
    n = math.sqrt(n2)
    x_unit = x / n
    psi = x_unit[:dim] + 1j * x_unit[dim:]
    rho_AB = reduced_rho_AB(psi, N)
    return -pair_cpsi(rho_AB) + 1e-4 * (n - 1.0) ** 2


def refine_full(psi, N, maxiter=300):
    dim = 2 ** N
    x0 = np.concatenate([psi.real, psi.imag])
    res = minimize(neg_cpsi_real_im, x0, args=(N, dim),
                    method="L-BFGS-B",
                    options={"ftol": 1e-14, "gtol": 1e-10, "maxiter": maxiter})
    x = res.x / np.linalg.norm(res.x)
    return x[:dim] + 1j * x[dim:]


def hessian_negative_check(psi, N, n_perturb=20, eps=1e-4, rng=None):
    """Random perturbation test: at a local max, all small perturbations
    should DECREASE pair-CPsi (negative second derivative)."""
    if rng is None:
        rng = np.random.default_rng(42)
    dim = 2 ** N
    cpsi_0 = pair_cpsi(reduced_rho_AB(psi, N))
    n_increases = 0
    n_decreases = 0
    for _ in range(n_perturb):
        delta_re = rng.normal(size=dim)
        delta_im = rng.normal(size=dim)
        delta = delta_re + 1j * delta_im
        # Project out parallel component
        delta = delta - (psi.conj() @ delta) * psi
        delta_norm = np.linalg.norm(delta)
        if delta_norm < 1e-12:
            continue
        delta = delta / delta_norm
        psi_pert = psi + eps * delta
        psi_pert = psi_pert / np.linalg.norm(psi_pert)
        cpsi_p = pair_cpsi(reduced_rho_AB(psi_pert, N))
        if cpsi_p > cpsi_0 + 1e-9:
            n_increases += 1
        else:
            n_decreases += 1
    return n_increases, n_decreases


def test_B_full_hilbert(N, n_random=50000, n_refine_top=200, seed=42):
    print("=" * 80)
    print(f"Test B: full Hilbert space random scan at N={N}")
    print("=" * 80)
    print(f"  CP^{2**N - 1} (state dim = {2**N}, real dim = {2 * 2**N - 1})")
    print(f"  Random samples: {n_random}, refine top: {n_refine_top}")
    print()

    dim = 2 ** N
    rng = np.random.default_rng(seed)

    # Step 1: random scan
    t0 = time.time()
    psi_samples = random_complex_unit(dim, n_random, rng)
    cpsis = np.empty(n_random)
    purities = np.empty(n_random)
    for i, psi in enumerate(psi_samples):
        rho_AB = reduced_rho_AB(psi, N)
        cpsis[i] = pair_cpsi(rho_AB)
        purities[i] = single_qubit_purity(rho_AB)
    print(f"  Random scan: {time.time() - t0:.1f} s")
    print(f"  Initial random pair-CPsi range: [{cpsis.min():.4f}, {cpsis.max():.4f}]")
    n_above_initial = int(np.sum(cpsis > FOLD))
    print(f"  Random samples above 1/4: {n_above_initial}/{n_random} "
          f"({100*n_above_initial/n_random:.1f}%)")
    n_above_nonproduct = int(np.sum((cpsis > FOLD)
                                       & (purities < PURITY_PRODUCT_THRESHOLD)))
    print(f"  ...and not near-product (purity < {PURITY_PRODUCT_THRESHOLD}): "
          f"{n_above_nonproduct}")
    print()

    # Step 2: refine top candidates (regardless of purity, to see what they
    # become — most should ascend to product manifold)
    print(f"  Refining top {n_refine_top} candidates...")
    sorted_idx = np.argsort(-cpsis)
    top_candidates = sorted_idx[:n_refine_top]
    refined = []
    t0 = time.time()
    for i, idx in enumerate(top_candidates):
        psi_init = psi_samples[idx]
        psi_opt = refine_full(psi_init, N)
        rho_AB = reduced_rho_AB(psi_opt, N)
        cpsi_o = pair_cpsi(rho_AB)
        purity_o = single_qubit_purity(rho_AB)
        refined.append({
            "init_idx": int(idx),
            "init_cpsi": float(cpsis[idx]),
            "init_purity": float(purities[idx]),
            "opt_cpsi": float(cpsi_o),
            "opt_purity": float(purity_o),
            "psi_opt": psi_opt,
        })
    print(f"  Refinement: {time.time() - t0:.1f} s")
    print()

    # Step 3: filter for non-product optima above 1/4
    nonproduct_above = [r for r in refined
                          if r["opt_purity"] < PURITY_PRODUCT_THRESHOLD
                          and r["opt_cpsi"] > FOLD]
    print(f"  After refinement:")
    n_above_refined = sum(1 for r in refined if r["opt_cpsi"] > FOLD)
    n_product_refined = sum(1 for r in refined
                              if r["opt_purity"] >= PURITY_PRODUCT_THRESHOLD)
    print(f"    refined above 1/4 (any purity): {n_above_refined}/{n_refine_top}")
    print(f"    refined to product manifold:    {n_product_refined}/{n_refine_top}")
    print(f"    NON-product above 1/4:          {len(nonproduct_above)}")
    print()

    # Step 4: Hessian / stability test on non-product survivors
    if nonproduct_above:
        print(f"  Stability check on {len(nonproduct_above)} non-product survivors:")
        print(f"  {'rank':>5} {'cpsi':>10} {'purity':>9} {'increases':>10} "
              f"{'decreases':>10}  verdict")
        survivors_with_stability = []
        for k, r in enumerate(nonproduct_above[:20]):
            n_inc, n_dec = hessian_negative_check(
                r["psi_opt"], N, n_perturb=50, eps=1e-4, rng=rng,
            )
            verdict = ("LOCAL MAX" if n_inc == 0
                        else "SADDLE" if n_inc < n_dec
                        else "UNSTABLE")
            print(f"  {k+1:>5d} {r['opt_cpsi']:>10.6f} "
                  f"{r['opt_purity']:>9.6f} {n_inc:>10d} {n_dec:>10d}  "
                  f"{verdict}")
            survivors_with_stability.append({
                **{k: v for k, v in r.items() if k != "psi_opt"},
                "n_increases": n_inc, "n_decreases": n_dec,
                "verdict": verdict,
            })
    else:
        survivors_with_stability = []
        print(f"  No non-product candidates above 1/4 to test.")
    print()

    return {
        "n_random": n_random,
        "n_refine_top": n_refine_top,
        "n_above_initial": n_above_initial,
        "n_above_nonproduct_initial": n_above_nonproduct,
        "n_above_refined": n_above_refined,
        "n_product_refined": n_product_refined,
        "nonproduct_above_count": len(nonproduct_above),
        "survivors": survivors_with_stability,
    }


def main():
    all_A = {}
    all_C = {}
    for N in [3, 4, 5, 6]:
        print("\n" + "#" * 80)
        print(f"# N = {N}")
        print("#" * 80 + "\n")
        all_A[f"N={N}"] = test_A_binary_dicke(N)
        all_C[f"N={N}"] = test_C_triple_dicke(N)

    # Cross-N summary
    print("\n" + "=" * 80)
    print("CROSS-N SUMMARY")
    print("=" * 80)
    print(f"\n  {'N':>3}  {'binary above 1/4':>18}  {'binary max':>12}  "
          f"{'triple above 1/4':>18}  {'triple max':>12}")
    print("  " + "-" * 80)
    for N in [3, 4, 5, 6]:
        bin_res = all_A[f"N={N}"]
        tri_res = all_C[f"N={N}"]
        n_bin = sum(1 for r in bin_res if r["cpsi_opt"] > FOLD)
        n_tri = sum(1 for r in tri_res if r["cpsi_opt"] > FOLD)
        max_bin = max(r["cpsi_opt"] for r in bin_res)
        max_tri = max(r["cpsi_opt"] for r in tri_res)
        print(f"  {N:>3d}  {n_bin:>10d}/{len(bin_res):<6}  {max_bin:>12.4f}  "
              f"{n_tri:>10d}/{len(tri_res):<6}  {max_tri:>12.4f}")

    p = RESULTS_DIR / "eq016_dicke_saddle_landscape.json"
    with open(p, "w") as f:
        json.dump({"binary": all_A, "triple": all_C}, f, indent=2, default=str)
    print(f"\nSaved: {p}")


if __name__ == "__main__":
    main()
