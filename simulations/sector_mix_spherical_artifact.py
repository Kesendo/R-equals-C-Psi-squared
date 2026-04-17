#!/usr/bin/env python3
"""
Diagnosis: the 3-state spherical scan 'optimum' is a product state.

When cpsi_sector_mix_optimization.py Part 2 scans the 3-state spherical
family

    |psi(theta, phi)> = sin(theta) cos(phi) |GHZ_3>
                      + sin(theta) sin(phi) |W_3>
                      + cos(theta) |W_bar_3>

it finds a peak at approximately (theta*, phi*) = (0.9121, 0.8856) with
min pair-CPsi(0) ~ 0.9998. This is not a new entangled state; it is a
disguised copy of |+>^3 that happens to lie in the 3-state subspace
spanned by {|GHZ>, |W>, |W_bar>}. A future Claude reading the optimizer
output without this context will mistake the artifact for a genuine
finding. This script exists to catch that mistake in 30 seconds.

Checks:
    1. Grid-scan to reproduce the peak.
    2. Extract amplitudes of |psi*>  => all ~ 1/sqrt(8).
    3. Single-qubit purities Tr(rho_q^2) for q = 0, 1, 2.
       Product state iff all three are 1.
    4. 3-tangle tau_ABC. Product state has tau = 0.
    5. Overlap |<+^3 | psi*>|^2. Product state match iff 1.
"""
from __future__ import annotations
import itertools
import math
import sys
import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def ket_ghz(n):
    v = np.zeros(2**n, dtype=complex); v[0] = 1/np.sqrt(2); v[-1] = 1/np.sqrt(2); return v

def ket_w(n):
    v = np.zeros(2**n, dtype=complex)
    for i in range(n): v[1 << i] = 1/np.sqrt(n)
    return v

def ket_w_bar(n):
    v = np.zeros(2**n, dtype=complex); mask = (1 << n) - 1
    for i in range(n): v[mask ^ (1 << i)] = 1/np.sqrt(n)
    return v

def ket_plus_N(n):
    return np.ones(2**n, dtype=complex) / np.sqrt(2**n)

def partial_trace_pair(rho, n, keep):
    a, b = keep
    cur = list(range(n)); t = rho.reshape([2]*(2*n))
    while len(cur) > 2:
        rm = next(q for q in cur if q not in keep)
        idx = cur.index(rm); nc = len(cur)
        t = np.trace(t, axis1=nc-1-idx, axis2=2*nc-1-idx)
        cur.pop(idx)
    if cur == [a, b]: return t.reshape(4, 4)
    return t.transpose(1, 0, 3, 2).reshape(4, 4)

def partial_trace_single(rho, n, keep):
    cur = list(range(n)); t = rho.reshape([2]*(2*n))
    while len(cur) > 1:
        rm = next(q for q in cur if q != keep)
        idx = cur.index(rm); nc = len(cur)
        t = np.trace(t, axis1=nc-1-idx, axis2=2*nc-1-idx)
        cur.pop(idx)
    return t.reshape(2, 2)


def cpsi_pair(rho2):
    C = float(np.real(np.trace(rho2 @ rho2)))
    diag = np.diag(np.diag(rho2))
    L1 = float(np.sum(np.abs(rho2 - diag)))
    return C * L1 / 3.0

def min_pair_cpsi(psi, n):
    rho = np.outer(psi, psi.conj())
    return min(cpsi_pair(partial_trace_pair(rho, n, (a, b)))
               for a, b in itertools.combinations(range(n), 2))

def concurrence_2q(rho2):
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex); YY = np.kron(Y, Y)
    eigs = np.sort(np.real(np.linalg.eigvals(rho2 @ (YY @ rho2.conj() @ YY))))[::-1]
    eigs = np.clip(eigs, 0.0, None); s = np.sqrt(eigs)
    return float(max(0.0, s[0] - s[1] - s[2] - s[3]))

def three_tangle(psi):
    n = 3
    rho = np.outer(psi, psi.conj())
    rho_A = partial_trace_single(rho, n, 0)
    lin_entropy = 2.0 * (1.0 - np.real(np.trace(rho_A @ rho_A)))
    c_AB = concurrence_2q(partial_trace_pair(rho, n, (0, 1)))
    c_AC = concurrence_2q(partial_trace_pair(rho, n, (0, 2)))
    return max(0.0, float(lin_entropy - c_AB**2 - c_AC**2))


def find_peak():
    """Two-stage scan: coarse then fine around the coarse optimum."""
    n = 3
    ghz, w, w_bar = ket_ghz(n), ket_w(n), ket_w_bar(n)

    best = {"min": -1.0}
    for th in np.linspace(0, math.pi, 201):
        for ph in np.linspace(0, math.pi / 2, 101):
            psi = (math.sin(th)*math.cos(ph)*ghz
                   + math.sin(th)*math.sin(ph)*w
                   + math.cos(th)*w_bar)
            nrm = np.linalg.norm(psi)
            if nrm < 1e-9: continue
            psi /= nrm
            m = min_pair_cpsi(psi, n)
            if m > best["min"]: best = {"theta": th, "phi": ph, "min": m}

    th0, ph0 = best["theta"], best["phi"]
    fine = dict(best)
    for dth in np.linspace(-0.05, 0.05, 101):
        for dph in np.linspace(-0.05, 0.05, 101):
            th, ph = th0 + dth, ph0 + dph
            psi = (math.sin(th)*math.cos(ph)*ghz
                   + math.sin(th)*math.sin(ph)*w
                   + math.cos(th)*w_bar)
            nrm = np.linalg.norm(psi)
            if nrm < 1e-9: continue
            psi /= nrm
            m = min_pair_cpsi(psi, n)
            if m > fine["min"]: fine = {"theta": th, "phi": ph, "min": m}

    th_star, ph_star = fine["theta"], fine["phi"]
    psi_star = (math.sin(th_star)*math.cos(ph_star)*ghz
                + math.sin(th_star)*math.sin(ph_star)*w
                + math.cos(th_star)*w_bar)
    psi_star /= np.linalg.norm(psi_star)
    return th_star, ph_star, fine["min"], psi_star


def diagnose(psi_star):
    """Run the three independent product-state tests on psi_star."""
    n = 3
    rho = np.outer(psi_star, psi_star.conj())

    # Test 1: single-qubit purities
    purities = []
    for q in range(n):
        rq = partial_trace_single(rho, n, q)
        purities.append(float(np.real(np.trace(rq @ rq))))

    # Test 2: 3-tangle
    tau = three_tangle(psi_star)

    # Test 3: overlap with |+>^3
    plus3 = ket_plus_N(3)
    overlap_plus = float(abs(np.vdot(plus3, psi_star))**2)

    return {
        "purities": purities,
        "tau": tau,
        "overlap_plus": overlap_plus,
        "max_comp_basis_overlap": max(
            float(abs(psi_star[i])**2) for i in range(8)
        ),
    }


def main():
    from pathlib import Path
    RESULTS = Path(__file__).parent / "results"
    RESULTS.mkdir(exist_ok=True)
    OUT = RESULTS / "sector_mix_spherical_artifact.txt"

    with open(OUT, "w", encoding="utf-8") as fh:
        def log(msg=""):
            print(msg, flush=True)
            fh.write(msg + "\n")

        log("=" * 72)
        log("  SPHERICAL SCAN ARTIFACT DIAGNOSIS")
        log("=" * 72)
        log()
        log("  The 3-state spherical scan in cpsi_sector_mix_optimization.py")
        log("  Part 2 finds a peak at min pair-CPsi(0) ~ 0.9998. This is NOT")
        log("  a new entangled state; it is a disguised |+>^3 product state.")
        log()


        log("-- Reproducing the peak --")
        th, ph, peak, psi_star = find_peak()
        log(f"  theta* = {th:.6f}, phi* = {ph:.6f}")
        log(f"  min pair-CPsi(0) = {peak:.8f}")
        log(f"  (cpsi_sector_mix_optimization.py reports 0.999929)")
        log()

        log("-- Amplitudes of |psi*> in the computational basis --")
        for i, amp in enumerate(psi_star):
            log(f"    |{i:03b}>:  {amp.real:+.4f}  |.|^2 = {abs(amp)**2:.4f}")
        log(f"  uniform-amplitude reference: 1/sqrt(8) = {1/math.sqrt(8):.4f}")
        log()

        log("-- Test 1: single-qubit purities (product state iff all = 1) --")
        d = diagnose(psi_star)
        for q, p in enumerate(d["purities"]):
            tag = "PURE (product factor)" if p > 0.99 else ""
            log(f"    qubit {q}: Tr(rho^2) = {p:.6f}   {tag}")
        log()

        log("-- Test 2: CKW 3-tangle (product state iff tau = 0) --")
        log(f"    tau_ABC = {d['tau']:.6f}")
        log()

        log("-- Test 3: overlap with |+>^3 (product state iff overlap = 1) --")
        log(f"    |<+^3 | psi*>|^2 = {d['overlap_plus']:.6f}")
        log()

        log("=" * 72)
        log("  VERDICT")
        log("=" * 72)
        verdict_product = (
            all(p > 0.99 for p in d["purities"])
            and d["tau"] < 0.01
            and d["overlap_plus"] > 0.99
        )
        if verdict_product:
            log("  All three tests confirm: psi* is |+>^3, a product state.")
            log("  The 0.9998 peak is the coordinate artifact of the 3-state")
            log("  spherical family containing |+>^3 as a specific angle.")
            log("  NOT a new tripartite-entangled state above the fold.")
            log("  The genuine sector-mix optimum is F69 (GHZ+W at N=3, min")
            log("  pair-CPsi = 0.3204, tau = 0.80, cross-checked in")
            log("  ghz_w_optimum_n3.py).")
        else:
            log("  Tests disagree. Inspect purities / tangle / overlap above.")
        log()
        log(f"  Output: {OUT}")


if __name__ == "__main__":
    main()
