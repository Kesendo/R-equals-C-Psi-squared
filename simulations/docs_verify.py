"""
SCOPED NUMERICAL DOCUMENTATION SMOKE TEST
==========================================

This script recomputes only the explicitly implemented checks below. It is not
a completeness or cross-document consistency gate.

Authors: Tom Wicht, Claude
Date: March 16, 2026
"""

import numpy as np
import sys
from itertools import product as iprod
from framework import max_f1_pairing_distance

I2 = np.eye(2)
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])

def tensor(*ops):
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def site_op(op, s, N):
    ops = [I2]*N; ops[s] = op
    return tensor(*ops)

def build_H(N, J, topo, J_SA=None, J_SB=None):
    d = 2**N; H = np.zeros((d,d), dtype=complex)
    if topo == "chain":
        pairs = [(i,i+1,J) for i in range(N-1)]
    elif topo == "ring":
        pairs = [(i,(i+1)%N,J) for i in range(N)]
    elif topo == "star":
        J_SA = J_SA or J; J_SB = J_SB or J
        pairs = [(0,1,J_SA)] + [(0,i,J_SB) for i in range(2,N)]
    else: raise ValueError(topo)
    for i,j,Jij in pairs:
        for p in [sx,sy,sz]:
            H += Jij * site_op(p,i,N) @ site_op(p,j,N)
    return H

def build_L(H, gammas, N):
    d = 2**N; d2 = d*d
    L = -1j*(np.kron(H,np.eye(d)) - np.kron(np.eye(d),H.T))
    for k in range(N):
        Zk = site_op(sz,k,N)
        L += gammas[k]*(np.kron(Zk,Zk.conj()) - np.eye(d2))
    return L

def make_ghz(N):
    d=2**N; psi=np.zeros(d,complex); psi[0]=psi[-1]=1/np.sqrt(2)
    return np.outer(psi, psi.conj())

def make_w(N):
    d=2**N; psi=np.zeros(d,complex)
    for i in range(N): psi[1<<(N-1-i)] = 1/np.sqrt(N)
    return np.outer(psi, psi.conj())

def concurrence_2q(rho):
    """Wootters concurrence for 2-qubit state."""
    sy2 = np.kron(sy, sy)
    rho_tilde = sy2 @ rho.conj() @ sy2
    R = rho @ rho_tilde
    eigs = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, eigs[0] - eigs[1] - eigs[2] - eigs[3])

def partial_trace_keep01(rho, N):
    """Trace out all but first 2 qubits."""
    d = 2**N
    rho_r = rho.reshape([2]*2*N)
    for i in range(N-1, 1, -1):
        rho_r = np.trace(rho_r, axis1=i, axis2=i+N-(N-1-i+1)+(N-1-i))
    # Simpler: reshape and trace manually
    rho_r = rho.reshape(4, 2**(N-2), 4, 2**(N-2))
    return np.trace(rho_r, axis1=1, axis2=3)

PASS = 0; FAIL = 0

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1; print(f"  PASS: {name}")
    else:
        FAIL += 1; print(f"  FAIL: {name} -- {detail}")

if __name__ == "__main__":
    print("=" * 70)
    print("SCOPED NUMERICAL DOCUMENTATION SMOKE TEST")
    print("=" * 70)

    if "--force-failure" in sys.argv:
        check("forced failure proves the process exits nonzero", False,
              "intentional mutation regression")

    # =============================================================
    # GLOSSARY.md claims
    # =============================================================
    print("\n--- GLOSSARY: Boundary values ---")

    # Claim: "2gamma = decay rate of c+ supermode (N=3)"
    # Claim: "8gamma/3 = concurrence envelope (N=3)"
    # Claim: "10gamma/3 = decay rate of c- supermode (N=3)"
    gamma = 0.05
    N = 3
    H = build_H(N, 1.0, "chain")
    L = build_L(H, [gamma]*N, N)
    evals = np.linalg.eigvals(L)
    reals = np.real(evals)
    # Get unique non-zero real parts
    nonzero_reals = reals[np.abs(reals) > 1e-10]
    unique_rates = np.sort(np.unique(np.round(nonzero_reals, 8)))

    expected_rates = [-2*gamma, -8*gamma/3, -10*gamma/3, -2*3*gamma]
    print(f"  gamma = {gamma}, N = 3, chain")
    print(f"  Expected rates: {[f'{r:.6f}' for r in sorted(expected_rates)]}")
    print(f"  Unique rates found: {[f'{r:.6f}' for r in unique_rates[:8]]}")

    for rate_name, expected in [("2gamma", -2*gamma),
                                 ("8gamma/3", -8*gamma/3),
                                 ("10gamma/3", -10*gamma/3)]:
        found = any(abs(r - expected) < 1e-3 for r in unique_rates)
        check(f"{rate_name} = {expected:.6f} exists in spectrum (tol=1e-3)", found,
              f"expected {expected:.6f}")

    # These are shared approximate anchors, not equality of the rate multisets.
    print("\n--- GLOSSARY: Shared N=3 rate anchors across sampled topologies ---")
    for topo in ["chain", "ring", "star"]:
        H = build_H(3, 1.0, topo)
        L = build_L(H, [0.05]*3, 3)
        evals = np.linalg.eigvals(L)
        reals = np.real(evals)
        for rate_name, expected in [("2gamma", -0.10), ("8gamma/3", -8*0.05/3), ("10gamma/3", -10*0.05/3)]:
            found = any(abs(r - expected) < 1e-3 for r in reals)
            check(f"{topo} contains anchor {rate_name} (tol=1e-3)", found)

    wrong_anchor = -0.115
    check("mutated shared-anchor value is rejected",
          not any(abs(r - wrong_anchor) < 1e-3 for r in reals),
          f"unexpected anchor near {wrong_anchor}")

    # Claim: "-2*sum_gamma = location of XOR modes"
    print("\n--- GLOSSARY: XOR modes at -2*sum_gamma ---")
    for N in [2, 3, 4]:
        gammas = [0.05]*N
        H = build_H(N, 1.0, "chain")
        L = build_L(H, gammas, N)
        evals = np.linalg.eigvals(L)
        reals = np.real(evals)
        target = -2*sum(gammas)
        count = np.sum(np.abs(reals - target) < 1e-6)
        check(f"N={N}: {count} modes at -2*sum_gamma={target:.3f} (expect N+1={N+1})",
              count == N+1, f"got {count}")

    # Claim: "0.039/gamma = approximate crossing time"
    print("\n--- GLOSSARY: t_cross = 0.039/gamma ---")
    # Simulate Bell+ under Heisenberg chain N=2 with dephasing
    from scipy.linalg import expm
    gamma = 0.05
    H2 = build_H(2, 1.0, "chain")
    L2 = build_L(H2, [gamma]*2, 2)
    bell = np.array([1,0,0,1], dtype=complex)/np.sqrt(2)
    rho0 = np.outer(bell, bell.conj())
    rho_vec = rho0.flatten()

    t_cross_found = None
    for t in np.linspace(0.01, 2.0, 2000):
        rho_t = (expm(L2 * t) @ rho_vec).reshape(4,4)
        C = concurrence_2q(rho_t)
        # l1 coherence normalized
        l1 = np.sum(np.abs(rho_t)) - np.sum(np.abs(np.diag(rho_t)))
        Psi = l1 / 3  # d-1 = 3 for d=4
        CPsi = C * Psi
        if CPsi < 0.25 and t_cross_found is None:
            t_cross_found = t
            break

    expected_tcross = 0.039 / gamma
    if t_cross_found:
        ratio = t_cross_found * gamma  # Should be ~0.039
        check(f"t_cross = {t_cross_found:.3f}, gamma*t_cross = {ratio:.4f} (expect ~0.039, tol 10%)",
              abs(ratio - 0.039) / 0.039 < 0.10,
              f"gamma*t_cross = {ratio:.4f}, expect ~0.039")
    else:
        check("t_cross found", False, "never crossed 1/4")

    # =============================================================
    # MIRROR_SYMMETRY_PROOF claims
    # =============================================================
    print("\n--- MIRROR_SYMMETRY_PROOF: Palindromic pairing ---")

    # Claim: "every decay rate d paired with 2*sum_gamma - d"
    for N in [2, 3, 4, 5]:
        gammas = [0.05]*N
        sg = sum(gammas)
        H = build_H(N, 1.0, "chain")
        L = build_L(H, gammas, N)
        evals = np.linalg.eigvals(L)
        distance = max_f1_pairing_distance(evals, sg)
        tolerance = 1e-5
        check(f"N={N}: full eigenvalue multiset is palindromically paired",
              distance < tolerance, f"max multiplicity-aware distance={distance:.3e}")

    # Mutation-sensitive controls for the global matcher. The fixed point is
    # lambda=-Sigma, not -2*Sigma; and deleting the zero endpoint must expose
    # the unmatched -2*Sigma endpoint.
    control_sigma = 0.25
    fixed_point = np.array([-control_sigma + 0.2j, -control_sigma - 0.2j])
    check("F1 fixed-point control is centered at lambda=-Sigma",
          max_f1_pairing_distance(fixed_point, control_sigma) < 1e-14)
    endpoint_without_zero = np.array([-2 * control_sigma, -control_sigma, -control_sigma])
    endpoint_distance = max_f1_pairing_distance(endpoint_without_zero, control_sigma)
    check("F1 endpoint-without-zero negative control fails pairing",
          endpoint_distance > 0.1,
          f"broken multiset residual unexpectedly {endpoint_distance:.3e}")

    # Claim: the full F1 multiset closure holds on each supported topology.
    print("\n--- MIRROR_SYMMETRY_PROOF: Full multiset closure by topology ---")
    for topo in ["chain", "ring", "star"]:
        for N in [3, 4]:
            gammas = [0.05]*N
            sg = sum(gammas)
            H = build_H(N, 1.0, topo)
            L = build_L(H, gammas, N)
            evals = np.linalg.eigvals(L)
            distance = max_f1_pairing_distance(evals, sg)
            check(f"N={N} {topo}: full eigenvalue multiset has F1 closure",
                  distance < 1e-5,
                  f"max multiplicity-aware distance={distance:.3e}")

    # A Hamiltonian outside the supported Pi-odd class must not make the
    # topology loop false-green merely because all rates stay in range.
    bad_n = 3
    bad_gammas = [0.05] * bad_n
    bad_h = build_H(bad_n, 1.0, "chain") + 0.37 * site_op(sz, 0, bad_n)
    bad_spectrum = np.linalg.eigvals(build_L(bad_h, bad_gammas, bad_n))
    bad_distance = max_f1_pairing_distance(bad_spectrum, sum(bad_gammas))
    check("non-Pi-compatible Hamiltonian negative control breaks F1 closure",
          bad_distance > 1e-3,
          f"broken-control distance unexpectedly {bad_distance:.3e}")

    print("\n--- XOR_SPACE: coordinate claims retired ---")
    print("  Right-eigenvector coordinate squares are not invariant state weights.")
    print("  Run f22_operator_charge.py for the operator-level F22 gate.")

    # =============================================================
    # CORE_ALGEBRA claims
    # =============================================================
    print("\n--- CORE_ALGEBRA: 1/4 boundary ---")

    # Claim: "discriminant of R = C(Psi+R)^2 changes sign at C*Psi = 1/4"
    # R = C(Psi+R)^2 -> R = C*Psi^2 + 2C*Psi*R + C*R^2
    # -> C*R^2 + (2C*Psi - 1)*R + C*Psi^2 = 0
    # discriminant = (2C*Psi - 1)^2 - 4*C^2*Psi^2 = 1 - 4C*Psi
    for CPsi in [0.1, 0.24, 0.25, 0.26, 0.5]:
        disc = 1 - 4*CPsi
        if CPsi < 0.25:
            check(f"CPsi={CPsi}: disc={disc:.4f} > 0 (real fixed points)", disc > 0)
        elif CPsi == 0.25:
            check(f"CPsi=0.25: disc={disc:.4f} = 0 (boundary)", abs(disc) < 1e-10)
        else:
            check(f"CPsi={CPsi}: disc={disc:.4f} < 0 (complex, no real FP)", disc < 0)

    # Claim: "Mandelbrot equivalence: u = C(Psi+R) maps to z^2 + c with c = C*Psi"
    print("\n--- CORE_ALGEBRA: Mandelbrot equivalence ---")
    C, Psi = 0.8, 0.2
    R0 = 0.1
    # R iteration: R_{n+1} = C*(Psi + R_n)^2
    R1 = C * (Psi + R0)**2
    # u iteration: u = C(Psi+R), u_{n+1} = u_n^2 + c where c = C*Psi
    u0 = C * (Psi + R0)
    u1_direct = C * (Psi + R1)
    u1_mandelbrot = u0**2 + C*Psi
    check(f"Mandelbrot map: u1_direct={u1_direct:.8f} = u1_mandelbrot={u1_mandelbrot:.8f}",
          abs(u1_direct - u1_mandelbrot) < 1e-10,
          f"diff = {abs(u1_direct - u1_mandelbrot):.2e}")

    # =============================================================
    # WEAKNESSES claims (spot check)
    # =============================================================
    print("\n--- WEAKNESSES: J threshold: run verify_star_topology.py ---")

    # The star-topology threshold used to be "checked" here by passing the
    # literal True, with a note saying verify_star_topology.py had verified it.
    # A check that cannot fail is not a check, and the number it stood behind
    # (1.466) was itself a sampling artifact; the converged value is 1.46295
    # (F29). Run verify_star_topology.py for the gate that can fail.

    # =============================================================
    # CONSISTENCY CHECKS between documents
    # =============================================================
    print("\n--- CONSISTENCY: retired mixed-XY correlation ---")
    print("  The reported correlations described a non-invariant coordinate diagnostic.")

    # Claim: Echo peak C_SB = 0.598 for N=3
    # This would require a full star simulation which we've done in qst_bridge.py
    # Reference check only
    print("\n--- REFERENCE: Echo peak C_SB = 0.598 ---")
    print("  Verified in simulations/qst_bridge.py. Not recomputed here.")
    print("  Cross-reference: QST_BRIDGE.md and GLOSSARY.md both cite 0.598.")

    # Claim: F_avg = 0.888 for star 2:1
    print("\n--- REFERENCE: F_avg = 0.888 ---")
    print("  Verified in simulations/qst_bridge.py. Not recomputed here.")
    print("  Would require average over random input states (expensive).")
    print("  Cross-reference: QST_BRIDGE.md, GLOSSARY.md, WEAKNESSES.md all cite 0.888.")

    # =============================================================
    # SUMMARY
    # =============================================================
    print("\n" + "=" * 70)
    print(f"DOCS VERIFICATION SUMMARY: {PASS} PASS, {FAIL} FAIL")
    print("=" * 70)

    if FAIL > 0:
        print(f"\n{FAIL} FAILURES. Fix before any release.")
    else:
        print("\nAll checks implemented in this script passed.")
        print("This is a scoped numerical smoke test, not a documentation-consistency gate.")

    print("\nNOT RECOMPUTED (expensive, verified in dedicated scripts):")
    print("  - F_avg = 0.888 (qst_bridge.py)")
    print("  - Echo peak C_SB = 0.598 (qst_bridge.py)")
    print("  - Holevo capacity 0.534 bits (verify_channel.py)")
    print("  - IBM Torino results (hardware, not reproducible locally)")
    raise SystemExit(1 if FAIL else 0)
