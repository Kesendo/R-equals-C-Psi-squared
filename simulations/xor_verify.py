"""
XOR Space Verification: Pre-release checks
============================================

Before we publish, every claim must be verified across:
- Multiple system sizes (N=2,3,4)
- Multiple topologies (chain, ring, star)
- The correlation r=0.976 for other N

Claims to verify:
1. Always N+1 center modes
2. GHZ -> 100% XOR for all N
3. W -> 100% palindrome for all N  
4. Mixed XY correlation holds for all N
5. XOR modes are off-diagonal for all N
"""

import numpy as np
from itertools import product as iprod

I2 = np.eye(2)
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])

def tensor(*ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r

def site_op(op, s, N):
    ops = [I2]*N; ops[s] = op
    return tensor(*ops)

def build_H(N, J, topo):
    d = 2**N; H = np.zeros((d,d), dtype=complex)
    if topo == "chain": pairs = [(i,i+1) for i in range(N-1)]
    elif topo == "ring": pairs = [(i,(i+1)%N) for i in range(N)]
    elif topo == "star": pairs = [(0,i) for i in range(1,N)]
    else: raise ValueError(topo)
    for i,j in pairs:
        for p in [sx,sy,sz]:
            H += J * site_op(p,i,N) @ site_op(p,j,N)
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
    return np.outer(psi,psi.conj())

def make_w(N):
    d=2**N; psi=np.zeros(d,complex)
    for i in range(N): psi[1<<(N-1-i)] = 1/np.sqrt(N)
    return np.outer(psi,psi.conj())

def make_product(N, s):
    m = {'0':np.array([1,0],complex), '1':np.array([0,1],complex),
         '+': np.array([1,1],complex)/np.sqrt(2), '-': np.array([1,-1],complex)/np.sqrt(2)}
    psi = m[s[0]]
    for c in s[1:]: psi = np.kron(psi, m[c])
    return np.outer(psi, psi.conj())

def pauli_xy_weight(rho, N):
    paulis = [I2, sx, sy, sz]
    d = 2**N; mixed_w = 0.0; total_w = 0.0
    for indices in iprod(range(4), repeat=N):
        P = paulis[indices[0]]
        for idx in indices[1:]: P = np.kron(P, paulis[idx])
        c = np.trace(P @ rho) / d
        w = np.abs(c)**2; total_w += w
        has_x = any(i==1 for i in indices)
        has_y = any(i==2 for i in indices)
        if has_x and has_y: mixed_w += w
    return mixed_w / total_w if total_w > 0 else 0

def decompose_and_measure(L, rho0, gammas, tol=1e-6):
    """Returns (xor_fraction, pal_fraction, xor_modes_are_offdiag)."""
    sg = sum(gammas); ts = -2*sg
    evals, rvecs = np.linalg.eig(L)
    lvecs = np.linalg.inv(rvecs)
    coeffs = lvecs @ rho0.flatten()
    weights = np.abs(coeffs)**2
    reals = np.real(evals)
    
    # Pair finding
    used = set(); paired = {}
    for i in range(len(evals)):
        if i in used or np.abs(evals[i])<1e-12: continue
        target = ts - reals[i]
        bj, bd = None, tol
        for j in range(len(evals)):
            if j in used or j==i or np.abs(evals[j])<1e-12: continue
            diff = abs(reals[j]-target)
            if diff < bd: bd=diff; bj=j
        if bj is not None:
            paired[i]=bj; paired[bj]=i; used.add(i); used.add(bj)
    
    pw=0; xw=0
    xor_diag_weights = []
    d = int(np.sqrt(L.shape[0]))
    for i in range(len(evals)):
        if np.abs(evals[i])<1e-12: continue
        if weights[i]<1e-15: continue
        if i in paired:
            pw += weights[i]
        else:
            xw += weights[i]
            # Check if mode is off-diagonal
            mode_op = rvecs[:,i].reshape(d,d)
            diag_w = np.sum(np.abs(np.diag(mode_op))**2) / np.sum(np.abs(mode_op)**2)
            xor_diag_weights.append(diag_w)
    
    total = pw+xw
    xf = xw/total if total>0 else 0
    pf = pw/total if total>0 else 0
    avg_diag = np.mean(xor_diag_weights) if xor_diag_weights else -1
    return xf, pf, avg_diag

if __name__ == "__main__":
    print("=" * 70)
    print("XOR SPACE VERIFICATION: Pre-release checks")
    print("=" * 70)
    
    PASS = 0; FAIL = 0
    
    def check(name, condition, detail=""):
        global PASS, FAIL
        if condition:
            PASS += 1
            print(f"  PASS: {name}")
        else:
            FAIL += 1
            print(f"  FAIL: {name} -- {detail}")
    
    # =========================================================
    # CHECK 1: N+1 center modes for N=2,3,4,5
    # =========================================================
    print("\n--- CHECK 1: Always N+1 center modes ---")
    
    for N in [2, 3, 4, 5]:
        gammas = [0.05]*N
        H = build_H(N, 1.0, "chain")
        L = build_L(H, gammas, N)
        evals = np.linalg.eigvals(L)
        sg = sum(gammas); ts = -2*sg
        
        nonzero = evals[np.abs(evals) > 1e-10]
        reals = np.real(nonzero)
        
        # Count modes at the center (within tolerance)
        center_count = np.sum(np.abs(reals - ts) < 1e-6)
        expected = N + 1
        check(f"N={N} chain: center modes = {center_count} (expect {expected})",
              center_count == expected,
              f"got {center_count}")

    # Also check different topologies
    for topo in ["ring", "star"]:
        N = 3; gammas = [0.05]*N
        H = build_H(N, 1.0, topo)
        L = build_L(H, gammas, N)
        evals = np.linalg.eigvals(L)
        nonzero = evals[np.abs(evals) > 1e-10]
        reals = np.real(nonzero)
        ts = -2*sum(gammas)
        center_count = np.sum(np.abs(reals - ts) < 1e-6)
        check(f"N=3 {topo}: center modes = {center_count} (expect {N+1})",
              center_count == N+1, f"got {center_count}")
    
    # Non-uniform gammas
    for gammas in [[0.05, 0.10, 0.05], [0.01, 0.05, 0.20]]:
        N = 3; sg = sum(gammas); ts = -2*sg
        H = build_H(N, 1.0, "chain")
        L = build_L(H, gammas, N)
        evals = np.linalg.eigvals(L)
        nonzero = evals[np.abs(evals) > 1e-10]
        reals = np.real(nonzero)
        center_count = np.sum(np.abs(reals - ts) < 1e-6)
        check(f"N=3 gammas={gammas}: center={center_count} (expect {N+1})",
              center_count == N+1, f"got {center_count}")
    
    # =========================================================
    # CHECK 2: GHZ -> 100% XOR for all N and topologies
    # =========================================================
    print("\n--- CHECK 2: GHZ always 100% XOR ---")
    
    for N in [2, 3, 4]:
        for topo in ["chain", "star"]:
            gammas = [0.05]*N
            H = build_H(N, 1.0, topo)
            L = build_L(H, gammas, N)
            xf, pf, diag = decompose_and_measure(L, make_ghz(N), gammas)
            check(f"GHZ N={N} {topo}: XOR={xf:.3f} (expect 1.0)",
                  abs(xf - 1.0) < 0.01 or pf < 0.01,
                  f"XOR={xf:.3f}, PAL={pf:.3f}")

    # =========================================================
    # CHECK 3: W -> 100% palindromic for all N and topologies
    # =========================================================
    print("\n--- CHECK 3: W always 100% palindromic ---")
    
    for N in [2, 3, 4]:
        for topo in ["chain", "star"]:
            gammas = [0.05]*N
            H = build_H(N, 1.0, topo)
            L = build_L(H, gammas, N)
            xf, pf, diag = decompose_and_measure(L, make_w(N), gammas)
            check(f"W N={N} {topo}: PAL={pf:.3f} (expect 1.0)",
                  abs(pf - 1.0) < 0.01 or xf < 0.01,
                  f"XOR={xf:.3f}, PAL={pf:.3f}")
    
    # =========================================================
    # CHECK 4: XOR modes are off-diagonal (coherences)
    # =========================================================
    print("\n--- CHECK 4: XOR modes are off-diagonal ---")
    
    for N in [2, 3, 4]:
        gammas = [0.05]*N
        H = build_H(N, 1.0, "chain")
        L = build_L(H, gammas, N)
        xf, pf, avg_diag = decompose_and_measure(L, make_ghz(N), gammas)
        check(f"GHZ N={N}: XOR mode diag_weight = {avg_diag:.4f} (expect ~0)",
              avg_diag < 0.01 or avg_diag < 0,
              f"diag_weight={avg_diag:.4f}")

    # =========================================================
    # CHECK 5: Mixed XY correlation holds for N=2 and N=3
    # =========================================================
    print("\n--- CHECK 5: Mixed XY Pauli weight predicts XOR ---")
    
    for N in [2, 3]:
        gammas = [0.05]*N
        H = build_H(N, 1.0, "chain")
        L = build_L(H, gammas, N)
        
        # Build state set
        if N == 2:
            states = {
                "00": make_product(2,"00"), "01": make_product(2,"01"),
                "+0": make_product(2,"+0"), "++": make_product(2,"++"),
                "+-": make_product(2,"+-"), "GHZ": make_ghz(2), "W": make_w(2),
            }
        else:
            states = {
                "000": make_product(3,"000"), "010": make_product(3,"010"),
                "+00": make_product(3,"+00"), "+++": make_product(3,"+++"),
                "+-+": make_product(3,"+-+"), "GHZ": make_ghz(3), "W": make_w(3),
            }
        
        xor_fracs = []
        xy_weights = []
        for name, rho in states.items():
            xf, pf, _ = decompose_and_measure(L, rho, gammas)
            xyw = pauli_xy_weight(rho, N)
            xor_fracs.append(xf)
            xy_weights.append(xyw)
        
        xor_arr = np.array(xor_fracs)
        xy_arr = np.array(xy_weights)
        if np.std(xor_arr) > 1e-10 and np.std(xy_arr) > 1e-10:
            corr = np.corrcoef(xor_arr, xy_arr)[0,1]
        else:
            corr = 0.0
        
        check(f"N={N}: XY-XOR correlation r={corr:.3f} (expect > 0.9)",
              corr > 0.9,
              f"r={corr:.3f}")

    # =========================================================
    # CHECK 6: Center modes are at EXACTLY -2*sum_gamma
    # =========================================================
    print("\n--- CHECK 6: XOR modes sit at -2*sum_gamma (max decay rate) ---")
    print("    NOTE: This is the maximum decay, NOT the midpoint of palindrome.")
    print("    Midpoint would be -sum_gamma. These modes pair with steady state (0).")
    
    for N in [2, 3, 4]:
        gammas = [0.05]*N
        H = build_H(N, 1.0, "chain")
        L = build_L(H, gammas, N)
        evals = np.linalg.eigvals(L)
        sg = sum(gammas)
        
        # Find the unpaired modes
        nonzero = evals[np.abs(evals) > 1e-10]
        reals = np.real(nonzero)
        target = -2*sg
        
        center_modes = nonzero[np.abs(reals - target) < 1e-6]
        if len(center_modes) > 0:
            actual_rate = np.mean(np.real(center_modes))
            check(f"N={N}: center modes at {actual_rate:.6f} (expect {target:.6f})",
                  abs(actual_rate - target) < 1e-5,
                  f"actual={actual_rate:.6f}, expected={target:.6f}")
        
        # Also check: is -2*sum_gamma the MAXIMUM decay rate?
        max_decay = np.min(reals)  # Most negative = fastest decay
        check(f"N={N}: max decay rate = {max_decay:.6f} (expect {target:.6f})",
              abs(max_decay - target) < 1e-5,
              f"max_decay={max_decay:.6f}, expected={target:.6f}")
    
    # =========================================================
    # CHECK 7: Document accuracy check
    # =========================================================
    print("\n--- CHECK 7: Claims in XOR_SPACE.md ---")
    
    # The doc says "They sit at the exact midpoint: lambda = -2*sum_gamma"
    # But -2*sum_gamma is NOT the midpoint, it's the MAXIMUM DECAY RATE.
    # The midpoint of the palindrome is -sum_gamma.
    # These modes pair with steady state (lambda=0), so they're at the EDGE.
    print("  WARNING: XOR_SPACE.md says 'midpoint' but -2*sum_gamma is the")
    print("  MAXIMUM DECAY RATE, not the midpoint. The midpoint is -sum_gamma.")
    print("  These modes pair with the STEADY STATE (lambda=0).")
    print("  The document needs correction: 'midpoint' -> 'maximum decay rate'")
    FAIL += 1  # Flag for document correction

    # =========================================================
    # SUMMARY
    # =========================================================
    print("\n" + "=" * 70)
    print(f"VERIFICATION SUMMARY: {PASS} PASS, {FAIL} FAIL")
    print("=" * 70)
    
    if FAIL > 0:
        print("\nFAILURES MUST BE RESOLVED BEFORE RELEASE.")
        print("DO NOT publish until all checks pass.")
    else:
        print("\nAll checks passed. Ready for release.")
    
    print("\nDOCUMENT CORRECTIONS NEEDED:")
    print("1. XOR_SPACE.md says 'midpoint' but should say 'maximum decay rate'")
    print("   The N+1 modes are at -2*sum_gamma (fastest decay).")
    print("   Their palindromic partner is the steady state (lambda=0).")
    print("   They are at the EDGE of the palindrome, not the CENTER.")
    print("   The palindrome analogy should be: in ABCBA, these are NOT C.")
    print("   They are like A paired with a silent letter outside the word.")
