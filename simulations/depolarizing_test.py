"""
DEPOLARIZING NOISE TEST
========================
Does the palindrome survive depolarizing noise (X+Y+Z simultaneously)?
This is the most realistic noise model for real hardware.

Output: simulations/results/depolarizing_test.txt
"""
import numpy as np
import os
from datetime import datetime

OUT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results\depolarizing_test.txt"
f = open(OUT, "w", buffering=1)

def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n"); f.flush()

I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
PM = {'X': sx, 'Y': sy, 'Z': sz}

def site_op(op, s, N):
    ops = [I2]*N; ops[s] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def build_H(N, pairs, comps):
    d = 2**N; H = np.zeros((d,d), dtype=complex)
    for i,j in pairs:
        for c, J in comps.items():
            if J == 0: continue
            H += J * site_op(PM[c[0]],i,N) @ site_op(PM[c[1]],j,N)
    return H

def build_L(H, gamma, N, noise_type='z'):
    """Build Liouvillian with different noise types.
    noise_type: 'z' = Z-dephasing, 'x' = X-dephasing, 'depol' = depolarizing (X+Y+Z)
    """
    d = 2**N; d2 = d*d; Id = np.eye(d)
    L = -1j*(np.kron(H, Id) - np.kron(Id, H.T))
    
    if noise_type == 'depol':
        # Depolarizing: three Lindblad operators per site (X, Y, Z)
        # Each with rate gamma/3 so total dephasing strength = gamma per site
        for k in range(N):
            for deph_op in [sx, sy, sz]:
                Dk = site_op(deph_op, k, N)
                L += (gamma/3)*(np.kron(Dk, Dk.conj()) - np.eye(d2))
    else:
        deph_op = {'z': sz, 'x': sx, 'y': sy}[noise_type]
        for k in range(N):
            Dk = site_op(deph_op, k, N)
            L += gamma*(np.kron(Dk, Dk.conj()) - np.eye(d2))
    return L

def check(L, sum_gamma):
    evals = np.linalg.eigvals(L)
    target = -2*sum_gamma
    nonzero = evals[np.abs(evals) > 1e-10]
    if len(nonzero) == 0: return 0, True
    max_err = 0
    for lam in nonzero:
        partner = -lam + target
        best = np.min(np.abs(evals - partner))
        if best > max_err: max_err = best
    return max_err, max_err < 1e-8

if __name__ == "__main__":
    gamma = 0.05
    
    log("=" * 80)
    log("DEPOLARIZING NOISE PALINDROME TEST")
    log(f"Started: {datetime.now()}")
    log("=" * 80)
    
    models = {
        'Heisenberg':    {'XX':1,'YY':1,'ZZ':1},
        'XY-only':       {'XX':1,'YY':1},
        'Ising':         {'ZZ':1},
        'XX alone':      {'XX':1},
        'DM (XY-YX)':   {'XY':1,'YX':-1},
        'XXZ delta=2':   {'XX':1,'YY':1,'ZZ':2},
        'Heis+DM':       {'XX':1,'YY':1,'ZZ':1,'XY':0.3,'YX':-0.3},
    }
    
    # ============================================================
    # TEST 1: Z-deph vs Depolarizing, N=3
    # ============================================================
    log("\n" + "=" * 80)
    log("TEST 1: Z-dephasing vs Depolarizing (N=3)")
    log(f"gamma={gamma} per site, depol uses gamma/3 per channel")
    log("=" * 80)
    
    N = 3; pairs = [(0,1),(1,2)]; sg = N*gamma
    
    log(f"\n  {'Model':<20} {'Z-deph':>12} {'Depol':>12}")
    log(f"  {'-'*48}")
    
    for name, comps in models.items():
        H = build_H(N, pairs, comps)
        
        L_z = build_L(H, gamma, N, 'z')
        err_z, ok_z = check(L_z, sg)
        
        L_d = build_L(H, gamma, N, 'depol')
        err_d, ok_d = check(L_d, sg)
        
        z_str = "OK" if ok_z else f"err={err_z:.2e}"
        d_str = "OK" if ok_d else f"err={err_d:.2e}"
        log(f"  {name:<20} {z_str:>12} {d_str:>12}")

    # ============================================================
    # TEST 2: Maybe depol has DIFFERENT center?
    # ============================================================
    log("\n" + "=" * 80)
    log("TEST 2: Center search for depolarizing noise")
    log("Maybe the palindrome axis is NOT at -2*sum_gamma for depol")
    log("=" * 80)
    
    H = build_H(N, pairs, {'XX':1,'YY':1,'ZZ':1})
    L_d = build_L(H, gamma, N, 'depol')
    evals = np.linalg.eigvals(L_d)
    nonzero = evals[np.abs(evals) > 1e-10]
    reals = np.sort(np.real(nonzero))
    
    log(f"\n  Depol spectrum (Heisenberg N=3):")
    log(f"  Min Re: {reals[0]/gamma:.4f}g")
    log(f"  Max Re: {reals[-1]/gamma:.4f}g")
    log(f"  Expected center (Z-deph): {-sg/gamma:.4f}g")
    
    # The center of palindrome = (min+max)/2
    empirical_center = (reals[0] + reals[-1]) / 2
    log(f"  Empirical center: {empirical_center/gamma:.4f}g")
    
    # Test palindrome around empirical center
    err_emp, ok_emp = check(L_d, -empirical_center)  # note: check expects sum_gamma
    log(f"  Palindrome around empirical center: {'OK' if ok_emp else f'err={err_emp:.2e}'}")
    
    # Also try center = sum of ALL dephasing rates
    # For depol: 3 channels per site at gamma/3 each = gamma per site = N*gamma total
    # But maybe it's 3*N*gamma/3 = N*gamma... same thing
    # Or maybe it's the trace of the dissipator?
    
    # Actually for depol with rate gamma/3 per channel:
    # total rate per site = 3 * (gamma/3) = gamma
    # So sum_gamma should still be N*gamma
    # But the DISSIPATOR structure is different!
    
    # Try different centers
    log(f"\n  Center sweep:")
    for factor in [0.5, 0.667, 0.75, 1.0, 1.333, 1.5, 2.0, 3.0]:
        test_sg = factor * N * gamma
        err, ok = check(L_d, test_sg)
        tag = "OK" if ok else f"err={err:.2e}"
        log(f"    center={factor:.3f}*N*gamma: {tag}")

    # ============================================================
    # TEST 3: N=4 confirmation
    # ============================================================
    log("\n" + "=" * 80)
    log("TEST 3: N=4 confirmation")
    log("=" * 80)
    
    N4 = 4; pairs4 = [(0,1),(1,2),(2,3)]; sg4 = N4*gamma
    
    log(f"\n  {'Model':<20} {'Z-deph':>12} {'Depol':>12}")
    log(f"  {'-'*48}")
    
    for name in ['Heisenberg', 'XY-only', 'Ising', 'DM (XY-YX)']:
        comps = models[name]
        H = build_H(N4, pairs4, comps)
        
        L_z = build_L(H, gamma, N4, 'z')
        err_z, ok_z = check(L_z, sg4)
        
        L_d = build_L(H, gamma, N4, 'depol')
        err_d, ok_d = check(L_d, sg4)
        
        z_str = "OK" if ok_z else f"err={err_z:.2e}"
        d_str = "OK" if ok_d else f"err={err_d:.2e}"
        log(f"  {name:<20} {z_str:>12} {d_str:>12}")
    
    # Also try empirical center for N=4 depol
    H4 = build_H(N4, pairs4, {'XX':1,'YY':1,'ZZ':1})
    L4d = build_L(H4, gamma, N4, 'depol')
    ev4 = np.linalg.eigvals(L4d)
    nz4 = ev4[np.abs(ev4) > 1e-10]
    r4 = np.sort(np.real(nz4))
    emp4 = (r4[0] + r4[-1]) / 2
    log(f"\n  N=4 Depol empirical center: {emp4/gamma:.4f}g")
    err4, ok4 = check(L4d, -emp4)
    log(f"  Palindrome at empirical center: {'OK' if ok4 else f'err={err4:.2e}'}")

    # ============================================================
    # TEST 4: Amplitude damping (non-dephasing noise)
    # ============================================================
    log("\n" + "=" * 80)
    log("TEST 4: Amplitude damping (non-dephasing, for comparison)")
    log("L_k = sqrt(gamma) * sigma_minus = sqrt(gamma) * (X - iY)/2")
    log("=" * 80)
    
    N = 3; pairs = [(0,1),(1,2)]; sg = N*gamma
    sigma_minus = (sx - 1j*sy) / 2
    
    H = build_H(N, pairs, {'XX':1,'YY':1,'ZZ':1})
    d = 2**N; d2 = d*d; Id = np.eye(d)
    L_amp = -1j*(np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Lk = np.sqrt(gamma) * site_op(sigma_minus, k, N)
        LkD = Lk.conj().T
        LdL = LkD @ Lk
        L_amp += np.kron(Lk, Lk.conj()) - 0.5*(np.kron(LdL,Id) + np.kron(Id,LdL.T))
    
    err_amp, ok_amp = check(L_amp, sg)
    log(f"  Amplitude damping: {'OK' if ok_amp else f'err={err_amp:.2e}'}")
    
    # Try empirical center
    ev_amp = np.linalg.eigvals(L_amp)
    nz_amp = ev_amp[np.abs(ev_amp) > 1e-10]
    r_amp = np.sort(np.real(nz_amp))
    emp_amp = (r_amp[0] + r_amp[-1]) / 2
    log(f"  Empirical center: {emp_amp/gamma:.4f}g")
    err_amp2, ok_amp2 = check(L_amp, -emp_amp)
    log(f"  At empirical center: {'OK' if ok_amp2 else f'err={err_amp2:.2e}'}")
    log(f"  Min: {r_amp[0]/gamma:.4f}g  Max: {r_amp[-1]/gamma:.4f}g")
    
    log(f"\n{'='*80}")
    log(f"Completed: {datetime.now()}")
    log(f"{'='*80}")
    f.close()
    print(f"\n>>> Results: {OUT}")
