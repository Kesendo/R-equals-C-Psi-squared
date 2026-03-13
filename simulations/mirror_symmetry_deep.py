"""
MIRROR SYMMETRY DEEP DIVE
==========================
The Liouvillian decay spectrum is exactly symmetric around Ng.
WHY? This script investigates:

1. Does it hold for non-uniform dephasing (different gamma per qubit)?
2. Does it hold for non-Z noise (X or Y dephasing)?
3. Does it hold for mixed noise types?
4. What is the conjugation operator that maps rate d -> 2Ng - d?
5. What do mirror pairs correspond to physically (Pauli weight)?
6. Does it hold for non-Heisenberg coupling (XXZ, XY models)?

Results -> results_mirror_symmetry.txt
"""
import numpy as np
import math
import time

I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)

def op_at(op, qubit, n_q):
    ops = [I2]*n_q
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result

def build_liouvillian_general(n_q, bonds, jump_ops):
    """General Liouvillian with arbitrary jump operators."""
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for i, j, J, paulis in bonds:
        for p in paulis:
            H += J * op_at(p,i,n_q) @ op_at(p,j,n_q)
    I_d = np.eye(d, dtype=complex)
    L_mat = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for Lk in jump_ops:
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L_mat += np.kron(Lk, Lk.conj()) - 0.5*(np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    return L_mat

def heisenberg_bonds(n_q, couplings):
    """Star topology with isotropic Heisenberg coupling."""
    return [(0, i, couplings[i-1], [X,Y,Z]) for i in range(1, n_q)]

def xxz_bonds(n_q, couplings, delta=0.5):
    """Star with XXZ anisotropy: J(XX+YY) + delta*J*ZZ."""
    bonds = []
    for i in range(1, n_q):
        J = couplings[i-1]
        bonds.append((0, i, J, [X, Y]))
        bonds.append((0, i, J*delta, [Z]))
    return bonds

def xy_bonds(n_q, couplings):
    """Star with XY coupling only (no ZZ)."""
    return [(0, i, couplings[i-1], [X, Y]) for i in range(1, n_q)]

def get_osc_rates_raw(L_mat, threshold=0.05):
    """Get oscillatory decay rates (raw, not normalized)."""
    evals = np.linalg.eigvals(L_mat)
    rates = []
    for ev in evals:
        if abs(np.imag(ev)) > threshold:
            rate = -np.real(ev)
            if rate > 0.001:
                rates.append(round(rate, 6))
    return sorted(rates)

def check_mirror_symmetry(rates, center, tol=0.005):
    """Check if rates are symmetric around center. Returns score and details."""
    below = [r for r in rates if r < center - tol]
    above = [r for r in rates if r > center + tol]
    at_center = [r for r in rates if abs(r - center) < tol]
    
    matched = 0
    unmatched = []
    for r in below:
        mirror = 2*center - r
        closest = min(above, key=lambda a: abs(a - mirror)) if above else 999
        if abs(closest - mirror) < tol:
            matched += 1
        else:
            unmatched.append(r)
    
    score = matched / max(len(below), 1) if below else 1.0
    return {
        'score': score,
        'below': len(below),
        'above': len(above),
        'at_center': len(at_center),
        'matched': matched,
        'unmatched': len(unmatched),
        'center': center,
    }

outfile = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results_mirror_symmetry.txt"
f = open(outfile, "w")

def log(msg):
    print(msg)
    f.write(msg + "\n")
    f.flush()

log("=" * 80)
log("MIRROR SYMMETRY DEEP DIVE")
log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 80)

# ============================================================
# TEST 1: Uniform Z dephasing (baseline - should be 100%)
# ============================================================
log("\n### TEST 1: Uniform Z dephasing (baseline)")
log("All qubits same gamma, Z noise, Heisenberg coupling")

for n_q in [3, 4, 5]:
    gamma = 0.05
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = [np.sqrt(gamma)*op_at(Z,k,n_q) for k in range(n_q)]
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    center = n_q * gamma
    sym = check_mirror_symmetry(rates, center)
    log(f"  N={n_q}: center={center:.3f}, symmetry={sym['score']:.1%}, "
        f"matched={sym['matched']}/{sym['below']}")

# ============================================================
# TEST 2: Non-uniform dephasing (different gamma per qubit)
# ============================================================
log("\n### TEST 2: Non-uniform Z dephasing")
log("Each qubit has different gamma. Does mirror symmetry survive?")

for n_q in [3, 4, 5]:
    gammas = [0.03, 0.05, 0.07, 0.10, 0.15][:n_q]
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = [np.sqrt(gammas[k])*op_at(Z,k,n_q) for k in range(n_q)]
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    
    # What should the center be now?
    sum_gamma = sum(gammas)
    center_guess = sum_gamma  # N * avg(gamma)?
    avg_gamma = sum_gamma / n_q
    center_Navg = n_q * avg_gamma
    
    # Test different center hypotheses
    for label, c in [("sum(gamma)", sum_gamma), ("N*avg(gamma)", center_Navg),
                      ("sum(gamma)/2", sum_gamma/2)]:
        sym = check_mirror_symmetry(rates, c)
        status = "EXACT" if sym['score'] > 0.99 else f"{sym['score']:.1%}"
        log(f"  N={n_q}, gammas={gammas[:n_q]}: center={label}={c:.4f} -> {status}")
    
    log(f"    Rate range: [{min(rates):.4f}, {max(rates):.4f}]")
    log(f"    Midpoint of range: {(min(rates)+max(rates))/2:.4f}")
    mid = (min(rates)+max(rates))/2
    sym_mid = check_mirror_symmetry(rates, mid)
    log(f"    Symmetry around midpoint: {sym_mid['score']:.1%}")

# ============================================================
# TEST 3: X dephasing instead of Z
# ============================================================
log("\n### TEST 3: X dephasing (instead of Z)")
log("Replace Z jump operators with X. Does symmetry hold?")

for n_q in [3, 4, 5]:
    gamma = 0.05
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = [np.sqrt(gamma)*op_at(X,k,n_q) for k in range(n_q)]
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    center = n_q * gamma
    sym = check_mirror_symmetry(rates, center)
    log(f"  N={n_q}: center={center:.3f}, symmetry={sym['score']:.1%}, "
        f"rates=[{min(rates):.4f}, {max(rates):.4f}]")

# ============================================================
# TEST 4: Y dephasing
# ============================================================
log("\n### TEST 4: Y dephasing")

for n_q in [3, 4, 5]:
    gamma = 0.05
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = [np.sqrt(gamma)*op_at(Y,k,n_q) for k in range(n_q)]
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    center = n_q * gamma
    sym = check_mirror_symmetry(rates, center)
    log(f"  N={n_q}: center={center:.3f}, symmetry={sym['score']:.1%}, "
        f"rates=[{min(rates):.4f}, {max(rates):.4f}]")

# ============================================================
# TEST 5: Mixed noise (Z on some qubits, X on others)
# ============================================================
log("\n### TEST 5: Mixed noise (Z on qubit 0, X on others)")
log("This breaks the noise uniformity. Does symmetry survive?")

for n_q in [3, 4, 5]:
    gamma = 0.05
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = [np.sqrt(gamma)*op_at(Z,0,n_q)]  # Z on center
    for k in range(1, n_q):
        jumps.append(np.sqrt(gamma)*op_at(X,k,n_q))  # X on leaves
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    center = n_q * gamma
    sym = check_mirror_symmetry(rates, center)
    mid = (min(rates)+max(rates))/2
    sym_mid = check_mirror_symmetry(rates, mid)
    log(f"  N={n_q}: center Ng={sym['score']:.1%}, center midpoint={sym_mid['score']:.1%}")
    log(f"    rates=[{min(rates):.4f}, {max(rates):.4f}], mid={mid:.4f}")

# ============================================================
# TEST 6: XXZ anisotropy (breaks Heisenberg symmetry)
# ============================================================
log("\n### TEST 6: XXZ coupling (J_xy=1, J_zz=delta*J)")
log("Does symmetry depend on isotropic Heisenberg?")

for delta in [0.0, 0.25, 0.5, 1.0, 2.0]:
    n_q = 4
    gamma = 0.05
    bonds = xxz_bonds(n_q, [1.0]*(n_q-1), delta=delta)
    jumps = [np.sqrt(gamma)*op_at(Z,k,n_q) for k in range(n_q)]
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    if rates:
        center = n_q * gamma
        sym = check_mirror_symmetry(rates, center)
        mid = (min(rates)+max(rates))/2
        sym_mid = check_mirror_symmetry(rates, mid)
        log(f"  delta={delta:.2f}: Ng symmetry={sym['score']:.1%}, "
            f"midpoint symmetry={sym_mid['score']:.1%}, "
            f"center_mid={mid:.4f}")
    else:
        log(f"  delta={delta:.2f}: no oscillatory rates")

# ============================================================
# TEST 7: XY coupling (no ZZ term)
# ============================================================
log("\n### TEST 7: XY coupling only (no ZZ)")

for n_q in [3, 4, 5]:
    gamma = 0.05
    bonds = xy_bonds(n_q, [1.0]*(n_q-1))
    jumps = [np.sqrt(gamma)*op_at(Z,k,n_q) for k in range(n_q)]
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    if rates:
        center = n_q * gamma
        sym = check_mirror_symmetry(rates, center)
        log(f"  N={n_q}: symmetry={sym['score']:.1%}, "
            f"rates=[{min(rates):.4f}, {max(rates):.4f}]")
    else:
        log(f"  N={n_q}: no oscillatory rates")

# ============================================================
# TEST 8: Pauli weight analysis of eigenvectors
# ============================================================
log("\n### TEST 8: Pauli weight of mirror-paired eigenmodes")
log("What distinguishes a mode from its mirror partner?")
log("For Z dephasing, a Pauli string with k non-commuting ops decays at 2kg.")
log("Mirror pairs should have complementary Pauli weights.")

n_q = 3
gamma = 0.05
bonds = heisenberg_bonds(n_q, [1.0, 2.0])
jumps = [np.sqrt(gamma)*op_at(Z,k,n_q) for k in range(n_q)]
L = build_liouvillian_general(n_q, bonds, jumps)
evals, evecs = np.linalg.eig(L)

# Build Pauli basis for 3 qubits
paulis_1q = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
pauli_basis = {}
for a in 'IXYZ':
    for b in 'IXYZ':
        for c in 'IXYZ':
            label = a+b+c
            op = np.kron(np.kron(paulis_1q[a], paulis_1q[b]), paulis_1q[c])
            pauli_basis[label] = op

# For each oscillatory eigenmode, decompose eigenvector into Pauli basis
log(f"\n  Oscillatory modes with Pauli decomposition:")
log(f"  {'rate/g':>8} | {'freq':>8} | {'top Pauli components':>40}")
log("-" * 65)

center = n_q * gamma
osc_modes = []
for i in range(len(evals)):
    freq = abs(np.imag(evals[i]))
    rate = -np.real(evals[i])
    if freq > 0.1 and rate > 0.001:
        # Reshape eigenvector to density matrix form
        v = evecs[:, i]
        rho_v = v.reshape(2**n_q, 2**n_q)
        
        # Project onto Pauli basis
        pauli_weights = {}
        for label, P in pauli_basis.items():
            weight = abs(np.trace(rho_v @ P)) / (2**n_q)
            if weight > 0.01:
                pauli_weights[label] = weight
        
        # Count X/Y content (non-commuting with Z)
        xy_count = 0
        total_weight = 0
        for label, w in pauli_weights.items():
            n_xy = sum(1 for c in label if c in 'XY')
            xy_count += n_xy * w
            total_weight += w
        avg_xy = xy_count / total_weight if total_weight > 0 else 0
        
        top = sorted(pauli_weights.items(), key=lambda x: -x[1])[:3]
        top_str = ", ".join(f"{l}:{w:.3f}" for l,w in top)
        
        osc_modes.append((rate/gamma, freq/(2*np.pi), avg_xy, top_str))

# Sort and print
for rate_g, freq, avg_xy, top_str in sorted(set((round(r,3), round(fr,3), round(a,3), t) 
                                                   for r,fr,a,t in osc_modes)):
    mirror_rate = 2*n_q - rate_g
    log(f"  {rate_g:>8.3f} | {freq:>8.4f} | avg_XY={avg_xy:.3f} | {top_str}")

# ============================================================
# TEST 9: The conjugation operator
# ============================================================
log("\n### TEST 9: Searching for the conjugation operator")
log("If L has eigenvalue lambda, the mirror symmetry means")
log("there exists an operator C such that C*L*C^-1 has eigenvalue -conj(lambda) + 2Ng")
log("")
log("For Z dephasing, the candidate is the 'spin flip' superoperator")
log("that maps X->X, Y->Y, Z->-Z on each qubit (or equivalently X on each qubit).")

n_q = 3
gamma = 0.05
bonds = heisenberg_bonds(n_q, [1.0, 2.0])
jumps = [np.sqrt(gamma)*op_at(Z,k,n_q) for k in range(n_q)]
L = build_liouvillian_general(n_q, bonds, jumps)

# Candidate: X^{\otimes n} conjugation in superoperator space
# This maps |rho>> to |X^n rho X^n>>
d = 2**n_q
Xn = np.eye(d, dtype=complex)
for k in range(n_q):
    Xn = Xn @ op_at(X, k, n_q)

# Superoperator version: S_X(rho) = Xn rho Xn^dag
S_X = np.kron(Xn, Xn.conj())  # since Xn is hermitian, Xn^dag = Xn

# Check: does S_X commute with the Hamiltonian part?
H = np.zeros((d,d), dtype=complex)
for i, j, J, paulis in bonds:
    for p in paulis:
        H += J * op_at(p,i,n_q) @ op_at(p,j,n_q)

I_d = np.eye(d, dtype=complex)
L_H = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
comm_H = np.linalg.norm(S_X @ L_H - L_H @ S_X)
log(f"\n  [S_X, L_Hamiltonian] = {comm_H:.2e}  {'COMMUTES' if comm_H < 1e-8 else 'BROKEN'}")

# Check: how does S_X interact with the dissipator?
L_D = np.zeros((d*d, d*d), dtype=complex)
for Lk in jumps:
    Lk_dag = Lk.conj().T
    LdL = Lk_dag @ Lk
    L_D += np.kron(Lk, Lk.conj()) - 0.5*(np.kron(LdL, I_d) + np.kron(I_d, LdL.T))

comm_D = np.linalg.norm(S_X @ L_D - L_D @ S_X)
anti_D = np.linalg.norm(S_X @ L_D + L_D @ S_X)
log(f"  [S_X, L_Dissipator] = {comm_D:.2e}  {'COMMUTES' if comm_D < 1e-8 else 'BROKEN'}")
log(f"  {{S_X, L_Dissipator}} = {anti_D:.2e}  {'ANTI-COMMUTES' if anti_D < 1e-8 else 'NO'}")

# The full relation: S_X L S_X = ?
SLS = S_X @ L @ S_X  # since S_X^2 = I (involution)
# If S_X L S_X = -L + 2*N*gamma*I, then eigenvalues map as lambda -> -lambda + 2Ng
# Which means Re(lambda) -> -Re(lambda) + 2Ng and Im(lambda) -> -Im(lambda)
# Decay rate d = -Re(lambda) -> -(-d + 2Ng) = d - 2Ng... no
# Let me check numerically what the relationship is

# Compare eigenvalues of L and S_X L S_X
evals_L = np.sort(np.linalg.eigvals(L))
evals_SLS = np.sort(np.linalg.eigvals(SLS))

log(f"\n  Comparing eigenvalues of L vs S_X*L*S_X:")
log(f"  Are they equal? (S_X commutes): {np.allclose(evals_L, evals_SLS, atol=1e-8)}")

# Check: L + S_X*L*S_X = ?
L_plus = L + SLS
L_minus = L - SLS
log(f"  ||L + S_X*L*S_X|| / ||L|| = {np.linalg.norm(L_plus)/np.linalg.norm(L):.6f}")
log(f"  ||L - S_X*L*S_X|| / ||L|| = {np.linalg.norm(L_minus)/np.linalg.norm(L):.6f}")

# Check if L + SLS is proportional to identity or to L_D
# If L_H commutes with S_X and L_D anti-commutes, then SLS = L_H - L_D
# So L + SLS = 2*L_H, and the symmetry center comes from L_D
log(f"\n  Check decomposition:")
log(f"  ||SLS - (L_H - L_D)|| = {np.linalg.norm(SLS - (L_H - L_D)):.2e}")
log(f"  ||SLS - L_H + L_D|| = {np.linalg.norm(SLS - L_H + L_D):.2e}")

# If SLS = L_H - L_D, then:
# L = L_H + L_D
# SLS = L_H - L_D
# eigenvalues of L: lambda = lambda_H + lambda_D
# eigenvalues of SLS: lambda' = lambda_H - lambda_D
# For the mirror: we need Re(lambda) + Re(lambda') = 2*something
# Re(lambda) + Re(lambda') = 2*Re(lambda_H)
# But L_H is anti-hermitian (unitary dynamics), so Re(lambda_H) = 0
# Therefore: Re(lambda) + Re(lambda') = 0
# But that means decay(lambda) + decay(lambda') = 0, which is wrong...
# Unless the constant shift from Tr(L_D) matters

# Let's check the actual trace
tr_L = np.trace(L)
tr_LD = np.trace(L_D)
tr_LH = np.trace(L_H)
log(f"\n  Tr(L) = {tr_L:.4f}")
log(f"  Tr(L_H) = {tr_LH:.4f}")
log(f"  Tr(L_D) = {tr_LD:.4f}")
log(f"  Tr(L_D) / d^2 = {tr_LD / d**2:.6f}")
log(f"  N * gamma = {n_q * gamma:.6f}")
log(f"  2 * N * gamma = {2 * n_q * gamma:.6f}")

# ============================================================
# TEST 10: Amplitude damping (not dephasing)
# ============================================================
log("\n### TEST 10: Amplitude damping (T1 decay)")
log("Jump operator = sqrt(gamma) * sigma_minus (not Z)")
log("This is physically different from dephasing. Does symmetry hold?")

sigma_minus = np.array([[0,0],[1,0]], dtype=complex)  # |0><1|

for n_q in [3, 4]:
    gamma = 0.05
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = [np.sqrt(gamma)*op_at(sigma_minus,k,n_q) for k in range(n_q)]
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    if rates:
        center = n_q * gamma
        sym = check_mirror_symmetry(rates, center)
        mid = (min(rates)+max(rates))/2
        sym_mid = check_mirror_symmetry(rates, mid)
        log(f"  N={n_q}: Ng symmetry={sym['score']:.1%}, "
            f"midpoint symmetry={sym_mid['score']:.1%}")
        log(f"    rates=[{min(rates):.4f}, {max(rates):.4f}], mid={mid:.4f}")
    else:
        log(f"  N={n_q}: no oscillatory rates")

# ============================================================
# TEST 11: Depolarizing noise (X+Y+Z on each qubit)
# ============================================================
log("\n### TEST 11: Depolarizing noise (X+Y+Z on each qubit)")

for n_q in [3, 4]:
    gamma = 0.05
    bonds = heisenberg_bonds(n_q, [1.0]*(n_q-1))
    jumps = []
    for k in range(n_q):
        for P in [X, Y, Z]:
            jumps.append(np.sqrt(gamma/3)*op_at(P,k,n_q))
    L = build_liouvillian_general(n_q, bonds, jumps)
    rates = get_osc_rates_raw(L)
    if rates:
        center = n_q * gamma  # for depolarizing: each qubit contributes gamma
        sym = check_mirror_symmetry(rates, center)
        mid = (min(rates)+max(rates))/2
        sym_mid = check_mirror_symmetry(rates, mid)
        log(f"  N={n_q}: Ng symmetry={sym['score']:.1%}, "
            f"midpoint symmetry={sym_mid['score']:.1%}")
        log(f"    rates=[{min(rates):.4f}, {max(rates):.4f}], mid={mid:.4f}")
        # For depolarizing, the effective rate per Pauli string should be different
        # Each non-identity Pauli on qubit k contributes 2*gamma*(2/3) = 4gamma/3
        # Hmm, let's just check numerically
    else:
        log(f"  N={n_q}: no oscillatory rates")

# ============================================================
# SUMMARY
# ============================================================
log("\n" + "=" * 80)
log("SUMMARY")
log("=" * 80)
log("""
Questions answered:

1. Non-uniform dephasing: Does mirror symmetry survive?
   -> Check results above. If center shifts, what determines it?

2. X/Y dephasing: Same symmetry as Z?
   -> For isotropic Heisenberg, X/Y/Z noise should all give symmetry
      (since H commutes with global rotations)

3. Mixed noise types: Z on some, X on others?
   -> This explicitly breaks the uniform structure

4. XXZ anisotropy: Does coupling type matter?
   -> If symmetry breaks with delta != 1, it's Heisenberg-specific

5. XY coupling: No ZZ term?
   -> Tests whether ZZ in H is essential

6. Amplitude damping: Different physics entirely?
   -> Dephasing preserves energy, amplitude damping doesn't

7. Depolarizing: All three noise channels?
   -> Maximum symmetry of noise

8. Conjugation operator: What maps rate d -> 2Ng - d?
   -> S_X = X^n superoperator is the candidate.
      If [S_X, L_H] = 0 and {S_X, L_D} = 0, then
      S_X L S_X = L_H - L_D, and the symmetry follows.

KEY INSIGHT: The mirror symmetry likely comes from:
  - Hamiltonian part L_H COMMUTES with X^n (Heisenberg is isotropic)
  - Dissipator part L_D ANTI-COMMUTES with X^n (Z dephasing flips under X)
  - This means S_X*L*S_X = L_H - L_D
  - Combined with L = L_H + L_D, the eigenvalue spectrum has the
    reflection symmetry around the mean dissipation rate.
""")

log(f"Completed: {time.strftime('%Y-%m-%d %H:%M:%S')}")
log("=" * 80)
f.close()
print(f"\n>>> Results saved to: {outfile}")
