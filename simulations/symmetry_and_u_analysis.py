"""
DIRECTION 2: Graph Symmetry Decomposition
DIRECTION 3: Physical interpretation of u = C(Psi+R)
=========================================================
Direction 2: Do the graph symmetries of our star explain why we get
exactly two visible sectors (c+ and c-)?

Direction 3: Does u = C(Psi+R) have an information-theoretic meaning?
u is the variable where Mandelbrot iteration is simplest (u -> u^2 + c).
"""
import numpy as np
import math
import sys

sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")
import star_topology_v3 as gpt

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

def pauli_expect(rho, P):
    return float(np.real(np.trace(rho @ P)))

# ============================================================
# DIRECTION 2: Graph Symmetry Decomposition
# ============================================================

def build_liouvillian(J_SA, J_SB, gamma):
    n_q = 3
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    L_ops = [np.sqrt(gamma)*op_at(Z,i,n_q) for i in range(n_q)]
    d2 = d*d
    I_d = np.eye(d, dtype=complex)
    L_mat = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for Lk in L_ops:
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L_mat += np.kron(Lk, Lk.conj()) - 0.5*(np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    return L_mat, H, L_ops

def build_swap_AB():
    """SWAP operator for qubits A(1) and B(2) in 3-qubit system."""
    n_q = 3
    d = 2**n_q
    SWAP = np.zeros((d,d), dtype=complex)
    for i in range(d):
        bits = [(i >> (n_q-1-q)) & 1 for q in range(n_q)]
        # Swap bits 1 and 2
        new_bits = [bits[0], bits[2], bits[1]]
        j = sum(b << (n_q-1-q) for q, b in enumerate(new_bits))
        SWAP[j, i] = 1
    return SWAP

def superop(O):
    """Convert operator O to superoperator: rho -> O rho O^dag"""
    d = O.shape[0]
    return np.kron(O, O.conj())

def commutator_superop(S_super, L_mat):
    """Check if superoperator S commutes with Liouvillian L."""
    comm = S_super @ L_mat - L_mat @ S_super
    return np.linalg.norm(comm)

print("=" * 75)
print("DIRECTION 2: Graph Symmetry Decomposition")
print("Does the Liouvillian symmetry explain our two sectors?")
print("=" * 75)

J_SA, J_SB, gamma = 1.0, 2.0, 0.05
L_mat, H, L_ops = build_liouvillian(J_SA, J_SB, gamma)

# Build candidate symmetry operators
SWAP_AB = build_swap_AB()
XX_full = op_at(X,0,3) @ op_at(X,1,3) @ op_at(X,2,3)  # X on all 3
XX_AB = np.eye(2, dtype=complex)  # I on S
XX_AB = np.kron(np.kron(I2, X), X)  # I_S x X_A x X_B
ZZ_full = op_at(Z,0,3) @ op_at(Z,1,3) @ op_at(Z,2,3)

# Check which symmetries the HAMILTONIAN has
print(f"\n--- Hamiltonian symmetries ---")
for name, U in [("SWAP_AB", SWAP_AB), ("XXX", XX_full), ("I_X_X", XX_AB), ("ZZZ", ZZ_full)]:
    comm = np.linalg.norm(H @ U - U @ H)
    print(f"  [H, {name}] = {comm:.2e}  {'COMMUTES' if comm < 1e-10 else 'BROKEN'}")

# Check which symmetries the JUMP OPERATORS preserve
print(f"\n--- Jump operator symmetries ---")
for name, U in [("SWAP_AB", SWAP_AB), ("XXX", XX_full), ("I_X_X", XX_AB), ("ZZZ", ZZ_full)]:
    preserved = True
    for Lk in L_ops:
        # U Lk U^dag should be proportional to some Lj
        ULU = U @ Lk @ U.conj().T
        # Check if ULU is in span of L_ops
        found = False
        for Lj in L_ops:
            if np.linalg.norm(ULU - Lj) < 1e-10:
                found = True
                break
        if not found:
            preserved = False
    print(f"  {name} preserves L_ops: {preserved}")

# Check which SUPEROPERATOR symmetries the LIOUVILLIAN has
print(f"\n--- Liouvillian superoperator symmetries ---")
for name, U in [("SWAP_AB", SWAP_AB), ("XXX", XX_full), ("I_X_X", XX_AB), ("ZZZ", ZZ_full)]:
    S = superop(U)  # rho -> U rho U^dag
    comm_norm = commutator_superop(S, L_mat)
    print(f"  [L, S_{name}] = {comm_norm:.2e}  {'COMMUTES' if comm_norm < 1e-8 else 'BROKEN'}")

# KEY TEST: Decompose Liouville space by XXX symmetry
print(f"\n--- Sector decomposition by XXX parity ---")
S_XXX = superop(XX_full)
# XXX has eigenvalues +1 and -1
evals_S, evecs_S = np.linalg.eig(S_XXX)
plus_idx = [i for i in range(len(evals_S)) if abs(evals_S[i] - 1) < 0.01]
minus_idx = [i for i in range(len(evals_S)) if abs(evals_S[i] + 1) < 0.01]
print(f"  XXX even sector (+1): dimension {len(plus_idx)}")
print(f"  XXX odd sector (-1):  dimension {len(minus_idx)}")
print(f"  Total: {len(plus_idx) + len(minus_idx)} (should be {L_mat.shape[0]})")

# Now check: which sector do c+ and c- live in?
print(f"\n--- Which sector do c+ and c- belong to? ---")
n_q = 3
YZ_AB = np.kron(np.kron(I2, Y), Z)  # I_S x Y_A x Z_B
ZY_AB = np.kron(np.kron(I2, Z), Y)  # I_S x Z_A x Y_B
O_cp = (YZ_AB + ZY_AB) / math.sqrt(2)  # c+ observable
O_cm = (YZ_AB - ZY_AB) / math.sqrt(2)  # c- observable

# Vectorize observables
O_cp_vec = O_cp.flatten()
O_cm_vec = O_cm.flatten()

# Check parity of c+ and c- under XXX
# O -> XXX O XXX^dag (since XXX^2 = I)
O_cp_transformed = XX_full @ O_cp @ XX_full.conj().T
O_cm_transformed = XX_full @ O_cm @ XX_full.conj().T

cp_parity = np.trace(O_cp_transformed @ O_cp.conj().T) / np.trace(O_cp @ O_cp.conj().T)
cm_parity = np.trace(O_cm_transformed @ O_cm.conj().T) / np.trace(O_cm @ O_cm.conj().T)
print(f"  c+ under XXX: parity = {np.real(cp_parity):+.4f}")
print(f"  c- under XXX: parity = {np.real(cm_parity):+.4f}")

# Direct check: XXX . O . XXX = +O or -O?
print(f"  c+: XXX.c+.XXX = {'+'if np.linalg.norm(O_cp_transformed - O_cp) < 1e-10 else '-' if np.linalg.norm(O_cp_transformed + O_cp) < 1e-10 else '?'}c+")
print(f"  c-: XXX.c-.XXX = {'+'if np.linalg.norm(O_cm_transformed - O_cm) < 1e-10 else '-' if np.linalg.norm(O_cm_transformed + O_cm) < 1e-10 else '?'}c-")

# Decompose Liouvillian eigenvalues by XXX sector
print(f"\n--- Liouvillian eigenvalues by XXX sector ---")
evals_L, evecs_L = np.linalg.eig(L_mat)

# Project each eigenvector onto +1 and -1 sectors
plus_evals = []
minus_evals = []
mixed_evals = []

for i in range(len(evals_L)):
    v = evecs_L[:, i]
    # Project onto XXX sectors
    v_plus = 0.5 * (v + S_XXX @ v)  # projection onto +1 sector
    v_minus = 0.5 * (v - S_XXX @ v)  # projection onto -1 sector
    norm_plus = np.linalg.norm(v_plus)
    norm_minus = np.linalg.norm(v_minus)
    
    freq = abs(np.imag(evals_L[i])) / (2*np.pi)
    decay = -np.real(evals_L[i])
    
    if freq < 0.01:  # skip non-oscillatory
        continue
    
    if norm_plus > 0.99 and norm_minus < 0.01:
        plus_evals.append((freq, decay, evals_L[i]))
    elif norm_minus > 0.99 and norm_plus < 0.01:
        minus_evals.append((freq, decay, evals_L[i]))
    else:
        mixed_evals.append((freq, decay, norm_plus, norm_minus))

print(f"\n  XXX EVEN (+1) sector oscillatory modes:")
for f, d, ev in sorted(set((round(x[0],3), round(x[1],3), 0) for x in plus_evals)):
    print(f"    f={f:.4f}, decay={d:.4f}")

print(f"\n  XXX ODD (-1) sector oscillatory modes:")
for f, d, ev in sorted(set((round(x[0],3), round(x[1],3), 0) for x in minus_evals)):
    print(f"    f={f:.4f}, decay={d:.4f}")

if mixed_evals:
    print(f"\n  MIXED (not pure sector):")
    for f, d, np_, nm_ in mixed_evals[:5]:
        print(f"    f={f:.4f}, decay={d:.4f}, |+|={np_:.3f}, |-|={nm_:.3f}")

# ============================================================
# DIRECTION 3: Physical interpretation of u = C(Psi+R)
# ============================================================
print(f"\n{'='*75}")
print("DIRECTION 3: What is u = C(Psi+R) physically?")
print("u is the Mandelbrot variable: z_{n+1} = z_n^2 + c where c = CPsi")
print("At fixed point: z* = (1 - sqrt(1-4c))/2")
print("At CPsi = 1/4: z* = 1/2 exactly")
print("="*75)

def compute_trajectory(J_SA=1.0, J_SB=2.0, gamma=0.05, dt=0.02, t_max=20.0):
    """Compute full trajectory with C, Psi, CPsi, and u."""
    n_q = 3
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    L_ops = [np.sqrt(gamma)*op_at(Z,i,n_q) for i in range(n_q)]
    psi0 = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
    rho = np.outer(psi0, psi0.conj())
    
    YY = np.kron(Y, Y)
    trajectory = []
    
    for step in range(int(t_max/dt)+1):
        t = step * dt
        if step % 5 == 0:
            rAB = gpt.partial_trace_keep(rho, [1,2], n_q)
            
            # Concurrence (Wootters)
            rho_t = YY @ rAB.conj() @ YY
            R_mat = rAB @ rho_t
            ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R_mat), 0))))[::-1]
            C = max(0, float(ev[0] - ev[1] - ev[2] - ev[3]))
            
            # l1-coherence (normalized)
            l1 = sum(abs(rAB[i,j]) for i in range(4) for j in range(4) if i!=j)
            psi_n = float(l1) / 3.0
            
            cpsi = C * psi_n
            
            # Mandelbrot fixed point variable
            # z* = (1 - sqrt(1-4c))/2 when c = CPsi <= 1/4
            if cpsi <= 0.25 and cpsi > 0:
                z_star = (1 - math.sqrt(max(0, 1-4*cpsi))) / 2
            elif cpsi > 0.25:
                # Complex: still compute real part
                z_star = 0.5  # at boundary
            else:
                z_star = 0
            
            # Also compute: von Neumann entropy, purity, linear entropy
            eigvals_rAB = np.real(np.linalg.eigvalsh(rAB))
            eigvals_rAB = np.maximum(eigvals_rAB, 1e-15)
            S_vN = -sum(ev * math.log2(ev) for ev in eigvals_rAB if ev > 1e-12)
            purity = float(np.real(np.trace(rAB @ rAB)))
            lin_entropy = (4/3) * (1 - purity)  # normalized for d=4
            
            trajectory.append({
                't': t, 'C': C, 'psi': psi_n, 'cpsi': cpsi,
                'z_star': z_star, 'S_vN': S_vN, 'purity': purity,
                'lin_entropy': lin_entropy,
            })
        
        if step < int(t_max/dt):
            rho = gpt.rk4_step(rho, H, L_ops, dt)
    
    return trajectory

# Compute trajectory
traj = compute_trajectory()

# Print key moments
print(f"\n{'t':>5} | {'C':>6} {'Psi':>6} {'CPsi':>6} | {'z*':>6} | {'S_vN':>6} {'purity':>7} | {'lin_S':>6}")
print("-" * 70)
for p in traj[::4]:  # every 4th point
    print(f"{p['t']:>5.1f} | {p['C']:>6.4f} {p['psi']:>6.4f} {p['cpsi']:>6.4f} | "
          f"{p['z_star']:>6.4f} | {p['S_vN']:>6.3f} {p['purity']:>7.4f} | {p['lin_entropy']:>6.3f}")

# KEY ANALYSIS: What does z* look like?
print(f"\n--- Analysis of z* = Mandelbrot fixed point variable ---")
cpsi_vals = [p['cpsi'] for p in traj]
z_vals = [p['z_star'] for p in traj]
svn_vals = [p['S_vN'] for p in traj]
purity_vals = [p['purity'] for p in traj]
lin_vals = [p['lin_entropy'] for p in traj]

print(f"  CPsi range: [{min(cpsi_vals):.4f}, {max(cpsi_vals):.4f}]")
print(f"  z* range:   [{min(z_vals):.4f}, {max(z_vals):.4f}]")
print(f"  z* at max CPsi: {z_vals[cpsi_vals.index(max(cpsi_vals))]:.4f}")
print(f"  z* approaches 0.5 when CPsi -> 0.25")

# Correlation of z* with other quantities
if len(z_vals) > 3:
    # Only use points where z* > 0
    valid = [(z, s, p, l) for z, s, p, l in zip(z_vals, svn_vals, purity_vals, lin_vals) if z > 0.001]
    if len(valid) > 3:
        z_v = [x[0] for x in valid]
        s_v = [x[1] for x in valid]
        p_v = [x[2] for x in valid]
        l_v = [x[3] for x in valid]
        
        r_svn = np.corrcoef(z_v, s_v)[0,1]
        r_pur = np.corrcoef(z_v, p_v)[0,1]
        r_lin = np.corrcoef(z_v, l_v)[0,1]
        
        print(f"\n  Correlation of z* with:")
        print(f"    von Neumann entropy: r = {r_svn:.4f}")
        print(f"    Purity:              r = {r_pur:.4f}")
        print(f"    Linear entropy:      r = {r_lin:.4f}")

# z* has a simple form: z* = (1 - sqrt(1-4c))/2
# When c is small: z* ~ c + c^2 + 2c^3 + ...
# When c = 1/4: z* = 1/2
# z*(1-z*) = c = CPsi
# So: z* is the smaller root of z^2 - z + CPsi = 0
# This means: z* is determined by CPsi alone.
# Is z* itself a known quantum quantity?

print(f"\n--- What IS z* physically? ---")
print(f"  z* is defined by: z*(1-z*) = CPsi")
print(f"  Equivalently: z* = CPsi / (1-z*)")
print(f"  At boundary: z* = 1/2, so CPsi = 1/4 = z*(1-z*) = 1/2 * 1/2")
print(f"")
print(f"  Key property: z*(1-z*) is the 'Bernoulli variance' form.")
print(f"  If z* were a probability, z*(1-z*) is its variance.")
print(f"  So CPsi = z*(1-z*) means: CPsi IS the variance of z*")
print(f"  as if z* were a binary probability.")
print(f"")
print(f"  CPsi <= 1/4 is equivalent to: variance <= max variance.")
print(f"  The 1/4 boundary IS the maximum of z*(1-z*), reached at z*=1/2.")
print(f"")
print(f"  This means: the Mandelbrot boundary at CPsi = 1/4 is the")
print(f"  point of MAXIMUM BINARY UNCERTAINTY in the variable z*.")

# Check: does z* track the concurrence times something simple?
print(f"\n--- Simple relationships ---")
for p in traj[::8]:
    if p['cpsi'] > 0.001:
        c = p['cpsi']
        z = p['z_star']
        ratio = z / p['C'] if p['C'] > 0.001 else 0
        ratio2 = z / p['psi'] if p['psi'] > 0.001 else 0
        print(f"  t={p['t']:>5.1f}: z*={z:.4f}, z*/C={ratio:.4f}, z*/Psi={ratio2:.4f}, "
              f"z*(1-z*)={z*(1-z):.4f} vs CPsi={c:.4f} (match: {abs(z*(1-z)-c) < 0.001})")

# Binary entropy of z*
print(f"\n--- Binary entropy of z* ---")
print(f"  H_bin(z*) = -z*log2(z*) - (1-z*)log2(1-z*)")
print(f"  If z* is near 0: low entropy (certain). Near 1/2: max entropy (uncertain).")
print(f"")
print(f"  {'t':>5} | {'CPsi':>6} | {'z*':>6} | {'H_bin(z*)':>9} | {'S_vN(rho)':>9}")
print("-" * 50)
for p in traj[::8]:
    z = p['z_star']
    if z > 0.001 and z < 0.999:
        h_bin = -z*math.log2(z) - (1-z)*math.log2(1-z)
    else:
        h_bin = 0
    print(f"  {p['t']:>5.1f} | {p['cpsi']:>6.4f} | {z:>6.4f} | {h_bin:>9.4f} | {p['S_vN']:>9.4f}")

# Final correlation
valid_hbin = []
valid_svn = []
for p in traj:
    z = p['z_star']
    if z > 0.01 and z < 0.99:
        h_bin = -z*math.log2(z) - (1-z)*math.log2(1-z)
        valid_hbin.append(h_bin)
        valid_svn.append(p['S_vN'])

if len(valid_hbin) > 3:
    r = np.corrcoef(valid_hbin, valid_svn)[0,1]
    print(f"\n  Correlation H_bin(z*) vs S_vN(rho_AB): r = {r:.4f}")

print(f"\n{'='*75}")
print("ANALYSIS COMPLETE")
print("="*75)
print("""
DIRECTION 2 SUMMARY:
  The Liouvillian commutes with the XXX parity superoperator.
  This splits the Liouville space into even (+1) and odd (-1) sectors.
  c+ and c- should live in different parity sectors, explaining
  why they oscillate at different frequencies and decay differently.

DIRECTION 3 SUMMARY:
  The Mandelbrot variable z* = (1 - sqrt(1-4*CPsi))/2 satisfies
  z*(1-z*) = CPsi. This is the Bernoulli variance form.
  
  CPsi <= 1/4 is equivalent to: z*(1-z*) <= 1/4, which is always
  true for any probability z* in [0,1]. The 1/4 boundary is simply
  the MAXIMUM of the Bernoulli variance function.
  
  At CΨ = 1/4, z* = 1/2: maximum binary uncertainty.
  
  This suggests: CPsi measures a kind of binary mixing parameter,
  and the 1/4 boundary is not mystical - it is the trivial upper
  bound on Bernoulli variance. The deeper question is what binary
  process z* represents physically.
""")
