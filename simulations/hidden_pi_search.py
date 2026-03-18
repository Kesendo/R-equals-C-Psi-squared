"""
HIDDEN PI OPERATOR SEARCH
===========================
For models that are palindromic but don't anti-commute with our known Pi,
numerically construct the conjugation operator Q from eigenvector pairing.

If lambda_j = -lambda_i - 2*sum_gamma, then Q maps |r_i> to |r_j>.
We construct Q, verify Q L Q^-1 = -L - 2Sg*I, and analyze Q's structure.

Output: simulations/results/hidden_pi_search.txt
"""
import numpy as np
from datetime import datetime
import os

RESULTS_DIR = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results"
OUT = os.path.join(RESULTS_DIR, "hidden_pi_search.txt")
f = open(OUT, "w", buffering=1)

def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n")
    f.flush()

I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
PM = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}
LABELS = ['I','X','Y','Z']

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

def build_L(H, gammas, N):
    d = 2**N; d2 = d*d; Id = np.eye(d)
    L = -1j*(np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gammas[k]*(np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L

def build_known_pi(N):
    """Build our known Pi operator as a d^2 x d^2 matrix."""
    from itertools import product as iprod
    d = 2**N
    pi_map = {'I': ('X',1), 'X': ('I',1), 'Y': ('Z',1j), 'Z': ('Y',1j)}
    
    basis_labels = list(iprod(LABELS, repeat=N))
    pauli_mats = []
    for bl in basis_labels:
        mat = PM[bl[0]]
        for i in range(1,N): mat = np.kron(mat, PM[bl[i]])
        pauli_mats.append(mat)
    
    Pi = np.zeros((d*d, d*d), dtype=complex)
    for idx, bl in enumerate(basis_labels):
        new_labels = []
        phase = 1.0
        for site_label in bl:
            nl, sp = pi_map[site_label]
            new_labels.append(nl); phase *= sp
        new_idx = basis_labels.index(tuple(new_labels))
        vec_old = pauli_mats[idx].flatten()
        vec_new = pauli_mats[new_idx].flatten()
        norm_sq = np.dot(vec_old.conj(), vec_old).real
        Pi += phase * np.outer(vec_new, vec_old.conj()) / norm_sq
    return Pi

def construct_Q_from_eigenvectors(L, sum_gamma, tol=1e-6):
    """
    Construct Q such that Q L Q^-1 = -L - 2*sum_gamma*I
    by pairing eigenvectors of partner eigenvalues.
    
    If lambda is eigenvalue with right eigenvector |r>, then
    -lambda - 2*sum_gamma should also be an eigenvalue with eigenvector |r'>.
    Q maps |r> -> |r'>.
    """
    d2 = L.shape[0]
    target = -2*sum_gamma
    
    evals, R = np.linalg.eig(L)  # R[:,i] is right eigenvector for evals[i]
    
    # Pair eigenvalues: lambda_i <-> -lambda_i - 2*sum_gamma
    pairs = []
    used = set()
    unpaired = []
    
    for i in range(len(evals)):
        if i in used: continue
        partner = -evals[i] + target
        
        # Find closest match
        best_j = -1
        best_err = 999
        for j in range(len(evals)):
            if j in used: continue
            err = abs(evals[j] - partner)
            if err < best_err:
                best_err = err; best_j = j
        
        if best_err < tol and best_j != i:
            pairs.append((i, best_j))
            used.add(i); used.add(best_j)
        elif abs(evals[i] - target/2) < tol:
            # Self-paired (at center)
            pairs.append((i, i))
            used.add(i)
        else:
            unpaired.append(i)
    
    if unpaired:
        return None, f"Failed: {len(unpaired)} unpaired eigenvalues"
    
    # Construct Q: Q |r_i> = |r_j> for each pair (i,j)
    # In matrix form: Q = sum_pairs |r_j><l_i| where <l_i| is left eigenvector
    # Left eigenvectors: L_left = inv(R), rows are left eigenvectors
    try:
        L_left = np.linalg.inv(R)
    except np.linalg.LinAlgError:
        L_left = np.linalg.pinv(R)
    
    Q = np.zeros((d2, d2), dtype=complex)
    for i, j in pairs:
        Q += np.outer(R[:,j], L_left[i,:])
    
    return Q, f"OK: {len(pairs)} pairs"

def analyze_Q_structure(Q, N, name, known_Pi=None):
    """Analyze the structure of a numerically constructed Q operator."""
    d2 = Q.shape[0]
    d = int(np.sqrt(d2))
    log(f"\n  === {name} ===")
    
    # 1. Basic properties
    det = np.linalg.det(Q)
    log(f"  det(Q) = {det:.6f} (magnitude {abs(det):.6f})")
    
    # Is Q unitary? Q^dag Q = I?
    QdQ = Q.conj().T @ Q
    unitarity = np.linalg.norm(QdQ - np.eye(d2))
    log(f"  ||Q^dag Q - I|| = {unitarity:.4e} ({'unitary' if unitarity < 1e-8 else 'NOT unitary'})")
    
    # Is Q involutory? Q^2 = cI?
    Q2 = Q @ Q
    # Check if Q^2 is proportional to identity
    if abs(Q2[0,0]) > 1e-10:
        Q2_normalized = Q2 / Q2[0,0]
        involutory_err = np.linalg.norm(Q2_normalized - np.eye(d2))
        log(f"  ||Q^2/Q^2[0,0] - I|| = {involutory_err:.4e} (Q^2 = {Q2[0,0]:.4f}*I? {'YES' if involutory_err < 1e-8 else 'NO'})")
    else:
        log(f"  Q^2[0,0] ~ 0, not involutory")
    
    # 2. Compare with known Pi
    if known_Pi is not None:
        diff = np.linalg.norm(Q - known_Pi)
        log(f"  ||Q - Pi_known|| = {diff:.4e} ({'SAME' if diff < 1e-8 else 'DIFFERENT'})")
        
        # Check if Q = phase * Pi
        if abs(Q[0,0]) > 1e-10 and abs(known_Pi[0,0]) > 1e-10:
            phase = Q[0,0] / known_Pi[0,0]
            diff_phase = np.linalg.norm(Q - phase * known_Pi)
            log(f"  ||Q - ({phase:.4f})*Pi|| = {diff_phase:.4e} ({'Q = phase*Pi' if diff_phase < 1e-8 else 'NOT proportional'})")
    
    # 3. Sparsity
    nonzero = np.sum(np.abs(Q) > 1e-10)
    log(f"  Nonzero entries: {nonzero}/{d2*d2} ({100*nonzero/(d2*d2):.1f}%)")
    
    # 4. Pauli decomposition of Q
    # Q acts on d^2-dim space (superoperator). Decompose in Pauli-Pauli basis.
    from itertools import product as iprod
    basis_labels = list(iprod(LABELS, repeat=N))
    pauli_mats = []
    for bl in basis_labels:
        mat = PM[bl[0]]
        for i in range(1,N): mat = np.kron(mat, PM[bl[i]])
        pauli_mats.append(mat)
    
    # Q maps vec(P_a) -> sum_b Q_ba vec(P_b)
    # Q_ba = Tr(P_b^dag Q(P_a)) / d
    # where Q(P_a) means: reshape Q @ vec(P_a) back to d x d
    log(f"\n  Pauli-to-Pauli mapping (top entries):")
    mapping = {}
    for a_idx, a_label in enumerate(basis_labels):
        vec_a = pauli_mats[a_idx].flatten()
        Qvec = Q @ vec_a
        Qmat = Qvec.reshape(d, d)
        
        for b_idx, b_label in enumerate(basis_labels):
            coeff = np.trace(pauli_mats[b_idx].conj().T @ Qmat) / d
            if abs(coeff) > 1e-8:
                a_str = ''.join(a_label)
                b_str = ''.join(b_label)
                mapping[(a_str, b_str)] = coeff
    
    # Sort by magnitude and show top mappings
    sorted_map = sorted(mapping.items(), key=lambda x: -abs(x[1]))
    shown = 0
    for (a_str, b_str), coeff in sorted_map:
        if shown >= 20: break
        phase_str = f"{coeff:.4f}"
        if abs(abs(coeff) - 1) < 1e-8:
            angle = np.angle(coeff)
            if abs(angle) < 1e-8: phase_str = "+1"
            elif abs(angle - np.pi/2) < 1e-8: phase_str = "+i"
            elif abs(angle + np.pi/2) < 1e-8: phase_str = "-i"
            elif abs(abs(angle) - np.pi) < 1e-8: phase_str = "-1"
        log(f"    {a_str} -> {phase_str} * {b_str}")
        shown += 1
    
    if len(sorted_map) > 20:
        log(f"    ... ({len(sorted_map)} total nonzero mappings)")
    
    # 5. Check if Q is a tensor product of per-site operators
    # For N=3, if Q = Q1 x Q2 x Q3, then Q acting on vec(P1xP2xP3) 
    # should map each site independently
    log(f"\n  Per-site factorization test:")
    # Check: does Q map single-site Paulis to single-site Paulis?
    factorizable = True
    for site in range(N):
        for p_name in ['I','X','Y','Z']:
            # Build P with p_name at site, I elsewhere
            test_labels = ['I']*N
            test_labels[site] = p_name
            a_str = ''.join(test_labels)
            
            # Find what Q maps it to
            targets = [(b, c) for (a,b),c in mapping.items() if a == a_str]
            if len(targets) == 1:
                b_str, coeff = targets[0]
                # Check if target also has I at all other sites
                other_sites_identity = all(b_str[s] == 'I' for s in range(N) if s != site)
                if not other_sites_identity:
                    factorizable = False
                    log(f"    Site {site}, {p_name}: maps to {b_str} (NOT per-site)")
            elif len(targets) == 0:
                log(f"    Site {site}, {p_name}: maps to NOTHING")
                factorizable = False
            else:
                factorizable = False
                log(f"    Site {site}, {p_name}: maps to {len(targets)} targets (NOT clean)")
    
    if factorizable:
        log(f"  -> Q IS a per-site tensor product!")
        # Extract per-site action
        for site in range(N):
            log(f"    Site {site}:")
            for p_name in ['I','X','Y','Z']:
                test_labels = ['I']*N
                test_labels[site] = p_name
                a_str = ''.join(test_labels)
                targets = [(b, c) for (a,b),c in mapping.items() if a == a_str]
                if targets:
                    b_str, coeff = targets[0]
                    target_pauli = b_str[site]
                    phase_str = f"{coeff:.4f}"
                    if abs(abs(coeff)-1) < 1e-8:
                        angle = np.angle(coeff)
                        if abs(angle) < 1e-8: phase_str = "+1"
                        elif abs(angle - np.pi/2) < 1e-8: phase_str = "+i"
                        elif abs(angle + np.pi/2) < 1e-8: phase_str = "-i"
                        elif abs(abs(angle) - np.pi) < 1e-8: phase_str = "-1"
                    log(f"      {p_name} -> {phase_str} * {target_pauli}")
    else:
        log(f"  -> Q is NOT a simple per-site tensor product")
    
    return mapping

if __name__ == "__main__":
    N = 3
    gamma = 0.05
    gammas = [gamma]*N
    sg = sum(gammas)
    pairs = [(0,1),(1,2)]
    
    log("=" * 90)
    log("HIDDEN PI OPERATOR SEARCH")
    log(f"N={N}, gamma={gamma}, Z-dephasing")
    log(f"Started: {datetime.now()}")
    log("=" * 90)
    
    known_Pi = build_known_pi(N)
    
    # Verify known Pi works for Heisenberg
    log("\n" + "=" * 90)
    log("CONTROL: Known Pi on Heisenberg")
    log("=" * 90)
    H_heis = build_H(N, pairs, {'XX':1,'YY':1,'ZZ':1})
    L_heis = build_L(H_heis, gammas, N)
    
    d2 = L_heis.shape[0]
    residual = known_Pi @ L_heis @ np.linalg.inv(known_Pi) + L_heis + 2*sg*np.eye(d2)
    log(f"  ||Pi L Pi^-1 + L + 2Sg|| = {np.linalg.norm(residual):.4e}")
    
    Q_heis, status = construct_Q_from_eigenvectors(L_heis, sg)
    log(f"  Eigenvector Q construction: {status}")
    if Q_heis is not None:
        analyze_Q_structure(Q_heis, N, "Heisenberg (control)", known_Pi)
    
    # ================================================================
    # Now the interesting cases: palindromic but Pi doesn't work
    # ================================================================
    models = {
        'DM (XY-YX)':      {'XY':1, 'YX':-1},
        'XY alone':         {'XY':1},
        'YX alone':         {'YX':1},
        'XZ alone':         {'XZ':1},
        'ZX alone':         {'ZX':1},
        'YZ alone':         {'YZ':1},
        'ZY alone':         {'ZY':1},
        'XY+YX symmetric':  {'XY':1, 'YX':1},
        'All9 uniform':     {t:1 for t in ['XX','XY','XZ','YX','YY','YZ','ZX','ZY','ZZ']},
        'ZZ+XY':            {'ZZ':1, 'XY':1},
    }
    
    for name, comps in models.items():
        log(f"\n{'=' * 90}")
        log(f"MODEL: {name}")
        log(f"Comps: {comps}")
        log(f"{'=' * 90}")
        
        H = build_H(N, pairs, comps)
        L = build_L(H, gammas, N)
        
        # Check if Pi works
        residual = known_Pi @ L @ np.linalg.inv(known_Pi) + L + 2*sg*np.eye(d2)
        pi_norm = np.linalg.norm(residual)
        log(f"  Known Pi residual: {pi_norm:.4e} ({'WORKS' if pi_norm < 1e-8 else 'BROKEN'})")
        
        # Construct Q from eigenvectors
        Q, status = construct_Q_from_eigenvectors(L, sg)
        log(f"  Q construction: {status}")
        
        if Q is not None:
            # Verify Q works (use pinv for robustness with degenerate eigenvalues)
            try:
                Q_inv = np.linalg.inv(Q)
            except np.linalg.LinAlgError:
                log(f"  Q is singular, using pseudoinverse")
                Q_inv = np.linalg.pinv(Q)
            residual_Q = Q @ L @ Q_inv + L + 2*sg*np.eye(d2)
            log(f"  ||Q L Q^-1 + L + 2Sg|| = {np.linalg.norm(residual_Q):.4e}")
            
            analyze_Q_structure(Q, N, name, known_Pi)

    # ================================================================
    # COMPARISON: all per-site actions side by side
    # ================================================================
    log(f"\n{'=' * 90}")
    log("SUMMARY: Per-site Pauli maps for each model")
    log("Format: I->?, X->?, Y->?, Z->?  (per site, if factorizable)")
    log(f"{'=' * 90}")
    
    log(f"\n  Known Pi (Z-dephasing):  I->+1*X, X->+1*I, Y->+i*Z, Z->+i*Y")
    log(f"  This is the operator that works for XX, YY, ZZ terms.\n")
    
    log(f"  If the hidden Q operators are also per-site tensor products,")
    log(f"  we can read off the per-site transformation rule for each model.")
    log(f"  If they all share the same per-site rule, there's a universal Pi")
    log(f"  we haven't recognized. If they differ, the family is genuinely diverse.")
    
    # ================================================================
    # BONUS: Try to find Pi for X-dephasing and Y-dephasing
    # ================================================================
    log(f"\n{'=' * 90}")
    log("BONUS: Construct Q for Heisenberg under X-dephasing and Y-dephasing")
    log("(Should reveal the rotated Pi operators)")
    log(f"{'=' * 90}")
    
    for deph_name, deph_op in [('X-dephasing', sx), ('Y-dephasing', sy)]:
        log(f"\n  --- {deph_name} ---")
        H = build_H(N, pairs, {'XX':1,'YY':1,'ZZ':1})
        d = 2**N; Id = np.eye(d)
        L = -1j*(np.kron(H,Id) - np.kron(Id,H.T))
        for k in range(N):
            Dk = site_op(deph_op, k, N)
            L += gamma*(np.kron(Dk, Dk.conj()) - np.eye(d*d))
        
        Q, status = construct_Q_from_eigenvectors(L, sg)
        log(f"  Q construction: {status}")
        if Q is not None:
            try:
                Q_inv = np.linalg.inv(Q)
            except np.linalg.LinAlgError:
                log(f"  Q singular, using pseudoinverse")
                Q_inv = np.linalg.pinv(Q)
            res = Q @ L @ Q_inv + L + 2*sg*np.eye(d*d)
            log(f"  ||Q L Q^-1 + L + 2Sg|| = {np.linalg.norm(res):.4e}")
            analyze_Q_structure(Q, N, f"Heisenberg + {deph_name}", known_Pi)
    
    log(f"\n{'=' * 90}")
    log(f"Completed: {datetime.now()}")
    log(f"{'=' * 90}")
    f.close()
    print(f"\n>>> Results saved to: {OUT}")
