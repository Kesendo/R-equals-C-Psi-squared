"""
NON-HEISENBERG DEEP ANALYSIS
==============================
Systematic investigation of which Hamiltonians preserve palindromic symmetry.

Writes results live to: simulations/results/non_heisenberg_deep.txt

Key questions:
1. Which single Pauli-pair terms are palindromic? (all 9)
2. Which 2-term combinations break? WHY?
3. Does the breaking grow with N?
4. Does dephasing axis matter?
5. What is the EXACT condition for palindrome preservation?

Uses COMPLEX eigenvalue pairing: lambda -> -lambda - 2*sum_gamma
"""
import numpy as np
import sys
import os
from datetime import datetime
from itertools import combinations

# Output file
RESULTS_DIR = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results"
OUT = os.path.join(RESULTS_DIR, "non_heisenberg_deep.txt")

f = open(OUT, "w", buffering=1)  # line-buffered for live reading

def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n")
    f.flush()

# Pauli matrices
I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)

PM = {'X': sx, 'Y': sy, 'Z': sz}
PAULI_NAMES = ['X', 'Y', 'Z']
ALL_PAIRS = [a+b for a in PAULI_NAMES for b in PAULI_NAMES]  # XX,XY,XZ,YX,...

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

def build_L(H, gammas, N, deph_op=None):
    """Build Liouvillian. deph_op=None means Z-dephasing."""
    if deph_op is None: deph_op = sz
    d = 2**N; d2 = d*d; Id = np.eye(d)
    L = -1j*(np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Dk = site_op(deph_op, k, N)
        L += gammas[k]*(np.kron(Dk, Dk.conj()) - np.eye(d2))
    return L

def palindrome_check(L, sum_gamma):
    """Check complex palindromic pairing. Returns max_error and details."""
    evals = np.linalg.eigvals(L)
    target = -2*sum_gamma
    nonzero = evals[np.abs(evals) > 1e-10]
    
    if len(nonzero) == 0:
        return {'max_err': 0, 'n_evals': 0, 'perfect': True}
    
    max_err = 0
    errors = []
    for lam in nonzero:
        partner = -lam + target  # -lambda - 2*sum_gamma
        dists = np.abs(evals - partner)
        best = np.min(dists)
        if best > max_err: max_err = best
        errors.append(best)
    
    return {
        'max_err': max_err,
        'mean_err': np.mean(errors),
        'median_err': np.median(errors),
        'n_evals': len(nonzero),
        'perfect': max_err < 1e-8,
        'errors': np.array(errors)
    }

def run_test(name, comps, N, gamma, deph_op=None, chain=True):
    pairs = [(i,i+1) for i in range(N-1)] if chain else [(0,i) for i in range(1,N)]
    gammas = [gamma]*N; sg = sum(gammas)
    H = build_H(N, pairs, comps)
    L = build_L(H, gammas, N, deph_op)
    return palindrome_check(L, sg)

# ================================================================
if __name__ == "__main__":
    gamma = 0.05
    
    log("=" * 90)
    log("NON-HEISENBERG DEEP ANALYSIS")
    log(f"Started: {datetime.now()}")
    log(f"Output: {OUT}")
    log("=" * 90)

    # ============================================================
    # SECTION 1: All 9 single terms, N=3,4,5
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 1: Single Pauli-pair terms")
    log("Each term tested alone with Z-dephasing")
    log("=" * 90)
    
    for N in [3, 4, 5]:
        log(f"\n  N={N}:")
        log(f"  {'Term':<8} {'Perfect':>8} {'MaxErr':>12} {'MeanErr':>12}")
        log(f"  {'-'*44}")
        for term in ALL_PAIRS:
            r = run_test(term, {term: 1.0}, N, gamma)
            tag = "YES" if r['perfect'] else "NO"
            log(f"  {term:<8} {tag:>8} {r['max_err']:>12.2e} {r['mean_err']:>12.2e}")

    # ============================================================
    # SECTION 2: All 36 two-term combinations, N=3
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 2: All 2-term combinations (N=3, Z-deph)")
    log("=" * 90)
    
    N = 3
    results_2term = {}
    log(f"\n  {'Combo':<12} {'Perfect':>8} {'MaxErr':>12}")
    log(f"  {'-'*36}")
    for t1, t2 in combinations(ALL_PAIRS, 2):
        name = f"{t1}+{t2}"
        r = run_test(name, {t1:1, t2:1}, N, gamma)
        results_2term[name] = r
        tag = "YES" if r['perfect'] else "NO"
        log(f"  {name:<12} {tag:>8} {r['max_err']:>12.2e}")
    
    # Summarize
    perfect = [k for k,v in results_2term.items() if v['perfect']]
    broken = [(k,v['max_err']) for k,v in results_2term.items() if not v['perfect']]
    log(f"\n  PALINDROMIC: {len(perfect)}/36")
    log(f"  BROKEN: {len(broken)}/36")
    if broken:
        log(f"\n  Broken combinations (sorted by error):")
        for name, err in sorted(broken, key=lambda x: -x[1]):
            log(f"    {name:<12} err={err:.4e}")

    # ============================================================
    # SECTION 3: Does breaking scale with N?
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 3: N-scaling of broken combinations")
    log("Testing the worst broken cases at N=3,4,5,6")
    log("=" * 90)
    
    # Pick the broken ones from section 2
    test_combos = {}
    if broken:
        # Top 5 worst broken + a few palindromic for control
        worst = sorted(broken, key=lambda x: -x[1])[:5]
        for name, _ in worst:
            parts = name.split('+')
            test_combos[name] = {parts[0]:1, parts[1]:1}
    # Always include controls
    test_combos['XX+YY (ctrl)'] = {'XX':1, 'YY':1}
    test_combos['Heisenberg'] = {'XX':1, 'YY':1, 'ZZ':1}
    
    log(f"\n  {'Combo':<18} {'N=3':>10} {'N=4':>10} {'N=5':>10} {'N=6':>10}")
    log(f"  {'-'*62}")
    for name, comps in test_combos.items():
        row = f"  {name:<18}"
        for N in [3, 4, 5, 6]:
            r = run_test(name, comps, N, gamma)
            if r['perfect']:
                row += f" {'OK':>10}"
            else:
                row += f" {r['max_err']:>10.2e}"
        log(row)

    # ============================================================
    # SECTION 4: Dephasing axis dependence
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 4: Dephasing axis dependence")
    log("Same Hamiltonian, three dephasing axes (N=3)")
    log("=" * 90)
    
    N = 3
    test_models = {
        'Heisenberg': {'XX':1,'YY':1,'ZZ':1},
        'XY-only': {'XX':1,'YY':1},
        'Ising': {'ZZ':1},
        'XX alone': {'XX':1},
        'DM': {'XY':1,'YX':-1},
        'XX+XY': {'XX':1,'XY':1},  # broken under Z-deph
        'XY+YZ': {'XY':1,'YZ':1},  # broken under Z-deph
    }
    
    log(f"\n  {'Model':<18} {'Z-deph':>10} {'X-deph':>10} {'Y-deph':>10}")
    log(f"  {'-'*52}")
    for name, comps in test_models.items():
        row = f"  {name:<18}"
        for deph in [sz, sx, sy]:
            r = run_test(name, comps, N, gamma, deph_op=deph)
            if r['perfect']:
                row += f" {'OK':>10}"
            else:
                row += f" {r['max_err']:>10.2e}"
        log(row)

    # ============================================================
    # SECTION 5: Coefficient dependence for broken combos
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 5: Does the error depend on coefficient ratio?")
    log("Testing XX+XY with varying J_XY/J_XX (N=3, Z-deph)")
    log("=" * 90)
    
    N = 3
    log(f"\n  {'J_XY/J_XX':>12} {'MaxErr':>12} {'Perfect':>8}")
    log(f"  {'-'*36}")
    for ratio in [0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]:
        r = run_test(f"ratio={ratio}", {'XX':1, 'XY':ratio}, N, gamma)
        tag = "YES" if r['perfect'] else "NO"
        log(f"  {ratio:>12.3f} {r['max_err']:>12.2e} {tag:>8}")

    # Same for XY+YZ
    log(f"\n  XY+YZ with varying ratio:")
    log(f"  {'J_YZ/J_XY':>12} {'MaxErr':>12} {'Perfect':>8}")
    log(f"  {'-'*36}")
    for ratio in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 10.0]:
        r = run_test(f"ratio={ratio}", {'XY':1, 'YZ':ratio}, N, gamma)
        tag = "YES" if r['perfect'] else "NO"
        log(f"  {ratio:>12.3f} {r['max_err']:>12.2e} {tag:>8}")

    # ============================================================
    # SECTION 6: Gamma dependence - does error scale with gamma?
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 6: Gamma dependence of the breaking")
    log("XX+XY at different gamma values (N=3)")
    log("=" * 90)
    
    N = 3
    log(f"\n  {'gamma':>12} {'MaxErr':>12} {'Err/gamma':>12}")
    log(f"  {'-'*40}")
    for g in [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]:
        r = run_test("XX+XY", {'XX':1,'XY':1}, N, g)
        ratio = r['max_err']/g if g > 0 else 0
        log(f"  {g:>12.3f} {r['max_err']:>12.2e} {ratio:>12.4f}")

    # ============================================================
    # SECTION 7: The Pi operator test
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 7: Does Pi L Pi^-1 = -L - 2Sg*I hold?")
    log("Direct check of the conjugation relation for each model")
    log("=" * 90)
    
    N = 3; d = 8; d2 = 64
    gammas = [gamma]*N; sg = sum(gammas)
    
    # Build Pi as a d2 x d2 matrix in the Pauli basis
    # Per-site: I->X, X->I, Y->iZ, Z->iY
    pi_single = {
        'I': ('X', 1),
        'X': ('I', 1),
        'Y': ('Z', 1j),
        'Z': ('Y', 1j),
    }
    
    # Build Pi in the computational basis
    # Pi acts on operators (d x d matrices), so it's a d^2 x d^2 superoperator
    # We represent operators as vectors: |rho>> = vec(rho)
    # Pi(rho) means: decompose rho in Pauli basis, transform each Pauli string, recompose
    
    def build_pi(N):
        d = 2**N
        pauli_map = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}
        labels = ['I','X','Y','Z']
        
        # All N-qubit Pauli strings
        from itertools import product as iprod
        basis_labels = list(iprod(labels, repeat=N))
        
        # Build each Pauli string as d x d matrix
        pauli_strings = []
        for bl in basis_labels:
            mat = pauli_map[bl[0]]
            for i in range(1, N):
                mat = np.kron(mat, pauli_map[bl[i]])
            pauli_strings.append(mat)
        
        # Pi matrix in Pauli basis: Pi maps |P_i>> to phase * |P_j>>
        # Then convert to computational vectorized basis
        Pi = np.zeros((d*d, d*d), dtype=complex)
        
        for idx, bl in enumerate(basis_labels):
            # Apply per-site transformation
            new_labels = []
            phase = 1.0
            for site_label in bl:
                new_label, site_phase = pi_single[site_label]
                new_labels.append(new_label)
                phase *= site_phase
            
            # Find index of new Pauli string
            new_idx = basis_labels.index(tuple(new_labels))
            
            # In vectorized basis: |P_i>> maps to phase * |P_j>>
            # |P>> = vec(P) = P.flatten() (column-major implied by kron structure)
            vec_old = pauli_strings[idx].flatten()
            vec_new = pauli_strings[new_idx].flatten()
            
            # Pi |vec_old> = phase * |vec_new>
            # This is a rank-1 contribution to Pi
            # But we need to be careful: Pauli strings are orthogonal under Frobenius
            # Tr(P_i^dag P_j) = d * delta_ij
            # So: Pi = sum_i phase_i |vec(P_new_i)><vec(P_old_i)| / (vec_old . vec_old)
            norm_sq = np.dot(vec_old.conj(), vec_old).real  # = d for normalized Paulis
            Pi += phase * np.outer(vec_new, vec_old.conj()) / norm_sq
        
        return Pi
    
    log("  Building Pi operator...")
    Pi = build_pi(N)
    Pi_inv = np.linalg.inv(Pi)
    
    # Test for each model
    test_models_pi = {
        'Heisenberg': {'XX':1,'YY':1,'ZZ':1},
        'XY-only': {'XX':1,'YY':1},
        'Ising': {'ZZ':1},
        'DM': {'XY':1,'YX':-1},
        'XX+XY': {'XX':1,'XY':1},
        'XY+YZ': {'XY':1,'YZ':1},
        'XX+YY(1,2)': {'XX':1,'YY':2},
        'All9 uniform': {t:1 for t in ALL_PAIRS},
        'Random': {},  # will be filled
    }
    
    # Generate a random Hermitian Hamiltonian
    np.random.seed(42)
    rand_comps = {}
    for t in ALL_PAIRS:
        rand_comps[t] = np.random.randn()
    test_models_pi['Random'] = rand_comps
    
    pairs = [(0,1),(1,2)]
    log(f"\n  {'Model':<20} {'||Pi L Pi^-1 + L + 2Sg||':>28} {'Palindrome':>12}")
    log(f"  {'-'*64}")
    
    for name, comps in test_models_pi.items():
        H = build_H(N, pairs, comps)
        L = build_L(H, gammas, N)
        
        # Check: Pi L Pi^-1 + L + 2*sg*I should be zero for palindromic
        PLP = Pi @ L @ Pi_inv
        residual = PLP + L + 2*sg*np.eye(d2)
        norm = np.linalg.norm(residual)
        
        # Also check palindrome via eigenvalues
        pr = palindrome_check(L, sg)
        tag = "YES" if pr['perfect'] else f"err={pr['max_err']:.2e}"
        
        log(f"  {name:<20} {norm:>28.6e} {tag:>12}")

    # ============================================================
    # SECTION 8: Decompose H into Pi-compatible and Pi-breaking
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 8: Analytical decomposition")
    log("Which part of the Hamiltonian breaks Pi L Pi^-1 = -L - 2Sg?")
    log("The Hamiltonian commutator: [H,.] should anti-commute with Pi")
    log("Check: Pi [H,rho] Pi^-1 = -[H, Pi rho Pi^-1] ?")
    log("=" * 90)
    
    # For each single-bond Pauli term, check if it anti-commutes with Pi
    # in the superoperator sense
    N = 3; d = 8; d2 = 64
    Id = np.eye(d)
    
    log(f"\n  {'Term':<8} {'||Pi[H,.]+[H,.]Pi||':>24} {'Anti-commutes':>16}")
    log(f"  {'-'*52}")
    
    for term in ALL_PAIRS:
        # Build H with just this one term on bond (0,1)
        H_term = site_op(PM[term[0]], 0, N) @ site_op(PM[term[1]], 1, N)
        # Commutator superoperator: [H_term, .] = H kron I - I kron H^T
        comm = np.kron(H_term, Id) - np.kron(Id, H_term.T)
        
        # Check: Pi @ comm + comm @ Pi should be zero for anti-commutation
        anticomm = Pi @ comm @ Pi_inv + comm  # Pi [H,.] Pi^-1 + [H,.] 
        norm = np.linalg.norm(anticomm)
        ok = "YES" if norm < 1e-10 else f"NO ({norm:.4e})"
        log(f"  {term:<8} {norm:>24.6e} {ok:>16}")
    
    # Same for bond (1,2) to check topology independence
    log(f"\n  Same analysis for bond (1,2):")
    log(f"  {'Term':<8} {'||Pi[H,.]+[H,.]Pi||':>24} {'Anti-commutes':>16}")
    log(f"  {'-'*52}")
    
    for term in ALL_PAIRS:
        H_term = site_op(PM[term[0]], 1, N) @ site_op(PM[term[1]], 2, N)
        comm = np.kron(H_term, Id) - np.kron(Id, H_term.T)
        anticomm = Pi @ comm @ Pi_inv + comm
        norm = np.linalg.norm(anticomm)
        ok = "YES" if norm < 1e-10 else f"NO ({norm:.4e})"
        log(f"  {term:<8} {norm:>24.6e} {ok:>16}")

    # ============================================================
    # SECTION 9: The real question - is there a DIFFERENT Pi 
    # for each Hamiltonian?
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 9: Could there be a DIFFERENT Pi for broken models?")
    log("If the standard Pi doesn't work for XX+XY, maybe a modified one does")
    log("The eigenvalue check showed palindrome even when Pi doesn't anti-commute")
    log("This would mean there exists ANOTHER symmetry operator")
    log("=" * 90)
    
    # Key insight: if eigenvalues are palindromic but Pi doesn't work,
    # there must be another conjugation operator Q such that Q L Q^-1 = -L - 2Sg
    
    # Let's check: for the "broken" cases from Section 2,
    # are the eigenvalues ACTUALLY palindromic at very tight tolerance?
    
    N = 3
    test_cases = {
        'Heisenberg': {'XX':1,'YY':1,'ZZ':1},
        'XX+XY': {'XX':1,'XY':1},
        'XY+YZ': {'XY':1,'YZ':1},
        'XZ+ZY': {'XZ':1,'ZY':1},
    }
    
    log(f"\n  Ultra-precise palindrome check (tol sweep):")
    for name, comps in test_cases.items():
        pairs = [(0,1),(1,2)]
        gammas = [gamma]*N; sg = sum(gammas)
        H = build_H(N, pairs, comps)
        L = build_L(H, gammas, N)
        evals = np.linalg.eigvals(L)
        target = -2*sg
        nonzero = evals[np.abs(evals) > 1e-10]
        
        # For each eigenvalue, find its best partner
        partner_errs = []
        for lam in nonzero:
            partner = -lam + target
            dists = np.abs(evals - partner)
            partner_errs.append(np.min(dists))
        
        partner_errs = np.array(partner_errs)
        log(f"\n  {name}:")
        log(f"    Eigenvalues: {len(nonzero)}")
        log(f"    Partner error - min:    {partner_errs.min():.2e}")
        log(f"    Partner error - max:    {partner_errs.max():.2e}")
        log(f"    Partner error - median: {np.median(partner_errs):.2e}")
        log(f"    Partner error - mean:   {partner_errs.mean():.2e}")
        log(f"    Errors > 1e-8:  {np.sum(partner_errs > 1e-8)}")
        log(f"    Errors > 1e-4:  {np.sum(partner_errs > 1e-4)}")
        log(f"    Errors > 1e-2:  {np.sum(partner_errs > 1e-2)}")

    # ============================================================
    # SECTION 10: Non-uniform gamma with broken combos
    # ============================================================
    log("\n" + "=" * 90)
    log("SECTION 10: Non-uniform gamma")
    log("Does the palindrome survive non-uniform dephasing for broken combos?")
    log("=" * 90)
    
    N = 3
    gamma_sets = [
        [0.05, 0.05, 0.05],
        [0.01, 0.05, 0.10],
        [0.02, 0.08, 0.15],
    ]
    
    for comps_name, comps in [('Heisenberg', {'XX':1,'YY':1,'ZZ':1}), 
                               ('XX+XY', {'XX':1,'XY':1})]:
        log(f"\n  {comps_name}:")
        for gammas in gamma_sets:
            sg = sum(gammas)
            H = build_H(N, [(0,1),(1,2)], comps)
            L = build_L(H, gammas, N)
            r = palindrome_check(L, sg)
            tag = "YES" if r['perfect'] else f"err={r['max_err']:.2e}"
            log(f"    gamma={gammas} -> {tag}")

    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    log("\n" + "=" * 90)
    log("SUMMARY")
    log("=" * 90)
    log(f"\nCompleted: {datetime.now()}")
    log("See detailed analysis above for conclusions.")
    log("Key sections: 7 (Pi operator), 8 (anti-commutation), 9 (ultra-precise check)")
    
    f.close()
    print(f"\n>>> Results saved to: {OUT}")
