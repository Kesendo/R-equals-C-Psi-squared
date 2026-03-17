"""
NON-HEISENBERG PALINDROME TEST v2 - Fixed palindrome checker
=============================================================
Bug in v1: tolerance and center matching was off.
Using the same approach as docs_verify.py which got 100% for Heisenberg.
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)

def tensor(*ops):
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def site_op(op, s, N):
    ops = [I2]*N; ops[s] = op
    return tensor(*ops)

def build_H(N, pairs, comps):
    d = 2**N; H = np.zeros((d,d), dtype=complex)
    pm = {'X': sx, 'Y': sy, 'Z': sz}
    for i,j in pairs:
        for c, J in comps.items():
            if J == 0: continue
            H += J * site_op(pm[c[0]],i,N) @ site_op(pm[c[1]],j,N)
    return H

def build_L(H, gammas, N):
    d = 2**N; d2 = d*d; Id = np.eye(d)
    L = -1j*(np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gammas[k]*(np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L

def check_palindrome(L, sum_gamma, tol=1e-4):
    """Match approach from docs_verify.py: pair every nonzero eigenvalue."""
    evals = np.linalg.eigvals(L)
    nonzero = evals[np.abs(evals) > 1e-10]
    reals = np.real(nonzero)
    target_sum = -2*sum_gamma

    # For each rate, find partner: Re(lambda) + Re(lambda') = -2*sum_gamma
    used = set()
    paired = 0
    at_center = 0
    unpaired_rates = []

    for i in range(len(reals)):
        if i in used: continue
        partner_target = target_sum - reals[i]

        # Check if self-paired (at center)
        if abs(reals[i] - target_sum/2) < tol:
            # Close to center, but still try to find distinct partner first
            pass

        found = False
        for j in range(len(reals)):
            if j in used or j == i: continue
            if abs(reals[j] - partner_target) < tol:
                used.add(i); used.add(j)
                paired += 1; found = True; break

        if not found:
            if abs(reals[i] - target_sum) < tol:
                # At the maximum decay rate (-2*sum_gamma), self-paired with steady state
                used.add(i); at_center += 1
            else:
                unpaired_rates.append(reals[i])

    total = len(reals)
    paired_count = paired * 2 + at_center
    score = paired_count / total if total > 0 else 1.0
    
    return {
        'score': score, 'paired': paired, 'at_max': at_center,
        'unpaired': len(unpaired_rates), 'total': total,
        'min_re': min(reals), 'max_re': max(reals),
        'unpaired_rates': unpaired_rates[:5]  # first 5 for debug
    }

if __name__ == "__main__":
    gamma = 0.05

    models = {
        'Heisenberg (XX+YY+ZZ)': {'XX':1, 'YY':1, 'ZZ':1},
        'XY-only (XX+YY)':       {'XX':1, 'YY':1, 'ZZ':0},
        'Ising (ZZ only)':       {'XX':0, 'YY':0, 'ZZ':1},
        'XX-only':               {'XX':1, 'YY':0, 'ZZ':0},
        'YY-only':               {'XX':0, 'YY':1, 'ZZ':0},
        'XXZ delta=0.5':         {'XX':1, 'YY':1, 'ZZ':0.5},
        'XXZ delta=2.0':         {'XX':1, 'YY':1, 'ZZ':2.0},
        'DM (XY-YX)':           {'XY':1, 'YX':-1},
        'Heisenberg+DM':        {'XX':1, 'YY':1, 'ZZ':1, 'XY':0.3, 'YX':-0.3},
    }

    for N in [3, 4]:
        pairs = [(i,i+1) for i in range(N-1)]
        gammas = [gamma]*N
        sg = sum(gammas)

        print(f"\n{'='*90}")
        print(f"N={N}, chain, gamma={gamma}, sum_gamma={sg}, palindrome center={-sg}")
        print(f"{'='*90}")
        print(f"{'Model':<28} {'Score':>6} {'Pairs':>6} {'AtMax':>6} {'Unp':>5} "
              f"{'Total':>6} {'Min/g':>8} {'Max/g':>8}")
        print("-"*90)

        for name, comps in models.items():
            H = build_H(N, pairs, comps)
            L = build_L(H, gammas, N)
            r = check_palindrome(L, sg)
            print(f"{name:<28} {r['score']:>6.1%} {r['paired']:>6} {r['at_max']:>6} "
                  f"{r['unpaired']:>5} {r['total']:>6} {r['min_re']/gamma:>8.3f} "
                  f"{r['max_re']/gamma:>8.3f}")
            if r['unpaired'] > 0:
                print(f"  -> unpaired rates: {[f'{x/gamma:.4f}g' for x in r['unpaired_rates']]}")

    # BONUS: Does dephasing BASIS matter?
    print(f"\n{'='*90}")
    print("DEPHASING BASIS TEST (N=3, Heisenberg chain)")
    print(f"{'='*90}")
    print(f"{'Dephasing':<28} {'Score':>6} {'Pairs':>6} {'AtMax':>6} {'Unp':>5}")
    print("-"*60)

    N = 3; pairs = [(0,1),(1,2)]; gammas = [gamma]*N; sg = sum(gammas)
    comps = {'XX':1, 'YY':1, 'ZZ':1}
    H = build_H(N, pairs, comps)
    d = 2**N; d2 = d*d; Id = np.eye(d)

    for deph_name, deph_op in [('Z-dephasing', sz), ('X-dephasing', sx), ('Y-dephasing', sy)]:
        L = -1j*(np.kron(H,Id) - np.kron(Id,H.T))
        for k in range(N):
            Dk = site_op(deph_op, k, N)
            L += gamma*(np.kron(Dk, Dk.conj()) - np.eye(d2))
        r = check_palindrome(L, sg)
        print(f"{deph_name:<28} {r['score']:>6.1%} {r['paired']:>6} {r['at_max']:>6} {r['unpaired']:>5}")

    # BONUS 2: Random Hamiltonian with Z-dephasing
    print(f"\n{'='*90}")
    print("RANDOM HAMILTONIAN + Z-DEPHASING (N=3)")
    print(f"{'='*90}")
    np.random.seed(42)
    for trial in range(5):
        # Random 2-body Hamiltonian: sum of random Pauli pairs
        H_rand = np.zeros((8,8), dtype=complex)
        for i,j in [(0,1),(1,2)]:
            for p1 in [sx,sy,sz]:
                for p2 in [sx,sy,sz]:
                    coeff = np.random.randn()
                    H_rand += coeff * site_op(p1,i,3) @ site_op(p2,j,3)
        H_rand = (H_rand + H_rand.conj().T) / 2  # ensure Hermitian
        L = build_L(H_rand, gammas, 3)
        r = check_palindrome(L, sg)
        print(f"  Random trial {trial+1}: {r['score']:.1%} paired, {r['unpaired']} unpaired")

    # BONUS 3: Which cross-terms break the palindrome?
    print(f"\n{'='*90}")
    print("CROSS-TERM ANALYSIS: Which Pauli pairs preserve palindrome? (N=3, Z-deph)")
    print(f"{'='*90}")
    N = 3; pairs = [(0,1),(1,2)]; gammas = [gamma]*N; sg = sum(gammas)

    single_terms = {
        'XX': {'XX':1}, 'YY': {'YY':1}, 'ZZ': {'ZZ':1},
        'XY': {'XY':1}, 'YX': {'YX':1},
        'XZ': {'XZ':1}, 'ZX': {'ZX':1},
        'YZ': {'YZ':1}, 'ZY': {'ZY':1},
    }

    print(f"{'Term':<12} {'Score':>6} {'Unpaired':>8}")
    print("-"*30)
    for name, comps in single_terms.items():
        H = build_H(N, pairs, comps)
        L = build_L(H, gammas, N)
        r = check_palindrome(L, sg)
        mark = "OK" if r['score'] == 1.0 else "BROKEN"
        print(f"{name:<12} {r['score']:>6.1%} {r['unpaired']:>8}  {mark}")

    # BONUS 4: Pairwise combinations
    print(f"\n{'='*90}")
    print("PAIR COMBINATIONS (N=3, Z-deph)")
    print(f"{'='*90}")
    import itertools
    pauli_names = ['XX','YY','ZZ','XY','YX','XZ','ZX','YZ','ZY']
    
    works = []
    breaks = []
    for pair in itertools.combinations(pauli_names, 2):
        comps = {pair[0]: 1, pair[1]: 1}
        H = build_H(N, pairs, comps)
        L = build_L(H, gammas, N)
        r = check_palindrome(L, sg)
        if r['score'] == 1.0:
            works.append(pair)
        else:
            breaks.append((pair, r['score']))
    
    print(f"\nPalinDROMIC ({len(works)} pairs): {', '.join(['+'.join(p) for p in works])}")
    print(f"\nBROKEN ({len(breaks)} pairs):")
    for pair, score in breaks:
        print(f"  {'+'.join(pair)}: {score:.1%}")

    # BONUS 5: Tolerance sensitivity for the "broken" pairs
    print(f"\n{'='*90}")
    print("TOLERANCE SENSITIVITY: XX+ZZ at different tolerances")
    print(f"{'='*90}")
    comps_xxzz = {'XX':1, 'ZZ':1}
    H = build_H(3, [(0,1),(1,2)], comps_xxzz)
    L = build_L(H, [gamma]*3, 3)
    for tol in [1e-6, 1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2]:
        r = check_palindrome(L, 3*gamma, tol=tol)
        print(f"  tol={tol:.0e}: {r['score']:.1%} ({r['unpaired']} unpaired)")

    # BONUS 6: Random Hamiltonians with tight tolerance
    print(f"\n{'='*90}")
    print("RANDOM HAMILTONIANS RETESTED (tol=1e-6 and tol=1e-3)")
    print(f"{'='*90}")
    np.random.seed(42)
    for trial in range(5):
        H_rand = np.zeros((8,8), dtype=complex)
        for i,j in [(0,1),(1,2)]:
            for p1 in [sx,sy,sz]:
                for p2 in [sx,sy,sz]:
                    H_rand += np.random.randn() * site_op(p1,i,3) @ site_op(p2,j,3)
        H_rand = (H_rand + H_rand.conj().T) / 2
        L = build_L(H_rand, [gamma]*3, 3)
        r6 = check_palindrome(L, 3*gamma, tol=1e-6)
        r3 = check_palindrome(L, 3*gamma, tol=1e-3)
        print(f"  Trial {trial+1}: tol=1e-6 -> {r6['score']:.1%} ({r6['unpaired']} unp)  "
              f"tol=1e-3 -> {r3['score']:.1%} ({r3['unpaired']} unp)")

    # BONUS 7: Systematic - sum of ALL 9 terms with uniform coefficient
    print(f"\n{'='*90}")
    print("SYSTEMATIC: All 9 Pauli pairs with various coefficient patterns (N=3, Z-deph)")
    print(f"{'='*90}")
    N=3; pairs=[(0,1),(1,2)]; gammas=[gamma]*3; sg=3*gamma

    # All 9 terms, coefficient 1
    all9 = {'XX':1,'XY':1,'XZ':1,'YX':1,'YY':1,'YZ':1,'ZX':1,'ZY':1,'ZZ':1}
    H = build_H(N, pairs, all9)
    L = build_L(H, gammas, N)
    r = check_palindrome(L, sg, tol=1e-6)
    print(f"  All 9 terms, coeff=1:       {r['score']:.1%} ({r['unpaired']} unpaired)")

    # All 9 with different fixed coefficients per term (same across bonds)
    fixed = {'XX':0.5,'XY':1.3,'XZ':-0.7,'YX':2.1,'YY':-0.3,'YZ':0.8,'ZX':-1.1,'ZY':0.4,'ZZ':1.7}
    H = build_H(N, pairs, fixed)
    L = build_L(H, gammas, N)
    r = check_palindrome(L, sg, tol=1e-6)
    print(f"  All 9 terms, fixed coeffs:  {r['score']:.1%} ({r['unpaired']} unpaired)")

    # Now the key: DIFFERENT coefficients per BOND
    print(f"\n  Key test: different coefficients per bond")
    d = 2**N
    H_diff = np.zeros((d,d), dtype=complex)
    np.random.seed(42)
    pm = {'X': sx, 'Y': sy, 'Z': sz}
    for i,j in [(0,1),(1,2)]:
        for n1 in ['X','Y','Z']:
            for n2 in ['X','Y','Z']:
                c = np.random.randn()
                H_diff += c * site_op(pm[n1],i,N) @ site_op(pm[n2],j,N)
    H_diff = (H_diff + H_diff.conj().T)/2
    L = build_L(H_diff, gammas, N)
    r = check_palindrome(L, sg, tol=1e-6)
    print(f"  Random per bond (seed=42):  {r['score']:.1%} ({r['unpaired']} unpaired)")
    r2 = check_palindrome(L, sg, tol=1e-3)
    print(f"  Same, tol=1e-3:             {r2['score']:.1%} ({r2['unpaired']} unpaired)")
    r3 = check_palindrome(L, sg, tol=1e-1)
    print(f"  Same, tol=1e-1:             {r3['score']:.1%} ({r3['unpaired']} unpaired)")

    # BONUS 8: Symmetric vs asymmetric coupling tensor
    print(f"\n{'='*90}")
    print("J-MATRIX SYMMETRY TEST (N=3, Z-deph)")
    print("J_ab is the coupling tensor: H = sum J_ab * Pa_i Pb_j")
    print(f"{'='*90}")
    N=3; pairs=[(0,1),(1,2)]; gammas=[gamma]*3; sg=3*gamma

    tests = {
        'J=I (Heisenberg)':     {'XX':1,'YY':1,'ZZ':1},
        'J=symmetric':          {'XX':1,'XY':0.5,'YX':0.5,'YY':1,'ZZ':1},
        'J=antisymmetric':      {'XY':1,'YX':-1},
        'J=sym+antisym':        {'XX':1,'YY':1,'ZZ':1,'XY':0.5,'YX':-0.5},
        'J=asymmetric(XY!=YX)': {'XX':1,'YY':1,'ZZ':1,'XY':1.3,'YX':0.7},
        'J=full symmetric':     {'XX':0.5,'XY':1.3,'XZ':-0.7,
                                 'YX':1.3,'YY':-0.3,'YZ':0.8,
                                 'ZX':-0.7,'ZY':0.8,'ZZ':1.7},
        'J=full asymmetric':    {'XX':0.5,'XY':1.3,'XZ':-0.7,
                                 'YX':2.1,'YY':-0.3,'YZ':0.8,
                                 'ZX':-1.1,'ZY':0.4,'ZZ':1.7},
    }

    print(f"{'Model':<28} {'Score(1e-6)':>12} {'Score(1e-3)':>12}")
    print("-"*55)
    for name, comps in tests.items():
        H = build_H(N, pairs, comps)
        L = build_L(H, gammas, N)
        r6 = check_palindrome(L, sg, tol=1e-6)
        r3 = check_palindrome(L, sg, tol=1e-3)
        print(f"{name:<28} {r6['score']:>12.1%} {r3['score']:>12.1%}")
