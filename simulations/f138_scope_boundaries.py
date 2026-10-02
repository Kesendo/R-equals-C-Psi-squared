"""Two boundaries of F138's sufficient direction, measured.

1. The three-axis lift needs ONE bond letter per component. Field-free, one dephasing axis per site, all three
   axes present, every bond a single term P⊗P: on K3 and P4 every assignment of bond letters and axes, read by
   F158's two ends (nullities of the commutator with H on the dark and the lit Pauli strings, ranks modulo two
   primes, the smaller nullity kept); on K3 also by the spectrum's optimal pairing about -sum(gamma) (dense
   Liouvillian, gamma = 0.05), which must agree, its two verdicts at least six decades apart.
2. The clauses need bonds built of P⊗P terms. P3 with bonds XX + XY under Z on every site, no field; and the
   Dzyaloshinskii-Moriya chain XY - YX at N = 3 under Z on every site, palindromic with no field and broken
   by X fields the clauses admit.
Run: python simulations/f138_scope_boundaries.py   (output: simulations/results/f138_scope_boundaries.txt)"""
import itertools
import numpy as np
from scipy.optimize import linear_sum_assignment

I2 = np.eye(2)
P = {'X': np.array([[0, 1], [1, 0]], complex), 'Y': np.array([[0, -1j], [1j, 0]]), 'Z': np.diag([1., -1]).astype(complex)}
L_IDX = {'X': 1, 'Y': 2, 'Z': 3}
PRIMES = [2147483629, 2147483587]

def site(op, l, n):
    r = np.eye(1)
    for k in range(n): r = np.kron(r, op if k == l else I2)
    return r

def pairing_distance(H, jumps, n, g=0.05):
    d = 2 ** n; Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for A in jumps: L += g * (np.kron(A, A.conj()) - np.kron(Id, Id))
    ev = np.linalg.eigvals(L); s = len(jumps) * g
    C = np.abs(ev[:, None] - (-ev[None, :] - 2 * s)); r, c = linear_sum_assignment(C)
    return C[r, c].max()

# --- F158's two ends by ranks modulo primes (strings commuting / anticommuting with every jump) ---
def mul1(a, b):
    if a == 0: return 0, b
    if b == 0: return 0, a
    if a == b: return 0, 0
    return (1 if (a, b) in ((1, 2), (2, 3), (3, 1)) else 3), 6 - a - b

def mul(s, t):
    ph, out = 0, []
    for a, b in zip(s, t):
        q, c = mul1(a, b); ph += q; out.append(c)
    return ph % 4, tuple(out)

def anti(s, t): return sum(1 for a, b in zip(s, t) if a and b and a != b) % 2 == 1

def rank_mod(M, p):
    A = M.copy() % p; r = 0; rows, cols = A.shape
    for c in range(cols):
        nz = np.nonzero(A[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: A[[r, i]] = A[[i, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): A[nzr] = (A[nzr] - (col[nzr, None] * A[r][None, :]) % p) % p
        r += 1
        if r == rows: break
    return r

def ends(n, terms, jumps):
    out = []
    for far in (False, True):
        cols = [s for s in itertools.product(range(4), repeat=n) if all(anti(s, j) == far for j in jumps)]
        idx, ent = {}, []
        for w, t in terms:
            for jx, s in enumerate(cols):
                if anti(t, s):
                    k, u = mul(t, s)
                    if u not in idx: idx[u] = len(idx)
                    ent.append((idx[u], jx, (1 if k == 1 else -1) * w))
        M = np.zeros((max(1, len(idx)), len(cols)), dtype=np.int64)
        for i, j, v in ent: M[i, j] += v
        out.append(min(len(cols) - rank_mod(M, p) for p in PRIMES))
    return out

def letter_terms(n, bonds, letters):
    """Single-term bonds P(x)P as weighted Pauli strings for ends()."""
    terms = []
    for (a, b), c in zip(bonds, letters):
        t = [0] * n; t[a] = L_IDX[c]; t[b] = L_IDX[c]; terms.append((1, tuple(t)))
    return terms

def site_jumps(n, axes):
    return [tuple(L_IDX[axes[l]] if k == l else 0 for k in range(n)) for l in range(n)]

def string(n, letters_at):
    t = [0] * n
    for k, c in letters_at.items(): t[k] = L_IDX[c]
    return tuple(t)

def three_axis_census(n, edges, spectral):
    """Every bond-letter set times every three-axis assignment, read by F158's two ranks; on K3 also by the
    spectrum, whose verdicts must agree with the ranks and whose threshold is the measured gap."""
    res = {'one letter': [0, 0], 'mixed letters': [0, 0]}
    paired, broke, disagree = [], [], 0
    axsets = [a for a in itertools.product('XYZ', repeat=n) if set(a) == {'X', 'Y', 'Z'}]
    for letters in itertools.product('XYZ', repeat=len(edges)):
        kind = 'one letter' if len(set(letters)) == 1 else 'mixed letters'
        for ax in axsets:
            nr, fr = ends(n, letter_terms(n, edges, letters), site_jumps(n, ax))
            broken = not (nr == fr >= 1)
            res[kind][0] += 1; res[kind][1] += broken
            if spectral:
                H = sum(site(P[c], a, n) @ site(P[c], b, n) for (a, b), c in zip(edges, letters))
                dist = pairing_distance(H, [site(P[x], l, n) for l, x in enumerate(ax)], n)
                (broke if broken else paired).append(dist)
    return res, paired, broke

def main():
    print(__doc__.split('\n\n')[0])
    print()
    k3, paired, broke = three_axis_census(3, [(0, 1), (0, 2), (1, 2)], True)
    p4, _, _ = three_axis_census(4, [(0, 1), (1, 2), (2, 3)], False)
    for name, r in (('K3', k3), ('P4', p4)):
        print(f'1. {name} (two ranks mod two primes): one bond letter {r["one letter"][1]} of {r["one letter"][0]} break, '
              f'mixed letters {r["mixed letters"][1]} of {r["mixed letters"][0]} break')
    gap = np.log10(min(broke) / max(paired))
    print('   K3 spectrum (gamma = 0.05): rank-paired rows pair to at most %.1e, rank-broken rows miss by at least %.1e, '
          '%.1f decades apart' % (max(paired), min(broke), gap))
    rows = []
    rows.append(('K3 bonds X0X1, X1X2, Y0Y2, axes Y, X, Z, no field',
                 ends(3, letter_terms(3, [(0, 1), (1, 2), (0, 2)], 'XXY'), site_jumps(3, 'YXZ')), [2, 0]))
    xxxy = [(1, string(3, {a: 'X', b: 'X'})) for a, b in [(0, 1), (1, 2)]] +            [(1, string(3, {a: 'X', b: 'Y'})) for a, b in [(0, 1), (1, 2)]]
    rows.append(('2. P3 bonds XX + XY, Z on every site, no field', ends(3, xxxy, site_jumps(3, 'ZZZ')), [2, 0]))
    dm = [(100, string(3, {a: 'X', b: 'Y'})) for a, b in [(0, 1), (1, 2)]] +          [(-100, string(3, {a: 'Y', b: 'X'})) for a, b in [(0, 1), (1, 2)]]
    rows.append(('   DM chain XY - YX, N = 3, Z on every site, no field', ends(3, dm, site_jumps(3, 'ZZZ')), [4, 4]))
    dmf = dm + [(m, string(3, {l: 'X'})) for l, m in enumerate((30, 22, 41))]
    rows.append(('   DM chain with X fields 0.30, 0.22, 0.41 (bond weight 1)', ends(3, dmf, site_jumps(3, 'ZZZ')), [1, 0]))
    for name, got, want in rows:
        print(f'{name}: near {got[0]}, far {got[1]} ({"pairs" if got[0] == got[1] else "breaks"})')
    ok = (k3 == {'one letter': [18, 0], 'mixed letters': [144, 42]}
          and p4 == {'one letter': [108, 0], 'mixed letters': [864, 96]}
          and gap >= 6 and all(list(got) == want for _, got, want in rows))
    print('BOUNDARIES', 'CONFIRMED' if ok else 'NOT AS STATED')

if __name__ == '__main__':
    main()
