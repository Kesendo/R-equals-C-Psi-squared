"""SKEPTIC (gate-first): the INDEPENDENT-CONVENTION adversarial companion to chain_gap_dominance.py.

PROOF_CHAIN_GAP_DOMINANCE is Tier1Derived (it graduated ClockHandLadderClaim + TopologyBandEdgeClaim). This
verifier re-derives it FROM SCRATCH with deliberately DIFFERENT conventions than chain_gap_dominance.py, so a
convention-dependent error in either would surface as a disagreement:

  - column-stack vectorisation vec(A)=A.flatten('F'), so L = -i(I⊗H − Hᵀ⊗I) + γΣ(Zᵀ⊗Z − I)
    (the OTHER stacking than the verifier's row-stack kron(H,I)−kron(I,Hᵀ));
  - σ⁻ and the Jordan-Wigner string built independently;
  - n_XY read as a Pauli-weight distribution (Tr(P†A)), NOT borrowed from the Absorption Theorem.

It also adds two checks the row-stack verifier does not carry, both physics-first-review moves:
  - PURITY (Attack A): the exact-(−2γ) eigenspace lies INSIDE the pure n_XY=1 Pauli subspace for N≥4 (the
    'L_D is scalar ⟹ free' step); N=3 LEAKS (its 4 extra {0,2} modes — the documented N=3 special).
  - REGIME (Attack B): the subspace max|Im|@−2γ = E1 holds at γ≪J, γ~J AND γ≫J — the −2γ band's internal
    structure is γ-independent (eigenvalues −2γ±iE_k exactly), so gap-dominance is not a γ≪J artifact.

(Attack D — the N=3-vanishing MECHANISM — is deliberately NOT gated here: its probe is grading-sensitive
(Pauli n_XY weight vs the doc's n_diff/coherence grading, total weight = n_diff + Z-shadow) and too loose to
be a reliable gate. The N=3 special itself is captured by the PURITY leak below + the n3 small-N catalogue.)

Run: python simulations/chain_gap_dominance_skeptic.py
"""
import sys
import itertools
import numpy as np
from math import cos, pi

sys.stdout.reconfigure(encoding="utf-8")
np.set_printoptions(precision=5, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
SM = 0.5 * (X + 1j * Y)            # [[0,1],[0,0]] = σ⁻ (lowers |1>→|0>)
PAULI = {0: I2, 1: X, 2: Y, 3: Z}
TOL = 1e-7


def kron(ops):
    o = np.array([[1.0 + 0j]])
    for m in ops:
        o = np.kron(o, m)
    return o


def op_at(M, l, N):
    return kron([M if k == l else I2 for k in range(N)])


def H_XY(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), complex)
    for i in range(N - 1):
        H += (J / 2) * (op_at(X, i, N) @ op_at(X, i + 1, N) + op_at(Y, i, N) @ op_at(Y, i + 1, N))
    return H


def jw(j, N):
    """c_j = (∏_{l<j} Z_l) σ⁻_j (independent Jordan-Wigner construction)."""
    return kron([Z if l < j else SM if l == j else I2 for l in range(N)])


def vecF(A):
    return A.flatten("F")               # column-stack (the OTHER convention)


def unvecF(v, d):
    return v.reshape((d, d), order="F")


def superL_col(N, g, J=1.0):
    """Column-stack Liouvillian: L = -i(I⊗H − Hᵀ⊗I) + γΣ(Zᵀ⊗Z − I)."""
    H = H_XY(N, J)
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for l in range(N):
        Zl = op_at(Z, l, N)
        L += g * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    return L


def band(N, J=1.0):
    h = np.zeros((N, N))
    for i in range(N - 1):
        h[i, i + 1] = h[i + 1, i] = J
    return np.linalg.eigh(h)


def nxy_mean(A, N):
    """HS-weighted mean Pauli n_XY of an operator A (mixture-aware)."""
    d = 2 ** N
    num = den = 0.0
    for idx in itertools.product(range(4), repeat=N):
        c = np.trace(kron([PAULI[t] for t in idx]).conj().T @ A) / d
        w = abs(c) ** 2
        if w < 1e-18:
            continue
        num += w * sum(1 for t in idx if t in (1, 2))
        den += w
    return num / den if den > 1e-15 else 0.0


def nxy1_basis(N):
    """Orthonormal vecF basis of the n_XY=1 Pauli operators."""
    cols = []
    for idx in itertools.product(range(4), repeat=N):
        if sum(1 for t in idx if t in (1, 2)) == 1:
            cols.append(vecF(kron([PAULI[t] for t in idx])) / (2 ** (N / 2)))
    return np.array(cols).T


def exact_minus2g_eigenspace_leakage(N, g):
    """Max leakage of the exact-(−2γ) eigenspace outside the pure n_XY=1 subspace (basis-independent)."""
    d = 2 ** N
    ev, evec = np.linalg.eig(superL_col(N, g))
    idxs = np.where(np.abs(ev.real - (-2 * g)) < TOL)[0]
    Q, _ = np.linalg.qr(nxy1_basis(N))
    leak = 0.0
    for i in idxs:
        v = evec[:, i] / np.linalg.norm(evec[:, i])
        leak = max(leak, np.linalg.norm(v - Q @ (Q.conj().T @ v)))
    return len(idxs), leak


g = 0.05
print("=" * 92)
print("SKEPTIC — independent (column-stack) re-derivation + attacks on PROOF_CHAIN_GAP_DOMINANCE")
print("=" * 92)

# ---- Stage 1: JW single-particle modes (independent conventions) ----
print("\n[1] JW modes: band = 2cos(πk/(N+1)), [H,c_k] = −E_k c_k, every c_k pure n_XY=1")
for N in (3, 4, 5):
    H = H_XY(N)
    w, v = band(N)
    cj = [jw(j, N) for j in range(N)]
    ck = [sum(v[j, k] * cj[j] for j in range(N)) for k in range(N)]
    Ek = sorted(2 * cos(pi * (k + 1) / (N + 1)) for k in range(N))
    assert np.allclose(sorted(w), Ek), f"[1] FIRED N={N}: band != 2cos"
    assert all(np.allclose(H @ ck[k] - ck[k] @ H, -w[k] * ck[k], atol=1e-9) for k in range(N)), \
        f"[1] FIRED N={N}: [H,c_k] != -E_k c_k"
    assert all(abs(nxy_mean(ck[k], N) - 1.0) < 1e-9 for k in range(N)), f"[1] FIRED N={N}: c_k not pure n_XY=1"
    print(f"    N={N}: band==2cos ✓  [H,c_k]=−E_k c_k ✓  all c_k n_XY=1 ✓  E1={max(Ek):.6f}")

# ---- Stage 2: c_k·P_m and c_k†·P_m are EXACT −2γ∓iE_k modes (col-stack) ----
print("\n[2] c_k·P_m, c_k†·P_m are exact L-eigenmodes at −2γ±iE_k (0 failures)")
for N in (4, 5):
    w, v = band(N)
    cj = [jw(j, N) for j in range(N)]
    ck = [sum(v[j, k] * cj[j] for j in range(N)) for k in range(N)]
    Ntot = sum(c.conj().T @ c for c in cj)
    dN = np.round(np.real(np.diag(Ntot))).astype(int)
    Pm = [np.diag((dN == m).astype(complex)) for m in range(N + 1) if (dN == m).any()]
    L = superL_col(N, g)
    bad = tot = 0
    for k in range(N):
        for P in Pm:
            for A, lam in ((ck[k] @ P, -2 * g + 1j * w[k]), (ck[k].conj().T @ P, -2 * g - 1j * w[k])):
                if np.linalg.norm(A) < 1e-12:
                    continue
                vv = vecF(A)
                tot += 1
                if np.linalg.norm(L @ vv - lam * vv) / max(1, np.linalg.norm(vv)) > 1e-8:
                    bad += 1
    assert bad == 0, f"[2] FIRED N={N}: {bad}/{tot} generators not exact −2γ±iE_k modes"
    print(f"    N={N}: {tot} generators, {bad} failures ✓")

# ---- Stage 3: spanning — dim(exact −2γ eigenspace) == rank{vec(c_k^(†) P_m)} ----
print("\n[3] spanning: dim(exact −2γ subspace) == rank{vec generators} (32/50 at N=4/5)")
for N, expect in ((4, 32), (5, 50)):
    w, v = band(N)
    cj = [jw(j, N) for j in range(N)]
    ck = [sum(v[j, k] * cj[j] for j in range(N)) for k in range(N)]
    Ntot = sum(c.conj().T @ c for c in cj)
    dN = np.round(np.real(np.diag(Ntot))).astype(int)
    Pm = [np.diag((dN == m).astype(complex)) for m in range(N + 1) if (dN == m).any()]
    L = superL_col(N, g)
    gens = [vecF(A) for k in range(N) for P in Pm for A in (ck[k] @ P, ck[k].conj().T @ P)
            if np.linalg.norm(A) > 1e-12]
    dim_span = np.linalg.matrix_rank(np.array(gens).T, tol=1e-9)
    ev = np.linalg.eigvals(L)
    dim_sub = int(np.sum(np.abs(ev.real - (-2 * g)) < TOL))
    assert dim_span == dim_sub == expect, f"[3] FIRED N={N}: span={dim_span} sub={dim_sub} expect={expect}"
    print(f"    N={N}: dim_sub={dim_sub} == dim_span={dim_span} == {expect} ✓")

# ---- Stage 4: max|Im| over exact −2γ modes == E1 ----
print("\n[4] max|Im| over the exact −2γ subspace == E1")
for N in (3, 4, 5):
    E1 = 2 * cos(pi / (N + 1))
    ev = np.linalg.eigvals(superL_col(N, g))
    mx = np.abs(ev.imag[np.abs(ev.real - (-2 * g)) < TOL]).max()
    assert abs(mx - E1) < 1e-6, f"[4] FIRED N={N}: max|Im|={mx} != E1={E1}"
    print(f"    N={N}: max|Im|={mx:.6f} == E1={E1:.6f} ✓")

# ---- Stage 5 (Attack A): purity — N≥4 eigenspace inside n_XY=1; N=3 leaks (its 4 extras) ----
print("\n[5] PURITY: exact −2γ eigenspace ⊂ pure n_XY=1 for N≥4; N=3 LEAKS (the documented special)")
dim3, leak3 = exact_minus2g_eigenspace_leakage(3, g)
assert leak3 > 0.5, f"[5] FIRED: N=3 should leak (the 4 {{0,2}} extras), got {leak3:.2e}"
print(f"    N=3: dim={dim3}, leakage={leak3:.2e}  -> LEAKS ✓ (the 4 extra {{0,2}} modes — N=3 special)")
for N in (4, 5):
    dim, leak = exact_minus2g_eigenspace_leakage(N, g)
    assert leak < 1e-7, f"[5] FIRED N={N}: eigenspace leaks outside n_XY=1 by {leak:.2e}"
    print(f"    N={N}: dim={dim}, leakage={leak:.2e}  -> PURE n_XY=1 ✓ (the 'L_D scalar ⟹ free' step holds)")

# ---- Stage 6 (Attack B): regime — subspace max|Im|@−2γ == E1 at γ≪J, γ~J, γ≫J ----
print("\n[6] REGIME: subspace max|Im|@−2γ == E1 across γ≪J, γ~J, γ≫J (γ-independent, not a small-γ artifact)")
for (gg, J, tag) in [(0.05, 1.0, "γ≪J"), (1.0, 1.0, "γ~J"), (5.0, 1.0, "γ≫J")]:
    for N in (4, 5):
        E1 = 2 * J * cos(pi / (N + 1))
        ev = np.linalg.eigvals(superL_col(N, gg, J))
        mx = np.abs(ev.imag[np.abs(ev.real - (-2 * gg)) < 1e-7]).max()
        assert abs(mx - E1) < 1e-6, f"[6] FIRED {tag} N={N}: max|Im|={mx} != E1={E1}"
    print(f"    {tag:4s}: N=4,5 subspace max|Im|@−2γ == E1 ✓")

print("\n" + "=" * 92)
print("ALL STAGES PASS — PROOF_CHAIN_GAP_DOMINANCE confirmed by an INDEPENDENT (column-stack) re-derivation,")
print("a basis-independent purity gate (N≥4 pure n_XY=1, N=3 the documented leak), and a γ-regime sweep.")
print("Two conventions agree (row-stack chain_gap_dominance.py + this col-stack), so it is not a stacking")
print("artifact; the −2γ band's max|Im|=E1 is γ-independent. DONE.")
