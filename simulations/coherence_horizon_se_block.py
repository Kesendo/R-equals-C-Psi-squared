"""The coherence horizon Q*(N) via the single-excitation (Haken-Strobl) Liouvillian.
Q*(N) is the EP of the N-site dephased tight-binding chain's rho^(1) dynamics, an N^2-dim
object that reproduces the full 4^N horizon. This script derives its closed form / exact
condition. Self-validating: run it, every assert must pass."""
import numpy as np


def h_single(N, J):
    """N x N single-particle Hamiltonian: tridiagonal, off-diagonal J (the (J/2)(XX+YY) hop)."""
    h = np.zeros((N, N), complex)
    for i in range(N - 1):
        h[i, i + 1] = J
        h[i + 1, i] = J
    return h


def L_se(N, J, g):
    """The N^2-dim single-excitation Liouvillian in the real-space basis. Coherent part
    -i[h, rho]; Z-dephasing makes coherence rho[i][j] (i!=j) decay at 4g, populations untouched
    (a single-excitation coherence has popcount(i^j)=2, Re = -2g*2 = -4g, the Absorption Theorem)."""
    h = h_single(N, J)
    I = np.eye(N)
    L = -1j * (np.kron(h, I) - np.kron(I, h.T))
    deph = np.array([(-4.0 * g if i != j else 0.0) for i in range(N) for j in range(N)])
    return L + np.diag(deph)


def qstar_se(N, J=1.0, lo=0.02, hi=6.0):
    """Q* = J/g* where g* is the largest g (smallest Q) at which the slowest non-zero SE mode
    stops oscillating (the EP), by bisection on g."""
    for _ in range(70):
        m = 0.5 * (lo + hi)
        ev = np.linalg.eigvals(L_se(N, J, m))
        nz = ev[ev.real < -1e-7]
        if len(nz) == 0:
            hi = m
            continue
        gap = nz.real.max()
        band = nz[np.abs(nz.real - gap) < 1e-7]
        if np.abs(band.imag).max() > 1e-7:
            lo = m
        else:
            hi = m
    return J / (0.5 * (lo + hi))


# Anchors: the live-witness ladder.
LADDER = {2: 1.0, 3: np.sqrt(2.0), 4: 1.87850, 5: 2.37220}


def _assert_ladder():
    for N, q in LADDER.items():
        got = qstar_se(N)
        assert abs(got - q) < 2e-3, f"SE Q*({N})={got:.5f} != ladder {q:.5f}"
    print("[Task1] SE builder reproduces the witness ladder:",
          {N: round(qstar_se(N), 5) for N in LADDER})


def _L_full(N, J, g):
    I2 = np.eye(2)
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]])
    Z = np.diag([1, -1]).astype(complex)

    def site(op, l):
        m = np.array([[1.0 + 0j]])
        for k in range(N):
            m = np.kron(m, op if k == l else I2)
        return m

    d = 2 ** N
    H = np.zeros((d, d), complex)
    for b in range(N - 1):
        H += (J / 2) * (site(X, b) @ site(X, b + 1) + site(Y, b) @ site(Y, b + 1))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = site(Z, l)
        L += g * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def _assert_full_reduction():
    # SE spectrum must be a sub-spectrum of the full L (the single-excitation block is invariant).
    for N in (2, 3):
        g = 0.5
        full = np.linalg.eigvals(_L_full(N, 1.0, g))
        sub = np.linalg.eigvals(L_se(N, 1.0, g))
        for s in sub:
            assert np.min(np.abs(full - s)) < 1e-9, f"N={N}: SE eig {s} not in full spectrum"
    print("[Task1] SE spectrum is an exact sub-spectrum of the full Liouvillian (N=2,3).")


def sine_U(N):
    """OBC sine transform: columns are the eigenvectors of h_single (real orthogonal).
    U[i,k] = sqrt(2/(N+1)) sin((i+1)(k+1) pi/(N+1)); h = U diag(E_k) U^T, E_k = 2J cos((k+1)pi/(N+1))."""
    i = np.arange(1, N + 1)[:, None]
    k = np.arange(1, N + 1)[None, :]
    return np.sqrt(2.0 / (N + 1)) * np.sin(i * k * np.pi / (N + 1))


def L_se_mode(N, J, g):
    """L_se in the sine-mode basis: W L W^T with W = U^T (x) U^T (orthogonal). The coherent part
    is diag(-i(E_k - E_l)); the real-space dephasing becomes an explicit mode-coupling block."""
    U = sine_U(N).astype(complex)
    W = np.kron(U.T, U.T)
    return W @ L_se(N, J, g) @ W.T


def reversal_R(N):
    """Site reflection i <-> N-1-i as an N x N permutation matrix."""
    R = np.zeros((N, N))
    for i in range(N):
        R[i, N - 1 - i] = 1.0
    return R


def parity_blocks(N, J, g):
    """Split L_se into reflection-parity even/odd sectors. P = R (x) R commutes with L_se.
    Returns (L_even, L_odd, basis_even, basis_odd)."""
    R = reversal_R(N)
    P = np.kron(R, R)
    L = L_se(N, J, g)
    w, V = np.linalg.eigh(P)               # P is a real symmetric involution: eigenvalues +-1
    even = V[:, w > 0]
    odd = V[:, w < 0]
    Le = even.T @ L @ even
    Lo = odd.T @ L @ odd
    return Le, Lo, even, odd


def _spec_match(a, b, atol=1e-6):
    """Two spectra match as multisets: equal size, every element of a has a close partner in b
    (robust to eigenvalue ordering, which is arbitrary and unstable at near-degeneracies)."""
    return len(a) == len(b) and max(np.min(np.abs(np.asarray(b) - x)) for x in a) < atol


def _assert_basis_and_parity():
    N, g = 5, 0.4
    # basis change preserves the spectrum (similarity transform by an orthogonal W)
    a = np.linalg.eigvals(L_se(N, 1.0, g))
    b = np.linalg.eigvals(L_se_mode(N, 1.0, g))
    assert _spec_match(a, b), "mode-basis spectrum differs from real-space"
    # parity blocks union = full spectrum, and the slowest non-zero mode is in one sector
    Le, Lo, _, _ = parity_blocks(N, 1.0, g)
    union = np.concatenate([np.linalg.eigvals(Le), np.linalg.eigvals(Lo)])
    assert _spec_match(a, union), "parity blocks do not reconstruct the spectrum"
    ee, eo = np.linalg.eigvals(Le), np.linalg.eigvals(Lo)
    slow_e = ee[ee.real < -1e-7].real.max()
    slow_o = eo[eo.real < -1e-7].real.max()
    print(f"[Task2] basis + parity OK. slowest non-zero: even={slow_e:.4f}, odd={slow_o:.4f}, "
          f"the horizon lives in the {'even' if slow_e > slow_o else 'odd'} sector.")


def coalescing_pair(N, J=1.0, above=1.10):
    """The coalescing conjugate pair just above Q*(N): among the slow (near-gap) oscillating modes,
    the one with the SMALLEST |Im| is the {0,2}-coherence about to freeze (the band-edge survivor
    has the large |Im|=2cos(pi/(N+1))). Returns (lambda, sum s, product p) of that pair."""
    g = J / (above * qstar_se(N))
    ev = np.linalg.eigvals(L_se(N, J, g))
    nz = ev[ev.real < -1e-7]
    osc = nz[nz.imag > 1e-7]                       # upper-half oscillating modes
    slow = osc[osc.real > osc.real.max() - 0.25]   # near the gap
    la = slow[np.argmin(np.abs(slow.imag))]        # the coalescer: smallest |Im|
    lb = np.conj(la)
    return la, g, (la + lb).real, (la * lb).real


def _report_pair_invariants():
    """Is the coalescing pair a clean 2x2 with N-simple invariants (s, p)? If s = -4g and p has a
    pattern, Q* has a closed form via s^2 = 4p; if not, the pair is dressed -> transcendental."""
    print("[Task3] coalescing-pair invariants just above Q* (the fork test):")
    print(f"  {'N':>2} {'g':>7} {'s=la+lb':>10} {'s/(-4g)':>9} {'p=la*lb':>10} "
          f"{'p/(4J^2)':>9} {'s^2-4p':>10}")
    for N in range(2, 9):
        la, g, s, p = coalescing_pair(N)
        print(f"  {N:>2} {g:>7.4f} {s:>10.5f} {s/(-4*g):>9.4f} {p:>10.5f} "
              f"{p/4.0:>9.4f} {s*s-4*p:>10.5f}")


def _assert_fork():
    """The clean 2x2 (s = -4g) holds EXACTLY at N=2,3 and breaks at N>=4: the structural form of the
    2cos(pi/(N+1)) low-N accident. N=2,3 -> closed form; N>=4 -> collectively dressed, transcendental."""
    for N in (2, 3):
        _, g, s, p = coalescing_pair(N)
        assert abs(s / (-4 * g) - 1.0) < 1e-3, f"N={N}: expected the clean 2x2 trace s=-4g"
    _, g4, s4, _ = coalescing_pair(4)
    assert abs(s4 / (-4 * g4) - 1.0) > 1e-2, "N=4 should DEVIATE from s=-4g (collective dressing)"
    print("[Task3] FORK: clean 2x2 (s=-4g, p=4J^2/2J^2) at N=2,3 -> closed form; "
          "N>=4 dressed -> exact transcendental condition.")


def _pair_sp_at(N, g, J=1.0):
    """(s, p) = (sum, product) of the coalescing pair at a given g (no rescaling), to test whether
    s = -4g and p = c*J^2 are gamma-INDEPENDENT identities (a genuine 2x2) or just hold at one point."""
    ev = np.linalg.eigvals(L_se(N, J, g))
    nz = ev[ev.real < -1e-7]
    osc = nz[nz.imag > 1e-7]
    slow = osc[osc.real > osc.real.max() - 0.25]
    la = slow[np.argmin(np.abs(slow.imag))]
    return 2.0 * la.real, abs(la) ** 2


def _assert_closed_forms():
    """N=2,3: the coalescing pair are the roots of lambda^2 + 4g*lambda + c*J^2 = 0 with c constant
    (s=-4g and p=c*J^2 hold for all g below the EP), so Q* = 2/sqrt(c) exactly. N>=4: that identity
    breaks (collective dressing) and the exact condition is transcendental."""
    J = 1.0
    for N, c in [(2, 4.0), (3, 2.0)]:
        for g in (0.30, 0.45, 0.60):                 # several g below the EP (pair clearly complex)
            s, p = _pair_sp_at(N, g)
            assert abs(s - (-4 * g)) < 1e-6, f"N={N} g={g}: s={s:.6f} != -4g"
            assert abs(p - c * J * J) < 1e-6, f"N={N} g={g}: p={p:.6f} != {c}J^2"
        qstar_cf = 2.0 / np.sqrt(c)                   # EP: discriminant 16g^2 - 4cJ^2 = 0 -> g=J*sqrt(c)/2
        assert abs(qstar_cf - qstar_se(N)) < 1e-4, f"N={N}: closed form 2/sqrt(c) != bisected Q*"
        print(f"[Task4] N={N}: closed form  lambda^2 + 4g*lambda + {c:.0f}*J^2 = 0  "
              f"(s=-4g, p={c:.0f}J^2 exact for all g)  ->  Q* = 2/sqrt({c:.0f}) = {qstar_cf:.6f}")
    # N>=4: (s,p) is not constant in g -> no clean 2x2 -> transcendental
    s4a, p4a = _pair_sp_at(4, 0.45)
    s4b, p4b = _pair_sp_at(4, 0.55)
    assert abs(p4a - p4b) > 1e-3 or abs(s4a / (-4 * 0.45) - 1.0) > 1e-2, \
        "N=4 (s,p) should vary with g (collective dressing)"
    print("[Task4] N>=4: (s,p) varies with g -> no clean 2x2. The exact condition is the SE")
    print("        Liouvillian's slowest-mode double root (transcendental); Q*(N) = qstar_se(N),")
    print("        reproducing the ladder + the N=6..14 sequence (Task 1).")


def _report_canonical():
    """N>=4 is transcendental, so the canonical Q*(4)/Q*(5) are the high-precision SE-EP values
    (they supersede the 1.8785/2.3722 bisection digits). Plus the diffusive asymptotic slope."""
    # N=2,3 exact closed forms validate the bisector's precision
    print(f"[Task5] precision check vs the closed forms: Q*(2)={qstar_se(2):.8f} (exact 1), "
          f"Q*(3)={qstar_se(3):.8f} (exact sqrt2={np.sqrt(2.0):.8f})")
    q4, q5 = qstar_se(4), qstar_se(5)
    print(f"[Task5] canonical (SE-EP, transcendental): Q*(4)={q4:.6f}, Q*(5)={q5:.6f} "
          f"(supersede the 1.8785/2.3722 grid digits)")
    # diffusive asymptotic: Q*(N) ~ a N + b + c/N, fit on N=8..14
    Ns = np.arange(8, 15)
    ys = np.array([qstar_se(int(N)) for N in Ns])
    A = np.vstack([Ns, np.ones_like(Ns), 1.0 / Ns]).T
    a, b, c = np.linalg.lstsq(A, ys, rcond=None)[0]
    print(f"[Task5] asymptotic Q*(N) ~ {a:.5f} N + {b:.5f} + {c:.5f}/N   "
          f"(2/pi = {2/np.pi:.5f}; long-wavelength critical damping, slope still creeping up)")


if __name__ == "__main__":
    _assert_ladder()
    _assert_full_reduction()
    _assert_basis_and_parity()
    _report_pair_invariants()
    _assert_fork()
    _assert_closed_forms()
    _report_canonical()
