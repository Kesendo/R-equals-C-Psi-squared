"""Gate from below: what IS the star's frozen slow mode, and is it a [H,A]=0 commutant?

The lens-agent (outer) claimed the committed docs' '[H,A]=0 commutant, kernel of ad_H' mechanism is
WRONG (||[H,rho_slow]||~2g != 0), the real reason graph-spectral (flat a=0 band). The skill says:
a contradiction the agent surfaces is a FLAG to gate, not a fact. So we compute, from inside.

RESOLVED 2026-06-20 (see star_survivor_heisenberg.py + THE_STAR_FROZEN_SEAM.md 'Model scope'): the docs
are NOT wrong. The committed mechanism is correct -- the survivor IS the [H,rho]=0 commutant, but only in
the high-Q LIMIT (||[H,rho]|| ~ 1/Q -> 0), frozen at every Q. This gate is HEISENBERG (H has the ZZ term),
so its survivor darkness is 4/N; the docs' 4/(N-1) is the XY (hopping-only) value, by design. The
||[H,rho]||~2g the lens-agent saw is the finite-Q (non-limit) residual, not a refutation of the commutant.
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)


def op(P, l, N):
    m = np.array([[1.0 + 0j]])
    for k in range(N):
        m = np.kron(m, P if k == l else I2)
    return m


def run(N, J, g):
    bonds = [(0, l) for l in range(1, N)]            # star: center 0 -> leaves
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        H += J * (op(X, i, N) @ op(X, j, N) + op(Y, i, N) @ op(Y, j, N) + op(Z, i, N) @ op(Z, j, N))
    Id = np.eye(d, dtype=complex)
    Lh = -1j * (np.kron(Id, H) - np.kron(H.T, Id))   # -i[H,.], column-stacking vec
    Ld = np.zeros((d * d, d * d), dtype=complex)
    for l in range(N):
        Zl = op(Z, l, N)
        Ld += g * (np.kron(Zl.T, Zl) - np.kron(Id, Id))   # Z_l rho Z_l - rho
    L = Lh + Ld
    w, V = np.linalg.eig(L)
    re = w.real
    cand = [k for k in range(len(w)) if re[k] < -1e-9]
    k = max(cand, key=lambda k: re[k])
    lam = w[k]
    rho = V[:, k].reshape(d, d, order="F")            # column-stacking unvec

    herm = np.linalg.norm(rho - rho.conj().T) / np.linalg.norm(rho)
    comm = np.linalg.norm(H @ rho - rho @ H) / np.linalg.norm(rho)
    diagwt = np.linalg.norm(np.diag(rho)) / np.linalg.norm(rho)

    print(f"\n==== STAR N={N}  J={J} g={g}  (Q={J/g:.3g}) ====", flush=True)
    print(f"slowest decaying: Re={lam.real:+.5f}  |Im|={abs(lam.imag):.5f}   "
          f"[frozen if |Im|=0]", flush=True)
    print(f"rho Hermitian residual: {herm:.2e}   (small => reshape/convention ok)", flush=True)
    print(f"||[H,rho]||/||rho|| = {comm:.4f}    (2g = {2*g}; ~0 => commutant/ad_H-kernel, "
          f"!=0 => NOT commutant)", flush=True)
    print(f"diagonal (population) weight: {diagwt:.3f}", flush=True)

    # sector content: where does |rho[a,b]|^2 live in (popcount_bra, popcount_ket)?
    pc = np.array([bin(a).count("1") for a in range(d)])
    P2 = np.abs(rho) ** 2
    secw = {}
    for a in range(d):
        for b in range(d):
            if P2[a, b] > 1e-10:
                key = (pc[a], pc[b])
                secw[key] = secw.get(key, 0.0) + P2[a, b]
    tot = sum(secw.values())
    top = sorted(secw.items(), key=lambda kv: -kv[1])[:5]
    print("top (p_bra,p_ket) sectors of rho:  " +
          "  ".join(f"{s}:{v/tot:.2f}" for s, v in top), flush=True)

    # Z-profile of the mode (hub vs leaf): Tr(rho Z_l)
    zp = [float(np.real(np.trace(rho @ op(Z, l, N)))) for l in range(N)]
    mx = max(abs(x) for x in zp) or 1.0
    print("Tr(rho Z_l), l=0(hub)..leaves (norm):  [" +
          " ".join(f"{x/mx:+.3f}" for x in zp) + "]", flush=True)

    # graph-spectral check: star single-excitation adjacency spectrum (the claimed flat band)
    A = np.zeros((N, N))
    for (i, j) in bonds:
        A[i, j] = 1.0
        A[j, i] = 1.0
    aval = np.round(np.linalg.eigvalsh(A), 4)
    print(f"star adjacency eigenvalues a: {aval}   (flat a=0 band mult {sum(abs(aval)<1e-9)}, "
          f"|a|_max={max(abs(aval)):.4f}=sqrt(N-1)={np.sqrt(N-1):.4f})", flush=True)
    # Haken-Strobl single-excitation rate per a: -g +- sqrt(g^2-(J a)^2); real(frozen) iff J|a|<g
    for a in sorted(set(np.round(np.abs(aval), 4))):
        disc = g * g - (J * a) ** 2
        if disc >= 0:
            print(f"   a={a:.3f}: real rates {-g+np.sqrt(disc):+.4f}, {-g-np.sqrt(disc):+.4f}  (FROZEN, J|a|<g)", flush=True)
        else:
            print(f"   a={a:.3f}: complex rate -{g} +- i{np.sqrt(-disc):.4f}  (OSCILLATES, J|a|>g)", flush=True)


if __name__ == "__main__":
    for g in (1.0, 0.5, 2.0):     # Q = 1, 2, 0.5
        run(5, 1.0, g)
