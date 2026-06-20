"""Gate-first: is the half-filling SURVIVOR a fixed point of a particle-hole self-mirror,
and is THAT why it is dark?  (a lens-generated conjecture handed to the physics from below)

THE DIRECTION (from above, borrowing-a-discipline / condensed matter):
  the longest-lived survivor coherence at half-filling is the PARTICLE-HOLE SELF-MIRROR state
  (half-filling = the bipartite particle-hole fixed point), and being mirror-fixed is WHY it is dark.

THE SURVEY (from below, 3-agent) sharpened WHICH mirror -- and pre-flagged one branch as a trap:
  CONTROL  X = X^(x)N  (uniform spin-flip = Pi^2):  X-self-paired modes are BRIGHT (<n_XY>=N/2),
           the survivor is DARK (<n_XY>->0).  X commutes with L, so the survivor has a definite
           X-parity -- PREDICTED ODD (-1): an ANTI-fixed point, NOT the lens's "+1 fixed".
  LIVE     U = X^(x)N . prod_{l odd} Z_l  (the staggered/bipartite half-filling particle-hole
           c_i -> (-1)^i c_i^dag):  UNBUILT in the repo (no staggered operator exists anywhere).
           The real condensed-matter tool.  Maps XX+YY -> -(XX+YY): an E->-E symmetry of XY,
           NOT a clean symmetry of Heisenberg.  Does it fix the survivor?  (open ground)
  BONUS    R = spatial reflection (site l -> N-1-l):  the survivor's ACTUAL typed symmetry
           (density gradient mirror-symmetric, SurvivorDiffusionGradientClaim).  Self-mirror, spatial.

GATES (gate-first; a FIRING gate is the FIND -- diagnose, do not loosen):
  G0  the survivor here IS the dark half-filling (N/2,N/2) interior mode (Heisenberg chain, low Q
      below the coherence horizon).  If not, the conjecture's premise is void in this regime.
  G1  the INTERTWINER: for each operator S, how does the conjugation superoperator relate to L
      (commute / anticommute)?  The relationship sets the MEANING of G2.
  G2  the FIXED-POINT PROBE: is the survivor an eigenvector of S (parity +-1)?
      PREDICTION: X -> -1 (dark anti-fixed, control), R -> +1 (the real self-mirror),
      U -> the open question.  A clean separation X(-1)/R(+1) IS the find.
  G3  CAUSE vs CORRELATION: across the (p,p) spectrum, does S-parity sort DARK from BRIGHT
      (<n_XY> = -Re(lambda)/2gamma, the Absorption Theorem)?  Is darkness == a definite parity?
  MATH-LENS GUARD: U and X must COMPLEMENT (|a> -> |~a>), non-trivial on diagonals;
      a pure Z-string would fix ANY diagonal mode trivially.  Assert they actually complement.
"""
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------- lattice + sectors
def bonds(N, topo):
    if topo == "chain":
        return [(i, i + 1) for i in range(N - 1)]
    if topo == "ring":
        return [(i, (i + 1) % N) for i in range(N)]
    if topo == "star":
        return [(0, i) for i in range(1, N)]
    raise ValueError(topo)


def pbasis(N, p):
    return [a for a in range(1 << N) if bin(a).count("1") == p]


def zval(a, l):                                    # Z eigenvalue at site l: +1 empty, -1 occupied
    return 1.0 - 2.0 * ((a >> l) & 1)


def sector_H(N, bs, J, bnds, model):
    """p-excitation-sector Hamiltonian on the basis bs.  XX+YY hopping (amp 2J) + ZZ diagonal (heisenberg)."""
    d = len(bs)
    idx = {a: i for i, a in enumerate(bs)}
    H = np.zeros((d, d))
    for i, a in enumerate(bs):
        if model == "heisenberg":
            H[i, i] += J * sum(zval(a, u) * zval(a, v) for (u, v) in bnds)
        for (u, v) in bnds:
            if ((a >> u) & 1) != ((a >> v) & 1):
                a2 = a ^ (1 << u) ^ (1 << v)
                H[idx[a2], i] += 2.0 * J
    return H


def block_L(N, prow, pcol, J, g, bnds, model):
    """The (prow,pcol) Liouvillian block on coherences rho_{ab}, a in basis(prow), b in basis(pcol).
    L = -i(H_r rho - rho H_c) + sum_l g (Z_l rho Z_l - rho).  Row-major vec: v = i*dc + j."""
    br, bc = pbasis(N, prow), pbasis(N, pcol)
    Hr, Hc = sector_H(N, br, J, bnds, model), sector_H(N, bc, J, bnds, model)
    dr, dc = len(br), len(bc)
    D = dr * dc
    L = np.zeros((D, D), dtype=complex)
    for i in range(dr):
        for j in range(dc):
            v = i * dc + j
            for k in range(dr):                    # -i (H_r rho)
                L[k * dc + j, v] += -1j * Hr[k, i]
            for k in range(dc):                    # +i (rho H_c)
                L[i * dc + k, v] += 1j * Hc[k, j]
            L[v, v] += sum(g * (zval(br[i], l) * zval(bc[j], l) - 1.0) for l in range(N))
    return L, br, bc


# ---------------------------------------------------------------- the mirror operators
def perm_phase(N, bs, kind):
    """U|bs[i]> = phi[i] |bs[pi[i]]> for a complementing/reflecting operator on the p-sector basis bs."""
    idx = {a: i for i, a in enumerate(bs)}
    mask = (1 << N) - 1
    pi = np.zeros(len(bs), dtype=int)
    phi = np.ones(len(bs), dtype=complex)
    for i, a in enumerate(bs):
        if kind == "X":                            # X^(x)N : complement
            pi[i] = idx[a ^ mask]
        elif kind == "UPH":                        # X^(x)N . prod_{odd} Z : complement w/ staggered phase
            ph = 1.0
            for l in range(N):
                if l % 2 == 1:
                    ph *= zval(a, l)
            pi[i] = idx[a ^ mask]
            phi[i] = ph
        elif kind == "R":                          # spatial reflection l -> N-1-l
            a2 = 0
            for l in range(N):
                if (a >> l) & 1:
                    a2 |= 1 << (N - 1 - l)
            pi[i] = idx[a2]
        else:
            raise ValueError(kind)
    return pi, phi


def superop(bs, kind, N):
    """Conjugation superoperator S: vec(U rho U^dag) = S vec(rho) (row-major), for a (p,p) block."""
    pi, phi = perm_phase(N, bs, kind)
    d = len(bs)
    S = np.zeros((d * d, d * d), dtype=complex)
    for i in range(d):
        for j in range(d):
            S[pi[i] * d + pi[j], i * d + j] = phi[i] * np.conj(phi[j])
    return S


# ---------------------------------------------------------------- spectral helpers
def slowest_eig(L):
    """Slowest strictly-DECAYING mode (largest Re below 0). Returns (None, None) if the block is
    entirely stationary (e.g. the 1-dim fully-polarized sector) -- it then has no survivor coherence."""
    w, V = np.linalg.eig(L)
    cand = [k for k in range(len(w)) if w[k].real < -1e-9]      # exclude ~stationary steady states
    if not cand:
        return None, None
    k = max(cand, key=lambda k: w[k].real)
    return w[k], V[:, k]


def parity(v, S):
    Sv = S @ v
    fid = abs(np.vdot(v, Sv)) / (np.linalg.norm(v) * np.linalg.norm(Sv) + 1e-300)
    sign = float(np.sign(np.real(np.vdot(v, Sv))))
    return fid, sign


def rel(A, B):
    return np.linalg.norm(A - B) / (np.linalg.norm(B) + 1e-300)


# ---------------------------------------------------------------- the gates
def g0_survivor_sector(N, J, g, bnds, model):
    """Scan (p,p) diagonal sectors + the (0,1) band edge; return (Re, sector) of the GLOBAL survivor."""
    best = None
    for p in range(1, N + 1):
        L, _, _ = block_L(N, p, p, J, g, bnds, model)
        lam, _ = slowest_eig(L)
        if lam is None:                                        # entirely stationary block (no survivor)
            continue
        if best is None or lam.real > best[0]:
            best = (lam.real, (p, p))
    L01, _, _ = block_L(N, 0, 1, J, g, bnds, model)
    lam01, _ = slowest_eig(L01)
    if lam01 is not None and (best is None or lam01.real > best[0]):
        best = (lam01.real, (0, 1))
    return best


def run(N, J, g, topo, model):
    bnds = bonds(N, topo)
    Q = J / g
    half = N // 2
    print(f"\n========  N={N}  {topo}  {model}  J={J} g={g}  (Q=J/g={Q:.3g})  ========", flush=True)

    # G0 -- is the global survivor the dark half-filling (half,half) interior mode?
    re_glob, sec = g0_survivor_sector(N, J, g, bnds, model)
    nxy_glob = -re_glob / (2 * g)
    half_ok = (sec == (half, half))
    print(f"G0  global survivor: sector {sec}, Re={re_glob:+.5f}, <n_XY>={nxy_glob:.4f}  "
          f"[half-filling ({half},{half})? {'YES' if half_ok else 'NO -> premise void here'}]", flush=True)
    if N % 2 == 1:
        print("    (odd N: no (N/2,N/2) sector exists; particle-hole has no fixed sector -- skipping fixed-point test)", flush=True)
        return
    if not half_ok:
        print("    (survivor is NOT at half-filling in this regime; G2 would test the wrong mode -- skipping)", flush=True)
        return

    # build the half-filling block + its survivor
    L, bs, _ = block_L(N, half, half, J, g, bnds, model)
    lam_s, v_s = slowest_eig(L)
    print(f"    half-filling block survivor: Re={lam_s.real:+.5f}, |Im|={abs(lam_s.imag):.5f}", flush=True)

    ops = {}
    for kind in ("X", "UPH", "R"):
        ops[kind] = superop(bs, kind, N)

    # MATH-LENS GUARD: X and UPH must COMPLEMENT (move a localized diagonal population), not act trivially
    d = len(bs)
    test = np.zeros(d * d, dtype=complex)
    test[0 * d + 0] = 1.0                                   # |bs[0]><bs[0]|, a localized population
    for kind in ("X", "UPH"):
        moved = np.linalg.norm(ops[kind] @ test - test) > 1e-9
        assert moved, f"GUARD FAIL: {kind} acts trivially on a diagonal mode (it is not complementing)"
    print(f"    [math-lens guard ok: X and UPH genuinely complement diagonals; "
          f"UPH != X (staggered phase present: {np.linalg.norm(ops['UPH']-ops['X'])>1e-9})]", flush=True)

    # G1 -- the intertwiner: how does each S relate to L?
    print("G1  intertwiner  (||SL-LS||/||L|| commute,  ||SL+LS||/||L|| anticommute):", flush=True)
    for kind, S in ops.items():
        comm = rel(S @ L, L @ S)
        anti = np.linalg.norm(S @ L + L @ S) / (np.linalg.norm(L) + 1e-300)
        tag = "COMMUTES" if comm < 1e-8 else ("ANTICOMMUTES" if anti < 1e-8 else "neither (mixed)")
        print(f"      {kind:4s}: comm={comm:.2e}  anti={anti:.2e}   -> {tag}", flush=True)

    # G2 -- the fixed-point probe on the survivor
    print("G2  survivor parity under each mirror (fid~1 => eigenvector; sign = +1 fixed / -1 anti-fixed):", flush=True)
    verdicts = {}
    for kind, S in ops.items():
        fid, sign = parity(v_s, S)
        eig = "+1 FIXED" if (fid > 0.99 and sign > 0) else ("-1 ANTI-FIXED" if (fid > 0.99 and sign < 0) else "not an eigenvector")
        verdicts[kind] = eig
        print(f"      {kind:4s}: fidelity={fid:.4f}  sign={sign:+.0f}   -> {eig}", flush=True)

    # G3 -- cause vs correlation: does parity sort dark from bright across the block?
    print("G3  parity vs darkness across the half-filling block (<n_XY> = -Re/2g):", flush=True)
    w, V = np.linalg.eig(L)
    for kind, S in ops.items():
        if rel(S @ L, L @ S) > 1e-8:
            print(f"      {kind:4s}: S does not commute with L -> parity undefined across modes (skip)", flush=True)
            continue
        plus_nxy, minus_nxy = [], []
        for k in range(len(w)):
            vk = V[:, k]
            fid, sign = parity(vk, S)
            if fid < 0.99:
                continue                                    # degenerate-mixed, no definite parity
            nxy = -w[k].real / (2 * g)
            (plus_nxy if sign > 0 else minus_nxy).append(nxy)
        mp = np.mean(plus_nxy) if plus_nxy else float("nan")
        mm = np.mean(minus_nxy) if minus_nxy else float("nan")
        print(f"      {kind:4s}: <n_XY|+1>={mp:.3f} (n={len(plus_nxy)})   <n_XY|-1>={mm:.3f} (n={len(minus_nxy)})", flush=True)

    return verdicts


def main():
    print("=== SURVIVOR vs PARTICLE-HOLE SELF-MIRROR (gate-first; a firing gate is the find) ===", flush=True)
    J = 1.0
    # half-filling survivor lives BELOW the coherence horizon Q*(N)~0.59N: use strong dephasing (low Q).
    for g in (1.0, 0.5):                                    # Q = 1, 2  (below Q*(6)~3.5)
        run(6, J, g, "chain", "heisenberg")
    # XY control: U_PH is the clean E->-E symmetry of XY (no ZZ)
    run(6, J, 1.0, "chain", "xy")
    # counter-case: the star has NO interior half-filling survivor (G0 should fire gracefully)
    run(6, J, 1.0, "star", "heisenberg")
    # weak dephasing (above the horizon): survivor should LEAVE half-filling (G0 fires gracefully)
    run(6, J, 0.05, "chain", "heisenberg")
    print("\n=== done ===", flush=True)


if __name__ == "__main__":
    main()
