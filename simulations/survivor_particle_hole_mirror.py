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
      below the chain's own handover Q*_gap(N), G5).  If not, the conjecture's premise is void in this regime.
  G1  the INTERTWINER: for each operator S, how does the conjugation superoperator relate to L
      (commute / anticommute)?  The relationship sets the MEANING of G2.
  G2  the FIXED-POINT PROBE: is the survivor an eigenvector of S (parity +-1)?
      PREDICTION: X -> -1 (dark anti-fixed, control), R -> +1 (the real self-mirror),
      U -> the open question.  A clean separation X(-1)/R(+1) IS the find.
  G3  CAUSE vs CORRELATION: across the (p,p) spectrum, does S-parity sort DARK from BRIGHT
      (<n_XY> = -Re(lambda)/2gamma, the Absorption Theorem)?  Is light content == a definite parity?
  G4  THE LEVEL, not the solver's pick: X and R commute with L and with each other, so the half-filling block
      splits into four joint sectors (X = +-1, R = +-1); the slowest rate of each, with how many eigenvalues of
      that sector lie within the eigensolver's rounding (64 eps ||L||_1) of it, says which parities the block's
      slowest level holds and how often, and that mode's weight by the Hamming distance of its cells (the number
      of X/Y letters).  Then the survivor's populations (the diagonal of the slowest
      eigenvector of the first sector holding the level, scaled to 1 at its largest entry) and their Rayleigh
      quotient under the block's exclusion Laplacian, beside that Laplacian's smallest nonzero and largest
      eigenvalue: near the smallest for a density wave, at the largest for the mode a hop changes most; and their
      overlap with the block's diagonal of H, its mean removed, the populations of the energy mode; and the site
      profile n(j), the populations summed over the configurations with j occupied, scaled to 1 at its largest, with
      n . n_reflected / |n|^2 (at half filling n . n_complement / |n|^2 = -1 for every traceless population, whatever
      its shape, so only the reflection reads the shape).
  G5  the Heisenberg chain's handover Q*_gap(N) (N = 4, 6), read by bisection: the Q where the diagonal
      blocks' slowest <n_XY> reaches the band edge's 1, the edge of the regime this note reads.
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


def sector_slowest(L, S_list, signs, tol):
    """the slowest strictly-decaying eigenvalue of L on the joint eigenspace of the commuting involutions S_list (real
    permutation superoperators commuting with L) with the given signs, how many of that sector's eigenvalues lie within
    tol of it, and its eigenvector in the block's coordinates"""
    d2 = L.shape[0]
    P = np.eye(d2)
    for S, s in zip(S_list, signs):
        P = P @ (np.eye(d2) + s * S.real) / 2
    w, U = np.linalg.eigh((P + P.T) / 2)
    B = U[:, w > 0.5]                                       # an orthonormal basis of the projector's range
    if B.shape[1] == 0:
        return None, 0, None
    Ls = B.T @ L @ B
    ev, V = np.linalg.eig(Ls)
    cand = [k for k in range(len(ev)) if ev[k].real < -1e-9]
    if not cand:
        return None, 0, None
    k = max(cand, key=lambda k: ev[k].real)
    mult = sum(1 for e in ev if abs(e - ev[k]) <= tol)
    return ev[k], mult, B @ V[:, k]


def exclusion_laplacian(N, bs, bnds):
    """the configuration graph of the p-sector: one edge per hop of one excitation across one bond"""
    idx = {a: i for i, a in enumerate(bs)}
    Lex = np.zeros((len(bs), len(bs)))
    for i, a in enumerate(bs):
        for (u, v) in bnds:
            if ((a >> u) & 1) != ((a >> v) & 1):
                Lex[i, i] += 1
                Lex[idx[a ^ (1 << u) ^ (1 << v)], i] -= 1
    return Lex


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
        Lh, _, _ = block_L(N, half, half, J, g, bnds, model)
        lam_h, _ = slowest_eig(Lh)
        # a tie across fillings, to the eigensolver's rounding (eps ||L||_1 per block): on the XY chain the ladder puts
        # each block's spectrum into the next one's up to half filling, and every filling ties for the slowest rate
        if lam_h is not None and abs(lam_h.real - re_glob) <= 64 * np.finfo(float).eps * np.abs(Lh).sum(axis=0).max():
            print("    (every filling ties for the slowest rate here, to rounding: the survivor is not unique -- skipping)",
                  flush=True)
        else:
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
    print("G3  parity vs light content across the half-filling block (<n_XY> = -Re/2g):", flush=True)
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

    # G4 -- the level in the joint sectors of X and R, and the survivor's populations
    tol = 64 * np.finfo(float).eps * np.abs(L).sum(axis=0).max()
    print(f"G4  the block's slowest level by sector (rounding 64 eps ||L||_1 = {tol:.1e}):", flush=True)
    sectors = {}
    for sx in (+1, -1):
        for sr in (+1, -1):
            lam, mult, vec = sector_slowest(L, [ops["X"], ops["R"]], [sx, sr], tol)
            sectors[(sx, sr)] = (lam, mult, vec)
            txt = "no decaying mode" if lam is None else f"Re={lam.real:+.6f} |Im|={abs(lam.imag):.6f}, {mult} eigenvalue(s) there"
            print(f"      X={sx:+d} R={sr:+d}: {txt}", flush=True)
    top = max(lam.real for lam, _, _ in sectors.values() if lam is not None)
    holding = [key for key, (lam, _, _) in sectors.items() if lam is not None and lam.real >= top - tol]
    mult = sum(sectors[key][1] for key in holding)
    print(f"      the slowest level, Re={top:+.6f}: {mult}-fold, in the sector(s) "
          f"{', '.join(f'X={sx:+d} R={sr:+d}' for sx, sr in holding)}", flush=True)
    vec = sectors[holding[0]][2]
    pops = np.array([vec[i * d + i] for i in range(d)])
    pops = pops / pops[np.argmax(np.abs(pops))]
    Lex = exclusion_laplacian(N, bs, bnds)
    ew = np.linalg.eigvalsh(Lex)
    pr = pops.real
    rq = pr @ Lex @ pr / (pr @ pr)
    label = {a: "".join(str((a >> l) & 1) for l in range(N)) for a in bs}
    print(f"      populations of its slowest mode in X={holding[0][0]:+d} R={holding[0][1]:+d} (largest = 1; imaginary "
          f"part {np.abs(pops.imag).max():.1e}): " + ", ".join(f"{label[a]} {pops[i].real:+.3f}" for i, a in enumerate(bs)),
          flush=True)
    print(f"      their Rayleigh quotient under the block's exclusion Laplacian {rq:.4f}, its smallest nonzero eigenvalue "
          f"{min(e for e in ew if e > 1e-9):.4f} and largest {ew.max():.4f}", flush=True)
    hdist = np.array([bin(a ^ b).count("1") for a in bs for b in bs])
    wt = np.abs(vec) ** 2 / np.sum(np.abs(vec) ** 2)
    print("      its weight by the Hamming distance of the cells |a><b|, the number of X/Y letters: " +
          ", ".join(f"{h}: {wt[hdist == h].sum():.3g}" for h in sorted(set(hdist.tolist()))), flush=True)
    prof = np.array([sum(pr[i] for i, a in enumerate(bs) if (a >> j) & 1) for j in range(N)])
    if np.linalg.norm(prof) <= 1e-9 * np.linalg.norm(pr):
        # every site equally occupied: the populations read a correlation, not a density
        print(f"      its site profile n(j) vanishes (|n| / |populations| = {np.linalg.norm(prof) / np.linalg.norm(pr):.1e}): "
              f"no density wave", flush=True)
    else:
        prof = prof / prof[np.argmax(np.abs(prof))]
        print(f"      its site profile n(j), j = 0..{N - 1}: [" + ", ".join(f"{v:+.3f}" for v in prof) + "], n . n_reflected / |n|^2 = "
              f"{prof @ prof[::-1] / (prof @ prof):+.4f}", flush=True)
    hd = np.diag(sector_H(N, bs, J, bnds, model)).copy()
    hd -= hd.mean()
    if np.linalg.norm(hd) > 0:
        print(f"      their overlap with the block's diagonal of H, its mean removed (the energy mode's populations): "
              f"|cos| = {abs(pr @ hd) / (np.linalg.norm(pr) * np.linalg.norm(hd)):.4f}", flush=True)

    return verdicts


def main():
    print("=== SURVIVOR vs PARTICLE-HOLE SELF-MIRROR (gate-first; a firing gate is the find) ===", flush=True)
    J = 1.0
    # the half-filling survivor lives below the Heisenberg chain's own handover Q*_gap(N) (G5): strong dephasing.
    for g in (1.0, 0.5):                                    # Q = 1, 2  (below Q*_gap(6) = 2.448, G5)
        run(6, J, g, "chain", "heisenberg")
    # the chain N=4 at Q = 1: half filling, odd under both; the rings N=4, 6 at Q = 1: the density wave as a twofold
    # level holding both R-parities (the R printed is the solver's pick); at Q = 2 the chain N=4 and the ring N=6 hand
    # over to the (0,1) band edge, and the ring N=4's survivor is the energy mode's Neel pattern, even under both,
    # the fastest population mode at the Zeno end become the slowest
    for N, topo in ((4, "chain"), (4, "ring"), (6, "ring")):
        for g in (1.0, 0.5):
            run(N, J, g, topo, "heisenberg")
    # XY control: U_PH is the clean E->-E symmetry of XY (no ZZ)
    run(6, J, 1.0, "chain", "xy")
    # counter-case: the star has NO interior half-filling survivor (G0 should fire gracefully)
    run(6, J, 1.0, "star", "heisenberg")
    # weak dephasing (above the handover): survivor should LEAVE half-filling (G0 fires gracefully)
    run(6, J, 0.05, "chain", "heisenberg")
    # G5 -- the Heisenberg chain's own handover Q*_gap(N), where the diagonal blocks' slowest light content reaches the
    # band edge's 1; below it the half-filling wave is the survivor this note reads (not the XY chain's horizon Q*(N))
    for N in (4, 6):
        bnds = bonds(N, "chain")

        def kmin(Q):
            g = J / Q
            return min(-slowest_eig(block_L(N, p, p, J, g, bnds, "heisenberg")[0])[0].real / (2 * g) for p in range(1, N))
        lo, hi = 0.5, 4.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if kmin(mid) < 1 else (lo, mid)
        print(f"G5  Heisenberg chain N={N}: the diagonal blocks' slowest <n_XY> reaches the band edge's 1 at "
              f"Q*_gap = {0.5 * (lo + hi):.4f}", flush=True)
    print("\n=== done ===", flush=True)


if __name__ == "__main__":
    main()
