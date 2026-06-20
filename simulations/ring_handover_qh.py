"""Gate-first: the RING half-filling double-excitation HANDOVER Q_h(N) -- the hard open piece of the
clock_hand_ladder arc (keeps CoherenceHorizonClaim + SecondClockRegimeClaim at Tier1Candidate).

THE OBJECT (operational, from handover_q.py, reframed via the Absorption Theorem):
  On the closed XY ring under Z-dephasing, below a crossover Q the longest-lived interior mode is a
  diagonal (k,k)-sector coherence (the "double-excitation seam"); above it, the single-excitation band
  edge takes over. The handover Q_h is where the (k,k) survivor's darkness <n_XY> reaches the F50 floor
  EXACTLY 1. By the Absorption Theorem Re(lam) = -2g*<n_XY>, so:

     Q_h(N) = the Q at which the SLOWEST interior mode of the (k,k) Liouvillian sector has Re = -2g.

  Below Q_h the (k,k) survivor is DARKER (<n_XY> < 1) and out-survives the band edge; above, it dips
  below the floor and the band edge (at exactly -2g) wins. NOT the second-clock floor-frequency
  sqrt(B^2-(2g)^2) (that is the on-floor {0,2} frequency at ring N=4 ONLY); this is the off-floor handover.

THE OPEN QUESTION THIS GATES: which sector k? Benzene N=6 survivor is the 2-excitation (2,2). But a FIXED
2-particle band saturates (4J cos(pi/N) -> 4J), so a linear Q_h~0.29N would suggest a GROWING half-filling
(N/2,N/2) sector. STAGE 1 computes BOTH k=2 and k=N/2 and lets the data decide which grows ~0.29N and
matches the full-L survivor.

  STAGE 0  GROUND TRUTH (full 4^N L, N=6): the global interior survivor below Q_h, its (ket#,bra#) sector
           and darkness; the (k,k)-block reproduces it. Plus Hk(N=4,k=2) = 2sqrt2 J (anti-periodic gate).
  STAGE 1  SECTOR DISCRIMINATOR: Q_h(N) for k=2 AND k=N/2 (even N); which is linear ~0.29N, which saturates.
  STAGE 2  THE LAW: c_eff = (N/Q_h)^2 -- flat (not sqrt-N)?  closed form (c_eff=12 -> Q_h=N/(2sqrt3))?

A firing gate is the find: if k=2 saturates while k=N/2 grows linearly, the survivor is half-filling, not
the fixed double-excitation -- diagnose, do not loosen.
"""
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from itertools import combinations
from math import sqrt, cos, pi

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

J = 1.0
TOL = 1e-7
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)


# ----------------------------------------------------------------- (k,k) sector machinery
def k_states(N, k):
    return [frozenset(c) for c in combinations(range(N), k)]


def Hk_ring(N, k):
    """XY hopping in the k-excitation sector: sign-free hardcore-boson hopping on the ring, wrap bond
    included. Sign-free in the spin model; for even k on a periodic ring this is the anti-periodic
    free-fermion sector (the 2sqrt2 source at N=4, k=2)."""
    states = k_states(N, k)
    idx = {s: a for a, s in enumerate(states)}
    D = len(states)
    H = np.zeros((D, D))
    bonds = [(i, (i + 1) % N) for i in range(N)]
    for a, s in enumerate(states):
        for (i, j) in bonds:
            si, sj = i in s, j in s
            if si ^ sj:
                moved = (s - {i}) | {j} if si else (s - {j}) | {i}
                H[idx[moved], a] += J
    return H, states


def hamming_pair(s1, s2):
    return len(s1 ^ s2)


def Lkk_sparse(N, k, gamma):
    H, states = Hk_ring(N, k)
    D = len(states)
    Hs = sp.csr_matrix(H)
    Id = sp.identity(D, format="csr")
    L = -1j * (sp.kron(Hs, Id) - sp.kron(Id, Hs.T))
    deph = np.array([hamming_pair(states[a], states[b]) for a in range(D) for b in range(D)], float)
    return (L - 2 * gamma * sp.diags(deph)).tocsc()


def darkness(N, k, Q):
    """<n_XY> = -Re/(2g) of the slowest interior (Re<0) mode of the (k,k) sector at J=1, gamma=1/Q.
    Dense eigvals for small blocks; sparse SHIFT-INVERT (sigma just right of 0 -> the rightmost = slowest
    modes) for large -- ARPACK 'LR' fails to converge on the near-anti-Hermitian L, shift-invert does not."""
    gamma = 1.0 / Q
    L = Lkk_sparse(N, k, gamma)
    Dsq = L.shape[0]
    if Dsq <= 6000:                                     # dense beats shift-invert overhead below ~6000
        ev = np.linalg.eigvals(L.toarray())
    else:
        kk = min(30, Dsq - 2)
        # sigma = +1e-3 is to the RIGHT of every eigenvalue (all Re <= 0); shift-invert 'LM' finds those
        # nearest sigma = the slowest (max-Re) modes. The handover mode is real (frozen) so this is exact there.
        ev = spla.eigs(L, k=kk, sigma=1e-3, which="LM", return_eigenvectors=False, maxiter=8000, tol=1e-10)
    interior = ev[ev.real < -1e-9]
    return -interior.real.max() / (2 * gamma) if interior.size else 0.0


def qh_sector(N, k, lo=1.2, hi=5.0):
    """Q where darkness(N,k,Q) = 1 (bisection; nan if no bracket in [lo,hi]). Coarse width break
    (~8 evals): Q_h to +-0.01 is plenty for the slope fit, and each darkness() eval is a full eig."""
    flo = darkness(N, k, lo) - 1
    fhi = darkness(N, k, hi) - 1
    if flo * fhi > 0:
        return float('nan')
    for _ in range(30):
        if hi - lo < 0.02:
            break
        mid = 0.5 * (lo + hi)
        fm = darkness(N, k, mid) - 1
        if abs(fm) < 1e-7:
            return mid
        if flo * fm < 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


# ----------------------------------------------------------------- full Liouvillian (ground truth, N<=6)
def site(op, l, N):
    out = np.array([[1.0 + 0j]])
    for kk in range(N):
        out = np.kron(out, op if kk == l else I2)
    return out


def full_L(N, gamma):
    d = 2 ** N
    H = np.zeros((d, d), complex)
    for i in range(N):
        j = (i + 1) % N
        H += (J / 2) * (site(X, i, N) @ site(X, j, N) + site(Y, i, N) @ site(Y, j, N))
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for l in range(N):
        Zl = site(Z, l, N)
        L += gamma * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    return L


def popcount(x):
    return bin(x).count("1")


def full_survivor(N, Q):
    """Global slowest interior mode of the full 4^N L: darkness, dominant (ket#,bra#) sector, its weight."""
    gamma = 1.0 / Q
    d = 2 ** N
    ev, evec = np.linalg.eig(full_L(N, gamma))
    interior = np.where(ev.real < -1e-6)[0]
    idx = interior[np.argmax(ev.real[interior])]
    v = evec[:, idx].reshape(d, d)
    w = np.abs(v) ** 2
    kn = np.array([popcount(i) for i in range(d)])
    sect = {}
    for i in range(d):
        row = w[i]
        for jj in range(d):
            if row[jj] > 1e-12:
                key = (kn[i], kn[jj])
                sect[key] = sect.get(key, 0.0) + row[jj]
    dom = max(sect, key=sect.get)
    return -ev[idx].real / (2 * gamma), dom, sect[dom] / sum(sect.values())


def main():
    GATE = {"fired": []}

    def gate(name, ok, detail=""):
        flag = "ok " if ok else "GATE-FIRE"
        if not ok:
            GATE["fired"].append(name)
        print(f"   [{flag}] {name}" + (f"   {detail}" if detail else ""))

    # ===================================================================== STAGE 0
    print("=" * 100)
    print("STAGE 0 -- GROUND TRUTH (full 4^N L, N=6) + the anti-periodic Hk gate")
    print("=" * 100)
    top4 = np.linalg.eigvalsh(Hk_ring(4, 2)[0]).max()
    gate("Hk(N=4,k=2) two-particle top = 2sqrt2 J (sign-free hardcore-boson = anti-periodic fermions)",
         abs(top4 - 2 * sqrt(2) * J) < 1e-9, f"top={top4:.6f}")

    # N=6 ground truth: find Q_h from the (2,2) block, confirm the full-L survivor there is (2,2) and matches
    qh26 = qh_sector(6, 2)
    Qchk = max(qh26 - 0.3, 1.0)
    dk_full, dom, frac = full_survivor(6, Qchk)
    dk_blk = darkness(6, 2, Qchk)
    print(f"   N=6 at Q={Qchk:.2f} (below Q_h~{qh26:.3f}): full-L survivor sector (ket#,bra#)={tuple(int(x) for x in dom)} "
          f"({frac*100:.0f}%), darkness={dk_full:.4f}; (2,2)-block darkness={dk_blk:.4f}")
    # the survivor is the 2-EXCITATION DOUBLET {(2,2),(N-2,N-2)} -- particle-hole partners, isospectral, so
    # the (2,2) block computes its darkness/Q_h exactly even when the full-L dominant component is (N-2,N-2).
    gate("N=6: the full-L survivor is the 2-excitation doublet {(2,2),(4,4)} (particle-hole partners)",
         dom in {(2, 2), (4, 4)}, f"sector {tuple(int(x) for x in dom)}")
    gate("N=6: the (2,2)-block darkness equals the full-L survivor darkness (PH-isospectral, the right object)",
         abs(dk_full - dk_blk) < 1e-3, f"full {dk_full:.5f} vs block {dk_blk:.5f}")

    # ===================================================================== STAGE 1: sector discriminator
    print("\n" + "=" * 100)
    print("STAGE 1 -- SECTOR DISCRIMINATOR: Q_h for k=2 (fixed double-excitation) vs k=N/2 (half-filling)")
    print("=" * 100)
    print(f"{'N':>3} {'Q_h(k=2)':>10} {'slope/N':>9} {'Q_h(k=N/2)':>12} {'slope/N':>9}  (k=2 matches full-L at N=6)")
    qh2, qhh = {}, {}
    for N in (6, 8, 10, 12):
        q2 = qh_sector(N, 2)
        qh2[N] = q2
        # half-filling only as the small-N discriminator (N=6,8): at N=6, k=2 matches the full-L Q_h while
        # k=N/2 does NOT, so the survivor is the 2-excitation doublet, NOT half-filling. Skip the costly N>=10.
        if N <= 8:
            qh = qh_sector(N, N // 2)
            qhh[N] = qh
            hh = f"{qh:>12.5f} {qh/N:>9.5f}"
        else:
            hh = f"{'(skip)':>12} {'':>9}"
        print(f"{N:>3} {q2:>10.5f} {q2/N:>9.5f}  {hh}")

    # ===================================================================== STAGE 2: the law (on the winning sector)
    print("\n" + "=" * 100)
    print("STAGE 2 -- THE LAW on the linear sector: c_eff=(N/Q_h)^2 flat? closed form (c_eff=12 -> N/(2sqrt3))?")
    print("=" * 100)

    def report(label, qh):
        Ns = [N for N in sorted(qh) if np.isfinite(qh[N])]
        print(f"\n   [{label}]")
        print(f"   {'N':>4} {'Q_h':>10} {'Q_h/N':>10} {'c_eff=(N/Q_h)^2':>16} {'N/(2sqrt3)':>12}")
        ce = []
        for N in Ns:
            q = qh[N]
            ce.append((N, (N / q) ** 2))
            print(f"   {N:>4} {q:>10.5f} {q/N:>10.5f} {(N/q)**2:>16.4f} {N/(2*sqrt(3)):>12.5f}")
        if len(ce) >= 2:
            flat = abs(ce[-1][1] - ce[-2][1]) / ce[-2][1] < 0.05
            print(f"   c_eff flat (last two within 5%): {flat}  "
                  f"[c_eff({ce[-2][0]})={ce[-2][1]:.3f}, c_eff({ce[-1][0]})={ce[-1][1]:.3f}]")
            return Ns, flat, ce
        return Ns, False, ce

    Ns2, flat2, ce2 = report("k = 2 (fixed double-excitation)", qh2)
    if qhh:
        report("k = N/2 (half-filling)", qhh)

    # decide the linear sector and fit the slope
    big = [(N, qh2[N]) for N in Ns2 if N >= 8 and np.isfinite(qh2[N])]
    if len(big) >= 2:
        A = np.array([[N, 1.0] for N, _ in big]); b = np.array([q for _, q in big])
        slope, intc = np.linalg.lstsq(A, b, rcond=None)[0]
        print(f"\n   k=2 large-N fit (N>=8): Q_h = {slope:.5f}*N + {intc:.5f}   "
              f"(1/(2sqrt3)={1/(2*sqrt(3)):.5f}, 1/pi={1/pi:.5f})")

    print("\n" + "=" * 100)
    print(f"GATES: {'ALL PASS' if not GATE['fired'] else str(len(GATE['fired'])) + ' FIRED -> ' + str(GATE['fired'])}")
    print("Read the Stage-1 table: which sector's Q_h grows linearly ~0.29N (the real survivor) vs saturates.")
    print("=" * 100)
    if GATE["fired"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
