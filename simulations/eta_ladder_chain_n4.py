"""The N=4 sideways-ladder walk: band and fold freight share the four sectors. GATED.

Sibling of eta_ladder_chain.py (the odd-N gate, N=5/7); the even-N case the arc left open,
run 2026-08-10 at the four real-q defective loci of the N=4 (1,2) block (ReferenceDefectiveLoci;
their lambda measured the same day, all on Re lambda = -4, the fold-partner map's fixed line).

What N=4 changes, measured here and gated:

  THE SAME FREIGHT, FOLDED INTO SHARED SECTORS. The walk's freight is, as at odd N, the
  band value lambda_A plus its fold partner s(lambda) = -lambda - 2N, the holomorphic fold
  of PROOF_CODIM1_BY_ADDITIVITY section 7 consequence (b). The loci sit on Re lambda = -4,
  which is the fixed LINE of the antiholomorphic partner map -conj(lambda) - 2N (F125's
  fold partner in the F-registry) -- and that placement is DERIVED, not accidental: both s
  and conjugation (T is diagonal, so the F89_SEED_EXISTENCE_REDUCTION bipartite similarity
  runs in every (p,q) block at every N) permute the block's defective set, so a SINGLE
  defective pair must be fixed by their composite, forcing Re lambda = -4 (gated below).
  There the fold value -lambda* - 2N equals conj(lambda*): band and fold freight are a
  conjugate pair, and foldset = bandset = the confined orbit {(1,2),(2,1),(2,3),(3,2)} by
  the walk's own definitions. Measured: EACH S+ chain transports BOTH values, so four
  l = 1/2 doublets in all (two S+ chains x two conjugate values) live on the 4N-12 = 4
  shared sectors. The odd-N seat identity 4(N-2) = |chains| x interior does NOT
  re-instantiate here: at N=4 the two eta/Phi chain interiors are the SAME four blocks
  ({(1,2),(2,3)} and {(2,1),(3,2)}), the four families of the odd-N accounting landing on
  one block set; the eta side is, as at odd N, not separately norm-measured. What even N
  adds is that its loci are complex-lambda (the real-lambda species is excluded there, the
  even-N check), so an even-N walk always carries the conjugate pair split open, while at
  the real odd-N seeds the pair is invisible, not absent. The block-level sharing itself
  (overlap 4, band-and-fold double role) is section 7's already; what this run adds is the
  VALUE-level pair in every shared sector, and its transport.

  Per doublet, the walk of eta_ladder_chain.py reproduces exactly: the interior rung
  transports at the CG norm sqrt(l(l+1) - m(m+1)) = 1 (measured 1.000000000; the split
  pair's two members and their mix agree, though near the EP the members are almost
  parallel, 1 - |overlap| between 2.6e-7 and 3.0e-6 across the loci, printed below, so
  that agreement is a smoke check here, not the N=9 mix-robustness); the image is an eigenvector of the target at the same value (relative
  residual ~1e-14); the terminal step into the boundary block dies at the eigensolver floor
  (ratio ~1e-15, judged as a ratio, never an absolute); survivors = 4 = dim(target); and
  the Phi + S+ intertwining CONTROL is exactly 0.0 on every rung of both chains (operator
  identity, pairwise cancellation, the same == 0.0 gate as the odd-N script).

  A counter-reading died in this run (recorded here, its only record): "conjugation splits
  the chain" assumed the orbit's conjugate-relatedness separates the two values across
  blocks; it forgot the block-local T-conjugacy that closes every block spectrum under
  conjugation, so both values live in every orbit block and nothing splits.

  Labels inherited, not re-measured: "defective" for the four loci is ReferenceDefectiveLoci's
  EpCharacter classification; this script gates walk properties only (no gap exponent, no
  loop). The scan behind that list shows ONE member per locus because it reads the R-even
  sector; the conjugate twin lives in R-odd (sigma_even = conj sigma_odd at even N, the
  sector-swap witness), while this walk's full blocks see both.

Caveat carried from the sigma_min fence: the CG value 1 at l = 1/2 is the value sigma_min
already coincides with at N=4, so this run's confirmation weight lives in the STRUCTURE
(both values in every sector, image-eigenvector, terminal death, survivors = dim(target)), not
in the single norm.

Run:  python simulations/eta_ladder_chain_n4.py     (gates everything; exit != 0 on failure)
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from eta_ladder_blocks import build_L, build_ladder

N = 4
GAMMA = [1.0] * N
ORBIT = [(1, 2), (2, 1), (2, 3), (3, 2)]
SEEDS = [0.460212, 0.854438, 0.857458, 1.738181]   # octic-book q; J = 2 q* (the repo convention)

failures = 0


def check(ok, label, detail=""):
    global failures
    print(("PASS " if ok else "FAIL ") + label + (f"  [{detail}]" if detail else ""))
    if not ok:
        failures += 1


for qstar in SEEDS:
    J = [2.0 * qstar] * (N - 1)
    print(f"\n=== seed q* = {qstar} ===")

    L12, _, _ = build_L(N, 1, 2, J, GAMMA)
    w12 = np.linalg.eigvals(L12)
    dist = np.abs(w12[:, None] - w12[None, :]) + np.eye(len(w12)) * 1e9
    i, j = np.unravel_index(np.argmin(dist), dist.shape)
    lam = (w12[i] + w12[j]) / 2
    half = dist[i, j] / 2
    tol = max(10 * half, 1e-9)
    print(f"lambda* = {lam:.6f}  (pair split {2*half:.2e}, tol {tol:.2e})")

    # the self-fold pin: the locus value sits on Re lambda = -4 = -N (the fixed line of
    # lambda -> -conj(lambda) - 2N), the mean of the split pair recovering it. Error model:
    # near a defective pair the eigenvalue error is O(sqrt(eps)) ~ 1e-8; the gate is 100x that.
    check(abs(lam.real + 4.0) < 1e-6, "self-fold pin Re lambda* = -4", f"{lam.real}")

    # BOTH VALUES: every orbit block carries lambda* AND its fold partner conj(lambda*)
    blocks = {}
    for (p, q) in ORBIT:
        L, _, _ = build_L(N, p, q, J, GAMMA)
        w, V = np.linalg.eig(L)
        blocks[(p, q)] = (L, w, V)
        d1 = np.min(np.abs(w - lam))
        d2 = np.min(np.abs(w - np.conj(lam)))
        check(d1 < tol and d2 < tol, f"block {p},{q} carries lambda* AND conj",
              f"d={d1:.2e}/{d2:.2e}")

    # the four doublets: two chains x two conjugate values
    for (b0, b1, b2) in (((1, 2), (2, 1), (3, 0)), ((2, 3), (3, 2), (4, 1))):
        for tag, val in (("lam*", lam), ("conj", np.conj(lam))):
            _, ws, Vs = blocks[b0]
            Lt, wt, Vt = blocks[b1]
            S01 = build_ladder(N, b0[0], b0[1], "splus")
            order = np.argsort(np.abs(ws - val))
            k1, k2 = int(order[0]), int(order[1])
            name = f"{b0}->{b1} {tag}"
            v1n = Vs[:, k1] / np.linalg.norm(Vs[:, k1])
            v2n = Vs[:, k2] / np.linalg.norm(Vs[:, k2])
            print(f"      {name}: pair-member 1-|overlap| = {1 - abs(np.vdot(v1n, v2n)):.2e}")
            # CG norm 1 for members and mix (the closed form, no free parameter)
            for vtag, v in (("m1", Vs[:, k1]), ("m2", Vs[:, k2]), ("mix", Vs[:, k1] + Vs[:, k2])):
                n_ = np.linalg.norm(S01 @ (v / np.linalg.norm(v)))
                check(abs(n_ - 1.0) < 1e-6, f"{name} {vtag} CG norm 1", f"{n_:.9f}")
            # image is an eigenvector of the target at the same value
            v1 = Vs[:, k1] / np.linalg.norm(Vs[:, k1])
            img = S01 @ v1
            res = np.linalg.norm(Lt @ img - ws[k1] * img) / np.linalg.norm(img)
            check(res < 1e-12, f"{name} image-eigenvector residual", f"{res:.2e}")
            # terminal: highest-weight death at the floor, as a ratio; survivors saturate rank
            S12 = build_ladder(N, b1[0], b1[1], "splus")
            img_n = img / np.linalg.norm(img)
            norms_all = np.linalg.norm(S12 @ (Vt / np.linalg.norm(Vt, axis=0)), axis=0)
            ratio = np.linalg.norm(S12 @ img_n) / norms_all.max()
            from math import comb
            dim_tgt = comb(N, b2[0]) * comb(N, b2[1])   # (3,0) and (4,1) are both 4-dimensional
            check(ratio < 1e-12, f"{b1}->{b2} {tag} terminal ratio", f"{ratio:.2e}")
            check(int(np.sum(norms_all > 0.1)) == dim_tgt,
                  f"{b1}->{b2} {tag} survivors = dim(target) = {dim_tgt}",
                  f"{int(np.sum(norms_all > 0.1))}")

    # CONTROL: the intertwining residual is EXACTLY 0.0, both ladders, every rung, both chains
    worst = 0.0
    for total in (N - 1, N + 1):
        for p in range(0, N):
            q = total - p
            for kind, dq in (("phi", 1), ("splus", -1)):
                q2 = q + dq
                if not (0 <= q <= N and 0 <= q2 <= N and p + 1 <= N):
                    continue
                Ls, _, _ = build_L(N, p, q, J, GAMMA)
                Lt2, _, _ = build_L(N, p + 1, q2, J, GAMMA)
                M = build_ladder(N, p, q, kind)
                worst = max(worst, np.abs(Lt2 @ M - M @ Ls).max())
    check(worst == 0.0, "CONTROL both ladders intertwine exactly", f"{worst}")

print(f"\n{'ALL GATES PASS' if failures == 0 else f'{failures} FAILURE(S)'}")
sys.exit(0 if failures == 0 else 1)
