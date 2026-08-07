"""The sideways SPIN ladder S+ carries F125's cross-fold family, and it is one SU(2) multiplet.

S+ is the SPIN raising, not the eta one: PROOF_FROZEN_BAND_SO4 section 2 gives Phi the
eta-pairing and S+ the spin. An earlier draft of this work had it backwards.

Gates the multiplet reading of the `sideways_spin_ladder` open arc at N = 5 and N = 7: the
script exits non-zero if any check fails. It does NOT gate everything that arc quotes; the rung
sweep lives in eta_ladder_blocks.py and the break-input table in eta_ladder_breakinput.py, and
both of those print rather than assert.

What it checks, in the order that makes the later checks readable:

  CONTROL   Phi and S+ both intertwine L exactly on every rung of the chain, residual == 0.0.
            Phi's half is what both arcs prove, so a nonzero residual there means the build is
            wrong and nothing below may be read.

  SHAPE     The interior blocks carrying the cross-fold partner -lambda_A - 2N are exactly the
            interiors of the two S+ chains at p+q = N-1 and p+q = N+1.  Written as a PREDICTION
            derived at N=5 and first tested at N=7.

  COUNT     Those fold sectors plus the band sectors number 4N-8, which is the orbit size
            PROOF_CODIM1_BY_ADDITIVITY states independently for odd N.

  LADDER    The transport norms along a chain are the Clebsch-Gordan coefficients
            sqrt(l(l+1) - m(m+1)) for l = (N-3)/2, with no free parameter, and the chain
            terminates at the EIGENSOLVER FLOOR at the step into a boundary block, which is the
            highest-weight condition S+|l,l> = 0.

The seeds are the repo's own, from compute/RCPsiSquared.Core/F89PathK/RealDefectiveSeeds.cs;
gamma is uniform 1 and J_b = 2q*, the convention of PROOF_CODIM1_BY_ADDITIVITY section 1.

The recorded lambda_A is a 4-decimal value, so the acceptance tolerance is DERIVED from the run
(ten times the smallest offset actually found) rather than set: an earlier pass mistook that
rounding offset for a physical spread.

Usage:  python simulations/eta_ladder_chain.py [N ...]      (default: 5 7)
"""
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("\\", 1)[0].rsplit("/", 1)[0])
from eta_ladder_blocks import build_L, build_ladder, block_basis  # noqa: E402

# (q*, lambda_A) from RealDefectiveSeeds.All, the +1 R-parity entries used by the census probes.
SEEDS = {5: (0.620878, -4.6189), 7: (1.514833, -4.8846)}

failures = []


def check(ok, label, detail=""):
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    if not ok:
        failures.append(label)


def run(N):
    if N not in SEEDS:
        print(f"  N = {N}: no recorded defective seed (RealDefectiveSeeds lists odd N only)")
        return
    qstar, lam_a = SEEDS[N]
    fold = -lam_a - 2 * N
    J = [2 * qstar] * (N - 1)
    gamma = [1.0] * N
    interior = lambda p, q: 1 <= p <= N - 1 and 1 <= q <= N - 1
    chains = [N - 1, N + 1]

    print(f"\n  N = {N}, gamma = 1 uniform, q* = {qstar}, J_b = {2 * qstar}")
    print(f"    lambda_A = {lam_a}, cross-fold partner = {fold:.4f}, "
          f"centre {(lam_a + fold) / 2:.4f} = -N")

    # CONTROL: both ladders intertwine exactly, on every rung of both chains.
    worst = 0.0
    for kind in ("phi", "splus"):
        dq = 1 if kind == "phi" else -1
        for total in chains:
            for p in range(N + 1):
                q = total - p
                if not (0 <= q <= N and 0 <= p + 1 <= N and 0 <= q + dq <= N):
                    continue
                Ls, _, _ = build_L(N, p, q, J, gamma)
                Lt, _, _ = build_L(N, p + 1, q + dq, J, gamma)
                M = build_ladder(N, p, q, kind)
                if M.size:
                    worst = max(worst, np.linalg.norm(Lt @ M - M @ Ls))
    check(worst == 0.0, "CONTROL both ladders intertwine exactly", f"worst residual {worst:.1e}")

    # Spectra once per block, and a tolerance derived from the run.
    spec, best = {}, np.inf
    for p in range(N + 1):
        for q in range(N + 1):
            L, _, _ = build_L(N, p, q, J, gamma)
            if L.shape[0]:
                spec[(p, q)] = np.linalg.eigvals(L)
                best = min(best, np.abs(spec[(p, q)] - fold).min(),
                           np.abs(spec[(p, q)] - lam_a).min())
    tol = max(best * 10, 1e-9)

    foldset = sorted(k for k, e in spec.items()
                     if interior(*k) and np.abs(e - fold).min() < tol)
    bandset = sorted(k for k, e in spec.items()
                     if interior(*k) and np.abs(e - lam_a).min() < tol)
    pred = sorted([(p, s - p) for s in chains for p in range(N + 1)
                   if 0 <= s - p <= N and interior(p, s - p)])
    check(foldset == pred, "SHAPE fold sectors are the two chain interiors",
          f"{len(foldset)} sectors, tolerance {tol:.1e}")
    check(len(foldset) + len(bandset) == 4 * N - 8, "COUNT fold + band = 4N-8",
          f"{len(foldset)} + {len(bandset)} = {len(foldset) + len(bandset)}")

    # LADDER: transport norms against the Clebsch-Gordan coefficients, no free parameter.
    l = (N - 3) / 2.0
    cg = [np.sqrt(l * (l + 1) - m * (m + 1)) for m in np.arange(-l, l)]
    for total in chains:
        chain = [(p, total - p) for p in range(N + 1) if 0 <= total - p <= N]
        norms, terminal = [], None
        for (p, q) in chain[:-1]:
            Ls, _, _ = build_L(N, p, q, J, gamma)
            Lt, _, _ = build_L(N, p + 1, q - 1, J, gamma)
            if not Ls.shape[0] or not Lt.shape[0]:
                continue
            w, V = np.linalg.eig(Ls)
            j = int(np.argmin(np.abs(w - fold)))
            if abs(w[j] - fold) > tol:
                continue                      # source does not carry the fold
            v = V[:, j] / np.linalg.norm(V[:, j])
            Sblk = build_ladder(N, p, q, "splus")
            img = Sblk @ v
            n = float(np.linalg.norm(img))
            src = block_basis(N, p, q)[0]
            tgt = block_basis(N, p + 1, q - 1)[0]
            if interior(p + 1, q - 1):
                res = np.linalg.norm(Lt @ img - w[j] * img) / n
                check(res < 1e-12, f"    image is an eigenvector at the same lambda "
                                   f"({p},{q})->({p+1},{q-1})", f"residual {res:.1e}")
                norms.append(n)
            else:
                terminal = n
                term_Ls, term_S, term_tgt = Ls, Sblk, tgt
        check(len(norms) == len(cg) and all(abs(a - b) < 1e-6 for a, b in zip(norms, cg)),
              f"LADDER p+q={total} norms are the CG coefficients for l = {l:g}",
              f"{[f'{x:.6f}' for x in norms]} vs {[f'{x:.6f}' for x in cg]}")
        # NOT an exact zero, and it cannot be one.  The intertwining residual above IS exactly
        # 0.0 because it is an operator identity that cancels pairwise, but v here comes from an
        # eigensolver and carries its own error, so ||S+ v|| at the terminal step measures the
        # EIGENVECTOR's accuracy, not S+.  The error model is ||S+|| * ||v_computed - v_true||,
        # so the quantity with meaning is the RATIO to the interior steps, which is at the
        # floor and stays there across N: gate that, and print it so a reader sees the size.
        ratio = terminal / max(norms) if norms and terminal is not None else float("nan")
        # THE TERMINAL CHECK IS WEAK ON ITS OWN and must not be read as the evidence.  At the
        # last rung S+ maps a large block onto a small one, so most of the source block is in
        # its kernel and a randomly chosen eigenvector would "die" too.  Measured and printed
        # below: the surviving fraction.  What carries the reading is the CONTRAST, every
        # interior step alive at the CG norm and this one dead, plus the negative control that
        # survivors exist at all.
        allv = np.linalg.eig(term_Ls)[1]
        allv = allv / np.linalg.norm(allv, axis=0)
        alln = np.linalg.norm(term_S @ allv, axis=0)
        survivors = int((alln > 0.1).sum())
        check(survivors >= len(term_tgt), f"LADDER p+q={total} terminal rung has survivors "
                                     f"(negative control: the map is not the zero map)",
              f"{survivors} of {len(allv[0])} eigenvectors transport at norm > 0.1, "
              f"target dim {len(term_tgt)}")
        check(ratio < 1e-12, f"LADDER p+q={total} the fold vector is among the dying "
                             f"(highest weight S+|l,l> = 0)",
              f"||S+ v|| = {terminal:.2e} against interior {max(norms):.6f}, ratio {ratio:.1e}; "
              f"{len(allv[0]) - survivors}/{len(allv[0])} of that block dies, so this check alone is "
              f"weak and the contrast with the interior steps is what counts")


def sigma_min_is_not_the_cg_coefficient():
    """sigma_min of a rung map is NOT the multiplet's CG coefficient, and the difference matters.

    A review round proposed that F125's pinned sigma_min(W) values (1 at N=4, sqrt(2) at N=5)
    ARE the eta-side CG coefficients, i.e. a free two-point confirmation of the multiplet
    reading. They are not. sigma_min is the smallest singular value of the WHOLE rung map, and
    a block carries several multiplets, so it reads the weakest direction present; it happens
    to coincide where the block is small. Gated here because the open arc quotes the
    disagreement in order to REJECT that proposal, and a rejection needs a route back too.
    """
    print()
    print("  sigma_min(W) on an eta rung against the CG coefficient of the spin-(N-3)/2 chain")
    agree, differ = [], []
    for N in (4, 5, 7):
        l = (N - 3) / 2.0
        m = -l
        for p in range(1, N - 1):
            if not (1 <= p + 1 <= N - 1 and 1 <= p + 2 <= N - 1):
                continue
            sv = np.linalg.svd(build_ladder(N, p, p + 1, "phi"), compute_uv=False)
            cg = np.sqrt(l * (l + 1) - m * (m + 1))
            same = abs(sv[-1] - cg) < 1e-9
            (agree if same else differ).append((N, p, float(sv[-1]), float(cg)))
            print(f"    N={N} ({p},{p+1})->({p+1},{p+2}): sigma_min = {sv[-1]:.6f}  "
                  f"CG = {cg:.6f}  {'equal' if same else 'DIFFERENT'}")
            m += 1
    check(len(differ) > 0, "sigma_min is NOT the CG coefficient in general",
          f"{len(differ)} rung(s) differ, e.g. "
          + ", ".join(f"N={n} rung {p}: {s:.6f} vs {c:.6f}" for n, p, s, c in differ[:2]))
    # the last rung starts at p = N-3, not N-2: its target block is (N-2, N-1)
    check(all(N in (4, 5) or p in (1, N - 3) for (N, p, _, _) in agree),
          "the agreements are the small blocks and the outer rungs, not a structural match",
          f"{len(agree)} agreeing rung(s)")


if __name__ == "__main__":
    for N in [int(a) for a in sys.argv[1:]] or [5, 7]:
        run(N)

    sigma_min_is_not_the_cg_coefficient()

    print()
    if failures:
        print(f"FAILED {len(failures)} check(s): {failures}")
        sys.exit(1)
    print("all checks passed")
