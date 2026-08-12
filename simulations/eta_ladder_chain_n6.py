"""The even-N chain walk from a (1,3) seed at N=6 -- GATED.

The arc sideways_spin_ladder's next step after the (1,3)@N=6 census: the census
restored the even-N walk input one block over. From (1,3) the S+ ladder reads
into (2,2) (p+q = 4 preserved) and Phi into (2,4) (p-q = -2 preserved). Neither
chain is one of the four canonical odd-N chains (p+q = N-/+1, d = +-1), and the
Phi ladder had never carried freight anywhere in the arc before this run.

Chains through (1,3)@N=6 (block dims in brackets):
  S+  : (0,4)[15] -> (1,3)[120] -> (2,2)[225] -> (3,1)[120] -> (4,0)[15]
  Phi : (0,2)[15] -> (1,3)[120] -> (2,4)[225] -> (3,5)[120] -> (4,6)[15]

WHAT THE LADDER OWNS (coupling-independent, gated once at an off-locus
non-uniform control): each ladder on (1,3) has EXACTLY TWO singular values,
sqrt2 on a 105-dim family and sqrt6 on a 15-dim family; the terminal rung has
the single value 2.0 at rank 15. Under each chain's su(2) the block decomposes
as 105 x (l=1, seat m=-1) + 15 x (l=2, seat m=-1): the su(2) weights
(p-q)/2 and (p+q-N)/2 fix m = -1, never l. So the CG rung norm sqrt2, the
highest-weight death and survivors = 15 are FAMILY properties: 105 of 120
eigenvectors pass the full walk gate set at ANY coupling (gated 105/120 at the
control). This is the N=4 sigma_min caveat's stronger sibling: the
confirmation weight lives in the structure, not in the norm. Each rate
value's exact 3-dim kappa=0 AT eigenspace sits inside one chain's 15-dim
l=2 family (rides sqrt6, sqrt6 over the walked rungs, no terminal death) and
inside the other chain's l=1 family (rides sqrt2, dies); which chain hosts
which rate value is gated in all four cases.

WHAT THE SEEDS OWN (the per-locus gated content): the pair at the census
lambda with truncation-sized split (two-sided gate, error model below); the
fold-line / axis pins (scale-relative, eps * ||L||); the twin s(lambda*) =
-2N - lambda* = -12 - lambda* present as a SECOND near-defective pair at the
same coupling (conj on the line, the distinct s-partner on the axis); BOTH
pairs defective-grade non-normal (the departure from normality s of the
restricted 2x2 exceeds 100x the split; measured 1.4e3..2.5e4x, against ~0 for
a semisimple degeneracy: this is the gate that a mere near-crossing would
fail); both freights lowest-weight in the l=1 family (entering-rung projector
at the floor, the AT strands at ~1 where they carry l=2); presence of both
freight values in every seat spectrum; member images exact eigenvectors of
every seat (Rayleigh floor). Mix vectors ride the same CG norm (in-plane
robustness); their Rayleigh residual is NOT gated against the split, which is
no bound for a non-normal pair (the true in-plane worst is s itself).

SPLIT ERROR MODEL (the two-sided gate): at a sqrt-branch locus the pair split
is 2*sqrt(|a|*|dq|) with |dq| <= 1e-7 (census q* recorded to 7 decimals gives
|dq| <= 5e-8; the two Sturm-only q* are isolated to 1e-7) and a the branch
coefficient; measured splits span 1.8e-5..1.6e-3. The gate is
1e-8 < split < 1e-2: the lower edge excludes the exact AT degeneracies
(split ~1e-15, seven orders below the edge), the upper edge excludes
non-locus pairs (O(0.1..1) spacing, one order above the edge). Not a bare
threshold: both edges sit >= 6x from every measured value.

Seeds: the 13 pinned (q*, lambda) of simulations/disc13_n6_census.py plus the
TWO loci the exact Sturm count added past the census floors (ef8f77d, doc
section "The exact real-root count" in experiments/F89_PATH_K_DIABOLIC.md):
q* = 0.6133162 R-even fold (inside the window, both census channels blind) and
q* = 0.4627729 R-odd fold (7.8e-4 beside the axis locus 0.4635508, below
channel A's grid step; the mixed-species close pair). Their lambda were read
with the census's sector + pair-filter device at the Sturm-exact q* (NOT the
census's detection chain: no trisection refinement, no exponent
classification), and that read is GATED below (sturm_only_sector_gate): each
locus's own R-sector holds exactly two truncation-split pairs, the
conjugate-doubled fold pair on Re lambda = -6, the +Im member at the SEEDS
lambda. The exact-device certificate |M*| = 4y^2 via M* = -s10/psc1 at the two
new A-roots stays open (psc1/s10 are landed in the slow C# pipeline, not
exported).
FIFTEEN seeds total (J = 2 q*, gamma = 1/site). Runtime ~1 min. Exit 0 iff
every gate passes.
"""
import sys, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eta_ladder_blocks import build_L, build_ladder

N = 6
GAMMA = [1.0] * N
SQRT2 = np.sqrt(2.0)
SQRT6 = np.sqrt(6.0)
EPS = np.finfo(float).eps

failures = 0


def check(ok, label, detail=""):
    global failures
    tag = "PASS" if ok else "FAIL"
    if not ok:
        failures += 1
    print(f"  {tag} {label}" + (f"  [{detail}]" if detail else ""))


# (q*, species, sector, lambda_census) from the census table
# (experiments/F89_PATH_K_DIABOLIC.md). species 'fold' = fold pair on
# Re lambda = -6, conj-doubled (lambda = -6 + y i, twin = conj);
# species 'axis' = conjugate pair on the real axis, s-doubled, PT-breaking
# (lambda real, twin = -12 - lambda is the other printed member).
SEEDS = [
    (0.4630013, "fold", "even", -6.0 + 0.3864821j),
    (0.6133162, "fold", "even", -6.0 + 1.7171734j),   # Sturm-only: both census channels blind
    (0.6400121, "fold", "even", -6.0 + 1.2687899j),
    (0.7026159, "axis", "even", -6.5624176 + 0.0j),
    (0.7450837, "fold", "even", -6.0 + 1.4621529j),
    (0.9184132, "fold", "even", -6.0 + 7.7721940j),
    (1.5109958, "fold", "even", -6.0 + 0.9288088j),
    (3.0922209, "fold", "even", -6.0 + 1.5169893j),
    (3.3536850, "axis", "even", -6.0787877 + 0.0j),
    (0.4627729, "fold", "odd", -6.0 + 3.0292680j),    # Sturm-only: below channel A's grid step,
    (0.4635508, "axis", "odd", -6.2680671 + 0.0j),    # ... the mixed-species close pair with this axis locus
    (0.8029541, "axis", "odd", -6.6051577 + 0.0j),
    (1.2273677, "axis", "odd", -6.3760926 + 0.0j),
    (1.9567664, "fold", "odd", -6.0 + 7.6439946j),
    (3.2608458, "fold", "odd", -6.0 + 25.1047962j),
]

# the two loci the exact Sturm count added past the census floors: their
# lambda is not a census pin, so the sector-level read is gated here
STURM_ONLY = {0.6133162, 0.4627729}

SPLUS_CHAIN = [(1, 3), (2, 2), (3, 1)]   # then boundary (4, 0)
PHI_CHAIN = [(1, 3), (2, 4), (3, 5)]     # then boundary (4, 6)
BOUNDARY = {"splus": (4, 0), "phi": (4, 6)}
ENTERING = {"splus": (0, 4), "phi": (0, 2)}   # the rung INTO the seed block


def eig_block(p, q, J):
    L, _, _ = build_L(N, p, q, J, GAMMA)
    w, V = np.linalg.eig(L)
    return L, w, V


def refl_site(x):
    y = 0
    for l in range(N):
        if x & (1 << (N - 1 - l)):
            y |= 1 << l
    return y


def sector_eigs(L13, sector):
    """Eigenvalues of the R-sector of the given (1,3) block (the census's
    sector device, disc13_n6_census.py: reflection permutation, eigh-split).
    Gates [L, P] == 0.0 first: B.T L B is a RESTRICTION only if L commutes
    with the reflection; without that gate a broken build would silently
    project and the downstream reads would be fabricated."""
    from eta_ladder_blocks import block_basis
    pairs, index = block_basis(N, 1, 3)
    P = np.zeros((len(pairs), len(pairs)))
    for col, (a, b) in enumerate(pairs):
        P[index[(refl_site(a), refl_site(b))], col] = 1.0
    check(np.abs(L13 @ P - P @ L13).max() == 0.0,
          "[L, P] == 0.0 exactly (the sector read restricts, not projects)")
    wP, VP = np.linalg.eigh((P + P.T) / 2)
    B = VP[:, wP > 0.5] if sector == "even" else VP[:, wP < -0.5]
    return np.linalg.eigvals(B.T @ L13 @ B)


def sturm_only_sector_gate(L13, sector, lam_census):
    """The lambda read of a Sturm-only locus, gated at sector level: the
    locus's own R-sector holds EXACTLY TWO truncation-split pairs, they are
    the conjugate-doubled fold pair on Re lambda = -6, and the +Im member
    sits at the SEEDS lambda. FOLD-ONLY by construction (pins Re = -6); an
    axis member of STURM_ONLY would need its own pin."""
    ws = sector_eigs(L13, sector)
    found = []
    for i in range(len(ws)):
        for j in range(i + 1, len(ws)):
            d = abs(ws[i] - ws[j])
            if 1e-8 < d < 1e-2:
                found.append((ws[i] + ws[j]) / 2)
    check(len(found) == 2,
          f"sector holds exactly two truncation-split pairs (got {len(found)})")
    if len(found) == 2:
        m_plus = max(found, key=lambda m: m.imag)
        m_minus = min(found, key=lambda m: m.imag)
        scale = np.abs(L13).sum(axis=1).max()
        check(abs(m_plus - np.conj(m_minus)) < 100 * EPS * scale,
              f"the two sector pairs are conjugates "
              f"(dev {abs(m_plus - np.conj(m_minus)):.1e})")
        check(abs(m_plus.real + 6.0) < 100 * EPS * scale,
              f"sector pair means on the fold line Re = -6 "
              f"({abs(m_plus.real + 6.0):.1e})")
        check(abs(m_plus - lam_census) < 5e-6,
              f"sector +Im pair mean at the SEEDS lambda "
              f"(dist {abs(m_plus - lam_census):.1e})")


def entering_image_basis(kind):
    """Orthonormal basis of the entering rung's image inside (1,3)."""
    p, q = ENTERING[kind]
    M = build_ladder(N, p, q, kind)
    Q, _ = np.linalg.qr(M)
    return Q[:, :np.linalg.matrix_rank(M)]


def pair_data(w, V, val):
    """The two eigenvalues closest to val: indices, mean, split, and the
    departure from normality s of the restricted 2x2 (s^2 = ||A||_F^2 -
    |l1|^2 - |l2|^2 for A the compression of L onto the pair's orthonormalized
    plane; s ~ 0 for a semisimple degeneracy, s = O(1) for a defective pair)."""
    order = np.argsort(np.abs(w - val))
    k1, k2 = order[0], order[1]
    return k1, k2, (w[k1] + w[k2]) / 2, abs(w[k1] - w[k2])


def departure_from_normality(L, V, k1, k2, w):
    # QR of two near-parallel eigenvectors is the ill-conditioned case; the
    # measured s is stable across both chains and carries a >= 14x margin
    # over the 100x gate, so the verdict is safe. A Schur-ordered invariant
    # subspace would be the fully robust route if the margin ever thins.
    Q, _ = np.linalg.qr(np.column_stack([V[:, k1], V[:, k2]]))
    A = Q.conj().T @ L @ Q
    s2 = np.linalg.norm(A, "fro") ** 2 - abs(w[k1]) ** 2 - abs(w[k2]) ** 2
    return np.sqrt(max(s2, 0.0))


def walk_gate_set(v, kind, chain, blocks, cg_value=SQRT2):
    """Run one vector up the chain; return (rungs_ok, term_ratio)."""
    src = chain[0]
    for tgt in chain[1:]:
        M = build_ladder(N, src[0], src[1], kind)
        img = M @ v
        nrm = np.linalg.norm(img)
        if abs(nrm - cg_value) > 1e-9:
            return False, None
        v = img / nrm
        src = tgt
    M_term = build_ladder(N, src[0], src[1], kind)
    term = np.linalg.norm(M_term @ v)
    _, ws, Vs = blocks[src]
    Vn = Vs / np.linalg.norm(Vs, axis=0, keepdims=True)
    return True, term / np.linalg.norm(M_term @ Vn, axis=0).max()


def walk_chain(chain, kind, val, blocks):
    """Transport the two closest eigenvectors to val (and their mix) up the
    chain. Members are eigenvectors, transported exactly (Rayleigh floor);
    the mix is kept for the CG in-plane robustness only."""
    p0, q0 = chain[0]
    L0, w0, V0 = blocks[(p0, q0)]
    k1, k2, _, _ = pair_data(w0, V0, val)
    mix = V0[:, k1] + V0[:, k2]
    vecs = {
        "m1": V0[:, k1] / np.linalg.norm(V0[:, k1]),
        "m2": V0[:, k2] / np.linalg.norm(V0[:, k2]),
        "mix": mix / np.linalg.norm(mix),
    }
    worst_cg = 0.0
    worst_res_member = 0.0
    worst_term = 0.0
    survivors_seen = set()
    presence = 0.0   # worst over seats of min distance of val to the spectrum
    for tgt in chain[1:]:
        presence = max(presence, np.min(np.abs(blocks[tgt][1] - val)))
    for name, v in vecs.items():
        v = v.copy()
        src = chain[0]
        for tgt in chain[1:]:
            M = build_ladder(N, src[0], src[1], kind)
            img = M @ v
            nrm = np.linalg.norm(img)
            worst_cg = max(worst_cg, abs(nrm - SQRT2))
            img /= nrm
            Lt = blocks[tgt][0]
            Limg = Lt @ img
            theta = np.vdot(img, Limg)
            if name != "mix":
                worst_res_member = max(
                    worst_res_member, np.linalg.norm(Limg - theta * img))
            v = img
            src = tgt
        M_term = build_ladder(N, src[0], src[1], kind)
        term_norm = np.linalg.norm(M_term @ v)
        _, ws, Vs = blocks[src]
        Vn = Vs / np.linalg.norm(Vs, axis=0, keepdims=True)
        all_norms = np.linalg.norm(M_term @ Vn, axis=0)
        worst_term = max(worst_term, term_norm / all_norms.max())
        survivors_seen.add(int(np.sum(all_norms > 0.1)))
    return (presence, worst_cg, worst_res_member, worst_term, survivors_seen)


def control(J):
    """Intertwining residual over every rung of both chains, exact-zero."""
    worst = 0.0
    for chain, kind in ((SPLUS_CHAIN, "splus"), (PHI_CHAIN, "phi")):
        seats = chain + [BOUNDARY[kind]]
        for src, tgt in zip(seats[:-1], seats[1:]):
            Ls, _, _ = build_L(N, src[0], src[1], J, GAMMA)
            Lt, _, _ = build_L(N, tgt[0], tgt[1], J, GAMMA)
            M = build_ladder(N, src[0], src[1], kind)
            worst = max(worst, np.abs(Lt @ M - M @ Ls).max())
    return worst


def at_control(blocks, ein_bases):
    """The kappa=0 AT eigenspaces at the two rate values, ALL FOUR cases.

    Each rate value's 3-dim exact eigenspace sits in the l=2 family of ONE
    chain (rides sqrt6, survives) and in the l=1 family of the OTHER (rides
    sqrt2, dies): the death sorts by multiplet family membership, not by
    defectiveness and not by rate value. Every member of each eigenspace is
    walked; the entering-rung projection ~1 marks the AT strands as NOT
    lowest-weight where they carry l=2.
    """
    for at_val, six_kind in ((-8.0, "splus"), (-4.0, "phi")):
        two_kind = "phi" if six_kind == "splus" else "splus"
        _, w0, V0 = blocks[(1, 3)]
        idx = [k for k in range(len(w0)) if abs(w0[k] - at_val) < 1e-8]
        check(len(idx) == 3, f"AT lambda={at_val}: exact eigenspace dim 3")
        for kind, cg, dies in ((six_kind, SQRT6, False), (two_kind, SQRT2, True)):
            chain = SPLUS_CHAIN if kind == "splus" else PHI_CHAIN
            ok_norms, ok_death = True, True
            for k in idx:
                v = V0[:, k] / np.linalg.norm(V0[:, k])
                rungs_ok, ratio = walk_gate_set(v, kind, chain, blocks, cg)
                ok_norms = ok_norms and rungs_ok
                if ratio is None:
                    ok_death = False
                elif dies:
                    ok_death = ok_death and ratio < 1e-12
                else:
                    ok_death = ok_death and ratio > 0.5
            check(ok_norms,
                  f"AT lambda={at_val} rides {kind} at "
                  f"{'sqrt6' if cg == SQRT6 else 'sqrt2'} (all 3 members)")
            check(ok_death,
                  f"AT lambda={at_val} on {kind}: "
                  f"{'dies at floor' if dies else 'no death'} (l="
                  f"{'1' if dies else '2'} family)")
        proj = max(np.linalg.norm(ein_bases[six_kind].conj().T @
                                  (V0[:, k] / np.linalg.norm(V0[:, k])))
                   for k in idx)
        check(proj > 0.9,
              f"AT lambda={at_val} NOT lowest-weight on {six_kind} "
              f"(entering-rung projection {proj:.3f})")


def structure_control():
    """The ladder's own coupling-independent structure, gated once at an
    off-locus non-uniform coupling: two singular values on (1,3), the
    terminal rank, and the 105/120 pass rate of the full walk gate set."""
    print("\n== structure control (off-locus, non-uniform J) ==")
    for kind in ("splus", "phi"):
        M = build_ladder(N, 1, 3, kind)
        sv = np.sort(np.linalg.svd(M, compute_uv=False))
        ok = (np.all(np.abs(sv[:105] - SQRT2) < 1e-12) and
              np.all(np.abs(sv[105:] - SQRT6) < 1e-12))
        check(ok, f"{kind}(1,3) singular values exactly {{sqrt2 x105, sqrt6 x15}}")
        p, q = SPLUS_CHAIN[-1] if kind == "splus" else PHI_CHAIN[-1]
        sv_t = np.sort(np.linalg.svd(build_ladder(N, p, q, kind),
                                     compute_uv=False))
        check(np.all(np.abs(sv_t - 2.0) < 1e-12) and len(sv_t) == 15,
              f"{kind} terminal rung singular values exactly {{2.0 x15}}")
    J = [1.55, 2.09, 1.85, 0.75, 0.90]
    blocks = {}
    for pq in set(SPLUS_CHAIN + PHI_CHAIN + list(BOUNDARY.values())):
        blocks[pq] = eig_block(pq[0], pq[1], J)
    _, w0, V0 = blocks[(1, 3)]
    for chain, kind in ((SPLUS_CHAIN, "splus"), (PHI_CHAIN, "phi")):
        n_pass = 0
        for k in range(len(w0)):
            v = V0[:, k] / np.linalg.norm(V0[:, k])
            rungs_ok, ratio = walk_gate_set(v, kind, chain, blocks)
            if rungs_ok and ratio is not None and ratio < 1e-12:
                n_pass += 1
        check(n_pass == 105,
              f"{kind}: 105/120 eigenvectors pass the full walk at the "
              f"off-locus control (selectivity is the family's, not the "
              f"loci's; got {n_pass})")
    check(control(J) == 0.0,
          "CONTROL intertwines exactly at the off-locus coupling too")


def main():
    structure_control()
    ein_bases = {k: entering_image_basis(k) for k in ("splus", "phi")}
    for qstar, species, sector, lam_census in SEEDS:
        J = [2.0 * qstar] * (N - 1)
        print(f"\n== seed q* = {qstar} ({species}, R-{sector}) ==")
        blocks = {}
        for pq in set(SPLUS_CHAIN + PHI_CHAIN + list(BOUNDARY.values())):
            blocks[pq] = eig_block(pq[0], pq[1], J)
        L13, w13, V13 = blocks[(1, 3)]
        if qstar in STURM_ONLY:
            sturm_only_sector_gate(L13, sector, lam_census)
        scale = np.abs(L13).sum(axis=1).max()   # inf-norm of the block
        k1, k2, lam, split = pair_data(w13, V13, lam_census)
        twin = -2.0 * N - lam
        check(abs(lam - lam_census) < 5e-6,
              f"seed pair at census lambda (dist {abs(lam - lam_census):.1e})")
        if species == "fold":
            check(abs(lam.real + 6.0) < 100 * EPS * scale,
                  f"fold-line pin Re lambda* = -6 "
                  f"({abs(lam.real + 6):.1e} vs eps scale {100 * EPS * scale:.1e})")
        else:
            check(abs(lam.imag) < 100 * EPS * scale,
                  f"axis pin Im lambda* = 0 "
                  f"({abs(lam.imag):.1e} vs eps scale {100 * EPS * scale:.1e})")
        # both freight pairs are defective-grade non-normal: the departure
        # from normality s of the restricted 2x2 dwarfs the split (a
        # semisimple degeneracy, e.g. the AT triples, has s ~ 0)
        for val, vname in ((lam, "lambda*"), (twin, "twin")):
            j1, j2, _, vsplit = pair_data(w13, V13, val)
            check(1e-8 < vsplit < 1e-2,
                  f"{vname} pair split truncation-sized, two-sided "
                  f"({vsplit:.1e}; model: 2*sqrt(|a|*|dq|), |dq| <= 1e-7)")
            s_dep = departure_from_normality(L13, V13, j1, j2, w13)
            check(s_dep > 100 * vsplit,
                  f"{vname} pair defective-grade non-normal "
                  f"(s {s_dep:.2e} = {s_dep / max(vsplit, 1e-300):.0f}x split)")
            for kind in ("splus", "phi"):
                proj = max(
                    np.linalg.norm(ein_bases[kind].conj().T @
                                   (V13[:, j] / np.linalg.norm(V13[:, j])))
                    for j in (j1, j2))
                check(proj < 1e-6,
                      f"{vname} lowest-weight on {kind} "
                      f"(entering-rung projection {proj:.1e})")
        for val, vname in ((lam, "lambda*"), (twin, "twin")):
            for chain, kind in ((SPLUS_CHAIN, "splus"), (PHI_CHAIN, "phi")):
                presence, cg, res_m, term, surv = walk_chain(
                    chain, kind, val, blocks)
                check(presence < 2e-3,
                      f"{vname} present in every {kind} seat spectrum "
                      f"(worst {presence:.1e})")
                check(cg < 1e-9,
                      f"{vname} {kind} rungs at CG sqrt2 (worst dev {cg:.1e})")
                check(res_m < 1e-9,
                      f"{vname} {kind} member images are eigenvectors "
                      f"(worst Rayleigh residual {res_m:.1e})")
                check(term < 1e-12,
                      f"{vname} {kind} terminal death at floor ({term:.1e})")
                check(surv == {15},
                      f"{vname} {kind} terminal rung surjective onto the "
                      f"15-dim boundary ({surv})")
        w_ctl = control(J)
        check(w_ctl == 0.0, "CONTROL both ladders intertwine exactly",
              repr(w_ctl))
        at_control(blocks, ein_bases)
    print(f"\n{'ALL GATES PASS' if failures == 0 else f'{failures} FAILURES'}")
    sys.exit(0 if failures == 0 else 1)


if __name__ == "__main__":
    main()
