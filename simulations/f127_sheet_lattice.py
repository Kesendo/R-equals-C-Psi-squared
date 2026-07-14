"""The fragile-thing hunt, scout 1: the residue lattice of cross_form.

DISCRIMINATING QUESTION (HANDOVER 2026-07-14): fix ONE pole location and sum the
residues of cross_form there over the lattice of atoms sharing it; does the
constraint Sum cos = 0 make that residue-sum vanish family-wise (telescoping),
or only globally?

STRUCTURE (verified backbone): cross_form = (1/4) Sum_{i,j,xi,up,e} (-1)^(i+j) (-e)
c_xi c_up (1+cos phi) [ 1/(cos phi - cos(a_i+b_j)) - 1/(cos phi - cos(a_i-b_j)) ],
phi = xi + e*up.  72 atoms, 144 simple poles, 288 sheet-events.  Every pole sheet
is {L = 0 mod 2pi} for a linear form L with ALL SIX coefficients in {+-1}:
at most 32 sheets up to sign.

Near a sheet with canonical coordinate ell:  cos phi - cos mu_s ~ -eps sin(phi) ell,
so the event residue in ell is  w = -A c_s eps / sin(phi),
A = (1/4)(-1)^(i+j)(-e) c_xi c_up (1+cos phi),  c_s = +1 (s=+1) / -1 (s=-1).

GATES first (fail loud), then the survey:
  G1: Sum(atoms) == cross_form at random generic points (S5 transcription correct).
  G2: residue bookkeeping == numeric lim ell * cross_form along a transversal.
Then, per sheet, the residue-sum R on: sheet only / sheet+Ca / sheet+Cb / sheet+Ca+Cb,
plus sub-family partial sums (by (i,j), by i, by j, by e) and a pairwise-matching hunt.
"""
import itertools
import math
import random
import sys

sys.path.insert(0, "simulations")
from cross_triple_orthogonality import cross_form  # the committed object

TWO_PI = 2 * math.pi


# ---------------------------------------------------------------- atoms and events
def pieces(ang):
    """Committed _pieces convention (cross_triple_orthogonality.py:454-463)."""
    P = []
    for i in range(3):
        j, l = [t for t in range(3) if t != i]
        P.append(dict(
            # coefficient 6-vectors (a1..a3 or b1..b3 slots) for psum/pdif
            psum_vec={j: 1, l: 1}, pdif_vec={j: -1, l: 1},
        ))
    return P


def angval(vec, ang, off):
    """Value of a coefficient dict {slot: c} over ang, slot offset off (0=a,3=b)."""
    return sum(c * ang[k] for k, c in vec.items())


def coeff(ang, i, which):
    """c_xi: alpha for psum, beta for pdif (committed convention)."""
    j, l = [t for t in range(3) if t != i]
    if which == "psum":
        return math.sin(ang[l]) - math.sin(ang[j])
    return -(math.sin(ang[l]) + math.sin(ang[j]))


def enum_atoms():
    """All 72 atoms as index tuples (i, j, xi, up, e)."""
    return [(i, j, xi, up, e)
            for i in range(3) for j in range(3)
            for xi in ("psum", "pdif") for up in ("psum", "pdif")
            for e in (1, -1)]


def atom_value(atom, a, b):
    i, j, xi, up, e = atom
    Pa, Pb = pieces(a), pieces(b)
    xi_val = angval(Pa[i][xi + "_vec"], a, 0)
    up_val = angval(Pb[j][up + "_vec"], b, 0)
    phi = xi_val + e * up_val
    A = 0.25 * ((-1) ** (i + j)) * (-e) * coeff(a, i, xi) * coeff(b, j, up) * (1 + math.cos(phi))
    mp, mm = a[i] + b[j], a[i] - b[j]
    return A * (1.0 / (math.cos(phi) - math.cos(mp)) - 1.0 / (math.cos(phi) - math.cos(mm)))


def atom_sum(a, b):
    return sum(atom_value(t, a, b) for t in enum_atoms())


def enum_events():
    """288 events: (atom, s, tau) with sheet form L = phi - tau*(a_i + s*b_j).

    Returns list of (atom, s, tau, Lvec) with Lvec a 6-tuple of {+-1} coefficients
    over (a1,a2,a3,b1,b2,b3), NOT yet canonicalized.
    """
    out = []
    Pa_v, Pb_v = pieces([0, 0, 0]), pieces([0, 0, 0])   # only the vecs matter
    for atom in enum_atoms():
        i, j, xi, up, e = atom
        vec = [0] * 6
        for k, c in Pa_v[i][xi + "_vec"].items():
            vec[k] += c
        for k, c in Pb_v[j][up + "_vec"].items():
            vec[3 + k] += e * c
        for s in (1, -1):
            for tau in (1, -1):
                L = list(vec)
                L[i] -= tau
                L[3 + j] -= tau * s
                assert all(x in (-1, 1) for x in L), (atom, s, tau, L)
                out.append((atom, s, tau, tuple(L)))
    return out


def canon(L):
    """Canonical sheet key: global sign so first coefficient is +1; eps = that sign."""
    eps = 1 if L[0] > 0 else -1
    return tuple(eps * x for x in L), eps


def event_weight(atom, s, tau, eps, a, b):
    """Residue in the canonical sheet coordinate ell of this event, at a sheet point."""
    i, j, xi, up, e = atom
    Pa, Pb = pieces(a), pieces(b)
    xi_val = angval(Pa[i][xi + "_vec"], a, 0)
    up_val = angval(Pb[j][up + "_vec"], b, 0)
    phi = xi_val + e * up_val
    A = 0.25 * ((-1) ** (i + j)) * (-e) * coeff(a, i, xi) * coeff(b, j, up) * (1 + math.cos(phi))
    c_s = 1 if s == 1 else -1
    return -A * c_s * eps / math.sin(phi)


# ---------------------------------------------------------------- point sampling
def rand_ang():
    return random.uniform(-math.pi, math.pi)


def solve_ca():
    """Random (a1,a2,a3) with cos a1 + cos a2 + cos a3 = 0."""
    while True:
        a1, a2 = rand_ang(), rand_ang()
        c = -math.cos(a1) - math.cos(a2)
        if abs(c) <= 0.98:
            a3 = random.choice([1, -1]) * math.acos(c)
            return [a1, a2, a3]


def on_sheet_solve(L, fixed, free_slot, branch=0):
    """Solve L . x = 2 pi * branch for x[free_slot]; L coeffs are +-1."""
    rest = sum(L[k] * fixed[k] for k in range(6) if k != free_slot)
    return (TWO_PI * branch - rest) / L[free_slot]


def sample_sheet(L, mode, branch=0, tries=400):
    """A random point on sheet {L=0 mod 2pi}, with constraints per mode:
    'none' | 'ca' | 'cb' | 'cacb'."""
    for _ in range(tries):
        if mode == "none":
            x = [rand_ang() for _ in range(6)]
            x[5] = on_sheet_solve(L, x, 5, branch)
        elif mode == "ca":
            a = solve_ca()
            x = a + [rand_ang(), rand_ang(), 0.0]
            x[5] = on_sheet_solve(L, x, 5, branch)
        elif mode == "cb":
            b = solve_ca()
            x = [rand_ang(), rand_ang(), 0.0] + b
            x[2] = on_sheet_solve(L, x, 2, branch)
        elif mode == "cacb":
            # a on Ca; b1 free; solve b2 by 1-d root-find so that with
            # b3 = sheet-solve(b2), cos b1 + cos b2 + cos b3 = 0.
            a = solve_ca()
            b1 = rand_ang()

            def g(b2):
                x_ = a + [b1, b2, 0.0]
                b3 = on_sheet_solve(L, x_, 5, branch)
                return math.cos(b1) + math.cos(b2) + math.cos(b3)

            # scan for a sign change
            grid = [(-math.pi) + TWO_PI * k / 720 for k in range(721)]
            vals = [g(t) for t in grid]
            iv = [(grid[k], grid[k + 1]) for k in range(720)
                  if vals[k] == 0 or (vals[k] < 0) != (vals[k + 1] < 0)]
            if not iv:
                continue
            lo, hi = random.choice(iv)
            for _ in range(80):                      # bisection
                mid = 0.5 * (lo + hi)
                if (g(lo) < 0) != (g(mid) < 0):
                    hi = mid
                else:
                    lo = mid
            b2 = 0.5 * (lo + hi)
            x = a + [b1, b2, 0.0]
            x[5] = on_sheet_solve(L, x, 5, branch)
        else:
            raise ValueError(mode)
        if guard_ok(x, L):
            return x
    return None


def guard_ok(x, L, tol=5e-2):
    """Reject points where an UNRELATED pole factor or sin(phi) is nearly singular."""
    a, b = x[:3], x[3:]
    for atom, s, tau, Le in EVENTS:
        i, j, xi, up, e = atom
        Pa, Pb = pieces(a), pieces(b)
        phi = angval(Pa[i][xi + "_vec"], a, 0) + e * angval(Pb[j][up + "_vec"], b, 0)
        if abs(math.sin(phi)) < tol:
            return False
        key, _ = canon(Le)
        if key != canon(L)[0]:
            # unrelated sheet must be far
            val = sum(Le[k] * x[k] for k in range(6))
            if abs(math.remainder(val, TWO_PI)) < tol:
                return False
    # mu_m removable-singularity guard for cross_form itself is irrelevant here
    # (we never evaluate cross_form at sheet points), but keep sin(phi) guard above.
    return True


# ---------------------------------------------------------------- gates
EVENTS = enum_events()


def gate_atoms(n=6, tol=1e-9):
    worst = 0.0
    for _ in range(n):
        a = [rand_ang() for _ in range(3)]
        b = [rand_ang() for _ in range(3)]
        cf, at = cross_form(a, b), atom_sum(a, b)
        den = max(1.0, abs(cf))
        worst = max(worst, abs(cf - at) / den)
    assert worst < tol, f"G1 FAIL: atom transcription off by {worst:.2e}"
    return worst


def gate_residue(tol=2e-4):
    """G2: bookkeeping residue == numeric lim ell*cross_form on one random sheet."""
    random.seed(20260714)
    atom, s, tau, L = random.choice(EVENTS)
    key, _ = canon(L)
    x0 = sample_sheet(key, "none")
    assert x0 is not None
    fam = [(at, ss, tt, ee) for (at, ss, tt, Le) in EVENTS
           for (kk, ee) in [canon(Le)] if kk == key]
    a0, b0 = x0[:3], x0[3:]
    R = sum(event_weight(at, ss, tt, ee, a0, b0) for at, ss, tt, ee in fam)
    # numeric: move along +e_0 direction scaled so that d ell = 1 (L[0] = +1)
    est = []
    for h in (1e-5, 1e-6):
        x = list(x0)
        x[0] += h        # canonical L has coefficient +1 on slot 0
        est.append(h * cross_form(x[:3], x[3:]))
    err = abs(est[1] - R) / max(1.0, abs(R))
    assert err < tol, f"G2 FAIL: bookkeeping {R:.6e} vs numeric {est} (err {err:.2e})"
    return R, est


# ---------------------------------------------------------------- the survey
def survey(pts=6, seed=1):
    random.seed(seed)
    sheets = {}
    for atom, s, tau, L in EVENTS:
        key, eps = canon(L)
        sheets.setdefault(key, []).append((atom, s, tau, eps))
    print(f"[survey] {len(EVENTS)} events on {len(sheets)} distinct sheets "
          f"(multiplicities: {sorted(set(len(v) for v in sheets.values()))})")

    modes = ["none", "ca", "cb", "cacb"]
    table = {}
    for key, fam in sorted(sheets.items()):
        row = {}
        for mode in modes:
            vals, got = [], 0
            for _ in range(pts):
                x = sample_sheet(key, mode)
                if x is None:
                    continue
                got += 1
                a, b = x[:3], x[3:]
                ws = [(evt, event_weight(*evt[:3], evt[3], a, b)) for evt in fam]
                vals.append(abs(sum(w for _, w in ws)))
            if not got:
                row[mode] = None
                continue
            # partition-level verdicts live in survey_granularity()
            row[mode] = dict(total=max(vals) if vals else None, n=got)
        table[key] = (fam, row)

    # print the headline: total residue per sheet per mode
    print("\nsheet (canonical L over a1 a2 a3 b1 b2 b3) | #evt | max|R| per mode")
    print("mode order: none / +Ca / +Cb / +Ca+Cb")
    for key, (fam, row) in sorted(table.items()):
        cells = []
        for mode in modes:
            r = row[mode]
            cells.append("     -    " if (r is None or r["total"] is None)
                         else f"{r['total']:9.2e}")
        print("  ".join(f"{c:+d}" for c in key), f"| {len(fam):2d} |", "  ".join(cells))
    return sheets


def survey_granularity(key, fam, pts=5, seed=7):
    """For ONE sheet: per-point partial sums over candidate partitions + pair hunt."""
    random.seed(seed)
    partitions = {
        "(i,j)": lambda atom, s, tau: (atom[0], atom[1]),
        "i": lambda atom, s, tau: atom[0],
        "j": lambda atom, s, tau: atom[1],
        "e": lambda atom, s, tau: atom[4],
        "s": lambda atom, s, tau: s,
        "(xi,up)": lambda atom, s, tau: (atom[2], atom[3]),
        "(i,j,e)": lambda atom, s, tau: (atom[0], atom[1], atom[4]),
    }
    results = {p: 0.0 for p in partitions}
    W = []           # per-point weight vectors for the pair hunt
    for _ in range(pts):
        x = sample_sheet(key, "cacb")
        if x is None:
            continue
        a, b = x[:3], x[3:]
        ws = [((atom, s, tau), event_weight(atom, s, tau, eps, a, b))
              for atom, s, tau, eps in fam]
        W.append([w for _, w in ws])
        scale = max(abs(w) for _, w in ws)
        for pname, kf in partitions.items():
            groups = {}
            for (atom, s, tau), w in ws:
                groups.setdefault(kf(atom, s, tau), 0.0)
                groups[kf(atom, s, tau)] += w
            results[pname] = max(results[pname],
                                 max(abs(v) for v in groups.values()) / scale)
    # pair hunt: events k,l with w_k + w_l ~ 0 at EVERY point
    n = len(fam)
    pairs = []
    if W:
        for k in range(n):
            for l in range(k + 1, n):
                if all(abs(Wp[k] + Wp[l]) < 1e-9 * max(abs(Wp[k]), abs(Wp[l]), 1e-30)
                       for Wp in W):
                    pairs.append((fam[k][:3], fam[l][:3]))
    return results, pairs, len(W)


def main():
    random.seed(20260714)
    g1 = gate_atoms()
    print(f"[G1 PASS] atom decomposition == cross_form, worst rel err {g1:.2e}")
    R, est = gate_residue()
    print(f"[G2 PASS] residue bookkeeping {R:.6e} vs numeric {est[1]:.6e}")

    sheets = survey(pts=5, seed=3)

    # granularity on every sheet, on V (Ca+Cb)
    print("\n[granularity on sheet cap V]  max_normalized |partial sum| per partition")
    header = None
    for key, fam in sorted(sheets.items()):
        res, pairs, npts = survey_granularity(key, fam, pts=4, seed=11)
        if header is None:
            header = list(res.keys())
            print("sheet | npts | " + " | ".join(f"{h:>8s}" for h in header)
                  + " | perfect pairs")
        print("".join(f"{c:+d}" for c in key), f"| {npts} | "
              + " | ".join(f"{res[h]:8.1e}" for h in header)
              + f" | {len(pairs)}")
        if pairs:
            for p in pairs[:6]:
                print("      pair:", p)


if __name__ == "__main__":
    main()
