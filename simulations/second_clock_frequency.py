"""Gate-first: the SECOND-CLOCK FLOOR FREQUENCY unifies the chain-N=3 and ring-N=4 specials.

The chain gap-dominance proof (chain N=3) and the ring gap-dominance proof (ring N=4) each carry a lone
small-N exception: a {0,2}-coherence sqrt-EP mode PINNED on the Re=−2γ floor (above its own EP, i.e. for
Q > Q* = 2J/B, because that sector's block closes there). Treated as two separate specials. This verifier
shows they are ONE law:

  THE SECOND-CLOCK FLOOR FREQUENCY:  |Im| = sqrt(B^2 − (2γ)^2)  for B > 2γ,
    the 2x2 {0,2}-block [[0, B], [−B, −4γ]]: a population (Hamming weight 0, dissipator cell 0) coupled to
    a 2-Hamming coherence (dissipator cell −4γ) by the bare, γ-INDEPENDENT amplitude B. The two cells average
    to −2γ, which is what puts the PAIR on the floor; the eigenvalues are −2γ ± sqrt((2γ)^2 − B^2), so the
    pair oscillates at sqrt(B^2 − (2γ)^2) above the EP at B = 2γ (i.e. Q* = 2J/B), is DEFECTIVE there, and
    below it splits along the real axis and leaves the floor. B = the FREE-FERMION BAND TOP of the sector
    the {0,2} block lives in (the γ→0 coupling). One envelope; only B is topology/sector-specific.

  Three on-floor instances (the {0,2} pins to −2γ only at these special N):
    chain N=3 (1,1):  B = E1   = 2J cos(π/4) = √2 J   (single-particle band edge; floor freq BELOW E1)
    ring  N=4 (1,1):  B = 2J             (periodic single-particle band top)
    ring  N=4 (2,2):  B = 2√2 J = 2·(2J cos(π/4))  (anti-periodic TWO-fermion top; floor freq ABOVE 2J)

  "Second clock exceeds the first (band edge J·ρ)" reduces to B vs J·ρ: ring (2,2) 2√2 > 2 (exceeds);
  chain N=3 B=E1=J·ρ_chain and ring (1,1) B=2J=J·ρ_ring (equal at γ→0, the √-shift puts the floor freq just below).

  STAGE 0  the envelope: |Im| at the floor = sqrt(B^2−(2γ)^2), B γ-INDEPENDENT (γ-sweep), at all 3 instances.
  STAGE 1  B = the free-fermion sector band top (closed forms): E1 / 2J / 2√2 J.
  STAGE 2  exceeds-the-first-clock iff B > J·ρ (ring (2,2) yes; the (1,1) cases marginal).
  STAGE 3  the block itself: the two diagonal cells (0 and −4γ) read off the full Liouvillian, the defect at
           B = 2γ, and the crossing (below the EP the pair leaves the floor) -- the discriminator against the
           equal-diagonal form [[−2γ, ω], [−ω, −2γ]], which would keep the pair on the floor at every γ.
  Scope note: the {0,2} pins to the −2γ floor only at chain N=3 / ring N=4 (other N: it drifts off-floor);
  this is the unification of those two specials, not an all-N floor law.

Run: python simulations/second_clock_frequency.py
"""
import sys
import numpy as np
from itertools import combinations
from math import cos, pi, sqrt

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

J = 1.0
TOL = 1e-7


def basis(N, p):
    return [sum(1 << i for i in c) for c in combinations(range(N), p)]


def H_p(N, p, ring):
    states = basis(N, p)
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)), complex)
    bonds = [(b, b + 1) for b in range(N - 1)] + ([(N - 1, 0)] if ring else [])
    for a, b in bonds:
        for s in states:
            if (s >> a) & 1 and not (s >> b) & 1:
                H[idx[(s & ~(1 << a)) | (1 << b)], idx[s]] += J
            if (s >> b) & 1 and not (s >> a) & 1:
                H[idx[(s & ~(1 << b)) | (1 << a)], idx[s]] += J
    return H, states


def floor_freq(N, p, g, ring):
    """max|Im| over the (p,p)-block modes pinned at Re = −2γ (the n_XY=1 floor); None if none."""
    H, s = H_p(N, p, ring)
    n = len(s)
    L = -1j * (np.kron(H, np.eye(n)) - np.kron(np.eye(n), H.T))
    deph = np.array([-2.0 * g * bin(s[a] ^ s[b]).count('1') for a in range(n) for b in range(n)])
    ev = np.linalg.eigvals(L + np.diag(deph))
    at = np.abs(ev.real - (-2 * g)) < 1e-6
    ims = np.abs(ev.imag[at]); ims = ims[ims > TOL]
    return float(ims.max()) if ims.size else None


# the three on-floor instances and their predicted free-fermion band-top B
E1_N3 = 2 * J * cos(pi / 4)          # = sqrt2 J, the chain N=3 single-particle band edge
INSTANCES = [
    ("chain N=3 (1,1)", False, 3, 1, E1_N3,            "E1 = 2J cos(pi/4) = sqrt2 J (single-particle)"),
    ("ring  N=4 (1,1)", True,  4, 1, 2 * J,            "2J (periodic single-particle band top)"),
    ("ring  N=4 (2,2)", True,  4, 2, 2 * sqrt(2) * J,  "2sqrt2 J (anti-periodic TWO-fermion top)"),
]

# ====================================================================================================
# STAGE 0 -- THE ENVELOPE: |Im| at the floor = sqrt(B^2 - (2g)^2), B gamma-INDEPENDENT
# ====================================================================================================
print("=" * 104)
print("STAGE 0 -- the second-clock floor frequency is the envelope |Im| = sqrt(B^2 - (2g)^2), B gamma-indep")
print("=" * 104)
print(f"{'instance':18} {'g':>5} {'|Im|@floor':>12} {'B=sqrt(Im^2+4g^2)':>18} {'B match across g?':>18}")
B_fit = {}
for name, ring, N, p, B_pred, _ in INSTANCES:
    Bs = []
    for g in (0.05, 0.10, 0.20):
        im = floor_freq(N, p, g, ring)
        assert im is not None, f"STAGE 0 GATE FIRED: {name} has no on-floor {{0,2}} sqrt-EP at g={g}"
        B = sqrt(im ** 2 + (2 * g) ** 2)
        Bs.append(B)
        print(f"{name:18} {g:>5.2f} {im:>12.6f} {B:>18.6f} {'':>18}")
    spread = max(Bs) - min(Bs)
    B_fit[name] = Bs[0]
    assert spread < 1e-6, f"STAGE 0 GATE FIRED: {name} B not gamma-independent (spread {spread:.2e})"
    print(f"{'':18} {'-->':>5} {'':>12} {Bs[0]:>18.6f} {f'spread {spread:.1e} OK':>18}")
print("\nSTAGE 0 PASS: every on-floor {0,2} mode obeys |Im| = sqrt(B^2 - (2g)^2) with B gamma-INDEPENDENT")
print("  (the 2x2 [[0, B], [-B, -4g]] envelope; B = the bare gamma->0 coupling, the cells 0 / -4g average to -2g).")

# ====================================================================================================
# STAGE 1 -- B = the FREE-FERMION BAND TOP of the sector
# ====================================================================================================
print("\n" + "=" * 104)
print("STAGE 1 -- B = the free-fermion band top of the {0,2} block's sector (closed forms)")
print("=" * 104)
print(f"{'instance':18} {'B (fitted)':>12} {'B (FF top)':>12} {'closed form':>46} {'match':>6}")
for name, ring, N, p, B_pred, desc in INSTANCES:
    fit = B_fit[name]
    ok = abs(fit - B_pred) < 1e-6
    print(f"{name:18} {fit:>12.6f} {B_pred:>12.6f} {desc:>46} {('YES' if ok else 'NO'):>6}")
    assert ok, f"STAGE 1 GATE FIRED: {name} B {fit} != FF top {B_pred} ({desc})"
print("\nSTAGE 1 PASS: B is the free-fermion sector band top -- single-particle E1 for the (1,1) sectors,")
print("  the anti-periodic TWO-fermion top 2sqrt2 J for the ring (2,2) half-filling. One envelope, B per sector.")

# ====================================================================================================
# STAGE 2 -- "second clock exceeds the first (band edge J*rho)" reduces to B vs J*rho
# ====================================================================================================
print("\n" + "=" * 104)
print("STAGE 2 -- exceeds the first clock (band edge J*rho) iff B > J*rho")
print("=" * 104)
rho = {("chain", 3): 2 * cos(pi / 4), ("ring", 4): 2.0}     # chain rho = 2cos(pi/(N+1)); ring rho = 2
print(f"{'instance':18} {'B':>10} {'J*rho (band edge)':>18} {'exceeds?':>9}")
for name, ring, N, p, B_pred, _ in INSTANCES:
    r = J * rho[("ring" if ring else "chain", N)]
    exceeds = B_pred > r + 1e-9
    print(f"{name:18} {B_pred:>10.6f} {r:>18.6f} {('YES (above)' if exceeds else 'no (<= band edge)'):>9}")
# the ring (2,2) is the ONLY one that exceeds -- the genuine "second clock overtakes first"
assert (2 * sqrt(2) * J) > 2 * J + 1e-9, "ring (2,2) must exceed the ring band edge 2J"
assert E1_N3 <= 2 * cos(pi / 4) * J + 1e-9, "chain N=3 B should equal its own band edge (not exceed)"
print("\nSTAGE 2 PASS: the ring (2,2) is the lone exceed-case (B=2sqrt2 > 2J = the anti-periodic two-fermion top")
print("  overtaking the periodic band edge); the (1,1) instances have B = J*rho (the √-shift puts the floor freq")
print("  just below). So 'second clock overtakes first' = 'a multi-fermion sector top exceeds J*rho'.")

# ====================================================================================================
# STAGE 3 -- THE BLOCK ITSELF: the cells 0 / -4g, the defect at B = 2g, the crossing off the floor.
#   The equal-diagonal form [[-2g, w], [-w, -2g]] has eigenvalues -2g +/- i*w at EVERY g: it would hold the
#   pair on the floor for all Q and leave no EP anywhere. The real block has UNEQUAL cells (0 and -4g,
#   averaging to -2g), and the discriminator is the crossing: past B = 2g its eigenvalues turn real and the
#   pair disappears from the floor.
# ====================================================================================================
print("\n" + "=" * 104)
print("STAGE 3 -- the {0,2} block from below: the cells 0 / -4g, the defect at B = 2g, the crossing")
print("=" * 104)


def block(B, g):
    return np.array([[0.0, B], [-B, -4.0 * g]], dtype=complex)


g0 = 0.05
H4, s4 = H_p(4, 2, True)
n4 = len(s4)
pop_cell = -2.0 * g0 * bin(s4[0] ^ s4[0]).count('1')
coh = [(a, b) for a in range(n4) for b in range(n4) if bin(s4[a] ^ s4[b]).count('1') == 2]
coh_cell = -2.0 * g0 * bin(s4[coh[0][0]] ^ s4[coh[0][1]]).count('1')
print(f"  ring N=4 (2,2), g={g0}: dissipator cell of a population {pop_cell:+.4f}, of a 2-Hamming coherence "
      f"{coh_cell:+.4f}, mean {0.5*(pop_cell+coh_cell):+.4f} = -2g")
assert abs(pop_cell) < 1e-12 and abs(coh_cell + 4 * g0) < 1e-12, "STAGE 3 GATE FIRED: the cells are not 0 / -4g"
assert abs(0.5 * (pop_cell + coh_cell) + 2 * g0) < 1e-12, "STAGE 3 GATE FIRED: the cells do not average to -2g"

print(f"{'instance':18} {'B':>10} {'block |Im|':>12} {'measured |Im|':>14} {'match':>6}")
for name, ring, N, p, B_pred, _ in INSTANCES:
    im_meas = floor_freq(N, p, g0, ring)
    w = np.linalg.eigvals(block(B_pred, g0))
    im_block = float(np.abs(w.imag).max())
    ok = abs(im_block - im_meas) < 1e-9
    print(f"{name:18} {B_pred:>10.6f} {im_block:>12.6f} {im_meas:>14.6f} {('YES' if ok else 'NO'):>6}")
    assert ok, f"STAGE 3 GATE FIRED: {name} block frequency {im_block} != measured {im_meas}"
    # the block is defective exactly at B = 2g, and the equal-diagonal form is not
    gEP = B_pred / 2
    M = block(B_pred, gEP)
    rank = np.linalg.matrix_rank(M + 2 * gEP * np.eye(2), tol=1e-9)
    equal_diag = np.array([[-2 * gEP, B_pred], [-B_pred, -2 * gEP]], dtype=complex)
    rank_eq = np.linalg.matrix_rank(equal_diag + 2 * gEP * np.eye(2), tol=1e-9)
    assert rank == 1, f"STAGE 3 GATE FIRED: {name} block not defective at B = 2g"
    assert rank_eq == 2, f"STAGE 3 GATE FIRED: the equal-diagonal form should stay semisimple at B = 2g"
    print(f"{'':18} {'EP at g = B/2 =':>10} {gEP:.6f}: rank(M+2g) = {rank} (defective) vs {rank_eq} for "
          f"[[-2g,w],[-w,-2g]] (semisimple, no EP)")

# the crossing, measured on the real Liouvillian: below Q* the pair is gone from the floor
print(f"\n{'instance':18} {'g':>8} {'2g vs B':>12} {'on-floor |Im|':>14}")
for name, ring, N, p, B_pred, _ in INSTANCES:
    for g in (0.9 * B_pred / 2, 1.1 * B_pred / 2):
        im = floor_freq(N, p, g, ring)
        rel = "2g < B (above EP)" if 2 * g < B_pred else "2g > B (below EP)"
        print(f"{name:18} {g:>8.4f} {rel:>12} {(f'{im:.6f}' if im is not None else 'GONE from the floor'):>14}")
        if 2 * g < B_pred:
            assert im is not None and abs(im - sqrt(B_pred ** 2 - (2 * g) ** 2)) < 1e-6, \
                f"STAGE 3 GATE FIRED: {name} should still ride the floor at g={g}"
        else:
            assert im is None, f"STAGE 3 GATE FIRED: {name} should have left the floor at g={g} (past its EP)"
print("\nSTAGE 3 PASS: the block is [[0, B], [-B, -4g]] (cells 0 and -4g, mean -2g), defective exactly at")
print("  B = 2g, and past its EP the pair leaves the floor -- which the equal-diagonal form cannot produce.")

print("\n" + "=" * 104)
print("UNIFIED: the chain-N=3 and ring-N=4 specials are ONE law -- the second-clock {0,2} sqrt-EP floor")
print("frequency = sqrt(B^2 - (2g)^2), B = the free-fermion band top of its sector. The EP is at B=2g (Q*=2J/B).")
print("Scope: the {0,2} pins to the -2g floor only at chain N=3 / ring N=4 (the unification of those specials).")
print("DONE.")
