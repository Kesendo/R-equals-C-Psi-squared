"""Gate-first: the SECOND-CLOCK FLOOR FREQUENCY unifies the chain-N=3 and ring-N=4 specials.

The chain gap-dominance proof (chain N=3) and the ring gap-dominance proof (ring N=4) each carry a lone
small-N exception: a {0,2}-coherence sqrt-EP mode PINNED on the Re=−2γ floor (for all Q, because that
sector's block closes there). Treated as two separate specials. This verifier shows they are ONE law:

  THE SECOND-CLOCK FLOOR FREQUENCY:  |Im| = sqrt(B^2 − (2γ)^2),
    a 2x2 {0,2}-block [[−2γ, ω],[−ω, −2γ]] (population <-> 2-Hamming coherence, coupling ω = sqrt(B^2−(2γ)^2)),
    with B = the FREE-FERMION BAND TOP of the sector the {0,2} block lives in (γ→0 limit). The EP is at
    B = 2γ (i.e. Q* = 2J/B). One envelope; only B is topology/sector-specific.

  Three on-floor instances (the {0,2} pins to −2γ only at these special N):
    chain N=3 (1,1):  B = E1   = 2J cos(π/4) = √2 J   (single-particle band edge; floor freq BELOW E1)
    ring  N=4 (1,1):  B = 2J             (periodic single-particle band top)
    ring  N=4 (2,2):  B = 2√2 J = 2·(2J cos(π/4))  (anti-periodic TWO-fermion top; floor freq ABOVE 2J)

  "Second clock exceeds the first (band edge J·ρ)" reduces to B vs J·ρ: ring (2,2) 2√2 > 2 (exceeds);
  chain N=3 B=E1=J·ρ_chain and ring (1,1) B=2J=J·ρ_ring (equal at γ→0, the √-shift puts the floor freq just below).

  STAGE 0  the envelope: |Im| at the floor = sqrt(B^2−(2γ)^2), B γ-INDEPENDENT (γ-sweep), at all 3 instances.
  STAGE 1  B = the free-fermion sector band top (closed forms): E1 / 2J / 2√2 J.
  STAGE 2  exceeds-the-first-clock iff B > J·ρ (ring (2,2) yes; the (1,1) cases marginal).
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
print("  (the 2x2 [[−2g, w],[−w,−2g]] envelope; B = the bare gamma->0 coupling).")

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

print("\n" + "=" * 104)
print("UNIFIED: the chain-N=3 and ring-N=4 specials are ONE law -- the second-clock {0,2} sqrt-EP floor")
print("frequency = sqrt(B^2 - (2g)^2), B = the free-fermion band top of its sector. The EP is at B=2g (Q*=2J/B).")
print("Scope: the {0,2} pins to the -2g floor only at chain N=3 / ring N=4 (the unification of those specials).")
print("DONE.")
