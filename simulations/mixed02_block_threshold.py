"""The {0,2} extras on the floor: a low-N threshold whose rung moves with the model.

Companion verifier for PROOF_CHAIN_GAP_DOMINANCE section 4.

Section 4 of that proof records four extra modes at N = 3 that are not free
fermions: equal-particle-number (n,n) coherences with n_XY histogram
{0: 1/2, 2: 1/2}, sitting at exactly Re = -2*gamma because their two-level
{n_XY = 0 <-> n_XY = 2} block closes. The proof called this "one more entry in
the chain's list of N = 3 accidents". It is not an accident and it is not
N = 3's. Two things this script establishes.

  THE RUNG MOVES WITH THE MODEL. XY closes the block at N = 2 and N = 3,
  Heisenberg only at N = 2. The ZZ term costs exactly one rung, and the reason
  is visible on the single-excitation sector's ZZ diagonal: CONSTANT at N = 2
  (so ZZ commutes with everything on the sector and cannot break the block),
  NOT constant at N = 3 (so it breaks it).

  AT N = 2 THE EXTRAS EXCEED THE BAND EDGE, and only above a gamma threshold.
  The proof's headline statement, max|Im| = E1 on the floor, is FALSE at N = 2
  for Q = J/gamma > 2/sqrt(3) and TRUE at or below it. The document never
  treats N = 2 and asserts the statement unqualified in N.

Gates (all must pass). Most are two-sided, i.e. they can fail in both
directions. Two are honestly not independent checks and are labelled as
such below: "G6 consequence" follows from "G6 diagonal" by algebra (if the
ZZ diagonal is constant on the sector then ZZ acts as a multiple of the
identity there, so the two models MUST agree), and it is kept as a
derivation check, not as evidence.

  G1  the extras exist exactly at XY N = 2, XY N = 3 and Heisenberg N = 2, and
      nowhere else in N = 2..5 for either model, at EVERY gamma swept.
  G2  every extra found has the {0: 1/2, 2: 1/2} histogram of the {0,2} family.
  G3  the closed form |Im| = sqrt(c^2 - (2 gamma)^2) holds, gamma-swept, with
      c = 2J at N = 2 and c = E1 = 2J cos(pi/(N+1)) at N = 3.
  G4  THE CROSSOVER. At N = 2 the extras sit above E1 exactly for
      Q > 2/sqrt(3) ~ 1.1547. Swept from Q = 0.5 to Q = 20, the measured facing
      must flip at that Q and nowhere else. (This replaces an earlier "the N = 2
      extras are above the band edge" check, which was run only at Q = 20 and
      was therefore unfalsifiable in the regime where it is false.)
  G5  UNIFORM BONDS ARE REQUIRED. Detuning one bond by delta takes the extras
      off the floor. The gate also records the smallest delta a given tolerance
      can still see, because the earlier verifier ran at TOL = 1e-8 and passed
      on detuned input where the property it tests is false.
  G6  THE MECHANISM. The single-excitation ZZ diagonal is constant at N = 2 and
      NOT constant at N = 3; consequently Heisenberg N = 2 reproduces the XY
      N = 2 frequency to machine precision, and Heisenberg has no extras at N = 3.
      (At N >= 4 the diagonal is not constant either, but it is not degenerate-free
      -- N = 4 reads [0.5, -0.5, -0.5, 0.5]. What matters is only "constant or not".)
  G7  BOND DISORDER DOES NOT BREAK THE MAXIMALITY. With random bond weights at
      N = 4 and 5 the floor's max|Im| still equals the largest single-particle
      energy of the WEIGHTED hopping matrix. This is the half of the bond
      question that section 4.3 says is NOT fenced, and it needs its own
      measurement because G5 measures the other half.

Conventions. The builder takes H = J*sum(XX + YY [+ ZZ]) over the bonds. The
proof's convention is H = (J/2)*sum(XX + YY), so the proof's J = 1 is this
builder's J = 0.5. The two conventions are NOT related by halving the
frequencies: gamma does not scale with J, so sqrt((4J)^2 - (2g)^2) at builder
J = 1 is 3.998750 while the proof-convention value is
sqrt((2J)^2 - (2g)^2) = 1.997498. Halve J, never the answer.

Run: python simulations/mixed02_block_threshold.py
Output: simulations/results/mixed02_block_threshold.txt
"""

import os
import sys

import numpy as np

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "mixed02_block_threshold.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")


X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

# The floor test is an EXACT membership test (Re = -2*gamma to machine
# precision), so the tolerance is a machine-precision one, not a physical
# window.  The earlier version of this script used 1e-8, which is wide enough
# to admit modes that a bond detuning of 1e-4 had already moved off the floor.
TOL = 1e-11


def op_at(n, site, p):
    m = np.array([[1.0 + 0j]])
    for k in range(n):
        m = np.kron(m, p if k == site else I2)
    return m


def hamiltonian(n, J, model, bonds=None):
    """H = sum_b J_b (XX + YY [+ ZZ]).  `bonds` overrides the per-bond J."""
    d = 2 ** n
    H = np.zeros((d, d), dtype=complex)
    Js = bonds if bonds is not None else [J] * (n - 1)
    for b in range(n - 1):
        H += Js[b] * (op_at(n, b, X) @ op_at(n, b + 1, X)
                      + op_at(n, b, Y) @ op_at(n, b + 1, Y))
        if model == "heisenberg":
            H += Js[b] * (op_at(n, b, Z) @ op_at(n, b + 1, Z))
    return H


def liouvillian(n, J, gamma, model, bonds=None):
    d = 2 ** n
    H = hamiltonian(n, J, model, bonds)
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for s in range(n):
        Zs = op_at(n, s, Z)
        L += gamma * (np.kron(Zs, Zs.T) - np.kron(Id, Id))
    return L


def floor_modes(n, J, gamma, model, bonds=None, tol=TOL):
    """Eigenmodes at Re lambda = -2 gamma, each tagged pure-n_XY=1 or mixed."""
    d = 2 ** n
    vals, vecs = np.linalg.eig(liouvillian(n, J, gamma, model, bonds))
    out = []
    for idx in np.where(np.abs(vals.real + 2 * gamma) < tol)[0]:
        w = np.abs(vecs[:, idx].reshape(d, d)) ** 2
        tot = w.sum()
        if tot < 1e-14:
            continue
        hist = {}
        for i in range(d):
            for j in range(d):
                if w[i, j] > 1e-12 * tot:
                    k = bin(i ^ j).count("1")
                    hist[k] = hist.get(k, 0.0) + w[i, j] / tot
        out.append((vals[idx], set(hist) == {1}, dict(sorted(hist.items()))))
    return out


def extras(n, J, gamma, model, bonds=None, tol=TOL):
    return [m for m in floor_modes(n, J, gamma, model, bonds, tol) if not m[1]]


def se_zz_diagonal(n, J):
    """ZZ energies of the single-excitation states, in the BUILDER's J."""
    out = []
    for i in range(n):
        e = 0.0
        for b in range(n - 1):
            same = (b != i) and (b + 1 != i)
            e += J * (1.0 if same else -1.0)
        out.append(e)
    return np.array(out)


failures = []


def gate(name, ok, detail):
    log(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


# Work in the proof's convention throughout: builder J = 0.5 is the proof's J = 1.
JP = 1.0           # the proof's J
JB = JP / 2        # the builder's J
GAMMA = 0.05

log("=" * 74)
log("THE {0,2} EXTRAS: A LOW-N THRESHOLD WHOSE RUNG MOVES WITH THE MODEL")
log("=" * 74)
log()
log("G1 + G2  where the extras are, and what they are (gamma-swept)")
log("-" * 74)

expected = {("xy", 2), ("xy", 3), ("heisenberg", 2)}
g1_ok = True
for g in (0.02, 0.05, 0.20):
    found = set()
    for model in ("xy", "heisenberg"):
        for n in (2, 3, 4, 5):
            ex = extras(n, JB, g, model)
            ims = sorted({round(abs(m[0].imag), 6) for m in ex
                          if abs(m[0].imag) > 1e-9})
            log(f"    gamma={g:<5} N={n} {model:10s} extras={len(ex):2d}  |Im| = {ims}")
            if ex:
                found.add((model, n))
                for _, _, hist in ex:
                    shape_ok = (set(hist) == {0, 2}
                                and abs(hist[0] - 0.5) < 1e-6
                                and abs(hist[2] - 0.5) < 1e-6)
                    if not shape_ok:
                        gate("G2 histogram", False, f"{model} N={n}: {hist}")
    if found != expected:
        g1_ok = False
        gate(f"G1 threshold gamma={g}", False,
             f"expected {sorted(expected)}, found {sorted(found)}")
gate("G1 threshold", g1_ok,
     f"extras exactly at {sorted(expected)} at every gamma swept (Q = 50, 20, 5). "
     f"SCOPE: this is a HIGH-Q statement. Each block has its own Q* below which the pair "
     f"coalesces off the floor (Q*(2)=1, Q*(3)=sqrt2), so at gamma=0.8 the N=3 extras are "
     f"already gone and at gamma=1.5 all of them are. The counts are functions of (N, Q), "
     f"not of N alone.")
if not any(f.startswith("G2") for f in failures):
    gate("G2 histogram", True, "every extra has n_XY histogram {0: 1/2, 2: 1/2}")

log()
log("G3  the closed form sqrt(c^2 - (2 gamma)^2), gamma-swept")
log("-" * 74)
for model, n in sorted(expected):
    E1 = 2 * JP * np.cos(np.pi / (n + 1))
    c = 2 * JP if n == 2 else E1
    for g in (0.05, 0.10, 0.20):
        ex = extras(n, JB, g, model)
        meas = max(abs(m[0].imag) for m in ex)
        pred = np.sqrt(c ** 2 - (2 * g) ** 2)
        gate(f"G3 {model} N={n} gamma={g}", abs(meas - pred) < 1e-9,
             f"c={c:.6f}  measured {meas:.9f}  closed form {pred:.9f}")

log()
log("G4  THE CROSSOVER at N = 2: above E1 exactly for Q > 2/sqrt(3)")
log("-" * 74)
Q_STAR = 2 / np.sqrt(3)
E1_N2 = 2 * JP * np.cos(np.pi / 3)          # = JP
log(f"    E1(N=2) = {E1_N2:.6f}, predicted crossover at Q = {Q_STAR:.6f}"
    f"  (gamma = {JP / Q_STAR:.6f})")
log("    Below Q = 1 the extras do not merely drop below E1: the two-level")
log("    block's square root turns real, the pair coalesces and LEAVES the")
log("    floor. So the statement is true below Q = 2/sqrt(3) for two different")
log("    reasons, and the floor max is E1 in both.")
cross_ok = True
for Q in (0.5, 0.8, 0.95, 1.05, 1.1, 1.15, 1.16, 1.2, 1.5, 2.0, 5.0, 20.0):
    g = JP / Q
    fm = floor_modes(2, JB, g, "xy")
    ex = [m for m in fm if not m[1]]
    # The claim is about the WHOLE floor, so measure the whole floor.
    meas = max((abs(m[0].imag) for m in fm), default=0.0)
    ex_meas = max((abs(m[0].imag) for m in ex), default=0.0)
    above = meas > E1_N2 + 1e-12
    want = Q > Q_STAR
    ok = (above == want)
    if not ok:
        cross_ok = False
    where = "none" if not ex else f"{ex_meas:.6f}"
    log(f"    Q={Q:<6} gamma={g:.6f}  floor dim={len(fm):2d}  extras |Im|={where:<9}  "
        f"FLOOR max|Im|={meas:.6f}  "
        f"{'ABOVE' if above else 'at or below'} E1   predicted "
        f"{'ABOVE' if want else 'at or below'}   {'ok' if ok else 'MISMATCH'}")
gate("G4 crossover", cross_ok,
     "the facing flips at Q = 2/sqrt(3) and nowhere else in the swept range")
dim_hi = len(floor_modes(2, JB, JP / 1.05, "xy"))
dim_lo = len(floor_modes(2, JB, JP / 0.8, "xy"))
gate("G4 coalescence", len(extras(2, JB, JP / 0.8, "xy")) == 0
     and len(extras(2, JB, JP / 1.05, "xy")) > 0
     and (dim_hi, dim_lo) == (10, 8),
     f"the extras leave the floor entirely below Q = 1 (the Q*(2) EP): "
     f"floor dim {dim_hi} above it, {dim_lo} below")

log()
log("G5  THE N=3 EXTRAS NEED UNIFORM BONDS: the drift law, which is tolerance-free")
log("-" * 74)
log("    Detuning one bond by delta moves the would-be extras OFF the floor. Any")
log("    fixed membership tolerance is therefore a window on delta, and quoting")
log("    'the extras are gone' without the law just moves the blindness one decade")
log("    down. So measure the drift itself and fit its power.")
drift = []
for delta in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    bonds = [JB, JB * (1 + delta)]
    # Select the {0,2} modes with a generous window, then read how far off the
    # floor they actually sit.  This does not depend on TOL.
    near = [m for m in floor_modes(3, JB, GAMMA, "xy", bonds=bonds, tol=1e-3)
            if not m[1]]
    off = max((abs(m[0].real + 2 * GAMMA) for m in near), default=float("nan"))
    drift.append((delta, off))
    log(f"    delta={delta:<8} max |Re + 2gamma| among the {{0,2}} modes = {off:.3e}"
        f"   (admitted by a tolerance of {off:.0e} or looser)")
# Fit the power: log(off) vs log(delta) over the clean decades.
xs = np.log10([d for d, o in drift[:4]])
ys = np.log10([o for d, o in drift[:4]])
power = float(np.polyfit(xs, ys, 1)[0])
gate("G5 drift law", abs(power - 2.0) < 0.05,
     f"the drift is O(delta^{power:.3f}), i.e. QUADRATIC in the detuning")
gate("G5 detuning", len(extras(3, JB, GAMMA, "xy",
                               bonds=[JB, JB * (1 + 1e-4)], tol=TOL)) == 0,
     f"so at delta=1e-4 the drift is {drift[2][1]:.1e} and this script's "
     f"TOL={TOL:.0e} correctly rejects all four")
gate("G5 tolerance window", len(extras(3, JB, GAMMA, "xy",
                                       bonds=[JB, JB * (1 + 1e-4)], tol=1e-8)) > 0,
     "while the earlier TOL=1e-8 admitted them, which is why THAT gate could not "
     f"fail. Ours fails in turn at delta below ~1e-5, where the drift falls under "
     f"{TOL:.0e}: the quadratic law above is the statement without a window")

log()
log("G6  THE MECHANISM: the single-excitation ZZ diagonal")
log("-" * 74)
for n in (2, 3, 4):
    dg = se_zz_diagonal(n, JB)
    const = bool(np.allclose(dg, dg[0], atol=1e-12))
    log(f"    N={n}: ZZ diagonal on the SE sector (proof units) = "
        f"{np.round(dg / JB * JP / 2, 6)}   {'CONSTANT' if const else 'not constant'}")
    if n == 2 and not const:
        gate("G6 N=2 constant", False, "expected constant")
    if n == 3 and const:
        gate("G6 N=3 split", False, "expected a non-constant diagonal")
gate("G6 diagonal", not any(f.startswith("G6") for f in failures),
     "constant at N=2 (ZZ cannot break the block), not constant at N=3 (it does). NOTE: not constant is not the same as non-degenerate -- (0,-J,0) repeats 0, and N=4 reads (J/2,-J/2,-J/2,J/2). Constancy is what matters and what is tested.")

xy2 = max(abs(m[0].imag) for m in extras(2, JB, GAMMA, "xy"))
he2 = max(abs(m[0].imag) for m in extras(2, JB, GAMMA, "heisenberg"))
gate("G6 consequence (derivation check, follows from G6 diagonal)", abs(xy2 - he2) < 1e-12,
     f"Heisenberg N=2 |Im|={he2:.12f} equals XY N=2 |Im|={xy2:.12f}")
gate("G6 one rung", len(extras(3, JB, GAMMA, "heisenberg")) == 0,
     "and Heisenberg has no extras at N=3: the ZZ term costs exactly one rung")

log()
log("G7  BOND DISORDER: the maximality survives it (only the extras do not)")
log("-" * 74)
rng = np.random.default_rng(20260802)
disorder_ok = True
for n in (4, 5):
    for trial in range(3):
        # Random positive bond weights, in the BUILDER's J.
        wb = rng.uniform(0.3, 1.7, size=n - 1) * JB
        # Single-particle hopping matrix: XX+YY hops with amplitude 2 per bond.
        hop = np.zeros((n, n))
        for b in range(n - 1):
            hop[b, b + 1] = hop[b + 1, b] = 2 * wb[b]
        top = max(abs(np.linalg.eigvalsh(hop)))
        fm = floor_modes(n, JB, GAMMA, "xy", bonds=list(wb))
        meas = max((abs(m[0].imag) for m in fm), default=0.0)
        ok = abs(meas - top) < 1e-8
        if not ok:
            disorder_ok = False
        log(f"    N={n} trial {trial}: floor max|Im|={meas:.9f}  "
            f"weighted band top={top:.9f}  {'ok' if ok else 'MISMATCH'}")
gate("G7 bond disorder", disorder_ok,
     "floor max|Im| = the weighted hopping matrix's band top at N=4,5, "
     "so section 4.3's fence is on the extras only")

log()
log("=" * 74)
if failures:
    log(f"FAILED: {failures}")
    raise SystemExit(1)
log("ALL GATES PASS")
log()
log("Reading. The {0,2} extras are not an N = 3 accident. They are the low-N end")
log("of a closure condition on the (n,n) sectors, met by XY at N = 2 and 3 and by")
log("Heisenberg at N = 2 only, the difference being the SE ZZ diagonal (G6). And")
log("at N = 2 they sit ABOVE the band edge whenever Q > 2/sqrt(3) (G4), so the")
log("chain has a second band-edge counterexample beside the ring's 4-cycle, in")
log("that regime. Uniform bonds are part of the statement, not a convenience (G5).")
log("=" * 74)
