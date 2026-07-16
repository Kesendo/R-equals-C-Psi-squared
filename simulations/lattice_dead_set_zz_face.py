"""Gate for the zz face of the dead-set law (experiments/LATTICE_DEAD_SET_ZZ_FACE.md):
how the dead set lifts as the interaction grows from zero.

Context: F132 (experiments/LATTICE_DEAD_SET_RULE.md) proves the free-world
dead set and gates its zz boundary only at zz = 0.7 (V and F die, layer K
persists, the lone fermionic survivor is the conserved parity Z^N). The open
question carried by the handover: jump or crossover at small zz? Answer,
gated here: BOTH, cleanly split by level.

  T1  ZZ = 0 REFERENCE. The F132 law is exact at the generic h direction used
      throughout (the parent gate's +-0.4 scan), all three configs.

  T2  THE SET JUMPS. For every zz in {1e-4, 1e-3, 1e-2, 0.1, 0.7} the alive
      set is already the FULL zz-large set: everything K-readable revives
      except the conserved parity Z^N at odd N (which starts at 0 for the
      population pair and stays there at every zz); nothing newly dies. The
      symmetry breaking is derived (any zz != 0 breaks V and F exactly);
      that no intermediate set appears is the gated observation.

  T3  THE MAGNITUDES CROSS OVER, WITH EXACT INTEGER ORDERS. Every revived
      readout obeys the revival-order law

          max_t |<O>(t)| ~ zz^m,   m = min |d(O) - d0| / 2  over populated
                                       prep degrees d0 of the same parity,

      i.e. m quartic zz vertices bridge the Majorana-degree gap (a
      contributing vertex moves d by exactly +-2, the shared-generator count
      must be odd; parity is preserved). Checked per readout: the
      local log-log slopes on (1e-4 -> 1e-3) and (1e-3 -> 1e-2) match m to
      < 0.02. At N = 5 s = 3 coherence the orders split three ways: the 255
      d = 2 mod 4 and the 240 odd-d strings at distance 2 from the coherence
      degree ride zz^1; the 20 F-killed strings (d in {1, 9}, distance 4
      from d = 5) ride zz^2. The revived pool is exactly accounted:
      20 + 255 + 240 = 515.

  T4  THE JUMP IS DETECTION-RELATIVE, GATED. At zz = 1e-5 the m = 2 strings
      sit near 7e-11 < the 1e-10 aliveness threshold while every m = 1
      string is already alive: the still-dead set at zz = 1e-5 is exactly
      {Z^N} + the 20 m = 2 strings. So any apparent growth of the SET with
      zz is a threshold artifact of the zz^m magnitudes; the identically-
      zero property itself is lost at every zz != 0.

Setup identical to the F132 gate (one-sided watching, generic fixed h
direction +-0.4 * [1.0, -0.6, 0.3, 0.5, ...], RK4 window [0, 2], dt = 0.05),
zz added as zz * sum ZZ over the chain bonds. Machinery imported from the
committed simulations/lattice_dead_set_h_zero.py. Runtime ~3 min.
"""
import sys

import numpy as np

sys.path.insert(0, "simulations")
from lattice_dead_set_h_zero import (  # noqa: E402
    Config, block_readable, majorana_degree, sz,
)

ALIVE_TOL = 1e-10
FAILURES = []


def check(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILURES.append(name)


def make_config(N, seed, coh_kind, gammas, zz, label):
    cfg = Config(N, seed, coh_kind, gammas, label)
    for b in range(N - 1):
        cfg.H0 = cfg.H0 + zz * (cfg.op(sz, b) @ cfg.op(sz, b + 1))
    return cfg


def populated_degrees(cfg):
    """Majorana degrees carrying a nonzero Pauli coefficient of rho(0),
    including the identity's degree 0."""
    degs = {0}
    for name, O in cfg.readouts:
        if abs(np.trace(cfg.rho0 @ O)) > 1e-12:
            degs.add(majorana_degree(name))
    return degs


def revival_order(name, degs):
    """m = min |d - d0| / 2 over populated degrees of the same parity."""
    d = majorana_degree(name)
    gaps = [abs(d - d0) // 2 for d0 in degs if (d - d0) % 2 == 0]
    return min(gaps) if gaps else None


class Scan:
    """One config swept over the zz values, magnitudes cached."""

    ZZ = [0.0, 1e-4, 1e-3, 1e-2, 0.1, 0.7]

    def __init__(self, N, seed, coh_kind, gammas, label, extra_zz=()):
        self.label = label
        self.N = N
        hdir = ([1.0, -0.6, 0.3] + [0.5] * N)[:N]
        self.hp = [0.4 * v for v in hdir]
        self.hm = [-0.4 * v for v in hdir]
        self.mags = {}
        self.ref = None
        for zz in list(self.ZZ) + list(extra_zz):
            cfg = make_config(N, seed, coh_kind, gammas, zz, label)
            if self.ref is None:
                self.ref = cfg
            mp = cfg.alive_magnitudes(self.hp)
            mm = cfg.alive_magnitudes(self.hm)
            self.mags[zz] = {n: max(mp[n], mm[n]) for n in mp}
        self.law = self.ref.law_alive
        self.a0 = self.alive(0.0)
        self.k_dead = set()
        for name, _ in self.ref.readouts:
            w = sum(1 for ch in name if ch in "XY")
            k_any = block_readable(N, w, self.ref.sup_pop) or \
                (self.ref.coh and block_readable(N, w, self.ref.sup_coh))
            if not k_any:
                self.k_dead.add(name)
        self.pool = {n for n, _ in self.ref.readouts} - self.k_dead - self.a0
        self.zn = "Z" * N

    def alive(self, zz):
        return {n for n, m in self.mags[zz].items() if m > ALIVE_TOL}


def main():
    scans = [
        Scan(3, 1, None, [0.0, 0.0, 0.5], "N=3 s=1 population"),
        Scan(4, 1, "real", [0.0] * 3 + [0.5], "N=4 s=1 real coherence"),
        Scan(5, 3, "real", [0.0] * 4 + [0.5], "N=5 s=3 real coherence",
             extra_zz=(1e-5,)),
    ]

    print("T1  zz = 0 reference: the F132 law is exact at the scan direction")
    for s in scans:
        check(f"{s.label}: law at zz=0", s.a0 == s.law,
              f"alive {len(s.a0)}/{len(s.ref.readouts)}, law {len(s.law)}")

    print("T2  the set jumps: every zz != 0 revives the whole pool except "
          "the conserved parity Z^N (odd N)")
    for s in scans:
        survivor = {s.zn} if s.N % 2 == 1 else set()
        expected = s.a0 | (s.pool - survivor)
        ok = True
        detail_counts = []
        for zz in Scan.ZZ[1:]:
            a = s.alive(zz)
            ok = ok and (a == expected)
            detail_counts.append(len(a))
        check(f"{s.label}: full jump at every zz",
              ok and all(s.mags[zz][n] < ALIVE_TOL
                         for zz in Scan.ZZ[1:] for n in survivor),
              f"alive at zz=1e-4..0.7: {detail_counts} (expected "
              f"{len(expected)} each); never revived beyond K: "
              f"{sorted(s.pool - (s.alive(0.7) - s.a0)) or 'none'}")

    print("T3  the revival-order law: max|<O>| ~ zz^m with m = min |d-d0|/2 "
          "over same-parity populated prep degrees, per readout")
    for s in scans:
        degs = populated_degrees(s.ref)
        revived = s.alive(0.7) - s.a0
        worst = 0.0
        orders = {}
        for n in revived:
            m_pred = revival_order(n, degs)
            for z1, z2 in ((1e-4, 1e-3), (1e-3, 1e-2)):
                slope = np.log10(s.mags[z2][n]) - np.log10(s.mags[z1][n])
                worst = max(worst, abs(slope - m_pred))
            orders[m_pred] = orders.get(m_pred, 0) + 1
        check(f"{s.label}: per-readout slopes match m",
              worst < 0.02,
              f"populated degrees {sorted(degs)}; revived by order "
              f"{dict(sorted(orders.items()))}; worst |slope - m| = {worst:.4f}")
    s5 = scans[2]
    degs5 = populated_degrees(s5.ref)
    revived5 = s5.alive(0.7) - s5.a0
    m2 = {n for n in revived5 if revival_order(n, degs5) == 2}
    check("N=5 s=3: the m=2 set is exactly the 20 F-killed strings",
          len(m2) == 20 and all(majorana_degree(n) in (1, 9) for n in m2),
          f"{len(m2)} strings, degrees "
          f"{sorted({majorana_degree(n) for n in m2})} = the a_k and their "
          f"Gamma-duals, two quartic vertices from the coherence degree 5")
    m1 = revived5 - m2
    n_pop = sum(1 for n in m1 if majorana_degree(n) % 4 == 2)
    n_coh = sum(1 for n in m1 if majorana_degree(n) % 2 == 1)
    check("N=5 s=3: the m=1 set splits 255 (d=2 mod 4) + 240 (odd d)",
          len(m1) == 495 and n_pop == 255 and n_coh == 240
          and n_pop + n_coh == len(m1),
          f"{len(m1)} m=1 strings = {n_pop} population-channel (d=2 mod 4) "
          f"+ {n_coh} coherence-channel (d in {{3, 7}})")

    print("T4  the jump is detection-relative: at zz = 1e-5 exactly the m=2 "
          "strings (and Z^N) sit below the threshold")
    still_dead = s5.pool - (s5.alive(1e-5) - s5.a0)
    peak_m2 = max(s5.mags[1e-5][n] for n in m2)
    check("N=5 s=3 at zz=1e-5: still dead = {Z^N} + the 20 m=2 strings",
          still_dead == m2 | {s5.zn} and peak_m2 < ALIVE_TOL,
          f"still dead {len(still_dead)}; m=2 peak magnitude {peak_m2:.1e} "
          f"< {ALIVE_TOL:g}: apparent set growth with zz is a threshold "
          f"artifact of the zz^m magnitudes")

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: the dead set's zz lift is a jump in the SET (any "
          "zz != 0 breaks V and F exactly; only the conserved parity Z^N "
          "stays silent at odd N) and a crossover in the MAGNITUDES, with "
          "exact integer orders m = min |d - d0|/2: the interaction revives "
          "each silence at the perturbative order set by its Majorana-degree "
          "distance from the prep.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
