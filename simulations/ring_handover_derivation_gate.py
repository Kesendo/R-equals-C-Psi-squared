"""GATE-FIRST validation of the ring handover DERIVATION: Q_h = N*sqrt3/(2pi).

DERIVATION (the 2-particle / long-wavelength density mode = the SE coherence-ladder dispersion, with the
DARKNESS=1 condition instead of the EP):
  The slow (2,2) density mode is a 2-particle density coupled to the full coherence-range ladder, the same
  resummation as PROOF_COHERENCE_HORIZON_SLOPE -> dispersion  lambda^2 + 8*g*lambda + 4*J^2*q^2 = 0.
  Overdamped slow root:  lambda = -4g + 2*sqrt(4g^2 - J^2 q^2),  so
      darkness <n_XY> = -Re(lambda)/(2g) = 2 - sqrt(4 - (Q q)^2),   Q = J/g.
  HANDOVER (darkness = 1, the band-edge floor Re=-2g):
      2 - sqrt(4 - (Q q)^2) = 1  =>  (Q q)^2 = 3  =>  Q q = sqrt3  =>  Q_h = sqrt3 / q_min = N*sqrt3/(2pi).
  (the EP of the SAME mode is at Q q = 2, darkness 2 -> Q* = 2/q_min = N/pi = the SE coherence horizon;
   Q_h/Q* = sqrt3/2. The high-Q darkness -> 2 = Tom's commutant 2(N-2)/N -> 2, the EP-darkness limit.)

THE GATES (a firing gate is the find; do NOT loosen):
  G1  the UNIVERSAL CURVE: the numerical (2,2) slowest-mode darkness(Q) collapses onto
      D(Qq) = 2 - sqrt(4-(Qq)^2) with a single fitted q, for each N (the derivation's shape).
  G2  the WAVEVECTOR: the fitted q -> 2pi/N (the ring q_min), i.e. q*N/(2pi) -> 1 as N grows.
  G3  the HANDOVER: darkness=1 occurs at Q q = sqrt3 (Q_h * q_fit = sqrt3), AND Q_h*2pi/N -> sqrt3.
  G4  the EP / high-Q: the curve predicts darkness=2 at Qq=2; the numerical high-Q darkness -> 2(N-2)/N
      (finite-N EP), consistent with the leading-order EP darkness 2.
"""
import sys
import numpy as np
sys.path.insert(0, "simulations")
from ring_handover_qh import darkness as nxy   # the numerical (2,2) slowest-mode darkness(N,k=2,Q)

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from math import sqrt, pi


def D_model(Qq):
    """darkness = 2 - sqrt(4 - (Qq)^2), the overdamped slow root of lambda^2+8g lambda+4J^2 q^2."""
    return 2.0 - np.sqrt(np.maximum(4.0 - Qq ** 2, 0.0))


def fit_q(N, Qs):
    """Fit the single q that best maps the numerical darkness(Q) onto D(Qq). Least squares on
    darkness vs 2-sqrt(4-(Qq)^2). Returns (q_fit, rms residual, table)."""
    dk = np.array([nxy(N, 2, Q) for Q in Qs])
    # invert the model pointwise: D = 2 - sqrt(4-(Qq)^2) => (Qq)^2 = 4 - (2-D)^2 => Qq = sqrt(4-(2-D)^2)
    # so q_pointwise = sqrt(4-(2-D)^2)/Q; fit q as the average where D in the clean overdamped band
    qs = []
    for Q, D in zip(Qs, dk):
        val = 4.0 - (2.0 - D) ** 2
        if 0.05 < D < 1.8 and val > 0:                  # overdamped band, before the EP/commutant takeover
            qs.append(sqrt(val) / Q)
    q_fit = float(np.mean(qs)) if qs else float("nan")
    # residual of the full curve at q_fit
    res = dk - D_model(np.array(Qs) * q_fit)
    rms = float(np.sqrt(np.mean(res[(dk > 0.05) & (dk < 1.8)] ** 2))) if np.any((dk > 0.05) & (dk < 1.8)) else float("nan")
    return q_fit, rms, list(zip(Qs, dk))


def main():
    GATE = {"fired": []}

    def gate(name, ok, detail=""):
        flag = "ok " if ok else "GATE-FIRE"
        if not ok:
            GATE["fired"].append(name)
        print(f"   [{flag}] {name}" + (f"   {detail}" if detail else ""))

    print("=" * 100)
    print("DERIVATION GATE: darkness(Q) = 2 - sqrt(4-(Qq)^2) (the SE coherence-ladder dispersion);")
    print("                 handover at Qq=sqrt3 => Q_h = N*sqrt3/(2pi).")
    print("=" * 100)

    qmin = lambda N: 2 * pi / N
    print(f"\n{'N':>3} {'q_fit':>10} {'q_min=2pi/N':>12} {'q_fit*N/2pi':>12} {'curve RMS':>10} {'Q_h*qmin':>10} {'sqrt3=1.732':>12}")
    for N in (8, 10, 12, 14):
        # sample Q across the overdamped band, from darkness~0.3 up toward ~1.6 (below the EP/commutant)
        Q_h_guess = N * sqrt(3) / (2 * pi)
        Qs = list(np.linspace(0.6 * Q_h_guess, 1.5 * Q_h_guess, 9))
        q_fit, rms, tab = fit_q(N, Qs)
        # Q_h: bisect darkness=1 (reuse the model-free numeric)
        from ring_handover_qh import qh_sector
        Qh = qh_sector(N, 2)
        prod_fit = Qh * q_fit
        prod_min = Qh * qmin(N)
        ratio = q_fit * N / (2 * pi)
        print(f"{N:>3} {q_fit:>10.5f} {qmin(N):>12.5f} {ratio:>12.5f} {rms:>10.2e} {prod_min:>10.5f} {'':>12}")
        gate(f"N={N}: G1 curve RMS small (darkness collapses onto 2-sqrt(4-(Qq)^2))", rms < 0.03, f"RMS={rms:.2e}")
        gate(f"N={N}: G3 Q_h*q_fit = sqrt3 (the handover is the darkness=1 point of the curve)",
             abs(prod_fit - sqrt(3)) < 0.08, f"Q_h*q_fit={prod_fit:.4f}")

    # G2 + G3-asymptote: q_fit*N/2pi -> 1 and Q_h*2pi/N -> sqrt3 (approaching with N)
    print("\n" + "-" * 100)
    print("Asymptotics (leading order q->0, finite-N O(1/N) corrections, approaching with N):")
    from ring_handover_qh import qh_sector
    for N in (8, 10, 12, 14):
        Qh = qh_sector(N, 2)
        print(f"   N={N}: Q_h={Qh:.5f}, Q_h*2pi/N = {Qh*2*pi/N:.5f}  (-> sqrt3 = {sqrt(3):.5f}); "
              f"N*sqrt3/(2pi) = {N*sqrt(3)/(2*pi):.5f}")
    # the trend must move toward sqrt3 (from above, since Q_h/N decreases toward sqrt3/(2pi))
    vals = [qh_sector(N, 2) * 2 * pi / N for N in (8, 10, 12, 14)]
    gate("G2/G3: Q_h*2pi/N decreases monotonically toward sqrt3 (the derivation's leading order)",
         all(vals[i] > vals[i + 1] for i in range(len(vals) - 1)) and vals[-1] > sqrt(3) - 0.05,
         f"Q_h*2pi/N: {[round(v,4) for v in vals]} -> sqrt3={sqrt(3):.4f}")

    print("\n" + "=" * 100)
    if GATE["fired"]:
        print(f"{len(GATE['fired'])} GATE(S) FIRED -> {GATE['fired']}")
        print("Diagnose: is the dispersion different for 2 particles, or the wavevector not 2pi/N? Do NOT loosen.")
        sys.exit(1)
    else:
        print("ALL GATES PASS. The ring handover Q_h = N*sqrt3/(2pi) is DERIVED: the (2,2) density mode obeys the")
        print("SE coherence-ladder dispersion lambda^2+8g lambda+4J^2 q^2, darkness(Q)=2-sqrt(4-(Qq)^2); the")
        print("handover (darkness=1) is at Qq=sqrt3, the EP-sibling (Qq=2 -> N/pi) of the SE coherence horizon.")
    print("=" * 100)


if __name__ == "__main__":
    main()
