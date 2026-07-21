"""Reproduce the IBM Run 3 (Torino q80, 2026-03-18) crossing predictions.

Pins the convention question for experiments/IBM_RUN3_PALINDROME.md: both
locked predictions (9.47 us from the stale T2* = 11.0 us, 15.01 us from the
same-day Ramsey T2* = 17.36 us) fall out of the STANDARD decay convention,
|rho_01|(t) = (1/2) e^{-t/T2*} with the T1 population term included
(repo rate convention gamma = 1/(2 T2*), coherence e^{-2 gamma t}).

The doc's earlier background sentence "gamma = 1/(pi T2*)" was a prose slip
(that pi belongs to the Ramsey linewidth Delta_f = 1/(pi T2*), not to the
decay rate); under it the same-day prediction would be 49.9 us against the
measured 15.29 us.

Model: single qubit prepared in |+>, free decoherence.
  |rho_01|(t) = (1/2) e^{-t/T2*}
  p_1(t) = (1/2) e^{-t/T1},  p_0 = 1 - p_1        (relaxation toward |0>)
  purity P = p_0^2 + p_1^2 + 2 |rho_01|^2
  L1/(d-1) = 2 |rho_01|
  CPsi = P * L1/(d-1),  crossing at CPsi = 1/4.

Run: python simulations/ibm_run3_crossing_check.py
"""
import numpy as np
from scipy.optimize import brentq

T1 = 143.1  # us, day-of calibration


def cpsi(t, T2s, rate_of_T2s=lambda T2s: 1.0 / T2s):
    rate = rate_of_T2s(T2s)
    c = 0.5 * np.exp(-rate * t)
    p1 = 0.5 * np.exp(-t / T1)
    p0 = 1.0 - p1
    purity = p0 ** 2 + p1 ** 2 + 2 * c ** 2
    return purity * 2 * c


def crossing(T2s, rate_of_T2s=lambda T2s: 1.0 / T2s):
    return brentq(lambda t: cpsi(t, T2s, rate_of_T2s) - 0.25, 0.1, 300.0)


if __name__ == "__main__":
    print("IBM Run 3 crossing predictions, standard convention "
          "(|rho01| ~ e^{-t/T2*}, T1 = %.1f us):" % T1)
    rows = [
        ("stale T2* (March 12)", 11.0, 9.47),
        ("same-day Ramsey T2* (March 18)", 17.36, 15.01),
    ]
    ok = True
    for label, T2s, doc in rows:
        t = crossing(T2s)
        match = abs(t - doc) < 0.005
        ok &= match
        print("  %-32s T2*=%6.2f us  t* = %6.2f us   doc: %5.2f   %s"
              % (label, T2s, t, doc, "MATCH" if match else "MISMATCH"))

    t_pi = crossing(17.36, lambda T2s: 1.0 / (np.pi * T2s))
    print("\nCounterfactual pi convention (gamma = 1/(pi T2*)): "
          "t* = %.2f us, measured 15.29 us -- refuted." % t_pi)
    print("\nMeasured (8-point tomography): t* = 15.29 us; "
          "same-day prediction deviation 1.9%.")
    print("\nAll doc predictions reproduce: %s" % ok)
    assert ok
