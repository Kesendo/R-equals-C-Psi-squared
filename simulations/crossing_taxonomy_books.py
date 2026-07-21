"""Crossing-taxonomy books probe: the five bridges, two models, all crossings.

Resolves the February crossing numbers (t = 0.652 / 0.773 / 1.437 at
gamma = 0.05; K = 0.033 / 0.039 / 0.072) of experiments/CROSSING_TAXONOMY.md,
experiments/OBSERVER_DEPENDENT_CROSSING.md and experiments/METRIC_DISCRIMINATION.md.

The trajectory: Bell+ under local Z-dephasing is J-independent (Heisenberg
eigenstate); rho(t) = diag(1/2,0,0,1/2) + (f/2)(|00><11| + h.c.). Psi = f/3.

Two models for f:
  CLEAN     df/dt = -4*gamma*f            (standard Lindblad; f = e^(-4*gamma*t))
  FEEDBACK  df/dt = -4*gamma*C(f)*f       (the retired February delta_calc tool's
            gamma_eff = gamma_base*C(t), disclosed in OBSERVER_DEPENDENT_CROSSING
            sections 3.2/5.1; for C = f it integrates to f(t) = 1/(1+4*gamma*t))

The five bridges C(f), as pinned from the surviving prose + reproduced behavior
(the tool itself was never committed):
  concurrence   C = f                       (Wootters concurrence of rho(t))
  mutual_info   C = [2 - h2((1+f)/2)]/2     (von Neumann MI in bits / its t=0 value)
  correlation   C = 1 exactly               (connected Z-basis correlation
                <Z1Z2> - <Z1><Z2> = 1 for all t: diagonals are frozen under pure
                dephasing. The old prose label "excess purity" does not reproduce
                the recorded C = 1.000: excess purity decays.)
  mutual_purity C = 1/2                     (single-subsystem purity; the prose
                "product of subsystem purities" would be 1/4)
  overlap       C = 1/4                     (= 1/d; the prose "fidelity with the
                initial state" would start at 1 and decay)

Crossing: C(f)*f/3 = 1/4. correlation is identical in both models (C = 1 makes
the feedback inert). mutual_purity and overlap have CPsi(0) < 1/4: never cross.

Expected output (verified 2026-07-21):
  bridge         clean K      clean t(g=.05)  feedback K   feedback t   Feb doc t
  mutual_info    0.029657     0.5931          0.032645     0.6529       0.652
  concurrence    0.035960     0.7192          0.038675     0.7735       0.773
  correlation    0.071921     1.4384          0.071921     1.4384       1.437
  (concurrence feedback crossing is analytic: t = (2/sqrt(3)-1)/(4*gamma))

K-invariance is exact in BOTH models (f depends on tau = gamma*t only), which is
why the February gamma-sweeps saw constant K: the surviving Lindblad-scaling
result is model-independent.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

GAMMA = 0.05
F_START = 1.0


def h2(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * np.log2(p) - (1 - p) * np.log2(1 - p)


BRIDGES = {
    "mutual_info": lambda f: (2.0 - h2((1.0 + f) / 2.0)) / 2.0,
    "concurrence": lambda f: f,
    "correlation": lambda f: 1.0,
    "mutual_purity": lambda f: 0.5,
    "overlap": lambda f: 0.25,
}


def cpsi(bridge, f):
    return BRIDGES[bridge](f) * f / 3.0


def clean_crossing(bridge):
    """f = e^(-4*gamma*t); solve C(f)*f/3 = 1/4 for f, then t."""
    g = lambda f: cpsi(bridge, f) - 0.25
    if g(1.0) <= 0:
        return None, None
    f_cross = brentq(g, 1e-9, 1.0, xtol=1e-14)
    K = -np.log(f_cross) / 4.0
    return K, K / GAMMA


def feedback_crossing(bridge):
    """df/dt = -4*gamma*C(f)*f (the February tool's gamma_eff = gamma*C)."""
    if cpsi(bridge, 1.0) <= 0.25:
        return None, None
    hit = lambda t, y: cpsi(bridge, y[0]) - 0.25
    hit.terminal = True
    hit.direction = -1
    sol = solve_ivp(lambda t, y: [-4.0 * GAMMA * BRIDGES[bridge](y[0]) * y[0]],
                    (0.0, 100.0), [F_START], events=hit,
                    rtol=1e-12, atol=1e-14, dense_output=True)
    t_cross = sol.t_events[0][0]
    return GAMMA * t_cross, t_cross


def main():
    print("Crossing-taxonomy books (Bell+ Z-dephasing, gamma = 0.05, Psi = f/3)")
    print(f"{'bridge':<14s} {'CPsi(0)':>8s} {'clean K':>10s} {'clean t':>9s} "
          f"{'feedbk K':>10s} {'feedbk t':>9s}")
    for name in BRIDGES:
        c0 = cpsi(name, 1.0)
        cK, ct = clean_crossing(name)
        fK, ft = feedback_crossing(name)
        fmt = lambda x, w: ("-" * 5).rjust(w) if x is None else f"{x:{w}.6f}"
        print(f"{name:<14s} {c0:8.4f} {fmt(cK,10)} "
              f"{('never' if ct is None else f'{ct:9.4f}')} {fmt(fK,10)} "
              f"{('never' if ft is None else f'{ft:9.4f}')}")
    t_analytic = (2.0 / np.sqrt(3.0) - 1.0) / (4.0 * GAMMA)
    print(f"\nconcurrence feedback crossing, analytic (2/sqrt(3)-1)/(4*gamma) "
          f"= {t_analytic:.6f}")
    print("correlation is identical in both models (C = 1: feedback inert).")
    print("K-invariance is exact in both models: f depends on tau = gamma*t only.")


if __name__ == "__main__":
    main()
