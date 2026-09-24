"""Crossing-taxonomy probe: five scalar readouts in two books.

Resolves the February crossing numbers (t = 0.652 / 0.773 / 1.437 at
gamma = 0.05; K = 0.033 / 0.039 / 0.072) of experiments/CROSSING_TAXONOMY.md,
experiments/OBSERVER_DEPENDENT_CROSSING.md and experiments/METRIC_DISCRIMINATION.md.

The trajectory: Bell+ under local Z-dephasing is J-independent (Heisenberg
eigenstate); rho(t) = diag(1/2,0,0,1/2) + (f/2)(|00><11| + h.c.). Psi = f/3.

Two named evolution-law families for f:
  CLEAN     df/dt = -4*gamma*f            (standard Lindblad; f = e^(-4*gamma*t))
  FEEDBACK  df/dt = -4*gamma*C(f)*f       (the retired February delta_calc tool's
            gamma_eff = gamma_base*C(t), disclosed in OBSERVER_DEPENDENT_CROSSING
            section 3.2; for C = f it integrates to f(t) = 1/(1+4*gamma*t))

The five bridges C(f), as the retired tool's bridge function (bridge_type; source kept
outside the repo) defines them, read on Bell+:
  concurrence   C = f                       (Wootters concurrence of rho(t))
  mutual_info   C = [2 - h2((1+f)/2)]/2     (von Neumann MI in bits / its t=0 value)
  correlation   C = min(1, 2(P_AB - P_A P_B)) = min(1, 1/2 + f^2), an excess
                purity doubled and capped; it is 1 for f >= 1/sqrt(2), so through
                every crossing here.
  mutual_purity C = sqrt(P_A P_B) = 1/2      (the geometric mean of the subsystem
                purities, not their product 1/4)
  overlap       C = |Tr(rho_A rho_B)|^2 = 1/4 (not fidelity with the initial state,
                which would start at 1 and decay)

Crossing: C(f)*f/3 = 1/4. correlation is identical in both models (C = 1 through the crossing
makes the feedback inert there). mutual_purity and overlap have CPsi(0) < 1/4: never cross.

Expected output (verified 2026-07-21):
  bridge         clean K      clean t(g=.05)  feedback K   feedback t   Feb doc t
  mutual_info    0.029657     0.5931          0.032645     0.6529       0.652
  concurrence    0.035960     0.7192          0.038675     0.7735       0.773
  correlation    0.071921     1.4384          0.071921     1.4384       1.437
  (concurrence feedback crossing is analytic: t = (2/sqrt(3)-1)/(4*gamma))

K = gamma*t is exact for these five equations in each book, with the bridge
fixed on this Hamiltonian-dead Bell+ family while gamma is swept. The feedback
law is state-dependent/nonlinear for mutual information and concurrence;
correlation (on its cap, which lasts past the crossing), mutual purity
and overlap give constant-rate linear scalar decay.
The feedback family is not one linear Lindblad generator. A scalar readout
choice does not establish a physical measurement event or experienced time.
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
    "correlation": lambda f: min(1.0, 0.5 + f * f),
    "mutual_purity": lambda f: 0.5,
    "overlap": lambda f: 0.25,
}


def cpsi(bridge, f, bridges=None):
    return (BRIDGES if bridges is None else bridges)[bridge](f) * f / 3.0


def clean_crossing(bridge, bridges=BRIDGES):
    """f = e^(-4*gamma*t); solve C(f)*f/3 = 1/4 for f, then t."""
    g = lambda f: cpsi(bridge, f, bridges) - 0.25
    if g(1.0) <= 0:
        return None, None
    f_cross = brentq(g, 1e-9, 1.0, xtol=1e-14)
    K = -np.log(f_cross) / 4.0
    return K, K / GAMMA


def feedback_crossing(bridge, bridges=None, *, rtol=1e-12, atol=1e-14):
    """df/dt = -4*gamma*C(f)*f (the February tool's gamma_eff = gamma*C)."""
    if cpsi(bridge, 1.0, bridges) <= 0.25:
        return None, None
    hit = lambda t, y: cpsi(bridge, y[0], bridges) - 0.25
    hit.terminal = True
    hit.direction = -1
    sol = solve_ivp(lambda t, y: [-4.0 * GAMMA * (BRIDGES if bridges is None else bridges)[bridge](y[0]) * y[0]],
                    (0.0, 100.0), [F_START], events=hit,
                    rtol=rtol, atol=atol, dense_output=True)
    t_cross = sol.t_events[0][0]
    return GAMMA * t_cross, t_cross


def crossings_for_book(book, bridges=BRIDGES):
    """Read every supplied scalar bridge through one named evolution law."""
    if book not in ("clean", "feedback"):
        raise ValueError(f"Unknown crossing book: {book!r}")
    crossing = clean_crossing if book == "clean" else feedback_crossing
    return {name: crossing(name, bridges=bridges) for name in bridges}


def main():
    print("Two books: Bell+ local Z-dephasing, gamma = 0.05, Psi = f/3")
    clean = crossings_for_book("clean")
    feedback = crossings_for_book("feedback")
    print(f"{'bridge':<14s} {'CPsi(0)':>8s} {'clean K':>10s} {'clean t':>9s} "
          f"{'feedbk K':>10s} {'feedbk t':>9s}")
    for name in BRIDGES:
        c0 = cpsi(name, 1.0)
        cK, ct = clean[name]
        fK, ft = feedback[name]
        fmt = lambda x, w: ("-" * 5).rjust(w) if x is None else f"{x:{w}.6f}"
        print(f"{name:<14s} {c0:8.4f} {fmt(cK,10)} "
              f"{('never' if ct is None else f'{ct:9.4f}')} {fmt(fK,10)} "
              f"{('never' if ft is None else f'{ft:9.4f}')}")
    t_analytic = (2.0 / np.sqrt(3.0) - 1.0) / (4.0 * GAMMA)
    print(f"\nconcurrence feedback crossing, analytic (2/sqrt(3)-1)/(4*gamma) "
          f"= {t_analytic:.6f}")
    print("correlation is identical in both models (C = 1 through the crossing: feedback inert there).")
    print("K = gamma*t is exact for the five listed equations in each fixed-bridge Bell+ book.")
    print("Feedback: mutual information and concurrence are state-dependent/nonlinear;")
    print("correlation (on its cap through the crossing) and the two constant bridges give constant-rate linear scalar decay.")
    print("Readout crossings do not establish physical observer events or experienced time.")


if __name__ == "__main__":
    main()
