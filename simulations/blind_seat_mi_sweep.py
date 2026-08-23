"""Section 7's prediction, run: the site-by-site gamma sweep on a blind state.

`experiments/MEDIATOR_NOISE_GATE_LEVEL_THREE.md` swept gamma site by site on an
eleven-site chain and read the mediator seat as unremarkable, fourth of the nine
interior sites.  `experiments/THE_BLIND_SITE.md` section 7 argued that the seat's
node structure is not inert but invisible to that sweep's preparation, a
Bell-on-vacuum state living in popcount {0, 2} where the blind dimension is 0 at
every seat, and predicted that the same sweep on a reflection-odd
single-excitation preparation would read exactly zero at the mediator.

This script runs it.  Six parts:

  gate      reproduce the committed eleven-row table, so the runner is the same
            instrument before it is pointed anywhere new
  sweep     the same sweep on both preparations under ONE functional
  algebra   the exact operator facts, including the one that says what the
            mechanism is NOT, and the one trajectory reading beside them
  support   what the residual at the centre is actually made of
  zeno      why the sign of that residual flips
  converge  which of these numbers is a law and which is a window

The prediction HOLDS: with nothing on the other sites the centre's response is a
null difference between two identical computations.  Everything after that is
what the run added and the prediction did not have.

WHAT THE MECHANISM IS, AND WHAT IT IS NOT.  It is a SUPPORT statement, not a
parity statement.  The single-excitation action of the jump is
Z_k = I - 2|k><k|, and every reflection-odd state has a node at the centre, so
Z_5 acts as the identity there: max|Z_5 P_odd - P_odd| = 0.0 exactly.  For any
other seat |k> is not in the odd span and the jump carries odd states out of it.
Reflection PARITY is not what protects.  In the configurations where the
blindness is exact, the ones with gamma on the centre, the profile is itself
reflection-symmetric and the Liouvillian is exactly reflection-covariant
(measured 0.0), while the preparation rho0 = |psi><psi| with R psi = -psi is
Ad_R-EVEN: the grading is conserved and the state sits in its trivial half.  The
off-centre rows of the sweep are not reflection-covariant at all, which is a
second reason covariance cannot be the explanation.  That distinction is the one
`experiments/THE_BLIND_SITE.md` section 3 was written to make, and this script
measures both sides of it in the `algebra` part.

THE RESIDUAL IS NOT A LEAK MAGNITUDE.  At the sweep's own uniform baseline of
0.05 the centre reads 0.081 percent rather than zero, and it is tempting to call
that the state leaking out of the odd subspace.  The `support` part refuses the
temptation: the ends alone give -0.26 percent, the interior alone +0.31, the
centre's own two neighbours +2.10.  The 0.081 is what survives a cancellation
between contributions of both signs and an order of magnitude larger.

Usage: python simulations/blind_seat_mi_sweep.py
           [gate | sweep | algebra | support | zeno | converge | all]
"""
import sys

import numpy as np

sys.path.insert(0, "simulations")
from bridge_sector import (SectorPropagator, bell_on_vacuum, build_h, build_mask,
                           mutual_information, sector_basis)

N = 11
A = [0, 1, 2, 3, 4]
B = [6, 7, 8, 9, 10]
CENTRE = 5
BONDS = [(i, i + 1, 1.0) for i in range(N - 1)]
T_MAX = 20.0

# The gate reproduces a committed run and must use that run's step.  Everything
# else uses the finer step, because at dt = 0.05 the PEAK of the pure arm exceeds
# its own exact upper bound of 2 bits by 1.3e-3, the same order as the difference
# the centre's span is computed from.  See `converge` part (a), which prints it.
DT_GATE = 0.05
DT = 0.02

# experiments/MEDIATOR_NOISE_GATE_LEVEL_THREE.md, the eleven-row table
PUBLISHED = [39.67, 36.33, 28.04, 19.47, 12.30,
             15.57, 23.96, 18.04, 9.97, 4.59, 5.02]


def se_basis(n):
    """The single-excitation basis, one state per site."""
    states = sorted(1 << (n - 1 - a) for a in range(n))
    return states, {s: i for i, s in enumerate(states)}


def reflection_odd(n, states, index):
    """(|1_0> - |1_{n-1}>)/sqrt(2): the excitation on the first site minus the
    excitation on the last, as a density matrix in the single-excitation block."""
    psi = np.zeros(len(states), dtype=complex)
    psi[index[1 << (n - 1)]] = 1.0 / np.sqrt(2.0)
    psi[index[1]] = -1.0 / np.sqrt(2.0)
    return np.outer(psi, psi.conj())


def reflection_odd_projector(n, index):
    """Projector onto the reflection-odd part of the single-excitation block."""
    p = np.zeros((n, n), dtype=complex)
    for a in range(n // 2):
        v = np.zeros(n, dtype=complex)
        v[index[1 << (n - 1 - a)]] = 1.0 / np.sqrt(2.0)
        v[index[1 << a]] = -1.0 / np.sqrt(2.0)
        p += np.outer(v, v.conj())
    return p


def site_reversal(n, index):
    """The site reversal a -> n-1-a on the single-excitation basis."""
    r = np.zeros((n, n))
    for a in range(n):
        r[index[1 << a], index[1 << (n - 1 - a)]] = 1.0
    return r


def mi_curve(states, index, rho0, gammas, dt=DT, t_max=T_MAX):
    """I(A:B) after every RK4 step, including t = 0."""
    prop = SectorPropagator(N, BONDS, gammas, states, index)
    vals = []
    prop.propagate_every_step(
        rho0.copy(), t_max, prop.stable_dt(dt),
        lambda t, rho: vals.append(mutual_information(rho, N, states, A, B)))
    return np.array(vals)


def peak(curve):
    return float(curve.max())


def mean(curve):
    return float(curve.mean())


def end(curve):
    return float(curve[-1])


def profiles(seat, base, support=None):
    """The sweep's two gamma profiles: the seat at 0, then at 0.5.

    `support`, when given, is the set of OTHER sites that carry the baseline;
    the default is every site except the swept one.
    """
    others = range(N) if support is None else support
    lo = [0.0] * N
    hi = [0.0] * N
    for k in others:
        if k != seat:
            lo[k] = base
            hi[k] = base
    hi[seat] = 0.5
    return lo, hi


def span(states, index, rho0, seat, functional, base=0.05, dt=DT, t_max=T_MAX,
         support=None):
    lo_g, hi_g = profiles(seat, base, support)
    lo = functional(mi_curve(states, index, rho0, lo_g, dt, t_max))
    hi = functional(mi_curve(states, index, rho0, hi_g, dt, t_max))
    return lo, hi, 100.0 * (lo - hi) / lo


def setups():
    """The two preparations, each with its own basis."""
    b_states, b_index = sector_basis(N)
    s_states, s_index = se_basis(N)
    return ((b_states, b_index, bell_on_vacuum(N, 0, 1, b_states, b_index)),
            (s_states, s_index, reflection_odd(N, s_states, s_index)))


def run_gate():
    """The instrument check: reproduce the committed table before using it."""
    (b_states, b_index, bell), _ = setups()
    print(f"GATE: Bell-on-vacuum, peak I(A:B) over t in [0, {T_MAX:.0f}],")
    print(f"baseline gamma = 0.05, RK4 dt = {DT_GATE} (the committed run's step),")
    print("against the eleven rows of MEDIATOR_NOISE_GATE_LEVEL_THREE.md.")
    print()
    print("A failure here would be a row off by more than the 0.005 that two")
    print("printed decimals can hide; that is what this gate can catch.")
    print()
    print(f"{'seat':>4} {'peak @0':>10} {'peak @0.5':>10} {'span%':>8} "
          f"{'published':>10} {'delta':>8}")
    worst = 0.0
    for seat in range(N):
        lo, hi, s = span(b_states, b_index, bell, seat, peak, dt=DT_GATE)
        d = s - PUBLISHED[seat]
        worst = max(worst, abs(d))
        print(f"{seat:>4} {lo:10.6f} {hi:10.6f} {s:8.2f} "
              f"{PUBLISHED[seat]:10.2f} {d:8.2f}")
    print()
    print(f"worst deviation from the published two decimals: {worst:.3f}")


def run_sweep():
    """The same sweep on both preparations, under ONE functional."""
    (b_states, b_index, bell), (s_states, s_index, odd) = setups()
    print(f"SWEEP: N = {N}, uniform chain J = 1, A = {A}, B = {B},")
    print(f"centre = {CENTRE}. I(A:B) in log base 2 over t in [0, {T_MAX:.0f}],")
    print(f"RK4 dt = {DT}, measured after every step; the listed seat is swept")
    print("gamma 0 -> 0.5 while every other site holds the baseline 0.05.")
    print()
    print("THE COMPARISON IS LIKE FOR LIKE.  The committed sweep's PEAK")
    print("functional cannot carry the reflection-odd state: that state starts")
    print("at the maximum I(A:B) and never rises above it, so its peak is the")
    print("t = 0 value at every seat and every gamma, printed below to show the")
    print("saturation.  With no dephasing at all it does not even fall: p_centre")
    print("is 0 for all t, so p_A = p_B = 1/2 and I(A:B) is exactly 2.")
    print("The window MEAN is therefore run on BOTH preparations, so the two")
    print("columns differ in the preparation and in nothing else.")
    print()
    print(f"{'seat':>4} | {'Bell mean%':>11} | {'odd peak':>9} "
          f"{'odd mean%':>13} {'odd end%':>12}")
    bell_spans, odd_spans = [], []
    for seat in range(N):
        _, _, b_m = span(b_states, b_index, bell, seat, mean)
        p_lo, _, _ = span(s_states, s_index, odd, seat, peak)
        _, _, o_m = span(s_states, s_index, odd, seat, mean)
        _, _, o_e = span(s_states, s_index, odd, seat, end)
        bell_spans.append(b_m)
        odd_spans.append(o_m)
        mark = "  <- centre" if seat == CENTRE else ""
        print(f"{seat:>4} | {b_m:11.3f} | {p_lo:9.6f} {o_m:13.6e} "
              f"{o_e:12.6e}{mark}")
    b_other = [v for j, v in enumerate(bell_spans) if j != CENTRE]
    o_other = [v for j, v in enumerate(odd_spans) if j != CENTRE]
    print()
    print(f"Bell:  centre {bell_spans[CENTRE]:.3f} %, "
          f"others {min(b_other):.2f} to {max(b_other):.2f} %")
    print(f"odd:   centre {odd_spans[CENTRE]:.4f} %, "
          f"others {min(o_other):.2f} to {max(o_other):.2f} %")
    print(f"across preparations at the centre: "
          f"{bell_spans[CENTRE] / odd_spans[CENTRE]:.0f} to 1")
    print(f"within the odd column, smallest other seat to centre: "
          f"{min(o_other) / odd_spans[CENTRE]:.0f} to 1")
    print()
    print("The odd column is exactly reflection-symmetric, seat j against seat")
    print(f"{N - 1} - j.  That is FORCED, not observed, and not by any single")
    print("run being covariant: R carries the seat-j profile to the")
    print(f"seat-({N - 1} - j) profile, rho0 is Ad_R-even (see `algebra`), and R")
    print("swaps A with B while I(A:B) is symmetric, so the two runs are")
    print("Ad_R-conjugate TO EACH OTHER.")
    print("The Bell column is not symmetric because its preparation sits on")
    print("sites 0 and 1; its spread is distance to that pair.")


def run_algebra():
    """The exact operator facts, including the one that says what this is NOT."""
    _, (s_states, s_index, odd) = setups()
    p = reflection_odd_projector(N, s_index)
    r = site_reversal(N, s_index)
    h = build_h(N, BONDS, s_states, s_index)
    print("ALGEBRA: three exact facts, no propagation, no eigensolver.")
    print()
    print("(1) THE MECHANISM IS SUPPORT.  Z_k = I - 2|k><k| on the")
    print("single-excitation block; every reflection-odd state has a node at the")
    print("centre, so Z_centre acts as the identity on that subspace.")
    for k in (CENTRE, CENTRE - 1, 0):
        z = np.eye(N)
        z[s_index[1 << (N - 1 - k)], s_index[1 << (N - 1 - k)]] = -1.0
        print(f"    max|Z_{k} P_odd - P_odd| = {np.abs(z @ p - p).max():.3e}")
    print()
    print("(2) THE MECHANISM IS NOT PARITY.  With gamma on the CENTRE the")
    print("profile [b, ..., b, gamma_centre, b, ..., b] is itself")
    print("reflection-symmetric, so R permutes the jump set and the Liouvillian")
    print("is exactly reflection-covariant.  Those are the configurations where")
    print("the blindness is exact.  The sweep's off-centre rows are NOT")
    print("reflection-covariant, which is a second reason covariance cannot be")
    print("the explanation.")
    print(f"    R is an involution: max|R R - I| = "
          f"{np.abs(r @ r - np.eye(N)).max():.1e}")
    print(f"    max|R H R - H| = {np.abs(r @ h @ r - h).max():.1e}")
    for base, centre in ((0.05, 0.0), (0.05, 0.5), (0.2, 0.5)):
        g = [base] * N
        g[CENTRE] = centre
        m = build_mask(N, g, s_states)
        print(f"    baseline {base}, centre {centre}: "
              f"max|R mask R - mask| = {np.abs(r @ m @ r - m).max():.1e}")
    print()
    print("(3) AND THE PREPARATION IS Ad_R-EVEN, so the conserved reflection")
    print("parity gives it nothing: R psi = -psi makes |psi><psi| invariant.")
    print(f"    max|R rho0 R - rho0| = {np.abs(r @ odd @ r - odd).max():.1e}")
    print()
    print("Together: reflection covariance holds and is NOT the protecting")
    print("thing, which is what THE_BLIND_SITE.md section 3 says in prose.  The")
    print("protecting thing is that the jump at the centre is the identity on")
    print("the subspace the state lives in.")
    print()
    print("(4) The trajectory consequence, the one thing here that is")
    print("propagated: Tr(P_odd rho), with P_odd of dimension "
          f"{int(round(np.trace(p).real))} of {N}.")
    print()
    print(f"{'baseline':>10} {'centre':>8} {'t = 0':>16} {'t = 5':>16} "
          f"{'t = 20':>16}")
    for base in (0.0, 0.05):
        for centre in (0.0, 0.5):
            g = [base] * N
            g[CENTRE] = centre
            prop = SectorPropagator(N, BONDS, g, s_states, s_index)
            rec = []
            prop.propagate_every_step(
                odd.copy(), T_MAX, prop.stable_dt(DT),
                lambda t, rho: rec.append((t, float(np.trace(p @ rho).real))))
            got = [next(v for t, v in rec if abs(t - w) < 1e-9)
                   for w in (0.0, 5.0, T_MAX)]
            print(f"{base:10.3f} {centre:8.2f} " + " ".join(f"{v:16.12f}"
                                                            for v in got))
    print()
    print("With nothing on the other sites the population is 1.000000000000")
    print("across the whole window whether or not the centre carries gamma, so")
    print("the centre is a node for the whole trajectory and not only at t = 0.")
    print("With the baseline it is down to 0.469 by t = 20, on its way to the")
    print(f"unbiased share {N // 2}/{N} = {(N // 2) / N:.6f}.")


def run_support():
    """What the residual at the centre is made of.  It is not one leak."""
    _, (s_states, s_index, odd) = setups()
    print("SUPPORT DECOMPOSITION: baseline 0.05 on the NAMED SET only, the")
    print("centre swept 0 -> 0.5, window mean of I(A:B).")
    print()
    print("It is tempting to read the centre's 0.081 percent as the size of the")
    print("leak out of the odd subspace.  This part refuses that reading.")
    print()
    named = [("all ten (the sweep itself)", [j for j in range(N) if j != CENTRE]),
             ("the two ends {0, 10}", [0, N - 1]),
             ("the interior, ends excluded", [j for j in range(1, N - 1)
                                              if j != CENTRE]),
             ("one end {0}", [0]),
             ("the centre's neighbours {4, 6}", [CENTRE - 1, CENTRE + 1])]
    print(f"{'baseline support':<30} {'centre span%':>14}")
    for label, s in named:
        _, _, v = span(s_states, s_index, odd, CENTRE, mean, support=s)
        print(f"{label:<30} {v:+14.4f}")
    print()
    print("The contributions carry BOTH signs and the neighbours alone are an")
    print("order of magnitude larger than the total.  So 0.081 percent is a")
    print("cancellation residue, not a leak magnitude, and the only thing the")
    print("run establishes about it is that it vanishes when the rest of the")
    print("chain is quiet.  Why the cancellation lands where it does is open.")


def run_zeno():
    """Why the residual changes sign.  The two-point span hides a crossover."""
    _, (s_states, s_index, odd) = setups()
    print("ZENO CROSSOVER: the sweep compares gamma_centre = 0 against 0.5 and")
    print("reports one number.  Scanning gamma_centre instead shows that number")
    print("is a two-point sample of a NON-MONOTONE function.")
    print()
    print("p_centre is the window-mean occupation of the centre site, the")
    print("bottleneck every A-to-B path crosses.")
    print()
    for base in (0.05, 0.2):
        print(f"baseline {base}")
        print(f"{'gamma_centre':>14} {'span% vs 0':>14} {'mean p_centre':>16}")
        g0 = [base] * N
        g0[CENTRE] = 0.0
        ref = mean(mi_curve(s_states, s_index, odd, g0))
        for gc in (0.0, 0.5, 2.0, 5.0, 20.0, 100.0):
            g = [base] * N
            g[CENTRE] = gc
            v = mean(mi_curve(s_states, s_index, odd, g))
            prop = SectorPropagator(N, BONDS, g, s_states, s_index)
            occ = []
            prop.propagate_every_step(
                odd.copy(), T_MAX, prop.stable_dt(DT),
                lambda t, rho: occ.append(
                    float(rho[s_index[1 << (N - 1 - CENTRE)],
                              s_index[1 << (N - 1 - CENTRE)]].real)))
            print(f"{gc:14.1f} {100.0 * (ref - v) / ref:14.3f} "
                  f"{float(np.mean(occ)):16.6f}")
        print()
    print("Strong dephasing at the centre BLOCKS transport through it, the")
    print("occupation there falls, and the preparation keeps more of its initial")
    print("maximal I(A:B).  That is the Zeno / ENAQT regime the repo already")
    print("names (docs/GLOSSARY.md, D06_SPECTRAL_GAP, PROOF_ABSORPTION_THEOREM).")
    print("The crossover exists at BOTH baselines; gamma_centre = 0.5 simply")
    print("sits on different sides of it.  So the sign is not unexplained, and")
    print("it is a counterexample, on THIS preparation, to the reading that the")
    print("mediator's own noise always harms (experiments/GAMMA_CONTROL.md).")


def run_converge():
    """Which of these numbers is a law and which is a window."""
    _, (s_states, s_index, odd) = setups()
    print("CONVERGENCE.")
    print()
    print("(a) THE EXACT BOUND, and why the step was refined.  For a single")
    print("excitation in a pure state, I(A:B) = h(p_A) + h(p_B) - h(p_centre)")
    print("<= 2 bits, and with no dephasing anywhere it is exactly 2 for all t")
    print("because p_centre stays 0.  So anything above 2 below is integrator")
    print("error and nothing else.")
    print(f"{'dt':>10} {'max I(A:B)':>18} {'excess over 2':>16}")
    zero = [0.0] * N
    for dt in (0.05, 0.02, 0.01, 0.005):
        m = peak(mi_curve(s_states, s_index, odd, zero, dt=dt))
        print(f"{dt:10.4f} {m:18.9f} {m - 2.0:+16.2e}")
    print()
    print(f"The committed run's dt = {DT_GATE} violates the bound by 1.3e-3, the")
    print("same order as the difference the centre's span is computed from, so")
    print(f"everything but the gate is run at dt = {DT}.")
    print()
    print("(b) THE SPANS ARE dt-CONVERGED ANYWAY, because the error is common")
    print("to both arms and cancels in the difference.")
    print(f"{'dt':>10} {'base 0.001':>14} {'base 0.01':>14} {'base 0.05':>14} "
          f"{'base 0.2':>14}")
    for dt in (0.05, 0.02, 0.01, 0.005):
        row = [span(s_states, s_index, odd, CENTRE, mean, base=b, dt=dt)[2]
               for b in (0.001, 0.01, 0.05, 0.2)]
        print(f"{dt:10.4f} " + " ".join(f"{v:14.4e}" for v in row))
    print()
    print("The baseline-0.001 column moves by 6 percent between the coarsest and")
    print("the finest step; the other three are stable in the fourth digit, and")
    print("the negative sign at baseline 0.2 survives every step.")
    print()
    print("(c) THE WINDOW IS NOT CONVERGED, and cannot be.  Past the transient")
    print("every profile decays toward the same limit, so a longer window")
    print("dilutes every span alike.  The committed PEAK is window-converged")
    print("(that page reports the same 15.566 percent at t_max = 20, 40 and 80);")
    print("the window MEAN is not.  What holds still is the RATIO.")
    print(f"{'t_max':>10} {'centre span%':>16} {'min other':>12} "
          f"{'max other':>12} {'min other / centre':>20}")
    for t_max in (20.0, 40.0, 80.0):
        c = span(s_states, s_index, odd, CENTRE, mean, t_max=t_max)[2]
        others = [span(s_states, s_index, odd, j, mean, t_max=t_max)[2]
                  for j in range(N) if j != CENTRE]
        print(f"{t_max:10.0f} {c:16.6e} {min(others):12.4f} "
              f"{max(others):12.4f} {min(others) / c:20.1f}")
    print()
    print("So the quotable number is the RATIO, stable to about one percent")
    print("across a factor of four in the window.  The absolute percentages")
    print(f"elsewhere in this run are read at t_max = {T_MAX:.0f} and only there.")
    print()
    print("(d) THE ZERO-BASELINE ANSWER is not an integrator floor and should")
    print("not be called one: with nothing on the other sites the two arms are")
    print("the SAME computation, because Z_centre is the identity on the")
    print("subspace, so what is printed is a null difference between identical")
    print("runs.  It certifies exact blindness and says nothing about accuracy.")
    print(f"{'t_max':>10} {'baseline 0':>16} {'baseline 0.2':>16}")
    for t_max in (20.0, 40.0, 80.0):
        z = span(s_states, s_index, odd, CENTRE, mean, base=0.0, t_max=t_max)[2]
        n = span(s_states, s_index, odd, CENTRE, mean, base=0.2, t_max=t_max)[2]
        print(f"{t_max:10.0f} {z:16.3e} {n:16.6e}")


PARTS = [("gate", run_gate), ("sweep", run_sweep), ("algebra", run_algebra),
         ("support", run_support), ("zeno", run_zeno),
         ("converge", run_converge)]


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    names = [n for n, _ in PARTS]
    if mode not in names + ["all"]:
        print(f"unknown part {mode!r}; known: {' | '.join(names)} | all")
        return
    for i, (name, fn) in enumerate(PARTS):
        if mode in (name, "all"):
            fn()
            if mode == "all" and i < len(PARTS) - 1:
                print()
                print("-" * 72)
                print()


if __name__ == "__main__":
    main()
