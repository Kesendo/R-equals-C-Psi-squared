"""F130 hardware infeasibility: the beat protocol is dead from below.

THE QUESTION. F130 says two equal-level triples never hybridize under dephasing:
the effective coupling on the (1,2) Liouville coherence block is 5.9e-5 q^4 for the
collision pair vs 1.9e-2 q^2 for a generic pair (gamma units, q = J/gamma; committed
gate simulations/f130_time_domain_decoupling.py, G2). The naive hardware protocol
("the beat") would prepare coherence weight in one Slater-lift multiplet and watch
the q^2 beat appear (generic control) or stay absent (collision pair). This script
pins WHY that protocol is not flyable: even an everywhere-OPTIMISTIC upper bound on
the beat signal sits ~2 orders below shot noise at 8192 shots.

Part A (Trotter): does the first-order Trotter step (odd/even Givens layering) kill
the exact degeneracy? Answer: no faster than 0.13 theta^3 per step; the collision
survives Trotterization at safe angles. (So the protocol dies on AMPLITUDE, not on
the degeneracy.)

Part B (the amplitude budget, optimistic in every factor): after M steps at hop
angle theta = q * gamma_step per step, the transferred-coherence amplitude is upper-
bounded by

    A(q, M) <= 0.5 * (M * g_step) * exp(-6 M gamma_step) * (1 - p2)^CX(M),

with g_step = c2 q^2 gamma_step (generic; c2 = 1.9e-2 from the committed gate),
0.5 = the maximal initial (1,2) coherence weight, exp(-6 gamma_step) = the rung -6
decay per step (A = -2 n_diff, n_diff = 3), and CX(M) the two-qubit count of prep
network + M steps + readout network (the 2/3-magnon networks the Ramsey design also
needs; counted at their theoretical minimum here). OPTIMISM LIST (every item favors
the protocol): linear-in-t growth M*g_step with no detuning suppression (G5 says the
self-blocks split at q^2, which only reduces transfer); no T1/T2* idle decay; no
readout error; no Trotter-splitting dephasing; unit fringe contrast; the coupling law
extrapolated at face value to q = 2 (fit window q <= 0.32). The scan maximizes over
theta <= 0.32 (safe Trotter, part A), gamma_step free, M free.

VERDICT CRITERION (pre-registered): the protocol is infeasible if even the optimistic
bound gives SNR = A * sqrt(shots) < 2 at 8192 shots everywhere on the scan; the
discriminating measurement (generic beat present vs collision beat absent) needs the
GENERIC arm visible at >= 5 sigma, so SNR < 2 on the optimistic bound is dead with
margin. The protected arm is smaller by another (q^2 * 5.9e-5/1.9e-2), i.e. invisible
a fortiori.

Doc: experiments/F130_HW_INFEASIBILITY.md. Constants: p2 = 0.3% (median Heron r2
2q error, Confirmations 21 noise model), shots 8192.
"""

import numpy as np

N = 11
n = N + 1
TAU = (1, 2, 10)
SIGMA = (3, 5, 6)
CTRL = (3, 4, 9)

C2_GENERIC = 1.9e-2   # gamma-units q^2 coefficient, committed gate G2
C4_PROTECTED = 5.9e-5  # gamma-units q^4 coefficient, committed gate G2
P2 = 0.003
SHOTS = 8192
THETA_MAX = 0.32


def one_particle_step(theta, nn=N):
    U = np.eye(nn, dtype=complex)
    for parity in (0, 1):
        L = np.eye(nn, dtype=complex)
        for i in range(parity, nn - 1, 2):
            B = np.eye(nn, dtype=complex)
            B[i, i] = B[i + 1, i + 1] = np.cos(theta)
            B[i, i + 1] = B[i + 1, i] = 1j * np.sin(theta)
            L = B @ L
        U = L @ U
    return U


def floquet_phases(theta):
    ev = np.linalg.eigvals(one_particle_step(theta))
    return np.sort(np.angle(ev))[::-1]


def level(triple, phases):
    return sum(phases[k - 1] for k in triple)


def part_a():
    print("Part A: Trotter splitting of the exact degeneracy (per step)")
    print("theta    split(tau,sigma)      detune(ctrl)     split/theta^3")
    for theta in (0.02, 0.04, 0.08, 0.16, 0.32):
        ph = floquet_phases(theta)
        split = abs(level(TAU, ph) - level(SIGMA, ph))
        det = abs(level(TAU, ph) - level(CTRL, ph))
        print(f"{theta:5.2f}   {split:.6e}      {det:.6e}     {split / theta ** 3:8.4f}")
    print()


def cx_count(m_steps):
    """Theoretical minimum: 3-magnon prep network (55 Givens x 2 CX) + seed (5) +
    M steps (10 bonds x 2 CX) + inverse network + un-seed. Cross-sector readout
    would cost MORE (it does not exist yet); counted as one more network here."""
    return 5 + 110 + m_steps * 20 + 110 + 110


def best_point(qs):
    """Maximize the bound over gamma_step AND M at each q. The bound
    A = 0.5 M c2 q^2 gamma e^(-6 M gamma) (1-p2)^CX has its interior optimum at
    gamma = 1/(6M); the Trotter-safety constraint theta = q*gamma <= THETA_MAX
    caps gamma at THETA_MAX/q, so gamma* = min(THETA_MAX/q, 1/(6M))."""
    best = (0.0, None)
    for q in qs:
        for m in range(1, 61):
            gamma_step = min(THETA_MAX / q, 1.0 / (6 * m))
            g_step = C2_GENERIC * q ** 2 * gamma_step
            amp = 0.5 * m * g_step * np.exp(-6 * m * gamma_step) * (1 - P2) ** cx_count(m)
            snr = amp * np.sqrt(SHOTS)
            if snr > best[0]:
                best = (snr, (q, gamma_step, m, amp))
    return best


def part_b():
    """The bound amp ~ q * exp(-1.92/q) at M ~ 1 GROWS monotonically in q, so there
    is no interior maximum: the scan must stop at the law's validity edge, and that
    edge is physics, not a knob. The q^2 coupling law is second-order perturbation
    theory in J/gamma, fit over q <= 0.32 and meaningless for q >> 1: there the
    chain is hopping-dominated, the dephasing-induced beat B is no longer the
    object, and the surviving first-order physics is the Ramsey fringe (the pivot).
    Verdict domain = q <= 1 (the perturbative edge); q = 2 is printed separately as
    a generous over-extrapolation that STILL fails. gamma_step and M are truly
    maximized over (see best_point); nothing is pinned except theta <= THETA_MAX."""
    print("Part B: optimistic beat-amplitude bound inside the law's validity domain")
    snr, (q, gs, m, amp) = best_point((0.05, 0.1, 0.2, 0.32, 0.5, 1.0))
    print(f"  best in-domain point (q <= 1): q = {q}, gamma_step = {gs:.3f}, M = {m}")
    print(f"  optimistic amplitude = {amp:.3e}, SNR at {SHOTS} shots = {snr:.4f} "
          f"({5 / snr:.0f}x short of 5 sigma)")
    snr2, (q2, gs2, m2, amp2) = best_point((2.0,))
    print(f"  over-extrapolation q = 2 (past the fit window q <= 0.32): "
          f"amplitude = {amp2:.3e}, SNR = {snr2:.3f} ({5 / snr2:.0f}x short) -- still dead")
    ratio = C4_PROTECTED * q2 ** 2 / C2_GENERIC
    print(f"  protected arm smaller by another x{ratio:.1e} at q = 2 (invisible a fortiori)")
    print(f"  VERDICT: {'INFEASIBLE (optimistic SNR < 2 everywhere in and beyond the domain)' if max(snr, snr2) < 2 else 'feasible?!'}")


if __name__ == "__main__":
    part_a()
    part_b()
