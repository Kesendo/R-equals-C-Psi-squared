"""
The operator-feedback runs of February 2026, regenerated and read against exact propagation
===========================================================================================
experiments/OPERATOR_FEEDBACK.md reports runs of the retired delta_calc MCP tool (its source is kept
outside the repo and is transcribed here, not imported): simulate_dynamic_lindblad with
noise_type = operator_feedback. The tool

  - builds H = J * sum over the chain bonds of (XX + YY + ZZ) + h * sum_k X_k (the field along x);
  - builds the jumps from jump_operator: sigma_z on every site (the default, and the fall-through for
    a setting it does not know, unless it ends in '_pairs' or '_all', where an unknown letter becomes
    sigma_x), with O_int = X_0 X_1; one sigma_a x sigma_a on qubits 0 and 1 ('xx', 'yy', 'zz'), with
    O_int = that jump itself; sigma_a x sigma_a on ALL pairs i < j ('x_pairs', 'y_pairs', 'z_pairs'),
    with O_int = their mean; one global sigma_a^N ('x_all', ...);
  - sets the rate gamma_eff = max(0, gamma_0 * (1 - kappa * <O_int>)) at every step;
  - integrates by an Euler step of dt 0.01, then takes the Hermitian part, sets negative
    eigenvalues to zero and renormalizes the trace;
  - records every 10th step, rounded to six digits: the purity Tr(rho^2) of the full state, its
    mutual-purity bridge C = (prod_k Tr(rho_k^2))^(1/n), Psi = l1(rho)/(d - 1), <O_int>, and
    delta = purity minus (diag + offdiag * exp(-2 n gamma_0 t)), a fixed function of rho_0 and t.

The Claude Desktop log of the chat's MCP holds nine calls of 2026-02-20 (17:22 to 17:29 UTC) on
Bell+, and the one operator-feedback x_pairs call on GHZ_4 the same evening (19:07 UTC, h = 0.2,
kappa 0.99, t_max 30);
the transcription reproduces every <O_int> value the log shows of each (the log truncates each
response to its tail) and regenerates the §8.2 table. No GHZ call at the §8.1 settings is in the
log; the transcription regenerates the page's GHZ rows under the §8.2 settings with 'x_pairs'.

The numbers experiments/OPERATOR_FEEDBACK.md derives from runs are compared here, at the page's
precision, with the values computed here (check: a literal from the page against a computation); the
page's quotes of the log are compared with this file's transcription of the log (quote); statements
that compare two computations, or a computation with a bound from its error model, are claims
(claim); the §4 table and the Simulation Evidence runs are regenerated in
simulations/delta_calc_feedback_runs.py. Beside the tool's Euler numbers stand exact routes:
  - the identities behind the readings, compared exactly (== 0) on integer or exactly representable
    matrices (products of permutations, sign changes and integer sums); a CONTROL is a check through
    the same function that must come out nonzero: [H, X_0 X_1] = 0 for the field along x and not
    along z; on the X_0 X_1 = +1 sector YY = -ZZ; S = Sum_{i<j} X_i X_j commutes with the
    Heisenberg chain, the x field and every X_i X_j jump, so also with the nearest-neighbour jumps,
    while the nearest-neighbour sum does not commute with the bonds and S does not commute with the
    XY chain; Z_k X_i X_j Z_k = -X_i X_j for k in {i, j}, so the sigma_z dissipator maps S to -4 S;
  - propagation by expm of the Liouvillian (linear runs) and by DOP853 at rtol 1e-12, 1e-11 for the
    GHZ runs of §8.1 (feedback runs); each closed form is gated against propagation by a one-sided
    bound: at every rtol from 1e-5 to 1e-11 the deviation must stay within rtol (converges), which a
    closed form off by more than the tightest tolerance cannot do; where the deviation scales with
    rtol the printed ratios show it, and where it sits at rounding it shows no scaling. Where the closed
    form is for x = <O_int>, the dynamics closes on x (dx/dt is a function of x alone), so the Runge-Kutta
    stages, projected on O_int, are the same method applied to that scalar equation, at the steps the
    faster components set; at the default rtol 1e-12 the deviation then sits at rounding, and those
    gates are bounds of 100 eps. The other
    tolerances here (10 rtol, plus or minus 100 eps, 500 eps, and departures of more than 0.01 for the
    controls that must leave a closed form or zero) are one-sided bounds too, each named in its label.

Conventions (the tool's): Bell+ = (|00> + |11>)/sqrt(2); GHZ_N = (|0..0> + |1..1>)/sqrt(2); W_N the
even superposition of the single excitations; qubit 0 the most significant bit; gamma_0 = 0.1,
J = 1, h = 0.5, t_max 5 unless stated. Exits with status 1 if any check fails. The summary counts
the values printed on experiments/OPERATOR_FEEDBACK.md and on the pages whose numbers it quotes.

Import-inert; prints only.
"""

from functools import reduce
from math import comb
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
PAULI = {"x": X, "y": Y, "z": Z}
EPS = np.finfo(float).eps


def site(P, k, n):
    return reduce(np.kron, [P if i == k else I2 for i in range(n)])


def pair(P, i, j, n):
    return site(P, i, n) @ site(P, j, n)


def heisenberg(n, xy=False, ring=False):
    """The bond part with J = 1: an integer matrix."""
    H = np.zeros((2**n, 2**n), dtype=complex)
    bonds = [(i, (i + 1) % n) for i in range(n if ring else n - 1)]
    for i, j in bonds:
        for P in ((X, Y) if xy else (X, Y, Z)):
            H += pair(P, i, j, n)
    return H


def field(n, P=X):
    """The field part with h = 1: an integer matrix."""
    return sum(site(P, k, n) for k in range(n))


def hamiltonian(n, J=1.0, h=0.5, axis="x", xy=False, ring=False):
    return J * heisenberg(n, xy, ring) + h * field(n, PAULI[axis])


def jumps_and_observable(n, jump):
    """The tool's _build_jump_operators for operator_feedback, transcribed (the branches used here)."""
    if jump.endswith("_pairs"):
        L = [pair(PAULI[jump[0]], i, j, n) for i in range(n) for j in range(i + 1, n)]
        return L, sum(L) / len(L)
    if jump in ("xx", "yy", "zz"):
        L = pair(PAULI[jump[0]], 0, 1, n)
        return [L], L
    return [site(Z, k, n) for k in range(n)], pair(X, 0, 1, n)


def one_site(rho, n, k):
    t = rho.reshape([2] * (2 * n))
    for i in sorted((q for q in range(n) if q != k), reverse=True):
        t = np.trace(t, axis1=i, axis2=i + t.ndim // 2)
    return t.reshape(2, 2)


def purity(r):
    return float(np.real(np.trace(r @ r)))


def l1(r):
    return float(np.abs(r).sum() - np.abs(np.diag(r)).sum())


def mutual_purity(rho, n):
    return float(np.prod([purity(one_site(rho, n, k)) for k in range(n)]) ** (1.0 / n))


def cpsi(rho, n):
    """The tool's reading: mutual-purity bridge times l1/(d - 1) of the full state."""
    return mutual_purity(rho, n) * l1(rho) / (2**n - 1)


def concurrence(rho):
    YY = np.kron(Y, Y)
    R = rho @ YY @ rho.conj() @ YY
    ev = np.sort(np.real(np.sqrt(np.maximum(0, np.abs(np.linalg.eigvals(R))))))[::-1]
    return float(max(0.0, ev[0] - ev[1] - ev[2] - ev[3]))


def pure(amplitudes):
    v = np.asarray(amplitudes, dtype=complex)
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


def bell():
    return pure([1, 0, 0, 1])


def ghz(n):
    a = np.zeros(2**n)
    a[0] = a[-1] = 1
    return pure(a)


def w_state(n):
    return pure([1 if bin(i).count("1") == 1 else 0 for i in range(2**n)])


def tool_run(rho0, n, H, gamma0, kappa, jump, dt=0.01, t_max=5.0, clip=True, record_every=10):
    """The tool's Euler run, transcribed. Returns the recorded rows (t, purity, C, Psi, <O_int>,
    delta), rounded as the tool rounds them."""
    L_ops, O_int = jumps_and_observable(n, jump)
    d = 2**n
    diag = float(np.sum(np.abs(np.diag(rho0)) ** 2))
    off = purity(rho0) - diag
    rho = rho0.copy()
    n_steps = int(t_max / dt)
    rows = []
    for step in range(n_steps + 1):
        if step % record_every == 0:
            t = step * dt
            o = float(np.clip(np.real(np.trace(rho @ O_int)), -1.0, 1.0))
            p = purity(rho)
            rows.append((round(t, 4), round(p, 6), round(mutual_purity(rho, n), 6),
                         round(l1(rho) / (d - 1), 6), round(o, 6),
                         round(round(p, 6) - (diag + off * np.exp(-2 * n * gamma0 * t)), 6)))
        if step == n_steps:
            break
        o = float(np.clip(np.real(np.trace(rho @ O_int)), -1.0, 1.0))
        g = max(0.0, gamma0 * (1 - kappa * o))
        drho = -1j * (H @ rho - rho @ H)
        for L in L_ops:
            LdL = L.conj().T @ L
            drho += g * (L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL))
        rho = rho + dt * drho
        if clip:
            rho = 0.5 * (rho + rho.conj().T)
            w, v = np.linalg.eigh(rho)
            rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
            trace = np.real(np.trace(rho))
            if trace > 0:
                rho /= trace
    return rows


def euler_final(rho0, H, gamma0, kappa, clip, t_max=10.0, dt=0.01):
    """The tool's sigma_z feedback step on two qubits (as in tool_run), returning the final state."""
    L_ops, O = jumps_and_observable(2, "sigma_z")
    rho = rho0.copy()
    for _ in range(int(round(t_max / dt))):
        o = float(np.clip(np.real(np.trace(rho @ O)), -1.0, 1.0))
        g = max(0.0, gamma0 * (1 - kappa * o))
        drho = -1j * (H @ rho - rho @ H)
        for L in L_ops:
            drho += g * (L @ rho @ L - rho)
        rho = rho + dt * drho
        if clip:
            rho = 0.5 * (rho + rho.conj().T)
            w, v = np.linalg.eigh(rho)
            rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
            rho /= np.real(np.trace(rho))
    return rho


def plain_extremes(rho0, n, H, gamma0, kappa, t_max, dt=0.01, jump="sigma_z"):
    """The tool's feedback step without the clipping: the largest purity and the lowest eigenvalue
    over every step up to t_max."""
    L_ops, O = jumps_and_observable(n, jump)
    rho, top, low = rho0.copy(), purity(rho0), 0.0
    for _ in range(int(t_max / dt)):
        o = float(np.clip(np.real(np.trace(rho @ O)), -1.0, 1.0))
        g = max(0.0, gamma0 * (1 - kappa * o))
        drho = -1j * (H @ rho - rho @ H)
        for L in L_ops:
            drho += g * (L @ rho @ L - rho)
        rho = rho + dt * drho
        top = max(top, purity(rho))
        low = min(low, float(np.min(np.linalg.eigvalsh(0.5 * (rho + rho.conj().T)))))
    return top, low


def euler_lowest(H, gamma0, kappa=0.5, t_max=10.0, dt=0.01):
    """The tool's clipped Euler run from Bell+ (sigma_z, two qubits): the lowest concurrence * l1/3
    over every step, the reading of the §4 sweep with the density matrix's own Psi."""
    L_ops, O = jumps_and_observable(2, "sigma_z")
    rho = bell()
    low = concurrence(rho) * l1(rho) / 3
    for _ in range(int(round(t_max / dt))):
        o = float(np.clip(np.real(np.trace(rho @ O)), -1.0, 1.0))
        g = max(0.0, gamma0 * (1 - kappa * o))
        drho = -1j * (H @ rho - rho @ H)
        for L in L_ops:
            drho += g * (L @ rho @ L - rho)
        rho = rho + dt * drho
        rho = 0.5 * (rho + rho.conj().T)
        w, v = np.linalg.eigh(rho)
        rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
        rho /= np.real(np.trace(rho))
        low = min(low, concurrence(rho) * l1(rho) / 3)
    return low


def liouvillian(H, L_ops, g):
    """Column-stacking vec: vec(A rho B) = (B^T kron A) vec(rho)."""
    d = H.shape[0]
    Id = np.eye(d)
    Lv = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for L in L_ops:
        LdL = L.conj().T @ L
        Lv += g * (np.kron(L.conj(), L) - 0.5 * np.kron(Id, LdL) - 0.5 * np.kron(LdL.T, Id))
    return Lv


class Linear:
    """Exact propagation of a run with a constant rate: rho(t) = expm(L t) rho0."""

    def __init__(self, rho0, H, L_ops, g):
        self.d = H.shape[0]
        self.Lv = liouvillian(H, L_ops, g)
        self.v0 = rho0.reshape(-1, order="F")

    def rho(self, t):
        return (expm(self.Lv * t) @ self.v0).reshape(self.d, self.d, order="F")


def feedback_exact(rho0, H, L_ops, O, g0, kappa, t_max, rtol=1e-12):
    d = H.shape[0]

    def f(_, v):
        rho = v.reshape(d, d)
        g = g0 * (1 - kappa * np.real(np.trace(rho @ O)))
        dr = -1j * (H @ rho - rho @ H)
        for L in L_ops:
            dr += g * (L @ rho @ L.conj().T - rho)  # every jump here is Hermitian with L^2 = I
        return dr.reshape(-1)

    return solve_ivp(f, (0, t_max), rho0.reshape(-1).astype(complex), method="DOP853",
                     rtol=rtol, atol=rtol * 1e-2, dense_output=True)


def scalar_exact(rho0, H, g0, t_max, rtol=1e-12):
    """The scalar law with the mutual-purity bridge on two qubits, sigma_z on both sites: the rate
    g0 * C(rho) recomputed from the state, propagated by DOP853."""
    L_ops, _ = jumps_and_observable(2, "sigma_z")

    def f(_, v):
        rho = v.reshape(4, 4)
        g = g0 * mutual_purity(rho, 2)
        dr = -1j * (H @ rho - rho @ H)
        for L in L_ops:
            dr += g * (L @ rho @ L - rho)
        return dr.reshape(-1)

    return solve_ivp(f, (0, t_max), rho0.reshape(-1).astype(complex), method="DOP853",
                     rtol=rtol, atol=rtol * 1e-2, dense_output=True)


def closed_form_x(x0, g0, kappa, t):
    """x/(1 - kappa x) = x0/(1 - kappa x0) * exp(-4 g0 t), solved for x."""
    q = x0 / (1 - kappa * x0) * np.exp(-4 * g0 * t)
    return q / (1 + kappa * q)


VALUES = []
LOGGED = []
CLAIMS = []


def check(label, page, value, digits=3):
    """A number the page prints (a literal) against the value computed here, at the page's precision."""
    ok = round(value, digits) == round(page, digits)
    print(f"  {label:<56s} page {page:.{digits}f}  computed {value:.{digits + 3}f}  {'match' if ok else 'DIFFERS'}")
    VALUES.append((label, ok))


def claim(label, ok, detail):
    print(f"  {label:<56s} {detail}  {'holds' if ok else 'FAILS'}")
    CLAIMS.append((label, ok))


QUOTES = []


def quote(label, page, record, digits):
    """A number the page quotes from the log, against the log's record (both literals)."""
    ok = round(record, digits) == round(page, digits)
    print(f"  {label:<56s} page {page:.{digits}f}  log {record:.{digits}f}  {'match' if ok else 'DIFFERS'}")
    QUOTES.append((label, ok))


def converges(label, closed, run):
    """A closed form against propagation by DOP853 at rtol 1e-5, 1e-7, 1e-9, 1e-11, at the times run
    returns: at every tolerance the largest deviation must stay within rtol. This is a one-sided
    threshold, not an error law (rtol bounds the solver's local error, not the global one); a closed form
    off by more than about 1e-11 at any of those times fails it."""
    rtols = (1e-5, 1e-7, 1e-9, 1e-11)
    devs = [float(np.max(np.abs(np.asarray(run(rtol)) - np.asarray(closed)))) for rtol in rtols]
    ok = all(d <= r for d, r in zip(devs, rtols))
    claim(label, ok, "deviation at rtol 1e-5..1e-11: " + " ".join(f"{d:.1e}" for d in devs)
          + "; ratio to rtol " + " ".join(f"{d / r:.3f}" for d, r in zip(devs, rtols)))


def exact(label, M, zero):
    """An identity compared exactly on matrices whose entries are integers or exactly representable,
    so the products are exact; zero=False is a control through the same door, which must be nonzero."""
    largest = float(np.max(np.abs(M)))
    ok = (largest == 0.0) if zero else (largest > 0.0)
    claim(label, ok, f"largest entry {largest:g} ({'== 0 required' if zero else 'nonzero required'})")


def logged(label, record, computed):
    hits = [a == b for a, b in zip(record, computed)]
    print(f"    {label:<56s} {sum(hits)} of {len(hits)} logged <O_int> match, last {record[-1]:.6f}")
    LOGGED.extend((f"{label} #{k}", ok) for k, ok in enumerate(hits))


def crossings(fun, t_max, grid=0.005):
    """Every time in (0, t_max] where fun(t) - 1/4 changes sign, refined by brentq."""
    ts = np.arange(0.0, t_max + grid / 2, grid)
    vals = np.array([fun(t) - 0.25 for t in ts])
    out = []
    for a, b, va, vb in zip(ts, ts[1:], vals, vals[1:]):
        if va == 0.0 or va * vb < 0:
            out.append(brentq(lambda s: fun(s) - 0.25, a, b, xtol=1e-12))
    return out, ts, vals + 0.25


def refined_extremum(fun, t_max, sign=1.0, grid=0.005):
    """The largest (sign 1) or smallest (sign -1) value of fun on [0, t_max]: a grid, then a bounded
    search on the bracket around the best grid point (the readings have cusps, which a grid misses)."""
    ts = np.arange(0.0, t_max + grid / 2, grid)
    vals = [sign * fun(t) for t in ts]
    k = int(np.argmax(vals))
    a, b = ts[max(k - 1, 0)], ts[min(k + 1, len(ts) - 1)]
    res = minimize_scalar(lambda s: -sign * fun(s), bounds=(a, b), method="bounded", options={"xatol": 1e-10})
    best = max((vals[k], ts[k]), (-res.fun, res.x))
    return sign * best[0], best[1]


# The <O_int> records the MCP log shows of each call of 2026-02-20 (the Bell+ calls: t = 3.9 to 5.0, the
# xx tail t = 3.3 to 5.0, the two zz tails losing their first record to the truncation; the GHZ_4 call:
# t = 28.9 to 30.0).
LOG_TAILS = {
    ("sigma_z", 0.99): [0.964748, 0.963022, 0.961235, 0.959383, 0.957466, 0.955481, 0.953425, 0.951298,
                        0.949096, 0.946817, 0.94446, 0.942021],
    ("sigma_z", 0.5): [0.347347, 0.336001, 0.324952, 0.314198, 0.303735, 0.29356, 0.283668, 0.274055,
                       0.264719, 0.255652, 0.246852, 0.238314],
    ("xx", 0.95): [1.0] * 18,
    ("sigma_z", 0.95): [0.845355, 0.838746, 0.831979, 0.825055, 0.817974, 0.810735, 0.80334, 0.795791,
                        0.788087, 0.780232, 0.772227, 0.764071],
    ("yy", 0.95): [-0.101073, 0.040669, 0.178618, 0.307717, 0.423264, 0.521049, 0.597499, 0.649808,
                   0.676053, 0.675299, 0.647664, 0.594354],
    ("zz", 0.95): [-0.091401, -0.233652, -0.363834, -0.476744, -0.568093, -0.634668, -0.674424,
                   -0.686512, -0.671239, -0.629981, -0.565059],
    ("sigma_z", 0.0): [0.209579, 0.201345, 0.193435, 0.185835, 0.178534, 0.17152, 0.164781, 0.158307,
                       0.152088, 0.146113, 0.140372, 0.134857],
    ("yy", 0.0): [-0.078082, 0.067301, 0.207711, 0.337669, 0.452197, 0.547007, 0.618658, 0.664684,
                  0.683669, 0.675291, 0.640316, 0.580553],
    ("zz", 0.0): [-0.067301, -0.207711, -0.337669, -0.452197, -0.547007, -0.618658, -0.664684,
                  -0.683669, -0.675291, -0.640316, -0.580553],
}
LOG_GHZ4 = [0.048332, 0.048507, 0.048681, 0.048856, 0.04903, 0.049205, 0.049379, 0.049554, 0.049729,
            0.049903, 0.050078, 0.050253]


def section_identities():
    print("Identities, compared exactly")
    n = 2
    XX, YY, ZZ = pair(X, 0, 1, n), pair(Y, 0, 1, n), pair(Z, 0, 1, n)
    Hb, Fx, Fz = heisenberg(n), field(n, X), field(n, Z)
    for name, M in (("bonds", Hb), ("x field", Fx)):
        exact(f"[{name}, X0X1] = 0", M @ XX - XX @ M, True)
    exact("control: [z field, X0X1] = 0 fails", Fz @ XX - XX @ Fz, False)
    exact("z field: [z field, Z0Z1] = 0", Fz @ ZZ - ZZ @ Fz, True)
    P = (np.eye(4) + XX) / 2
    exact("on the X0X1 = +1 sector, YY = -ZZ", P @ YY @ P + P @ ZZ @ P, True)
    exact("control: on that sector YY = +ZZ fails", P @ YY @ P - P @ ZZ @ P, False)
    phi = np.array([1, 0, 0, 1], dtype=complex)
    exact("X0X1 Bell+ = Bell+ (the xx jump fixes it)", XX @ phi - phi, True)
    exact("YY Bell+ = -Bell+, ZZ Bell+ = +Bell+", np.concatenate([YY @ phi + phi, ZZ @ phi - phi]), True)
    exact("control: [x field, Y0Y1] = 0 fails", Fx @ YY - YY @ Fx, False)
    # The mutual-purity bridge stays 1/2 on Bell+: X0X1 fixes H and Bell+ and maps each jump to +-itself
    # and every Y_k, Z_k to its negative, so <Y_k> = <Z_k> = 0; X0 + X1 commutes with H, and the adjoint
    # dissipator of each channel maps it to a multiple of itself, so <X0 + X1> stays 0; SWAP gives
    # <X0> = <X1>. These hold for any rate, so also under either feedback law.
    S = field(n, X)
    exact("[bonds, X0 + X1] = 0 (the x field is X0 + X1 itself)", Hb @ S - S @ Hb, True)
    for jump in ("sigma_z", "xx", "yy", "zz"):
        L_ops, _ = jumps_and_observable(n, jump)
        adj = sum(L @ S @ L - S for L in L_ops)
        c = {"sigma_z": -2, "xx": 0, "yy": -2, "zz": -2}[jump]
        exact(f"{jump}: D'(X0 + X1) = {c}(X0 + X1) at unit rate", adj - c * S, True)
        for L in L_ops:
            plus, minus = XX @ L @ XX - L, XX @ L @ XX + L
            exact(f"{jump}: X0X1 L X0X1 = +-L", plus if np.max(np.abs(plus)) <= np.max(np.abs(minus)) else minus, True)
    for k in range(2):
        exact(f"Z{k} X0X1 Z{k} = -X0X1 (the sigma_z channel)", site(Z, k, n) @ XX @ site(Z, k, n) + XX, True)
    print("S = sum_{i<j} X_i X_j, N = 3..6")
    for n in range(3, 7):
        Sp = sum(pair(X, i, j, n) for i in range(n) for j in range(i + 1, n))
        Snn = sum(pair(X, i, i + 1, n) for i in range(n - 1))
        Hb, Fx = heisenberg(n), field(n, X)
        exact(f"N={n}: [bonds, S] = 0", Hb @ Sp - Sp @ Hb, True)
        exact(f"N={n}: [ring bonds, S] = 0", heisenberg(n, ring=True) @ Sp - Sp @ heisenberg(n, ring=True), True)
        exact(f"N={n}: [x field, S] = 0", Fx @ Sp - Sp @ Fx, True)
        exact(f"N={n}: [X_iX_j, S] = 0 for every pair (nearest-neighbour ones included)",
              np.array([np.max(np.abs(pair(X, i, j, n) @ Sp - Sp @ pair(X, i, j, n)))
                        for i in range(n) for j in range(i + 1, n)]), True)
        exact(f"N={n}: sigma_z on every site: sum_k (Z_k S Z_k - S) = -4 S",
              sum(site(Z, k, n) @ Sp @ site(Z, k, n) - Sp for k in range(n)) + 4 * Sp, True)
        exact(f"N={n}: control: [bonds, nearest-neighbour sum] = 0 fails", Hb @ Snn - Snn @ Hb, False)
        Hxy = heisenberg(n, xy=True)
        exact(f"N={n}: control: [XY bonds, S] = 0 fails", Hxy @ Sp - Sp @ Hxy, False)
        g = np.zeros(2**n)
        g[0] = g[-1] = 1
        exact(f"N={n}: <GHZ|S|GHZ> = 0", np.array([g @ Sp.real @ g]), True)
    g2 = np.array([1, 0, 0, 1.0])
    exact("control: N=2: <Bell+|X0X1|Bell+> = 0 fails", np.array([g2 @ pair(X, 0, 1, 2).real @ g2]), False)


def section_2():
    print("\n§2, sigma_z feedback: x = <X0X1> obeys x/(1 - kappa x) = x0/(1 - kappa x0) exp(-4 gamma_0 t)")
    L_ops, O = jumps_and_observable(2, "sigma_z")
    for h in (0.0, 0.5, 0.9):
        sol = scalar_exact(bell(), hamiltonian(2, h=h), 0.005, 10.0)
        cdev = max(abs(mutual_purity(sol.sol(t).reshape(4, 4), 2) - 0.5) for t in np.linspace(0, 10, 101))
        claim(f"the scalar law on Bell+, h {h}: C stays 1/2 (|C - 1/2| <= 100 eps, a bound on rounding)",
              cdev <= 100 * EPS, f"largest |C - 1/2| {cdev:.1e}")
        if h == 0.0:
            check("the scalar law on Bell+, gamma_base 0.005, h = 0: Psi(10)", 0.302, l1(sol.sol(10.0).reshape(4, 4)) / 3)
    for name, bond in (("XY", heisenberg(2, xy=True)), ("Ising", pair(PAULI["z"], 0, 1, 2))):
        sol = scalar_exact(bell(), bond + 0.5 * field(2), 0.1, 10.0)
        cdev = max(abs(mutual_purity(sol.sol(t).reshape(4, 4), 2) - 0.5) for t in np.linspace(0, 10, 101))
        claim(f"control, the scalar law on Bell+ under the {name} bond, h 0.5, gamma_base 0.1: C leaves 1/2 "
              "(by more than 0.01)", cdev > 0.01, f"largest |C - 1/2| {cdev:.3f}")
    sol = scalar_exact(bell(), hamiltonian(2, h=0.0), 0.1, 5.0)
    cdev = max(abs(mutual_purity(sol.sol(t).reshape(4, 4), 2) - 0.5) for t in np.linspace(0, 5, 51))
    claim("the scalar law on Bell+, h 0, gamma_base 0.1: C stays 1/2, so the rate is gamma_base/2 (<= 100 eps)",
          cdev <= 100 * EPS, f"largest |C - 1/2| {cdev:.1e}")
    half = Linear(bell(), hamiltonian(2, h=0.0), jumps_and_observable(2, "sigma_z")[0], 0.05)
    ddev = max(abs(purity(half.rho(t)) - (0.5 + 0.5 * np.exp(-2 * 2 * 0.1 * t))) for t in (0.5, 2.0, 5.0))
    claim("  at that rate the model's purity is the delta curve: delta = 0 at t = 0.5, 2, 5 (expm, <= 100 eps)",
          ddev <= 100 * EPS, f"largest |delta| {ddev:.1e}")
    x = closed_form_x(1.0, 0.005, 0.5, 10.0)
    check("gamma_0 0.005, kappa 0.5: x at t = 10", 0.900, x)
    check("  its rate at t = 10", 0.00275, 0.005 * (1 - 0.5 * x), digits=5)
    check("  late-time shift ln(1/(1 - kappa))/(4 gamma_0)", 34.7, np.log(2) / 0.02, digits=1)
    check("  lag behind the plain curve at t = 10", 4.75, 10.0 + np.log(x) / 0.02, digits=2)
    x = closed_form_x(1.0, 0.1, 0.99, 5.0)
    check("gamma_0 0.1, kappa 0.99: x at t = 5", 0.940, x)
    check("  late-time shift", 11.5, np.log(100) / 0.4, digits=1)
    check("  lag at t = 5", 4.85, 5.0 + np.log(x) / 0.4, digits=2)
    check("kappa 0: exp(-2)", 0.135, closed_form_x(1.0, 0.1, 0.0, 5.0))
    check("late-time shift ln(1/(1 - kappa x0))/(4 gamma_0), W3: x0 = 2/3, gamma_0 0.005, kappa 0.5",
          20.3, np.log(1 / (1 - 0.5 * 2 / 3)) / 0.02, digits=1)
    quote("the page's quote of the log, kappa 0.99 at t = 5", 0.942, LOG_TAILS[("sigma_z", 0.99)][-1], 3)
    quote("the page's quote of the log, kappa 0 at t = 5", 0.135, LOG_TAILS[("sigma_z", 0.0)][-1], 3)
    times = (1.0, 2.0, 5.0)
    converges("the closed form against DOP853: Bell+, kappa 0.99, x at t = 1, 2, 5",
              [closed_form_x(1.0, 0.1, 0.99, t) for t in times],
              lambda r: [np.real(np.trace(feedback_exact(bell(), hamiltonian(2), L_ops, O, 0.1, 0.99, 5.0, rtol=r)
                                          .sol(t).reshape(4, 4) @ O)) for t in times])
    L3, O3 = jumps_and_observable(3, "sigma_z")
    converges("the closed form against DOP853: W3, ring, h 0.9, x at t = 1, 2, 5",
              [closed_form_x(2 / 3, 0.005, 0.5, t) for t in times],
              lambda r: [np.real(np.trace(feedback_exact(w_state(3), hamiltonian(3, h=0.9, ring=True), L3, O3, 0.005,
                                                         0.5, 5.0, rtol=r).sol(t).reshape(8, 8) @ O3)) for t in times])
    print("  the closed form against DOP853 on Bell+ (J, h vary; the form does not see them), at five times:")
    for g0, kappa, J, h, T in ((0.005, 0.5, 1.0, 0.5, 10.0), (0.005, 0.5, 1.0, 0.9, 10.0),
                               (0.1, 0.99, 1.0, 0.5, 5.0), (0.1, 0.5, 2.0, 1.3, 5.0)):
        sol = feedback_exact(bell(), hamiltonian(2, J=J, h=h), L_ops, O, g0, kappa, T)
        dev = max(abs(np.real(np.trace(sol.sol(t).reshape(4, 4) @ O)) - closed_form_x(1.0, g0, kappa, t))
                  for t in np.linspace(T / 5, T, 5))
        claim(f"gamma_0 {g0}, kappa {kappa}, J {J}, h {h}, t <= {T:g}: |DOP853 - closed form| <= 100 eps (rounding)",
              dev <= 100 * EPS, f"{dev:.1e}")
    sol = feedback_exact(bell(), hamiltonian(2, axis="z"), L_ops, O, 0.1, 0.5, 5.0)
    xz = np.real(np.trace(sol.sol(5.0).reshape(4, 4) @ O))
    check("control, field along z (J 1, h 0.5, gamma_0 0.1, kappa 0.5)", -0.115, xz)
    check("  against the closed form", 0.238, closed_form_x(1.0, 0.1, 0.5, 5.0))
    print("  states invariant under every permutation of the sites, sigma_z on every site, kappa 0.5:")
    for label, rho0, n, H, g0, T, page_x in (
            ("W3, ring, h 0.9, gamma_0 0.005", w_state(3), 3, hamiltonian(3, h=0.9, ring=True), 0.005, 5.0, 0.623),
            ("W4, chain, h 0.5, gamma_0 0.1", w_state(4), 4, hamiltonian(4), 0.1, 5.0, None),
            ("W4, chain, h 1.3, gamma_0 0.1", w_state(4), 4, hamiltonian(4, h=1.3), 0.1, 5.0, None),
            ("W5, chain, h 0.5, gamma_0 0.1", w_state(5), 5, hamiltonian(5), 0.1, 5.0, None),
            ("GHZ4, chain, h 0.5, gamma_0 0.1", ghz(4), 4, hamiltonian(4), 0.1, 5.0, None)):
        L_n, O_n = jumps_and_observable(n, "sigma_z")
        sol = feedback_exact(rho0, H, L_n, O_n, g0, 0.5, T)
        x0 = np.real(np.trace(rho0 @ O_n))
        xT = np.real(np.trace(sol.sol(T).reshape(2**n, 2**n) @ O_n))
        dev = max(abs(np.real(np.trace(sol.sol(t).reshape(2**n, 2**n) @ O_n)) - closed_form_x(x0, g0, 0.5, t))
                  for t in (1.0, 2.0, 5.0))
        claim(f"{label}: x from {x0:.4f} on the closed form at t = 1, 2, 5 (|deviation| <= 100 eps, rounding)",
              dev <= 100 * EPS, f"x(5) {xT:.6f}, deviation {dev:.1e}")
        if page_x is not None:
            check(f"{label}: x at t = 5", page_x, xT)
    n = 3
    L_n, O_n = jumps_and_observable(n, "sigma_z")
    rho0 = np.kron(bell(), pure([1, 0]))
    sol = feedback_exact(rho0, hamiltonian(3), L_n, O_n, 0.1, 0.5, 5.0)
    check("control, Bell+ x |0> on the N = 3 chain: x at t = 5", 0.039,
          np.real(np.trace(sol.sol(5.0).reshape(8, 8) @ O_n)))
    sol = feedback_exact(w_state(3), hamiltonian(3, xy=True), L_n, O_n, 0.1, 0.5, 5.0)
    x_xy = np.real(np.trace(sol.sol(5.0).reshape(8, 8) @ O_n))
    x_cf = closed_form_x(2 / 3, 0.1, 0.5, 5.0)
    claim("control, W3 on the XY chain (h 0.5, gamma_0 0.1, kappa 0.5): the closed form fails (by more than 0.01)",
          abs(x_xy - x_cf) > 0.01, f"x(5) {x_xy:.4f} against the form's {x_cf:.4f}")
    L4z, O4z = jumps_and_observable(4, "sigma_z")
    sol = feedback_exact(w_state(4), hamiltonian(4, h=0.0, xy=True), L4z, O4z, 0.1, 0.5, 5.0)
    x_xy4 = np.real(np.trace(sol.sol(5.0).reshape(16, 16) @ O4z))
    x_cf4 = closed_form_x(0.5, 0.1, 0.5, 5.0)
    claim("control, W4 on the XY chain without the field: the closed form fails (by more than 0.01)",
          abs(x_xy4 - x_cf4) > 0.01, f"x(5) {x_xy4:.4f} against the form's {x_cf4:.4f}")
    sol = feedback_exact(ghz(3), hamiltonian(3, xy=True), L_n, O_n, 0.1, 0.5, 5.0)
    top = max(abs(np.real(np.trace(sol.sol(t).reshape(8, 8) @ O_n))) for t in np.linspace(0, 5, 251))
    claim("control, GHZ3 on the XY chain under sigma_z jumps: <X0X1> leaves zero (by more than 0.01)", top > 0.01,
          f"largest {top:.3f}")
    rows = tool_run(bell(), 2, hamiltonian(2, h=0.0), 0.005, 0.5, "sigma_z", t_max=10.0)
    check("Dynamic Fixed Points §6, the operator law at h = 0: the tool's C*Psi at t = 10", 0.150,
          rows[-1][2] * rows[-1][3])
    for kappa, page_d in ((0.0, -0.074), (0.5, -0.004)):
        rows = tool_run(bell(), 2, hamiltonian(2, h=0.0), 0.1, kappa, "sigma_z")
        check(f"the tool's delta at t = 0.5, h 0, gamma_0 0.1, kappa {kappa}", page_d,
              next(r[5] for r in rows if abs(r[0] - 0.5) < 1e-9))


def section_82():
    print("\n§8.2, Bell+, heisenberg, J 1, h 0.5, gamma_0 0.1, t_max 5: the tool's Euler run")
    n, H = 2, hamiltonian(2)
    finals = {}
    for (jump, kappa), record in LOG_TAILS.items():
        rows = tool_run(bell(), n, H, 0.1, kappa, jump)
        finals[(jump, kappa)] = rows[-1]
        logged(f"{jump} kappa {kappa}", record, [r[4] for r in rows[-len(record):]])
    page = {"sigma_z": (0.124, 0.336), "xx": (0.348, 1.000), "yy": (0.287, 0.734), "zz": (0.287, 0.734)}
    for jump, (p_cpsi, p_pur) in page.items():
        rows = tool_run(bell(), n, H, 0.1, 0.0, jump)
        _, p, c, psi, _, delta = rows[-1]
        check(f"tool, {jump}, kappa 0: C*Psi", p_cpsi, c * psi)
        check(f"tool, {jump}, kappa 0: purity", p_pur, p)
        print(f"    delta at t = 5 is purity - {p - delta:.6f}")
    _, p, c, psi, _, _ = finals[("sigma_z", 0.95)]
    check("tool, sigma_z, kappa 0.95: C*Psi", 0.287, c * psi)
    check("tool, sigma_z, kappa 0.95: purity", 0.785, p)
    quote("the page's quote of the log, yy tail at kappa 0.95", 0.594354, LOG_TAILS[("yy", 0.95)][-1], 6)
    quote("the page's quote of the log, zz tail at kappa 0.95", -0.565059, LOG_TAILS[("zz", 0.95)][-1], 6)
    quote("the page's quote of the log, yy tail at kappa 0", 0.580553, LOG_TAILS[("yy", 0.0)][-1], 6)
    quote("the page's quote of the log, zz tail at kappa 0", -0.580553, LOG_TAILS[("zz", 0.0)][-1], 6)

    print("§8.2 exact (expm), kappa 0, read with the tool's C*Psi; C stays 1/2 (identities above)")
    expect = {"sigma_z": (0.112, 0.323, 0.413, 0.695, 1.66, 2.54),
              "xx": (0.348, 1.000, 0.500, None, 4.12, None),
              "yy": (0.275, 0.693, 0.476, 0.761, 3.81, 13.50),
              "zz": (0.275, 0.693, 0.476, 0.761, 3.81, 13.50)}
    states = {}
    for jump, (p_cpsi, p_pur, p_max, p_tmax, p_above, p_last) in expect.items():
        L_ops, _ = jumps_and_observable(n, jump)
        run = Linear(bell(), H, L_ops, 0.1)
        rho5 = run.rho(5.0)
        states[jump] = rho5
        fun = lambda t, run=run: cpsi(run.rho(t), n)
        check(f"exact, {jump}: C*Psi at t = 5", p_cpsi, cpsi(rho5, n))
        check(f"exact, {jump}: purity at t = 5", p_pur, purity(rho5))
        top, t_top = refined_extremum(fun, 5.0)
        check(f"exact, {jump}: largest C*Psi on [0, 5]", p_max, top)
        if p_tmax is not None:
            check(f"exact, {jump}: time of the largest", p_tmax, t_top)
        cross, ts, _ = crossings(fun, 60.0)
        inside = [c for c in cross if c <= 5.0]
        edges = [0.0] + inside + [5.0]
        above = sum(b - a for a, b in zip(edges, edges[1:]) if fun(0.5 * (a + b)) > 0.25)
        check(f"exact, {jump}: time above 1/4 in [0, 5]", p_above, above, digits=2)
        if p_last is not None:
            check(f"exact, {jump}: last crossing of 1/4 (t <= 60)", p_last, cross[-1], digits=2)
        else:
            print(f"    {jump}: {len(cross)} crossings up to t = 60, the last at {cross[-1]:.3f}")
        if jump in ("yy", "zz"):
            check(f"exact, {jump}: C*Psi at t = 60 (the limit 1/6)", 1 / 6, fun(60.0))
        cdev = max(abs(mutual_purity(run.rho(t), n) - 0.5) for t in np.arange(0, 20.0001, 0.5))
        print(f"    {jump}: largest |C - 1/2| over t <= 20 {cdev / EPS:.1f} eps")
    print(f"  yy against zz at t = 5: largest entry difference {np.max(np.abs(states['yy'] - states['zz'])) / EPS:.1f} eps")
    run = Linear(bell(), H, [pair(X, 0, 1, 2)], 0.1)
    ts = np.arange(0, 10.0001, 0.05)
    worst = max(abs(cpsi(run.rho(t), 2) - (1 + 2 * abs(np.sin(2.0 * t))) / 6) for t in ts)
    print(f"  xx: C*Psi = (1 + 2|sin 4ht|)/6, largest deviation over t <= 10: {worst / EPS:.1f} eps")
    pur_book = [purity(run.rho(t)) * l1(run.rho(t)) / 3 for t in ts]
    conc_book = [concurrence(run.rho(t)) * l1(run.rho(t)) / 3 for t in ts]
    check("xx, purity x l1/3: smallest over t <= 10", 1 / 3, min(pur_book))
    check("xx, purity x l1/3: largest over t <= 10", 1.0, max(pur_book))
    check("xx, concurrence x l1/3: smallest over t <= 10", 1 / 3, min(conc_book))
    check("xx, concurrence x l1/3: largest over t <= 10", 1.0, max(conc_book))
    check("Bell+, where every channel starts: purity x l1/3", 1 / 3, purity(bell()) * l1(bell()) / 3)
    check("Bell+, where every channel starts: concurrence x l1/3", 1 / 3, concurrence(bell()) * l1(bell()) / 3)
    Lx, Ox = jumps_and_observable(n, "xx")
    converges("(1 + 2|sin 4ht|)/6 against DOP853: xx, t = 3.3", (1 + 2 * abs(np.sin(2.0 * 3.3))) / 6,
              lambda r: cpsi(feedback_exact(bell(), H, Lx, Ox, 0.1, 0.0, 3.3, rtol=r).sol(3.3).reshape(4, 4), 2))
    print("  the same runs with the field along z:")
    Hz = hamiltonian(2, axis="z")
    for jump, page_p in (("xx", 0.693), ("yy", 0.693), ("zz", 1.000)):
        L_ops, _ = jumps_and_observable(n, jump)
        check(f"field along z, {jump}: purity at t = 5", page_p, purity(Linear(bell(), Hz, L_ops, 0.1).rho(5.0)))
    for jump, page_p in (("yy", 0.6918), ("zz", 0.6929)):
        L_ops, O = jumps_and_observable(n, jump)
        sol = feedback_exact(bell(), H, L_ops, O, 0.1, 0.5, 5.0)
        check(f"exact, {jump}, kappa 0.5: purity at t = 5", page_p, purity(sol.sol(5.0).reshape(4, 4)), digits=4)
    L_ops, O = jumps_and_observable(n, "sigma_z")
    sol = feedback_exact(bell(), H, L_ops, O, 0.1, 0.95, 30.0)
    r5 = sol.sol(5.0).reshape(4, 4)
    check("exact, sigma_z, kappa 0.95: C*Psi at t = 5", 0.284, cpsi(r5, 2))
    check("exact, sigma_z, kappa 0.95: purity at t = 5", 0.773, purity(r5))
    cross, _, _ = crossings(lambda t: cpsi(sol.sol(t).reshape(4, 4), 2), 30.0)
    check("exact, sigma_z, kappa 0.95: crossings of 1/4 (t <= 30)", 12, len(cross), digits=0)
    dip = [c for c in cross if 4.0 < c < 5.0]
    check("  it dips below 1/4 at", 4.52, dip[0], digits=2)
    check("  and is above again from", 4.92, dip[1], digits=2)
    check("  the last crossing", 8.85, cross[-1], digits=2)
    print("  the tool against exact at t = 5, kappa 0, gap in per cent:")
    gaps = {}
    for jump in ("sigma_z", "xx", "yy", "zz"):
        L_ops, _ = jumps_and_observable(n, jump)
        ex = Linear(bell(), H, L_ops, 0.1).rho(5.0)
        for clip in (True, False):
            r = tool_run(bell(), n, H, 0.1, 0.0, jump, clip=clip)[-1]
            g_c, g_p = 100 * abs(r[2] * r[3] / cpsi(ex, 2) - 1), 100 * abs(r[1] / purity(ex) - 1)
            if clip:
                gaps[(jump, "C*Psi")], gaps[(jump, "purity")] = g_c, g_p
            print(f"    {jump}, {'with' if clip else 'without'} clipping: C*Psi {g_c:.1f}, purity {g_p:.1f}")
    worst = max(gaps, key=gaps.get)
    check(f"the largest gap in the table, per cent ({worst[0]}, {worst[1]})", 11, gaps[worst], digits=0)
    claim("the largest gap is the C*Psi of the sigma_z run", worst == ("sigma_z", "C*Psi"), f"{worst}")
    L_ops, _ = jumps_and_observable(n, "sigma_z")
    ex = cpsi(Linear(bell(), H, L_ops, 0.1).rho(5.0), 2)
    for dt, clip, page_gap in ((0.01, False, 11.6), (0.01, True, 10.6), (0.005, True, 5.3), (0.0025, True, 2.7)):
        r = tool_run(bell(), n, H, 0.1, 0.0, "sigma_z", dt=dt, clip=clip)[-1]
        check(f"sigma_z C*Psi gap, dt {dt}, {'with' if clip else 'without'} clipping, per cent",
              page_gap, 100 * (r[2] * r[3] / ex - 1), digits=1)


def section_81():
    print("\n§8.1, GHZ_N, x_pairs, heisenberg chain, J 1, h 0.5, gamma_0 0.1, t = 5")
    page_exact = {3: (0.2637, 4), 4: (0.2503, 4), 5: (0.06261, 5), 6: (0.06250, 5)}
    page_max = {3: (0.339, 3), 4: (0.25004, 5), 5: (0.164, 3), 6: (0.115, 3)}
    for n in range(3, 7):
        H = hamiltonian(n)
        closed = 2.0 ** -(n - 1) * sum(comb(n, a) * np.exp(-4 * 0.1 * a * (n - a) * 5.0)
                                       for a in range(0, n + 1, 2))
        check(f"GHZ{n}: purity at t = 5, closed form", page_exact[n][0], closed, digits=page_exact[n][1])
        L_ops, O = jumps_and_observable(n, "x_pairs")
        sol = feedback_exact(ghz(n), H, L_ops, O, 0.1, 0.99, 10.0, rtol=1e-11)
        ode = purity(sol.sol(5.0).reshape(2**n, 2**n))
        worst_o = max(abs(np.real(np.trace(sol.sol(t).reshape(2**n, 2**n) @ O))) for t in np.linspace(0, 10, 101))
        print(f"    DOP853 at kappa 0.99: {ode:.6f} (difference {abs(ode - closed):.1e}), largest |<O_int>| {worst_o:.1e}")
        top, t_top = refined_extremum(lambda t: cpsi(sol.sol(t).reshape(2**n, 2**n), n), 5.0)
        check(f"GHZ{n}: largest C*Psi (tool's reading, exact) on [0, 5]", page_max[n][0], top, digits=page_max[n][1])
        print(f"    at t = {t_top:.3f}, value {top:.6f}")
    for n in range(3, 7):
        Ln, On = jumps_and_observable(n, "x_pairs")
        Pn = 2.0 ** -(n - 1) * sum(comb(n, a) * np.exp(-4 * 0.1 * a * (n - a) * 5.0) for a in range(0, n + 1, 2))
        converges(f"P(t) against DOP853: GHZ{n}, x_pairs, kappa 0.99, t = 5", Pn,
                  lambda r, n=n, Ln=Ln, On=On: purity(feedback_exact(ghz(n), hamiltonian(n), Ln, On, 0.1, 0.99, 5.0,
                                                                     rtol=r).sol(5.0).reshape(2**n, 2**n)))
    for h, page_c in ((0.5, 0.127), (0.45, 0.103), (0.55, 0.144)):
        r = tool_run(ghz(4), 4, hamiltonian(4, h=h), 0.1, 0.0, "x_pairs")[-1]
        check(f"the tool's GHZ4 C*Psi at t = 5, kappa 0, h {h} (the page: {page_c})", page_c, r[2] * r[3])
    L4, O4 = jumps_and_observable(4, "x_pairs")
    sol = feedback_exact(ghz(4), hamiltonian(4, h=0.45), L4, O4, 0.1, 0.99, 5.0, rtol=1e-11)
    check("GHZ4 at h = 0.45: largest C*Psi (tool's reading, exact) on [0, 5]", 0.238,
          refined_extremum(lambda t: cpsi(sol.sol(t).reshape(16, 16), 4), 5.0)[0])
    february = {4: 0.252, 5: 0.063, 6: 0.063}
    for n, feb in february.items():
        at05 = tool_run(ghz(n), n, hamiltonian(n, h=0.5), 0.1, 0.5, "x_pairs")[-1][1]
        at02 = tool_run(ghz(n), n, hamiltonian(n, h=0.2), 0.1, 0.5, "x_pairs")[-1][1]
        check(f"the February GHZ{n} purity, regenerated at h = 0.5", feb, at05)
        at045 = tool_run(ghz(n), n, hamiltonian(n, h=0.45), 0.1, 0.5, "x_pairs")[-1][1]
        check(f"the February GHZ{n} purity, regenerated at h = 0.45 too", feb, at045)
        if n == 4:
            check("the tool's GHZ4 purity at h = 0.2 (the page's 0.250)", 0.250, at02)
            claim("  so h = 0.2 does not regenerate the February GHZ4 row", round(at02, 3) != feb, f"{at02:.4f}")
        else:
            claim(f"the February GHZ{n} purity regenerates at h = 0.2 as well (no page states it)",
                  round(at02, 3) == feb, f"{at02:.4f}")
    for n, lo, hi, digits in ((4, 0.2523, 0.2524, 4), (5, 0.06263, 0.06263, 5), (6, 0.06333, 0.06333, 5)):
        ps = [tool_run(ghz(n), n, hamiltonian(n), 0.1, k, "x_pairs")[-1][1] for k in (0.0, 0.5, 0.99, 1.0)]
        check(f"tool, GHZ{n}: smallest purity over kappa 0, 0.5, 0.99, 1", lo, min(ps), digits=digits)
        check(f"tool, GHZ{n}: largest purity over kappa 0, 0.5, 0.99, 1", hi, max(ps), digits=digits)
        if n == 4:
            claim("tool, GHZ4: spread over kappa below 6e-5", max(ps) - min(ps) < 6e-5, f"spread {max(ps) - min(ps):.2e}")
    for clip, page_o, digits, page_p in ((True, 0.05, 2, 0.2523), (False, 0.0, 6, 0.2885)):
        rows = tool_run(ghz(4), 4, hamiltonian(4), 0.1, 0.5, "x_pairs", clip=clip)
        check(f"tool, GHZ4, kappa 0.5, Euler {'with' if clip else 'without'} clipping: largest |<O_int>|",
              page_o, max(abs(r[4]) for r in rows), digits=digits)
        check(f"tool, GHZ4, kappa 0.5, Euler {'with' if clip else 'without'} clipping: purity at t = 5",
              page_p, rows[-1][1], digits=4)
    rows = tool_run(ghz(4), 4, hamiltonian(4, h=0.2), 0.1, 0.99, "x_pairs", t_max=30.0)
    logged("GHZ4 x_pairs, h 0.2, kappa 0.99, t_max 30", LOG_GHZ4, [r[4] for r in rows[-len(LOG_GHZ4):]])
    print("  scope, N = 4, kappa 0 against 0.99 (DOP853):")
    n = 4
    S = sum(pair(X, i, j, n) for i in range(n) for j in range(i + 1, n))
    for label, pairs, xy, p_o, p_p0, p_p1 in (
            ("nearest-neighbour pairs, Heisenberg", [(i, i + 1) for i in range(n - 1)], False, 0.05, 0.2755, 0.2784),
            ("all pairs, XY chain", [(i, j) for i in range(n) for j in range(i + 1, n)], True, 0.30, 0.1723, 0.1738)):
        L_ops = [pair(X, i, j, n) for i, j in pairs]
        O = sum(L_ops) / len(L_ops)
        H = hamiltonian(n, xy=xy)
        out = []
        for kappa in (0.0, 0.99):
            sol = feedback_exact(ghz(n), H, L_ops, O, 0.1, kappa, 5.0, rtol=1e-11)
            states = [sol.sol(t).reshape(2**n, 2**n) for t in np.linspace(0, 5, 101)]
            out.append((purity(states[-1]), max(abs(np.real(np.trace(r @ O))) for r in states),
                        max(abs(np.real(np.trace(r @ S))) for r in states)))
        check(f"{label}: largest |<O_int>| at kappa 0", p_o, out[0][1], digits=2)
        check(f"{label}: purity at kappa 0", p_p0, out[0][0], digits=4)
        check(f"{label}: purity at kappa 0.99", p_p1, out[1][0], digits=4)
        print(f"    largest |<S>|: {out[0][2]:.1e} (kappa 0), {out[1][2]:.1e} (kappa 0.99)")


def section_4():
    print("\n§4, the agents' sweep (Bell+, sigma_z, kappa 0.5, t_max 10, concurrence) and Simulation Evidence §2")
    L_ops, O = jumps_and_observable(2, "sigma_z")
    lowest = 1.0
    for g0, hs, page_c in ((0.005, (0.7, 0.8, 0.9, 1.0), 0.899), (0.006, (0.7, 0.9), 0.879)):
        for h in hs:
            H = hamiltonian(2, h=h)
            sol = feedback_exact(bell(), H, L_ops, O, g0, 0.5, 10.0)
            f = lambda t, sol=sol: concurrence(sol.sol(t).reshape(4, 4)) * l1(sol.sol(t).reshape(4, 4)) / 3
            low, _ = refined_extremum(f, 10.0, sign=-1.0, grid=0.002)
            lowest = min(lowest, low)
            check(f"gamma_0 {g0}, h {h}: exact concurrence at t = 10", page_c, concurrence(sol.sol(10.0).reshape(4, 4)))
            plain = euler_final(bell(), H, g0, 0.5, False)
            if g0 == 0.005 and h in (0.7, 1.0):
                check(f"gamma_0 {g0}, h {h}: plain Euler purity at t = 10 (no state)", {0.7: 1.44, 1.0: 2.68}[h],
                      purity(plain), digits=2)
            neg = float(np.min(np.linalg.eigvalsh(0.5 * (plain + plain.conj().T))))
            claim(f"gamma_0 {g0}, h {h}: the plain Euler matrix at t = 10 has a negative eigenvalue", neg < 0,
                  f"lowest eigenvalue {neg:.3f}")
            print(f"    the tool (Euler with clipping): concurrence {concurrence(euler_final(bell(), H, g0, 0.5, True)):.6f}")
            print(f"    exact: lowest concurrence * l1/3 on t <= 10 {low:.6f}")
    check("lowest exact concurrence * l1/3 over the sweep", 0.277, lowest)
    tool_low = min(euler_lowest(hamiltonian(2, h=h), g0) for g0, hs in
                   ((0.005, (0.7, 0.8, 0.9, 1.0)), (0.006, (0.7, 0.9))) for h in hs)
    check("lowest concurrence * l1/3 over every step of the tool's six runs", 0.289, tool_low)
    check("the historical copy's row gamma_0 0.006, h 1.0: C_final (Euler with clipping)", 0.901,
          concurrence(euler_final(bell(), hamiltonian(2, h=1.0), 0.006, 0.5, True)))
    print("  Simulation Evidence §2 (sigma_z, kappa 0.5, gamma_0 0.005, h 0.9, t_max 5, the tool's reading), exact:")
    for label, rho0, n, H, p_end, p_max, p_above, pur_end in (
            ("Bell+, chain", bell(), 2, hamiltonian(2, h=0.9), 0.402, None, None, 0.951),
            ("GHZ3, ring", ghz(3), 3, hamiltonian(3, h=0.9, ring=True), 0.248, 0.493, 3.89, 0.864),
            ("W3, ring", w_state(3), 3, hamiltonian(3, h=0.9, ring=True), 0.421, None, None, 0.909)):
        L_n, O_n = jumps_and_observable(n, "sigma_z")
        sol = feedback_exact(rho0, H, L_n, O_n, 0.005, 0.5, 5.0)
        fun = lambda t, sol=sol, n=n: cpsi(sol.sol(t).reshape(2**n, 2**n), n)
        check(f"{label}: C*Psi at t = 5", p_end, fun(5.0))
        # the tool's purity of the same run is gated in delta_calc_feedback_runs.py
        check(f"{label}: purity at t = 5", pur_end, purity(sol.sol(5.0).reshape(2**n, 2**n)))
        xs = [np.real(np.trace(sol.sol(t).reshape(2**n, 2**n) @ O_n)) for t in (0.0, 5.0)]
        print(f"    <X0X1> from {xs[0]:.6f} to {xs[1]:.6f}")
        if p_max is not None:
            check(f"{label}: largest C*Psi on [0, 5]", p_max, refined_extremum(fun, 5.0)[0])
            cross, _, _ = crossings(fun, 5.0)
            edges = [0.0] + cross + [5.0]
            above = sum(b - a for a, b in zip(edges, edges[1:]) if fun(0.5 * (a + b)) > 0.25)
            check(f"{label}: time above 1/4 in [0, 5]", p_above, above, digits=2)
            check(f"{label}: <X0X1> at t = 5", 0.0, xs[1], digits=6)
            sol6 = feedback_exact(rho0, H, L_n, O_n, 0.005, 0.5, 6.0)
            near = [c for c in crossings(lambda t: cpsi(sol6.sol(t).reshape(8, 8), 3), 6.0, grid=0.002)[0]
                    if 4.5 < c < 6.0]
            check(f"{label}: exact downward crossing of 1/4 near t = 5", 4.998, near[0])
            check(f"{label}: exact upward crossing after it", 5.48, near[1], digits=2)
            rows = tool_run(rho0, n, H, 0.005, 0.5, "sigma_z", t_max=5.2, record_every=1)
            after = next(i for i, r in enumerate(rows) if r[0] > 4.5 and r[2] * r[3] < 0.25)
            check(f"{label}: the tool's last record above 1/4", 5.01, rows[after - 1][0], digits=2)
            check(f"{label}: the tool's first record below 1/4", 5.02, rows[after][0], digits=2)
        if n == 3:
            if p_max is None:
                check(f"{label}: <X0X1> at t = 5", 0.623, xs[1])
            plain = tool_run(rho0, n, H, 0.005, 0.5, "sigma_z", clip=False)[-1][4]
            clipped = tool_run(rho0, n, H, 0.005, 0.5, "sigma_z")[-1][4]
            check(f"{label}: the tool's <X0X1> at t = 5, with clipping", 0.07 if p_max else 0.70, clipped, digits=2)
            if p_max is not None:
                check(f"{label}: the same Euler step without clipping, as the tool records it", 0.0, plain, digits=6)
            else:
                check(f"{label}: the same Euler step without clipping", 0.623, plain)


def section_5():
    print("\n§5, the Euler step")
    H = hamiltonian(2)
    rho = bell()
    dt = 0.01
    comm = H @ rho - rho @ H
    step = rho - 1j * dt * comm
    gain = purity(step) - purity(rho)
    check("one step at gamma 0 from Bell+: purity gain (units of 1e-4)", 2.0, gain * 1e4, digits=1)
    first = np.min(np.linalg.eigvalsh(0.5 * (step + step.conj().T)))
    check("  and its lowest eigenvalue (units of 1e-4)", -1.0, first * 1e4, digits=1)
    print(f"    dt^2 ||[H, rho]||_F^2 = {dt**2 * np.real(np.trace(comm.conj().T @ comm)):.6e}")
    for t_max, page in ((5.0, 1.111), (10.0, 1.246)):
        rows = tool_run(bell(), 2, H, 0.0, 0.0, "sigma_z", t_max=t_max, clip=False)
        check(f"gamma 0, no clipping: purity at t = {t_max:g}", page, rows[-1][1])
    for g0, kappa, h, t_max, page_top in ((0.1, 0.99, 0.5, 5.0, 1.047),):
        rho, top = bell(), 1.0
        Lz, Oz = jumps_and_observable(2, "sigma_z")
        Hh = hamiltonian(2, h=h)
        for _ in range(int(round(t_max / dt))):
            o = float(np.clip(np.real(np.trace(rho @ Oz)), -1.0, 1.0))
            g = max(0.0, g0 * (1 - kappa * o))
            drho = -1j * (Hh @ rho - rho @ Hh)
            for L in Lz:
                drho += g * (L @ rho @ L - rho)
            rho = rho + dt * drho
            top = max(top, purity(rho))
        check(f"plain Euler, gamma_0 {g0}, kappa {kappa}, h {h}: largest purity, t <= {t_max:g}", page_top, top,
              digits=3 if page_top < 1.5 else 2)
    below = [round(0.05 * k, 2) for k in range(19)] + [0.925, 0.94, 0.949]
    tops = [plain_extremes(bell(), 2, H, 0.1, kappa, 5.0)[0] for kappa in below]
    claim("plain Euler, gamma_0 0.1: purity never exceeds 1 (t <= 5) at kappa 0 to 0.9 in steps of 0.05 and at "
          "0.925, 0.94, 0.949", max(tops) <= 1.0, f"largest purity - 1 {max(tops) - 1:.1e}")
    for kappa in (0.95, 0.99):
        top = plain_extremes(bell(), 2, H, 0.1, kappa, 5.0)[0]
        claim(f"plain Euler, gamma_0 0.1, kappa {kappa}: purity passes 1 (t <= 5)", top > 1.0, f"{top - 1:.1e}")
    for kappa in (0.0, 0.5):
        rho, top = bell(), 1.0
        for _ in range(500):
            o = float(np.clip(np.real(np.trace(rho @ Oz)), -1.0, 1.0))
            g = 0.1 * (1 - kappa * o)
            drho = -1j * (H @ rho - rho @ H)
            for L in Lz:
                drho += g * (L @ rho @ L - rho)
            rho = rho + dt * drho
            top = max(top, purity(rho))
        claim(f"plain Euler, gamma_0 0.1, kappa {kappa}: purity never exceeds its start, 1 (t <= 5; necessary only)",
              top <= 1.0, f"{top - 1:.1e}")
    rho, low = bell(), 0.0
    for _ in range(500):
        drho = -1j * (H @ rho - rho @ H)
        for L in Lz:
            drho += 0.1 * (L @ rho @ L - rho)
        rho = rho + dt * drho
        low = min(low, np.min(np.linalg.eigvalsh(0.5 * (rho + rho.conj().T))))
    check("plain Euler, sigma_z, gamma_0 0.1, kappa 0: lowest eigenvalue, t <= 5", -0.004, low)
    print("  the plain step in the runs with an active Hamiltonian, kappa 0, 0.5, 1:")
    active = [(f"sweep row gamma_0 {g0}, h {h}", bell(), 2, hamiltonian(2, h=h), g0, 10.0)
              for g0, hs in ((0.005, (0.7, 0.8, 0.9, 1.0)), (0.006, (0.7, 0.9))) for h in hs]
    active += [("SE §2 Bell+, chain, h 0.9", bell(), 2, hamiltonian(2, h=0.9), 0.005, 5.0),
               ("SE §2 GHZ3, ring, h 0.9", ghz(3), 3, hamiltonian(3, h=0.9, ring=True), 0.005, 5.0),
               ("SE §2 W3, ring, h 0.9", w_state(3), 3, hamiltonian(3, h=0.9, ring=True), 0.005, 5.0)]
    for label, rho0, n, Hh, g0, t_max in active:
        tops = [plain_extremes(rho0, n, Hh, g0, k, t_max)[0] for k in (0.0, 0.5, 1.0)]
        claim(f"{label}: the plain step's purity passes 1 at kappa 0, 0.5 and 1", min(tops) > 1.0,
              " ".join(f"{t:.3f}" for t in tops))
    print("  where H does not turn the state (h = 0):")
    for label, rho0, n, Hh, g0 in (("Bell+, J 1, h 0, gamma_0 0.1", bell(), 2, hamiltonian(2, h=0.0), 0.1),
                                   ("SE §2 Bell+, J 0, h 0, gamma_0 0.005", bell(), 2, hamiltonian(2, J=0.0, h=0.0), 0.005),
                                   ("Mediator §3.1 GHZ3, ring, h 0, gamma_0 0.05", ghz(3), 3,
                                    hamiltonian(3, h=0.0, ring=True), 0.05)):
        tops, lows = zip(*[plain_extremes(rho0, n, Hh, g0, k, 5.0) for k in (0.0, 0.5, 1.0)])
        claim(f"{label}: the plain step stays a state at kappa 0, 0.5, 1 (purity <= 1; eigenvalues >= -100 eps, a bound)",
              max(tops) <= 1.0 and min(lows) >= -100 * EPS, f"purity {max(tops):.6f}, lowest {min(lows):.1e}")
    for kappa in (0.0, 0.95):
        rows = tool_run(bell(), 2, H, 0.1, kappa, "xx", clip=False)
        check(f"plain Euler, xx, kappa {kappa}: purity at t = 5", 1.111, rows[-1][1])
    print("  the first plain step in the other runs where H turns the state (kappa 0.5):")
    for label, rho0, n, jump, Hh in (("sigma_z, Bell+", bell(), 2, "sigma_z", H),
                                     ("yy, Bell+", bell(), 2, "yy", H), ("zz, Bell+", bell(), 2, "zz", H),
                                     ("x_pairs, GHZ4", ghz(4), 4, "x_pairs", hamiltonian(4)),
                                     ("x_pairs, GHZ4, h 0.2 (the 19:07 call)", ghz(4), 4, "x_pairs",
                                      hamiltonian(4, h=0.2)),
                                     ("x_pairs, GHZ5", ghz(5), 5, "x_pairs", hamiltonian(5)),
                                     ("x_pairs, GHZ6", ghz(6), 6, "x_pairs", hamiltonian(6))):
        L_j, O_j = jumps_and_observable(n, jump)
        g = 0.1 * (1 - 0.5 * float(np.real(np.trace(rho0 @ O_j))))
        drho = -1j * (Hh @ rho0 - rho0 @ Hh)
        for L in L_j:
            drho += g * (L @ rho0 @ L - rho0)
        step1 = rho0 + dt * drho
        first = float(np.min(np.linalg.eigvalsh(0.5 * (step1 + step1.conj().T))))
        claim(f"{label}: the first plain step opens a negative eigenvalue (below -100 eps)", first < -100 * EPS,
              f"{first:.1e}")
    check("<X0X1> at t = 5, gamma_0 0.005, kappa 0.5: closed form", 0.9500416, closed_form_x(1.0, 0.005, 0.5, 5.0), digits=7)
    x = 1.0
    for _ in range(500):
        x = x * (1 - 4 * 0.005 * (1 - 0.5 * x) * 0.01)
    check("  the scalar Euler recursion x <- x (1 - 4 gamma_eff dt)", 0.9500415, x, digits=7)
    L_ops, O = jumps_and_observable(2, "sigma_z")
    for J, h, page_clip in ((0.0, 0.0, None), (1.0, 0.5, 0.951246), (1.0, 0.9, 0.953819)):
        Hh = hamiltonian(2, J=J, h=h)
        plain = np.real(np.trace(euler_final(bell(), Hh, 0.005, 0.5, False, t_max=5.0) @ O))
        claim(f"J {J:g}, h {h:g}: plain Euler <X0X1> is the recursion (|difference| <= 500 eps, its rounding)",
              abs(plain - x) <= 500 * EPS, f"difference {plain - x:.1e}")
        clipped = tool_run(bell(), 2, Hh, 0.005, 0.5, "sigma_z")[-1][4]
        if page_clip is None:
            claim("J 0, h 0: with nothing to turn, the clipping leaves <X0X1> on the recursion (6 digits)",
                  round(clipped, 6) == round(x, 6), f"{clipped:.6f}")
        else:
            check(f"J {J:g}, h {h:g}: Euler with clipping", page_clip, clipped, digits=6)
    print("  the adaptive solver near kappa = 1 (Bell+, sigma_z, gamma_0 0.1, h 0.5, t <= 5):")
    L_ops, O = jumps_and_observable(2, "sigma_z")
    for kappa in (0.99, 1.0):
        ratios = []
        for rtol in (1e-10, 1e-12):
            sol = feedback_exact(bell(), H, L_ops, O, 0.1, kappa, 5.0, rtol=rtol)
            states = [sol.sol(t).reshape(4, 4) for t in np.linspace(0, 5, 501)]
            low = min(np.min(np.linalg.eigvalsh(0.5 * (r + r.conj().T))) for r in states)
            high = max(purity(r) for r in states)
            ratios.append((abs(min(low, 0.0)) / rtol, abs(max(high - 1, 0.0)) / rtol))
        claim(f"DOP853, kappa {kappa}: any departure from the states stays within 10 rtol at rtol 1e-10 and 1e-12",
              all(a <= 10 and b <= 10 for a, b in ratios),
              " ".join(f"({a:.2f}, {b:.2f})" for a, b in ratios))


def main():
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    section_identities()
    section_2()
    section_4()
    section_5()
    section_81()
    section_82()
    matched = sum(ok for _, ok in VALUES)
    quoted = sum(ok for _, ok in QUOTES)
    logged_ok = sum(ok for _, ok in LOGGED)
    held = sum(ok for _, ok in CLAIMS)
    print(f"\n{matched} of {len(VALUES)} printed values match at their pages' precision; {quoted} of {len(QUOTES)}"
          f" quotes of the log; {logged_ok} of {len(LOGGED)} logged values regenerate;"
          f" {held} of {len(CLAIMS)} identities, controls, convergences and claims hold.")
    failed = [label for label, ok in VALUES + QUOTES + LOGGED + CLAIMS if not ok]
    for label in failed:
        print(f"  not reproduced: {label}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
