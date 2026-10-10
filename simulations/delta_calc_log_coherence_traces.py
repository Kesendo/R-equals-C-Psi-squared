"""
The log-coherence traces of February 19 in the retired tool's own reading, and exactly
=====================================================================================
experiments/ALGEBRAIC_EXPLORATION.md (Finding 1, and the slope spread of lambda among its rejected
findings) and docs/historical/CORE_ALGEBRA.md §11 report
rate ratios gamma_eff/gamma_base of xi = ln(Psi) from four runs of the retired delta_calc MCP tool
(its source is kept outside the repo and is transcribed here, not imported), and a fifth run under
its memory-kernel feedback. The Claude Desktop log of the chat's MCP holds the five calls
(2026-02-19, 03:31 to 03:47 UTC), all simulate_dynamic_lindblad with gamma_base 0.1, dt 0.05,
t_max 5 and the mutual-purity bridge:

  - Bell+ on the open Heisenberg chain, noise_type local;
  - W_3 on the Heisenberg ring, noise_type local;
  - Bell+, noise_type collective;
  - Bell+, noise_type local with jump_operator sigma_x;
  - Bell+, noise_type memory_kernel_feedback, kappa 0.5, tau 1.

No call passes h, so the field is off (the tool's default h = 0, J = 1). The tool source
read here is dated after the calls (evolution.py and simulations.py 2026-03-01, states.py
2026-02-19 in the afternoon); that every value
the log shows regenerates is the evidence that these paths did not change. The tool, outside its
two feedback modes, runs the scalar law: the rate gamma_base * C(rho) with C
the bridge, here the mutual purity (the geometric mean of the single-site purities), in front of
sigma_z on every site (local, whatever jump_operator says) or the single jump sum_k sigma_z^(k)
(collective). The memory kernel puts gamma_base * (1 - kappa * M / tau) in front of sigma_z on every
site, M <- M * exp(-dt/tau) + <X_0 X_1> * dt. Each step is an Euler step, then the Hermitian part,
the negative eigenvalues set to zero and the trace renormalized; it records every max(1, int(0.1/dt))
steps, rounded to six digits, Psi = l1(rho)/(d - 1).

The transcription regenerates the values the log shows of each call (the log keeps only the tail of
the last array of each response). This script then reads the rate ratios three ways: in the tool's run, in exact
propagation of the same law (DOP853), and under the plain Lindblad equation at the fixed rate
gamma_base (expm), and states why each is what it is. In all four constant-rate runs H does not
act: [H, rho] = 0 on the initial state and on every dephased image of it, checked exactly on
integer matrices. A real sigma_x channel, which the tool never ran, is propagated as a control.

Import-inert; prints only, and exits with status 1 if any check fails.
"""

from functools import reduce
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
EPS = np.finfo(float).eps
DT, GAMMA, T_MAX = 0.05, 0.1, 5.0

# The tails the log shows of each call (2026-02-19 UTC): Psi for the first four, <O_int> for the
# memory kernel, and its memory_integral_final as printed.
LOG = {
    "local Bell+ (03:31)": [0.140445, 0.13765, 0.134911, 0.132226, 0.129595, 0.127016, 0.124488, 0.122011],
    "local W3, ring (03:32)": [0.099953, 0.097744, 0.095584, 0.093472],
    "collective Bell+ (03:34)": [0.058657, 0.056334, 0.054104, 0.051961, 0.049903, 0.047927, 0.046029,
                                 0.044207],
    "local Bell+, jump_operator sigma_x (03:35)": [0.140445, 0.13765, 0.134911, 0.132226, 0.129595,
                                                   0.127016, 0.124488, 0.122011],
}
LOG_MEMORY_O = [0.259785, 0.251347, 0.243133, 0.235141, 0.227367, 0.219807, 0.212458]
LOG_MEMORY_INTEGRAL = 0.3118303352358058

CLAIMS = []
CHECKS = []


def site(P, k, n):
    return reduce(np.kron, [P if i == k else I2 for i in range(n)])


def heisenberg(n, ring):
    """J = 1: an integer matrix."""
    H = np.zeros((2**n, 2**n), dtype=complex)
    for i in range(n if ring else n - 1):
        j = (i + 1) % n
        for P in (X, Y, Z):
            H += site(P, i, n) @ site(P, j, n)
    return H


def one_site(rho, n, k):
    t = rho.reshape([2] * (2 * n))
    keep = [k, n + k]
    rest = [a for a in range(2 * n) if a not in keep]
    t = np.transpose(t, keep + rest).reshape(2, 2, -1)
    m = 2 ** (n - 1)
    t = t.reshape(2, 2, m, m)
    return np.einsum("abii->ab", t)


def purity(r):
    return float(np.real(np.trace(r @ r)))


def mutual_purity(rho, n):
    return float(np.prod([purity(one_site(rho, n, k)) for k in range(n)]) ** (1.0 / n))


def psi(r):
    return float((np.abs(r).sum() - np.abs(np.diag(r)).sum()) / (r.shape[0] - 1))


def pure(v):
    v = np.asarray(v, dtype=complex)
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


def bell():
    return pure([1, 0, 0, 1])


def w_state(n):
    return pure([1 if bin(i).count("1") == 1 else 0 for i in range(2**n)])


def dissipate(rho, L_ops, g):
    out = np.zeros_like(rho)
    for L in L_ops:
        LdL = L.conj().T @ L
        out += g * (L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL))
    return out


def tool_run(rho0, n, H, L_ops, memory=None):
    """The tool's run, transcribed: the recorded Psi, the recorded <X_0 X_1>, and M at the end."""
    O = site(X, 0, n) @ site(X, 1, n)
    rho, M = rho0.copy(), 0.0
    every, steps = max(1, int(0.1 / DT)), int(T_MAX / DT)
    psis, os_ = [], []
    for step in range(steps + 1):
        if step % every == 0:
            psis.append(round(psi(rho), 6))
            os_.append(round(float(np.clip(np.real(np.trace(rho @ O)), -1, 1)), 6))
        if step == steps:
            break
        if memory:
            kappa, tau = memory
            M = M * np.exp(-DT / tau) + float(np.clip(np.real(np.trace(rho @ O)), -1, 1)) * DT
            g = max(0.0, GAMMA * (1 - kappa * M / tau))
        else:
            g = GAMMA * mutual_purity(rho, n)
        rho = rho + DT * (-1j * (H @ rho - rho @ H) + dissipate(rho, L_ops, g))
        rho = 0.5 * (rho + rho.conj().T)
        w, v = np.linalg.eigh(rho)
        rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
        rho /= np.real(np.trace(rho))
    return psis, os_, M


def exact_run(rho0, n, H, L_ops, law, memory=None, rtol=1e-12):
    """DOP853 on rho (and M for the memory kernel): law 'bridge' is gamma * C(rho), 'memory' the
    kernel's rate."""
    d = rho0.shape[0]
    O = site(X, 0, n) @ site(X, 1, n)

    def f(_, y):
        r = y[:-1].reshape(d, d)
        if law == "bridge":
            g = GAMMA * mutual_purity(r, n)
            dM = 0.0
        else:
            kappa, tau = memory
            M = float(np.real(y[-1]))
            g = max(0.0, GAMMA * (1 - kappa * M / tau))
            dM = float(np.real(np.trace(r @ O))) - M / tau
        dr = -1j * (H @ r - r @ H) + dissipate(r, L_ops, g)
        return np.concatenate([dr.reshape(-1), [dM]])

    y0 = np.concatenate([rho0.reshape(-1), [0.0]]).astype(complex)
    return solve_ivp(f, (0, T_MAX), y0, method="DOP853", rtol=rtol, atol=rtol * 1e-2, dense_output=True)


def claim(label, ok, detail=""):
    CLAIMS.append((label, ok))
    print(f"  {label} {detail}  {'holds' if ok else 'FAILS'}")


def check(label, page, value, digits=3):
    ok = round(value, digits) == round(page, digits)
    CHECKS.append((label, ok))
    print(f"  {label:<70s} page {page:.{digits}f}  computed {value:.{digits + 3}f}  {'match' if ok else 'DIFFERS'}")


def ratio(psi0, psi5):
    return -(np.log(psi5) - np.log(psi0)) / (T_MAX * GAMMA)


def main():
    configs = [
        ("local Bell+ (03:31)", bell(), 2, heisenberg(2, False), [site(Z, k, 2) for k in range(2)], 4, 0.5,
         2.010),
        ("local W3, ring (03:32)", w_state(3), 3, heisenberg(3, True), [site(Z, k, 3) for k in range(3)], 4,
         5 / 9, 2.235),
        ("collective Bell+ (03:34)", bell(), 2, heisenberg(2, False), [site(Z, 0, 2) + site(Z, 1, 2)], 8, 0.5,
         4.041),
    ]
    print("The logged calls, regenerated (every value the log shows):")
    for label, rho0, n, H, L_ops, *_ in configs:
        psis, _, _ = tool_run(rho0, n, H, L_ops)
        tail = LOG[label]
        claim(f"{label}: the {len(tail)} logged Psi values regenerate", psis[-len(tail):] == tail,
              f"(end {psis[-1]})")
    claim("the logged tails of the sigma_x-labelled call and the sigma_z call, as transcribed, are identical",
          LOG["local Bell+, jump_operator sigma_x (03:35)"] == LOG["local Bell+ (03:31)"])

    print("\nH does nothing in these runs: the dynamics stays in span{rho0, D(rho0)}, D the dissipator at unit"
          " rate, and H commutes with both (integer matrices, exact):")
    for label, rho0, n, H, L_ops, *_ in configs:
        scale = 2 if n == 2 else 3
        base = np.rint(np.real(rho0 * scale)).astype(int)
        Li = [np.rint(np.real(L)).astype(int) for L in L_ops]

        def D(m):
            out = np.zeros_like(m)
            for L in Li:
                LdL = L.T @ L
                out = out + 2 * (L @ m @ L.T) - (LdL @ m + m @ LdL)
            return out  # twice the dissipator, to stay in the integers

        d1, d2 = D(base), D(D(base))
        Hi = np.rint(np.real(H)).astype(int)
        worst = max(int(np.abs(Hi @ m - m @ Hi).max()) for m in (base, d1))
        a_, b_ = d1.reshape(-1), d2.reshape(-1)
        closed = int(np.abs(np.outer(a_, b_) - np.outer(b_, a_)).max())  # every 2x2 minor of the pair
        turned = np.zeros_like(base)
        turned[1, 1] = 1  # |0..01><0..01|, which the bonds move
        control = int(np.abs(Hi @ turned - turned @ Hi).max())
        claim(f"{label}: [H, rho0] = [H, D rho0] = 0 and D(D rho0) is a multiple of D rho0 "
              f"(control [H, |0..01><0..01|] has {control})", worst == 0 and closed == 0 and control != 0,
              f"(largest entries {worst}, {closed})")

    print("\nThe February slope variation, read as the spread (max - min)/mean of xi's slope over ten windows"
          " of 0.5:")
    for label, rho0, n, H, L_ops, *_ in configs:
        psis, _, _ = tool_run(rho0, n, H, L_ops)
        s = -np.diff(np.log(psis[::5])) / 0.5
        print(f"    {label}: {100 * np.ptp(s) / s.mean():.4f}%")
        if "Bell+" in label:
            check(f"{label}: the spread, in per cent (February: 0.009)", 0.009, 100 * np.ptp(s) / s.mean())
        else:
            print(f"    (W3: February printed 0.002; this reading gives {100 * np.ptp(s) / s.mean():.3f})")
    print("  lambda = -ln(P Psi) with P the purity, on Bell+ P = (1 + 9 Psi^2)/2, the same windows:")
    for label, rho0, n, H, L_ops, *_ in configs:
        if "Bell+" not in label:
            continue
        p = np.array(tool_run(rho0, n, H, L_ops)[0])
        lam = -np.log((1 + 9 * p * p) / 2 * p)
        s = np.diff(lam[::5]) / 0.5
        page = 44 if "local" in label else 64
        check(f"{label}: lambda's spread, in per cent", page, 100 * np.ptp(s) / s.mean(), digits=0)

    print("\nThe rate ratio gamma_eff/gamma_base of xi = ln Psi, three ways:")
    for label, rho0, n, H, L_ops, k, C, page in configs:
        psis, _, _ = tool_run(rho0, n, H, L_ops)
        check(f"{label}: the tool's ratio (Euler, dt 0.05)", page, ratio(psis[0], psis[-1]))
        euler = -np.log(1 - k * GAMMA * C * DT) / (DT * GAMMA)
        print(f"    its closed Euler form -ln(1 - {k} gamma C dt)/(dt gamma) = {euler:.6f}")
        sol = exact_run(rho0, n, H, L_ops, "bridge")
        r5 = sol.sol(T_MAX)[:-1].reshape(rho0.shape)
        cdev = max(abs(mutual_purity(sol.sol(t)[:-1].reshape(rho0.shape), n) - C) for t in np.linspace(0, 5, 26))
        claim(f"{label}: exactly, C stays {C:.6f} (<= 100 eps)", cdev <= 100 * EPS, f"(largest |C - C0| {cdev:.1e})")
        ex = ratio(psi(rho0), psi(r5))
        claim(f"{label}: exactly, the ratio is {k} C = {k * C:.6f} (|deviation| <= 1e-9)",
              abs(ex - k * C) <= 1e-9, f"({ex:.12f})")
        plain = expm  # the plain Lindblad equation at the fixed rate gamma
        d = rho0.shape[0]
        Lv = -1j * (np.kron(np.eye(d), H) - np.kron(H.T, np.eye(d)))
        for L in L_ops:
            LdL = L.conj().T @ L
            Lv += GAMMA * (np.kron(L.conj(), L) - 0.5 * np.kron(np.eye(d), LdL) - 0.5 * np.kron(LdL.T, np.eye(d)))
        r5p = (plain(Lv * T_MAX) @ rho0.reshape(-1, order="F")).reshape(d, d, order="F")
        claim(f"{label}: the plain Lindblad equation at gamma gives {k} (|deviation| <= 1e-9)",
              abs(ratio(psi(rho0), psi(r5p)) - k) <= 1e-9, f"({ratio(psi(rho0), psi(r5p)):.12f})")

    print("\nA real sigma_x channel, which the tool never ran (the scalar law, Bell+, Heisenberg chain):")
    sol = exact_run(bell(), 2, heisenberg(2, False), [site(X, k, 2) for k in range(2)], "bridge")
    dev = max(abs(psi(sol.sol(t)[:-1].reshape(4, 4)) - 1 / 3) for t in np.linspace(0, 5, 26))
    claim("Psi stays 1/3 (<= 100 eps)", dev <= 100 * EPS, f"(largest |Psi - 1/3| {dev:.1e})")
    pur = purity(sol.sol(T_MAX)[:-1].reshape(4, 4))
    claim("  while the state dephases: purity at t = 5 is (1 + e^-2)/2 (|deviation| <= 1e-9)",
          abs(pur - (1 + np.exp(-2.0)) / 2) <= 1e-9, f"({pur:.6f})")

    print("\nThe memory kernel (kappa 0.5, tau 1), Bell+ on the chain:")
    psis, os_, M = tool_run(bell(), 2, heisenberg(2, False), [site(Z, k, 2) for k in range(2)], memory=(0.5, 1.0))
    claim(f"the {len(LOG_MEMORY_O)} logged <X_0 X_1> values regenerate", os_[-len(LOG_MEMORY_O):] == LOG_MEMORY_O,
          f"(end {os_[-1]})")
    claim("memory_integral_final regenerates bit for bit", M == LOG_MEMORY_INTEGRAL, f"({M!r})")
    slopes = -np.diff(np.log(psis)) / 0.1
    print(f"    the tool's run: xi's slope over the recorded intervals {slopes[0]:.3f} at the start, "
          f"{slopes.min():.3f} at its smallest, {slopes[-1]:.3f} at the end (no page prints these)")
    s5 = -np.diff(np.log(psis[::5])) / 0.5
    check("the spread of the slope over ten windows of 0.5, in per cent (February: 24.5)", 24.5,
          100 * np.ptp(s5) / s5.mean(), digits=1)
    sol = exact_run(bell(), 2, heisenberg(2, False), [site(Z, k, 2) for k in range(2)], "memory", memory=(0.5, 1.0))
    rate = lambda t: 4 * GAMMA * (1 - 0.5 * float(np.real(sol.sol(t)[-1])) / 1.0)
    fine = np.linspace(0, 5, 5001)
    rates = [rate(t) for t in fine]
    check("exactly, xi's decay rate at t = 0", 0.400, rates[0])
    check("exactly, xi's decay rate at its smallest", 0.284, min(rates))
    claim("  gamma(t) = rate/4 never falls below 0.071 (so it stays positive: plain time-dependent dephasing)",
          min(rates) / 4 > 0.071, f"(smallest {min(rates) / 4:.7f})")
    check("  where it falls (t)", 1.7, fine[int(np.argmin(rates))], digits=1)
    check("  at t = 5", 0.340, rates[-1])
    worst = 0.0
    L2 = [site(Z, k, 2) for k in range(2)]
    H2 = heisenberg(2, False)
    for t in np.linspace(0, 5, 11):
        r = sol.sol(t)[:-1].reshape(4, 4)
        dr = -1j * (H2 @ r - r @ H2) + dissipate(r, L2, rate(t) / 4)
        worst = max(worst, abs(np.real(dr[0, 3]) / np.real(r[0, 3]) + rate(t)) / rate(t))
    claim("  and it is xi's decay rate: d ln rho_03/dt from the generator equals -rate, as H does nothing "
          "(relative, <= 100 eps)",
          worst <= 100 * EPS, f"(largest relative difference {worst:.1e})")
    ts = np.arange(0, 5.0001, 0.1)
    ps = [psi(sol.sol(t)[:-1].reshape(4, 4)) for t in ts]
    s5e = -np.diff(np.log(ps[::5])) / 0.5
    check("  the spread over ten windows of 0.5, in per cent", 24.3, 100 * np.ptp(s5e) / s5e.mean(), digits=1)

    held = sum(ok for _, ok in CLAIMS)
    matched = sum(ok for _, ok in CHECKS)
    print(f"\n{matched} of {len(CHECKS)} printed values match at the page's precision; "
          f"{held} of {len(CLAIMS)} identities, regenerations and claims hold.")
    if held < len(CLAIMS) or matched < len(CHECKS):
        sys.exit(1)


if __name__ == "__main__":
    main()
