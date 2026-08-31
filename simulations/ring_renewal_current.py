"""Gate-first: F126's renewal ladder on the cracked-ring CURRENT (the warble's dressed clock).

Verifier for the renewal section of experiments/THE_CRACKED_BELL.md (2026-08-31).
Plays the open point of THE_CRACKED_BELL x F126: the (1,1)
warble equation rho-dot = -i[h,rho] - Gamma(rho - diag rho), Gamma = 4*gamma, IS
F126's generator for ANY h (the derivation never uses the chain). Splitting at the
diagonal refill and resumming gives, for the circulation I(t) = Tr[M rho],
M = the antisymmetric bond-imaginary reader:

    I(t) = e^{-Gamma t} * [ I_free(t) + Gamma * int_0^t ds  sum_m C_m(t-s) * ntil_m(s) ]

    ntil_m(t) = e^{+Gamma t} n_m(t)  solves the PLAIN-convolution renewal ladder
    ntil = n_free + Gamma * (K * ntil),   K_{ab}(tau) = |G_{ab}(tau)|^2,  G = e^{-i h tau}

    C_m(tau) = <m| G(tau)^dag M G(tau) |m>   (the REBIRTH KERNEL: the current a walker
               reborn at seat m re-develops after free time tau; odd in tau, zero at 0)

Structure facts checked analytically first:
  - perfect ring: C_m == 0 would follow from site symmetry; in any case the ladder
    then gives I = I(0) e^{-Gamma t}, the known exact answer (populations frozen).
  - crack reflection R: j -> N-1-j preserves the cracked h, flips M: C_{R(m)} = -C_m.
  - I^{(0)} = e^{-Gamma t} I_free has GAMMA-INDEPENDENT zeros: the never-caught wave
    cannot advance the clock; every bit of the dressing must sit in the reborn terms.

And the FIRST-ORDER DRESSED CLOCK, closed in gamma-free objects: expanding the ladder
to one rebirth and locating the zero,

    T_zero(Gamma) = T_0 - Gamma * V(T_0) / I_free'(T_0) + O(Gamma^2),
    V(T_0) = int_0^{T_0} ds  sum_m C_m(T_0 - s) * n_free_m(s),

the first rebirth's kernel integrated against the free populations. At N=8, m=1,
delta=0.15: V(T_0) = 1.2743, I_free'(T_0) = 0.05081, so dT_zero/dGamma = -25.08;
the advance coefficient c = 25.08*DeltaE/T_0 = 0.0956 in units of Gamma/DeltaE with
the EXACT pair split DeltaE = 0.077392 (the first-order book 4*delta*J/N = 0.075
reads 0.0927: pin the book). The prediction's residual against the full ladder is
2e-3 at gamma = 0.002 and falls as O(Gamma^2) (STAGE E).

Gates: route equality vs the exact expm superoperator (with a dt-halving order check),
the committed T_zero pins (20.295 / 19.347 / 16.477), C-kernel symmetry, delta=0 control,
and the first-order law with an O(Gamma^2)-modeled band across a Gamma-halving.
Reads: generation-resolved clock (truncate the ladder at j), reborn share at the deepest
reversal, small-Gamma sweep of the advance.
"""

import numpy as np
from scipy.linalg import expm
from scipy.signal import fftconvolve

J = 1.0
N, M_SEED, DELTA = 8, 1, 0.15
FAIL = []


def gate(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


def ring_h(n, delta):
    h = np.zeros((n, n))
    for a in range(n - 1):
        h[a, a + 1] = h[a + 1, a] = J
    h[n - 1, 0] = h[0, n - 1] = J * (1.0 - delta)
    return h


def plane_wave(n, m):
    return np.exp(2j * np.pi * m * np.arange(n) / n) / np.sqrt(n)


def pair_split(n, delta, m):
    """Exact splitting of the m-pair, matched by proximity (as in cracked_bell_gate.py)."""
    ev_ = np.linalg.eigvalsh(ring_h(n, delta))
    e0_ = 2.0 * J * np.cos(2.0 * np.pi * m / n)
    idx_ = np.argsort(np.abs(ev_ - e0_))[:2]
    pr = np.sort(ev_[idx_])
    return pr[1] - pr[0], pr - e0_


def reader(n):
    m = np.zeros((n, n), dtype=complex)
    for a in range(n):
        m[(a + 1) % n, a] += 1.0 / (2j)
        m[a, (a + 1) % n] -= 1.0 / (2j)
    return m


def kernels(n, delta, ts):
    """G-derived arrays on the grid: K[t,a,b] = |G_ab|^2, C[t,m], n_free[t,a], I_free[t]."""
    ev, vec = np.linalg.eigh(ring_h(n, delta))
    Mr = reader(n)
    psi0 = plane_wave(n, M_SEED)
    K = np.empty((len(ts), n, n))
    C = np.empty((len(ts), n))
    n_free = np.empty((len(ts), n))
    I_free = np.empty(len(ts))
    for k, t in enumerate(ts):
        G = (vec * np.exp(-1j * ev * t)) @ vec.T.conj()
        K[k] = np.abs(G) ** 2
        GMG = G.conj().T @ Mr @ G
        C[k] = np.real(np.diag(GMG))
        psi_t = G @ psi0
        n_free[k] = np.abs(psi_t) ** 2
        I_free[k] = np.real(psi_t.conj() @ Mr @ psi_t)
    return K, C, n_free, I_free


def ladder(n, delta, gamma, dt, t_max, j_keep=6, tol=1e-13):
    """Generation-summed ladder. Returns ts, I_total(t), I_gen[j<=j_keep](t), refill share fn."""
    Gam = 4.0 * gamma
    ts = np.arange(0.0, t_max + dt / 2, dt)
    K, C, n_free, I_free = kernels(n, delta, ts)

    def volt(sig2d):     # Gamma * int K(t-s) sig(s) ds, per target site (trapezoid)
        acc = np.zeros_like(sig2d)
        for b in range(n):
            acc += fftconvolve(K[:, :, b], sig2d[:, b:b + 1], axes=0)[: len(ts)]
        end0 = 0.5 * np.einsum('tab,b->ta', K, sig2d[0])
        endt = 0.5 * np.einsum('ab,tb->ta', K[0], sig2d)
        return Gam * dt * (acc - end0 - endt)

    def volt_current(sig2d):  # Gamma * int sum_m C_m(t-s) sig_m(s) ds
        f = np.zeros(len(ts))
        for b in range(n):
            f += fftconvolve(C[:, b], sig2d[:, b])[: len(ts)]
        end0 = 0.5 * (C @ sig2d[0])
        endt = 0.5 * (C[0] @ sig2d.T)      # C[0] = 0, kept for form
        return Gam * dt * (f - end0 - endt)

    gen = n_free.copy()                    # ntil^{(0)}
    ntil = gen.copy()
    I_gen = [I_free.copy()]                # Itil^{(0)}
    j = 0
    while True:
        j += 1
        i_this = volt_current(gen)         # walkers reborn exactly j times, kernel on gen j-1
        gen = volt(gen)                    # ntil^{(j)}
        ntil += gen
        if j <= j_keep:
            I_gen.append(i_this)
        if np.max(np.abs(gen)) < tol * max(1.0, np.max(np.abs(ntil))) and j > 3:
            break
        if j > 200:
            raise RuntimeError("ladder did not converge")
    Itil = I_free + volt_current(ntil)
    env = np.exp(-Gam * ts)
    return ts, env * Itil, [env * g for g in I_gen], env[:, None] * ntil, j


def expm_route(n, delta, gamma, dt, t_max):
    ts = np.arange(0.0, t_max + dt / 2, dt)
    h = ring_h(n, delta)
    eye = np.eye(n)
    L = -1j * (np.kron(h, eye) - np.kron(eye, h.T))
    for a in range(n):
        for b in range(n):
            if a != b:
                L[a * n + b, a * n + b] -= 4.0 * gamma
    P = expm(L * dt)
    psi = plane_wave(n, M_SEED)
    rho = np.outer(psi, psi.conj()).reshape(-1)
    Mr = reader(n)
    out = np.empty(len(ts))
    pops = np.empty((len(ts), n))
    for k in range(len(ts)):
        r = rho.reshape(n, n)
        out[k] = np.real(np.trace(Mr @ r))
        pops[k] = np.real(np.diag(r))
        rho = P @ rho
    return ts, out, pops


def t_zero(ts, cur):
    r = cur / cur[0]
    for k in range(1, len(ts)):
        if r[k - 1] > 0 >= r[k]:
            return ts[k - 1] + (ts[k] - ts[k - 1]) * r[k - 1] / (r[k - 1] - r[k])
    return np.nan


T_MAX, DT = 58.0, 0.01
print("=" * 78)
print("F126 renewal ladder on the cracked-ring current: N=8, m=1, delta=0.15")
print("=" * 78)

# structure gates first -- the EXACT route (case 1 of the no-rounding convention):
# both kernel theorems are operator identities, gated as entry-wise float equality.
R_ref = np.eye(N)[::-1]                      # j -> N-1-j, reflection through the cracked bond
T_shift = np.roll(np.eye(N), 1, axis=0)      # the cyclic shift of the perfect ring
h_c, h_0, Mr_ = ring_h(N, DELTA), ring_h(N, 0.0), reader(N)
gate("exact: R h R == h and R M R == -M (so C_{N-1-a} = -C_a is a theorem)",
     np.array_equal(R_ref @ h_c @ R_ref, h_c) and np.array_equal(R_ref @ Mr_ @ R_ref, -Mr_))
gate("exact: [h0, T] == 0 and [M, T] == 0 (so C == 0 on the perfect ring is a theorem)",
     np.array_equal(h_0 @ T_shift, T_shift @ h_0)
     and np.array_equal(Mr_ @ T_shift, T_shift @ Mr_))
ts_s = np.linspace(0.0, 10.0, 401)
K_s, C_s, _, _ = kernels(N, DELTA, ts_s)
_, C0_s, _, _ = kernels(N, 0.0, ts_s)
print(f"  (the eigh-route reads of the two theorems: max|C + R(C)| = "
      f"{np.max(np.abs(C_s + C_s[:, ::-1])):.1e}, max|C(delta=0)| = "
      f"{np.max(np.abs(C0_s)):.1e} -- float noise on exact statements)")

for gam, pin in [(0.01, 19.3470), (0.05, 16.4765)]:
    ts, I_lad, I_gens, n_lad, jmax = ladder(N, DELTA, gam, DT, T_MAX)
    ts_e, I_ex, pops_ex = expm_route(N, DELTA, gam, DT, T_MAX)
    err = np.max(np.abs(I_lad - I_ex))
    errp = np.max(np.abs(n_lad - pops_ex))
    # dt-halving order checks on BOTH route equalities
    ts2, I_lad2, _, n_lad2, _ = ladder(N, DELTA, gam, DT * 2, T_MAX)
    _, I_ex2, pops_ex2 = expm_route(N, DELTA, gam, DT * 2, T_MAX)
    err2 = np.max(np.abs(I_lad2 - I_ex2))
    errp2 = np.max(np.abs(n_lad2 - pops_ex2))
    gate(f"gamma={gam}: ladder == expm on the current (order check err(2dt)/err(dt))",
         err < 5e-5 and 2.5 < err2 / err < 6.5,
         f"err={err:.2e}, err(2dt)={err2:.2e}, ratio={err2 / err:.2f}, generations={jmax}")
    gate(f"gamma={gam}: ladder == expm on the populations (order check)",
         errp < 1e-6 and 2.5 < errp2 / errp < 6.5,
         f"err={errp:.2e}, err(2dt)={errp2:.2e}, ratio={errp2 / errp:.2f}")
    tz = t_zero(ts, I_lad)
    gate(f"gamma={gam}: ladder T_zero hits the committed pin {pin}",
         abs(tz - pin) < 0.02, f"T_zero(ladder) = {tz:.4f}")
    tz0 = t_zero(ts, I_gens[0])
    gate(f"gamma={gam}: the never-caught clock does NOT move (20.295, gamma-free)",
         abs(tz0 - 20.2953) < 0.02, f"T_zero(gen 0) = {tz0:.4f}")
    # generation-truncated clocks: where does the advance saturate?
    partial = I_gens[0].copy()
    row = [f"j=0: {tz0:.3f}"]
    for j in range(1, len(I_gens)):
        partial = partial + I_gens[j]
        row.append(f"<= {j}: {t_zero(ts, partial):.3f}")
    print(f"    generation-truncated T_zero (gamma={gam}): " + ", ".join(row))
    # refill share at the deepest reversal
    k_rev = int(np.argmin(I_lad / I_lad[0]))
    share = 1.0 - I_gens[0][k_rev] / I_lad[k_rev]
    print(f"    deepest reversal at t={ts[k_rev]:.2f}: reborn share = {share:.3f}")

# delta=0 control through the ladder
ts, I_lad, I_gens, _, _ = ladder(N, 0.0, 0.05, 0.02, 20.0)
ana = I_gens[0][0] * np.exp(-4 * 0.05 * ts)
gate("delta=0 control: ladder gives exactly I(0) e^{-Gamma t}",
     np.max(np.abs(I_lad - ana)) < 1e-12, f"max dev = {np.max(np.abs(I_lad - ana)):.1e}")

# the small-Gamma read: how the advance grows
print()
print("the dressing vs Gamma (full ladder, dt=0.02):")
print(f"{'gamma':>8} {'Gamma/DeltaE':>12} {'T_zero':>8} {'advance %':>10}")
tz_sweep = {}
for gam in [0.002, 0.005, 0.01, 0.02, 0.05]:
    ts, I_lad, _, _, _ = ladder(N, DELTA, gam, 0.02, T_MAX)
    tz_sweep[gam] = t_zero(ts, I_lad)
    print(f"{gam:>8} {4 * gam / (4 * DELTA * J / N):>12.3f} {tz_sweep[gam]:>8.3f} "
          f"{100 * (20.2953 - tz_sweep[gam]) / 20.2953:>10.2f}")

print()
print("=" * 78)
print("STAGE E: the first-order dressed clock, closed in gamma-free objects")
print("=" * 78)
dtE = 0.005
tsE = np.arange(0.0, 22.0 + dtE / 2, dtE)
_, C_E, nfree_E, Ifree_E = kernels(N, DELTA, tsE)
kz = np.where((Ifree_E[:-1] / Ifree_E[0] > 0) & (Ifree_E[1:] / Ifree_E[0] <= 0))[0][0]
T0 = tsE[kz] + dtE * Ifree_E[kz] / (Ifree_E[kz] - Ifree_E[kz + 1])
fV = np.zeros(len(tsE))
for b in range(N):
    fV += fftconvolve(C_E[:, b], nfree_E[:, b])[: len(tsE)]
kk = int(round(T0 / dtE))
V0 = dtE * (fV[kk] - 0.5 * (C_E[kk] @ nfree_E[0]) - 0.5 * (C_E[0] @ nfree_E[kk]))
dI0 = (Ifree_E[kk + 1] - Ifree_E[kk - 1]) / (2 * dtE)
slope = -V0 / dI0                                     # dT_zero / dGamma
split_E, _ = pair_split(N, DELTA, M_SEED)
c_adv = -slope * split_E / T0
print(f"  T_0 = {T0:.4f}, V(T_0) = {V0:.4f}, I_free'(T_0) = {dI0:.5f}, "
      f"dT_zero/dGamma = {slope:.3f}, c = {c_adv:.4f}")
# the LAW gate: the O(Gamma^2) residual of the prediction must fall ~4x when Gamma halves
res = {}
for gam in [0.002, 0.001]:
    if gam not in tz_sweep:
        ts, I_lad, _, _, _ = ladder(N, DELTA, gam, 0.02, T_MAX)
        tz_sweep[gam] = t_zero(ts, I_lad)
    res[gam] = abs(tz_sweep[gam] - (T0 + slope * 4.0 * gam))
gate("E1 first-order law: residual is O(Gamma^2) (Gamma-halving ratio ~ 4)",
     res[0.002] < 0.02 and 2.5 < res[0.002] / res[0.001] < 6.5,
     f"res(0.002) = {res[0.002]:.2e}, res(0.001) = {res[0.001]:.2e}, "
     f"ratio = {res[0.002] / res[0.001]:.2f}")
gate("E2 the coefficient pin: dT_zero/dGamma = -25.08 (N=8, m=1, delta=0.15)",
     abs(slope + 25.08) < 0.05, f"slope = {slope:.3f}")

print()
if FAIL:
    print(f"{len(FAIL)} GATE(S) FAILED:", *FAIL, sep="\n  ")
    raise SystemExit(1)
print("ALL GATES PASS.")
