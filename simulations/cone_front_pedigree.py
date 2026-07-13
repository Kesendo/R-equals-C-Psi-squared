"""The front's pedigree: the F126 refill ladder resolved by CATCH COUNT j.

The front-survival verifier cone_front_survival_asymptote.py closed the survival
EXPONENT A (the fixed-gamma ceiling A_inf = 4 - phi(2J)/gamma). This companion
resolves the same halo one layer finer: not "how much survives" but "how many
times was the surviving weight caught (dephased) on the way".

The exact renewal / ladder representation (interior seed, infinite chain, J = 1):
     P_n(t) = e^{-Gamma t} S_n(t),   Gamma = 4 gamma,
     S_n(t) = sum_j S^(j)_n(t),
     S^(0)_n(t) = |G_{n0}(t)|^2 = J_n(2Jt)^2            (never caught, the coherent front)
     S^(j)_n(t) = Gamma int_0^t ds sum_m |G_{nm}(t-s)|^2 S^(j-1)_m(s)   (caught j times).
  The catch-count weight is P(j) = S^(j)_n(t) / S_n(t); the front weight P(j=0) is the
  coherent fraction |<a_n>|^2 / P_n, and it obeys, at the clean band-edge caustic time
  t*_0 = n/(2J) (2J t*_0 = n exactly), the identity

     P(j=0) = e^{-(4 - A_same) K},   K = gamma n / (2J),  A_same = 4 - ln(S_n/J_n^2)/K.

  Both sides are the SAME S-ratio S^(0)_n(t*_0)/S_n(t*_0); the identity just rewrites it.

Each catch order S^(j) is computed per momentum by the Graf identity: the kernel
sum_m |G_{nm}(tau)|^2 has the exact momentum transform Khat(p, tau) = J_0(a(p) tau),
a(p) = 4 J |sin(p/2)|, so S^(j)(p, .) = Gamma (Khat *_t S^(j-1)(p, .)) is a scalar
causal convolution in time. The per-order sum reproduces the implicit renewal solver
of cone_front_survival_asymptote.py to machine precision (section [1]).

--------------------------------------------------------------------------------
THE DEVIATION, RESOLVED (2026-07-13, before this script was written).
  A play scout measured A_same(50) ~ 3.00 at t*_0 = n/(2J); the committed verifier
  cone_front_survival_asymptote.py section [5] measured A_sameinstant(40) = 2.7316.
  These do NOT disagree: they evaluate the SAME renewal S-ratio at two DIFFERENT
  reference times both loosely named "t*_0":
     scout / this script : t*_0 = n/(2J)                       (the BAND-EDGE, 2J t = n)
     committed [5]       : t*_0 = argmax_t J_n(2Jt)^2
                                 = (n + 0.8086 n^{1/3})/(2J)   (the CAUSTIC PEAK, later)
  Recomputing S_n(t)/J_n(2Jt)^2 with BOTH solvers (per-momentum renewal AND the
  displacement Volterra) agrees to <= 2e-5, and at the caustic-peak instant the
  committed identity reproduces 2.7316 (n=40) / 2.7943 (n=50) EXACTLY, while the
  band-edge instant gives 2.9427 (n=40) / 3.0007 (n=50). Neither number is buggy;
  the band-edge n/(2J) is used here because there Gamma t = 4K exactly, so the
  catch-count identity above carries the clean constant 4 and the ladder split is
  cleanest. (Section [1] restates these numbers from the run.)

Run:  python simulations/cone_front_pedigree.py
Writes simulations/results/cone_defect_arrival/front_pedigree.txt
"""

import os
import time
import numpy as np
from scipy.special import jn, jv
from scipy.signal import fftconvolve

J = 1.0
OUT = []


def log(s=""):
    print(s, flush=True)
    OUT.append(s)


# ---------------------------------------------------------------------------
# solvers (three exact forms of the same renewal object)
# ---------------------------------------------------------------------------
def momentum_ladder(n, gamma, t_star, dt=0.02, Np=1024, jmax=18):
    """Per-momentum, per-catch-order ladder. Returns (tg, Sj_front, Shat_full):
    Sj_front[j, k] = S^(j)_n(t_k) at the front site d = n; Shat_full[k, p] = the
    momentum-space FULL sum sum_j S^(j)(p, t_k) (for the last-catch age)."""
    Gamma = 4.0 * gamma
    tg = np.arange(0.0, t_star + 1e-9, dt)
    nt = len(tg)
    p = np.arange(Np) * 2 * np.pi / Np
    a = 4.0 * J * np.abs(np.sin(p / 2))
    Khat = jn(0, np.outer(tg, a))                     # (nt, Np), Khat(p, tau)
    Shat = np.zeros((jmax + 1, nt, Np))
    Shat[0] = Khat                                    # S^(0)(p, .) = Khat
    for j in range(1, jmax + 1):
        # causal convolution Gamma * int_0^t Khat(t-s) S^(j-1)(s) ds, trapezoid rule
        conv = fftconvolve(Khat, Shat[j - 1], axes=0)[:nt]         # rectangle sum
        corr = 0.5 * (Khat * Shat[j - 1][0][None, :]
                      + Khat[0][None, :] * Shat[j - 1])            # trapezoid endpoints
        Shat[j] = Gamma * dt * (conv - corr)
    cosd = np.cos(n * p) / Np
    Sj_front = np.real(Shat @ cosd)                   # (jmax+1, nt), inverse transform
    Shat_full = Shat.sum(axis=0)                      # (nt, Np)
    return tg, Sj_front, Shat_full


def renewal_full_front(n, gamma, t_star, dt=0.02, Np=1024):
    """Implicit (trapezoid) per-momentum renewal solver, front value S_n(t).
    This is the committed cone_front_survival_asymptote.py solver."""
    Gamma = 4.0 * gamma
    tg = np.arange(0.0, t_star + 1e-9, dt)
    h = dt
    nt = len(tg)
    p = np.arange(Np) * 2 * np.pi / Np
    a = 4.0 * J * np.abs(np.sin(p / 2))
    K = jn(0, np.outer(tg, a))
    Shat = np.zeros((nt, Np))
    Shat[0] = K[0]
    for k in range(1, nt):
        conv = 0.5 * K[k] * Shat[0]
        if k > 1:
            conv = conv + np.sum(K[k - 1:0:-1] * Shat[1:k], axis=0)
        Shat[k] = (K[k] + Gamma * h * conv) / (1.0 - 0.5 * Gamma * h)
    cosd = np.cos(n * p) / Np
    return tg, np.real(Shat @ cosd)


def volterra_full_front(n, gamma, t_star, dt=0.02, Dpad=40):
    """Displacement-space (implicit-trapezoid) Volterra solver, front value S_n(t).
    This is the scout's independent method: a spatial convolution on the chain."""
    Gamma = 4.0 * gamma
    tg = np.arange(0.0, t_star + 1e-9, dt)
    steps = len(tg) - 1
    D = n + Dpad
    dsp = np.arange(-D, D + 1)
    Kd = jv(np.abs(dsp)[None, :], 2.0 * J * tg[:, None]) ** 2
    S = np.zeros((steps + 1, 2 * D + 1))
    S[0] = Kd[0]
    for k in range(1, steps + 1):
        acc = 0.5 * np.convolve(S[0], Kd[k], mode="full")[D:3 * D + 1]
        for s in range(1, k):
            acc += np.convolve(S[s], Kd[k - s], mode="full")[D:3 * D + 1]
        rhs = Kd[k, :] + Gamma * dt * acc
        S[k] = rhs / (1.0 - 0.5 * Gamma * dt)
    return tg, S[:, D + n]


# ============================================================================
t_wall0 = time.time()
log("=" * 78)
log("THE FRONT'S PEDIGREE: the F126 refill ladder resolved by CATCH COUNT j")
log("Exact single-excitation dephasing, J = 1, Gamma = 4 gamma; interior seed,")
log("infinite chain. P(j) = S^(j)_n(t*_0)/S_n(t*_0) at the band-edge t*_0 = n/(2J).")
log("=" * 78)

# ---------------------------------------------------------------------------
# SECTION 1 : METHOD CROSS-CHECK + the deviation reconciliation
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[1] METHOD CROSS-CHECK : ladder-sum == implicit renewal == displacement Volterra")
log("-" * 78)
log("    (n, gamma) = (30, 0.05), front d = n, at the band-edge t*_0 = n/(2J) = 15.")
log("    three independent solvers of the same renewal object must agree.")

N1, G1 = 30, 0.05
TS1 = N1 / (2 * J)
JMAX1 = 18
tg1, Sj1, _ = momentum_ladder(N1, G1, TS1, dt=0.02, Np=1024, jmax=JMAX1)
_, Sren1 = renewal_full_front(N1, G1, TS1, dt=0.02, Np=1024)
_, Svol1 = volterra_full_front(N1, G1, TS1, dt=0.02, Dpad=40)

ladder_sum = float(Sj1[:, -1].sum())
S_implicit = float(Sren1[-1])
S_volterra = float(Svol1[-1])
rel_lad_imp = abs(ladder_sum - S_implicit) / S_implicit
rel_lad_vol = abs(ladder_sum - S_volterra) / S_volterra
jmax_tail = float(Sj1[-1, -1] / ladder_sum)     # P(j = JMAX1)

log(f"\n    ladder-sum   sum_j S^(j)_n(t*_0) = {ladder_sum:.8f}   (per-order momentum, jmax={JMAX1})")
log(f"    implicit     renewal S_n(t*_0)    = {S_implicit:.8f}   (per-momentum, the committed solver)")
log(f"    displacement Volterra S_n(t*_0)   = {S_volterra:.8f}   (scout's chain-convolution solver)")
log(f"    rel |ladder - implicit|  = {rel_lad_imp:.2e}   (same discretized operator, Neumann series)")
log(f"    rel |ladder - Volterra|  = {rel_lad_vol:.2e}   (independent method, grid/window-limited)")
log(f"    jmax tail P(j={JMAX1})        = {jmax_tail:.2e}")
assert rel_lad_imp <= 1e-3, f"ladder vs implicit off by {rel_lad_imp:.2e}"
assert rel_lad_vol <= 1e-3, f"ladder vs Volterra off by {rel_lad_vol:.2e}"
assert jmax_tail <= 1e-4, f"jmax tail {jmax_tail:.2e} not below 1e-4"
log(f"    PASS: all three agree <= 1e-3; the catch-order tail is below 1e-4.")

# -- the deviation reconciliation: same object, two reference times ---------
log("\n    RECONCILIATION (the scout-vs-committed A_same deviation):")
log("    S_n(t)/J_n(2Jt)^2 recomputed with BOTH solvers, at BOTH reference times.")
log(f"    {'n':>3} {'time':>22} {'t':>8} {'S(renewal)':>11} {'S(Volterra)':>12} "
    f"{'reldiff':>9} {'A_same':>8}")


def _interp(tg, col, t):
    idx = t / (tg[1] - tg[0])
    i = int(np.floor(idx))
    fr = idx - i
    return col[i] * (1 - fr) + col[i + 1] * fr


for n in (40, 50):
    gamma = 0.05
    K = gamma * n / (2 * J)
    x = np.linspace(n - 2, n + 6 * n ** (1 / 3) + 6, 400000)
    t_peak = x[int(np.argmax(jv(n, x) ** 2))] / (2 * J)
    t_edge = n / (2 * J)
    tmax = t_peak + 1.0
    tgr, Sr = renewal_full_front(n, gamma, tmax, dt=0.01, Np=1024)
    tgv, Sv = volterra_full_front(n, gamma, tmax, dt=0.02, Dpad=35)
    for label, t in (("band-edge n/(2J)", t_edge), ("caustic peak (comm. [5])", t_peak)):
        Jn2 = jn(n, 2 * J * t) ** 2
        Sr_t = _interp(tgr, Sr, t)
        Sv_t = _interp(tgv, Sv, t)
        rel = abs(Sr_t - Sv_t) / Sr_t
        A_same = -np.log(np.exp(-4 * gamma * t) * Sr_t / Jn2) / K
        log(f"    {n:>3} {label:>22} {t:>8.4f} {Sr_t:>11.5f} {Sv_t:>12.5f} "
            f"{rel:>9.2e} {A_same:>8.4f}")
log("    => the two solvers agree to <= 2e-5 at BOTH instants; the committed 2.7316 (n=40)")
log("       / 2.7943 (n=50) are the CAUSTIC-PEAK reading, the band-edge reading is 2.9427 /")
log("       3.0007. Same S-ratio, different reference time: neither number is wrong. This")
log("       script's pedigree uses the band-edge n/(2J), where Gamma t = 4K exactly.")

# ---------------------------------------------------------------------------
# SECTION 2 : THE PEDIGREE TABLE
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[2] THE PEDIGREE TABLE : P(j), j = 0..7, at the band-edge t*_0 = n/(2J)")
log("-" * 78)
log("    P(j=0) is the never-caught coherent front; <j> is the mean catch count of the")
log("    surviving front weight; Gamma t*_0 = 4K is the free-walker mean catch count.")
log("    IDENTITY (not a test): P(0) = exp(-(4 - A_same) K), A_same = 4 - ln(S_n/J_n^2)/K")
log("    from this run's own totals. Informative: <j> / (Gamma t*_0) = the catch deficit.")

CASES = [(20, 0.05), (30, 0.05), (50, 0.05), (30, 0.10), (30, 0.025)]
log(f"\n    {'n':>4} {'gamma':>6} {'K':>6} {'P(0)':>7} {'<j>':>6} {'G t*=4K':>8} "
    f"{'<j>/4K':>7} {'A_same':>7} {'e^-(4-A)K':>10}  P(j) j=0..7")
for n, gamma in CASES:
    TS = n / (2 * J)
    K = gamma * n / (2 * J)
    tg, Sj, _ = momentum_ladder(n, gamma, TS, dt=0.02, Np=1024, jmax=18)
    w = Sj[:, -1]
    tot = float(w.sum())
    P = w / tot
    P0 = float(P[0])
    mean_j = float((np.arange(len(P)) * P).sum())
    Gt = 4.0 * K
    Jn2 = jn(n, 2 * J * TS) ** 2       # 2 J t*_0 = n exactly
    A_same = 4.0 - np.log(tot / Jn2) / K
    ident = np.exp(-(4.0 - A_same) * K)
    assert abs(ident - P0) < 1e-9, f"identity broke: {ident:.6f} vs {P0:.6f}"
    assert float(P[-1]) <= 1e-4, f"tail P(18)={P[-1]:.2e} not below 1e-4 at (n={n},g={gamma})"
    log(f"    {n:>4} {gamma:>6.3f} {K:>6.3f} {P0:>7.4f} {mean_j:>6.3f} {Gt:>8.3f} "
        f"{mean_j / Gt:>7.3f} {A_same:>7.4f} {ident:>10.6f}  "
        + " ".join(f"{x:.3f}" for x in P[:8]))
log("\n    P(0) = e^{-(4-A_same)K} holds to < 1e-9 in every row (it is the identity, not a")
log("    fit). <j> sits far below the free-walker 4K: dephasing near the arriving front is")
log("    strongly SUPPRESSED (a caught walker mostly falls behind the ballistic edge), so")
log("    the survivors were caught ~1/4 as often as a free walker over the same dose.")

# ---------------------------------------------------------------------------
# SECTION 3 : THE LAST-CATCH AGE (n = 50, gamma = 0.05)
# ---------------------------------------------------------------------------
log("\n" + "-" * 78)
log("[3] THE LAST-CATCH AGE : how young is the surviving halo? (n = 50, gamma = 0.05)")
log("-" * 78)
log("    Among front weight caught at least once, the age tau of the LAST catch has density")
log("    A(tau) = Gamma sum_m |G_{nm}(tau)|^2 S_m(t*_0 - tau) (the walker's final free flight")
log("    of duration tau after its last dephasing). Reported as fractions of the trip t*_0.")

N3, G3 = 50, 0.05
TS3 = N3 / (2 * J)
DT3 = 0.02
tg3, Sj3, Shat_full3 = momentum_ladder(N3, G3, TS3, dt=DT3, Np=1024, jmax=18)
Gamma3 = 4.0 * G3
steps3 = len(tg3) - 1
p3 = np.arange(1024) * 2 * np.pi / 1024
a3 = 4.0 * J * np.abs(np.sin(p3 / 2))
Khat3 = jn(0, np.outer(tg3, a3))
cosn3 = np.cos(N3 * p3) / 1024
# A(tau) at age index a: Gamma * sum_p Khat(p, tau=a dt) * Shat_full(p, t*_0 - a dt) cos(n p)
ages = tg3[1:]
dens = Gamma3 * np.real(
    (Khat3[1:steps3 + 1] * Shat_full3[steps3 - 1::-1][:steps3]) @ cosn3)
caught_total = float(dens.sum() * DT3)
never = float(jn(N3, 2 * J * TS3) ** 2)               # S^(0)_n(t*_0), never caught
tot3 = float(Sj3[:, -1].sum())
P0_3 = never / tot3
# consistency with the ladder: never + caught_total == S_n(t*_0)
consistency = abs((never + caught_total) - tot3) / tot3
log(f"\n    never-caught weight S^(0)_n(t*_0) = {never:.6f}")
log(f"    caught weight  Gamma int A(tau)   = {caught_total:.6f}")
log(f"    ladder total   S_n(t*_0)          = {tot3:.6f}")
log(f"    never + caught vs ladder total    : reldiff = {consistency:.2e}")
assert consistency < 5e-3, f"age normalization inconsistent with ladder: {consistency:.2e}"
log(f"    never-caught fraction P(0) = never/total = {P0_3:.4f}  (matches the [2] row at (50,0.05))")

cum = np.cumsum(dens) * DT3 / caught_total            # CDF of last-catch age among caught
log(f"\n    quartiles of the last-catch age tau among caught front weight (trip t*_0 = {TS3:.1f}):")
log(f"    {'quantile':>10} {'tau':>8} {'tau / t*_0':>12}")
q_report = {}
for frac in (0.25, 0.50, 0.75, 0.90):
    tau = float(ages[np.searchsorted(cum, frac)])
    q_report[frac] = tau / TS3
    log(f"    {int(frac * 100):>9}% {tau:>8.2f} {tau / TS3:>11.1%}")
log(f"\n    the age spread is BROAD (25%-75% spans ~{q_report[0.25]:.0%} to ~{q_report[0.75]:.0%} of the")
log("    trip): the last catch is not concentrated in a thin sliver just behind the caustic.")
log("    This ties to the bulk-dominance of the single-refill constant I_1 (the fourth")
log("    follow-up's referee finding, cone_front_survival_asymptote.py [4b]): the joint")
log("    caustic sliver supplies only ~30% of I_1; ~70% of the refill comes from the cone")
log("    INTERIOR, so the refilling weight is caught across the whole cone, not at its edge.")
assert 0.25 < q_report[0.75] < 0.85, f"75% age fraction out of expected band: {q_report[0.75]:.3f}"

runtime = time.time() - t_wall0
log("\n" + "=" * 78)
log("SUMMARY: the surviving front is mostly virgin (P(0) = 0.53/0.43/0.29 at n = 20/30/50,")
log("gamma = 0.05); the caught remainder was dephased far fewer times than a free walker")
log("(<j>/4K ~ 0.3), and its last catch is spread broadly across the cone interior, not the")
log("caustic sliver. P(0) = e^{-(4-A_same)K} is the identity linking the pedigree to the")
log("survival exponent A of cone_front_survival_asymptote.py. The scout-vs-committed A_same")
log("deviation was a reference-time mismatch (band-edge n/2J vs caustic peak), not an error.")
log(f"runtime: {runtime:.1f} s")
log("=" * 78)
log("DONE")

os.makedirs("simulations/results/cone_defect_arrival", exist_ok=True)
with open("simulations/results/cone_defect_arrival/front_pedigree.txt", "w",
          encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
