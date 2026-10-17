"""
The shared field in a crystal: the coherence-variance law against multiple-quantum NMR data
===========================================================================================
experiments/THE_SHARED_FIELD_IN_A_CRYSTAL.md. Under Z noise with Kossakowski (covariance) matrix a
over the sites, a coherence |i><j| collects the quadratic form

    q(i, j) = 2 v^T a v,   v_k = (z_k(i) - z_k(j)) / 2  in {0, +1, -1},

which for a field shared by every pair with correlation eta reads 2 gamma [(1 - eta) h + eta p^2]
(h = Hamming distance, p = Delta popcount = the NMR coherence order; PROOF_ABSORPTION_THEOREM §2).
Read Markovian (Lindblad), q is a decay RATE; read quasi-static (a frozen Gaussian field), q is a
phase VARIANCE and the decay is a Gaussian of width 1/sqrt(q).

Parts:
  A  the quadratic form against the full Lindblad superoperator, exactly (integer covariances);
     the mean of q over all sign patterns of a K-site disagreement is 2 tr(a_DD), exactly;
  B  Markovian reading: the 1/e rate of a K-spin cluster signal (uniform over sign patterns) grows
     about linearly in K for every covariance family tried, never near the measured K^0.43..0.48
     (Krojanski & Suter, PRA 74, 062319 (2006));
  C  quasi-static reading against Fig. 4 of Cho, Cappellaro, Cory, Ramanathan, arXiv:cond-mat/0608620
     (CaF2, Gaussian decay times Td(n, N)): Td ∝ 1/sqrt((1 - eta) c N + eta n^2) with h = c N; only
     kappa = eta / ((1 - eta) c) is identified, ONE kappa for all points, compared with the authors'
     linear law (their constants, and their form refit) and the Fedorov-Fedichkin two-term form;
  D  (only with --rebirth, ~5 min) renewal on a small chain: an 8-site chain with hopping plus local
     dephasing gives a log-log exponent that moves with gamma; the rates are an offset plus a slope that
     approaches 2 gamma per added spin.

The Fig. 4 points were read off a screenshot of the published figure by two readers independently
(DIGITIZED_A by eye from pixel positions, DIGITIZED_B by a separate agent); part C runs on both.
Table 1 of the paper gives the cluster sizes N per MQ growth time.

Import-inert; prints only, and exits with status 1 if any check fails.
"""

import itertools
import sys

import numpy as np
from scipy.optimize import brentq

CLAIMS = []


def claim(label, ok, detail=""):
    CLAIMS.append(bool(ok))
    print(f"  [{'ok' if ok else 'FAIL'}] {label} {detail}")


# ---------------------------------------------------------------- part A
Z1 = np.diag([1, -1])
I1 = np.eye(2, dtype=int)


def zop(N, k):
    out = np.array([[1]])
    for l in range(N):
        out = np.kron(out, Z1 if l == k else I1)
    return out


def lindblad_diag_times_two(a):
    """2·L for jumps on Z with Kossakowski matrix a (integer), row-major vec; returns -diag as d×d.

    D(rho) = sum_kl a_kl (Z_k rho Z_l - 1/2 {Z_l Z_k, rho}); multiplied by 2 to stay in integers."""
    N = a.shape[0]
    d = 2 ** N
    Zs = [zop(N, k) for k in range(N)]
    Id = np.eye(d, dtype=int)
    L2 = np.zeros((d * d, d * d), dtype=np.int64)
    for k in range(N):
        for l in range(N):
            if a[k, l] == 0:
                continue
            ZlZk = Zs[l] @ Zs[k]
            L2 += a[k, l] * (2 * np.kron(Zs[k], Zs[l].T) - np.kron(ZlZk, Id) - np.kron(Id, ZlZk.T))
    offdiag = L2 - np.diag(np.diag(L2))
    return -np.diag(L2).reshape(d, d), int(np.abs(offdiag).max())


def quadratic_form(a):
    N = a.shape[0]
    d = 2 ** N
    z = np.array([[1 - 2 * ((i >> (N - 1 - k)) & 1) for k in range(N)] for i in range(d)])
    q = np.zeros((d, d), dtype=np.int64)
    for i in range(d):
        for j in range(d):
            v = (z[i] - z[j]) // 2
            q[i, j] = 2 * v @ a @ v
    return q


def part_a():
    print("A. the quadratic form against the full Lindblad generator (exact integers)")
    rng = np.random.default_rng(11)
    for N in (3, 4, 5):
        for name, a in [("local 3·I", 3 * np.eye(N, dtype=int)),
                        ("random B·Bᵀ", (lambda B: B @ B.T)(rng.integers(-2, 3, size=(N, N))))]:
            rate2, off = lindblad_diag_times_two(a)
            q = quadratic_form(a)
            claim(f"N={N} {name}: generator diagonal on every |i><j|, rate = 2 vᵀav exactly",
                  off == 0 and np.array_equal(rate2, 2 * q),
                  f"(max offdiag {off}, max |2L - 2q| {int(np.abs(rate2 - 2 * q).max())})")
    for K in (3, 5, 7):
        B = rng.integers(-3, 4, size=(K, K))
        a = B @ B.T
        V = np.array(list(itertools.product((-1, 1), repeat=K)))
        mean_q_times_count = int((2 * np.einsum("pk,kl,pl->p", V, a, V)).sum())
        claim(f"K={K}: mean of 2vᵀav over all 2^K sign patterns = 2·tr(a), exactly",
              mean_q_times_count == 2 * int(np.trace(a)) * len(V),
              f"({mean_q_times_count / len(V):.0f} vs {2 * np.trace(a)})")


# ---------------------------------------------------------------- part B
def cluster_rate(a, ts):
    K = a.shape[0]
    V = np.array(list(itertools.product((-1.0, 1.0), repeat=K)))
    r = 2 * np.einsum("pk,kl,pl->p", V, a, V)
    S = np.array([np.exp(-r * t).mean() for t in ts])
    k = int(np.argmax(S < np.exp(-1)))
    t0, t1, s0, s1 = ts[k - 1], ts[k], S[k - 1], S[k]
    return 1 / (t0 + (np.log(s0) + 1) / (np.log(s0) - np.log(s1)) * (t1 - t0))


def part_b():
    print("\nB. Markovian reading: 1/e rate of a K-spin cluster signal vs K (gamma = 1, K = 7..16)")
    dist = lambda K: np.abs(np.subtract.outer(range(K), range(K)))
    models = {
        "local (eta = 0)": lambda K: np.eye(K),
        "uniform eta = 0.5": lambda K: 0.5 * np.eye(K) + 0.5 * np.ones((K, K)),
        "fully shared (eta = 1)": lambda K: np.ones((K, K)),
        "exponential rho = 0.9": lambda K: 0.9 ** dist(K),
        "alternating rho = -0.9": lambda K: (-0.9) ** dist(K),
        "power law (1+d)^-1/2": lambda K: 1 / (1 + dist(K)) ** 0.5,
    }
    Ks = np.arange(7, 17)
    ts = np.linspace(0, 2.0, 40001)
    exps, jensen, ratio16 = {}, True, {}
    for name, mk in models.items():
        rates = [cluster_rate(mk(K), ts / K) for K in Ks]
        exps[name] = np.polyfit(np.log(Ks), np.log(rates), 1)[0]
        # Jensen: S(t) = E exp(-r t) >= exp(-E[r] t), so the 1/e rate never exceeds the mean 2·tr(a) = 2K
        jensen &= all(r <= 2 * K * (1 + 1e-9) for r, K in zip(rates, Ks))
        ratios = [r / (2 * K) for r, K in zip(rates, Ks)]
        ratio16[name] = ratios[-1]
        print(f"    {name:24s} exponent {exps[name]:.3f}   1/e rate over 2K: {ratios[0]:.3f} at K=7 .. "
              f"{ratios[-1]:.3f} at K=16 (range {min(ratios):.3f}..{max(ratios):.3f})")
    print("    (the alternating family is the exponential one under the gauge v_k -> (-1)^k v_k: five families, not six)")
    claim("in the five families tried the exponent stays at or above K^0.9 (measured: K^0.43..0.48)",
          min(exps.values()) > 0.9, f"(smallest {min(exps.values()):.3f})")
    claim("correlations only lower the 1/e rate below the mean 2K (Jensen), the fully shared field to under half",
          jensen and ratio16["fully shared (eta = 1)"] < 0.5, f"(fully shared: {ratio16['fully shared (eta = 1)']:.2f} of 2K)")
    claim("the local reading is exactly linear", abs(exps["local (eta = 0)"] - 1) < 1e-6,
          f"({exps['local (eta = 0)']:.6f}; the local rates are 2K exactly, the exponent is read off a 1/e interpolation)")


# ---------------------------------------------------------------- part C
# keys: cluster size N from Table 1 (MQ growth times 43.3, 86.8, 130.2, 173.6, 217, 260.4, 303.8 us)
DIGITIZED_A = {  # Td in microseconds, by eye from pixel positions, ~±0.05-0.1 us
    6: {2: 15.72},
    8: {2: 13.50, 4: 10.57},
    12: {2: 10.53, 4: 9.05, 6: 7.99, 8: 6.64},
    28: {2: 7.71, 4: 6.97, 6: 6.43, 8: 5.73, 10: 5.19, 12: 4.50},
    36: {2: 6.00, 4: 5.60, 6: 5.28, 8: 4.93, 10: 4.60, 12: 4.32, 14: 4.02},
    66: {2: 4.58, 4: 4.41, 6: 4.28, 8: 4.08, 10: 3.82, 12: 3.58, 14: 3.48, 16: 3.30, 18: 2.98, 20: 2.65,
         22: 2.32, 24: 2.11},
    96: {2: 3.63, 4: 3.58, 6: 3.48, 8: 3.30, 10: 3.22, 12: 3.07, 14: 3.02, 16: 2.80, 18: 2.67},
}
DIGITIZED_B = {  # independent reading: pixel analysis (PIL), axis fit residual < 0.003 us, ±0.02 us
    # isolated markers; flagged: N=28 n=12 and N=36 n=12 merged/partly hidden (±0.03-0.04),
    # N=36 n=16 mostly hidden (±0.05), N=66 n=24 cut by the axis, N=96 n=20 hidden, from error-bar caps only
    6: {2: 15.72},
    8: {2: 13.51, 4: 10.57},
    12: {2: 10.51, 4: 9.02, 6: 7.96, 8: 6.60},
    28: {2: 7.72, 4: 6.98, 6: 6.43, 8: 5.74, 10: 5.19, 12: 4.50},
    36: {2: 6.01, 4: 5.62, 6: 5.30, 8: 4.93, 10: 4.62, 12: 4.33, 14: 4.02, 16: 3.43},
    66: {2: 4.59, 4: 4.44, 6: 4.29, 8: 4.10, 10: 3.83, 12: 3.58, 14: 3.48, 16: 3.30, 18: 2.99, 20: 2.67,
         22: 2.33, 24: 2.10},
    96: {2: 3.62, 4: 3.59, 6: 3.49, 8: 3.28, 10: 3.22, 12: 3.07, 14: 3.00, 16: 2.84, 18: 2.67, 20: 2.63,
         22: 2.53, 24: 2.28},
}

FLAGGED_B = {(28, 12), (36, 12), (36, 16), (66, 24), (96, 20)}  # merged, hidden or cut by the axis


def flatten(data):
    N = np.array([k for k, v in data.items() for _ in v], float)
    n = np.array([m for v in data.values() for m in v], float)
    T = np.array([t for v in data.values() for t in v.values()])
    return N, n, T


def best_scale_fit(model, grid, T):
    best = None
    for p in grid:
        g = model(p)
        s = np.exp(np.mean(np.log(T) - np.log(g)))
        rel = s * g / T - 1
        rms = float(np.sqrt(np.mean(rel ** 2)))
        if best is None or rms < best[0]:
            best = (rms, p, rel)
    return best


def part_c(data, label):
    print(f"\nC. quasi-static reading against Cho et al. Fig. 4, {label}")
    N, n, T = flatten(data)
    collapse = {k: round(float(v[min(v)] * np.sqrt(k)), 1) for k, v in data.items()}
    print(f"    Td·sqrt(N) at the lowest order of each N: {collapse}")
    # the law sees (1 - eta)·c·N + eta·n², c the unknown constant in h = c·N: only kappa = eta / ((1 - eta)·c)
    # is identified, the shared variance per unit order² over the local variance per correlated spin
    kappa_model = lambda k, NN, nn: 1 / np.sqrt(NN + k * nn ** 2)
    kgrid = np.linspace(0.0, 2.0, 4001)
    h2 = best_scale_fit(lambda k: kappa_model(k, N, n), kgrid, T)
    e1 = best_scale_fit(lambda eta: 1 / np.sqrt((1 - eta) * N + eta * n ** 2), np.linspace(0.001, 0.999, 999), T)
    e2 = best_scale_fit(lambda eta: 1 / np.sqrt((1 - eta) * N / 2 + eta * n ** 2), np.linspace(0.001, 0.999, 999), T)
    lin = best_scale_fit(lambda r: np.clip(1 / np.sqrt(N) - r * n / N, 1e-9, None), np.linspace(0, 1, 1001), T)
    # the authors' own fit: ln Td(0) = 3.94 - 0.57 ln N, ln(-kappa) = 2.75 - 1.21 ln N, Td = Td(0) - kappa_C·n,
    # read with Td decreasing in n (their Eq. 10); A = e^3.94 = 51.4, B = e^2.75 = 15.6
    pub = (51.4 * N ** -0.57 - 15.6 * n * N ** -1.21) / T - 1
    pub_rms = float(np.sqrt(np.mean(pub ** 2)))
    pub_round = (51.4 / np.sqrt(N) - 15.6 * n / N) / T - 1
    pub_round_rms = float(np.sqrt(np.mean(pub_round ** 2)))

    def ff(p):
        out = []
        for Ni, ni in zip(N, n):
            f = lambda t: p * np.exp(-ni ** 2 * t ** 2) + (1 - p) * np.exp(-Ni / 2 * t ** 2) - np.exp(-0.5)
            out.append(brentq(f, 1e-9, 50.0))
        return np.array(out)

    fed = best_scale_fit(ff, np.linspace(0.0, 0.6, 301), T)
    print(f"    one kappa for all {len(T)} points: kappa = {h2[1]:.3f}, rms {100 * h2[0]:.2f} %")
    print(f"      the same fit read as eta: {e1[1]:.3f} if h = N (rms {100 * e1[0]:.2f} %), "
          f"{e2[1]:.3f} if h = N/2 (rms {100 * e2[0]:.2f} %): eta itself is not identified")
    print(f"    the authors' own fit (A = 51.4, B = 15.6 us, exponents -0.57, -1.21; four constants fitted by them): "
          f"rms {100 * pub_rms:.2f} %  (exponents rounded to -1/2, -1: {100 * pub_round_rms:.1f} %)")
    print(f"    the authors' functional form A/sqrt(N) - B n/N refit (2 params): rms {100 * lin[0]:.2f} %")
    print(f"    Fedorov-Fedichkin two-term form, Td read at S = e^(-1/2) (2 params): rms {100 * fed[0]:.2f} %")
    kap, within, signs = {}, {}, {}
    for k, v in data.items():
        if len(v) < 3:
            continue
        m = N == k
        b = best_scale_fit(lambda kk: kappa_model(kk, k, n[m]), kgrid, T[m])
        kap[k], within[k] = b[1], b[0]
    allres = h2[2]
    for k in data:
        m = N == k
        signs[k] = "".join("+" if r > 0 else "-" for r in allres[m])
    print("    kappa per cluster size: " + ", ".join(f"N={k} {e:.3f} ({100 * within[k]:.1f} %)" for k, e in kap.items()))
    print("    sign of model/data - 1 per N (one-kappa fit): " + ", ".join(f"N={k} {s}" for k, s in signs.items()))
    # the pre-registered shape criterion: at fixed N, 1/Td² must be affine in n²; against it, a straight
    # line Td in n (the authors' form at fixed N), both with 2 parameters per size, unweighted (the
    # figure's error bars are not digitized)
    shape = {}
    for k, v in data.items():
        if len(v) < 4:
            continue
        m = N == k
        lin_k = np.polyfit(n[m], T[m], 1)
        lin_rms = float(np.sqrt(np.mean((np.polyval(lin_k, n[m]) / T[m] - 1) ** 2)))
        shape[k] = (within[k], lin_rms)
    print("    per-size rms, one-kappa law vs straight line in n: " +
          ", ".join(f"N={k} {100 * a:.1f} vs {100 * b:.1f} %" for k, (a, b) in shape.items()))
    Nk = np.array(list(kap)); Kv = np.array(list(kap.values()))
    trend = np.polyfit(np.log(Nk), np.log(Kv), 1)[0]
    claim("one kappa fits all points within 6 % rms", h2[0] < 0.06, f"({100 * h2[0]:.2f} %)")
    claim("kappa does not drift with cluster size (|d ln kappa / d ln N| < 0.15 over N = 12..96)",
          abs(trend) < 0.15, f"(slope {trend:+.3f}, values {Kv.min():.2f}..{Kv.max():.2f})")
    claim("eta is not identified: h = N and h = N/2 give the same rms (to the grid)",
          abs(e1[0] - e2[0]) < 2e-4, f"({100 * e1[0]:.3f} vs {100 * e2[0]:.3f} %)")
    claim("the one-kappa law beats the two-term form", h2[0] < fed[0], f"({100 * h2[0]:.2f} vs {100 * fed[0]:.2f} %)")
    claim("recorded negative: the one-kappa law is slightly WORSE than the authors' form refit and their own fit",
          h2[0] > lin[0] and h2[0] > pub_rms, f"({100 * h2[0]:.2f} vs {100 * lin[0]:.2f} and {100 * pub_rms:.2f} %)")
    wins = [k for k, (a, b) in shape.items() if b < a - 0.002]
    claim("recorded negative: unweighted, a straight line in n fits better than the law at four or more of the "
          "five sizes with >= 4 points (margin 0.2 % set after the fact)", len(wins) >= 4,
          f"(line better by > 0.2 % at N = {wins})")
    return h2[1]


# ---------------------------------------------------------------- part D
def part_d():
    print("\nD. the rebirth toy: N = 8 chain, XX+YY hopping (J = 1) plus local Z dephasing, P = X on sites 0..K-1")
    N = 8
    d = 2 ** N
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]])
    I2 = np.eye(2)

    def op(ops):
        out = np.array([[1.0 + 0j]])
        for l in range(N):
            out = np.kron(out, ops.get(l, I2))
        return out

    Zm = np.diag([1.0 + 0j, -1.0])
    H_xy = sum(op({l: X, l + 1: X}) + op({l: Y, l + 1: Y}) for l in range(N - 1))
    H_dip = sum(op({l: Zm, l + 1: Zm}) - 0.5 * (op({l: X, l + 1: X}) + op({l: Y, l + 1: Y})) for l in range(N - 1))
    H = H_xy
    idx = np.arange(d)
    ham = np.array([[bin(i ^ j).count("1") for j in idx] for i in idx])

    def rate(P, g, dt=0.004, tmax=40.0, Hm=None):
        Hm = H if Hm is None else Hm
        mask = -2 * g * ham
        f = lambda r: -1j * (Hm @ r - r @ Hm) + mask * r
        rho, t, prev = P.copy(), 0.0, (0.0, 1.0)
        nrm = np.real(np.trace(P @ P))
        while t < tmax:
            k1 = f(rho); k2 = f(rho + dt / 2 * k1); k3 = f(rho + dt / 2 * k2); k4 = f(rho + dt * k3)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
            F = np.real(np.trace(P @ rho)) / nrm
            if F < np.exp(-1):
                return 1 / (prev[0] + (prev[1] - np.exp(-1)) / (prev[1] - F) * dt)
            prev = (t, F)
        return np.nan

    exps, slopes = {}, {}
    K = np.arange(1, 8)
    for g in (0.05, 0.2, 0.5, 1.0):
        r = np.array([rate(op({l: X for l in range(k)}), g) for k in K])
        exps[g] = np.polyfit(np.log(K[1:]), np.log(r[1:]), 1)[0]
        slopes[g] = np.polyfit(K[3:], r[3:], 1)[0]
        print(f"    gamma {g}: rates " + " ".join(f"{x:.3f}" for x in r) +
              f" | log-log exponent {exps[g]:.3f} | slope K=4..7 {slopes[g]:.3f} (2 gamma = {2 * g})")
    r = np.array([rate(op({l: X for l in range(k)}), 0.2, Hm=H_dip) for k in K])
    print("    dipolar-like ZZ - (XX+YY)/2, gamma 0.2: rates " + " ".join(f"{x:.3f}" for x in r) +
          f" | log-log exponent {np.polyfit(np.log(K[1:]), np.log(r[1:]), 1)[0]:.3f}")
    claim("the log-log exponent moves with gamma (an offset artifact, not a scaling law)",
          max(exps.values()) - min(exps.values()) > 0.5, f"({min(exps.values()):.2f}..{max(exps.values()):.2f})")
    claim("the slope per added spin approaches 2 gamma (within 20 % over K = 4..7 at every gamma)",
          all(abs(slopes[g] / (2 * g) - 1) < 0.2 for g in slopes),
          "(" + ", ".join(f"{slopes[g] / (2 * g):.2f}" for g in slopes) + " of 2 gamma)")


def main():
    part_a()
    part_b()
    part_c(DIGITIZED_A, "reading A (by eye)")
    part_c(DIGITIZED_B, "reading B (independent)")
    reduced = {k: {m: t for m, t in v.items() if (k, m) not in FLAGGED_B} for k, v in DIGITIZED_B.items()}
    part_c(reduced, "reading B without its five flagged points")
    if "--rebirth" in sys.argv:
        part_d()
    held = sum(CLAIMS)
    print(f"\n{held} of {len(CLAIMS)} claims hold.")
    if held < len(CLAIMS):
        sys.exit(1)


if __name__ == "__main__":
    main()
