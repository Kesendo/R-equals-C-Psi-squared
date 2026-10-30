"""The decoder and the front: one bond defect, read through two kinds of carrier
(experiments/THE_DECODER_AND_THE_FRONT.md). The second route beside the C# readings in compute/MirrorWorld
(WalkTime.Letters / Decode / DecodeOne / TrajectoryLetters, run mode `walk decode`, pinned in WalkTimeTests):
this script propagates the single-excitation block exactly (one matrix exponential per time step) where the
C# runs the Cone's RK4, so the two meet only through the physics.

The block: open XY chain, hopping J on the bonds (the Pauli convention H = J * sum (XX + YY)/2), every coherence
|a><b| (a != b) charged 4 gamma, populations p_i = rho_ii. J = 1, gamma = 0.05 (Q = 20, the DefectDecoder's
point). A letter is the derivative of a reading with respect to one bond's delta (J_b = J (1 + delta)).

  [1] the arrival letters of a site-0 seed (relative threshold 0.2 of the clean per-site peak, as in
      experiments/COUPLING_DEFECT_WALK_TIME_STEP.md), N = 3..10: worst signed cos, smallest singular value;
      the DefectDecoder's confuser at N = 5; every bond at once at N = 5 and 8, and the half-delta law.
  [2] the same observable for every carrier, the whole population trajectory (every site, t = 0.25..20):
      stationary eigenmodes, a mirror-broken mode superposition, mirror-even localized carriers, the site seed.
  [3] the decoder's own letter (the linear letter, the delta -> 0 limit of the painters' alpha fit: per site
      <g, dP>/<g, g>, g = t * dP/dt, P the purity, t on [0, 1/gamma], 60 points) for the same carriers. For
      the bonding mode it reproduces DecoderAntiCollinearitySeesawClaim's linear-letter row; for a localized
      carrier the alpha rescaling picture does not describe the data (hypotheses/PERSPECTIVAL_TIME_FIELD.md,
      "Where it degrades"), so that column reads a projection, not a fit the painters would accept.
  [2b] stationary carriers (mixtures of modes, [H, rho] = 0) against moving ones, both observables.
  [5] at equal cost: the per-shot Fisher information of the full Z-position readout for several carriers and
      schedules, two handshakes with the house's Fisher apparatus, and the efficiency of the decoder's own estimate.
  [4] random carriers at N = 7: built on m modes, or smooth site packets; the worst signed cos of [3]
      against the mode participation PR = 1 / sum_k |<psi_k|phi>|^4.
Output: simulations/results/decoder_carrier/decoder_carrier_dictionary.txt
"""
from pathlib import Path

import numpy as np
from scipy.linalg import expm

J, GAMMA = 1.0, 0.05
REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "simulations" / "results" / "decoder_carrier" / "decoder_carrier_dictionary.txt"
LINES = []


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def generator(n, deltas):
    h = np.zeros((n, n))
    for a in range(n - 1):
        h[a, a + 1] = h[a + 1, a] = J * (1.0 + deltas[a])
    eye = np.eye(n)
    return -1j * (np.kron(h, eye) - np.kron(eye, h.T)) - np.diag((4 * GAMMA * (1 - eye)).reshape(-1))


def populations(n, deltas, rho0, times):
    """p_i(t) on the given times, exactly through the generator's eigenbasis."""
    w, V = np.linalg.eig(generator(n, deltas))
    c = np.linalg.solve(V, rho0.reshape(-1).astype(complex))
    diag = [i * n + i for i in range(n)]
    E = np.exp(np.outer(times, w)) * c
    return np.real(E @ V[diag].T), np.real((E * w) @ V[diag].T)


def unit(n, k):
    v = np.zeros(n)
    v[k] = 1.0
    return v


def mode(n, k):
    v = np.sin(k * np.pi * np.arange(1, n + 1) / (n + 1))
    return v / np.linalg.norm(v)


def worst_cos(F):
    best = (2.0, -1, -1)
    for a in range(len(F)):
        for b in range(a + 1, len(F)):
            best = min(best, (float(F[a] @ F[b] / np.sqrt((F[a] @ F[a]) * (F[b] @ F[b]))), a, b))
    return best


def conditioning(F):
    """The smallest singular value of the row-normalised dictionary (0: some combination of letters vanishes, a
    direction the reading cannot see) and the largest pairwise cos (near +1: two letters nearly the same, the
    location blindness that worst_cos does not see)."""
    Fn = F / np.linalg.norm(F, axis=1, keepdims=True)
    C = Fn @ Fn.T
    iu = np.triu_indices(len(F), 1)
    return float(np.linalg.svd(Fn, compute_uv=False).min()), float(C[iu].max())


# ---- [1] the arrival letters ----
DT_ARR = 0.002


def arrival(P, times, ref):
    out = np.full(P.shape[1], np.nan)
    for i in range(1, P.shape[1]):
        thr = 0.2 * ref[i]
        idx = np.nonzero(P[:, i] >= thr)[0]
        if len(idx) and idx[0] > 0:
            k = idx[0]
            out[i] = times[k - 1] + (thr - P[k - 1, i]) / (P[k, i] - P[k - 1, i]) * (times[k] - times[k - 1])
    return out[1:]


def arrival_setup(n, theta=0.2):
    times = np.arange(0.0, 0.75 * n + 3.0, DT_ARR)
    seed = np.outer(unit(n, 0), unit(n, 0))
    P0 = populations(n, np.zeros(n - 1), seed, times)[0]
    ref = theta / 0.2 * P0.max(axis=0)
    return times, seed, ref, arrival(P0, times, ref)


def arrival_profile(n, deltas, setup):
    times, seed, ref, t0 = setup
    return arrival(populations(n, deltas, seed, times)[0], times, ref) - t0


def arrival_letters(n, setup, d=0.02):
    L = []
    for b in range(n - 1):
        e = np.zeros(n - 1)
        e[b] = d
        L.append((arrival_profile(n, e, setup) - arrival_profile(n, -e, setup)) / (2 * d))
    return np.array(L)


def part_one():
    say("[1] the arrival letters of a site-0 seed (exact propagation, dt = %.3f)" % DT_ARR)
    keep = {}
    for n in range(3, 11):
        setup = arrival_setup(n)
        L = arrival_letters(n, setup)
        keep[n] = (setup, L)
        c, a, b = worst_cos(L)
        sv = np.linalg.svd(L, compute_uv=False)
        cn, mx = conditioning(L)
        say("  N = %2d: worst signed cos %+.3f on (%d,%d); largest cos %+.3f; singular values %.3f .. %.3f, normalised smallest %.3f"
            % (n, c, a, b, mx, sv.min(), sv.max(), cn))
        if n == 5:
            for k in range(n - 1):
                say("      bond %d: %s" % (k, " ".join("%+.3f" % (0.0 if abs(x) < 5e-4 else x) for x in L[k])))
    say("  past N = 10 the worst signed cos drifts down (the bump before each bond grows with its distance from the seed):")
    for n in (12, 16, 20):
        L = arrival_letters(n, arrival_setup(n))
        c, a, b = worst_cos(L)
        cn, mx = conditioning(L)
        say("    N = %2d: worst signed cos %+.3f on (%d,%d); normalised smallest singular value %.3f" % (n, c, a, b, cn))
    say("  and it depends on the threshold theta (N = 12): " + ", ".join(
        "theta %.1f: %+.3f" % (th, worst_cos(arrival_letters(12, arrival_setup(12, th)))[0]) for th in (0.1, 0.2, 0.5)))
    setup, L = keep[5]
    for bond, delta in ((3, -0.02), (1, +0.02)):
        e = np.zeros(4)
        e[bond] = delta
        obs = arrival_profile(5, e, setup)
        res = sorted((float(np.sum((obs - (obs @ L[k] / (L[k] @ L[k])) * L[k]) ** 2)), k, float(obs @ L[k] / (L[k] @ L[k]))) for k in range(4))
        say("  the decoder's confuser at N = 5, bond %d at %+.2f: read bond %d at %+.4f, runner-up bond %d, residual ratio %.0f"
            % (bond, delta, res[0][1], res[0][2], res[1][1], res[1][0] / res[0][0]))
    truths = {5: [+0.0076, -0.0108, -0.0035, -0.0129], 8: [+0.0102, -0.0070, +0.0048, -0.0095, -0.0136, -0.0165, +0.0181]}
    for n, truth in truths.items():
        setup, L = keep[n]
        errs = []
        for scale in (1.0, 0.5):
            t = np.array(truth) * scale
            x = np.linalg.lstsq(L.T, arrival_profile(n, t, setup), rcond=None)[0]
            errs.append(float(np.abs(x - t).max()))
        say("  every bond at once, N = %d: largest error %.5f; every delta halved %.5f, ratio %.3f (the second order's 1/4)"
            % (n, errs[0], errs[1], errs[1] / errs[0]))


# ---- [2] and [3] the carriers ----
def trajectory_letters(n, rho0, d=1e-4):
    times = np.arange(1, 81) * 0.25
    L = []
    for b in range(n - 1):
        e = np.zeros(n - 1)
        e[b] = d
        L.append(((populations(n, e, rho0, times)[0] - populations(n, -e, rho0, times)[0]) / (2 * d)).reshape(-1))
    return np.array(L)


def linear_letters(n, rho0, d=1e-4):
    times = np.linspace(0.0, 1.0 / GAMMA, 60)

    def purity(deltas):
        p, pd = populations(n, deltas, rho0, times)
        return 0.5 * (1 + (1 - 2 * p) ** 2), -2 * (1 - 2 * p) * pd
    P0, Pd0 = purity(np.zeros(n - 1))
    g = times[:, None] * Pd0
    L = []
    for b in range(n - 1):
        e = np.zeros(n - 1)
        e[b] = d
        dP = (purity(e)[0] - purity(-e)[0]) / (2 * d)
        L.append((g * dP).sum(0) / (g * g).sum(0))
    return np.array(L)


def carriers(n):
    P = np.outer
    s = lambda a, b: (a + b) / np.sqrt(2)
    out = [("bonding mode psi1 (stationary, mirror-even)", P(mode(n, 1), mode(n, 1))),
           ("mode psi2 (stationary, mirror-even rho)", P(mode(n, 2), mode(n, 2))),
           ("(psi1 + psi2)/sqrt2 (moving, mirror-broken)", P(s(mode(n, 1), mode(n, 2)), s(mode(n, 1), mode(n, 2)))),
           ("end pair (|0> + |N-1>)/sqrt2 (mirror-even)", P(s(unit(n, 0), unit(n, n - 1)), s(unit(n, 0), unit(n, n - 1)))),
           ("end mixture (mirror-even)", 0.5 * (P(unit(n, 0), unit(n, 0)) + P(unit(n, n - 1), unit(n, n - 1))))]
    if n % 2 == 1:
        out.append(("centre site (mirror-even)", P(unit(n, n // 2), unit(n, n // 2))))
    out.append(("site 0 (mirror-broken)", P(unit(n, 0), unit(n, 0))))
    return out


def part_two_three():
    say("")
    say("[2] the whole population trajectory and [3] the decoder's linear letter, per carrier: worst signed cos,")
    say("    and [normalised smallest singular value, largest cos] of each dictionary")
    for n in (5, 6, 7, 8, 9, 12, 14):
        for name, rho0 in carriers(n):
            if n > 9 and not name.startswith(("bonding", "site 0")):
                continue
            Ft = trajectory_letters(n, rho0)
            ct = worst_cos(Ft)
            F = linear_letters(n, rho0)
            cl = worst_cos(F)
            mirror = ""
            if "mirror-even" in name:
                res = 0.0
                for b in range(n - 1):
                    f, g = F[b], F[n - 2 - b]
                    w = np.sum(((f - f[::-1]) / 2) ** 2) / np.sum(f * f)
                    res = max(res, abs(f @ g / np.sqrt((f @ f) * (g @ g)) - (1 - 2 * w)))
                mirror = "   cos = 1 - 2w on mirror pairs to %.0e" % res
            kt, kl = conditioning(Ft), conditioning(F)
            say("  N = %d  %-46s trajectory %+.3f (%d,%d) [smallest %.3f, largest cos %+.3f]   linear letter %+.3f (%d,%d) [smallest %.3f, largest cos %+.3f]%s"
                % (n, name, ct[0], ct[1], ct[2], kt[0], kt[1], cl[0], cl[1], cl[2], kl[0], kl[1], mirror))
    row = []
    for n in range(3, 13):
        row.append("%+.3f" % worst_cos(linear_letters(n, np.outer(mode(n, 1), mode(n, 1))))[0])
    say("  the bonding mode's linear letter, worst signed cos over all sites at N = 3..12: " + ", ".join(row))
    say("    (DecoderAntiCollinearitySeesawClaim's linear-letter row: -0.965, -0.536, -0.983, -0.973, -0.928, -0.960, -0.968, -0.931, -0.905, -0.915)")


# ---- [2b] stationary against moving carriers ----
def part_two_b():
    say("")
    say("[2b] stationary carriers ([H, rho] = 0, a mixture of modes) against moving ones: worst signed cos")
    P = np.outer
    for n in (5, 7, 8):
        h = np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
        s13 = (mode(n, 1) + mode(n, 3)) / np.sqrt(2)
        w = np.ones(n) / np.sqrt(n)
        for name, rho0 in (("psi1", P(mode(n, 1), mode(n, 1))), ("psi3", P(mode(n, 3), mode(n, 3))),
                           ("(P1 + P3)/2", (P(mode(n, 1), mode(n, 1)) + P(mode(n, 3), mode(n, 3))) / 2),
                           ("(P1 + P2)/2", (P(mode(n, 1), mode(n, 1)) + P(mode(n, 2), mode(n, 2))) / 2),
                           ("(psi1 + psi3)/sqrt2", P(s13, s13)), ("W state", P(w, w))):
            comm = np.linalg.norm(h @ rho0 - rho0 @ h) / np.linalg.norm(rho0)
            say("  N = %d  %-22s |[H, rho]|/|rho| %.2f   trajectory %+.3f   linear letter %+.3f"
                % (n, name, comm, worst_cos(trajectory_letters(n, rho0))[0], worst_cos(linear_letters(n, rho0))[0]))


# ---- [4] random carriers ----
def part_four():
    say("")
    say("[4] random carriers at N = 7, the worst signed cos of both observables against the mode participation PR")
    n = 7
    rng = np.random.default_rng(20261030)
    M = np.array([mode(n, k) for k in range(1, n + 1)])
    groups = {}

    def both(phi):
        r = np.outer(phi, phi)
        return 1 / np.sum((M @ phi) ** 4), worst_cos(linear_letters(n, r))[0], worst_cos(trajectory_letters(n, r))[0]
    for m in range(1, n + 1):
        for _ in range(12):
            c = np.zeros(n)
            c[rng.choice(n, m, replace=False)] = rng.normal(size=m)
            phi = c @ M
            phi /= np.linalg.norm(phi)
            groups.setdefault("%d mode%s" % (m, "" if m == 1 else "s"), []).append(both(phi))
    for width in (0.5, 1.0, 2.0, 4.0):
        for _ in range(12):
            x0 = rng.uniform(0, n - 1)
            phi = np.exp(-((np.arange(n) - x0) ** 2) / (4 * width ** 2))
            phi /= np.linalg.norm(phi)
            groups.setdefault("packet width %.1f" % width, []).append(both(phi))
    allv = []
    for lab, v in groups.items():
        v = np.array(v)
        allv.append(v)
        say("  %-18s PR mean %.2f   linear letter: mean %+.3f, range %+.3f .. %+.3f, below -0.5 in %2d of %d   trajectory: mean %+.3f, range %+.3f .. %+.3f, below -0.5 in %2d of %d"
            % (lab, v[:, 0].mean(), v[:, 1].mean(), v[:, 1].min(), v[:, 1].max(), int((v[:, 1] < -0.5).sum()), len(v),
               v[:, 2].mean(), v[:, 2].min(), v[:, 2].max(), int((v[:, 2] < -0.5).sum()), len(v)))
    allv = np.vstack(allv)
    say("  correlation of PR with the worst cos over all %d carriers: linear letter %+.3f, trajectory %+.3f"
        % (len(allv), np.corrcoef(allv[:, 0], allv[:, 1])[0, 1], np.corrcoef(allv[:, 0], allv[:, 2])[0, 1]))


# ---- [5] at equal cost: the Fisher information of the shots ----
def position_fisher(n, rho0, times, eps=1e-5):
    """Per-time categorical Fisher matrices of one Z shot (the excitation's position), every bond's relative delta
    a parameter: F_ab(t) = sum_j dp_j/da dp_j/db / p_j. Returns (F per time, p, dp)."""
    p = populations(n, np.zeros(n - 1), rho0, times)[0]
    dp = []
    for b in range(n - 1):
        e = np.zeros(n - 1)
        e[b] = eps
        dp.append((populations(n, e, rho0, times)[0] - populations(n, -e, rho0, times)[0]) / (2 * eps))
    dp = np.array(dp)
    w = np.where(p > 1e-13, 1.0 / np.maximum(p, 1e-300), 0.0)
    return np.einsum("ats,bts,ts->tab", dp, dp, w), p, dp


def decoder_efficiency(n, rho0, times, bond):
    """The decoder's own estimate of one known bond's delta from equal-split shots on the schedule: the projection
    of the linear letters, linear in the measured populations; its per-shot variance by propagation through the
    multinomial covariance, set against the single-bond Cramer-Rao bound 1/F_bb of the same shots."""
    Ft, p, dp = position_fisher(n, rho0, times)
    F = Ft.mean(0)
    pdot = populations(n, np.zeros(n - 1), rho0, times)[1]
    sl = -2 * (1 - 2 * p)
    g = times[:, None] * sl * pdot
    den = (g * g).sum(0)
    letter = np.array([(g * sl * dp[k]).sum(0) / den for k in range(n - 1)])
    fb = letter[bond]
    c = fb[None, :] * g * sl / den[None, :] / (fb @ fb)
    var = sum(c[t] @ ((np.diag(p[t]) - np.outer(p[t], p[t])) * len(times)) @ c[t] for t in range(len(times)))
    return (1 / F[bond, bond]) / var, float((c * dp[bond]).sum())


def part_five():
    global J
    say("")
    say("[5] at equal cost: per-shot Fisher information of the full Z-position readout (one shot, one position),")
    say("    every bond's relative delta a parameter; lambda_min = the worst direction, every bond jointly")
    # the two handshakes with the house's Fisher apparatus
    J = 1.0
    times = (1 / GAMMA) * np.arange(1, 81) / 80
    rho = np.outer(mode(4, 1), mode(4, 1))
    p0 = populations(4, np.zeros(3), rho, times)[0]
    p1 = populations(4, np.array([0.02, 0, 0]), rho, times)[0]
    say("  handshake 1: max over K in (0,1] of the bond-0 FI at N = 4, Q = 20, forward difference at 0.02: %.4f"
        " (ReadingPowerWitness reads 1.5934)" % (((p1 - p0) / 0.02) ** 2 / p0).sum(1).max())
    J = 0.075
    sched = np.array([2, 4, 6, 8, 10, 12, 15, 18, 24, 30], float)
    vals = []
    for v in (mode(7, 1), unit(7, 3), unit(7, 0)):
        Ft = position_fisher(7, np.outer(v, v), sched)[0]
        vals.append(np.linalg.eigvalsh(Ft.sum(0) * 1e5 / J ** 2).min())
    say("  handshake 2: N = 7, Q = 1.5, THE_BOND_AND_THE_SEAT_WITH_FINITE_SHOTS's baseline schedule (uniform profile), 1e5 shots"
        " per time, summed over its 1e6 shots, additive units: psi1 %.0f (that page: 190,443); on the same schedule the centre"
        " site %.0f, site 0 %.0f" % tuple(vals))
    sel = np.array([8.5, 9, 9.5, 10, 10.5, 24, 25.5, 27, 28.5, 30])
    vals = []
    for v in (mode(7, 1), unit(7, 3), unit(7, 0)):
        Ft = position_fisher(7, np.outer(v, v), sel)[0]
        vals.append(np.linalg.eigvalsh(Ft.sum(0) * 1e5 / J ** 2).min())
    say("  the same on that page's selected schedule (its exchange search for psi1): psi1 %.0f (that page: 344,084), centre"
        " site %.0f, site 0 %.0f" % tuple(vals))
    J = 1.0
    say("  Q = 20, the decoder's window (60 equal-split times on (0, 1/gamma]): lambda_min, the overlap of its eigenvector with"
        " the uniform rescaling u = (1, .., 1)/sqrt(N-1) of every bond, lambda_min on the directions orthogonal to u, and the"
        " worst single-bond variance bound max_b (F^-1)_bb, all per shot")
    win = np.linspace(1 / GAMMA / 60, 1 / GAMMA, 60)
    for n in (5, 7, 9):
        m = n - 1
        u = np.ones(m) / np.sqrt(m)
        Vc = np.linalg.svd(np.eye(m) - np.outer(u, u))[0][:, :m - 1]
        for name, rho0 in (("bonding mode psi1", np.outer(mode(n, 1), mode(n, 1))), ("site 0", np.outer(unit(n, 0), unit(n, 0))),
                           ("centre site", np.outer(unit(n, n // 2), unit(n, n // 2))),
                           ("end mixture", 0.5 * (np.outer(unit(n, 0), unit(n, 0)) + np.outer(unit(n, n - 1), unit(n, n - 1))))):
            F = position_fisher(n, rho0, win)[0].mean(0)
            w, V = np.linalg.eigh(F)
            say("    N = %d  %-18s lambda_min %.3e, overlap with u %.3f, lambda_min orthogonal to u %.3e, max_b (F^-1)_bb %.3e"
                % (n, name, w[0], abs(V[:, 0] @ u), np.linalg.eigvalsh(Vc.T @ F @ Vc).min(), np.diag(np.linalg.inv(F)).max()))
    for Jq, label in ((1.0, "Q = 20"), (0.075, "Q = 1.5")):
        J = Jq
        say("  %s: lambda_min per shot on equal-split schedules of 60 times over K = gamma t in (0, 1], (0, 3], (0, 10],"
            " and at the best single time of the three grids; mean single-bond F_bb at the best of these" % label)
        for n in (5, 7, 9):
            car = [("bonding mode psi1", np.outer(mode(n, 1), mode(n, 1))), ("mode psi2", np.outer(mode(n, 2), mode(n, 2))),
                   ("site 0", np.outer(unit(n, 0), unit(n, 0))), ("centre site", np.outer(unit(n, n // 2), unit(n, n // 2))),
                   ("end mixture", 0.5 * (np.outer(unit(n, 0), unit(n, 0)) + np.outer(unit(n, n - 1), unit(n, n - 1))))]
            for name, rho0 in car:
                cols, best, best_f, ls = [], 0.0, None, 0.0
                for K in (1, 3, 10):
                    Ft = position_fisher(n, rho0, np.linspace(K / GAMMA / 60, K / GAMMA, 60))[0]
                    F = Ft.mean(0)
                    lm = np.linalg.eigvalsh(F).min()
                    cols.append(lm)
                    if lm > best:
                        best, best_f = lm, F
                    for f in Ft:
                        lt = np.linalg.eigvalsh(f).min()
                        ls = max(ls, lt)
                        if lt > best:
                            best, best_f = lt, f
                eff = np.mean([decoder_efficiency(n, rho0, np.linspace(1 / GAMMA / 60, 1 / GAMMA, 60), b)[0] for b in range(n - 1)])
                say("    N = %d  %-18s lambda_min K<=1 %.2e, K<=3 %.2e, K<=10 %.2e, best single time %.2e; best %.2e;"
                    " mean F_bb at the best %.2e; decoder's single-bond efficiency (K<=1, mean over bonds) %.2f"
                    % (n, name, cols[0], cols[1], cols[2], ls, best, np.mean(np.diag(best_f)), eff))
    J = 1.0


def main():
    say("The decoder and the front: J = 1, gamma = 0.05 (Q = 20), open XY chain, single-excitation block")
    part_one()
    part_two_three()
    part_two_b()
    part_four()
    part_five()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(("\n".join(LINES) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main()
