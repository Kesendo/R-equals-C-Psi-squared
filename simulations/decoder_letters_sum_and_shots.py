"""The decoder and the front, the open items (experiments/THE_DECODER_AND_THE_FRONT.md §4 and §5): why a carrier that
holds few modes reads nearly opposite letters, an actual decode at finite shots, and the decoder's de-lossed path with
other carriers. Built on simulations/decoder_carrier_dictionary.py (the same exact propagation of the single-excitation
block, J = 1, gamma = 0.05 unless stated), so the two scripts meet on the same letters.

  [1] the sum rule: h = sum_b V_b, so the letters of all bonds add up to the letter of the uniform rescaling U of every
      bond; at gamma = 0 that rescaling is a rescaling of time and U = t dP/dt, the carrier's own motion. A carrier at
      rest under H has U = 0 in its clean evolution: its letters cancel. N = 3 (two letters, exactly opposite) and
      N = 4 (three letters, their cosines fixed by their norms); the bonding mode's weakest direction against the sum;
      the cancellation kappa = |U| / sum_b |l_b| against the worst signed cos over the 132 random carriers of
      decoder_carrier_dictionary.py [4], beside the mode participation PR, on both observables; the dominance of one
      direction (the rank-one share) as the null control.
  [2] an actual decode at finite shots on the decoder's window (60 equal-split times on (0, 1/gamma]): (A) one defect,
      maximum likelihood over (bond, delta) on a grid, the probability of the right bond and sign, and where the wrong
      answers go; (B) every bond at once, maximum likelihood by Fisher scoring run to convergence, the covariance of
      the error against the Cramer-Rao bound at the true deltas, with a bootstrap spread.
  [3] the de-lossed path (DefectDecoder.DecodeDeviation: per site the mean over the painters' grid of the purity's
      deviation, one letter per bond calibrated at delta = +0.02, the single-letter rule) for other carriers.
Output: simulations/results/decoder_carrier/decoder_letters_sum_and_shots.txt
"""
from collections import Counter
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import decoder_carrier_dictionary as D

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "simulations" / "results" / "decoder_carrier" / "decoder_letters_sum_and_shots.txt"
LINES = []
TIMES = np.arange(1, 81) * 0.25


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def uniform_letter(n, rho0, d=1e-4):
    e = np.full(n - 1, d)
    return ((D.populations(n, e, rho0, TIMES)[0] - D.populations(n, -e, rho0, TIMES)[0]) / (2 * d)).reshape(-1)


def kappa(F, U):
    return np.linalg.norm(U) / np.linalg.norm(F, axis=1).sum()


def cosines(F):
    nr = np.linalg.norm(F, axis=1)
    return (F @ F.T) / np.outer(nr, nr), nr


def random_carriers(n=7):
    """The 132 carriers of decoder_carrier_dictionary.py [4], drawn in the same order from the same seed."""
    rng = np.random.default_rng(20261030)
    M = np.array([D.mode(n, k) for k in range(1, n + 1)])
    out = []
    for m in range(1, n + 1):
        for _ in range(12):
            c = np.zeros(n)
            c[rng.choice(n, m, replace=False)] = rng.normal(size=m)
            phi = c @ M
            out.append(phi / np.linalg.norm(phi))
    for width in (0.5, 1.0, 2.0, 4.0):
        for _ in range(12):
            x0 = rng.uniform(0, n - 1)
            phi = np.exp(-((np.arange(n) - x0) ** 2) / (4 * width ** 2))
            out.append(phi / np.linalg.norm(phi))
    return M, out


def part_one():
    say("[1] the sum rule: the letters of all bonds add up to the letter U of the uniform rescaling")
    say("  the bare location dictionary M[b, k] = <psi_k|V_b|psi1>, k = 2..N (rank N - 2, the K-partner column k = N")
    say("  zero): every column sums over the bonds to <psi_k|h|psi1> = 0, so its one left null vector is (1, .., 1)")
    for n in range(3, 11):
        M = np.zeros((n - 1, n - 1))
        for b in range(n - 1):
            V = np.zeros((n, n))
            V[b, b + 1] = V[b + 1, b] = 1.0
            for k in range(2, n + 1):
                M[b, k - 2] = D.mode(n, k) @ V @ D.mode(n, 1)
        U_, s, Vt = np.linalg.svd(M.T)
        left = Vt[-1]
        say("    N = %2d  rank %d, max |column sum| %.1e, |K-partner column| %.1e, overlap of the left null vector with"
            " (1, .., 1)/sqrt(N-1) %.15f" % (n, int((s > 1e-12).sum()), np.abs(M.sum(0)).max(), np.abs(M[:, -1]).max(),
                                             abs(left @ np.ones(n - 1)) / np.sqrt(n - 1)))
    for g in (0.0, 0.05):
        D.GAMMA = g
        say("  gamma = %.2f" % g)
        for n in (5, 7):
            for name, rho0 in D.carriers(n):
                F = D.trajectory_letters(n, rho0)
                U = uniform_letter(n, rho0)
                res = np.linalg.norm(F.sum(0) - U) / np.linalg.norm(F, axis=1).sum()
                line = "    N = %d  %-46s |sum_b l_b - U| / sum_b |l_b| %.1e   kappa %.3f" % (n, name, res, kappa(F, U))
                pd = D.populations(n, np.zeros(n - 1), rho0, TIMES)[1]
                tP = (TIMES[:, None] * pd).reshape(-1)
                if g > 0.0:
                    # the generator is homogeneous of degree one in (h, gamma): U = t dp/dt - gamma dp/dgamma
                    eg = 1e-6
                    D.GAMMA = g + eg
                    up = D.populations(n, np.zeros(n - 1), rho0, TIMES)[0]
                    D.GAMMA = g - eg
                    dn = D.populations(n, np.zeros(n - 1), rho0, TIMES)[0]
                    D.GAMMA = g
                    tP = tP - g * ((up - dn) / (2 * eg)).reshape(-1)
                    line += "   |U - (t dP/dt - gamma dP/dgamma)| / |U| %.1e" % (np.linalg.norm(U - tP) / np.linalg.norm(U))
                elif np.linalg.norm(U) > 1e-6:
                    line += "   |U - t dP/dt| / |U| %.1e" % (np.linalg.norm(U - tP) / np.linalg.norm(U))
                else:
                    line += "   |U| %.1e, |t dP/dt| %.1e" % (np.linalg.norm(U), np.linalg.norm(tP))
                say(line)
    say("  every mode psi_k a carrier at rest: N = 3, the two letters' cos; N = 4, the three cosines against the triangle")
    say("  law cos_ab = (|l_c|^2 - |l_a|^2 - |l_b|^2) / (2 |l_a| |l_b|) of three vectors adding up to zero; exact without")
    say("  dephasing, broken at Q = 20 by the halo")
    for g, n, k in [(g, n, k) for g in (0.0, 0.05) for n in (3, 4) for k in range(1, n + 1)]:
        D.GAMMA = g
        if True:
            C, nr = cosines(D.trajectory_letters(n, np.outer(D.mode(n, k), D.mode(n, k))))
            if n == 3:
                say("    gamma = %.2f  N = 3  psi%d  cos %+.12f" % (g, k, C[0, 1]))
            else:
                tri = lambda a, b, c: (nr[c] ** 2 - nr[a] ** 2 - nr[b] ** 2) / (2 * nr[a] * nr[b])
                say("    gamma = %.2f  N = 4  psi%d  cos01 %+.8f (law %+.8f)  cos02 %+.8f (law %+.8f)  cos12 %+.8f (law %+.8f)"
                    % (g, k, C[0, 1], tri(0, 1, 2), C[0, 2], tri(0, 2, 1), C[1, 2], tri(1, 2, 0)))
    say("  the dictionary's weakest direction (the smallest singular vector of the normalised letters, read over the")
    say("  bonds) against the norm vector (|l_0|, .., |l_N-2|), along which the normalised letters add up to U")
    for g in (0.0, 0.05):
        D.GAMMA = g
        for n in (5, 7, 9):
            for name, rho0 in (("bonding mode psi1", np.outer(D.mode(n, 1), D.mode(n, 1))),
                               ("site 0", np.outer(D.unit(n, 0), D.unit(n, 0)))):
                F = D.trajectory_letters(n, rho0)
                C, nr = cosines(F)
                w, V = np.linalg.eigh(C)
                u = nr / np.linalg.norm(nr)
                say("    gamma = %.2f  N = %d  %-18s conditioning %.3f, overlap with the norm vector %.4f, Rayleigh quotient"
                    " of the norm vector over the smallest eigenvalue of the cosine matrix %.3f"
                    % (g, n, name, np.sqrt(max(w[0], 0.0)), abs(V[:, 0] @ u), (u @ C @ u) / w[0] if w[0] > 1e-12 else np.nan))
                off = ~np.eye(n - 1, dtype=bool)
                W = np.outer(nr, nr)
                say("      its worst cos %+.3f, its norm-weighted mean cos %+.3f, the equal-norm floor -1/(N-2) %+.3f"
                    % (C[off].min(), (W * C)[off].sum() / W[off].sum(), -1 / (n - 2)))
    D.GAMMA = 0.05
    say("  the linear letter at Q = 20 per carrier: kappa (the letters' sum over their summed norms), worst and largest cos")
    for n in (5, 7):
        for name, rho0 in D.carriers(n):
            Fl = D.linear_letters(n, rho0)
            Cl, _ = cosines(Fl)
            iu = np.triu_indices(n - 1, 1)
            say("    N = %d  %-46s kappa %.3f, worst cos %+.3f, largest cos %+.3f"
                % (n, name, kappa(Fl, Fl.sum(0)), Cl[iu].min(), Cl[iu].max()))
    say("  the 132 random carriers at N = 7: correlation of the worst signed cos with kappa, with PR, and with the rank-one")
    say("  share s1^2 / sum s^2 of the normalised dictionary (one dominant direction)")
    M, phis = random_carriers()
    for g in (0.0, 0.05):
        D.GAMMA = g
        rows = []
        for phi in phis:
            rho0 = np.outer(phi, phi)
            F = D.trajectory_letters(7, rho0)
            U = uniform_letter(7, rho0)
            C, nr = cosines(F)
            iu = np.triu_indices(6, 1)
            s = np.linalg.svd(F / nr[:, None], compute_uv=False)
            row = [1 / np.sum((M @ phi) ** 4), C[iu].min(), kappa(F, U), s[0] ** 2 / (s ** 2).sum()]
            if g > 0:
                Fl = D.linear_letters(7, rho0)
                Cl, _ = cosines(Fl)
                row += [Cl[iu].min(), Cl[iu].max(), kappa(Fl, Fl.sum(0))]
            rows.append(row)
        R = np.array(rows)
        cc = lambda i, j: np.corrcoef(R[:, i], R[:, j])[0, 1]
        say("    gamma = %.2f trajectory: worst cos with kappa %+.3f, with PR %+.3f, with the rank-one share %+.3f;"
            " kappa %.3f..%.3f where worst < -0.5 (%d), %.3f..%.3f where not"
            % (g, cc(1, 2), cc(1, 0), cc(1, 3), R[R[:, 1] < -0.5, 2].min(), R[R[:, 1] < -0.5, 2].max(),
               int((R[:, 1] < -0.5).sum()), R[R[:, 1] >= -0.5, 2].min(), R[R[:, 1] >= -0.5, 2].max()))
        say("      kappa with PR %+.3f" % cc(2, 0))
        if g > 0:
            say("    gamma = 0.05 linear letter: worst cos with kappa %+.3f, with PR %+.3f; largest cos with kappa %+.3f,"
                " with PR %+.3f" % (cc(4, 6), cc(4, 0), cc(5, 6), cc(5, 0)))
    D.GAMMA = 0.05


WIN = np.linspace(1 / 0.05 / 60, 1 / 0.05, 60)
RNG = np.random.default_rng(31)


def pops(n, deltas, rho0):
    p = np.clip(D.populations(n, deltas, rho0, WIN)[0], 1e-15, None)
    return p / p.sum(1, keepdims=True)


def draw(p, shots):
    return np.array([RNG.multinomial(shots // len(WIN), pt) for pt in p])


def score_fisher(n, rho0, delta, c, eps=1e-6):
    p = pops(n, delta, rho0)
    dp = []
    for b in range(n - 1):
        e = np.zeros(n - 1)
        e[b] = eps
        dp.append((pops(n, delta + e, rho0) - pops(n, delta - e, rho0)) / (2 * eps))
    dp = np.array(dp)
    return np.einsum("btj,tj->b", dp, c / p), np.einsum("atj,btj,tj->ab", dp, dp, c.sum(1)[:, None] / p)


def part_two():
    say("")
    say("[2] an actual decode at finite shots, the decoder's window at Q = 20, shots split equally over the 60 times")
    grid = np.concatenate([-np.linspace(0.06, 0.0005, 120), np.linspace(0.0005, 0.06, 120)])
    shots = [600, 2400, 6000, 24000, 60000, 240000]
    say("  (A) one defect, maximum likelihood over (bond, delta), delta on a grid of 240 values in [-0.06, 0.06]; 300 draws")
    say("      per cell; P(right bond and sign) at %s shots; where the answers go at 60000" % shots)
    for n, tb, td in ((5, 3, -0.02), (5, 1, 0.02), (7, 5, -0.02), (7, 1, 0.02)):
        for name, rho0 in (("bonding mode psi1", np.outer(D.mode(n, 1), D.mode(n, 1))),
                           ("site 0", np.outer(D.unit(n, 0), D.unit(n, 0))),
                           ("centre site", np.outer(D.unit(n, n // 2), D.unit(n, n // 2)))):
            table = np.array([[np.log(pops(n, np.eye(n - 1)[b] * g, rho0)) for g in grid] for b in range(n - 1)])
            p_true = pops(n, np.eye(n - 1)[tb] * td, rho0)
            probs, where = [], Counter()
            for S in shots:
                right = 0
                for _ in range(300):
                    LL = np.einsum("bgtj,tj->bg", table, draw(p_true, S))
                    b, gi = np.unravel_index(np.argmax(LL), LL.shape)
                    right += (b == tb) and (np.sign(grid[gi]) == np.sign(td))
                    if S == 60000:
                        where["bond %d %s" % (b, "+" if grid[gi] > 0 else "-")] += 1
                probs.append(right / 300)
            say("    N = %d, bond %d at %+.2f, %-18s %s   at 60000: %s"
                % (n, tb, td, name, "  ".join("%.3f" % x for x in probs), dict(where.most_common())))
    say("  (B) every bond at once, maximum likelihood by Fisher scoring from zero with step halving until the step is below")
    say("      1e-9 (at most 80 steps); one true configuration per cell, uniform in [-0.01, 0.01] per bond; 300 draws: the")
    say("      covariance of the error against the Cramer-Rao bound F^-1/S at the true deltas, the bootstrap spread of that")
    say("      ratio (500 resamples), the mean squared error over the bound, and the draws not converged")
    for n in (5, 7):
        for S in (60000, 600000):
            for name, rho0 in (("bonding mode psi1", np.outer(D.mode(n, 1), D.mode(n, 1))),
                               ("site 0", np.outer(D.unit(n, 0), D.unit(n, 0)))):
                true = RNG.uniform(-0.01, 0.01, n - 1)
                pt = pops(n, true, rho0)
                crb = np.linalg.inv(score_fisher(n, rho0, true, pt / len(WIN))[1]) / S
                errs, stuck = [], 0
                for _ in range(300):
                    c = draw(pt, S)
                    d = np.zeros(n - 1)
                    ll = float((c * np.log(pops(n, d, rho0))).sum())
                    done = False
                    for _ in range(80):
                        sc, F = score_fisher(n, rho0, d, c)
                        step = np.linalg.solve(F, sc)
                        t = 1.0
                        while True:
                            dn = d + t * step
                            lln = float((c * np.log(pops(n, dn, rho0))).sum())
                            if lln >= ll - 1e-9 or t < 1e-6:
                                break
                            t /= 2
                        d, ll = dn, lln
                        if np.abs(t * step).max() < 1e-9:
                            done = True
                            break
                    stuck += not done
                    errs.append(d - true)
                E = np.array(errs)
                tr = lambda e: np.trace(np.cov(e.T)) / np.trace(crb)
                boot = [tr(E[RNG.integers(0, len(E), len(E))]) for _ in range(500)]
                say("    N = %d  %7d shots  %-18s spread of the error %.2e (bound %.2e), covariance trace over the bound's"
                    " %.2f (bootstrap %.2f..%.2f), mean squared error over the bound %.2f, not converged %d of 300"
                    % (n, S, name, np.sqrt(np.trace(np.cov(E.T))), np.sqrt(np.trace(crb)), tr(E),
                       np.percentile(boot, 2.5), np.percentile(boot, 97.5), (E ** 2).sum(1).mean() / np.trace(crb), stuck))


GRID = np.linspace(0.0, 1.0 / 0.05, 60)


def purity_mean(n, deltas, rho0):
    p = D.populations(n, deltas, rho0, GRID)[0]
    return (0.5 * (1 + (1 - 2 * p) ** 2)).mean(0)


def part_three():
    say("")
    say("[3] the de-lossed path (DecodeDeviation) with other carriers: the deviation dictionary's worst signed cos and")
    say("    conditioning, and the single-letter decode of planted defects (decoded bond and strength, squared residual")
    say("    ratio runner-up over winner); psi1 at N = 5, bond 3 at -0.02, is the live inspect --root decoder case (near 516)")
    for n in (4, 5):
        cars = [("bonding mode psi1", np.outer(D.mode(n, 1), D.mode(n, 1))), ("site 0", np.outer(D.unit(n, 0), D.unit(n, 0))),
                ("end mixture", 0.5 * (np.outer(D.unit(n, 0), D.unit(n, 0)) + np.outer(D.unit(n, n - 1), D.unit(n, n - 1))))]
        if n % 2:
            cars.append(("centre site", np.outer(D.unit(n, n // 2), D.unit(n, n // 2))))
        cases = [(3, -0.02), (0, -0.02), (1, 0.01)] if n == 5 else [(2, -0.02), (1, -0.02), (0, 0.01)]
        for name, rho0 in cars:
            base = purity_mean(n, np.zeros(n - 1), rho0)
            dic = np.array([(purity_mean(n, np.eye(n - 1)[b] * 0.02, rho0) - base) / 0.02 for b in range(n - 1)])
            C, nr = cosines(dic)
            iu = np.triu_indices(n - 1, 1)
            cond = np.sqrt(max(np.linalg.eigvalsh(C)[0], 0.0))
            outs = []
            for tb, td in cases:
                obs = purity_mean(n, np.eye(n - 1)[tb] * td, rho0) - base
                res = sorted((np.sum((obs - (obs @ f / (f @ f)) * f) ** 2), b, obs @ f / (f @ f)) for b, f in enumerate(dic))
                outs.append("bond %d at %+.2f -> bond %d at %+.4f, ratio %.2f" % (tb, td, res[0][1], res[0][2], res[1][0] / res[0][0]))
            say("    N = %d  %-18s worst cos %+.3f, conditioning %.3f;  %s" % (n, name, C[iu].min(), cond, ";  ".join(outs)))


def main():
    say("The decoder and the front, the open items: J = 1, open XY chain, single-excitation block")
    part_one()
    part_two()
    part_three()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(("\n".join(LINES) + "\n").encode("utf-8"))


if __name__ == "__main__":
    main()
