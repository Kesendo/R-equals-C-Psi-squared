"""The comb test on the road (arc the_forced_and_the_met, item 1): F129's collisions under F160's crack.

F129's collisions live on the chain comb E_k = 2 cos(k pi/n), n = N+1, as EXACT equalities of two
clean-triple sums S(tau) = S(sigma). F160's road H(u) = H_chain + u*V, V = |0><N-1| + |N-1><0|, leaves the
comb at u = 0 (the wrap bond of the ring, u = J'/J). First order in u, nondegenerate perturbation theory
(every chain level is simple):

    dE_k/du |_{u=0} = <psi_k|V|psi_k> = 2 psi_k(0) psi_k(N-1) = (-1)^{k+1} (4/n) sin^2(k pi/n)

which is F65's endpoint rate comb alpha_k/gamma_0 with an alternating sign (psi_k(N-1) = (-1)^{k+1} psi_k(0)
by the chain reflection). So a collision stands at first order on the road iff the SIGNED F65 sums agree:

    D(tau, sigma) = sum_tau (-1)^{k+1} w_k - sum_sigma (-1)^{k+1} w_k,   w_k = (4/n) sin^2(k pi/n) = (2/n)(1 - cos(2k pi/n)).

D is an element of Q(zeta_2n) and is decided EXACTLY (integer vectors mod Phi_2n, the F129 script's own layer):
D != 0 is a theorem that the pair separates linearly on the road; D == 0 leaves the pair to second order.
(The word "survivor" is not used for those pairs: it is MirrorWorld's Survivor, the slowest mode, and F123's noun.)

Gates (exit 0 iff all pass):
  [R1] the velocity identity, exact for small n by sympy on the closed-form eigenvectors, and numerically
       (finite-difference eigenvalues vs the closed form) at every firing n.
  [R2] the chiral-pair control: a balanced pair k + k' = n (the pairs CLEAN excludes) has D = 0 at odd n
       (even N, the cracked ring stays bipartite) and D != 0 at even n (odd N, bipartiteness lost). Exact.
  [R3] every F129 collision pair at every firing n <= N_HI is enumerated exactly, D computed exactly, and
       the count of pairs left to second order (D == 0) reported; the verdict per n. R3b: at odd n
       D = (4/n)(o_tau - o_sigma) exactly (Galois k -> k(n+2)). R3c: at even n every pair left to second
       order is a collision of two ZERO-SUM triples, parity-uniform of one class, with equal doubled sums.
  [R4] the decade law for the derivative, on EVERY pair with D != 0 (the exact test, not a float): [S_tau(u) -
       S_sigma(u)]/u -> D with a residual whose leading power is 1 or 2 (one decade per decade of u, or two where
       the next coefficient vanishes), read from 40-digit eigenvalues at u = 1e-2, 1e-3, 1e-4. The class is read
       at the FINER decade pair (1e-3 -> 1e-4), where the leading power dominates, and the finer pair's ratio must
       be no farther from a pure power than the coarser one's (the residual converges to its
       power law; at u = 1e-2 the u^3 term still competes with the u^4 one, and four pairs at
       n = 30 show it, their c2 being EXACTLY zero by F161's dM_3 = 0). A first version selected the
       pairs with a float D == 0.0 and let five exact zeros through as "dissolving"; the exact test was three lines
       above. A second version read the class at the coarse pair and failed on those four.
  [R4b] EVERY pair left to second order leaves there: gap/u^2 -> c2 != 0 by a decade law.
  [R4c] the parity of the gap in u for EVERY pair left to second order: for a Theta-mirror pair sigma = n - tau
       at even n the gap is EVEN in u (K H(u) K = -H(-u) at odd N, ChiralK), odd part at the 40-digit floor; for
       the others the odd part is read and its leading power reported (u^5 wherever it was read).
  [R5] the finite-u READING (not a certificate): at u = 1/2 the smallest gap between distinct clean-triple
       sums, over all pairs, at 40 digits; and the former collision pairs' gaps.

Run: python simulations/comb_road_f129.py [--fast]
"""
import itertools
import sys

import numpy as np
import sympy as sp
import mpmath as mp

sys.path.insert(0, "simulations")
import f129_level_collision_law as f129

FAST = "--fast" in sys.argv
N_HI = 18 if FAST else 30
FIRING = [n for n in range(5, N_HI + 1) if (n % 3 == 0 and n >= 9) or (n % 10 == 0 and n >= 20)]
FAIL = []


def gate(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


def sign(k):
    return 1 if k % 2 == 1 else -1


def n_times_D_vec(n, tau, sigma):
    """n * D(tau, sigma) as an exact cyclotomic vector: n*w_k = 2 - 2cos(2k pi/n), 2cos(2k pi/n) = zeta^{2k} + zeta^{-2k}."""
    m = 2 * n
    vec = np.zeros(len(f129.root_sum_vec(n, [0])), dtype=np.int64)
    const = 0
    for k in tau:
        const += 2 * sign(k)
        vec -= sign(k) * f129.root_sum_vec(n, [2 * k, m - 2 * k])
    for k in sigma:
        const -= 2 * sign(k)
        vec += sign(k) * f129.root_sum_vec(n, [2 * k, m - 2 * k])
    vec[0] += const
    return vec


def D_float(n, tau, sigma):
    def w(k):
        return (4.0 / n) * np.sin(k * np.pi / n) ** 2
    return sum(sign(k) * w(k) for k in tau) - sum(sign(k) * w(k) for k in sigma)


def chain_h(N, u, prec=None):
    if prec is None:
        h = np.zeros((N, N))
        for a in range(N - 1):
            h[a, a + 1] = h[a + 1, a] = 1.0
        h[0, N - 1] = h[N - 1, 0] = u
        return h
    h = mp.zeros(N, N)
    for a in range(N - 1):
        h[a, a + 1] = h[a + 1, a] = mp.mpf(1)
    h[0, N - 1] = h[N - 1, 0] = mp.mpf(u)
    return h


def road_levels_mp(N, u):
    """The road's N levels at 40 digits, labelled k = 1..N by descending order (the chain comb's order)."""
    ev = mp.eigsy(chain_h(N, u, prec=True), eigvals_only=True)
    ev = sorted([mp.mpf(x) for x in ev], reverse=True)
    return {k + 1: ev[k] for k in range(N)}


def exact_collision_pairs(n):
    """All exact F129 collision pairs (tau, sigma) of clean triples at n: the F129 script's mod-p groups, refuted exactly."""
    pairs = []
    for grp in f129.collision_groups_mod_p(n):
        exact = {}
        for t in grp:
            exact.setdefault(f129.level_vec(n, t).tobytes(), []).append(t)
        for ts in exact.values():
            for a, b in itertools.combinations(ts, 2):
                pairs.append((a, b))
    return pairs


# ---------------------------------------------------------------- R1 the velocity identity
def gate_r1():
    ok_exact = True
    for n in (5, 6, 7, 9):
        N = n - 1
        for k in range(1, N + 1):
            psi = [sp.sqrt(sp.Rational(2, n)) * sp.sin(sp.pi * k * (i + 1) / n) for i in range(N)]
            lhs = 2 * psi[0] * psi[N - 1]
            rhs = (-1) ** (k + 1) * sp.Rational(4, n) * sp.sin(k * sp.pi / n) ** 2
            if sp.simplify(sp.expand_trig(lhs - rhs)) != 0:      # exact; a first version carried a dead 1e-25 fallback
                ok_exact = False
    gate("R1a velocity identity 2 psi_k(0) psi_k(N-1) = (-1)^(k+1) (4/n) sin^2(k pi/n), exact (sympy) at n = 5, 6, 7, 9", ok_exact)
    worst = 0.0
    for n in FIRING:
        N = n - 1
        du = 1e-6
        hp = np.linalg.eigvalsh(chain_h(N, du))[::-1]
        hm = np.linalg.eigvalsh(chain_h(N, -du))[::-1]
        vel = (hp - hm) / (2 * du)
        pred = np.array([sign(k) * (4.0 / n) * np.sin(k * np.pi / n) ** 2 for k in range(1, N + 1)])
        worst = max(worst, np.max(np.abs(vel - pred)))
    gate("R1b the same identity from finite-difference eigenvalues at every firing n (central difference)",
         worst < 1e-8, f"worst |dE/du - (-1)^(k+1) alpha_k/gamma_0| = {worst:.1e}")


# ---------------------------------------------------------------- R2 chiral-pair control
def gate_r2():
    ok = True
    for n in range(5, N_HI + 1):
        for k in range(1, n):
            kp = n - k
            if kp <= k:
                continue
            zero = not np.any(n_times_D_vec(n, (k, kp), ()))
            if n % 2 == 1 and not zero:
                ok = False
            if n % 2 == 0 and zero:
                ok = False
    gate("R2 balanced pair k + k' = n: first-order sum derivative exactly 0 at odd n (even N, bipartite kept) and nonzero at even n (odd N)", ok)


# ---------------------------------------------------------------- R3 the exact census on the road
def gate_r3():
    table = {}
    for n in FIRING:
        pairs = exact_collision_pairs(n)
        standing = [(a, b) for (a, b) in pairs if not np.any(n_times_D_vec(n, a, b))]
        table[n] = (pairs, standing)
        print(f"    n = {n:3d} (N = {n-1:3d}): {len(pairs):5d} exact collision pairs, "
              f"{len(pairs) - len(standing):5d} dissolve at first order (D != 0, a theorem), "
              f"{len(standing):3d} stand at first order (D == 0), left to second order")
        if standing and len(standing) <= 12:
            for a, b in standing:
                print(f"         left to second order: {a} ~ {b}")
    leak = sum(1 for n in FIRING for (a, b) in table[n][1] if D_float(n, a, b) != 0.0)
    total_standing = sum(len(table[n][1]) for n in FIRING)
    mirror_odd = sum(1 for n in FIRING if n % 2 == 1 for (a, b) in table[n][0] if tuple(sorted(n - k for k in a)) == tuple(b))
    mirror_even = sum(1 for n in FIRING if n % 2 == 0 for (a, b) in table[n][0] if tuple(sorted(n - k for k in a)) == tuple(b))
    print(f"    of the {total_standing} pairs with D == 0 exactly, {leak} have a float D != 0.0 (what a float test would have let through)")
    print(f"    Theta-mirror collision pairs sigma = n - tau: {mirror_odd} at odd n (none of them stands), {mirror_even} at even n (all of them stand)")
    gate("R3 every firing n <= N_HI has its exact collision pairs enumerated and D decided exactly",
         all(len(v[0]) > 0 for v in table.values()))
    # R3b, the standing' law at odd n, a THEOREM. (-1)^{k+1} cos(2k pi/n) = -cos(k(n+2) pi/n), and for odd n
    # the map k -> k(n+2) mod 2n is a Galois automorphism of Q(zeta_2n) (gcd(n+2, 2n) = 1), so the signed cosine
    # part of D is minus the Galois conjugate of S(tau) - S(sigma) = 0 and vanishes for EVERY collision. What is
    # left is the constant part: D = (4/n) (o_tau - o_sigma), o = the number of ODD labels in the triple. So at
    # odd n a collision stands at first order on the road iff its two triples carry the same number of odd labels.
    ok = True
    for n in FIRING:
        if n % 2 == 0:
            continue
        pairs, standing = table[n]
        def odd_count(t):
            return sum(1 for k in t if k % 2 == 1)
        for a, b in pairs:
            vec = n_times_D_vec(n, a, b)                 # n * D, exact; predicted (4)(o_a - o_b) in the constant slot
            pred = np.zeros_like(vec)
            pred[0] = 4 * (odd_count(a) - odd_count(b))
            if np.any(vec - pred):
                ok = False
        predicted = [(a, b) for (a, b) in pairs if odd_count(a) == odd_count(b)]
        if set(predicted) != set(standing):
            ok = False
        hist = {}
        for a, b in pairs:
            dd = abs(odd_count(a) - odd_count(b))
            hist[dd] = hist.get(dd, 0) + 1
        print(f"    n = {n}: |o_tau - o_sigma| over the pairs: {dict(sorted(hist.items()))} (|D| = that times 4/n)")
        print(f"    n = {n}: D == (4/n)(o_tau - o_sigma) on all {len(pairs)} pairs; pairs left to second order == equal-odd-count pairs ({len(standing)})")
    gate("R3b odd n: D = (4/n)(o_tau - o_sigma) exactly on every collision pair (Galois k -> k(n+2)); the pairs left to second order are the equal-odd-count pairs", ok)
    # even n: k -> k(n+2) is not an automorphism (n+2 even); the standing are read directly. Every one found is
    # parity-uniform of the same class, where D = -/+ (2/n)(sum_sigma c - sum_tau c), so D = 0 iff the DOUBLED
    # sums agree; for the twelve non-mirror ones F161 names the shape that forces it, the doubled
    # labels being ROT3; the eleven Theta-mirror ones are term-by-term and need no shape.
    ok_even = True
    n_zero = n_mirror = n_total = 0
    for n in FIRING:
        if n % 2 == 1:
            continue
        pairs, standing = table[n]
        m = 2 * n
        for a, b in standing:
            da = f129.root_sum_vec(n, [e for k in a for e in (2 * k, m - 2 * k)])
            db = f129.root_sum_vec(n, [e for k in b for e in (2 * k, m - 2 * k)])
            pa = {k % 2 for k in a}
            pb = {k % 2 for k in b}
            uniform = len(pa) == 1 and pa == pb
            zero_sum = not np.any(f129.level_vec(n, a))          # exact: S(tau) = 0 in Z[zeta_2n]
            mirror = tuple(sorted(n - k for k in a)) == tuple(b)   # sigma = n - tau, F129's Theta-mirror pair
            n_total += 1
            n_zero += zero_sum
            n_mirror += mirror
            if not uniform or np.any(da - db) or not zero_sum:
                ok_even = False
            print(f"    n = {n} left to second order {a} ~ {b}: zero-sum collision = {zero_sum}, Theta-mirror pair = {mirror}, "
                  f"parity-uniform same class = {uniform}, doubled sums equal = {not np.any(da - db)}, doubled sum zero = {not np.any(da)}")
    gate("R3c even n: every pair left to second order is a collision of two ZERO-SUM triples (exact), parity-uniform of one class, with equal doubled sums",
         ok_even, f"{n_total} pairs, {n_zero} zero-sum, {n_mirror} of them Theta-mirror pairs sigma = n - tau")
    table["_mirror_flags"] = True
    return table


# ---------------------------------------------------------------- R4 decade law, R5 finite-u reading
def decade_class(r):
    if abs(r / 10.0 - 1.0) <= 0.15:
        return 1
    if abs(r / 100.0 - 1.0) <= 0.15:
        return 2
    return 0


def gate_r4_r5(table):
    mp.mp.dps = 40
    worst_ratio_lo, worst_ratio_hi = 1e9, 0.0
    nrows = 0
    ratios, not_converging, coarse_between = [], [], []
    for n in FIRING:
        N = n - 1
        pairs, standing = table[n]
        lv = {u: road_levels_mp(N, u) for u in (1e-2, 1e-3, 1e-4)}
        for a, b in pairs:
            if not np.any(n_times_D_vec(n, a, b)):      # exact: D == 0 in Z[zeta_2n]; a float test let five through
                continue
            D = D_float(n, a, b)
            res = []
            for u in (1e-2, 1e-3, 1e-4):
                gap = sum(lv[u][k] for k in a) - sum(lv[u][k] for k in b)
                res.append(abs(float(gap / u) - D))
            r1, r2 = res[0] / res[1], res[1] / res[2]
            worst_ratio_lo = min(worst_ratio_lo, r1, r2)
            worst_ratio_hi = max(worst_ratio_hi, r1, r2)
            ratios.append(r2)
            dist = lambda r: min(abs(np.log10(r) - pw) for pw in (1, 2))
            if dist(r2) > dist(r1) + 1e-9:                 # the finer pair must sit no farther from a pure power
                not_converging.append((n, a, b, r1, r2))
            if decade_class(r1) == 0:
                coarse_between.append((n, a, b, r1, r2))
            nrows += 1
    classes = [decade_class(r) for r in ratios]
    gate("R4 decade law on every pair with D != 0: at the finer decade pair (1e-3 -> 1e-4) the residual's ratio is within 15% of 10 or of 100, and the finer pair's ratio is never farther from a pure power than the coarser one's",
         all(c > 0 for c in classes) and not not_converging,
         f"{nrows} pairs; finer-pair ratios: {sum(1 for c in classes if c == 1)} at 10, {sum(1 for c in classes if c == 2)} at 100; "
         f"{len(coarse_between)} coarse ratios sit between the classes (all converging); all ratios in [{worst_ratio_lo:.2f}, {worst_ratio_hi:.2f}]")
    for n, a, b, r1, r2 in coarse_between:
        print(f"    coarse pair between classes: n = {n} {a} ~ {b}: 1e-2/1e-3 = {r1:.2f}, 1e-3/1e-4 = {r2:.2f}")
    # R4b the standing leave at second order: gap(u)/u^2 tends to a nonzero constant, its residual falling one
    # decade per decade of u; and no standing pair has a zero second-order coefficient (which would ask for FOURTH order: F161
    # proves the third vanishes for every standing pair).
    lo2, hi2, nstand, minc2 = 1e9, 0.0, 0, 1e9
    ratios2 = []
    for n in FIRING:
        N = n - 1
        pairs, standing = table[n]
        if not standing:
            continue
        lv = {u: road_levels_mp(N, u) for u in (1e-2, 1e-3, 1e-4, 1e-5)}
        for a, b in standing:
            c2 = []
            for u in (1e-2, 1e-3, 1e-4, 1e-5):
                gap = sum(lv[u][k] for k in a) - sum(lv[u][k] for k in b)
                c2.append(float(gap / (u * u)))
            c2_lim = c2[-1]
            minc2 = min(minc2, abs(c2_lim))
            res = [abs(c - c2_lim) for c in c2[:-1]]
            if res[1] > 0 and res[2] > 0:
                r1, r2 = res[0] / res[1], res[1] / res[2]
                lo2, hi2 = min(lo2, r1, r2), max(hi2, r1, r2)
                ratios2 += [r1, r2]
            nstand += 1
    if nstand:
        classes2 = [decade_class(r) for r in ratios2]
        gate("R4b the pairs left to second order leave THERE: gap(u)/u^2 -> c2 != 0, the residual falling one or two decades per decade (within 15% of 10 or 100); smallest |c2| far above the 40-digit floor",
             all(c > 0 for c in classes2) and minc2 > 1e-6,
             f"{nstand} pairs; ratios in [{lo2:.2f}, {hi2:.2f}]; {sum(1 for c in classes2 if c == 1)} at 10, {sum(1 for c in classes2 if c == 2)} at 100; smallest |c2| = {minc2:.3e}")
    # R4c the parity of the gap in u. K = diag((-1)^site) flips every chain bond; the wrap bond joins sites 0 and
    # N-1, so K V K = (-1)^(N-1) V: at odd N (even n) K H(u) K = -H(-u), hence E_k(u) = -E_{n-k}(-u) with the levels
    # labelled by descending order, and for a Theta-mirror pair sigma = n - tau the gap S_tau(u) - S_sigma(u) =
    # S_tau(u) + S_tau(-u) is EVEN in u to all orders (the shape of PROOF_ZETA2_ANTI_PROTECTION, F131 Theorem B,
    # with u for zeta). For the other pairs left to second order the odd part is read: it starts at u^5.
    ok_par, n_mir, n_other, worst_mirror, ratios5 = True, 0, 0, 0.0, []
    powers = {}
    for n in FIRING:
        N = n - 1
        pairs, standing = table[n]
        if not standing:
            continue
        lv = {u: road_levels_mp(N, u) for u in (1e-2, -1e-2, 1e-3, -1e-3)}
        for a, b in standing:
            def gap(u):
                return sum(lv[u][k] for k in a) - sum(lv[u][k] for k in b)
            odd1 = (gap(1e-2) - gap(-1e-2)) / 2
            odd2 = (gap(1e-3) - gap(-1e-3)) / 2
            mirror = n % 2 == 0 and tuple(sorted(n - k for k in a)) == tuple(b)
            if mirror:
                n_mir += 1
                worst_mirror = max(worst_mirror, float(abs(odd1)), float(abs(odd2)))
                if abs(odd1) > mp.mpf(10) ** -35 or abs(odd2) > mp.mpf(10) ** -35:
                    ok_par = False
            else:
                n_other += 1
                r = float(abs(odd1 / odd2)) if odd2 != 0 else float('nan')
                ratios5.append(r)
                pw = round(np.log10(r)) if r == r and r > 0 else None
                powers[pw] = powers.get(pw, 0) + 1
                if pw is None or abs(r / 10.0 ** pw - 1.0) > 0.15:
                    ok_par = False
    gate("R4c every pair left to second order: a Theta-mirror pair (even n) has a gap EVEN in u (odd part at the 40-digit floor); every other pair's odd part has an integer leading power, reported",
         ok_par, f"{n_mir} mirror pairs, worst odd part {worst_mirror:.1e}; {n_other} others, odd part ~ u^p with p counted as {dict(sorted(powers.items()))}")
    print("  R5 the finite-u reading (u = 1/2, 40 digits; a reading, not a certificate):")
    for n in FIRING:
        N = n - 1
        lv = road_levels_mp(N, mp.mpf(1) / 2)
        clean = [t for t in itertools.combinations(range(1, N + 1), 3) if f129.is_clean(n, t)]
        sums = sorted((sum(lv[k] for k in t), t) for t in clean)
        gaps = [(sums[i + 1][0] - sums[i][0], sums[i][1], sums[i + 1][1]) for i in range(len(sums) - 1)]
        gmin = min(gaps, key=lambda g: g[0])
        pairs, _ = table[n]
        former = [abs(sum(lv[k] for k in a) - sum(lv[k] for k in b)) for a, b in pairs]
        print(f"    n = {n:3d}: {len(clean)} clean triples; smallest gap between distinct sums {mp.nstr(gmin[0], 6)} "
              f"({gmin[1]} vs {gmin[2]}); former collision pairs' gaps in [{mp.nstr(min(former), 4)}, {mp.nstr(max(former), 4)}]")


if __name__ == "__main__":
    print(f"The comb test on the road: F129 under F160, firing n = {FIRING}")
    gate_r1()
    gate_r2()
    tab = gate_r3()
    gate_r4_r5(tab)
    print("ALL GATES PASS" if not FAIL else f"FAILED: {FAIL}")
    sys.exit(0 if not FAIL else 1)
