"""Gate for the muteness axis (experiments/LATTICE_DEAD_SET_H_ZERO.md): the
dead-set law's sufficiency is GENERIC in h, and at h = 0 the law refines to a
finer conservation law that collapses to one kinematic line.

Context: the dead-set law F132 (experiments/LATTICE_DEAD_SET_RULE.md,
simulations/lattice_dead_set_rule.py),

    alive  iff  (K_pop and d = 0 mod 4)  or  (coherence and K_coh and d = N),

has its necessity face derived; its sufficiency face (allowed => alive) is the
gated observation, and the F132 gate's censuses scan exactly ONE fixed generic
h direction per config. This gate closes that gap and maps the h axis:

  T1  GENERIC SUFFICIENCY. Across 5 configs (N = 3..5, population and real
      coherence, seeds s = 1 and s = 3) x 6 seeded random h directions and
      scales, the law is exact at EVERY sample, and the separation is not
      marginal: alive readouts sit at >= 5e-3, dead ones at machine zero
      (~2e-16, asserted < 1e-12 vs an alive floor > 1e-3).

  T2  THE STRUCTURED POINTS. At staggered and single-site h the law is
      unchanged (both already generic). At h = 0 and at uniform h the alive
      set SHRINKS by exact, reproducible counts (extra deaths are sightings
      of new structure, not law violations; necessity holds everywhere).

  T3  THE BARE CHIRAL-CONJ SYMMETRY AT h = 0. Without the field, the
      sublattice gauge U_g flips the hopping alone (every bond has exactly
      one site in g), so V'_g(rho) = U_g conj(rho) U_g is a flow symmetry
      WITHOUT the X^N factor (pinned dynamically: V' commutes with the h = 0
      flow). Its kill sign is eps'_g = (-1)^(n_Y + xy_g); the
      per-channel sign rule reproduces the h = 0 alive set exactly at N = 3
      and N = 4, and OVER-predicts at N = 5 s = 3 coherence by exactly the
      10 strings the finer layer below kills: the same discovery shape as
      F132 itself (F first bit at N = 5 with the popcount-2 seed).

  T4  THE BI-DEGREE CONSERVATION LAW. In Majorana language the h = 0 hopping
      graph is DISCONNECTED: the edges (2l+1, 2l+2) and (2l, 2l+3) split the
      2N indices into E = {k mod 4 in {0,3}} and O = {k mod 4 in {1,2}}, and
      the field edges (2l, 2l+1) are the ONLY bridge. So at h = 0 the flow
      conserves the BI-degree (d_E, d_O), refining F132's total degree d
      (pinned: [H0, O] leaks exactly 0.0 outside O's bi-sector). The h = 0
      dead-set rule collapses to ONE kinematic line, no sign condition left:

          alive at h = 0  iff  (d_E, d_O) populated in rho(0)  and  K-readable

      (population sectors are (k, k) with k even, where eps'_g = +1 and
      d = 2k = 0 mod 4 hold identically; the sign rule and F132's mod-4
      condition are both automatic per populated sector). Exact on all 6
      configs, including the 10 strings of T3.

  T5  THE SHADOW IDENTITY, ONE LEVEL DOWN. eps'_odd = (-1)^(d(d-1)/2 + d_O)
      and eps'_even = (-1)^(d(d-1)/2 + d_E): the bare V' is the mod-4 shadow
      of the bi-degree exactly as F132's V is the mod-4 shadow of the total
      degree; the whole tower repeats one level finer. All 4^N strings x
      both gauges, N = 2..6.

  T6  UNIFORM h. The uniform field commutes with the hopping (it is a
      function of total popcount), so on the population channel the dynamics
      are IDENTICAL to h = 0 (the block phases cancel on diagonal blocks)
      while the coherence blocks pick up a rotating phase that revives every
      sign-killed coherence readout: alive(uniform) = alive(h=0) union
      {coherence on, K_coh, d = N}. Exact on all 5 configs.

Honesty split, as everywhere in this thread: the necessity directions of T3,
T4, T6 are conservation/symmetry arguments; every sufficiency (allowed =>
alive) is the gated observation. Alive = max_t |<O>(t)| > 1e-10 over an RK4
window [0, 2], dt = 0.05, one-sided watching as in the F132 gate.

Runtime ~4 min (the N = 5 configs dominate). Scout history: the probe-D sign
rule was the first mechanism guess and its 10-string failure at N = 5 s = 3
was the crack the bi-degree came through; the scout's discovery order
(A generic, B/C structured, D sign rule, E bi-degree, F shadow) is preserved
as T1..T5 here, with the uniform-h union rule added as T6.
"""
import sys
from itertools import product

import numpy as np

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
s0 = np.eye(2, dtype=complex)
PAULI = {"I": s0, "X": sx, "Y": sy, "Z": sz}

ALIVE_TOL = 1e-10
FAILURES = []


def check(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILURES.append(name)


def majorana_degree(name):
    """Left-JW Majorana degree (verbatim from lattice_dead_set_rule.py)."""
    tau = 0
    d = 0
    for l in range(len(name) - 1, -1, -1):
        ch = name[l]
        if ch in "XY":
            n = 1
        elif ch == "Z":
            n = 2 if tau % 2 == 0 else 0
        else:
            n = 0 if tau % 2 == 0 else 2
        d += n
        tau += n
    return d


def majorana_indices(name):
    """The Majorana index set S of a Pauli string (O ~ prod_{k in S} a_k),
    left-JW, computed right-to-left via the tail parity tau: with tau odd the
    site displays Z times its own content, so X reads as a_2l+1, Y as a_2l,
    Z as no Majorana and I as both."""
    tau = 0
    S = []
    for l in range(len(name) - 1, -1, -1):
        ch = name[l]
        even = (tau % 2 == 0)
        if ch == "X":
            S.append(2 * l if even else 2 * l + 1)
            tau += 1
        elif ch == "Y":
            S.append(2 * l + 1 if even else 2 * l)
            tau += 1
        elif ch == "Z":
            if even:
                S.extend([2 * l, 2 * l + 1])
                tau += 2
        else:
            if not even:
                S.extend([2 * l, 2 * l + 1])
                tau += 2
    return S


def bi_degree(name):
    """(d_E, d_O): Majorana count per component of the h = 0 hopping graph.
    XX+YY gives the Majorana edges (2l+1, 2l+2) and (2l, 2l+3); walking them
    from any index stays inside E = {k : k mod 4 in {0, 3}} respectively
    O = {k : k mod 4 in {1, 2}}; only the field edges (2l, 2l+1) bridge."""
    S = majorana_indices(name)
    dE = sum(1 for k in S if k % 4 in (0, 3))
    return dE, len(S) - dE


def block_readable(N, w, support):
    for p, q in support:
        if (p + w - q) % 2:
            continue
        a = (p + w - q) // 2
        if 0 <= a <= w and a <= p and p - a <= N - w:
            return True
    return False


class Config:
    """One census configuration with h as a free parameter."""

    def __init__(self, N, seed, coh_kind, gammas, label):
        self.N, self.seed, self.coh_kind, self.gammas = N, seed, coh_kind, gammas
        self.label = label
        D = 2 ** N
        self.D = D

        def op(single, site):
            ops = [s0] * N
            ops[site] = single
            out = ops[0]
            for o in ops[1:]:
                out = np.kron(out, o)
            return out

        self.op = op

        def bit(i, l):
            return (i >> (N - 1 - l)) & 1

        self.bit = bit

        self.H0 = np.zeros((D, D), dtype=complex)
        for b in range(N - 1):
            self.H0 += op(sx, b) @ op(sx, b + 1) + op(sy, b) @ op(sy, b + 1)
        self.Zdiag = [np.array([1 - 2 * bit(i, l) for i in range(D)], dtype=float)
                      for l in range(N)]

        self.mask = np.zeros((D, D))
        for i in range(D):
            for j in range(D):
                self.mask[i, j] = -2 * sum(g * (bit(i, l) ^ bit(j, l))
                                           for l, g in enumerate(gammas))

        rb = np.zeros((D, D), dtype=complex)
        rb[seed, seed] = 0.5
        rb[D - 1 - seed, D - 1 - seed] = 0.5
        if coh_kind == "real":
            rb[seed, D - 1 - seed] = 0.3
            rb[D - 1 - seed, seed] = 0.3
        self.rho0 = rb
        ps, pt = bin(seed).count("1"), N - bin(seed).count("1")
        self.sup_pop = {(ps, ps), (pt, pt)}
        self.sup_coh = {(ps, pt), (pt, ps)}
        self.coh = coh_kind is not None

        self.readouts = []
        for letters in product("IXYZ", repeat=N):
            name = "".join(letters)
            if name == "I" * N:
                continue
            O = op(PAULI[letters[0]], 0)
            for site in range(1, N):
                O = O @ op(PAULI[letters[site]], site)
            self.readouts.append((name, O))

        # F132's predicted alive set (h-independent by construction)
        self.law_alive = set()
        for name, _ in self.readouts:
            w = sum(1 for ch in name if ch in "XY")
            d = majorana_degree(name)
            pred = (block_readable(N, w, self.sup_pop) and d % 4 == 0) or \
                   (self.coh and block_readable(N, w, self.sup_coh) and d == N)
            if pred:
                self.law_alive.add(name)

        self._mag_cache = {}

    def hamiltonian(self, hs):
        H = self.H0.copy()
        for l in range(self.N):
            H += hs[l] * np.diag(self.Zdiag[l]).astype(complex)
        return H

    def propagate(self, rho, hs, dt=0.05, ticks=40):
        H = self.hamiltonian(hs)
        mask = self.mask

        def rhs(r):
            return -1j * (H @ r - r @ H) + mask * r

        states = [rho.copy()]
        for _ in range(ticks):
            k1 = rhs(rho)
            k2 = rhs(rho + dt / 2 * k1)
            k3 = rhs(rho + dt / 2 * k2)
            k4 = rhs(rho + dt * k3)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            states.append(rho.copy())
        return states

    def alive_magnitudes(self, hs):
        """{name: max_t |<O>(t)|} at fixed h, cached per h tuple."""
        key = tuple(np.round(hs, 12))
        if key not in self._mag_cache:
            states = self.propagate(self.rho0.copy(), hs)
            mags = {}
            for name, O in self.readouts:
                mags[name] = max(abs(float(np.real(np.trace(r @ O))))
                                 for r in states)
            self._mag_cache[key] = mags
        return self._mag_cache[key]

    def alive_set(self, hs):
        return {n for n, m in self.alive_magnitudes(hs).items() if m > ALIVE_TOL}


def sign_rule_alive(cfg):
    """The bare chiral-conj per-channel sign rule at h = 0 (T3):
    eps'_g = (-1)^(n_Y + xy_g); population channel needs eps'_g = +1 for both
    gauges, the coherence channel needs eps'_g to match the sign U_g puts on
    the coherence pair; the K and d conditions are F132's, unchanged."""
    N, D = cfg.N, cfg.D
    Ug = {}
    for gname, par in (("even", 0), ("odd", 1)):
        u = np.ones(D)
        for l in range(N):
            if l % 2 == par:
                u = u * cfg.Zdiag[l]
        Ug[gname] = u
    s_pair = {g: float(u[cfg.seed] * u[D - 1 - cfg.seed]) for g, u in Ug.items()}
    pred = set()
    for name, _ in cfg.readouts:
        w = sum(1 for ch in name if ch in "XY")
        d = majorana_degree(name)
        ny = sum(1 for ch in name if ch == "Y")
        eps = {}
        for gname, par in (("even", 0), ("odd", 1)):
            xy_g = sum(1 for l, ch in enumerate(name) if ch in "XY" and l % 2 == par)
            eps[gname] = (-1) ** (ny + xy_g)
        pop_ok = (block_readable(N, w, cfg.sup_pop) and d % 4 == 0
                  and eps["even"] == 1 and eps["odd"] == 1)
        coh_ok = (cfg.coh and block_readable(N, w, cfg.sup_coh) and d == N
                  and eps["even"] == s_pair["even"] and eps["odd"] == s_pair["odd"])
        if pop_ok or coh_ok:
            pred.add(name)
    return pred


def populated_bisectors(cfg):
    """The (d_E, d_O) sectors carrying a nonzero Pauli coefficient of rho(0),
    including the identity's (0, 0)."""
    populated = {(0, 0)}
    for name, O in cfg.readouts:
        if abs(np.trace(cfg.rho0 @ O)) > 1e-12:
            populated.add(bi_degree(name))
    return populated


def collapsed_h0_alive(cfg):
    """The collapsed h = 0 rule: alive iff bi-sector populated and K-readable
    (through either populated channel). No sign condition, no mod-4 condition:
    both are automatic per populated sector."""
    populated = populated_bisectors(cfg)
    pred = set()
    for name, _ in cfg.readouts:
        w = sum(1 for ch in name if ch in "XY")
        k_any = block_readable(cfg.N, w, cfg.sup_pop) or \
            (cfg.coh and block_readable(cfg.N, w, cfg.sup_coh))
        if k_any and bi_degree(name) in populated:
            pred.add(name)
    return pred


def main():
    configs = [
        Config(3, 1, None, [0.0, 0.0, 0.5], "N=3 s=1 population"),
        Config(3, 1, "real", [0.0, 0.0, 0.5], "N=3 s=1 real coherence"),
        Config(4, 1, None, [0.0] * 3 + [0.5], "N=4 s=1 population"),
        Config(4, 1, "real", [0.0] * 3 + [0.5], "N=4 s=1 real coherence"),
        Config(5, 3, "real", [0.0] * 4 + [0.5], "N=5 s=3 real coherence"),
    ]
    cfg5s1 = Config(5, 1, "real", [0.0] * 4 + [0.5], "N=5 s=1 real coherence")

    print("T1  generic sufficiency: the law is exact at every random h sample")
    for cfg in configs:
        rng = np.random.default_rng(20260716 + cfg.N * 100 + cfg.seed)
        exact = True
        floor, ceil = np.inf, 0.0
        for _ in range(6):
            hs = rng.uniform(-1.0, 1.0, cfg.N) * rng.choice([0.3, 0.7, 1.3])
            mags = cfg.alive_magnitudes(hs)
            alive = {n for n, m in mags.items() if m > ALIVE_TOL}
            exact = exact and (alive == cfg.law_alive)
            floor = min(floor, min(mags[n] for n in cfg.law_alive))
            ceil = max(ceil, max(mags[n] for n in set(mags) - cfg.law_alive))
        check(f"{cfg.label}: 6 random h samples",
              exact and floor > 1e-3 and ceil < 1e-12,
              f"law exact at every sample = {exact}; alive floor {floor:.2e}, "
              f"dead ceiling {ceil:.2e} (machine zero, not near-threshold)")

    print("T2  the structured points: staggered and single-site are generic; "
          "h = 0 and uniform h shrink the alive set by exact counts")
    extra_h0 = {}
    expect_h0 = {"N=3 s=1 population": 6, "N=3 s=1 real coherence": 16,
                 "N=4 s=1 population": 24, "N=4 s=1 real coherence": 32,
                 "N=5 s=3 real coherence": 256}
    expect_uni = {"N=3 s=1 population": 6, "N=3 s=1 real coherence": 6,
                  "N=4 s=1 population": 24, "N=4 s=1 real coherence": 0,
                  "N=5 s=3 real coherence": 130}
    for cfg in configs:
        N = cfg.N
        ok_generic = True
        for hs in ([0.7 * (-1) ** l for l in range(N)], [0.7] + [0.0] * (N - 1)):
            ok_generic = ok_generic and (cfg.alive_set(hs) == cfg.law_alive)
        a0 = cfg.alive_set([0.0] * N)
        au = cfg.alive_set([0.7] * N)
        extra_h0[cfg.label] = cfg.law_alive - a0
        check(f"{cfg.label}: structured points",
              ok_generic and not (a0 - cfg.law_alive) and not (au - cfg.law_alive)
              and len(cfg.law_alive - a0) == expect_h0[cfg.label]
              and len(cfg.law_alive - au) == expect_uni[cfg.label],
              f"staggered/single-site generic = {ok_generic}; necessity holds at "
              f"h=0 and uniform; extra deaths h=0: {len(cfg.law_alive - a0)} "
              f"(expect {expect_h0[cfg.label]}), uniform: "
              f"{len(cfg.law_alive - au)} (expect {expect_uni[cfg.label]})")

    print("T3  the bare chiral-conj symmetry at h = 0: V'_g = U_g conj(.) U_g, "
          "eps'_g = (-1)^(n_Y + xy_g)")
    # dynamic pin: V' commutes with the h = 0 flow (apply V' to rho0, evolve
    # both, compare), on a NON-V'-invariant coherence prep too
    worst = 0.0
    for cfg in (configs[0], configs[1]):
        for par in (0, 1):
            U = np.ones(cfg.D)
            for l in range(cfg.N):
                if l % 2 == par:
                    U = U * cfg.Zdiag[l]
            Vrho0 = (U[:, None] * cfg.rho0.conj()) * U[None, :]
            path_a = cfg.propagate(Vrho0.copy(), [0.0] * cfg.N)
            path_b = cfg.propagate(cfg.rho0.copy(), [0.0] * cfg.N)
            worst = max(worst, max(float(np.max(np.abs(
                a - (U[:, None] * b.conj()) * U[None, :])))
                for a, b in zip(path_a, path_b)))
    check("V' commutes with the h = 0 flow (dynamic, both gauges, both preps)",
          worst < 1e-12, f"worst residual V'(rho(t)) vs evolve(V'(rho0)): {worst:.1e}")
    for cfg in configs[:4]:
        pred = sign_rule_alive(cfg)
        a0 = cfg.alive_set([0.0] * cfg.N)
        check(f"{cfg.label}: sign rule exact at h = 0",
              pred == a0, f"predicts {len(pred)}, actual {len(a0)}")
    cfg = configs[4]
    pred = sign_rule_alive(cfg)
    a0 = cfg.alive_set([0.0] * cfg.N)
    over = pred - a0
    populated = populated_bisectors(cfg)
    check("N=5 s=3: the sign rule over-predicts by exactly the bi-degree kills",
          not (a0 - pred) and len(over) == 10
          and all(bi_degree(n) not in populated for n in over),
          f"over-prediction {len(over)} strings (never under), all with "
          f"unpopulated (d_E, d_O): {sorted(over)}")

    print("T4  the bi-degree conservation law and the collapsed h = 0 rule")
    all6 = configs + [cfg5s1]
    import random
    rnd = random.Random(7)
    worst_leak = 0.0
    for cfg in all6:
        sample = cfg.readouts if cfg.N == 3 else rnd.sample(cfg.readouts, 40)
        for name, O in sample:
            C = cfg.H0 @ O - O @ cfg.H0
            bd = bi_degree(name)
            for name2, O2 in cfg.readouts:
                coef = np.trace(C @ O2) / cfg.D
                if abs(coef) > 1e-12 and bi_degree(name2) != bd:
                    worst_leak = max(worst_leak, abs(coef))
    check("[H0, O] never leaves O's (d_E, d_O) sector",
          worst_leak == 0.0,
          f"worst leak {worst_leak:.1e} (all 4^N strings at N=3, 40 sampled "
          f"per config above)")
    for cfg in all6:
        pred = collapsed_h0_alive(cfg)
        a0 = cfg.alive_set([0.0] * cfg.N)
        populated = populated_bisectors(cfg)
        pop_ok = all(dE == dO and dE % 2 == 0 for dE, dO in populated
                     if dE + dO != cfg.N) if not cfg.coh else True
        check(f"{cfg.label}: collapsed h = 0 rule exact",
              pred == a0 and pop_ok,
              f"bi-sector populated AND K-readable predicts {len(pred)}, actual "
              f"{len(a0)}; populated sectors {sorted(populated)}")

    print("T5  the shadow identity: the bare V' is the mod-4 shadow of the "
          "bi-degree, one level below F132's eps_g(d)")
    bad = 0
    total = 0
    for N in range(2, 7):
        for letters in product("IXYZ", repeat=N):
            name = "".join(letters)
            d = majorana_degree(name)
            dE, dO = bi_degree(name)
            ny = sum(1 for ch in name if ch == "Y")
            for par, dg in ((0, dE), (1, dO)):
                xy_g = sum(1 for l, ch in enumerate(name)
                           if ch in "XY" and l % 2 == par)
                total += 1
                if (-1) ** (ny + xy_g) != (-1) ** ((d * (d - 1)) // 2 + dg):
                    bad += 1
    check("eps'_even = (-1)^(d(d-1)/2 + d_E), eps'_odd = (-1)^(d(d-1)/2 + d_O)",
          bad == 0, f"{total} (string, gauge) pairs N=2..6, {bad} mismatches")

    print("T6  uniform h: population channel identical to h = 0, coherence "
          "channel revived by the rotating block phase")
    for cfg in configs:
        a0 = cfg.alive_set([0.0] * cfg.N)
        au = cfg.alive_set([0.7] * cfg.N)
        coh_channel = set()
        for name, _ in cfg.readouts:
            w = sum(1 for ch in name if ch in "XY")
            if cfg.coh and block_readable(cfg.N, w, cfg.sup_coh) \
                    and majorana_degree(name) == cfg.N:
                coh_channel.add(name)
        check(f"{cfg.label}: alive(uniform) = alive(h=0) union coherence channel",
              au == a0 | coh_channel,
              f"h=0 {len(a0)} + coherence channel {len(coh_channel)} "
              f"(overlap {len(a0 & coh_channel)}) = {len(a0 | coh_channel)}, "
              f"actual {len(au)}")

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: F132's sufficiency is generic in h (machine-zero dead "
          "sets at every random sample); at h = 0 the Majorana hopping graph "
          "disconnects, the bi-degree (d_E, d_O) is conserved, and the dead-set "
          "rule collapses to one kinematic line (bi-sector populated and "
          "K-readable); the bare chiral-conj symmetry is its mod-4 shadow, one "
          "level below F132's; uniform h keeps the h = 0 population channel "
          "and revives the coherence channel.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
