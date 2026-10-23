"""One density matrix, two jump counts: the single-excitation sector under uniform Z-dephasing, unravelled twice.

Check for docs/quantum/QUANTUM_JUMP_TRANSLATED.md. The sector of F126 (docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md):
one excitation on N sites of a chain with hopping h (amplitude J), local Z-dephasing at rate γ per site, so

    ρ̇ = −i[h, ρ] + γ·Σ_k (Z_k ρ Z_k − ρ) = −i[h, ρ] − Γ·(ρ − diag ρ),   Γ = 4γ,

with Z_k = I − 2|k⟩⟨k| on the sector. Two jump sets give this one generator:

    the house's own jumps   c_k = √γ·Z_k            (N of them, Σ c†c = Nγ·I),
    the site projectors      c_n = √Γ·|n⟩⟨n|        (N of them, Σ c†c = Γ·I),

the second being the jump set whose refill order F126's ladder counts (experiments/FRONT_PEDIGREE.md). A counting field s
on the jump term, L_s = L₀ + s·𝒥, generates the jump-count distribution: the coefficient of s^j in tr e^{t·L_s}ρ₀ is the
probability of j jumps by time t, and the coefficient of s^j in its diagonal is the population carried by j-jump histories.

What is checked, exactly, over the Gaussian rationals (Taylor coefficients in t, polynomials in s, no float anywhere):
(1) the two jump sets give the same generator at s = 1, entry for entry;
(2) both no-jump parts are a number times the identity (Σ c†c ∝ I), so the no-jump evolution is the clean unitary times
    one scalar decay; the control: a single site projector, Σ c†c = Γ|0⟩⟨0|, is not;
(3) hence the count is Poisson whatever the state: tr e^{t·L_s}ρ₀ = e^{(s−1)·c·t} with c = Nγ for the Z jumps and c = Γ
    for the projectors, read on two different initial states (a site, a coherent superposition), order by order in t;
(4) the populations agree (the same ρ), and the j-resolved pieces are related exactly: on the sector
    γ·Σ_k Z_k ρ Z_k = (N − 4)γ·ρ + Γ·diag ρ, so L_s(Z) = L_s(projectors) + (s − 1)(N − 4)γ and
    e^{t·L_s(Z)} = e^{(s−1)(N−4)γt}·e^{t·L_s(projectors)}: the two counts differ by an independent Poisson number of jumps
    the ensemble does not register, the pieces coincide at N = 4 and differ at every other N;
(4') the record: each Z_k jump fires at rate γ·tr(Z_k ρ Z_k) = γ·tr ρ whatever ρ, so the whole record of the house's own
    jumps is state-independent, while a projector catch at site n fires at Γ·ρ_nn (the control: it depends on ρ);
(5) the control again: with the single projector the count depends on the initial state, the kind of no-jump term that
    lets a record announce a jump in advance.

Run: python simulations/quantum_jump_unravelings.py   (a few seconds)
"""
from __future__ import annotations

import sys

import sympy as sp

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

RESULTS = []


def check(name, ok):
    RESULTS.append((name, bool(ok)))
    print(f"{'PASS' if ok else 'FAIL'}  {name}", flush=True)


s = sp.symbols('s')
I = sp.I
ORDER = 6                       # Taylor order in t
GAMMA = sp.Rational(1, 3)       # γ per site
J = 1                           # hopping amplitude


def hopping(n):
    return sp.Matrix(n, n, lambda a, b: J if abs(a - b) == 1 else 0)


def projector(n, k):
    return sp.Matrix(n, n, lambda a, b: 1 if a == b == k else 0)


def zk(n, k):
    return sp.eye(n) - 2 * projector(n, k)


def generator(n, jumps, field):
    """ρ ↦ −i[h, ρ] + Σ_c (field·c ρ c† − ½{c†c, ρ}), the counting field on the jump term."""
    h = hopping(n)
    k = sum((c.H * c for c in jumps), sp.zeros(n, n))

    def apply(rho):
        out = -I * (h * rho - rho * h) - (k * rho + rho * k) / 2
        for c in jumps:
            out += field * c * rho * c.H
        return out.applyfunc(sp.expand)
    return apply, k


def taylor(apply, rho0):
    """The coefficients X_k = L^k ρ₀ / k! of e^{tL}ρ₀ up to ORDER."""
    xs, x = [rho0], rho0
    for k in range(1, ORDER + 1):
        x = (apply(x) / k).applyfunc(sp.expand)
        xs.append(x)
    return xs


def z_jumps(n):
    return [sp.sqrt(GAMMA) * zk(n, k) for k in range(n)]


def projector_jumps(n):
    return [sp.sqrt(4 * GAMMA) * projector(n, k) for k in range(n)]


def seeds(n):
    site = projector(n, 0)
    sup = sp.zeros(n, n)                       # (|0⟩ + i|2⟩)/√2
    sup[0, 0], sup[2, 2], sup[0, 2], sup[2, 0] = sp.Rational(1, 2), sp.Rational(1, 2), -I / 2, I / 2
    return {'site 0': site, '(|0⟩ + i|2⟩)/√2': sup}


def basis(n):
    for a in range(n):
        for b in range(n):
            e = sp.zeros(n, n)
            e[a, b] = 1
            yield e


for n in (3, 4, 5, 6):
    gZ, kZ = generator(n, z_jumps(n), s)
    gP, kP = generator(n, projector_jumps(n), s)
    gZ1, _ = generator(n, z_jumps(n), 1)
    gP1, _ = generator(n, projector_jumps(n), 1)
    print(f'N = {n}, γ = {GAMMA}, J = {J}')

    # (1) one generator
    check(f'(1) N = {n}: the Z jumps and the site projectors give the same generator, on every matrix unit',
          all((gZ1(e) - gP1(e)).applyfunc(sp.simplify) == sp.zeros(n, n) for e in basis(n)))
    check(f'(1) N = {n}: and that generator damps every coherence at Γ = 4γ and keeps the diagonal (a check of the build)',
          all(sp.simplify(generator(n, z_jumps(n), 1)[0](e)[a, b] + (-I * (hopping(n) * e - e * hopping(n)))[a, b] * -1
                          - (0 if a == b else -4 * GAMMA * e[a, b])) == 0
              for e in basis(n) for a in range(n) for b in range(n)))

    # (2) the no-jump parts
    check(f'(2) N = {n}: Σ c†c = Nγ·I for the Z jumps and Γ·I for the projectors',
          kZ == n * GAMMA * sp.eye(n) and kP == 4 * GAMMA * sp.eye(n))

    # (3) Poisson whatever the state
    for label, rho0 in seeds(n).items():
        for name, g, c in (('Z jumps', gZ, n * GAMMA), ('projectors', gP, 4 * GAMMA)):
            xs = taylor(g, rho0)
            ok = all(sp.expand(xs[k].trace() - (c * (s - 1)) ** k / sp.factorial(k)) == 0 for k in range(ORDER + 1))
            check(f'(3) N = {n}, {name}, seed {label}: tr e^(tL_s)ρ₀ = e^((s−1)·{c}·t) to order t^{ORDER}', ok)

    # (4) one ρ, two pedigrees
    rho0 = seeds(n)['site 0']
    xz, xp = taylor(gZ, rho0), taylor(gP, rho0)
    same = all((xz[k].subs(s, 1) - xp[k].subs(s, 1)).applyfunc(sp.expand) == sp.zeros(n, n) for k in range(ORDER + 1))
    check(f'(4) N = {n}: at s = 1 the two expansions are the same ρ(t), order by order', same)
    jump_map = all((sum((GAMMA * zk(n, k) * e * zk(n, k) for k in range(n)), sp.zeros(n, n))
                    - (n - 4) * GAMMA * e - 4 * GAMMA * sp.diag(*[e[a, a] for a in range(n)])).applyfunc(sp.expand)
                   == sp.zeros(n, n) for e in basis(n))
    check(f'(4) N = {n}: γ·Σ_k Z_k ρ Z_k = (N − 4)γ·ρ + Γ·diag ρ on every matrix unit', jump_map)
    tw = sp.symbols('tw')
    factor = sp.series(sp.exp((s - 1) * (n - 4) * GAMMA * tw), tw, 0, ORDER + 1).removeO()
    fk = [sp.expand(factor.coeff(tw, k)) for k in range(ORDER + 1)]
    related = all((xz[k] - sum((fk[i] * xp[k - i] for i in range(k + 1)), sp.zeros(n, n))).applyfunc(sp.expand)
                  == sp.zeros(n, n) for k in range(ORDER + 1))
    check(f'(4) N = {n}: e^(tL_s)ρ₀ for the Z jumps = e^((s−1)(N−4)γt) × the same for the projectors, to order t^{ORDER}',
          related)
    differ = any(sp.expand(xz[k][a, a] - xp[k][a, a]) != 0 for k in range(1, ORDER + 1) for a in range(n))
    check(f'(4) N = {n}: the populations carried by j-jump histories ' + ('coincide (N = 4)' if n == 4 else 'differ'),
          differ == (n != 4))

    # (4') the record
    rec_z = all(sp.expand((GAMMA * zk(n, k) * e * zk(n, k)).trace() - GAMMA * e.trace()) == 0
                for k in range(n) for e in basis(n))
    rec_p = any(sp.expand((4 * GAMMA * projector(n, k) * e * projector(n, k)).trace() - 4 * GAMMA * e.trace() / n) != 0
                for k in range(n) for e in basis(n))
    check(f"(4') N = {n}: every Z_k jump fires at γ·tr ρ on every matrix unit (the record carries no state), "
          f"while a projector catch does depend on ρ (the control)", rec_z and rec_p)

# (5) the control: one site projector, Σ c†c = Γ|0⟩⟨0|
n = 5
gC, kC = generator(n, [sp.sqrt(4 * GAMMA) * projector(n, 0)], s)
check('(5) the control: a single site projector has Σ c†c = Γ·|0⟩⟨0|, not a multiple of the identity',
      kC != kC[0, 0] * sp.eye(n) and kC != sp.zeros(n, n))
tr = {label: [sp.expand(x.trace()) for x in taylor(gC, rho0)] for label, rho0 in seeds(n).items()}
check('(5) the control: its jump count depends on the initial state (the two seeds differ at some order)',
      any(a != b for a, b in zip(*tr.values())))
check('(5) the control: and from site 0 it is not Poisson at rate Γ (some order up to t^6 departs from Γ^k(s − 1)^k/k!)',
      any(sp.expand(tr['site 0'][k] - (4 * GAMMA * (s - 1)) ** k / sp.factorial(k)) != 0 for k in range(ORDER + 1)))

print()
failed = [n_ for n_, ok in RESULTS if not ok]
print(f'{len(RESULTS) - len(failed)}/{len(RESULTS)} checks passed' + (f'; FAILED: {failed}' if failed else ''))
sys.exit(1 if failed else 0)
