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

What is checked, exactly, in sympy arithmetic (the jumps carry exact square roots that cancel in c ρ c†; Taylor
coefficients in t, polynomials in s, no float anywhere):
(1) the two jump sets give the same generator at s = 1, entry for entry;
(2) both sums Σ c†c are a number times the identity, so the no-jump evolution is the clean unitary times
    one scalar decay; the control: a single site projector, Σ c†c = Γ|0⟩⟨0|, is not;
(3) hence the count is Poisson whatever the state: tr e^{t·L_s}ρ₀ = e^{(s−1)·κ·t} with κ = Nγ for the Z jumps and κ = Γ
    for the projectors, read on two different initial states (a site, a coherent superposition), order by order in t;
(4) the populations agree (the same ρ), and the j-resolved pieces are related exactly: on the sector
    γ·Σ_k Z_k ρ Z_k = (N − 4)γ·ρ + Γ·diag ρ, so L_s(Z) = L_s(projectors) + (s − 1)(N − 4)γ and
    e^{t·L_s(Z)} = e^{(s−1)(N−4)γt}·e^{t·L_s(projectors)}: for N > 4 the Z-jump pieces are the projector pieces convolved
    with an independent Poisson count at rate (N − 4)γ, for N < 4 the roles swap, the pieces coincide at N = 4 and
    differ at N = 3, 5, 6;
(4') the record: each Z_k jump fires at rate γ·tr(Z_k ρ Z_k) = γ·tr ρ whatever ρ, so the whole record of the house's own
    jumps is state-independent, while a projector catch at site n fires at Γ·ρ_nn (the control: it depends on ρ);
(5) the control again: with the single projector the count depends on the initial state, the kind of no-jump term that
    lets a record announce a jump in advance;
(6) the drift identity behind the front's jump-free name, computed from the jump list c_k = √γ·U_k as γ·Σ_k (U_k − 1):
    for the house's Z jumps it is −(Γ/2)·1 on the sector. With the state-independent rates of (2) and (4'), the mean
    amplitude of the Z-kick trajectories then obeys d/dt E[ψ] = (−ih − Γ/2)·E[ψ] at every N, so |E[ψ]|² = e^{−Γt}|G|² is
    the coherent front |⟨a_n⟩|² of the renewal proof, in a reading where nothing is caught (that step is argued, not
    gated here). The name belongs to the phase convention Z_k = 1 − 2|k⟩⟨k| on the sector: the gauge control computes the
    drift of the jumps −√γ·Z_k, the same channel since c ρ c† and c†c do not see the sign, as −(2N − 2)γ, not −Γ/2;
(7) the front's number with no convention: on the full 2^N register (N = 3, 4) the vacuum-excitation coherence
    |1_n⟩⟨vac| runs with the clean wave at Γ/2, so 4|ρ_{n,vac}|² = e^{−Γt}|G_{n0}|² for the seed (|vac⟩ + |1_0⟩)/√2
    (F73's vac-SE coherence); a generator cannot see the sign of a jump operator, so no convention enters.

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

    # (4') the record, one rate function tr(c ρ c†) for both jump sets
    def rate(c, rho):
        return sp.expand((c * rho * c.H).trace())
    rec_z = all(rate(c, e) == sp.expand(GAMMA * e.trace()) for c in z_jumps(n) for e in basis(n))
    # every projector catch fires at Γ times the population of its site, on every matrix unit; a constant rate fails
    rec_p = all(rate(c, e) == sp.expand(4 * GAMMA * e[k, k])
                for k, c in enumerate(projector_jumps(n)) for e in basis(n))
    check(f"(4') N = {n}: every Z_k jump fires at γ·tr ρ on every matrix unit (with (2), the record carries no state), "
          f"while the catch at site n fires at Γ·ρ_nn (the control: it depends on ρ)", rec_z and rec_p)

    # (6) the mean amplitude of a unitary-kick unravelling c_k = √γ·U_k: between kicks the normalised state runs by −ih
    # (the no-jump part is a scalar), each U_k kicks at rate γ, so d/dt E[ψ] = (−ih + γ·Σ_k (U_k − 1))·E[ψ]. The drift
    # is computed from the jump list itself; with U_k = Z_k it is −(Γ/2) at every N, the noise-averaged amplitude of the
    # renewal proof, whose square is the front e^{−Γt}|G|²
    def drift_of(jumps):
        return sum((c / sp.sqrt(GAMMA) - sp.eye(n) for c in jumps), sp.zeros(n, n)).applyfunc(sp.expand) * GAMMA
    check(f'(6) N = {n}: the drift identity behind the mean amplitude of the Z kicks, γ·Σ_k (Z_k − 1) = −(Γ/2)·1, N-free',
          (drift_of(z_jumps(n)) - (-2 * GAMMA) * sp.eye(n)).applyfunc(sp.expand) == sp.zeros(n, n))
    # the gauge control: the jumps −√γ·Z_k give the same channel (c ρ c† and c†c do not see the sign) but, computed from
    # that jump list, the drift −(2N − 2)γ, which for N ≥ 3 is not −Γ/2: the front's second name is a phase convention
    check(f'(6) N = {n}: the gauge control: the jumps −√γ·Z_k drift the mean amplitude at −(2N − 2)γ = '
          f'{-(2 * n - 2) * GAMMA}, not −Γ/2',
          (drift_of([-c for c in z_jumps(n)]) - (-(2 * n - 2) * GAMMA) * sp.eye(n)).applyfunc(sp.expand) == sp.zeros(n, n)
          and (drift_of([-c for c in z_jumps(n)]) - (-2 * GAMMA) * sp.eye(n)).applyfunc(sp.expand) != sp.zeros(n, n))

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


# (7) the front's second name on the register: the full 2^N Lindbladian with the Pauli Z_k (the excitation is |1⟩) maps
# the vacuum-excitation coherence |1_n⟩⟨vac| to −i Σ_m h_mn |1_m⟩⟨vac| − (Γ/2)·|1_n⟩⟨vac|, for the jumps +√γ·Z_k and
# −√γ·Z_k alike (argued, not gated: c ρ c† and c†c do not see the sign): an entry of ρ that runs with the clean
# wave at Γ/2, so
# 4|ρ_{n,vac}|² = e^{−Γt}|G_{n0}|² for the seed (|vac⟩ + |1_0⟩)/√2, F73's vac-SE coherence
def kron_all(ops):
    out = ops[0]
    for o in ops[1:]:
        out = sp.kronecker_product(out, o)
    return out


PX, PY, PZ, P1 = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.Matrix([[1, 0], [0, -1]]), sp.eye(2)
for n in (3, 4):
    def site(o, k):
        return kron_all([o if j == k else P1 for j in range(n)])
    H = sum((J * (site(PX, k) * site(PX, k + 1) + site(PY, k) * site(PY, k + 1)) / 2 for k in range(n - 1)),
            sp.zeros(2 ** n, 2 ** n))
    d = 2 ** n
    vac = sp.zeros(d, 1); vac[0] = 1
    ket = [sp.zeros(d, 1) for _ in range(n)]
    for m in range(n):
        ket[m][1 << (n - 1 - m)] = 1
    h = hopping(n)
    cs = [sp.sqrt(GAMMA) * site(PZ, k) for k in range(n)]
    kk = sum((c.H * c for c in cs), sp.zeros(d, d))
    ok = True
    for m in range(n):
        E = ket[m] * vac.T
        LE = -I * (H * E - E * H) - (kk * E + E * kk) / 2 + sum((c * E * c.H for c in cs), sp.zeros(d, d))
        want = sum((-I * h[mm, m] * ket[mm] * vac.T for mm in range(n)), sp.zeros(d, d)) - 2 * GAMMA * E
        ok = ok and (LE - want).applyfunc(sp.expand) == sp.zeros(d, d)
    check(f'(7) N = {n}: on the full register the vacuum-excitation coherence runs with the clean wave at Γ/2, '
          f'so 4|ρ_n,vac|² = e^(−Γt)|G_n0|² for the seed (|vac⟩ + |1_0⟩)/√2', ok)

print()
failed = [n_ for n_, ok in RESULTS if not ok]
print(f'{len(RESULTS) - len(failed)}/{len(RESULTS)} checks passed' + (f'; FAILED: {failed}' if failed else ''))
sys.exit(1 if failed else 0)
