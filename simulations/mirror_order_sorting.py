"""Gate for the mirror's order-sorting law (docs/proofs/PROOF_MIRROR_ORDER_SORTING.md).

THE LAW (Theorem A, the unitary column; derived 2026-07-16): let R be the F71
site reversal and split the parameters of the Lindblad chain (bonds J, local
dephasing gamma, longitudinal h) into an R-even base and an R-odd scan
direction with scale t. The conjugation identity

    (R (x) R) . L(base + t*dir) . (R (x) R) = L(base - t*dir)     (exact)

forces, for an OPERATOR-R-even preparation (R rho0 R = rho0) and a readout O
of definite R-parity q (R O R = q O):

    <O>(t; time) = q * <O>(-t; time)   for every evolution time.

The four (q, sigma_eff) cells: (+,+) generic; (+,-) EVEN response; (-,-) ODD
response; (-,+) IDENTICALLY ZERO. sigma_eff is the sign mirror-conjugation
puts on the scan parameter (here = the operator parity of the direction;
chi_M = +1, R linear). The antiunitary column (chi_M = -1, the zeta^2 law,
Theorem B) is gated separately by simulations/zeta2_anti_protection.py.

GATES (exit 0 iff all PASS; default N = 4 and 5, --deep adds N = 6 ~ minutes):
  T0  the conjugation identity entry-wise (machine zero) + prep/readout parity
  T1  full Liouville spectrum EVEN in t (optimal-assignment multiset matching;
      naive complex sort mispairs the near-degenerate -2*gamma groups)
  T2  trajectory parities: <O_even> even, <O_odd> odd in t, several times,
      machine precision; odd cell verified non-vacuous
  T3  the ZERO cell: R-even generator + R-even prep + R-odd readout ->
      <O>(time) = 0 identically (THE cell witness; NOT the Pi-protected-
      observables mechanism, which is a distinct cluster cancellation)
  T4  the leak (hypothesis violation): prep = R-even + eps * (R-odd) opens an
      odd-in-t component EXACTLY affine in eps (master equation linear in
      rho0; ratio 2.000 on halving eps) with leading term O(eps*t)
      (t-halving check: the odd component's slope converges as t -> 0)

Fences: the identity requires reversing ALL scanned parameters together (a
mixed J/gamma scan is R-odd only if every component negates); the preparation
hypothesis is OPERATOR-even, expectation-even is provably insufficient; the
general vec condition is U (x) conj(U), and kron(R, R) is correct here
because R is a real permutation.

Runtime ~40 s default. Scout history: _reception_parity_scout.py (local).
"""

import sys

import numpy as np
from scipy.linalg import expm
from scipy.optimize import linear_sum_assignment

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
s0 = np.eye(2, dtype=complex)

FAILURES = []


def check(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        FAILURES.append(name)


def run_level(N):
    D = 2 ** N
    print(f"\n===== N = {N} (Liouville {D * D} x {D * D}) =====")

    def op(single, site):
        ops = [s0] * N
        ops[site] = single
        out = ops[0]
        for o in ops[1:]:
            out = np.kron(out, o)
        return out

    def liouvillian(J, gammas):
        H = np.zeros((D, D), dtype=complex)
        for b in range(N - 1):
            H += J[b] * (op(sx, b) @ op(sx, b + 1) + op(sy, b) @ op(sy, b + 1))
        Id = np.eye(D)
        L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
        for i in range(N):
            Z = op(sz, i)
            L += gammas[i] * (np.kron(Z, Z.T)
                              - 0.5 * np.kron(Z @ Z, Id)
                              - 0.5 * np.kron(Id, (Z @ Z).T))
        return L

    perm = np.zeros(D, dtype=int)
    for idx in range(D):
        bits = [(idx >> k) & 1 for k in range(N)]
        perm[idx] = sum(b << (N - 1 - k) for k, b in enumerate(bits))
    R = np.zeros((D, D))
    R[perm, np.arange(D)] = 1.0
    RR = np.kron(R, R)          # correct because R is real (general: U (x) conj U)

    NB = N - 1
    J_SYM = np.array([1.0 - 0.15 * min(b, NB - 1 - b) for b in range(NB)])
    J_ANTI = np.array([(2 * b - (NB - 1)) / max(NB - 1, 1) for b in range(NB)])
    G_SYM = np.array([0.05 + 0.03 * min(i, N - 1 - i) for i in range(N)])

    O_EVEN = op(sz, 0) + op(sz, N - 1)
    O_ODD = op(sz, 0) - op(sz, N - 1)

    psi = np.zeros(D, dtype=complex)
    psi[1 << 0] = 1.0
    psi[1 << (N - 1)] = 1.0
    psi /= np.linalg.norm(psi)
    RHO0 = np.outer(psi, psi.conj())
    psi_a = np.zeros(D, dtype=complex)
    psi_a[1 << 0] = 1.0
    psi_a[1 << (N - 1)] = -1.0
    psi_a /= np.linalg.norm(psi_a)
    RHO_ODD = 0.5 * (np.outer(psi, psi_a.conj()) + np.outer(psi_a, psi.conj()))

    def traj(J, gam, rho0, O, T):
        v = expm(liouvillian(J, gam) * T) @ rho0.reshape(-1)
        return float(np.real(np.trace(v.reshape(D, D) @ O)))

    print("T0  the conjugation identity + parities")
    t = 0.37
    res = np.max(np.abs(RR @ liouvillian(J_SYM + t * J_ANTI, G_SYM) @ RR
                        - liouvillian(J_SYM - t * J_ANTI, G_SYM)))
    check("(R x R) L(+t) (R x R) = L(-t)", res < 1e-12, f"residual {res:.2e}")
    check("prep operator-R-even", np.max(np.abs(R @ RHO0 @ R - RHO0)) < 1e-14, "rho0")
    check("readout parities q = +1 / -1",
          np.max(np.abs(R @ O_EVEN @ R - O_EVEN)) < 1e-14
          and np.max(np.abs(R @ O_ODD @ R + O_ODD)) < 1e-14, "O_even / O_odd")
    # gamma-direction corollary (F91's axis): R-odd gamma scan, same identity
    g_anti = np.array([(2 * i - (N - 1)) / max(N - 1, 1) for i in range(N)]) * 0.02
    resg = np.max(np.abs(RR @ liouvillian(J_SYM, G_SYM + g_anti) @ RR
                         - liouvillian(J_SYM, G_SYM - g_anti)))
    check("gamma-scan corollary (F91 axis)", resg < 1e-12, f"residual {resg:.2e}")

    print("T1  full spectrum even in t (multiset matching)")
    worst = 0.0
    for tv in (0.3, 0.9):
        ep = np.linalg.eigvals(liouvillian(J_SYM + tv * J_ANTI, G_SYM))
        em = np.linalg.eigvals(liouvillian(J_SYM - tv * J_ANTI, G_SYM))
        cost = np.abs(ep[:, None] - em[None, :])
        ri, ci = linear_sum_assignment(cost)
        worst = max(worst, float(cost[ri, ci].max()))
    check("spec L(+t) = spec L(-t)", worst < 1e-9, f"worst matched diff {worst:.2e}")

    print("T2  trajectory parities (cells (+,-) and (-,-))")
    worst_e, worst_o, odd_mag = 0.0, 0.0, 0.0
    for tv in (0.25, 0.6):
        for T in (0.7, 2.3):
            fe = (traj(J_SYM + tv * J_ANTI, G_SYM, RHO0, O_EVEN, T),
                  traj(J_SYM - tv * J_ANTI, G_SYM, RHO0, O_EVEN, T))
            fo = (traj(J_SYM + tv * J_ANTI, G_SYM, RHO0, O_ODD, T),
                  traj(J_SYM - tv * J_ANTI, G_SYM, RHO0, O_ODD, T))
            worst_e = max(worst_e, abs(fe[0] - fe[1]))
            worst_o = max(worst_o, abs(fo[0] + fo[1]))
            odd_mag = max(odd_mag, abs(fo[0]))
    check("<O_even> even in t", worst_e < 1e-12, f"worst {worst_e:.2e}")
    check("<O_odd> odd in t", worst_o < 1e-12, f"worst {worst_o:.2e}")
    check("odd cell non-vacuous", odd_mag > 1e-3, f"max |<O_odd>| = {odd_mag:.4f}")

    print("T3  the zero cell (-,+)")
    worst = max(abs(traj(J_SYM, G_SYM, RHO0, O_ODD, T))
                for T in (0.4, 1.1, 2.9, 5.0))
    check("<O_odd>(time) = 0 identically", worst < 1e-12, f"worst {worst:.2e}")

    print("T4  the leak: exactly affine in eps, leading O(eps*t)")
    def odd_component(eps, tv):
        rho = RHO0 + eps * RHO_ODD
        return (traj(J_SYM + tv * J_ANTI, G_SYM, rho, O_EVEN, 0.7)
                - traj(J_SYM - tv * J_ANTI, G_SYM, rho, O_EVEN, 0.7)) / (2 * tv)
    s1, s2 = odd_component(0.02, 0.05), odd_component(0.01, 0.05)
    check("leak exists", abs(s1) > 1e-6, f"slope(eps=0.02) = {s1:+.3e}")
    check("leak affine in eps (ratio 2)", abs(s1 / s2 - 2) < 1e-6,
          f"ratio {s1 / s2:.6f}")
    h1, h2 = odd_component(0.02, 0.05), odd_component(0.02, 0.025)
    check("leading term O(eps*t) (t-halving converges)",
          abs(h1 - h2) < 0.1 * abs(h1), f"slope t=0.05: {h1:+.4e}, "
          f"t=0.025: {h2:+.4e} (drift {abs(h1 - h2) / abs(h1) * 100:.1f}%)")


def main():
    levels = [4, 5] + ([6] if "--deep" in sys.argv else [])
    for N in levels:
        run_level(N)
    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print(f"GATE PASS: the order-sorting law holds from below at N = {levels}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
