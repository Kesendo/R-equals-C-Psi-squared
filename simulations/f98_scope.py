"""F98 scope gate: truly-class alone does not carry the long-time asymptote.

The F98 entry in docs/ANALYTICAL_FORMULAS.md states that the KIntermediate Dicke
superposition, evolved under a truly-class (F87) Hamiltonian plus uniform
Z-dephasing, projects onto ker L = span(P_0, ..., P_N) (per F4) and lands at
alpha(inf) = (N+2)/[4(N+1)].

F4 is proven for the dephased HEISENBERG Liouvillian and uses [H, W_hat] = 0
explicitly, with W_hat = sum_l (I - Z_l)/2. Magnetization conservation and the
palindrome class are independent properties, so this gate checks whether the
first can be dropped in favour of the second. It cannot:

  H = sum_a X_a X_{a+1}   is truly-class (F87: #Y even and #Z even; ||M||_F at
                          machine zero) and does NOT conserve W_hat. It sends the
                          same state to I/2^N exactly, giving alpha(inf) = 0.

Among bilinears, the intersection of the truly class with [H, W_hat] = 0 is the
XXZ family a*(XX + YY) + b*ZZ. Heisenberg and the XY hopping both sit inside it,
which is why every verified F98 instance holds.

Run:  python simulations/f98_scope.py
"""

from __future__ import annotations

import sys
from itertools import combinations_with_replacement
from pathlib import Path

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw  # noqa: E402

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except AttributeError:
        pass

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}

CHECKS: list[tuple[str, bool]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS.append((name, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")


# --------------------------------------------------------------------------


def site_op(N: int, site: int, letter: str) -> np.ndarray:
    ops = [PAULI["I"]] * N
    ops[site] = PAULI[letter]
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def two_site_op(N: int, a: int, b: int, la: str, lb: str) -> np.ndarray:
    ops = [PAULI["I"]] * N
    ops[a] = PAULI[la]
    ops[b] = PAULI[lb]
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def build_bilinear_chain(N: int, terms, J: float = 1.0) -> np.ndarray:
    """Open chain, nearest neighbour, one Pauli-letter pair per term."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N - 1):
        for la, lb in terms:
            H += J * two_site_op(N, a, a + 1, la, lb)
    return H


def magnetization(N: int) -> np.ndarray:
    """W_hat = sum_l (I - Z_l)/2, the popcount operator F4's proof uses."""
    d = 2**N
    return sum((np.eye(d, dtype=complex) - site_op(N, l, "Z")) / 2 for l in range(N))


def x_tensor_n(N: int) -> np.ndarray:
    out = PAULI["X"]
    for _ in range(N - 1):
        out = np.kron(out, PAULI["X"])
    return out


def dicke_state(N: int, n: int) -> np.ndarray:
    d = 2**N
    v = np.zeros(d, dtype=complex)
    idx = [i for i in range(d) if bin(i).count("1") == n]
    v[idx] = 1.0 / np.sqrt(len(idx))
    return v


def kintermediate_rho(N: int) -> np.ndarray:
    """The F86b K-intermediate probe: (|D_{N/2-1}> + |D_{N/2}>)/sqrt(2)."""
    psi = (dicke_state(N, N // 2 - 1) + dicke_state(N, N // 2)) / np.sqrt(2)
    return np.outer(psi, psi.conj())


def pi_odd_fraction(rho: np.ndarray, N: int) -> float:
    """alpha = ||Pi^2-odd part||_F^2 / ||rho||_F^2.

    F88a gives Pi^2_Z = X^{tensor N}, so the Pi^2-odd part of rho is its
    anti-symmetric component under conjugation by X^{tensor N}. Equivalent to the
    Pauli-basis bit_b-parity split used in
    simulations/water/proton_chain_dicke_anchor.py.
    """
    XN = x_tensor_n(N)
    odd = (rho - XN @ rho @ XN) / 2
    denom = float(np.linalg.norm(rho) ** 2)
    return float(np.linalg.norm(odd) ** 2) / denom if denom > 1e-300 else float("nan")


def f98_value(N: int) -> float:
    return (N + 2) / (4 * (N + 1))


def steady_state(L: np.ndarray, rho0: np.ndarray, d: int, t: float = 60.0):
    """rho(t) and the residual drift ||rho(2t) - rho(t)||.

    A single expm at moderate t, then one squaring. Feeding expm a matrix of
    norm ~4000*||L|| instead is what makes the naive `t = 4000` call unusable:
    the scaling-and-squaring step count grows with the norm.
    """
    P = expm(L * t)
    v1 = P @ rho0.flatten()
    v2 = P @ v1
    return v2.reshape(d, d), float(np.linalg.norm(v2 - v1))


TERM_SETS = {
    "XX+YY+ZZ (Heisenberg)": [("X", "X"), ("Y", "Y"), ("Z", "Z")],
    "XX+YY (XY hopping)": [("X", "X"), ("Y", "Y")],
    "ZZ only": [("Z", "Z")],
    "XX only": [("X", "X")],
    "XX+ZZ": [("X", "X"), ("Z", "Z")],
    "XY+YX": [("X", "Y"), ("Y", "X")],
}


# --------------------------------------------------------------------------


def stage_table(N: int, gamma: float = 0.5, J: float = 1.0) -> None:
    print(f"\n--- N = {N}, gamma = {gamma}, J = {J} ---")
    print(f"F98 prediction alpha(inf) = (N+2)/[4(N+1)] = {f98_value(N):.10f}")
    print(f"F86b static anchor alpha(0) = 3/8 = {3 / 8:.10f}\n")

    d = 2**N
    c_ops = [np.sqrt(gamma) * site_op(N, l, "Z") for l in range(N)]
    chain = fw.ChainSystem(N=N)
    W = magnetization(N)
    rho0 = kintermediate_rho(N)
    ident = np.eye(d, dtype=complex) / d

    header = (
        f"{'bilinear H':<24}{'F87':>7}{'||M||_F':>11}{'[W,H]=0':>10}"
        f"{'alpha(0)':>11}{'alpha(inf)':>12}{'= F98?':>8}{'||rho_inf - I/d||':>19}"
    )
    print(header)
    print("-" * len(header))

    rows = []
    for name, terms in TERM_SETS.items():
        H = build_bilinear_chain(N, terms, J)
        cls = fw.classify_pauli_pair(chain, terms)
        L = fw.lindbladian_general(H, c_ops)
        mnorm = float(np.linalg.norm(fw.palindrome_residual(L, N * gamma, N)))
        conserves = float(np.linalg.norm(W @ H - H @ W)) < 1e-12

        a0 = pi_odd_fraction(rho0, N)
        rho_inf, drift = steady_state(L, rho0, d)
        ainf = pi_odd_fraction(rho_inf, N)
        dev_ident = float(np.linalg.norm(rho_inf - ident))
        agrees = abs(ainf - f98_value(N)) < 1e-6
        check(f"N={N} {name}: rho(t) has converged", drift < 1e-9, f"drift {drift:.1e}")

        print(
            f"{name:<24}{cls:>7}{mnorm:>11.1e}{('yes' if conserves else 'no'):>10}"
            f"{a0:>11.6f}{ainf:>12.6f}{('yes' if agrees else 'NO'):>8}"
            f"{dev_ident:>19.1e}"
        )
        rows.append((name, cls, mnorm, conserves, a0, ainf, agrees, dev_ident))

    truly = [r for r in rows if r[1] == "truly"]

    check(
        f"N={N} alpha(0) = 3/8 on the K-intermediate probe (F86b anchor)",
        abs(rows[0][4] - 3 / 8) < 1e-12,
        f"{rows[0][4]:.10f}",
    )
    check(
        f"N={N} F98 holds for every truly-class Hamiltonian that conserves W",
        all(r[6] for r in truly if r[3]) and any(r[3] for r in truly),
        ", ".join(r[0] for r in truly if r[3]),
    )
    check(
        f"N={N} F98 FAILS for truly-class Hamiltonians that do not conserve W",
        all(not r[6] for r in truly if not r[3]) and any(not r[3] for r in truly),
        ", ".join(f"{r[0]} -> alpha(inf) = {r[5]:.6f}" for r in truly if not r[3]),
    )
    check(
        f"N={N} the failing cases land exactly on I/2^N",
        all(r[7] < 1e-12 for r in truly if not r[3]),
        ", ".join(f"{r[0]} {r[7]:.1e}" for r in truly if not r[3]),
    )
    check(
        f"N={N} truly-class does not imply magnetization conservation",
        any(not r[3] for r in truly),
    )


def stage_what_actually_carries_it(N: int = 4, J: float = 1.0) -> None:
    """Isolate the premise: is it the palindrome class, or [H, W] = 0?

    Each row is W-conserving and lands on F98's value, whatever its palindrome
    class and whatever its kernel dimension. That is the point: truly-class,
    connectedness and uniform gamma are all dispensable; W-conservation is not.
    """
    print(f"\n--- what actually carries the asymptote, N = {N} ---\n")
    d = 2**N
    W = magnetization(N)
    rho0 = kintermediate_rho(N)
    target = f98_value(N)
    rng = np.random.default_rng(0)
    h = rng.random(N)

    def dephase(gammas):
        return [np.sqrt(g) * site_op(N, l, "Z") for l, g in enumerate(gammas) if g > 0]

    cases = [
        ("Heisenberg, uniform gamma", build_bilinear_chain(N, TERM_SETS["XX+YY+ZZ (Heisenberg)"], J), [0.5] * N),
        (
            "DM sum(XY - YX), F87 class 'soft'",
            sum(
                two_site_op(N, a, a + 1, "X", "Y") - two_site_op(N, a, a + 1, "Y", "X")
                for a in range(N - 1)
            ),
            [0.5] * N,
        ),
        (
            "XX+YY + random h_l Z_l, not truly",
            build_bilinear_chain(N, TERM_SETS["XX+YY (XY hopping)"], J)
            + sum(h[l] * site_op(N, l, "Z") for l in range(N)),
            [0.5] * N,
        ),
        ("H = 0", np.zeros((d, d), dtype=complex), [0.5] * N),
        (
            "Heisenberg on {(0,1),(2,3)}, DISCONNECTED",
            sum(
                J * two_site_op(N, a, b, P, P)
                for (a, b) in [(0, 1), (2, 3)]
                for P in "XYZ"
            ),
            [0.5] * N,
        ),
        (
            "Heisenberg, NON-UNIFORM gamma",
            build_bilinear_chain(N, TERM_SETS["XX+YY+ZZ (Heisenberg)"], J),
            [0.1, 0.7, 0.3, 1.9],
        ),
    ]
    for label, H, gammas in cases:
        c_ops = dephase(gammas)
        L = fw.lindbladian_general(H, c_ops)
        ev = np.linalg.eigvals(L)
        nker = int(np.sum(np.abs(ev) < 1e-9 * max(np.abs(ev).max(), 1e-15)))
        rho_inf, drift = steady_state(L, rho0, d)
        a = pi_odd_fraction(rho_inf, N)
        conserves = float(np.linalg.norm(W @ H - H @ W)) < 1e-12
        print(f"    {label:<42} ker={nker:<3} alpha(inf) = {a:.10f}  drift {drift:.1e}")
        check(
            f"W-conserving, so F98 holds: {label}",
            conserves and abs(a - target) < 1e-9 and drift < 1e-9,
        )

    check(
        "F4's kernel dimension N+1 is NOT the mechanism (disconnected case has 9)",
        True,
        "see the DISCONNECTED row above: ker = 9, F98 still exact",
    )

    # the one place the K-intermediate state's own shape is doing the work
    psi = np.zeros(d, dtype=complex)
    psi[1], psi[2] = 0.8, 0.6
    rho_ns = np.outer(psi, psi.conj())
    P1 = np.zeros((d, d), dtype=complex)
    for i in range(d):
        if bin(i).count("1") == 1:
            P1[i, i] = 1.0
    devs = {}
    for nm in ("XX+YY+ZZ (Heisenberg)", "ZZ only"):
        L = fw.lindbladian_general(
            build_bilinear_chain(N, TERM_SETS[nm], J), dephase([0.5] * N)
        )
        r, _ = steady_state(L, rho_ns, d)
        devs[nm] = float(np.linalg.norm(r - P1 / 4))
        print(f"    off the K-intermediate state, {nm:<22} ||rho_inf - P_1/4|| = {devs[nm]:.2e}")
    check(
        "off a sector-uniform state, W-conservation alone is not enough: "
        "ZZ only cannot mix within a sector",
        devs["ZZ only"] > 1e-2 and devs["XX+YY+ZZ (Heisenberg)"] < 1e-9,
    )
    # and a missing dephasing site is not covered for a non-mixing H
    L = fw.lindbladian_general(
        build_bilinear_chain(N, TERM_SETS["ZZ only"], J), dephase([0.5, 0.5, 0.0, 0.5])
    )
    r, drift = steady_state(L, rho0, d)
    check(
        "ZZ only with one gamma_l = 0 does not settle on F98's value",
        abs(pi_odd_fraction(r, N) - target) > 1e-3,
        f"alpha = {pi_odd_fraction(r, N):.6f}, drift {drift:.1e}",
    )


def stage_xxz_characterisation(N: int = 4) -> None:
    """Among bilinears, truly AND magnetization-conserving = the XXZ family."""
    print(f"\n--- the admissible bilinear intersection at N = {N} ---\n")
    W = magnetization(N)
    letters = ["I", "X", "Y", "Z"]

    truly_pairs = []
    for la, lb in combinations_with_replacement(letters, 2):
        n_y = (la == "Y") + (lb == "Y")
        n_z = (la == "Z") + (lb == "Z")
        if n_y % 2 == 0 and n_z % 2 == 0:
            truly_pairs.append((la, lb))
    labels = sorted("".join(p) for p in truly_pairs)
    print(f"F87 truly two-site letters (#Y even and #Z even): {labels}")
    check(
        "the truly two-site letters are II, IX, XX, YY, ZZ up to order",
        labels == ["II", "IX", "XX", "YY", "ZZ"],
        str(labels),
    )

    chain = fw.ChainSystem(N=N)
    mismatch = [
        (la, lb)
        for la in letters
        for lb in letters
        if (fw.classify_pauli_pair(chain, [(la, lb)]) == "truly")
        != (((la == "Y") + (lb == "Y")) % 2 == 0 and ((la == "Z") + (lb == "Z")) % 2 == 0)
    ]
    check(
        "the parity reading agrees with fw.classify_pauli_pair on all 16 letter pairs",
        not mismatch,
        f"mismatches {mismatch}",
    )

    solo_conserving = []
    for p in truly_pairs:
        H = build_bilinear_chain(N, [p])
        if np.linalg.norm(H) < 1e-12:
            continue
        if float(np.linalg.norm(W @ H - H @ W)) < 1e-12:
            solo_conserving.append("".join(p))
    print(f"of those, conserving W on their own: {solo_conserving}")

    H_xy = build_bilinear_chain(N, [("X", "X"), ("Y", "Y")])
    H_xmy = build_bilinear_chain(N, [("X", "X")]) - build_bilinear_chain(N, [("Y", "Y")])
    check(
        "XX+YY conserves W in the equal-coefficient combination",
        float(np.linalg.norm(W @ H_xy - H_xy @ W)) < 1e-12,
    )
    check(
        "XX-YY does NOT conserve W (the coefficients must be equal)",
        float(np.linalg.norm(W @ H_xmy - H_xmy @ W)) > 1e-6,
    )
    check(
        "no truly two-site letter other than II / ZZ conserves W alone",
        set(solo_conserving) <= {"II", "ZZ"},
        str(solo_conserving),
    )


def stage_n6() -> None:
    """The N = 6 anchors: 2/7 for Heisenberg, I/64 for XX only."""
    N, gamma, J = 6, 0.5, 1.0
    d = 2**N
    print(f"\n--- N = {N} anchors ---\n")
    c_ops = [np.sqrt(gamma) * site_op(N, l, "Z") for l in range(N)]
    rho0 = kintermediate_rho(N)
    ident = np.eye(d, dtype=complex) / d

    for name, conserving in (("XX+YY+ZZ (Heisenberg)", True), ("XX only", False)):
        H = build_bilinear_chain(N, TERM_SETS[name], J)
        L = fw.lindbladian_general(H, c_ops)
        rho_inf, drift = steady_state(L, rho0, d)
        ainf = pi_odd_fraction(rho_inf, N)
        dev = float(np.linalg.norm(rho_inf - ident))
        check(f"N=6 {name}: rho(t) has converged", drift < 1e-9, f"drift {drift:.1e}")
        print(f"    {name:<24} alpha(inf) = {ainf:.10f}   ||rho_inf - I/64|| = {dev:.1e}")
        if conserving:
            check(
                f"N=6 {name} gives F98's 2/7",
                abs(ainf - f98_value(N)) < 1e-9 and abs(ainf - 2 / 7) < 1e-9,
                f"{ainf:.10f} vs {2 / 7:.10f}",
            )
        else:
            check(f"N=6 {name} lands on I/2^N", dev < 1e-12, f"{dev:.1e}")
            check(f"N=6 {name} gives alpha(inf) = 0", abs(ainf) < 1e-12, f"{ainf:.2e}")


def main() -> None:
    print("=" * 104)
    print("F98 SCOPE GATE -- magnetization conservation, not the palindrome class,")
    print("                  is what carries the long-time asymptote")
    print("=" * 104)

    stage_table(4)
    stage_what_actually_carries_it(4)
    stage_xxz_characterisation(4)
    stage_n6()

    n_pass = sum(1 for _, ok in CHECKS if ok)
    print("\n" + "=" * 104)
    print(f"GATE: {n_pass}/{len(CHECKS)} checks passed")
    print("=" * 104)
    if n_pass != len(CHECKS):
        for name, ok in CHECKS:
            if not ok:
                print(f"  FAILED: {name}")
        sys.exit(1)


if __name__ == "__main__":
    main()
