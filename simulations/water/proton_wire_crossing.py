"""The water-wire crossing pass: what the framework's grading IS on a proton wire.

The carbon pass (docs/carbon/BENZENE_THREE_DEPHASE_LETTERS.md) filtered its
Hamiltonians by the substrate's conserved quantity, pi-electron number
N_hat = sum_l (I - Z_l)/2. The carbon qubit is the OCCUPATION of a pi site, so
that operator IS the electron count.

The water qubit is different: it is the POSITION of one proton inside its own
hydrogen bond (|L> donor / |R> acceptor, docs/water/HYDROGEN_BOND_QUBIT.md:40-42).
The proton is always there, so the same operator cannot be a particle count. This
script asks what it is instead, and the answer is exact:

    popcount = the wire's electric dipole moment, in units of e * (O-O spacing).

Set s_l = 1 when the proton in bond l sits on its right-hand oxygen. A wire of N
bonds has N+1 oxygens; against the neutral reference s = 0...0 the wire charges are

    q_0 = -s_0,   q_j = s_{j-1} - s_j  (0 < j < N),   q_N = +s_{N-1}

which sum to zero in every configuration, so mu = sum_j j * q_j is origin-free.
It equals sum_l s_l identically. The framework's grading therefore survives the
crossing into water and changes its physical meaning: particle number in carbon,
transported charge in water.

Scope: this models a single-file water WIRE (carbon nanotube, gramicidin A,
aquaporin), not bulk water. Oxygen coordination is 4 in ice Ih and ~3.5 in the
liquid; a 1D chain of 2-coordinate oxygens is the confined case only.

Run:  python simulations/water/proton_wire_crossing.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
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
# primitives
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


def bit_string(i: int, N: int) -> list[int]:
    """s_0 .. s_{N-1} for basis index i, site 0 = most significant bit."""
    return [(i >> (N - 1 - l)) & 1 for l in range(N)]


def wire_charges(s: list[int]) -> list[int]:
    """q_0 .. q_N on the N+1 oxygens, against the neutral reference s = 0...0."""
    N = len(s)
    q = [0] * (N + 1)
    q[0] = -s[0]
    for j in range(1, N):
        q[j] = s[j - 1] - s[j]
    q[N] = s[N - 1]
    return q


# --------------------------------------------------------------------------
# the Hamiltonian inventory of docs/water/
# --------------------------------------------------------------------------


def H_heisenberg(N: int, J: float = 1.0) -> np.ndarray:
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N - 1):
        for P in ("X", "Y", "Z"):
            H += J * two_site_op(N, a, a + 1, P, P)
    return H


def H_tunnel(N: int, J: float = 1.0) -> np.ndarray:
    """-J sum_l X_l: the proton crossing its own double well."""
    return sum(-J * site_op(N, l, "X") for l in range(N))


def H_ising(N: int, K: float = 0.5) -> np.ndarray:
    """K sum_b Z_a Z_b: electrostatics between neighbouring proton positions."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N - 1):
        H += K * two_site_op(N, a, a + 1, "Z", "Z")
    return H


def H_tfi(N: int, J: float = 1.0, K: float = 0.5) -> np.ndarray:
    """The doc's own words: 'the physical proton model'."""
    return H_tunnel(N, J) + H_ising(N, K)


def H_field(N: int, delta: float = 0.3) -> np.ndarray:
    """Delta sum_l Z_l: the double-well bias. This IS the dipole coupled to a
    uniform field along the wire, since sum_l Z_l = N*I - 2*mu."""
    return sum(delta * site_op(N, l, "Z") for l in range(N))


def magnetization(N: int) -> np.ndarray:
    d = 2**N
    return sum((np.eye(d, dtype=complex) - site_op(N, l, "Z")) / 2 for l in range(N))


def z_bath(N: int, gamma: float = 0.5) -> list[np.ndarray]:
    return [np.sqrt(gamma) * site_op(N, l, "Z") for l in range(N)]


def rate_pairing_error(evals: np.ndarray, sigma: float) -> float:
    """Is the multiset of decay RATES palindromic about -sigma?"""
    re = np.sort(evals.real)
    return float(np.max(np.abs(re + re[::-1] + 2 * sigma)))


# --------------------------------------------------------------------------
# W1: popcount is the dipole moment
# --------------------------------------------------------------------------


def w1_popcount_is_dipole() -> None:
    print("\nW1  Is popcount the wire's dipole moment?")
    for N in (3, 4, 5, 6):
        worst_mu = 0.0
        worst_neutral = 0
        for i in range(2**N):
            s = bit_string(i, N)
            q = wire_charges(s)
            worst_neutral = max(worst_neutral, abs(sum(q)))
            mu = sum(j * q[j] for j in range(N + 1))
            worst_mu = max(worst_mu, abs(mu - sum(s)))
        check(
            f"N={N} the wire is neutral in every configuration",
            worst_neutral == 0,
        )
        check(
            f"N={N} mu = sum_j j*q_j equals popcount on every basis state",
            worst_mu == 0,
            f"max|mu - popcount| = {worst_mu}",
        )
    # and the field term is the dipole, up to a constant
    N = 4
    d = 2**N
    mu_op = magnetization(N)
    field = sum(site_op(N, l, "Z") for l in range(N))
    check(
        "sum_l Z_l = N*I - 2*mu (the bias field IS the dipole coupling)",
        float(np.linalg.norm(field - (N * np.eye(d, dtype=complex) - 2 * mu_op))) < 1e-12,
    )

    # zero net transport is not 'nothing moves': an XX+YY swap moves two protons
    # in opposite directions and leaves the dipole where it was.
    before, after = [1, 1, 0, 0], [1, 0, 1, 0]
    dq = [b - a for a, b in zip(wire_charges(before), wire_charges(after))]
    check(
        "an XX+YY swap moves protons but not the dipole",
        dq == [0, 1, -2, 1, 0] and sum(before) == sum(after),
        f"delta q = {dq}, popcount {sum(before)} -> {sum(after)}",
    )


# --------------------------------------------------------------------------
# W2: the domain-wall count is NOT the charged-oxygen count
# --------------------------------------------------------------------------


def w2_domain_walls_are_not_the_charge() -> None:
    print("\nW2  Is the domain-wall count the number of charged oxygens?")
    for N in (4, 5):
        disagree = 0
        for i in range(2**N):
            s = bit_string(i, N)
            walls = sum(1 for l in range(N - 1) if s[l] != s[l + 1])
            charged = sum(1 for x in wire_charges(s) if x != 0)
            if walls != charged:
                disagree += 1
        check(
            f"N={N} domain walls disagree with the charged-oxygen count",
            disagree > 0,
            f"{disagree}/{2**N} basis states",
        )
    # the sharp counterexample: a proton pushed the whole length of the wire
    N = 4
    s = [1, 1, 1, 1]
    q = wire_charges(s)
    walls = sum(1 for l in range(N - 1) if s[l] != s[l + 1])
    check(
        "s = 1111 has charged termini but zero domain walls",
        walls == 0 and q[0] == -1 and q[N] == 1,
        f"walls = {walls}, q = {q} (OH- at one end, H3O+ at the other)",
    )
    # and it is blind to the sign
    s_a, s_b = [1, 0, 0, 0], [0, 0, 0, 1]
    check(
        "the domain-wall count cannot tell H3O+ from OH-",
        sum(1 for l in range(3) if s_a[l] != s_a[l + 1])
        == sum(1 for l in range(3) if s_b[l] != s_b[l + 1])
        and wire_charges(s_a)[0] != wire_charges(s_b)[0],
        f"{wire_charges(s_a)} vs {wire_charges(s_b)}, both 1 wall",
    )


# --------------------------------------------------------------------------
# W3: the currency table
# --------------------------------------------------------------------------


def w3_currency_table() -> None:
    print("\nW3  Which operators move the dipole?")
    print("    ||[mu, H]||_F / ||H||_F   (0 = the dipole is conserved)\n")
    for N in (4, 5):
        mu = magnetization(N)
        inventory = {
            "Heisenberg J(XX+YY+ZZ)": H_heisenberg(N),
            "tunneling -J sum X": H_tunnel(N),
            "Ising K sum ZZ": H_ising(N),
            "TFI (physical proton model)": H_tfi(N),
            "bias Delta sum Z": H_field(N),
        }
        for name, H in inventory.items():
            r = float(np.linalg.norm(mu @ H - H @ mu) / np.linalg.norm(H))
            print(f"    N={N}  {name:<30} {r:.6f}")
        # the ratio 1 for the tunneling term is an identity, not a maximum
        Hx = H_tunnel(N)
        Hy = sum(-site_op(N, l, "Y") for l in range(N))
        check(
            f"N={N} [mu, -J sum X] = -i J sum Y, so the ratio 1 is an identity",
            float(np.linalg.norm((mu @ Hx - Hx @ mu) - (-1j) * Hy)) < 1e-12,
        )
        # something else exceeds it, so 1.0 is not a ceiling
        H3 = sum(two_site_op(N, a, a + 1, "X", "X") for a in range(N - 1))
        r3 = float(np.linalg.norm(mu @ H3 - H3 @ mu) / np.linalg.norm(H3))
        check(
            f"N={N} the ratio is not bounded by 1 (XX chain exceeds it)",
            r3 > 1.0 + 1e-9,
            f"XX chain gives {r3:.6f}",
        )
        # the charge to the right of the cut at bond k is s_k, NOT the popcount:
        # the popcount is the sum over all N cuts, not a single-cut integral
        bad = 0
        for i in range(2**N):
            s = bit_string(i, N)
            q = wire_charges(s)
            for k in range(N):
                if sum(q[k + 1 :]) != s[k]:
                    bad += 1
        check(
            f"N={N} charge right of the cut at bond k is s_k, so the popcount is the "
            "sum over all cuts and not a single-cut reading",
            bad == 0,
        )


# --------------------------------------------------------------------------
# W4: the kernel, read as the zero-transport sector
# --------------------------------------------------------------------------


def w4_kernel() -> None:
    print("\nW4  The kernel of L: what survives, and on which model")
    gamma = 0.5
    for N in (3, 4, 5):
        c_ops = z_bath(N, gamma)
        for name, H in (
            ("Heisenberg", H_heisenberg(N)),
            ("TFI (physical)", H_tfi(N)),
        ):
            L = fw.lindbladian_general(H, c_ops)
            ev = np.linalg.eigvals(L)
            nker = int(np.sum(np.abs(ev) < 1e-9 * max(np.abs(ev).max(), 1e-15)))
            expected = N + 1 if name == "Heisenberg" else 1
            check(
                f"N={N} {name}: dim ker L = {expected}",
                nker == expected,
                f"got {nker}",
            )


# --------------------------------------------------------------------------
# W5: what survives the crossing -- the F1 palindrome
# --------------------------------------------------------------------------


def w5_palindrome() -> None:
    print("\nW5  The F1 palindrome across the inventory")
    gamma = 0.5
    for N in (3, 4, 5):
        c_ops = z_bath(N, gamma)
        sigma = N * gamma
        for name, H in (
            ("tunneling -J sum X", H_tunnel(N)),
            ("Ising K sum ZZ", H_ising(N)),
            ("TFI (physical proton model)", H_tfi(N)),
        ):
            L = fw.lindbladian_general(H, c_ops)
            m = float(np.linalg.norm(fw.palindrome_residual(L, sigma, N)))
            check(f"N={N} {name}: palindrome bit-exact", m < 1e-10, f"||M|| = {m:.2e}")

    # the bias breaks it, linearly, and no profile repairs it
    N, gamma = 4, 0.5
    c_ops, sigma = z_bath(N, gamma), N * gamma
    base = H_tfi(N)
    norms = {}
    for delta in (0.0, 0.1, 0.3, 1.0):
        L = fw.lindbladian_general(base + H_field(N, delta), c_ops)
        norms[delta] = float(np.linalg.norm(fw.palindrome_residual(L, sigma, N)))
    print("    bias sweep: " + ", ".join(f"D={d} -> {v:.6f}" for d, v in norms.items()))
    check("Delta = 0 restores the palindrome", norms[0.0] < 1e-10)
    check(
        "||M|| is exactly linear in Delta",
        abs(norms[0.3] / norms[0.1] - 3.0) < 1e-9,
        f"ratio {norms[0.3] / norms[0.1]:.9f}",
    )

    # per-site: the map delta -> M is injective, so no bias profile cancels the break.
    # Injectivity alone carries that; the equal singular values say the extra thing,
    # that the map is a scaled isometry, so ||M|| depends only on ||delta||.
    for n in (3, 4, 5):
        cn, sn = z_bath(n, gamma), n * gamma
        basen = H_tfi(n)
        M0 = fw.palindrome_residual(fw.lindbladian_general(basen, cn), sn, n)
        cols = [
            np.asarray(
                fw.palindrome_residual(
                    fw.lindbladian_general(basen + site_op(n, l, "Z"), cn), sn, n
                )
                - M0
            ).flatten()
            for l in range(n)
        ]
        sv = np.linalg.svd(np.stack(cols, axis=1), compute_uv=False)
        check(
            f"N={n} the bias -> M map is injective, so no bias profile restores it",
            int(np.sum(sv > 1e-9 * sv[0])) == n,
            f"singular values {np.round(sv, 6).tolist()}",
        )
        check(
            f"N={n} the map is a scaled isometry with scale 2^(N+1)",
            float(sv.max() - sv.min()) < 1e-9 and abs(sv[0] - 2 ** (n + 1)) < 1e-9,
        )

    # 'only the bias breaks it' is false, and the rule is F87's letter parity:
    # a term survives iff #Y and #Z are both even, since Z-dephasing gives Y and Z
    # the same bit_b. Tested STANDALONE, so no interference with the tunneling term.
    # YZ+ZY is the discriminating case: #Y and #Z are each ODD (so not truly, and it
    # must break) while their SUM is even. A rule phrased on the bit_b parity alone
    # would wrongly let it through. It is the repo's canonical soft case.
    standalone = {
        "sum X_l": (sum(site_op(N, l, "X") for l in range(N)), True),
        "sum Y_l": (sum(site_op(N, l, "Y") for l in range(N)), False),
        "sum Z_l": (sum(site_op(N, l, "Z") for l in range(N)), False),
        "sum (XZ + ZX)": (
            sum(
                two_site_op(N, a, a + 1, "X", "Z") + two_site_op(N, a, a + 1, "Z", "X")
                for a in range(N - 1)
            ),
            False,
        ),
        "sum (YZ + ZY), #Y and #Z odd but their sum even": (
            sum(
                two_site_op(N, a, a + 1, "Y", "Z") + two_site_op(N, a, a + 1, "Z", "Y")
                for a in range(N - 1)
            ),
            False,
        ),
    }
    for name, (Hx, survives) in standalone.items():
        L = fw.lindbladian_general(Hx, c_ops)
        m = float(np.linalg.norm(fw.palindrome_residual(L, sigma, N)))
        check(
            f"F87 letter parity predicts {name} "
            + ("survives" if survives else "breaks"),
            (m < 1e-10) if survives else (m > 1e-6),
            f"||M|| = {m:.3f}",
        )


# --------------------------------------------------------------------------
# W6: ||M|| is not the severity meter
# --------------------------------------------------------------------------


def w6_rate_vs_frequency_break() -> None:
    print("\nW6  Does the bias break the RATES, or only the frequencies?")
    N, gamma, delta = 4, 0.5, 0.3
    c_ops, sigma = z_bath(N, gamma), N * gamma
    for name, base in (("Heisenberg", H_heisenberg(N)), ("TFI", H_tfi(N))):
        L = fw.lindbladian_general(base + H_field(N, delta), c_ops)
        ev = np.linalg.eigvals(L)
        err = rate_pairing_error(ev, sigma)
        m = float(np.linalg.norm(fw.palindrome_residual(L, sigma, N)))
        print(f"    {name:<12} ||M|| = {m:8.3f}   rate-pairing error = {err:.2e}")
        if name == "Heisenberg":
            check(
                "on Heisenberg the bias leaves the decay-rate palindrome exact",
                err < 1e-9,
                f"{err:.2e}",
            )
        else:
            check(
                "on the physical proton model the bias breaks the rates too",
                err > 1e-3,
                f"{err:.2e}",
            )


def main() -> None:
    print("=" * 88)
    print("THE WATER-WIRE CROSSING PASS")
    print("=" * 88)
    w1_popcount_is_dipole()
    w2_domain_walls_are_not_the_charge()
    w3_currency_table()
    w4_kernel()
    w5_palindrome()
    w6_rate_vs_frequency_break()

    n_pass = sum(1 for _, ok in CHECKS if ok)
    print("\n" + "=" * 88)
    print(f"GATE: {n_pass}/{len(CHECKS)} checks passed")
    print("=" * 88)
    if n_pass != len(CHECKS):
        for name, ok in CHECKS:
            if not ok:
                print(f"  FAILED: {name}")
        sys.exit(1)


if __name__ == "__main__":
    main()
