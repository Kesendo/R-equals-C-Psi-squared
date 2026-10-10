"""
The runs of Dynamic Fixed Points §4 in the retired tool's own reading, and exactly
=================================================================================
experiments/DYNAMIC_FIXED_POINTS.md §4 reports runs of the retired delta_calc MCP tool (its source
is kept outside the repo and is transcribed here, not imported): simulate_dynamic_lindblad with
noise_type local and the mutual-purity bridge, J = 1, t_max = 10, dt = 0.01. The Claude Desktop log
of the chat's MCP holds the six local runs among the eight calls of 2026-02-08, 13:57 to 14:01 UTC
(Bell+ at h = 0 was called twice, with identical output; the operator-feedback call belongs to §6);
the log keeps the head and only the tail of the last array of each response, Psi, and the
transcription regenerates every value of that tail.

Outside its feedback modes the tool runs the scalar law: the rate gamma * C, with C the mutual
purity (the geometric mean of the single-site purities), in front of sigma_z on every site; each
Euler step is followed by the Hermitian part, the negative eigenvalues set to zero and the trace
renormalized; it records every tenth step, rounded to six digits, Psi = l1/(d - 1), and reports
C_psi_dynamic = C_final * Psi_final and R_inf_dynamic = C_final * Psi_final^2, the final values at
t_max: the first iterate of R = C (Psi + R)^2 from R = 0; above C Psi = 1/4 the recurrence has no
real fixed point and the iteration diverges. The 0.327 the page once called R_inf matches the
larger, unstable root of that recurrence at C near 0.917 and Psi = 0.27 (0.329 at C = 0.917, 0.327
at 0.9175); the tool never forms that root, and how the number arose is not recorded.
H = J (XX + YY + ZZ) per bond plus h sum_k X_k; 'product' is |0..0>.

At h = 0 H does nothing on these states (Bell+ and its dephased image lie in the bond's triplet,
GHZ_3's branches are degenerate eigenstates of the ring, W_3 is an eigenstate of the ring and
dephasing mixes it only with the maximally mixed state P_1/3 of the one-excitation block, |00> is an
eigenstate), and the bridge C stays at its start in every run, so the rate is fixed; at gamma 0.005
the Euler step is close to exact propagation (DOP853), at gamma 0.5 it departs. The
field-driven run of regime 2 is the clipped Euler image and is read beside its exact value.

Import-inert; prints only, and exits with status 1 if any check fails.
"""

from functools import reduce
import sys

import numpy as np
from scipy.integrate import solve_ivp

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)

# The tails the log shows of each call's Psi (2026-02-08 UTC).
LOG = {
    "Bell+, h 0": [0.30373, 0.303426, 0.303123, 0.30282, 0.302517, 0.302215, 0.301913, 0.301611],
    "GHZ_3 ring, h 0": [0.124255, 0.124069, 0.123883, 0.123697, 0.123511, 0.123326, 0.123141, 0.122957],
    "W_3 ring, h 0": [0.257663, 0.257377, 0.257091, 0.256806, 0.25652, 0.256236, 0.255951, 0.255667],
    "product, h 0": [0.0] * 13,
    "Bell+, h 0.9": [0.880716, 0.739132, 0.542666, 0.32841, 0.528985, 0.727228, 0.871794, 0.944177],
    "Bell+, gamma 0.5": [3.2e-05, 2.9e-05, 2.6e-05, 2.4e-05, 2.2e-05, 1.9e-05, 1.8e-05, 1.6e-05, 1.4e-05],
}

CLAIMS, CHECKS = [], []


def site(P, k, n):
    return reduce(np.kron, [P if i == k else I2 for i in range(n)])


def hamiltonian(n, ring, h):
    H = np.zeros((2**n, 2**n), dtype=complex)
    bonds = [(i, (i + 1) % n) for i in range(n if ring else n - 1)]
    for i, j in bonds:
        for P in (X, Y, Z):
            H += site(P, i, n) @ site(P, j, n)
    return H + h * sum(site(X, k, n) for k in range(n))


def pure(v):
    v = np.asarray(v, dtype=complex)
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


def purity(r):
    return float(np.real(np.trace(r @ r)))


def one_site(rho, n, k):
    t = rho.reshape([2] * (2 * n))
    keep = [k, n + k]
    rest = [a for a in range(2 * n) if a not in keep]
    m = 2 ** (n - 1)
    return np.einsum("abii->ab", np.transpose(t, keep + rest).reshape(2, 2, m, m))


def mutual_purity(rho, n):
    return float(np.prod([purity(one_site(rho, n, k)) for k in range(n)]) ** (1.0 / n))


def psi(r):
    return float((np.abs(r).sum() - np.abs(np.diag(r)).sum()) / (r.shape[0] - 1))


def states(name, n):
    if name == "Bell+":
        return pure([1, 0, 0, 1])
    if name == "GHZ":
        v = np.zeros(2**n)
        v[0] = v[-1] = 1
        return pure(v)
    if name == "W":
        return pure([1 if bin(i).count("1") == 1 else 0 for i in range(2**n)])
    v = np.zeros(2**n)
    v[0] = 1
    return pure(v)


def tool_run(rho0, n, H, gamma, t_max=10.0, dt=0.01):
    """simulate_dynamic_lindblad, noise_type local, transcribed: rows (C, Psi), and the final state."""
    L = [site(Z, k, n) for k in range(n)]
    rho, rows = rho0.copy(), []
    steps = int(t_max / dt)
    for step in range(steps + 1):
        if step % 10 == 0:
            rows.append((round(mutual_purity(rho, n), 6), round(psi(rho), 6)))
        if step == steps:
            break
        g = gamma * mutual_purity(rho, n)
        d = -1j * (H @ rho - rho @ H)
        for Lk in L:
            d += g * (Lk @ rho @ Lk - rho)
        rho = rho + dt * d
        rho = 0.5 * (rho + rho.conj().T)
        w, v = np.linalg.eigh(rho)
        rho = v @ np.diag(np.maximum(w, 0.0)) @ v.conj().T
        rho /= np.real(np.trace(rho))
    return rows, rho


def exact_run(rho0, n, H, gamma, t_max=10.0):
    L = [site(Z, k, n) for k in range(n)]
    d = rho0.shape[0]

    def f(_, y):
        r = y.reshape(d, d)
        g = gamma * mutual_purity(r, n)
        dr = -1j * (H @ r - r @ H)
        for Lk in L:
            dr += g * (Lk @ r @ Lk - r)
        return dr.reshape(-1)

    return solve_ivp(f, (0, t_max), rho0.reshape(-1).astype(complex), method="DOP853", rtol=1e-12,
                     atol=1e-14, dense_output=True)


def claim(label, ok, detail=""):
    CLAIMS.append((label, ok))
    print(f"  {label} {detail}  {'holds' if ok else 'FAILS'}")


def check(label, page, value, digits=3):
    ok = round(value, digits) == round(page, digits)
    CHECKS.append((label, ok))
    print(f"  {label:<60s} page {page:.{digits}f}  computed {value:.{digits + 3}f}  {'match' if ok else 'DIFFERS'}")


def main():
    runs = [
        ("Bell+, h 0", "Bell+", 2, False, 0.0, 0.005, (0.500, 0.333, 0.302, 0.151, 0.167)),
        ("GHZ_3 ring, h 0", "GHZ", 3, True, 0.0, 0.005, (0.500, 0.143, 0.123, 0.061, 0.071)),
        ("W_3 ring, h 0", "W", 3, True, 0.0, 0.005, (0.556, 0.286, 0.256, 0.142, 0.159)),
        ("product, h 0", "product", 2, False, 0.0, 0.005, (1.000, 0.000, 0.000, 0.000, 0.000)),
    ]
    print("Regime 1 (h = 0, gamma 0.005): the tool's run, and exactly:")
    for label, st, n, ring, h, g, page in runs:
        rho0, H = states(st, n), hamiltonian(n, ring, h)
        rows, _ = tool_run(rho0, n, H, g)
        tail = LOG[label]
        claim(f"{label}: the {len(tail)} logged Psi values regenerate", [p for _, p in rows][-len(tail):] == tail)
        cpsi = [c * p for c, p in rows]
        for name, pv, v in zip(("C", "Psi(0)", "Psi(10)", "C Psi final", "C Psi max"), page,
                               (rows[-1][0], rows[0][1], rows[-1][1], cpsi[-1], max(cpsi))):
            check(f"{label}: {name} (the tool)", pv, v)
        sol = exact_run(rho0, n, H, g)
        r10 = sol.sol(10.0).reshape(rho0.shape)
        claim(f"{label}: exactly, Psi at t = 10 agrees with the tool to within 2e-6",
              abs(psi(r10) - rows[-1][1]) < 2e-6, f"(|difference| {abs(psi(r10) - rows[-1][1]):.1e})")
        _, rho_end = tool_run(rho0, n, H, g)
        claim(f"{label}: the tool's C stays at its start (every record, and the unrounded final state to 1e-12)",
              all(c == rows[0][0] for c, _ in rows) and abs(mutual_purity(rho_end, n) - mutual_purity(rho0, n)) < 1e-12)

    print("\nRegime 2 (Bell+, h 0.9, gamma 0.0045): the tool's clipped Euler image, and exactly:")
    rho0, H = states("Bell+", 2), hamiltonian(2, False, 0.9)
    rows, _ = tool_run(rho0, 2, H, 0.0045)
    claim(f"the {len(LOG['Bell+, h 0.9'])} logged Psi values regenerate",
          [p for _, p in rows][-len(LOG["Bell+, h 0.9"]):] == LOG["Bell+, h 0.9"])
    c, p = rows[-1]
    check("C Psi at t = 10 (the tool's C_psi_dynamic)", 0.472, c * p)
    check("C Psi^2 at t = 10 (the tool's R_inf_dynamic)", 0.446, c * p * p)
    check("Psi over the records, smallest (the tool)", 0.33, min(q for _, q in rows), digits=2)
    check("Psi over the records, largest (the tool)", 0.99, max(q for _, q in rows), digits=2)
    sol = exact_run(rho0, 2, H, 0.0045)
    ts = np.linspace(0, 10, 101)
    ex = [(mutual_purity(sol.sol(t).reshape(4, 4), 2), psi(sol.sol(t).reshape(4, 4))) for t in ts]
    ce, pe = ex[-1]
    check("exactly, C Psi at t = 10", 0.468, ce * pe)
    check("exactly, C Psi^2 at t = 10", 0.439, ce * pe * pe)
    check("exactly, Psi over the same records, smallest", 0.32, min(q for _, q in ex), digits=2)
    check("exactly, Psi over the same records, largest", 0.99, max(q for _, q in ex), digits=2)
    fine = [psi(sol.sol(t).reshape(4, 4)) for t in np.linspace(0, 10, 20001)]
    check("exactly, Psi on a fine grid, smallest", 0.32, min(fine), digits=2)
    check("exactly, Psi on a fine grid, largest", 0.997, max(fine), digits=3)
    claim("the bridge C stays 1/2 along the exact run and at every tool record (the uniform field with the "
          "computed; Bell+'s symmetries keep its single sites maximally mixed)",
          max(abs(c - 0.5) for c, _ in ex) < 1e-12 and all(c == 0.5 for c, _ in rows))
    it = [0.0]
    for _ in range(6):
        it.append(c * (p + it[-1]) ** 2)
    claim("at the tool's final C, Psi the recurrence has no real fixed point (1 - 4 C Psi < 0), so from R = 0 it "
          "diverges; its first iterate is R_inf_dynamic", 1 - 4 * c * p < 0 and round(it[1], 6) == round(c * p * p, 6),
          "(" + ", ".join(f"{v:.3f}" for v in it[1:]) + ")")
    print("  the match noted by the historical copy: the larger root of R = C (Psi + R)^2 at Psi = 0.27")
    for cc, page in ((0.917, 0.329), (0.9175, 0.327)):
        cp = cc * 0.27
        rp = (1 - 2 * cp + np.sqrt(1 - 4 * cp)) / (2 * cc)
        check(f"  at C = {cc} (the historical copy's values)", page, rp)
        slope = 2 * cc * (0.27 + rp)
        claim(f"    slope 2 C (Psi + R) = 1 + sqrt(1 - 4 C Psi) = {slope:.3f} > 1: unstable (<= 100 eps)",
              abs(slope - (1 + np.sqrt(1 - 4 * cp))) <= 100 * np.finfo(float).eps and slope > 1)

    print("\nRegime 3 (Bell+, h 0, gamma 0.5):")
    rho0, H = states("Bell+", 2), hamiltonian(2, False, 0.0)
    rows, rho = tool_run(rho0, 2, H, 0.5)
    rho_tool_end = rho
    claim(f"the {len(LOG['Bell+, gamma 0.5'])} logged Psi values regenerate",
          [p for _, p in rows][-len(LOG["Bell+, gamma 0.5"]):] == LOG["Bell+, gamma 0.5"])
    check("purity at t = 10", 0.5, purity(rho))
    check("C Psi at t = 10", 0.0, rows[-1][0] * rows[-1][1])
    sol = exact_run(rho0, 2, H, 0.5)
    pex = psi(sol.sol(10.0).reshape(4, 4))
    claim("the Euler step departs here (gamma C dt = 0.0025): Psi(10) differs from exact by more than 1 %",
          abs(psi(rho_tool_end) - pex) / pex > 0.01, f"(tool {psi(rho_tool_end):.3e}, exact {pex:.3e})")

    held = sum(ok for _, ok in CLAIMS)
    matched = sum(ok for _, ok in CHECKS)
    print(f"\n{matched} of {len(CHECKS)} printed values match at the page's precision; "
          f"{held} of {len(CLAIMS)} regenerations and claims hold.")
    if held < len(CLAIMS) or matched < len(CHECKS):
        sys.exit(1)


if __name__ == "__main__":
    main()
