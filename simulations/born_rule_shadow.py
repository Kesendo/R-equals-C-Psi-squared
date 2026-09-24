"""
The Born-rule shadow: a slow/fast spectral split of rho at the CΨ = 1/4 fold
============================================================================
Splits the density matrix of a Z-dephased Heisenberg chain at the rate Σγ into
a slow part (the modes that decay slower than Σγ) and a fast part (all others,
the modes that sit exactly on Σγ included) and reads two things off the split:
the Born probabilities P(i) = <i|rho|i>, which are linear in rho and so carry
no cross term, and the purity Tr(rho^2), which is quadratic and can carry one.
Standard Born probabilities are assumed, not derived.

The cut is decided on rate levels. At N=2 (J above γ/2) ten of the sixteen
modes sit exactly on Σγ: the Hamiltonian anticommutes with the dissipator
centred at Σγ there (F48), which puts every eigenvalue off the real axis on
the line Re λ = -Σγ, and the four real ones among the ten (XI+IX, YI+IY,
XZ+ZX, YZ+ZY) are swap-symmetric, so the Hamiltonian cannot touch them, and dephasing drains depth 1 at
exactly Σγ. The script prints how far the eigensolver places the ten from
the line and how far the nearest mode off it sits; every band inside that gap gives
the same split. |++> is Hamiltonian-dead (the Heisenberg bond is 2·SWAP - 1 and
uniform dephasing keeps the state swap-symmetric), so its split has a closed
form, which is checked.

CΨ = Tr(rho^2) * L1 / (d - 1), with L1 the sum of the off-diagonal moduli.

Import-inert; writes only through main(output_path).
Output: simulations/results/born_rule_shadow.txt
"""

import argparse
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm, schur, solve_sylvester

DEFAULT_OUTPUT = Path(__file__).parent / "results" / "born_rule_shadow.txt"
J = 1.0
GAMMA = 0.05
CUT_BAND = 1e-9       # any band inside the printed gap gives the same split
SEED = 20260404
N_RANDOM = 50

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
ZERO = np.array([1, 0], dtype=complex)
ONE = np.array([0, 1], dtype=complex)
PLUS = np.array([1, 1], dtype=complex) / np.sqrt(2)


def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_liouvillian(N, gammas):
    d = 2**N
    Id = np.eye(d, dtype=complex)
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N; ops[i] = P; ops[i + 1] = P
            H += J * kron_chain(ops)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        ops = [I2] * N; ops[k] = Zm
        Lk = np.sqrt(gammas[k]) * kron_chain(ops)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T))
    return L


def cpsi(rho):
    d = rho.shape[0]
    purity = np.real(np.trace(rho @ rho))
    l1 = np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))
    return purity * l1 / (d - 1)


def rate_levels(L, sigma):
    """Rates -Re(lambda), the modes on the cut, and the gap that makes the cut band-free."""
    rates = -np.linalg.eigvals(L).real
    on_cut = np.abs(rates - sigma) <= CUT_BAND
    off = rates[~on_cut]
    residual = np.max(np.abs(rates[on_cut] - sigma)) if on_cut.any() else 0.0
    gap = np.min(np.abs(off - sigma))
    slow = int(np.sum((rates < sigma) & ~on_cut))
    return slow, int(on_cut.sum()), len(rates) - slow - int(on_cut.sum()), residual, gap


def spectral_projector(L, select):
    """Spectral projector onto the modes select(λ) picks, along all others."""
    T, Q, k = schur(L, output="complex", sort=select)
    n = T.shape[0]
    X = solve_sylvester(T[:k, :k], -T[k:, k:], -T[:k, k:])
    S = np.eye(n, dtype=complex); S[:k, k:] = X
    S_inv = np.eye(n, dtype=complex); S_inv[:k, k:] = -X
    E = np.zeros((n, n), dtype=complex); E[:k, :k] = np.eye(k)
    return Q @ S @ E @ S_inv @ Q.conj().T


def slow_projector(L, sigma):
    """The modes strictly slower than the cut; the modes on it go to the fast part."""
    return spectral_projector(L, lambda z: (-z.real) < sigma - CUT_BAND)


def on_cut_projector(L, sigma):
    return spectral_projector(L, lambda z: abs(-z.real - sigma) <= CUT_BAND)


def split(L, P_slow, rho0, t):
    d = rho0.shape[0]
    v = expm(L * t) @ rho0.ravel()
    vp = P_slow @ v
    rho = v.reshape(d, d); rho = (rho + rho.conj().T) / 2
    rp = vp.reshape(d, d); rp = (rp + rp.conj().T) / 2
    rf = rho - rp
    pp = np.real(np.trace(rp @ rp))
    pf = np.real(np.trace(rf @ rf))
    pc = 2 * np.real(np.trace(rp @ rf))
    return rho, rp, rf, pp, pf, pc


def first_downward_quarter(L, rho0, t_max=60.0, steps=12000):
    """First time CΨ falls through 1/4 from above, refined by bisection; None if it never does."""
    d = rho0.shape[0]
    dt = t_max / steps
    U = expm(L * dt)
    v = rho0.ravel().copy()
    c_prev = cpsi(rho0)
    for k in range(1, steps + 1):
        v = U @ v
        c = cpsi(v.reshape(d, d))
        if c_prev > 0.25 >= c:
            lo, hi = (k - 1) * dt, k * dt
            for _ in range(60):
                mid = (lo + hi) / 2
                if cpsi((expm(L * mid) @ rho0.ravel()).reshape(d, d)) > 0.25:
                    lo = mid
                else:
                    hi = mid
            return hi
        c_prev = c
    return None


def product_rho(kets):
    psi = kron_chain(kets)
    return np.outer(psi, psi.conj())


def main(output_path=DEFAULT_OUTPUT):
    destination = Path(output_path)
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    out = []

    def log(msg=""):
        print(msg)
        out.append(msg)

    log("=" * 75)
    log("THE BORN-RULE SHADOW: SLOW/FAST SPLIT AT THE CΨ = 1/4 FOLD")
    log(f"Heisenberg chain, J = {J}, uniform Z-dephasing γ = {GAMMA} (Lindblad book)")
    log("=" * 75)
    log()

    # ─── N = 2: the rate levels and the cut ───
    N = 2
    d = 2**N
    sigma = N * GAMMA
    L = build_liouvillian(N, [GAMMA] * N)
    slow, on_cut, fast, residual, gap = rate_levels(L, sigma)
    log(f"N={N}: {d * d} modes, Σγ = {sigma}")
    log(f"  slower than Σγ: {slow}   on Σγ: {on_cut}   faster: {fast}")
    log(f"  modes on the cut sit within {residual:.1e} of Σγ; the nearest mode off it is {gap:.4f} away")
    log(f"  the modes on the cut go to the fast part (|Re λ| ≥ Σγ)")
    log()
    P_slow = slow_projector(L, sigma)

    # ─── |++> at its fold ───
    rho0 = product_rho([PLUS, PLUS])
    t_fold = first_downward_quarter(L, rho0)
    rho, rp, rf, pp, pf, pc = split(L, P_slow, rho0, t_fold)
    total = pp + pf + pc
    log("=" * 75)
    log("|++> AT ITS CΨ = 1/4 FOLD")
    log("=" * 75)
    log(f"  fold time t = {t_fold:.4f} (K = γt = {GAMMA * t_fold:.5f}), CΨ = {cpsi(rho):.6f}")
    log(f"  ||[H, rho(t)]|| = {np.linalg.norm(build_liouvillian(N, [0.0] * N) @ rho.ravel()):.1e}")
    log()
    log("  Born probabilities, slow + fast (linear in rho, no cross term):")
    log(f"  {'basis':>6s} {'P(i)':>10s} {'slow':>10s} {'fast':>10s}")
    for i, label in enumerate(["|00>", "|01>", "|10>", "|11>"]):
        log(f"  {label:>6s} {rho[i, i].real:10.6f} {rp[i, i].real:10.6f} {rf[i, i].real:10.6f}")
    log()
    log("  Purity, slow + fast + cross:")
    log(f"    Tr(rho_slow²)        = {pp:.6f}  ({pp / total * 100:.2f}%)")
    log(f"    Tr(rho_fast²)        = {pf:.6f}  ({pf / total * 100:.2f}%)")
    log(f"    2 Tr(rho_slow·rho_fast) = {pc:+.1e}")
    log(f"    total                = {total:.6f}")
    e = np.exp(-4 * GAMMA * t_fold)
    closed = ((1 + e) / 2) ** 2
    log(f"  closed form (pure dephasing of |++>): purity ((1 + e^(-4γt))/2)² = {closed:.6f},"
        f" slow part I/4 with purity 1/4, cross term 0")
    log(f"    deviations: purity {total - closed:+.1e}, slow {pp - 0.25:+.1e}, rho_slow - I/4 {np.linalg.norm(rp - np.eye(d) / 4):.1e}")
    vc = (on_cut_projector(L, sigma) @ (expm(L * t_fold) @ rho0.ravel())).reshape(d, d)
    vc = (vc + vc.conj().T) / 2
    log(f"  of the fast part's purity, the modes on the cut hold {np.real(np.trace(vc @ vc)) / pf * 100:.1f}%")
    log()

    # ─── |01> at its fold: a state the Hamiltonian moves ───
    r01 = product_rho([ZERO, ONE])
    t01 = first_downward_quarter(L, r01)
    rho, rp, rf, pp, pf, pc = split(L, P_slow, r01, t01)
    log(f"|01> at its fold t = {t01:.4f} (CΨ first rises through 1/4, then falls):")
    log(f"  {'basis':>6s} {'P(i)':>10s} {'slow':>10s} {'fast':>10s}")
    for i, label in enumerate(["|00>", "|01>", "|10>", "|11>"]):
        log(f"  {label:>6s} {rho[i, i].real:10.6f} {rp[i, i].real:10.6f} {rf[i, i].real:+10.6f}")
    log(f"  purity cross term {pc:+.1e}")
    log()

    # ─── late-time image: C_i = d * P(i) ───
    log("=" * 75)
    log("LATE-TIME IMAGE, t = 50: C_i = d × P(i)")
    log("=" * 75)
    psi_bell = np.zeros(d, dtype=complex); psi_bell[0] = psi_bell[3] = 1 / np.sqrt(2)
    states = {
        "Bell+": np.outer(psi_bell, psi_bell.conj()),
        "|01>": product_rho([ZERO, ONE]),
        "|++>": product_rho([PLUS, PLUS]),
    }
    U50 = expm(L * 50.0)
    for name, r0 in states.items():
        r = (U50 @ r0.ravel()).reshape(d, d)
        C = [d * r[i, i].real for i in range(d)]
        log(f"  {name:>6s}: " + "  ".join(f"C({lab}) = {c:.2f}" for lab, c in zip(["|00>", "|01>", "|10>", "|11>"], C)))
    log()

    # ─── N = 3: where a cross term can live ───
    N3 = 3
    d3 = 2**N3
    sigma3 = N3 * GAMMA
    L3 = build_liouvillian(N3, [GAMMA] * N3)
    slow3, on_cut3, fast3, residual3, gap3 = rate_levels(L3, sigma3)
    P3 = slow_projector(L3, sigma3)
    log("=" * 75)
    log("N = 3: WHERE A CROSS TERM CAN LIVE")
    log("=" * 75)
    log(f"  slower than Σγ: {slow3}   on Σγ: {on_cut3}   faster: {fast3}; nearest mode to the cut {gap3:.4f} away")
    for name, kets in (("|+++>", [PLUS, PLUS, PLUS]), ("|++0>", [PLUS, PLUS, ZERO])):
        r0 = product_rho(kets)
        tf = first_downward_quarter(L3, r0)
        rho, rp, rf, pp, pf, pc = split(L3, P3, r0, tf)
        tot = pp + pf + pc
        log(f"  {name} at its fold t = {tf:.4f}: slow {pp / tot * 100:.2f}%  fast {pf / tot * 100:.2f}%"
            f"  cross {pc / tot * 100:+.4f}% of the purity {tot:.6f}")
    rng = np.random.default_rng(SEED)
    shares = []
    for _ in range(N_RANDOM):
        v = rng.normal(size=d3) + 1j * rng.normal(size=d3)
        v /= np.linalg.norm(v)
        _, _, _, pp, pf, pc = split(L3, P3, np.outer(v, v.conj()), 3.0)
        shares.append(abs(pc) / (pp + pf + pc))
    shares = np.array(shares) * 100
    log(f"  {N_RANDOM} random pure states (seed {SEED}) at t = 3: |cross| share median {np.median(shares):.3f}%,"
        f" max {np.max(shares):.3f}%")
    log()

    log("=" * 75)
    log("SUMMARY")
    log("=" * 75)
    log("1. The Born probabilities split into slow + fast with no cross term: the diagonal of a")
    log("   sum is the sum of the diagonals. For |++> the fast part is pure coherence, so the whole")
    log("   image comes from the slow part, the steady state I/4; for |01> the fast part carries")
    log("   part of the image.")
    log("2. The purity can carry a cross term only where slow and fast modes overlap. At N = 2")
    log("   they never do: every mode off the real axis sits on Σγ (F48), the real ones sit at 0,")
    log("   Σγ and 2Σγ, the slow side is the steady states alone, and those are orthogonal to the")
    log("   rest, so the cross term is 0 for every state. At N = 3 |+++> (Hamiltonian-dead) still")
    log("   has none, |++0> carries a small one, and random states a fraction of a percent.")
    log("3. Photography and shutter language is a reading of this split, not a measurement mechanism.")

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    print(f"\n>>> Results saved to: {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()
    main(args.output)
