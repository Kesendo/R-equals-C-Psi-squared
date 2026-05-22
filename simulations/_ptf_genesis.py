"""PTF on the genesis system: does the closure law survive the on-site clocks?

docs/carbon master-question thread (2026-05-22, Tom + Claude).

PTF (Perspectival Time Field, hypotheses/PERSPECTIVAL_TIME_FIELD.md): under a
local J-defect, each site's single-qubit purity is a time-rescaling of the
unperturbed chain, P_B(i,t) ~ P_A(i, alpha_i * t); the per-site alpha_i satisfy
a closure law Sum_i ln(alpha_i) ~ 0. Standard PTF (n7_coupling_defect_overlay.py
+ observer_time_rescale.py) used a pure XY chain, no on-site terms.

The genesis arc concluded the qubits are SOURCES, each carrying an on-site clock
(its level splitting). This re-runs the PTF analysis WITH the on-site clocks in:
H = Sum_l h_l Z_l + Sum_bonds (J/2)(XX+YY). Two cases run side by side:
  h = 0   : standard PTF, a machinery check against the published alpha_i.
  h != 0  : the genesis version, distinct incommensurate on-site fields.
The question: does the closure Sum_i ln(alpha_i) ~ 0 survive the on-site clocks?

Machinery reused from n7_coupling_defect_overlay.py (propagation) and
observer_time_rescale.py (the alpha fit). Investigation only.
"""
import sys

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

sys.stdout.reconfigure(encoding="utf-8")

N = 7
GAMMA0 = 0.05
T_FINAL = 80.0
N_STEPS = 400
DT = T_FINAL / N_STEPS
T_FIT = 20.0
J_UNIFORM = 1.0
J_MOD = 1.1                        # delta J = +0.1 on bond (0,1)
ALPHA_BOUNDS = (0.1, 10.0)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
_LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    factors = [I2] * N
    factors[site] = op
    return kron_chain(*factors)


def onsite_fields(N, scale):
    """Distinct, incommensurate per-qubit on-site fields (the qubits' clocks)."""
    primes = np.array([2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0])[:N]
    return scale * np.sqrt(primes)


def build_H(J_list, h_field, N):
    """H = sum_l h_l Z_l + sum_i (J_i/2)(X_i X_{i+1} + Y_i Y_{i+1})."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for l in range(N):
        if h_field[l] != 0.0:
            H += h_field[l] * site_op(Z, l, N)
    for i in range(N - 1):
        H += (J_list[i] / 2.0) * (site_op(X, i, N) @ site_op(X, i + 1, N) +
                                  site_op(Y, i, N) @ site_op(Y, i + 1, N))
    return H


def build_hamming_matrix(N):
    d = 2 ** N
    idx = np.arange(d, dtype=np.uint32)
    xor = idx[:, None] ^ idx[None, :]
    h = np.zeros((d, d), dtype=np.int32)
    for i in range(N):
        h += ((xor >> i) & 1).astype(np.int32)
    return h


def lindblad_rhs(rho, H, hamming, gamma_0):
    commutator = H @ rho - rho @ H
    dephasing = -2.0 * gamma_0 * hamming * rho
    return -1j * commutator + dephasing


def rk4_step(rho, H, hamming, gamma_0, dt):
    k1 = lindblad_rhs(rho, H, hamming, gamma_0)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, hamming, gamma_0)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, hamming, gamma_0)
    k4 = lindblad_rhs(rho + dt * k3, H, hamming, gamma_0)
    return rho + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def single_excitation_mode(N, k=1):
    psi = np.zeros(2 ** N, dtype=complex)
    norm = np.sqrt(2.0 / (N + 1))
    for i in range(N):
        psi[2 ** (N - 1 - i)] = norm * np.sin(np.pi * k * (i + 1) / (N + 1))
    return psi


def initial_rho_bonding(N):
    """rho_0 = |phi><phi|, phi = (|vac> + |psi_1>)/sqrt(2)."""
    phi = np.zeros(2 ** N, dtype=complex)
    phi[0] = 1.0
    phi = phi + single_excitation_mode(N, 1)
    phi /= np.linalg.norm(phi)
    return np.outer(phi, phi.conj())


def reduced_single(rho, i, N):
    row = [_LETTERS[q] for q in range(N)]
    col = [_LETTERS[q + N] for q in range(N)]
    for q in range(N):
        if q != i:
            col[q] = row[q]
    spec = ''.join(row) + ''.join(col) + '->' + row[i] + col[i]
    return np.einsum(spec, rho.reshape([2] * (2 * N)))


def site_purity_all(rho, N):
    out = np.empty(N)
    for i in range(N):
        r = reduced_single(rho, i, N)
        out[i] = float(np.real(np.trace(r @ r)))
    return out


def propagate(J_list, h_field, hamming):
    H = build_H(J_list, h_field, N)
    rho = initial_rho_bonding(N)
    purity = np.empty((N_STEPS + 1, N))
    purity[0] = site_purity_all(rho, N)
    for step in range(1, N_STEPS + 1):
        rho = rk4_step(rho, H, hamming, GAMMA0, DT)
        purity[step] = site_purity_all(rho, N)
    return purity


def alpha_fit(t, p_A_i, p_B_i, t_max=T_FIT):
    """One-parameter time-rescale fit (observer_time_rescale.py convention)."""
    interp = interp1d(t, p_A_i, bounds_error=False,
                      fill_value=(float(p_A_i[0]), float(p_A_i[-1])),
                      kind='cubic')
    mask = t <= t_max
    t_eval = t[mask]
    b = p_B_i[mask]

    def mse(alpha):
        d = interp(alpha * t_eval) - b
        return float(np.mean(d * d))

    res = minimize_scalar(mse, bounds=ALPHA_BOUNDS, method='bounded',
                          options={'xatol': 1e-6})
    alpha = float(res.x)
    rmse = float(np.sqrt(res.fun))
    lo, hi = ALPHA_BOUNDS
    boundary = bool(abs(alpha - lo) < 1e-4 or abs(alpha - hi) < 1e-4)
    return alpha, rmse, boundary


def run_ptf(h_field, label):
    hamming = build_hamming_matrix(N)
    times = np.arange(N_STEPS + 1) * DT
    J_A = [J_UNIFORM] * (N - 1)
    J_B = [J_MOD] + [J_UNIFORM] * (N - 2)
    p_A = propagate(J_A, h_field, hamming)
    p_B = propagate(J_B, h_field, hamming)

    alphas, rmses, bds = [], [], []
    for i in range(N):
        a, rmse, bd = alpha_fit(times, p_A[:, i], p_B[:, i])
        alphas.append(a)
        rmses.append(rmse)
        bds.append(bd)
    closure = float(np.sum(np.log(np.clip(alphas, 1e-9, None))))

    print(f"--- {label} ---")
    print(f"  h_l    = [{', '.join(f'{h:.3f}' for h in h_field)}]")
    print(f"  alpha_i= [{', '.join(f'{a:.3f}' for a in alphas)}]")
    print(f"  RMSE   = [{', '.join(f'{r:.4f}' for r in rmses)}]"
          f"{'   BOUNDARY HIT' if any(bds) else ''}")
    print(f"  Sum ln(alpha_i) = {closure:+.4f}   "
          f"({'closed' if abs(closure) < 0.2 else 'broken'})")
    print()
    return alphas, closure


if __name__ == "__main__":
    print(f"PTF on the genesis system: N={N}, gamma_0={GAMMA0}, "
          f"J-defect delta J = {J_MOD - J_UNIFORM:+.1f} on bond (0,1)\n")
    run_ptf(np.zeros(N), "h = 0   standard PTF (machinery check)")
    run_ptf(onsite_fields(N, 0.15),
            "h != 0   genesis: qubits with on-site clocks")
