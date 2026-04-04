"""
Beer-Lambert test for quantum cavities (v2)
============================================
The original K_DOSIMETRY.md found K_system = N × K_qubit at uniform γ.
That ratio is trivially N (arithmetic cancellation of t).

The REAL Beer-Lambert question: with non-uniform γ (sacrifice zone),
does each qubit absorb proportionally to its share γ_k/Σγ of the total
dose? Or do correlations (strong J) redistribute absorption away from
the Beer-Lambert prediction?

We measure per-qubit decoherence rates from the Liouvillian eigenvalues
and compare the ACTUAL absorption share of each qubit against the
naive prediction γ_k/Σγ. We scan J/γ from 0.01 to 100 to test whether
correlations break the proportional sharing.

Additionally: does K_qubit (per-qubit dose to reach purity target)
depend systematically on J, and how does the Schwarzschild effect
interact with non-uniform dephasing?

Output: simulations/results/k_dosimetry_beer_lambert.txt
"""

import numpy as np
from pathlib import Path
import sys, os, time

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)

def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r

def build_hamiltonian(N, J):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N; ops[i] = P; ops[i + 1] = P
            H += J * kron_chain(ops)
    return H

def build_liouvillian(N, J, gammas):
    d = 2**N
    Id = np.eye(d, dtype=complex)
    H = build_hamiltonian(N, J)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        ops = [I2] * N; ops[k] = Zm
        Lk = np.sqrt(gammas[k]) * kron_chain(ops)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T))
    return L

def purity(rho):
    return np.real(np.trace(rho @ rho))

def partial_trace_qubit(rho_full, N, keep):
    """Trace out all qubits except 'keep', return 2×2 reduced density matrix."""
    a = keep
    b = N - keep - 1
    da = 2**a
    db = 2**b
    rho_r = rho_full.reshape(da, 2, db, da, 2, db)
    rho_reduced = np.zeros((2, 2), dtype=complex)
    for i in range(da):
        for j in range(db):
            rho_reduced += rho_r[i, :, j, i, :, j]
    return rho_reduced


def subsystem_purity(rho_full, N, qubit_k):
    """Purity of the reduced state of qubit k."""
    rho_k = partial_trace_qubit(rho_full, N, qubit_k)
    return np.real(np.trace(rho_k @ rho_k))


def propagate_and_measure(N, J, gammas, rho0, times):
    """Propagate via eigendecomposition, return per-qubit purities at each time."""
    d = 2**N
    L = build_liouvillian(N, J, gammas)
    eigvals, R = np.linalg.eig(L)
    _, Lf = np.linalg.eig(L.T)
    ov = Lf.conj().T @ R
    for j in range(len(eigvals)):
        Lf[:, j] /= ov[j, j]
    coeffs = Lf.conj().T @ rho0.ravel()

    results = []
    for t in times:
        exp_l = np.exp(eigvals * t)
        rv = R @ (coeffs * exp_l)
        rho = rv.reshape(d, d)
        rho = (rho + rho.conj().T) / 2
        tr = np.trace(rho).real
        if tr > 1e-15:
            rho /= tr
        sys_pur = purity(rho)
        qubit_purs = [subsystem_purity(rho, N, k) for k in range(N)]
        results.append((t, sys_pur, qubit_purs))
    return results


def measure_decoherence_rates(N, J, gammas, rho0, t_window=5.0):
    """Measure effective decoherence rate per qubit from early-time purity decay.

    Rate_k = -d(purity_k)/dt at t→0, normalized.
    Returns list of rates per qubit.
    """
    times = np.linspace(0.01, t_window, 100)
    data = propagate_and_measure(N, J, gammas, rho0, times)

    rates = []
    for k in range(N):
        # Extract qubit-k purity time series
        purs = [d[2][k] for d in data]
        ts = [d[0] for d in data]

        # Linear fit to early part (first 20 points) of log(purity - 0.5)
        # Since single-qubit purity → 0.5 for Z-dephasing
        early = min(20, len(purs))
        p_arr = np.array(purs[:early])
        t_arr = np.array(ts[:early])

        # Rate from finite difference at earliest time
        if len(p_arr) >= 2:
            dp_dt = (p_arr[-1] - p_arr[0]) / (t_arr[-1] - t_arr[0])
            rates.append(-dp_dt)  # positive = losing purity
        else:
            rates.append(0.0)
    return rates


# ═══════════════════════════════════════════════════════════
# Main experiment
# ═══════════════════════════════════════════════════════════

out = []
def log(msg=""):
    print(msg)
    out.append(msg)

log("=" * 80)
log("BEER-LAMBERT v2: NON-UNIFORM γ — DO CORRELATIONS REDISTRIBUTE ABSORPTION?")
log("=" * 80)
log()
log("Question: With sacrifice-zone dephasing (γ_edge ≫ γ_interior),")
log("does each qubit lose purity proportionally to γ_k/Σγ?")
log("Or do Hamiltonian correlations redistribute the 'light'?")
log()

plus = np.array([1, 1], dtype=complex) / np.sqrt(2)

# ─────────────────────────────────────────────────────────
# TEST 1: Absorption share at fixed J, varying asymmetry
# ─────────────────────────────────────────────────────────

log("=" * 80)
log("TEST 1: ABSORPTION SHARE vs ASYMMETRY (N=4, J=1.0)")
log("=" * 80)
log()
log("γ profiles from uniform to extreme sacrifice zone.")
log("Beer-Lambert prediction: rate_k/Σ(rates) = γ_k/Σγ")
log()

N = 4
J = 1.0
gamma_total = 0.2  # fixed Σγ = 0.2 (= 4 × 0.05)

plus_n = plus
for _ in range(N - 1):
    plus_n = np.kron(plus_n, plus)
rho0 = np.outer(plus_n, plus_n.conj())

# Asymmetry profiles: ratio of edge γ to interior γ
asymmetries = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0]

log(f"  {'Asym':>6s} {'γ_edge':>8s} {'γ_int':>8s} | {'Share_edge':>10s} {'Pred_edge':>10s} {'Dev%':>8s} | {'Share_int':>10s} {'Pred_int':>10s} {'Dev%':>8s}")
log(f"  {'-'*6} {'-'*8} {'-'*8}   {'-'*10} {'-'*10} {'-'*8}   {'-'*10} {'-'*10} {'-'*8}")

for asym in asymmetries:
    # γ_edge = asym × γ_int, and γ_edge + 3×γ_int = gamma_total
    gamma_int = gamma_total / (asym + (N - 1))
    gamma_edge = asym * gamma_int
    gammas = [gamma_edge] + [gamma_int] * (N - 1)

    # Beer-Lambert prediction
    pred_edge = gamma_edge / gamma_total
    pred_int = gamma_int / gamma_total

    # Measure actual rates
    rates = measure_decoherence_rates(N, J, gammas, rho0)
    total_rate = sum(rates)

    if total_rate > 1e-15:
        share_edge = rates[0] / total_rate
        share_int = rates[1] / total_rate  # any interior qubit

        dev_edge = (share_edge - pred_edge) / pred_edge * 100
        dev_int = (share_int - pred_int) / pred_int * 100 if pred_int > 1e-10 else 0

        log(f"  {asym:6.0f} {gamma_edge:8.5f} {gamma_int:8.5f} | {share_edge:10.4f} {pred_edge:10.4f} {dev_edge:+7.2f}% | {share_int:10.4f} {pred_int:10.4f} {dev_int:+7.2f}%")
    else:
        log(f"  {asym:6.0f} {gamma_edge:8.5f} {gamma_int:8.5f} | (rates too small)")

log()

# ─────────────────────────────────────────────────────────
# TEST 2: Same asymmetry, scan J/γ
# ─────────────────────────────────────────────────────────

log("=" * 80)
log("TEST 2: ABSORPTION SHARE vs COUPLING STRENGTH (N=4, asym=20)")
log("=" * 80)
log()
log("Fixed sacrifice zone (γ_edge/γ_int = 20). Scan J.")
log("At J→0: qubits independent, Beer-Lambert trivially holds.")
log("At J→∞: qubits strongly entangled, what happens?")
log()

asym_fixed = 20.0
gamma_int_fixed = gamma_total / (asym_fixed + (N - 1))
gamma_edge_fixed = asym_fixed * gamma_int_fixed
gammas_fixed = [gamma_edge_fixed] + [gamma_int_fixed] * (N - 1)
pred_edge_fixed = gamma_edge_fixed / gamma_total

J_over_gamma_edge = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]

log(f"  γ_edge = {gamma_edge_fixed:.5f}, γ_int = {gamma_int_fixed:.5f}")
log(f"  Beer-Lambert: edge share = {pred_edge_fixed:.4f}")
log()
log(f"  {'J/γ_e':>8s} {'J':>8s} {'Share_edge':>10s} {'Pred':>10s} {'Dev%':>8s} {'Rate_edge':>10s} {'Rate_int':>10s}")
log(f"  {'-'*8} {'-'*8} {'-'*10} {'-'*10} {'-'*8} {'-'*10} {'-'*10}")

for jg_ratio in J_over_gamma_edge:
    J_val = jg_ratio * gamma_edge_fixed

    rates = measure_decoherence_rates(N, J_val, gammas_fixed, rho0)
    total_rate = sum(rates)

    if total_rate > 1e-15:
        share_edge = rates[0] / total_rate
        dev = (share_edge - pred_edge_fixed) / pred_edge_fixed * 100
        log(f"  {jg_ratio:8.2f} {J_val:8.5f} {share_edge:10.4f} {pred_edge_fixed:10.4f} {dev:+7.2f}% {rates[0]:10.6f} {rates[1]:10.6f}")
    else:
        log(f"  {jg_ratio:8.2f} {J_val:8.5f} (rates too small)")

log()

# ─────────────────────────────────────────────────────────
# TEST 3: Per-qubit purity evolution (visualize redistribution)
# ─────────────────────────────────────────────────────────

log("=" * 80)
log("TEST 3: PER-QUBIT PURITY TRAJECTORIES (N=4, asym=20, J=1.0)")
log("=" * 80)
log()
log("Track how each qubit's purity decays over time.")
log("Beer-Lambert predicts: edge qubit decays fastest, proportional to γ_k.")
log()

J_test = 1.0
times_traj = np.linspace(0.1, 30.0, 50)
data_traj = propagate_and_measure(N, J_test, gammas_fixed, rho0, times_traj)

log(f"  {'t':>6s} {'Sys_pur':>8s} {'Q0(edge)':>9s} {'Q1(int)':>9s} {'Q2(int)':>9s} {'Q3(int)':>9s}")
log(f"  {'-'*6} {'-'*8} {'-'*9} {'-'*9} {'-'*9} {'-'*9}")

for i in range(0, len(data_traj), 5):  # every 5th point
    t, sp, qp = data_traj[i]
    log(f"  {t:6.1f} {sp:8.4f} {qp[0]:9.4f} {qp[1]:9.4f} {qp[2]:9.4f} {qp[3]:9.4f}")

log()

# Check: at what time does each qubit reach 0.5 + epsilon?
log("  Per-qubit time to reach purity 0.55:")
for k in range(N):
    for i in range(len(data_traj)-1):
        if data_traj[i][2][k] >= 0.55 and data_traj[i+1][2][k] < 0.55:
            t_cross = data_traj[i][0]
            label = "EDGE" if k == 0 else "interior"
            gamma_k = gammas_fixed[k]
            log(f"    Q{k} ({label}, γ={gamma_k:.5f}): t ≈ {t_cross:.1f}, K_k = γ_k × t = {gamma_k * t_cross:.4f}")
            break
    else:
        label = "EDGE" if k == 0 else "interior"
        log(f"    Q{k} ({label}): not reached in window")

log()

# ─────────────────────────────────────────────────────────
# TEST 4: N-scaling of the deviation
# ─────────────────────────────────────────────────────────

log("=" * 80)
log("TEST 4: DOES DEVIATION GROW WITH N? (asym=20, J=1.0)")
log("=" * 80)
log()

for N_test in range(2, 6):
    gamma_int_t = gamma_total / (asym_fixed + (N_test - 1))
    gamma_edge_t = asym_fixed * gamma_int_t
    gammas_t = [gamma_edge_t] + [gamma_int_t] * (N_test - 1)
    pred_edge_t = gamma_edge_t / sum(gammas_t)

    plus_t = plus
    for _ in range(N_test - 1):
        plus_t = np.kron(plus_t, plus)
    rho0_t = np.outer(plus_t, plus_t.conj())

    rates_t = measure_decoherence_rates(N_test, J_test, gammas_t, rho0_t)
    total_rate_t = sum(rates_t)

    if total_rate_t > 1e-15:
        share_edge_t = rates_t[0] / total_rate_t
        dev_t = (share_edge_t - pred_edge_t) / pred_edge_t * 100
        log(f"  N={N_test}: edge share = {share_edge_t:.4f}, predicted = {pred_edge_t:.4f}, deviation = {dev_t:+.2f}%")
    else:
        log(f"  N={N_test}: rates too small")

log()

# ─────────────────────────────────────────────────────────
# TEST 5: Absorption Theorem cross-check
# ─────────────────────────────────────────────────────────

log("=" * 80)
log("TEST 5: ABSORPTION THEOREM — Re(λ) = −2γ⟨n_XY⟩")
log("=" * 80)
log()
log("The Absorption Theorem says Re(λ) depends on γ × ⟨n_XY⟩.")
log("For non-uniform γ, does Σ_k γ_k × ⟨n_XY⟩_k give the system rate?")
log()

N_at = 4
for J_at in [0.1, 1.0, 5.0]:
    L = build_liouvillian(N_at, J_at, gammas_fixed)
    eigvals = np.linalg.eigvals(L)
    # Slowest decaying mode (smallest |Re(λ)| excluding zero)
    real_parts = np.sort(np.abs(eigvals.real))
    nonzero = real_parts[real_parts > 1e-10]
    if len(nonzero) > 0:
        slowest = nonzero[0]
        fastest = nonzero[-1]
        center = np.sum(gammas_fixed)
        log(f"  J={J_at:.1f}: slowest |Re(λ)| = {slowest:.6f}, fastest = {fastest:.6f}, Σγ = {center:.4f}")

log()

# ─────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────

log("=" * 80)
log("SUMMARY")
log("=" * 80)
log()
log("CONTEXT: The original K_DOSIMETRY.md found K_system/K_qubit = N exactly.")
log("This is arithmetically trivial for uniform γ (ratio = Σγ/γ = N).")
log()
log("THE REAL QUESTION: With non-uniform γ, does each qubit absorb")
log("proportional to its γ_k? (Beer-Lambert for correlated absorbers)")
log()
log("If deviations exist, they quantify how much Hamiltonian correlations")
log("redistribute absorption away from the naive local prediction.")
log("This would be a genuine correction to Beer-Lambert from entanglement.")

out_path = RESULTS_DIR / "k_dosimetry_beer_lambert.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(out))
log(f"\n>>> Results saved to: {out_path}")
