"""
Dwell-time prefactor generalization: verify the formula

    prefactor = (4/k) * (W_0 + W_k) / (W_0 + 3*W_k)

on three test cases under Z-dephasing.

Verification A: GHZ_3, GHZ_4 born below the fold (CPsi(0) < 1/4)
Verification B: Bell+ smoke test (regression against known values)
Verification C: W_3 out-of-sample test (the real test)

Task: TASK_DWELL_PREFACTOR_GENERALIZATION_V1.md
Date: 2026-04-06
"""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
import sympy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT_DIR, exist_ok=True)

lines = []


def pr(s=""):
    print(s)
    lines.append(s)


# ============================================================
# Utility functions (self-contained, no repo imports)
# ============================================================

# Pauli matrices
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS_1Q = [I2, X, Y, Z]
PAULI_LABELS_1Q = ["I", "X", "Y", "Z"]


def pauli_basis(n):
    """Return list of (label, matrix, xy_weight) for all 4^n Pauli strings."""
    if n == 1:
        weights = [0, 1, 1, 0]
        return list(zip(PAULI_LABELS_1Q, PAULIS_1Q, weights))
    sub = pauli_basis(n - 1)
    result = []
    for sl, sm, sw in sub:
        for j, (pl, pm, pw) in enumerate(zip(PAULI_LABELS_1Q, PAULIS_1Q, [0, 1, 1, 0])):
            result.append((sl + pl, np.kron(sm, pm), sw + pw))
    return result


def build_lindblad_superop(n, gamma=1.0):
    """Build the Z-dephasing Lindbladian superoperator for n qubits.

    L[rho] = sum_k gamma * (Z_k rho Z_k - rho)
    No Hamiltonian (pure dephasing).

    Returns d^2 x d^2 matrix acting on vec(rho).
    """
    d = 2**n
    d2 = d * d
    L_super = np.zeros((d2, d2), dtype=complex)
    eye_d = np.eye(d, dtype=complex)

    for k in range(n):
        # Build Z_k = I^{k} tensor Z tensor I^{n-k-1}
        Zk = np.eye(1, dtype=complex)
        for j in range(n):
            Zk = np.kron(Zk, Z if j == k else I2)

        # Dissipator: gamma * (Z_k rho Z_k - rho)
        # In superoperator form: gamma * (Z_k tensor Z_k^T - I tensor I)
        # vec(A rho B) = (B^T tensor A) vec(rho)
        L_super += gamma * (np.kron(Zk.T, Zk) - np.eye(d2, dtype=complex))

    return L_super


def propagate(L_super, rho0, t_array):
    """Propagate rho0 under superoperator L for each time in t_array."""
    d2 = L_super.shape[0]
    d = int(np.sqrt(d2))
    vec0 = rho0.flatten()
    states = []
    for t in t_array:
        vec_t = expm(L_super * t) @ vec0
        rho_t = vec_t.reshape(d, d)
        states.append(rho_t)
    return states


def propagate_at(L_super, rho0, t):
    """Propagate to a single time point."""
    vec0 = rho0.flatten()
    vec_t = expm(L_super * t) @ vec0
    d = int(np.sqrt(len(vec_t)))
    return vec_t.reshape(d, d)


def purity(rho):
    return np.real(np.trace(rho @ rho))


def l1_coherence(rho):
    d = rho.shape[0]
    return np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho)))


def psi_norm(rho):
    d = rho.shape[0]
    return l1_coherence(rho) / (d - 1)


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def sector_weights(rho, n):
    """Compute sector weights W_k for k=0..n from Pauli decomposition."""
    d = 2**n
    basis = pauli_basis(n)
    weights_by_k = {}
    for label, P, w in basis:
        a_P = np.trace(rho @ P)
        w_contrib = np.abs(a_P)**2 / d
        weights_by_k[w] = weights_by_k.get(w, 0.0) + np.real(w_contrib)
    return weights_by_k


def pauli_decomposition(rho, n):
    """Full Pauli decomposition: returns dict label -> a_P."""
    basis = pauli_basis(n)
    coeffs = {}
    for label, P, w in basis:
        coeffs[label] = np.trace(rho @ P)
    return coeffs


# ============================================================
# State constructors
# ============================================================

def ghz_state(n):
    """GHZ_n = (|00...0> + |11...1>) / sqrt(2)."""
    d = 2**n
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1.0 / np.sqrt(2)
    psi[d - 1] = 1.0 / np.sqrt(2)
    return np.outer(psi, psi.conj())


def bell_plus():
    """Bell+ = (|00> + |11>) / sqrt(2)."""
    return ghz_state(2)


def w_state_3():
    """W_3 = (|100> + |010> + |001>) / sqrt(3)."""
    psi = np.zeros(8, dtype=complex)
    psi[4] = 1.0 / np.sqrt(3)  # |100>
    psi[2] = 1.0 / np.sqrt(3)  # |010>
    psi[1] = 1.0 / np.sqrt(3)  # |001>
    return np.outer(psi, psi.conj())


# ============================================================
# VERIFICATION A: GHZ born below the fold
# ============================================================

pr("=" * 72)
pr("VERIFICATION A: GHZ_N born below the fold (N=3, N=4)")
pr("=" * 72)
pr()

results_a = {}
all_pass = True

for N in [3, 4]:
    pr(f"--- GHZ_{N} ---")
    rho0 = ghz_state(N)
    d = 2**N

    # Analytical CPsi(0) = 1/(2^N - 1)
    cpsi_analytical = 1.0 / (d - 1)
    cpsi_measured = cpsi(rho0)

    pr(f"  CPsi(0) analytical: {cpsi_analytical:.10f}")
    pr(f"  CPsi(0) measured:   {cpsi_measured:.10f}")
    err = abs(cpsi_measured - cpsi_analytical)
    pr(f"  Error: {err:.2e}")
    p1 = err < 1e-10
    pr(f"  CPsi(0) matches 1/{d-1}: {'PASS' if p1 else 'FAIL'}")
    if not p1:
        all_pass = False

    # Propagate
    L = build_lindblad_superop(N, gamma=1.0)
    t_arr = np.linspace(0, 5.0, 2000)
    cpsi_traj = []
    for t in t_arr:
        rho_t = propagate_at(L, rho0, t)
        cpsi_traj.append(cpsi(rho_t))
    cpsi_traj = np.array(cpsi_traj)

    # Monotonic decay check
    diffs = np.diff(cpsi_traj)
    p2 = np.all(diffs <= 1e-10)
    pr(f"  Monotonically decreasing: {'PASS' if p2 else 'FAIL'}")
    if not p2:
        all_pass = False

    # Never crosses 1/4
    p3 = np.all(cpsi_traj < 0.25 + 1e-10)
    pr(f"  Always below 1/4: {'PASS' if p3 else 'FAIL'}")
    if not p3:
        all_pass = False

    # Final value near zero
    pr(f"  CPsi(K=5) = {cpsi_traj[-1]:.2e}")
    pr()

    results_a[N] = {"t": t_arr, "cpsi": cpsi_traj}

pr(f"Verification A overall: {'PASS' if all_pass else 'FAIL'}")
pr()
verif_a_pass = all_pass


# ============================================================
# VERIFICATION B: Bell+ smoke test
# ============================================================

pr("=" * 72)
pr("VERIFICATION B: Bell+ smoke test (regression)")
pr("=" * 72)
pr()

all_pass = True
rho0_bell = bell_plus()
L_bell = build_lindblad_superop(2, gamma=1.0)

# Fine time grid for crossing detection
t_fine = np.linspace(0, 2.0, 10000)
cpsi_bell_traj = []
for t in t_fine:
    rho_t = propagate_at(L_bell, rho0_bell, t)
    cpsi_bell_traj.append(cpsi(rho_t))
cpsi_bell_traj = np.array(cpsi_bell_traj)

# Find crossing at CPsi = 1/4 via bisection on the propagator
def cpsi_bell_at_t(t):
    rho_t = propagate_at(L_bell, rho0_bell, t)
    return cpsi(rho_t) - 0.25

# Bracket: CPsi(0) = 1/3 > 1/4, CPsi(2) << 1/4
t_cross_bell = brentq(cpsi_bell_at_t, 0.01, 2.0, xtol=1e-12)
pr(f"Bell+ crossing time t_cross = {t_cross_bell:.10f}")

rho_cross = propagate_at(L_bell, rho0_bell, t_cross_bell)
f_cross_bell = np.exp(-4 * t_cross_bell)
pr(f"  f_cross = exp(-4*t_cross) = {f_cross_bell:.7f}")

# Reference: 0.8612241
err_f = abs(f_cross_bell - 0.8612241)
p_f = err_f < 1e-4
pr(f"  Reference f_cross = 0.8612241, error = {err_f:.2e}: {'PASS' if p_f else 'FAIL'}")
if not p_f:
    all_pass = False

# Sector weights at crossing
sw_cross = sector_weights(rho_cross, 2)
W2_cross_bell = sw_cross.get(2, 0.0)
pr(f"  W_2 at crossing = {W2_cross_bell:.6f}")
err_w2 = abs(W2_cross_bell - 0.3709) / 0.3709
p_w2 = err_w2 < 0.005
pr(f"  Reference W_2 = 0.3709, rel error = {err_w2:.4e}: {'PASS' if p_w2 else 'FAIL'}")
if not p_w2:
    all_pass = False

# |dCPsi/dt| via central difference
dt_cd = 1e-7
rho_p = propagate_at(L_bell, rho0_bell, t_cross_bell + dt_cd)
rho_m = propagate_at(L_bell, rho0_bell, t_cross_bell - dt_cd)
dcpsi_num = (cpsi(rho_p) - cpsi(rho_m)) / (2 * dt_cd)
abs_dcpsi_bell = abs(dcpsi_num)
pr(f"  |dCPsi/dt|_cross = {abs_dcpsi_bell:.6f}")
err_dc = abs(abs_dcpsi_bell - 1.851701) / 1.851701
p_dc = err_dc < 0.005
pr(f"  Reference |dCPsi/dt|/gamma = 1.851701, rel error = {err_dc:.4e}: {'PASS' if p_dc else 'FAIL'}")
if not p_dc:
    all_pass = False

# Prefactor from direct dwell measurement
pr()
pr("  Direct dwell-time measurement:")
deltas = [0.01, 0.005, 0.002, 0.001]
prefactors_direct_bell = []
for delta in deltas:
    t_upper = brentq(lambda t: cpsi_bell_at_t(t) - delta, 0.001, t_cross_bell, xtol=1e-12)
    t_lower = brentq(lambda t: cpsi_bell_at_t(t) + delta, t_cross_bell, 2.0, xtol=1e-12)
    t_dwell = t_lower - t_upper
    pf = t_dwell / delta  # t_dwell spans full 2*delta band, prefactor = 2*gamma/|dCPsi/dt|
    prefactors_direct_bell.append(pf)
    pr(f"    delta={delta:.4f}: t_dwell={t_dwell:.8f}, prefactor={pf:.6f}")

pf_direct_bell = prefactors_direct_bell[-1]  # smallest delta
err_pf_direct = abs(pf_direct_bell - 1.080088) / 1.080088
p_pf_direct = err_pf_direct < 0.001
pr(f"  Prefactor (direct, delta=0.001) = {pf_direct_bell:.6f}")
pr(f"  Reference = 1.080088, rel error = {err_pf_direct:.4e}: {'PASS' if p_pf_direct else 'FAIL'}")
if not p_pf_direct:
    all_pass = False

# Prefactor from generalized formula
W0_bell = sw_cross.get(0, 0.0)
k_bell = 2
pf_formula_bell = (4.0 / k_bell) * (W0_bell + W2_cross_bell) / (W0_bell + 3 * W2_cross_bell)
err_pf_formula = abs(pf_formula_bell - 1.080088) / 1.080088
p_pf_formula = err_pf_formula < 0.001
pr(f"  Prefactor (generalized formula) = {pf_formula_bell:.6f}")
pr(f"  Reference = 1.080088, rel error = {err_pf_formula:.4e}: {'PASS' if p_pf_formula else 'FAIL'}")
if not p_pf_formula:
    all_pass = False

pr()
pr(f"Verification B overall: {'PASS' if all_pass else 'FAIL'}")
pr()
verif_b_pass = all_pass

results_b = {"t": t_fine, "cpsi": cpsi_bell_traj, "t_cross": t_cross_bell}


# ============================================================
# VERIFICATION C: W_3 out-of-sample
# ============================================================

pr("=" * 72)
pr("VERIFICATION C: W_3 out-of-sample test")
pr("=" * 72)
pr()

all_pass = True
rho0_w3 = w_state_3()

# --- Pauli decomposition ---
pr("Step C.1: Pauli decomposition of W_3")
coeffs = pauli_decomposition(rho0_w3, 3)
sw0 = sector_weights(rho0_w3, 3)

# Print nonzero coefficients
pr("  Nonzero Pauli coefficients a_P = Tr(rho*P):")
for label, val in sorted(coeffs.items()):
    if abs(val) > 1e-12:
        pr(f"    {label}: {val.real:+.6f}{val.imag:+.6f}j" if abs(val.imag) > 1e-12
           else f"    {label}: {val.real:+.10f}")

# Check W_0 and W_2
W0_w3 = sw0.get(0, 0.0)
W2_w3 = sw0.get(2, 0.0)
pr(f"\n  W_0 = {W0_w3:.12f} (expected 1/3 = {1/3:.12f})")
err_w0 = abs(W0_w3 - 1/3)
p_w0 = err_w0 < 1e-12
pr(f"  W_0 = 1/3: {'PASS' if p_w0 else 'FAIL'} (error {err_w0:.2e})")
if not p_w0:
    all_pass = False

pr(f"  W_2 = {W2_w3:.12f} (expected 2/3 = {2/3:.12f})")
err_w2_w3 = abs(W2_w3 - 2/3)
p_w2_w3 = err_w2_w3 < 1e-12
pr(f"  W_2 = 2/3: {'PASS' if p_w2_w3 else 'FAIL'} (error {err_w2_w3:.2e})")
if not p_w2_w3:
    all_pass = False

# Check no other sectors
for k in range(4):
    if k not in [0, 2]:
        wk = sw0.get(k, 0.0)
        if abs(wk) > 1e-12:
            pr(f"  WARNING: W_{k} = {wk:.2e} (expected 0)")

# CPsi(0)
cpsi0_w3 = cpsi(rho0_w3)
cpsi0_expected = 2.0 / 7.0
pr(f"\n  CPsi(0) = {cpsi0_w3:.12f} (expected 2/7 = {cpsi0_expected:.12f})")
err_cpsi0 = abs(cpsi0_w3 - cpsi0_expected)
p_cpsi0 = err_cpsi0 < 1e-10
pr(f"  CPsi(0) = 2/7: {'PASS' if p_cpsi0 else 'FAIL'} (error {err_cpsi0:.2e})")
if not p_cpsi0:
    all_pass = False

# --- Propagation ---
pr("\nStep C.2: Lindblad propagation of W_3")
L_w3 = build_lindblad_superop(3, gamma=1.0)

# Verify trace preservation at a test time
rho_test = propagate_at(L_w3, rho0_w3, 0.5)
tr_test = np.real(np.trace(rho_test))
pr(f"  Trace at t=0.5: {tr_test:.15f} (should be 1)")
assert abs(tr_test - 1.0) < 1e-12, f"Trace not preserved: {tr_test}"

# Fine trajectory
t_fine_w3 = np.linspace(0, 3.0, 10000)
cpsi_w3_traj = []
for t in t_fine_w3:
    rho_t = propagate_at(L_w3, rho0_w3, t)
    cpsi_w3_traj.append(cpsi(rho_t))
cpsi_w3_traj = np.array(cpsi_w3_traj)

# --- Crossing ---
pr("\nStep C.3: Find crossing at CPsi = 1/4")

def cpsi_w3_at_t(t):
    rho_t = propagate_at(L_w3, rho0_w3, t)
    return cpsi(rho_t) - 0.25

t_cross_w3 = brentq(cpsi_w3_at_t, 0.001, 2.0, xtol=1e-12)
f_cross_w3 = np.exp(-4 * t_cross_w3)
pr(f"  t_cross = {t_cross_w3:.10f}")
pr(f"  f_cross = exp(-4*t_cross) = {f_cross_w3:.10f}")

# Verify against cubic root: 16f^3 + 8f - 21 = 0
cubic_residual = 16 * f_cross_w3**3 + 8 * f_cross_w3 - 21
pr(f"  Cubic residual 16f^3 + 8f - 21 = {cubic_residual:.2e}")

# Sympy exact solution
pr("\n  Sympy exact solution of 16f^3 + 8f - 21 = 0:")
f_sym = sympy.Symbol("f", positive=True)
cubic_eq = 16 * f_sym**3 + 8 * f_sym - 21
roots = sympy.solve(cubic_eq, f_sym)
real_roots = [r for r in roots if r.is_real]
if real_roots:
    f_exact = real_roots[0]
    f_exact_float = float(f_exact.evalf(20))
    pr(f"  f_cross (sympy) = {f_exact}")
    pr(f"  f_cross (float) = {f_exact_float:.15f}")
    err_f_w3 = abs(f_cross_w3 - f_exact_float)
    p_f_w3 = err_f_w3 < 1e-6
    pr(f"  Error vs sympy: {err_f_w3:.2e}: {'PASS' if p_f_w3 else 'FAIL'}")
    if not p_f_w3:
        all_pass = False
else:
    pr("  WARNING: No real positive root found by sympy")
    p_f_w3 = False
    all_pass = False

# Sector weights at crossing
rho_cross_w3 = propagate_at(L_w3, rho0_w3, t_cross_w3)
sw_cross_w3 = sector_weights(rho_cross_w3, 3)
W2_cross_w3 = sw_cross_w3.get(2, 0.0)
W2_cross_expected = (2.0 / 3.0) * f_cross_w3**2
pr(f"\n  W_2 at crossing = {W2_cross_w3:.10f}")
pr(f"  Expected (2/3)*f^2 = {W2_cross_expected:.10f}")
err_w2c = abs(W2_cross_w3 - W2_cross_expected)
p_w2c = err_w2c < 1e-6
pr(f"  Error: {err_w2c:.2e}: {'PASS' if p_w2c else 'FAIL'}")
if not p_w2c:
    all_pass = False

# |dCPsi/dt| at crossing via central difference
dt_cd = 1e-7
rho_p = propagate_at(L_w3, rho0_w3, t_cross_w3 + dt_cd)
rho_m = propagate_at(L_w3, rho0_w3, t_cross_w3 - dt_cd)
dcpsi_num_w3 = (cpsi(rho_p) - cpsi(rho_m)) / (2 * dt_cd)
abs_dcpsi_w3 = abs(dcpsi_num_w3)

# Analytical: (8/21)*f*(1+6f^2) at f_cross
dcpsi_analytical_w3 = (8.0 / 21.0) * f_cross_w3 * (1 + 6 * f_cross_w3**2)
pr(f"\n  |dCPsi/dt| (numerical) = {abs_dcpsi_w3:.10f}")
pr(f"  |dCPsi/dt| (analytical) = {dcpsi_analytical_w3:.10f}")
err_dc_w3 = abs(abs_dcpsi_w3 - dcpsi_analytical_w3) / dcpsi_analytical_w3
p_dc_w3 = err_dc_w3 < 0.005
pr(f"  Rel error: {err_dc_w3:.4e}: {'PASS' if p_dc_w3 else 'FAIL'}")
if not p_dc_w3:
    all_pass = False

# Prefactor from generalized formula
W0_cross_w3 = sw_cross_w3.get(0, 0.0)
k_w3 = 2
pf_formula_w3 = (4.0 / k_w3) * (W0_cross_w3 + W2_cross_w3) / (W0_cross_w3 + 3 * W2_cross_w3)
pr(f"\n  Prefactor (generalized formula) = {pf_formula_w3:.6f}")

# Prefactor from direct dwell measurement
pr("  Direct dwell-time measurement:")
deltas_w3 = [0.01, 0.005, 0.002, 0.001]
prefactors_direct_w3 = []
for delta in deltas_w3:
    t_upper = brentq(lambda t: cpsi_w3_at_t(t) - delta, 0.001, t_cross_w3, xtol=1e-12)
    t_lower = brentq(lambda t: cpsi_w3_at_t(t) + delta, t_cross_w3, 2.5, xtol=1e-12)
    t_dwell = t_lower - t_upper
    pf = t_dwell / delta  # t_dwell spans full 2*delta band, prefactor = 2*gamma/|dCPsi/dt|
    prefactors_direct_w3.append(pf)
    pr(f"    delta={delta:.4f}: t_dwell={t_dwell:.8f}, prefactor={pf:.6f}")

pf_direct_w3 = prefactors_direct_w3[-1]
pr(f"  Prefactor (direct, delta=0.001) = {pf_direct_w3:.6f}")

# Agreement checks
err_pf_agree = abs(pf_formula_w3 - pf_direct_w3) / pf_direct_w3
p_pf_agree = err_pf_agree < 0.005
pr(f"\n  Formula vs direct rel error: {err_pf_agree:.4e}: {'PASS' if p_pf_agree else 'FAIL'}")
if not p_pf_agree:
    all_pass = False

# Both near 0.877
pf_expected_approx = 0.877
err_pf_val_formula = abs(pf_formula_w3 - pf_expected_approx) / pf_expected_approx
err_pf_val_direct = abs(pf_direct_w3 - pf_expected_approx) / pf_expected_approx
p_pf_val = err_pf_val_formula < 0.01 and err_pf_val_direct < 0.01
pr(f"  Formula value ~0.877: rel error {err_pf_val_formula:.4e}")
pr(f"  Direct value ~0.877: rel error {err_pf_val_direct:.4e}")
pr(f"  Both near 0.877: {'PASS' if p_pf_val else 'FAIL'}")
if not p_pf_val:
    all_pass = False

# Sympy exact prefactor
if real_roots:
    f_ex = real_roots[0]
    W2_ex = sympy.Rational(2, 3) * f_ex**2
    W0_ex = sympy.Rational(1, 3)
    pf_exact = 2 * (W0_ex + W2_ex) / (W0_ex + 3 * W2_ex)
    pf_exact_simplified = sympy.simplify(pf_exact)
    pf_exact_float = float(pf_exact.evalf(20))
    pr(f"\n  Exact prefactor (sympy): {pf_exact_simplified}")
    pr(f"  Exact prefactor (float): {pf_exact_float:.15f}")

pr(f"\nVerification C overall: {'PASS' if all_pass else 'FAIL'}")
pr()
verif_c_pass = all_pass


# ============================================================
# OVERALL RESULT
# ============================================================

overall = verif_a_pass and verif_b_pass and verif_c_pass
pr("=" * 72)
pr(f"OVERALL: {'ALL PASS' if overall else 'FAIL'}")
pr(f"  Verification A (GHZ below fold): {'PASS' if verif_a_pass else 'FAIL'}")
pr(f"  Verification B (Bell+ regression): {'PASS' if verif_b_pass else 'FAIL'}")
pr(f"  Verification C (W_3 out-of-sample): {'PASS' if verif_c_pass else 'FAIL'}")
pr("=" * 72)


# ============================================================
# Write text output
# ============================================================

out_txt = os.path.join(OUT_DIR, "dwell_prefactor_generalization_v1.txt")
with open(out_txt, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
print(f"\nText output written to {out_txt}")


# ============================================================
# 4-panel plot
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel 1: GHZ_3 and GHZ_4
ax = axes[0, 0]
for N in [3, 4]:
    ax.plot(results_a[N]["t"], results_a[N]["cpsi"], label=f"GHZ$_{N}$", linewidth=1.5)
ax.axhline(0.25, color="red", linestyle="--", alpha=0.7, label="C$\\Psi$ = 1/4")
ax.set_xlabel("K = $\\gamma t$")
ax.set_ylabel("C$\\Psi$")
ax.set_title("A: GHZ states born below the fold")
ax.legend()
ax.set_ylim(-0.01, 0.35)

# Panel 2: Bell+
ax = axes[0, 1]
ax.plot(results_b["t"], results_b["cpsi"], "b-", linewidth=1.5, label="Bell$^+$")
ax.axhline(0.25, color="red", linestyle="--", alpha=0.7, label="C$\\Psi$ = 1/4")
ax.axvline(results_b["t_cross"], color="green", linestyle=":", alpha=0.7, label=f"crossing at K={results_b['t_cross']:.4f}")
ax.set_xlabel("K = $\\gamma t$")
ax.set_ylabel("C$\\Psi$")
ax.set_title("B: Bell$^+$ smoke test")
ax.legend(fontsize=8)

# Panel 3: W_3
ax = axes[1, 0]
ax.plot(t_fine_w3, cpsi_w3_traj, "g-", linewidth=1.5, label="$W_3$")
ax.axhline(0.25, color="red", linestyle="--", alpha=0.7, label="C$\\Psi$ = 1/4")
ax.axvline(t_cross_w3, color="orange", linestyle=":", alpha=0.7, label=f"crossing at K={t_cross_w3:.4f}")
ax.set_xlabel("K = $\\gamma t$")
ax.set_ylabel("C$\\Psi$")
ax.set_title("C: $W_3$ out-of-sample test")
ax.legend(fontsize=8)

# Panel 4: Predicted vs measured prefactors
ax = axes[1, 1]
predicted = [pf_formula_bell, pf_formula_w3]
measured = [pf_direct_bell, pf_direct_w3]
labels_scatter = ["Bell$^+$", "$W_3$"]
colors_scatter = ["blue", "green"]
for i in range(2):
    ax.scatter(predicted[i], measured[i], c=colors_scatter[i], s=100, zorder=5, label=labels_scatter[i])
# y=x line
mn = min(min(predicted), min(measured)) - 0.05
mx = max(max(predicted), max(measured)) + 0.05
ax.plot([mn, mx], [mn, mx], "k--", alpha=0.5, label="y = x")
ax.set_xlabel("Predicted prefactor (formula)")
ax.set_ylabel("Measured prefactor (dwell time)")
ax.set_title("D: Formula vs measurement")
ax.legend()
ax.set_aspect("equal")

plt.tight_layout()
out_png = os.path.join(OUT_DIR, "dwell_prefactor_generalization_v1.png")
plt.savefig(out_png, dpi=150)
print(f"Plot written to {out_png}")
print("\nDONE")
