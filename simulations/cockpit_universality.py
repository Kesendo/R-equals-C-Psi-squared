"""
Cockpit Universality Test
=========================
Is the low-dimensional structure of decoherence universal?

Sweep:
  - Topologies: Chain, Star, Ring, Complete (N=2,3,4)
  - Noise: Z-dephasing, Depolarizing
  - Observed pair: always (0,1) with Bell+ initial entanglement

Question: Is the PCA dimensionality (~4 for 95%) stable across
topologies, or topology-specific?

April 2, 2026
"""
import numpy as np
from scipy import stats
import sys, os, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
os.makedirs(results_dir, exist_ok=True)
results_path = os.path.join(results_dir, 'cockpit_universality.txt')
_lines = []

def out(s=""):
    print(s)
    _lines.append(s)


# ================================================================
# HELPERS
# ================================================================
I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)


def build_hamiltonian(n_qubits, bonds, couplings):
    dim = 2**n_qubits
    H = np.zeros((dim, dim), dtype=complex)
    for (q1, q2), J in zip(bonds, couplings):
        H += gpt.two_qubit_heisenberg_term(q1, q2, n_qubits, J)
    return H


def build_initial_state(n_qubits):
    """Bell+ on (0,1), |+> on rest."""
    bell = np.array([1,0,0,1], dtype=complex)/np.sqrt(2)
    plus = np.array([1,1], dtype=complex)/np.sqrt(2)
    psi = bell
    for _ in range(2, n_qubits):
        psi = np.kron(psi, plus)
    return np.outer(psi, psi.conj())


def z_dephasing_ops(gamma, n_qubits):
    """sqrt(gamma) * Z_k for each qubit."""
    ops = []
    for q in range(n_qubits):
        paulis = [I2]*n_qubits
        paulis[q] = Z
        op = paulis[0]
        for p in paulis[1:]:
            op = np.kron(op, p)
        ops.append(np.sqrt(gamma) * op)
    return ops


def depolarizing_ops(gamma, n_qubits):
    """sqrt(gamma/3) * P_k for P in {X,Y,Z}, each qubit."""
    ops = []
    for q in range(n_qubits):
        for P in [X, Y, Z]:
            paulis = [I2]*n_qubits
            paulis[q] = P
            op = paulis[0]
            for p in paulis[1:]:
                op = np.kron(op, p)
            ops.append(np.sqrt(gamma/3) * op)
    return ops


def bell_fids(rho):
    bells = {
        "Phi+": np.array([1,0,0,1], dtype=complex)/np.sqrt(2),
        "Phi-": np.array([1,0,0,-1], dtype=complex)/np.sqrt(2),
        "Psi+": np.array([0,1,1,0], dtype=complex)/np.sqrt(2),
        "Psi-": np.array([0,1,-1,0], dtype=complex)/np.sqrt(2),
    }
    return {n: float(np.real(np.trace(rho @ np.outer(b, b.conj()))))
            for n, b in bells.items()}


def von_neumann(rho):
    ev = np.real(np.linalg.eigvalsh(rho))
    ev = ev[ev > 1e-14]
    return float(-np.sum(ev * np.log2(ev)))


def cpsi_true(rho):
    pur = float(np.real(np.trace(rho @ rho)))
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    psi = l1 / (d - 1) if d > 1 else 0.0
    return pur * psi, pur, psi


def theta_deg(cpsi_val):
    return float(np.degrees(np.arctan(np.sqrt(4*cpsi_val - 1)))) if cpsi_val > 0.25 else 0.0


feat_names = ['Phi+', 'Phi-', 'Psi+', 'Psi-', 'Pur', 'SvN', 'C', 'Psi', 'ph03']


def extract_features(rho_pair):
    f = bell_fids(rho_pair)
    pur = float(np.real(np.trace(rho_pair @ rho_pair)))
    svn = von_neumann(rho_pair)
    c = gpt.concurrence_two_qubit(rho_pair)
    psi = gpt.psi_norm(rho_pair)
    ph03 = np.angle(rho_pair[0, 3]) / np.pi
    return [f['Phi+'], f['Phi-'], f['Psi+'], f['Psi-'],
            pur, svn, c, psi, ph03]


def run_analysis(name, n_qubits, bonds, couplings, noise_type, gamma=0.05):
    """Run full trajectory + PCA for one configuration."""
    H = build_hamiltonian(n_qubits, bonds, couplings)
    if noise_type == 'z_dephasing':
        L_ops = z_dephasing_ops(gamma, n_qubits)
    else:
        L_ops = depolarizing_ops(gamma, n_qubits)

    rho = build_initial_state(n_qubits)
    dt, sample_every = 0.005, 4
    n_steps = 2000

    features, cpsis, thetas, concs, purs, psis = [], [], [], [], [], []

    for step in range(n_steps + 1):
        if step % sample_every == 0:
            if n_qubits == 2:
                rho_pair = rho
            else:
                rho_pair = gpt.partial_trace_keep(rho, keep=[0, 1], n_qubits=n_qubits)

            feat = extract_features(rho_pair)
            cp, pu, ps = cpsi_true(rho_pair)

            features.append(feat)
            cpsis.append(cp)
            thetas.append(theta_deg(cp))
            concs.append(gpt.concurrence_two_qubit(rho_pair))
            purs.append(pu)
            psis.append(ps)

        if step < n_steps:
            rho = gpt.rk4_step(rho, H, L_ops, dt)

    X = np.array(features)
    cpsis = np.array(cpsis)
    thetas = np.array(thetas)
    concs = np.array(concs)
    purs = np.array(purs)
    psis = np.array(psis)
    n_pts = len(cpsis)

    # PCA
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0) + 1e-10
    X_norm = (X - X_mean) / X_std
    U, S, Vt = np.linalg.svd(X_norm, full_matrices=False)
    var_exp = S**2 / np.sum(S**2)
    cum_var = np.cumsum(var_exp)

    n90 = int(np.searchsorted(cum_var, 0.90)) + 1
    n95 = int(np.searchsorted(cum_var, 0.95)) + 1
    n99 = int(np.searchsorted(cum_var, 0.99)) + 1

    # PC scores
    scores = X_norm @ Vt.T

    # Best proxy for PC1, PC2, PC3
    observables = {'Concurrence': concs, 'Purity': purs, 'Psi-norm': psis,
                   'CPsi': cpsis, 'theta': thetas}
    proxies = []
    for k in range(min(3, scores.shape[1])):
        best_r, best_name = 0, "?"
        for oname, ovals in observables.items():
            if np.std(ovals) > 1e-10 and np.std(scores[:, k]) > 1e-10:
                r = abs(stats.pearsonr(ovals, scores[:, k])[0])
                if r > best_r:
                    best_r = r
                    best_name = oname
        proxies.append((best_name, best_r))

    # Top loading for PC1
    pc1_top = feat_names[np.argmax(np.abs(Vt[0]))]

    return {
        'name': name, 'n_qubits': n_qubits, 'noise': noise_type,
        'n_pts': n_pts,
        'var_exp': var_exp, 'cum_var': cum_var,
        'n90': n90, 'n95': n95, 'n99': n99,
        'proxies': proxies, 'pc1_top': pc1_top,
        'Vt': Vt,
        'theta_range': (thetas.min(), thetas.max()),
        'cpsi_range': (cpsis.min(), cpsis.max()),
        'n_quantum': int(np.sum(thetas > 0)),
    }


# ================================================================
# TOPOLOGY SWEEP
# ================================================================
out("=" * 70)
out("COCKPIT UNIVERSALITY TEST")
out("Is the ~4D decoherence structure universal?")
out("=" * 70)

configs = [
    # N=2 baseline
    ('Pair N=2', 2, [(0,1)], [1.0]),

    # N=3 topologies
    ('Chain N=3', 3, [(0,1),(1,2)], [1.0,1.0]),
    ('Star N=3 sym', 3, [(0,1),(0,2)], [1.0,1.0]),
    ('Star N=3 asym', 3, [(0,1),(0,2)], [1.0,2.0]),
    ('Ring N=3', 3, [(0,1),(1,2),(2,0)], [1.0,1.0,1.0]),
    ('Complete N=3', 3, [(0,1),(0,2),(1,2)], [1.0,1.0,1.0]),

    # N=4 topologies
    ('Chain N=4', 4, [(0,1),(1,2),(2,3)], [1.0]*3),
    ('Star N=4', 4, [(0,1),(0,2),(0,3)], [1.0]*3),
    ('Ring N=4', 4, [(0,1),(1,2),(2,3),(3,0)], [1.0]*4),
]

results_z = []
results_d = []

for topo_name, n, bonds, J in configs:
    out(f"\n  Running {topo_name}...")

    # Z-dephasing
    r = run_analysis(topo_name, n, bonds, J, 'z_dephasing')
    results_z.append(r)

    # Depolarizing
    r_d = run_analysis(topo_name + ' (depol)', n, bonds, J, 'depolarizing')
    results_d.append(r_d)


# ================================================================
# RESULTS: Z-DEPHASING
# ================================================================
out(f"\n{'=' * 70}")
out("RESULTS: Z-DEPHASING (gamma=0.05)")
out("=" * 70)

out(f"\n  {'Topology':>20} | {'n90':>3} {'n95':>3} {'n99':>3} | {'PC1%':>5} {'PC2%':>5} {'PC3%':>5} | "
    f"{'PC1~':>12} {'|r|':>5} | {'theta>0':>8}")
out(f"  {'-'*90}")

for r in results_z:
    out(f"  {r['name']:>20} | {r['n90']:>3} {r['n95']:>3} {r['n99']:>3} | "
        f"{r['var_exp'][0]*100:>5.1f} {r['var_exp'][1]*100:>5.1f} {r['var_exp'][2]*100:>5.1f} | "
        f"{r['proxies'][0][0]:>12} {r['proxies'][0][1]:>5.2f} | "
        f"{r['n_quantum']:>4}/{r['n_pts']}")

# PC1 loading comparison
out(f"\n  PC1 TOP LOADING (which feature dominates PC1?):")
out(f"  {'Topology':>20} | {'Top feat':>8} | {'Load':>6} | {'2nd feat':>8} {'Load':>6}")
out(f"  {'-'*65}")
for r in results_z:
    Vt = r['Vt']
    idx1 = np.argsort(np.abs(Vt[0]))[::-1]
    out(f"  {r['name']:>20} | {feat_names[idx1[0]]:>8} | {Vt[0,idx1[0]]:>+6.3f} | "
        f"{feat_names[idx1[1]]:>8} {Vt[0,idx1[1]]:>+6.3f}")


# ================================================================
# RESULTS: DEPOLARIZING
# ================================================================
out(f"\n{'=' * 70}")
out("RESULTS: DEPOLARIZING NOISE (gamma=0.05)")
out("=" * 70)

out(f"\n  {'Topology':>20} | {'n90':>3} {'n95':>3} {'n99':>3} | {'PC1%':>5} {'PC2%':>5} {'PC3%':>5} | "
    f"{'PC1~':>12} {'|r|':>5}")
out(f"  {'-'*80}")

for r in results_d:
    out(f"  {r['name']:>20} | {r['n90']:>3} {r['n95']:>3} {r['n99']:>3} | "
        f"{r['var_exp'][0]*100:>5.1f} {r['var_exp'][1]*100:>5.1f} {r['var_exp'][2]*100:>5.1f} | "
        f"{r['proxies'][0][0]:>12} {r['proxies'][0][1]:>5.2f}")


# ================================================================
# COMPARISON: Z-DEPHASING vs DEPOLARIZING
# ================================================================
out(f"\n{'=' * 70}")
out("COMPARISON: Z-DEPHASING vs DEPOLARIZING")
out("=" * 70)

out(f"\n  {'Topology':>20} | {'n95(Z)':>6} {'n95(D)':>6} | {'PC1(Z)':>12} {'PC1(D)':>12}")
out(f"  {'-'*65}")
for rz, rd in zip(results_z, results_d):
    out(f"  {rz['name']:>20} | {rz['n95']:>6} {rd['n95']:>6} | "
        f"{rz['proxies'][0][0]:>12} {rd['proxies'][0][0]:>12}")


# ================================================================
# UNIVERSALITY ANALYSIS
# ================================================================
out(f"\n{'=' * 70}")
out("UNIVERSALITY ANALYSIS")
out("=" * 70)

# Dimensionality statistics
dims_z = [r['n95'] for r in results_z]
dims_d = [r['n95'] for r in results_d]

out(f"\n  Dimensionality for 95% variance:")
out(f"    Z-dephasing: {dims_z}")
out(f"    Mean = {np.mean(dims_z):.1f}, Std = {np.std(dims_z):.1f}, Range = [{min(dims_z)}, {max(dims_z)}]")
out(f"    Depolarizing: {dims_d}")
out(f"    Mean = {np.mean(dims_d):.1f}, Std = {np.std(dims_d):.1f}, Range = [{min(dims_d)}, {max(dims_d)}]")

# Is PC1 always Concurrence?
pc1_z = [r['proxies'][0][0] for r in results_z]
pc1_d = [r['proxies'][0][0] for r in results_d]
conc_count_z = sum(1 for p in pc1_z if p == 'Concurrence')
conc_count_d = sum(1 for p in pc1_d if p == 'Concurrence')

out(f"\n  PC1 ~ Concurrence?")
out(f"    Z-dephasing: {conc_count_z}/{len(pc1_z)} topologies")
out(f"    Depolarizing: {conc_count_d}/{len(pc1_d)} topologies")
out(f"    PC1 proxies (Z): {pc1_z}")
out(f"    PC1 proxies (D): {pc1_d}")

# PC1 variance fraction
pc1_vars_z = [r['var_exp'][0]*100 for r in results_z]
pc1_vars_d = [r['var_exp'][0]*100 for r in results_d]
out(f"\n  PC1 variance fraction:")
out(f"    Z-dephasing: {[f'{v:.0f}%' for v in pc1_vars_z]}")
out(f"    Mean = {np.mean(pc1_vars_z):.1f}%, Range = [{min(pc1_vars_z):.0f}%, {max(pc1_vars_z):.0f}%]")
out(f"    Depolarizing: {[f'{v:.0f}%' for v in pc1_vars_d]}")
out(f"    Mean = {np.mean(pc1_vars_d):.1f}%, Range = [{min(pc1_vars_d):.0f}%, {max(pc1_vars_d):.0f}%]")

# Verdict
out(f"\n  VERDICT:")
dim_stable = (max(dims_z) - min(dims_z)) <= 2 and (max(dims_d) - min(dims_d)) <= 2
pc1_stable = conc_count_z >= len(pc1_z) * 0.7

if dim_stable and pc1_stable:
    out(f"    UNIVERSAL: Dimensionality and PC1 are stable across topologies.")
    out(f"    The cockpit instruments have the same meaning regardless of topology.")
elif dim_stable:
    out(f"    PARTIALLY UNIVERSAL: Dimensionality is stable, but PC1 identity varies.")
    out(f"    The number of instruments is universal, but their meaning is topology-specific.")
else:
    out(f"    TOPOLOGY-SPECIFIC: Both dimensionality and PC1 vary with topology.")
    out(f"    Each architecture needs its own cockpit calibration.")

# Additional detail
if not pc1_stable:
    out(f"\n  When PC1 != Concurrence, what is it?")
    for r in results_z:
        if r['proxies'][0][0] != 'Concurrence':
            out(f"    {r['name']}: PC1 ~ {r['proxies'][0][0]} (|r|={r['proxies'][0][1]:.2f})")

out(f"\n{'=' * 70}")
out("ANALYSIS COMPLETE")
out("=" * 70)

# Save
with open(results_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(_lines) + '\n')
print(f"\nResults saved to {results_path}")
