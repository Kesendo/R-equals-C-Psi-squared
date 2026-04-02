"""
Cockpit at N=5: The Missing Sweet Spot
=======================================
N=5 Heisenberg chain, the sweet spot throughout the repo.
PCA dimensionality, all 7 instruments, sacrifice zone comparison.

V(5) = 1.81 = 90% of max. 2+2=104 frequencies. 360x sacrifice improvement.
Liouvillian: 1024 x 1024. Pair (0,1)-(4,3) edge, (1,2)-(2,3) center.

April 2, 2026
"""
import numpy as np
from scipy import linalg, stats
from scipy.optimize import curve_fit
import sys, os, io, time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
os.makedirs(results_dir, exist_ok=True)
results_path = os.path.join(results_dir, 'cockpit_n5.txt')
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


def build_chain_H(n_qubits, J=1.0):
    dim = 2**n_qubits
    H = np.zeros((dim, dim), dtype=complex)
    for i in range(n_qubits - 1):
        H += gpt.two_qubit_heisenberg_term(i, i+1, n_qubits, J)
    return H


def build_initial_state(n_qubits):
    bell = np.array([1,0,0,1], dtype=complex)/np.sqrt(2)
    plus = np.array([1,1], dtype=complex)/np.sqrt(2)
    psi = bell
    for _ in range(2, n_qubits):
        psi = np.kron(psi, plus)
    return np.outer(psi, psi.conj())


def z_dephasing_ops(gammas, n_qubits):
    ops = []
    for q in range(n_qubits):
        if gammas[q] <= 0:
            continue
        paulis = [I2]*n_qubits
        paulis[q] = Z
        op = paulis[0]
        for p in paulis[1:]:
            op = np.kron(op, p)
        ops.append(np.sqrt(gammas[q]) * op)
    return ops


def build_liouvillian(H, gammas, n_qubits):
    dim = 2**n_qubits
    dim2 = dim**2
    Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for q in range(n_qubits):
        if gammas[q] <= 0:
            continue
        paulis = [I2]*n_qubits
        paulis[q] = Z
        Zq = paulis[0]
        for p in paulis[1:]:
            Zq = np.kron(Zq, p)
        L += gammas[q] * (np.kron(Zq, Zq.conj()) - np.eye(dim2))
    return L


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


def bures_distance(rho, sigma):
    try:
        sqrt_rho = linalg.sqrtm(rho)
        prod = sqrt_rho @ sigma @ sqrt_rho
        ev = np.real(np.linalg.eigvalsh(prod))
        fid = float(np.sum(np.sqrt(np.maximum(ev, 0))))**2
        return float(np.sqrt(max(0, 2*(1 - np.sqrt(max(0, min(1, fid)))))))
    except Exception:
        return 0.0


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


def simulate_and_analyze(name, n_qubits, gammas, H, focus_pair=(1,2)):
    """Full simulation + PCA + instruments for one noise profile."""
    L_ops = z_dephasing_ops(gammas, n_qubits)
    rho = build_initial_state(n_qubits)

    dt, sample_every = 0.005, 4
    n_steps = 2000
    nn_pairs = [(i, i+1) for i in range(n_qubits-1)]
    all_pairs = [(i, j) for i in range(n_qubits) for j in range(i+1, n_qubits)]

    times, feats = [], []
    pair_cpsi = {p: [] for p in all_pairs}
    pair_theta = {p: [] for p in all_pairs}
    pair_conc = {p: [] for p in all_pairs}
    focus_pur, focus_psi, focus_bures = [], [], []
    rho_prev_focus = None

    for step in range(n_steps + 1):
        t = step * dt
        if step % sample_every == 0:
            # Focus pair: full features + Bures
            rho_f = gpt.partial_trace_keep(rho, keep=list(focus_pair), n_qubits=n_qubits)
            feat = extract_features(rho_f)
            cp_f, pu_f, ps_f = cpsi_true(rho_f)

            if rho_prev_focus is not None:
                dB = bures_distance(rho_f, rho_prev_focus)
            else:
                dB = 0.0
            rho_prev_focus = rho_f.copy()

            times.append(t)
            feats.append(feat)
            focus_pur.append(pu_f)
            focus_psi.append(ps_f)
            focus_bures.append(dB)

            # ALL 10 pairs: CPsi, theta, Concurrence
            for p in all_pairs:
                rho_p = gpt.partial_trace_keep(rho, keep=list(p), n_qubits=n_qubits)
                cp, _, _ = cpsi_true(rho_p)
                pair_cpsi[p].append(cp)
                pair_theta[p].append(theta_deg(cp))
                pair_conc[p].append(gpt.concurrence_two_qubit(rho_p))

        if step < n_steps:
            rho = gpt.rk4_step(rho, H, L_ops, dt)

    # Convert
    times = np.array(times)
    X = np.array(feats)
    for p in all_pairs:
        pair_cpsi[p] = np.array(pair_cpsi[p])
        pair_theta[p] = np.array(pair_theta[p])
        pair_conc[p] = np.array(pair_conc[p])

    # PCA
    X_mean = X.mean(axis=0)
    X_std = X.std(axis=0) + 1e-10
    X_norm = (X - X_mean) / X_std
    U, S, Vt = np.linalg.svd(X_norm, full_matrices=False)
    var_exp = S**2 / np.sum(S**2)
    cum_var = np.cumsum(var_exp)
    scores = X_norm @ Vt.T

    n95 = int(np.searchsorted(cum_var, 0.95)) + 1

    # PC proxies
    cpsi_focus = pair_cpsi[focus_pair]
    theta_focus = pair_theta[focus_pair]
    conc_focus = pair_conc[focus_pair]
    pur_arr = np.array(focus_pur)
    psi_arr = np.array(focus_psi)

    observables = {'Concurrence': conc_focus, 'Purity': pur_arr,
                   'Psi-norm': psi_arr, 'CPsi': cpsi_focus, 'theta': theta_focus}
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

    return {
        'name': name, 'times': times, 'n_pts': len(times),
        'var_exp': var_exp, 'cum_var': cum_var, 'n95': n95,
        'Vt': Vt, 'scores': scores, 'proxies': proxies,
        'pair_cpsi': pair_cpsi, 'pair_theta': pair_theta,
        'pair_conc': pair_conc, 'nn_pairs': nn_pairs, 'all_pairs': all_pairs,
        'focus_pur': pur_arr, 'focus_psi': psi_arr,
        'focus_bures': np.array(focus_bures),
        'cpsi_focus': cpsi_focus, 'theta_focus': theta_focus,
        'conc_focus': conc_focus,
    }


# ================================================================
# PHASE 1: PCA AT N=5
# ================================================================
out("=" * 70)
out("COCKPIT AT N=5: THE MISSING SWEET SPOT")
out("Heisenberg chain [0-1-2-3-4], J=1.0, gamma=0.05")
out("Bell+(0,1) x |+>^3, focus pair: (1,2) center")
out("=" * 70)

N = 5
gamma = 0.05
H5 = build_chain_H(N, J=1.0)
gammas_uniform = [gamma]*N

t0 = time.time()
r = simulate_and_analyze("Chain N=5 uniform", N, gammas_uniform, H5, focus_pair=(1,2))
t_sim = time.time() - t0
out(f"\n  Simulation: {t_sim:.1f}s, {r['n_pts']} snapshots")

out(f"\n  PCA DIMENSIONALITY:")
out(f"  {'PC':>3} {'Var%':>7} {'Cum%':>7}")
out(f"  {'-'*20}")
for i in range(min(6, len(r['var_exp']))):
    out(f"  {i+1:>3} {r['var_exp'][i]*100:>7.1f} {r['cum_var'][i]*100:>7.1f}")
out(f"\n  PCs for 95% variance: {r['n95']}")

out(f"\n  PC LOADINGS:")
out(f"  {'':>8} | {'PC1':>7} {'PC2':>7} {'PC3':>7}")
out(f"  {'-'*35}")
for j in range(len(feat_names)):
    out(f"  {feat_names[j]:>8} | {r['Vt'][0,j]:>+7.3f} {r['Vt'][1,j]:>+7.3f} {r['Vt'][2,j]:>+7.3f}")

out(f"\n  PC PROXIES:")
for k, (pname, pr) in enumerate(r['proxies']):
    out(f"    PC{k+1} ~ {pname} (|r|={pr:.3f})")


# ================================================================
# PHASE 2: ALL 7 INSTRUMENTS
# ================================================================
out(f"\n{'=' * 70}")
out("PHASE 2: DASHBOARD (N=5, uniform noise)")
out("=" * 70)

# CPsi and theta per pair -- ALL 10 PAIRS
out(f"\n  CPsi RANGES PER PAIR (all {len(r['all_pairs'])} pairs):")
out(f"  {'Pair':>6} | {'min CPsi':>8} {'max CPsi':>8} {'max theta':>9} {'theta>0':>8} | {'Type':>8}")
out(f"  {'-'*60}")
for p in r['all_pairs']:
    cp = r['pair_cpsi'][p]
    th = r['pair_theta'][p]
    dist = abs(p[1] - p[0])
    if dist == 1:
        ptype = "NN-edge" if 0 in p or N-1 in p else "NN-cent"
    else:
        ptype = f"dist-{dist}"
    out(f"  {str(p):>6} | {cp.min():>8.4f} {cp.max():>8.4f} {th.max():>9.1f} "
        f"{int(np.sum(th>0)):>4}/{r['n_pts']} | {ptype:>8}")

# Dashboard at key times
key_t = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
out(f"\n  DASHBOARD (center pair (1,2)):")
out(f"  {'t':>5} | {'theta':>6} {'CPsi':>7} {'C':>6} {'Pur':>6} {'Psi':>6} | "
    f"{'PC1':>7} {'PC2':>7} {'PC3':>7} | {'v_B':>7}")
out(f"  {'-'*80}")
dt_sample = 0.005 * 4
for tk in key_t:
    idx = int(np.argmin(np.abs(r['times'] - tk)))
    vB = r['focus_bures'][idx] / dt_sample if idx > 0 else 0
    out(f"  {r['times'][idx]:>5.2f} | {r['theta_focus'][idx]:>6.1f} "
        f"{r['cpsi_focus'][idx]:>7.4f} {r['conc_focus'][idx]:>6.3f} "
        f"{r['focus_pur'][idx]:>6.3f} {r['focus_psi'][idx]:>6.3f} | "
        f"{r['scores'][idx,0]:>+7.2f} {r['scores'][idx,1]:>+7.2f} "
        f"{r['scores'][idx,2]:>+7.2f} | {vB:>7.4f}")

# Liouvillian eigenvalues
out(f"\n  LIOUVILLIAN SPECTRUM (1024x1024):")
t0 = time.time()
L_super = build_liouvillian(H5, gammas_uniform, N)
eigvals = linalg.eigvals(L_super)
t_eig = time.time() - t0
out(f"  Eigendecomposition: {t_eig:.1f}s")

unique_rates = sorted(set(np.round(-eigvals.real, 5)))
gap = unique_rates[1] if len(unique_rates) > 1 else 0
fastest = unique_rates[-1]
out(f"  Spectral gap (slowest decay): {gap:.5f}")
out(f"  Fastest decay rate: {fastest:.4f}")
out(f"  Number of distinct rates: {len(unique_rates)}")
out(f"  Predicted gap: 2*gamma = {2*gamma:.3f}")
out(f"  Predicted fastest: 2*(N-1)*gamma = {2*(N-1)*gamma:.3f}")

# Concurrence: edge vs center
out(f"\n  ENTANGLEMENT DYNAMICS: ALL 10 PAIRS")
out(f"  {'Pair':>6} {'Type':>8} | {'max C':>7} {'C(t=5)':>8} {'death':>10}")
out(f"  {'-'*50}")
for p in r['all_pairs']:
    conc = r['pair_conc'][p]
    t_death = None
    for i in range(1, len(conc)):
        if conc[i-1] > 0.01 and conc[i] < 0.01:
            t_death = r['times'][i]
            break
    dist = abs(p[1] - p[0])
    ptype = f"NN" if dist == 1 else f"d={dist}"
    if dist == 1 and (0 in p or N-1 in p):
        ptype = "NN-edge"
    elif dist == 1:
        ptype = "NN-cent"
    out(f"  {str(p):>6} {ptype:>8} | {conc.max():>7.3f} "
        f"{conc[min(250,len(conc)-1)]:>8.4f} "
        f"{'t='+f'{t_death:.2f}' if t_death else '  alive':>10}")

# Bures curvature at N=5
out(f"\n  BURES CURVATURE (center pair (1,2)):")
bv = r['focus_bures'].copy()
bv[0] = bv[1] if len(bv) > 1 else 0
dt_s = 0.005 * 4
bures_vel = bv / dt_s
dcpsi = np.gradient(r['cpsi_focus'], r['times'])
g_n5 = np.zeros(len(r['times']))
for i in range(len(r['times'])):
    if abs(dcpsi[i]) > 1e-12 and bures_vel[i] > 0:
        g_n5[i] = (bures_vel[i] / abs(dcpsi[i]))**2

valid_g = g_n5 > 1e-8
if np.sum(valid_g) > 6:
    cpsi_g = r['cpsi_focus'][valid_g]
    g_g = g_n5[valid_g]
    lng = np.log(g_g + 1e-30)
    d2lng = np.gradient(np.gradient(lng, cpsi_g), cpsi_g)
    K_n5 = -d2lng / (2 * g_g + 1e-30)

    out(f"  {'CPsi':>7} | {'g(CPsi)':>8} {'K_Gauss':>10}")
    out(f"  {'-'*30}")
    for cpsi_t in [0.28, 0.25, 0.22, 0.20, 0.18, 0.15]:
        idx = int(np.argmin(np.abs(cpsi_g - cpsi_t)))
        if abs(cpsi_g[idx] - cpsi_t) < 0.03:
            out(f"  {cpsi_g[idx]:>7.3f} | {g_g[idx]:>8.2f} {K_n5[idx]:>+10.1f}")

    idx_fold = int(np.argmin(np.abs(cpsi_g - 0.25)))
    out(f"\n  K at fold (CPsi~0.25): {K_n5[idx_fold]:.1f}")
    out(f"  N=2 simulation: K = -25 at fold")
else:
    out(f"  Insufficient valid g(CPsi) points for curvature")


# ================================================================
# PHASE 3: SACRIFICE ZONE
# ================================================================
out(f"\n{'=' * 70}")
out("PHASE 3: SACRIFICE ZONE (N=5)")
out("Same total gamma = 0.25, different distributions")
out("=" * 70)

total_gamma = N * gamma  # 0.25
profiles = {
    'Uniform':      [gamma]*N,
    'Edge sacrifice':  [0.20, 0.0125, 0.0125, 0.0125, 0.0125],
    'Double edge':  [0.10, 0.0167, 0.0167, 0.0167, 0.10],
}

sacrifice_results = {}
for prof_name, gammas in profiles.items():
    out(f"\n  Running {prof_name}: gamma={[f'{g:.4f}' for g in gammas]}...")
    t0 = time.time()
    sr = simulate_and_analyze(prof_name, N, gammas, H5, focus_pair=(1,2))
    out(f"    ({time.time()-t0:.1f}s)")
    sacrifice_results[prof_name] = sr

# Comparison table
out(f"\n  SACRIFICE ZONE COMPARISON (center pair (1,2)):")
out(f"\n  {'Instrument':>25} | {'Uniform':>10} {'Edge Sacr':>10} {'Double':>10} | {'Best':>12}")
out(f"  {'-'*75}")

instruments = [
    ('Max theta (center)', lambda r: r['theta_focus'].max()),
    ('Max CPsi (center)', lambda r: r['cpsi_focus'].max()),
    ('Max Concurrence (center)', lambda r: r['conc_focus'].max()),
    ('CPsi(t=5, center)', lambda r: r['cpsi_focus'][min(250, r['n_pts']-1)]),
    ('Conc(t=5, center)', lambda r: r['conc_focus'][min(250, r['n_pts']-1)]),
    ('CPsi(t=10, center)', lambda r: r['cpsi_focus'][-1]),
    ('PCs for 95%', lambda r: float(r['n95'])),
    ('PC1 variance %', lambda r: r['var_exp'][0]*100),
]

for inst_name, inst_fn in instruments:
    vals = {}
    for pn, sr in sacrifice_results.items():
        vals[pn] = inst_fn(sr)
    best = max(vals, key=vals.get)
    out(f"  {inst_name:>25} | {vals['Uniform']:>10.4f} {vals['Edge sacrifice']:>10.4f} "
        f"{vals['Double edge']:>10.4f} | {best:>12}")

# Most sensitive instrument
out(f"\n  MOST SENSITIVE INSTRUMENT (largest ratio Edge/Uniform):")
for inst_name, inst_fn in instruments:
    v_uni = inst_fn(sacrifice_results['Uniform'])
    v_edge = inst_fn(sacrifice_results['Edge sacrifice'])
    if v_uni > 1e-6:
        ratio = v_edge / v_uni
        out(f"    {inst_name:>25}: Edge/Uniform = {ratio:.2f}x")


# ================================================================
# PHASE 4: SUMMARY TABLE N=2-5
# ================================================================
out(f"\n{'=' * 70}")
out("PHASE 4: SCALING TABLE N=2 THROUGH N=5")
out("=" * 70)

# Run N=2, N=3, N=4 for comparison
scaling_results = []
for n in [2, 3, 4]:
    H_n = build_chain_H(n, J=1.0)
    gammas_n = [gamma]*n
    fp = (0, 1)
    sr = simulate_and_analyze(f"Chain N={n}", n, gammas_n, H_n, focus_pair=fp)
    scaling_results.append(sr)

# Add N=5
scaling_results.append(r)

out(f"\n  {'N':>3} | {'n95':>3} | {'PC1%':>5} {'PC2%':>5} {'PC3%':>5} | "
    f"{'PC1~':>12} {'|r|':>5} | {'max_th':>6} {'th>0':>6}")
out(f"  {'-'*70}")

for i, sr in enumerate(scaling_results):
    n = [2, 3, 4, 5][i]
    nq = sr['n_pts']
    n_theta = int(np.sum(sr['theta_focus'] > 0))
    out(f"  {n:>3} | {sr['n95']:>3} | {sr['var_exp'][0]*100:>5.1f} "
        f"{sr['var_exp'][1]*100:>5.1f} {sr['var_exp'][2]*100:>5.1f} | "
        f"{sr['proxies'][0][0]:>12} {sr['proxies'][0][1]:>5.2f} | "
        f"{sr['theta_focus'].max():>6.1f} {n_theta:>4}/{nq}")

# Dimensionality trend
dims = [sr['n95'] for sr in scaling_results]
out(f"\n  Dimensionality trend: {dims}")
if len(set(dims[1:])) == 1:
    out(f"  STABLE at {dims[1]} for N >= 3")
elif max(dims[1:]) - min(dims[1:]) <= 1:
    out(f"  STABLE within +/-1: {min(dims[1:])}-{max(dims[1:])} for N >= 3")
else:
    out(f"  GROWING: dimensionality increases with N")

# PC1 trend
pc1s = [sr['proxies'][0][0] for sr in scaling_results]
out(f"  PC1 identity: {pc1s}")


# ================================================================
# PHASE 5: FRAMEWORK IMPLICATIONS
# ================================================================
out(f"\n{'=' * 70}")
out("PHASE 5: WHAT N=5 REVEALS ABOUT THE FRAMEWORK")
out("=" * 70)

n5_n95 = r['n95']
n5_pc1 = r['proxies'][0][0]
n5_pc1_r = r['proxies'][0][1]

out(f"""
  1. DIMENSIONALITY AT N=5: {n5_n95} PCs for 95% variance.
     This {'confirms' if n5_n95 <= 4 else 'extends beyond'} the 3-4D pattern from N=3,4.
     The cockpit {'scales' if n5_n95 <= 4 else 'needs more instruments'} to N=5.

  2. PC1 AT N=5: {n5_pc1} (|r|={n5_pc1_r:.3f}).
     {'Consistent with N=4 chain (also Purity).' if n5_pc1 == 'Purity' else 'Different from N=4.'}

  3. THE 3-OBSERVABLE COCKPIT: For N=5 Heisenberg chain,
     monitoring {r['proxies'][0][0]}, {r['proxies'][1][0]}, and {r['proxies'][2][0]}
     captures {r['cum_var'][2]*100:.0f}% of the decoherence dynamics.
     This is 3 measurements instead of 4^5 = 1024 tomographic bases.

  4. SACRIFICE ZONE: Edge sacrifice {'enhances' if inst_fn(sacrifice_results['Edge sacrifice']) > inst_fn(sacrifice_results['Uniform']) else 'does not help'} center-pair coherence.
     The most sensitive instrument identifies where to look.

  5. SPECTRAL GAP: {gap:.4f} (predicted 2*gamma = {2*gamma:.3f}).
     {'Matches prediction.' if abs(gap - 2*gamma) < 0.001 else f'Deviates by {abs(gap-2*gamma):.4f}.'}
""")

out(f"\n  OPEN QUESTIONS (CORRECTED):")
out(f"  1. N=5 sweet spot: DONE (this analysis)")
out(f"  2. N=5 -> N=7: Does dimensionality hold?")
out(f"     (N=7 Liouvillian: 16384x16384, needs C# engine)")
out(f"  3. Non-Markovian noise: Does the cockpit hold?")
out(f"  4. Universal optimization: Which 3 observables +")
out(f"     which noise strategy optimizes any N?")

out(f"\n{'=' * 70}")
out("ANALYSIS COMPLETE")
out("=" * 70)

# Save
with open(results_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(_lines) + '\n')
print(f"\nResults saved to {results_path}")
