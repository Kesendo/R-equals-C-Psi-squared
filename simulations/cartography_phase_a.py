"""
STRUCTURAL CARTOGRAPHY - Phase A: Full Feature Charts
All 4 layers for all CΨ_AB visibility windows.
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

def bell_fids(rho):
    bells = {
        "Phi+": np.array([1,0,0,1], dtype=complex)/np.sqrt(2),
        "Phi-": np.array([1,0,0,-1], dtype=complex)/np.sqrt(2),
        "Psi+": np.array([0,1,1,0], dtype=complex)/np.sqrt(2),
        "Psi-": np.array([0,1,-1,0], dtype=complex)/np.sqrt(2),
    }
    return {n: float(np.real(np.trace(rho @ np.outer(b,b.conj())))) for n,b in bells.items()}

def von_neumann(rho):
    ev = np.real(np.linalg.eigvalsh(rho))
    ev = ev[ev > 1e-14]
    return float(-np.sum(ev * np.log2(ev)))

def bell_entropy(fids):
    p = np.array(list(fids.values()))
    p = p[p > 1e-14]
    return float(-np.sum(p * np.log2(p)))

def trace_dist(r1, r2):
    return 0.5 * np.sum(np.abs(np.linalg.eigvalsh(r1 - r2)))

def state_fidelity(r1, r2):
    from scipy.linalg import sqrtm
    sr = sqrtm(r1)
    prod = sr @ r2 @ sr
    ev = np.real(np.linalg.eigvalsh(prod))
    return float(np.sum(np.sqrt(np.maximum(ev, 0))))**2

def check_xx_sym(rho):
    XX = np.kron(np.array([[0,1],[1,0]]), np.array([[0,1],[1,0]]))
    return float(np.linalg.norm(rho @ XX - XX @ rho))

# ============================================================
# EVOLVE: noisy + unitary baseline
# ============================================================
print("=" * 70)
print("STRUCTURAL CARTOGRAPHY - Phase A")
print("=" * 70)

H = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
L_ops = gpt.dephasing_ops_n([0.05]*3)
L_zero = gpt.dephasing_ops_n([0.0]*3)
psi = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho_init = gpt.density_from_statevector(psi)

dt = 0.005
rho_n = rho_init.copy()
rho_u = rho_init.copy()

# Collect all windows (CΨ peaks)
peaks = []
prev_cpsi = 0
rising = True

for step in range(2001):
    t = step * dt
    if step % 4 == 0:
        rAB = gpt.partial_trace_keep(rho_n, keep=[1,2], n_qubits=3)
        rS  = gpt.partial_trace_keep(rho_n, keep=[0], n_qubits=3)
        rAB_u = gpt.partial_trace_keep(rho_u, keep=[1,2], n_qubits=3)
        c = gpt.concurrence_two_qubit(rAB)
        p = gpt.psi_norm(rAB)
        cpsi = c * p

        if prev_cpsi > 0.03 and cpsi < prev_cpsi and rising:
            peaks.append({
                't': t - dt*4, 'cpsi': prev_cpsi,
                'rAB': gpt.partial_trace_keep(rho_n, [1,2], 3).copy(),
                'rAB_u': gpt.partial_trace_keep(rho_u, [1,2], 3).copy(),
                'rS': gpt.partial_trace_keep(rho_n, [0], 3).copy(),
            })
            rising = False
        if cpsi > prev_cpsi:
            rising = True
        prev_cpsi = cpsi

    if step < 2000:
        rho_n = gpt.rk4_step(rho_n, H, L_ops, dt)
        rho_u = gpt.rk4_step(rho_u, H, L_zero, dt)

print(f"Found {len(peaks)} CPsi_AB peaks\n")

# ============================================================
# LAYER 1: Per-window feature chart
# ============================================================
print("LAYER 1: PER-WINDOW FEATURE CHART")
print("=" * 70)
hdr = (f"{'#':>2} {'t':>5} {'CPsi':>6} | {'F+':>5} {'F-':>5} {'P+':>5} {'P-':>5} "
       f"| {'H_B':>5} {'Pur':>5} {'SvN':>5} {'C':>5} {'Psi':>5} "
       f"| {'ph03':>6} {'ph12':>6} {'S_c':>5} {'XX':>8}")
print(hdr)
print("-" * len(hdr))

charts = []
for i, pk in enumerate(peaks[:12]):
    r = pk['rAB']
    f = bell_fids(r)
    pur = float(np.real(np.trace(r @ r)))
    svn = von_neumann(r)
    hb = bell_entropy(f)
    c_val = gpt.concurrence_two_qubit(r)
    psi_val = gpt.psi_norm(r)
    ph03 = np.angle(r[0,3]) / np.pi
    ph12 = np.angle(r[1,2]) / np.pi
    s_coh = gpt.l1_coherence(pk['rS'])
    xx = check_xx_sym(r)

    charts.append({
        'i': i, 't': pk['t'], 'cpsi': pk['cpsi'],
        'fids': f, 'pur': pur, 'svn': svn, 'hb': hb,
        'c': c_val, 'psi': psi_val,
        'ph03': ph03, 'ph12': ph12, 's_coh': s_coh, 'xx': xx,
        'rAB': r, 'rAB_u': pk['rAB_u'],
    })

    dom = max(f, key=f.get)
    print(f"{i:>2} {pk['t']:>5.2f} {pk['cpsi']:>6.3f} | "
          f"{f['Phi+']:>5.3f} {f['Phi-']:>5.3f} {f['Psi+']:>5.3f} {f['Psi-']:>5.3f} "
          f"| {hb:>5.2f} {pur:>5.3f} {svn:>5.2f} {c_val:>5.3f} {psi_val:>5.3f} "
          f"| {ph03:>+6.2f} {ph12:>+6.2f} {s_coh:>5.3f} {xx:>8.1e}")

# ============================================================
# LAYER 2: Inter-window transitions
# ============================================================
print(f"\nLAYER 2: INTER-WINDOW TRANSITIONS")
print("=" * 70)
print(f"{'n->n+1':>6} | {'TD':>6} {'Fid':>6} {'BellDr':>6} {'dph03':>6} {'dS_c':>6} | {'dt':>5}")
print("-" * 55)

for i in range(len(charts)-1):
    c1, c2 = charts[i], charts[i+1]
    td = trace_dist(c1['rAB'], c2['rAB'])
    fid = state_fidelity(c1['rAB'], c2['rAB'])
    # Bell vector drift
    bv1 = np.array([c1['fids'][k] for k in ['Phi+','Phi-','Psi+','Psi-']])
    bv2 = np.array([c2['fids'][k] for k in ['Phi+','Phi-','Psi+','Psi-']])
    bell_drift = np.linalg.norm(bv2 - bv1)
    # Phase advance
    dph = c2['ph03'] - c1['ph03']
    # S-coherence drift
    ds = c2['s_coh'] - c1['s_coh']
    dt_gap = c2['t'] - c1['t']

    print(f" {i:>2}->{i+1:<2} | {td:>6.3f} {fid:>6.3f} {bell_drift:>6.3f} "
          f"{dph:>+6.2f} {ds:>+6.3f} | {dt_gap:>5.2f}")

# ============================================================
# LAYER 4: Coherent vs noisy split
# ============================================================
print(f"\nLAYER 4: COHERENT vs NOISY SPLIT")
print("=" * 70)
print(f"{'#':>2} {'t':>5} | {'TD_nu':>6} {'dPur':>6} {'dph03':>7} | Verdict")
print("-" * 55)

for ch in charts:
    td_nu = trace_dist(ch['rAB'], ch['rAB_u'])
    pur_n = ch['pur']
    pur_u = float(np.real(np.trace(ch['rAB_u'] @ ch['rAB_u'])))
    dpur = pur_u - pur_n
    ph_n = np.angle(ch['rAB'][0,3])
    ph_u = np.angle(ch['rAB_u'][0,3])
    dph = abs(ph_n - ph_u)
    if dph > np.pi: dph = 2*np.pi - dph

    verdict = "phase OK" if dph < 0.01 else f"phase drift {dph:.3f}"
    print(f"{ch['i']:>2} {ch['t']:>5.2f} | {td_nu:>6.3f} {dpur:>6.3f} {dph:>7.4f} | {verdict}")

# ============================================================
# DIMENSIONALITY: How many independent features?
# ============================================================
print(f"\nDIMENSIONALITY ANALYSIS")
print("=" * 70)

# Build feature matrix: each row = one window, columns = features
features = []
feat_names = ['Phi+','Phi-','Psi+','Psi-','Pur','SvN','C','Psi','ph03','S_coh']
for ch in charts:
    row = [
        ch['fids']['Phi+'], ch['fids']['Phi-'],
        ch['fids']['Psi+'], ch['fids']['Psi-'],
        ch['pur'], ch['svn'], ch['c'], ch['psi'],
        ch['ph03'], ch['s_coh'],
    ]
    features.append(row)

X = np.array(features)
# Standardize
X_std = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-10)

# PCA via SVD
U, S, Vt = np.linalg.svd(X_std, full_matrices=False)
explained = S**2 / np.sum(S**2)
cumulative = np.cumsum(explained)

print(f"{'PC':>3} {'SingVal':>8} {'Var%':>7} {'Cum%':>7}")
print("-" * 30)
for i in range(min(len(S), 10)):
    print(f"{i+1:>3} {S[i]:>8.3f} {explained[i]*100:>7.1f} {cumulative[i]*100:>7.1f}")

# How many PCs for 95%?
n95 = int(np.searchsorted(cumulative, 0.95)) + 1
n99 = int(np.searchsorted(cumulative, 0.99)) + 1
print(f"\nDimensions for 95% variance: {n95}")
print(f"Dimensions for 99% variance: {n99}")
print(f"X*X symmetry says: 7 free parameters (out of 15)")
print(f"PCA says: {n95} effective dimensions (out of {len(feat_names)} features)")

# Top loadings for PC1 and PC2
print(f"\nPC1 loadings (strongest feature directions):")
for j in np.argsort(np.abs(Vt[0]))[::-1][:5]:
    print(f"  {feat_names[j]:>8}: {Vt[0,j]:>+.3f}")
print(f"PC2 loadings:")
for j in np.argsort(np.abs(Vt[1]))[::-1][:5]:
    print(f"  {feat_names[j]:>8}: {Vt[1,j]:>+.3f}")

print(f"\n{'=' * 70}")
print("PHASE A COMPLETE")
print("=" * 70)
