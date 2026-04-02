"""
Theta-PC Correlation Analysis
==============================
Correlates theta = arctan(sqrt(4*CPsi - 1)) with PCA components
from structural cartography. Tests whether theta is a function of
one PC or a combination.

Setup: Star topology (S=0, A=1, B=2)
  Initial: Bell_SA x |+>_B
  Coupling: J_SA=1.0, J_SB=2.0
  Noise: Z-dephasing gamma=0.05

CPsi definition: Purity x Psi-norm = Tr(rho^2) x L1/(d-1)
theta definition: arctan(sqrt(4*CPsi - 1)) for CPsi > 1/4, else 0

April 2, 2026
"""
import numpy as np
from scipy import stats
import sys, os, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

# ================================================================
# Output: print to stdout AND collect for results file
# ================================================================
results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
os.makedirs(results_dir, exist_ok=True)
results_path = os.path.join(results_dir, 'theta_pc_correlation.txt')
_lines = []

def out(s=""):
    print(s)
    _lines.append(s)


# ================================================================
# Helpers
# ================================================================
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
    """CPsi = Purity x Psi-norm (the R=CPsi^2 definition)."""
    pur = float(np.real(np.trace(rho @ rho)))
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    psi = l1 / (d - 1) if d > 1 else 0.0
    return pur * psi, pur, psi


def theta_deg(cpsi_val):
    """theta in degrees. 0 if CPsi <= 1/4."""
    return float(np.degrees(np.arctan(np.sqrt(4*cpsi_val - 1)))) if cpsi_val > 0.25 else 0.0


def extract_features(rAB, rS=None):
    """10 features matching cartography_phase_a.py."""
    f = bell_fids(rAB)
    pur = float(np.real(np.trace(rAB @ rAB)))
    svn = von_neumann(rAB)
    c = gpt.concurrence_two_qubit(rAB)
    psi = gpt.psi_norm(rAB)
    ph03 = np.angle(rAB[0, 3]) / np.pi
    s_coh = gpt.l1_coherence(rS) if rS is not None else 0.0
    return [f['Phi+'], f['Phi-'], f['Psi+'], f['Psi-'],
            pur, svn, c, psi, ph03, s_coh]


feat_names = ['Phi+', 'Phi-', 'Psi+', 'Psi-', 'Pur', 'SvN', 'C', 'Psi', 'ph03', 'S_coh']


# ================================================================
# PHASE 1: SIMULATION
# ================================================================
out("=" * 70)
out("THETA-PC CORRELATION ANALYSIS")
out("Star: Bell_SA x |+>_B, J_SA=1.0, J_SB=2.0, gamma=0.05")
out("CPsi = Purity x Psi-norm (R=CPsi^2 definition)")
out("=" * 70)

H = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
L_ops = gpt.dephasing_ops_n([0.05] * 3)
psi_init = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho = gpt.density_from_statevector(psi_init)

dt, sample_every = 0.005, 4
n_steps = 2000

# Storage
times, feats_AB = [], []
cpsi_AB, theta_AB_arr, pur_AB, psinorm_AB, conc_AB = [], [], [], [], []
cpsi_SA, theta_SA_arr = [], []
cpsi_SB, theta_SB_arr = [], []
# Also track cartography-style CPsi (Concurrence x Psi-norm) for comparison
cpsi_cart_AB = []

out("\nPhase 1: Simulating trajectory...")

for step in range(n_steps + 1):
    t = step * dt
    if step % sample_every == 0:
        rAB = gpt.partial_trace_keep(rho, keep=[1, 2], n_qubits=3)
        rSA = gpt.partial_trace_keep(rho, keep=[0, 1], n_qubits=3)
        rSB = gpt.partial_trace_keep(rho, keep=[0, 2], n_qubits=3)
        rS  = gpt.partial_trace_keep(rho, keep=[0], n_qubits=3)

        # AB pair: true CPsi + features
        cp_ab, pu_ab, ps_ab = cpsi_true(rAB)
        cn_ab = gpt.concurrence_two_qubit(rAB)
        feat = extract_features(rAB, rS)

        times.append(t)
        feats_AB.append(feat)
        cpsi_AB.append(cp_ab)
        theta_AB_arr.append(theta_deg(cp_ab))
        pur_AB.append(pu_ab)
        psinorm_AB.append(ps_ab)
        conc_AB.append(cn_ab)
        cpsi_cart_AB.append(cn_ab * ps_ab)

        # SA pair
        cp_sa, _, _ = cpsi_true(rSA)
        cpsi_SA.append(cp_sa)
        theta_SA_arr.append(theta_deg(cp_sa))

        # SB pair
        cp_sb, _, _ = cpsi_true(rSB)
        cpsi_SB.append(cp_sb)
        theta_SB_arr.append(theta_deg(cp_sb))

    if step < n_steps:
        rho = gpt.rk4_step(rho, H, L_ops, dt)

# Convert to arrays
times = np.array(times)
X = np.array(feats_AB)
cpsi_AB = np.array(cpsi_AB);       theta_AB = np.array(theta_AB_arr)
cpsi_SA = np.array(cpsi_SA);       theta_SA = np.array(theta_SA_arr)
cpsi_SB = np.array(cpsi_SB);       theta_SB = np.array(theta_SB_arr)
pur_AB  = np.array(pur_AB);        psinorm_AB = np.array(psinorm_AB)
conc_AB = np.array(conc_AB);       cpsi_cart_AB = np.array(cpsi_cart_AB)

N = len(times)

out(f"  {N} snapshots, t in [0, {times[-1]:.1f}]")
out(f"\n  CPsi ranges (Purity x Psi-norm, the R=CPsi^2 definition):")
out(f"    AB: [{cpsi_AB.min():.4f}, {cpsi_AB.max():.4f}]  "
    f"theta>0: {np.sum(theta_AB > 0)} pts ({100*np.sum(theta_AB > 0)/N:.1f}%)")
out(f"    SA: [{cpsi_SA.min():.4f}, {cpsi_SA.max():.4f}]  "
    f"theta>0: {np.sum(theta_SA > 0)} pts ({100*np.sum(theta_SA > 0)/N:.1f}%)")
out(f"    SB: [{cpsi_SB.min():.4f}, {cpsi_SB.max():.4f}]  "
    f"theta>0: {np.sum(theta_SB > 0)} pts ({100*np.sum(theta_SB > 0)/N:.1f}%)")

out(f"\n  Cartography CPsi (Concurrence x Psi-norm) for comparison:")
out(f"    AB: [{cpsi_cart_AB.min():.4f}, {cpsi_cart_AB.max():.4f}]  "
    f"(used for peak detection in cartography)")

out(f"\n  NOTE: Cartography used Concurrence x Psi-norm, this analysis")
out(f"  uses Purity x Psi-norm. They differ for mixed states.")


# ================================================================
# PCA ON AB FEATURES
# ================================================================
out(f"\n{'=' * 70}")
out("PCA ON AB-PAIR FEATURES (full trajectory)")
out("=" * 70)

X_mean = X.mean(axis=0)
X_std = X.std(axis=0) + 1e-10
X_norm = (X - X_mean) / X_std

U, S_vals, Vt = np.linalg.svd(X_norm, full_matrices=False)
var_exp = S_vals**2 / np.sum(S_vals**2)
cum_var = np.cumsum(var_exp)

out(f"\n{'PC':>3} {'SingVal':>8} {'Var%':>7} {'Cum%':>7}")
out("-" * 30)
for i in range(min(len(S_vals), 6)):
    out(f"{i+1:>3} {S_vals[i]:>8.3f} {var_exp[i]*100:>7.1f} {cum_var[i]*100:>7.1f}")

n95 = int(np.searchsorted(cum_var, 0.95)) + 1
n99 = int(np.searchsorted(cum_var, 0.99)) + 1
out(f"\n  Dimensions for 95% variance: {n95}")
out(f"  Dimensions for 99% variance: {n99}")

scores = X_norm @ Vt.T
PC1, PC2, PC3 = scores[:, 0], scores[:, 1], scores[:, 2]

for k in range(3):
    out(f"\n  PC{k+1} loadings (top 5):")
    for j in np.argsort(np.abs(Vt[k]))[::-1][:5]:
        out(f"    {feat_names[j]:>8}: {Vt[k, j]:>+.3f}")


# ================================================================
# PHASE 2: CORRELATION ANALYSIS
# ================================================================
out(f"\n{'=' * 70}")
out("PHASE 2: CORRELATION ANALYSIS")
out("=" * 70)


def run_correlations(pair_name, theta, cpsi_arr, pc1, pc2, pc3,
                     pur=None, psi=None, conc=None):
    """Run Pearson/Spearman correlations and linear regression."""
    masks = [
        ("All points", np.ones(len(theta), dtype=bool)),
        ("theta > 0 only", theta > 0),
    ]

    for sub_label, mask in masks:
        n = int(np.sum(mask))
        if n < 5:
            out(f"\n  [{sub_label}]: {n} points, too few, skipping")
            continue

        th = theta[mask]
        if np.std(th) < 1e-10:
            out(f"\n  [{sub_label}] ({n} pts): theta constant ({th[0]:.1f} deg), skipping")
            continue

        cp = cpsi_arr[mask]
        p1, p2, p3 = pc1[mask], pc2[mask], pc3[mask]

        out(f"\n  [{sub_label}] ({n} points)")

        # --- Pearson ---
        out(f"\n    Pearson r:")
        out(f"    {'':>12} | {'PC1':>8} {'PC2':>8} {'PC3':>8}")
        out(f"    {'-' * 47}")

        rows = [("theta", th), ("CPsi", cp)]
        if pur is not None:
            rows += [("Concurrence", conc[mask]), ("Purity", pur[mask]),
                     ("Psi-norm", psi[mask])]
        for name, vals in rows:
            if np.std(vals) < 1e-10:
                out(f"    r({name:>10}, PCk) |      n/a      n/a      n/a")
                continue
            r1 = stats.pearsonr(vals, p1)[0]
            r2 = stats.pearsonr(vals, p2)[0]
            r3 = stats.pearsonr(vals, p3)[0]
            out(f"    r({name:>10}, PCk) | {r1:>+8.4f} {r2:>+8.4f} {r3:>+8.4f}")

        # --- Spearman ---
        out(f"\n    Spearman rho:")
        out(f"    {'':>12} | {'PC1':>8} {'PC2':>8} {'PC3':>8}")
        out(f"    {'-' * 47}")
        for name, vals in [("theta", th), ("CPsi", cp)]:
            if np.std(vals) < 1e-10:
                continue
            r1 = stats.spearmanr(vals, p1).statistic
            r2 = stats.spearmanr(vals, p2).statistic
            r3 = stats.spearmanr(vals, p3).statistic
            out(f"    rho({name:>8}, PCk) | {r1:>+8.4f} {r2:>+8.4f} {r3:>+8.4f}")

        # --- Regression ---
        ss_tot = np.sum((th - np.mean(th))**2)
        if ss_tot < 1e-14:
            continue

        out(f"\n    Linear regression: theta = a1*PC1 + a2*PC2 + a3*PC3 + c")
        A = np.column_stack([p1, p2, p3, np.ones(n)])
        coeffs, _, _, _ = np.linalg.lstsq(A, th, rcond=None)
        r2 = 1 - np.sum((th - A @ coeffs)**2) / ss_tot
        out(f"      a1={coeffs[0]:>+8.4f}  a2={coeffs[1]:>+8.4f}  "
            f"a3={coeffs[2]:>+8.4f}  c={coeffs[3]:>+8.4f}")
        out(f"      R^2 = {r2:.6f}")

        for k, pc in enumerate([p1, p2, p3]):
            Ak = np.column_stack([pc, np.ones(n)])
            ck, _, _, _ = np.linalg.lstsq(Ak, th, rcond=None)
            r2k = 1 - np.sum((th - Ak @ ck)**2) / ss_tot
            out(f"      R^2(theta ~ PC{k+1} only) = {r2k:.6f}")


# --- AB pair (primary, matches cartography) ---
out(f"\n--- AB pair (matching cartography topology) ---")
run_correlations("AB", theta_AB, cpsi_AB, PC1, PC2, PC3,
                 pur=pur_AB, psi=psinorm_AB, conc=conc_AB)

# --- If AB lacks theta variation, also try SA and SB ---
for pair_name, theta_pair, cpsi_pair in [("SA", theta_SA, cpsi_SA),
                                          ("SB", theta_SB, cpsi_SB)]:
    n_q = int(np.sum(theta_pair > 0))
    if n_q >= 10:
        out(f"\n--- {pair_name} pair (theta>0 in {n_q} snapshots) ---")
        out(f"  CPsi range: [{cpsi_pair.min():.4f}, {cpsi_pair.max():.4f}]")
        out(f"  theta range: [{theta_pair.min():.1f} deg, {theta_pair.max():.1f} deg]")
        run_correlations(pair_name, theta_pair, cpsi_pair, PC1, PC2, PC3)


# ================================================================
# PHASE 3: MISSING DIMENSIONS
# ================================================================
out(f"\n{'=' * 70}")
out("PHASE 3: MISSING DIMENSIONS")
out("=" * 70)

# Find pair with most theta > 0 variation
best_pair = "AB"
best_theta = theta_AB
best_cpsi = cpsi_AB
for name, th, cp in [("SA", theta_SA, cpsi_SA), ("SB", theta_SB, cpsi_SB)]:
    if np.sum(th > 0) > np.sum(best_theta > 0):
        best_pair = name
        best_theta = th
        best_cpsi = cp

mask_q = best_theta > 0
n_q = int(np.sum(mask_q))

if n_q >= 5 and np.std(best_theta[mask_q]) > 1e-10:
    th_q = best_theta[mask_q]

    r_abs = []
    for k in range(3):
        r = abs(stats.pearsonr(th_q, scores[:, k][mask_q])[0])
        r_abs.append(r)

    out(f"\n  Best pair for theta analysis: {best_pair} ({n_q} quantum-regime points)")
    for k in range(3):
        out(f"  |r(theta_{best_pair}, PC{k+1})| = {r_abs[k]:.4f}  "
            f"(PC{k+1}: {var_exp[k]*100:.1f}% var)")

    dominant = int(np.argmax(r_abs))
    missing = [k for k in range(3) if k != dominant]

    out(f"\n  theta aligns most with PC{dominant+1} (|r|={r_abs[dominant]:.4f})")
    out(f"  Missing dimensions for complete 3D coordinate:")
    for m in missing:
        out(f"    PC{m+1} ({var_exp[m]*100:.1f}% variance)")

    out(f"\n  theta is a 1D scalar. The manifold is {n95}D (for 95% var).")
    out(f"  theta covers at most 1 direction. The other {n95-1} are invisible to theta.")
else:
    out(f"\n  No pair has sufficient theta > 0 variation for this analysis.")
    out(f"  The manifold lies entirely in the classical regime (CPsi < 1/4).")
    out(f"  The compass has no direction here: theta = 0 throughout.")

    # Even without theta, report which features the PCs capture
    out(f"\n  The 3D manifold structure exists independently of theta:")
    for k in range(3):
        top_feat = np.argsort(np.abs(Vt[k]))[::-1][0]
        out(f"    PC{k+1} ({var_exp[k]*100:.1f}% var): dominated by {feat_names[top_feat]}")


# ================================================================
# PHASE 4: CONNECTION TO theta-FIDELITY r=0.87
# ================================================================
out(f"\n{'=' * 70}")
out("PHASE 4: CONNECTION TO theta-FIDELITY r=0.87")
out("=" * 70)

out(f"\n  Reference: THETA_PALINDROME_ECHO measured r(theta_SB, F_channel) = 0.87")
out(f"  for the SB pair in channel mode. r^2 = 0.757, unexplained = 24.3%.")

if n_q >= 5 and np.std(best_theta[mask_q]) > 1e-10:
    th_q = best_theta[mask_q]
    ss_tot_q = np.sum((th_q - np.mean(th_q))**2)

    # 3-PC regression in quantum regime
    p1q = PC1[mask_q]; p2q = PC2[mask_q]; p3q = PC3[mask_q]
    A3 = np.column_stack([p1q, p2q, p3q, np.ones(n_q)])
    c3, _, _, _ = np.linalg.lstsq(A3, th_q, rcond=None)
    r2_3pc = 1 - np.sum((th_q - A3 @ c3)**2) / ss_tot_q

    out(f"\n  theta-PC(3D) regression (quantum regime): R^2 = {r2_3pc:.4f}")

    for k in range(3):
        Ak = np.column_stack([scores[:, k][mask_q], np.ones(n_q)])
        ck, _, _, _ = np.linalg.lstsq(Ak, th_q, rcond=None)
        r2k = 1 - np.sum((th_q - Ak @ ck)**2) / ss_tot_q
        out(f"  R^2(theta ~ PC{k+1}): {r2k:.4f}")

    dominant = int(np.argmax(r_abs))
    if r_abs[dominant] > 0.9:
        out(f"\n  theta is essentially a PC{dominant+1} proxy (|r|={r_abs[dominant]:.3f}).")
        out(f"  The 24% fidelity gap likely comes from the PCs that theta ignores.")
    elif r_abs[dominant] > 0.7:
        out(f"\n  theta primarily tracks PC{dominant+1} (|r|={r_abs[dominant]:.3f}).")
        out(f"  Some fidelity gap may stem from the non-PC{dominant+1} components.")
    else:
        out(f"\n  theta is not well-described by any single PC.")
        out(f"  The fidelity gap has a more complex origin.")
else:
    out(f"\n  The theta-Fidelity correlation was measured on CPsi_SB (channel mode).")
    out(f"  This PCA uses AB-pair features. The spaces are different.")
    out(f"  A direct quantitative connection requires PCA on SB features,")
    out(f"  which is beyond the scope of this analysis.")


# ================================================================
# BONUS: NONLINEAR ANALYSIS (Phase 2, Step 3)
# ================================================================
out(f"\n{'=' * 70}")
out("NONLINEAR ANALYSIS")
out("=" * 70)

if n_q >= 10 and np.std(best_theta[mask_q]) > 1e-10:
    th_q = best_theta[mask_q]
    ss_tot_q = np.sum((th_q - np.mean(th_q))**2)

    # theta vs sqrt(PC2^2 + PC3^2)
    p2q = PC2[mask_q]; p3q = PC3[mask_q]
    r_circ = np.sqrt(p2q**2 + p3q**2)
    if np.std(r_circ) > 1e-10:
        Ar = np.column_stack([r_circ, np.ones(n_q)])
        cr, _, _, _ = np.linalg.lstsq(Ar, th_q, rcond=None)
        r2_circ = 1 - np.sum((th_q - Ar @ cr)**2) / ss_tot_q
        r_pearson = stats.pearsonr(th_q, r_circ)[0]
        out(f"\n  theta vs sqrt(PC2^2 + PC3^2): r={r_pearson:+.4f}, R^2={r2_circ:.4f}")

    # theta vs PC_k^2 (quadratic)
    for k in range(3):
        pk = scores[:, k][mask_q]
        Aq = np.column_stack([pk, pk**2, np.ones(n_q)])
        cq, _, _, _ = np.linalg.lstsq(Aq, th_q, rcond=None)
        r2q = 1 - np.sum((th_q - Aq @ cq)**2) / ss_tot_q
        out(f"  theta ~ PC{k+1} + PC{k+1}^2 (quadratic): R^2 = {r2q:.4f}")
else:
    out(f"\n  Insufficient theta>0 data for nonlinear analysis.")


# ================================================================
# CAVEATS
# ================================================================
out(f"\n{'=' * 70}")
out("CAVEATS")
out("=" * 70)
out(f"")
out(f"  1. CPsi definition: This analysis uses CPsi = Purity x Psi-norm")
out(f"     (the R=CPsi^2 definition). The cartography used")
out(f"     Concurrence x Psi-norm for peak detection. These differ")
out(f"     for mixed states.")
out(f"")
out(f"  2. PCA basis: Computed on full trajectory ({N} snapshots),")
out(f"     not just the 9 peak windows from the cartography.")
out(f"     Loadings may differ from the cartography PCA.")
out(f"")
out(f"  3. Topology: Star (N=3) with asymmetric coupling (J_SB=2*J_SA).")
out(f"     Results may not transfer to other topologies.")
out(f"")
out(f"  4. Correlation != causation. Co-variation along a decoherence")
out(f"     trajectory does not imply a functional relationship.")

out(f"\n{'=' * 70}")
out("ANALYSIS COMPLETE")
out("=" * 70)

# Save results
with open(results_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(_lines) + '\n')
print(f"\nResults saved to {results_path}")
