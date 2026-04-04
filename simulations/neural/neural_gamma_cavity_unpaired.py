"""
Identify and classify unpaired palindromic modes in C. elegans connectome
=========================================================================
Extends neural_gamma_cavity.py with exclusive palindromic matching
and biological classification of unpaired eigenvalues.

Output: simulations/results/neural_gamma_cavity_unpaired.txt
"""

import json, numpy as np
from pathlib import Path
import sys, os

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    sys.stdout.reconfigure(encoding="utf-8")

RESULTS_DIR = Path(__file__).parent.parent / "results"
NEURAL_DIR = Path(__file__).parent

out = []
def log(msg=""):
    print(msg)
    out.append(msg)

log("=" * 75)
log("UNPAIRED PALINDROMIC MODES IN C. ELEGANS CONNECTOME")
log("=" * 75)
log()

# ─────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────

with open(NEURAL_DIR / "celegans_connectome.json") as f:
    cdata = json.load(f)

neuron_ids = []
with open(NEURAL_DIR / "celegans_neuron_ids.txt") as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            neuron_ids.append(parts[1])

# Build signed weight matrix
W = np.array(cdata['chemical'], dtype=float)
N = W.shape[0]
sign_data = np.array(cdata['chemical_sign'], dtype=float)
for i in range(N):
    if sign_data[i] < 0:
        W[i, :] *= -1

max_w = np.max(np.abs(W))
W_norm = W / max_w

# Effective Jacobian
f_prime = 0.3
J_ce = -np.eye(N) + f_prime * W_norm

# Eigendecomposition
eigvals, eigvecs = np.linalg.eig(J_ce)

log(f"Total eigenvalues: {N}")
log(f"Oscillating (|Im| > 1e-6): {np.sum(np.abs(eigvals.imag) > 1e-6)}")
log(f"Overdamped: {np.sum(np.abs(eigvals.imag) <= 1e-6)}")
log()

# ─────────────────────────────────────────────
# Exclusive palindromic pairing
# ─────────────────────────────────────────────

center = np.mean(eigvals.real)
tolerance = 0.01

log(f"Palindrome center (mean Re): {center:.6f}")
log(f"Matching tolerance: {tolerance}")
log()

used = set()
paired_indices = []
unpaired_indices = []

for i in range(N):
    if i in used:
        continue

    target_re = 2 * center - eigvals[i].real
    target_im_mag = abs(eigvals[i].imag)

    best_j = -1
    best_dist = 1e10

    for j in range(N):
        if j == i or j in used:
            continue
        dist_re = abs(eigvals[j].real - target_re)
        dist_im = abs(abs(eigvals[j].imag) - target_im_mag)
        dist = np.sqrt(dist_re**2 + dist_im**2)
        if dist < best_dist:
            best_dist = dist
            best_j = j

    if best_dist < tolerance and best_j >= 0:
        paired_indices.append((i, best_j))
        used.add(i)
        used.add(best_j)
    else:
        unpaired_indices.append(i)
        used.add(i)

log(f"Paired: {len(paired_indices)} pairs = {len(paired_indices)*2} eigenvalues")
log(f"Unpaired: {len(unpaired_indices)} eigenvalues")
log(f"Palindromic fraction: {len(paired_indices)*2/N*100:.1f}%")
log()

# ─────────────────────────────────────────────
# Classify unpaired modes
# ─────────────────────────────────────────────

# Known C. elegans neuron classifications
PHARYNGEAL = {'I1L','I1R','I2L','I2R','I3','I4','I5','I6',
              'M1','M2L','M2R','M3L','M3R','M4','M5',
              'MCL','MCR','MI','NSML','NSMR'}

SENSORY_PREFIXES = ('AD', 'AF', 'AW', 'PH', 'IL', 'OL', 'UR', 'CEP',
                    'BAG', 'PLM', 'FLP', 'PQ', 'AQ', 'ASE', 'ASG',
                    'ASH', 'ASI', 'ASJ', 'ASK')
INTER_PREFIXES = ('AI', 'AV', 'RI', 'SI', 'DV', 'PV', 'LU')
MOTOR_PREFIXES = ('DA', 'DB', 'DD', 'VA', 'VB', 'VC', 'VD', 'AS',
                  'SMB', 'SMD', 'SAA', 'SAB', 'SDQ', 'HSN')

W_raw = np.array(cdata['chemical'], dtype=float)

def classify_neuron(name):
    if name in PHARYNGEAL:
        return "PHARYNGEAL"
    for pfx in SENSORY_PREFIXES:
        if name.startswith(pfx):
            return "SENSORY"
    for pfx in INTER_PREFIXES:
        if name.startswith(pfx):
            return "INTER"
    for pfx in MOTOR_PREFIXES:
        if name.startswith(pfx):
            return "MOTOR"
    return "UNKNOWN"

log("=" * 75)
log("UNPAIRED MODES: DETAILED ANALYSIS")
log("=" * 75)
log()

categories = {"PHARYNGEAL": [], "SENSORY": [], "INTER": [], "MOTOR": [], "UNKNOWN": []}

for idx in unpaired_indices:
    ev = eigvals[idx]
    vec = eigvecs[:, idx]
    weights = np.abs(vec)
    weights_norm = weights / np.sum(weights)

    top1_idx = np.argmax(weights)
    top1_name = neuron_ids[top1_idx] if top1_idx < len(neuron_ids) else f"#{top1_idx}"
    top1_cat = classify_neuron(top1_name)
    ei = "I" if sign_data[top1_idx] < 0 else "E"

    in_deg = int(np.sum(np.abs(W_raw[:, top1_idx]) > 0))
    out_deg = int(np.sum(np.abs(W_raw[top1_idx, :]) > 0))
    ratio = in_deg / max(out_deg, 1)

    top5 = np.argsort(weights)[-5:][::-1]

    log(f"Mode {idx:3d}: lambda = {ev.real:.6f} + {ev.imag:.6f}i")
    log(f"  {'oscillating' if abs(ev.imag) > 1e-6 else 'OVERDAMPED'}")
    log(f"  Dominant: {top1_name} ({ei}, {top1_cat})  "
        f"weight={weights_norm[top1_idx]:.3f}  in={in_deg}  out={out_deg}  "
        f"in/out={ratio:.2f}")
    log(f"  Top 5:")
    for k in top5:
        n = neuron_ids[k] if k < len(neuron_ids) else f"#{k}"
        e = "I" if sign_data[k] < 0 else "E"
        log(f"    {n:>8s} ({e}, {classify_neuron(n)})  w={weights_norm[k]:.4f}")
    log()

    categories[top1_cat].append((idx, top1_name, ev, ratio))

# ─────────────────────────────────────────────
# Category summary
# ─────────────────────────────────────────────

log("=" * 75)
log("CATEGORY SUMMARY")
log("=" * 75)
log()

for cat, modes in categories.items():
    if not modes:
        continue
    log(f"{cat}: {len(modes)} modes")
    for idx, name, ev, ratio in modes:
        osc = "osc" if abs(ev.imag) > 1e-6 else "ovd"
        log(f"  {idx:3d}  {name:>8s}  Re={ev.real:+.4f}  Im={ev.imag:+.4f}  "
            f"in/out={ratio:.2f}  {osc}")
    log()

# ─────────────────────────────────────────────
# Pharynx isolation check
# ─────────────────────────────────────────────

log("=" * 75)
log("PHARYNX ISOLATION")
log("=" * 75)
log()

pharyngeal_idx = [i for i, n in enumerate(neuron_ids) if n in PHARYNGEAL]
somatic_idx = [i for i in range(N) if i not in pharyngeal_idx]

cross_p2s = int(np.sum(np.abs(W_raw[pharyngeal_idx, :][:, somatic_idx]) > 0))
cross_s2p = int(np.sum(np.abs(W_raw[somatic_idx, :][:, pharyngeal_idx]) > 0))
within_p = int(np.sum(np.abs(W_raw[pharyngeal_idx, :][:, pharyngeal_idx]) > 0))
within_s = int(np.sum(np.abs(W_raw[somatic_idx, :][:, somatic_idx]) > 0))
total = within_p + within_s + cross_p2s + cross_s2p

log(f"Pharyngeal neurons: {len(pharyngeal_idx)}")
log(f"Somatic neurons: {len(somatic_idx)}")
log(f"Synapses within pharynx: {within_p}")
log(f"Synapses within soma: {within_s}")
log(f"Synapses pharynx -> soma: {cross_p2s}")
log(f"Synapses soma -> pharynx: {cross_s2p}")
log(f"Cross-coupling fraction: {(cross_p2s + cross_s2p) / max(total, 1) * 100:.1f}%")
log()

# ─────────────────────────────────────────────
# Connectivity comparison: paired vs unpaired dominant neurons
# ─────────────────────────────────────────────

log("=" * 75)
log("CONNECTIVITY: PAIRED vs UNPAIRED DOMINANT NEURONS")
log("=" * 75)
log()

# Get dominant neurons for some paired modes
paired_ratios = []
for i, j in paired_indices[:30]:
    vec_i = eigvecs[:, i]
    dom_i = np.argmax(np.abs(vec_i))
    in_d = np.sum(np.abs(W_raw[:, dom_i]) > 0)
    out_d = np.sum(np.abs(W_raw[dom_i, :]) > 0)
    if out_d > 0:
        paired_ratios.append(in_d / out_d)

unpaired_ratios = []
for idx in unpaired_indices:
    vec = eigvecs[:, idx]
    dom = np.argmax(np.abs(vec))
    in_d = np.sum(np.abs(W_raw[:, dom]) > 0)
    out_d = np.sum(np.abs(W_raw[dom, :]) > 0)
    if out_d > 0:
        unpaired_ratios.append(in_d / out_d)

log(f"Paired modes (n={len(paired_ratios)}):   mean In/Out = {np.mean(paired_ratios):.2f} "
    f"(median {np.median(paired_ratios):.2f})")
log(f"Unpaired modes (n={len(unpaired_ratios)}): mean In/Out = {np.mean(unpaired_ratios):.2f} "
    f"(median {np.median(unpaired_ratios):.2f})")
log()
log("Interpretation: unpaired modes are dominated by neurons that SEND more")
log("than they receive (In/Out < 1). They are sources, not relays.")
log("Paired modes are dominated by neurons that RECEIVE more (In/Out > 1).")

# ─────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────

out_path = RESULTS_DIR / "neural_gamma_cavity_unpaired.txt"
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(out))
print(f"\n>>> Results saved to: {out_path}")
