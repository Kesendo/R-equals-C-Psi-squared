"""
BRIGHT-TRANSITION MAP: Exact diagonalization + visibility analysis
For 4-qubit star S-A-B-C, sweep J_SC, compute which Bohr frequencies
are visible to AB and with what weight.
"""
import numpy as np
import sys
sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")
import star_topology_v3 as gpt

I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)

def op_at(op, qubit, n_q):
    ops = [I2]*n_q
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result

def star_H(J_list, n_q):
    H = np.zeros((2**n_q, 2**n_q), dtype=complex)
    for i, J in enumerate(J_list):
        for p in [X, Y, Z]:
            H += J * op_at(p, 0, n_q) @ op_at(p, i+1, n_q)
    return H

def bright_lines(H, rho0, O_AB, threshold=0.001):
    """
    Exact diagonalization + visibility analysis.
    Returns list of (frequency, weight) for all bright Bohr transitions.
    """
    # Diagonalize
    evals, V = np.linalg.eigh(H)
    
    # Transform to energy basis
    rho_tilde = V.conj().T @ rho0 @ V
    O_tilde = V.conj().T @ O_AB @ V
    
    # Compute visibility weight for each transition
    d = len(evals)
    transitions = {}
    
    for m in range(d):
        for n in range(d):
            freq = (evals[m] - evals[n]) / (2 * np.pi)
            if freq < 0.01:  # skip zero and negative
                continue
            
            weight = abs(rho_tilde[m, n] * O_tilde[n, m])
            
            if weight < threshold * 0.01:  # skip truly negligible
                continue
            
            # Bin to nearest 0.001
            freq_bin = round(freq, 3)
            if freq_bin in transitions:
                transitions[freq_bin] += weight
            else:
                transitions[freq_bin] = weight
    
    # Sort by weight
    result = [(f, w) for f, w in transitions.items()]
    result.sort(key=lambda x: -x[1])
    return result

# Setup: 4-qubit star, S(0)-A(1)-B(2)-C(3)
n_q = 4
J_SA, J_SB = 1.0, 2.0

# Observable: c+ on AB = I_S x (YZ + ZY)/sqrt(2) x I_C
O_AB = (op_at(Y, 1, n_q) @ op_at(Z, 2, n_q) + 
        op_at(Z, 1, n_q) @ op_at(Y, 2, n_q)) / np.sqrt(2)

# Initial state: Bell_SA x |+>_B x |+>_C
bell = gpt.bell_phi_plus()
plus = gpt.plus_state()
psi = np.kron(np.kron(bell, plus), plus)
rho0 = np.outer(psi, psi.conj())

print("=" * 70)
print("BRIGHT-TRANSITION MAP")
print("4-qubit star: S-A(1.0)-B(2.0)-C(J_SC)")
print("Observable: c+ on AB pair")
print("=" * 70)

# Sweep J_SC
jsc_values = [0.0, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]

# Collect all bright frequencies across sweep for the map
all_freqs = set()
sweep_data = {}

for jsc in jsc_values:
    H = star_H([J_SA, J_SB, jsc], n_q)
    lines = bright_lines(H, rho0, O_AB)
    sweep_data[jsc] = lines
    for f, w in lines:
        if w > 0.001:
            all_freqs.add(round(f, 2))

# Print the map
print(f"\n--- Bright lines per J_SC (weight > 0.001) ---\n")
print(f"{'J_SC':>5} | {'#lines':>6} | Top 6 bright lines (freq: weight)")
print("-" * 75)

for jsc in jsc_values:
    lines = [(f,w) for f,w in sweep_data[jsc] if w > 0.001]
    top = lines[:6]
    desc = "  ".join(f"{f:.3f}:{w:.4f}" for f,w in top)
    print(f"{jsc:>5.1f} | {len(lines):>6} | {desc}")

# Print the full transition map: rows=frequencies, columns=J_SC
sorted_freqs = sorted(all_freqs)
print(f"\n--- FULL BRIGHT-TRANSITION MAP ---")
print(f"Rows = Bohr frequencies, Columns = J_SC")
print(f"Values = visibility weight (blank = invisible)\n")

# Header
header = f"{'freq':>7} |"
for jsc in jsc_values:
    header += f" {jsc:>5.1f}"
print(header)
print("-" * (8 + 6 * len(jsc_values)))

for freq in sorted_freqs:
    row = f"{freq:>7.3f} |"
    for jsc in jsc_values:
        # Find weight at this freq for this J_SC
        w = 0
        for f, wt in sweep_data[jsc]:
            if abs(f - freq) < 0.02:
                w = max(w, wt)
        if w > 0.001:
            row += f" {w:>5.3f}"
        else:
            row += f"     ."
    print(row)

# Key analysis: what changes with J_SC?
print(f"\n--- ANALYSIS ---\n")

# 1. Are new frequencies genuinely new eigenvalue gaps?
baseline = set(round(f,2) for f,w in sweep_data[0.0] if w > 0.001)
print(f"Baseline (J_SC=0): {len(baseline)} bright lines: {sorted(baseline)}")

for jsc in [0.5, 1.0, 2.0, 3.0]:
    current = set(round(f,2) for f,w in sweep_data[jsc] if w > 0.001)
    new = current - baseline
    lost = baseline - current
    print(f"J_SC={jsc}: {len(current)} lines. NEW: {sorted(new)}  LOST: {sorted(lost)}")

# 2. Which lines are most sensitive to J_SC?
print(f"\n--- Sensitivity: dW/dJ_SC for each frequency ---")
print(f"{'freq':>7} | {'W at 0':>6} {'W at 1':>6} {'W at 3':>6} | {'Change':>7} | Type")
print("-" * 60)

for freq in sorted_freqs:
    w0 = max([w for f,w in sweep_data[0.0] if abs(f-freq)<0.02], default=0)
    w1 = max([w for f,w in sweep_data[1.0] if abs(f-freq)<0.02], default=0)
    w3 = max([w for f,w in sweep_data[3.0] if abs(f-freq)<0.02], default=0)
    
    change = w3 - w0
    if abs(change) < 0.001 and w0 < 0.001:
        continue  # skip invisible throughout
    
    if w0 < 0.001 and w3 > 0.001:
        typ = "APPEARS"
    elif w0 > 0.001 and w3 < 0.001:
        typ = "DISAPPEARS"
    elif abs(change) > 0.01:
        typ = "SHIFTS" if change > 0 else "FADES"
    else:
        typ = "stable"
    
    print(f"{freq:>7.3f} | {w0:>6.4f} {w1:>6.4f} {w3:>6.4f} | {change:>+7.4f} | {typ}")

# 3. Can AB invert? How many independent sensitive lines?
print(f"\n--- Invertibility: Can AB recover J_SC? ---")
sensitive_lines = []
for freq in sorted_freqs:
    weights = []
    for jsc in jsc_values:
        w = max([w for f,w in sweep_data[jsc] if abs(f-freq)<0.02], default=0)
        weights.append(w)
    # Check if this line varies monotonically with J_SC
    diffs = [weights[i+1] - weights[i] for i in range(len(weights)-1)]
    varies = any(abs(d) > 0.005 for d in diffs)
    if varies:
        sensitive_lines.append(freq)

print(f"Lines sensitive to J_SC: {len(sensitive_lines)} of {len(sorted_freqs)}")
print(f"Sensitive frequencies: {[f'{f:.3f}' for f in sensitive_lines]}")
if len(sensitive_lines) >= 1:
    print(f"-> Detection: YES (at least one line changes)")
if len(sensitive_lines) >= 2:
    print(f"-> Characterization: MAYBE (multiple independent lines)")
else:
    print(f"-> Characterization: UNLIKELY (not enough independent info)")

print(f"\n{'='*70}")
print("BRIGHT-TRANSITION MAP COMPLETE")
print("="*70)
