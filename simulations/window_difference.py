"""What's in the difference between the two sharpest windows?"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

H = gpt.star_hamiltonian_n(n_observers=2, J_SA=1.0, J_SB=2.0)
L_ops = gpt.dephasing_ops_n([0.05]*3)
psi = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho = gpt.density_from_statevector(psi)

dt = 0.005
peaks = []
prev_cpsi = 0
rising = True

for step in range(2001):
    t = step * dt
    if step % 4 == 0:
        rAB = gpt.partial_trace_keep(rho, [1,2], 3)
        c = gpt.concurrence_two_qubit(rAB)
        p = gpt.psi_norm(rAB)
        cpsi = c * p
        if prev_cpsi > 0.03 and cpsi < prev_cpsi and rising:
            peaks.append({'t': t-dt*4, 'cpsi': prev_cpsi, 'rAB': rAB.copy()})
            rising = False
        if cpsi > prev_cpsi: rising = True
        prev_cpsi = cpsi
    if step < 2000:
        rho = gpt.rk4_step(rho, H, L_ops, dt)

r0 = peaks[0]['rAB']  # sharpest window 0
r1 = peaks[1]['rAB']  # sharpest window 1

diff = r1 - r0

print("=" * 60)
print("DIFFERENCE BETWEEN THE TWO SHARPEST WINDOWS")
print(f"Window 0: t={peaks[0]['t']:.2f}, CPsi={peaks[0]['cpsi']:.3f}")
print(f"Window 1: t={peaks[1]['t']:.2f}, CPsi={peaks[1]['cpsi']:.3f}")
print("=" * 60)

print(f"\nWindow 0 (rho_AB):")
for i in range(4):
    row = "  ".join(f"{r0[i,j].real:+.4f}{r0[i,j].imag:+.4f}j" for j in range(4))
    print(f"  [{row}]")

print(f"\nWindow 1 (rho_AB):")
for i in range(4):
    row = "  ".join(f"{r1[i,j].real:+.4f}{r1[i,j].imag:+.4f}j" for j in range(4))
    print(f"  [{row}]")

print(f"\nDIFFERENCE (Window 1 - Window 0):")
for i in range(4):
    row = "  ".join(f"{diff[i,j].real:+.5f}{diff[i,j].imag:+.5f}j" for j in range(4))
    print(f"  [{row}]")

# Eigendecomposition of the difference
evals, evecs = np.linalg.eigh(diff)
print(f"\nEigenvalues of the difference:")
for i, ev in enumerate(evals):
    print(f"  lambda_{i} = {ev:+.6f}")

# Pauli decomposition of the difference
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
si = np.eye(2, dtype=complex)

paulis = {'II': np.kron(si,si), 'IX': np.kron(si,sx), 'IY': np.kron(si,sy), 'IZ': np.kron(si,sz),
          'XI': np.kron(sx,si), 'XX': np.kron(sx,sx), 'XY': np.kron(sx,sy), 'XZ': np.kron(sx,sz),
          'YI': np.kron(sy,si), 'YX': np.kron(sy,sx), 'YY': np.kron(sy,sy), 'YZ': np.kron(sy,sz),
          'ZI': np.kron(sz,si), 'ZX': np.kron(sz,sx), 'ZY': np.kron(sz,sy), 'ZZ': np.kron(sz,sz)}

print(f"\nPauli decomposition of difference (non-zero components):")
coeffs = {}
for name, P in paulis.items():
    c = np.real(np.trace(diff @ P)) / 4
    coeffs[name] = c

for name in sorted(coeffs, key=lambda x: abs(coeffs[x]), reverse=True):
    c = coeffs[name]
    if abs(c) > 1e-6:
        bar = "#" * int(abs(c) * 200)
        print(f"  {name}: {c:+.6f}  {bar}")

# What fraction of the difference is in each Pauli direction?
total = sum(c**2 for c in coeffs.values())
print(f"\nDominant directions (fraction of ||diff||^2):")
for name in sorted(coeffs, key=lambda x: abs(coeffs[x]), reverse=True)[:5]:
    c = coeffs[name]
    frac = c**2 / total if total > 0 else 0
    print(f"  {name}: {frac*100:.1f}%")

# Now compare: difference between OTHER adjacent pairs
print(f"\n{'=' * 60}")
print("COMPARISON: ALL ADJACENT DIFFERENCES")
print("=" * 60)

print(f"\n{'Pair':>6} | {'||diff||':>8} | {'Top Pauli':>10} | {'coeff':>8} | {'frac%':>6} | Mode")
print("-" * 65)

for k in range(min(len(peaks)-1, 8)):
    d = peaks[k+1]['rAB'] - peaks[k]['rAB']
    norm = np.linalg.norm(d)
    
    cs = {}
    for name, P in paulis.items():
        cs[name] = np.real(np.trace(d @ P)) / 4
    
    tot = sum(c**2 for c in cs.values())
    top_name = max(cs, key=lambda x: abs(cs[x]))
    top_c = cs[top_name]
    top_frac = top_c**2 / tot if tot > 0 else 0
    
    # Detect mode
    td = 0.5 * np.sum(np.abs(np.linalg.eigvalsh(d)))
    mode = "SWITCH" if td > 0.5 else "glide"
    
    print(f" {k}->{k+1} | {norm:>8.4f} | {top_name:>10} | {top_c:>+8.5f} | {top_frac*100:>5.1f} | {mode}")

# The key question: is the glide-mode difference always the same direction?
print(f"\n{'=' * 60}")
print("GLIDE-MODE DIFFERENCES: Same direction?")
print("=" * 60)

glide_diffs = []
for k in range(min(len(peaks)-1, 8)):
    d = peaks[k+1]['rAB'] - peaks[k]['rAB']
    td = 0.5 * np.sum(np.abs(np.linalg.eigvalsh(d)))
    if td < 0.5:  # glide mode
        cs = {}
        for name, P in paulis.items():
            cs[name] = np.real(np.trace(d @ P)) / 4
        vec = np.array([cs[name] for name in sorted(paulis.keys())])
        glide_diffs.append((k, k+1, vec, cs))

if len(glide_diffs) >= 2:
    print(f"\nCosine similarity between glide-mode Pauli vectors:")
    for i in range(len(glide_diffs)):
        for j in range(i+1, len(glide_diffs)):
            v1 = glide_diffs[i][2]
            v2 = glide_diffs[j][2]
            cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-15)
            p1 = f"{glide_diffs[i][0]}->{glide_diffs[i][1]}"
            p2 = f"{glide_diffs[j][0]}->{glide_diffs[j][1]}"
            print(f"  {p1} vs {p2}: cos = {cos:+.4f}")

print(f"\n{'=' * 60}")
print("DONE")
print("=" * 60)
