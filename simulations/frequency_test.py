"""
FREQUENCY TEST: Does f_dom = (J_SA + J_SB + sqrt(J_SA^2 - J_SA*J_SB + J_SB^2)) / pi
match the actual oscillation frequency of the CΨ windows?
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import star_topology_v3 as gpt

J_SA, J_SB = 1.0, 2.0
gamma = 0.05

# GPT's formula
f_predicted = (J_SA + J_SB + np.sqrt(J_SA**2 - J_SA*J_SB + J_SB**2)) / np.pi
# Our earlier finding
f_jtotal = (J_SA + J_SB) / 2

print("=" * 60)
print("FREQUENCY TEST")
print(f"J_SA={J_SA}, J_SB={J_SB}, gamma={gamma}")
print("=" * 60)
print(f"\nPredictions:")
print(f"  GPT formula: f_dom = (J_SA+J_SB+sqrt(J^2-JJ+J^2))/pi = {f_predicted:.4f}")
print(f"  Our earlier: f = J_total/2 = {f_jtotal:.4f}")

# Evolve and collect CΨ_AB time series
H = gpt.star_hamiltonian_n(n_observers=2, J_SA=J_SA, J_SB=J_SB)
L_ops = gpt.dephasing_ops_n([gamma]*3)
L_zero = gpt.dephasing_ops_n([0.0]*3)
psi = np.kron(gpt.bell_phi_plus(), gpt.plus_state())
rho_n = gpt.density_from_statevector(psi)
rho_u = rho_n.copy()

dt = 0.005
times = []
cpsi_noisy = []
cpsi_unitary = []

for step in range(4001):  # t=0 to 20
    t = step * dt
    if step % 2 == 0:  # sample every 0.01
        rAB_n = gpt.partial_trace_keep(rho_n, [1,2], 3)
        rAB_u = gpt.partial_trace_keep(rho_u, [1,2], 3)
        cn = gpt.concurrence_two_qubit(rAB_n)
        pn = gpt.psi_norm(rAB_n)
        cu = gpt.concurrence_two_qubit(rAB_u)
        pu = gpt.psi_norm(rAB_u)
        times.append(t)
        cpsi_noisy.append(cn * pn)
        cpsi_unitary.append(cu * pu)
    if step < 4000:
        rho_n = gpt.rk4_step(rho_n, H, L_ops, dt)
        rho_u = gpt.rk4_step(rho_u, H, L_zero, dt)

times = np.array(times)
cpsi_n = np.array(cpsi_noisy)
cpsi_u = np.array(cpsi_unitary)

# FFT on both noisy and unitary
dt_sample = times[1] - times[0]
N = len(times)

# Remove mean
cpsi_u_centered = cpsi_u - np.mean(cpsi_u)
cpsi_n_centered = cpsi_n - np.mean(cpsi_n)

freqs = np.fft.rfftfreq(N, d=dt_sample)
fft_u = np.abs(np.fft.rfft(cpsi_u_centered))
fft_n = np.abs(np.fft.rfft(cpsi_n_centered))

# Find peaks
# Unitary
mask = freqs > 0.1  # ignore DC
top_u_idx = np.argsort(fft_u[mask])[::-1][:5]
top_u_freqs = freqs[mask][top_u_idx]
top_u_amps = fft_u[mask][top_u_idx]

# Noisy
top_n_idx = np.argsort(fft_n[mask])[::-1][:5]
top_n_freqs = freqs[mask][top_n_idx]
top_n_amps = fft_n[mask][top_n_idx]

print(f"\n--- FFT Results (unitary, gamma=0) ---")
print(f"  Top 5 frequencies:")
for i in range(5):
    match_gpt = abs(top_u_freqs[i] - f_predicted) / f_predicted * 100
    match_our = abs(top_u_freqs[i] - f_jtotal) / f_jtotal * 100
    print(f"  f={top_u_freqs[i]:.4f}  amp={top_u_amps[i]:.3f}"
          f"  (GPT:{match_gpt:.1f}% off, Ours:{match_our:.1f}% off)")

print(f"\n--- FFT Results (noisy, gamma={gamma}) ---")
print(f"  Top 5 frequencies:")
for i in range(5):
    match_gpt = abs(top_n_freqs[i] - f_predicted) / f_predicted * 100
    match_our = abs(top_n_freqs[i] - f_jtotal) / f_jtotal * 100
    print(f"  f={top_n_freqs[i]:.4f}  amp={top_n_amps[i]:.3f}"
          f"  (GPT:{match_gpt:.1f}% off, Ours:{match_our:.1f}% off)")

# Direct comparison
f_dominant_u = top_u_freqs[0]
f_dominant_n = top_n_freqs[0]
print(f"\n--- VERDICT ---")
print(f"  Dominant unitary freq: {f_dominant_u:.4f}")
print(f"  Dominant noisy freq:   {f_dominant_n:.4f}")
print(f"  GPT prediction:        {f_predicted:.4f}")
print(f"  Our prediction:        {f_jtotal:.4f}")
print(f"")
err_gpt_u = abs(f_dominant_u - f_predicted) / f_predicted * 100
err_our_u = abs(f_dominant_u - f_jtotal) / f_jtotal * 100
err_gpt_n = abs(f_dominant_n - f_predicted) / f_predicted * 100
err_our_n = abs(f_dominant_n - f_jtotal) / f_jtotal * 100
print(f"  Unitary: GPT error={err_gpt_u:.1f}%, Our error={err_our_u:.1f}%")
print(f"  Noisy:   GPT error={err_gpt_n:.1f}%, Our error={err_our_n:.1f}%")
print(f"")
if err_gpt_u < err_our_u:
    print(f"  -> GPT formula WINS for unitary")
else:
    print(f"  -> Our formula WINS for unitary")
if err_gpt_n < err_our_n:
    print(f"  -> GPT formula WINS for noisy")
else:
    print(f"  -> Our formula WINS for noisy")

# Test with OTHER J values to see which formula generalizes
print(f"\n{'=' * 60}")
print("GENERALIZATION TEST: Multiple J values")
print("=" * 60)

test_cases = [
    (0.5, 1.0),
    (1.0, 1.0),
    (1.0, 2.0),
    (1.0, 3.0),
    (2.0, 4.0),
]

print(f"\n{'J_SA':>5} {'J_SB':>5} | {'f_meas':>7} {'f_GPT':>7} {'f_ours':>7} | {'GPT%':>6} {'ours%':>6} | winner")
print("-" * 70)

for jsa, jsb in test_cases:
    H2 = gpt.star_hamiltonian_n(n_observers=2, J_SA=jsa, J_SB=jsb)
    L2 = gpt.dephasing_ops_n([0.0]*3)  # unitary for clean frequency
    rho2 = gpt.density_from_statevector(np.kron(gpt.bell_phi_plus(), gpt.plus_state()))
    
    t2, c2 = [], []
    for step in range(4001):
        t = step * dt
        if step % 2 == 0:
            rAB2 = gpt.partial_trace_keep(rho2, [1,2], 3)
            cc = gpt.concurrence_two_qubit(rAB2)
            pp = gpt.psi_norm(rAB2)
            t2.append(t)
            c2.append(cc * pp)
        if step < 4000:
            rho2 = gpt.rk4_step(rho2, H2, L2, dt)
    
    c2 = np.array(c2) - np.mean(c2)
    fft2 = np.abs(np.fft.rfft(c2))
    freqs2 = np.fft.rfftfreq(len(c2), d=dt_sample)
    mask2 = freqs2 > 0.1
    f_meas = freqs2[mask2][np.argmax(fft2[mask2])]
    
    f_gpt = (jsa + jsb + np.sqrt(jsa**2 - jsa*jsb + jsb**2)) / np.pi
    f_our = (jsa + jsb) / 2
    
    err_g = abs(f_meas - f_gpt) / f_meas * 100 if f_meas > 0 else 999
    err_o = abs(f_meas - f_our) / f_meas * 100 if f_meas > 0 else 999
    winner = "GPT" if err_g < err_o else "OURS" if err_o < err_g else "TIE"
    
    print(f"{jsa:>5.1f} {jsb:>5.1f} | {f_meas:>7.4f} {f_gpt:>7.4f} {f_our:>7.4f} | {err_g:>5.1f}% {err_o:>5.1f}% | {winner}")

print(f"\n{'=' * 60}")
print("FREQUENCY TEST COMPLETE")
print("=" * 60)
