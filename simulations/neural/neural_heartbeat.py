#!/usr/bin/env python3
"""
Neural Heartbeat: Wilson-Cowan dynamics around the palindromic midpoint.

Does the E/I balance oscillate? Does the crossing frequency map to
known neural rhythms (Delta 1-4 Hz, Theta 4-8, Alpha 8-13, Gamma 30-80)?

Usage: python neural_heartbeat.py
"""
import numpy as np

def sigmoid(x, a=1.3, theta=4.0):
    return 1.0 / (1.0 + np.exp(np.clip(-a * (x - theta), -500, 500)))


def wilson_cowan_dynamics(N=3, tau_E=8.0, tau_I=18.0,
                           w_EE=16.0, w_EI=12.0, w_IE=15.0, w_II=3.0,
                           w_lat=2.0, P=1.5, Q=0.0,
                           a_E=1.3, a_I=2.0, theta_E=4.0, theta_I=3.7,
                           dt=0.1, t_max=1000.0, noise=0.01):
    """
    Simulate Wilson-Cowan chain with full nonlinear dynamics.
    Returns time series of E and I activity per node.
    Time in milliseconds, dt in ms.
    """
    n_steps = int(t_max / dt)
    E = np.ones(N) * 0.2 + np.random.randn(N) * 0.01
    I = np.ones(N) * 0.1 + np.random.randn(N) * 0.01

    # Store mean E and I over time
    t_series = []
    E_mean_series = []
    I_mean_series = []
    balance_series = []  # E_mean - I_mean

    for step in range(n_steps):
        # Compute derivatives
        dE = np.zeros(N)
        dI = np.zeros(N)

        for i in range(N):
            # Local input
            inp_E = w_EE * E[i] - w_EI * I[i] + P
            inp_I = w_IE * E[i] - w_II * I[i] + Q

            # Lateral coupling
            if i > 0:
                inp_E += w_lat * E[i-1]
            if i < N - 1:
                inp_E += w_lat * E[i+1]

            dE[i] = (-E[i] + sigmoid(inp_E, a_E, theta_E)) / tau_E
            dI[i] = (-I[i] + sigmoid(inp_I, a_I, theta_I)) / tau_I

        # Add small noise (biological systems are driven)
        dE += np.random.randn(N) * noise / tau_E
        dI += np.random.randn(N) * noise / tau_I

        # Euler step
        E = np.clip(E + dt * dE, 0, 1)
        I = np.clip(I + dt * dI, 0, 1)

        # Record every 1 ms
        if step % max(1, int(1.0 / dt)) == 0:
            t_series.append(step * dt)
            E_mean_series.append(np.mean(E))
            I_mean_series.append(np.mean(I))
            balance_series.append(np.mean(E) - np.mean(I))

    return (np.array(t_series), np.array(E_mean_series),
            np.array(I_mean_series), np.array(balance_series))


def count_crossings(series, threshold=0.0):
    """Count zero-crossings (or threshold crossings) in a time series."""
    crossings_up = 0
    crossings_down = 0
    for i in range(1, len(series)):
        if series[i-1] < threshold and series[i] >= threshold:
            crossings_up += 1
        elif series[i-1] >= threshold and series[i] < threshold:
            crossings_down += 1
    return crossings_up, crossings_down


def analyze_frequency(t, series, threshold=0.0):
    """Estimate oscillation frequency from crossing rate."""
    total_time_s = (t[-1] - t[0]) / 1000.0  # ms to seconds
    up, down = count_crossings(series, threshold)
    total = up + down
    if total_time_s > 0 and total > 0:
        freq = total / (2 * total_time_s)  # crossings/2 = full cycles
    else:
        freq = 0
    return freq, up, down


print("=" * 60)
print("NEURAL HEARTBEAT: Wilson-Cowan Dynamics")
print("=" * 60)

# ============================================================
# Step 1: Basic dynamics - does it oscillate?
# ============================================================
print("\n### Step 1: Basic dynamics (N=3, default parameters)")

t, E_m, I_m, balance = wilson_cowan_dynamics(
    N=3, tau_E=8.0, tau_I=18.0, t_max=2000.0, noise=0.02)

# Find the equilibrium balance
eq_balance = np.mean(balance[len(balance)//2:])
print(f"Mean E: {np.mean(E_m):.4f}, Mean I: {np.mean(I_m):.4f}")
print(f"E-I balance mean: {eq_balance:.4f}")
print(f"Balance range: [{np.min(balance):.4f}, {np.max(balance):.4f}]")

freq, up, down = analyze_frequency(t, balance, threshold=eq_balance)
print(f"Crossings of mean: {up} up + {down} down = {up+down}")
print(f"Estimated frequency: {freq:.1f} Hz")

# Show first 100ms of balance
print(f"\nFirst 200ms (E-I balance):")
for i in range(min(200, len(t))):
    if i % 10 == 0:
        b = balance[i]
        bar_len = int(abs(b - eq_balance) * 200)
        if b > eq_balance:
            bar = "+" * min(bar_len, 40)
            print(f"  t={t[i]:6.0f}ms  bal={b:.4f}  E>{bar}")
        else:
            bar = "-" * min(bar_len, 40)
            print(f"  t={t[i]:6.0f}ms  bal={b:.4f}  I>{bar}")

# ============================================================
# Step 2: Coupling strength sweep
# ============================================================
print("\n\n" + "=" * 60)
print("### Step 2: Coupling sweep (frequency vs w_EI)")
print("=" * 60)

print(f"\n{'w_EI':>6s} {'w_IE':>6s} {'Freq(Hz)':>8s} {'Crossings':>10s} {'Amplitude':>10s} {'Band':>10s}")

for w_scale in [0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0]:
    w_EI = 12.0 * w_scale
    w_IE = 15.0 * w_scale

    t, E_m, I_m, balance = wilson_cowan_dynamics(
        N=3, tau_E=8.0, tau_I=18.0,
        w_EI=w_EI, w_IE=w_IE,
        t_max=5000.0, noise=0.02)

    eq = np.mean(balance[len(balance)//2:])
    freq, up, down = analyze_frequency(t, balance, threshold=eq)
    amp = np.std(balance[len(balance)//2:])

    # Classify frequency band
    if freq < 1: band = "sub-Delta"
    elif freq < 4: band = "Delta"
    elif freq < 8: band = "Theta"
    elif freq < 13: band = "Alpha"
    elif freq < 30: band = "Beta"
    elif freq < 80: band = "Gamma"
    else: band = "high-Gamma"

    print(f"{w_EI:6.1f} {w_IE:6.1f} {freq:8.1f} {up+down:10d} {amp:10.4f} {band:>10s}")

# ============================================================
# Step 3: Damping test - sustained or decaying?
# ============================================================
print("\n\n" + "=" * 60)
print("### Step 3: Damping (with and without noise)")
print("=" * 60)

for noise_level in [0.0, 0.01, 0.05]:
    t, E_m, I_m, balance = wilson_cowan_dynamics(
        N=3, tau_E=8.0, tau_I=18.0, t_max=5000.0, noise=noise_level)

    eq = np.mean(balance[len(balance)//2:])

    # Measure amplitude in first 1000ms vs last 1000ms
    first_idx = len(t) // 5
    last_start = 4 * len(t) // 5
    amp_first = np.std(balance[:first_idx])
    amp_last = np.std(balance[last_start:])
    freq, up, down = analyze_frequency(t, balance, threshold=eq)

    damping = "DAMPED" if amp_last < amp_first * 0.5 else "SUSTAINED"

    print(f"  noise={noise_level}: freq={freq:.1f}Hz, "
          f"amp_first={amp_first:.4f}, amp_last={amp_last:.4f}, {damping}")

# ============================================================
# Step 4: tau_E/tau_I ratio effect
# ============================================================
print("\n\n" + "=" * 60)
print("### Step 4: tau ratio effect on frequency")
print("=" * 60)

for tau_I in [8.0, 12.0, 18.0, 30.0, 50.0]:
    t, E_m, I_m, balance = wilson_cowan_dynamics(
        N=3, tau_E=8.0, tau_I=tau_I, t_max=5000.0, noise=0.02)

    eq = np.mean(balance[len(balance)//2:])
    freq, up, down = analyze_frequency(t, balance, threshold=eq)
    amp = np.std(balance[len(balance)//2:])

    print(f"  tau_E=8.0, tau_I={tau_I:5.1f} (ratio={tau_I/8:.1f}): "
          f"freq={freq:.1f}Hz, amp={amp:.4f}, crossings={up+down}")

print("\nDone.")
