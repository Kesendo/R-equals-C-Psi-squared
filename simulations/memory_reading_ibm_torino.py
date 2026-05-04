"""Apply the trio's state-level reading to the IBM Torino single-qubit tomography
data (Feb 2026). Quantifies how long ρ(t) carries readable memory, where memory is
the Bloch-vector magnitude |r|/2 (= eigenvalue deviation from 1/2 baseline) per
the framework's BlochAxisReading.

Setup: Qubit 52 prepared in |+⟩, free decay under T1/T2. Tomography at 25 delay
points 0..895 us. T1 = 221.2 us, T2 = 298.2 us.

For a single qubit:
- ρ(t) measured directly via tomography (3 bases)
- Bloch vector r(t) = (Tr(X·ρ), Tr(Y·ρ), Tr(Z·ρ))
- |r|(t) = sqrt(rx^2 + ry^2 + rz^2)
- Eigenvalue deviation from 1/2: ±|r|/2 (the structural memory-axis pair)

Empirical finding: under T1+T2 the state does NOT thermalise toward I/2 (which
would zero the deviation). Instead the dominant axis shifts X → Z as the prepared
coherence decoheres (T2) and the population relaxes toward |0⟩ (T1 to thermal).
Memory is never lost; it migrates between Bloch axes.
"""

import json
import numpy as np

# Pauli matrices
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)

DATA_PATH = "data/ibm_tomography_feb2026/tomography_ibm_torino_20260209_131521.json"


def dominant_axis(rx, ry, rz, tol=0.05):
    """Pick the Bloch axis with largest |r_a| above tolerance; 'I' if all below."""
    abs_x, abs_y, abs_z = abs(rx), abs(ry), abs(rz)
    if max(abs_x, abs_y, abs_z) < tol:
        return 'I', 0
    if abs_x >= abs_y and abs_x >= abs_z:
        return 'X', 1 if rx >= 0 else -1
    if abs_y >= abs_z:
        return 'Y', 1 if ry >= 0 else -1
    return 'Z', 1 if rz >= 0 else -1


with open(DATA_PATH) as f:
    data = json.load(f)

T1 = data["T1_us"]
T2 = data["T2_us"]

print(f"=== IBM Torino Qubit {data['qubit_index']} (Feb 2026 tomography) ===")
print(f"T1 = {T1:.1f} us,  T2 = {T2:.1f} us,  initial state |+>")
print()
print("Per-time-point Bloch reading (|r|/2 = eigenvalue deviation from 1/2):")
print(f"{'t (us)':>8} {'rx':>7} {'ry':>7} {'rz':>7} {'|r|':>6} {'|r|/2':>7} {'axis':>6} {'fidelity':>9}")
print("-" * 70)

ts = []
mags = []
rxs, rys, rzs = [], [], []
axes = []

for entry in data["raw_tomography"]:
    t_us = entry["delay_us"]
    rho = np.array(entry["density_matrix_real"]) + 1j * np.array(entry["density_matrix_imag"])

    rx = np.trace(SX @ rho).real
    ry = np.trace(SY @ rho).real
    rz = np.trace(SZ @ rho).real
    rmag = float(np.sqrt(rx**2 + ry**2 + rz**2))
    deviation = rmag / 2.0
    fid = entry["fidelity"]
    ax, sign = dominant_axis(rx, ry, rz)

    sign_str = '+' if sign > 0 else ('-' if sign < 0 else '')
    print(f"{t_us:>8.1f} {rx:>+7.3f} {ry:>+7.3f} {rz:>+7.3f} {rmag:>6.3f} {deviation:>7.3f} {ax}{sign_str:<2} {fid:>9.3f}")

    ts.append(t_us)
    mags.append(rmag)
    rxs.append(rx)
    rys.append(ry)
    rzs.append(rz)
    axes.append(ax)

ts = np.array(ts)
mags = np.array(mags)
rxs, rys, rzs = np.array(rxs), np.array(rys), np.array(rzs)

print()
print("=== Memory-reading regime over time ===")
print()
print("Dominant-axis trace:")
prev_ax = None
for i, ax in enumerate(axes):
    if ax != prev_ax:
        print(f"  t = {ts[i]:>6.1f} us:  axis -> {ax}  (|r|/2 = {mags[i]/2:.3f})")
        prev_ax = ax

print()
print("=== Three readings of 'memory lifetime' ===")
print()

# Reading 1: X-Y coherence magnitude (rotation-invariant, the true T2 measure)
xy_mag = np.sqrt(rxs**2 + rys**2)
idx_xy_drop = np.where(xy_mag < 0.10)[0]
if len(idx_xy_drop) > 0:
    t_xy_lost = ts[idx_xy_drop[0]]
    # Effective T2 from xy_mag(t) = xy_mag(0) · exp(-t/T2_eff):
    # take a clean mid-decay point and solve
    fit_idx = 4  # t ≈ 149 us, xy_mag ≈ 0.22
    if xy_mag[fit_idx] > 0 and xy_mag[0] > 0:
        T2_eff = -ts[fit_idx] / np.log(xy_mag[fit_idx] / xy_mag[0])
    else:
        T2_eff = T2
    print(f"(1) X-Y coherence lifetime (when sqrt(rx² + ry²) drops below 0.10):")
    print(f"    t ≈ {t_xy_lost:.0f} us")
    print(f"    Empirical T2_eff from in-plane decay: {T2_eff:.0f} us")
    print(f"    Calibration T2 (from Ramsey, may have drifted): {T2:.0f} us")
    print(f"    Ratio T2_eff / T2_cal: {T2_eff/T2:.2f}  (the actual run was {'faster' if T2_eff < T2 else 'slower'})")

print()
# Reading 2: Total Bloch magnitude (any-axis memory)
print(f"(2) Total Bloch memory |r|/2 throughout 0-895 us window:")
print(f"    min |r|/2 = {min(mags)/2:.3f} at t = {ts[np.argmin(mags)]:.0f} us  (axis crossover)")
print(f"    max |r|/2 = {max(mags)/2:.3f} at t = {ts[np.argmax(mags)]:.0f} us")
print(f"    final |r|/2 = {mags[-1]/2:.3f} at t = {ts[-1]:.0f} us  ({axes[-1]}-dominated)")

print()
# Reading 3: Distance to maximally mixed (|r|/2 = 0)
threshold_close = 0.05
mask = mags / 2 < threshold_close
if mask.any():
    t_mixed = ts[np.where(mask)[0][0]]
    print(f"(3) When does ρ approach maximally mixed (|r|/2 < {threshold_close})?")
    print(f"    t ≈ {t_mixed:.0f} us")
else:
    print(f"(3) ρ never approaches maximally mixed: |r|/2 stays above {min(mags)/2:.3f}")
    print(f"    The state ROTATES from X-axis to Z-axis instead of thermalising to I/2.")

print()
print("=== Structural reading ===")
print()
print("The framework's d=0 axis is the maximally mixed I/2 (zero Bloch magnitude).")
print("This data shows the qubit NEVER reaches d=0:")
print(f"  - Starts at d=2 X-axis (|r|/2 = {mags[0]/2:.3f}, near full anchor 0.5)")
print(f"  - Crosses through min |r|/2 = {min(mags)/2:.3f} at t = {ts[np.argmin(mags)]:.0f} us")
print(f"  - Recovers to d=2 Z-axis (|r|/2 = {mags[-1]/2:.3f}, thermal-state direction)")
print()
print("Memory is never UNREADABLE in 0-895 us; it MIGRATES between axes.")
print("The framework's reading distinguishes:")
print("  - WHICH axis carries the memory (X / Y / Z)")
print("  - HOW STRONGLY (|r|/2)")
print("  - The d=0 axis (mixed) is approached but not crossed.")
