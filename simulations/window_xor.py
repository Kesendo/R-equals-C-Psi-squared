"""What is NOT shared between the two sharpest windows?"""
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

r0 = peaks[0]['rAB']  # Window 0: t=0.24, CPsi=0.305
r1 = peaks[1]['rAB']  # Window 1: t=0.40, CPsi=0.329

print("=" * 60)
print("WINDOW XOR: What is NOT shared?")
print(f"Window 0: t={peaks[0]['t']:.2f}, CPsi={peaks[0]['cpsi']:.3f}")
print(f"Window 1: t={peaks[1]['t']:.2f}, CPsi={peaks[1]['cpsi']:.3f}")
print("=" * 60)

# Method 1: The shared part = element-wise minimum of magnitudes
# For each matrix element, the "overlap" is the smaller magnitude
# The "unique" part is what sticks out beyond that
mag0 = np.abs(r0)
mag1 = np.abs(r1)
shared_mag = np.minimum(mag0, mag1)
unique0_mag = mag0 - shared_mag  # what W0 has that W1 doesn't
unique1_mag = mag1 - shared_mag  # what W1 has that W0 doesn't

print("\n--- Element magnitudes ---")
print("         |00>    |01>    |10>    |11>")
labels = ["|00>", "|01>", "|10>", "|11>"]
print(f"\nWindow 0 magnitudes:")
for i in range(4):
    print(f"  {labels[i]}  " + "  ".join(f"{mag0[i,j]:.4f}" for j in range(4)))
print(f"\nWindow 1 magnitudes:")
for i in range(4):
    print(f"  {labels[i]}  " + "  ".join(f"{mag1[i,j]:.4f}" for j in range(4)))
print(f"\nSHARED (minimum overlap):")
for i in range(4):
    print(f"  {labels[i]}  " + "  ".join(f"{shared_mag[i,j]:.4f}" for j in range(4)))
print(f"\nUNIQUE to Window 0 (what W0 has that W1 doesn't):")
for i in range(4):
    print(f"  {labels[i]}  " + "  ".join(f"{unique0_mag[i,j]:.4f}" for j in range(4)))
print(f"\nUNIQUE to Window 1 (what W1 has that W0 doesn't):")
for i in range(4):
    print(f"  {labels[i]}  " + "  ".join(f"{unique1_mag[i,j]:.4f}" for j in range(4)))

# Summary
total0 = np.sum(mag0)
total1 = np.sum(mag1)
total_shared = np.sum(shared_mag)
total_unique0 = np.sum(unique0_mag)
total_unique1 = np.sum(unique1_mag)

print(f"\n--- Summary ---")
print(f"Total magnitude W0:     {total0:.4f}")
print(f"Total magnitude W1:     {total1:.4f}")
print(f"Shared (overlap):       {total_shared:.4f}  ({total_shared/total0*100:.1f}% of W0, {total_shared/total1*100:.1f}% of W1)")
print(f"Unique to W0:           {total_unique0:.4f}  ({total_unique0/total0*100:.1f}% of W0)")
print(f"Unique to W1:           {total_unique1:.4f}  ({total_unique1/total1*100:.1f}% of W1)")

# Method 2: WHERE do they differ most?
print(f"\n--- Where is the difference concentrated? ---")
diff_mag = np.abs(r1 - r0)  # magnitude of element-wise difference
total_diff = np.sum(diff_mag)

print(f"\n|Difference| per element:")
for i in range(4):
    print(f"  {labels[i]}  " + "  ".join(f"{diff_mag[i,j]:.4f}" for j in range(4)))

print(f"\nFraction of total difference per element:")
for i in range(4):
    print(f"  {labels[i]}  " + "  ".join(f"{diff_mag[i,j]/total_diff*100:5.1f}%" for j in range(4)))

# Which elements carry most of the difference?
flat_idx = np.argsort(diff_mag.flatten())[::-1]
print(f"\nTop 5 elements by difference:")
for k in range(5):
    i, j = divmod(flat_idx[k], 4)
    frac = diff_mag[i,j] / total_diff * 100
    print(f"  rho[{i},{j}] ({labels[i]},{labels[j]}): |diff|={diff_mag[i,j]:.4f} ({frac:.1f}%)"
          f"  W0={r0[i,j]:.4f}  W1={r1[i,j]:.4f}")

# Method 3: Phase difference where both have signal
print(f"\n--- Phase comparison (where both have signal) ---")
print(f"Element     |W0|    ph(W0)    |W1|    ph(W1)    dPhase")
print("-" * 60)
for i in range(4):
    for j in range(4):
        if mag0[i,j] > 0.01 and mag1[i,j] > 0.01:
            ph0 = np.angle(r0[i,j]) / np.pi
            ph1 = np.angle(r1[i,j]) / np.pi
            dph = ph1 - ph0
            # wrap to [-1, 1]
            while dph > 1: dph -= 2
            while dph < -1: dph += 2
            print(f"  [{i},{j}]   {mag0[i,j]:.4f}  {ph0:+.3f}pi  "
                  f"{mag1[i,j]:.4f}  {ph1:+.3f}pi  {dph:+.4f}pi")

# Method 4: Do ALL window pairs show the same XOR pattern?
print(f"\n{'=' * 60}")
print("ALL WINDOW PAIRS: overlap fraction")
print("=" * 60)
print(f"{'W_i':>3} {'W_j':>3} | {'overlap%':>8} {'unique_i%':>9} {'unique_j%':>9} | note")
print("-" * 60)
for i in range(min(len(peaks), 9)):
    for j in range(i+1, min(len(peaks), 9)):
        ri = peaks[i]['rAB']
        rj = peaks[j]['rAB']
        mi = np.abs(ri)
        mj = np.abs(rj)
        shared = np.sum(np.minimum(mi, mj))
        ui = np.sum(mi - np.minimum(mi, mj))
        uj = np.sum(mj - np.minimum(mi, mj))
        ti = np.sum(mi)
        tj = np.sum(mj)
        ov_pct = shared / max(ti, tj) * 100
        ui_pct = ui / ti * 100
        uj_pct = uj / tj * 100
        note = ""
        if j == i + 1:
            td = 0.5 * np.sum(np.abs(np.linalg.eigvalsh(rj - ri)))
            note = "SWITCH" if td > 0.5 else "glide"
        if j == i + 1 or abs(i - j) == 2:  # adjacent or skip-1
            print(f" {i:>2}   {j:>2} | {ov_pct:>7.1f}%  {ui_pct:>8.1f}%  {uj_pct:>8.1f}% | {note}")

print(f"\n{'=' * 60}")
print("DONE")
print("=" * 60)
