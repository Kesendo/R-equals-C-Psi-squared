#!/usr/bin/env python3
"""What happens AFTER the EP collapse? Span the EP, below and above.

At the EP the two modes coalesce (defective, Jordan block, angle -> 0, Petermann K -> inf). The
question is what lies on the far side. Read the F86a 2-level across x = Q/Q_EP from below to above:
the angle between the two eigenvectors, the Petermann factor (sensitivity), the rotation omega,
and the clock angle theta.

Reading: below the EP, overdamped , two real decay modes (two forgetting channels), angle
collapsing toward 0, K climbing. At the EP, the pinch. ABOVE the EP, the system HEALS: the two
modes re-separate, but now as a complex-conjugate ROTATING pair (omega opens, theta climbs), the
Petermann sensitivity comes back DOWN off the singularity, and the eigenvectors re-separate toward
orthogonal. The lost degree of freedom is restored, transformed , decay-split becomes
rotation-split. In the memory reading: after the collapse the surplus is remembered (rotation), the
coherent phase persisting on top of the now-fixed decay 4*gamma0*k.
"""
import numpy as np

G0, K_CH, G_EFF = 1.0, 1, 4.0 / 3.0
Q_EP = 2.0 / G_EFF


def read(x):
    Q = x * Q_EP
    J = Q * G0
    L = np.array([[-2 * G0 * (2 * K_CH - 1), 1j * J * G_EFF],
                  [1j * J * G_EFF,           -2 * G0 * (2 * K_CH + 1)]], dtype=complex)
    w, V = np.linalg.eig(L)
    v0, v1 = V[:, 0], V[:, 1]
    cosang = abs(np.vdot(v0, v1)) / (np.linalg.norm(v0) * np.linalg.norm(v1))
    angle = np.degrees(np.arccos(min(cosang, 1.0)))
    Vinv = np.linalg.inv(V)
    petermann = max(float(np.linalg.norm(Vinv[n, :]) ** 2) for n in range(2))
    decay, omega = -w.real, np.abs(w.imag)
    i = int(np.argmin(decay))
    theta = np.degrees(np.arctan2(omega[i], decay[i]))
    return angle, petermann, omega[i], theta


def main():
    print(f"Q_EP = {Q_EP:.3f} (g_eff={G_EFF:.4f}, gamma0={G0}, k={K_CH});  the EP is at x=1")
    print(f"  {'x=Q/Q_EP':>9}  {'angle':>7}  {'Petermann':>10}  {'omega':>7}  {'theta':>7}   regime")
    for x in [0.5, 0.9, 0.99, 1.001, 1.1, 1.5, 2.197, 3.0, 6.0]:
        a, p, om, th = read(x)
        reg = "overdamped  (forgetting)" if x < 1.0 else "underdamped (remembering)"
        tag = "  <- x_peak (resonance peak)" if abs(x - 2.197) < 1e-2 else ""
        print(f"  {x:9.3f}  {a:6.2f}°  {p:10.4g}  {om:7.3f}  {th:6.1f}°   {reg}{tag}")
    print("\n  Below x=1: two real decay modes, angle -> 0, K climbing (approaching the pinch).")
    print("  At x=1: the EP , modes coalesce, defective, K -> inf.")
    print("  Above x=1: HEALED into a rotating conjugate pair , angle re-grows, K comes back down,")
    print("  omega/theta open. The resonance peaks at x_peak=2.197, past the collapse, in the")
    print("  remembering regime; the decay stays pinned at 4*gamma0*k throughout.")


if __name__ == "__main__":
    main()
