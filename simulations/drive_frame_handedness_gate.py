"""The handedness of the virtual phase in the crosstalk letter (experiments/DRIVE_FRAME_HARDWARE_PREDICTION.md, the opening paragraph and §2).

Two constructions of one two-qubit circuit with random virtual rz angles on both qubits and sx pulses, plus a
crosstalk kick on qubit j at the physical drive phase of qubit k's pulse:
  PHYSICAL: no rz unitaries; every pulse on qubit q is a rotation about the axis at angle a_q = -TH_q (the deferred
            virtual phase: rz(theta) = diag(e^{-i theta/2}, e^{i theta/2}) deferred past a pulse turns the pulse's
            axis by -theta), and the kick on j is a small rotation about the axis at angle -TH_k + delta;
  REPLAY:   rz applied as unitaries at their instants, pulses about X, and the kick on j about the axis at angle
            s * (TH_k - TH_j) + delta with s = +1 (the fit script's convention) or s = -1.
Both end before a Z measurement on both qubits (the trailing rz are inert), so the Z-basis probabilities must agree
exactly for the right sign and not for the wrong one, at every random draw. Prints the sign that agrees.
"""
import numpy as np

rng = np.random.default_rng(20261029)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2)


def rot(axis_angle, theta):
    n = np.cos(axis_angle) * X + np.sin(axis_angle) * Y
    return np.cos(theta / 2) * I2 - 1j * np.sin(theta / 2) * n


def rz(t):
    return np.diag([np.exp(-1j * t / 2), np.exp(1j * t / 2)])


def on(q, M):
    return np.kron(M, I2) if q == 0 else np.kron(I2, M)


def run(ops, sign, eps, delta, physical):
    """ops: list of ('rz', q, angle) | ('sx', q). The kick accompanies every sx: on the other qubit, angle eps * pi/2."""
    psi = np.zeros(4, complex)
    psi[0] = 1
    TH = [0.0, 0.0]
    for op in ops:
        if op[0] == "rz":
            TH[op[1]] += op[2]
            if not physical:
                psi = on(op[1], rz(op[2])) @ psi
        else:
            k = op[1]
            j = 1 - k
            if physical:
                U = on(k, rot(-TH[k], np.pi / 2)) @ on(j, rot(-TH[k] + delta, eps * np.pi / 2))
            else:
                U = on(k, rot(0.0, np.pi / 2)) @ on(j, rot(sign * (TH[k] - TH[j]) + delta, eps * np.pi / 2))
            psi = U @ psi
    return np.abs(psi) ** 2


def main():
    verdict = {+1: 0, -1: 0}
    for trial in range(200):
        ops = []
        for _ in range(rng.integers(3, 9)):
            q = int(rng.integers(0, 2))
            if rng.random() < 0.5:
                ops.append(("rz", q, float(rng.uniform(-np.pi, np.pi))))
            else:
                ops.append(("sx", q))
        eps = float(rng.uniform(0.005, 0.05))
        delta = float(rng.uniform(-np.pi, np.pi))
        p_phys = run(ops, 0, eps, delta, True)
        for s in (+1, -1):
            if np.allclose(run(ops, s, eps, delta, False), p_phys, atol=1e-12, rtol=0):
                verdict[s] += 1
    print("trials where the replay construction equals the physical one: sign +1 (TH_k - TH_j):", verdict[+1],
          " sign -1 (TH_j - TH_k):", verdict[-1], " of 200")
    right = [s for s, v in verdict.items() if v == 200]
    print("handedness gate:", ("the crosstalk angle in the replay frame is %s + delta" % ("TH_k - TH_j" if right == [1] else "TH_j - TH_k")) if len(right) == 1 else "UNDECIDED")


if __name__ == "__main__":
    main()
