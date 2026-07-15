"""The Givens compiler: the 8x8 Floquet unitary as a NN matchgate schedule, per arm.

WHAT IT DOES. For each arm (a pair of 3-magnon branches a, b), chooses a column
permutation P putting the arm's DIFFERING modes (a triangle b, the cat block) on
contiguous qubits starting at 0 (a-private first, then b-private), the SHARED
modes next, spectators last; then decomposes the permuted single-particle
Floquet-mode matrix into NN SU(2) Givens rotations plus a diagonal phase layer,

    G_m ... G_1 (V P) = D   =>   Gamma(V P) = Gamma(G_1^dag) ... Gamma(G_m^dag) Gamma(D),

m = N(N-1)/2 = 28 for any P. The seed per arm: shared-mode qubits set to |1>,
GHZ over the cat block + X layer, i.e. (|a-bits> + |b-bits>)/sqrt2. The collision
arm (disjoint) has a 6-qubit cat (5 CX); overlapping control arms have SMALLER
cats (4 and 2 qubits: 3 and 1 CX), so the doc's GHZ-6 budget is conservative.

SIGNS. The circuit produces Slater states with columns in QUBIT order; reordering
columns to mode order costs the permutation sign per branch. The compiler reports
sign(a), sign(b); their product shifts the fringe by a constant 0-or-pi offset
(intercept, never slope). Certificates check exactness INCLUDING the sign.

CERTIFICATES (machine zero required):
  C1  circuit|bitstring> = sign * direct Slater state, both branches, all arms.
  C2  full chain (seed + network + M steps + inverse network) reads fringe phase
      M*dPhi + offset on the cat block, zero leakage, all arms, M = 1..8.
Doc: experiments/IBM_F129_RAMSEY_FRINGE.md section 3; consumed by f129_ramsey_7a_gate.py."""

import numpy as np
from itertools import combinations, permutations

N, THETA = 8, 0.5
ARMS = [("A0", (1, 5, 7), (2, 4, 8)),
        ("A1", (1, 5, 7), (2, 5, 6)),
        ("A2", (1, 5, 7), (1, 6, 7))]
DIM = 2 ** N


def one_particle_step(theta):
    U = np.eye(N, dtype=complex)
    for parity in (0, 1):
        L = np.eye(N, dtype=complex)
        for i in range(parity, N - 1, 2):
            B = np.eye(N, dtype=complex)
            B[i, i] = B[i + 1, i + 1] = np.cos(theta)
            B[i, i + 1] = B[i + 1, i] = 1j * np.sin(theta)
            L = B @ L
        U = L @ U
    return U


def floquet(theta):
    w, V = np.linalg.eig(one_particle_step(theta))
    ph = np.angle(w)
    order = np.argsort(-ph)
    return ph[order], V[:, order]


def givens_decompose(A):
    """(rotations, phases) with rotations = [(top row p, 2x2 SU(2) G)] in kill
    order such that G_m ... G_1 A = diag(phases)."""
    A = A.copy().astype(complex)
    n = A.shape[0]
    rots = []
    for j in range(n - 1):
        for i in range(n - 1, j, -1):
            a, b = A[i - 1, j], A[i, j]
            if abs(b) < 1e-15:
                continue
            r = np.hypot(abs(a), abs(b))
            G = np.array([[np.conj(a) / r, np.conj(b) / r],
                          [-b / r, a / r]], dtype=complex)
            A[i - 1:i + 1, :] = G @ A[i - 1:i + 1, :]
            rots.append((i - 1, G))
    phases = np.diag(A).copy()
    assert np.max(np.abs(A - np.diag(phases))) < 1e-12, "not diagonal after sweep"
    return rots, phases


def perm_sign(seq):
    """Parity sign of the permutation sorting seq ascending."""
    seq = list(seq)
    sign = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                sign = -sign
    return sign


class ArmNetwork:
    """Per-arm compiled network: permutation, Givens schedule, seed bitstrings,
    branch signs."""

    def __init__(self, V, a, b):
        self.a, self.b = a, b
        a_priv = [k for k in a if k not in b]
        b_priv = [k for k in b if k not in a]
        shared = [k for k in a if k in b]
        rest = [k for k in range(1, N + 1) if k not in set(a) | set(b)]
        self.mode_of_qubit = a_priv + b_priv + shared + rest   # 1-based modes
        self.cat_size = len(a_priv) + len(b_priv)
        perm = [m - 1 for m in self.mode_of_qubit]
        self.rots, self.phases = givens_decompose(V[:, perm])
        self.qubits_a = sorted(self.mode_of_qubit.index(k) for k in a)
        self.qubits_b = sorted(self.mode_of_qubit.index(k) for k in b)
        self.bits_a = sum(1 << q for q in self.qubits_a)
        self.bits_b = sum(1 << q for q in self.qubits_b)
        # circuit column order (by qubit) vs mode order (sorted triple):
        self.sign_a = perm_sign([self.mode_of_qubit[q] for q in self.qubits_a])
        self.sign_b = perm_sign([self.mode_of_qubit[q] for q in self.qubits_b])


def apply_2x2_number_conserving(state, p, B):
    out = state.copy()
    bp, bq = 1 << p, 1 << (p + 1)
    for x in range(DIM):
        if (x & bp) and not (x & bq):
            y = x ^ bp ^ bq
            out[x] = B[0, 0] * state[x] + B[0, 1] * state[y]
            out[y] = B[1, 0] * state[x] + B[1, 1] * state[y]
    return out


def apply_phase_layer(state, phases):
    out = state.copy()
    for x in range(DIM):
        ph = 1.0 + 0j
        for k in range(N):
            if x & (1 << k):
                ph *= phases[k]
        out[x] = ph * state[x]
    return out


def network_apply(state, net, inverse=False):
    if not inverse:
        state = apply_phase_layer(state, net.phases)
        for p, G in reversed(net.rots):
            state = apply_2x2_number_conserving(state, p, G.conj().T)
    else:
        for p, G in net.rots:
            state = apply_2x2_number_conserving(state, p, G)
        state = apply_phase_layer(state, np.conj(net.phases))
    return state


def step_circuit(s, theta):
    B = np.array([[np.cos(theta), 1j * np.sin(theta)],
                  [1j * np.sin(theta), np.cos(theta)]])
    for parity in (0, 1):
        for i in range(parity, N - 1, 2):
            s = apply_2x2_number_conserving(s, i, B)
    return s


def slater_vec(V, triple):
    cols = V[:, [k - 1 for k in triple]]
    vec = np.zeros(DIM, dtype=complex)
    for c in combinations(range(N), 3):
        idx = sum(1 << i for i in c)
        vec[idx] = np.linalg.det(cols[list(c), :])
    return vec / np.linalg.norm(vec)


def main():
    ph, V = floquet(THETA)
    print("certificate C1: circuit|bitstring> = sign * direct Slater, per branch")
    nets = {}
    for name, a, b in ARMS:
        net = ArmNetwork(V, a, b)
        nets[name] = net
        print(f"  {name} {a}~{b}: qubit->mode {net.mode_of_qubit}, "
              f"cat = {net.cat_size} qubits ({net.cat_size - 1} CX GHZ), "
              f"{len(net.rots)} Givens")
        for triple, bits, sign in ((a, net.bits_a, net.sign_a),
                                   (b, net.bits_b, net.sign_b)):
            state = np.zeros(DIM, dtype=complex)
            state[bits] = 1.0
            state = network_apply(state, net)
            err = np.linalg.norm(state - sign * slater_vec(V, triple))
            print(f"    {triple} sign {sign:+d}: |circuit - sign*direct| = {err:.2e}")

    print("\ncertificate C2: seed + network + M steps + inverse, fringe = M*dPhi + offset")
    for name, a, b in ARMS:
        net = nets[name]
        seed = np.zeros(DIM, dtype=complex)
        seed[net.bits_a] = 1 / np.sqrt(2)
        seed[net.bits_b] = 1 / np.sqrt(2)
        state = network_apply(seed, net)
        dPhi = sum(ph[k - 1] for k in a) - sum(ph[k - 1] for k in b)
        offset = 0.0 if net.sign_a * net.sign_b > 0 else np.pi
        worst_ph, worst_leak = 0.0, 0.0
        for m in range(1, 9):
            state = step_circuit(state, THETA)
            back = network_apply(state, net, inverse=True)
            aA, aB = back[net.bits_a], back[net.bits_b]
            fringe = np.angle(aA * np.conj(aB))
            pred = ((m * dPhi + offset + np.pi) % (2 * np.pi)) - np.pi
            worst_ph = max(worst_ph, abs(((fringe - pred + np.pi) % (2 * np.pi)) - np.pi))
            worst_leak = max(worst_leak, 1.0 - abs(aA) ** 2 - abs(aB) ** 2)
        print(f"  {name}: dPhi = {dPhi:+.6f}, sign offset = {offset:.0f}, "
              f"worst |fringe err| = {worst_ph:.2e}, worst leakage = {worst_leak:.2e}")


if __name__ == "__main__":
    main()
