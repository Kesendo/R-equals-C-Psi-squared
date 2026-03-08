"""
STRUCTURAL CALCULATIONS (GPT original, March 2026)
Verifies: X*X symmetry, Bell-sector structure, S-coherence gating
Path adjusted from GPT sandbox to local repo structure.
"""
import math
import os
import sys
import numpy as np

# Adjust path for local use
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

import star_topology_v3 as st


def simulate_store(gamma=0.05, J_SA=1.0, J_SB=2.0, J_AB=0.0, dt=0.005, t_max=5.0):
    n_observers = 2
    n_qubits = 3
    H = st.star_hamiltonian_n(n_observers=n_observers, J_SA=J_SA, J_SB=J_SB, J_AB=J_AB)
    L_ops = st.dephasing_ops_n([gamma] * n_qubits)
    rho = st.bell_sa_plus_rest(n_observers)
    times, rho_ABs, rho_Ss, cpsis = [], [], [], []
    steps = int(round(t_max / dt))
    for step in range(steps + 1):
        t = step * dt
        rho_AB = st.partial_trace_keep(rho, [1, 2], n_qubits)
        rho_S = st.partial_trace_keep(rho, [0], n_qubits)
        times.append(t)
        rho_ABs.append(rho_AB.copy())
        rho_Ss.append(rho_S.copy())
        cpsis.append(st.pair_metrics(rho_AB)["cpsi"])
        if step < steps:
            rho = st.rk4_step(rho, H, L_ops, dt)
    return np.array(times), rho_ABs, rho_Ss, np.array(cpsis)


def local_maxima(vals):
    idx = []
    for i in range(1, len(vals) - 1):
        if vals[i - 1] < vals[i] and vals[i] >= vals[i + 1]:
            idx.append(i)
    return idx


BELL_STATES = {
    "Phi+": np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2),
    "Phi-": np.array([1, 0, 0, -1], dtype=complex) / np.sqrt(2),
    "Psi+": np.array([0, 1, 1, 0], dtype=complex) / np.sqrt(2),
    "Psi-": np.array([0, 1, -1, 0], dtype=complex) / np.sqrt(2),
}


def bell_fidelities(rho):
    out = {}
    for name, ket in BELL_STATES.items():
        out[name] = float(np.real(np.vdot(ket, rho @ ket)))
    return out


def shannon_entropy(probs_dict):
    p = np.array(list(probs_dict.values()), dtype=float)
    p = np.clip(p, 1e-15, 1.0)
    return float(-(p * np.log2(p)).sum())


def mutual_information_binned(x, y, bins=8):
    x = np.asarray(x)
    y = np.asarray(y)
    x_edges = np.quantile(x, np.linspace(0, 1, bins + 1))
    y_edges = np.quantile(y, np.linspace(0, 1, bins + 1))
    def strictly_increasing(edges):
        out = [edges[0]]
        for e in edges[1:]:
            if e <= out[-1]:
                e = out[-1] + 1e-12
            out.append(e)
        return np.array(out)
    x_edges = strictly_increasing(x_edges)
    y_edges = strictly_increasing(y_edges)
    xi = np.clip(np.digitize(x, x_edges[1:-1], right=False), 0, bins - 1)
    yi = np.clip(np.digitize(y, y_edges[1:-1], right=False), 0, bins - 1)
    joint = np.zeros((bins, bins))
    for a, b in zip(xi, yi):
        joint[a, b] += 1
    joint /= joint.sum()
    px = joint.sum(axis=1, keepdims=True)
    py = joint.sum(axis=0, keepdims=True)
    mi = 0.0
    for i in range(bins):
        for j in range(bins):
            p = joint[i, j]
            if p > 0 and px[i, 0] > 0 and py[0, j] > 0:
                mi += p * np.log2(p / (px[i, 0] * py[0, j]))
    return float(mi)


def main():
    times, rho_ABs, rho_Ss, cpsis = simulate_store(gamma=0.05)
    times_u, rho_ABs_u, rho_Ss_u, cpsis_u = simulate_store(gamma=0.0)
    peaks = local_maxima(cpsis)
    rows = []
    for w, idx in enumerate(peaks[:8]):
        rho = rho_ABs[idx]
        f = bell_fidelities(rho)
        rows.append({
            "window": w,
            "t_peak": round(float(times[idx]), 3),
            "cpsi_AB": float(cpsis[idx]),
            "bell_dominant": max(f, key=f.get),
            "Phi+": f["Phi+"],
            "Phi-": f["Phi-"],
            "Psi+": f["Psi+"],
            "Psi-": f["Psi-"],
            "bell_entropy_bits": shannon_entropy(f),
            "ab_purity": st.purity(rho),
            "S_coherence": st.psi_norm(rho_Ss[idx]),
            "rho03_phase_rad": float(np.angle(rho[0, 3])),
            "rho12_phase_rad": float(np.angle(rho[1, 2])),
        })
    for r in rows:
        print(r)
    s_coh = np.array([st.psi_norm(r) for r in rho_Ss])
    print("\nPearson corr(CPsi_AB, S coherence):", np.corrcoef(cpsis, s_coh)[0, 1])
    print("MI_8bins(CPsi_AB, S coherence):", mutual_information_binned(cpsis, s_coh, bins=8))


if __name__ == "__main__":
    main()
