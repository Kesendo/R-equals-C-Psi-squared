"""F64 on non-chain topologies: test gamma_eff = gamma_B |a_B|^2 beyond chains.

F64 states: for single-excitation modes under single-site Z-dephasing at site B,
the Liouvillian decay rate is alpha_k = 2 gamma_B |a_B(psi_k)|^2, where a_B is the
B-site amplitude of the Hamiltonian eigenvector psi_k.

Currently verified: N=3 chain (1.8% max err), N=4 chain (9 configs, 1.000 +/- 0.0003).

This script tests F64 on:
  - chain (control, must match)
  - star (one center, four leaves)
  - ring (five-cycle)
  - complete (all 10 pairs connected)
  - Y-junction / tree (center with three branches)
  - at N = 5
  - for XX+YY (pure hopping) AND Heisenberg (XX+YY+ZZ)
  - with Z-dephasing at each of the N possible sites B

For each (topology, Hamiltonian, B) combination:
  1. Build single-excitation Hamiltonian H^(1) (N x N, tridiagonal for chain,
     general symmetric for other graphs).
  2. Diagonalize: eigenvalues E_k, eigenvectors psi_k.
  3. Compute a_B(psi_k) = psi_k[B] and predicted alpha_k^F64 = 2 gamma |a_B|^2.
  4. Build single-excitation coherence Liouvillian L^(coh) = -i H^(1) - 2 gamma P_B,
     where P_B is the projector onto site B (|1_B><vac| in superoperator sense).
  5. Diagonalize L^(coh): measured alpha_k^meas = -Re(eigenvalues).
  6. Match measured to predicted pairwise, compute max relative error.

If F64 generalises cleanly to non-chain topologies: max relative error should be
at the perturbative-correction scale O(gamma/J), same as for chains.
"""
import sys
import itertools
import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


# ======================================================================
# Topology bond generators
# ======================================================================

def chain_bonds(N):
    return [(i, i + 1) for i in range(N - 1)]


def star_bonds(N):
    # center = 0, leaves = 1..N-1
    return [(0, i) for i in range(1, N)]


def ring_bonds(N):
    return [(i, (i + 1) % N) for i in range(N)]


def complete_bonds(N):
    return [(i, j) for i, j in itertools.combinations(range(N), 2)]


def y_tree_bonds(N):
    # N=5: center = 0, branches at 1, 2, 3 (each length 1), plus one extra edge
    # Use: 0-1, 0-2, 0-3, 1-4 for Y with one extended branch
    if N == 5:
        return [(0, 1), (0, 2), (0, 3), (1, 4)]
    # N=7: similar but longer
    if N == 7:
        return [(0, 1), (0, 2), (0, 3), (1, 4), (2, 5), (3, 6)]
    raise NotImplementedError(f"y_tree_bonds: N={N} not implemented")


TOPOLOGIES = {
    "chain": chain_bonds,
    "star": star_bonds,
    "ring": ring_bonds,
    "complete": complete_bonds,
    "y_tree": y_tree_bonds,
}


# ======================================================================
# Single-excitation Hamiltonian
# ======================================================================

def single_excitation_H(N, bonds, J=1.0, ham_type="xy"):
    """Build NxN single-excitation Hamiltonian for given bond list.

    J may be scalar (uniform) or dict {(a,b): J_value} (non-uniform).
    Non-uniform J values should be given for each bond; missing bonds take J=1.

    ham_type = "xy"   -> XX+YY only: H_ij = J_bond if (i,j) in bonds, 0 else
    ham_type = "heis" -> Heisenberg (XX+YY+ZZ): adds diagonal from ZZ terms
    """
    def bond_J(a, b):
        if isinstance(J, dict):
            key = (min(a, b), max(a, b))
            return J.get(key, 1.0)
        return float(J)

    H = np.zeros((N, N), dtype=complex)
    # Off-diagonal: XX+YY hopping, per-bond J
    for (a, b) in bonds:
        jb = bond_J(a, b)
        H[a, b] += jb
        H[b, a] += jb
    if ham_type == "xy":
        return H
    if ham_type == "heis":
        # Diagonal drift: sum over bonds (a,b) of -2 J_bond if j in {a,b} else 0
        # = -2 * sum_{b: j adjacent to b} J_b relative to vacuum energy
        for j in range(N):
            drift = 0.0
            for (a, b) in bonds:
                if j == a or j == b:
                    drift += -2.0 * bond_J(a, b)
            H[j, j] += drift
        return H
    raise ValueError(f"ham_type must be 'xy' or 'heis', got {ham_type}")


# ======================================================================
# Single-excitation coherence Liouvillian
# ======================================================================

def single_excitation_coherence_L(N, bonds, B, gamma, J=1.0, ham_type="xy"):  # J may be scalar or dict
    """NxN Liouvillian acting on single-excitation coherences |1_j><vac|.

    L(|1_j><vac|)_i = -i H^(1)_ij + gamma (<1_i|Z_B|1_j><vac|Z_B|vac> - <ij>)
                    = -i H^(1)_ij - 2 gamma delta_{iB} delta_{jB}

    So L = -i H^(1) - 2 gamma P_B, where (P_B)_ij = delta_iB delta_jB.
    """
    H = single_excitation_H(N, bonds, J, ham_type)
    L = -1j * H
    L[B, B] += -2.0 * gamma
    return L


# ======================================================================
# F64 prediction
# ======================================================================

def f64_predict(H, B, gamma, deg_tol=1e-8):
    """Predicted alpha_k^F64 using degenerate first-order perturbation theory.

    For each H-eigenvalue group (degenerate subspace), diagonalise the local
    projector P_B restricted to that subspace to get the "correct" a_B values.
    Then alpha_k = 2 gamma * (eigenvalue of P_B in that basis).

    For non-degenerate H, reduces to the naive formula.
    """
    evals, evecs = np.linalg.eigh(H)
    N = len(evals)
    predicted = np.zeros(N)
    idx = 0
    while idx < N:
        # Find degenerate group starting at idx
        jdx = idx + 1
        while jdx < N and abs(evals[jdx] - evals[idx]) < deg_tol:
            jdx += 1
        # Subspace spanned by evecs[:, idx:jdx]
        V_sub = evecs[:, idx:jdx]  # N x dim_sub
        # Restricted P_B: P_B has only entry at (B, B). In subspace basis:
        # (V_sub^† P_B V_sub)_ab = V_sub[B, a].conj() * V_sub[B, b]
        P_restricted = np.outer(V_sub[B, :].conj(), V_sub[B, :])
        # Eigenvalues of P_restricted are the degenerate-PT-correct |a_B|^2
        pt_evals = np.linalg.eigvalsh(P_restricted)
        pt_evals = np.clip(np.real(pt_evals), 0.0, None)
        predicted[idx:jdx] = 2.0 * gamma * pt_evals
        idx = jdx
    return predicted, evals


def measure_decay_rates(L):
    """Extract alpha_k^meas = -Re(eigenvalues) from coherence Liouvillian."""
    evals = np.linalg.eigvals(L)
    alpha = -np.real(evals)
    alpha_sorted = np.sort(alpha)
    return alpha_sorted


# ======================================================================
# Matching predicted to measured
# ======================================================================

def compare(predicted, measured):
    """Sort both, compute pairwise relative errors.

    Predicted has N entries (one per H eigenstate); measured has N entries
    (one per L eigenvalue). They should match in sorted order for small gamma.
    """
    pred_sorted = np.sort(predicted)
    meas_sorted = np.sort(measured)
    abs_err = np.abs(pred_sorted - meas_sorted)
    rel_err = np.where(pred_sorted > 1e-12, abs_err / pred_sorted, abs_err)
    return pred_sorted, meas_sorted, rel_err


# ======================================================================
# Main scan
# ======================================================================

def main():
    gamma = 0.01  # small enough that first-order F64 should be accurate
    ham_types = ["xy", "heis"]

    # =========================================
    # Part 1: uniform J (control, verified earlier today)
    # =========================================
    print("=" * 74)
    print("PART 1: UNIFORM J (control, already verified)")
    print("=" * 74)
    for N in (5, 7):
        J = 1.0
        print(f"\nF64 test: alpha_k = 2 gamma |a_B(psi_k)|^2, uniform J = {J}")
        print(f"N = {N}, gamma = {gamma}, gamma/J = {gamma/J}")
        _run_scan(N, J, gamma, ham_types)

    # =========================================
    # Part 2: non-uniform J (new: EQ-015 final sub-case)
    # =========================================
    print()
    print("=" * 74)
    print("PART 2: NON-UNIFORM J (F64 non-uniform extension)")
    print("=" * 74)
    np.random.seed(42)  # reproducible
    for N in (5, 7):
        for trial in range(3):
            # Build random J per bond in [0.5, 1.5]
            print(f"\nN = {N}, trial {trial+1}:  random J per bond in [0.5, 1.5]")
            for topo_name, bond_fn in TOPOLOGIES.items():
                bonds = bond_fn(N)
                j_dict = {(min(a, b), max(a, b)): float(np.random.uniform(0.5, 1.5))
                          for (a, b) in bonds}
                # Report range of this trial's J values
                j_vals = sorted(j_dict.values())
                print(f"  {topo_name:10s} J range [{j_vals[0]:.3f}, {j_vals[-1]:.3f}]")
                _run_scan_single_topo(N, j_dict, gamma, ham_types, topo_name, bond_fn)


def _run_scan_single_topo(N, J, gamma, ham_types, topo_name, bond_fn):
    """Scan F64 for a single (topology, J-profile) across all B and both H types."""
    bonds = bond_fn(N)
    j_eff = max(J.values()) if isinstance(J, dict) else J
    for ham_type in ham_types:
        worst_rel = 0.0
        worst_abs = 0.0
        for B in range(N):
            H = single_excitation_H(N, bonds, J, ham_type)
            predicted, energies = f64_predict(H, B, gamma)
            L = single_excitation_coherence_L(N, bonds, B, gamma, J, ham_type)
            measured = measure_decay_rates(L)
            pred_sorted, meas_sorted, rel_err_arr = compare(predicted, measured)
            abs_err = np.abs(pred_sorted - meas_sorted)
            max_abs = np.max(abs_err)
            threshold = gamma * 0.01
            nonzero_mask = pred_sorted > threshold
            max_rel = np.max(rel_err_arr[nonzero_mask]) if np.any(nonzero_mask) else 0.0
            if max_rel > worst_rel:
                worst_rel = max_rel
            if max_abs > worst_abs:
                worst_abs = max_abs
        status = "OK" if worst_rel < 0.05 and worst_abs < gamma * 0.05 else \
                 ("CHECK" if worst_rel < 0.20 else "FAIL")
        print(f"    {ham_type:5s}  worst over B: rel = {worst_rel:.4f}, abs = {worst_abs:.4g}  [{status}]")


def _run_scan(N, J, gamma, ham_types):

    for ham_type in ham_types:
        print(f"=" * 72)
        print(f"Hamiltonian type: {ham_type}")
        print(f"=" * 72)
        for topo_name, bond_fn in TOPOLOGIES.items():
            bonds = bond_fn(N)
            for B in range(N):
                H = single_excitation_H(N, bonds, J, ham_type)
                predicted, energies = f64_predict(H, B, gamma)
                L = single_excitation_coherence_L(N, bonds, B, gamma, J, ham_type)
                measured = measure_decay_rates(L)
                pred_sorted, meas_sorted, rel_err_arr = compare(predicted, measured)

                # For non-matching entries (zero predicted), use absolute instead of relative
                abs_err = np.abs(pred_sorted - meas_sorted)
                max_abs = np.max(abs_err)
                # Relative tolerance: max_rel where predicted > threshold
                threshold = gamma * 0.01  # 1% of gamma scale
                nonzero_mask = pred_sorted > threshold
                max_rel = np.max(rel_err_arr[nonzero_mask]) if np.any(nonzero_mask) else 0.0

                status = "OK" if max_rel < 0.05 and max_abs < gamma * 0.05 else \
                         ("CHECK" if max_rel < 0.20 else "FAIL")
                print(f"  {topo_name:10s} B={B}  max rel err = {max_rel:.4f}  abs err = {max_abs:.4g}  [{status}]")
            print()


if __name__ == "__main__":
    main()
