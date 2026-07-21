"""
Pauli Weight Conjugation Proof - March 14, 2026
=================================================
Proves the mirror symmetry palindrome analytically.

The conjugation operator Π acts per site on Pauli indices:
  I → X (+1),  X → I (+1),  Y → iZ (+i),  Z → iY (+i)

Satisfies: Π · L · Π⁻¹ = -L - 2(Σγᵢ)·I
Therefore: decay rates d and 2Σγᵢ - d are always paired.

Scope: bond terms only. Adding a longitudinal on-site field h·ΣZ_i breaks the
identity (the field does not anti-commute with Π); the sweep below stays inside
the bond-only family the theorem covers.

Tested: N=2..5 (binary_tree from N=4) × star/chain/ring/complete/binary_tree;
        the δ and non-uniform-γ sweeps run at N=3,4; the dephasing-axis table
        (Z/Y/X/mixed/depolarizing) and the scope boundaries at N=3,4,5.
This script regenerates the N=3,4,5 rows and every scope-boundary number in the
"Numerical verification" section of docs/proofs/MIRROR_SYMMETRY_PROOF.md, and
produces simulations/results/conjugation_proof.txt in full. The N=6,7,8 rows of
that document come from the C# engine (RCPsiSquared.Compute), not from here.

Note on topology labels: several of the small graphs below coincide. At N=3,
star = path 1-0-2 is the chain relabelled, and ring and complete have the same
edge set; at N=4 the binary tree {(0,1),(0,2),(1,3)} is the path 2-0-1-3. Their
spectra are bit-identical, so the run counts below count configurations run,
not distinct physical systems.
"""

import numpy as np
from itertools import permutations, product as iproduct

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, Xm, Ym, Zm]
PAULI_NAMES = ['I', 'X', 'Y', 'Z']

# Π per-site map
PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}   # I↔X, Y↔Z
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}  # phases


def xy_weight(indices):
    return sum(1 for i in indices if i in (1, 2))


def build_hamiltonian_xxz(N, bonds, J=1.0, delta=1.0):
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    for (i, j) in bonds:
        for pidx, pauli in enumerate([Xm, Ym, Zm]):
            term = np.eye(1, dtype=complex)
            for k in range(N):
                term = np.kron(term, pauli if k in (i, j) else I2)
            H += J * (delta if pidx == 2 else 1.0) * term
    return H


def dissipator_rate(letter_idx, axis_idx):
    """Decay contribution of one site, in units of 2γ.

    A dephasing jump operator A (Hermitian, A² = I) acts on a Pauli letter P as
    A P A = ±P: the sign is −1 exactly when P anti-commutes with A, i.e. when P
    is neither the identity nor A itself. Each anti-commuting site contributes
    2γ to the decay rate, so the dissipator stays diagonal in the Pauli basis
    for every single-axis dephasing, not only for Z.
    """
    return 0 if letter_idx in (0, axis_idx) else 1


def build_liouvillian_pauli(N, H, gamma_per_site, axes=(3,), axis_scale=1.0):
    dim = 2**N
    num = 4**N
    all_idx = list(iproduct(range(4), repeat=N))
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for i in idx[1:]:
            m = np.kron(m, PAULIS[i])
        pmats.append(m)

    L_H = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)
        for a in range(num):
            L_H[a, b] = np.trace(pmats[a] @ comm) / dim

    # L_D diagonal: every site whose Pauli letter anti-commutes with a dephasing
    # axis contributes 2γᵢ per axis. For the default Z-dephasing (axes=(3,)) this
    # reduces to "each site carrying X or Y contributes 2γᵢ", the XY-weight rule.
    L_D = np.zeros((num, num), dtype=complex)
    for a, idx in enumerate(all_idx):
        rate = 0.0
        for site in range(N):
            # axes is either one global tuple of axis indices, or one tuple per
            # site (used by the mixed-axis rows, e.g. ZX = Z on even sites, X on odd)
            site_axes = axes[site] if len(axes) == N and isinstance(axes[0], tuple) else axes
            for axis in site_axes:
                rate += (2 * gamma_per_site[site] * axis_scale
                         * dissipator_rate(idx[site], axis))
        L_D[a, a] = -rate

    Pi = np.zeros((num, num), dtype=complex)
    for b, idx_b in enumerate(all_idx):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b:
            sign *= PI_SIGN[i]
        a = all_idx.index(mapped)
        Pi[a, b] = sign

    return L_H, L_D, L_H + L_D, Pi, all_idx


def get_bonds(N, topo):
    if isinstance(topo, (list, tuple)):
        return list(topo)
    if topo == 'star':
        return [(0, i) for i in range(1, N)]
    elif topo == 'chain':
        return [(i, i + 1) for i in range(N - 1)]
    elif topo == 'ring':
        return [(i, (i + 1) % N) for i in range(N)]
    elif topo == 'complete':
        return [(i, j) for i in range(N) for j in range(i + 1, N)]
    elif topo == 'binary_tree':
        bonds = []
        for i in range(N):
            for c in [2 * i + 1, 2 * i + 2]:
                if c < N:
                    bonds.append((i, c))
        return bonds
    return []


def amplitude_damping_spectrum(N, topo, gammas, delta=1.0, dephasing=None,
                               dephasing_axis=3, mixed_axes=None,
                               field=None, field_axis=3):
    """Liouvillian spectrum for Heisenberg + per-site amplitude damping (T1).

    Built in the vectorized (computational) basis because the T1 jump operator
    σ⁻ is not Hermitian, so the Pauli-basis shortcut used for dephasing (where
    A P A = ±P) does not apply.
    """
    dim = 2**N
    ident = np.eye(dim, dtype=complex)
    H = build_hamiltonian_xxz(N, get_bonds(N, topo), delta=delta)
    L = -1j * (np.kron(H, ident) - np.kron(ident, H.T))
    if field is not None:
        for site in range(N):
            term = np.eye(1, dtype=complex)
            for k in range(N):
                term = np.kron(term, AXIS_MATRIX[field_axis] if k == site else I2)
            H += field[site] * term
        L = -1j * (np.kron(H, ident) - np.kron(ident, H.T))
    lower = np.array([[0, 1], [0, 0]], dtype=complex)
    for site in range(N):
        A = np.eye(1, dtype=complex)
        for k in range(N):
            A = np.kron(A, lower if k == site else I2)
        AdA = A.conj().T @ A
        L += gammas[site] * (np.kron(A, A.conj())
                             - 0.5 * np.kron(AdA, ident)
                             - 0.5 * np.kron(ident, AdA.T))
    if dephasing is not None:
        for site in range(N):
            Aop = np.eye(1, dtype=complex)
            for k in range(N):
                ax = (mixed_axes[site] if mixed_axes is not None
                      else dephasing_axis)
                Aop = np.kron(Aop, AXIS_MATRIX[ax] if k == site else I2)
            L += dephasing[site] * (np.kron(Aop, Aop.T)
                                    - np.kron(ident, ident))
    return np.linalg.eigvals(L)


AXIS_MATRIX = {1: Xm, 2: Ym, 3: Zm}


def arbitrary_graph_spectrum(N, edges, axes_per_site, gammas, delta=1.0,
                             field=None, field_axis=3, J=None,
                             bond_terms=(1, 2, 3)):
    """Spectrum for an arbitrary edge set, per-site axis sets, and on-site field.

    Returns (eigenvalues, shift) where shift = 2·Σ over sites and axes of γ, the
    center the palindrome is tested about. Used for the scope-boundary rows that
    the topology helpers cannot express: disconnected graphs, several dephasing
    axes on one site, and transverse or longitudinal on-site fields.
    """
    dim = 2**N
    ident = np.eye(dim, dtype=complex)
    H = np.zeros((dim, dim), dtype=complex)
    for bidx, (i, j) in enumerate(edges):
        jj = J[bidx] if J is not None else 1.0
        # bond_terms selects which of XX, YY, ZZ the bond carries; a bond with
        # only ONE term relaxes clause 1 of the scope rule, so it is a knob
        for tidx in bond_terms:
            pauli = AXIS_MATRIX[tidx]
            term = np.eye(1, dtype=complex)
            for k in range(N):
                term = np.kron(term, pauli if k in (i, j) else I2)
            H += jj * (delta if tidx == 3 else 1.0) * term
    if field is not None:
        for site in range(N):
            term = np.eye(1, dtype=complex)
            for k in range(N):
                term = np.kron(term, AXIS_MATRIX[field_axis] if k == site else I2)
            H += field[site] * term
    L = -1j * (np.kron(H, ident) - np.kron(ident, H.T))
    shift = 0.0
    for site in range(N):
        for axis in axes_per_site[site]:
            A = np.eye(1, dtype=complex)
            for k in range(N):
                A = np.kron(A, AXIS_MATRIX[axis] if k == site else I2)
            L += gammas[site] * (np.kron(A, A.T) - np.kron(ident, ident))
            shift += gammas[site]
    return np.linalg.eigvals(L), 2 * shift


def field_identity_residual(N, edges, field, field_axis, gamma=0.05,
                            mix_angle=None):
    """max|Π·L·Π⁻¹ + L + 2Σγ·I| for a bond Hamiltonian plus an on-site field.

    Separate from the pairing counts on purpose: the identity and the spectral
    palindrome are different claims, and the transverse-Y field is exactly the
    case where one holds and the other does not.
    """
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    for (i, j) in edges:
        for pidx, pauli in enumerate([Xm, Ym, Zm]):
            term = np.eye(1, dtype=complex)
            for k in range(N):
                term = np.kron(term, pauli if k in (i, j) else I2)
            H += term
    # mix_angle rotates the field inside the XY plane, X at 0 and Y at π/2,
    # which is how the gauge nature of the X/Y asymmetry becomes visible
    axis_mat = (AXIS_MATRIX[field_axis] if mix_angle is None
                else np.cos(mix_angle) * Xm + np.sin(mix_angle) * Ym)
    for site in range(N):
        term = np.eye(1, dtype=complex)
        for k in range(N):
            term = np.kron(term, axis_mat if k == site else I2)
        H += field[site] * term
    _, _, L, Pi, _ = build_liouvillian_pauli(N, H, [gamma] * N)
    c = 2 * N * gamma
    return np.max(np.abs(Pi @ L @ np.linalg.inv(Pi) + L + c * np.eye(4**N)))


def mixed_direction_field(N, edges, field, angles, gamma=0.05, deph_axes=None):
    """Bond Hamiltonian plus a transverse field whose DIRECTION varies per site.

    angles[k] is site k's direction inside the XY plane, 0 = X, π/2 = Y. A
    single common angle is a global R_z rotation away from the X case and stays
    exact; several different angles are not, and that is where it breaks.
    Returns (eigenvalues, shift, identity residual).
    """
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    for (i, j) in edges:
        for pauli in (Xm, Ym, Zm):
            term = np.eye(1, dtype=complex)
            for k in range(N):
                term = np.kron(term, pauli if k in (i, j) else I2)
            H += term
    for site in range(N):
        axis = np.cos(angles[site]) * Xm + np.sin(angles[site]) * Ym
        term = np.eye(1, dtype=complex)
        for k in range(N):
            term = np.kron(term, axis if k == site else I2)
        H += field[site] * term
    axes = deph_axes if deph_axes is not None else (3,)
    _, _, L, Pi, _ = build_liouvillian_pauli(N, H, [gamma] * N, axes)
    c = 2 * N * gamma
    residual = np.max(np.abs(Pi @ L @ np.linalg.inv(Pi) + L + c * np.eye(4**N)))
    return np.linalg.eigvals(L), c, residual


def multiset_pairs(values, partner_of, tol=1e-7):
    """Count values that pair up as a MULTISET: each partner is consumed once.

    The naive counter ("does some partner exist within tol") double-counts and
    ignores multiplicity. Inside the theorem's scope the two agree at 100%, so
    the distinction never mattered for the passing rows; outside the scope the
    naive count inflates badly (depolarizing γ/3 at N=3: 29 naive vs 14 here),
    and those are exactly the numbers quoted as evidence of breakage.

    The matching is greedy nearest-first, so in general it is a LOWER bound on
    the maximum matching rather than the maximum itself. For the fully paired
    rows that makes no difference, and the 14 above was checked against an
    optimal assignment. Any new partial count taken from here should be
    re-checked the same way before being quoted.
    """
    pool = list(values)
    used = [False] * len(pool)
    matched = 0
    for v in values:
        target = partner_of(v)
        best, best_d = -1, tol
        for k, w in enumerate(pool):
            if used[k]:
                continue
            d = abs(w - target)
            if d < best_d:
                best_d, best = d, k
        if best >= 0:
            used[best] = True
            matched += 1
    return matched


def pauli_to_vec_basis(N):
    """Unitary taking the Pauli-string basis to the vectorized (computational)
    basis. Needed because a MAX ENTRY is basis-dependent: the same residual
    operator reads (2/3)Σγ in the Pauli basis and (2/9)Σγ here."""
    dim = 2**N
    cols = []
    for idx in iproduct(range(4), repeat=N):
        mat = PAULIS[idx[0]]
        for i in idx[1:]:
            mat = np.kron(mat, PAULIS[i])
        cols.append(mat.reshape(-1) / np.sqrt(dim))
    return np.column_stack(cols)


def permutation_sweep(N=3, topo='chain'):
    """All 24 label permutations × 256 phase assignments in {±1, ±i}.

    Returns [(label, working_phase_count)] for the permutations that admit any
    working phase. The survivors are exactly those swapping {I,Z} with {X,Y},
    which is Step 1 (w → N−w) stated as a permutation.
    """
    H = build_hamiltonian_xxz(N, get_bonds(N, topo))
    _, _, L, _, all_idx = build_liouvillian_pauli(N, H, [0.05] * N)
    c = 2 * N * 0.05
    num = 4**N
    pos = {t: i for i, t in enumerate(all_idx)}
    out = []
    for perm_t in permutations(range(4)):
        perm = {i: perm_t[i] for i in range(4)}
        good = 0
        for ph in iproduct([1, -1, 1j, -1j], repeat=4):
            Pi = np.zeros((num, num), dtype=complex)
            for b, idx_b in enumerate(all_idx):
                mapped = tuple(perm[i] for i in idx_b)
                sign = 1
                for i in idx_b:
                    sign *= ph[i]
                Pi[pos[mapped], b] = sign
            if np.max(np.abs(Pi @ L @ np.linalg.inv(Pi) + L
                             + c * np.eye(num))) < 1e-9:
                good += 1
        if good:
            label = " ".join(f"{PAULI_NAMES[i]}→{PAULI_NAMES[perm[i]]}"
                             for i in range(4))
            out.append((label, good))
    return out


def verify(N, topo, delta=1.0, gammas=None, axes=(3,), axis_scale=1.0,
           J=1.0, field_h=0.0):
    bonds = get_bonds(N, topo)
    if not bonds:
        return None
    H = build_hamiltonian_xxz(N, bonds, J=J, delta=delta)
    if field_h:
        # longitudinal on-site field h·ΣZ_i: outside the theorem's scope,
        # included so the scope-boundary section can measure the breakage
        for site in range(N):
            term = np.eye(1, dtype=complex)
            for k in range(N):
                term = np.kron(term, Zm if k == site else I2)
            H = H + field_h * term
    if gammas is None:
        gammas = [0.05] * N
    L_H, L_D, L, Pi, _ = build_liouvillian_pauli(N, H, gammas, axes, axis_scale)
    num = 4**N
    Pi_inv = np.linalg.inv(Pi)
    c = 2 * sum(gammas)

    err_H = np.max(np.abs(Pi @ L_H @ Pi_inv + L_H))
    err_D = np.max(np.abs(Pi @ L_D @ Pi_inv + L_D + c * np.eye(num)))
    residual_M = Pi @ L @ Pi_inv + L + c * np.eye(num)
    err_L = np.max(np.abs(residual_M))
    # the same operator, read in the computational basis: a max entry is not a
    # norm, so this number differs from err_L by a basis change alone
    V = pauli_to_vec_basis(N)
    err_L_vec = np.max(np.abs(V @ residual_M @ V.conj().T))

    eigs = np.linalg.eigvals(L)
    rates = -eigs.real
    center = c / 2
    paired = multiset_pairs(list(rates), lambda d: 2 * center - d)
    # The theorem is stated on the complex eigenvalue λ ↦ −λ − 2Σγ. Rate-only
    # pairing is the weaker reading: a U(1)-conserving on-site field leaves the
    # rates paired while destroying the complex pairing, so the two are tracked
    # separately.
    paired_complex = multiset_pairs(list(eigs), lambda lam: -lam - c)

    return {
        'err_H': err_H, 'err_D': err_D, 'err_L': err_L,
        'err_L_vec': err_L_vec,
        'shift': c, 'paired': paired, 'total': len(rates),
        'palindrome': f"{paired}/{len(rates)}",
        'palindrome_complex': f"{paired_complex}/{len(eigs)}",
        'palindromic': paired == len(rates),
        'min_re': float(eigs.real.min()),
        'oscillatory': int((np.abs(eigs.imag) > 0.05).sum()),
        'distinct_rates': len(set(np.round(eigs.real, 10))),
        'ok': err_L < 1e-10
    }


def topologies(N):
    if N == 2:
        return ['chain']          # every graph on 2 sites is the single bond
    topos = ['star', 'chain', 'ring', 'complete']
    if N >= 4:
        topos.append('binary_tree')
    return topos


def term_residual(N, topo, coeffs):
    """Max entry of Π·L_H·Π⁻¹ + L_H for a bond built from selected terms.

    coeffs multiplies (XX, YY, ZZ). Used to show that each term anti-commutes
    with Π on its own, which is what carries the identity across the whole XXZ
    family instead of a cancellation between terms at one particular δ.
    """
    dim = 2**N
    H = np.zeros((dim, dim), dtype=complex)
    for (i, j) in get_bonds(N, topo):
        for pidx, pauli in enumerate([Xm, Ym, Zm]):
            if coeffs[pidx] == 0:
                continue
            term = np.eye(1, dtype=complex)
            for k in range(N):
                term = np.kron(term, pauli if k in (i, j) else I2)
            H += coeffs[pidx] * term
    L_H, _, _, Pi, _ = build_liouvillian_pauli(N, H, [0.05] * N)
    return np.max(np.abs(Pi @ L_H @ np.linalg.inv(Pi) + L_H))


if __name__ == '__main__':
    print("Conjugation Proof: Π·L·Π⁻¹ = -L - 2Σγ·I")
    print("=" * 80)
    print("Scope: bond-only Hamiltonians (no longitudinal on-site field).")
    grand_total = 0
    grand_passed = 0

    # Section 1: topologies × N, Heisenberg δ=1, uniform γ=0.05
    print("\n## Topologies and N (Heisenberg δ=1, uniform γ=0.05, J=1)\n")
    print("Pairing tolerance: a partner counts if |d − (2Σγ − d')| < 1e-7;")
    print("the identity passes if max|Π·L·Π⁻¹ + L + 2Σγ·I| < 1e-10.\n")
    print(f"{'N':>3} {'Topo':>12} {'err_H':>10} {'err_L':>10} "
          f"{'Palindrome (λ)':>16} {'(rates)':>12}")
    print("-" * 70)
    n_sec = n_ok = 0
    for N in [2, 3, 4, 5]:
        # N=5 is where this construction stops being practical: projecting L_H
        # onto the Pauli basis is a double loop over 4^N strings, so the cost
        # grows as 16^N. N=6 is 256x the N=5 work and runs for hours here.
        for topo in topologies(N):
            r = verify(N, topo)
            n_sec += 1
            n_ok += r['ok'] and r['palindromic']
            print(f"{N:>3} {topo:>12} {r['err_H']:>10.1e} {r['err_L']:>10.1e} "
                  f"{r['palindrome_complex']:>16} {r['palindrome']:>12}")
    print(f"\n{n_ok}/{n_sec} passed.")
    grand_total += n_sec; grand_passed += n_ok

    # Section 1b: each bond term anti-commutes with Π on its own
    print("\n## Per-term anti-commutation (N=3 chain, max|Π·L_H·Π⁻¹ + L_H|)\n")
    print("The 16-row table in the proof document is δ=1. The XXZ family follows")
    print("because no cancellation BETWEEN terms is doing the work:\n")
    print(f"{'Bond term':>34} {'residual':>12}")
    print("-" * 48)
    for lbl, co in [("XX only", (1, 0, 0)), ("YY only", (0, 1, 0)),
                    ("ZZ only (the Ising bond)", (0, 0, 1)),
                    ("XX + YY (δ=0, the XY bond)", (1, 1, 0))]:
        print(f"{lbl:>34} {term_residual(3, 'chain', co):>12.3e}")
    for d in [-0.5, 0.3, 1.5, 2.0]:
        print(f"{'XX + YY + δZZ, δ=' + str(d):>34} "
              f"{term_residual(3, 'chain', (1, 1, d)):>12.3e}")

    # Section 2: XXZ δ sweep, all topologies at N=3,4
    DELTAS = [-0.5, 0.0, 0.3, 0.5, 1.0, 1.5, 2.0]
    print("\n## XXZ coupling (H = XX + YY + δ·ZZ, Z-dephasing γ=0.05)\n")
    print(f"δ = {', '.join(str(d) for d in DELTAS)}")
    d_sec = d_ok = 0
    for N in [3, 4]:
        n_here = ok_here = 0
        for topo in topologies(N):
            for delta in DELTAS:
                r = verify(N, topo, delta=delta)
                n_here += 1
                ok_here += r['ok'] and r['palindromic']
        print(f"N={N} × {', '.join(topologies(N))}: {ok_here}/{n_here} passed")
        d_sec += n_here; d_ok += ok_here
    print(f"Total: {d_ok}/{d_sec} passed. δ is irrelevant to the proof.")
    grand_total += d_sec; grand_passed += d_ok

    # Section 3: non-uniform γ, all topologies at N=3,4
    PROFILES = {
        3: [[0.03, 0.05, 0.07], [0.01, 0.02, 0.03], [0.10, 0.01, 0.05]],
        4: [[0.03, 0.05, 0.07, 0.04], [0.01, 0.02, 0.03, 0.04],
            [0.10, 0.01, 0.05, 0.02]],
    }
    print("\n## Non-uniform γ (Heisenberg δ=1, all topologies)\n")
    print(f"{'N':>3} {'Topo':>12} {'γ_values':>28} {'err_L':>8} {'center':>9} {'Palindrome':>14}")
    print("-" * 78)
    g_sec = g_ok = 0
    for N in [3, 4]:
        for topo in topologies(N):
            for gammas in PROFILES[N]:
                r = verify(N, topo, gammas=gammas)
                g_sec += 1
                g_ok += r['ok'] and r['palindromic']
                print(f"{N:>3} {topo:>12} {str(gammas):>28} {'✓' if r['ok'] else '✗':>8} "
                      f"{r['shift'] / 2:>9.4f} {r['palindrome']:>14}")
    print(f"\n{g_ok}/{g_sec} passed. Center = Σγᵢ, shift = 2Σγᵢ.")
    grand_total += g_sec; grand_passed += g_ok

    # Section 4: dephasing axis dependence (N=3 chain), Z-dephasing Π throughout
    print("\n## Dephasing axis dependence (N=3 chain, γ=0.05, the Z-dephasing Π)\n")
    print("The residual is the max entry of Π·L·Π⁻¹ + L + 2Σγ·I IN THE PAULI-STRING")
    print("BASIS (a max entry is basis-dependent, so the basis is part of the number).")
    print("'complex' pairs λ ↦ −λ−2Σγ, the theorem's statement; 'rates' pairs d ↦ 2Σγ−d.")
    print(f"\n{'Axis':>26} {'Π·L_H':>7} {'Π·L_D':>7} {'residual':>10} {'complex':>10} {'rates':>10}")
    print("-" * 76)
    AXIS_ROWS = [
        ("Z", (3,), 1.0),
        ("Y", (2,), 1.0),
        ("X", (1,), 1.0),
        ("mixed ZX", ((3,), (1,), (3,)), 1.0),
        ("depolarizing (γ per axis)", (1, 2, 3), 1.0),
        ("depolarizing (γ/3 per axis)", (1, 2, 3), 1.0 / 3.0),
    ]
    for label, axes, scale in AXIS_ROWS:
        r = verify(3, 'chain', axes=axes, axis_scale=scale)
        print(f"{label:>26} {'✓' if r['err_H'] < 1e-10 else '✗':>7} "
              f"{'✓' if r['err_D'] < 1e-10 else '✗':>7} "
              f"{r['err_L']:>10.4f} {r['palindrome_complex']:>10} "
              f"{r['palindrome']:>10}")
    sg = 3 * 0.05
    print(f"\nΣγ = {sg:.2f}. The two depolarizing rows differ only in convention:")
    print(f"  γ per axis   → residual = 6·Σγ = {6 * sg:.2f}")
    print(f"  γ/3 per axis → residual = (2/3)·Σγ = {2 / 3 * sg:.4f}  (the registered F1 value)")
    dep3 = verify(3, 'chain', axes=(1, 2, 3), axis_scale=1.0 / 3.0)
    print(f"\nThe residual is a MAX ENTRY, which is basis-dependent, not a norm.")
    print(f"Same operator, γ/3 depolarizing row, two bases:")
    print(f"  Pauli-string basis: {dep3['err_L']:.4f} = (2/3)·Σγ = {2 / 3 * sg:.4f}")
    print(f"  computational basis: {dep3['err_L_vec']:.4f} = (2/9)·Σγ = {2 / 9 * sg:.4f}")
    print("\nThe residual is |−8Nγs + 2Nγ| for scale s, so the ratio 9 between the two")
    print("rows is NOT 'three times the rate': the +2Σγ offset is not tripled with it.")
    print("Both break the palindrome outright: the complex pairing is 0/64 in BOTH")
    print("conventions. The γ/3 row's surviving rate-only count is the weaker reading,")
    print("which tolerates breakage the theorem's statement does not.")

    # Section 5: the scope boundaries claimed in the proof document
    print("\n## Scope boundaries (N=3 chain unless noted, γ=0.05)\n")
    print("Each row is a claim the proof document makes about what lies OUTSIDE")
    print("the theorem. 'complex' pairs λ ↦ −λ−2Σγ; 'rates' pairs d ↦ 2Σγ−d.\n")
    print(f"{'Case':>34} {'residual':>10} {'complex':>10} {'rates':>10}")
    print("-" * 68)
    base = verify(3, 'chain')
    print(f"{'bond-only H (in scope)':>34} {base['err_L']:>10.4f} "
          f"{base['palindrome_complex']:>10} {base['palindrome']:>10}")
    # One field vector for every axis, so the columns are comparable: the
    # residual turns out to depend only on max|h| and NOT on the axis, so it
    # cannot grade the cases. Only the pairing columns can.
    H_UNI = [0.3, 0.3, 0.3]
    H_NON = [0.3, -0.2, 0.5]
    rows = [
        ("+ longitudinal Z field, uniform h=0.3", H_UNI, 3),
        ("+ longitudinal Z field, non-uniform", H_NON, 3),
        ("+ transverse X field, non-uniform", H_NON, 1),
        ("+ transverse Y field, non-uniform", H_NON, 2),
    ]
    edges3 = [(0, 1), (1, 2)]
    for lbl, hs, axis in rows:
        ev, shift = arbitrary_graph_spectrum(3, edges3, [(3,)] * 3, [0.05] * 3,
                                             field=hs, field_axis=axis)
        cx = multiset_pairs(list(ev), lambda lam: -lam - shift)
        rt = multiset_pairs(list(-ev.real), lambda d: shift - d)
        res = field_identity_residual(3, edges3, hs, axis)
        print(f"{lbl:>40} {res:>10.4f} {f'{cx}/64':>10} {f'{rt}/64':>10}")
    print("\nThe Y row's nonzero residual is a GAUGE artifact, not a property of the")
    print("Y direction. Rotating the field in the XY plane by an angle φ away from X,")
    print("the residual under this FIXED Π traces 2·max|h|·|sin φ|:\n")
    print(f"{'φ (deg)':>9} {'residual':>10} {'2·max|h|·|sin φ|':>18}")
    print("-" * 40)
    for deg in [0, 15, 30, 45, 60, 90]:
        r = np.deg2rad(deg)
        res = field_identity_residual(3, edges3, H_NON, 1, mix_angle=r)
        print(f"{deg:>9} {res:>10.6f} "
              f"{2 * max(abs(x) for x in H_NON) * abs(np.sin(r)):>18.6f}")
    print("\nA per-site R_z(π/2) maps the X-field Hamiltonian to the Y-field one and")
    print("is a symmetry of Z-dephasing, so the two Liouvillians are unitarily")
    print("equivalent and their spectra are identical. Conjugating Π by the same")
    print("rotation satisfies the identity exactly for the Y field. So the residual")
    print("under a fixed Π grades the ANGLE to that mirror's preferred axis.")
    print("\nBut that argument needs ONE common direction: a global rotation cannot")
    print("align several at once. Let the direction vary per site and it breaks.\n")
    print(f"{'per-site field directions':>34} {'residual':>10} {'pairing':>10}")
    print("-" * 58)
    for angles, lbl in [
        ([0.0, 0.0, 0.0], "all along X"),
        ([np.pi / 2] * 3, "all along Y"),
        ([0.7, 0.7, 0.7], "all at 40 deg"),
        ([0.0, np.pi / 2, 0.0], "X, Y, X"),
        ([0.0, 0.7, 1.9], "0, 40, 109 deg"),
    ]:
        ev, shift, res = mixed_direction_field(3, edges3, H_NON, angles)
        print(f"{lbl:>34} {res:>10.4f} "
              f"{f'{multiset_pairs(list(ev), lambda lam: -lam - shift)}/64':>10}")
    print("\nUniform direction, any angle: exact. Direction varying between sites:")
    print("gone. Magnitudes may vary freely in both cases. The general residual law")
    print("is 2·maxᵢ|hᵢ·sin φᵢ|, which collapses to 2·max|h|·|sin φ| for a common")
    print("angle. This rhymes with the dephasing-axis rule below, where mixing axes")
    print("across sites is also the thing that eventually breaks the palindrome.")
    print("\nONE rule covers both field directions: the field survives exactly when")
    print("its axis is ORTHOGONAL to every dephasing axis in play. Avoiding the")
    print("used axes is not enough, as the 45° rows below show. 'Transverse' and")
    print("'longitudinal' are not two facts, they are this one fact seen from a")
    print("Z-dephasing default. Every axis-set and field-direction combination:\n")
    print(f"{'dephasing axes present':>24} {'field X':>9} {'field Y':>9} {'field Z':>9}")
    print("-" * 54)
    for ax in [(3, 3, 3), (3, 1, 3), (3, 2, 3), (1, 2, 1), (1, 1, 2), (2, 3, 2)]:
        row = []
        for fld in (1, 2, 3):
            ev, shift = arbitrary_graph_spectrum(
                3, edges3, [(a,) for a in ax], [0.05] * 3,
                field=H_NON, field_axis=fld)
            row.append(multiset_pairs(list(ev), lambda lam: -lam - shift))
        present = "+".join(sorted({ {1: 'X', 2: 'Y', 3: 'Z'}[a] for a in ax }))
        print(f"{present:>24} {row[0]:>9} {row[1]:>9} {row[2]:>9}")
    print("\nEvery 64 sits orthogonal to the noise's axes, every 0 does not. Note")
    print("that AVOIDING the used axes is not sufficient, only orthogonality is:\n")
    print(f"{'dephasing':>14} {'field direction':>22} {'pairing':>10}")
    print("-" * 50)
    for ax, ang, lbl in [
        ((3, 3, 3), 0.0, "along X"),
        ((3, 3, 3), np.pi / 4, "45° in the XY plane"),
        ((3, 1, 3), np.pi / 2, "along Y"),
        ((3, 1, 3), np.pi / 4, "45° in the XY plane"),
    ]:
        ev, shift, _ = mixed_direction_field(3, edges3, H_NON, [ang] * 3,
                                             deph_axes=[(a,) for a in ax])
        present = "+".join(sorted({ {1: 'X', 2: 'Y', 3: 'Z'}[a] for a in ax }))
        print(f"{present:>14} {lbl:>22} "
              f"{f'{multiset_pairs(list(ev), lambda lam: -lam - shift)}/64':>10}")
    print("\nThe 45° field lies along neither Z nor X, yet dies once both are")
    print("dephasing: it is not orthogonal to X. With Z alone the whole plane is")
    print("orthogonal to the noise, so the same 45° field is exact.")

    print("\n### Per-component freedom, and T1's field intolerance\n")
    ev, shift, _ = mixed_direction_field(
        4, [(0, 1), (2, 3)], [0.3, -0.2, 0.5, 0.4],
        [0.0, 0.0, np.pi / 2, np.pi / 2], deph_axes=[(3,)] * 4)
    print(f"  two disjoint bonds, X field on one, Y on the other: "
          f"{multiset_pairs(list(ev), lambda lam: -lam - shift)}/256")
    ev, shift, _ = mixed_direction_field(
        4, [(0, 1), (1, 2), (2, 3)], [0.3, -0.2, 0.5, 0.4],
        [0.0, 0.0, np.pi / 2, np.pi / 2], deph_axes=[(3,)] * 4)
    print(f"  the same pattern on a connected chain:              "
          f"{multiset_pairs(list(ev), lambda lam: -lam - shift)}/256")
    ev = amplitude_damping_spectrum(3, 'chain', [0.05] * 3,
                                    dephasing=[0.05, 0.05, 0.05],
                                    dephasing_axis=1, mixed_axes=(1, 2, 1))
    print(f"  T1 + two-axis XYX dephasing:                        "
          f"{multiset_pairs(list(ev), lambda lam: -lam - (3 * 0.05 + 2 * 3 * 0.05))}/64")
    for fa, lbl in [(1, "transverse X"), (3, "longitudinal Z")]:
        ev = amplitude_damping_spectrum(3, 'chain', [0.05] * 3,
                                        field=H_NON, field_axis=fa)
        s = 3 * 0.05
        cx = multiset_pairs(list(ev), lambda lam: -lam - s)
        rt = multiset_pairs(list(-ev.real), lambda d: s - d)
        print(f"  T1 alone plus a {lbl:14s} field:             "
              f"{f'{cx}/64':>7} complex, {f'{rt}/64':>7} rates")
    print("\nT1 tolerates no on-site field at all, unlike the dephasing-only case.")

    print("\nBoth T1 restrictions are per component. Two disjoint Heisenberg")
    print("bonds (0-1),(2-3) at N=4, T1 on sites 0-1 at 0.05; pairing about the")
    print("combined shift 2Σγ_deph + Σγ_T1:\n")
    disjoint = [(0, 1), (2, 3)]
    for lbl, kwargs, s_mix in [
        ("co-axial Z-dephasing on the OTHER bond",
         dict(dephasing=[0, 0, 0.05, 0.05], dephasing_axis=3), 0.3),
        ("X field on the OTHER bond",
         dict(field=[0, 0, 0.3, -0.2], field_axis=1), 0.1),
        ("co-axial Z-dephasing on T1's OWN bond",
         dict(dephasing=[0.05, 0.05, 0, 0], dephasing_axis=3), 0.3),
        ("X field on T1's OWN bond",
         dict(field=[0.3, -0.2, 0, 0], field_axis=1), 0.1),
    ]:
        ev = amplitude_damping_spectrum(4, disjoint, [0.05, 0.05, 0, 0],
                                        **kwargs)
        n = multiset_pairs(list(ev), lambda lam: -lam - s_mix)
        print(f"  {lbl:42s} {n:>3}/256")
    print("\nOn the other component neither channel matters; on T1's own")
    print("component both bite.")

    print("\n### Clause 1 needs a two-term bond: single-term bonds take three axes\n")
    print("All eighteen three-distinct-axis assignments at N=4 chain, by bond type.")
    print("Two bond terms are what force the two-axis ceiling; one term does not.\n")
    print(f"{'bond terms':>26} {'3-axis: hold/fail':>20}")
    print("-" * 48)
    for terms, lbl in [((1, 2, 3), "XX+YY+ZZ (Heisenberg)"), ((1, 3), "XX+ZZ"),
                       ((1, 2), "XX+YY"), ((3,), "ZZ only (the Ising bond)"),
                       ((2,), "YY only"), ((1,), "XX only")]:
        hold = fail = 0
        for assign in permutations([1, 2, 3]):
            for extra in (1, 2, 3):
                ax = list(assign) + [extra]
                if len(set(ax)) < 3:
                    continue
                ev, shift = arbitrary_graph_spectrum(
                    4, [(0, 1), (1, 2), (2, 3)], [(a,) for a in ax],
                    [0.05] * 4, bond_terms=terms)
                ok = multiset_pairs(list(ev), lambda lam: -lam - shift) == 256
                hold += ok
                fail += not ok
        print(f"{lbl:>26} {f'{hold}/{fail}':>20}")
    ax = [1, 2, 3]
    ev, shift = arbitrary_graph_spectrum(
        3, edges3, [(1, 2, 3), (3,), (3,)], [0.05] * 3, bond_terms=(3,))
    print(f"\nThe escape has two limits. Three axes stacked on ONE site is depolarizing")
    print(f"there, and breaks the Ising bond too "
          f"({multiset_pairs(list(ev), lambda lam: -lam - shift)}/64). And a")
    print("three-axis component must be FIELD-FREE, which is clause 2 read carefully:")
    print("no direction is orthogonal to X, Y and Z at once, so only the empty field")
    print("passes it. The two clauses look like they conflict here; they do not:\n")
    print(f"{'Ising bond, axes X,Y,Z':>34} {'pairing':>10}")
    print("-" * 46)
    for fld, lbl in [(None, "no field"), (1, "+ X field"), (2, "+ Y field"),
                     (3, "+ Z field")]:
        ev2, sh2 = arbitrary_graph_spectrum(
            3, edges3, [(1,), (2,), (3,)], [0.05] * 3, bond_terms=(3,),
            field=(H_NON if fld else None), field_axis=(fld or 3))
        print(f"{lbl:>34} "
              f"{f'{multiset_pairs(list(ev2), lambda lam: -lam - sh2)}/64':>10}")

    print("\n### The three qualifiers in the rule, each with its counterexample\n")
    print("Each row is a case the rule would forbid if the qualifier were dropped,")
    print("and which in fact holds. This is what the words COMPONENT (nonzero")
    print("coupling), PRESENT (nonzero rate) and 'carries dephasing at all' buy.\n")
    print(f"{'case':>50} {'pairing':>10}")
    print("-" * 62)
    # component = nonzero coupling: a J=0 edge splits the path
    ev, shift = arbitrary_graph_spectrum(
        4, [(0, 1), (1, 2), (2, 3)], [(1,), (2,), (3,), (3,)], [0.05] * 4,
        J=[1.0, 0.0, 1.0])
    print(f"{'N=4 path, J=[1,0,1], axes X,Y,Z,Z (3 axes)':>50} "
          f"{f'{multiset_pairs(list(ev), lambda lam: -lam - shift)}/256':>10}")
    ev, shift = arbitrary_graph_spectrum(
        4, [(0, 1), (1, 2), (2, 3)], [(1,), (2,), (3,), (3,)], [0.05] * 4)
    print(f"{'   control, J=[1,1,1] (one component)':>50} "
          f"{f'{multiset_pairs(list(ev), lambda lam: -lam - shift)}/256':>10}")
    # present = nonzero rate
    ev, shift = arbitrary_graph_spectrum(
        3, edges3, [(3,), (3,), (1,)], [0.05, 0.05, 0.0],
        field=H_NON, field_axis=1)
    print(f"{'N=3 chain, axes Z,Z,X with γ_X=0, X field':>50} "
          f"{f'{multiset_pairs(list(ev), lambda lam: -lam - shift)}/64':>10}")
    ev, shift = arbitrary_graph_spectrum(
        3, edges3, [(3,), (3,), (1,)], [0.05] * 3, field=H_NON, field_axis=1)
    print(f"{'   control, γ_X=0.05':>50} "
          f"{f'{multiset_pairs(list(ev), lambda lam: -lam - shift)}/64':>10}")
    # a component with no dephasing constrains nothing
    ev, shift = arbitrary_graph_spectrum(
        2, [(0, 1)], [(), ()], [0.0, 0.0], field=[0.3, 0.4], field_axis=1)
    print(f"{'N=2 bond, NO dephasing, field present':>50} "
          f"{f'{multiset_pairs(list(ev), lambda lam: -lam - shift)}/16':>10}")

    print("\n### The boundary is sharp, the failure is not\n")
    print("Chain dephased along Z with an X field (legal), then a competing")
    print("X-dephasing switched on at rate γ_X on one site (forbidden). The")
    print("pairing degrades with γ_X rather than collapsing at the first")
    print("nonzero value, so a small violating term costs proportionally little:\n")
    print(f"{'γ_X':>10} {'pairs @1e-7':>13} {'pairs @1e-10':>14}")
    print("-" * 40)
    for gx in [0.0, 1e-12, 1e-9, 1e-6, 1e-3, 0.05]:
        axes = [(3,), (3,), (1,)]
        gam = [0.05, 0.05, gx]
        ev, shift = arbitrary_graph_spectrum(3, edges3, axes, gam,
                                             field=H_NON, field_axis=1)
        p7 = multiset_pairs(list(ev), lambda lam: -lam - shift, tol=1e-7)
        p10 = multiset_pairs(list(ev), lambda lam: -lam - shift, tol=1e-10)
        print(f"{gx:>10.0e} {f'{p7}/64':>13} {f'{p10}/64':>14}")
    print("\nA measurement can only see a violation larger than its own tolerance.")
    print("The rule says where the palindrome is EXACT, not that everything")
    print("outside is equally far away.")
    print("\nAnd the surviving rate pairing under the longitudinal field comes from")
    print("UNIFORMITY, not from U(1) conservation: a uniform h·ΣZ commutes with H")
    print("and with the dephasing, so it enters as a pure imaginary shift. A")
    print("non-uniform longitudinal field is just as U(1)-conserving and wrecks the")
    print("rate pairing too, as the second row shows.")

    print("\n### Amplitude damping (T1) alone: palindromic, at HALF the shift\n")
    print("T1 is not a dephasing channel and no site-local Π of the shape used in")
    print("this proof reproduces it. But the spectrum is palindromic anyway, about")
    print("−Σγ/2 instead of −Σγ: the one-site T1 dissipator has Pauli-basis")
    print("eigenvalues {0, −γ/2, −γ/2, −γ}, already mirrored about −γ/2, and N sites")
    print("are a tensor sum. So 'Π fails' and 'the palindrome fails' are different")
    print("claims, and only the first holds for pure T1.\n")
    print(f"{'N':>3} {'Topo':>10} {'δ':>6} {'γ':>22} {'pairs about −Σγ/2':>20}")
    print("-" * 66)
    t1_rows = []
    for N in [2, 3, 4, 5]:
        for topo in topologies(N):
            if topo == 'binary_tree':
                continue
            t1_rows.append((N, topo, 1.0, [0.05] * N))
    for d in [-0.5, 0.0, 2.0]:                    # anisotropy sweep at N=3
        t1_rows.append((3, 'chain', d, [0.05] * 3))
    t1_rows.append((3, 'chain', 1.0, [0.03, 0.07, 0.05]))   # non-uniform rates
    t1_rows.append((4, 'ring', 1.0, [0.02, 0.06, 0.04, 0.08]))
    t1_all = True
    for N, topo, d, gv in t1_rows:
        ev = amplitude_damping_spectrum(N, topo, gv, delta=d)
        s = sum(gv)
        pairs = multiset_pairs(list(ev), lambda lam: -lam - s)
        t1_all &= pairs == len(ev)
        print(f"{N:>3} {topo:>10} {d:>6.1f} {str(gv):>22} "
              f"{f'{pairs}/{len(ev)}':>20}")
    print(f"\nAll rows fully paired: {t1_all}")
    print("\nT1 WITH dephasing depends entirely on the axis, and the obvious guess is")
    print("wrong: dephasing TRANSVERSE to T1 composes with it exactly. Only co-axial")
    print("Z-dephasing breaks it (N=3 chain, both channels at 0.05, paired about the")
    print("combined shift 2Σγ_deph + Σγ_T1):\n")
    print(f"{'dephasing beside T1':>26} {'pairing':>10}")
    print("-" * 38)
    for ax, lbl in [(1, 'X (transverse)'), (2, 'Y (transverse)'),
                    (3, 'Z (co-axial)')]:
        ev = amplitude_damping_spectrum(3, 'chain', [0.05] * 3,
                                        dephasing=[0.05] * 3, dephasing_axis=ax)
        s = 3 * 0.05 + 2 * 3 * 0.05
        print(f"{lbl:>26} "
              f"{f'{multiset_pairs(list(ev), lambda lam: -lam - s)}/64':>10}")
    print("\nSo the realistic hardware case, T1 alongside Z-dephasing, is the one")
    print("combination that fails, and it fails because the two channels share an")
    print("axis, not because dephasing is present at all.\n")
    print("At the theorem's own center, ignoring T1's half-shift (N=3 chain):\n")
    ev_mix = amplitude_damping_spectrum(3, 'chain', [0.05] * 3, dephasing=[0.05] * 3)
    theorem_shift = 2 * 3 * 0.05
    best_c, best_n = None, -1
    for cand in np.linspace(0.0, 1.5, 1501):
        n = multiset_pairs(list(ev_mix), lambda lam: -lam - cand)
        if n > best_n:
            best_n, best_c = n, cand
    at_theorem = multiset_pairs(list(ev_mix), lambda lam: -lam - theorem_shift)
    print(f"  T1 alone:            64/64 about −Σγ/2")
    print(f"  T1 + dephasing:      {at_theorem}/64 at the theorem's own center,")
    print(f"                       {best_n}/64 at the most favourable center")
    print(f"                       found anywhere (the count is matcher-stable,")
    print(f"                       its location is not, so it is not quoted)")
    print("\nThe second number is the generous test. Applied to depolarizing it is")
    print("no kinder: judged the same way, γ-per-axis depolarizing also reaches a")
    print("comparable count, so neither channel degrades gracefully. Quote the two")
    print("yardsticks separately or not at all.")

    print("\n### Mixed dephasing axes: ≤ 2 axes per CONNECTED COMPONENT\n")
    print("Every per-site assignment of an axis in {X, Y, Z}, exhaustively. The")
    print("criterion is per connected component, not global: a disconnected graph")
    print("carrying three axes overall still holds, as the control rows show.\n")
    print(f"{'N':>3} {'Topo':>10} {'≤2 axes: hold/fail':>20} {'3 axes: hold/fail':>20}")
    print("-" * 56)
    for N, topo in [(3, 'chain'), (3, 'ring'), (3, 'complete'), (4, 'chain')]:
        h2 = f2 = h3 = f3 = 0
        for assign in iproduct([1, 2, 3], repeat=N):
            r = verify(N, topo, axes=tuple((a,) for a in assign))
            full = r['palindrome_complex'] == f"{4**N}/{4**N}"
            if len(set(assign)) <= 2:
                h2 += full
                f2 += not full
            else:
                h3 += full
                f3 += not full
        print(f"{N:>3} {topo:>10} {f'{h2}/{f2}':>20} {f'{h3}/{f3}':>20}")
    print("\nZero exceptions in either direction on these connected graphs. Control")
    print("rows, where the graph is NOT connected (pairing about 2Σγ):\n")
    for lbl, N, edges, per_site in [
        ("two disjoint bonds, axes X,Y,Z,Z (3 global, 2 per part)", 4,
         [(0, 1), (2, 3)], [(1,), (2,), (3,), (3,)]),
        ("one bond + isolated site, axes X,Y,Z (3 global)", 3,
         [(0, 1)], [(1,), (2,), (3,)]),
        ("connected chain, axes X,Y,Z,Z (3 in one component)", 4,
         [(0, 1), (1, 2), (2, 3)], [(1,), (2,), (3,), (3,)]),
    ]:
        ev, shift = arbitrary_graph_spectrum(N, edges, per_site, [0.05] * N)
        print(f"  {lbl:56s} "
              f"{multiset_pairs(list(ev), lambda lam: -lam - shift)}/{len(ev)}")
    print("\nSo single-axis dephasing, the scope this proof claims, is strictly")
    print("conservative: two axes per component also work. Three axes inside one")
    print("component is depolarizing's mechanism, and that is where it dies.")

    print("\n### Maximum decay rate is pinned by trace preservation\n")
    print("L is trace-preserving, so λ=0 is always in the spectrum, and Re λ ≤ 0")
    print("for any Lindblad generator. The palindrome forces the partner −2Σγ to")
    print("be an eigenvalue, fixing the fastest decay at exactly 2Σγ:\n")
    print(f"{'N':>3} {'Topo':>12} {'min Re λ':>24} {'−2Σγ':>10}")
    print("-" * 52)
    for N, topo in [(3, 'chain'), (4, 'chain'), (4, 'complete'), (5, 'ring')]:
        r = verify(N, topo)
        print(f"{N:>3} {topo:>12} {r['min_re']:>24.15f} {-r['shift']:>10.4f}")

    print("\n### Decay rates are NOT topology-independent (only the palindrome is)\n")
    print(f"{'N':>3} {'Topo':>12} {'distinct Re λ':>15} {'min Re λ':>12}")
    print("-" * 45)
    ref = None
    for topo in ['chain', 'star', 'ring', 'complete']:
        r = verify(4, topo)
        print(f"{4:>3} {topo:>12} {r['distinct_rates']:>15} {r['min_re']:>12.6f}")
        if topo == 'chain':
            ref = r
    print("\nThe centers and the pairing agree across all four; the spectra do not.")

    print("\n### The oscillatory count depends on J (the palindrome does not)\n")
    print(f"{'N':>3} {'Topo':>8} {'J':>8} {'|Im λ|>0.05':>13} {'Palindrome':>14}")
    print("-" * 50)
    for J in [1.0, 0.5, 0.25, 0.075]:
        r = verify(5, 'chain', J=J)
        print(f"{5:>3} {'chain':>8} {J:>8.3f} {r['oscillatory']:>13} "
              f"{r['palindrome']:>14}")

    print("\n### The dedicated X-dephasing Π (I↔Z, X↔Y) and its phase\n")
    print("The Z-dephasing Π fails on X-dephasing, but a dedicated Π works. Its")
    print("phase on the X↔Y swap only has to be imaginary; the sign is free:\n")
    print(f"{'phase on X↔Y':>16} {'residual':>12}")
    print("-" * 30)
    N_x = 3
    H_x = build_hamiltonian_xxz(N_x, get_bonds(N_x, 'chain'))
    _, _, L_x, _, idx_x = build_liouvillian_pauli(N_x, H_x, [0.05] * N_x, axes=(1,))
    c_x = 2 * N_x * 0.05
    pos_x = {t: i for i, t in enumerate(idx_x)}
    perm_x = {0: 3, 3: 0, 1: 2, 2: 1}          # I↔Z, X↔Y
    for phase, label in [(-1j, "−i"), (1j, "+i"), (1, "+1"), (-1, "−1")]:
        sgn = {0: 1, 3: 1, 1: phase, 2: phase}
        Pi_x = np.zeros((4**N_x, 4**N_x), dtype=complex)
        for b, idx_b in enumerate(idx_x):
            mapped = tuple(perm_x[i] for i in idx_b)
            s = 1
            for i in idx_b:
                s *= sgn[i]
            Pi_x[pos_x[mapped], b] = s
        r = np.max(np.abs(Pi_x @ L_x @ np.linalg.inv(Pi_x) + L_x
                          + c_x * np.eye(4**N_x)))
        print(f"{label:>16} {r:>12.4f}")

    print("\n### Which label permutations admit a working Π at all\n")
    print("All 24 permutations of {I,X,Y,Z} × all 256 phase assignments in {±1,±i}:\n")
    for label, count in permutation_sweep():
        print(f"  {label}    {count}/256 phases work")
    print("\nExactly these four, and they are exactly the four swapping the sets")
    print("{I,Z} and {X,Y}: two involutions and two 4-cycles. That swap IS the")
    print("XY-weight flip w → N−w of Step 1, so the condition is structural, not")
    print("a lucky choice of labels. No purely real phase assignment works.")

    print(f"\n{'=' * 80}")
    print(f"TOTAL (sections 1-3): {grand_passed}/{grand_total} configurations passed")
    print(f"{'=' * 80}")
