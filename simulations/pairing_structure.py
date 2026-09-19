"""Finite N=2/N=5 frequency-bin comparison and F1 multiset matching.

At gamma=0.05, oscillatory entries satisfy |Im(lambda)|>0.01. Absolute
frequencies are rounded to six decimals and compared at tolerance 1e-4.
OLD/NEW denotes tolerance-matched/unmatched frequency values across N.
This census does not transport eigenvectors or identify persistent modes.
The XY-weight histogram is a reading of the returned N=5 eigenbasis.
"""
from itertools import product as _product

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def site_op(op, site, N):
    ops = [I2] * N
    ops[site] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def heisenberg_H(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, N) @ site_op(P, i + 1, N)
    return H


def build_liouvillian(H, gamma, N):
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d * d, dtype=complex))
    return L


def f1_partner(value, sigma):
    """Holomorphic F1 map with centre -sigma."""
    return -value - 2.0 * sigma


def maximum_cardinality_f1_match(spectrum, sigma, tolerance):
    """Maximum number of one-use matches with absolute distance < tolerance.

    This is a bipartite multiset count, not a list of unordered mode pairs.
    The supplied tolerance is a diagnostic budget, not an eigensolver error bound.
    """
    values = np.asarray(spectrum, dtype=complex)
    rows = []
    columns = []
    for source_index, source in enumerate(values):
        hits = np.flatnonzero(
            np.abs(f1_partner(source, sigma) - values) < tolerance
        )
        rows.extend([source_index] * hits.size)
        columns.extend(hits.tolist())
    graph = csr_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, columns)),
        shape=(len(values), len(values)),
    )
    matching = maximum_bipartite_matching(graph, perm_type="column")
    return int(np.count_nonzero(matching != -1))


def osc_frequencies(evals, threshold=0.01):
    """Six-decimal absolute-frequency bins above the stated threshold."""
    return sorted(set(round(abs(ev.imag), 6) for ev in evals
                      if abs(ev.imag) > threshold))


def freq_is_old(f, old_set, tol=1e-4):
    """Whether a frequency value is within absolute tolerance of the reference."""
    return any(abs(f - fo) < tol for fo in old_set)


def main():
    # Step 1: N=2 spectrum
    print("=" * 60)
    print("STEP 1: N=2 SINGLE RESONATOR")
    print("=" * 60)

    N2 = 2
    gamma = 0.05
    Sg2 = N2 * gamma  # 0.10

    H2 = heisenberg_H(N2)
    L2 = build_liouvillian(H2, gamma, N2)
    evals2 = np.linalg.eigvals(L2)

    freqs2 = osc_frequencies(evals2)
    print(f"Eigenvalues: {len(evals2)}")
    print(f"Distinct frequencies (|Im|>0.01): {len(freqs2)}")
    print(f"Frequencies: {freqs2}")
    print(f"Sg = {Sg2}")

    n_paired2 = maximum_cardinality_f1_match(evals2, Sg2, 1e-4)
    print(f"F1 matched entries: {n_paired2}/{len(evals2)}")


    # Step 2: N=5 spectrum
    print("\n" + "=" * 60)
    print("STEP 2: N=5 COUPLED SYSTEM (MediatorBridge)")
    print("=" * 60)

    N5 = 5
    Sg5 = N5 * gamma  # 0.25

    H5 = heisenberg_H(N5)
    L5 = build_liouvillian(H5, gamma, N5)
    evals5 = np.linalg.eigvals(L5)

    freqs5 = osc_frequencies(evals5)
    print(f"Eigenvalues: {len(evals5)}")
    print(f"Distinct frequencies (|Im|>0.01): {len(freqs5)}")
    print(f"Sg = {Sg5}")

    n_paired5 = maximum_cardinality_f1_match(evals5, Sg5, 1e-4)
    print(f"F1 matched entries: {n_paired5}/{len(evals5)}")


    # Step 3: Classify each N=5 frequency as OLD or NEW
    print("\n" + "=" * 60)
    print("STEP 3: OLD vs NEW FREQUENCY CLASSIFICATION")
    print("=" * 60)

    freq_tol = 1e-4  # absolute frequency-value comparison budget

    old_freqs = []
    new_freqs = []
    for f in freqs5:
        is_old = any(abs(f - f2) < freq_tol for f2 in freqs2)
        if is_old:
            old_freqs.append(f)
        else:
            new_freqs.append(f)

    print(f"OLD frequency bins (within tolerance of N=2): {len(old_freqs)}")
    print(f"NEW frequency bins (outside tolerance of N=2): {len(new_freqs)}")


    # Step 4: Population counts, independently of any unordered pairing.
    print("\n" + "=" * 60)
    print("STEP 4: SPECTRAL POPULATIONS")
    print("=" * 60)
    osc_mask = np.abs(evals5.imag) > 0.01
    n_osc = int(osc_mask.sum())
    n_real_axis = len(evals5) - n_osc
    n_old_osc = sum(1 for ev in evals5 if abs(ev.imag) > 0.01
                    and freq_is_old(round(abs(ev.imag), 6), freqs2))
    print(f"Oscillatory population (|Im|>0.01): {n_osc}")
    print(f"Complement population (|Im|<=0.01): {n_real_axis}")
    print(f"Oscillatory entries at OLD frequency values: {n_old_osc}")
    print("A full F1 match counts entries; it does not identify unordered pairs.")
    print("Cross-N frequency mismatches do not establish mode destruction or creation.")


    # Step 5: XY-weight distribution of the oscillating modes (the balance histogram)
    # For each oscillating right eigenvector (|Im| > 0.01), project onto the
    # orthonormal Pauli basis vec(P)/sqrt(d) and accumulate the probability mass
    # per XY-weight (number of X/Y letters). The histogram is a mass SHARE over
    # modes, not an integer count.
    print("\n" + "=" * 60)
    print("STEP 5: XY-WEIGHT DISTRIBUTION (oscillating modes, N=5)")
    print("=" * 60)

    _evals5v, _evecs5 = np.linalg.eig(L5)
    _osc_idx = [k for k in range(len(_evals5v)) if abs(_evals5v[k].imag) > 0.01]
    print(f"Oscillating modes: {len(_osc_idx)}/{len(_evals5v)}")

    _d5 = 2 ** N5
    _letters = [I2, sx, sy, sz]
    _weights = np.zeros(N5 + 1)
    _P_cache = []
    for _combo in _product(range(4), repeat=N5):
        _P = _letters[_combo[0]]
        for _c in _combo[1:]:
            _P = np.kron(_P, _letters[_c])
        _w = sum(1 for _c in _combo if _c in (1, 2))  # X or Y letters
        _P_cache.append((_P.reshape(-1) / np.sqrt(_d5), _w))

    for _k in _osc_idx:
        _v = _evecs5[:, _k]
        _v = _v / np.linalg.norm(_v)
        for _pvec, _w in _P_cache:
            _weights[_w] += abs(np.vdot(_pvec, _v)) ** 2

    _weights_pct = 100 * _weights / _weights.sum()
    for _w in range(N5 + 1):
        print(f"  w={_w}: {_weights_pct[_w]:5.1f}%  (mass {_weights[_w]:.3f})")
    print(f"  interior (w=2,3): {_weights_pct[2]+_weights_pct[3]:.1f}%")
    print(f"  boundary (w=1,4): {_weights_pct[1]+_weights_pct[4]:.1f}%")
    print(f"  extremes (w=0,5): {_weights_pct[0]+_weights_pct[5]:.1f}%")


if __name__ == "__main__":
    main()
