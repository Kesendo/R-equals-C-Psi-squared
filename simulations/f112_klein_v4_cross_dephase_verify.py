"""F112 cross-dephase (X / Y) verifier: Welle 13.

For each dephase letter d in {X, Y, Z}, F112-d asks: for which class of
Hermitian H and which class of bit_?-homogeneous c does the
`polarity_coordinates_from_L` asymmetry (computed against Π_d) vanish
exactly?

The Welle 11 closure handles d = Z (and Hermitian H + bit_b-homogeneous c).
This script tests the natural cross-dephase conjectures:

  Conjecture F112-Y: Hermitian H + bit_b-homogeneous c gives Π_Y-polarity
    asymmetry = 0 (same bit_b hypothesis as F112-Z because Π_Y² uses bit_b
    too).

  Conjecture F112-X: Hermitian H + bit_a-homogeneous c gives Π_X-polarity
    asymmetry = 0 (bit_a hypothesis because Π_X² uses bit_a).

We also test the Hadamard-conjugation route for F112-X:

  Conjecture F112-X-Hadamard: given any (Hermitian H_Z, bit_b-homogeneous
    c_Z) F112-Z configuration, rotate H' = U^⊗N H U^⊗N and c' = U^⊗N c U^⊗N
    where U is the single-qubit Hadamard (Hilbert space). Then (H', c')
    should give Π_X-polarity asymmetry = 0 via Klein-V₄ equivariance.
    Verify the rotated c is bit_a-homogeneous (Hadamard swaps X ↔ Z at
    Pauli-letter level, but leaves Y and I; check that bit_b parity of
    original c equals bit_a parity of rotated c).

We test at N = 2, 3 with 10 random configurations each, plus a few
explicitly chosen edge cases (Y-only c, mixed c, complex coefficients on
Pauli strings of c).

Output: pass / fail summary table per conjecture + a couple of failing
configurations if any conjecture fails.
"""
from __future__ import annotations

import sys
from itertools import product
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from framework.pauli import _vec_to_pauli_basis_transform
from framework.symmetry import build_pi_full
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}
LETTERS = ("I", "X", "Y", "Z")
BIT_A = {"I": 0, "X": 1, "Y": 1, "Z": 0}
BIT_B = {"I": 0, "X": 0, "Y": 1, "Z": 1}

# Single-qubit Hadamard: U X U = Z, U Z U = X, U Y U = -Y.
U_HADAMARD = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

TOL = 1e-9


def pauli_string(letters):
    op = PAULI[letters[0]]
    for L in letters[1:]:
        op = np.kron(op, PAULI[L])
    return op


def kron_n(op, N):
    out = op
    for _ in range(N - 1):
        out = np.kron(out, op)
    return out


def random_hermitian_pauli_h(N, rng, n_terms=4):
    """Random Hermitian H as a real-coefficient sum of Pauli strings."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    seen = set()
    for _ in range(n_terms * 3):
        if len(seen) >= n_terms:
            break
        s = tuple(rng.choice(LETTERS) for _ in range(N))
        if s in seen or all(c == "I" for c in s):
            continue
        seen.add(s)
        c = rng.normal()
        H = H + c * pauli_string(s)
    return H


def random_non_hermitian_pauli_h(N, rng, n_terms=4):
    """Random non-Hermitian H as a COMPLEX-coefficient sum of Pauli strings.
    The Welle 11 non-Hermitian extension covers this case under bit_b homo c."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    seen = set()
    for _ in range(n_terms * 3):
        if len(seen) >= n_terms:
            break
        s = tuple(rng.choice(LETTERS) for _ in range(N))
        if s in seen or all(c == "I" for c in s):
            continue
        seen.add(s)
        coef = rng.normal() + 1j * rng.normal()
        H = H + coef * pauli_string(s)
    return H


def random_homogeneous_c(N, rng, parity_axis, parity_value, n_terms=3, complex_coefs=True):
    """Random c with all Pauli strings sharing the given parity (bit_a or bit_b)
    at the given value (0 or 1)."""
    assert parity_axis in ("bit_a", "bit_b")
    parity_fn = (lambda s: sum(BIT_A[c] for c in s) % 2) if parity_axis == "bit_a" \
        else (lambda s: sum(BIT_B[c] for c in s) % 2)
    d = 2 ** N
    c = np.zeros((d, d), dtype=complex)
    seen = set()
    attempts = 0
    while len(seen) < n_terms and attempts < 1000:
        attempts += 1
        s = tuple(rng.choice(LETTERS) for _ in range(N))
        if s in seen or all(c0 == "I" for c0 in s):
            continue
        if parity_fn(s) != parity_value:
            continue
        seen.add(s)
        coef = (rng.normal() + 1j * rng.normal()) if complex_coefs else rng.normal()
        c = c + coef * pauli_string(s)
    if len(seen) < 1:
        # Fallback: pick any non-identity string with correct parity.
        for s in product(LETTERS, repeat=N):
            if all(c0 == "I" for c0 in s):
                continue
            if parity_fn(s) == parity_value:
                c = pauli_string(s).astype(complex)
                seen.add(s)
                break
    return c, list(seen)


def build_lindblad_pauli_basis(H, c_ops, gammas, N):
    """Build the standard Lindblad L = -i[H,·] + Σγ_k D[c_k] in vec basis
    then transform to Pauli basis. Identical to the construction in
    polarity_coordinates_from_hc but kept inline for clarity."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for c, g in zip(c_ops, gammas):
        c_dag_c = c.conj().T @ c
        anti = 0.5 * (np.kron(c_dag_c, Id) + np.kron(Id, c_dag_c.T))
        L_vec = L_vec + g * (np.kron(c, c.conj()) - anti)
    T = _vec_to_pauli_basis_transform(N)
    L_pauli = (T.conj().T @ L_vec @ T) / (2 ** N)
    return L_pauli


def asymmetry_for_dephase(H, c_ops, gammas, N, dephase_letter):
    """Compute the polarity_coordinates_from_L asymmetry against Π_d."""
    L_pauli = build_lindblad_pauli_basis(H, c_ops, gammas, N)
    sigma = float(np.real(sum(gammas)))
    Pi = build_pi_full(N, dephase_letter=dephase_letter)
    result = polarity_coordinates_from_L(L_pauli, N, sigma, Pi=Pi)
    return result["asymmetry"], result["norm_sq"]


def rotate_by_hadamard(H, N):
    """Apply Hadamard ⊗N conjugation to a Hilbert-space operator."""
    U_n = kron_n(U_HADAMARD, N)
    return U_n @ H @ U_n.conj().T


def pauli_strings_for_op(op, N, tol=1e-10):
    """Decompose op in Pauli basis; return list of (letters, coefficient)
    for non-zero strings."""
    d = 2 ** N
    out = []
    for letters in product(LETTERS, repeat=N):
        P = pauli_string(letters)
        coef = np.trace(P.conj().T @ op) / d
        if abs(coef) > tol:
            out.append((letters, complex(coef)))
    return out


def parity_set(letters_list, parity_fn):
    return {parity_fn(s) for (s, _) in letters_list}


# ---------------------------------------------------------------------------
# Sub-task 1.1: F112-Y direct check
# ---------------------------------------------------------------------------

def verify_f112_y_direct(N, n_configs=10, seed=42, non_hermitian=False):
    """For random (Hermitian or non-Hermitian H, bit_b-homogeneous c)
    configurations, build L using standard Lindblad construction, then
    evaluate Π_Y-polarity asymmetry."""
    rng = np.random.default_rng(seed)
    results = []
    max_abs = 0.0
    for cfg in range(n_configs):
        H = (random_non_hermitian_pauli_h(N, rng, n_terms=3) if non_hermitian
             else random_hermitian_pauli_h(N, rng, n_terms=3))
        parity_value = int(rng.integers(0, 2))
        c, seen = random_homogeneous_c(N, rng, "bit_b", parity_value,
                                        n_terms=3, complex_coefs=True)
        gammas = [float(rng.uniform(0.1, 1.0))]
        asym, norms = asymmetry_for_dephase(H, [c], gammas, N, "Y")
        max_abs = max(max_abs, abs(asym))
        results.append({
            "cfg": cfg,
            "parity_value": parity_value,
            "c_terms": seen,
            "asymmetry": asym,
            "abs": abs(asym),
        })
    return results, max_abs


# ---------------------------------------------------------------------------
# Sub-task 1.2: F112-X direct check (bit_a-homogeneous c)
# ---------------------------------------------------------------------------

def verify_f112_x_direct(N, n_configs=10, seed=43, non_hermitian=False):
    """For random (Hermitian or non-Hermitian H, bit_a-homogeneous c)
    configurations, build L using standard Lindblad construction, then
    evaluate Π_X-polarity asymmetry."""
    rng = np.random.default_rng(seed)
    results = []
    max_abs = 0.0
    for cfg in range(n_configs):
        H = (random_non_hermitian_pauli_h(N, rng, n_terms=3) if non_hermitian
             else random_hermitian_pauli_h(N, rng, n_terms=3))
        parity_value = int(rng.integers(0, 2))
        c, seen = random_homogeneous_c(N, rng, "bit_a", parity_value,
                                        n_terms=3, complex_coefs=True)
        gammas = [float(rng.uniform(0.1, 1.0))]
        asym, norms = asymmetry_for_dephase(H, [c], gammas, N, "X")
        max_abs = max(max_abs, abs(asym))
        results.append({
            "cfg": cfg,
            "parity_value": parity_value,
            "c_terms": seen,
            "asymmetry": asym,
            "abs": abs(asym),
        })
    return results, max_abs


# ---------------------------------------------------------------------------
# Sub-task 1.3: F112-X via Hadamard check
# ---------------------------------------------------------------------------

def verify_f112_x_via_hadamard(N, n_configs=10, seed=44):
    """Start from (Hermitian H_Z, bit_b-homogeneous c_Z). Rotate to (H', c')
    by Hadamard^⊗N. Check that c' is bit_a-homogeneous (per the Hadamard
    swap X↔Z, leaving I and Y fixed up to sign). Then build L_X using the
    rotated operators and measure Π_X-polarity asymmetry."""
    rng = np.random.default_rng(seed)
    results = []
    max_abs = 0.0
    homog_failures = 0
    for cfg in range(n_configs):
        H = random_hermitian_pauli_h(N, rng, n_terms=3)
        parity_value = int(rng.integers(0, 2))
        c, seen = random_homogeneous_c(N, rng, "bit_b", parity_value,
                                        n_terms=3, complex_coefs=True)
        gammas = [float(rng.uniform(0.1, 1.0))]
        # Rotate to Hadamard-basis operators.
        H_prime = rotate_by_hadamard(H, N)
        c_prime = rotate_by_hadamard(c, N)
        # Check bit_a-homogeneity of c_prime (post-Hadamard).
        c_prime_strings = pauli_strings_for_op(c_prime, N)
        bit_a_parities = parity_set(c_prime_strings,
                                     lambda s: sum(BIT_A[L] for L in s) % 2)
        bit_b_parities_prime = parity_set(c_prime_strings,
                                           lambda s: sum(BIT_B[L] for L in s) % 2)
        is_bit_a_homo = len(bit_a_parities) <= 1
        if not is_bit_a_homo:
            homog_failures += 1
        # Measure asymmetry against Π_X.
        asym, norms = asymmetry_for_dephase(H_prime, [c_prime], gammas, N, "X")
        max_abs = max(max_abs, abs(asym))
        results.append({
            "cfg": cfg,
            "orig_bit_b": parity_value,
            "orig_c_terms": seen,
            "post_hadamard_bit_a_parities": sorted(bit_a_parities),
            "post_hadamard_bit_b_parities": sorted(bit_b_parities_prime),
            "is_bit_a_homo_post": is_bit_a_homo,
            "asymmetry": asym,
            "abs": abs(asym),
        })
    return results, max_abs, homog_failures


# ---------------------------------------------------------------------------
# Sub-task 1.4: Edge cases - what about MIXED bit_a/bit_b?
# ---------------------------------------------------------------------------

def edge_case_table(N):
    """Test a handful of explicit edge cases, using COMPLEX coefficients to
    expose the bit_b / bit_a homogeneity structure (per probe 11: real-only
    coefficients are accidentally always-balanced).
    """
    rng = np.random.default_rng(99)
    H = random_hermitian_pauli_h(N, rng, n_terms=3)
    gammas = [0.4]
    cases = []
    site = 0

    def site_op(letter):
        if N == 1:
            return PAULI[letter]
        ops = [I2] * N
        ops[site] = PAULI[letter]
        op = ops[0]
        for o in ops[1:]:
            op = np.kron(op, o)
        return op

    # Pure single-Pauli c (trivially homogeneous on every axis)
    for letter in ("X", "Y", "Z"):
        c = site_op(letter)
        row = {"c": f"site-{site} {letter}"}
        for dl in ("X", "Y", "Z"):
            asym, _ = asymmetry_for_dephase(H, [c], gammas, N, dl)
            row[f"asym_pi_{dl}"] = asym
        cases.append(row)

    # Mixed c, real coefficients (probe-11 "accidentally balanced" regime)
    c_real_mixed = site_op("X") + site_op("Z")
    row = {"c": f"(X+Z) on site-{site}, real coefs (probe-11 reg)"}
    for dl in ("X", "Y", "Z"):
        asym, _ = asymmetry_for_dephase(H, [c_real_mixed], gammas, N, dl)
        row[f"asym_pi_{dl}"] = asym
    cases.append(row)

    # Mixed bit_b, COMPLEX coefficients (the real structural axis)
    # c = (1+2j) X_site + (3+4j) Z_site mixes bit_b (X bit_b=0, Z bit_b=1)
    c_complex_mixed_bitb = (1 + 2j) * site_op("X") + (3 + 4j) * site_op("Z")
    row = {"c": f"(1+2j)X + (3+4j)Z, mixed bit_b complex"}
    for dl in ("X", "Y", "Z"):
        asym, _ = asymmetry_for_dephase(H, [c_complex_mixed_bitb], gammas, N, dl)
        row[f"asym_pi_{dl}"] = asym
    cases.append(row)

    # Mixed bit_a, complex coefficients
    # c = (1+2j) X_site + (3+4j) Z_site has bit_a(X)=1, bit_a(Z)=0 → also mixed bit_a
    # Use c = (1+2j) X + (3+4j) Y to keep bit_b mixed but bit_a homogeneous
    # X has bit_b=0, Y has bit_b=1 (mixed); X has bit_a=1, Y has bit_a=1 (homogeneous bit_a)
    c_complex_homo_bita = (1 + 2j) * site_op("X") + (3 + 4j) * site_op("Y")
    row = {"c": f"(1+2j)X + (3+4j)Y, mixed bit_b, homo bit_a"}
    for dl in ("X", "Y", "Z"):
        asym, _ = asymmetry_for_dephase(H, [c_complex_homo_bita], gammas, N, dl)
        row[f"asym_pi_{dl}"] = asym
    cases.append(row)

    # Homogeneous bit_b, mixed bit_a, complex coefficients
    # Y has bit_b=1, Z has bit_b=1 (homo bit_b); Y bit_a=1, Z bit_a=0 (mixed bit_a)
    c_complex_homo_bitb = (1 + 2j) * site_op("Y") + (3 + 4j) * site_op("Z")
    row = {"c": f"(1+2j)Y + (3+4j)Z, homo bit_b, mixed bit_a"}
    for dl in ("X", "Y", "Z"):
        asym, _ = asymmetry_for_dephase(H, [c_complex_homo_bitb], gammas, N, dl)
        row[f"asym_pi_{dl}"] = asym
    cases.append(row)

    return cases


# ---------------------------------------------------------------------------
# Structural per-pair F-identity verifier for cross-dephase Π
# ---------------------------------------------------------------------------
# Mirror the Welle-11 lemmas N-A, N-B for Π_X and Π_Y:
#   - For Π_Y (graded by bit_b like Π_Z, but with phase -i instead of +i),
#     check ||L_{σ,-i}||² = 4^N for bit_b-odd σ and <L_{α,-i}, L_{β,-i}> = 0
#     for α ≠ β both bit_b-odd.
#   - For Π_X (graded by bit_a, phase -i^bit_a), check the same identities
#     with bit_a-odd σ replacing bit_b-odd σ.


def build_L_sigma_pauli_basis(sigma, N):
    """L_σ = -i [σ, ·] in Pauli basis (same construction as Welle-11 verifier)."""
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    T = _vec_to_pauli_basis_transform(N)
    L_vec = -1j * (np.kron(sigma, Id) - np.kron(Id, sigma.T))
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def project_minus_i_via_pi(M, Pi):
    """Project operator-space matrix M onto Π-conjugation eigenvalue -i
    subspace via the standard 4-term spectral formula."""
    Pi_inv = Pi.conj().T
    out = np.zeros_like(M)
    cur_pi = np.eye(Pi.shape[0], dtype=complex)
    cur_pi_inv = np.eye(Pi.shape[0], dtype=complex)
    for k in range(4):
        coef = ((-1j) ** (-k)) / 4.0
        out = out + coef * (cur_pi @ M @ cur_pi_inv)
        cur_pi = cur_pi @ Pi
        cur_pi_inv = cur_pi_inv @ Pi_inv
    return out


def frob(A, B):
    return complex(np.sum(A.conj() * B))


def all_strings(N):
    return list(product(LETTERS, repeat=N))


def bit_parity_total(s, axis):
    """Sum bit_axis(L) mod 2 over a Pauli string."""
    tbl = BIT_A if axis == "bit_a" else BIT_B
    return sum(tbl[L] for L in s) % 2


def verify_structural_identity_for_dephase(N, dephase_letter):
    """Run the Welle-11 Lemma-N-A and N-B style structural checks for the
    given dephase letter.

    For d in {Z, Y}: Π_d² grades by bit_b → check L_σ with σ bit_b-odd.
    For d = X:        Π_X² grades by bit_a → check L_σ with σ bit_a-odd.

    Reports:
      - Step A.2 max deviation: Π² L_σ Π⁻² + L_σ for σ in the relevant odd sector.
      - Step A.3 max deviation: <L_σ, Π L_σ Π⁻¹>.
      - Lemma N-A: ||L_{σ,-i}||² = 4^N for σ in odd sector.
      - Lemma N-B: <L_{α,-i}, L_{β,-i}> = 0 for α ≠ β both odd-sector.
      - Main theorem (all pairs, including mixed-sector): max|Im F| = 0.
    """
    Pi = build_pi_full(N, dephase_letter=dephase_letter)
    Pi_inv = Pi.conj().T

    axis = "bit_a" if dephase_letter == "X" else "bit_b"
    odd_strings = [s for s in all_strings(N) if bit_parity_total(s, axis) == 1]
    all_str = all_strings(N)

    L_cache = {s: build_L_sigma_pauli_basis(pauli_string(s), N) for s in all_str}
    L_mi_cache = {s: project_minus_i_via_pi(L_cache[s], Pi) for s in all_str}

    expected_raw = 2 * (4 ** N)
    expected_proj = 4 ** N

    max_step_a2 = 0.0
    max_step_a3 = 0.0
    max_norm_raw_dev = 0.0
    max_norm_proj_dev = 0.0
    for s in odd_strings:
        L = L_cache[s]
        # Step A.1: ||L_σ||² = 2·4^N
        raw = frob(L, L).real
        max_norm_raw_dev = max(max_norm_raw_dev, abs(raw - expected_raw))
        # Step A.2: Π² L Π⁻² = -L (for sectors graded by axis with σ odd)
        Pi2 = Pi @ Pi
        L_pi2 = Pi2 @ L @ Pi2.conj().T
        max_step_a2 = max(max_step_a2, float(np.max(np.abs(L_pi2 + L))))
        # Step A.3: <L_σ, Π L_σ Π⁻¹> = 0
        L_pi = Pi @ L @ Pi_inv
        cross = frob(L, L_pi)
        max_step_a3 = max(max_step_a3, abs(cross))
        # Lemma N-A: ||L_{σ,-i}||² = 4^N
        L_mi = L_mi_cache[s]
        proj = frob(L_mi, L_mi).real
        max_norm_proj_dev = max(max_norm_proj_dev, abs(proj - expected_proj))

    # Lemma N-B: <L_{α,-i}, L_{β,-i}> = 0 for α ≠ β, both odd-sector
    max_b = 0.0
    pair_count = 0
    for i, sa in enumerate(odd_strings):
        for j, sb in enumerate(odd_strings):
            if sa == sb:
                continue
            ip = frob(L_mi_cache[sa], L_mi_cache[sb])
            max_b = max(max_b, abs(ip))
            pair_count += 1

    # Main theorem: max|Im F(α, β)| = 0 over all string pairs (including
    # bit-even sector where L_{σ,-i} = 0).
    max_im = 0.0
    for sa in all_str:
        for sb in all_str:
            ip = frob(L_mi_cache[sa], L_mi_cache[sb])
            max_im = max(max_im, abs(ip.imag))

    return {
        "dephase_letter": dephase_letter,
        "axis": axis,
        "n_odd_strings": len(odd_strings),
        "n_off_diag_pairs": pair_count,
        "n_all_pairs": len(all_str) ** 2,
        "max_norm_raw_dev": max_norm_raw_dev,
        "max_step_a2": max_step_a2,
        "max_step_a3": max_step_a3,
        "max_norm_proj_dev": max_norm_proj_dev,
        "max_lemma_b": max_b,
        "max_im_main": max_im,
    }


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("F112 cross-dephase verifier (Welle 13)")
    print("Conjectures:")
    print("  F112-Z (proven Welle 11): Hermitian H + bit_b-homo c → Π_Z asym = 0")
    print("  F112-Y (conjecture):      Hermitian H + bit_b-homo c → Π_Y asym = 0")
    print("                            (same bit_b axis because Π_Y² = Π_Z² grading)")
    print("  F112-X (conjecture):      Hermitian H + bit_a-homo c → Π_X asym = 0")
    print("                            (Π_X² grades by bit_a, not bit_b)")
    print("  F112-X via Hadamard:      Rotate F112-Z config → predict F112-X")
    print("  Non-Hermitian extensions: per Welle 11, the Π-eigenspace structural")
    print("                            lemmas extend the result to non-Hermitian H")
    print("=" * 78)

    overall_pass = True

    for N in (2, 3):
        print()
        print(f"--- N = {N} ---")

        # Sub-task 1.1a: F112-Y direct check (Hermitian H).
        res_y, max_y = verify_f112_y_direct(N, n_configs=10, seed=42 + N,
                                              non_hermitian=False)
        print(f"F112-Y direct (Hermitian H + bit_b-homo c, Π_Y polarity):")
        print(f"  10 random configs, max|asymmetry| = {max_y:.3e}",
              "  PASS" if max_y < TOL else "  FAIL")
        overall_pass = overall_pass and max_y < TOL

        # Sub-task 1.1b: F112-Y direct check (non-Hermitian H).
        res_y_nh, max_y_nh = verify_f112_y_direct(N, n_configs=10, seed=72 + N,
                                                    non_hermitian=True)
        print(f"F112-Y direct (non-Hermitian H + bit_b-homo c, Π_Y polarity):")
        print(f"  10 random configs, max|asymmetry| = {max_y_nh:.3e}",
              "  PASS" if max_y_nh < TOL else "  FAIL")
        overall_pass = overall_pass and max_y_nh < TOL

        # Sub-task 1.2a: F112-X direct check (Hermitian H).
        res_x, max_x = verify_f112_x_direct(N, n_configs=10, seed=143 + N,
                                              non_hermitian=False)
        print(f"F112-X direct (Hermitian H + bit_a-homo c, Π_X polarity):")
        print(f"  10 random configs, max|asymmetry| = {max_x:.3e}",
              "  PASS" if max_x < TOL else "  FAIL")
        overall_pass = overall_pass and max_x < TOL

        # Sub-task 1.2b: F112-X direct check (non-Hermitian H).
        res_x_nh, max_x_nh = verify_f112_x_direct(N, n_configs=10, seed=173 + N,
                                                    non_hermitian=True)
        print(f"F112-X direct (non-Hermitian H + bit_a-homo c, Π_X polarity):")
        print(f"  10 random configs, max|asymmetry| = {max_x_nh:.3e}",
              "  PASS" if max_x_nh < TOL else "  FAIL")
        overall_pass = overall_pass and max_x_nh < TOL

        # Sub-task 1.3: F112-X via Hadamard.
        res_h, max_h, homog_fail = verify_f112_x_via_hadamard(N, n_configs=10, seed=244 + N)
        print(f"F112-X via Hadamard (Z-config + rotation → Π_X):")
        print(f"  10 random configs, max|asymmetry| = {max_h:.3e}, "
              f"bit_a-homogeneity holds post-Hadamard = {10 - homog_fail}/10",
              "  PASS" if max_h < TOL else "  FAIL")
        overall_pass = overall_pass and max_h < TOL

        # Sub-task 1.4: Edge case table (witnesses the structural picture).
        print(f"Edge cases at N = {N} (witnesses structural axis):")
        edge = edge_case_table(N)
        for row in edge:
            label = row["c"]
            asyms = " ".join(f"Π_{k.split('_')[-1]}={v:+.3e}"
                              for k, v in row.items() if k.startswith("asym_"))
            print(f"  c = {label:50s}  {asyms}")

        # Sub-task 1.5: Structural per-pair F-identity (Welle-11 lemma analogs).
        print(f"Welle-11 lemma analogs at N = {N}:")
        for dl in ("Z", "Y", "X"):
            r = verify_structural_identity_for_dephase(N, dl)
            ok = (r["max_step_a2"] < TOL and r["max_step_a3"] < TOL
                  and r["max_norm_proj_dev"] < TOL and r["max_lemma_b"] < TOL
                  and r["max_im_main"] < TOL)
            print(f"  Π_{dl} (axis={r['axis']}, "
                  f"n_odd={r['n_odd_strings']}, "
                  f"n_offdiag={r['n_off_diag_pairs']}, "
                  f"n_all={r['n_all_pairs']}):")
            print(f"    Step A.2 |Π² L Π⁻² + L| max = {r['max_step_a2']:.3e}")
            print(f"    Step A.3 |<L, Π L Π⁻¹>|  max = {r['max_step_a3']:.3e}")
            print(f"    Lemma N-A |‖L_{{,-i}}‖² − 4^N| max = {r['max_norm_proj_dev']:.3e}")
            print(f"    Lemma N-B |<L_{{α,-i}}, L_{{β,-i}}>| max = {r['max_lemma_b']:.3e}")
            print(f"    Main theorem |Im F| max = {r['max_im_main']:.3e}",
                  "  PASS" if ok else "  FAIL")
            overall_pass = overall_pass and ok

    print()
    print("=" * 78)
    if overall_pass:
        print("ALL CHECKS PASSED at N = 2, 3 within tolerance 1e-9.")
        print("F112-Y and F112-X hold under their natural axis hypotheses,")
        print("for both Hermitian and non-Hermitian H, mirroring the F112-Z")
        print("Welle 11 closure. The structural Welle-11 lemmas (N-A, N-B)")
        print("transfer directly to Π_Y (via the same bit_b grading) and to")
        print("Π_X (via the bit_a grading, with the proof argument unchanged).")
    else:
        print("SOME CHECKS FAILED. Inspect per-step output above.")
        sys.exit(1)
    print("=" * 78)


if __name__ == "__main__":
    main()
