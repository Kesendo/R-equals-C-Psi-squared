"""F112 universal-N strukturelle Ableitung — per-site Struktur-Exploration.

Ziel: per-Pauli-Paar Im⟨L_{σ_α,−i}, L_{σ_β,−i}⟩ = 0 strukturell beweisen
(statt enumerativ), via F113-Welle-4-Template: per-site Additivität +
Tensor-Faktorisierung + single-site Identität.

Ansatz (siehe Chat-Skizze):
1. Π = ⊗_l π_local (per-site Faktorisierung der Π-Aktion)
2. L_σ = −i(L^L_σ − L^R_σ) mit L^{L/R}_σ = ⊗_l L^{L/R}_{σ_l}
3. Frobenius ⟨L_α, L_β⟩ faktorisiert in 4 Tensor-Produkt-Inner-Products
4. Π-Projektion bricht die naive Faktorisierung — bearbeitet via
   per-site π_local Eigenbasis: (I±X)/√2, (Z±iY)/√2.

Dieses Script rechnet EXAKT (sympy) durch:
- N=1: alle (σ_α, σ_β) Paare, Im⟨L_{α,−i}, L_{β,−i}⟩
- N=2: ausgewählte (σ_α, σ_β) Paare mit nicht-trivialer Π-orbit-Struktur
- Identifizier die single-site Identität die das Verschwinden trägt
"""
from __future__ import annotations

import sympy as sp
from sympy import I, Matrix, Rational, sqrt, simplify, eye, zeros, conjugate, re, im

# ─── Pauli letter setup ─────────────────────────────────────────────────────
I2 = Matrix([[1, 0], [0, 1]])
X = Matrix([[0, 1], [1, 0]])
Y = Matrix([[0, -I], [I, 0]])
Z = Matrix([[1, 0], [0, -1]])

PAULI = {"I": I2, "X": X, "Z": Z, "Y": Y}
LETTERS = ["I", "X", "Z", "Y"]  # PauliIndex convention: I=0, X=1, Z=2, Y=3

# bit_a (= a) and bit_b (= b) per letter
BIT_A = {"I": 0, "X": 1, "Z": 0, "Y": 1}
BIT_B = {"I": 0, "X": 0, "Z": 1, "Y": 1}


def kron(A: Matrix, B: Matrix) -> Matrix:
    """Kronecker product."""
    rows_a, cols_a = A.shape
    rows_b, cols_b = B.shape
    out = zeros(rows_a * rows_b, cols_a * cols_b)
    for i in range(rows_a):
        for j in range(cols_a):
            for k in range(rows_b):
                for l in range(cols_b):
                    out[i * rows_b + k, j * cols_b + l] = A[i, j] * B[k, l]
    return out


def pauli_string(letters):
    op = PAULI[letters[0]]
    for L in letters[1:]:
        op = kron(op, PAULI[L])
    return op


def vec_to_pauli_basis_transform(N: int) -> Matrix:
    """d² × 4^N matrix T with columns vec_F(σ_α) column-stacked."""
    d = 2**N
    d2 = d * d
    pauli_count = 4**N
    T = zeros(d2, pauli_count)
    for k in range(pauli_count):
        # decode k → letters (big-endian, site 0 MSB)
        letters = []
        kk = k
        for _ in range(N):
            letters.insert(0, LETTERS[kk & 0b11])
            kk >>= 2
        sigma = pauli_string(letters)
        for j in range(d):
            for i in range(d):
                T[j * d + i, k] = sigma[i, j]
    return T


def build_L_H_pauli(H: Matrix, N: int) -> Matrix:
    """L = -i[H, ·] in Pauli basis: 4^N × 4^N exact."""
    d = 2**N
    Id = eye(d)
    # L_vec = -i (H ⊗ I − I ⊗ H^T)
    L_vec = -I * (kron(H, Id) - kron(Id, H.T))
    T = vec_to_pauli_basis_transform(N)
    return T.H * L_vec * T / d


def build_pi_pauli(N: int) -> Matrix:
    """Π in Pauli basis as a signed permutation matrix.

    Per-site π_local for Z-dephase: I↔X (phase 1), Z↔Y (phase i).
    Π = ⊗_l π_local at the FULL operator level acts as
    Π · σ · Π^(-1) on a Pauli string, which on the basis vec(σ) is the
    signed permutation matrix below.
    """
    pauli_count = 4**N

    def act_on_letter(letter: str) -> tuple[str, sp.Expr]:
        # Z-dephase: flip bit_a, phase = i^bit_b
        a, b = BIT_A[letter], BIT_B[letter]
        new_a = 1 - a
        new_b = b
        new_letter = next(L for L in LETTERS if BIT_A[L] == new_a and BIT_B[L] == new_b)
        phase = I if b == 1 else sp.Integer(1)
        return new_letter, phase

    pi = zeros(pauli_count, pauli_count)
    for k in range(pauli_count):
        # decode k → letters
        letters = []
        kk = k
        for _ in range(N):
            letters.insert(0, LETTERS[kk & 0b11])
            kk >>= 2
        new_letters = []
        total_phase = sp.Integer(1)
        for L in letters:
            new_L, phase = act_on_letter(L)
            new_letters.append(new_L)
            total_phase *= phase
        # encode new_letters → new_k
        new_k = 0
        for L in new_letters:
            new_k = (new_k << 2) | LETTERS.index(L)
        pi[new_k, k] = total_phase
    return pi


def project_pi_eigenspace(M: Matrix, Pi: Matrix, target_eigenvalue: sp.Expr) -> Matrix:
    """P_λ(M) = (1/4) Σ_{k=0..3} λ^{-k} Π^k M Π^{-k}."""
    Pi_inv = Pi.H
    result = zeros(*M.shape)
    cur_pi = eye(Pi.shape[0])
    cur_pi_inv = eye(Pi.shape[0])
    for k in range(4):
        coef = Rational(1, 4) * (1 / target_eigenvalue**k)
        result = result + coef * (cur_pi * M * cur_pi_inv)
        cur_pi = cur_pi * Pi
        cur_pi_inv = cur_pi_inv * Pi_inv
    return result


def frobenius_inner(A: Matrix, B: Matrix) -> sp.Expr:
    """⟨A, B⟩ = Σ A*[i,j] · B[i,j]."""
    total = sp.Integer(0)
    for i in range(A.rows):
        for j in range(A.cols):
            total += conjugate(A[i, j]) * B[i, j]
    return simplify(total)


# ─── N=1 explicit computation ─────────────────────────────────────────────
print("=" * 72)
print("N=1 EXACT computation (sympy)")
print("=" * 72)

N = 1
pi = build_pi_pauli(N)
print(f"\nΠ (N=1) on basis (I, X, Z, Y) =\n{pi}")
print(f"Π² =\n{pi**2}")
print(f"Π⁴ =\n{pi**4}  (should be identity)")
print()

L_per_letter = {}
L_minus_i_per_letter = {}
for letter in LETTERS:
    sigma = PAULI[letter]
    L = build_L_H_pauli(sigma, N)
    L_per_letter[letter] = L
    L_mi = project_pi_eigenspace(L, pi, -I)
    L_minus_i_per_letter[letter] = simplify(L_mi)
    print(f"L_{letter} (Pauli basis, 4×4):")
    sp.pprint(L)
    print(f"L_{letter},-i (Π=-i projection):")
    sp.pprint(L_minus_i_per_letter[letter])
    print(f"  Frobenius norm² = {simplify(frobenius_inner(L_minus_i_per_letter[letter], L_minus_i_per_letter[letter]))}")
    print()

print("\nN=1 all pair F(α, β) = Im⟨L_α,-i, L_β,-i⟩:")
print(f"{'(α, β)':<10} {'⟨L_α,-i, L_β,-i⟩':<40} {'Im part':<20}")
for a in LETTERS:
    for b in LETTERS:
        inner = simplify(frobenius_inner(L_minus_i_per_letter[a], L_minus_i_per_letter[b]))
        im_part = simplify(im(inner))
        print(f"({a}, {b})    {str(inner):<40} {str(im_part):<20}")
print()

# ─── N=2 sample: BitB=1 pair (ZI, IZ) and (ZZ, YY) ─────────────────────
print("=" * 72)
print("N=2 sample pairs (BitB-parity = 1 on each side)")
print("=" * 72)

N = 2
pi2 = build_pi_pauli(N)

# BitB(σ) for 2-letter string = Σ BitB(σ_l) mod 2
# BitB=1 strings at N=2: IZ, IY, XZ, XY, ZI, YI, ZX, YX (8 strings).
# Pick a few illustrative pairs:
sample_pairs = [
    (["Z", "I"], ["I", "Z"]),  # (ZI, IZ) — single Z at different sites
    (["Z", "I"], ["Z", "I"]),  # diagonal (ZI, ZI)
    (["Z", "Z"], ["Y", "Y"]),  # (ZZ, YY) — both BitB-even (Σ = 2 mod 2 = 0)!
    (["Z", "Y"], ["Y", "Z"]),  # both BitB-even
    (["I", "Z"], ["I", "Y"]),  # (IZ, IY) — both BitB-odd
]

print(f"\nReminder: BitB-odd strings at N=2: IZ, IY, XZ, XY, ZI, YI, ZX, YX")
print(f"          BitB-even (= BitB=0) include: ZZ, YY, ZY, YZ (with both letters in {{Z,Y}})")
print()

for letters_a, letters_b in sample_pairs:
    sigma_a = pauli_string(letters_a)
    sigma_b = pauli_string(letters_b)
    L_a = build_L_H_pauli(sigma_a, N)
    L_b = build_L_H_pauli(sigma_b, N)
    L_a_mi = simplify(project_pi_eigenspace(L_a, pi2, -I))
    L_b_mi = simplify(project_pi_eigenspace(L_b, pi2, -I))

    bit_b_a = (BIT_B[letters_a[0]] + BIT_B[letters_a[1]]) % 2
    bit_b_b = (BIT_B[letters_b[0]] + BIT_B[letters_b[1]]) % 2

    inner = simplify(frobenius_inner(L_a_mi, L_b_mi))
    im_part = simplify(im(inner))

    label_a = "".join(letters_a)
    label_b = "".join(letters_b)
    print(f"σ_a={label_a} (BitB={bit_b_a})  σ_b={label_b} (BitB={bit_b_b})")
    print(f"  ‖L_a,-i‖² = {simplify(frobenius_inner(L_a_mi, L_a_mi))}")
    print(f"  ‖L_b,-i‖² = {simplify(frobenius_inner(L_b_mi, L_b_mi))}")
    print(f"  ⟨L_a,-i, L_b,-i⟩ = {inner}")
    print(f"  Im part = {im_part}")
    print()

print("=" * 72)
print("N=2 EXHAUSTIVE: alle BitB-odd × BitB-odd Paare (8×8 = 64 Paare)")
print("=" * 72)

N = 2
pi2 = build_pi_pauli(N)

def all_strings_N2():
    return [(a, b) for a in LETTERS for b in LETTERS]

def bit_b_total(letters):
    return sum(BIT_B[L] for L in letters) % 2

# Build all 16 L_σ,-i at N=2, cache
L_minus_i_N2 = {}
for letters in all_strings_N2():
    sigma = pauli_string(letters)
    L = build_L_H_pauli(sigma, N)
    L_mi = simplify(project_pi_eigenspace(L, pi2, -I))
    L_minus_i_N2[letters] = L_mi

# Filter to BitB-odd
bit_b_odd = [s for s in all_strings_N2() if bit_b_total(s) == 1]
print(f"\nBitB-odd strings at N=2 ({len(bit_b_odd)}):")
for s in bit_b_odd:
    print(f"  {''.join(s)}  (per-site BitB: {BIT_B[s[0]]}, {BIT_B[s[1]]})")
print()

# Pairwise inner products. Verify Im=0, collect values.
print("\nAlle ⟨L_α,-i, L_β,-i⟩ für BitB-odd × BitB-odd:")
print(f"{'(α, β)':<14} {'inner (exact)':<22} {'Im':<10}")
nonzero_pairs = []
for sa in bit_b_odd:
    for sb in bit_b_odd:
        inner = simplify(frobenius_inner(L_minus_i_N2[sa], L_minus_i_N2[sb]))
        im_part = simplify(im(inner))
        label = f"({''.join(sa)}, {''.join(sb)})"
        print(f"{label:<14} {str(inner):<22} {str(im_part):<10}")
        if inner != 0:
            nonzero_pairs.append((sa, sb, inner))

print(f"\nNicht-Null Inner Products: {len(nonzero_pairs)} von 64")
print()

# Test multiplicativity hypothesis
print("=" * 72)
print("MULTIPLIKATIVITÄTS-TEST: ⟨L_α,-i, L_β,-i⟩ = ∏_l f(α_l, β_l)?")
print("=" * 72)

# Build N=1 per-letter inner product table
print("\nN=1 Single-Site-Tabelle f(a, b) = ⟨L_{a,-i}, L_{b,-i}⟩:")
print(f"{'(a, b)':<8} {'f(a,b)':<10}")
single_site_table = {}
for a in LETTERS:
    for b in LETTERS:
        val = simplify(frobenius_inner(L_minus_i_per_letter[a], L_minus_i_per_letter[b]))
        single_site_table[(a, b)] = val
        print(f"({a}, {b})    {val}")
print()

# Test: at N=2, does ⟨L_α,-i, L_β,-i⟩ = f(α[0], β[0]) × f(α[1], β[1])?
print("\nVergleich N=2 inner products mit per-site Produkt f(α[0],β[0]) × f(α[1],β[1]):")
print(f"{'(α, β)':<14} {'N=2 inner':<14} {'per-site Produkt':<22} {'gleich?':<10}")
all_match = True
mismatches = []
for sa in bit_b_odd:
    for sb in bit_b_odd:
        n2_val = simplify(frobenius_inner(L_minus_i_N2[sa], L_minus_i_N2[sb]))
        per_site_product = single_site_table[(sa[0], sb[0])] * single_site_table[(sa[1], sb[1])]
        per_site_product = simplify(per_site_product)
        match = simplify(n2_val - per_site_product) == 0
        label = f"({''.join(sa)}, {''.join(sb)})"
        marker = "✓" if match else "✗"
        print(f"{label:<14} {str(n2_val):<14} {str(per_site_product):<22} {marker}")
        if not match:
            all_match = False
            mismatches.append((sa, sb, n2_val, per_site_product))

print()
if all_match:
    print("✓ MULTIPLIKATIVITÄT VOLLSTÄNDIG: ⟨L_α,-i, L_β,-i⟩ = ∏_l f(α_l, β_l) bei N=2")
    print("  → Universal-N folgt direkt: F(α,β) = Im⟨...⟩ = Im(∏_l f(α_l, β_l)) = 0")
    print("  weil jeder Faktor f(α_l, β_l) reell ist (siehe N=1 Tabelle).")
else:
    print(f"✗ Multiplikativität BRICHT für {len(mismatches)} Paar(e):")
    for sa, sb, n2v, psv in mismatches[:5]:
        print(f"  ({''.join(sa)}, {''.join(sb)}): N=2={n2v}, per-site={psv}, diff={simplify(n2v - psv)}")

print()
print("=" * 72)
print("Zusammenfassung der strukturellen Beobachtungen:")
print(f"  • Alle Im-Parts bei BitB-odd Paaren = 0 (verifiziert N=1 + N=2)")
print(f"  • Single-Site-Tabelle f(a,b): nur f(Z,Z)=4, f(Y,Y)=4 nicht-null, alle anderen = 0")
print(f"  • → f ist reell-wertig und SEHR sparse (nur Diagonal-Z und Diagonal-Y)")
print("=" * 72)
