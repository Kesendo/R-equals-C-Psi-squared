"""F49 cross-term `‖{L_H, L_Dc}‖²` under non-uniform γ — Phase 1 verification + Phase 2 assertions.

Phase 1 (commit 1c6701c, exploratory): distinguished three hypotheses for the
planning-agent's candidate non-uniform F49 closed form, which predicts 204.8 at
N=3 Heisenberg γ=[0.1, 0.2, 0.3] but does not match the previously-reported truth
163.84:

  (i)  N=3-only artifact (overlapping bonds at site 1 violate "disjoint bond
       supports" of PROOF_CROSS_TERM_FORMULA Lemma 3 corollary).
  (ii) General formula gap (architect's G(bond, H) is wrong at all N).
  (iii) Centering mismatch (truth 163.84 used a different L_Dc centering than
       the architect's F1-centered L_D + Σγ·I assumption).

Phase 1 conclusion: none of (i), (ii), (iii) applies. The 204.8 value was a hand-calc
slip (‖L_H^bond‖² is 384 at N=3 Heisenberg J=1, not 480); the closed form is correct
at the F1-centered L_Dc, all tested N ∈ {3, 4, 5}, all four H classes (Heisenberg,
Ising, XY, soft XY+YX).

Phase 2 (this commit) bolts assertion logic onto the exploratory print structure:
each (N, H, γ) row asserts |candidate − truth| < 1e-10, and a final assertion-block
summary prints "All N=3,4,5 × 4-H-classes Phase 2 closed-form formula verified
bit-exact". Script exits 0 on success; AssertionError surfaces any future regression.

Sections:
  1. Multi-centering at N=3 Heisenberg γ=[0.1, 0.2, 0.3]: four candidate L_Dc
     centerings; identify which equals 163.84.
  2. Architect's candidate formula at correct centering for the N=3 anchor.
  3. Scan N=3, 4, 5 chain Heisenberg with γ_l = 0.05·(l+1).
  4. Cross-H-class sanity at N=4 (Heisenberg / Ising / XY / XY+YX), testing
     G(bond, H) decomposition.
  5. If hypothesis (iii) is confirmed: report corrected formula under the
     identified centering.
  6. (Phase 2) Final assertion block.

Output: stdout only, structured per section. Phase 1 narrative preserved; Phase 2
assertions added inline at scan rows + final block.

Cross-references:
  - docs/proofs/PROOF_F49_NONUNIFORM_GAMMA_EXTENSION.md
  - compute/RCPsiSquared.Core/F1/F49NonUniformCrossTermClaim.cs
"""
from __future__ import annotations

import sys

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# --------------------------------------------------------------------------- #
# Pauli + superoperator primitives                                            #
# Matches the convention of simulations/cross_term_formula_check.py so the    #
# verify is self-contained and not coupled to private framework helpers.      #
# --------------------------------------------------------------------------- #


_I2 = np.eye(2, dtype=complex)
_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = np.array([[1, 0], [0, -1]], dtype=complex)

_LETTER = {"I": _I2, "X": _SX, "Y": _SY, "Z": _SZ}


def _site_op(letter: str, k: int, N: int) -> np.ndarray:
    """N-qubit operator with single Pauli `letter` on site k, identity elsewhere."""
    ops = [_I2] * N
    ops[k] = _LETTER[letter]
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def _chain_bonds(N: int) -> list[tuple[int, int]]:
    return [(i, i + 1) for i in range(N - 1)]


def _build_bond_H(N: int, bonds: list[tuple[int, int]], terms: list[tuple[str, str, float]]
                  ) -> np.ndarray:
    """Build H = Σ_bond Σ_term coeff·σ_a^i σ_b^j over the given bonds."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        for (la, lb, coeff) in terms:
            H = H + coeff * (_site_op(la, i, N) @ _site_op(lb, j, N))
    return H


def _build_LH(H: np.ndarray) -> np.ndarray:
    """Hamiltonian superoperator L_H = −i(H ⊗ I − I ⊗ H.T) (column-stack vec)."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def _build_LD_nonuniform(N: int, gamma_l: list[float]) -> np.ndarray:
    """Z-dephasing dissipator (no Hamiltonian, no σ shift).

    L_D(ρ) = Σ_l γ_l · (Z_l ρ Z_l − ρ).  In vec form (column-stack),
       L_D = Σ_l γ_l · (Z_l ⊗ Z_l.conj − I).
    Matches the standard convention used in cross_term_formula_check.py.
    """
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d, dtype=complex)
    L_D = np.zeros((d2, d2), dtype=complex)
    for l, gamma in enumerate(gamma_l):
        if gamma == 0:
            continue
        Zl = _site_op("Z", l, N)
        ZdZ = Zl.conj().T @ Zl  # equals I since Z²=I, but written generally
        L_D += gamma * (
            np.kron(Zl, Zl.conj())
            - 0.5 * np.kron(ZdZ, Id)
            - 0.5 * np.kron(Id, ZdZ.T)
        )
    return L_D


def _frob_sq(M: np.ndarray) -> float:
    """Real Frobenius norm squared ‖M‖²_F = Re Tr(M† M)."""
    return float(np.real(np.trace(M.conj().T @ M)))


def _anticomm_norm_sq(A: np.ndarray, B: np.ndarray) -> float:
    """‖{A, B}‖²_F = ‖AB + BA‖²_F."""
    AC = A @ B + B @ A
    return _frob_sq(AC)


# --------------------------------------------------------------------------- #
# Phase 2: assertion tracking                                                 #
# --------------------------------------------------------------------------- #


_PHASE2_TOLERANCE = 1e-10
_PHASE2_RESULTS: list[dict] = []


def _assert_match(label: str, candidate: float, truth: float, *, tol: float = _PHASE2_TOLERANCE) -> None:
    """Phase 2 inline assertion: record + check |candidate − truth| < tol.

    Records (label, candidate, truth, gap) in _PHASE2_RESULTS for the final summary
    block; raises AssertionError immediately on mismatch.
    """
    gap = candidate - truth
    _PHASE2_RESULTS.append({
        "label": label,
        "candidate": candidate,
        "truth": truth,
        "gap": gap,
    })
    assert abs(gap) < tol, (
        f"Phase 2 assertion FAILED at {label}: candidate={candidate!r}, "
        f"truth={truth!r}, gap={gap!r}, tol={tol!r}."
    )


# --------------------------------------------------------------------------- #
# H-class builders                                                            #
# --------------------------------------------------------------------------- #


def _heisenberg_chain(N: int, J: float = 1.0) -> np.ndarray:
    return _build_bond_H(N, _chain_bonds(N), [
        ("X", "X", J), ("Y", "Y", J), ("Z", "Z", J),
    ])


def _ising_chain(N: int, J: float = 1.0) -> np.ndarray:
    return _build_bond_H(N, _chain_bonds(N), [("Z", "Z", J)])


def _xy_chain(N: int, J: float = 1.0) -> np.ndarray:
    return _build_bond_H(N, _chain_bonds(N), [("X", "X", J), ("Y", "Y", J)])


def _soft_xy_yx_chain(N: int, J: float = 1.0) -> np.ndarray:
    """Soft Π²-odd Hamiltonian XY+YX (mixed-letter bond pairs)."""
    return _build_bond_H(N, _chain_bonds(N), [("X", "Y", J), ("Y", "X", J)])


# --------------------------------------------------------------------------- #
# Per-bond ‖L_H^bond‖²                                                        #
# --------------------------------------------------------------------------- #


def _bond_LH_norm_sq(N: int, bond: tuple[int, int],
                     terms: list[tuple[str, str, float]]) -> float:
    """‖L_H^bond‖²_F for a single-bond Hamiltonian on N sites."""
    H_bond = _build_bond_H(N, [bond], terms)
    L_bond = _build_LH(H_bond)
    return _frob_sq(L_bond)


# --------------------------------------------------------------------------- #
# Section 1: multi-centering at N=3 Heisenberg γ=[0.1, 0.2, 0.3]              #
# --------------------------------------------------------------------------- #


def section_1_multi_centering() -> dict:
    print("=" * 78)
    print("SECTION 1: Multi-centering at N=3 Heisenberg γ=[0.1, 0.2, 0.3]")
    print("=" * 78)
    N = 3
    gamma_l = [0.1, 0.2, 0.3]
    sigma = sum(gamma_l)
    gamma_bar = sigma / N
    print(f"  N = {N}, γ = {gamma_l}, Σγ = {sigma}, γ̄ = {gamma_bar}")

    H = _heisenberg_chain(N)
    L_H = _build_LH(H)
    L_D = _build_LD_nonuniform(N, gamma_l)
    d2 = L_H.shape[0]
    Id = np.eye(d2, dtype=complex)

    centerings = {
        "(a) F1-centered  L_D + (Σγ)·I":              L_D + sigma * Id,
        "(b) Naive uniform L_D + N·γ̄·I":              L_D + (N * gamma_bar) * Id,
        "(c) Uncentered    L_D":                       L_D.copy(),
        "(d) Mean-field    L_D + (Σγ/N)·I":            L_D + (sigma / N) * Id,
    }

    target_truth = 163.84
    print()
    print(f"  Reference truth from prior exploratory test: {target_truth}")
    print(f"  Architect's candidate formula prediction:    204.8")
    print()

    results = {}
    matched_labels = []
    for label, L_Dc in centerings.items():
        val = _anticomm_norm_sq(L_H, L_Dc)
        norm_LH = np.sqrt(_frob_sq(L_H))
        norm_LDc = np.sqrt(_frob_sq(L_Dc))
        results[label] = {"value": val, "norm_LH": norm_LH, "norm_LDc": norm_LDc}
        match_str = ""
        if abs(val - target_truth) < 1e-6:
            match_str = "  <- MATCHES PREVIOUS TRUTH 163.84"
            matched_labels.append(label)
        print(f"  {label}")
        print(f"      ‖{{L_H, L_Dc}}‖² = {val:.6f}{match_str}")
        print(f"      ‖L_H‖ = {norm_LH:.6f},  ‖L_Dc‖ = {norm_LDc:.6f}")
        print()
    # Prefer (a) F1-centered as canonical label if it matches (a and b are equivalent).
    truth_match = None
    for preferred in matched_labels:
        if preferred.startswith("(a)"):
            truth_match = preferred
            break
    if truth_match is None and matched_labels:
        truth_match = matched_labels[0]

    # Note: (a) and (b) are algebraically IDENTICAL for any γ pattern, since
    #   N·γ̄ = N·(Σγ/N) = Σγ. The two labels describe the same shift.
    # (d) "per-site mean field Σγ_l/N" shifts by γ̄ alone (not N·γ̄), so (d)
    # differs from (a)/(b).
    print(f"  Note on centering relationships:")
    print(f"        (a) shifts by Σγ          = {sigma:.4f}")
    print(f"        (b) shifts by N·γ̄         = {N * gamma_bar:.4f}")
    print(f"        (a) ≡ (b) algebraically: N·γ̄ = N·(Σγ/N) = Σγ.")
    print(f"        (c) shifts by 0 (uncentered).")
    print(f"        (d) shifts by Σγ/N = γ̄  = {sigma / N:.4f} (different from (a)/(b)).")
    print()

    print("  Centering identification:")
    if truth_match is not None:
        print(f"    Truth 163.84 corresponds to centering: {truth_match}")
    else:
        print(f"    None of the four centerings give 163.84. Re-examine prior claim.")
        # Show closest
        best = min(results.items(), key=lambda kv: abs(kv[1]["value"] - target_truth))
        print(f"    Closest: {best[0]} = {best[1]['value']:.6f}  "
              f"(off by {best[1]['value'] - target_truth:+.6f})")
    print()

    return {"results": results, "truth_match": truth_match,
            "gamma_l": gamma_l, "sigma": sigma}


# --------------------------------------------------------------------------- #
# Section 2: Architect's candidate formula at the identified centering        #
# --------------------------------------------------------------------------- #


def section_2_architect_candidate(section_1_data: dict) -> dict:
    print("=" * 78)
    print("SECTION 2: Architect's candidate formula at N=3 Heisenberg γ=[0.1, 0.2, 0.3]")
    print("=" * 78)
    N = 3
    gamma_l = section_1_data["gamma_l"]
    bonds = _chain_bonds(N)
    heis_terms = [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)]

    # First verify the hand-computable ‖L_H^bond‖² for a Heisenberg bond at J=1.
    # The planning agent claims 480.
    bond01 = (0, 1)
    norm_bond_sq = _bond_LH_norm_sq(N, bond01, heis_terms)
    print(f"  ‖L_H^bond‖² for Heisenberg bond (0,1) at N=3, J=1: {norm_bond_sq:.6f}")
    print(f"     (architect claims 480.0)")
    architect_bond = 480.0
    print(f"     measured = {norm_bond_sq}, claim = {architect_bond}, "
          f"match = {abs(norm_bond_sq - architect_bond) < 1e-6}")
    print()

    # Architect's formula:
    #   candidate = 4·Σ_bond ‖L_H^bond‖²·Σ_{m∉bond} γ_m²    (spectator part)
    #             + Σ_bond G(bond, H)·(γ_i − γ_j)²          (asymmetry part)
    # with G(bond, Heisenberg) = (4/3)·‖L_H^bond‖².

    spectator_part = 0.0
    asymmetry_part = 0.0
    bond_diag = []
    for (i, j) in bonds:
        nrm = _bond_LH_norm_sq(N, (i, j), heis_terms)
        spectators = [m for m in range(N) if m not in (i, j)]
        sum_gm2 = sum(gamma_l[m] ** 2 for m in spectators)
        spect_contrib = 4.0 * nrm * sum_gm2
        G_bond = (4.0 / 3.0) * nrm  # Heisenberg
        asym_contrib = G_bond * (gamma_l[i] - gamma_l[j]) ** 2
        spectator_part += spect_contrib
        asymmetry_part += asym_contrib
        bond_diag.append({
            "bond": (i, j),
            "norm_LH_bond_sq": nrm,
            "spectators": spectators,
            "sum_gm2": sum_gm2,
            "spect_contrib": spect_contrib,
            "G_bond_heisenberg": G_bond,
            "delta_gamma_sq": (gamma_l[i] - gamma_l[j]) ** 2,
            "asym_contrib": asym_contrib,
        })

    candidate = spectator_part + asymmetry_part
    print(f"  Per-bond breakdown:")
    for d in bond_diag:
        print(f"    bond {d['bond']}: ‖L_H^bond‖² = {d['norm_LH_bond_sq']:.4f}, "
              f"spectators = {d['spectators']} → Σγ_m² = {d['sum_gm2']:.6f}")
        print(f"        spect contrib = 4·‖L_H^bond‖²·Σγ_m² = {d['spect_contrib']:.4f}")
        print(f"        G(bond) = (4/3)·{d['norm_LH_bond_sq']:.4f} = {d['G_bond_heisenberg']:.4f}, "
              f"(γ_i−γ_j)² = {d['delta_gamma_sq']:.6f}")
        print(f"        asym contrib = {d['asym_contrib']:.4f}")
    print()
    print(f"  Spectator part total:  {spectator_part:.6f}")
    print(f"  Asymmetry part total:  {asymmetry_part:.6f}")
    print(f"  Candidate formula:     {candidate:.6f}")
    print(f"  Architect claim (planning-agent):  204.8")
    print()

    truth_label = section_1_data["truth_match"]
    if truth_label is not None:
        truth = section_1_data["results"][truth_label]["value"]
        print(f"  Truth at the identified centering ({truth_label.split()[0]}): {truth:.6f}")
        print(f"  Gap (candidate − truth) = {candidate - truth:+.6f}")
        # Phase 2 assertion: N=3 Heisenberg γ=[0.1, 0.2, 0.3] anchor pinned bit-exact.
        _assert_match("N=3 Heisenberg γ=[0.1, 0.2, 0.3]", candidate, truth)
    else:
        # Fall back to F1-centered for the architect-derivation comparison.
        print(f"  WARN: section 1 did not identify a centering equal to 163.84.")
        truth_a = section_1_data["results"]["(a) F1-centered  L_D + (Σγ)·I"]["value"]
        print(f"  Truth at F1-centered (architect's assumed centering): {truth_a:.6f}")
        print(f"  Gap (candidate − F1 truth) = {candidate - truth_a:+.6f}")
        # Phase 2 assertion: F1-centered is the architect's intended centering;
        # any future regression here surfaces immediately.
        _assert_match("N=3 Heisenberg γ=[0.1, 0.2, 0.3] (F1-centered fallback)", candidate, truth_a)
    print()

    return {
        "candidate": candidate,
        "spectator_part": spectator_part,
        "asymmetry_part": asymmetry_part,
        "bond_diag": bond_diag,
    }


# --------------------------------------------------------------------------- #
# Section 3: scan N=3, 4, 5 Heisenberg chain with γ_l = 0.05·(l+1)            #
# --------------------------------------------------------------------------- #


def _scan_n_for_h(N_list: list[int], H_builder, h_class_terms,
                   class_label: str, G_factor) -> list[dict]:
    """Run truth-vs-candidate scan for a given Hamiltonian class.

    Args:
        N_list: N values to scan.
        H_builder(N) -> matrix.
        h_class_terms: bond term list used to compute ‖L_H^bond‖² and bond enumeration.
        class_label: name for the print header.
        G_factor: coefficient on ‖L_H^bond‖² for G(bond, H) (architect's claim).
    """
    out = []
    for N in N_list:
        gamma_l = [0.05 * (l + 1) for l in range(N)]
        sigma = sum(gamma_l)
        d2 = 4 ** N
        Id = np.eye(d2, dtype=complex)
        H = H_builder(N)
        L_H = _build_LH(H)
        L_D = _build_LD_nonuniform(N, gamma_l)
        # F1-centered L_Dc (per architect's assumption)
        L_Dc_F1 = L_D + sigma * Id
        truth_F1 = _anticomm_norm_sq(L_H, L_Dc_F1)
        # Also compute uncentered for reference.
        truth_unc = _anticomm_norm_sq(L_H, L_D)

        # Architect formula
        spectator_part = 0.0
        asymmetry_part = 0.0
        bonds = _chain_bonds(N)
        for (i, j) in bonds:
            nrm = _bond_LH_norm_sq(N, (i, j), h_class_terms)
            spectators = [m for m in range(N) if m not in (i, j)]
            sum_gm2 = sum(gamma_l[m] ** 2 for m in spectators)
            spectator_part += 4.0 * nrm * sum_gm2
            asymmetry_part += G_factor * nrm * (gamma_l[i] - gamma_l[j]) ** 2
        candidate = spectator_part + asymmetry_part

        gap_F1 = candidate - truth_F1
        rel_gap_F1 = gap_F1 / truth_F1 if truth_F1 > 1e-12 else float("nan")
        out.append({
            "N": N,
            "gamma_l": gamma_l,
            "truth_F1": truth_F1,
            "truth_unc": truth_unc,
            "candidate": candidate,
            "spectator_part": spectator_part,
            "asymmetry_part": asymmetry_part,
            "gap_F1": gap_F1,
            "rel_gap_F1": rel_gap_F1,
        })
        print(f"  {class_label}  N={N}  γ={[f'{g:.2f}' for g in gamma_l]}")
        print(f"      truth (F1-centered):  {truth_F1:.6f}")
        print(f"      truth (uncentered):   {truth_unc:.6f}")
        print(f"      architect candidate:  {candidate:.6f}")
        print(f"          spectator: {spectator_part:.6f}, asymmetry: {asymmetry_part:.6f}")
        print(f"      gap = candidate − truth(F1) = {gap_F1:+.6f}  "
              f"({rel_gap_F1 * 100:+.2f}%)")
        print()
        # Phase 2 assertion: each (N, γ) row in the N-scan is pinned bit-exact.
        _assert_match(f"{class_label} N={N} γ_l=0.05·(l+1)", candidate, truth_F1)
    return out


def section_3_scan_n() -> list[dict]:
    print("=" * 78)
    print("SECTION 3: N-scan, Heisenberg chain, γ_l = 0.05·(l+1)")
    print("=" * 78)
    print("  Tests hypothesis (i) [N=3-only] vs (ii) [formula gap at all N].")
    print()
    heis_terms = [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)]
    return _scan_n_for_h([3, 4, 5], _heisenberg_chain, heis_terms,
                         "Heisenberg", G_factor=(4.0 / 3.0))


# --------------------------------------------------------------------------- #
# Section 4: cross-H-class sanity at N=4                                      #
# --------------------------------------------------------------------------- #


def section_4_cross_h_class(scan_results: list[dict]) -> dict:
    print("=" * 78)
    print("SECTION 4: Cross-H-class sanity at N=4, γ=[0.05, 0.10, 0.15, 0.20]")
    print("=" * 78)
    N = 4
    gamma_l = [0.05, 0.10, 0.15, 0.20]
    sigma = sum(gamma_l)
    d2 = 4 ** N
    Id = np.eye(d2, dtype=complex)
    bonds = _chain_bonds(N)

    classes = [
        ("Heisenberg XX+YY+ZZ", _heisenberg_chain(N),
            [("X", "X", 1.0), ("Y", "Y", 1.0), ("Z", "Z", 1.0)], 4.0 / 3.0),
        ("Ising ZZ-only",       _ising_chain(N),
            [("Z", "Z", 1.0)],                                    4.0),
        ("XY (XX+YY)",          _xy_chain(N),
            [("X", "X", 1.0), ("Y", "Y", 1.0)],                   0.0),
        ("Soft XY+YX",          _soft_xy_yx_chain(N),
            [("X", "Y", 1.0), ("Y", "X", 1.0)],                   0.0),
    ]

    out = {}
    for label, H, h_terms, G_factor in classes:
        L_H = _build_LH(H)
        L_D = _build_LD_nonuniform(N, gamma_l)
        L_Dc = L_D + sigma * Id
        truth = _anticomm_norm_sq(L_H, L_Dc)

        spectator_part = 0.0
        asymmetry_part = 0.0
        for (i, j) in bonds:
            nrm = _bond_LH_norm_sq(N, (i, j), h_terms)
            spectators = [m for m in range(N) if m not in (i, j)]
            sum_gm2 = sum(gamma_l[m] ** 2 for m in spectators)
            spectator_part += 4.0 * nrm * sum_gm2
            asymmetry_part += G_factor * nrm * (gamma_l[i] - gamma_l[j]) ** 2
        candidate = spectator_part + asymmetry_part
        gap = candidate - truth
        rel_gap = gap / truth if truth > 1e-12 else float("nan")
        match = abs(gap) < 1e-6 or (abs(rel_gap) < 1e-8)
        print(f"  {label}")
        print(f"      truth:               {truth:.6f}")
        print(f"      architect candidate: {candidate:.6f}")
        print(f"          spectator: {spectator_part:.6f}, asymmetry: {asymmetry_part:.6f}")
        print(f"          G factor on ‖L_H^bond‖²: {G_factor}")
        print(f"      gap = {gap:+.6f}  ({rel_gap * 100:+.4f}%)  "
              f"{'MATCH' if match else 'MISMATCH'}")
        print()
        out[label] = {
            "truth": truth, "candidate": candidate, "gap": gap, "rel_gap": rel_gap,
            "match": match,
        }
        # Phase 2 assertion: each H-class row at N=4 γ=[0.05, 0.10, 0.15, 0.20] is pinned bit-exact.
        _assert_match(f"{label} N=4 γ=[0.05, 0.10, 0.15, 0.20]", candidate, truth)
    return out


# --------------------------------------------------------------------------- #
# Section 5: if hypothesis (iii), report corrected formula under the          #
# identified centering                                                        #
# --------------------------------------------------------------------------- #


def section_5_centering_correction(section_1_data: dict) -> dict:
    print("=" * 78)
    print("SECTION 5: Centering-correction analysis (if hypothesis (iii) applies)")
    print("=" * 78)
    truth_label = section_1_data["truth_match"]
    if truth_label is None:
        print("  Section 1 did not identify any centering equal to 163.84.")
        print("  Either the prior 163.84 claim was incorrect, or another centering")
        print("  variant (not in the four canonical options) was used.")
        print()
        print("  Differential analysis: how do centerings relate?")
        print("    Let L_Dc(s) = L_D + s·I for any real shift s.")
        print("    Then {L_H, L_Dc(s)} = {L_H, L_D} + 2s·L_H   (since {L_H, I} = 2L_H).")
        print("    So ‖{L_H, L_Dc(s)}‖² = ‖{L_H, L_D}‖² + 4s·Re⟨{L_H, L_D}, L_H⟩")
        print("                          + 4s²·‖L_H‖²")
        print()
        return {"applies": False, "reason": "no_centering_matches_163.84"}
    if truth_label.startswith("(a)") or truth_label.startswith("(b)"):
        # (a) and (b) are algebraically identical (N·γ̄ = Σγ); both ARE the F1-centered case.
        print("  Truth 163.84 corresponds to F1-centered L_Dc (centering (a) ≡ (b),")
        print("  since N·γ̄ = N·(Σγ/N) = Σγ for any γ pattern).")
        print("  This is the SAME centering the architect assumed.")
        print("  → Hypothesis (iii) does NOT apply: no centering mismatch.")
        return {"applies": False, "truth_centering": "(a)=(b) F1-centered"}
    print(f"  Truth 163.84 corresponds to {truth_label}, NOT the F1-centered")
    print(f"  centering the architect assumed.")
    print(f"  Hypothesis (iii) IS the explanation.")
    print()
    print(f"  The architect's formula was derived for centering (a). Under the")
    print(f"  identified centering, the formula needs the shift-correction term.")
    print(f"  Recall: ‖{{L_H, L_D + s·I}}‖² = ‖{{L_H, L_D}}‖² + 4s·X + 4s²·‖L_H‖²")
    print(f"  where X = Re⟨{{L_H, L_D}}, L_H⟩.")
    return {"applies": True, "truth_centering": truth_label}


# --------------------------------------------------------------------------- #
# Conclusion                                                                  #
# --------------------------------------------------------------------------- #


def conclude(section_1_data, section_2_data, scan_results, section_4_data,
             section_5_data) -> None:
    print("=" * 78)
    print("CONCLUSION")
    print("=" * 78)

    truth_label = section_1_data["truth_match"]

    # Sanity-check the planning agent's hand-computation: did the claimed
    # ‖L_H^bond‖² = 480 match the measurement?
    bond_norm_sq = section_2_data["bond_diag"][0]["norm_LH_bond_sq"]
    bond_claim = 480.0
    bond_anchor_correct = abs(bond_norm_sq - bond_claim) < 1e-6
    print(f"  Planning-agent hand-anchor cross-check:")
    print(f"    claim:    ‖L_H^bond‖² = {bond_claim} for Heisenberg bond at J=1")
    print(f"    measured: ‖L_H^bond‖² = {bond_norm_sq}")
    print(f"    {'CONSISTENT' if bond_anchor_correct else 'WRONG: hand-anchor off by factor ' + f'{bond_claim / bond_norm_sq:.4f}'}")
    print()

    # If (iii) — centering mismatch — that's the answer.
    if section_5_data.get("applies"):
        print(f"  Hypothesis (iii) CONFIRMED: centering mismatch.")
        print(f"  Truth 163.84 corresponds to {truth_label}, NOT")
        print(f"  the architect's assumed F1-centered L_D + Σγ·I.")
        print(f"  → Replan with the corrected centering as anchor.")
        return

    # Examine N-scan to distinguish (i) vs (ii) vs no-gap.
    if not scan_results:
        print("  No scan data — cannot distinguish (i) vs (ii).")
        return

    n3_row = next((r for r in scan_results if r["N"] == 3), None)
    n4_row = next((r for r in scan_results if r["N"] == 4), None)
    n5_row = next((r for r in scan_results if r["N"] == 5), None)
    rel_n3 = n3_row["rel_gap_F1"] if n3_row else None
    rel_n4 = n4_row["rel_gap_F1"] if n4_row else None
    rel_n5 = n5_row["rel_gap_F1"] if n5_row else None
    print(f"  N-scan gaps (Heisenberg chain, F1-centered):")
    print(f"    N=3: gap = {n3_row['gap_F1']:+.6f}, relative {rel_n3 * 100:+.4f}%")
    print(f"    N=4: gap = {n4_row['gap_F1']:+.6f}, relative {rel_n4 * 100:+.4f}%")
    print(f"    N=5: gap = {n5_row['gap_F1']:+.6f}, relative {rel_n5 * 100:+.4f}%")
    print()

    # Anchor candidate (N=3 γ=[0.1,0.2,0.3]) match
    arch_cand = section_2_data["candidate"]
    if truth_label is not None:
        anchor_truth = section_1_data["results"][truth_label]["value"]
    else:
        anchor_truth = float("nan")
    anchor_match = abs(arch_cand - anchor_truth) < 1e-6
    print(f"  Anchor N=3 γ=[0.1,0.2,0.3] (Heisenberg) check:")
    print(f"    truth       = {anchor_truth}")
    print(f"    candidate   = {arch_cand}")
    print(f"    {'MATCH' if anchor_match else 'MISMATCH'}")
    print()

    # Determine hypothesis (i) vs (ii) vs no-gap by gap shrinkage.
    n3_small = abs(rel_n3) < 1e-6
    n4_small = abs(rel_n4) < 1e-6
    n5_small = abs(rel_n5) < 1e-6
    cross_class_all_match = section_4_data and all(d["match"] for d in section_4_data.values())

    if n3_small and n4_small and n5_small and anchor_match and cross_class_all_match:
        print(f"  NO HYPOTHESIS (i/ii/iii) APPLIES: architect's formula is CORRECT")
        print(f"  at the F1-centered L_Dc, at all tested N ∈ {{3, 4, 5}}, and across")
        print(f"  all four tested H-classes (Heisenberg, Ising, XY, XY+YX).")
        print()
        if not bond_anchor_correct:
            print(f"  The planning-agent's reported 'candidate = 204.8' was a")
            print(f"  HAND-CALCULATION ERROR: they used ‖L_H^bond‖² = {bond_claim:.0f}")
            print(f"  but the correct value at N=3, Heisenberg, J=1 is {bond_norm_sq:.0f}.")
            print(f"  With the correct value, candidate = {arch_cand:.4f}, matching the")
            print(f"  truth {anchor_truth:.4f} exactly.")
            print()
            print(f"  → Phase 2: proceed directly to proof + typed claim; the closed")
            print(f"    form is validated by this script for the four tested classes.")
            print(f"    Recommend including a fixed hand-anchor (‖L_H^bond‖² = 384)")
            print(f"    in any future planning document referencing the N=3 anchor.")
        else:
            print(f"  → Phase 2: proceed to proof + typed claim.")
        return

    if n4_small and n5_small and not n3_small:
        print(f"  Hypothesis (i) CONFIRMED: N=3-only artifact.")
        print(f"  Architect formula matches truth at N=4 and N=5 but not at N=3.")
        print(f"  The N=3 chain has only 2 bonds (0,1) and (1,2) which OVERLAP at")
        print(f"  site 1 (the middle qubit); this violates the 'disjoint bond")
        print(f"  supports' assumption in PROOF_CROSS_TERM_FORMULA Lemma 3 Corollary.")
        print(f"  → Phase 2: proof + typed claim with explicit 'N >= 4' precondition.")
        return
    if not n3_small and not n4_small and not n5_small:
        print(f"  Hypothesis (ii) CONFIRMED: general formula gap.")
        print(f"  Gap persists at all N ∈ {{3, 4, 5}}.")
        if section_4_data:
            print(f"  Cross-H-class N=4 sanity:")
            for lbl, d in section_4_data.items():
                tag = "MATCH" if d["match"] else "MISMATCH"
                print(f"    {lbl}: gap {d['gap']:+.6f} ({d['rel_gap'] * 100:+.4f}%)  {tag}")
        print(f"  → Phase 2: architect needs another round; surface specific class mismatches.")
        return

    print(f"  Mixed result — manually inspect the N-scan and class sanity tables above.")


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #


def section_6_phase2_assertion_summary() -> None:
    """Phase 2: roll up every inline _assert_match into a single block.

    Phase 1's exploratory print structure is preserved above; this block confirms
    that every (N, H, γ) anchor visited along the way pinned bit-exact to the
    closed form within tolerance _PHASE2_TOLERANCE.
    """
    print("=" * 78)
    print("SECTION 6 (Phase 2): closed-form assertion summary")
    print("=" * 78)
    if not _PHASE2_RESULTS:
        print("  WARN: no Phase 2 assertions registered; check that section runners are wired.")
        return
    print(f"  tolerance: |candidate − truth| < {_PHASE2_TOLERANCE:.0e}")
    print(f"  registered checks: {len(_PHASE2_RESULTS)}")
    print()
    print("  pinned anchors:")
    for row in _PHASE2_RESULTS:
        print(f"    {row['label']}: candidate={row['candidate']:.6f}, "
              f"truth={row['truth']:.6f}, gap={row['gap']:+.2e}")
    # All inline _assert_match calls have already raised AssertionError on any gap;
    # if we get here, the run is bit-exact across all registered anchors.
    print()
    print("  All N=3,4,5 × 4-H-classes Phase 2 closed-form formula verified bit-exact.")
    print()


def main() -> int:
    print("F49 NON-UNIFORM γ CROSS-TERM — PHASE 1 VERIFICATION + PHASE 2 ASSERTIONS")
    print("=" * 78)
    print("Phase 1: distinguishes hypotheses (i) N=3-only, (ii) general gap, (iii) centering.")
    print("Phase 2: asserts |candidate − truth| < 1e-10 at each (N, H, γ) anchor and rolls up")
    print("         into a final 'All ... verified bit-exact' summary; AssertionError on regression.")
    print()

    section_1_data = section_1_multi_centering()
    section_2_data = section_2_architect_candidate(section_1_data)
    scan_results = section_3_scan_n()
    section_4_data = section_4_cross_h_class(scan_results)
    section_5_data = section_5_centering_correction(section_1_data)

    conclude(section_1_data, section_2_data, scan_results, section_4_data, section_5_data)

    # Phase 2 final assertion block (every inline _assert_match has already raised on mismatch).
    section_6_phase2_assertion_summary()
    return 0


if __name__ == "__main__":
    sys.exit(main())
