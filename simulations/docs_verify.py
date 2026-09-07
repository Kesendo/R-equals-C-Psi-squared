"""
SCOPED ROUTE-B CURRENT-TRUTH AND NUMERICAL DOCUMENTATION GATE
==========================================================

This is a scoped Route-B/current-truth cross-document consistency gate plus the
explicitly implemented numerical smoke checks. It is not a repo-wide completeness proof.

Authors: Tom Wicht, Claude
Date: March 16, 2026
"""

import numpy as np
import json
import sys
import re
import subprocess
from itertools import product as iprod
from pathlib import Path
from framework import max_f1_pairing_distance

I2 = np.eye(2)
sx = np.array([[0,1],[1,0]])
sy = np.array([[0,-1j],[1j,0]])
sz = np.array([[1,0],[0,-1]])

def tensor(*ops):
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def site_op(op, s, N):
    ops = [I2]*N; ops[s] = op
    return tensor(*ops)

def build_H(N, J, topo, J_SA=None, J_SB=None):
    d = 2**N; H = np.zeros((d,d), dtype=complex)
    if topo == "chain":
        pairs = [(i,i+1,J) for i in range(N-1)]
    elif topo == "ring":
        pairs = [(i,(i+1)%N,J) for i in range(N)]
    elif topo == "star":
        J_SA = J_SA or J; J_SB = J_SB or J
        pairs = [(0,1,J_SA)] + [(0,i,J_SB) for i in range(2,N)]
    else: raise ValueError(topo)
    for i,j,Jij in pairs:
        for p in [sx,sy,sz]:
            H += Jij * site_op(p,i,N) @ site_op(p,j,N)
    return H

def build_L(H, gammas, N):
    d = 2**N; d2 = d*d
    L = -1j*(np.kron(H,np.eye(d)) - np.kron(np.eye(d),H.T))
    for k in range(N):
        Zk = site_op(sz,k,N)
        L += gammas[k]*(np.kron(Zk,Zk.conj()) - np.eye(d2))
    return L

def make_ghz(N):
    d=2**N; psi=np.zeros(d,complex); psi[0]=psi[-1]=1/np.sqrt(2)
    return np.outer(psi, psi.conj())

def make_w(N):
    d=2**N; psi=np.zeros(d,complex)
    for i in range(N): psi[1<<(N-1-i)] = 1/np.sqrt(N)
    return np.outer(psi, psi.conj())

def concurrence_2q(rho):
    """Wootters concurrence for 2-qubit state."""
    sy2 = np.kron(sy, sy)
    rho_tilde = sy2 @ rho.conj() @ sy2
    R = rho @ rho_tilde
    eigs = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, eigs[0] - eigs[1] - eigs[2] - eigs[3])

def partial_trace_keep01(rho, N):
    """Trace out all but first 2 qubits."""
    d = 2**N
    rho_r = rho.reshape([2]*2*N)
    for i in range(N-1, 1, -1):
        rho_r = np.trace(rho_r, axis1=i, axis2=i+N-(N-1-i+1)+(N-1-i))
    # Simpler: reshape and trace manually
    rho_r = rho.reshape(4, 2**(N-2), 4, 2**(N-2))
    return np.trace(rho_r, axis1=1, axis2=3)

PASS = 0; FAIL = 0
SUCCESS_SCOPE = ("This is a scoped Route-B/current-truth cross-document consistency gate with numerical smoke checks; "
                 "not a repo-wide completeness proof.")

def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1; print(f"  PASS: {name}")
    else:
        FAIL += 1; print(f"  FAIL: {name} -- {detail}")

def verify_docs_gate_scope_description():
    def matches_scope(text):
        normalized = " ".join(text.split())
        return ("scoped Route-B/current-truth cross-document consistency gate" in normalized
                and "not a repo-wide completeness proof" in normalized)

    for label, text in (("module description", __doc__), ("success footer", SUCCESS_SCOPE)):
        check(f"Docs verifier accurate scope: {label}", matches_scope(text))
        check(f"Docs verifier rejects missing completeness boundary: {label}",
              not matches_scope(text.replace("not a repo-wide completeness proof", "")))


def route_b_current_arc_slice(source):
    """The current NextStep before the explicit journal delimiter; history stays history."""
    current = source[source.index('NextStep: "CURRENT STATE'):]
    current, delimiter, _journal = current.partition("--- THE JOURNAL (how we got here) ---")
    if not delimiter:
        raise ValueError("OpenArcs current-truth slice requires its journal delimiter")
    return " ".join(current.replace('" +', "").replace('"', "").split())


def verify_route_b_mechanism_current_truth():
    """Prose scope gate, not a proof of residual integrability or of its all-N extension."""
    root = Path(__file__).resolve().parents[1]
    forbidden = (
        "ε_k = 2J cos(kπ/(N+1))", "Tier 1 for the integrability mechanism",
        "Tier 1 for the mechanism, Tier 2 for completeness",
        "integrability keeps the diabolic crossings at every N",
        "integrability keeps the crossings", "diabolics survive at every N",
        "it auto-satisfies the node conditions", "at all N, throughout the complex-q plane",
        "integrability governs existence", "integrability confirmed",
        "That is now verified, not just read", "integrability mechanism confirmed",
        "The Δ-test confirms integrability-protection",
        "integrability character of these real-q diabolics, by contrast, is now settled",
        "protected by the same free-fermion structure as every other diabolic",
        "residual SET is Δ-stable because the ZZ term is Hermitian",
        "the diabolic was integrability-protected",
        "Confirms integrability-protection generalizes",
        "ALSO integrability-protected", "the integrability signature",
        "diabolics themselves DIE under Δ (integrability-protected",
        "Move 2, integrability-protected", "refutes integrability-protection",
        "defects/lifts => integrability-protected",
        "additivity broken => an integrable diabolic dies",
        "kills the diabolics themselves", "the diabolics themselves die under Δ",
        "makes the diabolics defective", "an integrable diabolic flips",
        "free-fermion integrability's", "silence is integrability's",
        "an integrability fact", "a free-fermion integrability fact",
        "same integrability reason", "silence = free-fermion's",
        "the diabolic's silence is integrability's",
        "the R-1/R-2 generically-EXPECTED outcome at N>=5",
        "the N=4 self-fold that auto-",
        "diabolics defect instantly", "the recorded phenomenology follows at once",
        "rests on this hypothesis", "load-bearing there, inside a Tier-1 result",
        "the instant Δ ≠ 0", "generic XXZ gives the ordinary defective EP",
        "the Hamiltonian and the diagonal dissipator do not mix the occupation-number sectors",
        "the palindrome does fix the diabolic's position",
        "Jordan block (a √-branch)", "two sheets crossing ⟹ diabolic",
        "p≈1 diabolic / p≈½ defective",
        "1 = diabolic crossing", "EpCharacter and gap-exponent (≈ 1) read diabolic",
        "Residual semisimplicity follows from free-fermion integrability at every N",
        "exponent that classifies a coalescence",
        "is measured by the Δ-test (does an XXZ anisotropy flip them defective, as it does the complex-q ones?) is run",
    )
    surfaces = {
        "experiments/F89_PATH_K_DIABOLIC.md": (
            "ε_k = 4J cos(kπ/(N+1))", "AT-locked Hamiltonian eigenoperators",
            "Residual semisimplicity is conditional on the twin-scalar restriction",
            "Tier 2 for the residual integrability mechanism",
            "finite-N Δ controls do not prove all-N existence",
            "Last refreshed 2026-09-07",
            "analytic-defective collisions can split linearly",
            "analytic-split candidate", "full-sector geometric multiplicity decides character",
        ),
        "experiments/README.md": (
            "Tier 2 for the residual integrability mechanism",
            "conditional twin-scalar restriction", "finite-N Δ controls",
        ),
        "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs": (
            "finite-N Delta response", "conditional residual mechanism",
            "does not alone prove causality or all-N protection",
            "residual labels are chosen exactly at the base (q0=2, Delta=0)",
            "full-spectrum strands are tracked by nearest-neighbour continuation",
            "finite-path labeling, not a proof of a globally Delta-invariant AT/residual split",
            "Full-block character remains required",
            "proposal tracking can exchange identity near dense or degenerate encounters",
            "diagonal term -i*q*(fieldEnergy(ket)-fieldEnergy(bra))",
            "for complex q this need not be purely imaginary",
            "unchanged dissipator entries do not imply unchanged eigenmode real parts",
            "w=null reproduces BuildFull(n,q,Δ)",
        ),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzDeltaFlipTests.cs": (
            "finite-N Delta response", "conditional residual mechanism",
            "Delta-response signature",
        ),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzCoherenceBlockFieldTests.cs": (
            "These tests use real q=2",
            "diagonal addition -i*q*(fieldEnergy(ket)-fieldEnergy(bra))",
            "not eigenmode real parts", "matrix diagonal's imaginary component",
        ),
        "compute/RCPsiSquared.Diagnostics/Foundation/CrossFoldSimilarityWitness.cs": (
            "sampled finite-N Delta response", "does not establish all-N protection",
        ),
        "compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs": (
            "sampled finite-N Delta response", "does not establish all-N protection",
        ),
        "compute/RCPsiSquared.Cli/Commands/PathKMonodromyScanCommand.cs": (
            "finite-N Delta response", "conditional residual mechanism",
            "falsifies this sampled-locus Delta-death prediction",
            "not proof of causality or all-N protection",
            "bounded scan null", "not a replacement for the exact Route B inventory",
        ),
        "compute/RCPsiSquared.Core.Tests/F89PathK/WeightCoherenceBlockTests.cs": (
            "finite-N Delta response", "conditional residual mechanism",
        ),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/CrossFoldSimilarityWitnessTests.cs": (
            "finite-N Delta response", "conditional residual mechanism",
        ),
        "compute/RCPsiSquared.Cli/Program.cs": (
            "finite-N Delta response", "conditional residual mechanism",
        ),
        "compute/RCPsiSquared.Diagnostics/Foundation/BranchLocusPalindromeWitness.cs": (
            "twin-scalar restriction at the N=4 point", "not an all-N cause", "γ=1, J=q_EP",
        ),
        "compute/RCPsiSquared.Core/Symmetry/F89BranchLocusPalindromeClaim.cs": (
            "twin-scalar restriction at the N=4 point", "not an all-N cause",
        ),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": (
            "twin-scalar restriction at the N=4 point", "not an all-N cause",
        ),
        "docs/ANALYTICAL_FORMULAS.md": (
            "sampled finite-N Delta responses", "tested N=4 control",
            "general protection principle remains interpretive and conditional",
            "Loss of the W/H intertwining identity is an algebraic statement",
            "Finite-N defect-or-lift is a separate observed response",
        ),
        "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md": (
            "Twin-scalar lemma (conditional)",
            "Loss of the W/H intertwining identity is an algebraic statement",
            "Finite-N defect-or-lift is a separate observed response",
        ),
        "experiments/F89_BRANCH_LOCUS_PALINDROME.md": (
            "N=4 twin-scalar restriction", "sampled N=4 Delta control",
            "palindrome proof is independent of the Tier-2 hypothesis",
        ),
        "hypotheses/DIABOLIC_BY_INTEGRABILITY.md": (
            "palindrome proof is independent", "sampled N=4 Delta control",
            "mirror fixes the line and pairing", "N=4 overlap/twin-scalar restriction supplies the on-line position",
        ),
        "compute/RCPsiSquared.Diagnostics/Foundation/PathKMonodromyScout.cs": (
            "analytic-defective collisions can split linearly",
            "analytic-split candidate", "full-sector geometric multiplicity decides character",
        ),
        "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs": (
            "The finite-N Delta-flip controls support that residual mechanism but do not prove existence at every N",
            "only the H-half of the twin-scalar condition, while the D-half remains a locus-specific obligation",
            "MECHANISM BOUNDARY",
            "Free-fermion (Jordan-Wigner) integrability alone does not prove residual semisimplicity or a Delta response at every N",
            "exponent that reads the local eigenvalue split (~1 for a linear analytic split, ~0.5 for a sqrt branch)",
            "the exponent corroborates branch type but does not decide geometric multiplicity or semisimplicity",
            "held fixed as a proposal space",
            "Hermiticity leaves the dissipative rate diagonal untouched but does not by itself prove that complement invariant",
            "Character is therefore read on the full block",
        ),
    }

    def defects(source, required):
        normalized = " ".join(source.replace("**", "").replace("///", "").split()).casefold()
        scoped_forbidden = ()
        if required == surfaces["compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs"]:
            scoped_forbidden = ("integrability-protected", "delta-stable because")
        if required == surfaces["compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs"]:
            scoped_forbidden = ("split is delta-stable", "hermiticity proves the at complement invariant",
                                "the zz term is hermitian, so", "hermitian-implies-invariance",
                                "so no at capture", "at strands are in the full spectrum but never adopted",
                                "only shifts the imaginary frequency", "random w breaks integrability")
        if required == surfaces["compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzCoherenceBlockFieldTests.cs"]:
            scoped_forbidden = ("only shifts the imaginary frequency", "random field breaks integrability",
                                "absorption-theorem real rate untouched", "at real rate")
        delta_normalized = normalized.replace("δ", "delta")
        return ([f"missing {phrase}" for phrase in required if phrase.casefold() not in normalized]
                + [f"forbidden {phrase}" for phrase in forbidden if phrase.casefold() in normalized]
                + [f"forbidden {phrase}" for phrase in scoped_forbidden if phrase in delta_normalized])

    arc_path = "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs"
    arc_source = (root / arc_path).read_text(encoding="utf-8")
    wrong_arc = arc_source.replace("but do not prove existence at every N", "and prove existence at every N")
    check("Route B mechanism current OpenArcs mutation baseline",
          not defects(route_b_current_arc_slice(arc_source), surfaces.get(arc_path, ())),
          ascii("; ".join(defects(route_b_current_arc_slice(arc_source), surfaces.get(arc_path, ())))))
    check("Route B mechanism rejects current OpenArcs all-N promotion",
          bool(set(defects(route_b_current_arc_slice(wrong_arc), surfaces.get(arc_path, ())))
               - set(defects(route_b_current_arc_slice(arc_source), surfaces.get(arc_path, ())))))
    current_arc = route_b_current_arc_slice(arc_source)
    xxz_path = "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs"
    xxz_source = (root / xxz_path).read_text(encoding="utf-8")
    for mutation_name, original, required, mutated in (
        ("current OpenArcs causal protection", current_arc, surfaces[arc_path],
         current_arc + " The residual split is integrability-protected."),
        ("current OpenArcs removed proposal-space caveat", current_arc, surfaces[arc_path],
         current_arc.replace("held fixed as a proposal space", "globally invariant")),
        ("current OpenArcs removed full-block character caveat", current_arc, surfaces[arc_path],
         current_arc.replace("Character is therefore read on the full block", "The compressed block decides character")),
        ("Xxz Delta-stable split", xxz_source, surfaces[xxz_path],
         xxz_source + " The split is Delta-stable."),
        ("Xxz Hermitian-implies-invariance", xxz_source, surfaces[xxz_path],
         xxz_source + " Hermiticity proves the AT complement invariant."),
        ("Xxz removed tracking limitation", xxz_source, surfaces[xxz_path],
         xxz_source.replace("finite-path labeling, not a proof of a globally Delta-invariant AT/residual split", "globally invariant split")),
        ("Xxz removed full-block requirement", xxz_source, surfaces[xxz_path],
         xxz_source.replace("Full-block character remains required", "Compressed character is sufficient")),
        ("Xxz guaranteed no AT capture", xxz_source, surfaces[xxz_path],
         xxz_source + " AT strands are in the full spectrum but never adopted, so no AT capture."),
        ("Xxz removed identity-exchange caveat", xxz_source, surfaces[xxz_path],
         xxz_source.replace("proposal tracking can exchange identity near dense or degenerate encounters", "identities cannot exchange")),
        ("Xxz imaginary-only field shift", xxz_source, surfaces[xxz_path],
         xxz_source + " The field only shifts the imaginary frequency."),
        ("Xxz random-field integrability claim", xxz_source, surfaces[xxz_path],
         xxz_source + " A random w breaks integrability."),
        ("Xxz removed complex-q field boundary", xxz_source, surfaces[xxz_path],
         xxz_source.replace("for complex q this need not be purely imaginary", "the field shift is imaginary")),
        ("Xxz removed eigenmode-rate boundary", xxz_source, surfaces[xxz_path],
         xxz_source.replace("unchanged dissipator entries do not imply unchanged eigenmode real parts", "eigenmode real parts remain unchanged")),
    ):
        check(f"Route B mechanism rejects {mutation_name}",
              bool(set(defects(mutated, required)) - set(defects(original, required))))
    alternative_promotion = current_arc + " Residual semisimplicity follows from free-fermion integrability at every N."
    wrong_exponent = current_arc.replace(
        "exponent that reads the local eigenvalue split (~1 for a linear analytic split, ~0.5 for a sqrt branch)",
        "exponent that classifies a coalescence (~1 for a diabolic / semisimple, ~0.5 for a defective EP / sqrt-branch)")
    for mutation_name, mutated_current in (("alternative all-N residual promotion", alternative_promotion),
                                           ("old exponent character classifier", wrong_exponent),
                                           ("removed JW-alone boundary", current_arc.replace(
                                               "Free-fermion (Jordan-Wigner) integrability alone does not prove residual semisimplicity or a Delta response at every N", "")),
                                           ("removed exponent limitation", current_arc.replace(
                                               "the exponent corroborates branch type but does not decide geometric multiplicity or semisimplicity", ""))):
        check(f"Route B mechanism rejects current OpenArcs {mutation_name}",
              bool(set(defects(mutated_current, surfaces[arc_path]))
                   - set(defects(current_arc, surfaces[arc_path]))))
    journal_delimiter = "--- THE JOURNAL (how we got here) ---"
    wrong_journal = arc_source.replace(journal_delimiter,
                                      journal_delimiter + " Integrability governs existence at every N.")
    check("Route B mechanism keeps journal assertions outside current-truth slice",
          route_b_current_arc_slice(wrong_journal) == route_b_current_arc_slice(arc_source))
    try:
        route_b_current_arc_slice(arc_source.replace(journal_delimiter, ""))
        missing_delimiter_rejected = False
    except ValueError:
        missing_delimiter_rejected = True
    check("Route B mechanism rejects a missing OpenArcs journal delimiter", missing_delimiter_rejected)

    for relative_path, required in surfaces.items():
        source = (root / relative_path).read_text(encoding="utf-8")
        if relative_path == arc_path:
            source = route_b_current_arc_slice(source)
        if relative_path.endswith("README.md"):
            source = next(line for line in source.splitlines() if "[Path-K Diabolic]" in line)
        if relative_path == "docs/ANALYTICAL_FORMULAS.md":
            source = "\n".join(line for line in source.splitlines()
                               if ("zz(b̄) = zz(b)" in line or "WHY diabolic (Tier 2 hypothesis)" in line
                                   or line.startswith("**Δ-tightness.")))
        errors = defects(source, required)
        check(f"Route B mechanism current truth: {relative_path}", not errors, ascii("; ".join(errors)))
        # Feed wrong prose through the same door: missing scope and restored
        # all-N claims must both fail. This guards the words, not the physics.
        check(f"Route B mechanism rejects missing scope: {relative_path}",
              bool(defects(source.casefold().replace(required[0].casefold(), ""), required)))
        check(f"Route B mechanism rejects all-N promotion: {relative_path}",
              bool(defects(source + " Integrability governs existence at every N.", required)))
        if relative_path == "experiments/F89_PATH_K_DIABOLIC.md":
            check("Route B mechanism rejects the half-hopping energy normalization",
                  bool(defects(source.replace("ε_k = 4J cos", "ε_k = 2J cos"), required)))
        if relative_path.endswith("/BranchLocusPalindromeWitness.cs"):
            check("Route B mechanism rejects q_EP with contradictory J=gamma=1 label",
                  bool(defects(source.replace("γ=1, J=q_EP", "γ=J=1"), required)))
        if relative_path.endswith("/XxzCoherenceBlockFieldTests.cs"):
            for mutation_name, wrong_text in (
                ("imaginary-only field", "The field only shifts the imaginary frequency."),
                ("integrability break", "A random field breaks integrability."),
                ("eigenmode-rate inference", "The absorption-theorem real rate untouched."),
            ):
                check(f"Route B mechanism rejects field-test {mutation_name}",
                      bool(set(defects(source + wrong_text, required)) - set(defects(source, required))))


def verify_random_field_current_truth():
    """Guard Hamiltonian knob distinctions, not an integrability or CSR theorem."""
    root = Path(__file__).resolve().parents[1]
    distinction = (
        "uniform XXZ remains Bethe-integrable",
        "random longitudinal Z disorder at Delta=0 remains quadratic",
        "random field plus Delta!=0 is the interacting disordered nonintegrable test",
        "finite executed CSR evidence",
    )
    surfaces = {
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzDeltaFlipTests.cs": (
            "base-label nearest-neighbour continuation", "proposal tracking, not an invariant-set proof"),
        "compute/RCPsiSquared.Core/F89PathK/WeightCoherenceBlock.cs": (
            "diagonal addition -i*q*(fe(ket)-fe(bra))", "for complex q it need not be purely imaginary",
            "unchanged dissipator entries do not imply unchanged eigenmode real parts"),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": distinction,
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdWitness.cs": distinction,
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdCsr.cs": distinction,
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/FillingThresholdCsrTests.cs": distinction,
        "compute/RCPsiSquared.Diagnostics/Foundation/IntegrabilityBreakingCsr.cs": (
            "recomputed eigenvalues by their current real parts", "not an invariant AT/residual decomposition",
            "unchanged dissipator does not fix eigenmode real parts"),
        "experiments/FILLING_THRESHOLD_CHAOS.md": distinction,
        "experiments/RANDOM_MATRIX_THEORY.md": distinction + (
            "Poisson statistics do not prove integrability",
            "does not run back up the chain"),
    }
    forbidden = (
        "is Delta-stable (ZZ Hermitian", "every integrability-breaking knob",
        "integrability-/symmetry-breaking disorder knob",
        "Hamiltonian's integrability (a random Z-field), with or without interactions",
        "random field breaks integrability", "Delta is integrability-breaking",
        "random field is nonintegrable with or without interactions",
        "however hard you break integrability", "however hard you break it",
        "merge only at extensive filling", "cannot thermalize",
        "AT real-part split is unchanged under Delta", "which Delta leaves untouched",
        "spectrum is Poisson (integrable", "Poisson (Integrable)",
        "palindromic Liouvillian is integrable", "system is integrable in every band",
        "spectrum is integrable", "bands are integrable", "bands stay integrable",
        "chiral but integrable", "integrable chiral Lindbladian", "SFF) confirms integrability",
    )

    def errors(source, required):
        normalized = " ".join(source.replace("///", "").replace("//", "")
                              .replace('" +', "").replace('"', "").replace("**", "")
                              .replace("Δ", "Delta").split()).casefold()
        return ([f"missing {phrase}" for phrase in required if phrase.casefold() not in normalized]
                + [f"forbidden {phrase}" for phrase in forbidden if phrase.casefold() in normalized])

    for path, required in surfaces.items():
        source = (root / path).read_text(encoding="utf-8")
        baseline = set(errors(source, required))
        check(f"Random-field current truth: {path}", not baseline, ascii("; ".join(sorted(baseline))))
        for phrase in forbidden:
            check(f"Random-field rejects {ascii(phrase)}: {path}",
                  bool(set(errors(source + " " + phrase, required)) - baseline))
        for phrase in required:
            check(f"Random-field requires {ascii(phrase)}: {path}",
                  bool(set(errors(source.replace(phrase, ""), required)) - baseline))


def verify_poisson_inference_current_truth():
    """Poisson/no-ramp observations are not converse proofs of integrability."""
    root = Path(__file__).resolve().parents[1]
    surfaces = {
        "docs/outbound/SHIFTED_ORDER4_CHIRAL_SYMMETRY.md": (
            "does not prove integrability", "AIII is a symmetry label",
            "uniform XXZ remains Bethe-integrable", "Delta=0 remains quadratic",
            "Delta!=0 plus generic disorder", "finite-N"),
        "hypotheses/DIABOLIC_BY_INTEGRABILITY.md": (
            "uniform XXZ remains Bethe-integrable", "does not isolate integrability causality",
            "Tier-2 hypothesis"),
        "docs/ANALYTICAL_FORMULAS.md": ("does not prove integrability",),
        "docs/THE_INTERPRETATION.md": ("does not prove integrability", "AIII is a symmetry label"),
        "docs/KMS_DETAILED_BALANCE.md": (
            "statistics do not remove the symmetry class or prove integrability", "random-matrix universality"),
        "experiments/TOPOLOGICAL_EDGE_MODES.md": ("do not prove integrability or topology",),
        "experiments/SPECTRAL_FORM_FACTOR.md": ("does not prove integrability",),
        "experiments/RANDOM_MATRIX_THEORY.md": ("AIII is a symmetry label",),
    }
    forbidden = ("Poisson (integrable)", "Poisson, integrable", "SFF confirms integrability",
                 "Poisson statistics of the integrable system", "precisely because it is the integrable member",
                 "single integrability-breaking knob", "outside all 38", "does not fit any of the 38",
                 "integrable chiral Lindbladian", "Liouvillian is integrable",
                 "Poisson under integrability breaking", "statistics of none of the 38",
                 "integrable at every N")
    for path, required in surfaces.items():
        source = (root / path).read_text(encoding="utf-8")
        if path == "docs/ANALYTICAL_FORMULAS.md":
            source = source.split("### F42.", 1)[1].split("### F43.", 1)[0]

        def errors(text):
            normalized = " ".join(text.replace("**", "").replace("*", "").split()).casefold()
            return ([f"missing {phrase}" for phrase in required if phrase.casefold() not in normalized]
                    + [f"forbidden {phrase}" for phrase in forbidden if phrase.casefold() in normalized])

        baseline = set(errors(source))
        check(f"Poisson implication current truth: {path}", not baseline, ascii("; ".join(sorted(baseline))))
        for phrase in forbidden:
            check(f"Poisson implication rejects {ascii(phrase)}: {path}",
                  bool(set(errors(source + " " + phrase)) - baseline))
        for phrase in required:
            check(f"Poisson implication requires {ascii(phrase)}: {path}",
                  bool(set(errors(source.replace(phrase, ""))) - baseline))


def verify_csr_producer_and_crossover_current_truth():
    root = Path(__file__).resolve().parents[1]
    crossover = ("filling-associated crossover evidence", "N=6..8", "movement toward GinUE",
                 "not a causal or thermodynamic threshold theorem")
    surfaces = {
        "experiments/README.md": crossover + ("does not prove integrability",),
        "experiments/F89_MULTI_SECTOR_MONODROMY.md": crossover,
        "experiments/RANDOM_MATRIX_THEORY.md": crossover,
        "experiments/FILLING_THRESHOLD_CHAOS.md": crossover,
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdWitness.cs": crossover,
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/IntegrabilityBreakingCsrTests.cs": (
            "random field at Delta=0 remains Anderson/free-fermion",
            "random field plus Delta!=0 is interacting disordered nonintegrable"),
        "compute/RCPsiSquared.Core/F89PathK/WeightCoherenceBlock.cs": (
            "-i*q*Delta*(zz(ket)-zz(bra))", "for Im q != 0 it has a real component"),
        "simulations/spectral_form_factor.py": ("does not prove integrability", "not a spectral-class or no-ramp test"),
        "simulations/results/spectral_form_factor.txt": ("does not prove integrability", "not a spectral-class or no-ramp test"),
        "simulations/rmt_topology_csr.py": (
            "does not prove integrability", "simultaneous Hamiltonian/dissipator eigenoperators",
            "not a formula for mixed Liouvillian modes"),
    }
    forbidden = ("Liouvillian is integrable", "Poisson (integrable)", "is INTEGRABLE. No dip",
                 "-> Heisenberg integrability", "is chaotic", "chaos located as a FILLING threshold",
                 "dissipative chaos is a filling threshold", "turns on with filling",
                 "at the GinUE value already", "DENSE block IS chaotic",
                 "random field U[−W,W] breaks integrability",
                 "adds only the frequency", "Poisson-like SFF at", "GUE-like SFF at",
                 "Classification: POISSON-like", "Classification: GUE-like",
                 "Classification: intermediate / palindromic")
    for path, required in surfaces.items():
        source = (root / path).read_text(encoding="utf-8")
        if path == "experiments/README.md":
            source = "\n".join(line for line in source.splitlines()
                               if "[Random Matrix Theory]" in line or "[Filling Threshold Chaos]" in line)

        def errors(text):
            normalized = " ".join(text.replace("///", "").replace('" +', "").replace('"', "")
                                  .replace("**", "").split()).casefold()
            return ([f"missing {phrase}" for phrase in required if phrase.casefold() not in normalized]
                    + [f"forbidden {phrase}" for phrase in forbidden if phrase.casefold() in normalized])

        baseline = set(errors(source))
        check(f"CSR producer/crossover current truth: {path}", not baseline, ascii("; ".join(sorted(baseline))))
        for phrase in forbidden:
            check(f"CSR producer/crossover rejects {ascii(phrase)}: {path}",
                  bool(set(errors(source + " " + phrase)) - baseline))
        for phrase in required:
            check(f"CSR producer/crossover requires {ascii(phrase)}: {path}",
                  bool(set(errors(source.casefold().replace(phrase.casefold(), ""))) - baseline))

    # The topology producer is bounded and stdout-only. Check the denominator
    # against the actual CSV population, not a prose percentage or all eigenvalues.
    output = subprocess.run([sys.executable, str(root / "simulations/rmt_topology_csr.py")],
                            check=True, capture_output=True, text=True, encoding="utf-8").stdout
    complete = np.loadtxt(root / "simulations/results/rmt_eigenvalues_complete_N6.csv", delimiter="\t", skiprows=1)
    upper = int(np.count_nonzero(complete[:, 1] > 1e-6))

    def valid_collapse(text):
        match = re.search(r"complete N=6: distinct_upper=(\d+) total_upper=(\d+) collapse=([\d.eE+-]+)", text)
        table = re.search(r"^complete\s+6\s+\d+\s+(\d+)\s", text, re.MULTILINE)
        if match is None or table is None:
            return False
        distinct, reported_upper, collapse = int(match[1]), int(match[2]), float(match[3])
        return (reported_upper == upper and distinct == int(table[1])
                and collapse == 1 - distinct / reported_upper)

    check("Topology collapse uses the executed upper-half population", valid_collapse(output))
    wrong_denominator = re.sub(r"total_upper=\d+", f"total_upper={len(complete)}", output)
    reported_distinct = re.search(r"distinct_upper=(\d+)", output)
    if reported_distinct is not None:
        # Keep the wrong fraction internally consistent: only the population
        # boundary, not an arithmetic mismatch, must reject all-spectrum counting.
        wrong_denominator = re.sub(r"collapse=[\d.eE+-]+",
                                   f"collapse={1 - int(reported_distinct[1]) / len(complete):.17g}",
                                   wrong_denominator)
    check("Topology collapse rejects all-spectrum denominator", valid_collapse(output) and not valid_collapse(wrong_denominator))


def verify_sff_windows_and_current_boundaries():
    root = Path(__file__).resolve().parents[1]
    # The dispersion's Taylor coefficient is computed independently of prose:
    # 4J(1-cos x) = J * coefficient * x^2 + O(x^4).
    import sympy as sp
    angle = sp.Symbol("angle")
    coefficient = sp.diff(4 * (1 - sp.cos(angle)), angle, 2).subs(angle, 0) / 2
    sff_doc = (root / "experiments/SPECTRAL_FORM_FACTOR.md").read_text(encoding="utf-8")

    def asymptotic_matches_dispersion(text):
        matches = re.findall(r"ω_min\s*~\s*(\d*)J?π²/N²", text)
        return len(matches) == 1 and sp.Integer(matches[0] or "1") == coefficient

    check("SFF asymptotic coefficient equals exact dispersion Taylor coefficient",
          asymptotic_matches_dispersion(sff_doc))
    check("SFF asymptotic rejects half coefficient with J",
          not asymptotic_matches_dispersion("ω_min ~ Jπ²/N²"))
    check("SFF asymptotic rejects half coefficient at J=1",
          not asymptotic_matches_dispersion("ω_min ~ π²/N²"))
    check("SFF asymptotic rejects missing coefficient formula",
          not asymptotic_matches_dispersion("t_Π grows as N²"))
    surfaces = {
        "hypotheses/DIABOLIC_BY_INTEGRABILITY.md": (
            "N=4 twin-scalar restriction", "Tier-2 residual interpretation",
            "uniform XXZ remains Bethe-integrable", "N=6 inventory is complete", "N=7 remains the inventory boundary"),
        "experiments/SPECTRAL_FORM_FACTOR.md": (
            "raw non-unfolded frequency SFF", "1/M", "not sampled", "reached window",
            "not a universality-class proof", "overflow/NaN remains unresolved",
            "raw multiset density scale", "multiplicity-dependent", "N=4 third", "N=6 second", "N=7 has no match", "114%"),
        "experiments/README.md": ("raw non-unfolded frequency SFF", "reached window", "not a universality-class proof", "raw multiset density scale"),
        "experiments/RANDOM_MATRIX_THEORY.md": ("cannot classify universality",),
        "docs/outbound/SHIFTED_ORDER4_CHIRAL_SYMMETRY.md": ("measured finite-N contrast", "convergence is a proposed extension"),
        "simulations/spectral_form_factor.py": ("independent_phase_reference", "format_sample", "raw non-unfolded frequency SFF",
                                               "not a physical time-scale estimate", "second-half-of-current-grid mean",
                                               "raw multiset density scale", "multiplicity-dependent", "RCPSI_SFF_OUTPUT_PATH", "raw_multiset_density_scale(ev)"),
        "simulations/sff_window_summary.py": ("def sff_frequency", "def raw_multiset_density_scale", "multiplicity-dependent"),
        "simulations/results/spectral_form_factor.txt": (
            "independent-phase reference 1/M", "not sampled", "reached window", "not classifiable",
            "not a physical time-scale estimate", "second-half-of-current-grid mean",
            "raw multiset density scale", "multiplicity-dependent"),
    }
    arc = (root / "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs").read_text(encoding="utf-8")
    entry = arc.split('Name: "f89_galois_open_doors"', 1)[1].split("new OpenArc(", 1)[0]
    current = entry.split("NextStep:", 1)[1].split("Status: OpenArcStatus", 1)[0]
    current = " ".join(current.replace('" +', "").replace('"', "").split())
    arc_label = "OpenArcs f89_galois current NextStep"
    surfaces[arc_label] = ("finite N=6..8", "not a causal or thermodynamic threshold theorem",
                          "XXZ chain remains Bethe-integrable", "Delta=0 remains a quadratic Anderson/free-fermion Hamiltonian",
                          "not proof of convergence, thermalization, or a thermodynamic threshold")
    forbidden = ("Because the Magnons Don't Interact", "and that is the protection",
                 "mechanism probe shows exactly the predicted cause", "integrability survives by elimination",
                 "Generalize beyond N=4 path-3", "at path-5/6", "Generic XXZ gives the ordinary defective EP",
                 "prediction is GinUE convergence", "structurally pinned Poisson",
                 "K(t) = 1 for all t > 0", "plateau at K=1", "Poisson: 1.0", "plateau: 1.0",
                 "Neither Poisson nor GUE: unique palindromic signature", "UNIQUE SFF signature",
                 "a fingerprint unique to", "chaos is a filling threshold", "dense block thermalizes",
                 "Thouless time (where K_freq", "first rises above Poisson baseline",
                 "no clear Thouless time detected", "1.5× its late-time mean",
                 "t_H", "Heisenberg time", "t_Th", "Thouless-time estimate",
                 "structurally Poisson", "neither standard Poisson nor GUE",
                 "dominant FFT peak matches", "ω_min peak is present")
    for path, required in surfaces.items():
        source = current if path == arc_label else (root / path).read_text(encoding="utf-8")
        if path == "experiments/README.md":
            source = next(line for line in source.splitlines() if "[Spectral Form Factor]" in line)
        # Mutate the same whitespace-normalized surface the validator reads,
        # including required fragments that cross Markdown line boundaries.
        source = " ".join(source.replace("**", "").replace('" +', "").replace('"', "").split())

        def errors(text):
            normalized = " ".join(text.replace("**", "").replace('" +', "").replace('"', "").split()).casefold()
            return ([f"missing {phrase}" for phrase in required if phrase.casefold() not in normalized]
                    + [f"forbidden {phrase}" for phrase in forbidden
                       if (bool(re.search(r"(?<!\w)" + re.escape(phrase.casefold()) + r"(?!\w)", normalized))
                           if phrase in ("t_H", "t_Th") else phrase.casefold() in normalized)])

        baseline = set(errors(source))
        check(f"SFF windows/current boundary: {path}", not baseline, ascii("; ".join(sorted(baseline))))
        for phrase in forbidden:
            check(f"SFF windows/current boundary rejects {ascii(phrase)}: {path}",
                  bool(set(errors(source + " " + phrase)) - baseline))
        for phrase in required:
            check(f"SFF windows/current boundary requires {ascii(phrase)}: {path}",
                  bool(set(errors(source.casefold().replace(phrase.casefold(), ""))) - baseline))
    # Origin/ParkedAt are intentionally outside this second arc's current surface.
    check("Galois current gate excludes its dated origin",
          "the philosophical-zoom-out reading" not in current and "The perturbation controls distinguish three cases" in current)

    # Bind each executed Phase-4 row to its own mode count and reached window;
    # a correct label elsewhere must not hide a false zero or unit reference.
    output = (root / "simulations/results/spectral_form_factor.txt").read_text(encoding="utf-8")
    phase4 = output.split("PHASE 4:", 1)[1].split("PHASE 5:", 1)[0]

    def window_report_errors(text):
        problems = []
        blocks = re.findall(r"  N=(\d+):\n(.*?)(?=\n  N=|\Z)", text, re.S)
        if [int(n) for n, _ in blocks] != [3, 4, 5, 6, 7]:
            problems.append("missing or duplicate N=3..7 window reports")
        for n_text, block in blocks:
            n = int(n_text)
            header = re.search(r"Reached time: ([\d.]+); raw scale = ([\d.]+); independent-phase reference 1/M = ([\deE.+-]+)", block)
            if header is None:
                problems.append(f"N={n} missing normalization/window header")
                continue
            reached, raw_scale, reference = map(float, header.groups())
            data = np.loadtxt(root / f"simulations/results/rmt_eigenvalues_N{n}.csv", delimiter="\t", skiprows=1)
            frequencies = np.abs(data[:, 1][np.abs(data[:, 1]) > 1e-10])
            # Independent telescoping form of the mean-adjacent-gap denominator:
            # duplicated frequencies count, rather than treating a printed scale as evidence.
            expected_scale = 2 * np.pi * (len(frequencies) - 1) / np.ptp(frequencies)
            period = 2 * np.pi / (4 * (1 - np.cos(np.pi / n)))
            expected_end = min(3 * expected_scale, 50 * period, 200)
            if abs(raw_scale - expected_scale) > 0.0051 or abs(reached - expected_end) > 0.0051:
                problems.append(f"N={n} raw scale/grid disagrees with CSV multiset")
            # Producer prints eight significant digits: allow its last-place
            # rounding, not a normalization-scale change (unit control below).
            if not np.isclose(reference, 1 / len(data), rtol=5e-8, atol=0):
                problems.append(f"N={n} wrong M^-2 independent-phase reference")
            for label, reached_window in (("Beyond", expected_end > expected_scale), ("Intermediate", expected_end > 0.1 * expected_scale)):
                line = next((line for line in block.splitlines() if line.strip().startswith(label + " (")), "")
                if not line or (("<K> = not sampled" in line) == reached_window):
                    problems.append(f"N={n} {label} window availability mismatch")
                if label == "Intermediate" and not reached_window:
                    if "slope = not sampled" not in line or "not classifiable" not in block:
                        problems.append(f"N={n} unreached slope classified")
        return problems

    check("SFF executed per-N normalization/window provenance", not window_report_errors(phase4),
          ascii("; ".join(window_report_errors(phase4))))
    for name, mutated in (
        ("unit plateau", phase4.replace("1/M = 0.015625", "1/M = 1")),
        ("empty beyond bin as zero", phase4.replace("Beyond (t > raw scale):        <K> = not sampled", "Beyond (t > raw scale):        <K> = 0.0000")),
        ("empty intermediate bin as zero", phase4.replace("Intermediate (0.1-1.0 raw scale):    <K> = not sampled", "Intermediate (0.1-1.0 raw scale):    <K> = 0.0000")),
        ("absent slope as zero", phase4.replace("slope = not sampled", "slope = 0.0000")),
        ("printed raw scale without CSV support", phase4.replace("raw scale = 5810.10", "raw scale = 100.00")),
    ):
        check(f"SFF window report rejects {name}", mutated != phase4 and bool(window_report_errors(mutated)))


def route_b_n6_current_truth_errors(source, label):
    """Return N6 provenance/count defects; source has normalized Markdown/C# whitespace."""
    marker = "Route B is complete at N=6 as a parity-labelled direct-t inventory"
    if marker not in source:
        return ["missing N6 current-truth section"]
    source = source[source.index(marker):]
    if label == "arc":
        # Only the current pre-journal boundary may satisfy the contract.
        source = source.split("Outside Route B", 1)[0]
    required = ["ROUTE_B_A2_N6_RECONCILE", "0 unresolved", "0.06589058580248003"]
    if label == "primary":
        required += ["consumes all 266 locus IDs once in ordinal order", "ExactAlgebraic=266",
                     "the loci and algebraic multiplicity, while character has the separate sources",
                     "| R parity | direct-t loci | HermitianAxis | EpCharacterStable | ExactRankExecuted | verdict |",
                     "| even | 133 | 59 | 74 | 0 | 133 Diabolic, alg=geo=2 |",
                     "| odd | 133 | 59 | 74 | 0 | 133 Diabolic, alg=geo=2 |",
                     "| total | 266 | 118 | 148 | 0 | 266 Diabolic; 0 Defective |",
                     "EpCharacterStable records stable numerical character from the isolated two-dimensional restriction; it is not an exact rank proof"]
    elif label == "arc":
        required += ["ExactAlgebraic owns all 266 loci (133 per parity), consumed exactly once with 0 unresolved",
                     "character is Diabolic with alg=geo=2 at every locus",
                     "118 by executed full-sector HermitianAxis and 148 by stable EpCharacterStable readings, with ExactRankExecuted=0",
                     "The numerical rule uses all three intrinsic radii", "No exact fallback was used",
                     "neither determines Jordan character"]
    else:
        required += ["consumes 266 distinct loci once", "ExactAlgebraic=266",
                     "HermitianAxis=118, EpCharacterStable=148, ExactRankExecuted=0",
                     "parity contributes 59 HermitianAxis and 74 EpCharacterStable readings",
                     "All 266 are Diabolic with alg=geo=2; 0 are Defective",
                     "not an exact rank proof", "exactRankCertificates array is empty"]
        required += (["133 per parity", "EpCharacterStable is numerical evidence"] if label == "path experiment"
                     else ["133 in each parity", "is stable numerical evidence"])
    # The exact algebraic/structural sources do not license a positive exact-rank
    # claim for numerical character. Negated rank-proof wording is the scope fence.
    errors = [f"missing {phrase}" for phrase in required if phrase not in source]
    if "an exact rank proof" in source.casefold().replace("not an exact rank proof", ""):
        errors.append("numerical character promoted to an exact rank proof")
    return errors


def verify_route_b_a2_current_truth():
    """Guard completion prose and inventory shape; C# owns the local character verdicts."""
    verify_route_b_mechanism_current_truth()
    verify_random_field_current_truth()
    verify_poisson_inference_current_truth()
    verify_csr_producer_and_crossover_current_truth()
    verify_sff_windows_and_current_boundaries()
    root = Path(__file__).resolve().parents[1]
    surfaces = (
        ("primary", "docs/THE_DOUBLE_ROOT.md"),
        ("path experiment", "experiments/F89_PATH_K_DIABOLIC.md"),
        ("genericity", "experiments/F89_BETA_EXOTIC_GENERICITY.md"),
        ("arc", "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs"),
    )
    forbidden = (
        "already gated positive-real R-odd diabolic",
        "the arc is effectively complete",
        "Remaining work is pair isolation and character gates",
        "positive-real R-odd A2 locus remains character-uncertified",
        "Route B remains open",
        "N=6 Route B remains open",
        "What remains is local pair extraction for Riesz seeds",
    )
    for label, relative_path in surfaces:
        source = (root / relative_path).read_text(encoding="utf-8")
        # Join C# string literals and whitespace within this file only.
        source = " ".join(source.replace('" +', "").replace('"', "").replace("**", "").split())
        check(f"Route B {label}: N6 executed counts and provenance",
              not route_b_n6_current_truth_errors(source, label),
              ascii("; ".join(route_b_n6_current_truth_errors(source, label))))
        # Mutations retain the source labels while corrupting their counts or
        # promoting their evidential strength. No tracked document is rewritten.
        for mutation_name, wrong_source in (
            ("total count", source.replace("266", "200")),
            ("numerical source count", source.replace("148", "147")),
            ("exact rank promotion", source.replace("EpCharacterStable", "EpCharacterStable (an exact rank proof)")),
        ):
            check(f"Route B N6 {label} mutation rejects {mutation_name}",
                  bool(route_b_n6_current_truth_errors(wrong_source, label)))
        if label == "primary":
            wrong_counts = source.replace("266 Diabolic; 0 Defective", "200 Diabolic; 66 Defective")
            wrong_proof = source.replace(
                "stable numerical character from the isolated two-dimensional restriction; it is not an exact rank proof",
                "an exact rank proof")
            check("Route B N6 mutation rejects 200 Diabolic / 66 Defective",
                  bool(route_b_n6_current_truth_errors(wrong_counts, label)))
            check("Route B N6 mutation rejects numerical character promoted to exact rank",
                  bool(route_b_n6_current_truth_errors(wrong_proof, label)))
        required = ["Route B is complete at N=5 as a root-by-root inventory",
                    "PSC1", "S1", "all 34 stable EpCharacter",
                    "No exact fallback was used", "not an all-N theorem"]
        # Sensitive to deleting any N=6 completion/provenance boundary or restoring
        # the obsolete N=6-open statement; numerical character must stay numerical.
        required += ["Route B is complete at N=6 as a parity-labelled direct-t inventory",
                     "ExactAlgebraic", "HermitianAxis", "EpCharacterStable", "ExactRankExecuted",
                     "ROUTE_B_A2_N6_RECONCILE", "F_53 (N=7)"]
        if label != "arc":
            required += ["29 A2(w) roots / 58 q-loci", "24 imaginary-q loci",
                         "executed full-sector Hermiticity", "alg=geo=2",
                         "exactRankCertificates array is empty", "ROUTE_B_A2_RECONCILE"]
        if label == "primary":
            required += ["Route B is complete at N=5 as a root-by-root inventory: 29 A2(w) roots / 58 q-loci, all semisimple with alg=geo=2.",
                         "12 negative-real w roots", "2 real-q loci", "32 nonreal-q loci",
                         "three isolating radii", "neither is a local Jordan character test"]
        surface_forbidden = forbidden
        if label == "path experiment":
            required.append("N=6 Route B is complete as the parity-labelled direct-t A2 inventory")
            surface_forbidden += (
                "higher-N completeness remains open",
                "positive-real R-odd A2 locus; its Jordan character remains open",
                "positive-real R-odd A2 locus requires its own local character determination",
                "N=5 R-odd A2 character remains open",
            )
        if label == "genericity":
            surface_forbidden += (
                "5.100831, whose character remains uncertified",
            )
        if label == "arc":
            required.append("NEXT: the F_53 (N=7) doubled-layer character/completeness edge remains unmeasured; do not extrapolate the N=6 inventory.")
            surface_forbidden += (
                "this arc is COMPLETE but for ONE minor non-blocking edge",
                "what remains is ONLY a closed form for the threshold N",
                "Outside the still-open Route B inventory",
                "This does not close Route B or the positive R-odd A2 character.",
            )
        for phrase in required:
            check(f"Route B {label}: {phrase}", phrase in source,
                  f"{relative_path}: missing current-truth statement")
        for phrase in surface_forbidden:
            check(f"Route B {label}: excludes {phrase}",
                  phrase.casefold() not in source.casefold(),
                  f"{relative_path}: stale Route B claim")

    # The exact producer/export tests own root isolation and PSC1/S1 validity.
    # Do not infer no-NearEp or classification sources from this JSON: those
    # are executed on the full parity operators by ROUTE_B_A2_RECONCILE.
    try:
        inventory = json.loads((root / "simulations/results/route_b_a2_n5.json").read_text(encoding="utf-8"))
        roots = [entry for sector in inventory["sectors"] for entry in sector["a2Roots"]]
        loci = [locus for entry in roots for locus in entry["qLoci"]]
        check("Route B artifact: schema 1, N=5", inventory["schemaVersion"] == 1 and inventory["n"] == 5)
        check("Route B artifact: 29 distinct A2 roots", len(roots) == len({r["id"] for r in roots}) == 29)
        check("Route B artifact: 58 distinct q loci", len(loci) == len({q["id"] for q in loci}) == 58)
        check("Route B artifact: root-kind split 12/1/16",
              tuple(sum(r["rootKind"] == kind for r in roots) for kind in ("negativeReal", "positiveReal", "nonreal")) == (12, 1, 16))
        check("Route B artifact: no exact rank certificates used", inventory["exactRankCertificates"] == [])
    except (OSError, ValueError, KeyError, TypeError) as error:
        check("Route B artifact: readable schema", False, str(error))


def verify_f151_route_b_current_truth():
    """Keep F151's bounded N=5 scan distinct from the completed Route B inventory."""
    path = Path(__file__).resolve().parents[1] / "docs/ANALYTICAL_FORMULAS.md"
    try:
        source = path.read_text(encoding="utf-8")
        section = source[source.index("### F151."):].split("\n### F", 1)[0]
    except (OSError, ValueError) as error:
        check("F151 Route B: readable section", False, str(error))
        return
    section = " ".join(section.split())
    check("F151 Route B: excludes unqualified N=5 absence",
          "N=5 is odd and shows none" not in section,
          "docs/ANALYTICAL_FORMULAS.md F151: stale unqualified N=5 claim")
    for phrase in ("bounded R-even real-q scan at N=5 shows none",
                   "positive-real R-odd A2 pair"):
        check(f"F151 Route B: {phrase}", phrase in section,
              "docs/ANALYTICAL_FORMULAS.md F151: missing current-truth statement")


if __name__ == "__main__":
    print("=" * 70)
    print("SCOPED ROUTE-B CURRENT-TRUTH AND NUMERICAL DOCUMENTATION GATE")
    print("=" * 70)

    if "--force-failure" in sys.argv:
        check("forced failure proves the process exits nonzero", False,
              "intentional mutation regression")

    verify_route_b_a2_current_truth()
    verify_docs_gate_scope_description()
    verify_f151_route_b_current_truth()

    # =============================================================
    # GLOSSARY.md claims
    # =============================================================
    print("\n--- GLOSSARY: Boundary values ---")

    # Claim: "2gamma = decay rate of c+ supermode (N=3)"
    # Claim: "8gamma/3 = concurrence envelope (N=3)"
    # Claim: "10gamma/3 = decay rate of c- supermode (N=3)"
    gamma = 0.05
    N = 3
    H = build_H(N, 1.0, "chain")
    L = build_L(H, [gamma]*N, N)
    evals = np.linalg.eigvals(L)
    reals = np.real(evals)
    # Get unique non-zero real parts
    nonzero_reals = reals[np.abs(reals) > 1e-10]
    unique_rates = np.sort(np.unique(np.round(nonzero_reals, 8)))

    expected_rates = [-2*gamma, -8*gamma/3, -10*gamma/3, -2*3*gamma]
    print(f"  gamma = {gamma}, N = 3, chain")
    print(f"  Expected rates: {[f'{r:.6f}' for r in sorted(expected_rates)]}")
    print(f"  Unique rates found: {[f'{r:.6f}' for r in unique_rates[:8]]}")

    for rate_name, expected in [("2gamma", -2*gamma),
                                 ("8gamma/3", -8*gamma/3),
                                 ("10gamma/3", -10*gamma/3)]:
        found = any(abs(r - expected) < 1e-3 for r in unique_rates)
        check(f"{rate_name} = {expected:.6f} exists in spectrum (tol=1e-3)", found,
              f"expected {expected:.6f}")

    # These are shared approximate anchors, not equality of the rate multisets.
    print("\n--- GLOSSARY: Shared N=3 rate anchors across sampled topologies ---")
    for topo in ["chain", "ring", "star"]:
        H = build_H(3, 1.0, topo)
        L = build_L(H, [0.05]*3, 3)
        evals = np.linalg.eigvals(L)
        reals = np.real(evals)
        for rate_name, expected in [("2gamma", -0.10), ("8gamma/3", -8*0.05/3), ("10gamma/3", -10*0.05/3)]:
            found = any(abs(r - expected) < 1e-3 for r in reals)
            check(f"{topo} contains anchor {rate_name} (tol=1e-3)", found)

    wrong_anchor = -0.115
    check("mutated shared-anchor value is rejected",
          not any(abs(r - wrong_anchor) < 1e-3 for r in reals),
          f"unexpected anchor near {wrong_anchor}")

    # Claim: "-2*sum_gamma = location of XOR modes"
    print("\n--- GLOSSARY: XOR modes at -2*sum_gamma ---")
    for N in [2, 3, 4]:
        gammas = [0.05]*N
        H = build_H(N, 1.0, "chain")
        L = build_L(H, gammas, N)
        evals = np.linalg.eigvals(L)
        reals = np.real(evals)
        target = -2*sum(gammas)
        count = np.sum(np.abs(reals - target) < 1e-6)
        check(f"N={N}: {count} modes at -2*sum_gamma={target:.3f} (expect N+1={N+1})",
              count == N+1, f"got {count}")

    # Claim: "0.039/gamma = approximate crossing time"
    print("\n--- GLOSSARY: t_cross = 0.039/gamma ---")
    # Simulate Bell+ under Heisenberg chain N=2 with dephasing
    from scipy.linalg import expm
    gamma = 0.05
    H2 = build_H(2, 1.0, "chain")
    L2 = build_L(H2, [gamma]*2, 2)
    bell = np.array([1,0,0,1], dtype=complex)/np.sqrt(2)
    rho0 = np.outer(bell, bell.conj())
    rho_vec = rho0.flatten()

    t_cross_found = None
    for t in np.linspace(0.01, 2.0, 2000):
        rho_t = (expm(L2 * t) @ rho_vec).reshape(4,4)
        C = concurrence_2q(rho_t)
        # l1 coherence normalized
        l1 = np.sum(np.abs(rho_t)) - np.sum(np.abs(np.diag(rho_t)))
        Psi = l1 / 3  # d-1 = 3 for d=4
        CPsi = C * Psi
        if CPsi < 0.25 and t_cross_found is None:
            t_cross_found = t
            break

    expected_tcross = 0.039 / gamma
    if t_cross_found:
        ratio = t_cross_found * gamma  # Should be ~0.039
        check(f"t_cross = {t_cross_found:.3f}, gamma*t_cross = {ratio:.4f} (expect ~0.039, tol 10%)",
              abs(ratio - 0.039) / 0.039 < 0.10,
              f"gamma*t_cross = {ratio:.4f}, expect ~0.039")
    else:
        check("t_cross found", False, "never crossed 1/4")

    # =============================================================
    # MIRROR_SYMMETRY_PROOF claims
    # =============================================================
    print("\n--- MIRROR_SYMMETRY_PROOF: Palindromic pairing ---")

    # Claim: "every decay rate d paired with 2*sum_gamma - d"
    for N in [2, 3, 4, 5]:
        gammas = [0.05]*N
        sg = sum(gammas)
        H = build_H(N, 1.0, "chain")
        L = build_L(H, gammas, N)
        evals = np.linalg.eigvals(L)
        distance = max_f1_pairing_distance(evals, sg)
        tolerance = 1e-5
        check(f"N={N}: full eigenvalue multiset is palindromically paired",
              distance < tolerance, f"max multiplicity-aware distance={distance:.3e}")

    # Mutation-sensitive controls for the global matcher. The fixed point is
    # lambda=-Sigma, not -2*Sigma; and deleting the zero endpoint must expose
    # the unmatched -2*Sigma endpoint.
    control_sigma = 0.25
    fixed_point = np.array([-control_sigma + 0.2j, -control_sigma - 0.2j])
    check("F1 fixed-point control is centered at lambda=-Sigma",
          max_f1_pairing_distance(fixed_point, control_sigma) < 1e-14)
    endpoint_without_zero = np.array([-2 * control_sigma, -control_sigma, -control_sigma])
    endpoint_distance = max_f1_pairing_distance(endpoint_without_zero, control_sigma)
    check("F1 endpoint-without-zero negative control fails pairing",
          endpoint_distance > 0.1,
          f"broken multiset residual unexpectedly {endpoint_distance:.3e}")

    # Claim: the full F1 multiset closure holds on each supported topology.
    print("\n--- MIRROR_SYMMETRY_PROOF: Full multiset closure by topology ---")
    for topo in ["chain", "ring", "star"]:
        for N in [3, 4]:
            gammas = [0.05]*N
            sg = sum(gammas)
            H = build_H(N, 1.0, topo)
            L = build_L(H, gammas, N)
            evals = np.linalg.eigvals(L)
            distance = max_f1_pairing_distance(evals, sg)
            check(f"N={N} {topo}: full eigenvalue multiset has F1 closure",
                  distance < 1e-5,
                  f"max multiplicity-aware distance={distance:.3e}")

    # A Hamiltonian outside the supported Pi-odd class must not make the
    # topology loop false-green merely because all rates stay in range.
    bad_n = 3
    bad_gammas = [0.05] * bad_n
    bad_h = build_H(bad_n, 1.0, "chain") + 0.37 * site_op(sz, 0, bad_n)
    bad_spectrum = np.linalg.eigvals(build_L(bad_h, bad_gammas, bad_n))
    bad_distance = max_f1_pairing_distance(bad_spectrum, sum(bad_gammas))
    check("non-Pi-compatible Hamiltonian negative control breaks F1 closure",
          bad_distance > 1e-3,
          f"broken-control distance unexpectedly {bad_distance:.3e}")

    print("\n--- XOR_SPACE: coordinate claims retired ---")
    print("  Right-eigenvector coordinate squares are not invariant state weights.")
    print("  Run f22_operator_charge.py for the operator-level F22 gate.")

    # =============================================================
    # CORE_ALGEBRA claims
    # =============================================================
    print("\n--- CORE_ALGEBRA: 1/4 boundary ---")

    # Claim: "discriminant of R = C(Psi+R)^2 changes sign at C*Psi = 1/4"
    # R = C(Psi+R)^2 -> R = C*Psi^2 + 2C*Psi*R + C*R^2
    # -> C*R^2 + (2C*Psi - 1)*R + C*Psi^2 = 0
    # discriminant = (2C*Psi - 1)^2 - 4*C^2*Psi^2 = 1 - 4C*Psi
    for CPsi in [0.1, 0.24, 0.25, 0.26, 0.5]:
        disc = 1 - 4*CPsi
        if CPsi < 0.25:
            check(f"CPsi={CPsi}: disc={disc:.4f} > 0 (real fixed points)", disc > 0)
        elif CPsi == 0.25:
            check(f"CPsi=0.25: disc={disc:.4f} = 0 (boundary)", abs(disc) < 1e-10)
        else:
            check(f"CPsi={CPsi}: disc={disc:.4f} < 0 (complex, no real FP)", disc < 0)

    # Claim: "Mandelbrot equivalence: u = C(Psi+R) maps to z^2 + c with c = C*Psi"
    print("\n--- CORE_ALGEBRA: Mandelbrot equivalence ---")
    C, Psi = 0.8, 0.2
    R0 = 0.1
    # R iteration: R_{n+1} = C*(Psi + R_n)^2
    R1 = C * (Psi + R0)**2
    # u iteration: u = C(Psi+R), u_{n+1} = u_n^2 + c where c = C*Psi
    u0 = C * (Psi + R0)
    u1_direct = C * (Psi + R1)
    u1_mandelbrot = u0**2 + C*Psi
    check(f"Mandelbrot map: u1_direct={u1_direct:.8f} = u1_mandelbrot={u1_mandelbrot:.8f}",
          abs(u1_direct - u1_mandelbrot) < 1e-10,
          f"diff = {abs(u1_direct - u1_mandelbrot):.2e}")

    # =============================================================
    # WEAKNESSES claims (spot check)
    # =============================================================
    print("\n--- WEAKNESSES: J threshold: run verify_star_topology.py ---")

    # The star-topology threshold used to be "checked" here by passing the
    # literal True, with a note saying verify_star_topology.py had verified it.
    # A check that cannot fail is not a check, and the number it stood behind
    # (1.466) was itself a sampling artifact; the converged value is 1.46295
    # (F29). Run verify_star_topology.py for the gate that can fail.

    # =============================================================
    # CONSISTENCY CHECKS between documents
    # =============================================================
    print("\n--- CONSISTENCY: retired mixed-XY correlation ---")
    print("  The reported correlations described a non-invariant coordinate diagnostic.")

    # Claim: Echo peak C_SB = 0.598 for N=3
    # This would require a full star simulation which we've done in qst_bridge.py
    # Reference check only
    print("\n--- REFERENCE: Echo peak C_SB = 0.598 ---")
    print("  Verified in simulations/qst_bridge.py. Not recomputed here.")
    print("  Cross-reference: QST_BRIDGE.md and GLOSSARY.md both cite 0.598.")

    # Claim: F_avg = 0.888 for star 2:1
    print("\n--- REFERENCE: F_avg = 0.888 ---")
    print("  Verified in simulations/qst_bridge.py. Not recomputed here.")
    print("  Would require average over random input states (expensive).")
    print("  Cross-reference: QST_BRIDGE.md, GLOSSARY.md, WEAKNESSES.md all cite 0.888.")

    # =============================================================
    # SUMMARY
    # =============================================================
    print("\n" + "=" * 70)
    print(f"DOCS VERIFICATION SUMMARY: {PASS} PASS, {FAIL} FAIL")
    print("=" * 70)

    if FAIL > 0:
        print(f"\n{FAIL} FAILURES. Fix before any release.")
    else:
        print("\nAll checks implemented in this script passed.")
        print(SUCCESS_SCOPE)

    print("\nNOT RECOMPUTED (expensive, verified in dedicated scripts):")
    print("  - F_avg = 0.888 (qst_bridge.py)")
    print("  - Echo peak C_SB = 0.598 (qst_bridge.py)")
    print("  - Holevo capacity 0.534 bits (verify_channel.py)")
    print("  - IBM Torino results (hardware, not reproducible locally)")
    raise SystemExit(1 if FAIL else 0)
