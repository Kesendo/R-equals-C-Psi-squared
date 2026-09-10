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
import ast
import hashlib
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
            "They do not supply an all-N existence theorem",
            "Last refreshed 2026-09-10",
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
            "Path5_ResidualProposals_RequireIndependentFullPairCertification",
            "Path4_ProposalsAndDefectiveControl_RequireIndependentFullPairCertification",
            "Path5_TrackedNearSplit_IsUncertifiedRegardlessOfSearchTolerance",
        ),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzCoherenceBlockFieldTests.cs": (
            "These tests use real q=2",
            "diagonal addition -i*q*(fieldEnergy(ket)-fieldEnergy(bra))",
            "not eigenmode real parts", "matrix diagonal's imaginary component",
        ),
        "compute/RCPsiSquared.Diagnostics/Foundation/CrossFoldSimilarityWitness.cs": (
            "Delta proposals are Uncertified", "If a character is independently certified",
        ),
        "compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs": (
            "proposals are Uncertified", "if character is independently certified",
        ),
        "compute/RCPsiSquared.Cli/Commands/PathKMonodromyScanCommand.cs": (
            "finite-N Delta response", "conditional residual mechanism",
            "falsifies this sampled-locus Delta-death prediction",
            "not proof of causality or all-N protection",
            "bounded scan null", "not a replacement for the exact Route B inventory",
        ),
        "compute/RCPsiSquared.Core.Tests/F89PathK/WeightCoherenceBlockTests.cs": (
            "finite-N Delta response", "every sampled positive-Delta proposal at N=4/N=5/N=6 is Uncertified",
        ),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/CrossFoldSimilarityWitnessTests.cs": (
            "independently certified Jordan character", "without treating it as a character certificate",
        ),
        "compute/RCPsiSquared.Cli/Program.cs": (
            "finite-N Delta response", "conditional residual mechanism",
        ),
        "compute/RCPsiSquared.Diagnostics/Foundation/BranchLocusPalindromeWitness.cs": (
            "twin-scalar restriction at the N=4 point", "positive Delta proposals remain Uncertified", "γ=1, J=q_EP",
        ),
        "compute/RCPsiSquared.Core/Symmetry/F89BranchLocusPalindromeClaim.cs": (
            "twin-scalar restriction at the N=4 point", "positive Delta proposals remain Uncertified",
        ),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": (
            "twin-scalar restriction at the N=4 point", "N=5/N=6 positive Delta proposals remain Uncertified",
        ),
        "docs/ANALYTICAL_FORMULAS.md": (
            "At N=4 only the Delta=0 control is character-certified", "positive-Delta N=4 proposals are split and Uncertified",
            "general protection principle remains interpretive and conditional",
            "Loss of the W/H intertwining identity is an algebraic statement",
            "all sampled positive-Delta proposals are Uncertified",
        ),
        "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md": (
            "Twin-scalar lemma (conditional)",
            "Loss of the W/H intertwining identity is an algebraic statement",
            "every sampled positive-Delta proposal is Uncertified",
        ),
        "experiments/F89_BRANCH_LOCUS_PALINDROME.md": (
            "N=4 twin-scalar restriction", "sampled positive-Delta proposals remain split",
            "palindrome proof is independent of this unresolved local-character question",
        ),
        "hypotheses/DIABOLIC_BY_INTEGRABILITY.md": (
            "palindrome proof is independent", "sampled positive-Delta proposals are not certified coalescences",
            "mirror fixes the line and pairing", "N=4 overlap/twin-scalar restriction supplies the on-line position",
        ),
        "compute/RCPsiSquared.Diagnostics/Foundation/PathKMonodromyScout.cs": (
            "analytic-defective collisions can split linearly",
            "analytic-split candidate", "full-sector geometric multiplicity decides character",
        ),
        "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs": (
            "positive-Delta N=4/N=5/N=6 proposals are Uncertified after strict full-block checking",
            "the H-half of the twin-scalar condition only after single-multiplet descent is independently certified, while the D-half remains a locus-specific obligation",
            "MECHANISM BOUNDARY",
            "Free-fermion (Jordan-Wigner) integrability alone does not prove residual semisimplicity or a Delta response at every N",
            "exponent that reads the local eigenvalue split (~1 for a linear analytic split, ~0.5 for a sqrt branch)",
            "the exponent corroborates branch type but does not decide geometric multiplicity or semisimplicity",
            "proposal space, not an invariant space at Delta!=0",
            "ZZ couples it to removed directions",
            "first locate a genuine full-block degeneracy, then ask EpCharacter",
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
    wrong_arc = arc_source.replace("integrability alone does not prove residual semisimplicity",
                                   "integrability alone proves residual semisimplicity")
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
         current_arc.replace("proposal space, not an invariant space at Delta!=0", "globally invariant")),
        ("current OpenArcs removed full-block character caveat", current_arc, surfaces[arc_path],
         current_arc.replace("first locate a genuine full-block degeneracy, then ask EpCharacter", "The compressed block decides character")),
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
                                   or line.startswith("**Δ-tightness.") or line.startswith("**Delta-tightness.")))
        errors = defects(source, required)
        baseline = set(errors)
        check(f"Route B mechanism current truth: {relative_path}", not errors, ascii("; ".join(errors)))
        # Feed wrong prose through the same door: missing scope and restored
        # all-N claims must both fail. This guards the words, not the physics.
        check(f"Route B mechanism rejects missing scope: {relative_path}",
              bool(set(defects(source.casefold().replace(required[0].casefold(), ""), required)) - baseline))
        check(f"Route B mechanism rejects all-N promotion: {relative_path}",
              bool(set(defects(source + " Integrability governs existence at every N.", required)) - baseline))
        if relative_path == "experiments/F89_PATH_K_DIABOLIC.md":
            check("Route B mechanism rejects the half-hopping energy normalization",
                  bool(set(defects(source.replace("ε_k = 4J cos", "ε_k = 2J cos"), required)) - baseline))
        if relative_path.endswith("/BranchLocusPalindromeWitness.cs"):
            check("Route B mechanism rejects q_EP with contradictory J=gamma=1 label",
                  bool(set(defects(source.replace("γ=1, J=q_EP", "γ=J=1"), required)) - baseline))
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
        "random-field XY Hamiltonian remains quadratic",
        "does not classify the Z-dephasing Liouvillian",
        "random field plus Delta!=0 is the interacting disordered nonintegrable test",
        "finite executed CSR evidence",
    )
    surfaces = {
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzDeltaFlipTests.cs": (
            "base-label nearest-neighbour continuation", "proposal tracking, not an invariant-set proof"),
        "compute/RCPsiSquared.Core/F89PathK/WeightCoherenceBlock.cs": (
            "diagonal addition -i*q*(fe(ket)-fe(bra))", "for complex q it need not be purely imaginary",
            "unchanged dissipator entries do not imply unchanged eigenmode real parts"),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": distinction + (
            "finite-precision cluster representatives", "tolerance-dependent", "not an exact degeneracy census"),
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdWitness.cs": distinction + (
            "finite-precision cluster representatives", "tolerance-dependent", "not an exact degeneracy census"),
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdCsr.cs": distinction + (
            "finite-precision cluster representatives", "tolerance-dependent", "not an exact degeneracy census"),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/FillingThresholdCsrTests.cs": distinction,
        "compute/RCPsiSquared.Diagnostics/Foundation/IntegrabilityBreakingCsr.cs": (
            "recomputed eigenvalues by their current real parts", "not an invariant AT/residual decomposition",
            "unchanged dissipator does not fix eigenmode real parts",
            "finite-precision cluster representatives", "tolerance-dependent", "not an exact degeneracy census"),
        "experiments/FILLING_THRESHOLD_CHAOS.md": distinction + (
            "1e-9 finite-precision cluster representative", "tolerance-dependent", "not an exact degeneracy census"),
        "experiments/RANDOM_MATRIX_THEORY.md": distinction + (
            "nor proves integrability inside a block",
            "does not run backward to integrability"),
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
        "signature of an integrable / non-chaotic system",
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

    rmt = (root / "experiments/RANDOM_MATRIX_THEORY.md").read_text(encoding="utf-8")
    result5 = rmt.split("## Result 5:", 1)[1].split("## What This Does Not Answer", 1)[0]
    result5_required = (
        "1e-9 finite-precision cluster representative", "tolerance-dependent",
        "not an exact degeneracy census",
    )
    result5_errors = lambda text: [p for p in result5_required if p.casefold() not in text.casefold()]
    baseline = set(result5_errors(result5))
    check("Random-field RMT Result 5 declares CSR clustering scope", not baseline,
          ascii("; ".join(sorted(baseline))))
    for phrase in result5_required:
        check(f"Random-field RMT Result 5 requires {ascii(phrase)}",
              bool(set(result5_errors(result5.replace(phrase, ""))) - baseline))


def verify_poisson_inference_current_truth():
    """Poisson/no-ramp observations are not converse proofs of integrability."""
    root = Path(__file__).resolve().parents[1]
    surfaces = {
        "docs/outbound/SHIFTED_ORDER4_CHIRAL_SYMMETRY.md": (
            "does not prove integrability", "not by itself a class assignment",
            "uniform XXZ remains Bethe-integrable", "random-field XY Hamiltonian remains",
            "does not classify the Z-dephasing Liouvillian",
            "Delta!=0 plus generic disorder", "finite-N"),
        "hypotheses/DIABOLIC_BY_INTEGRABILITY.md": (
            "uniform XXZ chain remains Bethe-integrable", "does not measure a character flip",
            "Tier-2 hypothesis"),
        "docs/ANALYTICAL_FORMULAS.md": ("not a physical Heisenberg time", "do not establish a short/long-time boundary",),
        "docs/THE_INTERPRETATION.md": ("do not prove integrability", "remaining sectorwise algebra", "class label, is OPEN"),
        "docs/KMS_DETAILED_BALANCE.md": (
            "statistics do not remove the symmetry class or prove integrability", "random-matrix universality"),
        "experiments/TOPOLOGICAL_EDGE_MODES.md": ("do not prove integrability or topology",),
        "experiments/SPECTRAL_FORM_FACTOR.md": ("does not prove integrability",),
        "experiments/RANDOM_MATRIX_THEORY.md": ("does not assign the full irreducible-sector symmetry class",),
        "compute/RCPsiSquared.Diagnostics/Foundation/GaloisSpectralChaosWitness.cs": (
            "structured Delta=0 XY frequency lattice", "does not classify the Liouvillian as integrable"),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/GaloisSpectralChaosWitnessTests.cs": (
            "structured Delta=0 XY frequency lattice", "does not classify the Liouvillian as integrable"),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": (
            "structured Delta=0 XY frequency lattice", "does not classify the Liouvillian as integrable"),
    }
    forbidden = ("Poisson (integrable)", "Poisson, integrable", "SFF confirms integrability",
                 "Poisson statistics of the integrable system", "precisely because it is the integrable member",
                 "single integrability-breaking knob", "outside all 38", "does not fit any of the 38",
                 "integrable chiral Lindbladian", "Liouvillian is integrable",
                 "Poisson under integrability breaking", "statistics of none of the 38",
                 "integrable at every N", "integrable frequency lattice", "reads integrable-lattice")
    for path, required in surfaces.items():
        source = (root / path).read_text(encoding="utf-8")
        if path == "docs/ANALYTICAL_FORMULAS.md":
            source = source.split("### F42.", 1)[1].split("### F43.", 1)[0]
        source = " ".join(source.replace("**", "").replace("*", "").split())

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
        "experiments/README.md": crossover + ("No standard-ensemble calibration for this degenerate unresolved population",),
        "experiments/F89_MULTI_SECTOR_MONODROMY.md": crossover,
        "experiments/RANDOM_MATRIX_THEORY.md": crossover,
        "experiments/FILLING_THRESHOLD_CHAOS.md": crossover,
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdWitness.cs": crossover + (
            "finite-precision cluster representatives", "tolerance-dependent", "not an exact degeneracy census"),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/IntegrabilityBreakingCsrTests.cs": (
            "random-field XY Hamiltonian remains Anderson/free-fermion",
            "without making the dephasing Liouvillian a quadratic generator",
            "random field plus Delta!=0 is interacting disordered nonintegrable"),
        "compute/RCPsiSquared.Core/F89PathK/WeightCoherenceBlock.cs": (
            "-i*q*Delta*(zz(ket)-zz(bra))", "for Im q != 0 it has a real component"),
        "simulations/spectral_form_factor.py": ("does not prove integrability", "not a spectral-class or no-ramp test"),
        "simulations/results/spectral_form_factor.txt": ("does not prove integrability", "not a spectral-class or no-ramp test"),
        "simulations/rmt_topology_csr.py": (
            "does not prove integrability", "simultaneous Hamiltonian/dissipator eigenoperators",
            "not a formula for mixed Liouvillian modes", "CLUSTER_DECIMALS = 9",
            "finite-precision clusters", "tolerance-dependent", "cluster_decimals=",
            "uncorrelated reference", "not an integrability verdict"),
    }
    forbidden = ("Liouvillian is integrable", "Poisson (integrable)", "is INTEGRABLE. No dip",
                 "-> Heisenberg integrability", "is chaotic", "chaos located as a FILLING threshold",
                 "dissipative chaos is a filling threshold", "turns on with filling",
                 "at the GinUE value already", "DENSE block IS chaotic",
                 "random field U[−W,W] breaks integrability",
                 "adds only the frequency", "Poisson-like SFF at", "GUE-like SFF at",
                 "Classification: POISSON-like", "Classification: GUE-like",
                 "Classification: intermediate / palindromic", "physically correct diagnostic",
                 "Exact degeneracies are removed first", "DISTINCT Liouvillian eigenvalues",
                 "2D Poisson (integrable / symmetry-fragmented)", "(integrable/fragmented)")
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
        match = re.search(r"complete N=6: clustered_upper=(\d+) total_upper=(\d+) clustered_fraction=([\d.eE+-]+)", text)
        table = re.search(r"^complete\s+6\s+\d+\s+(\d+)\s", text, re.MULTILINE)
        if match is None or table is None:
            return False
        distinct, reported_upper, collapse = int(match[1]), int(match[2]), float(match[3])
        return (reported_upper == upper and distinct == int(table[1])
                and collapse == 1 - distinct / reported_upper)

    check("Topology collapse uses the executed upper-half population", valid_collapse(output))
    wrong_denominator = re.sub(r"total_upper=\d+", f"total_upper={len(complete)}", output)
    reported_clustered = re.search(r"clustered_upper=(\d+)", output)
    if reported_clustered is not None:
        # Keep the wrong fraction internally consistent: only the population
        # boundary, not an arithmetic mismatch, must reject all-spectrum counting.
        wrong_denominator = re.sub(r"clustered_fraction=[\d.eE+-]+",
                                   f"clustered_fraction={1 - int(reported_clustered[1]) / len(complete):.17g}",
                                   wrong_denominator)
    check("Topology collapse rejects all-spectrum denominator", valid_collapse(output) and not valid_collapse(wrong_denominator))

    declared_digits = re.search(r"cluster_decimals=(\d+)", output)
    declared_sensitivity = re.search(
        r"cluster_sensitivity_N6_complete="
        r"d7:(\d+),d8:(\d+),d9:(\d+),d10:(\d+)", output)
    complete_complex = complete[:, 0] + 1j * complete[:, 1]
    complete_upper = complete_complex[complete_complex.imag > 1e-6]
    expected_sensitivity = tuple(
        len(np.unique(np.round(complete_upper, digits))) for digits in range(7, 11))
    check("Topology CSR declares the protected 1e-9 clustering tolerance",
          declared_digits is not None and int(declared_digits[1]) == 9)
    check("Topology CSR exposes clustering sensitivity from the executed population",
          declared_sensitivity is not None
          and tuple(map(int, declared_sensitivity.groups())) == expected_sensitivity)

    shared_csr_surfaces = {
        "compute/RCPsiSquared.Core/Numerics/ComplexSpacingRatio.cs": (
            "ClusterDecimals = 9", "finite-precision cluster", "tolerance-dependent",
            "not an exact degeneracy", "canonically ordered rounded",
            "OrderBy(key => key.Item1).ThenBy(key => key.Item2)",
            "new Complex(key.Item1 / scale, key.Item2 / scale)"),
        "compute/RCPsiSquared.Diagnostics/Foundation/GaloisSpectralChaosWitness.cs": (
            "finite-precision cluster representatives", "tolerance-dependent",
            "not an exact degeneracy census"),
    }
    shared_forbidden = ("distinct spectrum", "distinct-point count", "distinct points/q",
                        "pts.Add(p)")
    for path, required in shared_csr_surfaces.items():
        source = (root / path).read_text(encoding="utf-8")

        def shared_errors(text):
            normalized = " ".join(text.replace("///", "").replace('" +', "")
                                  .replace('"', "").split()).casefold()
            return ([f"missing {phrase}" for phrase in required if phrase.casefold() not in normalized]
                    + [f"forbidden {phrase}" for phrase in shared_forbidden
                       if phrase.casefold() in normalized])

        baseline = set(shared_errors(source))
        check(f"Shared CSR clustering scope: {path}", not baseline,
              ascii("; ".join(sorted(baseline))))
        for phrase in required:
            check(f"Shared CSR clustering requires {ascii(phrase)}: {path}",
                  bool(set(shared_errors(source.casefold().replace(phrase.casefold(), ""))) - baseline))
        for phrase in shared_forbidden:
            check(f"Shared CSR clustering rejects {ascii(phrase)}: {path}",
                  bool(set(shared_errors(source + " " + phrase)) - baseline))


def verify_sff_windows_and_current_boundaries():
    root = Path(__file__).resolve().parents[1]
    # The dispersion's Taylor coefficient is computed independently of prose:
    # 4J(1-cos x) = J * coefficient * x^2 + O(x^4).
    import sympy as sp
    angle = sp.Symbol("angle")
    coefficient = sp.diff(4 * (1 - sp.cos(angle)), angle, 2).subs(angle, 0) / 2
    sff_doc = (root / "experiments/SPECTRAL_FORM_FACTOR.md").read_text(encoding="utf-8")

    def asymptotic_matches_dispersion(text):
        matches = re.findall(r"ω_F2\s*~\s*(\d*)J?π²/N²", text)
        return len(matches) == 1 and sp.Integer(matches[0] or "1") == coefficient

    check("SFF asymptotic coefficient equals exact dispersion Taylor coefficient",
          asymptotic_matches_dispersion(sff_doc))
    check("SFF asymptotic rejects half coefficient with J",
          not asymptotic_matches_dispersion("ω_F2 ~ Jπ²/N²"))
    check("SFF asymptotic rejects half coefficient at J=1",
          not asymptotic_matches_dispersion("ω_F2 ~ π²/N²"))
    check("SFF asymptotic rejects missing coefficient formula",
          not asymptotic_matches_dispersion("t_Π grows as N²"))
    surfaces = {
        "hypotheses/DIABOLIC_BY_INTEGRABILITY.md": (
            "N=4 twin-scalar restriction", "Tier-2 residual interpretation",
            "uniform XXZ chain remains Bethe-integrable", "N=6 inventory is complete", "N=7 remains the inventory boundary"),
        "experiments/SPECTRAL_FORM_FACTOR.md": (
            "raw non-unfolded frequency SFF", "1/M", "not sampled", "reached window",
            "not a universality-class proof", "only this raw oscillation-frequency SFF",
            "physical Liouvillian propagation would involve", "requires a specified ensemble or averaging prescription",
            "raw multiset density scale", "multiplicity-dependent", "N=4 third", "N=6 second", "N=7 has no match", "114%",
            "(0,1)-block k=1 reference", "not the smallest nonzero frequency of the full spectrum"),
        "experiments/README.md": ("raw non-unfolded frequency SFF", "reached window", "not a universality-class proof", "raw multiset density scale"),
        "experiments/RANDOM_MATRIX_THEORY.md": ("cannot classify universality",),
        "docs/outbound/SHIFTED_ORDER4_CHIRAL_SYMMETRY.md": ("measured finite-N contrast", "convergence is a proposed extension"),
        "simulations/spectral_form_factor.py": ("independent_phase_reference", "format_sample", "raw non-unfolded frequency SFF",
                                                "2*cos", "doubled and cross frequencies",
                                                "raw multiset density scale", "multiplicity-dependent", "RCPSI_SFF_OUTPUT_PATH", "raw_multiset_density_scale(ev)",
                                                "Only the raw oscillation-frequency SFF is computed", "no decay-weighted or connected estimator is defined",
                                                "(0,1)-block k=1 reference", "not asserted to be the smallest nonzero frequency of the full spectrum"),
        "simulations/sff_window_summary.py": ("def sff_frequency", "def raw_multiset_density_scale", "multiplicity-dependent"),
        "simulations/results/spectral_form_factor.txt": (
            "independent-phase reference 1/M", "not sampled", "reached window", "not classifiable",
            "2*cos", "doubled and cross frequencies",
            "raw multiset density scale", "multiplicity-dependent",
            "Only the raw oscillation-frequency SFF is computed", "no decay-weighted or connected estimator is defined",
            "(0,1)-block k=1 reference", "not asserted to be the smallest nonzero frequency of the full spectrum"),
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
                 "dominant FFT peak matches", "ω_min peak is present", "first sampled exceedance",
                 "ω_min (slowest)", "slowest pair's trace-amplitude term",
                 "def sff_dissipative", "def sff_connected", "K_diss", "overflow/NaN")
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


def task8_text_errors(source, required, forbidden):
    """One predicate for live current text and deliberately corrupted current text."""
    text = " ".join(source.replace("**", "").replace("///", "")
                    .replace('" +', "").replace('"', "").split()).casefold()
    return ([f"missing {p}" for p in required if p.casefold() not in text]
            + [f"forbidden {p}" for p in forbidden if p.casefold() in text])


def f44_ordered_pair_api_errors(source):
    """Pin both public F44 evaluators to the same finite ordered-pair domain."""
    errors = []
    if source.count("ValidateOrderedPairDomain(deltaD, totalGamma);") != 2:
        errors.append("exact and linear F44 APIs must both use shared domain validation")
    for token in ("!double.IsFinite(totalGamma)", "totalGamma <= 0.0",
                  "!double.IsFinite(deltaD)", "deltaD < 0.0",
                  "deltaD / totalGamma >= SumCoefficient"):
        if token not in source:
            errors.append(f"F44 domain guard missing {token}")
    return errors


def mirrorworld_f44_ordered_pair_api_errors(source):
    errors = []
    for token in ("ValidateF44PositivePalindromicPair(dFast, dSlow, sg);",
                  "!double.IsFinite(sg)", "sg <= 0.0",
                  "!double.IsFinite(dFast)", "dFast <= 0.0",
                  "!double.IsFinite(dSlow)", "dSlow <= 0.0",
                  "Math.Abs(dFast / scale + dSlow / scale - 2.0 * (sg / scale))",
                  "!double.IsFinite(argument)", "argument >= 1.0"):
        if token not in source:
            errors.append(f"MirrorWorld F44 domain guard missing {token}")
    return errors


def f95_positive_b_api_errors(source, mirrorworld=False):
    guard = "if (!double.IsFinite(b) || b <= 0.0)"
    return [] if guard in source else [
        f"{'MirrorWorld ' if mirrorworld else ''}F95 must reject non-finite or non-positive b"]


def f95_python_positive_b_api_errors(source):
    guard = "if not math.isfinite(b) or b <= 0:"
    return [] if guard in source else ["F95 Python producer must reject non-finite or non-positive b"]


def delta_tracker_path_errors(source):
    # Select declarations, not matching call text somewhere else in the file.
    # Comments/strings are masked before brace matching; offsets stay unchanged.
    masked = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',
                    lambda m: " " * len(m.group()), source)
    def method(name):
        declaration = re.search(r"\b(?:public|private|internal)\s+static\s+[^\n(]+\s+" + name + r"\s*\(", masked)
        if not declaration:
            return ""
        start = declaration.start()
        brace, arrow = masked.find("{", start), masked.find("=>", start)
        if arrow >= 0 and (brace < 0 or arrow < brace):
            return source[start:masked.index(";", arrow) + 1]
        depth = 0
        for end in range(brace, len(masked)):
            depth += (masked[end] == "{") - (masked[end] == "}")
            if depth == 0:
                return source[start:end + 1]
        return ""
    seeded = "return CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, mid, lambdaSeed, seedRadius);"
    contracts = {
        "CharacterAtDiabolicNear": ("public static DeltaTrackResult CharacterAtDiabolicNear",
            "=> TrackDiabolicUnderDelta(n, qSeed, lambdaSeed, delta, boxHalf: 0, boxCell: cell);"),
        "TrackDiabolicUnderDeltaCompressed": (seeded,),
        "TrackDiabolicUnderDeltaResidual": (
            "return CertifyFullBlockProposal(BuildSym(n, center, delta), center, mid, lambdaSeed, seedRadius);",),
        "TrackDiabolicUnderDelta": (seeded,
            "TrackDiabolicUnderDeltaCompressed(n, qSeed, lambdaSeed, seedRadius, delta, boxHalf, boxCell)",
            "if (n >= 7 && !exact)", "if (!(seedDistances[1] < seedRadius))"),
        "CertifyFullBlockProposal": ("fullGap > FullBlockCoincidenceTolerance",
            "lambdaSeed.HasValue && !((mid - lambdaSeed.Value).Magnitude < seedRadius)",
            "(full[0] - mid).Magnitude <= FullBlockCoincidenceTolerance",
            "(full[1] - mid).Magnitude <= FullBlockCoincidenceTolerance",
            "if (!(radius > 2 * ds[1]))",
            "rr.Algebraic == 2 ? CertifiedCharacterVerdict(rr.Kind) : DeltaFlipVerdict.Uncertified"),
    }
    return [f"{name}: {error}" for name, required in contracts.items()
            for error in task8_text_errors(method(name), required, ())]


def verify_task8_current_truth():
    root = Path(__file__).resolve().parents[1]
    surfaces = {
        "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md": (("ε_k = 4J·cos(πk/(N+1))",), ("ε_k = 2J·cos(πk/(N+1))",)),
        "compute/RCPsiSquared.Core/F86/JordanWigner/XyJordanWignerModes.cs": ((
            "H_XY = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1})",
            "h[i, i±1] = J", "ε_k = 2·J·cos"),
            ("H_XY = J·Σ_b (X_b X_{b+1} + Y_b Y_{b+1})",)),
        "docs/proofs/PROOF_F90_F86C2_BRIDGE.md": (("ε_k = 2J·cos(πk/(N+1))", "ε_k = 4J·cos(πk/(N+1))"), ("ε_k = J·cos", "from L_super")),
        "experiments/F89_BRANCH_LOCUS_PALINDROME.md": (("q_octic = Q_carrier/2", "Q_carrier = 1.5 corresponds to q_octic = 0.75"), ("On the real axis q = Q", "q and Q (same ratio, different role)")),
        "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md": (("single-multiplet descent remains conditional", "numerical twin-scalar character"), ("additivity's semisimplicity route extends verbatim", "no unexplained residual coincidence remains at N=5")),
        "experiments/F89_PATH_K_DIABOLIC.md": (("AT semisimple crossings are diabolics, not residual A2 diabolics", "3/16", "Uncertified", "full-pair gap", "nonzero-Delta response remains open"),
                                                  ("These are not diabolics but", "each path-5 diabolic reads cleanly", "three real-q diabolics become defective and one lifts", "ResidualRootsExactXxz",
                                                   "N=4, N=5, N=6 and N=7 diabolics flip defective or lift")),
        "docs/proofs/derivations/D09_SECTOR_SFF_PAIRING.md": (("palindrome-paired decay-rate bands", "average-light bins", "constant (N+1)²", "normalized SFF is 1"), ("restricted to XY-weight sector", "(N+1)² × δ(t)", "stationary + near-stationary")),
        "experiments/RANDOM_MATRIX_THEORY.md": (("executed N=6..8 fixed-parameter comparison",), ("question is resolved in Result 5",)),
        "compute/RCPsiSquared.Diagnostics/Foundation/PathKMonodromyScout.cs": (("ResidualRootsCompressedXxz", "proposal roots"), ("ResidualRootsExactXxz",)),
        "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs": ((
            "FullBlockCoincidenceTolerance = 1e-6", "fullGap > FullBlockCoincidenceTolerance",
            "(full[0] - mid).Magnitude <= FullBlockCoincidenceTolerance",
            "(full[1] - mid).Magnitude <= FullBlockCoincidenceTolerance",
            "TrackDiabolicUnderDeltaCompressed(n, qSeed, lambdaSeed, seedRadius, delta, boxHalf, boxCell)",
            "CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, mid, lambdaSeed, seedRadius)",
            "CertifyFullBlockProposal(BuildSym(n, center, delta), center, mid, lambdaSeed, seedRadius)",
            "if (n >= 7 && !exact) return new DeltaTrackResult(DeltaFlipVerdict.Uncertified",
            "lambdaSeed.HasValue && !((mid - lambdaSeed.Value).Magnitude < seedRadius)",
            "if (!(seedDistances[1] < seedRadius))",
            "EpCharacter.EpKind.Diabolic => DeltaFlipVerdict.Diabolic",
            "EpCharacter.EpKind.Defective => DeltaFlipVerdict.Defective",
            "_ => DeltaFlipVerdict.Uncertified", "rr.Algebraic == 2",
            "public bool IsCertifiedDiabolic", "public bool? Survived", "_ => null",
            "Complex QCandidate, Complex LambdaCandidate, double Gap",
            "It is neither survival nor death", "unavailable character, not measured multiplicities"),
            ("ResidualRootsExactXxz", "fullGap > coalesceTol",
             "return CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, mid);",
             "public bool Survived => Verdict == DeltaFlipVerdict.Diabolic;",
             "EpKind.Normal => DeltaFlipVerdict.Defective", "EpKind.Normal => DeltaFlipVerdict.Diabolic",
             "EpKind.NearEp => DeltaFlipVerdict.Defective", "EpKind.NearEp => DeltaFlipVerdict.Diabolic",
             "DeltaFlipVerdict.Uncertified => false", "CertifyCompressedProposal(",
             "Complex QStar, Complex LambdaStar, double Gap",
             "new DeltaTrackResult(DeltaFlipVerdict.Lifted")),
        "compute/RCPsiSquared.Cli/Commands/PathKMonodromyScanCommand.cs": ((
            "Uncertified", "full-block pair", "t.IsCertifiedDiabolic", "UNRESOLVED",
            "no defect/lift conclusion", "N/A N/A",
            "Only the N=4, Delta=0 control is character-certified",
            "all sampled positive-Delta proposals at N=4..7 are Uncertified",
            "track the candidate pair from q_seed=", "q_candidate", "t.QCandidate"),
            ("ResidualRootsExactXxz", "exact-residual: AT-flood excluded",
             "if (d > 0 && t.Survived) survivesAtPositive = true;",
             "Only N=4 is character-certified", "track the coalescence at q=", "t.QStar")),
        "compute/RCPsiSquared.Cli/Program.cs": (("Uncertified", "compressed proposals", "unknown survival", "independent 1e-6"), ()),
        "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs": ((
            "the H-half of the twin-scalar condition only after single-multiplet descent is independently certified",
            "the D-half remains a locus-specific obligation", "Moves 1, 3 and 4 are done",
            "Move 2 remains open at N=7 after the compression-locator retraction",
            "their nonzero-Delta response is OPEN", "not an invariant space at Delta!=0",
            "first locate a genuine full-block degeneracy",
            "one basis coherence with disagreement n_diff pays the positive cell rate",
            "an AT-locked eigenmode has Re lambda = -2*gamma*<n_XY>_v",
            "<n_XY>_v need not equal the integer n_diff"),
            ("additivity supplies the H-half unconditionally", "Moves 1, 2, 3 and 4 are done",
             "N=7 nonzero-Delta response is closed", "an AT eigenmode has the integer basis-cell disagreement count")),
        "compute/RCPsiSquared.Core/Symmetry/F41PalindromicTimePi2Inheritance.cs": ((
            "2cos(ω_min t)", "trace amplitude", "doubled and cross frequencies",
            "F41 is not by itself the period of the complete SFF",
            "finite-N descriptive ratio, not a physical short/long-time division"),
            ("physical Heisenberg time t_H ~ 4^N", "F41 is the period of the complete SFF",
             "t_Pi/t_H ~ N²/4^N", "slowest palindromic SFF modulation frequency")),
        "compute/RCPsiSquared.Core/Symmetry/F43BandSffPairingPi2Inheritance.cs": ((
            "reflected decay-rate bands", "average-light bands", "not assigned to an invariant fixed-integer-weight sector",
            "N+1 stationary modes", "frequency SFF is the constant 1 for every time",
            "not an impulse at time zero", "finite-width numerical band can contain more"),
            ("fixed XY-weight eigenvalue sectors", "normalized endpoint K(t) = delta(t)", "endpoint multiplicity 2^N",
             "canonical Pauli-string weight partition", "delta-spike at zero frequency")),
        "docs/ANALYTICAL_FORMULAS.md": ((
            "doubled and cross frequencies", "not a physical Heisenberg time",
            "reflected decay-rate bands", "normalized frequency SFF is the constant 1",
            "N+1 stationary modes", "not a delta spike"),
            ("physical Heisenberg time t_H ~ 4^N", "normalized endpoint K(t) = delta(t)", "endpoint multiplicity 2^N")),
        "experiments/DEGENERACY_HUNT.md": (("Re(λ) = −2γ·⟨n_XY⟩_v", "Only a basis coherence has the integer count", "neither equal"),
            ("every eigenmode has integer n_diff",)),
        "experiments/SYMMETRY_CENSUS.md": (("Re(λ) = −2γ·⟨n_XY⟩_v", "not necessarily one integer-weight basis sector",
            "not invariant fixed-XY-weight eigenspaces"), ("every eigenmode has integer n_diff",)),
        "docs/outbound/SHIFTED_ORDER4_CHIRAL_SYMMETRY.md": ((
            "Order 4 globally, involutive sectorwise", "P_px=√px·Π|px", "P_px²=I",
            "Spectral reflection alone cannot select AIII or any other class",
            "cannot confirm or refute a sector-level universality prediction",
            "Prior-art coverage and equivalence to published constructions remain OPEN",
            "No priority or novelty claim is made"),
            ("the global Liouvillian is class AIII", "the constant shift is outside SRP",
             "Priority bookkeeping", "does not state the operator or the spectral reflection",
             "The repository-specific object here")),
        "docs/KMS_DETAILED_BALANCE.md": ((
            "No global AIII, BDI, or CI label is assigned here", "irreducible-sector classification remains open",
            "SRP applies negative symmetries to the trace-shifted generator",
            "Prior-art coverage and equivalence for that interacting, locality-resolved scope remain OPEN",
            "not an ownership, independence, priority, or novelty claim",
            "Strict GNS/KMS symmetry of the full generator",
            "Alicki/standard QDB permits a separate Hamiltonian derivation",
            "only the dissipative part is self-adjoint",
            "does not force the full generator to have a real spectrum",
            "For a QDB dissipative semigroup with no unitary part",
            "equal to its Petz recovery map",
            "With an additional commuting unitary part",
            "reverses the unitary sign while retaining the same dissipative evolution"),
            ("the global Liouvillian is class AIII", "the constant shift is outside SRP",
             "Our Π does not fit cleanly into any of their 38 classes", "constant shift 2Sγ is not part of their formalism",
             "QDB relates L to its adjoint L† (real eigenvalues)",
             "What QDB implies spectrally: all eigenvalues of L are real",
             "what stays ours", "are the new part",
             "new variant of known framework", "We built this", "with no literature input",
             "reached independently", "recognized afterward", "Ours works",
             "showed that quantum detailed balance is equivalent to the Petz recovery map",
             "QDB = Petz recovery map being exact channel reversal")),
        "simulations/universal_carrier_demo.py": ((
            "fixed Q=J/gamma=20", "J=2.0, gamma_0=gamma_10",
            "γ-only sweep at fixed J", "not this scaling identity",
            "joint J,γ scaling"),
            ("J=1.0, gamma_0=gamma_10", "abs(closest - pred) < 1e-4",
             "What scales: γ alone", "dimensionless n_XY values: γ-invariant")),
        "docs/THE_INTERPRETATION.md": (("sqrt(p_x) Π", "remaining sectorwise algebra", "class label, is OPEN",
            "do not prove integrability"), ("the global Liouvillian is class AIII",)),
        "experiments/README.md": (("full irreducible-sector SRP class is OPEN",
            "involutive sectorwise P generator", "No standard-ensemble calibration for this degenerate unresolved population"),
            ("the global Liouvillian is class AIII",)),
        "compute/RCPsiSquared.Core/F86/ShiftedGeneratorSectorwisePClaim.cs": ((
            "public sealed class ShiftedGeneratorSectorwisePClaim : Claim",
            "public sealed class FullIrreducibleSrpClassQuestion : Claim", "Tier.OpenQuestion",
            "P_px=sqrt(p_x)Π obeys P_px²=I", "does not assign the full irreducible SRP class",
            "sectorwise P-type anticommutation alone does not fix an SRP class"),
            ("ChiralAiiiClassification", "global SRP class AIII is established")),
        "compute/RCPsiSquared.Runtime/F86Main/F86MainRegistration.cs": ((
            ".Register<ShiftedGeneratorSectorwisePClaim>(b => { _ = b.Get<F1PalindromeIdentity>(); return new ShiftedGeneratorSectorwisePClaim(); })",),
            ("ChiralAiiiClassification",)),
        "compute/RCPsiSquared.Runtime/F86Main/F86OpenQuestionsRegistration.cs": ((
            ".Register<FullIrreducibleSrpClassQuestion>(b => { _ = b.Get<ShiftedGeneratorSectorwisePClaim>(); return new FullIrreducibleSrpClassQuestion(); })",),
            ("ChiralAiiiClassification",)),
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdWitness.cs": ((
            "GinUE comparison only: the full sector symmetry algebra and irreducible SRP class remain open",
            "GinUE is a comparison ensemble only; other unitary or antiunitary relations have not been excluded",
            "not a causal or thermodynamic threshold theorem"),
            ("Class A licensed", "no residual antiunitary", "GinUE licenses the full symmetry class")),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzDeltaFlipTests.cs": ((
            "Path5_ResidualProposals_RequireIndependentFullPairCertification",
            "Path4_ProposalsAndDefectiveControl_RequireIndependentFullPairCertification",
            "Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, d0.Verdict)",
            "AssertUncertifiedSplit(6, .1, d)", "AssertUncertifiedSplit(5, .05, d)",
            "Assert.Equal(0, result.Algebraic)", "Assert.Null(result.Survived)"),
            ("Assert.False(d.Survived", "Path5_Diabolics_DieUnderDelta_ResidualOnly",
             "Path4_Diabolics_DieUnderDelta_ControlStaysPut")),
        "compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs": ((
            "SingleDisagreementCellCost(double gammaZero)", "BasisPairDissipatorCost(int nDiff, double gammaZero)",
            "EigenmodeDecayRate(double averageNXy, double gammaZero)",
            "EigenmodeDecayRateCeiling(int n, double gammaZero)",
            "AverageNXyFromEigenmodeDecayRate(double rate, double gammaZero)",
            "HammingComplementCellCostSum(int blockSize, double gammaZero)",
            "return DissipatorCoefficient * gammaZero * averageNXy;",
            "generally non-integer expectation", "not an eigenvalue step"),
            ("public double Rate(", "AbsorptionQuantum", "PerCoherenceRateComputationalBasis(",
             "public double MaxRate(", "public double NXyFromRate(", "HammingComplementPairSum(",
             "eigenmode decay rates are integer multiples of 2γ")),
    }
    required, forbidden = surfaces["docs/KMS_DETAILED_BALANCE.md"]
    # Round-3 connected current consumers: the same predicate owns their live
    # contract and each removed-boundary / restored-old-claim input mutation.
    surfaces.update({
        "compute/RCPsiSquared.Core/Symmetry/PiOperator.cs": ((
            "P_px = sqrt(px)·Π", "sectorwise P-type generator; the full irreducible SRP class remains open"),
            ("existing class-AIII-chiral classification",)),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": ((
            "N=5/N=6 positive Delta proposals remain Uncertified",
            "Sampled small gaps beyond N=4 are similarity reads, not defect, lift, or diabolic certificates",
            "GinUE comparison only", "irreducible SRP class, which remain open"),
            ("N=5/N=6 positive Delta proposals defect or lift", "FillingThresholdWitness class-A gate")),
        "compute/RCPsiSquared.Diagnostics/Foundation/CrossFoldSimilarityWitness.cs": ((
            "Delta proposals are Uncertified", "If a character is independently certified",
            "coincidence and Jordan character remain Uncertified for this proposal"),
            ("N=5/N=6 positive Delta proposals defect or lift", "sampled small gaps certify Jordan character")),
        "compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs": ((
            "N=5/N=6 positive Delta proposals are Uncertified",
            "Its character transport is conditional on an independent coincidence/Jordan",
            "The similarity itself does not turn a sampled small gap into that certificate"),
            ("N=5/N=6 positive Delta proposals defect or lift", "sampled small gaps certify Jordan character")),
        "compute/RCPsiSquared.Core.Tests/F89PathK/WeightCoherenceBlockTests.cs": ((
            "Only the N=4, Delta=0 control is", "every sampled positive-Delta proposal at N=4/N=5/N=6 is Uncertified"),
            ("Only N=4 is character-certified", "N=5/N=6 positive Delta proposals defect or lift")),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/CrossFoldSimilarityWitnessTests.cs": ((
            "without treating it as a character certificate", "Summary_LeavesN5N6PositiveDeltaCharacterUncertified"),
            ("sampled small gaps certify Jordan character",)),
        "simulations/pt_symmetry_analysis.py": ((
            "P_px=sqrt(p_x)·Π is involutive", "The full irreducible SRP class remains OPEN",
            "not a global class assignment", "within palindrome axis",
            "off palindrome axis; operator relation remains exact"),
            ("Class AIII", "chiral phase", "chiral breaking")),
        "compute/RCPsiSquared.Compute/README.md": ((
            "sectorwise P symmetry after resolving Π²", "full irreducible symmetry class remains OPEN"),
            ("chiral symmetry class AIII", "ChiralAiiiClassification")),
        "compute/RCPsiSquared.Core/F_FORMULA_CROSSWALK.md": ((
            "AbsorptionTheoremClaim.SingleDisagreementCellCost(γ₀)",
            "AbsorptionTheoremClaim.BasisPairDissipatorCost(nDiff, γ₀)",
            "AbsorptionTheoremClaim.EigenmodeDecayRate(averageNXy, γ₀)",
            "AbsorptionTheoremClaim.EigenmodeDecayRateCeiling(N, γ₀)",
            "AbsorptionTheoremClaim.AverageNXyFromEigenmodeDecayRate(rate, γ₀)",
            "AbsorptionTheoremClaim.HammingComplementCellCostSum(blockSize, γ₀)",
            "generally non-integer expectation", "not thereby quantized on a 2γ₀ grid"),
            ("AbsorptionTheoremClaim.AbsorptionQuantum", "AbsorptionTheoremClaim.Rate(",
             "AbsorptionTheoremClaim.MaxRate(", "AbsorptionTheoremClaim.NXyFromRate(",
             "AbsorptionTheoremClaim.PerCoherenceRateComputationalBasis(", "ChiralAiiiClassification")),
        "experiments/FILLING_THRESHOLD_CHAOS.md": ((
            "not a symmetry-class assignment", "does not exhaust possible residual antiunitary symmetries",
            "every strong/unitary sector, is OPEN"), ("FillingThresholdWitness class-A gate",)),
        "experiments/CONCENTRATOR_MAPPING.md": ((
            "FillingThresholdWitness` conjugation-match diagnostic", "alongside its GinUE comparison",
            "That diagnostic does not assign a symmetry class"), ("FillingThresholdWitness class-A gate",)),
        "hypotheses/ZERO_IS_THE_MIRROR.md": ((
            "generally non-integer expectation", "2γ is a cell-cost quantum, not a proved spacing of the interacting rate spectrum"),
            ("makes γ the rung spacing of the rate spectrum",)),
    })
    path = "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md"
    proof_required, proof_forbidden = surfaces[path]
    surfaces[path] = (proof_required + (
        "algebraic/geometric multiplicities and departure are unavailable, not measured character data",
        "The proposals therefore establish neither persistence, defectiveness, nor lifting"),
        proof_forbidden + ("Uncertified proposals supply measured departure",))
    path = "docs/THE_INTERPRETATION.md"
    interpretation_required, interpretation_forbidden = surfaces[path]
    surfaces[path] = (interpretation_required + (
        "It does not decide priority or novelty", "the structure is new rather than a rediscovery, remain OPEN"),
        interpretation_forbidden + ("the structure is new rather than a rediscovery is established",))
    surfaces["docs/KMS_DETAILED_BALANCE.md"] = (required + (
        "CONFIRMED algebra / OPEN class", "complete irreducible-sector SRP algebra is not yet computed"), forbidden)
    for path in ("experiments/F89_PATH_K_DIABOLIC.md", "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs"):
        required, forbidden = surfaces[path]
        surfaces[path] = (required, forbidden + ("Every positive-Delta N=5/N=6 proposal is Defective",))
    path = "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs"
    required, forbidden = surfaces[path]
    surfaces[path] = (required + (
        "while their Delta=0 controls remain Diabolic with alg=geo=2",), forbidden)
    path = "experiments/F89_PATH_K_DIABOLIC.md"
    required, forbidden = surfaces[path]
    surfaces[path] = (required + (
        "The Delta=0 N=4/N=5 controls are Diabolic (alg=geo=2)",
        "same holds for the N=6 seeds before perturbation",
        "fixed 1e-6 full-block certificate",
        "Every sampled positive-Delta N=4/N=5/N=6 proposal is Uncertified"),
        forbidden + ("full-pair gap exceeds coalesceTol=0.001",))
    required, forbidden = surfaces["docs/proofs/derivations/D09_SECTOR_SFF_PAIRING.md"]
    surfaces["docs/proofs/derivations/D09_SECTOR_SFF_PAIRING.md"] = (required, forbidden + (
        "fixed XY-weight eigenvalue sectors", "normalized endpoint K(t) = delta(t)"))
    required, forbidden = surfaces["experiments/RANDOM_MATRIX_THEORY.md"]
    surfaces["experiments/RANDOM_MATRIX_THEORY.md"] = (required + (
        "Finite-size computational study", "does not assign the full irreducible-sector symmetry class",
        "not identify Hamiltonian or Liouvillian integrability",
        "exactly invariant under a common affine rescaling",
        "not invariant under an arbitrary smooth nonlinear transformation",
        "only approximately robust when the mean density varies slowly",
        "does not establish complete positivity or trace preservation",
        "The Lindblad construction supplies those properties independently"), forbidden + (
        "the dense block is chaotic at every N", "the global Liouvillian is class AIII",
        "confirming that the Liouvillian is a proper generator",
        "spacing ratios are invariant under arbitrary smooth nonlinear transformations"))
    for path, (required, forbidden) in surfaces.items():
        source = (root / path).read_text(encoding="utf-8")
        if path.endswith("OpenArcsRegistry.cs"):
            source = route_b_current_arc_slice(source)
        if path == "docs/ANALYTICAL_FORMULAS.md":
            source = source.split("### F41.", 1)[1].split("### F44.", 1)[0]
        source = " ".join(source.replace("**", "").replace("///", "")
                          .replace('" +', "").replace('"', "").split()).casefold()
        def errors(text):
            return task8_text_errors(text, required, forbidden)
        baseline = set(errors(source))
        check(f"Task8 current truth: {path}", not baseline, ascii("; ".join(sorted(baseline))))
        for phrase in required:
            check(f"Task8 rejects removed boundary {ascii(phrase)}: {path}",
                  bool(set(errors(source.replace(phrase.casefold(), ""))) - baseline))
        for phrase in forbidden:
            check(f"Task8 rejects old variant {ascii(phrase)}: {path}",
                  bool(set(errors(source + " " + phrase.casefold())) - baseline))
    arc_source = (root / "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs").read_text(encoding="utf-8")
    current = " ".join(route_b_current_arc_slice(arc_source).replace('" +', "").replace('"', "").split()).casefold()
    arc_required = ("not an invariant space at Delta!=0", "first locate a genuine full-block degeneracy", "their nonzero-Delta response is open")
    for phrase in arc_required:
        check(f"Task8 current OpenArcs requires {phrase}", phrase.casefold() in current)
        check(f"Task8 current OpenArcs mutation rejects removal of {phrase}",
              bool(set(task8_text_errors(current.replace(phrase.casefold(), ""), arc_required, ()))
                   - set(task8_text_errors(current, arc_required, ()))))

    # Each stale claim is an INPUT mutation sent through the same current-text predicate.
    # Comparing new defects to baseline prevents an unrelated pre-existing failure from masking a blind gate.
    mutation_cases = {
        "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs": (
            "public bool Survived => Verdict == DeltaFlipVerdict.Diabolic;",
            "fullGap > coalesceTol", "EpKind.Normal => DeltaFlipVerdict.Defective",
            "EpKind.Normal => DeltaFlipVerdict.Diabolic", "EpKind.NearEp => DeltaFlipVerdict.Defective",
            "EpKind.NearEp => DeltaFlipVerdict.Diabolic", "DeltaFlipVerdict.Uncertified => false",
            "new DeltaTrackResult(DeltaFlipVerdict.Lifted, 0, 0, double.NaN, qd, Complex.Zero, refined)",
            "CertifyCompressedProposal(BuildSym(n, qd, delta), qd, mid)"),
        "compute/RCPsiSquared.Cli/Commands/PathKMonodromyScanCommand.cs": ("if (d > 0 && t.Survived) survivesAtPositive = true;",),
        "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs": (
            "additivity supplies the H-half unconditionally", "Moves 1, 2, 3 and 4 are done",
            "N=7 nonzero-Delta response is closed", "an AT eigenmode has the integer basis-cell disagreement count"),
        "compute/RCPsiSquared.Core/Symmetry/F41PalindromicTimePi2Inheritance.cs": (
            "physical Heisenberg time t_H ~ 4^N", "F41 is the period of the complete SFF",
            "t_Pi/t_H ~ N²/4^N", "slowest palindromic SFF modulation frequency"),
        "compute/RCPsiSquared.Core/Symmetry/F43BandSffPairingPi2Inheritance.cs": (
            "fixed XY-weight eigenvalue sectors", "normalized endpoint K(t) = delta(t)", "endpoint multiplicity 2^N",
            "canonical Pauli-string weight partition", "delta-spike at zero frequency"),
        "docs/proofs/derivations/D09_SECTOR_SFF_PAIRING.md": ("fixed XY-weight eigenvalue sectors", "normalized endpoint K(t) = delta(t)"),
        "experiments/DEGENERACY_HUNT.md": ("every eigenmode has integer n_diff",),
        "experiments/SYMMETRY_CENSUS.md": ("every eigenmode has integer n_diff",),
        "experiments/RANDOM_MATRIX_THEORY.md": ("the dense block is chaotic at every N", "the global Liouvillian is class AIII"),
        "docs/outbound/SHIFTED_ORDER4_CHIRAL_SYMMETRY.md": ("the global Liouvillian is class AIII", "the constant shift is outside SRP"),
        "docs/KMS_DETAILED_BALANCE.md": ("the global Liouvillian is class AIII", "the constant shift is outside SRP"),
        "docs/THE_INTERPRETATION.md": ("the global Liouvillian is class AIII",),
        "experiments/README.md": ("the global Liouvillian is class AIII",),
        "compute/RCPsiSquared.Core/F86/ShiftedGeneratorSectorwisePClaim.cs": (
            "ChiralAiiiClassification", "global SRP class AIII is established"),
        "compute/RCPsiSquared.Runtime/F86Main/F86MainRegistration.cs": ("ChiralAiiiClassification",),
        "compute/RCPsiSquared.Runtime/F86Main/F86OpenQuestionsRegistration.cs": ("ChiralAiiiClassification",),
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdWitness.cs": (
            "Class A licensed", "no residual antiunitary", "GinUE licenses the full symmetry class"),
        "compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs": (
            "public double Rate(int nXY, double gammaZero)", "public double AbsorptionQuantum(double gammaZero)",
            "public double PerCoherenceRateComputationalBasis(int nDiff, double gammaZero)",
            "public double NXyFromRate(double rate, double gammaZero)",
            "eigenmode decay rates are integer multiples of 2γ"),
    }
    mutation_cases["experiments/RANDOM_MATRIX_THEORY.md"] += (
        "spacing ratios are invariant under arbitrary smooth nonlinear transformations",
        "confirming that the Liouvillian is a proper generator",)
    mutation_cases["experiments/F89_PATH_K_DIABOLIC.md"] = (
        "Every positive-Delta N=5/N=6 proposal is Defective", "full-pair gap exceeds coalesceTol=0.001",)
    mutation_cases["compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs"] += (
        "Every positive-Delta N=5/N=6 proposal is Defective",)
    mutation_cases.update({
        "compute/RCPsiSquared.Core/Symmetry/PiOperator.cs": ("existing class-AIII-chiral classification",),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": (
            "N=5/N=6 positive Delta proposals defect or lift", "FillingThresholdWitness class-A gate"),
        "compute/RCPsiSquared.Diagnostics/Foundation/CrossFoldSimilarityWitness.cs": (
            "N=5/N=6 positive Delta proposals defect or lift", "sampled small gaps certify Jordan character"),
        "compute/RCPsiSquared.Core/Symmetry/F89CrossFoldSimilarityClaim.cs": (
            "N=5/N=6 positive Delta proposals defect or lift", "sampled small gaps certify Jordan character"),
        "compute/RCPsiSquared.Core.Tests/F89PathK/WeightCoherenceBlockTests.cs": (
            "N=5/N=6 positive Delta proposals defect or lift",),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/CrossFoldSimilarityWitnessTests.cs": (
            "sampled small gaps certify Jordan character",),
        "simulations/pt_symmetry_analysis.py": ("Class AIII (chiral unitary)", "chiral phase", "chiral breaking"),
        "compute/RCPsiSquared.Compute/README.md": ("chiral symmetry class AIII",),
        "compute/RCPsiSquared.Core/F_FORMULA_CROSSWALK.md": (
            "AbsorptionTheoremClaim.AbsorptionQuantum", "AbsorptionTheoremClaim.Rate(", "ChiralAiiiClassification"),
        "experiments/FILLING_THRESHOLD_CHAOS.md": ("FillingThresholdWitness class-A gate",),
        "experiments/CONCENTRATOR_MAPPING.md": ("FillingThresholdWitness class-A gate",),
        "hypotheses/ZERO_IS_THE_MIRROR.md": ("makes γ the rung spacing of the rate spectrum",),
        "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md": ("Uncertified proposals supply measured departure",),
    })
    mutation_cases["docs/THE_INTERPRETATION.md"] += ("the structure is new rather than a rediscovery is established",)
    for path, variants in mutation_cases.items():
        source = (root / path).read_text(encoding="utf-8")
        if path.endswith("OpenArcsRegistry.cs"):
            source = route_b_current_arc_slice(source)
        required, forbidden = surfaces.get(path, ((), ()))
        baseline = set(task8_text_errors(source, required, forbidden))
        for variant in variants:
            check(f"Task8 integration mutation rejects {ascii(variant)}: {path}",
                  bool(set(task8_text_errors(source + " " + variant, required, forbidden)) - baseline))

    # These mutate the actual edge, branch, threshold, or API signature in place;
    # unrelated valid prose elsewhere cannot satisfy the removed contract.
    replacements = (
        ("compressed dispatch discards caller seed", "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs",
         "TrackDiabolicUnderDeltaCompressed(n, qSeed, lambdaSeed, seedRadius, delta, boxHalf, boxCell)",
         "TrackDiabolicUnderDeltaCompressed(n, qSeed, new Complex(100, 0), seedRadius, delta, boxHalf, boxCell)"),
        ("compressed certificate discards seed isolation", "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs",
         "return CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, mid, lambdaSeed, seedRadius);",
         "return CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, mid);"),
        ("wrong F1 parent edge", "compute/RCPsiSquared.Runtime/F86Main/F86MainRegistration.cs",
         "_ = b.Get<F1PalindromeIdentity>();", "_ = b.Get<DressedModeWeightClaim>();"),
        ("closed SRP question", "compute/RCPsiSquared.Core/F86/ShiftedGeneratorSectorwisePClaim.cs",
         "Tier.OpenQuestion", "Tier.Tier1Derived"),
        ("N7 unsafe locator enabled", "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs",
         "if (n >= 7 && !exact)", "if (n >= 8 && !exact)"),
        ("full coincidence tolerance inflated", "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs",
         "FullBlockCoincidenceTolerance = 1e-6", "FullBlockCoincidenceTolerance = 1e-3"),
        ("residual path bypasses shared certificate", "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs",
         "return CertifyFullBlockProposal(BuildSym(n, center, delta), center, mid, lambdaSeed, seedRadius);",
         "return new DeltaTrackResult(DeltaFlipVerdict.Diabolic, 3, 3, 0, center, mid, 0);"),
        ("default path bypasses shared certificate", "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs",
         "return CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, mid, lambdaSeed, seedRadius);",
         "return new DeltaTrackResult(DeltaFlipVerdict.Defective, 2, 1, 1, qd, mid, .001);"),
        ("algebraic triple accepted", "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs",
         "rr.Algebraic == 2", "rr.Algebraic >= 2"),
        ("midpoint correspondence removed", "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs",
         "(full[1] - mid).Magnitude <= FullBlockCoincidenceTolerance", "true"),
        ("integer eigenmode API", "compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs",
         "EigenmodeDecayRate(double averageNXy, double gammaZero)", "EigenmodeDecayRate(int averageNXy, double gammaZero)"),
        ("N5/N6 positive controls erased", "compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzDeltaFlipTests.cs",
         "Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, d0.Verdict)",
         "Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, d0.Verdict)"),
        ("N6 positive-Delta false death assertion", "compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzDeltaFlipTests.cs",
         "AssertUncertifiedSplit(6, .1, d)", "Assert.False(d.Survived)"),
        ("nonlinear spacing invariance", "experiments/RANDOM_MATRIX_THEORY.md",
         "exactly invariant under a common affine rescaling", "exactly invariant under every smooth nonlinear rescaling"),
    )
    for name, path, original, replacement in replacements:
        source = (root / path).read_text(encoding="utf-8")
        required, forbidden = surfaces[path]
        if name in ("default path bypasses shared certificate", "compressed certificate discards seed isolation"):
            offset = source.rfind(original) if name.startswith("default") else source.find(original)
            changed = source[:offset] + replacement + source[offset + len(original):] if offset >= 0 else source
            # A matching call in another method cannot mask this mutation.
            extra_errors = set(delta_tracker_path_errors(changed)) - set(delta_tracker_path_errors(source))
        else:
            changed = source.replace(original, replacement)
            extra_errors = set()
        check(f"Task8 in-place mutation rejects {name}", changed != source and
              bool((set(task8_text_errors(changed, required, forbidden))
                    - set(task8_text_errors(source, required, forbidden))) | extra_errors))

    tracker = (root / "compute/RCPsiSquared.Diagnostics/Foundation/XxzCoherenceBlock.cs").read_text(encoding="utf-8")
    check("Task8 every public/proposal path shares strict certification", not delta_tracker_path_errors(tracker))
    call = "return CertifyFullBlockProposal(BuildSym(n, qd, delta), qd, mid, lambdaSeed, seedRadius);"
    occurrences = [m.start() for m in re.finditer(re.escape(call), tracker)]
    check("Task8 independently addresses both identical qd certificate calls", len(occurrences) == 2)
    for index, offset in enumerate(occurrences):
        changed = tracker[:offset] + call.replace("lambdaSeed, seedRadius", "null, double.PositiveInfinity") + tracker[offset + len(call):]
        check(f"Task8 single-path mutation rejects certificate bypass occurrence {index}",
              bool(set(delta_tracker_path_errors(changed)) - set(delta_tracker_path_errors(tracker))))
    for name, original, replacement in (
        ("public direct character bypass", "=> TrackDiabolicUnderDelta(n, qSeed, lambdaSeed, delta, boxHalf: 0, boxCell: cell);",
         "=> new DeltaTrackResult(DeltaFlipVerdict.Diabolic, 2, 2, 0, qSeed, lambdaSeed, 0);"),
        ("compressed caller-seed discarded", "TrackDiabolicUnderDeltaCompressed(n, qSeed, lambdaSeed, seedRadius, delta, boxHalf, boxCell)",
         "TrackDiabolicUnderDeltaCompressed(n, qSeed, new Complex(100, 0), seedRadius, delta, boxHalf, boxCell)"),
        ("isolated-pair guard removed", "if (!(radius > 2 * ds[1]))", "if (false)"),
    ):
        changed = tracker.replace(original, replacement, 1)
        check(f"Task8 single-path mutation rejects {name}", changed != tracker and
              bool(set(delta_tracker_path_errors(changed)) - set(delta_tracker_path_errors(tracker))))

    # Outreach owns a different current NextStep, outside the Route-B arc.
    # Its ParkedAt history must neither satisfy nor violate this live boundary.
    def outreach_slice(text):
        _, found, section = text.partition('Name: "outbound_label_adapters"')
        if not found:
            raise ValueError("missing outbound_label_adapters arc")
        section = section.split("new OpenArc(", 1)[0]
        _, found, current = section.partition("NextStep:")
        if not found or "S4 BRIDGE" not in current:
            raise ValueError("missing current outreach NextStep delimiters")
        return current.split("S4 BRIDGE", 1)[0]

    outreach_required = ("shifted-generator, sectorwise involutive P relation",
                         "complete irreducible SRP class remains open",
                         "movement toward the GinUE comparison")
    outreach_forbidden = ("class AIII", "global AIII relation")
    source = (root / "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs").read_text(encoding="utf-8")
    current = outreach_slice(source)
    baseline = set(task8_text_errors(current, outreach_required, outreach_forbidden))
    check("Task8 outreach current class boundary", not baseline, ascii("; ".join(baseline)))
    for name, changed in (
        ("AIII outreach promotion", current + " class AIII is established"),
        ("closed outreach SRP class", current.replace("complete irreducible SRP class remains open", "complete irreducible SRP class is established")),
        ("erased outreach sectorwise P", current.replace("shifted-generator, sectorwise involutive P relation", "global AIII relation")),
    ):
        check(f"Task8 outreach mutation rejects {name}", changed != current and
              bool(set(task8_text_errors(changed, outreach_required, outreach_forbidden)) - baseline))
    # Wrong-input control: a historical ParkedAt promotion is excluded, while
    # a removed structural delimiter must not silently select an empty slice.
    historical = source.replace('Name: "outbound_label_adapters"',
                                'Name: "outbound_label_adapters", ParkedAt: "class AIII is established"')
    check("Task8 outreach keeps ParkedAt history outside current slice", outreach_slice(historical) == current)
    for delimiter in ('Name: "outbound_label_adapters"', "S4 BRIDGE"):
        try:
            outreach_slice(source.replace(delimiter, "missing"))
        except ValueError:
            rejected = True
        else:
            rejected = False
        check(f"Task8 outreach rejects missing delimiter {ascii(delimiter)}", rejected)


def producer_artifact_path_errors(producer, producer_path, artifact_path):
    """Resolve the producer's actual file sink and compare it with the canonical artifact."""
    tree = ast.parse(producer)
    producer_path = Path(producer_path).resolve()
    artifact_path = Path(artifact_path).resolve()
    assignments = []
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name)):
            assignments.append((node.targets[0].id, node.value))

    environment = {"__file__": producer_path}

    def dotted_name(node):
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Attribute):
            prefix = dotted_name(node.value)
            return f"{prefix}.{node.attr}" if prefix else node.attr
        return None

    def evaluate_path(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            return node.value
        if isinstance(node, ast.Name):
            return environment.get(node.id)
        if isinstance(node, ast.Attribute) and node.attr == "parent":
            value = evaluate_path(node.value)
            return Path(value).parent if value is not None else None
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            left = evaluate_path(node.left)
            right = evaluate_path(node.right)
            return Path(left) / right if left is not None and right is not None else None
        if isinstance(node, ast.BoolOp) and isinstance(node.op, ast.Or):
            # Output paths may expose an environment override followed by a
            # canonical fallback. Static verification binds the fallback sink.
            for value_node in reversed(node.values):
                value = evaluate_path(value_node)
                if value is not None:
                    return value
            return None
        if isinstance(node, ast.Call):
            name = dotted_name(node.func)
            args = [evaluate_path(arg) for arg in node.args]
            if any(arg is None for arg in args):
                return None
            if name == "Path" and len(args) == 1:
                return Path(args[0])
            if name == "os.path.abspath" and len(args) == 1:
                return Path(args[0]).resolve()
            if name == "os.path.dirname" and len(args) == 1:
                return Path(args[0]).parent
            if name == "os.path.join" and args:
                return Path(args[0]).joinpath(*args[1:])
        return None

    # Resolve dependencies such as RESULTS = Path(__file__).parent / "results"
    # before the local out_path assignment, without executing the producer.
    for _ in range(len(assignments) + 1):
        changed = False
        for name, value_node in assignments:
            value = evaluate_path(value_node)
            if value is not None and environment.get(name) != value:
                environment[name] = value
                changed = True
        if not changed:
            break

    problems = []
    expected = str(artifact_path).casefold()

    def is_writable_open(call):
        if not (isinstance(call, ast.Call) and isinstance(call.func, ast.Name)
                and call.func.id == "open" and call.args):
            return False
        mode_node = call.args[1] if len(call.args) > 1 else next(
            (keyword.value for keyword in call.keywords if keyword.arg == "mode"), None)
        mode = "r" if mode_node is None else evaluate_path(mode_node)
        return isinstance(mode, str) and any(flag in mode for flag in "wax+")

    open_handles = {}
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and isinstance(node.targets[0], ast.Name) and is_writable_open(node.value)):
            open_handles[id(node.value)] = node.targets[0].id
        elif isinstance(node, ast.withitem) and is_writable_open(node.context_expr):
            if isinstance(node.optional_vars, ast.Name):
                open_handles[id(node.context_expr)] = node.optional_vars.id

    writer_handles = set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "write" and isinstance(node.func.value, ast.Name)):
            writer_handles.add(node.func.value.id)

    sinks = []
    for node in ast.walk(tree):
        if is_writable_open(node):
            handle = open_handles.get(id(node))
            label = handle or "unbound open"
            sinks.append((label, node.args[0], handle in writer_handles))
        elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
              and node.func.attr in ("write_text", "write_bytes")):
            sinks.append((dotted_name(node.func.value) or "direct writer",
                          node.func.value, True))

    actual_sinks = [sink for sink in sinks if sink[2]]
    if not actual_sinks:
        return ["producer has no recognized writer associated with a file sink"]
    if len(sinks) != 1:
        problems.append(f"producer has {len(sinks)} writable file sinks; expected exactly one")

    for label, path_node, is_used in sinks:
        if not is_used:
            problems.append(f"unused writable file sink {label}")
        value = evaluate_path(path_node)
        if value is None:
            problems.append(f"cannot resolve producer file sink {label}")
            continue
        observed = Path(value)
        if not observed.is_absolute():
            observed = producer_path.parent / observed
        if str(observed.resolve()).casefold() != expected:
            problems.append(f"{label} resolves to {observed.resolve()}, expected {artifact_path}")
    return problems


def phase3_artifact_integrity_errors(producer, artifact, matcher_source=""):
    """Bind the deterministic Phase-3 snapshot to its source bundle and every payload byte."""
    problems = []
    source_matches = re.findall(r"^Source SHA256: ([0-9a-f]{64})$", artifact, re.MULTILINE)
    source_bundle = producer + ("\0" + matcher_source if matcher_source else "")
    expected_source = hashlib.sha256(source_bundle.encode("utf-8")).hexdigest()
    if source_matches != [expected_source]:
        problems.append("Phase-3 artifact source hash differs from producer")

    marker = "Payload SHA256: "
    payload, separator, suffix = artifact.rpartition(marker)
    digest = suffix.strip()
    if not separator or not re.fullmatch(r"[0-9a-f]{64}", digest):
        problems.append("Phase-3 artifact has no unique terminal payload hash")
    elif hashlib.sha256(payload.encode("utf-8")).hexdigest() != digest:
        problems.append("Phase-3 artifact payload hash differs from its content")
    return problems


def emitted_artifact_claim_errors(producer, artifact, markers):
    """Semantic freshness of protected emitted sentences, not timestamp/hash freshness.

    No producer is imported or executed here: PT/entropy write their result at import.
    Numerical output remains separately checked where an independent input route exists.
    """
    def normalize(text):
        return " ".join(text.split()).casefold()
    emitted = []
    for node in ast.walk(ast.parse(producer)):
        if not isinstance(node, ast.Call) or not node.args:
            continue
        is_log = isinstance(node.func, ast.Name) and node.func.id == "log"
        is_append = (isinstance(node.func, ast.Attribute) and node.func.attr == "append"
                     and isinstance(node.func.value, ast.Name) and node.func.value.id == "out")
        if (is_log or is_append) and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str):
            emitted.append(node.args[0].value)
    problems = []
    for marker in markers:
        messages = [m for m in emitted if normalize(marker) in normalize(m)]
        if not messages:
            problems.append(f"producer does not emit protected claim {marker}")
        for message in messages:
            if normalize(message) not in normalize(artifact):
                problems.append(f"artifact stale/missing emitted claim {message.strip()}")
    return problems


def rmt_artifact_row_errors(artifact, eigenvalues_by_n):
    """Independent triple-by-triple oracle: never remove a gap or call the producer."""
    problems = []
    for n, eigenvalues in eigenvalues_by_n.items():
        section = re.search(rf"^N={n}: (\d+) eigenvalues(?P<body>[\s\S]*?)(?=^N=\d+:|\Z)", artifact, re.M)
        if not section:
            problems.append(f"N={n} missing printed input row")
            continue
        if int(section[1]) != len(eigenvalues):
            problems.append(f"N={n} population differs from CSV")
        all_rates = sorted(float(-value.real) for value in eigenvalues)
        for label, rates in (("a", all_rates), ("b", [x for x in all_rates if x < n * 0.05])):
            ratios, undefined = [], 0
            for left, middle, right in zip(rates, rates[1:], rates[2:]):
                s1, s2 = middle - left, right - middle
                if s1 == 0 and s2 == 0:
                    undefined += 1
                else:
                    ratios.append(min(s1, s2) / max(s1, s2))
            expected = dict(levels=len(rates), gaps=max(0, len(rates)-1),
                zero_gaps=sum(a == b for a, b in zip(rates, rates[1:])),
                adjacent_pairs=max(0, len(rates)-2), defined_ratios=len(ratios), undefined_zero_zero=undefined)
            row = re.search(rf"\({label}\) [^\n]+\n\s*([^\n]+)\n\s*<r> defined = ([\w.+-]+)", section['body'])
            if not row:
                problems.append(f"N={n} {label} missing population/ratio row")
                continue
            actual = dict((key, int(value)) for key, value in re.findall(r"(\w+)=(\d+)", row[1]))
            if actual != expected:
                problems.append(f"N={n} {label} multiset population differs from CSV")
            measured = float(row[2])
            expected_ratio = sum(ratios) / len(ratios) if ratios else np.nan
            # 17 printed significant digits; allow only accumulation-order rounding,
            # not deletion of small gaps or an eigensolver resolution tolerance.
            matches = (np.isnan(measured) if np.isnan(expected_ratio)
                       else np.isfinite(measured) and abs(measured - expected_ratio) <= 2e-14)
            if not matches:
                problems.append(f"N={n} {label} ratio differs from independent CSV calculation")
    return problems


def rmt_band_artifact_errors(artifact, eigenvalues_by_n):
    """Rebuild rate-window membership and every absolute-frequency gap from CSV."""
    problems = []
    for n, eigenvalues in eigenvalues_by_n.items():
        section = re.search(rf"^N={n}: (\d+) eigenvalues(?P<body>[\s\S]*?)(?=^N=\d+:|\Z)", artifact, re.M)
        if not section or int(section[1]) != len(eigenvalues):
            problems.append(f"band N={n} missing/wrong input population")
            continue
        for w in range(n + 1):
            # Independent scalar selection, not the producer's band helper.
            levels = sorted(abs(float(z.imag)) for z in eigenvalues
                            if abs(-float(z.real) - 2*w*0.05) < 0.3*0.05)
            gaps = [b-a for a, b in zip(levels, levels[1:])]
            ratios = [min(a,b)/max(a,b) for a,b in zip(gaps, gaps[1:]) if a != 0 or b != 0]
            expected = dict(levels=len(levels), gaps=len(gaps), zero_gaps=gaps.count(0),
                adjacent_pairs=max(0,len(gaps)-1), defined_ratios=len(ratios),
                undefined_zero_zero=sum(a == b == 0 for a,b in zip(gaps,gaps[1:])))
            row = re.search(rf"w={w}:\s*([^\n]+)\n\s*<r> defined = ([\w.+-]+)", section['body'])
            if not row:
                problems.append(f"band N={n} w={w} missing row")
                continue
            actual = dict((key,int(value)) for key,value in re.findall(r"(\w+)=(\d+)",row[1]))
            if actual != expected:
                problems.append(f"band N={n} w={w} population differs from CSV")
            mean = sum(ratios)/len(ratios) if ratios else np.nan
            printed = float(row[2])
            if not (np.isnan(printed) if np.isnan(mean) else np.isfinite(printed) and abs(printed-mean) <= 2e-14):
                problems.append(f"band N={n} w={w} ratio differs from CSV")
    return problems


def rmt_document_table_errors(document, artifact, band_artifact):
    """Doc tables consume the artifact rows that the independent CSV oracles verify."""
    problems = []
    result1 = document.split("## Result 1:")[-1].split("## Result 2:")[0]
    result2 = document.split("## Result 2:")[-1].split("## Result 3:")[0]
    result3 = document.split("## Result 3:")[-1].split("## Result 4:")[0]
    rows = [[cell.strip() for cell in line.strip().strip('|').split('|')]
            for line in result1.splitlines() if re.match(r"\|\s*\d+\s*\|", line)]
    band_rows = [[cell.strip() for cell in line.strip().strip('|').split('|')]
                 for line in result3.splitlines() if re.match(r"\|\s*\d+\s*\|", line)]
    reflection_rows = [[cell.strip() for cell in line.strip().strip('|').split('|')]
                       for line in result2.splitlines() if re.match(r"\|\s*\d+\s*\|", line)]
    expected_rows = []
    for population in ('a', 'b'):
        for n in range(2,8):
            section = re.search(rf"^N={n}: [\s\S]*?(?=^N=\d+:|\Z)", artifact, re.M)
            match = re.search(rf"\({population}\) [^\n]+\n\s*([^\n]+)\n\s*<r> defined = ([\w.+-]+)", section[0]) if section else None
            expected_rows.append(([n], match))
    expected_band_rows = []
    for n,w in ((4,2),(5,2),(5,3),(6,3),(7,3),(7,4)):
        section = re.search(rf"^N={n}: [\s\S]*?(?=^N=\d+:|\Z)", band_artifact, re.M)
        match = re.search(rf"w={w}:\s*([^\n]+)\n\s*<r> defined = ([\w.+-]+)", section[0]) if section else None
        expected_band_rows.append(([n,w],match))
    for label, printed, expected in (("pooled", rows, expected_rows), ("band", band_rows, expected_band_rows)):
        if len(printed) != len(expected):
            problems.append(f"RMT {label} document row count")
        for index, (prefix, match) in enumerate(expected):
            if not match or index >= len(printed):
                problems.append(f"RMT {label} missing row {prefix}")
                continue
            values = dict((k,int(v)) for k,v in re.findall(r"(\w+)=(\d+)",match[1]))
            counts = prefix + [values[k] for k in ('levels','zero_gaps','defined_ratios','undefined_zero_zero')]
            try:
                cells = printed[index]
                if ([int(x) for x in cells[:-1]] != counts or not np.isfinite(float(cells[-1]))
                        or abs(float(cells[-1])-float(match[2])) > 2e-14):
                    problems.append(f"RMT {label} document numeric drift {prefix}")
            except ValueError:
                problems.append(f"RMT {label} invalid row {prefix}")
    if len(reflection_rows) != 6:
        problems.append("RMT reflection document row count")
    for index, n in enumerate(range(2, 8)):
        section = re.search(rf"^N={n}: [\s\S]*?(?=^N=\d+:|\Z)", artifact, re.M)
        match = re.search(
            r"positive=(\d+), negative=(\d+), central=(\d+), "
            r"multiplicity bottleneck error: ([\w.+-]+)", section[0]) if section else None
        if not match or index >= len(reflection_rows):
            problems.append(f"RMT reflection missing row N={n}")
            continue
        try:
            cells = reflection_rows[index]
            integers = [int(cell.replace(',', '')) for cell in cells[:-1]]
            expected_integers = [n, int(match[1]), int(match[2]), int(match[3])]
            if integers != expected_integers or float(cells[-1]) != float(match[4]):
                problems.append(f"RMT reflection document numeric drift N={n}")
        except ValueError:
            problems.append(f"RMT reflection invalid row N={n}")
    return problems


def verify_round4_artifacts_and_consumers():
    root = Path(__file__).resolve().parents[1]
    entropy_required = ("chosen I/d reference", "largest sampled dS/dt",
        "algebraic decay-rate pair sum", "not a thermodynamic entropy-production scale",
        "formal gain-loss generator is not a physical Lindblad channel",
        "four-qubit GHZ/cat coherence", "all palindrome pairs")
    entropy_forbidden = ("maximum entropy production", "σ peak", "ρ → ρ_ss", "Jarzynski", "Crooks",
        "β_eff", "infinite-temperature bath", "T_eff", "Carnot", "η =", "No net entropy production", "heat engine",
        "Bell pair across the bridge")
    rmt_required = ("affine-invariant", "no irreducible-class inference", "no cutoff or deduplication",
        "0/0 = undefined", "mean over defined ratios only", "No standard Poisson/GOE/GUE calibration",
        "near-degenerate", "undefined_zero_zero", "decay-rate projection reflection",
        "multiplicity-preserving perfect matching of the full real-part projection",
        "does not match the full complex multiset")
    rmt_forbidden = ("def unfold", "NNSD", "Nearest:", "Class (all)", "Class (half)", "classify_r",
        "spacings = spacings[spacings >", "rates[rates > 1e-10]", "Poisson=0.386")
    n4_false = ("N=4 positive-Delta proposals are certified defective",)
    surfaces = {
        "compute/RCPsiSquared.Core/Symmetry/F44AlgebraicPairRateLogIdentityPi2Inheritance.cs": ((
            "F44 algebraic pair-rate log identity", "inverse-rate linear coefficient",
            "No probability/work ensemble or dimensionless thermodynamic exponent",
            "finite and satisfy 0 ≤ Δd < 2·Σγ"),
            ("Crooks-like", "β_eff", "EffectiveInverseTemperature", "EmpiricalJarzynskiMean",
             "IsCrooksFluctuationTheorem", "Jarzynski equality fails")),
        "compute/RCPsiSquared.Runtime/PolarityArchitecture/F44AlgebraicPairRateLogIdentityPi2InheritanceRegistration.cs": ((
            "F44AlgebraicPairRateLogIdentityPi2Inheritance", "no probability/work ensemble",
            "dimensionless thermodynamic exponent", "0≤Δd&lt;2Σγ"),
            ("F44CrooksLikeRateIdentityPi2Inheritance", "β_eff", "Jarzynski")),
        "compute/MirrorWorld/Formulas.cs": (("F44 (T1, D08): algebraic pair-rate log identity",
            "no probability/work ensemble is defined", "finite ordered positive palindromic pair",
            "finite b > 0", "positive decay variable z=-lambda"),
            ("F44 (T1, D08): Crooks-like", "β_eff", "Jarzynski")),
        "docs/proofs/derivations/D08_CROOKS_RATE_IDENTITY.md": ((
            "Algebraic Pair-Rate Log Identity", "inverse-rate scale", "No probability ratio",
            "exponent and weighting were not derived from a work protocol",
            "0 ≤ Δd < 2Σγ", "exact typed API and its linear approximation"),
            ("β_eff = 1/Σγ", "effective inverse temperature")),
        "docs/proofs/derivations/README.md": ((
            "D8 (algebraic pair-rate log identity",),
            ("D8 (Crooks rate identity",)),
        "docs/ANALYTICAL_FORMULAS.md": (("F44. Algebraic pair-rate log identity",
            "algebraic linear coefficient is 1/Σγ", "not a Jarzynski test",
            "0 <= Delta_d < 2*Σγ", "finite b > 0", "positive decay variable `z=−λ`",
            "F40. Fragile bridge axis-departure threshold",
            "does not execute a threshold coalescence or Jordan-rank test",
            "EP, Hopf, and Jordan character remain OPEN"),
            ("F44. Crooks-like", "β_eff = 1/Σγ", "F40. Fragile bridge gamma_crit",
             "coalesces and leaves the axis at a second-order exceptional", "diverges as 1/delta",
             "b=−γ₀")),
        "review/OPEN_QUESTIONS_INDEX.md": (("OQ-168", "pure Z-dephasing assigns no bath temperature",
            "is neither KMS nor quantum detailed balance", "KMS/operator connection open"),
            ("pure dephasing (infinite-temperature bath)", "pure Z-dephasing is an infinite-temperature bath")),
        "compute/RCPsiSquared.Core/Symmetry/F95AngleAtQuadraticZeroPi2Inheritance.cs": ((
            "Returns zero at c = 1/4", "if (c < Threshold)", "if (c < thresh)",
            "finite b > 0", "if (!double.IsFinite(b) || b <= 0.0)"),
            ("Returns NaN for c ≤ 1/4", "if (c <= Threshold)", "if (c <= thresh)",
             "with real (b, c)")),
        "compute/RCPsiSquared.Runtime/PolarityArchitecture/F95AngleAtQuadraticZeroPi2InheritanceRegistration.cs": ((
            "positive-b closed form", "finite b &gt; 0", "currently registered callers use b &gt; 0"),
            ("universal closed form",)),
        "docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md": ((
            "finite b > 0", "second-quadrant principal argument π−θ",
            "b = 0 gives π/2", "neither branch is part of the current F95 API"),
            ("real coefficients (b, c) ∈ ℝ²",)),
        "docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md": ((
            "finite $b>0$", "positive-b polynomial-foundation identity used here"),
            ("valid for any monic quadratic", "universal angle-emergence formula at any quadratic discriminant zero")),
        "experiments/BOUNDARY_NAVIGATION.md": ((
            "finite b > 0", "positive decay variable z=−λ", "angle from the negative λ-axis"),
            ("for any monic quadratic", "universal polynomial-foundation identity")),
        "compute/RCPsiSquared.Core/Symmetry/F97CardioidHalfFixedPointPi2Inheritance.cs": ((
            "real c > b² and finite b>0",), ()),
        "docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md": ((
            "real-c, finite-b>0", "framework's b=1/2 case"), ()),
        "docs/proofs/INCOMPLETENESS_PROOF.md": (("finite-positive-b quadratic angle formula",),
            ("universal quadratic angle formula",)),
        "simulations/angle_at_zero_tier1_candidate.py": ((
            "finite b > 0", "if not math.isfinite(b) or b <= 0:"), ()),
        "compute/RCPsiSquared.Core/Symmetry/TransitionBridgeF95SiblingClaim.cs": ((
            "zero at the EP", "if (q == QEp(gEff)) return 0.0", "if (disc > 0.0)",
            "FRAGILE_BRIDGE spectral-abscissa axis departure; EP/Hopf/Jordan character OPEN",
            "positive decay variable z=−λ", "b = 4γ₀ &gt; 0"),
            ("NaN at or below the EP", "if (disc >= 0.0)")),
        "compute/RCPsiSquared.Runtime/PolarityArchitecture/TransitionBridgeF95SiblingClaimRegistration.cs": ((
            "genuine F86 toy 2x2 exceptional point", "angle zero at each boundary",
            "FRAGILE_BRIDGE axis departure has OPEN EP/Hopf/Jordan character"),
            ("FRAGILE_BRIDGE exceptional point", "FRAGILE_BRIDGE axis departure is a Hopf")),
        "compute/RCPsiSquared.Diagnostics/Foundation/BranchLocusPalindromeWitness.cs": ((
            "Delta=0", "N=4/N=5/N=6 positive Delta proposals remain Uncertified"),
            ("The sampled XXZ Delta control stays on-line yet defects",)),
        "compute/RCPsiSquared.Core/Symmetry/F89BranchLocusPalindromeClaim.cs": ((
            "Delta=0", "N=4/N=5/N=6 positive Delta proposals remain Uncertified"),
            ("control stays on-line yet becomes defective",)),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": ((
            "N=4", "N=4/N=5/N=6 positive Delta proposals remain Uncertified"),
            ("the certified N=4 XXZ Delta control defects",)),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/XxzDeltaFlipTests.cs": ((
            "CharacterAtDiabolicNear_RejectsUnrelatedLambdaSeed", "new Complex(100, 0)",
            "N4_DeltaProposals_RequireStrictFullPairCoincidence", "AssertUncertifiedSplit(4, delta, result)",
            "Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, r0.Verdict)"),
            ("N4_DeltaFlip_ReproducesCommittedTable", "Assert.InRange(r2.Departure, 0.012, 0.030)")),
        "simulations/results/pt_symmetry_analysis.txt": (("sectorwise P", "full irreducible SRP class remains OPEN",
            "within palindrome axis", "off palindrome axis; operator relation remains exact",
            "reports the spectral-abscissa axis departure only",
            "does not report single-eigenvector Petermann factors or angles",
            "Degenerate eigenspaces make those values basis-dependent"),
            ("Class AIII", "chiral phase", "chiral breaking", "No classical EP", "No exceptional point", "Hopf",
             "Eigenvector coalescence", "MIRROR-PARTNER ANGLE CONTROL",
             "Petermann factor K and phase rigidity")),
        "simulations/pt_palindrome_breaking.py": (("trace(L)/dim(L)", "IMAGINARY-AXIS OCCUPANCY",
            "does not classify the threshold as Hopf, EP, or a Jordan defect",
            "this sweep does not isolate either change as the cause",
            "Finite-offset most-unstable-mode check", "No branch continuation is performed"),
            ("chiral phase", "still a HOPF bifurcation", "changed to SADDLE-NODE",
             "palindrome IS the protection mechanism", "added damping dominates",
             "instability-branch check", "unstable branch oscillates")),
        "simulations/results/pt_palindrome_breaking.txt": (("trace(L)/dim(L)", "IMAGINARY-AXIS OCCUPANCY",
            "does not classify the threshold as Hopf, EP, or a Jordan defect",
            "this sweep does not isolate either change as the cause",
            "Finite-offset most-unstable-mode check", "No branch continuation is performed"),
            ("chiral phase", "still a HOPF bifurcation", "changed to SADDLE-NODE",
             "palindrome IS the protection mechanism", "added damping dominates",
             "instability-branch check", "unstable branch oscillates")),
        "experiments/PT_SYMMETRY_ANALYSIS.md": (("the final SRP class remains OPEN", "sqrt(p_x) Pi",
            "Spectral reflection alone does not assign AIII", "m = trace(L)/dim(L)",
            "imaginary-axis regime", "does not isolate their causal contributions",
            "2γ - 4epsilon/3", "eigenmode property of the full coupled generator",
            "EP, Hopf, and Jordan character remain OPEN", "spectral-abscissa statement",
            "transient norm growth can occur", "not a physical Lindblad channel",
            "does not report single-eigenvector Petermann factors or eigenvector angles",
            "basis-dependent inside a degenerate eigenspace", "invariant subspace-level conditioning",
            "nonstationary spectrum", "trace-preserving stationary eigenvalue remains at lambda = 0",
            "the stationary zero does not move", "Fragile-bridge axis-confined regime",
            "Fragile-bridge off-axis regime", "Pi remains exact on both sides"),
            ("The full irreducible SRP class is AIII", "chiral phase",
             "not the protection mechanism", "-(gamma - epsilon/3)", "Effect (2) dominates",
             "rules out the prior Hopf", "is a Hopf bifurcation", "is not a Hopf",
             "the system oscillates but does not grow", "This is the explosion",
             "CHIRAL (sublattice) symmetry", "chiral symmetry of the full operator structure",
             "confirming the chiral symmetry", "MIRROR-PARTNER ANGLE CONTROL",
             "Our system: linear Pi, imaginary eigenvalues in symmetric phase",
             "| System | Symmetric phase | Broken phase |",
             "| 1.463 | 402.7", "Above γ_crit: cos(theta) ~ 0.09",
             "by ANY epsilon > 0: eigenvalues immediately develop nonzero Re parts",
             "The eigenvalues move off the imaginary axis")),
        "hypotheses/FRAGILE_BRIDGE.md": (("spectral-abscissa threshold γ_crit",
            "EP, Hopf, and Jordan character remain OPEN",
            "does not perform branch continuation or execute a strict threshold coalescence or Jordan-rank certificate",
            "If a limit exists, is it 1/2?", "have not converged",
            "three finite-coupling frequency samples do not determine a limiting frequency",
            "Which mode geometry produces the turnover remains open",
            "derived partner gap adds no independent evidence",
            "across-axis partner −λ* is not supplied by Π alone"),
            ("not a Hopf", "is a Hopf bifurcation", "genuine EP at real γ_crit",
             "γ_crit × J_bridge → 0.50", "Perfect coupling between gain and loss = immediate instability",
             "saddle-node in character", "the two chains merge into one",
             "Dimer formation, destabilization", "| Resonance |")),
        "simulations/fragile_bridge_ep_signature.py": ((
            "spectral-abscissa axis departure", "does not perform branch continuation",
            "single-vector K not reported for a degenerate eigenspace",
            "derived from max Re and supplies no independent EP evidence",
            "near-threshold simplicity gate: nearest gap > 1e-10",
            "if gap <= SIMPLE_MODE_GAP_TOL:", "raise RuntimeError("),
            ("tracked axis departure", "secondary bump", "nearest-neighbour gap")),
        "simulations/results/fragile_bridge_ep_signature.txt": ((
            "spectral-abscissa axis departure", "not branch continuation",
            "single-vector K not reported for a degenerate eigenspace",
            "derived from max Re and supplies no independent EP evidence",
            "near-threshold simplicity gate: nearest gap > 1e-10"),
            ("tracked axis departure", "secondary bump", "nearest-neighbour gap")),
        "hypotheses/ZERO_IS_THE_MIRROR.md": (("spectral-abscissa axis departure",
            "off-axis quartets appear", "No branch continuation"),
            ("tracked axis departure", "tracked pair")),
        "hypotheses/PAIR_BREAKING_AT_THE_HORIZON.md": (("spectral-abscissa axis departure",
            "performs no branch continuation", "turnover mechanism remains open"),
            ("tracked axis departure", "tracked mirror pair", "tracked gain-loss instability")),
        "docs/proofs/PROOF_F86A_EP_MECHANISM.md": (("spectral-abscissa axis departure",
            "EP character OPEN"), ("tracked axis crossing",)),
        "experiments/F86_EP_THROUGH_THE_CLOCK.md": (("separate spectral-abscissa axis departure",
            "EP, Hopf, and Jordan character remain OPEN", "Only the toy 2×2 end has the stated discriminant"),
            ("is a Hopf bifurcation", "SEPARATE genuine EP", "two SEPARATE genuine EPs",
             "the dynamics escape into self-sustained oscillation", "Hopf can push")),
        "simulations/memory_fate_across_the_takt.py": (("linear amplitude GROWS",
            "EP, Hopf, and Jordan character remain OPEN"),
            ("real gain-loss Hopf", "Hopf / runaway", "the Hopf, the runaway")),
        "simulations/entropy_production.py": (entropy_required, entropy_forbidden),
        "simulations/results/entropy_production.txt": (entropy_required, entropy_forbidden),
        "experiments/ENTROPY_PRODUCTION.md": (("no thermodynamic interpretation is established",
            "not an inverse temperature", "run does not establish a unique stationary state",
            "not a Jarzynski test", "largest sampled dS/dt", "not a continuum-time peak",
            "grows logarithmically as t approaches zero from above",
            "both a Z-dephasing pointer state and an eigenstate of this Hamiltonian",
            "Z-diagonality alone is not enough", "|01⟩⟨01| initially has no Z-basis coherence",
            "root refinement gives CΨ=1/4 at t ≈ 0.747",
            "both Z-dephasing and the Hamiltonian preserve this Bell+ population",
            "no variance-peak time is defined", "four-qubit GHZ/cat coherence",
            "not a Bell pair"),
            ("The rate sum is a thermodynamic entropy-production scale", "reports a peak dS/dt",
             "reported dS/dt peak", "Only initial states with coherence",
             "CΨ crosses 1/4 at t ≈ 0.51", "Max variance: 0.250000 at t = 0.01",
             "Bell pair across the bridge")),
        "compute/RCPsiSquared.Core.Tests/F86/F86PetermannProbe.cs": ((
            "spectral-abscissa axis departure", "EP/Hopf/Jordan character",
            "single-eigenvector K inside a degenerate eigenspace is basis-dependent"),
            ("genuine EP on the real γ axis",)),
        "compute/RCPsiSquared.Core/F86/PolarityInheritanceLink.cs": ((
            "FRAGILE_BRIDGE's spectral-abscissa axis departure", "EP character OPEN"),
            ("FRAGILE_BRIDGE's tracked axis departure",)),
        "compute/RCPsiSquared.Core/F86/F86KnowledgeBase.cs": ((
            "FRAGILE_BRIDGE spectral-abscissa axis departure", "EP character OPEN"),
            ("tracked real-gamma axis departure",)),
        "docs/proofs/PROOF_F86B_UNIVERSAL_SHAPE.md": ((
            "sampled spectral-abscissa axis departure", "EP character is OPEN"),
            ("tracked real-γ axis crossing",)),
        "docs/NAVIGATING_THE_DIMENSIONS.md": ((
            "sampled spectral-abscissa axis departure", "no branch continuation",
            "EP character remains OPEN"),
            ("tracked real-γ axis crossing",)),
        "compute/RCPsiSquared.Diagnostics/Knowledge/KnowledgeRegistryFactory.cs": ((
            "sampled spectral-abscissa axis departure", "no branch continuation was executed",
            "EP character OPEN"),
            ("tracked axis departure",)),
        "simulations/ep_transition.py": ((
            "two-level exceptional-point pinch", "spectral character remains open"),
            ("fragile-bridge pinch",)),
        "compute/RCPsiSquared.Core.Tests/F86/LocalGlobalEpLinkTests.cs": ((
            "spectral-abscissa axis departure", "tracked real-gamma axis departure"),
            ()),
        "compute/RCPsiSquared.Core.Tests/Symmetry/TransitionBridgeF95SiblingClaimTests.cs": ((
            "FRAGILE_BRIDGE spectral-abscissa axis departure", "tracked axis departure"),
            ()),
        "experiments/README.md": (("no temperature, heat/work, efficiency, or unique global stationary state is inferred",
            "full SRP class remains OPEN", "not a Jarzynski observable", "r = +0.987",
            "causal contributions are not isolated", "real-γ spectral-abscissa axis departure"),
            ("The rate sum is a thermodynamic entropy-production scale", "+0.988",
             "stabilizes the sampled system") + n4_false),
        "simulations/rmt_analysis.py": (rmt_required, rmt_forbidden),
        "simulations/results/rmt_analysis.txt": (rmt_required, rmt_forbidden),
        "experiments/RANDOM_MATRIX_THEORY.md": (("No explicit spectral unfolding was applied",
            "only approximately robust when the mean density varies slowly",
            "not invariant under general smooth reparameterizations", "0/0 is undefined",
            "no standard Poisson/GOE/GUE calibration", "no cutoff or deduplication",
            "The F1 palindrome (the centered spectrum has exact ± pairing)",
            "Result 2: Centered Decay-Rate Reflection Check",
            "multiplicity-preserving perfect matching", "does not match the imaginary parts",
            "does not isolate block mixing as the cause",
            "does not measure whether that pairing changes the pooled statistic"),
            ("rank unfolding supplies a valid NNSD", "spacing ratios are invariant under nonlinear maps",
             "The chiral symmetry (centered spectrum has exact ± pairing)",
             "chiral symmetry Liouvillian eigenvalues",
             "Result 2: Exact Centered Spectral Pairing", "What does the work is the conserved structure",
             "it is not the one that suppresses repulsion here")),
        "simulations/rmt_goe_hint_verdict.py": (("average-light bands", "zero frequencies and multiplicities retained",
            "no standard-ensemble calibration", "half-width 0.15", "0/0 undefined"),
            ("classify_r", "np.unique", "np.round", "poisson_band", "XY-weight-w sector")),
        "simulations/results/rmt_band_multiset.txt": (("average-light bands", "zero frequencies and multiplicities retained",
            "no standard-ensemble calibration", "half-width 0.15", "0/0 undefined"),
            ("Poisson <r>=", "WITHIN Poisson band", "GOE-like", "XY-weight-w sector")),
        "docs/WHAT_WE_FOUND.md": (("cell costs form a 2γ-spaced ladder",
            "Hamiltonian mixing can make non-integer", "need not lie on a 2γ grid",
            "J_bridge/J in [1.8, 2.0]", "reaches 0.508 at J_bridge=100",
            "does not determine its limiting value", "spectral-abscissa axis departure",
            "independently selects max Re at each γ"),
            ("interacting eigenmode rates form a 2γ-spaced ladder",
             "γ_crit × J_bridge = 0.50", "approaches a constant: 0.50",
             "optimal region (twice the internal coupling)", "tracked oscillatory axis departure",
             "tracked mirror pair")),
        "experiments/F89_BRANCH_LOCUS_PALINDROME.md": (("N=4 Delta=0 control", "is certified diabolic",
            "sampled positive-Delta proposals remain split", "neither persistence of an on-line degeneracy nor defective character"), n4_false),
        "experiments/F89_PATH_K_DIABOLIC.md": (("at N=4/N=5/N=6 only the Delta=0 controls are Diabolic (alg=geo=2)",
            "Every sampled positive-Delta proposal is Uncertified", "no positive-Delta character is strictly certified at N=4 through N=7"), n4_false),
        "review/EMERGING_QUESTIONS.md": (("PT-like phenomenology, full class OPEN",
            "not licensed instances of a computed global AIII class",
            "No analytic continuation connecting them and no global AIII classification has been established",
            "does not identify K₁ with the full Liouvillian Π"), ("Π class AIII chiral", "K_1 = Π on H_1 sine-mode basis")),
        "hypotheses/DIABOLIC_BY_INTEGRABILITY.md": (("Exact N=4 Delta=0 twin-scalar restriction and Tier-2 residual interpretation",
            "their character and departure are Uncertified", "Only Delta=0 passes the strict full-block double-root/correspondence test"), n4_false),
        "hypotheses/README.md": (("At Delta=0 the F89 path-3 octic crossing is certified diabolic",
            "sampled positive-Delta proposals are split and Uncertified",
            "no diabolic-to-defective flip or general integrability protection is established",
            "spectral abscissa leaves zero and off-axis quartets appear"),
            n4_false + ("tracked pair leaves the imaginary axis",)),
    }
    fragile_tail_surfaces = {
        "compute/RCPsiSquared.Core/F_FORMULA_CROSSWALK.md": ((
            "FRAGILE_BRIDGE spectral-abscissa axis departure whose EP character is OPEN",),
            ("FRAGILE_BRIDGE tracked axis crossing",)),
        "docs/ANALYTICAL_FORMULAS.md": (("F19. Fragile bridge sampled large-coupling tail",
            "do not establish that a limit exists or that its value is 0.50",
            "Does not replace: an asymptotic analysis"),
            ("F19. Fragile bridge asymptotic constant", "gamma_crit * J_bridge -> 0.50")),
        "docs/GLOSSARY.md": (("sampled maximum is bracketed", "0.578 at J_bridge=10",
            "0.508 at 100 and has not converged", "spectral abscissa leaves zero",
            "off-axis quartets appear"),
            ("optimal (2x internal coupling)", "Asymptotic constant γ_crit x J_bridge = 0.50",
             "tracked pair leaves the imaginary axis")),
        "review/OPEN_QUESTIONS_INDEX.md": (("if a limit exists, is it 1/2?",
            "still decreases from 0.578 at J_bridge=10 to 0.508 at 100"),
            ("proposed large-coupling product limit γ_crit × J_bridge → 1/2",)),
        "review/OPEN_QUESTIONS_INDEX_PROPOSAL_math-proof.md": (("if so, is its limit 1/2?",
            "have not converged"),
            ("computed quantum large-bridge asymptote γ_crit × J_bridge → 0.50",)),
    }
    for path, (required_add, forbidden_add) in fragile_tail_surfaces.items():
        required, forbidden = surfaces.get(path, ((), ()))
        surfaces[path] = (required + required_add, forbidden + forbidden_add)
    physics_scope_surfaces = {
        "compute/RCPsiSquared.Core/Symmetry/F2W1DispersionPi2Inheritance.cs": ((
            "zero-frequency decaying mode", "at λ = −2γ", "not stationary"),
            ("zero-frequency stationary mode",)),
        "docs/THE_ATMOSPHERE_AND_THE_CANCELLED_FORMULAS.md": ((
            "basis-cell costs", "eigenmode expectation", "eigenmode decay is not generally quantised"),
            ("decay quantised in units of γ₀", "clean discrete ladder exists only because γ is uniform")),
        "review/OPEN_QUESTIONS_INDEX_PROPOSAL_math-proof.md": ((
            "gap=2γ requires the D6 coupling threshold", "weak-coupling mixed modes lie below that floor"),
            ("spectral gap setting (one absorption quantum)",)),
        "compute/RCPsiSquared.Core/F86/DressedModeWeightClaim.cs": ((
            "Tier 2 empirical", "does not identify Q_peak with an exceptional point",
            "does not identify Q_peak with an EP or isolate a causal eigenvalue pair"),
            ("generalised exceptional-point resonance condition", "generalised EP resonance")),
        "docs/ANALYTICAL_FORMULAS.md": ((
            "Tier-2 empirical weight redistribution", "does not identify the full block's Q_peak",
            "any full-block EP mechanism remain open"),
            ("Q_peak is a generalised exceptional-point resonance condition",)),
        "experiments/COCKPIT_UNIVERSALITY.md": ((
            "old blanket Petermann null is refuted", "pure Z-dephasing Liouvillians can be strongly non-normal",
            "not basis-invariant inside a degenerate eigenspace"),
            ("K_P ~ 1", "relevant only in gain-loss")),
        "simulations/cockpit_validation.py": ((
            "old K_P ~ 1 pure-dephasing null is refuted",
            "single-vector K_P is meaningful only for a simple isolated mode",
            "invariant-subspace or Jordan diagnostics"),
            ("And it's K_P ~ 1 anyway",)),
        "simulations/results/cockpit_validation.txt": ((
            "old K_P ~ 1 pure-dephasing null is refuted",
            "single-vector K_P is meaningful only for a simple isolated mode",
            "invariant-subspace or Jordan diagnostics"),
            ("And it's K_P ~ 1 anyway",)),
        "review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md": ((
            "blanket Petermann null is refuted", "needs-human / open",
            "single-vector Petermann factors apply only to simple isolated modes"),
            ("confirming eigenvector near-orthogonality",)),
        "experiments/FILLING_THRESHOLD_CHAOS.md": ((
            "separately sampled realization ensembles", "not a within-realization intervention",
            "does not thereby make the full Liouvillian a quadratic free-fermion generator",
            "fixed interacting-disorder distribution, not one fixed disorder realization"),
            ("At the same Liouvillian", "same Liouvillian's dilute-vs-dense",
             "at fixed interacting disorder")),
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdWitness.cs": ((
            "separately sampled realization ensembles", "not a within-realization causal intervention",
            "does not classify the Z-dephasing Liouvillian itself"),
            ("The SAME Liouvillian's",)),
        "simulations/f89_zz_break_gate.py": ((
            "historical locator without a positive-Delta character verdict",
            "strict full-block coincidence/correspondence tolerance",
            "positive-Delta proposals are Uncertified",
            "neither persistence, defectiveness, nor lifting",
            "pair correspondence not independently certified"),
            ("Delta>0 off-axis/defective", "BECOMES DEFECTIVE: breaking free-fermion",
             "FREE-FERMION INTEGRABILITY WAS THE PROTECTION",
             "the EP has MOVED OFF the real axis", "genuine defective Jordan EP",
             "OFF-AXIS EP LOCATOR", "TRUE coalescence")),
        "hypotheses/DIABOLIC_BY_INTEGRABILITY.md": ((
            "historical locator without a positive-Delta character verdict",
            "positive-Delta sweep supplies numerical pair proposals only",
            "strict current verdict is owned by XxzDeltaFlipTests"),
            ("gate-first; Stage 0 reproduces", "The why probes (gate-first")),
    }
    for path, (required_add, forbidden_add) in physics_scope_surfaces.items():
        required, forbidden = surfaces.get(path, ((), ()))
        surfaces[path] = (required + required_add, forbidden + forbidden_add)
    pooled_consumers = ("experiments/RANDOM_MATRIX_THEORY.md", "experiments/README.md",
        "experiments/SPECTRAL_FORM_FACTOR.md", "compute/RCPsiSquared.Compute/README.md",
        "docs/outbound/SHIFTED_ORDER4_CHIRAL_SYMMETRY.md", "docs/KMS_DETAILED_BALANCE.md")
    for path in pooled_consumers:
        required, forbidden = surfaces.get(path, ((), ()))
        surfaces[path] = (required + ("raw-multiset", "zero gaps", "0/0"), forbidden + (
            "0.383", "0.36 to 0.39", "0.36-0.39", "pooled finite-N mean lies near the Poisson reference",
            "sampled statistics lie near the Poisson reference"))
    for path, (required, forbidden) in surfaces.items():
        source = (root / path).read_text(encoding="utf-8")
        baseline = set(task8_text_errors(source, required, forbidden))
        check(f"Round4 live artifact/consumer: {path}", not baseline, ascii("; ".join(sorted(baseline))))
        for phrase in required:
            # Normalize first, so multiline Markdown/C# literals remain mutation-addressable.
            normalized = " ".join(source.replace("**", "").replace("///", "")
                                  .replace('" +', "").replace('"', "").split()).casefold()
            changed = normalized.replace(phrase.casefold(), "MISSING")
            check(f"Round4 missing boundary {ascii(phrase)}: {path}", changed != normalized and
                  bool(set(task8_text_errors(changed, required, forbidden)) - baseline))
        for phrase in forbidden:
            check(f"Round4 forbidden variant {ascii(phrase)}: {path}",
                  bool(set(task8_text_errors(source + " " + phrase, required, forbidden)) - baseline))

    f44_path = "compute/RCPsiSquared.Core/Symmetry/F44AlgebraicPairRateLogIdentityPi2Inheritance.cs"
    f44_source = (root / f44_path).read_text(encoding="utf-8")
    f44_baseline = set(f44_ordered_pair_api_errors(f44_source))
    check("F44 exact and linear APIs enforce finite ordered-pair domain",
          not f44_baseline, ascii("; ".join(sorted(f44_baseline))))
    for name, old, new in (
        ("negative Delta_d accepted", "deltaD < 0.0 ||", ""),
        ("non-finite Delta_d accepted", "!double.IsFinite(deltaD) ||", ""),
        ("linear API bypasses shared domain", "ValidateOrderedPairDomain(deltaD, totalGamma);", "",),
    ):
        changed = f44_source.replace(old, new, 1)
        check(f"F44 domain mutation rejects {name}", changed != f44_source and
              bool(set(f44_ordered_pair_api_errors(changed)) - f44_baseline))

    mirrorworld_path = "compute/MirrorWorld/Formulas.cs"
    mirrorworld_source = (root / mirrorworld_path).read_text(encoding="utf-8")
    mirror_f44_baseline = set(mirrorworld_f44_ordered_pair_api_errors(mirrorworld_source))
    check("MirrorWorld F44 enforces finite ordered-pair domain", not mirror_f44_baseline,
          ascii("; ".join(sorted(mirror_f44_baseline))))
    for name, old, new in (
        ("positive d_fast", "dFast <= 0.0", "dFast < 0.0"),
        ("positive d_slow", "dSlow <= 0.0", "dSlow < 0.0"),
        ("palindromic pair sum", "Math.Abs(dFast / scale + dSlow / scale - 2.0 * (sg / scale))", "Math.Abs(dFast - dSlow)"),
        ("public validation call", "ValidateF44PositivePalindromicPair(dFast, dSlow, sg);", ""),
    ):
        changed = mirrorworld_source.replace(old, new, 1)
        check(f"MirrorWorld F44 mutation rejects missing {name}", changed != mirrorworld_source and
              bool(set(mirrorworld_f44_ordered_pair_api_errors(changed)) - mirror_f44_baseline))

    f95_path = "compute/RCPsiSquared.Core/Symmetry/F95AngleAtQuadraticZeroPi2Inheritance.cs"
    f95_source = (root / f95_path).read_text(encoding="utf-8")
    f95_baseline = set(f95_positive_b_api_errors(f95_source))
    check("F95 typed API enforces finite positive b", not f95_baseline,
          ascii("; ".join(sorted(f95_baseline))))
    changed = f95_source.replace("b <= 0.0", "b == 0.0", 1)
    check("F95 typed mutation rejects negative b", changed != f95_source and
          bool(set(f95_positive_b_api_errors(changed)) - f95_baseline))
    mirror_f95_baseline = set(f95_positive_b_api_errors(mirrorworld_source, mirrorworld=True))
    check("MirrorWorld F95 enforces finite positive b", not mirror_f95_baseline,
          ascii("; ".join(sorted(mirror_f95_baseline))))
    changed = mirrorworld_source.replace("b <= 0.0", "b == 0.0", 1)
    check("MirrorWorld F95 mutation rejects negative b", changed != mirrorworld_source and
          bool(set(f95_positive_b_api_errors(changed, mirrorworld=True)) - mirror_f95_baseline))
    f95_python = (root / "simulations/angle_at_zero_tier1_candidate.py").read_text(encoding="utf-8")
    f95_python_baseline = set(f95_python_positive_b_api_errors(f95_python))
    check("F95 Python producer enforces finite positive b", not f95_python_baseline,
          ascii("; ".join(sorted(f95_python_baseline))))
    changed = f95_python.replace("b <= 0:", "b == 0:", 1)
    check("F95 Python producer mutation rejects negative b", changed != f95_python and
          bool(set(f95_python_positive_b_api_errors(changed)) - f95_python_baseline))

    variants = {
        "compute/RCPsiSquared.Core/Symmetry/F44AlgebraicPairRateLogIdentityPi2Inheritance.cs": (
            "Crooks-like", "β_eff", "EffectiveInverseTemperature", "EmpiricalJarzynskiMean",
            "IsCrooksFluctuationTheorem", "Jarzynski equality fails"),
        "compute/RCPsiSquared.Runtime/PolarityArchitecture/F44AlgebraicPairRateLogIdentityPi2InheritanceRegistration.cs": (
            "F44CrooksLikeRateIdentityPi2Inheritance", "β_eff", "Jarzynski"),
        "compute/MirrorWorld/Formulas.cs": ("F44 (T1, D08): Crooks-like", "β_eff", "Jarzynski"),
        "docs/proofs/derivations/D08_CROOKS_RATE_IDENTITY.md": ("β_eff = 1/Σγ", "effective inverse temperature"),
        "docs/ANALYTICAL_FORMULAS.md": ("F44. Crooks-like", "β_eff = 1/Σγ"),
        "review/OPEN_QUESTIONS_INDEX.md": ("pure dephasing (infinite-temperature bath)",
            "pure Z-dephasing is an infinite-temperature bath"),
        "compute/RCPsiSquared.Core/Symmetry/F95AngleAtQuadraticZeroPi2Inheritance.cs": (
            "Returns NaN for c ≤ 1/4", "if (c <= Threshold)", "if (c <= thresh)"),
        "compute/RCPsiSquared.Core/Symmetry/TransitionBridgeF95SiblingClaim.cs": (
            "NaN at or below the EP", "if (disc >= 0.0)"),
        "compute/RCPsiSquared.Runtime/PolarityArchitecture/TransitionBridgeF95SiblingClaimRegistration.cs": (
            "FRAGILE_BRIDGE exceptional point", "FRAGILE_BRIDGE axis departure is a Hopf"),
        "compute/RCPsiSquared.Diagnostics/Foundation/BranchLocusPalindromeWitness.cs": ("The sampled XXZ Delta control stays on-line yet defects",),
        "compute/RCPsiSquared.Core/Symmetry/F89BranchLocusPalindromeClaim.cs": ("control stays on-line yet becomes defective",),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": ("the certified N=4 XXZ Delta control defects",),
        "simulations/results/pt_symmetry_analysis.txt": ("Class AIII", "chiral breaking", "No classical EP"),
        "simulations/pt_palindrome_breaking.py": ("chiral phase", "still a HOPF bifurcation", "changed to SADDLE-NODE", "palindrome IS the protection mechanism"),
        "simulations/results/pt_palindrome_breaking.txt": ("chiral phase", "still a HOPF bifurcation", "changed to SADDLE-NODE", "palindrome IS the protection mechanism"),
        "experiments/PT_SYMMETRY_ANALYSIS.md": ("The full irreducible SRP class is AIII",
            "rules out the prior Hopf", "is a Hopf bifurcation", "is not a Hopf"),
        "hypotheses/FRAGILE_BRIDGE.md": ("not a Hopf", "is a Hopf bifurcation", "genuine EP at real γ_crit"),
        "experiments/F86_EP_THROUGH_THE_CLOCK.md": ("is a Hopf bifurcation", "SEPARATE genuine EP",
            "two SEPARATE genuine EPs", "the dynamics escape into self-sustained oscillation", "Hopf can push"),
        "simulations/memory_fate_across_the_takt.py": ("real gain-loss Hopf", "Hopf / runaway",
            "the Hopf, the runaway"),
        "simulations/entropy_production.py": ("maximum entropy production", "No net entropy production", "Carnot"),
        "simulations/results/entropy_production.txt": ("maximum entropy production", "No net entropy production", "Carnot"),
        "experiments/ENTROPY_PRODUCTION.md": ("The rate sum is a thermodynamic entropy-production scale",),
        "experiments/README.md": ("The rate sum is a thermodynamic entropy-production scale", "N=4 positive-Delta proposals are certified defective", "+0.988", "stabilizes the sampled system"),
        "simulations/rmt_analysis.py": ("def unfold", "NNSD", "classify_r"),
        "simulations/results/rmt_analysis.txt": ("NNSD", "Class (all)", "Nearest: Poisson"),
        "experiments/RANDOM_MATRIX_THEORY.md": ("rank unfolding supplies a valid NNSD", "spacing ratios are invariant under nonlinear maps"),
        "docs/WHAT_WE_FOUND.md": ("interacting eigenmode rates form a 2γ-spaced ladder",),
        "experiments/F89_BRANCH_LOCUS_PALINDROME.md": ("N=4 positive-Delta proposals are certified defective",),
        "experiments/F89_PATH_K_DIABOLIC.md": ("N=4 positive-Delta proposals are certified defective",),
        "review/EMERGING_QUESTIONS.md": ("Π class AIII chiral", "K_1 = Π on H_1 sine-mode basis"),
        "hypotheses/DIABOLIC_BY_INTEGRABILITY.md": ("N=4 positive-Delta proposals are certified defective",),
        "hypotheses/README.md": ("N=4 positive-Delta proposals are certified defective",),
    }
    for path, old_claims in variants.items():
        source = (root / path).read_text(encoding="utf-8")
        required, forbidden = surfaces.get(path, ((), ()))
        baseline = set(task8_text_errors(source, required, forbidden))
        for old in old_claims:
            check(f"Round4 consumer/artifact rejects {ascii(old)}: {path}",
                  bool(set(task8_text_errors(source + " " + old, required, forbidden)) - baseline))
    for name, markers in (
        ("pt_symmetry_analysis", ("full irreducible SRP class remains OPEN", "off palindrome axis; operator relation remains exact")),
        ("pt_palindrome_breaking", ("one-parameter composite perturbation sweep", "does not classify the threshold as Hopf, EP, or a Jordan defect")),
        ("entropy_production", ("chosen I/d reference", "not a thermodynamic entropy-production scale", "not a physical Lindblad channel")),
        ("rmt_analysis", ("affine-invariant", "no irreducible-class inference")),
        ("spectral_form_factor", ("Only the raw oscillation-frequency SFF is computed", "no decay-weighted or connected estimator is defined")),
        ("fragile_bridge_ep_signature", ("spectral-abscissa axis departure", "not branch continuation", "A single-vector K is omitted wherever the selected eigenvalue is degenerate")),
    ):
        producer_path = root / "simulations" / (name + ".py")
        artifact_path = root / "simulations/results" / (name + ".txt")
        producer = producer_path.read_text(encoding="utf-8")
        artifact = artifact_path.read_text(encoding="utf-8")
        path_baseline = set(producer_artifact_path_errors(producer, producer_path, artifact_path))
        check(f"Round4 {name} producer writes canonical artifact path", not path_baseline,
              ascii("; ".join(sorted(path_baseline))))
        output_literal = f'"{name}.txt"'
        literal_index = producer.rfind(output_literal)
        changed_path = (producer[:literal_index] + f'"{name}_moved.txt"' +
                        producer[literal_index + len(output_literal):]) if literal_index >= 0 else producer
        check(f"Round4 output-path mutation rejects {name}", changed_path != producer and
              bool(set(producer_artifact_path_errors(changed_path, producer_path, artifact_path)) - path_baseline))
        baseline = set(emitted_artifact_claim_errors(producer, artifact, markers))
        check(f"Round4 protected emitted sentences fresh in {name} artifact", not baseline, ascii("; ".join(sorted(baseline))))
        for marker in markers:
            changed_artifact = artifact.replace(marker, "STALE_OUTPUT")
            changed_producer = producer.replace(marker, marker + " NEW_SCOPE_NOT_IN_ARTIFACT")
            check(f"Round4 artifact-only mutation rejects {name}: {ascii(marker)}", changed_artifact != artifact and
                  bool(set(emitted_artifact_claim_errors(producer, changed_artifact, markers)) - baseline))
            check(f"Round4 producer-only mutation requires regenerated {name}: {ascii(marker)}", changed_producer != producer and
                  bool(set(emitted_artifact_claim_errors(changed_producer, artifact, markers)) - baseline))
        if name == "pt_palindrome_breaking":
            matcher_source = (root / "simulations/pt_multiset_matching.py").read_text(encoding="utf-8")
            integrity_baseline = set(phase3_artifact_integrity_errors(producer, artifact, matcher_source))
            check("Round4 Phase-3 deterministic artifact binds source and payload", not integrity_baseline,
                  ascii("; ".join(sorted(integrity_baseline))))
            changed_source = producer + "\n# SOURCE_MUTATION"
            changed_numeric = artifact.replace("0.1873101", "0.9999999", 1)
            check("Round4 Phase-3 source mutation requires regenerated artifact",
                  bool(set(phase3_artifact_integrity_errors(changed_source, artifact, matcher_source)) - integrity_baseline))
            check("Round4 Phase-3 matcher mutation requires regenerated artifact",
                  bool(set(phase3_artifact_integrity_errors(producer, artifact, matcher_source + "\n# MATCHER_MUTATION")) - integrity_baseline))
            check("Round4 Phase-3 numerical artifact mutation is rejected", changed_numeric != artifact and
                  bool(set(phase3_artifact_integrity_errors(producer, changed_numeric, matcher_source)) - integrity_baseline))
    # Recompute from exported eigenvalues, independently of the producer's spacing helper.
    inputs = {}
    for n in range(2, 8):
        lines = (root / f"simulations/results/rmt_eigenvalues_N{n}.csv").read_text(encoding="utf-8").splitlines()[1:]
        pairs = [[float(x.replace(",", ".")) for x in line.split("\t")] for line in lines if line.strip()]
        inputs[n] = np.array([real + 1j * imaginary for real, imaginary in pairs])
    artifact = (root / "simulations/results/rmt_analysis.txt").read_text(encoding="utf-8")
    baseline = set(rmt_artifact_row_errors(artifact, inputs))
    check("Round4 RMT artifact rows agree with independent CSV calculation", not baseline, ascii("; ".join(baseline)))
    for name, changed in (
        ("wrong N3 ratio", artifact.replace("<r> defined = 0.22260363264213059", "<r> defined = 0.7777", 1)),
        ("zero gaps deleted", artifact.replace("zero_gaps=1340", "zero_gaps=0", 1)),
        ("undefined counted as defined", artifact.replace("defined_ratios=3107, undefined_zero_zero=987", "defined_ratios=4094, undefined_zero_zero=0", 1)),
        ("lower-half population drift", artifact.replace("levels=1973", "levels=1900", 1)),
        ("wrong N7 population", artifact.replace("N=7: 16384 eigenvalues", "N=7: 100 eigenvalues", 1)),
        ("missing N6 row", artifact.replace("N=6: 4096 eigenvalues", "MISSING N6", 1)),
    ):
        check(f"Round4 RMT CSV provenance mutation rejects {name}", changed != artifact and
              bool(set(rmt_artifact_row_errors(changed, inputs)) - baseline))
    band_artifact = (root / "simulations/results/rmt_band_multiset.txt").read_text(encoding="utf-8")
    baseline = set(rmt_band_artifact_errors(band_artifact, inputs))
    check("RMT band artifact agrees with independent CSV calculation", not baseline, ascii("; ".join(baseline)))
    for name, old, new in (("deduplicated band", "levels=1710", "levels=414"),
        ("rounded-away zero gap", "zero_gaps=17", "zero_gaps=0"),
        ("wrong band ratio", "0.19686572210722628", "0.283"),
        ("missing endpoint band", "w=0:", "MISSING:")):
        changed = band_artifact.replace(old, new, 1)
        check(f"RMT band CSV mutation rejects {name}", changed != band_artifact and
              bool(set(rmt_band_artifact_errors(changed, inputs)) - baseline))
    document = (root / "experiments/RANDOM_MATRIX_THEORY.md").read_text(encoding="utf-8")
    baseline = set(rmt_document_table_errors(document, artifact, band_artifact))
    check("RMT document tables consume executed artifact rows", not baseline, ascii("; ".join(baseline)))
    for old,new in (("| 16384 | 136 | 16366 | 16 |", "| 16384 | 136 | 16382 | 0 |"),
        ("0.2021456120489688 |", "0.383 |"), ("0.2021456120489688 |", "nan |"),
        ("| 7 | 3 | 1710 |", "| 7 | 3 | 414 |")):
        changed = document.replace(old,new,1)
        check(f"RMT document numeric mutation rejects {old}", changed != document and
              bool(set(rmt_document_table_errors(changed,artifact,band_artifact))-baseline))


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
                     "ROUTE_B_N6_A2_LOCUS_ATLAS.md", "their 96 conjugation/parity orbits",
                     "| even | 133 | 59 | 74 | 0 | 133 Diabolic, alg=geo=2 |",
                     "| odd | 133 | 59 | 74 | 0 | 133 Diabolic, alg=geo=2 |",
                     "| total | 266 | 118 | 148 | 0 | 266 Diabolic; 0 Defective |",
                     "EpCharacterStable records stable numerical character from the isolated two-dimensional restriction; it is not an exact rank proof",
                     "The semisimplicity itself does not rest on those readings"]
    elif label == "arc":
        required += ["ExactAlgebraic owns all 266 loci (133 per parity), consumed exactly once with 0 unresolved",
                     "character is Diabolic with alg=geo=2 at every locus",
                     "118 by executed full-sector HermitianAxis and 148 by stable EpCharacterStable readings, with ExactRankExecuted=0",
                     "The numerical rule uses all three intrinsic radii", "No exact fallback was used",
                     "neither determines Jordan character",
                     "F163 (typed as RouteBN6A2UnfoldingClaim, live at inspect --root n6unfolding) independently proves all 266 semisimple"]
    else:
        required += ["consumes 266 distinct loci once", "ExactAlgebraic=266",
                     "HermitianAxis=118, EpCharacterStable=148, ExactRankExecuted=0",
                     "parity contributes 59 HermitianAxis and 74 EpCharacterStable readings",
                     "All 266 are Diabolic with alg=geo=2; 0 are Defective",
                     "not an exact rank proof", "exactRankCertificates array is empty",
                     "Semisimplicity itself is settled algebraically for all 266"]
        required += (["133 per parity", "EpCharacterStable is numerical evidence"] if label == "path experiment"
                     else ["133 in each parity", "is stable numerical evidence"])
    # The exact algebraic/structural sources do not license a positive exact-rank
    # claim for numerical character. Negated rank-proof wording is the scope fence.
    errors = [f"missing {phrase}" for phrase in required if phrase not in source]
    if "an exact rank proof" in source.casefold().replace("not an exact rank proof", ""):
        errors.append("numerical character promoted to an exact rank proof")
    return errors


def _route_b_exact_integer(value):
    """Parse an integer JSON atom without accepting bools, floats, or decimal strings."""
    if type(value) is int:
        return value
    if isinstance(value, str) and re.fullmatch(r"[+-]?\d+", value):
        return int(value)
    raise ValueError("not an exact integer")


def route_b_n6_real_q_exclusion_errors(inventory):
    """Validate the schema-3 N=6 certificate and require every exact t-box to exclude Re(t)=0."""
    errors = []
    if not isinstance(inventory, dict):
        return ["N6 inventory root must be an object"]
    if type(inventory.get("schemaVersion")) is not int or inventory.get("schemaVersion") != 3:
        errors.append("N6 inventory must have integer schemaVersion=3")
    if type(inventory.get("n")) is not int or inventory.get("n") != 6:
        errors.append("N6 inventory must have integer n=6")
    if inventory.get("model") != "open-uniform-XY-delta0-field0-gamma1-SEket-DEbra":
        errors.append("N6 inventory must use the canonical model")
    conventions = inventory.get("conventions")
    if (not isinstance(conventions, dict)
            or conventions.get("parameter") != "t=i*qCSharp;qCSharp=-i*t"
            or conventions.get("eigenvalue") != "Lambda=2*lambda"
            or conventions.get("coefficientOrder") != "lambda-lowest-first,t-lowest-first"):
        errors.append("N6 inventory must use the canonical t/q/lambda conventions")

    degrees = inventory.get("a2Degrees")
    if (not isinstance(degrees, dict)
            or type(degrees.get("E")) is not int or degrees.get("E") != 133
            or type(degrees.get("O")) is not int or degrees.get("O") != 133):
        errors.append("N6 inventory must declare integer A2 degrees E=133 and O=133")

    loci = inventory.get("loci")
    if not isinstance(loci, list) or len(loci) != 266:
        errors.append("N6 inventory must contain exactly 266 certified loci")
        return errors

    ids = set()
    parity_counts = {"E": 0, "O": 0}
    for locus in loci:
        locus_id = locus.get("id", "<missing-id>") if isinstance(locus, dict) else "<invalid-locus>"
        if not isinstance(locus, dict):
            errors.append("<invalid-locus>: locus must be an object")
            continue
        if not isinstance(locus_id, str) or not locus_id.strip() or locus_id in ids:
            errors.append(f"{locus_id}: locus id must be nonempty and unique")
        else:
            ids.add(locus_id)
        parity = locus.get("parity")
        if parity not in parity_counts:
            errors.append(f"{locus_id}: parity must be E or O")
        else:
            parity_counts[parity] += 1
        if type(locus.get("algebraicMultiplicity")) is not int or locus.get("algebraicMultiplicity") != 2:
            errors.append(f"{locus_id}: certified A2 locus must have integer algebraicMultiplicity=2")
        try:
            interval = locus["tBox"]["real"]
            lo = interval["lower"]
            hi = interval["upper"]
            lo_num = _route_b_exact_integer(lo["numerator"])
            lo_den = _route_b_exact_integer(lo["denominator"])
            hi_num = _route_b_exact_integer(hi["numerator"])
            hi_den = _route_b_exact_integer(hi["denominator"])
            if lo_den <= 0 or hi_den <= 0:
                errors.append(f"{locus_id}: nonpositive rational denominator")
            elif lo_num * hi_den > hi_num * lo_den:
                errors.append(f"{locus_id}: inverted tBox.real interval")
            elif lo_num <= 0 <= hi_num:
                errors.append(f"{locus_id}: tBox.real does not exclude zero")
        except (KeyError, TypeError, ValueError):
            errors.append(f"{locus_id}: malformed exact tBox.real interval")
    if parity_counts != {"E": 133, "O": 133}:
        errors.append(f"N6 inventory parity counts must be E=133 and O=133, got {parity_counts}")
    semantic_digest = hashlib.sha256(json.dumps(
        inventory, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")).hexdigest()
    if semantic_digest != "6ae3c0fb76d3415177d8a9bec932ba7abfb9bf9be42b8a341e931cb672532eda":
        errors.append("N6 boxes are detached from the semantically certified canonical carrier")
    return errors


READING_POWER_MONOTONICITY_TEST = (
    "ReadingPowerWitnessTests.FiIsMonotoneOnTheSevenPointSampledQGrid_InEveryBasis"
)


def verify_final_five_repairs():
    """Guard the final F44/F86/flow/ledger/CSR repairs with same-door mutations."""
    root = Path(__file__).resolve().parents[1]
    text_surfaces = {
        "docs/F86_VALUES_INVENTORY.md": (
            ("SE-walk population handover", "Q_label", "Q_Lindblad = 2 Q_label",
             "spectral character remains open", "not an F86-block EP anchor"),
            ("Q_EP got its hardware anchor", "rotation born at the F86a exceptional point",
             "Q_EP onset (hardware)", "post-EP dynamics", "reborn mode")),
        "experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md": (
            ("sampled SE-walk population handover", "Q_label", "Q_Lindblad = 2 Q_label",
             "spectral character remains open", "fixed positive per-site profile",
             "not a physical CP dephasing channel", "not an isolated global sink",
             "global slowest non-kernel mode is not identified with the population-visible approach rate",
             "0.28 → 0.84", "0.8417853730254796"),
            ("two modes of L collide and merge", "is a Hopf bifurcation",
             "At Q_EP the two coalesce defectively", "verified facts are the two singularity types",
             "reborn mode", "post-EP dynamics", "equipartition floor", "natural Q≈30",
             "populations converge toward", "simple λ=0 kernel")),
        "README.md": (
            ("Q_label = 1.5→2.5", "Q_Lindblad = 3→5", "population handover",
             "spectral character open"),
            ("across the EP", "EP onset: revival")),
        "data/ibm_ep_onset_may2026/README.md": (
            ("historical run identifier", "Q_label", "Q_Lindblad = 2 Q_label",
             "population handover", "spectral character remains open", "0.28 → 0.84",
             "sum 1.103", "not a normalized", "0.8417853730254796"),
            ("the EP onset on a real chip", "reborn memory", "dephasing EP onset",
             "equipartition floor", "populations converge toward", "floor and the onset are clean")),
        "docs/Q_REGIME_ANCHORS.md": (
            ("not a spectral Q anchor", "Q_label=1.5→2.5", "Q_Lindblad=3→5",
             "spectral character remains open"),
            ("single-excitation flow exceptional point", "Birth Canal's birth singularity",
             "chosen to sit on the flow EP")),
        "compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs": (
            ("Q_label", "Q_Lindblad = 2 Q_label", "spectral character remains open",
             "population handover", "0.28 → 0.84", "sum 1.103",
             "does not isolate a local Z-dephasing rate", "0.8417853730254796"),
            ("equipartition floor", "Confirms the typed UniversalCarrierClaim.DefaultGammaZero", "hardware-anchored",
             "single-excitation walk's critical-damping transition", "the reborn memory",
             "ExceptionalPointClock (the toy 2×2 reduction")),
        "simulations/framework/confirmations.py": (
            ("Q_label_grid", "Q_lindblad_grid", "Q_Lindblad = 2 Q_label",
             "spectral character remains open", "population handover", "0.28 → 0.84",
             "sum 1.103", "does not isolate a local Z-dephasing rate", "0.8417853730254796"),
            ("equipartition floor", "Confirms the typed", "hardware-anchored",
             "floor_below_ep", "liftoff_above_ep", "critical-damping transition",
             "the reborn memory", "ExceptionalPointClock (the toy 2×2 reduction")),
        "compute/RCPsiSquared.Diagnostics/Foundation/EpField.cs": (
            ("Q_Lindblad = 2 Q_label", "spectral character remains open",
              "population handover", "F89's scattered seeds"),
            ("at the same Q_EP", "post-EP regime it opens", "critical-damping",
             "has NO defective EP there (eigenvalues simple)")),
        "compute/RCPsiSquared.Diagnostics/Foundation/ExceptionalPointClock.cs": (
            ("not a hardware calibration",),
            ("hardware-anchored", "BirthCanal of")),
        "compute/RCPsiSquared.Diagnostics/Foundation/PostEpFlowField.cs": (
            ("single-excitation population flow", "no spectral-transition claim",
              "H = (Q/2)·H_unit", "new Complex(q / 2.0, 0.0)",
              "strictly > 0", "not a flow-overlap rate",
              "HasBirthCanalClassification", "N=5 open-chain-only birth/sterile surface",
              "N == 5 && Topology == FlowTopology.Chain",
              "slow.Parity", "SlowClusterMeanRate", "AbsorptionResidual"),
            ("The post-EP flow as", "post-EP flow is", "every trajectory relaxes",
             "Always available, never throws", "depth mod 2", "equal to the rate, bit-exact")),
        "compute/RCPsiSquared.Cli/Commands/AssemblyCommand.cs": (
            ("1.0,2.0,3.0,5.0,10.0,40.0", "odd support", "even support",
             "mixed support", "global slow-rate drift", "SlowClusterMeanRate", "abs(mean)"),
            ("0.5,1.0,1.5,2.5,5.0,20.0", "BIRTH(odd)", "flow(even)",
             "depth (= light n_XY = rate)")),
        "compute/RCPsiSquared.Cli/Program.cs": (
            ("--q-list 1,2,3,5,10,40", "weighted Absorption rate",
             "projector parity support"),
            ("--q-list 0.5,1,1.5,2.5,5,20", "depth=light=rate")),
        "compute/RCPsiSquared.Diagnostics/Foundation/SectorReductionWitness.cs": (
            ("public Q is the canonical carrier Q", "H=(Q/2)",
             "new ChainSystem(n, canonicalQ", "canonical Q=3", "canonical Q=2000"),
            ("new ChainSystem(n, 2.0 *", "rate(Q=1000)", "rate(Q=1.5)")),
        "compute/RCPsiSquared.Diagnostics/Foundation/VacuumBlockReductionClaim.cs": (
            ("L_(1,0) = -i(Q/2)", "canonical Q=3", "Q=2000"),
            ("-iQh", "Q=1.5 ->", "Q=1000 ->")),
        "compute/RCPsiSquared.Diagnostics/Foundation/BirthCanalSurfaceWitness.cs": (
            ("rate(Q=2000) - rate(Q=3)", "light (Q=3)", "rate(Q=3)",
             "lo.AbsorptionResidual", "cluster-mean rate", "probes rate-sterility vs rate drift",
             "Deviation is the spectral-edge rate at Q=2000 minus",
             "selected tolerance-cluster's mean rate to 2·Σ_l γ_l·light_l",
             "light-freeze verdict requires", "NumericalConditioningSuspect",
             "only a numerical-conditioning diagnostic", "not an EP certificate"),
            ("rate(Q=1000)", "rate(Q=1.5)", "light (Q=1.5)", "reads the light-freeze",
             "): rate = 2·Σ_l", "EP detector", "EpSuspect")),
        "compute/RCPsiSquared.Diagnostics/Ptf/SlowLightDistribution.cs": (
            ("canonical Q=2000", "rate at Q=3", "at Q=3", "at Q=2000",
             "rate-sterility", "not sufficient", "N=5", "EvenParityFraction",
             "OddParityFraction", "SlowLightParity", "SpectralEdgeRate",
             "double ClusterMeanRate,", "ClusterDimension",
             "slow tolerance cluster. Membership does not assert exact rate degeneracy."),
            ("rate at Q=1000", "rate at Q=1.5",
             "sterile ⇔ this distribution is", "at N ≤ 6", "the slowest rate is")),
        "compute/RCPsiSquared.Diagnostics.Tests/Ptf/SlowLightDistributionTests.cs": (
            ("ToleranceCluster_SeparatesSpectralEdgeFromClusterMeanAbsorptionRate",
             "1.00000025", "ClusterMeanRate", "AbsorptionRate", "ClusterDimension"),
            ()),
        "compute/RCPsiSquared.Diagnostics/Foundation/TrichotomyWitness.cs": (
            ("Public canonical convention", "new ChainSystem(n, canonicalQ",
             "BirthCanalProbeQLow", "BirthCanalProbeQHigh", "CarbonHamiltonianQ = 1.0"),
            ("new ChainSystem(n, 2.0 *", "const double qLo = 1.5", "const double qHi = 1000")),
        "simulations/gamma_profile_isolation.py": (
            ("distinct finite-scan observables", "without using rate",
             "drift as an isolation classifier", "Q / 2.0", "rates(N, 3.0",
             "rates(N, 2000.0", "No isolation verdict"),
            ("rate(Q=1000) - rate(Q=1.5)", "-1j * Q *",
             "reveals spectral isolation", "0 = isolated/shortcut",
             "gap is the isolation itself", "isolated (shortcut)")),
        "simulations/gamma_profile_shortcut.py": (
            ("compares orthogonal projectors", "no exact equivalence", "Q / 2.0",
             "subspace_overlap", "commutator_residual", "not an iff theorem"),
            ("v(Q=1.5)", "v(Q=1000)", "-1j * Q *", "exactly when",
             "SHARED eigenvector", "slowest-mode overlap")),
        "simulations/light_content.py": (
            ("orthogonal projector", "CLUSTER-MEAN decay rate",
             "spectral-edge rate is reported separately", "Q / 2.0",
             "for Q in [3.0, 2000.0]", "slow_subspace", "projector_light",
             "max per-site distribution drift"),
            ("Q=1.5 and Q=1000", "-1j * Q *",
             "for the slowest mode, decompose its eigenvector", "FROZEN (sterile)",
             "Hamiltonian mixes light in -> birth canal")),
        "simulations/tests/test_task8_birth_canal_repairs.py": (
            ("test_tolerance_cluster_separates_edge_from_projector_mean",
             "assert mean != edge", "assert np.isclose(absorption, mean",
             "test_shortcut_diagnostics_are_invariant_under_subspace_basis_rotation",
             "shortcut.subspace_overlap", "shortcut.commutator_residual"),
            ()),
        "simulations/birth_canal_junction_nature.py": (
            ("rate(Q=2000) - rate(Q=3)", "Q / 2.0", "QLO, QHI, DEV_TOL = 3.0, 2000.0",
             "public canonical convention H=(Q/2)*H_unit"),
            ("rate(Q=1000) - rate(Q=1.5)", "-1j * Q *",
             "Q*H_unit + per-site profile dephasing")),
        "simulations/birth_canal_boundary_pathdependence.py": (
            ("rate@Q=3", "Q / 2.0", "slowest(N, 2000.0", "slowest(N, 3.0"),
            ("rate@1.5", "-1j * Q *")),
        "simulations/birth_canal_n6_mode_crossing.py": (
            ("Q / 2.0", "for Q in (3.0, 2000.0)"),
            ("for Q in (1.5, 1000.0)", "-1j * Q *")),
        "simulations/birth_canal_vacuum_block_verifier.py": (
            ("PROBE_LO, PROBE_HI = 3.0, 2000.0", "Q / 2.0", "-i(Q/2)"),
            ("PROBE_LO, PROBE_HI = 1.5, 1000.0", "-1j * Q *", "-iQh")),
        "compute/RCPsiSquared.Cli/Commands/InspectCommand.cs": (
            ("single-excitation population flow",),
            ("the post-EP single-excitation flow",)),
        "simulations/ep_transition.py": (
            ("population_scan()", "1/N reference level",
             "spectral character remains open", "population handover"),
            ("matching the hardware", "memory switches on as Q crosses the EP",
             "rotation born on a real chip", "confirmed on IBM Kingston")),
        "simulations/journey_between_singularities.py": (
            ("map(np.asarray, population_scan())", "H += (Q / 2.0)",
             "spectral character remains open", "rate-book corrected; spectral character remains open",
             "1/N reference"),
            ("reborn memory", "the chip sits on the EP", "EP is reached only by INJECTING noise")),
        "simulations/journey_control.py": (
            ("Illustrative dephasing control", "not a calibration of the Kingston runner",
             "no hardware or EP calibration", "control_domain(starts)",
             "tau = gamma_per_site * t"),
            ("Pushing the chip onto the EP", "birth-axis", "newborn rotation", "gamma_total")),
        "simulations/post_ep_dynamics_4d.py": (
            ("canonical Q = J/gamma", "Q / 2.0", "not a critical-damping or spectral-transition verdict"),
            ("post-EP dynamics", "underdamped (hops, remembers)")),
        "simulations/at_the_target.py": (
            ("common per-site profile scalar f", "not asserted to be population-visible",
             "non-CP inverse dephasing", "Q / 2.0"),
            ("Sigma-gamma (f)", "the Hopf", "SINK (all modes flow in)")),
        "simulations/the_flow_endpoints.py": (
            ("N+1-fold semisimple", "Q / 2.0", "no branch continuation"),
            ("lambda=0 SIMPLE", "TWO SINGULAR ENDPOINTS")),
        "docs/NAVIGATING_THE_DIMENSIONS.md": (
            ("Q_label=J/Γ", "canonical Q_Lindblad=3→5", "spectral character open",
             "no branch continuation or doorway"),
            ("1/N equipartition floor", "overdamped→revival handover", "post-EP flow")),
        "experiments/COUPLING_DEFECT_WALK_TIME_STEP.md": (
            ("Q_label=1.5→2.5", "canonical Q_Lindblad=3→5",
             "populations supply no Q*, critical-damping, or EP certificate"),
            ("two hardware-confirmed handovers",)),
        "experiments/GAMMA0_IS_ALWAYS_THERE.md": (
            ("finite-grid transfer-overshoot bracket", "not locate a critical point",
             "does not by itself identify a single local Z-dephasing coefficient"),
            ("reading the carrier off its only lever", "threshold is γ₀, read off",
             "chip telling us its true γ₀")),
        "compute/RCPsiSquared.Core/F86/PolarityPairQPeakDecompositionClaim.cs": (
            ("Illustrative γ₀ value", "finite-time", "does not isolate a dephasing rate",
             "unresolved finite-grid/finite-N bracket"),
            ("Hardware anchor for the Q", "critical-damping threshold", "hardware-anchored axis")),
        "experiments/RELAY_PROTOCOL.md": (
            ("baseline γ remains a model choice", "finite-time transfer overshoot",
             "no critical damping", "no calibrated error model", "does not measure γ₀"),
            ("turned out to be a measurement", "Trotterization-limited", "post-EP/reborn mode")),
        "experiments/README.md": (
            ("finite-time transfer overshoot", "does not measure γ₀", "certify critical damping"),
            ("Reading the carrier Q = J/γ₀ off its only lever",)),
        "docs/Q_BELONGS_TO_NO_SUBSTANCE.md": (
            ("gamma_book_enforcement_boundary", "finite-grid transfer-overshoot bracket",
             "does not isolate γ₀", "not a hardware calibration"),
            ("gamma_book_enforced_nowhere", "hardware read-off of γ₀")),
        "docs/proofs/PROOF_DIFFUSION_RAYLEIGH_CLOSURE.md": (
            ("chosen repository model point", "not a hardware calibration",
             "finite-grid population-overshoot change"),
            ("canonical hardware operating point", "hardware dephasing rate")),
        "experiments/CONCENTRATOR_AB_MECHANISM_TEST.md": (
            ("illustrative default", "named hardware records do not calibrate it"),
            ("critical damping exactly", "hardware-confirmed")),
        "docs/outbound/SELECTIVE_DECOUPLING_SELECTION_RULE.md": (
            ("illustrative carrier default", "not a hardware measurement"),
            ("hardware-measured carrier rate",)),
        "experiments/CORNER_BEAT_HARDWARE_PREDICTION.md": (
            ("illustrative model point", "gamma_book_enforcement_boundary",
             "remaining legacy/future-ingestion schema boundary"),
            ("canonical hardware regime", "canonical hardware-anchored point",
             "gamma_book_enforced_nowhere")),
        "experiments/ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION.md": (
            ("no dedicated `absorption_rung_ladder` arc exists",
             "gamma_book_enforcement_boundary", "not registered anywhere"),
            ("gamma_book_enforced_nowhere",
             "registered in the arc `absorption_rung_ladder`")),
        "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs": (
            ("gamma_book_enforcement_boundary",
             "corner beat's compliance with gamma_book_enforcement_boundary"),
            ("gamma_book_enforced_nowhere",)),
        "simulations/flight/README.md": (
            ("gamma_book_enforcement_boundary", "`corner_beat` arc now exists",
             "current-truth registry pointers", "2026-08-18 record"),
            ("gamma_book_enforced_nowhere", "does not exist yet")),
        "docs/CAUGHT_ERRORS.md": (
            ("K≈403 is a finite-grid, finite-offset simple-mode non-normality reading",
             "FRAGILE_BRIDGE measurement is a spectral-abscissa axis departure",
             "sectorwise involutive P-type anticommutation",
             "does not assign a global irreducible SRP/AIII class",
             "Q_Lindblad = 2 Q_label", "Populations alone certify no critical damping",
             "H=(Q/2)Σ(XX+YY)", "probe Q=3/2000", "0.8417853730254796",
             "finite-grid transfer-overshoot bracket", "distinct decompositions"),
            ()),
    }
    for path, (required, forbidden) in text_surfaces.items():
        source = (root / path).read_text(encoding="utf-8")
        baseline = set(task8_text_errors(source, required, forbidden))
        check(f"Final repair current truth: {path}", not baseline,
              ascii("; ".join(sorted(baseline))))
        for phrase in required:
            changed = source.replace(phrase, "")
            check(f"Final repair mutation rejects missing {ascii(phrase)}: {path}",
                  changed != source and
                  bool(set(task8_text_errors(changed, required, forbidden)) - baseline))
        for phrase in forbidden:
            changed = source + "\n" + phrase
            check(f"Final repair mutation rejects {ascii(phrase)}: {path}",
                  bool(set(task8_text_errors(changed, required, forbidden)) - baseline))

    # The arc names above are also prose, so substring presence cannot certify the
    # registry topology: a description may still mention a name after its actual
    # `Name:` field has disappeared. Parse the fields and gate the three relevant
    # membership claims directly.
    arc_registry_path = "compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs"
    arc_registry_source = (root / arc_registry_path).read_text(encoding="utf-8")

    def task8_arc_registry_errors(source):
        names = re.findall(r'\bName:\s*"([^"]+)"', source)
        errors = []
        for name in ("gamma_book_enforcement_boundary", "corner_beat"):
            if names.count(name) != 1:
                errors.append(f"registry must contain exactly one Name field for {name}")
        for name in ("absorption_rung_ladder", "gamma_book_enforced_nowhere"):
            if name in names:
                errors.append(f"registry must not contain a Name field for {name}")
        return errors

    arc_baseline = set(task8_arc_registry_errors(arc_registry_source))
    check("Final repair exact OpenArc registry topology", not arc_baseline,
          ascii("; ".join(sorted(arc_baseline))))
    for name in ("gamma_book_enforcement_boundary", "corner_beat"):
        changed = arc_registry_source.replace(
            f'Name: "{name}"', f'Name: "{name}_MUTATED"', 1)
        check(f"Final repair registry mutation rejects renamed {name}",
              changed != arc_registry_source and
              bool(set(task8_arc_registry_errors(changed)) - arc_baseline))
    for name in ("absorption_rung_ladder", "gamma_book_enforced_nowhere"):
        changed = arc_registry_source + f'\nName: "{name}"\n'
        check(f"Final repair registry mutation rejects inserted {name}",
              bool(set(task8_arc_registry_errors(changed)) - arc_baseline))

    # Sweep every active consumer involved in this repair. The generated menu and
    # test breadcrumb were the two places a registry-only check previously missed.
    arc_surface_requirements = {
        "compute/RCPsiSquared.Core.Tests/Inspection/OpenArcsInspectableNodeTests.cs":
            ("gamma_book_enforcement_boundary", "corner_beat"),
        "simulations/results/rescued_arcs_menu.txt":
            ("gamma_book_enforcement_boundary", "corner_beat"),
        "docs/Q_BELONGS_TO_NO_SUBSTANCE.md":
            ("gamma_book_enforcement_boundary",),
        "experiments/CORNER_BEAT_HARDWARE_PREDICTION.md":
            ("gamma_book_enforcement_boundary",),
        "experiments/ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION.md":
            ("gamma_book_enforcement_boundary",),
        "simulations/flight/README.md":
            ("gamma_book_enforcement_boundary", "corner_beat"),
    }
    deprecated_arc_name = "gamma_book_enforced_nowhere"

    def task8_arc_surface_errors(source, required_names):
        errors = []
        for name in required_names:
            if name not in source:
                errors.append(f"missing current arc name {name}")
        if deprecated_arc_name in source:
            errors.append(f"deprecated arc name remains: {deprecated_arc_name}")
        return errors

    for path, required_names in arc_surface_requirements.items():
        source = (root / path).read_text(encoding="utf-8")
        baseline = set(task8_arc_surface_errors(source, required_names))
        check(f"Final repair OpenArc propagation: {path}", not baseline,
              ascii("; ".join(sorted(baseline))))
        for name in required_names:
            changed = source.replace(name, "")
            check(f"Final repair OpenArc mutation rejects missing {name}: {path}",
                  changed != source and
                  bool(set(task8_arc_surface_errors(changed, required_names)) - baseline))
        changed = source + "\n" + deprecated_arc_name
        check(f"Final repair OpenArc mutation rejects deprecated name: {path}",
              bool(set(task8_arc_surface_errors(changed, required_names)) - baseline))

    # Raw records are the numerical source. Recompute the observable and rate-book map,
    # and mutate each mechanism so this gate demonstrably fails.
    data_dir = root / "data/ibm_ep_onset_may2026"
    hw_record = json.loads((data_dir / "ep_onset_hardware_ep_ibm_kingston_20260531_064022.json").read_text())
    sim_record = json.loads((data_dir / "ep_onset_simulate_twirl_20260531_063048.json").read_text())

    def raw_population_errors(record, stored_key, q_factor=2.0):
        rows = record["scan"]
        q_label = [float(row["Q"]) for row in rows]
        revival = [max(float(pops[0]) for t, pops in row["pops"].items() if float(t) >= 2.0)
                   for row in rows]
        errors = []
        if [q_factor * q for q in q_label] != [1.0, 2.0, 3.0, 5.0, 10.0, 40.0]:
            errors.append("canonical Q grid is not 2*Q_label")
        if any(abs(value - float(row[stored_key])) > 1e-12 for value, row in zip(revival, rows)):
            errors.append("stored revival disagrees with max n0(t>=2 us)")
        return errors, revival

    raw_errors, hw_revival = raw_population_errors(hw_record, "revival")
    check("F86 raw hardware populations define revival and factor-two Q map", not raw_errors,
          ascii("; ".join(raw_errors)))
    mutated_factor_errors, _ = raw_population_errors(hw_record, "revival", q_factor=1.0)
    check("F86 raw rate-book mutation rejects missing factor two", bool(mutated_factor_errors))
    mutated_hw = json.loads(json.dumps(hw_record))
    mutated_hw["scan"][0]["pops"]["3.0"][0] += 0.2
    mutated_pop_errors, _ = raw_population_errors(mutated_hw, "revival")
    check("F86 raw population mutation rejects stale stored revival", bool(mutated_pop_errors))

    sim_errors, sim_revival = raw_population_errors(sim_record, "revival_max_n0")
    check("F86 raw exact-twirl endpoints are 0.28 to 0.84",
          not sim_errors and round(sim_revival[0], 2) == 0.28 and round(sim_revival[-1], 2) == 0.84)
    exact_twirl_endpoint = repr(sim_revival[-1])
    exact_twirl_consumers = {
        "raw-data README": root / "data/ibm_ep_onset_may2026/README.md",
        "flow synthesis": root / "experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md",
        "C# confirmation registry": root / "compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs",
        "Python confirmation registry": root / "simulations/framework/confirmations.py",
    }
    for label, path in exact_twirl_consumers.items():
        source = path.read_text(encoding="utf-8")
        check(f"F86 {label} carries exact raw twirl endpoint", exact_twirl_endpoint in source)
        mutated = source.replace(exact_twirl_endpoint, "0.837", 1)
        check(f"F86 {label} exact-endpoint mutation is rejected",
              mutated != source and exact_twirl_endpoint not in mutated)
    mutated_sim = json.loads(json.dumps(sim_record))
    mutated_sim["scan"][-1]["revival_max_n0"] += 0.01
    mutated_sim_errors, _ = raw_population_errors(mutated_sim, "revival_max_n0")
    check("F86 raw exact-twirl stored-endpoint mutation is rejected", bool(mutated_sim_errors))

    ep_source = (root / "compute/RCPsiSquared.Diagnostics/Foundation/EpField.cs").read_text(encoding="utf-8")
    match = re.search(r"var hwRev = new double\[\] \{([^}]+)\}", ep_source, re.S)
    ep_values = [float(x) for x in re.findall(r"\d+\.\d+", match.group(1))] if match else []
    check("EpField hardware curve equals raw recomputed revival",
          len(ep_values) == len(hw_revival) and np.allclose(ep_values, hw_revival, atol=1e-15))

    exact_literal = ", ".join(str(value) for value in hw_revival)
    registry_sources = {
        "C# confirmation registry": (root / "compute/RCPsiSquared.Core/Confirmations/ConfirmationsRegistry.cs").read_text(encoding="utf-8"),
        "Python confirmation registry": (root / "simulations/framework/confirmations.py").read_text(encoding="utf-8"),
    }
    for label, source in registry_sources.items():
        check(f"{label} revival vector equals raw recomputed values", exact_literal in source)
        mutated = source.replace(str(hw_revival[0]), str(hw_revival[0] + 0.01), 1)
        check(f"{label} raw-vector mutation is rejected",
              mutated != source and exact_literal not in mutated)

    part_a = json.loads((data_dir / "ep_onset_hardware_ibm_kingston_20260531_060943.json").read_text())
    endpoint = [float(x) for x in part_a["pops"][-1]]
    check("F86 Part-A endpoint is finite-time and not a normalized one-excitation distribution",
          part_a["ts_us"][-1] == 20.0 and not np.isclose(sum(endpoint), 1.0))

    csr_surfaces = {
        "compute/RCPsiSquared.Diagnostics/Foundation/IntegrabilityBreakingCsr.cs": (
            ("SpectrumClusterBootstrap95", "NoneDeterministicGrid",
             "InsufficientIndependentSpectra", "ReduceIndependentSpectra",
             "rng.Next(nonempty.Length)", "ReduceDeterministicGrid",
             "z-values within one spectrum are not treated as independent",
             "return w == 0.0", "bootstraps = 400"),
            ("abs[r.Next(nz)]", "Reduce(IReadOnlyList<Complex> zs")),
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdCsr.cs": (
            ("cluster-bootstrap whole independent disorder realizations",
             "treat a fixed q-grid as descriptive with no sampling CI",
             "ReduceIndependentSpectra(spectra", "ReduceDeterministicGrid(spectra)",
             "return w == 0.0"),
            ("IntegrabilityBreakingCsr.Reduce(pool",)),
        "compute/RCPsiSquared.Diagnostics/Foundation/FillingThresholdWitness.cs": (
            ("whole-spectrum cluster bootstrap", "independent disorder realizations"),
            ("individual-z bootstrap",)),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/IntegrabilityBreakingCsrTests.cs": (
            ("ClusterBootstrap_ResamplesWholeSpectra_NotIndividualSpacingRatios",
             "ClusterBootstrap_DuplicateIdenticalSpectraDoNotManufacturePrecision",
             "DeterministicGridReading_DoesNotClaimSamplingConfidenceInterval"),
            ("with 95% bootstrap CI",)),
        "experiments/FILLING_THRESHOLD_CHAOS.md": (
            ("ReduceIndependentSpectra", "whole independent spectra",
             "SpectrumClusterBootstrap95", "seed=4001", "seed=4002",
             "400 bootstrap resamples", "[0.676,0.685]", "[0.716,0.719]",
             "13,440/15,680 pooled z-values"),
            ("class's `Reduce`", "bootstrap the CI",
             "[0.677,0.684]", "[0.715,0.721]")),
        "compute/RCPsiSquared.Diagnostics.Tests/Foundation/FillingThresholdCsrTests.cs": (
            ("SpectrumClusterBootstrap95", "Assert.Equal(60, dilute.IndependentSpectrumCount)",
             "Assert.Equal(4, dense.IndependentSpectrumCount)",
             "seed: 4001", "seed: 4002",
             "Assert.Equal(0.676, dilute.CiLo, 3)", "Assert.Equal(0.685, dilute.CiHi, 3)",
             "Assert.Equal(0.716, dense.CiLo, 3)", "Assert.Equal(0.719, dense.CiHi, 3)"),
            ("Assert.Equal(0.677, dilute.CiLo, 3)", "Assert.Equal(0.684, dilute.CiHi, 3)",
             "Assert.Equal(0.715, dense.CiLo, 3)", "Assert.Equal(0.721, dense.CiHi, 3)")),
    }
    for path, (required, forbidden) in csr_surfaces.items():
        source = (root / path).read_text(encoding="utf-8")
        baseline = set(task8_text_errors(source, required, forbidden))
        check(f"CSR uncertainty current semantics: {path}", not baseline,
              ascii("; ".join(sorted(baseline))))
        for phrase in required:
            changed = source.replace(phrase, "")
            check(f"CSR uncertainty mutation rejects missing {ascii(phrase)}: {path}",
                  changed != source and
                  bool(set(task8_text_errors(changed, required, forbidden)) - baseline))
        for phrase in forbidden:
            changed = source + "\n" + phrase
            check(f"CSR uncertainty mutation rejects {ascii(phrase)}: {path}",
                  bool(set(task8_text_errors(changed, required, forbidden)) - baseline))


def verify_route_b_a2_current_truth():
    """Guard completion prose and inventory shape; C# owns the local character verdicts."""
    verify_route_b_mechanism_current_truth()
    verify_random_field_current_truth()
    verify_poisson_inference_current_truth()
    verify_csr_producer_and_crossover_current_truth()
    verify_sff_windows_and_current_boundaries()
    verify_task8_current_truth()
    verify_round4_artifacts_and_consumers()
    verify_final_five_repairs()
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
        baseline = set(route_b_n6_current_truth_errors(source, label))
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
                  bool(set(route_b_n6_current_truth_errors(wrong_source, label)) - baseline))
        # The F163 pointer is what keeps the numerical routes from reading as the
        # strongest available statement; deleting it must be caught, not tolerated.
        f163_sentence = ("The semisimplicity itself does not rest on those readings" if label == "primary"
                         else "F163 (typed as RouteBN6A2UnfoldingClaim, live at inspect --root n6unfolding) independently proves all 266 semisimple" if label == "arc"
                         else "Semisimplicity itself is settled algebraically for all 266")
        check(f"Route B N6 {label} mutation rejects dropping the F163 semisimplicity pointer",
              bool(set(route_b_n6_current_truth_errors(source.replace(f163_sentence, ""), label)) - baseline))
        if label == "primary":
            wrong_counts = source.replace("266 Diabolic; 0 Defective", "200 Diabolic; 66 Defective")
            wrong_proof = source.replace(
                "stable numerical character from the isolated two-dimensional restriction; it is not an exact rank proof",
                "an exact rank proof")
            check("Route B N6 mutation rejects 200 Diabolic / 66 Defective",
                  bool(set(route_b_n6_current_truth_errors(wrong_counts, label)) - baseline))
            check("Route B N6 mutation rejects numerical character promoted to exact rank",
                  bool(set(route_b_n6_current_truth_errors(wrong_proof, label)) - baseline))
        required = ["Route B is complete at N=5 as a root-by-root inventory",
                    "PSC1", "S1", "24 imaginary-q", "2 real-q", "32 nonreal-q",
                    "26-member R-odd", "38 distinct", "20 nonreal-q R-even",
                    "F164", "not an all-N theorem"]
        # Sensitive to deleting any N=6 completion/provenance boundary or restoring
        # the obsolete N=6-open statement; numerical character must stay numerical.
        required += ["Route B is complete at N=6 as a parity-labelled direct-t inventory",
                     "ExactAlgebraic", "HermitianAxis", "EpCharacterStable", "ExactRankExecuted",
                     "ROUTE_B_A2_N6_RECONCILE", "F_53 (N=7)"]
        if label != "arc":
            required += ["29 A2(w) roots / 58 q-loci", "24 imaginary-q loci",
                         "executed full-sector Hermiticity", "alg=geo=2",
                         "exactRankCertificates", "ROUTE_B_A2_RECONCILE"]
        if label == "primary":
            required += ["Route B is complete at N=5 as a root-by-root inventory: 29 A2(w) roots / 58 q-loci, all semisimple with alg=geo=2.",
                         "12 negative-real w roots", "2 real-q loci", "32 nonreal-q loci",
                         "three isolating radii", "neither is a local Jordan character test"]
        surface_forbidden = forbidden
        surface_forbidden += (
            "32 nonreal-q loci remain numerical-only",
            "the other 32 keep only the reading",
        )
        if label == "path experiment":
            required.append("N=6 Route B is the parity-labelled direct-t A2 inventory (266 loci, 133 per parity, all semisimple by F163)")
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
        contract_baseline = set(task8_text_errors(source, required, surface_forbidden))
        check(f"Route B {label}: N5/F164 text contract", not contract_baseline,
              ascii("; ".join(sorted(contract_baseline))))
        for phrase in required:
            check(f"Route B {label}: {phrase}", phrase in source,
                  f"{relative_path}: missing current-truth statement")
            changed = re.sub(re.escape(phrase), "", source, flags=re.IGNORECASE)
            check(f"Route B {label} mutation rejects missing {phrase}",
                  changed != source and
                  bool(set(task8_text_errors(changed, required, surface_forbidden)) - contract_baseline),
                  f"{relative_path}: required phrase mutation was ineffective")
        for phrase in surface_forbidden:
            check(f"Route B {label}: excludes {phrase}",
                  phrase.casefold() not in source.casefold(),
                  f"{relative_path}: stale Route B claim")
            changed = source + " " + phrase
            check(f"Route B {label} mutation rejects forbidden {phrase}",
                  bool(set(task8_text_errors(changed, required, surface_forbidden)) - contract_baseline),
                  f"{relative_path}: forbidden phrase mutation was ineffective")

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

    try:
        n6_inventory = json.loads((root / "simulations/results/route_b_a2_n6.json").read_text(encoding="utf-8"))
        baseline = set(route_b_n6_real_q_exclusion_errors(n6_inventory))
        check("Route B N6 artifact: all 266 exact tBox.real intervals exclude zero",
              not baseline, ascii("; ".join(sorted(baseline))))
        mutated = json.loads(json.dumps(n6_inventory))
        mutated["loci"][0]["tBox"]["real"]["upper"]["numerator"] = "0"
        check("Route B N6 artifact mutation rejects a tBox.real touching zero",
              bool(set(route_b_n6_real_q_exclusion_errors(mutated)) - baseline))
        wrong_n = json.loads(json.dumps(n6_inventory))
        wrong_n["n"] = 5
        check("Route B N6 artifact mutation rejects n=5",
              bool(set(route_b_n6_real_q_exclusion_errors(wrong_n)) - baseline))
        float_inverted = json.loads(json.dumps(n6_inventory))
        real_interval = float_inverted["loci"][0]["tBox"]["real"]
        real_interval["lower"]["numerator"] = 1.25
        real_interval["upper"]["numerator"] = -1.25
        check("Route B N6 artifact mutation rejects float inverted interval",
              bool(set(route_b_n6_real_q_exclusion_errors(float_inverted)) - baseline))
        duplicate = json.loads(json.dumps(n6_inventory))
        duplicate["loci"][1]["id"] = duplicate["loci"][0]["id"]
        check("Route B N6 artifact mutation rejects duplicate locus id",
              bool(set(route_b_n6_real_q_exclusion_errors(duplicate)) - baseline))
    except (OSError, ValueError, KeyError, TypeError) as error:
        check("Route B N6 artifact: readable exact t boxes", False, str(error))

    current_surfaces = {
        "docs/proofs/PROOF_N5_REAL_Q_DIABOLIC.md": (
            ("RouteBN5RealQSemisimpleClaim", "RouteBN5RealQSemisimpleWitness",
             "inspect --root n5diabolic", "24 imaginary-q", "2 real-q", "32 nonreal-q",
             "26-member R-odd", "38 distinct", "20 nonreal-q R-even", "numerical-only",
             "every one of the 266 persisted exact rational", "tBox.real",
             "only real-q points in the N=5/N=6 A2 inventories",
             "q is the settable parameter", "lambda is the resulting spectral eigenvalue"),
            ("No `Claim` and no `IInspectable` owns", "float read of a field", "classified and uncertified")),
        "docs/ANALYTICAL_FORMULAS.md": (
            ("RouteBN5RealQSemisimpleClaim", "RouteBN5RealQSemisimpleWitness",
             "inspect --root n5diabolic", "24 imaginary-q", "2 real-q", "32 nonreal-q",
             "26-member R-odd", "38 distinct", "20 nonreal-q R-even", "numerical-only",
             "every one of the 266 persisted exact rational"),
            ("classified-not-certified", "float read of a field")),
        "docs/GLOSSARY.md": (
            ("outward-rounded interval/ball arithmetic", "artifact field describes that producer's route",
             "24 imaginary-q", "2 real-q", "32 nonreal-q", "26-member R-odd",
             "38 exact-certified", "20 nonreal-q R-even", "numerical-only"),
            ("empty means classified only",)),
    }
    for relative_path, (required, forbidden) in current_surfaces.items():
        source = (root / relative_path).read_text(encoding="utf-8")
        source = " ".join(source.split())
        baseline = set(task8_text_errors(source, required, forbidden))
        check(f"Route B living surface current truth: {relative_path}",
              not baseline, ascii("; ".join(sorted(baseline))))
        for phrase in required:
            changed = re.sub(re.escape(phrase), "", source, flags=re.IGNORECASE)
            check(f"Route B living surface mutation rejects missing {ascii(phrase)}: {relative_path}",
                  changed != source and bool(set(task8_text_errors(changed, required, forbidden)) - baseline))
        for phrase in forbidden:
            changed = source + "\n" + phrase
            check(f"Route B living surface mutation rejects {ascii(phrase)}: {relative_path}",
                  bool(set(task8_text_errors(changed, required, forbidden)) - baseline))

    readout_surfaces = {
        path: (
            ("ReadoutFisher", "static calculation helper", "ReadingPowerWitness",
             "live `IInspectable`", "inspect --root decoder", "lowest endpoint Q=1",
             "does not sample the N=4 coherence-horizon EP", "no EP-specific verdict",
             READING_POWER_MONOTONICITY_TEST),
            ("live readout lab", "NoEpPeak_InAnyBasis_FiMonotoneInQ",
             "exceptional point Q = 1", "Q = 1 (EP)"))
        for path in (
            "experiments/ROUTE_B_N4_VIRTUAL_READOUT.md",
            "experiments/ROUTE_B_N4_READOUT_SEARCH.md",
            "experiments/ROUTE_B_N4_HISTOGRAM_FILTER.md",
        )
    }
    readout_surfaces["hypotheses/HANDSHAKE_GEOMETRY.md"] = (
        ("Q=1 is the lowest endpoint", "Q*(4)=1.87874", "does not sample",
         "no EP-specific verdict", "sampled Q=1 endpoint", "1555.3×",
         "X/Y remain nonzero", "population basis remains much stronger"),
        ("Q = 1 (EP)", "exceptional point Q = 1", "at the EP every read",
         "from 1.65 at the EP", "1670×", "only the population basis still reads"))
    for relative_path, (required, forbidden) in readout_surfaces.items():
        source = (root / relative_path).read_text(encoding="utf-8")
        source = " ".join(source.split())
        baseline = set(task8_text_errors(source, required, forbidden))
        check(f"ReadingPower living surface current truth: {relative_path}",
              not baseline, ascii("; ".join(sorted(baseline))))
        for phrase in required:
            changed = re.sub(re.escape(phrase), "", source, flags=re.IGNORECASE)
            check(f"ReadingPower mutation rejects missing {ascii(phrase)}: {relative_path}",
                  changed != source and bool(set(task8_text_errors(changed, required, forbidden)) - baseline))
        for phrase in forbidden:
            changed = source + "\n" + phrase
            check(f"ReadingPower mutation rejects {ascii(phrase)}: {relative_path}",
                  bool(set(task8_text_errors(changed, required, forbidden)) - baseline))


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
