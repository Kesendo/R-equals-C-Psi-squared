"""Finite-offset diagnostics around the fragile bridge's spectral-abscissa axis departure.

Two chains of `n_per_chain` qubits, Heisenberg within each, joined by a Heisenberg
bridge; chain A dephases at +gamma, chain B at -gamma, so Sigma-gamma is exactly 0
and the palindrome is centred on zero. Pi forces inversion lambda <-> -lambda;
Hermiticity preservation separately supplies conjugate pairing.  Together they
give the generic off-axis quartet {lambda, lambda*, -lambda, -lambda*}.  Below
gamma_crit every eigenvalue instead lies on the imaginary axis.

This producer independently selects a max-real-part eigenvalue at each gamma; it
does not perform branch continuation.  It reports two finite-offset diagnostics
and one derived symmetry value above the numerically bisected spectral-abscissa threshold:

  1. Re lambda ~ sqrt(delta), delta = gamma/gamma_crit - 1, the coefficient constant
     across the sampled decades.
  2. The finite-offset gap to the across-axis partner -lambda* equals
     2*abs(Re lambda), a derived symmetry identity rather than independent evidence.
  3. The finite-offset Petermann reading grows approximately by a decade per
     sampled decade of delta.

The square-root-like onset and simple-mode Petermann sequence motivate an EP2
hypothesis; the derived partner gap adds no independent evidence. The producer
does not execute a strict threshold
coalescence or Jordan-rank certificate. The existence of off-axis pairs is established;
its EP character remains OPEN.

It also reads the max-real-part selected mode's K across a wide gamma range,
omitting degenerate selections for which a single-vector K is basis-dependent.
No finite-offset reading is a threshold Jordan-character certificate.

Anchors: hypotheses/FRAGILE_BRIDGE.md Section 3, experiments/PT_SYMMETRY_ANALYSIS.md,
docs/ANALYTICAL_FORMULAS.md F19.
"""
import os
import sys

import numpy as np
import scipy.sparse as sp
from scipy.linalg import eig, eigvals

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "fragile_bridge_ep_signature.txt")

_LINES = []


def log(msg=""):
    print(msg, flush=True)
    _LINES.append(msg)


I2 = sp.eye(2, format="csr", dtype=complex)
SX = sp.csr_matrix(np.array([[0, 1], [1, 0]], dtype=complex))
SY = sp.csr_matrix(np.array([[0, -1j], [1j, 0]], dtype=complex))
SZ = sp.csr_matrix(np.array([[1, 0], [0, -1]], dtype=complex))


def op_at(op, qubit, n_qubits):
    result = sp.eye(1, format="csr", dtype=complex)
    for k in range(n_qubits):
        result = sp.kron(result, op if k == qubit else I2, format="csr")
    return result


def build_L(n_per_chain, gamma, J=1.0, J_bridge=0.1):
    """The gain-loss Liouvillian: +gamma on chain A, -gamma on chain B."""
    n_total = 2 * n_per_chain
    d = 2 ** n_total
    d2 = d * d
    H = sp.csr_matrix((d, d), dtype=complex)
    for i in range(n_per_chain - 1):
        for P in (SX, SY, SZ):
            H = H + J * (op_at(P, i, n_total) @ op_at(P, i + 1, n_total))
    for i in range(n_per_chain, 2 * n_per_chain - 1):
        for P in (SX, SY, SZ):
            H = H + J * (op_at(P, i, n_total) @ op_at(P, i + 1, n_total))
    for P in (SX, SY, SZ):
        H = H + J_bridge * (op_at(P, n_per_chain - 1, n_total)
                            @ op_at(P, n_per_chain, n_total))
    Id = sp.eye(d, format="csr", dtype=complex)
    L = -1j * (sp.kron(H, Id, format="csr") - sp.kron(Id, H.T, format="csr"))
    Id2 = sp.eye(d2, format="csr", dtype=complex)
    for k in range(n_per_chain):
        Zk = op_at(SZ, k, n_total)
        L = L + gamma * (sp.kron(Zk, Zk.conj(), format="csr") - Id2)
    for k in range(n_per_chain, 2 * n_per_chain):
        Zk = op_at(SZ, k, n_total)
        L = L + (-gamma) * (sp.kron(Zk, Zk.conj(), format="csr") - Id2)
    return L


def max_re(n_per_chain, gamma, J_bridge):
    return float(np.max(eigvals(build_L(n_per_chain, gamma, J_bridge=J_bridge).toarray()).real))


def gamma_crit(n_per_chain, J_bridge, lo=0.0, hi=2.0, iterations=60):
    """Bisect on max Re lambda leaving zero. The axis below is clean to ~1e-14."""
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        if max_re(n_per_chain, mid, J_bridge) > 1e-12:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def leading_mode(n_per_chain, gamma, J_bridge):
    """Return (lambda, gap to the nearest other eigenvalue, Petermann K)."""
    A = build_L(n_per_chain, gamma, J_bridge=J_bridge).toarray()
    w, vr = eig(A)
    wl, vl = eig(A.conj().T)
    i = int(np.argmax(w.real))
    gap = float(np.sort(np.abs(w - w[i]))[1])
    j = int(np.argmin(np.abs(np.conj(wl) - w[i])))
    r = vr[:, i] / np.linalg.norm(vr[:, i])
    l = vl[:, j] / np.linalg.norm(vl[:, j])
    overlap = abs(np.vdot(l, r)) ** 2
    return w[i], gap, (float("inf") if overlap == 0.0 else 1.0 / overlap)


NPC = 2
DELTAS = (1e-2, 1e-3, 1e-4, 1e-5)
SIMPLE_MODE_GAP_TOL = 1e-10

log("=" * 78)
log("The fragile bridge spectral-abscissa axis departure: EP character remains OPEN")
log("=" * 78)
log()
log("Two chains of 2 qubits, Heisenberg, gain-loss dephasing +/- gamma, Sigma-gamma = 0.")
log()

for J_BRIDGE in (1.0, 1.9):
    GC = gamma_crit(NPC, J_BRIDGE)
    log("J_bridge = {0}   gamma_crit = {1:.9f}".format(J_BRIDGE, GC))
    log()
    log("  {0:>8} {1:>15} {2:>15} {3:>15} {4:>13} {5:>13}".format(
    "delta", "max Re lambda", "Re/sqrt(delta)", "gap to -lambda*", "nearest gap", "Petermann K"))
    coefficients = []
    ks = []
    last_gap = None
    for DELTA in DELTAS:
        lam, gap, K = leading_mode(NPC, GC * (1.0 + DELTA), J_BRIDGE)
        if gap <= SIMPLE_MODE_GAP_TOL:
            raise RuntimeError(
                "near-threshold Petermann row is not a simple mode: "
                "J_bridge={0}, delta={1}, nearest gap={2}".format(J_BRIDGE, DELTA, gap))
        coefficients.append(lam.real / np.sqrt(DELTA))
        ks.append(K)
        last_gap = 2.0 * abs(lam.real)
        log("  {0:>8.0e} {1:>15.8f} {2:>15.4f} {3:>15.3e} {4:>13.3e} {5:>13.4g}".format(
            DELTA, lam.real, lam.real / np.sqrt(DELTA), last_gap, gap, K))
    log("  near-threshold simplicity gate: nearest gap > 1e-10 for every reported K")
    spread = max(coefficients) / min(coefficients)
    decades = np.sqrt(DELTAS[0] / DELTAS[-1])
    log()
    log("  1. Re lambda / sqrt(delta) spans {0:.4f} over four decades of delta."
        " A transversal crossing (Re ~ delta) would make this ratio grow by"
        " {1:.0f} across the same range.".format(spread, decades))
    log("  2. the gap to the across-axis partner -lambda* is 2*abs(Re lambda) = {0:.2e};"
        " it is derived from max Re and supplies no independent EP evidence."
        .format(last_gap))
    ratios = [ks[i + 1] / ks[i] for i in range(len(ks) - 1)]
    log("  3. Petermann K multiplies by {0} per sampled decade of delta; this"
        " finite-offset trend is not a Jordan-character certificate.".format(
            ", ".join("{0:.2f}".format(r) for r in ratios)))
    for below in (-1e-3, -1e-2, -1e-1):
        log("     below threshold, delta = {0:>+.0e}: max Re lambda = {1:.2e}".format(
            below, max_re(NPC, GC * (1.0 + below), J_BRIDGE)))
    log()

log("=" * 78)
log("What a coarse scan sees instead")
log("=" * 78)
log()
log("A finite grid reads sampling-dependent Petermann values for simple selected modes.")
log("A single-vector K is omitted wherever the selected eigenvalue is degenerate:")
log()
GC = gamma_crit(NPC, 1.0)
log("  {0:>17} {1:>43}".format("gamma/gamma_crit", "Petermann reading"))
for RATIO in (1.001, 1.01, 1.05, 1.2, 1.46, 1.5, 2.0, 3.0):
    _, gap, K = leading_mode(NPC, GC * RATIO, 1.0)
    if gap <= SIMPLE_MODE_GAP_TOL:
        reading = "single-vector K not reported for a degenerate eigenspace"
    else:
        reading = "K={0:.4g} (nearest gap {1:.2e})".format(K, gap)
    log("  {0:>17.3f} {1:>43}".format(RATIO, reading))
log()
log("Degenerate rows do not support a basis-invariant single-eigenvector conditioning value.")
log("Each gamma independently selects max Re(lambda);")
log("this producer establishes off-axis-pair existence, not branch continuation.")
log("It does not execute a strict threshold coalescence or Jordan-rank certificate;")
log("the local EP character remains OPEN.")

with open(OUT_PATH, "w", encoding="utf-8") as handle:
    handle.write("\n".join(_LINES) + "\n")
print("\nResults: " + OUT_PATH)
