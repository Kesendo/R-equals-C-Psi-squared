"""The fragile bridge's threshold is a second-order exceptional point on the REAL gamma axis.

Two chains of `n_per_chain` qubits, Heisenberg within each, joined by a Heisenberg
bridge; chain A dephases at +gamma, chain B at -gamma, so Sigma-gamma is exactly 0
and the palindrome is centred on zero. Pi then forces lambda <-> -lambda, which
leaves an eigenvalue two options: sit on the imaginary axis, or carry a partner
mirrored across it. Below gamma_crit every eigenvalue takes the first.

This producer measures the three signatures that separate a coalescence from a
crossing, all at real gamma = gamma_crit:

  1. Re lambda ~ sqrt(delta), delta = gamma/gamma_crit - 1, the coefficient constant
     across decades. A Hopf bifurcation is a transversal crossing, Re lambda ~ delta
     with a finite nonzero slope.
  2. The gap to the mirror partner closes to zero: the pair merges.
  3. The Petermann factor of the leading mode diverges as 1/delta. At an ordinary
     non-normal point it is large and finite, so a quoted "peak height" is a
     statement about the scan's step size rather than about the system.

It also reads the leading mode's K across a wide gamma range, which is what a coarse
scan sees: a secondary bump well above threshold, unrelated to the divergence at it.

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

log("=" * 78)
log("The fragile bridge threshold: coalescence, not crossing")
log("=" * 78)
log()
log("Two chains of 2 qubits, Heisenberg, gain-loss dephasing +/- gamma, Sigma-gamma = 0.")
log()

for J_BRIDGE in (1.0, 1.9):
    GC = gamma_crit(NPC, J_BRIDGE)
    log("J_bridge = {0}   gamma_crit = {1:.9f}".format(J_BRIDGE, GC))
    log()
    log("  {0:>8} {1:>15} {2:>15} {3:>15} {4:>13}".format(
        "delta", "max Re lambda", "Re/sqrt(delta)", "gap to partner", "Petermann K"))
    coefficients = []
    ks = []
    last_gap = None
    for DELTA in DELTAS:
        lam, gap, K = leading_mode(NPC, GC * (1.0 + DELTA), J_BRIDGE)
        coefficients.append(lam.real / np.sqrt(DELTA))
        ks.append(K)
        last_gap = gap
        log("  {0:>8.0e} {1:>15.8f} {2:>15.4f} {3:>15.3e} {4:>13.4g}".format(
            DELTA, lam.real, lam.real / np.sqrt(DELTA), gap, K))
    spread = max(coefficients) / min(coefficients)
    decades = np.sqrt(DELTAS[0] / DELTAS[-1])
    log()
    log("  1. Re lambda / sqrt(delta) spans {0:.4f} over four decades of delta."
        " A transversal crossing (Re ~ delta) would make this ratio grow by"
        " {1:.0f} across the same range.".format(spread, decades))
    log("  2. the gap to the mirror partner falls to {0:.2e}: the pair merges."
        .format(last_gap))
    ratios = [ks[i + 1] / ks[i] for i in range(len(ks) - 1)]
    log("  3. Petermann K multiplies by {0} per decade of delta: the 1/delta"
        " divergence of an exceptional point.".format(
            ", ".join("{0:.2f}".format(r) for r in ratios)))
    for below in (-1e-3, -1e-2, -1e-1):
        log("     below threshold, delta = {0:>+.0e}: max Re lambda = {1:.2e}".format(
            below, max_re(NPC, GC * (1.0 + below), J_BRIDGE)))
    log()

log("=" * 78)
log("What a coarse scan sees instead")
log("=" * 78)
log()
log("K diverges only AT threshold, so a scan that never samples close to it reads a")
log("modest number somewhere else. The leading mode's K across a wide range:")
log()
GC = gamma_crit(NPC, 1.0)
log("  {0:>17} {1:>13}".format("gamma/gamma_crit", "Petermann K"))
for RATIO in (1.001, 1.01, 1.05, 1.2, 1.46, 1.5, 2.0, 3.0):
    _, _, K = leading_mode(NPC, GC * RATIO, 1.0)
    log("  {0:>17.3f} {1:>13.4g}".format(RATIO, K))
log()
log("There is a secondary bump above threshold, near gamma/gamma_crit = 1.46, an order")
log("of magnitude below the reading at delta = 1e-3. A height and a location quoted")
log("from different points are not one measurement.")

with open(OUT_PATH, "w", encoding="utf-8") as handle:
    handle.write("\n".join(_LINES) + "\n")
print("\nResults: " + OUT_PATH)
