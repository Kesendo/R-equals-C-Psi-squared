"""Producer for the rate-embedding claim behind MirrorAnalysis.CheckSymmetry.

`compute/RCPsiSquared.Compute/MirrorAnalysis.cs` scores the F1 palindrome of a
DECAY RATE multiset by handing it to the canonical complex check
`F1SpectrumStatistics.MaxF1PairingDistance`, which reflects lambda -> -2*sigma - lambda.
It does that by embedding lambda = -d. The claim this script exists to make
reproducible is that the embedding is not an approximation:

    with lambda = -d,  -2*sigma - lambda = d - 2*sigma = -(2*sigma - d)

which is the embedding of the mirrored rate. Since d -> -d is an isometry and
preserves the index order the greedy matcher walks, the embedded route and a
direct with-removal matcher on {d} vs {2*sigma - d} return the SAME number, not
merely a close one.

It also reports what the rate projection costs, by printing the full complex F1
distance beside the rate-only one. That ratio is a single-system reading, not a
law: it moves with N, topology and solver.

Convention: this uses the Python framework's Pauli-J book, which is the same
book `compute/RCPsiSquared.Compute/Topology.cs` builds in (H = sum J*(XX+YY+ZZ)
with plain Pauli operators, no 1/4). Absolute distances still differ from the
C# run because numpy's eigensolver is not MKL zgeev; the EQUALITY of the two
routes is what is exact here, and that is solver-independent.

Run:  python simulations/f1_rate_embedding_check.py
Out:  simulations/results/f1_rate_embedding_check.txt
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import framework as fw

GAMMA = 0.05
IM_THRESHOLD = 0.05          # the engine's oscillatory cut, Liouvillian.ExtractRates
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "results", "f1_rate_embedding_check.txt")


def oscillatory_rates(evals, threshold=IM_THRESHOLD):
    """The engine's projection: d = -Re(lambda) for the modes with |Im| > threshold.

    No rounding and no low-rate floor. Both were removed from ExtractRates on
    2026-08-06; the floor in particular is not F1-covariant, since it cuts the
    low-d end of a reflection that maps low d to high d.
    """
    return sorted(-ev.real for ev in evals if abs(ev.imag) > threshold)


def direct_rate_distance(rates, sigma):
    """Max nearest-neighbour distance WITH REMOVAL between {d} and {2*sigma - d}."""
    rates = np.asarray(rates, dtype=float)
    mirror = 2.0 * sigma - rates
    taken = np.zeros(len(mirror), dtype=bool)
    worst = 0.0
    for d in rates:
        dist = np.abs(mirror - d)
        dist[taken] = np.inf
        j = int(np.argmin(dist))
        taken[j] = True
        if dist[j] > worst:
            worst = float(dist[j])
    return worst


def embedded_rate_distance(rates, sigma):
    """The route the C# code takes: embed lambda = -d, call the canonical check."""
    return fw.max_f1_pairing_distance(np.asarray([-d for d in rates], dtype=complex), sigma)


def main():
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("F1 rate-embedding check")
    emit("=" * 78)
    emit(f"gamma = {GAMMA}, J = 1.0, |Im| > {IM_THRESHOLD}, Pauli-J convention")
    emit()
    emit(f"{'topology':<10} {'N':>2} {'rates':>6} {'direct':>13} {'embedded':>13} "
         f"{'equal':>6} {'full complex':>13} {'ratio':>7}")
    emit("-" * 78)

    all_equal = True
    for topology in ("chain", "star"):
        for N in (2, 3, 4, 5):
            system = fw.ChainSystem(N=N, gamma_0=GAMMA, J=1.0, topology=topology)
            evals = np.linalg.eigvals(system.L)
            rates = oscillatory_rates(evals)
            sigma = N * GAMMA

            direct = direct_rate_distance(rates, sigma)
            embedded = embedded_rate_distance(rates, sigma)
            equal = (direct == embedded)
            all_equal = all_equal and equal

            full = fw.max_f1_pairing_distance(evals, sigma)
            ratio = full / direct if direct > 0 else float("nan")

            emit(f"{topology:<10} {N:>2} {len(rates):>6} {direct:>13.6e} {embedded:>13.6e} "
                 f"{str(equal):>6} {full:>13.6e} {ratio:>7.2f}")

    emit()
    emit(f"embedded route equals direct route on every row: {all_equal}")
    emit()
    emit("Reading. The 'equal' column is an EXACT equality, not a tolerance, and that")
    emit("is the point: the embedding is an isometry, so there is nothing to round.")
    emit("The 'ratio' column is what the projection to rates costs on that one system.")
    emit("It varies across these rows, so do not quote a single value for it.")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    if not all_equal:
        raise SystemExit("EMBEDDING CHECK FAILED: the two routes disagree")


if __name__ == "__main__":
    main()
