"""C_block lens on IBM hardware framework_snapshots.

For each existing IBM Marrakesh / Kingston / Fez 2026-04-26 tomography snapshot
(plus the Aer simulator baseline at the same parameters), we reconstruct the
2-qubit reduced ρ on (q0, q2) from the 16-Pauli expectations and compute the
block-purity content C_block on the (popcount-0, popcount-1) and
(popcount-1, popcount-2) coherence blocks. We compare to the universal 1/4
ceiling per Theorem 2 of PROOF_BLOCK_CPSI_QUARTER and to the noiseless
continuous-time ideal evolution.

Question: how close to the 1/4 boundary do real hardware states sit, and does
the Theorem 2 lens differentiate the F87 trichotomy (truly / soft / hard) the
same way the F88 lens did at the operator level?

Snapshot configurations (all from |+−+⟩ at N=3, t=0.8, J=1.0):
  - heisenberg : H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1})
  - truly      : H = (J/4) Σ_b (X_b X_{b+1} + Y_b Y_{b+1})         (Π²-even)
  - soft       : H = (J/4) Σ_b (X_b Y_{b+1} + Y_b X_{b+1})         (Π²-odd)
  - hard       : H = (J/4) Σ_b (X_b X_{b+1} + X_b Y_{b+1})         (mixed)

Reduced ρ is on (q0, q2). On a 2-qubit Hilbert space the (popcount-n, popcount-(n+1))
blocks for n ∈ {0, 1} are the only non-trivial coherence-block options.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

sys.path.insert(0, str(Path(__file__).parent))
from _f88_lens_ibm_framework_snapshots import (  # noqa: E402
    I2,
    SX,
    SY,
    SZ,
    reconstruct_2qubit_rho,
)


SNAPSHOT_DIR = Path(
    r"D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI"
    r"\experiments\ibm_quantum_tomography\results"
)
SNAPSHOT_TIMESTAMP = "20260426_105948"
RESULTS_FILE = Path(__file__).parent / "results" / "_block_cpsi_lens_ibm_snapshots.txt"

QUARTER = 0.25  # universal C_block ceiling per Theorem 2

SCENARIOS = ["heisenberg", "truly", "soft", "hard"]


def c_block_2qubit(rho: np.ndarray, n: int) -> float:
    """Block-purity content on the (popcount-n, popcount-(n+1)) block of a 2-qubit ρ.

    For a 2-qubit ρ (4×4 in the {|00⟩, |01⟩, |10⟩, |11⟩} computational basis),
    sum |ρ_{ab}|² over (a, b) with popcount(a) = n, popcount(b) = n+1.
    Theorem 2 of PROOF_BLOCK_CPSI_QUARTER: this is ≤ 1/4 for any density matrix.
    """
    if n not in (0, 1):
        raise ValueError(f"For 2-qubit ρ only n in (0, 1) is admissible; got {n}")
    indices_a = [i for i in range(4) if bin(i).count("1") == n]
    indices_b = [j for j in range(4) if bin(j).count("1") == n + 1]
    s = 0.0
    for a in indices_a:
        for b in indices_b:
            s += abs(rho[a, b]) ** 2
    return float(s)


def ideal_evolution_q0q2(h_type: str, t: float = 0.8, J: float = 1.0) -> np.ndarray:
    """|+−+⟩ at N=3 evolves under bond Hamiltonian H_type; partial trace to (q0, q2)."""
    psi_plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi_minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(psi_plus, psi_minus), psi_plus)

    def site_op(P: np.ndarray, k: int, N: int = 3) -> np.ndarray:
        ops = [I2] * N
        ops[k] = P
        return np.kron(np.kron(ops[0], ops[1]), ops[2])

    bond_pairs = {
        "heisenberg": [(SX, SX), (SY, SY), (SZ, SZ)],
        "truly":      [(SX, SX), (SY, SY)],
        "soft":       [(SX, SY), (SY, SX)],
        "hard":       [(SX, SX), (SX, SY)],
    }
    if h_type not in bond_pairs:
        raise ValueError(f"unknown h_type: {h_type}")

    H = np.zeros((8, 8), dtype=complex)
    for b in [0, 1]:
        for Pa, Pb in bond_pairs[h_type]:
            H += (J / 4.0) * site_op(Pa, b) @ site_op(Pb, b + 1)

    from scipy.linalg import expm
    U = expm(-1j * H * t)
    psi_t = U @ psi
    rho_full = np.outer(psi_t, psi_t.conj())

    rho_q0q2 = np.zeros((4, 4), dtype=complex)
    for i0 in range(2):
        for i2 in range(2):
            for j0 in range(2):
                for j2 in range(2):
                    s = 0.0 + 0j
                    for q1 in range(2):
                        idx_i = (i0 << 2) | (q1 << 1) | i2
                        idx_j = (j0 << 2) | (q1 << 1) | j2
                        s += rho_full[idx_i, idx_j]
                    rho_q0q2[(i0 << 1) | i2, (j0 << 1) | j2] = s
    return rho_q0q2


def load_snapshot_rhos(snapshot_json: Path) -> dict[str, np.ndarray]:
    """Reconstruct 2-qubit ρ for each scenario present in the snapshot JSON."""
    with open(snapshot_json, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    rhos: dict[str, np.ndarray] = {}
    rhos["heisenberg"] = reconstruct_2qubit_rho(data["snapshot_c_heisenberg"]["expectations"])
    if "snapshot_d_softbreak_trichotomy" in data:
        for cat in ["truly", "soft", "hard"]:
            rhos[cat] = reconstruct_2qubit_rho(
                data["snapshot_d_softbreak_trichotomy"]["expectations_per_category"][cat]
            )
    return rhos


def find_snapshot_files() -> dict[str, Path]:
    paths: dict[str, Path] = {}
    paths["aer (Trotter noiseless)"] = SNAPSHOT_DIR / "framework_snapshots_aer_20260426_012516.json"
    for b in ["ibm_marrakesh", "ibm_kingston", "ibm_fez"]:
        paths[b] = SNAPSHOT_DIR / f"framework_snapshots_{b}_{SNAPSHOT_TIMESTAMP}.json"
    return paths


def main() -> None:
    out_lines: list[str] = []

    def emit(line: str = "") -> None:
        out_lines.append(line)
        print(line)

    emit("C_block lens on IBM framework_snapshots (2026-04-26)")
    emit("=" * 80)
    emit(f"Universal Theorem 2 ceiling: C_block ≤ 1/4 = {QUARTER:.4f} (any density matrix on 2^N).")
    emit("Block (n, n+1): coherence content Σ |ρ_{ab}|² with popcount(a)=n, popcount(b)=n+1.")
    emit("Reduced ρ on (q0, q2) of the N=3 |+−+⟩ chain at t=0.8, J=1.0.")
    emit()

    # Continuous-time noiseless ideal per scenario, for reference and comparison.
    ideal_c_block: dict[str, tuple[float, float]] = {}
    for scen in SCENARIOS:
        rho_id = ideal_evolution_q0q2(scen, t=0.8, J=1.0)
        ideal_c_block[scen] = (
            c_block_2qubit(rho_id, n=0),
            c_block_2qubit(rho_id, n=1),
        )

    emit("Ideal (continuous, noiseless) per scenario:")
    emit("  scenario     block(0,1) = c   % of 1/4   block(1,2) = c   % of 1/4")
    for scen in SCENARIOS:
        c01, c12 = ideal_c_block[scen]
        emit(
            f"  {scen:<11s}  {c01:.4f}           {c01/QUARTER*100:5.1f}%      "
            f"{c12:.4f}           {c12/QUARTER*100:5.1f}%"
        )
    emit()

    # Per-backend table (Aer first, then hardware).
    paths = find_snapshot_files()
    backend_rhos: dict[str, dict[str, np.ndarray]] = {}
    for backend, path in paths.items():
        if not path.exists():
            emit(f"=== {backend} ===  [SKIP: file not found] {path}")
            emit()
            continue
        rhos = load_snapshot_rhos(path)
        backend_rhos[backend] = rhos
        emit(f"=== {backend} ===")
        emit(f"  source: {path.name}")
        emit("  scenario     block   C_block    % of 1/4   ideal     Δ(this − ideal)")
        for scen in SCENARIOS:
            if scen not in rhos:
                continue
            rho = rhos[scen]
            for n in (0, 1):
                cb = c_block_2qubit(rho, n)
                cbi = ideal_c_block[scen][n]
                emit(
                    f"  {scen:<11s} ({n},{n+1})   {cb:.4f}     {cb/QUARTER*100:5.1f}%     "
                    f"{cbi:.4f}    {cb - cbi:+.4f}"
                )
        emit()

    # Cross-backend spread, hardware only.
    hw_rhos = {b: r for b, r in backend_rhos.items() if b.startswith("ibm_")}
    if hw_rhos:
        emit("Cross-backend spread (ibm_* hardware only):")
        for scen in SCENARIOS:
            present = [b for b, r in hw_rhos.items() if scen in r]
            if not present:
                continue
            emit(f"  [{scen}]")
            for n in (0, 1):
                vals = [(b, c_block_2qubit(hw_rhos[b][scen], n)) for b in present]
                cs = [v for _, v in vals]
                spread_str = "  ".join(
                    f"{b.replace('ibm_',''):>10s}={v:.4f}" for b, v in vals
                )
                emit(
                    f"    block({n},{n+1}):  mean={np.mean(cs):.4f}  "
                    f"std={np.std(cs):.4f}  range={max(cs) - min(cs):.4f}  | {spread_str}"
                )
            emit()

    # Aer (noiseless Trotter) vs hardware: pure-noise contribution per scenario.
    aer_label = "aer (Trotter noiseless)"
    if aer_label in backend_rhos and hw_rhos:
        emit("Hardware noise contribution (mean over hw backends − Aer baseline):")
        emit("  scenario     block   Aer C_block   ⟨hw⟩ C_block   noise Δ⟨hw⟩−Aer")
        aer_rhos = backend_rhos[aer_label]
        for scen in SCENARIOS:
            if scen not in aer_rhos:
                continue
            present = [b for b, r in hw_rhos.items() if scen in r]
            if not present:
                continue
            for n in (0, 1):
                aer_c = c_block_2qubit(aer_rhos[scen], n)
                hw_mean = np.mean([c_block_2qubit(hw_rhos[b][scen], n) for b in present])
                emit(
                    f"  {scen:<11s} ({n},{n+1})   {aer_c:.4f}        "
                    f"{hw_mean:.4f}         {hw_mean - aer_c:+.4f}"
                )
        emit()

    # Sanity: every state must satisfy C_block ≤ 1/4 (Theorem 2). Assert.
    for backend, rhos in backend_rhos.items():
        for scen, rho in rhos.items():
            for n in (0, 1):
                cb = c_block_2qubit(rho, n)
                if cb > QUARTER + 1e-6:
                    emit(
                        f"  WARNING: Theorem 2 violation candidate: "
                        f"{backend} {scen} block({n},{n+1}) = {cb:.6f} > 0.25"
                    )

    emit("=" * 80)
    emit("Reading legend (no claims yet, look at the numbers):")
    emit("  - C_block ≤ 1/4: hard universal ceiling per Theorem 2; no state can exceed.")
    emit("  - % of 1/4: how close to the Mandelbrot-cardioid boundary the state sits at this lens.")
    emit("  - Δ(this − ideal): deviation from the noiseless continuous-time evolution; mixes")
    emit("    Trotter discretisation error and hardware noise (Aer row isolates Trotter).")
    emit("  - Cross-backend std: how uniform the three ibm_* hardware platforms are at this lens.")

    RESULTS_FILE.parent.mkdir(parents=True, exist_ok=True)
    RESULTS_FILE.write_text("\n".join(out_lines), encoding="utf-8")
    print(f"\nResults written to {RESULTS_FILE}")


if __name__ == "__main__":
    main()
