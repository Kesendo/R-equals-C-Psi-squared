"""F99 anchor lens on IBM Marrakesh + Kingston soft-break hardware data.

Companion to `_f88b_lens_ibm_framework_snapshots.py` (F88b Π²-odd-of-memory lens).
This script applies the F99 canonical-Niven anchor lens — testing whether
state-level measurements on hardware cluster around the five Niven anchors
{0, 1/8, 1/4, 3/8, 1/2} (plus polarity complements {5/8, 3/4, 7/8, 1}) that
F99 / CanonicalTrigAnchorPi2Inheritance derives via α = sin²(θ)/2 at the
five canonical trig angles {0°, 30°, 45°, 60°, 90°}.

Data sources (in-repo, not external):
  - data/ibm_soft_break_april2026/soft_break_ibm_marrakesh_20260426_001101.json
    (3 categories: truly_unbroken, soft_broken, hard_broken; q-path [48,49,50])
  - data/ibm_soft_break_april2026/soft_break_ibm_kingston_20260505_102806.json
    (4 categories: truly_unbroken, pi2_odd_pure, pi2_even_nontruly,
     mixed_anti_one_sixth; q-path [43,56,63] uniform-quantum)

Both runs measure 16-Pauli tomography on (q0, q2) of the N=3 |+−+⟩ state
at t=0.8, J=1.0, 3 Trotter steps, 4096 shots.

Niven anchors are equispaced on [0,1] with spacing 1/8. Uniform-random null:
expected nearest-anchor distance E[|X - nearest|] = 1/16 = 0.0625. A
measurement closer than 1/16 to an anchor is "tighter than random"; the
strength of clustering is quantified per measurement.

Reading: do hardware Π²-odd / memory_frac / static_frac values pile up on
specific Niven anchors, or are they anchor-blind? With only 7 (backend, category)
combinations this is a structural pilot, not a statistical claim.

Run:
  PYTHONIOENCODING=utf-8 python simulations/_f99_anchor_lens_ibm_snapshots.py
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


REPO_ROOT = Path(__file__).resolve().parents[1]
MARRAKESH_JSON = REPO_ROOT / "data" / "ibm_soft_break_april2026" / "soft_break_ibm_marrakesh_20260426_001101.json"
KINGSTON_JSON = REPO_ROOT / "data" / "ibm_soft_break_april2026" / "soft_break_ibm_kingston_20260505_102806.json"


# F99 Niven anchors on [0,1]: {0, 1/8, 1/4, 3/8, 1/2} + polarity complements.
NIVEN_ANCHORS = [Fraction(n, 8) for n in range(9)]
NIVEN_NAMES = {
    Fraction(0, 8): "0     (vacuum / θ=0°)",
    Fraction(1, 8): "1/8   (F99 θ=30°)",
    Fraction(2, 8): "1/4   (Quarter / θ=45°)",
    Fraction(3, 8): "3/8   (F86b KIntermediate / θ=60°)",
    Fraction(4, 8): "1/2   (Half fixed point / θ=90°)",
    Fraction(5, 8): "5/8   (F99 θ=60° β-face)",
    Fraction(6, 8): "3/4   (Quarter complement)",
    Fraction(7, 8): "7/8   (F99 θ=30° β-face)",
    Fraction(8, 8): "1     (full / noble endpoint)",
}
SPACING = 1.0 / 8.0
RANDOM_NULL_EXPECTED = SPACING / 2.0  # = 1/16 = 0.0625


# ---------------------- 2-qubit ρ reconstruction + lens (lifted from F88b-lens script) ----------------------

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {"I": I2, "X": SX, "Y": SY, "Z": SZ}
BIT_B = {"I": 0, "X": 0, "Y": 1, "Z": 1}


def reconstruct_2qubit_rho(expectations: dict[str, float]) -> np.ndarray:
    rho = np.zeros((4, 4), dtype=complex)
    for key, expv in expectations.items():
        a, b = key.split(",")
        sigma = np.kron(PAULIS[a], PAULIS[b])
        rho += expv * sigma / 4.0
    return rho


def state_level_lens_2qubit(rho: np.ndarray) -> dict[str, float]:
    """Static / memory / Π²-odd-of-memory split at N=2."""
    P0 = np.diag([1, 0, 0, 0]).astype(complex)
    P1 = np.diag([0, 1, 1, 0]).astype(complex)
    P2 = np.diag([0, 0, 0, 1]).astype(complex)

    rho_d0 = (
        np.real(np.trace(P0 @ rho)) / 1.0 * P0
        + np.real(np.trace(P1 @ rho)) / 2.0 * P1
        + np.real(np.trace(P2 @ rho)) / 1.0 * P2
    )
    rho_d2 = rho - rho_d0

    odd_memory_sq = 0.0
    for a in "IXYZ":
        for b in "IXYZ":
            if (BIT_B[a] + BIT_B[b]) & 1 == 0:
                continue
            sigma = np.kron(PAULIS[a], PAULIS[b])
            coeff = np.trace(sigma @ rho_d2) / 4.0
            odd_memory_sq += abs(coeff) ** 2 * 4

    norm_total = np.linalg.norm(rho, "fro") ** 2
    norm_static = np.linalg.norm(rho_d0, "fro") ** 2
    norm_memory = np.linalg.norm(rho_d2, "fro") ** 2

    return {
        "static_frac": float(norm_static / norm_total) if norm_total > 0 else 0.0,
        "memory_frac": float(norm_memory / norm_total) if norm_total > 0 else 0.0,
        "pi2_odd_in_memory": float(odd_memory_sq / norm_memory) if norm_memory > 1e-14 else 0.0,
    }


# ---------------------- F99 anchor lens ----------------------

def nearest_niven_anchor(x: float) -> tuple[Fraction, float, float]:
    """Return (anchor, residual = x - anchor, relative_tightness = |residual| / SPACING).

    relative_tightness < 0.5 means the value sits inside its anchor's bin (the
    half-spacing neighbourhood). < 0.1 means within 10% of the bin spacing
    (= within 1.25% on the absolute [0,1] axis); strong anchor alignment.
    """
    distances = [(a, x - float(a)) for a in NIVEN_ANCHORS]
    anchor, residual = min(distances, key=lambda t: abs(t[1]))
    return anchor, residual, abs(residual) / SPACING


def render_anchor_row(label: str, value: float) -> str:
    anchor, residual, tightness = nearest_niven_anchor(value)
    tighter_than_random = abs(residual) < RANDOM_NULL_EXPECTED
    flag = "★" if tightness < 0.1 else ("·" if tighter_than_random else " ")
    return (
        f"  {label:<24} {value:+.4f}  →  nearest {str(anchor):>5}  "
        f"(Δ = {residual:+.4f}, tightness = {tightness:.2f}) {flag}  {NIVEN_NAMES[anchor]}"
    )


def render_category(backend_label: str, category: str, expectations: dict[str, float]) -> dict[str, float]:
    rho = reconstruct_2qubit_rho(expectations)
    lens = state_level_lens_2qubit(rho)
    print(f"--- {backend_label}  /  {category} ---")
    print(render_anchor_row("pi2_odd_in_memory", lens["pi2_odd_in_memory"]))
    print(render_anchor_row("memory_frac",        lens["memory_frac"]))
    print(render_anchor_row("static_frac",        lens["static_frac"]))
    print()
    return lens


def main():
    print("=" * 96)
    print("F99 Anchor Lens on IBM soft-break hardware (Marrakesh 2026-04-26 + Kingston 2026-05-05)")
    print("=" * 96)
    print()
    print(f"Niven anchors: {{0, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, 1}} (spacing 1/8)")
    print(f"Random-null E[|x - nearest anchor|] = 1/16 = {RANDOM_NULL_EXPECTED:.4f}")
    print(f"Legend: ★ tightness < 0.10 (within 10% of bin, very tight)")
    print(f"        · tightness < 0.50 (closer than uniform-random)")
    print(f"        (blank) tightness ≥ 0.50 (anchor-blind)")
    print()

    all_measurements: list[tuple[str, str, str, float]] = []

    print("=" * 96)
    print(f"MARRAKESH (Heron-r2, ibm_marrakesh; q-path={json.load(open(MARRAKESH_JSON))['path']}, t=0.8, 4096 shots)")
    print("=" * 96)
    with open(MARRAKESH_JSON, "r", encoding="utf-8") as fh:
        marr = json.load(fh)
    for cat, expects in marr["expectations"].items():
        lens = render_category("Marrakesh", cat, expects)
        for k, v in lens.items():
            all_measurements.append(("Marrakesh", cat, k, v))

    print("=" * 96)
    print(f"KINGSTON  (Heron-r2, ibm_kingston;  q-path={json.load(open(KINGSTON_JSON))['path']}, t=0.8, 4096 shots, uniform-quantum)")
    print("=" * 96)
    with open(KINGSTON_JSON, "r", encoding="utf-8") as fh:
        king = json.load(fh)
    for cat, expects in king["expectations"].items():
        lens = render_category("Kingston", cat, expects)
        for k, v in lens.items():
            all_measurements.append(("Kingston", cat, k, v))

    # ---------------------- Summary across all measurements ----------------------

    print("=" * 96)
    print("SUMMARY: do hardware measurements cluster on Niven anchors?")
    print("=" * 96)
    print()

    n_total = len(all_measurements)
    tightnesses = [nearest_niven_anchor(v)[2] for _, _, _, v in all_measurements]
    n_very_tight = sum(1 for t in tightnesses if t < 0.1)
    n_tight = sum(1 for t in tightnesses if t < 0.5)
    n_blind = n_total - n_tight
    mean_residual_abs = float(np.mean([abs(nearest_niven_anchor(v)[1]) for _, _, _, v in all_measurements]))

    print(f"  Total measurements:                {n_total}")
    print(f"  ★ Very tight (tightness < 0.10):   {n_very_tight}  ({100*n_very_tight/n_total:.0f}%)")
    print(f"  · Tight     (tightness < 0.50):   {n_tight}  ({100*n_tight/n_total:.0f}%)  [≥ random-null]")
    print(f"      blind     (tightness ≥ 0.50):   {n_blind}  ({100*n_blind/n_total:.0f}%)")
    print()
    print(f"  Mean |Δ| across all measurements:  {mean_residual_abs:.4f}")
    print(f"  Uniform-random null E[|Δ|]:        {RANDOM_NULL_EXPECTED:.4f}")
    ratio = mean_residual_abs / RANDOM_NULL_EXPECTED
    print(f"  Ratio (measured / null):           {ratio:.2f}  ({'tighter than random' if ratio < 1 else 'looser than random'})")
    print()

    # Per-anchor cluster count
    cluster_counts: dict[Fraction, int] = {a: 0 for a in NIVEN_ANCHORS}
    very_tight_per_anchor: dict[Fraction, list[str]] = {a: [] for a in NIVEN_ANCHORS}
    for backend, cat, name, v in all_measurements:
        anchor, residual, tightness = nearest_niven_anchor(v)
        cluster_counts[anchor] += 1
        if tightness < 0.1:
            very_tight_per_anchor[anchor].append(f"{backend}/{cat}/{name} ({v:+.4f}, Δ={residual:+.4f})")

    print("  Per-anchor measurement-count distribution:")
    for a in NIVEN_ANCHORS:
        n_at_anchor = cluster_counts[a]
        bar = "█" * n_at_anchor
        very_tight_hits = very_tight_per_anchor[a]
        marker = "  ★ " + "; ".join(very_tight_hits) if very_tight_hits else ""
        print(f"    {str(a):>5}  {bar:<8}  n = {n_at_anchor}{marker}")
    print()

    # ---------------------- Reading ----------------------

    print("=" * 96)
    print("READING")
    print("=" * 96)
    print()
    print("This is a structural pilot, not a statistical claim — 7 (backend,category)")
    print("combinations × 3 lens measurements = 21 data points. Honest framing:")
    print()
    print(" 1. The strongest 'tight on anchor' hits (★ marker, tightness < 0.10) name the")
    print("    specific anchors that hardware actually lands on. If the very-tight set is")
    print("    dominated by 0 (truly categories → near-zero Π²-odd) and 3/4 (soft → near")
    print("    3/4 polarity complement), the lens replicates the F88b-reading at anchor level.")
    print()
    print(" 2. The mean |Δ| vs uniform-random null is the bulk test. Ratio < 1 means the")
    print("    measurements as a whole sit closer to Niven anchors than uniform-random")
    print("    values would. Ratio ≈ 1 is anchor-blind. Ratio > 1 would be ANTI-clustering.")
    print()
    print(" 3. The per-anchor cluster count surfaces WHICH anchors carry the measurements.")
    print("    Anchors that pile up multiple hits across backends are reproducible structural")
    print("    targets; lonely anchors are coincidence candidates.")
    print()
    print("Open follow-ups (not addressed here):")
    print(" - Compare hardware measurements to noiseless-simulation (Trotter + closed-form)")
    print("   to separate F99-anchor alignment from F87-trichotomy alignment. F87 alone might")
    print("   account for some/all of the clustering without invoking the Niven anchor set.")
    print(" - Extend to block-CΨ saturation data (Kingston May 2026) and f83 signature")
    print("   data (Marrakesh Apr 30) for a larger N independent of the soft-break setup.")
    print(" - The 1/8 anchor (F99 θ=30°, γ=√3/2, depth-3) has no in-protocol predicted hit;")
    print("   would a state designed to target α=1/8 (non-uniform Dicke) actually land there?")


if __name__ == "__main__":
    main()
