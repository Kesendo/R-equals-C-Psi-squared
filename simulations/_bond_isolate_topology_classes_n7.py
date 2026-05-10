"""Multi-bond topology-class comparison at N=7.

Tests the S_N-orbit prediction: bond configurations in the same topology class
(same multiset of connected-path-lengths) give bit-identical S(t).

Topology classification: each active-bond set yields a graph on the 7 sites.
The connected components are paths; we classify by the sorted multiset of path
lengths. E.g. bonds {0,1,3} = bond 0 (sites 0-1) + bond 1 (sites 1-2) + bond 3
(sites 3-4) → path-of-length-2 (sites 0-1-2) + path-of-length-1 (sites 3-4)
→ topology class (1, 2).

Outputs:
  - Topology-class table: per class, list of representatives and max within-class
    pairwise |S_a(t) - S_b(t)|.
  - Cross-class S(t) at sample times (does the topology actually distinguish?).
  - Plot: one curve per topology class (representative S(t)), all overlaid.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

REPO_ROOT = Path(__file__).resolve().parents[1]
CSV_DIR = REPO_ROOT / "simulations" / "results" / "bond_isolate"
PARAMS = "J0.0750_gamma0.0500_probe-coherence"
N = 7


def classify_bonds(bonds: list[int]) -> tuple[int, ...]:
    """Multiset of connected-path-lengths in the bond graph."""
    if not bonds:
        return ()
    bs = sorted(bonds)
    components: list[int] = []
    current = [bs[0]]
    for b in bs[1:]:
        if b == current[-1] + 1:
            current.append(b)
        else:
            components.append(len(current))
            current = [b]
    components.append(len(current))
    return tuple(sorted(components))


def csv_path(bonds: list[int]) -> Path:
    bs = sorted(bonds)
    label = f"b{bs[0]}" if len(bs) == 1 else "b" + "-".join(str(b) for b in bs)
    return CSV_DIR / f"N{N}_{label}_{PARAMS}.csv"


def load_S(bonds: list[int]) -> tuple[np.ndarray, np.ndarray] | None:
    p = csv_path(bonds)
    if not p.exists():
        return None
    data = np.loadtxt(p, delimiter=",", skiprows=1)
    return data[:, 0], data[:, -1]


# All bond configurations to test
CONFIGS = [
    # 1-bond (already exist as single-bond runs)
    [0], [1], [2], [3], [4], [5],
    # 2-bond adjacent (path-2)
    [0, 1], [2, 3],
    # 2-bond disjoint (1+1)
    [0, 2], [0, 5],
    # 3-bond path-3
    [0, 1, 2], [2, 3, 4],
    # 3-bond all disjoint (1+1+1)
    [0, 2, 4],
    # 3-bond mixed: path-2 + isolated edge (1, 2)
    [0, 1, 3], [0, 1, 4],
    # 4-bond path-4
    [0, 1, 2, 3], [1, 2, 3, 4],
    # 4-bond mixed: two path-2 (2, 2)
    [0, 1, 3, 4], [0, 1, 4, 5],
    # 4-bond mixed: path-3 + isolated (1, 3)
    [0, 1, 2, 4],
    # 4-bond mixed: path-2 + 2 isolated edges (1, 1, 2)
    [0, 1, 3, 5],
    # 5-bond path-5
    [0, 1, 2, 3, 4], [1, 2, 3, 4, 5],
    # 5-bond mixed: path-4 + isolated edge (1, 4)
    [0, 2, 3, 4, 5], [0, 1, 2, 3, 5],
    # 5-bond mixed: path-2 + path-3 (2, 3)
    [0, 1, 3, 4, 5], [0, 1, 2, 4, 5],
    # full chain (6-bond path-6)
    [0, 1, 2, 3, 4, 5],
]


def main() -> None:
    series: dict[tuple[int, ...], tuple[np.ndarray, np.ndarray, list[int]]] = {}
    missing = []
    for cfg in CONFIGS:
        loaded = load_S(cfg)
        if loaded is None:
            missing.append(cfg)
            continue
        series[tuple(cfg)] = (*loaded, cfg)
    if missing:
        print(f"WARNING: missing CSVs for: {missing}")
        print("Run the bond-isolate compute for these first.")
        return

    # Group by topology class
    by_class: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
    for cfg_key in series:
        cls = classify_bonds(list(cfg_key))
        by_class[cls].append(cfg_key)

    # Within-class consistency table
    print()
    print("# Topology-class consistency: max |S_a(t) - S_b(t)| within each class")
    print()
    print("| Class      | k bonds | Representatives                                      | Within-class max diff | Pred (S_N) |")
    print("|------------|---------|------------------------------------------------------|-----------------------|------------|")
    class_repr: dict[tuple[int, ...], tuple[np.ndarray, np.ndarray]] = {}
    for cls in sorted(by_class.keys(), key=lambda c: (sum(c), len(c), c)):
        members = by_class[cls]
        reps_str = ", ".join("{" + ",".join(str(b) for b in m) + "}" for m in members)
        # Compute max pairwise within-class diff
        if len(members) >= 2:
            t_ref, S_ref = series[members[0]][0], series[members[0]][1]
            max_diff = 0.0
            for m in members[1:]:
                _, S_m, _ = series[m]
                d = float(np.max(np.abs(S_ref - S_m)))
                if d > max_diff:
                    max_diff = d
        else:
            max_diff = float("nan")
        k = sum(cls)
        cls_str = "(" + ", ".join(str(c) for c in cls) + ")"
        pred = "identical" if len(members) >= 2 else "n/a"
        diff_str = f"{max_diff:.2e}" if not np.isnan(max_diff) else "n/a"
        print(f"| {cls_str:10s} | {k}       | {reps_str:52s} | {diff_str:21s} | {pred:10s} |")
        # Store representative for cross-class plot
        t, S, _ = series[members[0]]
        class_repr[cls] = (t, S)
    print()

    # Cross-class comparison: how do classes differ at fixed sample times?
    sample_t_indices = [0, 30, 50, 100, 200, 300]
    print("# Cross-class S(t) at sample times")
    print()
    head = "| Class      | k | " + " | ".join(f"t={class_repr[next(iter(class_repr))][0][i]:.1f}" for i in sample_t_indices) + " |"
    print(head)
    print("|" + "------|" * (len(sample_t_indices) + 2))
    for cls in sorted(class_repr.keys(), key=lambda c: (sum(c), len(c), c)):
        t, S = class_repr[cls]
        k = sum(cls)
        cls_str = "(" + ", ".join(str(c) for c in cls) + ")"
        vals = " | ".join(f"{S[i]:.6f}" for i in sample_t_indices)
        print(f"| {cls_str:10s} | {k} | {vals} |")
    print()

    # Plot: one curve per class (linear + log)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    cmap = plt.cm.viridis
    classes_sorted = sorted(class_repr.keys(), key=lambda c: (sum(c), len(c), c))
    for i, cls in enumerate(classes_sorted):
        t, S = class_repr[cls]
        color = cmap(i / max(len(classes_sorted) - 1, 1))
        members = by_class[cls]
        rep = members[0]
        label = f"{cls}  k={sum(cls)}  rep=[{','.join(str(b) for b in rep)}]"
        ax1.plot(t, S, color=color, label=label, linewidth=1.5)
        ax2.semilogy(t, S, color=color, label=label, linewidth=1.5)
    ax1.set_ylabel("S(t) (spatial-sum coherence)")
    ax1.set_title(f"Multi-bond topology classes at N={N} (ρ_cc probe, J=0.075, γ=0.05)")
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(alpha=0.3)
    ax2.set_xlabel("t")
    ax2.set_ylabel("log S(t)")
    ax2.grid(alpha=0.3, which="both")
    plt.tight_layout()
    out = CSV_DIR / "topology_classes_N7_S_overlay.png"
    plt.savefig(out, dpi=120)
    print(f"Plot saved: {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
