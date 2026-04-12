"""
Cockpit Scaling Analysis
========================
Reads cockpit_scaling CSVs from the C# matrix-free propagator, applies
trajectory sanity gates, drops trivial configurations, runs PCA per
(N, topology, pair), and produces:
  - cockpit_scaling_results.txt   (human-readable scaling table with gate status)
  - cockpit_scaling_results.json  (machine-readable full dump)
  - cockpit_scaling_curve.png     (4-panel scaling plot, far_edge excluded)

Initial state: Bell+(c1,c2) at the chain center, with spectator qubits in
|+>. The sanity gates prevent trivial trajectories from inflating PCA
coverage numbers.

Note on naming: an earlier version of this experiment placed the Bell pair
on the boundary qubits (0,1) and was found to be an initial-state artifact.
That version is archived in ClaudeTasks/Archiv/ as
RESULT_COCKPIT_SCALING_V1_OBSOLETE_initial_state_artifact.md. The current
script and its outputs no longer carry a v1/v2 suffix because there is only
one live version.

Feature order (must match C# ExtractCockpitFeatures):
  phi_plus, phi_minus, psi_plus, psi_minus, purity, svn, concurrence, psi_norm, ph03
"""

import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Configuration ----------------------------------------------------------

FEAT_NAMES = ["Phi+", "Phi-", "Psi+", "Psi-", "Pur", "SvN", "C", "Psi", "ph03"]
PROXY_CANDIDATES = {"C": 6, "Pur": 4, "Psi": 7}  # name -> column index

# Pair classes that the cockpit framework is designed to monitor.
# far_edge and far_leaf pairs are valid classical decoherence trajectories
# but are not informative about the cockpit claim (see COCKPIT_SCALING.md
# Section 6). They are analyzed but flagged cockpit_relevant=False.
COCKPIT_RELEVANT_CLASSES = {"center_bell", "adjacent", "center_leaf"}
RESULTS_DIR = os.path.join(os.path.dirname(__file__), "results", "cockpit_scaling")

# ---- Data loading -----------------------------------------------------------

def load_all_csvs():
    pattern = os.path.join(RESULTS_DIR, "cockpit_scaling_N*_*.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"ERROR: No CSV files found matching {pattern}", file=sys.stderr)
        sys.exit(1)

    records = []
    for fpath in files:
        fname = os.path.basename(fpath)
        m = re.match(r"cockpit_scaling_N(\d+)_(\w+)\.csv", fname)
        if not m:
            continue
        N = int(m.group(1))
        topo = m.group(2)

        lines = open(fpath, "r").readlines()
        data = []
        for line in lines[1:]:
            parts = line.strip().split(",")
            t = float(parts[0])
            pair = parts[1]
            feats = [float(x) for x in parts[2:]]
            data.append((t, pair, feats))

        records.append({"N": N, "topo": topo, "file": fpath, "data": data})
        print(f"  Loaded {fname}: {len(data)} rows")

    return records


# ---- Sanity gates -----------------------------------------------------------

def classify_label(pair_label):
    """Extract pair class from label suffix."""
    if "_center_bell" in pair_label:
        return "center_bell"
    elif "_adjacent" in pair_label:
        return "adjacent"
    elif "_far_edge" in pair_label:
        return "far_edge"
    elif "_center_leaf" in pair_label:
        return "center_leaf"
    elif "_far_leaf" in pair_label:
        return "far_leaf"
    return "unknown"


def check_trajectory_sanity(X, pair_label, N, topo):
    """
    Apply 3 sanity gates to a feature matrix X (n_samples, 9).
    Returns dict with gate results and pass/fail status.
    """
    pclass = classify_label(pair_label)
    stds = X.std(axis=0)
    failures = []

    # Gate 1: Concurrence variation (only required for center_bell)
    conc_std = float(stds[6])
    gate1_required = (pclass == "center_bell")
    gate1_pass = conc_std >= 0.01
    if gate1_required and not gate1_pass:
        failures.append(
            f"SANITY FAIL: N={N} {topo} pair {pair_label}: std(concurrence) = {conc_std:.2e}\n"
            f"  Expected an entangled-observer pair to have rich concurrence dynamics.\n"
            f"  This (N, topology, pair) is DROPPED from the scaling table."
        )

    # Gate 2: Feature richness
    active_mask = stds > 1e-6
    n_active = int(active_mask.sum())
    active_list = [FEAT_NAMES[i] for i in range(9) if active_mask[i]]
    dropped_list = [FEAT_NAMES[i] for i in range(9) if not active_mask[i]]
    gate2_pass = n_active >= 4
    if not gate2_pass:
        failures.append(
            f"TRIVIAL TRAJECTORY: N={N} {topo} pair {pair_label}: only {n_active} of 9 features have std > 1e-6\n"
            f"  Active features: {active_list}\n"
            f"  Dropped features: {dropped_list}\n"
            f"  PCA on this trajectory would amplify numerical noise. DROPPED."
        )

    # Gate 3: Purity range
    pur_range = float(X[:, 4].max() - X[:, 4].min())
    gate3_pass = pur_range >= 0.05
    if not gate3_pass:
        failures.append(
            f"PURITY RANGE FAIL: N={N} {topo} pair {pair_label}: purity range = {pur_range:.4f}\n"
            f"  Purity barely changes over the trajectory. DROPPED."
        )

    overall_pass = (gate1_pass or not gate1_required) and gate2_pass and gate3_pass

    return {
        "pass": overall_pass,
        "gate1_concurrence_std": conc_std,
        "gate1_required": gate1_required,
        "gate1_pass": gate1_pass,
        "gate2_active_features": n_active,
        "gate2_active_list": active_list,
        "gate2_dropped_list": dropped_list,
        "gate2_pass": gate2_pass,
        "gate3_purity_range": pur_range,
        "gate3_pass": gate3_pass,
        "failures": failures,
        "active_mask": active_mask,
    }


# ---- PCA -------------------------------------------------------------------

def run_pca_on_active(X, active_mask):
    """
    Drop near-zero features, standardize, run PCA via SVD.
    Returns var_exp, cum_var, Vt, scores, active_feat_names.
    """
    X_active = X[:, active_mask]
    active_feat_names = [FEAT_NAMES[i] for i, k in enumerate(active_mask) if k]

    mu = X_active.mean(axis=0)
    std = X_active.std(axis=0)
    X_norm = (X_active - mu) / std

    U, S, Vt = np.linalg.svd(X_norm, full_matrices=False)
    var_exp = S**2 / np.sum(S**2)
    cum_var = np.cumsum(var_exp)
    scores = X_norm @ Vt.T

    return var_exp, cum_var, Vt, scores, active_feat_names


def analyze_one(X, active_mask):
    """Run gated PCA on feature matrix, return result dict."""
    var_exp, cum_var, Vt, scores, active_names = run_pca_on_active(X, active_mask)
    n_active = len(active_names)

    n90 = int(np.searchsorted(cum_var, 0.90) + 1)
    n95 = int(np.searchsorted(cum_var, 0.95) + 1)
    n99 = int(np.searchsorted(cum_var, 0.99) + 1)

    # Top loading feature on PC1 (within active features)
    top_loading_idx_active = int(np.argmax(np.abs(Vt[0])))
    top_loading_name = active_names[top_loading_idx_active]
    top_loading_val = float(Vt[0, top_loading_idx_active])

    # Best proxy for PC1 among Concurrence, Purity, Psi-norm (if active)
    pc1_scores = scores[:, 0]
    best_proxy_name = None
    best_proxy_corr = -1
    for pname, pidx in PROXY_CANDIDATES.items():
        if not active_mask[pidx]:
            continue
        feat_col = X[:, pidx]
        feat_std = (feat_col - feat_col.mean()) / (feat_col.std() + 1e-10)
        pc1_std = (pc1_scores - pc1_scores.mean()) / (pc1_scores.std() + 1e-10)
        corr = float(np.abs(np.mean(feat_std * pc1_std)))
        if corr > best_proxy_corr:
            best_proxy_corr = corr
            best_proxy_name = pname

    three_pc = float(cum_var[min(2, len(cum_var) - 1)] * 100)

    return {
        "n90": n90, "n95": n95, "n99": n99,
        "n_active": n_active,
        "active_features": active_names,
        "var_exp": var_exp.tolist(),
        "cum_var": cum_var.tolist(),
        "Vt": Vt.tolist(),
        "top_loading_name": top_loading_name,
        "top_loading_val": top_loading_val,
        "best_proxy_name": best_proxy_name,
        "best_proxy_corr": best_proxy_corr,
        "pc1_pct": float(var_exp[0] * 100),
        "pc2_pct": float(var_exp[1] * 100) if len(var_exp) > 1 else 0,
        "pc3_pct": float(var_exp[2] * 100) if len(var_exp) > 2 else 0,
        "three_pc_pct": three_pc,
    }


# ---- Main analysis ----------------------------------------------------------

def main():
    print("=" * 60)
    print("COCKPIT SCALING ANALYSIS")
    print("=" * 60)
    print("Initial state: Bell+(c1, c2) tensor |+>^(n-2), c1 = (n-1)/2")
    print()

    records = load_all_csvs()
    print()

    all_results = {}  # key -> result dict (analyzed configs)
    dropped = []      # list of drop-info dicts
    gate_summary = [] # per-config gate status for the report

    for rec in records:
        N, topo = rec["N"], rec["topo"]
        pairs_in_data = sorted(set(row[1] for row in rec["data"]))

        for pair in pairs_in_data:
            pair_rows = [row for row in rec["data"] if row[1] == pair]
            X = np.array([row[2] for row in pair_rows])
            pclass = classify_label(pair)

            # Run sanity gates
            sanity = check_trajectory_sanity(X, pair, N, topo)
            gate_entry = {
                "N": N, "topo": topo, "pair": pair, "pclass": pclass,
                "cockpit_relevant": pclass in COCKPIT_RELEVANT_CLASSES,
                "gate1_pass": sanity["gate1_pass"],
                "gate1_required": sanity["gate1_required"],
                "gate1_conc_std": sanity["gate1_concurrence_std"],
                "gate2_pass": sanity["gate2_pass"],
                "gate2_n_active": sanity["gate2_active_features"],
                "gate3_pass": sanity["gate3_pass"],
                "gate3_pur_range": sanity["gate3_purity_range"],
                "overall": sanity["pass"],
            }
            gate_summary.append(gate_entry)

            if not sanity["pass"]:
                for msg in sanity["failures"]:
                    print(msg)
                dropped.append({
                    "N": N, "topo": topo, "pair": pair, "pclass": pclass,
                    "cockpit_relevant": pclass in COCKPIT_RELEVANT_CLASSES,
                    "failures": sanity["failures"],
                    "gate1_conc_std": sanity["gate1_concurrence_std"],
                    "gate2_n_active": sanity["gate2_active_features"],
                    "gate3_pur_range": sanity["gate3_purity_range"],
                })
                continue

            # PCA with feature dropping
            result = analyze_one(X, sanity["active_mask"])
            result["N"] = N
            result["topo"] = topo
            result["pair"] = pair
            result["pclass"] = pclass
            result["cockpit_relevant"] = pclass in COCKPIT_RELEVANT_CLASSES
            result["gate1_conc_std"] = sanity["gate1_concurrence_std"]
            result["gate3_pur_range"] = sanity["gate3_purity_range"]
            key = f"N{N}_{topo}_{pair}"
            all_results[key] = result

    print()

    # ---- Output 1: text report ----
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("COCKPIT SCALING ANALYSIS")
    report_lines.append("=" * 60)
    report_lines.append("Initial state: Bell+(c1, c2) tensor |+>^(n-2), c1 = (n-1)/2")
    report_lines.append("")

    for topo_label in ["chain", "star"]:
        topo_keys = sorted(
            [k for k, v in all_results.items() if v["topo"] == topo_label],
            key=lambda k: (all_results[k]["N"], all_results[k]["pair"])
        )
        # Also gather dropped entries for this topology
        topo_dropped = [d for d in dropped if d["topo"] == topo_label]

        if not topo_keys and not topo_dropped:
            continue

        report_lines.append(f"{topo_label.upper()} TOPOLOGY")
        report_lines.append(f"{'':20s} n_active n95  3-PC% | gate status     | PC1 best proxy | top loading")

        # Merge analyzed + dropped, sorted by (N, pair)
        all_entries = []
        for k in topo_keys:
            r = all_results[k]
            all_entries.append((r["N"], r["pair"], "analyzed", r))
        for d in topo_dropped:
            all_entries.append((d["N"], d["pair"], "dropped", d))
        all_entries.sort(key=lambda x: (x[0], x[1]))

        for N_val, pair, status, data in all_entries:
            pclass = classify_label(pair)
            if status == "analyzed":
                r = data
                line = (f"N={r['N']:<3d} {pclass:15s}"
                        f"  {r['n_active']:3d}     {r['n95']:3d}  {r['three_pc_pct']:5.1f}"
                        f" | PASS            "
                        f"| {r['best_proxy_name'] or 'N/A':14s}"
                        f" | {r['top_loading_name']}")
            else:
                d = data
                which = []
                if d.get("gate2_n_active", 9) < 4:
                    which.append(f"Gate2: {d['gate2_n_active']} feats")
                if d.get("gate3_pur_range", 1) < 0.05:
                    which.append(f"Gate3: pur_range={d['gate3_pur_range']:.4f}")
                if d.get("gate1_conc_std", 1) < 0.01 and pclass == "center_bell":
                    which.append(f"Gate1: conc_std={d['gate1_conc_std']:.2e}")
                reason = ", ".join(which) if which else "see log"
                line = f"N={N_val:<3d} {pclass:15s}  ---     ---  ---   | DROPPED ({reason})"
            report_lines.append(line)
        report_lines.append("")

    # Dropped configurations section
    if dropped:
        report_lines.append("=" * 60)
        report_lines.append("DROPPED CONFIGURATIONS")
        report_lines.append("=" * 60)
        for d in dropped:
            for msg in d["failures"]:
                report_lines.append(msg)
        report_lines.append("")

    # Scaling trend (center_bell only)
    report_lines.append("=" * 60)
    report_lines.append("SCALING TREND (CENTER_BELL pairs only)")
    report_lines.append("=" * 60)
    report_lines.append(f"{'':16s}  Chain center_bell       Star center_bell")

    all_N = sorted(set(
        list(set(v["N"] for v in all_results.values()))
        + [d["N"] for d in dropped]
    ))
    center_bell_chain = {}
    center_bell_star = {}
    for k, v in all_results.items():
        if v["pclass"] != "center_bell":
            continue
        if v["topo"] == "chain":
            center_bell_chain[v["N"]] = v
        elif v["topo"] == "star":
            center_bell_star[v["N"]] = v

    for N_val in all_N:
        chain_str = "---"
        star_str = "---"
        if N_val in center_bell_chain:
            v = center_bell_chain[N_val]
            chain_str = f"{v['three_pc_pct']:.1f}% (n95={v['n95']})"
        if N_val in center_bell_star:
            v = center_bell_star[N_val]
            star_str = f"{v['three_pc_pct']:.1f}% (n95={v['n95']})"
        report_lines.append(f"N={N_val:<14d}  {chain_str:24s} {star_str}")

    # Constrained verdict
    report_lines.append("")
    report_lines.append("VERDICT (constrained):")
    center_bell_analyzed = [v for v in all_results.values() if v["pclass"] == "center_bell"]
    center_bell_dropped_check = [d for d in dropped if classify_label(d["pair"]) == "center_bell"]

    if center_bell_dropped_check:
        report_lines.append("  PIPELINE FAILURE OR REAL DEGRADATION: center_bell pair(s) were dropped by sanity gates.")
        for d in center_bell_dropped_check:
            report_lines.append(f"    N={d['N']} {d['topo']}: {d['failures'][0].split(chr(10))[0]}")
    elif center_bell_analyzed:
        coverages = [v["three_pc_pct"] for v in center_bell_analyzed]
        min_cov = min(coverages)
        max_cov = max(coverages)
        n_analyzed = len(center_bell_analyzed)
        n_dropped_total = len(dropped)
        if min_cov >= 85:
            verdict = "COCKPIT SCALES"
        elif min_cov >= 70:
            verdict = "GRADUAL DEGRADATION"
        else:
            verdict = "COCKPIT N-BOUNDED"
        report_lines.append(f"  {verdict}")
        report_lines.append(f"  center_bell: {n_analyzed} configs analyzed, coverage [{min_cov:.1f}%, {max_cov:.1f}%]")
        report_lines.append(f"  {n_dropped_total} configuration(s) dropped total (expected for far_edge/far_leaf)")
    else:
        report_lines.append("  NO DATA: no center_bell configurations were analyzed.")

    # Cockpit-relevant configurations table
    report_lines.append("")
    report_lines.append("=" * 60)
    report_lines.append("COCKPIT-RELEVANT CONFIGURATIONS (cockpit_relevant=True)")
    report_lines.append("Excluded classes (far_edge, far_leaf) remain in the JSON under cockpit_relevant=False.")
    report_lines.append("=" * 60)
    report_lines.append(f"{'':20s} n_active n95  3-PC%  PC1%  best proxy")
    cr_entries = sorted(
        [v for v in all_results.values() if v.get("cockpit_relevant")],
        key=lambda v: (v["topo"], v["N"], v["pclass"])
    )
    for v in cr_entries:
        report_lines.append(
            f"{v['topo']:5s} N={v['N']:<3d} {v['pclass']:15s}"
            f"  {v['n_active']:3d}   {v['n95']:3d}  {v['three_pc_pct']:5.1f}"
            f"  {v['pc1_pct']:5.1f}  {v['best_proxy_name'] or 'N/A'}"
        )

    report_text = "\n".join(report_lines)
    print(report_text)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    txt_path = os.path.join(RESULTS_DIR, "cockpit_scaling_results.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write(report_text + "\n")
    print(f"\nSaved: {txt_path}")

    # ---- Output 2: JSON dump ----
    json_data = {
        "analyzed": all_results,
        "dropped": dropped,
        "gate_summary": gate_summary,
    }
    json_path = os.path.join(RESULTS_DIR, "cockpit_scaling_results.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_data, f, indent=2, default=lambda x: x.tolist() if hasattr(x, 'tolist') else x)
    print(f"Saved: {json_path}")

    # ---- Output 3: scaling plot ----
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Cockpit Scaling: Bell Pair at Center", fontsize=14, fontweight="bold")

    pair_styles = {
        "center_bell": "-",
        "adjacent": "--",
        "center_leaf": "--",
    }
    topo_colors = {"chain": "tab:blue", "star": "tab:orange"}

    # Panel A: 3-PC cumulative variance vs N (center_bell + adjacent only)
    ax = axes[0, 0]
    for topo in ["chain", "star"]:
        for pclass, ls in pair_styles.items():
            Ns, vals = [], []
            for k, v in sorted(all_results.items(), key=lambda x: x[1]["N"]):
                if v["topo"] != topo or v["pclass"] != pclass:
                    continue
                Ns.append(v["N"])
                vals.append(v["three_pc_pct"])
            if Ns:
                ax.plot(Ns, vals, ls, color=topo_colors[topo],
                        label=f"{topo} {pclass}", marker="o", markersize=4)
    ax.axhline(85, color="gray", linestyle="-.", linewidth=0.8, label="85% threshold")
    ax.text(0.02, 0.02, "only cockpit_relevant classes shown (center_bell, adjacent, center_leaf)",
            transform=ax.transAxes, fontsize=7, fontstyle="italic", color="gray")
    ax.set_xlabel("N (qubits)")
    ax.set_ylabel("3-PC cumulative variance (%)")
    ax.set_title("Panel A: 3-PC Coverage vs N")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    # Panel B: n95 vs N
    ax = axes[0, 1]
    for topo in ["chain", "star"]:
        for pclass, ls in pair_styles.items():
            Ns, vals = [], []
            for k, v in sorted(all_results.items(), key=lambda x: x[1]["N"]):
                if v["topo"] != topo or v["pclass"] != pclass:
                    continue
                Ns.append(v["N"])
                vals.append(v["n95"])
            if Ns:
                ax.plot(Ns, vals, ls, color=topo_colors[topo],
                        label=f"{topo} {pclass}", marker="s", markersize=4)
    if all_N:
        all_N_arr = np.array(sorted(all_N))
        ax.plot(all_N_arr, all_N_arr, "k--", alpha=0.4, label="n95 = N")
    ax.set_xlabel("N (qubits)")
    ax.set_ylabel("n95")
    ax.set_title("Panel B: n95 vs N")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    # Panel C: PC1 variance fraction vs N
    ax = axes[1, 0]
    for topo in ["chain", "star"]:
        for pclass, ls in pair_styles.items():
            Ns, vals = [], []
            for k, v in sorted(all_results.items(), key=lambda x: x[1]["N"]):
                if v["topo"] != topo or v["pclass"] != pclass:
                    continue
                Ns.append(v["N"])
                vals.append(v["pc1_pct"])
            if Ns:
                ax.plot(Ns, vals, ls, color=topo_colors[topo],
                        label=f"{topo} {pclass}", marker="^", markersize=4)
    ax.set_xlabel("N (qubits)")
    ax.set_ylabel("PC1 variance (%)")
    ax.set_title("Panel C: PC1 Variance vs N")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    # Panel D: PC1 best proxy correlation vs N (center_bell only)
    ax = axes[1, 1]
    proxy_markers = {"C": "o", "Pur": "s", "Psi": "^"}
    proxy_colors_map = {"C": "tab:red", "Pur": "tab:green", "Psi": "tab:purple"}
    plotted_labels = set()
    for k, v in sorted(all_results.items(), key=lambda x: (x[1]["topo"], x[1]["N"])):
        if v["pclass"] != "center_bell":
            continue
        pname = v["best_proxy_name"]
        if pname is None:
            continue
        topo = v["topo"]
        marker = proxy_markers.get(pname, "x")
        facecolor = proxy_colors_map.get(pname, "gray")
        edgecolor = topo_colors[topo]
        label_str = f"{topo}/{pname}"
        label = label_str if label_str not in plotted_labels else None
        plotted_labels.add(label_str)
        ax.scatter(v["N"], v["best_proxy_corr"], marker=marker,
                   c=facecolor, edgecolors=edgecolor, s=60, linewidths=1.5,
                   label=label, zorder=3)
    ax.set_xlabel("N (qubits)")
    ax.set_ylabel("|correlation with PC1|")
    ax.set_title("Panel D: PC1 Best Proxy vs N (center_bell only)")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    png_path = os.path.join(RESULTS_DIR, "cockpit_scaling_curve.png")
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {png_path}")


if __name__ == "__main__":
    main()
