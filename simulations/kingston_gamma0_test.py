"""
Phase 1: Kingston gamma_0 consistency test (zero QPU cost).

Tests whether per-qubit dephasing on IBM Kingston is consistent with the
Primordial Gamma Constant hypothesis: gamma_0 is a framework constant, and
per-qubit gamma_phi is related to mode-exposure |a_B|^2 in the chip's cavity
structure.

Steps:
1. Parse 8 Kingston calibration CSVs -> gamma_phi(i) median and CV per qubit
2. Query ibm_kingston backend -> frequencies, coupling map, RZZ gate data
3. Extract or estimate J-coupling per edge
4. Select linear sub-chain of 5-7 qubits (prefer CPsi-crossing qubits as anchors)
5. Build tridiagonal single-excitation Hamiltonian, diagonalize
6. Compute gamma_0_inferred under multiple interpretations, report verdict

Outputs:
  simulations/results/kingston_backend_properties.json
  simulations/results/gamma0_consistency_test.txt

Rules from task: numbers only from script output, no head calculation.
"""

import os
import sys
import json
import csv
import math
from pathlib import Path
from collections import defaultdict

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

# --- Paths ---------------------------------------------------------------

REPO = Path(__file__).resolve().parent.parent
CALIB_DIR = REPO / "ClaudeTasks" / "IBM_R2_calibrations"
RESULTS_DIR = REPO / "simulations" / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
BACKEND_JSON = RESULTS_DIR / "kingston_backend_properties.json"
REPORT_TXT = RESULTS_DIR / "gamma0_consistency_test.txt"

# Anchor qubits from existing Kingston CPsi=1/4 crossing experiment (commit 20ef49e)
CPSI_CROSSING_QUBITS = [14, 15, 124, 125]


# --- Step 1: Parse calibration CSVs --------------------------------------

def parse_partner_field(s: str) -> dict:
    if not s or s.strip() in ("", '""'):
        return {}
    out = {}
    for part in s.split(";"):
        part = part.strip().strip('"')
        if ":" not in part:
            continue
        p_str, v_str = part.split(":", 1)
        try:
            out[int(p_str)] = float(v_str)
        except ValueError:
            continue
    return out


def load_calibrations():
    """Load all calibration CSVs. Returns per-qubit stats and per-edge stats."""
    snaps = sorted(CALIB_DIR.glob("ibm_kingston_calibrations_*.csv"))
    qubit_series = defaultdict(lambda: {"T1": [], "T2": [], "gamma_phi": [], "operational": []})
    edge_series = defaultdict(lambda: {"cz_err": [], "rzz_err": [], "gate_len_ns": []})

    for snap in snaps:
        with open(snap, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                q = int(row["Qubit"])
                t1_raw = row["T1 (us)"].strip()
                t2_raw = row["T2 (us)"].strip()
                if not t1_raw or not t2_raw:
                    qubit_series[q]["operational"].append(False)
                    continue
                t1 = float(t1_raw)
                t2 = float(t2_raw)
                gamma_phi = 1.0 / t2 - 1.0 / (2.0 * t1)
                qubit_series[q]["T1"].append(t1)
                qubit_series[q]["T2"].append(t2)
                qubit_series[q]["gamma_phi"].append(gamma_phi)
                qubit_series[q]["operational"].append(row["Operational"] == "Yes")

                cz = parse_partner_field(row["CZ error"])
                rzz = parse_partner_field(row["RZZ error"])
                glen = parse_partner_field(row["Gate length (ns)"])
                partners = set(cz.keys()) | set(rzz.keys()) | set(glen.keys())
                for partner in partners:
                    edge = tuple(sorted((q, partner)))
                    if partner in cz:
                        edge_series[edge]["cz_err"].append(cz[partner])
                    if partner in rzz:
                        edge_series[edge]["rzz_err"].append(rzz[partner])
                    if partner in glen:
                        edge_series[edge]["gate_len_ns"].append(glen[partner])

    qubit_stats = {}
    for q, series in qubit_series.items():
        if not series["gamma_phi"]:
            continue
        g = np.array(series["gamma_phi"])
        t2 = np.array(series["T2"])
        t1 = np.array(series["T1"])
        med_g = float(np.median(g))
        cv_g = float(np.std(g) / abs(med_g)) if med_g != 0 else float("inf")
        qubit_stats[q] = {
            "T1_median": float(np.median(t1)),
            "T2_median": float(np.median(t2)),
            "gamma_phi_median": med_g,
            "gamma_phi_min": float(g.min()),
            "gamma_phi_max": float(g.max()),
            "gamma_phi_std": float(np.std(g)),
            "gamma_phi_cv": cv_g,
            "n_snapshots": int(len(g)),
            "always_operational": bool(all(series["operational"])),
        }

    edge_stats = {}
    for edge, series in edge_series.items():
        edge_stats[edge] = {
            k: float(np.median(v)) if v else None
            for k, v in series.items()
        }

    return qubit_stats, edge_stats, len(snaps)


# --- Step 2: Query backend ------------------------------------------------

def query_backend(use_cache: bool = True) -> dict | None:
    """Query ibm_kingston backend. Returns dict, None on failure.

    Token is read from the IBM_QUANTUM_TOKEN environment variable. If not set,
    the backend query is skipped and the script proceeds with CSV-only data
    (or the cached JSON if present).
    """
    if use_cache and BACKEND_JSON.exists():
        print(f"  Using cached backend properties: {BACKEND_JSON.name}")
        with open(BACKEND_JSON, "r", encoding="utf-8") as f:
            return json.load(f)

    token = os.environ.get("IBM_QUANTUM_TOKEN")
    if not token:
        print("  IBM_QUANTUM_TOKEN not set; skipping backend query.")
        print("  To populate kingston_backend_properties.json, run with:")
        print("    IBM_QUANTUM_TOKEN=<your_token> python simulations/kingston_gamma0_test.py")
        return None

    try:
        from qiskit_ibm_runtime import QiskitRuntimeService
    except ImportError:
        print("  WARNING: qiskit-ibm-runtime not installed. Skipping backend query.")
        return None

    print(f"  Token: ...{token[-8:]}")
    service = None
    for channel in ("ibm_quantum_platform", "ibm_quantum"):
        try:
            service = QiskitRuntimeService(channel=channel, token=token)
            print(f"  Connected via channel={channel}")
            break
        except Exception as e:
            print(f"  channel={channel} failed: {e}")
            try:
                QiskitRuntimeService.save_account(channel=channel, token=token, overwrite=True)
                service = QiskitRuntimeService(channel=channel)
                print(f"  Connected via channel={channel} after save_account")
                break
            except Exception as e2:
                print(f"  save_account fallback also failed: {e2}")

    if service is None:
        return None

    try:
        backend = service.backend("ibm_kingston")
    except Exception as e:
        print(f"  Backend ibm_kingston not accessible: {e}")
        return None

    props = backend.properties()
    target = backend.target
    cmap = backend.coupling_map

    data = {
        "backend_name": backend.name,
        "num_qubits": backend.num_qubits,
        "coupling_map_directed": [list(e) for e in cmap.get_edges()],
        "qubits": {},
        "edges": {},
        "target_operations": sorted(list(target.operation_names)),
    }

    for q in range(backend.num_qubits):
        entry = {}
        try:
            entry["frequency_GHz"] = props.frequency(q) / 1e9
        except Exception:
            pass
        try:
            entry["T1_us"] = props.t1(q) * 1e6
        except Exception:
            pass
        try:
            entry["T2_us"] = props.t2(q) * 1e6
        except Exception:
            pass
        for prop_name in ("anharmonicity",):
            try:
                val = props.qubit_property(q, prop_name)
                if val is not None:
                    entry[prop_name] = val[0] if hasattr(val, "__getitem__") else val
            except Exception:
                pass
        data["qubits"][str(q)] = entry

    for op_name in ("rzz", "cz", "ecr"):
        if op_name in target.operation_names:
            try:
                props_map = target[op_name]
            except Exception:
                continue
            for qargs, inst_props in props_map.items():
                if qargs is None:
                    continue
                try:
                    edge_key = ",".join(str(q) for q in sorted(qargs))
                except TypeError:
                    continue
                data["edges"].setdefault(edge_key, {})
                if inst_props is not None:
                    duration = getattr(inst_props, "duration", None)
                    error = getattr(inst_props, "error", None)
                    if duration is not None:
                        data["edges"][edge_key][f"{op_name}_duration_s"] = duration
                    if error is not None:
                        data["edges"][edge_key][f"{op_name}_error"] = error

    with open(BACKEND_JSON, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
    print(f"  Saved {BACKEND_JSON.name} ({len(data['qubits'])} qubits, {len(data['edges'])} edges)")
    return data


# --- Step 3: Extract J couplings -----------------------------------------

def extract_J_per_edge(backend_data, edge_stats):
    """
    Estimate J per edge. Reports which route produced each estimate.

    Route (b): From RZZ native gate. Heron r2 native gate is RZZ(theta) with
    variable theta; the characterized gate time corresponds to a specific
    reference angle. Following IBM Heron r2 docs, the listed gate duration is
    for the calibrated native gate at a fixed nominal angle (Heron r2 typically
    calibrates at theta = pi/2). With H = J_ZZ * Z Z active during the gate,
    U = exp(-i * J_ZZ * t_gate * Z Z), so the effective angle per Z Z generator
    is phi = J_ZZ * t_gate; the standard RZZ convention is RZZ(theta) =
    exp(-i * theta/2 * Z Z), which gives theta = 2 * J_ZZ * t_gate, hence
    J_ZZ = theta / (2 * t_gate). Assume theta = pi/2.

    All J values reported in MHz (angular frequency / 2pi).

    Returns (J_dict, source_dict).
    """
    J_est = {}
    sources = {}

    for edge, stats in edge_stats.items():
        gate_len_ns = stats.get("gate_len_ns")
        if gate_len_ns and gate_len_ns > 0:
            t_s = gate_len_ns * 1e-9
            theta = math.pi / 2  # assumed nominal native angle
            J_rad_per_s = theta / (2.0 * t_s)
            J_MHz = J_rad_per_s / (2.0 * math.pi) / 1e6
            J_est[edge] = J_MHz
            sources[edge] = "route_b_gate_duration_CSV"

    if backend_data is not None:
        for edge_key, edata in backend_data.get("edges", {}).items():
            qs = tuple(int(x) for x in edge_key.split(","))
            if len(qs) != 2:
                continue
            edge = tuple(sorted(qs))
            for op in ("rzz", "cz", "ecr"):
                dur = edata.get(f"{op}_duration_s")
                if dur and dur > 0:
                    theta = math.pi / 2
                    J_rad_per_s = theta / (2.0 * dur)
                    J_MHz = J_rad_per_s / (2.0 * math.pi) / 1e6
                    if edge not in J_est:
                        J_est[edge] = J_MHz
                        sources[edge] = f"route_b_gate_duration_backend_{op}"
                    break

    return J_est, sources


# --- Step 4: Select linear sub-chain -------------------------------------

def find_linear_path(edges, qubit_stats, length=5, anchors=None, max_cv=None):
    """
    Find a linear path of `length` qubits. Linear means each interior qubit has
    exactly two neighbors in the path. Prefers high median T2 and anchor inclusion.
    If max_cv given, only qubits with gamma_phi_cv <= max_cv are considered.
    """
    anchors = set(anchors or [])
    adj = defaultdict(set)
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)

    t2_values = np.array([s["T2_median"] for s in qubit_stats.values()])
    median_T2 = float(np.median(t2_values))

    def ok(q):
        if q not in qubit_stats:
            return False
        s = qubit_stats[q]
        if not s.get("always_operational", False):
            return False
        if max_cv is not None and s["gamma_phi_cv"] > max_cv:
            return False
        return True

    best_path = None
    best_score = -np.inf

    starts = sorted(
        [q for q in qubit_stats.keys() if ok(q)],
        key=lambda q: (-int(q in anchors), -qubit_stats[q]["T2_median"]),
    )

    for start in starts:
        stack = [(start, (start,))]
        while stack:
            node, path = stack.pop()
            if len(path) == length:
                avg_t2 = float(np.mean([qubit_stats[q]["T2_median"] for q in path]))
                cv_penalty = max(qubit_stats[q]["gamma_phi_cv"] for q in path)
                anchor_bonus = 500.0 * len(anchors & set(path))
                above_median_bonus = 50.0 * sum(
                    1 for q in path if qubit_stats[q]["T2_median"] >= median_T2
                )
                score = avg_t2 + anchor_bonus + above_median_bonus - 100.0 * cv_penalty
                if score > best_score:
                    best_score = score
                    best_path = list(path)
                continue
            for nbr in adj[node]:
                if nbr in path:
                    continue
                if not ok(nbr):
                    continue
                stack.append((nbr, path + (nbr,)))

    return best_path, best_score, median_T2


# --- Step 5: Build tridiagonal Hamiltonian -------------------------------

def build_H(path, J_est, uniform=False):
    N = len(path)
    H = np.zeros((N, N))
    edge_J_used = []
    for i in range(N - 1):
        edge = tuple(sorted((path[i], path[i + 1])))
        J = 1.0 if uniform else J_est.get(edge, None)
        if J is None:
            J = 1.0
        H[i, i + 1] = J
        H[i + 1, i] = J
        edge_J_used.append((edge, J))
    evs, evecs = np.linalg.eigh(H)
    return H, evs, evecs, edge_J_used


# --- Step 6: gamma_0 consistency test ------------------------------------

def gamma0_test(path, eigvals, eigvecs, qubit_stats):
    """
    For each qubit i in path, compute gamma_0_inferred via several
    interpretations. See report section for physical meaning.
    """
    N = len(path)
    per_qubit = []
    for i, q in enumerate(path):
        gp = qubit_stats[q]["gamma_phi_median"]
        amp_sq = np.abs(eigvecs[i, :]) ** 2  # |psi_k(i)|^2 for all modes k
        max_amp = float(np.max(amp_sq))
        ipr = float(np.sum(amp_sq ** 2))
        # Participation number M_i = 1/IPR -> effective number of modes qubit i sits in
        participation = 1.0 / ipr if ipr > 0 else float("inf")
        # Slowest non-trivial mode amplitude at site i
        # (by convention: k=0 has smallest eigenvalue magnitude for our tridiagonal)
        # Use smallest |eigval| for "slowest" (most dangerous) mode
        slow_idx = int(np.argmin(np.abs(eigvals)))
        slow_amp = float(amp_sq[slow_idx])

        per_qubit.append({
            "qubit": q,
            "gamma_phi": gp,
            "T2_us": qubit_stats[q]["T2_median"],
            "T1_us": qubit_stats[q]["T1_median"],
            "gamma_phi_cv": qubit_stats[q]["gamma_phi_cv"],
            "max_amp_sq": max_amp,
            "IPR": ipr,
            "participation": participation,
            "slow_mode_amp_sq": slow_amp,
            "amp_sq_per_mode": amp_sq.tolist(),
            # Interpretations of gamma_0:
            # A: mode-exposure, pick max |psi_k(i)|^2 (most weighted mode)
            "gamma_0_A_max_mode": gp / max_amp if max_amp > 0 else None,
            # B: use slowest mode amplitude
            "gamma_0_B_slow_mode": gp / slow_amp if slow_amp > 0 else None,
            # C: uniform sum over modes (sum |psi_k(i)|^2 = 1 by orthonormality, trivial)
            "gamma_0_C_total_participation": gp,
            # D: IPR-weighted (gamma_phi * IPR), interpretation: if gamma_phi is
            # proportional to localization, high IPR means more localized
            "gamma_0_D_IPR_weighted": gp * ipr,
        })
    return per_qubit


# --- Reporting ------------------------------------------------------------

def write_report(report_lines, path_info, qubit_stats, edge_stats, backend_data,
                 J_est, J_src, path, result, uniform_result, undirected_edges,
                 sweep_results=None):
    g_all = np.array([s["gamma_phi_median"] for s in qubit_stats.values()])
    cv_all = np.array([s["gamma_phi_cv"] for s in qubit_stats.values()])

    out = []
    out.append("=" * 78)
    out.append("Phase 1: IBM Kingston gamma_0 Consistency Test")
    out.append("Task: TASK_GAMMA0_HARDWARE_TEST.md")
    out.append("=" * 78)
    out.append("")
    out.append("1) gamma_phi characterization (8 calibration snapshots, Apr 12-19)")
    out.append("-" * 78)
    out.append(f"  Qubits with data:           {len(qubit_stats)}")
    out.append(f"  gamma_phi range (/us):      [{g_all.min():.6f}, {g_all.max():.6f}]")
    out.append(f"  gamma_phi max/min ratio:    {g_all.max()/g_all.min():.1f}x")
    out.append(f"  gamma_phi median:           {np.median(g_all):.6f} /us")
    out.append(f"  gamma_phi mean +/- std:     {g_all.mean():.6f} +/- {g_all.std():.6f} /us")
    out.append(f"  Median CV (8 snapshots):    {np.median(cv_all):.1%}")
    out.append(f"  TLS outliers (CV > 50%):    {int((cv_all > 0.5).sum())} qubits")
    out.append("")
    out.append("2) Coupling map")
    out.append("-" * 78)
    out.append(f"  Edges from CSV partner fields: {len(edge_stats)} undirected")
    if backend_data is not None:
        bd_edges = set(tuple(sorted(e)) for e in backend_data.get("coupling_map_directed", []))
        out.append(f"  Edges from backend.coupling_map: {len(bd_edges)} undirected")
    else:
        out.append("  Backend query unavailable; relying on CSV-derived edges only.")
    out.append("")
    out.append("3) J-coupling extraction")
    out.append("-" * 78)
    src_counts = defaultdict(int)
    for s in J_src.values():
        src_counts[s] += 1
    for s, c in sorted(src_counts.items()):
        out.append(f"  Source={s}: {c} edges")
    J_values = np.array(list(J_est.values()))
    if len(J_values) > 0:
        out.append(f"  J [MHz] range: [{J_values.min():.3f}, {J_values.max():.3f}]")
        out.append(f"  J [MHz] median: {np.median(J_values):.3f}")
        out.append(f"  J [MHz] std: {J_values.std():.3f}")
    out.append("")
    out.append("  NOTE: J extracted from gate duration assuming native RZZ angle = pi/2")
    out.append("  (Heron r2 nominal calibration). J_ZZ = theta / (2 * t_gate).")
    out.append("  Actual coupling strengths during idle (between gates) are NOT accessible")
    out.append("  from calibration data; residual idle ZZ is typically << gate-active J.")
    out.append("  This is a limitation: the test uses gate-active J as a proxy for the")
    out.append("  coupling that determines chain mode structure during a dynamics experiment.")
    out.append("")
    out.append("4) Sub-chain selection")
    out.append("-" * 78)
    out.append(f"  Linear path length: {len(path)}")
    out.append(f"  Path: {path}")
    out.append(f"  Selection score: {path_info['score']:.2f}")
    out.append(f"  Anchor CPsi-crossing qubits: {CPSI_CROSSING_QUBITS}")
    out.append(f"  Anchors included in path: {sorted(set(CPSI_CROSSING_QUBITS) & set(path))}")
    out.append(f"  T2 median (all qubits):    {path_info['median_T2']:.1f} us")
    for i, q in enumerate(path):
        s = qubit_stats[q]
        edge = tuple(sorted((path[i], path[i+1]))) if i < len(path)-1 else None
        J = J_est.get(edge, None) if edge else None
        J_str = f"J_{path[i]}-{path[i+1]}={J:.3f} MHz" if (J is not None and edge) else ""
        out.append(f"    Q{q:3d}: T1={s['T1_median']:7.1f} us, T2={s['T2_median']:7.1f} us, "
                   f"gamma_phi={s['gamma_phi_median']:.5f} /us, CV={s['gamma_phi_cv']:5.1%}  {J_str}")
    out.append("")
    out.append("5) Single-excitation Hamiltonian eigenmodes")
    out.append("-" * 78)
    eigvals_ext = np.array(result["eigvals"])
    eigvals_uni = np.array(uniform_result["eigvals"])
    out.append("  (a) Extracted J per edge:")
    out.append(f"      Eigenvalues [MHz]: {[f'{e:.4f}' for e in eigvals_ext]}")
    out.append("  (b) Uniform J=1 (unitless):")
    out.append(f"      Eigenvalues:       {[f'{e:.4f}' for e in eigvals_uni]}")
    out.append("")
    out.append("  |psi_k(i)|^2 matrix (rows=qubits, cols=modes, extracted J):")
    evecs = np.array(result["eigvecs"])
    amp_sq = np.abs(evecs) ** 2
    header = "    " + "Q".ljust(6) + "".join([f"mode{k}".rjust(10) for k in range(len(path))])
    out.append(header)
    for i, q in enumerate(path):
        row = "    " + f"Q{q}".ljust(6) + "".join([f"{amp_sq[i,k]:10.4f}" for k in range(len(path))])
        out.append(row)
    out.append("")
    out.append("6) gamma_0 inference test")
    out.append("-" * 78)
    interpretations = [
        ("A_max_mode", "gamma_phi(i) / max_k |psi_k(i)|^2",
         "If qubit i is mostly exposed via its dominant mode."),
        ("B_slow_mode", "gamma_phi(i) / |psi_{k_slow}(i)|^2",
         "If the slowest chain mode (smallest |E|) dominates long-time dephasing."),
        ("C_total_participation", "gamma_phi(i)  (sum |psi_k(i)|^2 = 1)",
         "Trivial baseline: framework predicts uniform gamma_0 across all i."),
        ("D_IPR_weighted", "gamma_phi(i) * IPR(i) = gamma_phi * sum_k |psi_k(i)|^4",
         "Localized qubits (high IPR) amplified; delocalized ones suppressed."),
    ]
    per_qubit = result["per_qubit"]
    for key, formula, desc in interpretations:
        name = f"gamma_0_{key}"
        vals = np.array([r[name] for r in per_qubit], dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        out.append(f"  Interpretation {key}:")
        out.append(f"    Formula: {formula}")
        out.append(f"    Meaning: {desc}")
        out.append(f"    Values [/us]: {[f'{v:.5f}' for v in vals]}")
        out.append(f"    Median:       {np.median(vals):.6f} /us")
        out.append(f"    Mean +/- std: {vals.mean():.6f} +/- {vals.std():.6f} /us")
        cv = vals.std() / abs(vals.mean()) if vals.mean() != 0 else float("inf")
        ratio = vals.max() / vals.min() if vals.min() > 0 else float("inf")
        out.append(f"    CV:           {cv:.1%}")
        out.append(f"    max/min:      {ratio:.2f}x")
        out.append("")
    if sweep_results:
        out.append("6b) Multi-configuration path sweep")
        out.append("-" * 78)
        out.append("  Test repeated for 4 sub-chain configurations. Best-case max/min")
        out.append("  across configurations indicates how sensitive the test is to path choice.")
        out.append("")
        header_cols = ("config", "path", "best_interp", "median_gamma_0", "CV", "max/min")
        out.append(f"  {'config':<24}  {'path':<30}  {'interp':<28}  {'median /us':>11}  {'CV':>6}  {'max/min':>8}")
        out.append("  " + "-" * 96)
        for r in sweep_results:
            if r["path"] is None:
                out.append(f"  {r['label']:<24}  (no valid path)")
                continue
            p_str = "-".join(f"Q{q}" for q in r["path"])
            # Best-case: smallest max/min
            best = None
            for k, s in r["interp_summary"].items():
                if best is None or s["max_over_min"] < best[1]["max_over_min"]:
                    best = (k, s)
            if best is None:
                out.append(f"  {r['label']:<24}  {p_str:<30}  (no valid interpretation)")
                continue
            k, s = best
            out.append(f"  {r['label']:<24}  {p_str:<30}  {k:<28}  "
                       f"{s['median']:>11.6f}  {s['cv']:>6.1%}  {s['max_over_min']:>8.2f}x")
        out.append("")
        # Summary of best config
        finite = [r for r in sweep_results if r.get("interp_summary")]
        if finite:
            best_cfg = min(
                finite,
                key=lambda r: min((v["max_over_min"] for v in r["interp_summary"].values()), default=float("inf"))
            )
            best_k = min(best_cfg["interp_summary"], key=lambda k: best_cfg["interp_summary"][k]["max_over_min"])
            best_s = best_cfg["interp_summary"][best_k]
            out.append(f"  Best configuration overall: {best_cfg['label']}")
            out.append(f"    Path: {best_cfg['path']}")
            out.append(f"    Best interpretation: {best_k}")
            out.append(f"    gamma_0 median: {best_s['median']:.6f} /us,  CV: {best_s['cv']:.1%},  max/min: {best_s['max_over_min']:.2f}x")
            out.append("")

    out.append("7) Alternative reading: gamma_0 as a floor, not a multiplier")
    out.append("-" * 78)
    out.append("  EQ-017 reading (b): gamma_phi(i) = gamma_0 + device_noise(i)")
    out.append("  with gamma_0 constant and non-negative TLS/fabrication noise on top.")
    out.append("  Under this reading:")
    non_outliers = [s for s in qubit_stats.values() if s["gamma_phi_cv"] <= 0.5]
    g_non = np.array([s["gamma_phi_median"] for s in non_outliers])
    out.append(f"    Non-outlier qubits (CV <= 50%):   {len(non_outliers)}")
    out.append(f"    min gamma_phi (floor candidate):  {g_non.min():.6f} /us")
    out.append(f"    5th percentile:                   {np.percentile(g_non, 5):.6f} /us")
    out.append(f"    median gamma_phi:                 {np.median(g_non):.6f} /us")
    out.append(f"    Implied gamma_0 (floor):          {g_non.min():.6f} /us")
    out.append(f"    Implied T2_framework (2/gamma_0): {2.0/g_non.min():.1f} us")
    out.append("")
    out.append("8) Verdict")
    out.append("-" * 78)

    best_cv = None
    best_key = None
    for key, _, _ in interpretations:
        name = f"gamma_0_{key}"
        vals = np.array([r[name] for r in per_qubit], dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue
        cv = vals.std() / abs(vals.mean()) if vals.mean() != 0 else float("inf")
        if best_cv is None or cv < best_cv:
            best_cv = cv
            best_key = key
    out.append(f"  Best-case CV across interpretations: {best_cv:.1%} (interp={best_key})")
    max_ratio = None
    for key, _, _ in interpretations:
        name = f"gamma_0_{key}"
        vals = np.array([r[name] for r in per_qubit], dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0 or vals.min() <= 0:
            continue
        ratio = vals.max() / vals.min()
        if max_ratio is None or ratio < max_ratio:
            max_ratio = ratio
            best_ratio_key = key
    out.append(f"  Best-case max/min ratio: {max_ratio:.2f}x (interp={best_ratio_key})")
    out.append("")
    out.append("  Success criterion (task):   max/min < 3x across non-outlier qubits")
    out.append("  Failure criterion (task):   max/min > 10x after outlier removal")
    out.append(f"  Sub-chain observed max/min: {max_ratio:.2f}x")
    out.append("")
    if max_ratio is not None:
        if max_ratio < 3.0:
            verdict = "CONSISTENT: mode-exposure reading (a) supported on this sub-chain."
        elif max_ratio > 10.0:
            verdict = "INCONSISTENT: mode-exposure reading (a) does NOT explain hardware gamma_phi."
        else:
            verdict = "INCONCLUSIVE: falls between success and failure thresholds."
    else:
        verdict = "INCONCLUSIVE: no valid interpretation gave finite values."
    out.append(f"  VERDICT: {verdict}")
    out.append("")
    out.append("  Caveat: this test uses idle-Ramsey gamma_phi against gate-active J.")
    out.append("  Idle dynamics has H_idle ~ diagonal (no chain mode structure), so the")
    out.append("  per-qubit T2 directly measures local dephasing, NOT chain-mode dephasing.")
    out.append("  The mode-exposure formula applies to a chain running its coupling")
    out.append("  Hamiltonian (e.g. during a CPsi crossing experiment). A proper test")
    out.append("  requires Phase 2: run chain dynamics on the selected sub-chain and")
    out.append("  compare observed decoherence to gamma_0 * |a_B|^2.")
    out.append("")
    out.append("  EQ-017 reading (b) is easier to test from idle data alone:")
    out.append("  the minimum gamma_phi across non-outlier qubits gives a gamma_0 floor")
    out.append(f"  candidate of {g_non.min():.6f} /us, corresponding to T2_framework = {2.0/g_non.min():.1f} us.")
    out.append("")
    out.append("=" * 78)
    out.append("End of report")
    out.append("=" * 78)

    return "\n".join(out)


# --- Main ----------------------------------------------------------------

def main():
    print("=" * 72)
    print("Phase 1: Kingston gamma_0 consistency test")
    print("=" * 72)

    print("\n[1/6] Parsing 8 calibration CSVs...")
    qubit_stats, edge_stats, n_snaps = load_calibrations()
    g_all = np.array([s["gamma_phi_median"] for s in qubit_stats.values()])
    cv_all = np.array([s["gamma_phi_cv"] for s in qubit_stats.values()])
    print(f"  n_snapshots={n_snaps}, qubits={len(qubit_stats)}, edges={len(edge_stats)}")
    print(f"  gamma_phi /us: range [{g_all.min():.6f}, {g_all.max():.6f}], "
          f"ratio {g_all.max()/g_all.min():.1f}x, median {np.median(g_all):.6f}")
    print(f"  CV median {np.median(cv_all):.1%}, TLS outliers (CV>50%) = "
          f"{int((cv_all > 0.5).sum())}")

    print("\n[2/6] Query ibm_kingston backend (cached if available)...")
    backend_data = query_backend(use_cache=True)
    if backend_data is None:
        print("  Proceeding with CSV-only data.")

    print("\n[3/6] Extract J couplings...")
    J_est, J_src = extract_J_per_edge(backend_data, edge_stats)
    J_vals = np.array(list(J_est.values()))
    print(f"  {len(J_est)} edges with J estimate.")
    if len(J_vals) > 0:
        print(f"  J [MHz] range [{J_vals.min():.3f}, {J_vals.max():.3f}], "
              f"median {np.median(J_vals):.3f}")
    src_counts = defaultdict(int)
    for s in J_src.values():
        src_counts[s] += 1
    for src, c in sorted(src_counts.items()):
        print(f"  source={src}: {c} edges")

    print("\n[4/6] Select linear sub-chain (multi-configuration sweep)...")
    undirected = sorted(set(edge_stats.keys()))
    if backend_data is not None:
        undirected = sorted(set(tuple(sorted(e)) for e in backend_data["coupling_map_directed"]))

    configs = [
        dict(length=5, anchors=CPSI_CROSSING_QUBITS, max_cv=None, label="len5_anchor_any_CV"),
        dict(length=7, anchors=CPSI_CROSSING_QUBITS, max_cv=None, label="len7_anchor_any_CV"),
        dict(length=5, anchors=None,                 max_cv=0.30, label="len5_free_CV30"),
        dict(length=7, anchors=None,                 max_cv=0.30, label="len7_free_CV30"),
    ]

    sweep_results = []
    for cfg in configs:
        path, score, median_T2 = find_linear_path(
            undirected, qubit_stats,
            length=cfg["length"], anchors=cfg["anchors"], max_cv=cfg["max_cv"],
        )
        if path is None:
            print(f"  [{cfg['label']}] no path found")
            sweep_results.append({**cfg, "path": None})
            continue
        H_ext, evs_ext, evecs_ext, _ = build_H(path, J_est, uniform=False)
        H_uni, evs_uni, evecs_uni, _ = build_H(path, J_est, uniform=True)
        per_qubit = gamma0_test(path, evs_ext, evecs_ext, qubit_stats)

        interp_summary = {}
        for key in ("gamma_0_A_max_mode", "gamma_0_B_slow_mode",
                    "gamma_0_C_total_participation", "gamma_0_D_IPR_weighted"):
            vals = np.array([r[key] for r in per_qubit], dtype=float)
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0 or vals.min() <= 0:
                continue
            cv = float(vals.std() / abs(vals.mean()))
            interp_summary[key] = {
                "median": float(np.median(vals)),
                "cv": cv,
                "max_over_min": float(vals.max() / vals.min()),
            }

        sweep_results.append({
            **cfg,
            "path": path, "score": score, "median_T2_all": median_T2,
            "eigvals_ext": evs_ext.tolist(), "eigvals_uni": evs_uni.tolist(),
            "evecs_ext": evecs_ext.tolist(),
            "per_qubit": per_qubit,
            "interp_summary": interp_summary,
        })
        print(f"  [{cfg['label']}] path={path}")
        for key, s in interp_summary.items():
            print(f"     {key}: median={s['median']:.6f}  CV={s['cv']:.1%}  max/min={s['max_over_min']:.2f}x")

    # Use the length-5-anchor config as the primary path for the detailed report.
    primary = sweep_results[0]
    if primary["path"] is None:
        # Fall back to first config that found a path
        primary = next((r for r in sweep_results if r["path"] is not None), None)
    if primary is None:
        print("  ERROR: no valid linear path found for any configuration. Aborting.")
        return

    path = primary["path"]
    result = {
        "path": path,
        "eigvals": primary["eigvals_ext"],
        "eigvecs": primary["evecs_ext"],
        "per_qubit": primary["per_qubit"],
    }
    uniform_result = {
        "path": path,
        "eigvals": primary["eigvals_uni"],
        "eigvecs": None,
    }

    print("\n[5/6] Writing report...")
    report = write_report(
        report_lines=None,
        path_info={"score": primary["score"], "median_T2": primary["median_T2_all"]},
        qubit_stats=qubit_stats,
        edge_stats=edge_stats,
        backend_data=backend_data,
        J_est=J_est,
        J_src=J_src,
        path=path,
        result=result,
        uniform_result=uniform_result,
        undirected_edges=undirected,
        sweep_results=sweep_results,
    )
    REPORT_TXT.write_text(report, encoding="utf-8")
    print(f"  Wrote {REPORT_TXT}")


if __name__ == "__main__":
    main()
