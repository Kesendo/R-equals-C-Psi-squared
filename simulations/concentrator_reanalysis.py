"""Re-analysis of the March 2026 ibm_torino concentrator run (zero QPU).

Data: data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_{config}_*.json
(raw_counts: 45 five-qubit counts dicts per config = 5 time points x 9
tomography settings, time-major; settings in (a,b) row-major order over
BASES = [X,Y,Z], even qubits measured in a, odd in b).

The estimator is re-implemented VERBATIM from the original runner
(external pipeline, run_sacrifice_zone.py): linear-inversion pair
reconstruction, eigenvalue clip at 1e-15 + renormalize, entropies in
bits, MI = max(0, S_A + S_B - S_AB).

GATE-FIRST discipline:
  Stage 0 (the handshake): reproduce the PUBLISHED pair_mi / sum_mi from
          raw_counts. Any mismatch beyond float tolerance = abort.
  Stage 1: bootstrap error bars (multinomial resampling of the full
          5-qubit counts per setting, preserving cross-pair correlation).
  Stage 2: null bias floor per config/time: product-of-marginals per
          pair per setting (kills correlations, keeps singles), same
          estimator -> the MI the pipeline reports on a ZERO-true-MI
          state with these exact marginals and shot numbers.
  Stage 3: floor-aware verdicts: ordering significance (bootstrap),
          uniform-row-vs-floor, floor-corrected differences and ratios.

Stage 0 passed 2026-07-11 with worst deviation 0.0 over all 15 cells.
Results are folded into experiments/IBM_CONCENTRATOR.md (the 2026-07-11
re-analysis subsection). Deterministic under SEED; verdicts checked
seed-stable at SEED = 20260711 and 424242. The JSON summary and run logs
this script writes (leading-underscore siblings) are throwaway sidecars,
deliberately left uncommitted (note: the gitignore WIP rule covers only
simulations/_*.py, so these are untracked, not ignored); the printed
output is the record.

Two statistics notes (2026-07-11 empty review, statistician lens):
  - The ordering confidence in Stage 3a is reported BOTH raw and
    floor-adjusted. The floor is config-dependent (selective preserves
    more single-qubit coherence, hence carries the larger positive MI
    bias), so the raw comparison is anti-conservative; the adjusted one
    shifts the comparison by the floor-mean difference. The floor itself
    is a heuristic bias proxy (it conditions on the measured marginals
    only), so both numbers are reported.
  - P here is a bootstrap ordering confidence (the fraction of shot-noise
    resamples preserving the observed ordering), not a frequentist
    p-value against an equal-MI null; no such null is available from two
    distinct-state configurations.
"""

import json
import sys
from pathlib import Path

import numpy as np

# ================================================================
# Constants (verbatim from the original runner)
# ================================================================
N_QUBITS = 5
BASES = ['X', 'Y', 'Z']
TROTTER_STEPS = [2, 4, 6, 8, 10]
DT = 0.5
TIME_POINTS = [n * DT for n in TROTTER_STEPS]

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {'I': I2, 'X': SX, 'Y': SY, 'Z': SZ}

REPO = Path(__file__).resolve().parents[1]
DATA_DIR = REPO / "data" / "ibm_sacrifice_zone_march2026"
CONFIG_FILES = {
    "selective_dd": "sacrifice_zone_hw_selective_dd_20260324_191523.json",
    "uniform_dd": "sacrifice_zone_hw_uniform_dd_20260324_191614.json",
    "no_dd": "sacrifice_zone_hw_no_dd_20260324_191713.json",
}
PAIRS = [(i, i + 1) for i in range(N_QUBITS - 1)]

N_BOOT = 1000
N_NULL = 1000
SEED = 20260711


def build_basis_list():
    """Setting s -> per-qubit bases. Even qubit -> a, odd -> b; (a,b) row-major."""
    basis_list = []
    for a in BASES:
        for b in BASES:
            basis_list.append([a if q % 2 == 0 else b for q in range(N_QUBITS)])
    return basis_list


BASIS_LIST = build_basis_list()


# ================================================================
# Estimator, verbatim re-implementation
# ================================================================

def marginal_counts(counts, qi, qj):
    """2-qubit marginal; Qiskit bitstring reversed so index i = qubit i."""
    marg = {'00': 0, '01': 0, '10': 0, '11': 0}
    for bitstr, cnt in counts.items():
        bits = bitstr.replace(' ', '')[::-1]
        marg[bits[qi] + bits[qj]] += cnt
    return marg


def reconstruct_rho_from_marginals(margs, pair_bases):
    """margs: list of 9 marginal dicts; pair_bases: list of 9 (bi, bj)."""
    exps = {('I', 'I'): 1.0}
    single_i = {'X': [], 'Y': [], 'Z': []}
    single_j = {'X': [], 'Y': [], 'Z': []}
    for m, (bi, bj) in zip(margs, pair_bases):
        total = sum(m.values())
        if total == 0:
            continue
        exps[(bi, bj)] = (m['00'] - m['01'] - m['10'] + m['11']) / total
        single_i[bi].append((m['00'] + m['01'] - m['10'] - m['11']) / total)
        single_j[bj].append((m['00'] + m['10'] - m['01'] - m['11']) / total)
    for p in BASES:
        if single_i[p]:
            exps[(p, 'I')] = float(np.mean(single_i[p]))
        if single_j[p]:
            exps[('I', p)] = float(np.mean(single_j[p]))
    rho = np.zeros((4, 4), dtype=complex)
    for a in ['I', 'X', 'Y', 'Z']:
        for b in ['I', 'X', 'Y', 'Z']:
            rho += exps.get((a, b), 0.0) * np.kron(PAULIS[a], PAULIS[b])
    return rho / 4.0


def vn_entropy(rho):
    eigs = np.linalg.eigvalsh(rho)
    eigs = np.clip(eigs.real, 1e-15, None)
    eigs = eigs / eigs.sum()
    return float(-np.sum(eigs * np.log2(eigs)))


def mi(rho_ab):
    rho_a = np.array([
        [rho_ab[0, 0] + rho_ab[1, 1], rho_ab[0, 2] + rho_ab[1, 3]],
        [rho_ab[2, 0] + rho_ab[3, 1], rho_ab[2, 2] + rho_ab[3, 3]],
    ])
    rho_b = np.array([
        [rho_ab[0, 0] + rho_ab[2, 2], rho_ab[0, 1] + rho_ab[2, 3]],
        [rho_ab[1, 0] + rho_ab[3, 2], rho_ab[1, 1] + rho_ab[3, 3]],
    ])
    return max(0.0, vn_entropy(rho_a) + vn_entropy(rho_b) - vn_entropy(rho_ab))


def sum_mi_from_counts(tp_counts):
    """tp_counts: the 9 five-qubit counts dicts of one time point."""
    pair_mi = {}
    for qi, qj in PAIRS:
        margs = [marginal_counts(c, qi, qj) for c in tp_counts]
        pair_bases = [(BASIS_LIST[s][qi], BASIS_LIST[s][qj]) for s in range(9)]
        rho = reconstruct_rho_from_marginals(margs, pair_bases)
        pair_mi[f"({qi},{qj})"] = mi(rho)
    return sum(pair_mi.values()), pair_mi


# ================================================================
# Stage 0: the handshake gate
# ================================================================

def load_config(name):
    with open(DATA_DIR / CONFIG_FILES[name], 'r') as f:
        return json.load(f)


def stage0_gate(data_by_config):
    print("=" * 72)
    print("STAGE 0 - HANDSHAKE GATE: reproduce published pair_mi from raw_counts")
    print("=" * 72)
    worst = 0.0
    for cfg, data in data_by_config.items():
        raw = data["raw_counts"]
        assert len(raw) == 45, f"{cfg}: expected 45 counts dicts, got {len(raw)}"
        for tp_idx, t in enumerate(TIME_POINTS):
            tp_counts = raw[tp_idx * 9:(tp_idx + 1) * 9]
            s, pm = sum_mi_from_counts(tp_counts)
            s_pub = data["results"]["sum_mi"][tp_idx]
            pm_pub = data["results"]["pair_mi"][tp_idx]
            ds = abs(s - s_pub)
            dp = max(abs(pm[k] - pm_pub[k]) for k in pm)
            worst = max(worst, ds, dp)
            status = "OK" if max(ds, dp) < 1e-9 else "MISMATCH"
            print(f"  {cfg:>13} t={t:.1f}: sum_mi recomputed {s:.12f}"
                  f"  published {s_pub:.12f}  |d|={ds:.2e}  [{status}]")
    print(f"\n  worst deviation (sum_mi and per-pair): {worst:.3e}")
    if worst >= 1e-9:
        print("\n  GATE FAILED - handshake break. ABORTING before any statistics.")
        sys.exit(1)
    print("  GATE PASSED - estimator and data ordering verified.\n")


# ================================================================
# Stage 1: bootstrap error bars
# ================================================================

def resample_counts(counts, rng):
    """Multinomial resample of one 5-qubit counts dict (same total)."""
    keys = list(counts.keys())
    vals = np.array([counts[k] for k in keys], dtype=float)
    tot = int(vals.sum())
    draw = rng.multinomial(tot, vals / vals.sum())
    return {k: int(n) for k, n in zip(keys, draw) if n > 0}


def stage1_bootstrap(data_by_config, rng):
    print("=" * 72)
    print(f"STAGE 1 - BOOTSTRAP ({N_BOOT} replicates, full-counts multinomial)")
    print("=" * 72)
    boots = {}  # cfg -> array (N_BOOT, 5)
    for cfg, data in data_by_config.items():
        raw = data["raw_counts"]
        arr = np.zeros((N_BOOT, len(TIME_POINTS)))
        for b in range(N_BOOT):
            for tp_idx in range(len(TIME_POINTS)):
                tp_counts = [resample_counts(c, rng)
                             for c in raw[tp_idx * 9:(tp_idx + 1) * 9]]
                arr[b, tp_idx], _ = sum_mi_from_counts(tp_counts)
        boots[cfg] = arr
        pub = np.array(data["results"]["sum_mi"])
        print(f"\n  {cfg}: published / boot mean +- sd / [2.5%, 97.5%]")
        for i, t in enumerate(TIME_POINTS):
            lo, hi = np.percentile(arr[:, i], [2.5, 97.5])
            print(f"    t={t:.1f}: {pub[i]:.4f} / {arr[:, i].mean():.4f} "
                  f"+- {arr[:, i].std(ddof=1):.4f} / [{lo:.4f}, {hi:.4f}]")
    return boots


# ================================================================
# Stage 2: null bias floor (product-of-marginals per pair per setting)
# ================================================================

def null_floor(data, rng, n_null=N_NULL):
    """Per time point: Sum-MI distribution the estimator reports when the
    TRUE pair state is the product of the measured per-setting singles."""
    raw = data["raw_counts"]
    floors = np.zeros((n_null, len(TIME_POINTS)))
    # Precompute per (tp, pair, setting): total shots + product distribution
    prod_p = {}
    for tp_idx in range(len(TIME_POINTS)):
        tp_counts = raw[tp_idx * 9:(tp_idx + 1) * 9]
        for (qi, qj) in PAIRS:
            for s, c in enumerate(tp_counts):
                m = marginal_counts(c, qi, qj)
                tot = sum(m.values())
                pi1 = (m['10'] + m['11']) / tot  # P(bit_i = 1)
                pj1 = (m['01'] + m['11']) / tot  # P(bit_j = 1)
                p = np.array([(1 - pi1) * (1 - pj1), (1 - pi1) * pj1,
                              pi1 * (1 - pj1), pi1 * pj1])
                prod_p[(tp_idx, qi, qj, s)] = (tot, p)
    outcomes = ['00', '01', '10', '11']
    for n in range(n_null):
        for tp_idx in range(len(TIME_POINTS)):
            total_mi = 0.0
            for (qi, qj) in PAIRS:
                margs = []
                for s in range(9):
                    tot, p = prod_p[(tp_idx, qi, qj, s)]
                    draw = rng.multinomial(tot, p)
                    margs.append(dict(zip(outcomes, draw.tolist())))
                pair_bases = [(BASIS_LIST[s][qi], BASIS_LIST[s][qj])
                              for s in range(9)]
                rho = reconstruct_rho_from_marginals(margs, pair_bases)
                total_mi += mi(rho)
            floors[n, tp_idx] = total_mi
    return floors


def stage2_floor(data_by_config, rng):
    print("\n" + "=" * 72)
    print(f"STAGE 2 - NULL BIAS FLOOR ({N_NULL} draws, zero-true-MI product nulls)")
    print("=" * 72)
    floors = {}
    for cfg, data in data_by_config.items():
        fl = null_floor(data, rng)
        floors[cfg] = fl
        pub = np.array(data["results"]["sum_mi"])
        print(f"\n  {cfg}: floor mean +- sd / [2.5%, 97.5%] / published / excess")
        for i, t in enumerate(TIME_POINTS):
            lo, hi = np.percentile(fl[:, i], [2.5, 97.5])
            print(f"    t={t:.1f}: {fl[:, i].mean():.4f} +- {fl[:, i].std(ddof=1):.4f}"
                  f" / [{lo:.4f}, {hi:.4f}] / {pub[i]:.4f}"
                  f" / {pub[i] - fl[:, i].mean():+.4f}")
    return floors


# ================================================================
# Stage 3: floor-aware verdicts
# ================================================================

def stage3_verdicts(data_by_config, boots, floors):
    print("\n" + "=" * 72)
    print("STAGE 3 - VERDICTS")
    print("=" * 72)
    sel = boots["selective_dd"]
    uni = boots["uniform_dd"]
    nod = boots["no_dd"]
    pub = {c: np.array(d["results"]["sum_mi"]) for c, d in data_by_config.items()}

    print("\n  (a) Ordering selective > uniform, bootstrap ordering confidence.")
    print("      raw = P(sel_b > uni_b); floor-adjusted shifts by the")
    print("      config-dependent floor-mean difference (see docstring).")
    fdiff = floors["selective_dd"].mean(axis=0) - floors["uniform_dd"].mean(axis=0)
    for i, t in enumerate(TIME_POINTS):
        p_raw = float(np.mean(sel[:, i] > uni[:, i]))
        p_adj = float(np.mean(sel[:, i] - uni[:, i] > fdiff[i]))
        d = sel[:, i] - uni[:, i]
        d_hat = pub['selective_dd'][i] - pub['uniform_dd'][i]
        lo, hi = np.percentile(d, [2.5, 97.5])
        piv_lo, piv_hi = 2 * d_hat - hi, 2 * d_hat - lo
        print(f"    t={t:.1f}: P_raw = {p_raw:.3f}  P_adj = {p_adj:.3f};"
              f"  diff {d_hat:+.4f}"
              f"  percentile CI [{lo:+.4f}, {hi:+.4f}]"
              f"  pivotal CI [{piv_lo:+.4f}, {piv_hi:+.4f}]")
    joint_raw = float(np.mean(np.all(sel > uni, axis=1)))
    joint_adj = float(np.mean(np.all(sel - uni > fdiff[None, :], axis=1)))
    print(f"    joint (all 5 simultaneously): P_raw = {joint_raw:.3f}"
          f"  P_adj = {joint_adj:.3f}")

    print("\n  (b) Is the uniform row consistent with the ZERO-true-MI floor?")
    fl_u = floors["uniform_dd"]
    for i, t in enumerate(TIME_POINTS):
        p_ge = float(np.mean(fl_u[:, i] >= pub["uniform_dd"][i]))
        print(f"    t={t:.1f}: uniform {pub['uniform_dd'][i]:.4f}  floor "
              f"{fl_u[:, i].mean():.4f} +- {fl_u[:, i].std(ddof=1):.4f}"
              f"  P(floor >= measured) = {p_ge:.3f}")

    print("\n  (c) Floor-corrected Sum-MI (published - own floor mean), per config.")
    print("      Parenthesis = shot-noise SE of the measured value (boot sd);")
    print("      the floor MEAN's own SE is floor_sd/sqrt(N_NULL), negligible.")
    for cfg in ["selective_dd", "uniform_dd", "no_dd"]:
        exc = pub[cfg] - floors[cfg].mean(axis=0)
        sd = boots[cfg].std(axis=0, ddof=1)
        line = "  ".join(f"{e:+.4f}({s:.4f})" for e, s in zip(exc, sd))
        print(f"    {cfg:>13}: {line}")

    print("\n  (d) Ratio sel/uni raw vs floor-corrected (point estimates):")
    for i, t in enumerate(TIME_POINTS):
        raw_ratio = pub["selective_dd"][i] / pub["uniform_dd"][i]
        num = pub["selective_dd"][i] - floors["selective_dd"][:, i].mean()
        den = pub["uniform_dd"][i] - floors["uniform_dd"][:, i].mean()
        corr = num / den if den > 0 else float('nan')
        print(f"    t={t:.1f}: raw {raw_ratio:.2f}  floor-corrected "
              f"{corr:.2f}" + ("  (denominator <= 0: ratio undefined)" if den <= 0 else ""))

    print("\n  (e) Selective vs no_dd (the quieter comparison):")
    for i, t in enumerate(TIME_POINTS):
        p_gt = float(np.mean(sel[:, i] > nod[:, i]))
        print(f"    t={t:.1f}: sel {pub['selective_dd'][i]:.4f} vs no_dd "
              f"{pub['no_dd'][i]:.4f}  P(sel_b > no_dd_b) = {p_gt:.3f}")


def main():
    rng = np.random.default_rng(SEED)
    data_by_config = {c: load_config(c) for c in CONFIG_FILES}
    stage0_gate(data_by_config)
    boots = stage1_bootstrap(data_by_config, rng)
    floors = stage2_floor(data_by_config, rng)
    stage3_verdicts(data_by_config, boots, floors)
    out = {
        "seed": SEED, "n_boot": N_BOOT, "n_null": N_NULL,
        "boot_sd": {c: boots[c].std(axis=0, ddof=1).tolist() for c in boots},
        "floor_mean": {c: floors[c].mean(axis=0).tolist() for c in floors},
        "floor_sd": {c: floors[c].std(axis=0, ddof=1).tolist() for c in floors},
    }
    out_path = Path(__file__).with_name("_concentrator_reanalysis_out.json")
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\n  summary written: {out_path}")


if __name__ == "__main__":
    main()
