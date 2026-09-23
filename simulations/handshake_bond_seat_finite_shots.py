"""Finite-shot Z-position reading of the fixed-gamma N=7 bond response.

This is a local simulated protocol, not a hardware measurement. A shot in the
single-excitation sector has exactly one of seven site-position outcomes.
"""

from __future__ import annotations

import sys

import numpy as np
import scipy
from scipy.linalg import expm
from scipy.special import xlogy

from framework.sine_slater import sine_mode_matrix
from handshake_bond_seat_fixed_gamma import (
    bond_superoperator, liouvillian, response_tensor,
)


N = 7
GAMMA0 = 0.05
J_HOP = 0.075
SHOTS_PER_TIME = 100000
REPEATS = 256
BASELINE_TIMES = np.array([2., 4., 6., 8., 10., 12., 15., 18., 24., 30.])
PROFILES = {"uniform": tuple(range(N)), "single_centre": (N // 2,)}


def finite_probabilities(times, watched_sites: tuple[int, ...],
                         delta) -> np.ndarray:
    """Return [time, site] probabilities with nominal psi1 kept as preparation."""
    grid = np.asarray(times, dtype=float)
    changes = np.asarray(delta, dtype=float)
    if grid.ndim != 1 or not np.all(np.isfinite(grid)) or np.any(grid < 0):
        raise ValueError("times must be finite and nonnegative")
    if changes.shape != (N - 1,) or not np.all(np.isfinite(changes)):
        raise ValueError("delta must contain six finite additive hoppings")
    if np.any(J_HOP + changes <= 0):
        raise ValueError("all defected hoppings must remain positive")

    generator = liouvillian(N, J_HOP, GAMMA0, watched_sites)
    for bond, change in enumerate(changes):
        if change:
            generator = generator + change * bond_superoperator(N, bond)
    psi = sine_mode_matrix(N)[0]
    initial = np.outer(psi, psi).ravel()
    probabilities = np.empty((len(grid), N), dtype=float)
    for index, time in enumerate(grid):
        evolved = (expm(time * generator) @ initial).reshape(N, N)
        diagonal = np.diag(evolved)
        if np.max(np.abs(diagonal.imag)) > 1e-10:
            raise ArithmeticError("population trace lost Hermiticity")
        probabilities[index] = diagonal.real
    _validate_probabilities(probabilities)
    return probabilities


def fisher_matrix(probabilities: np.ndarray, slopes: np.ndarray,
                  shots_per_time: int, site: int | None = None) -> np.ndarray:
    """Shot-weighted local bond Fisher matrix for categorical or binary Z data.

    `slopes` has the existing response_tensor order [time, bond, site].
    """
    p = np.asarray(probabilities, dtype=float)
    s = np.asarray(slopes, dtype=float)
    if p.ndim != 2 or s.ndim != 3 or s.shape != (len(p), s.shape[1], p.shape[1]):
        raise ValueError("probabilities and slopes have incompatible shapes")
    if (not isinstance(shots_per_time, (int, np.integer)) or
            shots_per_time <= 0):
        raise ValueError("shots_per_time must be a positive integer")
    if not np.all(np.isfinite(p)) or not np.all(np.isfinite(s)):
        raise ValueError("probabilities and slopes must be finite")
    _validate_probabilities(p)
    if site is None:
        if np.any(p <= 0):
            raise ValueError("categorical Fisher requires positive probabilities")
        return shots_per_time * np.einsum("tbj,tcj,tj->bc", s, s, 1.0 / p)
    if not 0 <= site < p.shape[1] or np.any(p[:, site] <= 0) or np.any(p[:, site] >= 1):
        raise ValueError("binary Fisher requires an interior site probability")
    return shots_per_time * np.einsum(
        "tb,tc,t->bc", s[:, :, site], s[:, :, site],
        1.0 / (p[:, site] * (1.0 - p[:, site])),
    )


def select_schedule(per_time_fisher: np.ndarray,
                    baseline_indices: tuple[int, ...]) -> tuple[int, ...]:
    """Local E-optimal one-for-one time exchanges from a fixed baseline.

    Candidates are visited in index order. Changes within a roundoff budget
    count as ties and keep the first best schedule encountered.
    """
    matrices = np.asarray(per_time_fisher, dtype=float)
    if (matrices.ndim != 3 or matrices.shape[1] != matrices.shape[2] or
            not np.all(np.isfinite(matrices))):
        raise ValueError("per-time Fisher matrices must be finite square matrices")
    chosen = tuple(sorted(baseline_indices))
    if (not chosen or len(set(chosen)) != len(chosen) or
            any(index < 0 or index >= len(matrices) for index in chosen)):
        raise ValueError("baseline indices must be distinct candidate positions")

    while True:
        current_matrix = np.sum(matrices[list(chosen)], axis=0)
        current_min = float(np.linalg.eigvalsh(current_matrix)[0])
        budget = (128.0 * np.finfo(float).eps *
                  max(1.0, float(np.linalg.norm(current_matrix, ord=2))))
        best_indices = None
        best_min = current_min
        occupied = set(chosen)
        for removed in chosen:
            for inserted in range(len(matrices)):
                if inserted in occupied:
                    continue
                proposed_matrix = current_matrix - matrices[removed] + matrices[inserted]
                candidate_min = float(np.linalg.eigvalsh(proposed_matrix)[0])
                if candidate_min > best_min + budget:
                    best_indices = tuple(sorted((occupied - {removed}) | {inserted}))
                    best_min = candidate_min
        if best_indices is None:
            return chosen
        chosen = best_indices


def hypotheses() -> list[tuple[str, np.ndarray]]:
    """Fixed null and signed single-bond alternatives, in tie-break order."""
    cases = [("null", np.zeros(N - 1))]
    for bond in range(N - 1):
        for sign, label in ((-1, "-"), (+1, "+")):
            change = np.zeros(N - 1)
            change[bond] = sign * 0.1 * J_HOP
            cases.append((f"b{bond}{label}", change))
    return cases


def _validate_probabilities(probabilities: np.ndarray) -> None:
    if (not np.all(np.isfinite(probabilities)) or
            np.any(probabilities < 0) or np.any(probabilities > 1) or
            not np.allclose(probabilities.sum(axis=-1), 1.0, rtol=0, atol=1e-10)):
        raise ValueError("probability rows must be finite, nonnegative and sum to one")


def log_likelihoods(counts: np.ndarray, models: np.ndarray,
                    site: int | None = None) -> np.ndarray:
    """Exact categorical or site-Bernoulli scores; [repeat, hypothesis].

    Terms with zero count contribute zero even at a model-zero probability.
    A positive count at a model-zero probability contributes negative infinity.
    """
    observed = np.asarray(counts)
    predicted = np.asarray(models, dtype=float)
    if (observed.ndim != 3 or predicted.ndim != 3 or
            observed.shape[1:] != predicted.shape[1:] or
            not len(predicted) or not np.all(np.isfinite(observed)) or
            np.any(observed < 0)):
        raise ValueError("counts and models must have compatible nonnegative shapes")
    _validate_probabilities(predicted)
    if site is None:
        return np.sum(xlogy(observed[:, None, :, :], predicted[None, :, :, :]),
                      axis=(2, 3))
    if not 0 <= site < observed.shape[-1]:
        raise ValueError("site outside outcome range")
    yes = observed[:, :, site]
    no = observed.sum(axis=2) - yes
    p = predicted[:, :, site]
    return np.sum(xlogy(yes[:, None, :], p[None, :, :]) +
                  xlogy(no[:, None, :], (1.0 - p)[None, :, :]), axis=2)


def classify(counts: np.ndarray, models: np.ndarray,
             site: int | None = None) -> np.ndarray:
    """Maximum known-model likelihood, earliest hypothesis on exact ties."""
    return np.argmax(log_likelihoods(counts, models, site=site), axis=1)


def draw_counts(probabilities: np.ndarray, shots_per_time: int, repeats: int,
                seed_parts: tuple[int, ...]) -> np.ndarray:
    """Independent categorical shots by time, reproducible per hypothesis."""
    p = np.asarray(probabilities, dtype=float)
    if p.ndim != 2 or not len(p) or not p.shape[1]:
        raise ValueError("probabilities must have [time, outcome] shape")
    _validate_probabilities(p)
    if (not isinstance(shots_per_time, (int, np.integer)) or
            shots_per_time <= 0 or not isinstance(repeats, (int, np.integer)) or
            repeats <= 0):
        raise ValueError("shots and repeats must be positive integers")
    rng = np.random.default_rng(np.random.SeedSequence(seed_parts))
    return np.stack([rng.multinomial(shots_per_time, row, size=repeats)
                     for row in p], axis=1)


def wilson_interval(successes: int, trials: int) -> tuple[float, float]:
    """Two-sided 95% Wilson score interval for one fixed hypothesis stratum."""
    if (not isinstance(successes, (int, np.integer)) or
            not isinstance(trials, (int, np.integer)) or
            trials <= 0 or not 0 <= successes <= trials):
        raise ValueError("require 0 <= integer successes <= positive trials")
    z = 1.959963984540054
    fraction = successes / trials
    scale = 1.0 + z * z / trials
    centre = (fraction + z * z / (2 * trials)) / scale
    radius = z * np.sqrt(fraction * (1 - fraction) / trials +
                         z * z / (4 * trials * trials)) / scale
    return max(0.0, float(centre - radius)), min(1.0, float(centre + radius))


def evaluate_counts(models: np.ndarray, shots_per_time: int, repeats: int,
                    seed_prefix: tuple[int, ...]) -> dict:
    """Two readings of each synthetic record, with paired correctness counts."""
    predicted = np.asarray(models, dtype=float)
    if predicted.ndim != 3 or not len(predicted) or predicted.shape[2] < 3:
        raise ValueError("models need [hypothesis,time,outcome>=3] shape")
    _validate_probabilities(predicted)
    count = len(predicted)
    full = np.zeros((count, count), dtype=int)
    site = np.zeros((count, count), dtype=int)
    paired = {key: np.zeros(count, dtype=int) for key in
              ("paired_both", "paired_full_only", "paired_site_only",
               "paired_neither")}
    for truth in range(count):
        records = draw_counts(predicted[truth], shots_per_time, repeats,
                              seed_parts=(*seed_prefix, truth))
        if not np.all(records.sum(axis=2) == shots_per_time):
            raise ArithmeticError("multinomial shot conservation failed")
        full_guess = classify(records, predicted)
        site_guess = classify(records, predicted, site=2)
        full[truth] = np.bincount(full_guess, minlength=count)
        site[truth] = np.bincount(site_guess, minlength=count)
        full_correct = full_guess == truth
        site_correct = site_guess == truth
        paired["paired_both"][truth] = np.count_nonzero(full_correct & site_correct)
        paired["paired_full_only"][truth] = np.count_nonzero(full_correct & ~site_correct)
        paired["paired_site_only"][truth] = np.count_nonzero(~full_correct & site_correct)
        paired["paired_neither"][truth] = np.count_nonzero(~full_correct & ~site_correct)
    correct_full = np.diag(full).copy()
    correct_site = np.diag(site).copy()
    return {
        "confusion_full": full, "confusion_site": site,
        "correct_full": correct_full, "correct_site": correct_site,
        "wilson_full": [wilson_interval(int(value), repeats) for value in correct_full],
        "wilson_site": [wilson_interval(int(value), repeats) for value in correct_site],
        **paired,
    }


def _fisher_at_times(probabilities: np.ndarray, slopes: np.ndarray,
                     shots_per_time: int) -> dict:
    full = fisher_matrix(probabilities, slopes, shots_per_time)
    site = fisher_matrix(probabilities, slopes, shots_per_time, site=2)
    centre = fisher_matrix(probabilities, slopes, shots_per_time, site=3)
    centre_singular = np.linalg.svd(centre, compute_uv=False)
    tolerance = 64 * np.finfo(float).eps * len(centre_singular) * centre_singular[0]
    return {
        "fisher_full_eigenvalues": np.linalg.eigvalsh(full),
        "fisher_site_eigenvalues": np.linalg.eigvalsh(site),
        "fisher_centre_eigenvalues": np.linalg.eigvalsh(centre),
        "centre_fisher_rank": int(np.count_nonzero(centre_singular > tolerance)),
        "fisher_dominance_residual": float(np.linalg.eigvalsh(full - site)[0]),
        "fisher_dominance_budget": float(128 * np.finfo(float).eps *
                                         max(1.0, np.linalg.norm(full, ord=2))),
    }


def run_fixture() -> dict:
    """Predeclared N=7 protocol; no accuracy threshold enters its checks."""
    candidates = np.unique(np.concatenate([
        np.linspace(left, right, 5)
        for left, right in zip(BASELINE_TIMES[:-1], BASELINE_TIMES[1:])
    ]))
    baseline_indices = tuple(int(np.where(candidates == time)[0][0])
                             for time in BASELINE_TIMES)
    cases = hypotheses()
    nominal = {}
    # Freeze both time schedules using the unperturbed, full-position Fisher
    # before computing any finite alternative or drawing any counts.
    for name, watched_sites in PROFILES.items():
        p = finite_probabilities(candidates, watched_sites, np.zeros(N - 1))
        _validate_probabilities(p)
        slopes = response_tensor(N, candidates, J_HOP, GAMMA0, watched_sites)
        per_time = np.array([
            fisher_matrix(p[index:index + 1], slopes[index:index + 1],
                          SHOTS_PER_TIME)
            for index in range(len(candidates))
        ])
        nominal[name] = {
            "probabilities": p, "slopes": slopes,
            "indices": {
                "baseline": baseline_indices,
                "selected": select_schedule(per_time, baseline_indices),
            },
        }

    checks = {}
    output = {}
    checks["candidate_count_37"] = len(candidates) == 37
    for profile_index, (name, watched_sites) in enumerate(PROFILES.items()):
        p = nominal[name]["probabilities"]
        slopes = nominal[name]["slopes"]
        scale = max(1.0, float(np.max(np.abs(slopes))))
        checks[f"{name}_each_bond_population_slope_sums_zero"] = bool(
            np.max(np.abs(slopes.sum(axis=2))) <=
            128 * np.finfo(float).eps * scale)
        rows = {"watched_sites": watched_sites,
                "sum_gamma": len(watched_sites) * GAMMA0}
        for schedule_index, schedule in enumerate(("baseline", "selected")):
            chosen = nominal[name]["indices"][schedule]
            times = candidates[list(chosen)]
            statistics = _fisher_at_times(p[list(chosen)], slopes[list(chosen)],
                                          SHOTS_PER_TIME)
            models = np.array([
                finite_probabilities(times, watched_sites, delta)
                for _, delta in cases
            ])
            _validate_probabilities(models)
            observations = evaluate_counts(models, SHOTS_PER_TIME, REPEATS,
                                           seed_prefix=(20260923, profile_index,
                                                        schedule_index))
            run = {
                "times": times.tolist(),
                "shots_total": int(len(times) * SHOTS_PER_TIME),
                **statistics, **observations,
            }
            rows[schedule] = run
            key = f"{name}_{schedule}"
            checks[f"{key}_budget"] = run["shots_total"] == 1000000
            checks[f"{key}_full_dominates_site_fisher"] = (
                run["fisher_dominance_residual"] >=
                -run["fisher_dominance_budget"])
            checks[f"{key}_centre_rank_three"] = run["centre_fisher_rank"] == 3
            checks[f"{key}_confusion_accounting"] = bool(
                np.all(run["confusion_full"].sum(axis=1) == REPEATS) and
                np.all(run["confusion_site"].sum(axis=1) == REPEATS) and
                np.all(sum(run[paired_key] for paired_key in (
                    "paired_both", "paired_full_only", "paired_site_only",
                    "paired_neither")) == REPEATS))
            tolerance = 128 * np.finfo(float).eps * max(1.0, np.max(models))
            checks[f"{key}_centre_reflection"] = bool(all(
                np.max(np.abs(models[1 + 2 * bond + offset, :, 3] -
                              models[1 + 2 * (N - 2 - bond) + offset, :, 3]))
                <= tolerance
                for bond in range((N - 1) // 2) for offset in (0, 1)
            ))
            checks[f"{key}_offcentre_end_bonds_separate"] = bool(
                np.max(np.abs(models[2, :, 2] - models[12, :, 2])) >
                1000 * tolerance)
        checks[f"{name}_selected_not_worse"] = (
            rows["selected"]["fisher_full_eigenvalues"][0] >=
            rows["baseline"]["fisher_full_eigenvalues"][0])
        output[name] = rows

    return {
        "gamma0": GAMMA0, "hopping": J_HOP, "q": J_HOP / GAMMA0,
        "baseline_times": BASELINE_TIMES.tolist(),
        "candidate_times": candidates.tolist(),
        "shots_per_time": SHOTS_PER_TIME, "repeats": REPEATS,
        "seed_root": 20260923,
        "hypotheses": [label for label, _ in cases],
        "profiles": output, "checks": checks,
    }


def format_run(result: dict) -> str:
    """Stable, fully labeled text for the retained local simulation run."""
    def times(values):
        return "[" + ", ".join(f"{value:g}" for value in values) + "]"

    def spectrum(values):
        return "[" + ", ".join(f"{value:.9g}" for value in values) + "]"

    lines = [
        "=== N=7 finite-shot bond and seat reading ===",
        "MODEL: ideal single-excitation XY chain; simulated Z-position shots; no hardware data",
        f"gamma0={result['gamma0']}; J_hop={result['hopping']}; Q={result['q']:g}; "
        "H_SE=J_hop*A_path+deltaJ_b*V_b",
        "H=(J_hop/2)*sum(XX+YY); D_ab=-2(gamma_a+gamma_b) for a!=b",
        "PREPARATION: nominal psi1 held fixed under every finite bond defect",
        "PROTOCOL: null or one bond b=0..5 changed by deltaJ_b=+/-0.0075 "
        "(10% J_hop); candidate signs and magnitude fixed, realized sign unknown; "
        "gamma, time, preparation and detector map known",
        "HYPOTHESIS ORDER: " + ", ".join(result["hypotheses"]),
        f"SHOTS: {result['shots_per_time']} independent multinomial Z-position "
        f"draws per time; 10 times = {10 * result['shots_per_time']} "
        f"per experiment; {result['repeats']} experiments per hypothesis",
        f"SEED: SeedSequence([{result['seed_root']}, profile_index, "
        "schedule_index, truth_index]); profile_index=uniform:0,single_centre:1; "
        "schedule_index=baseline:0,selected:1",
        f"RUNTIME: numpy={np.__version__}; scipy={scipy.__version__}",
        "BASELINE TIMES: " + times(result["baseline_times"]),
        f"CANDIDATE POOL: {len(result['candidate_times'])} times from "
        "fourfold subdivision of baseline intervals",
        "TIME RULE: deterministic one-for-one exchange locally maximizes the smallest "
        "nominal full-position Fisher eigenvalue; schedules frozen before "
        "finite-defect models and counts",
        "SCORES: categorical seven-outcome log likelihood versus site-2 "
        "Bernoulli marginal of the same shots; exact tie takes first hypothesis",
        "INTERVALS: two-sided Wilson 95% for 256 repeats in each fixed truth row",
    ]
    for name in PROFILES:
        profile = result["profiles"][name]
        lines.append(f"PROFILE {name}: watched sites={list(profile['watched_sites'])}; "
                     f"sum_gamma={profile['sum_gamma']:g}")
        for schedule in ("baseline", "selected"):
            run = profile[schedule]
            lines.append(f"  SCHEDULE {schedule}: times={times(run['times'])}; "
                         f"shots={run['shots_total']}")
            lines.append("  nominal Fisher full eigenvalues=" +
                         spectrum(run["fisher_full_eigenvalues"]))
            lines.append("  nominal Fisher site2 eigenvalues=" +
                         spectrum(run["fisher_site_eigenvalues"]))
            lines.append("  nominal Fisher centre3 eigenvalues=" +
                         spectrum(run["fisher_centre_eigenvalues"]) +
                         f"; rank={run['centre_fisher_rank']}")
            lines.append(f"  Fisher full-site2 minimum eigenvalue="
                         f"{run['fisher_dominance_residual']:.9g}; "
                         f"roundoff budget={run['fisher_dominance_budget']:.9g}")
            for detector, confusion_key, correct_key, interval_key in (
                ("full", "confusion_full", "correct_full", "wilson_full"),
                ("site2", "confusion_site", "correct_site", "wilson_site"),
            ):
                lines.append(f"  {detector.upper()} CONFUSION: columns follow HYPOTHESIS ORDER")
                confusion = run[confusion_key]
                for truth, label in enumerate(result["hypotheses"]):
                    lo, hi = run[interval_key][truth]
                    row = ",".join(str(int(value)) for value in confusion[truth])
                    lines.append(f"    true={label:4s}: {row} | "
                                 f"correct={int(run[correct_key][truth])}/"
                                 f"{result['repeats']}; Wilson95=[{lo:.6f},{hi:.6f}]")
                lines.append(f"  {detector.upper()} TOTAL: "
                             f"{int(np.sum(run[correct_key]))}/"
                             f"{len(result['hypotheses']) * result['repeats']} "
                             "fixed-stratum trials (descriptive only)")
            lines.append("  PAIRED CORRECTNESS: both,full_only,site2_only,neither "
                         "from the same records")
            for truth, label in enumerate(result["hypotheses"]):
                counts = ",".join(str(int(run[key][truth])) for key in (
                    "paired_both", "paired_full_only", "paired_site_only",
                    "paired_neither"))
                lines.append(f"    true={label:4s}: {counts}")
    lines.append("CONSISTENCY CHECKS (no favorable accuracy threshold):")
    for key, value in result["checks"].items():
        lines.append(f"  {key}={'PASS' if value else 'FAIL'}")
    lines.append("SCOPE: local ideal known-model classification and nominal "
                 "six-parameter Fisher; not unrestricted coupling recovery, "
                 "T1/leakage robustness, or a hardware confirmation")
    lines.append("VERDICT: " + ("PASS" if all(result["checks"].values())
                                else "A CHECK FIRED"))
    return "\n".join(lines) + "\n"


def main() -> int:
    result = run_fixture()
    print(format_run(result), end="")
    return 0 if all(result["checks"].values()) else 1


if __name__ == "__main__":
    sys.exit(main())
