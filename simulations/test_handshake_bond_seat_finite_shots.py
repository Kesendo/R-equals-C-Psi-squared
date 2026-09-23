"""Independent physical and statistical controls for finite-shot bond reading."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


sys.path.insert(0, str(Path(__file__).parent))


def producer():
    return importlib.import_module("handshake_bond_seat_finite_shots")


def test_categorical_and_binary_fisher_on_three_outcome_two_bond_toy():
    probabilities = np.array([[0.2, 0.3, 0.5]])
    slopes = np.array([[[0.1, 0.0, -0.1],
                        [0.0, 0.2, -0.2]]])
    full = producer().fisher_matrix(probabilities, slopes, 1)
    site_2 = producer().fisher_matrix(probabilities, slopes, 1, site=2)
    np.testing.assert_allclose(
        full, [[0.07, 0.04], [0.04, 0.21333333333333335]],
        rtol=0, atol=4 * np.finfo(float).eps,
    )
    np.testing.assert_allclose(
        site_2, [[0.04, 0.08], [0.08, 0.16]],
        rtol=0, atol=4 * np.finfo(float).eps,
    )
    assert np.linalg.eigvalsh(full - site_2)[0] > -4 * np.finfo(float).eps


def test_finite_defect_populations_match_explicit_site_z_matrix_ode():
    n, hopping, gamma0 = 7, 0.075, 0.05
    times = np.array([2.0, 10.0])
    delta = np.zeros(n - 1)
    delta[0] = 0.0075
    psi = np.sqrt(2.0 / (n + 1)) * np.sin(
        np.pi * np.arange(1, n + 1) / (n + 1)
    )
    initial = np.outer(psi, psi).astype(complex)
    h = np.zeros((n, n))
    for bond in range(n - 1):
        h[bond, bond + 1] = h[bond + 1, bond] = hopping + delta[bond]

    for watched_sites in (tuple(range(n)), (3,)):
        z_operators = []
        for site in watched_sites:
            signs = np.ones(n)
            signs[site] = -1.0
            z_operators.append(np.diag(signs))

        def rhs(_, flat):
            rho = flat.reshape(n, n)
            derivative = -1j * (h @ rho - rho @ h)
            for z in z_operators:
                derivative += gamma0 * (z @ rho @ z - rho)
            return derivative.ravel()

        rtol, atol = 2e-12, 2e-14
        solved = solve_ivp(rhs, (0.0, times[-1]), initial.ravel(),
                           t_eval=times, method="DOP853", rtol=rtol, atol=atol)
        assert solved.success
        direct = np.diagonal(solved.y.T.reshape(len(times), n, n),
                             axis1=1, axis2=2).real
        actual = producer().finite_probabilities(times, watched_sites, delta)
        np.testing.assert_allclose(actual, direct, rtol=0,
                                   atol=8 * (atol + rtol * np.max(direct)))
        np.testing.assert_allclose(actual.sum(axis=1), 1.0, rtol=0,
                                   atol=64 * np.finfo(float).eps)
        assert np.min(actual) > 0.0


def test_halved_channel_rate_breaks_the_independent_physical_control(monkeypatch):
    module = producer()
    delta = np.zeros(6)
    delta[0] = 0.0075
    reference = module.finite_probabilities([10.0], tuple(range(7)), delta)
    monkeypatch.setattr(module, "GAMMA0", 0.025)
    wrong_rate = module.finite_probabilities([10.0], tuple(range(7)), delta)
    assert np.max(np.abs(wrong_rate - reference)) > 0.001


def test_exchange_schedule_improves_a_known_weak_direction_deterministically():
    per_time = np.array([np.diag([10.0, 0.0]),
                         np.diag([0.0, 1.0]),
                         np.diag([0.0, 10.0])])
    selected = producer().select_schedule(per_time, (0, 1))
    assert selected == (0, 2)
    assert producer().select_schedule(per_time, (0, 1)) == selected


def test_nominal_full_record_schedule_and_centre_reflection_controls():
    from handshake_bond_seat_fixed_gamma import response_tensor

    module = producer()
    baseline = np.array([2., 4., 6., 8., 10., 12., 15., 18., 24., 30.])
    candidates = np.unique(np.concatenate([
        np.linspace(left, right, 5)
        for left, right in zip(baseline[:-1], baseline[1:])
    ]))
    assert len(candidates) == 37
    baseline_indices = tuple(int(np.where(candidates == time)[0][0])
                             for time in baseline)
    watched_sites = tuple(range(7))
    p = module.finite_probabilities(candidates, watched_sites, np.zeros(6))
    slopes = response_tensor(7, candidates, 0.075, 0.05, watched_sites)
    per_time = np.array([
        module.fisher_matrix(p[t:t + 1], slopes[t:t + 1], 1)
        for t in range(len(candidates))
    ])
    selected = module.select_schedule(per_time, baseline_indices)
    assert len(selected) == len(set(selected)) == 10
    assert all(0 <= index < len(candidates) for index in selected)
    old_min = np.linalg.eigvalsh(np.sum(per_time[list(baseline_indices)], axis=0))[0]
    new_min = np.linalg.eigvalsh(np.sum(per_time[list(selected)], axis=0))[0]
    assert new_min >= old_min

    base_p = p[list(baseline_indices)]
    base_slopes = slopes[list(baseline_indices)]
    full = module.fisher_matrix(base_p, base_slopes, 100000)
    site_2 = module.fisher_matrix(base_p, base_slopes, 100000, site=2)
    centre = module.fisher_matrix(base_p, base_slopes, 100000, site=3)
    budget = (128 * np.finfo(float).eps *
              max(1.0, np.linalg.norm(full, ord=2)))
    assert np.linalg.eigvalsh(full - site_2)[0] >= -budget
    singular = np.linalg.svd(centre, compute_uv=False)
    centre_tol = 64 * np.finfo(float).eps * len(singular) * singular[0]
    assert int(np.count_nonzero(singular > centre_tol)) == 3

    left = np.zeros(6)
    right = np.zeros(6)
    left[0] = right[5] = 0.0075
    pl = module.finite_probabilities(baseline, watched_sites, left)
    pr = module.finite_probabilities(baseline, watched_sites, right)
    roundoff = 128 * np.finfo(float).eps * max(np.max(pl), np.max(pr))
    assert np.max(np.abs(pl[:, 3] - pr[:, 3])) <= roundoff
    assert abs(pl[4, 2] - pr[4, 2]) > 1000 * roundoff


def test_hypothesis_order_and_single_bond_strength_are_declared():
    cases = producer().hypotheses()
    assert [label for label, _ in cases] == [
        "null", "b0-", "b0+", "b1-", "b1+", "b2-", "b2+",
        "b3-", "b3+", "b4-", "b4+", "b5-", "b5+",
    ]
    np.testing.assert_array_equal(cases[0][1], np.zeros(6))
    for bond in range(6):
        for sign, offset in ((-1, 1), (+1, 2)):
            expected = np.zeros(6)
            expected[bond] = sign * 0.0075
            np.testing.assert_array_equal(cases[2 * bond + offset][1],
                                          expected)


def test_full_likelihood_sees_outcome_split_that_binary_site_loses():
    models = np.array([[[0.2, 0.3, 0.5]], [[0.3, 0.2, 0.5]]])
    counts = np.array([[[3, 2, 5]]])
    full = producer().log_likelihoods(counts, models)
    site = producer().log_likelihoods(counts, models, site=2)
    assert full.shape == site.shape == (1, 2)
    assert full[0, 1] > full[0, 0]
    assert site[0, 0] == site[0, 1]
    np.testing.assert_array_equal(producer().classify(counts, models), [1])
    np.testing.assert_array_equal(producer().classify(counts, models, site=2),
                                  [0])


def test_zero_bins_are_not_smoothed_and_impossible_observations_lose():
    models = np.array([[[1.0, 0.0, 0.0]], [[0.5, 0.5, 0.0]]])
    counts = np.array([[[0, 1, 0]]])
    scores = producer().log_likelihoods(counts, models)
    assert scores[0, 0] == -np.inf
    assert np.isfinite(scores[0, 1])
    np.testing.assert_array_equal(producer().classify(counts, models), [1])


def test_multinomial_draws_conserve_shots_and_share_single_site_record():
    probabilities = np.array([[0.2, 0.3, 0.5], [0.1, 0.6, 0.3]])
    counts = producer().draw_counts(probabilities, 100000, 256,
                                    seed_parts=(20260923, 0, 0, 0))
    assert counts.shape == (256, 2, 3)
    np.testing.assert_array_equal(counts.sum(axis=2),
                                  np.full((256, 2), 100000))
    assert np.array_equal(counts, producer().draw_counts(
        probabilities, 100000, 256, seed_parts=(20260923, 0, 0, 0)
    ))


def test_detector_only_reflection_swaps_the_physical_bond_label():
    module = producer()
    grid = np.array([2.0, 4.0, 10.0, 18.0, 30.0])
    left = np.zeros(6)
    right = np.zeros(6)
    left[0] = right[5] = 0.0075
    models = np.array([
        module.finite_probabilities(grid, tuple(range(7)), left),
        module.finite_probabilities(grid, tuple(range(7)), right),
    ])
    asimov = (1000000 * models[0])[None, :, :]
    np.testing.assert_array_equal(module.classify(asimov, models), [0])
    reflected_detector_only = asimov[:, :, ::-1]
    np.testing.assert_array_equal(
        module.classify(reflected_detector_only, models), [1]
    )
    np.testing.assert_array_equal(reflected_detector_only[:, :, 3],
                                  asimov[:, :, 3])


def test_wilson_intervals_cover_boundary_counts_without_normal_shortcut():
    module = producer()
    lower, upper = module.wilson_interval(0, 256)
    z = 1.959963984540054
    assert lower == 0.0
    np.testing.assert_allclose(upper, z * z / (256 + z * z), rtol=0,
                               atol=4 * np.finfo(float).eps)
    lower, upper = module.wilson_interval(128, 256)
    np.testing.assert_allclose(lower + upper, 1.0, rtol=0,
                               atol=4 * np.finfo(float).eps)


def test_confusions_and_paired_outcomes_use_identical_multinomial_records():
    module = producer()
    models = np.array([
        [[0.2, 0.3, 0.5], [0.3, 0.4, 0.3]],
        [[0.3, 0.2, 0.5], [0.4, 0.3, 0.3]],
    ])
    result = module.evaluate_counts(models, 100, 64, seed_prefix=(17, 2, 4))
    for truth in range(2):
        records = module.draw_counts(models[truth], 100, 64,
                                     seed_parts=(17, 2, 4, truth))
        full = module.classify(records, models)
        site = module.classify(records, models, site=2)
        np.testing.assert_array_equal(result["confusion_full"][truth],
                                      np.bincount(full, minlength=2))
        np.testing.assert_array_equal(result["confusion_site"][truth],
                                      np.bincount(site, minlength=2))
        assert result["paired_full_only"][truth] == np.count_nonzero(
            (full == truth) & (site != truth))
        assert result["paired_site_only"][truth] == np.count_nonzero(
            (full != truth) & (site == truth))
        assert sum(result[key][truth] for key in (
            "paired_both", "paired_full_only", "paired_site_only",
            "paired_neither")) == 64


def test_complete_fixture_fixes_budget_schedule_and_consistency_gates():
    result = producer().run_fixture()
    assert result["gamma0"] == 0.05
    assert result["hopping"] == 0.075
    assert result["shots_per_time"] == 100000
    assert result["repeats"] == 256
    assert len(result["hypotheses"]) == 13
    assert len(result["candidate_times"]) == 37
    assert all(result["checks"].values())
    for profile in ("uniform", "single_centre"):
        rows = result["profiles"][profile]
        assert len(rows["selected"]["times"]) == 10
        assert rows["selected"]["fisher_full_eigenvalues"][0] >= (
            rows["baseline"]["fisher_full_eigenvalues"][0])
        for schedule in ("baseline", "selected"):
            run = rows[schedule]
            assert len(run["times"]) == 10
            assert run["shots_total"] == 1000000
            assert run["confusion_full"].shape == (13, 13)
            assert run["confusion_site"].shape == (13, 13)
            np.testing.assert_array_equal(run["confusion_full"].sum(axis=1),
                                          np.full(13, 256))
            np.testing.assert_array_equal(run["confusion_site"].sum(axis=1),
                                          np.full(13, 256))
            assert len(run["wilson_full"]) == 13
            assert len(run["wilson_site"]) == 13
            assert run["centre_fisher_rank"] == 3
