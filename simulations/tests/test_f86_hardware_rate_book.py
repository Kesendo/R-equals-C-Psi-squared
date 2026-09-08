import sys
from pathlib import Path

import numpy as np


SIM = Path(__file__).resolve().parents[1]
if str(SIM) not in sys.path:
    sys.path.insert(0, str(SIM))

import at_the_target as target
import f86_hardware_rate_book as rate_book
import journey_control as control
import post_ep_dynamics_4d as flow


def test_hardware_population_scan_is_recomputed_from_raw_records():
    q_label, q_lindblad, revival = rate_book.population_scan()
    assert q_label == [0.5, 1.0, 1.5, 2.5, 5.0, 20.0]
    assert q_lindblad == [1.0, 2.0, 3.0, 5.0, 10.0, 40.0]
    assert np.allclose(revival, [0.2978515625, 0.3623809814453125, 0.34356689453125,
                                 0.4898834228515625, 0.560699462890625, 0.70330810546875])


def test_twirl_simulation_endpoint_is_028_not_truncated_subscan_031():
    _, _, revival = rate_book.population_scan(rate_book.TWIRL_SCAN, "revival_max_n0")
    assert round(revival[0], 2) == 0.28
    assert round(revival[-1], 2) == 0.84
    assert revival[-1] == 0.8417853730254796


def test_part_a_endpoint_is_finite_time_unnormalized_marginal_sighting():
    time_us, populations = rate_book.part_a_endpoint()
    assert time_us == 20.0
    assert np.allclose(populations, [0.339111328125, 0.426025390625, 0.338134765625])
    assert not np.isclose(sum(populations), 1.0)


def _expected_dimensionless_liouvillian(module, n: int, q: float, f: float = 1.0):
    d = 2 ** n
    ident = np.eye(d)
    h_unit = sum(module.bond_op(n, b, module.X, module.X) + module.bond_op(n, b, module.Y, module.Y)
                 for b in range(n - 1))
    expected = -1j * (q / 2.0) * (np.kron(ident, h_unit) - np.kron(h_unit.T, ident))
    for site in range(n):
        z_site = module.op_at(n, site, module.Z)
        expected += f * (np.kron(z_site, z_site) - np.kron(ident, ident))
    return expected


def test_flow_producer_uses_canonical_q_half_coefficient():
    assert np.allclose(flow.liouvillian_dimensionless(2, 3.0),
                       _expected_dimensionless_liouvillian(flow, 2, 3.0))


def test_target_producer_uses_canonical_q_half_coefficient():
    assert np.allclose(target.liouvillian(2, 3.0, 1.0),
                       _expected_dimensionless_liouvillian(target, 2, 3.0))


def test_control_domain_contains_handover_endpoints_and_uses_per_site_time_book():
    starts = [3.0, 5.0, 10.0, 20.0, control.Q_BASE]
    lo, hi = control.control_domain(starts)
    assert all(lo <= control.inj_for_q(q) <= hi for q in starts)
    assert hi == control.inj_for_q(3.0)
    source = Path(control.__file__).read_text(encoding="utf-8")
    assert "tau = gamma_per_site * t" in source
    assert "gamma_total" not in source
