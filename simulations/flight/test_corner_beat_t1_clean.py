"""Tests for the T1-CLEAN us->gbar bridge.

WHY THIS FILE EXISTS. Section 9 promotes T1-CLEAN to a registered gate
and gives its conversion in prose,

    gamma_l^hw = tau_step / (4 * T1_l * dt)

the Gamma_l/4-equivalent dephasing in the simulation's J units, read in
gbar as gamma_l^hw / gbar. The runner enforced that `tau_step_us` is
PRESENT in the day-of addendum and then never read it: as of 2026-08-18
no line in either repository computes that conversion. The scale face
would therefore have been evaluated by hand, from prose, across two
repositories, on flight day, once. This file makes the arithmetic code.

The bridge's own history is the argument for it. It has already eaten a
unit repair (round 14, the formula returns J units) and a missing
division by gbar (round 3); both were found because someone looked, not
because anything failed.

THE TWO DEFINITIONS THIS FILE PINS, neither of them invented here:
 - the conversion and its two anchor values are section 9's, quoted
   there as ~0.0033 J units at tau_step = 0.5 us and T1 = 250 us, and
   ~0.033 gbar at the flown Q = 10;
 - "scale" is the committed gate's, simulations/corner_beat_gate.py in
   the research repo: "'scale' = Gam_mean/4 in units gbar", drawn as
   uniform(0, 2*scale*gbar), whose MEAN is scale*gbar. It is a mean, not
   a spread; the G4 v2 rows that produced the <= 0.25 gbar bound were
   generated from that draw.
"""

import importlib.util as _ilu
from pathlib import Path

import pytest

HERE = Path(__file__).parent


def _load(name, path):
    spec = _ilu.spec_from_file_location(name, path)
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rcb = _load("rcb_t1_clean_test", HERE / "run_corner_beat.py")


# ------------------------------------------------- the conversion itself
def test_profile_reproduces_the_section_9_anchor():
    """Section 9's own worked value: tau_step = 0.5 us against T1 = 250 us
    is ~0.0033 in J units and ~0.033 in gbar at the flown Q = 10."""
    prof = rcb.t1_clean_profile(0.5, [250.0] * 6)
    assert prof["gamma_hw_j"] == pytest.approx([0.5 / 150.0] * 6)
    assert prof["gamma_hw_j"][0] == pytest.approx(0.003333, abs=5e-7)
    assert prof["gamma_hw_gbar"][0] == pytest.approx(0.03333, abs=5e-6)


def test_profile_is_per_site_and_inverse_in_T1():
    """A station with half the T1 carries twice the equivalent dephasing;
    the profile is per-site, not a chain scalar."""
    prof = rcb.t1_clean_profile(0.5, [250.0, 125.0, 500.0, 250.0, 250.0, 250.0])
    g = prof["gamma_hw_j"]
    assert g[1] == pytest.approx(2 * g[0])
    assert g[2] == pytest.approx(0.5 * g[0])


def test_profile_scales_linearly_in_tau_step():
    """tau_step is the wall duration per Strang step; doubling it doubles
    the T1 dose the same circuit pays."""
    a = rcb.t1_clean_profile(0.5, [250.0] * 6)["gamma_hw_j"][0]
    b = rcb.t1_clean_profile(1.0, [250.0] * 6)["gamma_hw_j"][0]
    assert b == pytest.approx(2 * a)


# ------------------------------------------------------- the scale face
def test_scale_is_the_mean_in_gbar_not_a_spread():
    """The committed gate draws the background as uniform(0, 2*scale*gbar),
    so 'scale' is the profile MEAN in gbar units. A uniform profile has a
    scale equal to its own per-site value and a spread of zero; if the
    implementation ever returns a spread, this fails."""
    prof = rcb.t1_clean_profile(0.5, [250.0] * 6)
    assert prof["scale_gbar"] == pytest.approx(prof["gamma_hw_gbar"][0])
    assert prof["scale_gbar"] == pytest.approx(0.03333, abs=5e-6)


def test_scale_averages_over_stations():
    """Two stations at 125 us and four at 500 us: the scale is the mean of
    the six equivalents, which a worst-station reading would overstate."""
    t1 = [125.0, 125.0, 500.0, 500.0, 500.0, 500.0]
    prof = rcb.t1_clean_profile(0.5, t1)
    expected = sum(0.5 / (4.0 * t * rcb.JDT) for t in t1) / 6.0 / rcb.GBAR
    assert prof["scale_gbar"] == pytest.approx(expected)
    assert prof["scale_gbar"] < max(prof["gamma_hw_gbar"])


# -------------------------------------------------------------- the gate
def test_scale_face_passes_under_the_registered_quarter():
    """G4 v2's registered bound: scattered Gamma_l/4-equivalent profile
    scale <= 0.25 gbar."""
    ok, prof = rcb.t1_clean_scale_face(0.5, [250.0] * 6)
    assert ok
    assert prof["scale_gbar"] < 0.25


def test_scale_face_fails_above_the_registered_quarter():
    """T1 an order of magnitude worse pushes the mean past the bound and
    the face must say so rather than round toward passing."""
    ok, prof = rcb.t1_clean_scale_face(0.5, [25.0] * 6)
    assert not ok
    assert prof["scale_gbar"] > 0.25


def test_scale_face_is_closed_at_the_bound():
    """A scale landing exactly on 0.25 passes: the registered wording is
    'scale <= 0.25 gbar', and a boundary that silently excludes its own
    edge is how a gate acquires an unregistered margin."""
    t1_at_bound = 0.5 / (4.0 * rcb.JDT * 0.25 * rcb.GBAR)
    ok, prof = rcb.t1_clean_scale_face(0.5, [t1_at_bound] * 6)
    assert prof["scale_gbar"] == pytest.approx(0.25)
    assert ok


# ------------------------------------------------------------- the book
def test_the_conversion_is_the_lindblad_book():
    """Gamma/4 is the Lindblad-book T1 contribution to Z-dephasing: a T1
    process at rate Gamma costs Gamma/2 of coherence rate, and the
    repository's canonical gamma is half a coherence rate (GLOSSARY, the
    T2 -> gamma conversion). Writing Gamma/2 here, the coherence book,
    would double every station. Pinned so that a future edit that
    'simplifies' the 4 has to argue with a test."""
    t1 = 200.0
    tau = 0.4
    prof = rcb.t1_clean_profile(tau, [t1] * 6)
    coherence_book = tau / (2.0 * t1 * rcb.JDT)
    assert prof["gamma_hw_j"][0] == pytest.approx(coherence_book / 2.0)


def test_rejects_a_station_count_that_is_not_the_chain():
    """Six stations. A calibration snapshot of the wrong length is an
    executor error on flight day, not something to broadcast over."""
    with pytest.raises(ValueError):
        rcb.t1_clean_profile(0.5, [250.0] * 5)


def test_rejects_a_nonpositive_t1():
    with pytest.raises(ValueError):
        rcb.t1_clean_profile(0.5, [250.0, 0.0, 250.0, 250.0, 250.0, 250.0])


# ------------------------------------- reading it out of the addendum
def test_reads_t1_from_the_snapshot_in_line_order():
    """The runner writes snap["qubits"][str(q)]["T1_us"]; the face must
    read the SAME structure, and in the flown line's order, since the
    profile is per-station and the stations are the line."""
    snap = {"qubits": {"13": {"T1_us": 250.0}, "14": {"T1_us": 125.0},
                       "15": {"T1_us": 500.0}, "16": {"T1_us": 250.0},
                       "17": {"T1_us": 250.0}, "18": {"T1_us": 250.0}}}
    t1 = rcb.t1_list_from_snapshot(snap, [13, 14, 15, 16, 17, 18])
    assert t1 == [250.0, 125.0, 500.0, 250.0, 250.0, 250.0]


def test_snapshot_order_follows_the_chain_not_the_dict():
    """Reversing the line reverses the profile: a dict has no order and
    reading one would silently mis-assign stations."""
    snap = {"qubits": {"13": {"T1_us": 100.0}, "14": {"T1_us": 200.0},
                       "15": {"T1_us": 300.0}, "16": {"T1_us": 400.0},
                       "17": {"T1_us": 500.0}, "18": {"T1_us": 600.0}}}
    fwd = rcb.t1_list_from_snapshot(snap, [13, 14, 15, 16, 17, 18])
    rev = rcb.t1_list_from_snapshot(snap, [18, 17, 16, 15, 14, 13])
    assert rev == list(reversed(fwd))


def test_missing_station_is_an_error_not_a_default():
    """A station absent from the snapshot must abort, never silently
    become a nominal T1: that is how a bad line passes a gate."""
    snap = {"qubits": {"13": {"T1_us": 250.0}}}
    with pytest.raises(KeyError):
        rcb.t1_list_from_snapshot(snap, [13, 14, 15, 16, 17, 18])
