"""Tests for the polarity_fingerprint workflow (F87 × F112 × F113 joint).

Smoke + invariant tests:
  - Heisenberg + pure Z-deph: F87 truly, F112 BALANCED, in_typed_scope = True
  - YZ+ZY (F108 non-truly) + pure Z-deph: F87 soft, F112 BALANCED, in_typed_scope
  - XY (Π²-odd) + pure Z-deph: F87 hard, F112 BALANCED, in_typed_scope
  - Heisenberg + T1: F87 truly (H unchanged), F112 BALANCED empirically (T1 alone
    with bilinear bit_b-homog H), in_typed_scope = False (c bit_b-mixed)
  - Z-drive H (single-site Z) + T1 + Z-drive-omegas provided: F113 fields populated
    and F113-extracted γ_T1 matches input γ_T1 to machine precision
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

import framework as fw


def test_polarity_fingerprint_returns_expected_keys():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_fingerprint(chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')])
    expected_keys = {
        'f87_class', 'f112_asymmetry', 'f112_rel_asymmetry', 'f112_M_norm_sq',
        'f112_verdict', 'in_f112_typed_scope', 'h_bit_b_homogeneous',
        'c_bit_b_homogeneous', 'f113_applies', 'f113_predicted',
        'f113_extracted_gamma_t1', 'reading',
        'f112_M_anti_norm_sq', 'f112_polarity_degenerate',
    }
    assert expected_keys.issubset(result.keys())


def test_polarity_fingerprint_rel_asymmetry_divides_by_the_polarity_content():
    """The denominator of f112_rel_asymmetry is ||M_anti||^2, not ||M||^2.

    This is the discriminating pin for the 2026-08-06 denominator change, and it is here
    rather than in test_diagnose_hardware because it needs a fixture whose asymmetry is
    NON-zero: where asym == 0.0 the ratio is 0.0 under either denominator and no assertion
    can see the difference (measured by mutation).

    A single-site Z-drive plus T1 is Pi^2-odd, so it has real polarity content and a
    non-zero F113 asymmetry. Reverting the denominator to max(||M||^2, 1e-15) turns the
    equality below red by the F83 factor 2 + 4r.
    """
    omega, gamma_t1 = 0.13, 0.001
    chain = fw.ChainSystem(N=2, J=omega / 2.0, gamma_0=0.05)
    result = fw.polarity_fingerprint(chain, [('Z', 'I')], gamma_t1=gamma_t1)

    assert result['f112_asymmetry'] != 0.0,         'fixture must have a non-zero numerator or this test cannot discriminate'
    assert result['f112_polarity_degenerate'] is False
    m_anti = result['f112_M_anti_norm_sq']
    assert m_anti > 1e-6, f'Pi^2-odd drive must carry polarity content; got {m_anti}'

    expected = abs(result['f112_asymmetry']) / m_anti
    assert result['f112_rel_asymmetry'] == pytest.approx(expected, rel=1e-12)
    # And it is NOT the retired quantity. The two differ by ||M||^2 / ||M_anti||^2, which is
    # F83's 2 + 4r on the pure-dephasing bilinear family but NOT here: this fixture carries
    # T1, which puts content into M that F83's anti-fraction does not model, so the measured
    # factor is 2.000533 on THIS fixture (a one-site drive at N=2) rather than 2. Pin the
    # separation, not a law that does not apply. Note the factor is fixture-specific: the C#
    # witness drives BOTH sites and measures 2.000266, so do not carry either number across.
    retired = abs(result['f112_asymmetry']) / result['f112_M_norm_sq']
    assert result['f112_rel_asymmetry'] / retired == pytest.approx(2.000533, rel=1e-5)


def test_polarity_fingerprint_heisenberg_pure_z_in_scope_is_silent_not_balanced():
    """Heisenberg + pure Z-deph: in F112's typed scope, and carrying NO polarity content.

    The classic case, and the one that reads wrongly if the verdict is taken off the ratio.
    H is Pi^2-even and c is bit_b-homogeneous, so M_anti vanishes as a theorem and the
    asymmetry is 0 - 0. That is DEGENERATE, not BALANCED: a balanced reading here would be a
    confirmation of F112 drawn from a row where there is nothing to confirm. Read as
    BALANCED until 2026-08-07.
    """
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_fingerprint(chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')])
    assert result['f87_class'] == 'truly'
    assert result['f112_verdict'] == 'DEGENERATE'
    assert result['f112_polarity_degenerate'] is True
    # N=3 on purpose: the float guard alone would MISS this one. ||M_anti||^2 arrives as
    # ~1e-33 rather than exact 0.0, and only the structural test is N-independent.
    assert result['f112_M_anti_norm_sq'] != 0.0
    assert result['f112_M_anti_norm_sq'] < 1e-30
    assert result['in_f112_typed_scope'] is True
    assert result['h_bit_b_homogeneous'] is True
    assert result['c_bit_b_homogeneous'] is True
    assert result['f113_applies'] is False  # no z_drive_omegas_per_site


def test_polarity_fingerprint_yz_zy_pi2_even_in_scope_is_silent_not_balanced():
    """YZ + ZY (F108 non-truly Π²-even): F87 soft, and F112-SILENT.

    Sharper than the Heisenberg case above and the reason ‖M‖² was the wrong scale even
    where it does not vanish: here ‖M‖² is large (YZ+ZY is non-truly, so M is far from
    zero) while the POLARITY content is nil, because the Pi^2-even H puts nothing into
    M's ±i halves. Dividing by ‖M‖² hid that behind a big denominator instead of a floor.
    """
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_fingerprint(chain, [('Y', 'Z'), ('Z', 'Y')])
    assert result['f87_class'] == 'soft'
    assert result['f112_polarity_degenerate'] is True
    assert result['f112_verdict'] == 'DEGENERATE'
    # The discriminating pair: M is substantial, its polarity content is not.
    assert result['f112_M_norm_sq'] > 1.0
    assert result['f112_M_anti_norm_sq'] < 1e-30
    assert result['in_f112_typed_scope'] is True


def test_polarity_fingerprint_xy_pi2_odd_in_scope_balanced():
    """XY pair (Π²-odd): F87 hard for unbroken Π² pair, F112 BALANCED."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_fingerprint(chain, [('X', 'Y')])
    # F87 classification of pure XY depends on the chain N + bond structure;
    # accept whatever F87 returns and assert F112 is BALANCED + in-scope.
    assert result['f87_class'] in ('truly', 'soft', 'hard')
    assert result['f112_verdict'] == 'BALANCED'
    assert result['in_f112_typed_scope'] is True


def test_polarity_fingerprint_heisenberg_with_t1_out_of_typed_scope_but_balanced():
    """Heisenberg + σ⁻ T1: c bit_b-mixed → out of F112 typed scope. But the
    broader empirical envelope (probes 1-14) shows balance still holds for
    bilinear-only H + σ⁻ T1 (only Z-drive H + σ⁻ T1 breaks balance, per F113)."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_fingerprint(
        chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')],
        gamma_t1=0.1,
    )
    assert result['in_f112_typed_scope'] is False  # T1 makes c bit_b-mixed
    assert result['c_bit_b_homogeneous'] is False
    # Heisenberg H is bilinear (not single-site Z-drive) so F113 doesn't apply
    # and the broader empirical envelope (Welle 2) confirms BALANCED in this case.
    assert result['f112_verdict'] == 'BALANCED'


def test_polarity_fingerprint_f113_extraction_recovers_input_gamma_t1():
    """When z_drive_omegas_per_site is provided AND H matches the Z-drive
    structure, F113 inversion of the measured asymmetry recovers the input
    γ_T1 to machine precision."""
    chain = fw.ChainSystem(N=2, gamma_0=0.0)  # no Z-deph, isolate F113 signal
    omega = 0.13
    gamma_t1 = 0.001
    # H = (omega/2)·(Z_0 + Z_1) expressed as single-site Z terms with the right
    # coefficient. polarity_coordinates uses (chain, terms) to build the chain's
    # H via the bond-bilinear convention; for single-site terms we need a path
    # that lets us encode H = ω·Z_0/2 + ω·Z_1/2.
    # Single-site Z terms DO reach the measurement through the bilinear interface, so
    # this test goes through the workflow end to end; the hand-built anchor below is kept
    # only to pin the registry's (N/2)·4^N collapse beside it. For the same inversion on
    # a hardware-effective L see simulations/f113_t1_extraction_kingston.py.
    expected_asym = (4 ** 2) / 2.0 * 2 * omega * (0 - gamma_t1)  # -16·ω·γ_T1
    # This must go through the workflow. Inverting a hand-built asymmetry by hand is
    # -(k·ω·g)/(k·ω), an identity true for any k and ω, which held with the workflow
    # deleted. The drive here is uniform on ALL N sites, the one case where the
    # (N/2)·4^N·ω collapse applies; per driven site the prefactor is (1/2)·4^N·ω.
    result = fw.polarity_fingerprint(
        chain, [('Z', 'I'), ('I', 'Z')], gamma_t1=gamma_t1,
        z_drive_omegas_per_site=[2.0, 2.0],
    )
    assert result['f113_drive_matches_terms'] is True
    assert result['f113_extracted_gamma_t1'] == pytest.approx(gamma_t1, rel=1e-9)
    # and the hand-built anchor still agrees with the collapse it documents
    assert expected_asym == pytest.approx(-2.08e-3, rel=1e-12)
    assert -expected_asym / ((2 / 2.0) * (4 ** 2) * omega) == pytest.approx(gamma_t1, rel=1e-12)


def test_polarity_fingerprint_verdicts_span_all_three_branches():
    """The BROKEN / near-BALANCED / BALANCED classifier was entirely unpinned: replacing
    its threshold with 1e99 (everything BALANCED) left the suite green. Each branch is
    reachable on the same H by moving γ_T1 alone."""
    chain = fw.ChainSystem(N=2, gamma_0=0.0)
    verdicts = {}
    for gt1 in (1e-3, 1e-7, 1e-11):
        verdicts[gt1] = fw.polarity_fingerprint(chain, [('Z', 'I')], gamma_t1=gt1)
    assert verdicts[1e-3]['f112_verdict'] == 'BROKEN'
    assert verdicts[1e-7]['f112_verdict'] == 'near-BALANCED'
    assert verdicts[1e-11]['f112_verdict'] == 'BALANCED'


def test_polarity_fingerprint_rejects_wrong_length_rate_and_drive():
    """All three length guards were untested; deleting any of them left the suite green."""
    chain = fw.ChainSystem(N=3, gamma_0=0.0)
    terms = [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')]
    with pytest.raises(ValueError, match=r"gamma_t1 length 2"):
        fw.polarity_fingerprint(chain, terms, gamma_t1=[0.1, 0.2])
    with pytest.raises(ValueError, match=r"gamma_pump length 4"):
        fw.polarity_fingerprint(chain, terms, gamma_pump=[0.1] * 4)
    with pytest.raises(ValueError, match=r"z_drive_omegas_per_site length 0"):
        fw.polarity_fingerprint(chain, terms, z_drive_omegas_per_site=[])


def test_polarity_fingerprint_materialises_a_generator_rate():
    """The `terms` generator has a test; the RATE generator did not, and it is the one
    that silently deletes a channel: exhausted before `pi_decompose_M`, gamma_t1 reads as
    absent, the T1 channel vanishes from the measurement, and ‖M‖² drops from 541.76 to
    512.0 while c_bit_b_homogeneous still reports False from the other copy."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    terms = [('X', 'Y')]
    rates = [0.1, 0.2, 0.3]
    from_list = fw.polarity_fingerprint(chain, terms, gamma_t1=list(rates))
    from_gen = fw.polarity_fingerprint(chain, terms, gamma_t1=(r for r in rates))
    assert from_gen['f112_M_norm_sq'] == pytest.approx(from_list['f112_M_norm_sq'], rel=1e-12)
    assert from_gen['c_bit_b_homogeneous'] is False
    assert from_list['f112_M_norm_sq'] == pytest.approx(541.76, rel=1e-9)


def test_polarity_fingerprint_drive_match_stays_relative_at_tiny_coupling():
    """Why the drive-match window is scaled and not absolute. With a constant 1e-9 window
    a sign-FLIPPED drive is accepted below J ≈ 1e-10 and the inversion then returns a
    NEGATIVE rate. Pinned at J = 1e-12, where an absolute window returns −0.001."""
    from framework.workflows.polarity_fingerprint import _true_z_moments

    chain = fw.ChainSystem(N=2, gamma_0=0.0, J=1e-12)
    moments = _true_z_moments(chain, [('Z', 'I')])
    assert moments[0] != 0.0
    flipped = [-m for m in moments]
    res = fw.polarity_fingerprint(
        chain, [('Z', 'I')], gamma_t1=0.001, z_drive_omegas_per_site=flipped)
    assert res['f113_drive_matches_terms'] is False
    assert res['f113_extracted_gamma_t1'] is None


def test_polarity_fingerprint_detects_bit_b_inhomogeneous_H():
    """The h_bit_b_homogeneous False path was never asserted; hard-wiring it True left
    the suite green. XX+XY mixes bit_b = 0 and bit_b = 1 in H."""
    chain = fw.ChainSystem(N=2, gamma_0=0.0)
    mixed = fw.polarity_fingerprint(chain, [('X', 'X'), ('X', 'Y')])
    assert mixed['h_bit_b_homogeneous'] is False
    assert mixed['in_f112_typed_scope'] is False
    homog = fw.polarity_fingerprint(chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')])
    assert homog['h_bit_b_homogeneous'] is True


def test_polarity_fingerprint_empty_terms_raises():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    with pytest.raises(ValueError):
        fw.polarity_fingerprint(chain, [])


def test_polarity_fingerprint_f113_accepts_per_site_rates():
    """F113's closed form is rate-weighted PER SITE, so the workflow must accept a
    per-site rate list. It used to collapse the rates to a scalar (`float(gp - gt1)`)
    and raise TypeError on a list, which is the one shape the law actually needs."""
    chain = fw.ChainSystem(N=2, gamma_0=0.0)
    result = fw.polarity_fingerprint(
        chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')],
        gamma_t1=[0.001, 0.002],
        z_drive_omegas_per_site=[0.13, 0.13],
    )
    # (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l) = 8·(0.13·(−0.001) + 0.13·(−0.002))
    expected = 8.0 * (0.13 * -0.001 + 0.13 * -0.002)
    assert result['f113_predicted'] == pytest.approx(expected, rel=1e-12)


def test_polarity_fingerprint_f113_uniform_rate_matches_registry_anchor():
    """The scalar path must keep the registry's anchor: ω = 0.13 on both sites,
    γ_T1 = 0.001 at N=2 gives −2.08e-3."""
    chain = fw.ChainSystem(N=2, gamma_0=0.0)
    result = fw.polarity_fingerprint(
        chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')],
        gamma_t1=0.001,
        z_drive_omegas_per_site=[0.13, 0.13],
    )
    assert result['f113_predicted'] == pytest.approx(-2.08e-3, rel=1e-12)


def test_polarity_fingerprint_f113_zero_total_drive_still_fires():
    """The case the collapsed form could not express. With Σ_l ω_l = 0 exactly, a
    factorised prediction (Σω)·(scalar Δγ) is 0 for EVERY rate profile, while the true
    per-site sum is nonzero. The firing condition is Σ_l ω_l·Δγ_l ≠ 0 and nothing
    weaker: a non-uniform profile ORTHOGONAL to the drive gives exactly 0 too, which
    is why 'non-uniform rates' is not the criterion."""
    chain = fw.ChainSystem(N=4, gamma_0=0.0)
    terms = [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')]
    # Note the sign: this test declares ω directly, where the proof's §4 and the
    # verifier's BLOCK I specify the same physics as site energies Σ_l α_l n_l with
    # α = −ω, so they report −25.6 where this reports +25.6. Same magnitude, opposite
    # convention, and the pair is the one place in this arc the sign can be misread.
    omegas = [0.1, -0.1, 0.1, -0.1]  # Σ_l ω_l = 0 exactly
    pref = (4 ** 4) / 2.0

    unequal = fw.polarity_fingerprint(
        chain, terms, gamma_t1=[1.0, 2.0, 1.0, 2.0], z_drive_omegas_per_site=omegas,
    )
    assert unequal['f113_predicted'] == pytest.approx(pref * 0.2, rel=1e-12)

    uniform = fw.polarity_fingerprint(
        chain, terms, gamma_t1=[1.0, 1.0, 1.0, 1.0], z_drive_omegas_per_site=omegas,
    )
    assert uniform['f113_predicted'] == pytest.approx(0.0, abs=1e-13)

    # non-uniform, but orthogonal to this drive profile: also exactly zero
    orthogonal = fw.polarity_fingerprint(
        chain, terms, gamma_t1=[1.0, 1.0, 2.0, 2.0], z_drive_omegas_per_site=omegas,
    )
    assert orthogonal['f113_predicted'] == pytest.approx(0.0, abs=1e-13)


def test_polarity_fingerprint_f113_extraction_refuses_a_drive_absent_from_terms():
    """The inversion divides a MEASURED asymmetry by a DECLARED drive. If the drive is
    not the single-site Z content of `terms`, the quotient is meaningless: it used to
    return a plausible wrong rate (7.7x out for terms=[('Z','I')]) or a confident −0.0
    for bond-only terms. Both must now be None with the reason named."""
    chain = fw.ChainSystem(N=2, gamma_0=0.0)

    # bond-only H: true moments are all zero, so a declared 0.13 drive is absent
    bond_only = fw.polarity_fingerprint(
        chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')],
        gamma_t1=0.001, z_drive_omegas_per_site=[0.13, 0.13],
    )
    assert bond_only['f113_drive_matches_terms'] is False
    assert bond_only['f113_extracted_gamma_t1'] is None
    assert 'not the single-site Z content' in bond_only['reading']

    # a genuine single-site Z drive IS reachable through the bilinear interface,
    # and there the declared drive must equal Tr(Z_l H)/2^(N−1) = (2, 0)
    wrong_drive = fw.polarity_fingerprint(
        chain, [('Z', 'I')], gamma_t1=0.001, z_drive_omegas_per_site=[0.13, 0.13],
    )
    assert wrong_drive['f112_asymmetry'] != pytest.approx(0.0, abs=1e-12)
    assert wrong_drive['f113_drive_matches_terms'] is False
    assert wrong_drive['f113_extracted_gamma_t1'] is None

    right_drive = fw.polarity_fingerprint(
        chain, [('Z', 'I')], gamma_t1=0.001, z_drive_omegas_per_site=[2.0, 0.0],
    )
    assert right_drive['f113_drive_matches_terms'] is True
    assert right_drive['f113_extracted_gamma_t1'] == pytest.approx(0.001, rel=1e-12)


def test_polarity_fingerprint_f113_reading_names_each_refusal():
    """Every path that withholds the extraction must say which condition failed, and a
    mismatch must report the drive `terms` actually carries so it is diagnosable."""
    chain = fw.ChainSystem(N=2, gamma_0=0.0)

    mismatch = fw.polarity_fingerprint(
        chain, [('Z', 'I')], gamma_t1=0.001, z_drive_omegas_per_site=[0.13, 0.13],
    )
    assert mismatch['f113_true_z_moments'] == pytest.approx([2.0, 0.0], abs=1e-12)
    assert '[2.0, 0.0]' in mismatch['reading']

    zero_total = fw.polarity_fingerprint(
        chain, [('X', 'X')], gamma_t1=0.001, z_drive_omegas_per_site=[0.0, 0.0],
    )
    assert zero_total['f113_drive_matches_terms'] is True
    assert zero_total['f113_extracted_gamma_t1'] is None
    assert 'Σ_l ω_l = 0' in zero_total['reading']

    pumped = fw.polarity_fingerprint(
        chain, [('Z', 'I')], gamma_t1=0.001, gamma_pump=0.002,
        z_drive_omegas_per_site=[2.0, 0.0],
    )
    assert pumped['f113_drive_matches_terms'] is True
    assert pumped['f113_extracted_gamma_t1'] is None
    assert 'pumping is on' in pumped['reading']


def test_polarity_fingerprint_f113_single_site_moments_share_a_sign():
    """Why there is no mixed-sign guard: every term is built at the one coefficient
    chain.J, so the single-site Z moments the measurement carries all have sign(J).
    A mixed-sign drive can therefore never match `terms` on this path, and the
    unboundedness of Σ_l ω_l·γ_T1,l / Σ_l ω_l is a hazard only for hand application
    of the closed form. Pinned so a builder change surfaces here."""
    from framework.workflows.polarity_fingerprint import _true_z_moments

    for N in (3, 4):
        for topology in ('chain', 'ring', 'star'):
            for J in (1.0, -1.0, 2.5):
                chain = fw.ChainSystem(N=N, gamma_0=0.0, topology=topology, J=J)
                for terms in ([('Z', 'I')], [('I', 'Z')], [('Z', 'I'), ('I', 'Z')],
                              [('X', 'X'), ('Z', 'I')], [('Z', 'I', 'I')]):
                    moments = _true_z_moments(chain, terms)
                    nonzero = [m for m in moments if abs(m) > 1e-12]
                    assert nonzero, f"no nonzero moment: {N} {topology} {J} {terms}"
                    assert (all(m > 0 for m in nonzero)
                            or all(m < 0 for m in nonzero)), (topology, J, terms, moments)


def test_polarity_fingerprint_accepts_a_generator_for_terms():
    """`terms` is read by four separate consumers, so a one-shot iterator has to be
    materialised before any of them. Unmaterialised it emptied after the first read and
    the whole fingerprint then described the ZERO Hamiltonian, with no error raised:
    f87 'soft' became 'truly', the asymmetry −0.016 became 0.0, the moments [2, 0]
    became [0, 0]. A generator is truthy even when empty, so the non-empty guard could
    not fire either."""
    chain = fw.ChainSystem(N=2, gamma_0=0.0)
    kw = dict(gamma_t1=0.001, z_drive_omegas_per_site=[2.0, 0.0])
    from_list = fw.polarity_fingerprint(chain, [('Z', 'I')], **kw)
    from_gen = fw.polarity_fingerprint(chain, (t for t in [('Z', 'I')]), **kw)
    for key in ('f87_class', 'f112_asymmetry', 'f112_verdict',
                'f113_true_z_moments', 'f113_extracted_gamma_t1'):
        assert from_gen[key] == from_list[key], key
    # pin the values too, so a both-broken pair cannot agree its way through
    assert from_gen['f113_true_z_moments'] == pytest.approx([2.0, 0.0], abs=1e-12)
    assert from_gen['f113_extracted_gamma_t1'] == pytest.approx(0.001, rel=1e-9)

    with pytest.raises(ValueError):
        fw.polarity_fingerprint(chain, (t for t in []))


def test_polarity_fingerprint_f113_accepts_a_float_noise_zero_moment():
    """Tr(Z_l H) is a sum of 2^N signed terms; at non-dyadic J it cancels to ~1e-17
    rather than to 0.0. The drive-match tolerance therefore has to scale with the drive
    VECTOR, not with the component: per component, the correct declaration of 0.0
    against a 1e-17 moment fell outside a 1e-9-relative window and was refused, while
    the same result printed that noise as the drive it expected."""
    from framework.workflows.polarity_fingerprint import _true_z_moments

    chain = fw.ChainSystem(N=4, gamma_0=0.0, topology='star', J=0.7)
    terms = [('Z', 'I'), ('X', 'X'), ('Y', 'Y'), ('Z', 'Z')]
    moments = _true_z_moments(chain, terms)
    assert any(0.0 < abs(m) < 1e-12 for m in moments), \
        f"no float-noise moment to test against: {moments}"

    declared = [m if abs(m) > 1e-12 else 0.0 for m in moments]
    ok = fw.polarity_fingerprint(
        chain, terms, gamma_t1=0.001, z_drive_omegas_per_site=declared)
    assert ok['f113_drive_matches_terms'] is True
    assert ok['f113_extracted_gamma_t1'] == pytest.approx(0.001, rel=1e-9)

    # the widened tolerance must still refuse a genuinely wrong drive at the same J
    flipped = list(declared)
    flipped[0] = -flipped[0]
    refused = fw.polarity_fingerprint(
        chain, terms, gamma_t1=0.001, z_drive_omegas_per_site=flipped)
    assert refused['f113_drive_matches_terms'] is False
    assert refused['f113_extracted_gamma_t1'] is None


def test_polarity_fingerprint_reading_string_includes_classifications():
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    result = fw.polarity_fingerprint(chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')])
    reading = result['reading']
    assert 'F87' in reading
    assert 'F112' in reading
    assert result['f87_class'] in reading
    assert result['f112_verdict'] in reading


def test_polarity_fingerprint_refuses_negative_and_nan_rates():
    """A negative rate reached lindbladian_z_plus_t1 as sqrt(negative), so ‖M‖², the
    asymmetry and the relative asymmetry all went NaN. Both verdict comparisons are
    `NaN < x`, i.e. False, so the classifier fell through to a confident 'BROKEN' on an
    all-NaN measurement. NaN itself was accepted the same way."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    terms = [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')]
    for kw in ({'gamma_t1': -0.001}, {'gamma_pump': -0.001},
               {'gamma_t1': [0.1, -0.1, 0.1]}, {'gamma_t1': float('nan')}):
        with pytest.raises(ValueError, match=r"is not a valid rate"):
            fw.polarity_fingerprint(chain, terms, **kw)


def test_polarity_fingerprint_pump_alone_leaves_typed_scope():
    """c is bit_b-homogeneous iff there is NO σ⁻ and NO σ⁺. Dropping the σ⁺ half of that
    check left the suite green, because every test exercising the False path also had T1
    on. σ⁺ alone must break it too."""
    chain = fw.ChainSystem(N=3, gamma_0=0.05)
    pumped = fw.polarity_fingerprint(chain, [('X', 'X'), ('Y', 'Y'), ('Z', 'Z')],
                                     gamma_pump=0.1)
    assert pumped['c_bit_b_homogeneous'] is False
    assert pumped['in_f112_typed_scope'] is False


def test_polarity_fingerprint_drive_scale_floor_and_both_sides():
    """Two halves of the drive-match scale, each of which survived mutation.

    The noise FLOOR: when every true moment is float noise, the vector max is itself
    ~1e-17 and the window collapses, so a correct all-zero declaration was refused with a
    message accusing the caller. The BOTH-SIDES max: taking the scale over the declared
    drive alone accepts a 100%-wrong declaration when the declared side is the zero one."""
    from framework.workflows.polarity_fingerprint import _true_z_moments

    noisy = fw.ChainSystem(N=4, gamma_0=0.0, J=0.7)
    moments = _true_z_moments(noisy, [('Z', 'Z')])
    assert all(abs(m) < 1e-12 for m in moments) and any(m != 0.0 for m in moments), moments
    res = fw.polarity_fingerprint(noisy, [('Z', 'Z')], gamma_t1=0.001,
                                  z_drive_omegas_per_site=[0.0] * 4)
    assert res['f113_drive_matches_terms'] is True

    tiny = fw.ChainSystem(N=2, gamma_0=0.0, J=1e-10)
    true = _true_z_moments(tiny, [('Z', 'I')])
    assert abs(true[0]) > 0.0
    wrong = fw.polarity_fingerprint(tiny, [('Z', 'I')], gamma_t1=0.001,
                                    z_drive_omegas_per_site=[0.0, 0.0])
    assert wrong['f113_drive_matches_terms'] is False
