"""Finite product-preparation/local-readout search, ideal equal-shot mean SNR.

Not a global measurement optimum, Fisher information, or EP certificate.
Run: OPENBLAS_NUM_THREADS=1 python simulations/route_b_n4_readout_search.py
"""
import json
import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import numpy as np
import scipy.linalg as la
from route_b_n4_virtual_readout import (
    ROOT, EP, TIMES, generator, verify_generator_contract,
    artifact_provenance, verify_artifact_provenance,
)


def preparations():
    labels, columns = [], []
    for excited in range(4):
        for plus in range(4):
            if plus == excited:
                continue
            for phase in [0, 1]:
                a, b = 1 << excited, (1 << excited) | (1 << plus)
                psi = np.zeros(16, complex)
                psi[a], psi[b] = 1/np.sqrt(2), (1j)**phase/np.sqrt(2)
                labels.append(dict(excited_site=excited, superposition_site=plus,
                                   phase_pi_over_2=phase))
                columns.append(np.outer(psi, psi.conj()).ravel())
    return labels, np.array(columns).T


def measurements():
    labels, means, seconds = [], [], []
    for site in range(4):
        edges = [(a, a | (1 << site)) for a in range(16) if not a & (1 << site)]
        for quadrature in ['X', 'Y']:
            for selected in [None]+list(range(8)):
                chosen = edges if selected is None else [edges[selected]]
                mean, second = np.zeros(256, complex), np.zeros(256)
                for a, b in chosen:
                    mean[16*a+b] = 1 if quadrature == 'X' else 1j
                    mean[16*b+a] = 1 if quadrature == 'X' else -1j
                    second[16*a+a] = second[16*b+b] = 1
                labels.append(dict(site=site, quadrature=quadrature,
                                   spectator_mask=None if selected is None else edges[selected][0]))
                means.append(mean); seconds.append(second)
    return labels, np.array(means), np.array(seconds)


def signals(epsilon, initial, means, seconds):
    g = generator(epsilon)
    step = la.expm(g*(TIMES[1]-TIMES[0]))
    state = initial.copy()
    values, moments = [], []
    for index, time in enumerate(TIMES):
        value, moment = means@state, seconds@state
        assert np.max(abs(value.imag)) < 2e-11
        assert np.max(abs(moment.imag)) < 2e-11
        assert np.min(moment.real-abs(value.real)) > -2e-11
        assert np.max(moment.real) < 1+2e-11
        values.append(value.real); moments.append(moment.real)
        if index in [75, 400]:
            assert la.norm(state-la.expm(g*time)@initial) < 2e-11
        state = step@state
    return np.array(values), np.array(moments)


def scores(first, second):
    signal, moment = first
    other, other_moment = second
    variance = moment-signal**2+other_moment-other**2
    result = abs(other-signal)/np.sqrt(np.maximum(variance, 1e-28))
    result[0] = 0
    return result


def main():
    verify_generator_contract(builder=generator)
    prep, initial = preparations()
    read, means, seconds = measurements()
    for pi, label in enumerate(prep):
        r = read.index(dict(site=label['superposition_site'],
                            quadrature='Y' if label['phase_pi_over_2'] else 'X', spectator_mask=None))
        assert abs(means[r]@initial[:, pi]-1) < 1e-14
    baseline_p = prep.index(dict(excited_site=0, superposition_site=1, phase_pi_over_2=0))
    baseline_r = read.index(dict(site=1, quadrature='X', spectator_mask=1))
    rows, stored = [], {}
    for offset in [-.2, 0., .2]:
        center = signals(EP+offset, initial, means, seconds)
        for shift in [-.05, .05]:
            other = signals(EP+offset+shift, initial, means, seconds)
            score = scores(center, other)
            assert np.max(scores(center, center)) == 0
            ti, ri, pi = np.unravel_index(np.argmax(score), score.shape)
            canonical_p = prep.index(dict(excited_site=2, superposition_site=1, phase_pi_over_2=1))
            canonical_r = read.index(dict(site=2, quadrature='X', spectator_mask=2))
            reflected_p = prep.index(dict(excited_site=1, superposition_site=2, phase_pi_over_2=1))
            reflected_r = read.index(dict(site=1, quadrature='X', spectator_mask=4))
            assert np.max(abs(score[:, canonical_r, canonical_p]-score[:, reflected_r, reflected_p])) < 2e-12
            assert abs(np.max(score[:, canonical_r, canonical_p])-np.max(score)) < 2e-12
            baseline = float(np.max(score[:, baseline_r, baseline_p]))
            row = dict(center_offset=offset, epsilon_shift=shift,
                       time=float(TIMES[ti]), preparation=prep[pi], readout=read[ri],
                       signal=float(center[0][ti, ri, pi]),
                       difference=float(other[0][ti, ri, pi]-center[0][ti, ri, pi]),
                       shots_each=float(25/score[ti, ri, pi]**2),
                       baseline_shots_each=float(25/baseline**2),
                       shot_reduction=float((score[ti, ri, pi]/baseline)**2))
            rows.append(row)
            stored[offset, shift] = (score, pi, ri)
            print(json.dumps(row), flush=True)
    # Transfer the EP-selected protocol to both off-EP centers. Reoptimize
    # time only, and also report exactly the EP time with everything fixed.
    transfers = []
    for shift in [-.05, .05]:
        ep_score, pi, ri = stored[0., shift]
        ep_t = int(np.argmax(ep_score[:, ri, pi]))
        for offset in [-.2, 0., .2]:
            score = stored[offset, shift][0][:, ri, pi]
            best = int(np.argmax(score))
            transfers.append(dict(epsilon_shift=shift, center_offset=offset,
                                  best_time=float(TIMES[best]),
                                  shots_each_best_time=float(25/score[best]**2),
                                  shots_each_fixed_ep_time=float(25/score[ep_t]**2)))
    # Reproduce the independently propagated baseline, not a copied number.
    from route_b_n4_virtual_readout import run
    _, original, original_second, _ = run(EP)
    check = signals(EP, initial[:, [baseline_p]], means[[baseline_r]], seconds[[baseline_r]])
    assert np.max(abs(check[0][:, 0, 0]-original)) < 2e-12
    assert np.max(abs(check[1][:, 0, 0]-original_second)) < 2e-12
    # Physical winner's selected Jordan overlap, not inferred from the score.
    from route_b_other_n_unfolding import parts
    d, c, left, right, _ = parts(4)
    block = d+2*(c+EP*(left+right)/2)
    pairs = [(a, b) for a in range(16) if a.bit_count() == 1
             for b in range(16) if b.bit_count() == 2]
    indices = [16*a+b for a, b in pairs]
    lam = -4+2j
    p = np.zeros_like(block)
    for theta in 2*np.pi*(np.arange(192)+.5)/192:
        z = .2*np.exp(1j*theta)
        p += z*la.solve((lam+z)*np.eye(24)-block, np.eye(24))/192
    pi = prep.index(dict(excited_site=2, superposition_site=1, phase_pi_over_2=1))
    ri = read.index(dict(site=2, quadrature='X', spectator_mask=2))
    a = means[ri, indices]@p@initial[indices, pi]
    b = means[ri, indices]@(block-lam*np.eye(24))@p@initial[indices, pi]
    nil = (block-lam*np.eye(24))@p
    assert la.norm(p@p-p) < 2e-11 and abs(np.trace(p)-2) < 2e-11
    assert la.norm(nil@nil) < 2e-11 and abs(b) > .04
    assert la.norm(la.expm(block*.345)@p-np.exp(lam*.345)*(p+.345*nil)) < 2e-11
    jordan = dict(a=[float(a.real), float(a.imag)], b=[float(b.real), float(b.imag)])
    output = dict(preparations=len(prep), readouts=len(read), times=len(TIMES),
                  q=2, gamma=1, rows=rows, transferred_ep_protocols=transfers,
                  canonical_winner_jordan=jordan,
                  scope='Finite ideal mean-SNR search; M each setting, 2M total; design search is numerical, not experimental shots.',
                  **artifact_provenance(__file__))
    verify_artifact_provenance(output, __file__)
    (ROOT/'simulations/results/route_b_n4_readout_search.json').write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
    print('ALL PASS')


if __name__ == '__main__':
    main()
