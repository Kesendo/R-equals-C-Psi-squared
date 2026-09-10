"""Covariance-weighted histogram contrast with the fixed N4 product protocol."""
import json
import os

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import expm_multiply
from route_b_n4_virtual_readout import (
    ROOT, EP, TIMES, generator, verify_generator_contract,
    artifact_provenance, verify_artifact_provenance,
)


def record(epsilon):
    psi = np.zeros(16, complex); psi[4] = 1/np.sqrt(2); psi[6] = 1j/np.sqrt(2)
    initial = np.outer(psi, psi.conj()).ravel()
    states = expm_multiply(generator(epsilon), initial, start=0, stop=2, num=len(TIMES)).reshape(-1, 16, 16)
    edges = [(a, a | 4) for a in range(16) if not a & 4]
    outcomes = []
    for a, b in edges:
        total = (states[:, a, a]+states[:, b, b]).real
        mean = 2*states[:, a, b].real
        outcomes.extend([(total+mean)/2, (total-mean)/2])
    p = np.array(outcomes).T
    assert p.min() > -2e-12 and np.max(abs(p.sum(axis=1)-1)) < 2e-12
    # Independently rotate into the actual local measurement basis.
    h = np.eye(16)
    for a, b in edges:
        h[np.ix_([a, b], [a, b])] = np.array([[1, 1], [1, -1]])/np.sqrt(2)
    for k in [0, 69, 200]:
        rotated = h@states[k]@h.T
        assert max(abs(p[k]-np.array([rotated[a, a].real for edge in edges for a in edge]))) < 2e-12
    return np.maximum(p, 0)


def compress(p, kind):
    # Edge a=2,b=6 is the previously selected conditional readout.
    if kind == 'ternary':
        return np.column_stack([p[:, 4], p[:, 5], p.sum(axis=1)-p[:, 4]-p[:, 5]])
    if kind == 'spectators':
        return p.reshape(-1, 8, 2).sum(axis=2)
    return p


def matched(p, q):
    s, d = p+q, q-p
    weights = np.divide(d, s, out=np.zeros_like(d), where=s > 0)
    g = (d*weights).sum(axis=1)
    score2 = np.divide(g, 1-g/2, out=np.full_like(g, np.inf), where=g < 2)
    return score2, weights


def main():
    verify_generator_contract(builder=generator)
    rows = []
    for offset in [-.2, 0., .2]:
        p = record(EP+offset)
        for shift in [-.05, .05]:
            q = record(EP+offset+shift)
            row = dict(center_offset=offset, shift=shift, methods={})
            pt, qt = compress(p, 'ternary'), compress(q, 'ternary')
            old = np.array([1., -1., 0.])
            difference = (qt-pt)@old
            variance = (pt+qt)@(old**2)-(pt@old)**2-(qt@old)**2
            old_score = difference**2/np.maximum(variance, 1e-28); old_score[0] = 0
            all_scores = {'original_mean': old_score}
            for kind in ['ternary', 'spectators', 'full']:
                pp, qq = compress(p, kind), compress(q, kind)
                score, weights = matched(pp, qq)
                assert np.max(matched(pp, pp)[0]) == 0
                for k in [69, 200]:
                    cov = np.diag(pp[k]+qq[k])-np.outer(pp[k], pp[k])-np.outer(qq[k], qq[k])
                    d = qq[k]-pp[k]
                    # Total probability is fixed: its covariance eigenvalue
                    # is zero, even when floating propagation leaves ~1e-15.
                    assert abs(d@la.pinvh(cov, atol=1e-12)@d-score[k]) < 2e-12
                    w = weights[k]
                    assert abs((d@w)**2/(w@cov@w)-score[k]) < 2e-12
                all_scores[kind] = score
            assert np.min(all_scores['ternary']-old_score) > -2e-12
            assert np.min(all_scores['full']-all_scores['ternary']) > -2e-12
            assert np.min(all_scores['full']-all_scores['spectators']) > -2e-12
            fixed = int(np.argmax(old_score))
            for name, score in all_scores.items():
                best = int(np.argmax(score))
                row['methods'][name] = dict(best_time=float(TIMES[best]),
                    shots_each=float(25/score[best]),
                    shots_each_at_original_time=float(25/score[fixed]))
            rows.append(row)
            print(json.dumps(row), flush=True)
    # A contrast wholly hidden by pooling must be recovered by full outcomes.
    p = np.array([[.25, .25, .25, .25]])
    q = np.array([[.3, .2, .2, .3]])
    assert matched(p, q)[0][0] > .02
    assert matched(p.reshape(1, 2, 2).sum(2), q.reshape(1, 2, 2).sum(2))[0][0] == 0
    assert np.isinf(matched(np.array([[1., 0]]), np.array([[0., 1]]))[0][0])
    output = dict(protocol='mask4+i mask6, X site2 and Z spectators; fixed across all rows',
        accounting='M shots each, 2M total, known-model weights; mean SNR5 not test power', rows=rows,
        **artifact_provenance(__file__))
    verify_artifact_provenance(output, __file__)
    (ROOT/'simulations/results/route_b_n4_histogram_filter.json').write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
    print('ALL PASS')


if __name__ == '__main__':
    main()
