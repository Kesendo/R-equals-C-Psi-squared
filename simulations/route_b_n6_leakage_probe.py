"""Numerical Route-B N6 two-plane leakage and full-block bond probe.

Exploration at atlas seeds, not algebraic certification. Units: gamma=1,
q=J/gamma in the octic book, t=i*q, H=J sum(XX+YY), delta=0.
"""
from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np
import scipy.linalg as la
import scipy
from scipy.optimize import root
from framework.weight_coherence_block import weight_block_build, weight_block_configs

ROOT = Path(__file__).resolve().parents[1]


def parts():
    d = weight_block_build(6, 1, 2, 0)
    c = weight_block_build(6, 1, 2, 1, gamma=0)
    # Isolate the existing primitive's entries for the end bond, sites 0 and 1.
    basis = [(a, b) for a in weight_block_configs(6, 1)
             for b in weight_block_configs(6, 2)]
    end = np.zeros_like(c)
    for i, (a, b) in enumerate(basis):
        for j, (aa, bb) in enumerate(basis):
            if (a == aa and b ^ bb == 3) or (b == bb and a ^ aa == 3):
                end[i, j] = c[i, j]
    return d, c, end


def seed(row, key):
    return complex(float(row[key]['real']), float(row[key]['imag']))


def choose(d, c):
    source = json.loads((ROOT / 'simulations/results/route_b_a2_n6.json').read_text())
    candidates = {'off_axis': [], 'hermitian': []}
    for row in source['loci']:
        q = seed(row, 'qPhysicalCSharpSeed')
        if row['parity'] != 'E' or q.real < 0 or q.imag < 0:
            continue
        lam = seed(row, 'lambdaPhysicalSeed')
        distances = np.sort(abs(la.eigvals(d + q*c) - lam))
        candidates['hermitian' if q.real == 0 else 'off_axis'].append(
            (float(distances[2]), row, distances[:3].tolist()))
    return {kind: max(rows, key=lambda item: item[0]) for kind, rows in candidates.items()}


def plane(l, center, radius):
    t, v, count = la.schur(l, output='complex', sort=lambda z: abs(z-center) < radius)
    if count != 2:
        raise ValueError(f'Expected isolated pair; found {count}')
    return v[:, :2], t[:2, :2]


def leakage(d, h, v):
    def outside(a):
        av = a @ v
        return av-v@(v.conj().T@av)
    x, y = outside(d), outside(h)
    return {'D': float(la.norm(x)), 'qC': float(la.norm(y)),
            'sum': float(la.norm(x+y)),
            'D_singular_values': la.svdvals(x).tolist(),
            'qC_singular_values': la.svdvals(y).tolist(),
            'D_compression': encode(v.conj().T@d@v),
            'qC_compression': encode(v.conj().T@h@v)}


def contour(l, center, radius, nodes):
    p = np.zeros_like(l)
    identity = np.eye(len(l))
    for theta in 2*np.pi*(np.arange(nodes)+0.5)/nodes:
        z = radius*np.exp(1j*theta)
        p += z*la.solve((center+z)*identity-l, identity)/nodes
    u, s, _ = la.svd(p)
    return u[:, :2], {'projector_singular_values': s[:4].tolist(),
                     'idempotency': float(la.norm(p@p-p))}


def encode(x):
    a = np.asarray(x)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def independent_hop(epsilon):
    """Pauli XY action on all 64 spin states, then restrict each Hilbert leg."""
    h = np.zeros((64, 64))
    for s in range(5):
        for a in range(64):
            if ((a >> s) & 1) != ((a >> (s+1)) & 1):
                h[a ^ (3 << s), a] = 2*(1+epsilon if s == 0 else 1)
    a, b = weight_block_configs(6, 1), weight_block_configs(6, 2)
    return -1j*(np.kron(h[np.ix_(a, a)], np.eye(15))
                - np.kron(np.eye(6), h[np.ix_(b, b)]))


def discriminant(m):
    return (m[0, 0]-m[1, 1])**2+4*m[0, 1]*m[1, 0]


def local_signature(d, hop, q, center, isolation, circle_radius=1e-4):
    """Gap exponents and discriminant winding; no threshold-based EP verdict."""
    radii = circle_radius*np.array([1, 0.1, 0.01, 0.001])
    gaps = []
    for radius in radii:
        _, m = plane(d+(q+radius)*hop, center, isolation*0.25)
        gaps.append(float(abs(discriminant(m))**0.5))
    values = []
    for theta in np.linspace(0, 2*np.pi, 65):
        _, m = plane(d+(q+circle_radius*np.exp(1j*theta))*hop,
                     center, isolation*0.25)
        values.append(discriminant(m))
    angles = np.unwrap(np.angle(values))
    return {'q_radii': radii.tolist(), 'gaps': gaps,
            'gap_exponents': (np.diff(np.log(gaps))/np.diff(np.log(radii))).tolist(),
            'discriminant_winding': float((angles[-1]-angles[0])/(2*np.pi)),
            'maximum_phase_step': float(max(abs(np.diff(angles)))),
            'minimum_discriminant_on_loop': float(min(abs(np.array(values))))}


def continuation(d, c, end, q0, lam, isolation, v0):
    """Track two roots of the local pair discriminant in the complex q plane.

    A zero discriminant is a numerical candidate only. The traceless Schur
    compression and singular values distinguish splitting from scalar crossing.
    The same fixed disk is checked to contain precisely two FULL-block roots.
    """
    radius = isolation*0.4
    # Complex symmetric L: transpose pairing supplies the left dual on the plane.
    left = la.solve(v0.T@v0, v0.T)
    a, b = left@c@v0, left@(q0*end)@v0
    f0, f1, fm = discriminant(b), discriminant(a+b), discriminant(b-a)
    coefficients = [(f1+fm)/2-f0, (f1-fm)/2, f0]
    slopes = np.roots(coefficients)
    results = []
    for direction, slope in enumerate(slopes):
        previous, previous_epsilon = q0, 0.0
        for epsilon in (0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.05):
            hop = c+epsilon*end
            start = previous+slope*(epsilon-previous_epsilon)
            def evaluate(q):
                return plane(d+q*hop, lam, radius)
            def objective(x):
                _, m = evaluate(complex(*x))
                value = discriminant(m)/epsilon**2
                return [value.real, value.imag]
            solved = root(objective, [start.real, start.imag], tol=1e-9)
            q = complex(*solved.x)
            v, m = evaluate(q)
            center = np.trace(m)/2
            departure = la.norm(m-center*np.eye(2))
            values = la.eigvals(d+q*hop)
            ds = np.sort(abs(values-center))
            fixed_v, fixed_m = evaluate(q0)
            sv = la.svdvals(d+q*hop-center*np.eye(90))
            results.append({'branch': direction, 'epsilon': epsilon, 'q': encode(q),
                            'lambda': encode(center), 'solver_success': bool(solved.success),
                            'scaled_discriminant_residual': float(la.norm(objective(solved.x))),
                            'eigenvalue_gap': float(abs(la.eigvals(m)[0]-la.eigvals(m)[1])),
                            'fixed_q_gap': float(abs(la.eigvals(fixed_m)[0]-la.eigvals(fixed_m)[1])),
                            'traceless_compression_norm': float(departure),
                            'smallest_singular_values': sv[-3:].tolist(),
                            'third_distance': float(ds[2]),
                            'initial_plane_overlap': la.svdvals(v0.conj().T@v).tolist(),
                            'invariance_residual': float(la.norm((d+q*hop)@v-v@m)),
                            'leakage': leakage(d, q*hop, v)})
            if epsilon in (0.01, 0.05):
                results[-1]['local_signature'] = local_signature(d, hop, q, center, ds[2])
                results[-1]['contours'] = []
                for fraction in (0.25, 0.5, 0.75):
                    vc, diagnostics = contour(d+q*hop, center, fraction*ds[2], 128)
                    mc = vc.conj().T@(d+q*hop)@vc
                    results[-1]['contours'].append({'fraction': fraction, **diagnostics,
                        'traceless_compression_norm': float(la.norm(mc-np.trace(mc)/2*np.eye(2))),
                        'plane_distance': float(la.norm(vc@vc.conj().T-v@v.conj().T))})
            slope = (q-previous)/(epsilon-previous_epsilon)
            previous, previous_epsilon = q, epsilon
    return {'linear_slopes': encode(slopes), 'disk_center': encode(lam),
            'disk_radius': radius, 'rows': results}


def main():
    d, c, end = parts()
    report = {'conventions': __doc__, 'dimension': len(d), 'selection':
              'Largest full-90D third-eigenvalue distance among E seeds with Re q, Im q >= 0',
              'probes': {}}
    source = ROOT / 'simulations/results/route_b_a2_n6.json'
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    atlas = json.loads((ROOT / 'simulations/results/route_b_a2_n6_atlas.json').read_text())
    if digest != atlas['sourceSha256']:
        raise ValueError('Atlas and seed source disagree')
    report['provenance'] = {'source_sha256': digest,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'numpy_version': np.__version__, 'scipy_version': scipy.__version__}
    report['builder_checks'] = {str(e): float(la.norm(c+e*end-independent_hop(e)))
                                for e in (0, 0.125, -0.25)}
    if any(report['builder_checks'].values()):
        raise ValueError('Weighted builder disagrees with full Hilbert-space construction')
    report['missing_bond_mutation_residual'] = float(la.norm(c-independent_hop(0.125)))
    if report['missing_bond_mutation_residual'] == 0:
        raise ValueError('Independent construction failed to see the bond mutation')
    for kind, (gap, row, distances) in choose(d, c).items():
        q, lam = seed(row, 'qPhysicalCSharpSeed'), seed(row, 'lambdaPhysicalSeed')
        l = d+q*c
        v, _ = plane(l, lam, gap/2)
        readings = []
        for fraction in (0.25, 0.5, 0.75):
            for nodes in (64, 128):
                vc, diagnostics = contour(l, lam, fraction*gap, nodes)
                readings.append({'fraction': fraction, 'nodes': nodes, **diagnostics,
                                 'plane_distance': float(la.norm(vc@vc.conj().T-v@v.conj().T)),
                                 **leakage(d, q*c, vc)})
        # Wrong plane: replace one direction with a deterministic external direction.
        e = np.eye(len(d), dtype=complex)[:, 0]
        e -= v@(v.conj().T@e)
        e /= la.norm(e)
        wrong = np.column_stack((v[:, 0], e))
        probe = {'id': row['id'], 'q': encode(q), 'lambda': encode(lam),
                 'nearest_distances': distances, 'schur': leakage(d, q*c, v),
                 'contours': readings, 'wrong_plane': leakage(d, q*c, wrong)}
        if kind == 'off_axis':
            dc = v.conj().T@d@v
            scalar = np.trace(dc).real/2
            outside = d@v-v@dc
            predicted = (scalar+6)*(-2-scalar)
            probe['two_rung_variance'] = {'d': float(scalar),
                'weight_n_diff_1': float((scalar+6)/4),
                'weight_n_diff_3': float((-2-scalar)/4),
                'predicted_leakage_singular_value': float(np.sqrt(predicted)),
                'gram_residual': float(la.norm(outside.conj().T@outside-predicted*np.eye(2)))}
            probe['local_signature'] = local_signature(d, c, q, lam, gap)
            probe['continuation'] = continuation(d, c, end, q, lam, gap, v)
        report['probes'][kind] = probe
        print(kind, row['id'], q, lam, 'isolation', gap, 'leakage', probe['schur'], flush=True)
    output = ROOT / 'simulations/results/route_b_n6_leakage_probe.json'
    output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(output)


if __name__ == '__main__':
    main()
