"""Compare end-bond profiles at the same Route-B N6 crossing in full 90D.

Run: python simulations/route_b_n6_end_profile_probe.py
Exploratory double-precision analytic continuation, not an exact EP certificate.
gamma=1, delta=0, q=qCSharp, H=q sum j_s(XX+YY), t=i*q.
Profiles (left end, right end): one=(1+e,1), even=(1+e/2,1+e/2),
odd=(1+e/2,1-e/2). All three internal bonds retain strength 1.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
import scipy.linalg as la
from scipy.optimize import root

import route_b_n6_leakage_probe as base

ROOT = Path(__file__).resolve().parents[1]
EPSILONS = (0.002, 0.003, 0.005, 0.01, 0.02, 0.03, 0.05)
DIRECTIONS = {'one': (1., 0.), 'even': (0.5, 0.5), 'odd': (0.5, -0.5)}


def reflected_end(left):
    basis = [(a, b) for a in base.weight_block_configs(6, 1)
             for b in base.weight_block_configs(6, 2)]
    index = {pair: i for i, pair in enumerate(basis)}
    reverse = lambda x: int(f'{x:06b}'[::-1], 2)
    permutation = [index[(reverse(a), reverse(b))] for a, b in basis]
    return left[np.ix_(permutation, permutation)], permutation


def weighted_hop(weights):
    """Independent full Hilbert-space XY construction followed by leg restriction."""
    h = np.zeros((64, 64))
    for s, weight in enumerate(weights):
        for a in range(64):
            if ((a >> s) & 1) != ((a >> (s+1)) & 1):
                h[a ^ (3 << s), a] = 2*weight
    a, b = base.weight_block_configs(6, 1), base.weight_block_configs(6, 2)
    return -1j*(np.kron(h[np.ix_(a, a)], np.eye(15))
                - np.kron(np.eye(6), h[np.ix_(b, b)]))


def leading_slopes(a, b):
    f0 = base.discriminant(b)
    fp, fm = base.discriminant(a+b), base.discriminant(b-a)
    coefficients = [(fp+fm)/2-f0, (fp-fm)/2, f0]
    return np.roots(coefficients)


def measure(d, hop, q, center, radius, initial):
    l = d+q*hop
    v, m = base.plane(l, center, radius)
    lam = np.trace(m)/2
    traceless = m-lam*np.eye(2)
    distances = np.sort(abs(la.eigvals(l)-lam))
    return v, m, {'q': base.encode(q), 'lambda': base.encode(lam),
        'eigenvalue_gap': float(abs(base.discriminant(m))**0.5),
        'discriminant_magnitude': float(abs(base.discriminant(m))),
        'traceless_compression_norm': float(la.norm(traceless)),
        'smallest_singular_values': la.svdvals(l-lam*np.eye(90))[-3:].tolist(),
        'third_distance': float(distances[2]),
        'invariance_residual': float(la.norm(l@v-v@m)),
        'initial_plane_overlap': la.svdvals(initial.conj().T@v).tolist()}


def trace_profile(d, c, perturbation, q0, lam0, isolation, v0, order, slopes):
    rows = []
    radius = isolation*0.4
    epsilons = ((0.0001, 0.0003, 0.001)+EPSILONS if order == 1 else EPSILONS)
    for branch, slope in enumerate(slopes):
        for epsilon in epsilons:
            scale = epsilon**order
            hop = c+epsilon*perturbation
            def objective(x):
                q = q0+scale*complex(*x)
                _, m = base.plane(d+q*hop, lam0, radius)
                value = base.discriminant(m)/scale**2
                return np.array([value.real, value.imag])
            def jacobian(x):
                # Differentiate in scaled displacement, not tiny absolute q steps.
                h = 1e-3
                return np.column_stack([(objective(x+h*np.eye(2)[i])
                    - objective(x-h*np.eye(2)[i]))/(2*h) for i in range(2)])
            solution = root(objective, [slope.real, slope.imag], jac=jacobian, tol=1e-7)
            if not solution.success or not np.isfinite(solution.x).all():
                raise ValueError(f'Continuation failed: branch={branch}, e={epsilon}: {solution.message}')
            slope = complex(*solution.x)
            q = q0+scale*slope
            v, m, row = measure(d, hop, q, lam0, radius, v0)
            _, _, fixed = measure(d, hop, q0, lam0, radius, v0)
            row.update({'epsilon': epsilon, 'branch': branch,
                'solver_success': bool(solution.success),
                'scaled_discriminant_residual': float(la.norm(objective(solution.x))),
                'scaled_q_displacement': base.encode(slope),
                'fixed_q_gap': fixed['eigenvalue_gap'],
                'leakage': base.leakage(d, q*hop, v)})
            rows.append(row)
    for epsilon in epsilons:
        pair = [r for r in rows if r['epsilon'] == epsilon]
        qs = [complex(r['q']['real'], r['q']['imag']) for r in pair]
        separation = abs(qs[0]-qs[1])
        if separation == 0 or not np.isfinite(separation):
            raise ValueError('Two predictors converged to the same point')
        for row, q in zip(pair, qs):
            row['other_q_branch_distance'] = float(separation)
            if epsilon in (0.01, 0.05):
                hop = c+epsilon*perturbation
                lam = complex(row['lambda']['real'], row['lambda']['imag'])
                v, _ = base.plane(d+q*hop, lam, row['third_distance']/2)
                row['local_signature'] = base.local_signature(d, hop, q, lam,
                    row['third_distance'], circle_radius=min(1e-4, separation/10))
                row['contours'] = []
                for fraction in (0.25, 0.5, 0.75):
                    vc, diagnostics = base.contour(d+q*hop, lam,
                        fraction*row['third_distance'], 128)
                    mc = vc.conj().T@(d+q*hop)@vc
                    row['contours'].append({'fraction': fraction, **diagnostics,
                        'traceless_compression_norm': float(la.norm(mc-np.trace(mc)/2*np.eye(2))),
                        'plane_distance': float(la.norm(vc@vc.conj().T-v@v.conj().T))})
    return {'order': order, 'leading_slopes': base.encode(slopes),
            'lambda_disk_center': base.encode(lam0), 'lambda_disk_radius': radius, 'rows': rows}


def scaling_summary(rows):
    first = [r for r in rows if r['branch'] == 0]
    epsilon = np.array([r['epsilon'] for r in first])
    separation = np.array([r['other_q_branch_distance'] for r in first])
    gap = np.array([r['fixed_q_gap'] for r in first])
    return {'epsilon': epsilon.tolist(), 'q_branch_separation': separation.tolist(),
        'separation_exponents': (np.diff(np.log(separation))/np.diff(np.log(epsilon))).tolist(),
        'fixed_q_gap_exponents': (np.diff(np.log(gap))/np.diff(np.log(epsilon))).tolist()}


def main():
    d, c, el = base.parts()
    er, reflection = reflected_end(el)
    source_path = ROOT / 'simulations/results/route_b_a2_n6.json'
    digest = hashlib.sha256(source_path.read_bytes()).hexdigest()
    atlas = json.loads((ROOT / 'simulations/results/route_b_a2_n6_atlas.json').read_text())
    if digest != atlas['sourceSha256']:
        raise ValueError('Atlas source digest mismatch')
    source = json.loads(source_path.read_text())
    seed = next(r for r in source['loci'] if r['id'] == 'N6-E-A2-T-007')
    q0, lam0 = base.seed(seed, 'qPhysicalCSharpSeed'), base.seed(seed, 'lambdaPhysicalSeed')
    l0 = d+q0*c
    isolation = float(np.sort(abs(la.eigvals(l0)-lam0))[2])
    v0, _ = base.plane(l0, lam0, isolation/2)
    dual = la.solve(v0.T@v0, v0.T)
    p = v0@dual
    complement = np.eye(90)-p
    # Reduced resolvent is inverse of lambda-L only on the complementary space.
    resolvent = complement@la.solve(lam0*np.eye(90)-l0+p, complement)
    a = dual@c@v0
    report = {'conventions': __doc__, 'seed_id': seed['id'], 'q0': base.encode(q0),
        'lambda0': base.encode(lam0), 'dimension': 90, 'initial_third_distance': isolation,
        'provenance': {'source_sha256': digest,
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'base_script_sha256': hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest(),
            'numpy_version': np.__version__, 'scipy_version': scipy.__version__},
        'initial_reflection_even_residual': float(la.norm(v0[reflection]-v0)),
        'reduced_resolvent_residual': float(la.norm((lam0*np.eye(90)-l0)@resolvent-complement)),
        'profiles': {}}
    for name, (left, right) in DIRECTIONS.items():
        perturbation = left*el+right*er
        b1 = dual@(q0*perturbation)@v0
        b2 = dual@(q0*perturbation)@resolvent@(q0*perturbation)@v0
        order = 2 if name == 'odd' else 1
        slopes = leading_slopes(a, b2 if order == 2 else b1)
        checks = {}
        for epsilon in (0.125, -0.25):
            weights = [1+epsilon*left, 1, 1, 1, 1+epsilon*right]
            checks[str(epsilon)] = float(la.norm(c+epsilon*perturbation-weighted_hop(weights)))
        if any(checks.values()):
            raise ValueError('Independent weighted Hamiltonian disagrees')
        character = -1 if name == 'odd' else 1
        parity_residual = float(la.norm(perturbation[np.ix_(reflection, reflection)]
                                       -character*perturbation))
        if name != 'one' and parity_residual != 0:
            raise ValueError('End perturbation has incorrect reflection character')
        result = trace_profile(d, c, perturbation, q0, lam0, isolation, v0, order, slopes)
        result.update({'end_direction': [left, right], 'builder_checks': checks,
            'reflection_even_residual': float(la.norm(perturbation[np.ix_(reflection, reflection)]-perturbation)),
            'reflection_odd_residual': float(la.norm(perturbation[np.ix_(reflection, reflection)]+perturbation)),
            'first_order_effective': base.encode(b1), 'second_order_effective': base.encode(b2)})
        result['scaling'] = scaling_summary(result['rows'])
        result['epsilon_sign_reflection_residual'] = float(la.norm(
            (d+q0*(c+0.01*perturbation))[np.ix_(reflection, reflection)]
            -(d+q0*(c-0.01*perturbation))))
        if name == 'odd':
            excursion = q0*perturbation@v0
            intermediate = resolvent@excursion
            result['opposite_parity_excursion'] = {
                'norm': float(la.norm(excursion)),
                'odd_parity_residual': float(la.norm(excursion[reflection]+excursion)),
                'resolvent_odd_parity_residual': float(la.norm(intermediate[reflection]+intermediate))}
        report['profiles'][name] = result
        print(name, 'leading slopes', slopes, flush=True)
        for row in result['rows']:
            if row['epsilon'] in (0.01, 0.05):
                print(row['branch'], row['epsilon'], row['q'], 'departure',
                    row['traceless_compression_norm'], 'signature', row['local_signature'], flush=True)
    path = ROOT / 'simulations/results/route_b_n6_end_profile_probe.json'
    path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(path)


if __name__ == '__main__':
    main()
