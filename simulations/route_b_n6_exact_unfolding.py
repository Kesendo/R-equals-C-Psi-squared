"""Exact arithmetic certificates for the Route-B N6 local unfolding theorem.

Run: python simulations/route_b_n6_exact_unfolding.py
The proof tying these certificates together is PROOF_ROUTE_B_N6_UNFOLDING.md.
All matrix arithmetic below is integer or finite-field arithmetic; no eigensolver.
L(t)=D+tT, t=i*qCSharp, gamma=1, delta=0, physical lambda (Lambda=2*lambda).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

from flint import fmpz_mat, fmpz_poly, fmpz_mod_poly_ctx, nmod_mat, nmod_poly
import sympy as sp

import route_b_a2_n6 as exact

ROOT = Path(__file__).resolve().parents[1]


def integer_pencil():
    basis = [(a, b) for a in range(64) if a.bit_count() == 1
             for b in range(64) if b.bit_count() == 2]
    index = {pair: j for j, pair in enumerate(basis)}
    n = len(basis)
    d = [[0]*n for _ in basis]
    bonds = [[[0]*n for _ in basis] for _ in range(5)]
    for j, (a, b) in enumerate(basis):
        d[j][j] = -2*(a ^ b).bit_count()
        for s in range(5):
            mask = 3 << s
            if (a & mask).bit_count() == 1:
                bonds[s][index[(a ^ mask, b)]][j] = -2
            if (b & mask).bit_count() == 1:
                bonds[s][index[(a, b ^ mask)]][j] = 2
    reverse = lambda a: int(f'{a:06b}'[::-1], 2)
    reflection = [index[(reverse(a), reverse(b))] for a, b in basis]
    stagger = [(-1)**sum(((a >> s) & 1)+((b >> s) & 1) for s in (1, 3, 5))
               for a, b in basis]
    return d, bonds, reflection, stagger


def evaluate_rows(rows, t):
    return fmpz_poly([sum(int(c)*t**j for j, c in enumerate(row)) for row in rows])


def exact_builder_certificate(d, bonds, reflection, stagger, fixture):
    n = len(d)
    hop = [[sum(b[i][j] for b in bonds) for j in range(n)] for i in range(n)]
    assert d == list(map(list, zip(*d)))
    assert hop == list(map(list, zip(*hop)))
    assert all(reflection[reflection[i]] == i and reflection[i] != i for i in range(n))
    for i in range(n):
        for j in range(n):
            assert d[reflection[i]][reflection[j]] == d[i][j]
            assert hop[reflection[i]][reflection[j]] == hop[i][j]
            assert stagger[i]*d[i][j]*stagger[j] == d[i][j]
            assert stagger[i]*hop[i][j]*stagger[j] == -hop[i][j]
            for s in range(5):
                assert bonds[s][reflection[i]][reflection[j]] == bonds[4-s][i][j]
                assert stagger[i]*bonds[s][i][j]*stagger[j] == -bonds[s][i][j]
        assert stagger[reflection[i]] == -stagger[i]
    representatives = [i for i in range(n) if i < reflection[i]]
    # No fixed reflection orbit: both sectors have constant Gram matrix 2I.
    # Every characteristic coefficient has t-degree <=45. 46 exact evaluations
    # therefore prove a polynomial identity, not a sampled numerical comparison.
    for sign, name in ((1, 'rEven'), (-1, 'rOdd')):
        assert max(len(row)-1 for row in fixture[name]['residualInT'])+max(
            len(row)-1 for row in fixture[name]['atFactorInT']) <= 45
        for t in range(46):
            sector = fmpz_mat([[2*(d[i][j]+t*hop[i][j]
                +sign*(d[i][reflection[j]]+t*hop[i][reflection[j]]))
                for j in representatives] for i in representatives])
            expected = evaluate_rows(fixture[name]['residualInT'], t)*evaluate_rows(
                fixture[name]['atFactorInT'], t)
            assert sector.charpoly() == expected, (name, t)
        print('integer characteristic identity', name, '46/46', flush=True)
    return hop, {'dimension': n, 'sector_dimension': len(representatives),
                 'degree_bound': 45, 'exact_evaluations_per_sector': 46,
                 'hermitian_real_t': True, 'stagger_swaps_reflection_parity': True}


def rabin_irreducibility(coefficients, p=367):
    context = fmpz_mod_poly_ctx(p)
    f = context(coefficients)
    assert f.degree() == 133
    x = context([0, 1])
    assert (pow(x, p**133, f)-x) % f == 0
    for divisor in (7, 19):
        assert (pow(x, p**(133//divisor), f)-x).gcd(f).degree() == 0
    # A primitive odd-degree irreducible polynomial has a real embedding.
    return {'prime': p, 'degree': 133, 'prime_divisors_of_degree': [7, 19],
            'frobenius_remainder_zero': True, 'proper_frobenius_gcd_degrees': [0, 0]}


def modular_matrix(rows, p):
    return nmod_mat(len(rows), len(rows[0]), [int(c) % p for row in rows for c in row], p)


def trace(matrix, p):
    return sum(int(matrix[i, i]) for i in range(matrix.nrows())) % p


def discriminant(matrix, p):
    return (2*trace(matrix*matrix, p)-trace(matrix, p)**2) % p


def matrix_polynomial(f, m, identity):
    out = identity*0
    for coefficient in reversed(f.coeffs()):
        out = out*m+identity*int(coefficient)
    return out


def modular_witness(d, hop, bonds, reflection, a2_coefficients, fixture, p=101):
    context = fmpz_mod_poly_ctx(p)
    a2 = context(a2_coefficients)
    assert a2.degree() == 133
    candidates = sorted(int(-f[0]/f[1]) for f, multiplicity in a2.factor()[1]
                        if f.degree() == 1 and multiplicity == 1)
    n = len(d)
    identity = modular_matrix([[int(i == j) for j in range(n)] for i in range(n)], p)
    dm, tm = modular_matrix(d, p), modular_matrix(hop, p)
    left, right = modular_matrix(bonds[0], p), modular_matrix(bonds[4], p)
    inv2 = pow(2, -1, p)
    reflection_matrix = modular_matrix([[int(reflection[j] == i) for j in range(n)] for i in range(n)], p)
    for t in candidates:
        f = nmod_poly([int(c) % p for c in evaluate_rows(fixture['rEven']['residualInT'], t).coeffs()], p)
        gcd = f.gcd(f.derivative())
        if gcd.degree() != 1 or f.degree() != 32:
            continue
        cleared_lambda = int(-gcd[0]/gcd[1])
        lam = cleared_lambda*inv2 % p
        m = dm+tm*t
        characteristic = m.charpoly()
        linear = nmod_poly([-lam, 1], p)
        quotient, remainder = divmod(characteristic, linear*linear)
        if remainder or int(quotient(lam)) == 0:
            continue
        projector = matrix_polynomial(quotient, m, identity)*pow(int(quotient(lam)), -1, p)
        zero = identity*0
        if projector*projector != projector or (m-identity*lam)*projector != zero:
            continue
        assert projector.rref()[1] == 2
        assert reflection_matrix*projector == projector
        complement = identity-projector
        resolvent = complement*(identity*lam-m+projector).inv()*complement
        assert (identity*lam-m)*resolvent == complement
        a = projector*tm*projector
        alpha = discriminant(a, p)
        profiles = {}
        for name, direction in (('one', left), ('even', (left+right)*inv2),
                                ('odd', (left-right)*inv2)):
            first = projector*(direction*t)*projector
            if name == 'odd':
                assert first == zero
                b = projector*(direction*t)*resolvent*(direction*t)*projector
            else:
                b = first
            beta = (4*trace(a*b, p)-2*trace(a, p)*trace(b, p)) % p
            gamma = discriminant(b, p)
            slope_discriminant = (beta*beta-4*alpha*gamma) % p
            profiles[name] = {'alpha': alpha, 'beta': beta, 'gamma': gamma,
                              'slope_discriminant': slope_discriminant}
        if not alpha or not all(r['slope_discriminant'] for r in profiles.values()):
            continue
        assert profiles['one'] == profiles['even']
        # Break-inputs: scalar and commuting effective responses do NOT unfold.
        controls = {}
        for name, b in (('scalar', projector), ('commuting', a)):
            beta = (4*trace(a*b, p)-2*trace(a, p)*trace(b, p)) % p
            controls[name] = (beta*beta-4*alpha*discriminant(b, p)) % p
            assert controls[name] == 0
        # Adding even admixture destroys the odd first-order suppression.
        assert projector*(left*t)*projector != zero
        return {'prime': p, 't': t, 'lambda_cleared': cleared_lambda,
                'lambda_physical': lam, 'a2_derivative': int(a2.derivative()(t)),
                'residual_gcd_degree': gcd.degree(), 'full_algebraic_multiplicity': 2,
                'projector_rank': 2, 'complement_characteristic_at_lambda': int(quotient(lam)),
                'profiles': profiles, 'non_unfolding_controls': controls,
                'even_admixture_breaks_first_order_suppression': True}
    raise ValueError('No valid witness among the modular linear factors')


def main():
    started = time.perf_counter()
    fixture_path = ROOT/'simulations/tests/fixtures/route_b_a2_n6_residual.json'
    fixture = json.loads(fixture_path.read_text())
    d, bonds, reflection, stagger = integer_pencil()
    hop, builder = exact_builder_certificate(d, bonds, reflection, stagger, fixture)
    # Reuse the existing bounded-CRT producer, including its independently
    # reconstructed discriminant and rigorous integer coefficient bound.
    layer = exact.prove_layer_identity(exact.load_exact_pencils(fixture_path)['E'], workers=1)
    coefficients = list(reversed([int(c) for c in layer.a2.all_coeffs()]))
    irreducibility = rabin_irreducibility(coefficients)
    try:
        rabin_irreducibility([1]+[0]*132+[1])  # t^133+1 has the factor t+1.
    except AssertionError:
        reducible_control_rejected = True
    else:
        raise AssertionError('Rabin check accepted a reducible polynomial')
    witness = modular_witness(d, hop, bonds, reflection, coefficients, fixture)
    report = {'conventions': __doc__, 'builder': builder, 'irreducibility': irreducibility,
        'reducible_control_rejected': reducible_control_rejected,
        'layer': {'a2_coefficients_lowest_first': list(map(str, coefficients)),
                  'proof_modulus': str(layer.proof_modulus), 'proof_bound': str(layer.proof_bound),
                  'proof_prime_count': len(layer.proof_primes), 'a2_degree': layer.a2.degree(),
                  'discriminant_degree': layer.deg_d, 'a2_exponent': 2},
        'modular_witness': witness,
        'provenance': {'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'fixture_sha256': hashlib.sha256(fixture_path.read_bytes()).hexdigest(),
            'source_pencil_digest': exact.N6_SOURCE_PENCIL_DIGEST},
        'runtime_seconds': time.perf_counter()-started}
    output = ROOT/'simulations/results/route_b_n6_exact_unfolding.json'
    output.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(witness, indent=2), flush=True)
    print(output, 'seconds', report['runtime_seconds'], flush=True)


if __name__ == '__main__':
    main()
