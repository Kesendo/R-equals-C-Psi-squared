"""Strict raw-pencil contract and falsifiable integer discriminant certificates."""
import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import sympy as sp

from simulations.route_b_a2_n6 import (
    Lambda, t, load_exact_pencils, prove_layer_identity,
    certify_layer_identity, coefficient_difference_bound, discriminant_degree_bound,
)

RAW = Path(__file__).parents[1] / "results/route_b_a2_n6_residual.tmp.json"


@pytest.fixture
def raw():
    if not RAW.exists():
        pytest.skip("Generate raw input with route-b-n6-residual --out first")
    return json.loads(RAW.read_text())


def write_raw(tmp_path, data):
    path = tmp_path / "raw.json"
    path.write_text(json.dumps(data))
    return path


def test_load_exact_axes_and_transport(raw, tmp_path):
    pencils = load_exact_pencils(write_raw(tmp_path, raw))
    e, o = pencils["E"], pencils["O"]
    assert e.residual.domain == sp.ZZ
    assert e.residual.degree(Lambda) == 32
    assert e.at.degree(Lambda) == 13
    assert e.residual.coeff_monomial(t**3) == int(raw["rEven"]["residualInT"][0][3])
    assert o.residual == sp.Poly(e.residual.as_expr().subs(t, -t), Lambda, t)


@pytest.mark.parametrize("mutation", [
    lambda d: d.update(schemaVersion=2), lambda d: d.update(n=True),
    lambda d: d.update(extra=0), lambda d: d.pop("rOdd"),
    lambda d: d.update(lambdaConvention="lambda=Lambda"),
    lambda d: d.update(parameterConvention="t=qCSharp"),
    lambda d: d.update(coefficientOrder="descending"),
    lambda d: d["rOdd"].update(rOdd=False),
    lambda d: d["rEven"].update(sectorDimension=44),
    lambda d: d["rEven"].update(imaginary=1),
    lambda d: d["rEven"]["residualInT"][0].__setitem__(0, 3),
    lambda d: d["rEven"]["residualInT"][0].__setitem__(0, "3.0"),
    lambda d: d["rEven"]["residualInT"][0].__setitem__(0, "1j"),
    lambda d: d["rEven"]["residualInT"].__setitem__(4, []),
    lambda d: d["rEven"]["residualInT"][0].append("0"),
    lambda d: d["rEven"]["residualInT"][-1].__setitem__(0, "2"),
    lambda d: d["rEven"]["atFactorInT"].pop(),
    lambda d: d["rOdd"]["residualInT"][0].__setitem__(0, "1"),
    lambda d: d["rOdd"]["atFactorInT"][0].__setitem__(0, "1"),
])
def test_load_rejects_mutation(raw, tmp_path, mutation):
    damaged = copy.deepcopy(raw)
    mutation(damaged)
    with pytest.raises(ValueError):
        load_exact_pencils(write_raw(tmp_path, damaged))


def test_load_rejects_duplicate_json_key(raw, tmp_path):
    path = write_raw(tmp_path, raw)
    path.write_text(path.read_text().replace('"n": 6', '"n": 6, "n": 6'))
    with pytest.raises(ValueError, match="duplicate"):
        load_exact_pencils(path)


def toy():
    # Quadratic roots collide at -1; a quadratic root meets 2 at t=3.
    f = sp.Poly((Lambda**2-t-1)*(Lambda-2), Lambda, t, domain=sp.ZZ)
    return f, sp.Poly(t+1,t), sp.Poly(t-3,t)


def test_layer_exact_toy_and_bound():
    f, a1, a2 = toy()
    cert = certify_layer_identity(f, 0, a1, a2, constant=4)
    exact = sp.Poly(sp.discriminant(f.as_expr(), Lambda), t)
    assert exact == 4*a1*a2**2
    assert cert.discriminant_lhs == exact.as_expr()
    assert sp.expand(cert.discriminant_rhs-cert.discriminant_lhs) == 0
    assert cert.proof_modulus > 2*cert.proof_bound
    assert cert.proof_bound == coefficient_difference_bound(f, 4, a1, a2)
    assert max(map(abs, exact.all_coeffs())) <= cert.proof_bound
    assert exact.degree() <= discriminant_degree_bound(f)


def test_layer_reconstructs_constant_from_independent_discriminant():
    f, a1, a2 = toy()
    cert = certify_layer_identity(f, 0, a1, a2)
    assert cert.constant == 4
    assert cert.proof_modulus > 2*cert.proof_bound


@pytest.mark.parametrize("expression", [
    Lambda**2+t**7+3*t+1,
    Lambda**3+t*Lambda**2+(t**2+1)*Lambda+t**3+2,
    Lambda**4+(t**5+1)*Lambda**3+t**2*Lambda+3,
    Lambda**5+Lambda+t,
])
def test_layer_sylvester_support_bound_against_exact_discriminant(expression):
    f = sp.Poly(expression,Lambda,t,domain=sp.ZZ)
    exact = sp.Poly(sp.discriminant(expression,Lambda),t)
    assert exact.degree() <= discriminant_degree_bound(f)
    # Derivative resultant has the opposite sign for degrees 2 and 3.
    from simulations.route_b_a2_n6 import _disc_stream
    from simulations.o2b_gcd_certificate import prime_stream
    stream = _disc_stream(f,prime_stream(),1)
    p, modular = next(stream)
    stream.close()
    assert list(map(int,modular)) == [int(c)%p for c in exact.all_coeffs()]


@pytest.mark.parametrize("field", ["residual", "valuation", "constant", "a1", "a2"])
def test_layer_identity_mutation(field):
    f, a1, a2 = toy()
    args = dict(valuation=0, a1=a1, a2=a2, constant=4)
    if field == "residual":
        f = sp.Poly(f.as_expr()+1, Lambda,t)
    elif field == "valuation":
        args[field] = 1
    elif field == "constant":
        args[field] = 5
    else:
        args[field] += sp.Poly(1,t)
    with pytest.raises(ValueError):
        certify_layer_identity(f, **args)


def test_layer_bound_prevents_sampled_false_identity():
    f, a1, a2 = toy()
    # A change invisible at the first proof prime must still be rejected.
    from simulations.o2b_gcd_certificate import prime_stream
    p = next(prime_stream())
    with pytest.raises(ValueError):
        certify_layer_identity(f, 0, a1, a2, constant=4+p)


@pytest.mark.parametrize("case", ["repeated", "shared"])
def test_layer_rejects_wrong_shape_without_unbounded_prime_search(case):
    script = '''
import sympy as sp
from simulations.route_b_a2_n6 import Lambda,t,certify_layer_identity
case = CASE
f = sp.Poly(Lambda**2-(t+1)**(2 if case == 'repeated' else 3),Lambda,t)
a1 = sp.Poly((t+1)**(2 if case == 'repeated' else 1),t)
a2 = sp.Poly(1 if case == 'repeated' else t+1,t)
try:
    certify_layer_identity(f,0,a1,a2,constant=4)
except ValueError:
    raise SystemExit(0)
raise SystemExit('incorrectly certified a non-squarefree or noncoprime decomposition')
'''.replace("CASE",repr(case))
    result = subprocess.run([sys.executable,"-c",script],capture_output=True,text=True,timeout=4)
    assert result.returncode == 0, result.stderr


@pytest.fixture(scope="module")
def canonical_proof():
    if os.environ.get("ROUTE_B_N6_FULL_PROOF") != "1":
        pytest.skip("full CRT proof opt-in: ROUTE_B_N6_FULL_PROOF=1")
    if not RAW.exists():
        pytest.skip("Generate raw input with route-b-n6-residual --out first")
    e = load_exact_pencils(RAW)["E"]
    cert = prove_layer_identity(e, workers=int(os.environ.get("ROUTE_B_N6_WORKERS", "1")))
    return e, cert


def test_layer_actual_canonical_e(canonical_proof):
    _, cert = canonical_proof
    assert (cert.deg_d, cert.valuation, cert.a1.degree(), cert.a2.degree()) == (926,536,124,133)
    assert cert.constant != 0
    assert cert.proof_modulus > 2*cert.proof_bound
    assert sp.expand(cert.discriminant_rhs-cert.discriminant_lhs) == 0


def test_benchmark_selects_workers_from_timings_and_memory():
    from simulations.route_b_a2_n6 import recommend_workers
    samples = [{"workers":1,"seconds_per_prime":0.3,"peak_tree_rss_bytes":100},
               {"workers":2,"seconds_per_prime":0.2,"peak_tree_rss_bytes":180},
               {"workers":4,"seconds_per_prime":0.4,"peak_tree_rss_bytes":300}]
    assert recommend_workers(samples,available_memory=1000) == 2
    assert recommend_workers(samples,available_memory=250) == 1


@pytest.mark.skipif(os.environ.get("ROUTE_B_N6_BENCHMARK_TEST") != "1",reason="fresh C# benchmark opt-in")
def test_benchmark_measures_export_and_n6_roots():
    from simulations.route_b_a2_n6 import benchmark
    report = benchmark(RAW)
    assert report["export"]["returncode"] == 0
    assert report["export"]["seconds"] > 0
    assert report["export"]["peak_tree_rss_bytes"] > 0
    assert report["s1"]["peak_tree_rss_bytes"] > 0
    assert report["roots"]["source"] == "fresh N6 E A2 CRT candidate; timing only"
    assert set(report["roots"]["probes"]) == {"real","nonreal"}
    assert report["recommended_workers"] in [sample["workers"] for sample in report["worker_samples"]]
    if report["s1"]["status"] == "timeout":
        assert report["task4_checkpoint"]["required_route"] == "CRT_WITH_RIGOROUS_COEFFICIENT_BOUND"


@pytest.mark.parametrize("field", ["residual", "valuation", "constant", "a1", "a2"])
def test_layer_actual_canonical_mutation(canonical_proof, field):
    e, cert = canonical_proof
    f = e.residual
    args = dict(valuation=cert.valuation, a1=cert.a1, a2=cert.a2, constant=cert.constant)
    if field == "residual":
        f = sp.Poly(f.as_expr()+1,Lambda,t,domain=sp.ZZ)
    elif field in ("valuation", "constant"):
        args[field] += 1
    else:
        args[field] += sp.Poly(1,t)
    with pytest.raises(ValueError):
        certify_layer_identity(f,**args)
