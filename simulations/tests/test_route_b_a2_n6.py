"""Strict raw-pencil contract and falsifiable integer discriminant certificates."""
import copy
from fractions import Fraction
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
    LayerCertificate, first_subresultant_linear_exact, prove_pair_uniqueness,
    prove_at_nonoverlap, isolate_a2_roots, _write_generic_inventory as write_inventory,
    fujiwara_integer_radius, verify_fujiwara_radius,
)

RAW = Path(__file__).parent / "fixtures/route_b_a2_n6_residual.json"


@pytest.fixture(autouse=True)
def count_capability(monkeypatch):
    """One local test-run credential, explicitly propagated, never production state."""
    from functools import wraps
    from simulations import route_b_a2_n6 as m
    capability=os.urandom(32)
    for name in ("certify_a2_boxes","certify_odd_transport","isolate_a2_roots",
                 "verify_box_certificate","_write_generic_inventory","write_inventory"):
        original=getattr(m,name)
        @wraps(original)
        def scoped(*args,_original=original,**kwargs):
            kwargs.setdefault("capability",capability)
            return _original(*args,**kwargs)
        monkeypatch.setattr(m,name,scoped)
    original_phase=m._bounded_exact_phase
    @wraps(original_phase)
    def phase(kind,*args,**kwargs):
        if kind in {"certify","certify_odd","isolate","write","write_generic"}:
            kwargs.setdefault("capability",capability)
        return original_phase(kind,*args,**kwargs)
    monkeypatch.setattr(m,"_bounded_exact_phase",phase)
    monkeypatch.setitem(globals(),"isolate_a2_roots",m.isolate_a2_roots)
    monkeypatch.setitem(globals(),"write_inventory",m._write_generic_inventory)
    return capability


@pytest.fixture
def raw():
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
    # Runs by default. These are the only tests in the module that EXECUTE the CRT
    # layer proof rather than reading the stored artifact, so without them the only
    # check on the layer identity is a comparison against literals held by the module
    # that wrote it. Module-scoped, so the proof is paid once for the seven tests that
    # consume it. The opt-out is for iterating on the cheap tests, not for reporting
    # a result.
    if os.environ.get("ROUTE_B_N6_SKIP_FULL_PROOF") == "1":
        pytest.skip("full CRT proof skipped by request: ROUTE_B_N6_SKIP_FULL_PROOF=1")
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
def test_benchmark_measures_export_and_n6_roots(tmp_path):
    from simulations.route_b_a2_n6 import benchmark
    (tmp_path/"residual.json").write_bytes(RAW.read_bytes())
    report = benchmark(tmp_path/"residual.json")
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


def test_subresultant_fixed_minor_toy_retains_parameter_content():
    f = sp.Poly(Lambda**3 + (t-1)*Lambda + (t-1), Lambda, t, domain=sp.ZZ)
    cert = first_subresultant_linear_exact(f)
    assert cert.s1 == sp.Poly((t-1)*(6*Lambda+9), Lambda,t,domain=sp.ZZ)
    assert cert.a_raw == sp.Poly(6*(t-1),t,domain=sp.ZZ)
    assert cert.b_raw == sp.Poly(9*(t-1),t,domain=sp.ZZ)
    assert cert.psc1_raw == cert.a_raw
    assert cert.verification_modulus > 2*cert.verification_bound
    with pytest.raises(ValueError, match="pair uniqueness"):
        prove_pair_uniqueness(sp.Poly(t-1,t), cert)


def test_pair_rejects_two_simultaneous_double_roots():
    f = sp.Poly((Lambda**2-1)**2+t, Lambda,t,domain=sp.ZZ)
    cert = first_subresultant_linear_exact(f)
    with pytest.raises(ValueError, match="pair uniqueness"):
        prove_pair_uniqueness(sp.Poly(t,t),cert)


def test_overlap_gate_rejects_at_strand_at_repeated_seed():
    f = sp.Poly((Lambda**2-t)*(Lambda-3),Lambda,t,domain=sp.ZZ)
    cert = first_subresultant_linear_exact(f)
    a2 = sp.Poly(t,t,domain=sp.ZZ)
    prove_pair_uniqueness(a2,cert)
    with pytest.raises(ValueError,match="AT overlap"):
        prove_at_nonoverlap(a2,sp.Poly(Lambda,Lambda,t,domain=sp.ZZ),cert)


def test_isolate_and_export_reject_structural_mutations(tmp_path):
    from simulations.route_b_a2_n6 import certify_odd_transport
    a2 = sp.Poly(t**3-2,t,domain=sp.ZZ)
    roots = isolate_a2_roots(a2)
    odd = certify_odd_transport(a2,roots)
    assert len(roots) == 3
    assert sum(box.root_count for box in roots) == 3
    # A tiny exact double-root family with seed Lambda=2 exercises all transports.
    f = source_bound_toy()[0].residual
    s1 = first_subresultant_linear_exact(f)
    layer = source_bound_toy()[1]
    out = tmp_path/"inventory.json"
    payload = write_inventory(out,layer,s1,a2,roots,odd,source_pencil=source_bound_toy()[0])
    assert len(payload["loci"]) == 6
    assert out.read_bytes().endswith(b"\n") and not out.read_bytes().endswith(b"\n\n")
    for mutation in ("remove","duplicate","overlap","conjugation","parity","plusminus"):
        damaged = copy.deepcopy(payload)
        if mutation == "remove": damaged["loci"].pop()
        elif mutation == "duplicate": damaged["loci"].append(copy.deepcopy(damaged["loci"][0]))
        elif mutation == "overlap": damaged["loci"][1]["tBox"] = damaged["loci"][0]["tBox"]
        elif mutation == "conjugation": damaged["loci"][0]["conjugationPartnerId"] = "bad"
        elif mutation == "parity": damaged["loci"][0]["parityPartnerId"] = "bad"
        else: damaged["loci"][0]["qPhysicalCSharpBox"] = damaged["loci"][0]["tBox"]
        with pytest.raises(ValueError):
            write_inventory(tmp_path/f"{mutation}.json",layer,s1,a2,roots,odd,payload=damaged,source_pencil=source_bound_toy()[0])


@pytest.mark.parametrize("poly,real_count",[
    (sp.Poly(t**3-2,t,domain=sp.ZZ),1),
    (sp.Poly(3*t**2-12,t,domain=sp.ZZ),2),
    (sp.Poly(7*t**4+0*t**3-28*t**2+0*t+7,t,domain=sp.ZZ),4),
])
def test_isolate_fujiwara_integer_bound_is_exact_and_mutation_sensitive(poly,real_count):
    radius=fujiwara_integer_radius(poly)
    assert type(radius) is int and radius>=1
    assert verify_fujiwara_radius(poly,radius)
    # Direct exact root counts, not floating roots, show the known examples fit.
    assert poly.count_roots(-radius,radius) == real_count
    if radius>1:
        with pytest.raises(ValueError,match="Fujiwara"):
            verify_fujiwara_radius(poly,radius-1)


def test_isolate_fujiwara_rejects_missing_polynomial():
    for bad in (sp.Poly(0,t,domain=sp.ZZ),sp.Poly(3,t,domain=sp.ZZ)):
        with pytest.raises(ValueError):
            fujiwara_integer_radius(bad)


def test_hybrid_requires_executed_count_and_outward_containment(monkeypatch):
    from fractions import Fraction
    from simulations import route_b_a2_n6 as m
    assert hasattr(m, "ProposedBall"), "hybrid proposal type is missing"
    ball=m.ProposedBall("0", "1", "0.001", "0.001", 128)
    lower=m.ProposedBall("0", "-1", "0.001", "0.001", 128)
    box=m.outward_dyadic_box(ball, 20)
    assert box.real[0] <= Fraction(-1,1000) and box.real[1] >= Fraction(1,1000)
    assert box.root_count is None
    p=sp.Poly(t*t+1,t,domain=sp.ZZ)
    calls=[]
    original=m.exact_count_in_box
    def count(poly, candidate):
        calls.append(candidate)
        return original(poly,candidate)
    monkeypatch.setattr(m,"exact_count_in_box",count)
    roots=m.certify_a2_boxes(p,[ball,lower],dyadic_bits=20)
    assert len(roots)==2 and len(calls)>=1
    monkeypatch.setattr(m,"exact_count_in_box",lambda *_:0)
    with pytest.raises(ValueError,match="count"):
        m.certify_a2_boxes(p,[ball,lower],dyadic_bits=20)


def test_local_real_certification_never_calls_global_isolators(monkeypatch):
    from simulations import route_b_a2_n6 as m
    import sympy.polys.rootisolation as isolation
    def forbidden(*_args,**_kwargs):
        pytest.fail("global root isolation was called")
    monkeypatch.setattr(sp.Poly,"intervals",forbidden)
    monkeypatch.setattr(isolation,"dup_isolate_real_roots_sqf",forbidden)
    monkeypatch.setattr(isolation,"dup_isolate_complex_roots_sqf",forbidden)
    p=sp.Poly((t-2)*(t*t+1),t,domain=sp.ZZ)
    balls=[m.ProposedBall("2","0","0.001","0",128),
           m.ProposedBall("0","1","0.001","0.001",128),
           m.ProposedBall("0","-1","0.001","0.001",128)]
    roots=m.certify_a2_boxes(p,balls,dyadic_bits=96)
    assert len(roots)==3 and sum(box.root_count for box in roots)==3
    assert sum(box.imag==(0,0) for box in roots)==1


def test_local_real_count_is_executed_and_rejects_boundary(monkeypatch):
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"exact_count_in_real_interval"), "local real counter missing"
    p=sp.Poly(t-2,t,domain=sp.ZZ)
    assert m.exact_count_in_real_interval(p,1,3)==1
    assert m.exact_count_in_real_interval(p,3,4)==0
    for ends in ((2,3),(1,2),(2,2),(3,1)):
        with pytest.raises(ValueError): m.exact_count_in_real_interval(p,*ends)
    monkeypatch.setattr(m,"exact_count_in_real_interval",lambda *_:0)
    with pytest.raises(ValueError,match="count"):
        m.certify_a2_boxes(p,[m.ProposedBall("2","0","0.001","0",128)])


def test_real_proposal_zero_axis_is_preserved_outward():
    from simulations import route_b_a2_n6 as m
    ball=m.ProposedBall("2","0","0.001","0",128)
    box=m.outward_dyadic_box(ball,96)
    assert box.imag==(0,0) and box.root_count is None
    m._check_containment(ball,box)


def test_no_global_intervals_call_remains_in_producer():
    import ast
    from simulations import route_b_a2_n6 as m
    tree=ast.parse(Path(m.__file__).read_text(encoding="utf-8"))
    forbidden=[node for node in ast.walk(tree) if isinstance(node,ast.Call)
               and isinstance(node.func,ast.Attribute) and node.func.attr=="intervals"]
    assert not forbidden, "global Poly.intervals call remains"


def test_exact_sturm_chain_counts_match_sympy_without_real_count_roots(monkeypatch):
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"build_exact_sturm_chain"), "exact Sturm chain missing"
    p=sp.Poly((t-2)*(t+1)*(t*t+1),t,domain=sp.ZZ)
    intervals=[(Fraction(-3),Fraction(-2)),(Fraction(-2),Fraction(0)),
               (Fraction(1),Fraction(3)),(Fraction(-3),Fraction(3))]
    expected=[int(p.count_roots(*ends)) for ends in intervals]
    chain=m.build_exact_sturm_chain(p)
    assert m.verify_exact_sturm_chain(p,chain)
    def forbidden(*_args,**_kwargs): pytest.fail("per-interval SymPy real count called")
    monkeypatch.setattr(sp.Poly,"count_roots",forbidden)
    assert [m.exact_count_in_real_interval(chain,*ends) for ends in intervals]==expected
    with pytest.raises(ValueError,match="zero|boundary"):
        m.sturm_variations_at(chain,Fraction(2))
    monkeypatch.setattr(m,"sturm_variations_at",lambda _chain,x: int(x))
    with pytest.raises(ValueError,match="negative|variation"):
        m.exact_count_in_real_interval(chain,Fraction(1),Fraction(3))


@pytest.mark.parametrize("mutation",["order","derivative","remainder_sign","member","terminal"])
def test_exact_sturm_chain_rejects_mutations(mutation):
    from dataclasses import replace
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"build_exact_sturm_chain"), "exact Sturm chain missing"
    p=sp.Poly(t**3-2*t+2,t,domain=sp.ZZ)
    chain=m.build_exact_sturm_chain(p)
    members=list(chain.members)
    index={"order":0,"derivative":1,"remainder_sign":2,"member":2,"terminal":len(members)-1}[mutation]
    member=list(members[index])
    if mutation=="order": member.reverse()
    elif mutation=="remainder_sign": member=[(-n,d) for n,d in member]
    elif mutation=="terminal": member=[(0,1)]
    else: member[0]=(member[0][0]+member[0][1],member[0][1])
    members[index]=tuple(member)
    with pytest.raises(ValueError):
        m.verify_exact_sturm_chain(p,replace(chain,members=tuple(members)))


def test_exact_sturm_chain_rejects_nonsquarefree_and_zero_member_endpoint():
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"build_exact_sturm_chain"), "exact Sturm chain missing"
    with pytest.raises(ValueError): m.build_exact_sturm_chain(sp.Poly((t-1)**2,t))
    chain=m.build_exact_sturm_chain(sp.Poly(t*t-2,t))
    with pytest.raises(ValueError,match="zero"):
        m.sturm_variations_at(chain,Fraction(0))


@pytest.mark.parametrize("mutation",["sign","value","missing"])
def test_exact_sturm_positive_scale_witness_rejects_mutations(mutation):
    from dataclasses import replace
    from simulations import route_b_a2_n6 as m
    chain=m.build_exact_sturm_chain(sp.Poly(t**3-2*t+2,t))
    assert hasattr(chain,"scales"), "positive scale witnesses missing"
    scales=list(chain.scales)
    if mutation=="missing": scales.pop()
    else:
        n,d=scales[0]
        scales[0]=(-n,d) if mutation=="sign" else (n+d,d)
    with pytest.raises(ValueError):
        m.verify_exact_sturm_chain(sp.Poly(t**3-2*t+2,t),replace(chain,scales=tuple(scales)))


def test_exact_sturm_positive_scaling_is_canonical_and_preserves_sign():
    from simulations import route_b_a2_n6 as m
    chain=m.build_exact_sturm_chain(sp.Poly(t**3-2*t+2,t))
    assert hasattr(chain,"scales"), "positive scale witnesses missing"
    assert any(Fraction(n,d)!=1 for n,d in chain.scales)
    for index,(n,d) in enumerate(chain.scales[:-1],1):
        scale=Fraction(n,d)
        assert scale>0
        prev,current=map(m._unpack_sturm_poly,chain.members[index-1:index+1])
        quotient,remainder=divmod(prev,current)
        following=m._unpack_sturm_poly(chain.members[index+1])
        assert following==-remainder*m._require_flint().fmpq(n,d)
        assert following.denom()==1 and following.numer().content()==1


def test_certification_reuses_one_verified_sturm_chain(monkeypatch):
    from simulations import route_b_a2_n6 as m
    p=sp.Poly((t-2)*(t-3),t)
    balls=[m.ProposedBall(str(x),"0","0.001","0",128) for x in (2,3)]
    original=m.build_exact_sturm_chain
    calls=[]
    def build(poly):
        calls.append(poly)
        return original(poly)
    monkeypatch.setattr(m,"build_exact_sturm_chain",build)
    def forbidden(*_args,**_kwargs): pytest.fail("per-interval real SymPy counter called")
    monkeypatch.setattr(sp.Poly,"count_roots",forbidden)
    roots=m.certify_a2_boxes(p,balls,dyadic_bits=96)
    assert [root.root_count for root in roots]==[1,1]
    assert len(calls)==1


def test_real_counter_has_no_sympy_count_roots_call():
    import ast,inspect
    from simulations import route_b_a2_n6 as m
    tree=ast.parse(inspect.getsource(m.exact_count_in_real_interval))
    assert not [node for node in ast.walk(tree) if isinstance(node,ast.Call)
                and isinstance(node.func,ast.Attribute) and node.func.attr=="count_roots"]


def test_real_proposal_rejects_unresolved_imaginary_interval():
    from simulations import route_b_a2_n6 as m
    with pytest.raises(ValueError,match="classif|axis"):
        m.certify_a2_boxes(sp.Poly(t-2,t),[m.ProposedBall("2","0","0.001","0.001",128)])


def test_local_real_proposals_reject_overlap_and_wrong_location():
    from simulations import route_b_a2_n6 as m
    p=sp.Poly((t-2)*(t-3),t,domain=sp.ZZ)
    with pytest.raises(ValueError,match="overlap"):
        m.certify_a2_boxes(p,[m.ProposedBall("2","0","1","0",128),
                              m.ProposedBall("3","0","1","0",128)])
    with pytest.raises(ValueError,match="count"):
        m.certify_a2_boxes(p,[m.ProposedBall("2","0","0.001","0",128),
                              m.ProposedBall("4","0","0.001","0",128)])


def test_proposals_preserve_exact_integer_input(monkeypatch):
    import flint
    from simulations import route_b_a2_n6 as m
    def forbidden(*_args,**_kwargs):
        pytest.fail("integer polynomial was converted to rounded ACB coefficients")
    monkeypatch.setattr(flint,"acb_poly",forbidden)
    proposals=m.propose_a2_boxes_flint(sp.Poly(t*t+1,t,domain=sp.ZZ),128)
    assert len(proposals)==2


def test_missing_lower_proposal_cannot_pass_completeness():
    from simulations import route_b_a2_n6 as m
    p=sp.Poly(t*t+1,t,domain=sp.ZZ)
    upper=m.ProposedBall("0","1","0.001","0.001",128)
    with pytest.raises(ValueError,match="proposal"):
        m.certify_a2_boxes(p,[upper])


def test_hybrid_rejects_missing_duplicate_and_shrunk_boxes(monkeypatch):
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"ProposedBall"), "hybrid proposal type is missing"
    p=sp.Poly(t*t+1,t,domain=sp.ZZ)
    ball=m.ProposedBall("0","1","0.001","0.001",128)
    lower=m.ProposedBall("0","-1","0.001","0.001",128)
    for proposals in ([],[ball,ball]):
        with pytest.raises(ValueError): m.certify_a2_boxes(p,proposals)
    monkeypatch.setattr(m,"outward_dyadic_box",lambda *_:m.ExactRootBox((0,0),(1,1)))
    with pytest.raises(ValueError,match="contain"):
        m.certify_a2_boxes(p,[ball,lower])


def test_exact_counter_rejects_boundary_root():
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"exact_count_in_box"), "hybrid exact counter is missing"
    with pytest.raises(ValueError,match="boundary"):
        m.exact_count_in_box(sp.Poly(t*t+1,t),m.ExactRootBox((0,1),(0,2)))


def test_decimal_point_seed_is_exact():
    from fractions import Fraction
    from simulations.route_b_a2_n6 import _decimal_in_interval
    for value in (Fraction(1,2),Fraction(-3,8),Fraction(7,20)):
        text,decoded=_decimal_in_interval((value,value))
        assert Fraction(text)==value==decoded


def test_quotient_enclosure_has_bounded_exact_dyadic_endpoints():
    from simulations import route_b_a2_n6 as m
    box=m._lambda_box(sp.Poly(3,t),sp.Poly(-1,t),m.ExactRootBox((0,0),(0,0)))
    assert box.real[0]<=Fraction(1,3)<=box.real[1]
    assert all(v.denominator & (v.denominator-1)==0 for v in box.real+box.imag)


def test_interval_horner_rounds_outward_at_every_step():
    from simulations import route_b_a2_n6 as m
    p=sp.Poly(t**7-3*t**4+2,t,domain=sp.ZZ)
    box=m.ExactRootBox((Fraction(1,3),Fraction(2,3)),(Fraction(1,5),Fraction(2,5)))
    result=m._poly_box(p,box)
    assert all(v.denominator.bit_length()<=257 for v in result.real+result.imag)
    assert all(v.denominator & (v.denominator-1)==0 for v in result.real+result.imag)
    for r in box.real:
        for i in box.imag:
            value=sp.expand_complex(p.as_expr().subs(t,sp.Rational(r)+sp.I*sp.Rational(i)))
            assert result.real[0]<=sp.re(value)<=result.real[1]
            assert result.imag[0]<=sp.im(value)<=result.imag[1]


def test_raw_monomial_content_is_evaluated_without_linear_wrapping():
    from simulations import route_b_a2_n6 as m
    width=Fraction(1,2**64)
    box=m.ExactRootBox((1-width,1+width),(1-width,1+width))
    raw=sp.Poly(t**512,t,domain=sp.ZZ)
    result=m._poly_box(raw,box)
    # (1+i)**512=2**256; |arg(z)-pi/4|<2**-62 on this box.
    # Thus the entire raw power has positive real part.
    assert result.real[0]>0
    assert result.real[0]<=2**256<=result.real[1]


def test_monomial_localization_preserves_raw_quotient_at_a2_root():
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"localize_monomial_quotient"), "approved localization is missing"
    a2=sp.Poly(t+1,t,domain=sp.ZZ)
    a=sp.Poly(t**2*(t+2),t,domain=sp.ZZ)
    b=sp.Poly(t**2*(3*t-1),t,domain=sp.ZZ)
    localized=m.localize_monomial_quotient(a2,a,b,expected_valuation=2)
    assert localized.valuation==2
    assert localized.a_unit==sp.Poly(t+2,t,domain=sp.ZZ)
    assert localized.b_unit==sp.Poly(3*t-1,t,domain=sp.ZZ)
    assert -b.eval(-1)/a.eval(-1)==-localized.b_unit.eval(-1)/localized.a_unit.eval(-1)==4
    box=m.ExactRootBox((Fraction(-1),Fraction(-1)),(Fraction(0),Fraction(0)))
    enclosure=m.localized_lambda_box(localized,box)
    assert enclosure.real==(Fraction(4),Fraction(4)) and enclosure.imag==(0,0)
    assert a==sp.Poly(t**2*(t+2),t,domain=sp.ZZ)
    assert b==sp.Poly(t**2*(3*t-1),t,domain=sp.ZZ)


@pytest.mark.parametrize("case",["unequal","a2_zero","canonical_valuation"])
def test_monomial_localization_rejects_invalid_domain(case):
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"localize_monomial_quotient"), "approved localization is missing"
    a2=sp.Poly(t if case=="a2_zero" else t+1,t,domain=sp.ZZ)
    a=sp.Poly(t**2*(t+2),t,domain=sp.ZZ)
    b=sp.Poly(t**(3 if case=="unequal" else 2)*(3*t-1),t,domain=sp.ZZ)
    with pytest.raises(ValueError):
        m.localize_monomial_quotient(a2,a,b,expected_valuation=496 if case=="canonical_valuation" else 2)


def test_monomial_division_rejects_nonzero_remainder():
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"divide_monomial_exact"), "exact monomial division is missing"
    with pytest.raises(ValueError,match="remainder"):
        m.divide_monomial_exact(sp.Poly(t**2+1,t,domain=sp.ZZ),2)


def test_monomial_localization_rejects_zero_touching_box():
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"localize_monomial_quotient"), "approved localization is missing"
    localized=m.localize_monomial_quotient(sp.Poly(t+1,t),sp.Poly(t**2*(t+2),t),sp.Poly(t**2*(3*t-1),t),expected_valuation=2)
    with pytest.raises(ValueError,match="zero"):
        m.localized_lambda_box(localized,m.ExactRootBox((Fraction(-1),Fraction(0)),(Fraction(0),Fraction(0))))


def test_monomial_localization_rejects_nonmonomial_cancellation():
    from dataclasses import replace
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"localize_monomial_quotient"), "approved localization is missing"
    a2=sp.Poly(t+1,t,domain=sp.ZZ)
    a=sp.Poly(t**2*(t-1)*(t+2),t,domain=sp.ZZ)
    b=sp.Poly(t**2*(t-1)*(3*t-1),t,domain=sp.ZZ)
    localized=m.localize_monomial_quotient(a2,a,b,expected_valuation=2)
    assert localized.a_unit.eval(1)==localized.b_unit.eval(1)==0
    damaged=replace(localized,a_unit=sp.Poly(t+2,t),b_unit=sp.Poly(3*t-1,t))
    with pytest.raises(ValueError,match="monomial"):
        m.verify_monomial_localization(a2,a,b,damaged,expected_valuation=2)


def quotient_ring_toy():
    from simulations import route_b_a2_n6 as m
    a2=sp.Poly(t+1,t,domain=sp.ZZ)
    localized=m.localize_monomial_quotient(a2,sp.Poly(t**2*(t+2),t),
                                          sp.Poly(t**2*(3*t-1),t),expected_valuation=2)
    return m,a2,localized


def test_quotient_ring_exact_toy_and_polynomial_only_evaluation(monkeypatch):
    m,a2,localized=quotient_ring_toy()
    assert hasattr(m,"build_quotient_ring_certificate"), "quotient-ring certificate missing"
    cert=m.build_quotient_ring_certificate(a2,localized)
    assert cert.H==sp.Poly(4,t,domain=sp.ZZ) and cert.D==1
    assert cert.u*localized.a_unit+cert.v*a2==sp.Poly(1,t,domain=sp.QQ)
    assert (localized.a_unit*cert.H+cert.D*localized.b_unit).rem(a2).is_zero
    monkeypatch.setattr(m,"_lambda_box",lambda *_:pytest.fail("parameter-dependent interval division"))
    box=m.ExactRootBox((Fraction(-1),Fraction(-1)),(Fraction(0),Fraction(0)))
    result=m.quotient_ring_lambda_box(cert,box)
    assert result.real==(4,4) and result.imag==(0,0)
    with pytest.raises(ValueError,match="zero"):
        m.quotient_ring_lambda_box(cert,m.ExactRootBox((Fraction(-1),Fraction(0)),(Fraction(0),Fraction(0))))


def test_quotient_ring_rejects_noninvertible_unit():
    from dataclasses import replace
    m,a2,localized=quotient_ring_toy()
    assert hasattr(m,"build_quotient_ring_certificate"), "quotient-ring certificate missing"
    with pytest.raises(ValueError,match="invertible"):
        m.build_quotient_ring_certificate(a2,replace(localized,a_unit=a2))


@pytest.mark.parametrize("field",["u","v","H","D","remainder","common_scale"])
def test_quotient_ring_rejects_mutated_certificate(field):
    from dataclasses import replace
    m,a2,localized=quotient_ring_toy()
    assert hasattr(m,"build_quotient_ring_certificate"), "quotient-ring certificate missing"
    cert=m.build_quotient_ring_certificate(a2,localized)
    if field=="remainder":
        localized=replace(localized,b_unit=localized.b_unit+1)
        damaged=cert
    elif field=="common_scale":
        damaged=replace(cert,H=2*cert.H,D=2*cert.D)
    else:
        damaged=replace(cert,**{field:getattr(cert,field)+(2 if field=="D" else 1)})
    with pytest.raises(ValueError):
        m.verify_quotient_ring_certificate(a2,localized,damaged)


def test_inventory_uses_only_polynomial_interval_evaluation(tmp_path,monkeypatch):
    from simulations import route_b_a2_n6 as m
    a2=sp.Poly(t**3-2,t,domain=sp.ZZ)
    roots=m.isolate_a2_roots(a2)
    odd=m.certify_odd_transport(a2,roots)
    raw=m.first_subresultant_linear_exact(source_bound_toy()[0].residual)
    layer=source_bound_toy()[1]
    monkeypatch.setattr(m,"_lambda_box",lambda *_:pytest.fail("inventory used parameter-dependent interval division"))
    assert len(m._write_generic_inventory(tmp_path/"ring.json",layer,raw,a2,roots,odd,source_pencil=source_bound_toy()[0])["loci"])==6


def phase_split_toy():
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"certify_odd_transport"), "separate O transport certifier missing"
    a2=sp.Poly(t**3-2,t,domain=sp.ZZ)
    e=m.isolate_a2_roots(a2)
    o=m.certify_odd_transport(a2,e)
    raw=m.first_subresultant_linear_exact(source_bound_toy()[0].residual)
    layer=source_bound_toy()[1]
    return m,a2,e,o,raw,layer


def test_phase_split_writer_is_count_and_isolation_free(tmp_path,monkeypatch):
    m,a2,e,o,raw,layer=phase_split_toy()
    def forbidden(*_args,**_kwargs): pytest.fail("writer called count/isolation API")
    for name in ("exact_count_in_box","exact_count_in_real_interval","certify_a2_boxes",
                 "certify_odd_transport","isolate_a2_roots","propose_a2_boxes_flint",
                 "sturm_variations_at","build_exact_sturm_chain"):
        monkeypatch.setattr(m,name,forbidden)
    monkeypatch.setattr(sp.Poly,"count_roots",forbidden)
    monkeypatch.setattr(sp.Poly,"intervals",forbidden)
    payload=m._write_generic_inventory(tmp_path/"split.json",layer,raw,a2,e,o,source_pencil=source_bound_toy()[0])
    assert len(payload["loci"])==6
    assert payload["n"] is None and payload["model"]=="generic-polynomial"
    assert all(row["id"].startswith("GENERIC-") for row in payload["loci"])


@pytest.mark.parametrize("field",["count","partner","box","provenance","seal"])
def test_phase_split_writer_rejects_mutated_carriers(tmp_path,field):
    from dataclasses import replace
    m,a2,e,o,raw,layer=phase_split_toy()
    if field=="count":
        boxes=list(e.boxes); boxes[0]=replace(boxes[0],root_count=2)
        e=replace(e,boxes=tuple(boxes))
    elif field=="partner": e=replace(e,conjugation_partners=(0,)*len(e))
    elif field=="box":
        boxes=list(e.boxes); boxes[0]=replace(boxes[0],real=(Fraction(4),Fraction(5)))
        e=replace(e,boxes=tuple(boxes))
    elif field=="provenance": e=replace(e,count_methods=("literal",)*len(e))
    else: e=replace(e,seal="forged")
    with pytest.raises(ValueError): m._write_generic_inventory(tmp_path/"bad.json",layer,raw,a2,e,o,source_pencil=source_bound_toy()[0])
    assert not (tmp_path/"bad.json").exists()


def test_phase_split_odd_counts_only_odd_boxes(monkeypatch):
    m,a2,e,_o,_raw,_layer=phase_split_toy()
    observed=[]
    original=m.exact_count_in_box
    def count(poly,box):
        observed.append((poly,box))
        return original(poly,box)
    monkeypatch.setattr(m,"exact_count_in_box",count)
    o=m.certify_odd_transport(a2,e)
    assert len(observed)==len(o)==3
    odd=sp.Poly(a2.as_expr().subs(t,-t),t)
    assert all(poly==odd for poly,box in observed if box.imag!=(0,0))


def test_phase_split_carriers_survive_bounded_child_handoff(tmp_path):
    m,a2,_e,_o,raw,layer=phase_split_toy()
    balls=m.propose_a2_boxes_flint(a2,128)
    e=m._bounded_exact_phase("certify",(a2,balls,96,1),timeout_seconds=20)
    o=m._bounded_exact_phase("certify_odd",(a2,e,1),timeout_seconds=20)
    payload=m._bounded_exact_phase("write_generic",(tmp_path/"handoff.json",layer,raw,a2,e,o,None,1,None,source_bound_toy()[0]),timeout_seconds=20)
    assert len(payload["loci"])==6


@pytest.mark.parametrize("field",["count","partner"])
def test_phase_split_cheap_checks_reject_resigned_inconsistent_fields(tmp_path,field,count_capability):
    from dataclasses import replace
    m,a2,e,o,raw,layer=phase_split_toy()
    if field=="count":
        boxes=list(e.boxes); boxes[0]=replace(boxes[0],root_count=2)
        e=replace(e,boxes=tuple(boxes))
    else: e=replace(e,conjugation_partners=(0,)*len(e))
    # Isolate the cheap mathematical shape check from the integrity check.
    e=replace(e,seal=__import__("hmac").new(count_capability,m._carrier_material(e),__import__("hashlib").sha256).hexdigest())
    with pytest.raises(ValueError,match="count|conjugation"):
        m._write_generic_inventory(tmp_path/"bad.json",layer,raw,a2,e,o,source_pencil=source_bound_toy()[0])


def test_exact_gaussian_projective_seed_has_the_correct_complex_sign():
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"_projective_seed_target"), "bounded-size scalar seed evaluation missing"
    assert m._projective_seed_target(sp.Poly(t+2,t),sp.Poly(3*t-1,t),Fraction(0),Fraction(1))==(
        Fraction(-1,5),Fraction(-7,5))


def test_writer_reverifies_supplied_quotient_ring_certificate(tmp_path):
    from dataclasses import replace
    from simulations import route_b_a2_n6 as m
    pencil,layer,raw=source_bound_toy()
    a2=layer.a2
    roots=m.isolate_a2_roots(a2)
    odd=m.certify_odd_transport(a2,roots)
    cert=m.build_quotient_ring_certificate(a2,raw)
    m._write_generic_inventory(tmp_path/"valid.json",layer,raw,a2,roots,odd,quotient_certificate=cert,source_pencil=pencil)
    with pytest.raises(ValueError):
        m._write_generic_inventory(tmp_path/"bad.json",layer,raw,a2,roots,odd,quotient_certificate=replace(cert,H=cert.H+1),source_pencil=source_bound_toy()[0])
    assert not (tmp_path/"bad.json").exists()


def test_bounded_proposal_and_count_phases():
    from simulations import route_b_a2_n6 as m
    p=sp.Poly(t*t+1,t,domain=sp.ZZ)
    proposals=m._bounded_exact_phase("propose",(p,128),timeout_seconds=20)
    boxes=[m.outward_dyadic_box(ball,20) for ball in proposals if Fraction(ball.imag_midpoint)>0]
    assert m._bounded_exact_phase("count",(p,boxes,1),timeout_seconds=20)==[1]
    result=m._bounded_exact_phase("quotient",(sp.Poly(1,t),sp.Poly(-t,t),boxes),timeout_seconds=20)
    assert len(result)==1 and result[0].imag[0]<1<result[0].imag[1]


def test_export_rejects_literal_count_without_execution(tmp_path,monkeypatch):
    from dataclasses import replace
    from simulations import route_b_a2_n6 as m
    a2=sp.Poly(t**3-2,t,domain=sp.ZZ)
    roots=m.isolate_a2_roots(a2)
    odd=m.certify_odd_transport(a2,roots)
    s1=m.first_subresultant_linear_exact(source_bound_toy()[0].residual)
    layer=source_bound_toy()[1]
    forged=replace(roots,boxes=tuple(replace(box,root_count=1) for box in roots),seal="")
    monkeypatch.setattr(m,"exact_count_in_box",lambda *_:pytest.fail("writer recounted forged carrier"))
    with pytest.raises(ValueError,match="count"):
        m._write_generic_inventory(tmp_path/"bad.json",layer,s1,a2,forged,odd,source_pencil=source_bound_toy()[0])
    assert not (tmp_path/"bad.json").exists()


def test_inventory_write_failure_preserves_previous_complete_file(tmp_path,monkeypatch):
    from simulations import route_b_a2_n6 as m
    a2=sp.Poly(t**3-2,t,domain=sp.ZZ)
    roots=m.isolate_a2_roots(a2)
    odd=m.certify_odd_transport(a2,roots)
    s1=m.first_subresultant_linear_exact(source_bound_toy()[0].residual)
    layer=source_bound_toy()[1]
    out=tmp_path/"inventory.json"
    out.write_bytes(b"previous-complete\n")
    def fail(path,*_args,**_kwargs):
        path.write_bytes(b"partial")
        raise OSError("simulated interrupted write")
    monkeypatch.setattr(Path,"write_text",fail)
    with pytest.raises(OSError):
        m._write_generic_inventory(out,layer,s1,a2,roots,odd,source_pencil=source_bound_toy()[0])
    assert out.read_bytes()==b"previous-complete\n"


def test_missing_flint_fails_before_output(tmp_path,monkeypatch):
    from simulations import route_b_a2_n6 as m
    import builtins
    original=builtins.__import__
    def missing(name,*args,**kwargs):
        if name=="flint": raise ImportError("deliberately absent")
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,"__import__",missing)
    monkeypatch.setattr(sys,"argv",["route_b_a2_n6","--residual-json",str(RAW),"--out",str(tmp_path/"absent.json")])
    monkeypatch.setattr(m,"prove_layer_identity",lambda *_:pytest.fail("expensive work before dependency preflight"))
    with pytest.raises(RuntimeError,match="python-flint"): m.main()
    assert not (tmp_path/"absent.json").exists()


@pytest.mark.parametrize("field",["a_raw","b_raw","psc1_raw","degree_bounds","coefficient_bounds","verification_modulus","verification_bound"])
def test_raw_subresultant_certificate_rejects_mutations(field):
    from dataclasses import replace
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"verify_subresultant_certificate"), "raw certificate verifier missing"
    f=sp.Poly(Lambda**3+(t-1)*Lambda+(t-1),Lambda,t,domain=sp.ZZ)
    cert=m.first_subresultant_linear_exact(f)
    assert m.verify_subresultant_certificate(f,cert)
    old=getattr(cert,field)
    new=tuple(x+1 for x in old) if isinstance(old,tuple) else old+1
    with pytest.raises(ValueError):
        m.verify_subresultant_certificate(f,replace(cert,**{field:new}))


@pytest.mark.parametrize("method",["flint","numpy"])
def test_modular_subresultant_specialization_node_rejection(monkeypatch,method):
    from simulations import route_b_a2_n6 as m
    # Degree seven exercises the modular worker, not the degree<=6 exact toy path.
    regular=sp.Poly(Lambda**7+t*Lambda+t,Lambda,t,domain=sp.ZZ)
    shifted=sp.Poly(Lambda**7+(t-1)*Lambda+(t-1),Lambda,t,domain=sp.ZZ)
    bounds=m._fixed_minor_data(shifted)[2]
    assert m._fixed_minor_values(m._poly_matrix(shifted),[1],101,method)==[[0],[0]]
    starts=[]
    original=m._interp_consecutive_modp
    def record(start,values,prime):
        starts.append(start)
        return original(start,values,prime)
    monkeypatch.setattr(m,"_interp_consecutive_modp",record)
    m._minor_prime_job((regular,101,bounds,method))
    assert starts==[1,1]
    starts.clear()
    prime,polynomials=m._minor_prime_job((shifted,101,bounds,method))
    assert starts==[2,2], "simultaneous-zero specialization entered interpolation"
    # A held-out specialization meets the fixed minors through the other backend.
    opposite="numpy" if method=="flint" else "flint"
    exact=m._fixed_minor_values(m._poly_matrix(shifted),[9],prime,opposite)
    assert [int(sp.Poly.from_list(coefficients,t).eval(9))%prime for coefficients in polynomials]==[row[0] for row in exact]


@pytest.mark.parametrize("method",["flint","numpy"])
def test_modular_subresultant_specialization_exhaustion(method):
    from simulations import route_b_a2_n6 as m
    everywhere_degenerate=sp.Poly(Lambda**7,Lambda,t,domain=sp.ZZ)
    with pytest.raises(ValueError,match="too many rejected subresultant specializations"):
        m._minor_prime_job((everywhere_degenerate,101,(0,0),method))


def source_bound_toy():
    from simulations import route_b_a2_n6 as m
    a2=sp.Poly(t**3-2,t,domain=sp.ZZ)
    a1=sp.Poly(t**3+2,t,domain=sp.ZZ)
    pencil=m.ParityPencil("E",sp.Poly((Lambda**2-t**3-2)*(Lambda-2),Lambda,t,domain=sp.ZZ),
                          sp.Poly(Lambda+7,Lambda,t,domain=sp.ZZ))
    return pencil,m.certify_layer_identity(pencil.residual,0,a1,a2,4),m.first_subresultant_linear_exact(pencil.residual)


@pytest.mark.parametrize("field",["valuation","constant","proof_modulus","proof_bound","proof_primes",
                                  "deg_d","discriminant_lhs","discriminant_rhs","b_raw","s1","residual","at","carrier"])
def test_writer_binds_certificates_to_source_pencil(tmp_path,field):
    from dataclasses import replace
    from simulations import route_b_a2_n6 as m
    pencil,layer,s1=source_bound_toy()
    e=m.isolate_a2_roots(layer.a2); o=m.certify_odd_transport(layer.a2,e)
    if field in {"valuation","constant","proof_modulus","proof_bound","deg_d","discriminant_lhs","discriminant_rhs"}:
        layer=replace(layer,**{field:getattr(layer,field)+1})
    elif field=="proof_primes": layer=replace(layer,proof_primes=(3,))
    elif field=="b_raw": s1=replace(s1,b_raw=s1.b_raw+1)
    elif field=="s1": s1=m.first_subresultant_linear_exact(sp.Poly((Lambda+123)**2*(Lambda-3),Lambda,t,domain=sp.ZZ))
    elif field=="residual": pencil=replace(pencil,residual=pencil.residual+1)
    elif field=="at": pencil=replace(pencil,at=sp.Poly(Lambda-2,Lambda,t,domain=sp.ZZ))
    else:
        e=m.isolate_a2_roots(sp.Poly(t**3-3,t)); o=m.certify_odd_transport(sp.Poly(t**3-3,t),e)
    with pytest.raises(ValueError):
        m._write_generic_inventory(tmp_path/"unbound.json",layer,s1,layer.a2,e,o,source_pencil=pencil)
    assert not (tmp_path/"unbound.json").exists()


def test_fujiwara_multicoefficient_counterexample():
    p=sp.Poly(t**4-t**3-2*t**2-4*t-8,t)
    assert p.eval(2)==-16 and p.eval(3)==16 and p.count_roots(2,3)==1
    with pytest.raises(ValueError): verify_fujiwara_radius(p,2)
    assert fujiwara_integer_radius(p)==4


def test_actual_inventory_structural_contract():
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"verify_inventory_artifact_structure")
    payload=m.verify_inventory_artifact_structure(Path(__file__).parents[1]/"results/route_b_a2_n6.json")
    assert len(payload["loci"])==266


def test_actual_inventory_layer_metadata_matches_independent_source_proof(canonical_proof):
    # Separate from the cheap JSON check: the fixture executes the CRT proof.
    from simulations import route_b_a2_n6 as m
    _,layer=canonical_proof
    proof=m.verify_inventory_artifact_structure(Path(__file__).parents[1]/"results/route_b_a2_n6.json")["layerIdentity"]
    assert proof=={"discriminantDegree":layer.deg_d,"valuation":layer.valuation,"a1Degree":layer.a1.degree(),
                   "a2Degree":layer.a2.degree(),"constant":str(layer.constant),
                   "proofModulus":str(layer.proof_modulus),"proofBound":str(layer.proof_bound)}


def test_writer_rejects_missing_source_even_with_valid_carriers(tmp_path):
    m,a2,e,o,raw,layer=phase_split_toy()
    with pytest.raises(ValueError,match="source pencil"):
        m._write_generic_inventory(tmp_path/"missing.json",layer,raw,a2,e,o)


@pytest.mark.parametrize("field",["schema","id","count","box","partner","seed","proof","encoding","source_digest","missing_digest"])
def test_actual_inventory_structural_mutations(tmp_path,field):
    from simulations import route_b_a2_n6 as m
    data=json.loads((Path(__file__).parents[1]/"results/route_b_a2_n6.json").read_text())
    # Mutate a schema-3 structural control in a TEST TEMPFILE only.
    data["schemaVersion"]=3
    data["sourcePencilDigest"]=load_exact_pencils(RAW)["E"].source_digest
    path=tmp_path/"damaged.json"
    path.write_bytes((json.dumps(data,indent=2,ensure_ascii=False)+"\n").encode())
    assert m.verify_inventory_artifact_structure(path)==data
    if field=="schema": data["schemaVersion"]=1
    elif field=="id": data["loci"][1]["id"]=data["loci"][0]["id"]
    elif field=="count": data["loci"].pop()
    elif field=="box": data["loci"][0]["qPhysicalCSharpBox"]=data["loci"][0]["tBox"]
    elif field=="partner": data["loci"][0]["parityPartnerId"]="missing"
    elif field=="seed": data["loci"][0]["lambdaPhysicalSeed"]["real"]="123456"
    elif field=="proof": data["layerIdentity"]["proofModulus"]="1"
    elif field=="source_digest": data["sourcePencilDigest"]="0"*64
    elif field=="missing_digest": del data["sourcePencilDigest"]
    path.write_bytes((json.dumps(data,indent=2,ensure_ascii=False)+("\n\n" if field=="encoding" else "\n")).encode())
    with pytest.raises(ValueError): m.verify_inventory_artifact_structure(path)


def test_global_resigning_cannot_turn_zero_counts_into_evidence():
    from simulations import route_b_a2_n6 as m
    p=sp.Poly(t**3-2,t)
    boxes=[m.ExactRootBox((Fraction(10),Fraction(11)),(Fraction(0),Fraction(0)),1),
           m.ExactRootBox((Fraction(20),Fraction(21)),(Fraction(1),Fraction(2)),1),
           m.ExactRootBox((Fraction(20),Fraction(21)),(Fraction(-2),Fraction(-1)),1)]
    assert [m.exact_count_in_box(p,b) for b in boxes]==[0,0,0]
    if hasattr(m,"_issue_box_certificate"):
        fake=m._issue_box_certificate(m.EBoxCertificate,p,boxes,["sturm","complex-direct","complex-direct"])
        with pytest.raises(ValueError): m.verify_box_certificate(p,fake,m.EBoxCertificate)
    assert not hasattr(m,"_CERTIFICATE_KEY") and not hasattr(m,"_carrier_seal")
    unsigned=m._assemble_unsigned_box_certificate(m.EBoxCertificate,p,boxes,["sturm","complex-direct","complex-direct"])
    with pytest.raises(ValueError,match="transport integrity"):
        m.verify_box_certificate(p,unsigned,m.EBoxCertificate)


def test_public_n6_writer_rejects_consistent_toy_source(tmp_path):
    m,a2,e,o,raw,layer=phase_split_toy()
    with pytest.raises(ValueError,match="canonical N6 source"):
        m.write_inventory(tmp_path/"not-n6.json",layer,raw,a2,e,o,source_pencil=source_bound_toy()[0])


def test_loader_rejects_consistent_noncanonical_e_o_coefficients(tmp_path,raw):
    for parity in ("rEven","rOdd"):
        raw[parity]["residualInT"][0][0]=str(int(raw[parity]["residualInT"][0][0])+1)
    with pytest.raises(ValueError,match="canonical N6 source"):
        load_exact_pencils(write_raw(tmp_path,raw))


@pytest.mark.parametrize("name",["ROUTE_B_N6_WORKERS","ROUTE_B_N6_PHASE_TIMEOUT_SECONDS","ROUTE_B_N6_PHASE_MEMORY_GIB"])
@pytest.mark.parametrize("value",["0","-1","garbage","1.5",""])
def test_runtime_environment_rejects_invalid_positive_integers(monkeypatch,name,value):
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"parse_runtime_arguments"), "environment defaults are not implemented"
    monkeypatch.setenv(name,value)
    with pytest.raises(SystemExit): m.parse_runtime_arguments(["--residual-json",str(RAW),"--prove-layers"])


def test_runtime_environment_defaults_and_cli_overrides(monkeypatch):
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"parse_runtime_arguments"), "environment defaults are not implemented"
    for name,value in (("WORKERS","4"),("PHASE_TIMEOUT_SECONDS","23"),("PHASE_MEMORY_GIB","2")):
        monkeypatch.setenv("ROUTE_B_N6_"+name,value)
    argv=["--residual-json",str(RAW),"--prove-layers"]
    args=m.parse_runtime_arguments(argv)
    assert (args.workers,args.phase_timeout_seconds,args.phase_memory_gib)==(4,23,2)
    args=m.parse_runtime_arguments(argv+["--workers","2","--phase-timeout-seconds","19","--phase-memory-gib","1"])
    assert (args.workers,args.phase_timeout_seconds,args.phase_memory_gib)==(2,19,1)


def test_count_capability_cannot_cross_runs_or_sign_unsigned_carriers(count_capability):
    from simulations import route_b_a2_n6 as m
    p=sp.Poly(t**3-2,t)
    e=m.isolate_a2_roots(p)
    assert m.verify_box_certificate(p,e,m.EBoxCertificate,capability=count_capability)
    for other in (None,b"",os.urandom(32)):
        with pytest.raises(ValueError): m.verify_box_certificate(p,e,m.EBoxCertificate,capability=other)
        with pytest.raises(ValueError): m.certify_odd_transport(p,e,capability=other)
    unsigned=m._assemble_unsigned_box_certificate(m.EBoxCertificate,p,e.boxes,e.count_methods)
    with pytest.raises(ValueError,match="transport integrity"):
        m.verify_box_certificate(p,unsigned,m.EBoxCertificate,capability=count_capability)


@pytest.mark.parametrize("field",["missing","wrong","residual","at"])
def test_n6_source_pin_recomputed_not_trusted_from_field(monkeypatch,field):
    from dataclasses import replace
    from simulations import route_b_a2_n6 as m
    e=load_exact_pencils(RAW)["E"]
    m._require_canonical_n6_source(e)
    if field=="missing": e=replace(e,source_digest="")
    elif field=="wrong": e=replace(e,source_digest="0"*64)
    else: e=replace(e,**{field:getattr(e,field)+1})
    monkeypatch.setattr(m,"_verify_generic_inventory_source",lambda *_:pytest.fail("bad source reached algebraic proof"))
    with pytest.raises(ValueError,match="canonical N6 source"):
        m.verify_inventory_source(e,None,None)


def test_semantic_source_digest_ignores_json_formatting(tmp_path,raw):
    path=tmp_path/"reformatted.json"
    path.write_text(json.dumps(raw,sort_keys=True,indent=4),encoding="utf-8")
    assert load_exact_pencils(path)==load_exact_pencils(RAW)


def test_phase_environment_caps_are_logged_and_enforced(monkeypatch,capsys):
    from simulations import route_b_a2_n6 as m
    monkeypatch.setenv("ROUTE_B_N6_PHASE_TIMEOUT_SECONDS","7")
    monkeypatch.setenv("ROUTE_B_N6_PHASE_MEMORY_GIB","1")
    assert len(m._bounded_exact_phase("propose",(sp.Poly(t*t+1,t),128)))==2
    log=capsys.readouterr().out
    assert "timeout_seconds=7 memory_bytes=1073741824" in log
    with pytest.raises(MemoryError,match="exceeded 1 bytes"):
        m._bounded_exact_phase("propose",(sp.Poly(t*t+1,t),128),memory_bytes=1)
    with pytest.raises(RuntimeError,match="exceeded 0.001s"):
        m._bounded_exact_phase("propose",(sp.Poly(t*t+1,t),128),timeout_seconds=0.001)


def test_cli_override_ignores_invalid_environment_for_overridden_fields(monkeypatch):
    from simulations import route_b_a2_n6 as m
    for suffix in ("WORKERS","PHASE_TIMEOUT_SECONDS","PHASE_MEMORY_GIB"):
        monkeypatch.setenv("ROUTE_B_N6_"+suffix,"invalid")
    args=m.parse_runtime_arguments(["--residual-json",str(RAW),"--prove-layers",
            "--workers","2","--phase-timeout-seconds","17","--phase-memory-gib","1"])
    assert (args.workers,args.phase_timeout_seconds,args.phase_memory_gib)==(2,17,1)


@pytest.mark.parametrize("override",[False,True])
def test_main_propagates_resolved_runtime_limits_to_phase(monkeypatch,tmp_path,override):
    from simulations import route_b_a2_n6 as m
    for suffix,value in (("WORKERS","4"),("PHASE_TIMEOUT_SECONDS","23"),("PHASE_MEMORY_GIB","2")):
        monkeypatch.setenv("ROUTE_B_N6_"+suffix,value)
    argv=["route_b_a2_n6","--residual-json",str(RAW),"--out",str(tmp_path/"never-written.json")]
    if override: argv += ["--workers","2","--phase-timeout-seconds","17","--phase-memory-gib","1"]
    monkeypatch.setattr(sys,"argv",argv)
    seen=[]
    class StopAfterDispatch(Exception): pass
    def phase(kind,arguments,timeout,memory,**kwargs):
        seen.append((kind,arguments[1],timeout,memory,kwargs["capability"]))
        raise StopAfterDispatch
    monkeypatch.setattr(m,"_bounded_exact_phase",phase)
    with pytest.raises(StopAfterDispatch): m.main()
    assert seen==[("layer",2,17,1024**3,None) if override else ("layer",4,23,2*1024**3,None)]
    assert not (tmp_path/"never-written.json").exists()


def test_verifier_degree_bound_is_row_max_not_a_second_assignment(monkeypatch):
    from simulations import route_b_a2_n6 as m
    import scipy.optimize
    rows,specs,candidate,_=m._fixed_minor_data(load_exact_pencils(RAW)["E"].residual)
    assert candidate==(930,931)
    assert tuple(m._row_degree_bound(rows,c) for c in specs)==(1919,1919)
    monkeypatch.setattr(scipy.optimize,"linear_sum_assignment",lambda *_args,**_kwargs: (_ for _ in ()).throw(AssertionError("assignment reused")))
    assert tuple(m._verification_degree_bound(rows,c) for c in specs)==(1919,1919)


def test_shared_too_small_assignment_cannot_certify_aliasing_minors(monkeypatch):
    from simulations import route_b_a2_n6 as m
    import numpy as np
    import scipy.optimize
    # At candidate nodes 1 and 2 g=1, hence false bound zero aliases both raw
    # minors to constants at EVERY prime. A transposed assignment repeats it.
    g=1+(t-1)*(t-2)
    f=sp.Poly(Lambda**7+g*Lambda+g,Lambda,t,domain=sp.ZZ)
    monkeypatch.setattr(scipy.optimize,"linear_sum_assignment",
                        lambda *_args,**_kwargs:(np.array([],dtype=int),np.array([],dtype=int)))
    with pytest.raises(ValueError,match="fixed-minor verification"):
        m.first_subresultant_linear_exact(f)


def test_subresultant_recheck_rejects_aliased_candidate_certificate(monkeypatch):
    from simulations import route_b_a2_n6 as m
    import numpy as np
    import scipy.optimize
    g=1+(t-1)*(t-2)
    f=sp.Poly(Lambda**7+g*Lambda+g,Lambda,t,domain=sp.ZZ)
    monkeypatch.setattr(scipy.optimize,"linear_sum_assignment",
                        lambda *_args,**_kwargs:(np.array([],dtype=int),np.array([],dtype=int)))
    # Manufacture the old defective certificate, then restore the independent
    # verifier while leaving the common candidate-assignment defect in place.
    with monkeypatch.context() as defective:
        defective.setattr(m,"_verification_degree_bound",lambda *_:0)
        certificate=m.first_subresultant_linear_exact(f)
    assert certificate.degree_bounds==(0,0) and certificate.a_raw.degree()==0
    with pytest.raises(ValueError,match="raw fixed-minor specialization mismatch"):
        m.verify_subresultant_certificate(f,certificate)


def test_verifier_padding_and_truncation_at_independent_degree_1919():
    from simulations import route_b_a2_n6 as m
    assert hasattr(m,"_minor_matches_degree_bound")
    p=101; f=sp.Poly(3*t**873-5,t)
    actual=[0]*(1919-873)+[3]+[0]*872+[96]
    assert len(actual)==1920 and m._minor_matches_degree_bound(actual,f,p,1919)
    assert not m._minor_matches_degree_bound(actual[1:],f,p,1919)
    assert not m._minor_matches_degree_bound(actual+[0],f,p,1919)
    assert not m._minor_matches_degree_bound(actual,sp.Poly(t**1920+f.as_expr(),t),p,1919)


@pytest.mark.parametrize("where",["env","cli"])
def test_oversized_workers_fail_before_any_phase(monkeypatch,where):
    from simulations import route_b_a2_n6 as m
    argv=["route_b_a2_n6","--residual-json",str(RAW),"--prove-layers"]
    if where=="env": monkeypatch.setenv("ROUTE_B_N6_WORKERS","1000000")
    else: argv += ["--workers","1000000"]
    monkeypatch.setattr(sys,"argv",argv)
    monkeypatch.setattr(m,"load_exact_pencils",lambda *_:pytest.fail("oversized workers passed preflight into source loading"))
    monkeypatch.setattr(m,"_bounded_exact_phase",lambda *_args,**_kwargs:pytest.fail("oversized worker pool reached phase"))
    with pytest.raises(SystemExit): m.main()


def test_main_has_no_unbounded_layer_proof_call():
    import ast,inspect
    from simulations import route_b_a2_n6 as m
    tree=ast.parse(inspect.getsource(m.main))
    assert not any(isinstance(node,ast.Call) and isinstance(node.func,ast.Name)
                   and node.func.id in {"prove_layer_identity","_lift_layers"} for node in ast.walk(tree))


def test_initial_layer_proof_obeys_real_child_timeout(monkeypatch,tmp_path,capsys):
    from simulations import route_b_a2_n6 as m
    # No layer mock: the genuine canonical proof must be killed by its own cap.
    monkeypatch.setattr(sys,"argv",["route_b_a2_n6","--residual-json",str(RAW),
        "--out",str(tmp_path/"absent.json"),"--workers","1","--phase-timeout-seconds","1"])
    with pytest.raises(RuntimeError,match="layer exact phase exceeded 1s"):
        m.main()
    log=capsys.readouterr().out
    assert "phase layer:" in log and "timeout_seconds=1" in log
    assert not (tmp_path/"absent.json").exists()


def test_worker_limit_respects_cpus_and_cli_precedence(monkeypatch):
    from simulations import route_b_a2_n6 as m
    monkeypatch.setattr(m.os,"cpu_count",lambda:2)
    monkeypatch.setenv("ROUTE_B_N6_WORKERS","1000000")
    argv=["--residual-json",str(RAW),"--prove-layers"]
    assert m.parse_runtime_arguments(argv+["--workers","2"]).workers==2
    with pytest.raises(SystemExit): m.parse_runtime_arguments(argv+["--workers","3"])
    with pytest.raises(ValueError): m._validate_workers(3)


@pytest.mark.parametrize("entry",["s1","verify_s1","certify_e","certify_o","disc"])
def test_direct_exact_entry_points_reject_oversized_pools(entry):
    from simulations import route_b_a2_n6 as m
    with pytest.raises(ValueError,match="workers"):
        if entry=="s1": m.first_subresultant_linear_exact(None,workers=1000000)
        elif entry=="verify_s1": m.verify_subresultant_certificate(None,None,workers=1000000)
        elif entry=="certify_e": m.certify_a2_boxes(None,None,workers=1000000)
        elif entry=="certify_o": m.certify_odd_transport(None,None,workers=1000000)
        else: next(m._disc_stream(None,iter(()),1000000))
