"""Exact N=6 exported pencils and bounded-CRT discriminant layer identities.

Works directly in t. Block reconstruction and odd-N orbit premises stay outside.
N6_TABLE is a proposal until the coefficient-difference bound proves the identity.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from fractions import Fraction
import json
import hashlib
import hmac
import multiprocessing as mp
import queue as queue_module
import os
from pathlib import Path
import re
import subprocess
import threading
import tempfile
import time

import numpy as np
import sympy as sp

if __package__:
    from . import o2b_gcd_certificate as arithmetic
else:
    import o2b_gcd_certificate as arithmetic

Lambda, t = sp.symbols("Lambda t")
N6_TABLE = {"degD": 926, "valuation": 536, "A1": 124, "A2": 133}
# Semantic digest of the exact C# source fixture, independent of JSON whitespace.
N6_SOURCE_PENCIL_DIGEST="cf1549c54a116132373e481d0ce7a7ea03409c07f75e3f6ab5b9fd0f738dc916"


def _validate_workers(workers):
    # Only 1/2/4/8-worker exact phases have been cost-piloted. Higher concurrency
    # has no measured memory/time justification; also respect available CPUs.
    limit=min(8,os.cpu_count() or 1)
    if type(workers) is not int or not 1<=workers<=limit:
        raise ValueError(f"workers must be between 1 and {limit} (CPU/piloted limit)")
    return workers


@dataclass(frozen=True)
class ParityPencil:
    parity: str
    residual: sp.Poly
    at: sp.Poly
    source_digest: str = ""


@dataclass(frozen=True)
class LayerCertificate:
    deg_d: int
    valuation: int
    a1: sp.Poly
    a2: sp.Poly
    constant: int
    proof_modulus: int
    proof_bound: int
    proof_primes: tuple[int, ...]
    discriminant_lhs: sp.Expr
    discriminant_rhs: sp.Expr


@dataclass(frozen=True)
class SubresultantCertificate:
    s1: sp.Poly
    a_raw: sp.Poly
    b_raw: sp.Poly
    psc1_raw: sp.Poly
    degree_bounds: tuple[int, int, int]
    coefficient_bounds: tuple[int, int, int]
    candidate_primes: tuple[int, ...]
    candidate_modulus: int
    verification_primes: tuple[int, ...]
    verification_modulus: int
    verification_bound: int
    convention: str = (
        "Sylvester1840: descending Lambda columns; F shifts high-to-low then dF shifts "
        "high-to-low; delete last F and last dF row; a/PSC columns 0..r-1; "
        "b columns 0..r-2 plus r; determinant sign unchanged"
    )


@dataclass(frozen=True)
class ExactRootBox:
    real: tuple[Fraction, Fraction]
    imag: tuple[Fraction, Fraction]
    root_count: int | None = None


@dataclass(frozen=True)
class ParityBoxCertificate:
    polynomial_coefficients: tuple
    boxes: tuple
    count_methods: tuple
    conjugation_partners: tuple
    e_partners: tuple = ()
    source_e_seal: str = ""
    seal: str = ""

    def __len__(self): return len(self.boxes)
    def __iter__(self): return iter(self.boxes)
    def __getitem__(self,index): return self.boxes[index]


@dataclass(frozen=True)
class EBoxCertificate(ParityBoxCertificate):
    pass


@dataclass(frozen=True)
class OBoxCertificate(ParityBoxCertificate):
    pass


def _carrier_material(certificate):
    material=(type(certificate).__name__,certificate.polynomial_coefficients,
              certificate.boxes,certificate.count_methods,certificate.conjugation_partners,
              certificate.e_partners,certificate.source_e_seal)
    return repr(material).encode("ascii")


def _assemble_unsigned_box_certificate(cls,poly,boxes,methods,e_certificate=None):
    """Shape assembly only: this helper cannot issue transport credentials."""
    pairs=sorted(zip(boxes,methods),key=lambda pair:(pair[0].real,pair[0].imag))
    boxes,methods=map(tuple,zip(*pairs))
    conjugation=tuple(_find_exact_partner(boxes,_box_conjugate(box)) for box in boxes)
    e_partners=(() if e_certificate is None else
                tuple(_find_exact_partner(e_certificate,_box_neg(box)) for box in boxes))
    result=cls(tuple(int(c) for c in reversed(poly.all_coeffs())),boxes,methods,conjugation,
               e_partners,"" if e_certificate is None else e_certificate.seal)
    return result


def verify_box_certificate(poly,certificate,cls,e_certificate=None,*,capability=None):
    """Check run-local transport integrity and shape, NOT mathematical roots.

    The explicit ephemeral capability belongs to the run orchestrator. Only the
    E/O certifiers sign, inline after executed exact counts; no re-sign API exists.
    This is not a security boundary against code possessing that run capability.
    """
    if type(capability) is not bytes or len(capability)!=32:
        raise ValueError("explicit run-local count transport capability required")
    if type(certificate) is not cls:
        raise ValueError("completed executed-count parity certificate required")
    expected=hmac.new(capability,_carrier_material(certificate),hashlib.sha256).hexdigest()
    if not hmac.compare_digest(certificate.seal,expected):
        raise ValueError("count carrier transport integrity failed")
    if certificate.polynomial_coefficients!=tuple(int(c) for c in reversed(poly.all_coeffs())):
        raise ValueError("certificate polynomial mismatch")
    boxes=certificate.boxes; degree=poly.degree()
    if len(boxes)!=degree or any(type(box.root_count) is not int or box.root_count!=1 for box in boxes):
        raise ValueError("certificate count completeness failed")
    real_count=sum(box.imag==(0,0) for box in boxes)
    if degree==133 and (real_count,degree-real_count)!=(59,74):
        raise ValueError("canonical N6 real/nonreal count split failed")
    if boxes!=tuple(sorted(boxes,key=lambda box:(box.real,box.imag))):
        raise ValueError("certificate order changed")
    if any(not _box_disjoint(a,b) for j,a in enumerate(boxes) for b in boxes[j+1:]):
        raise ValueError("certificate closures overlap")
    expected_methods=tuple(("sturm" if cls is EBoxCertificate else "sturm-pullback")
                           if box.imag==(0,0) else "complex-direct" for box in boxes)
    if certificate.count_methods!=expected_methods:
        raise ValueError("certificate execution provenance shape failed")
    partners=certificate.conjugation_partners
    if len(partners)!=degree: raise ValueError("incomplete conjugation partners")
    for index,partner in enumerate(partners):
        if (type(partner) is not int or not 0<=partner<degree or partners[partner]!=index
                or boxes[partner]!=_box_conjugate(boxes[index])):
            raise ValueError("conjugation partner involution failed")
    if cls is EBoxCertificate:
        if certificate.e_partners or certificate.source_e_seal:
            raise ValueError("invalid E provenance")
    else:
        if e_certificate is None or certificate.source_e_seal!=e_certificate.seal:
            raise ValueError("O certificate lacks its completed E source")
        if sorted(certificate.e_partners)!=list(range(degree)):
            raise ValueError("parity partner map is not total")
        for index,partner in enumerate(certificate.e_partners):
            if boxes[index]!=_box_neg(e_certificate[partner]):
                raise ValueError("exact parity transport failed")
    return True


@dataclass(frozen=True)
class MonomialLocalization:
    """Evaluation-only carrier; never a substitute for the raw S1 certificate."""
    valuation: int
    a_unit: sp.Poly
    b_unit: sp.Poly


def valuation_at_zero(poly):
    if (not isinstance(poly,sp.Poly) or poly.gens!=(t,) or poly.domain!=sp.ZZ
            or poly.is_zero):
        raise ValueError("valuation requires a nonzero ZZ[t] polynomial")
    return min(power[0] for power,coefficient in poly.terms())


def divide_monomial_exact(poly,valuation):
    if type(valuation) is not int or valuation<0:
        raise ValueError("invalid monomial valuation")
    divisor=sp.Poly(t**valuation,t,domain=sp.ZZ)
    quotient,remainder=poly.div(divisor,auto=False)
    if not remainder.is_zero or quotient*divisor!=poly:
        raise ValueError("monomial division has a nonzero remainder")
    return quotient


def localize_monomial_quotient(a2,a_raw,b_raw,expected_valuation=None):
    va,vb=valuation_at_zero(a_raw),valuation_at_zero(b_raw)
    if va!=vb or (expected_valuation is not None and va!=expected_valuation):
        raise ValueError("raw monomial valuations differ from the required common valuation")
    if (not isinstance(a2,sp.Poly) or a2.gens!=(t,) or a2.domain!=sp.ZZ
            or a2.is_zero or a2.eval(0)==0 or sp.gcd(a2,sp.Poly(t,t)).degree()!=0):
        raise ValueError("A2 does not exclude the localization point zero")
    return MonomialLocalization(va,divide_monomial_exact(a_raw,va),
                                divide_monomial_exact(b_raw,vb))


def verify_monomial_localization(a2,a_raw,b_raw,localized,expected_valuation=None):
    expected=localize_monomial_quotient(a2,a_raw,b_raw,expected_valuation)
    if localized!=expected:
        raise ValueError("evaluation polynomials do not equal the exact monomial quotients")
    return True


def localized_lambda_box(localized,box):
    if box.real[0]<=0<=box.real[1] and box.imag[0]<=0<=box.imag[1]:
        raise ValueError("localized t box contains or touches zero")
    return _lambda_box(localized.a_unit,localized.b_unit,box)


@dataclass(frozen=True)
class QuotientRingCertificate:
    """Evaluation-only Bezout witness and canonical polynomial H/D."""
    u: sp.Poly
    v: sp.Poly
    H: sp.Poly
    D: int


def _quotient_source_polynomials(source):
    # Small generic fixtures may use their unmodified raw pair directly. N6's
    # production caller supplies only its verified t**496 localization.
    if isinstance(source,SubresultantCertificate):
        return source.a_raw,source.b_raw
    return source.a_unit,source.b_unit


def verify_quotient_ring_certificate(a2,localized,certificate):
    modulus=a2.set_domain(sp.QQ)
    source_a,source_b=_quotient_source_polynomials(localized)
    a,b=source_a.set_domain(sp.QQ),source_b.set_domain(sp.QQ)
    if (modulus.degree()<1 or certificate.u.gens!=(t,) or certificate.v.gens!=(t,)
            or certificate.u.domain!=sp.QQ or certificate.v.domain!=sp.QQ):
        raise ValueError("invalid quotient-ring Bezout carrier")
    if certificate.u*a+certificate.v*modulus!=sp.Poly(1,t,domain=sp.QQ):
        raise ValueError("quotient-ring Bezout identity fails")
    H,D=certificate.H,certificate.D
    if (not isinstance(H,sp.Poly) or H.gens!=(t,) or H.domain!=sp.ZZ
            or H.degree()>=modulus.degree() or type(D) is not int or D<=0
            or int(sp.igcd(D,abs(int(H.content()))))!=1):
        raise ValueError("noncanonical quotient-ring H/D")
    if not (a*H+D*b).rem(modulus).is_zero:
        raise ValueError("quotient-ring remainder identity fails")
    return True


def build_quotient_ring_certificate(a2,localized):
    modulus=a2.set_domain(sp.QQ)
    source_a,source_b=_quotient_source_polynomials(localized)
    a,b=source_a.set_domain(sp.QQ),source_b.set_domain(sp.QQ)
    u,v,g=sp.gcdex(a,modulus)
    if g!=sp.Poly(1,t,domain=sp.QQ):
        raise ValueError("a_unit is not invertible modulo A2")
    h=(-b*u).rem(modulus)
    denominator,H=h.clear_denoms(convert=True)
    D=int(denominator)
    if D<0: D,H=-D,-H
    common=int(sp.igcd(D,abs(int(H.content()))))
    if common>1:
        D//=common
        H=H.exquo_ground(common)
    certificate=QuotientRingCertificate(u,v,H,D)
    verify_quotient_ring_certificate(a2,localized,certificate)
    return certificate


def quotient_ring_lambda_box(certificate,box):
    if box.real[0]<=0<=box.real[1] and box.imag[0]<=0<=box.imag[1]:
        raise ValueError("quotient-ring t box contains or touches zero")
    if type(certificate.D) is not int or certificate.D<=0:
        raise ValueError("quotient-ring denominator must be a positive integer")
    value=_poly_box(certificate.H,box)
    scale=2**96
    def outward(bounds):
        lo,hi=(endpoint/certificate.D for endpoint in bounds)
        return (Fraction((lo*scale).__floor__(),scale),
                Fraction((hi*scale).__ceil__(),scale))
    return ExactRootBox(outward(value.real),outward(value.imag))


@dataclass(frozen=True)
class ProposedBall:
    real_midpoint: str
    imag_midpoint: str
    real_radius: str
    imag_radius: str
    precision_bits: int


def _require_flint():
    try:
        import flint
    except ImportError as exc:
        raise RuntimeError("python-flint is required for proposals; no artifact written") from exc
    return flint


def propose_a2_boxes_flint(poly, precision_bits=128):
    """Numerical suggestions only. No root count or certificate is returned."""
    flint=_require_flint()
    with flint.ctx.workprec(precision_bits):
        # Retain exact input coefficients; the returned ACB balls are still
        # merely proposals and their multiplicities are deliberately ignored.
        values=[z for z,_multiplicity in flint.fmpz_poly(
            [int(c) for c in reversed(poly.all_coeffs())]).complex_roots()]
        def decimal(value):
            rational=Fraction(str(value.fmpq()))
            return _decimal_in_interval((rational,rational))[0]
        return [ProposedBall(decimal(z.real.mid()),decimal(z.imag.mid()),
                             decimal(z.real.rad()),decimal(z.imag.rad()),precision_bits)
                for z in values]


def _ball_bounds(ball):
    if not isinstance(ball,ProposedBall) or type(ball.precision_bits) is not int or ball.precision_bits<2:
        raise ValueError("invalid proposal")
    r,i,rr,ir=map(Fraction,(ball.real_midpoint,ball.imag_midpoint,ball.real_radius,ball.imag_radius))
    if rr<0 or ir<0: raise ValueError("negative proposal radius")
    return (r-rr,r+rr),(i-ir,i+ir)


def outward_dyadic_box(ball,dyadic_bits=64):
    if type(dyadic_bits) is not int or dyadic_bits<1:
        raise ValueError("invalid dyadic precision")
    scale=2**dyadic_bits
    def outward(bounds):
        lo,hi=bounds
        return (Fraction((lo*scale).__floor__()-1,scale),
                Fraction((hi*scale).__ceil__()+1,scale))
    real,imag=_ball_bounds(ball)
    return ExactRootBox(outward(real),(Fraction(0),Fraction(0)) if imag==(0,0) else outward(imag))


def _check_containment(ball,box):
    real,imag=_ball_bounds(ball)
    if not all(got[0]<=want[0]<=want[1]<=got[1] for got,want in ((box.real,real),(box.imag,imag))):
        raise ValueError("dyadic box does not contain proposal ball")


@dataclass(frozen=True)
class ExactSturmChain:
    # Ascending exact (numerator, positive denominator) coefficients are both
    # immutable and multiprocessing-safe; FLINT polynomial objects are not.
    members: tuple
    quotients: tuple
    scales: tuple


def _pack_sturm_poly(poly):
    return tuple((int(c.numerator),int(c.denominator)) for c in poly)


def _unpack_sturm_poly(coefficients):
    flint=_require_flint()
    return flint.fmpq_poly([flint.fmpq(n,d) for n,d in coefficients])


def _sturm_input(poly):
    if not isinstance(poly,sp.Poly) or poly.gens!=(t,) or poly.domain!=sp.ZZ or poly.degree()<1:
        raise ValueError("Sturm input must be nonconstant ZZ[t]")
    result=_require_flint().fmpq_poly([int(c) for c in reversed(poly.all_coeffs())])
    back=sp.Poly.from_list([sp.Rational(int(c.numerator),int(c.denominator))
                           for c in reversed(list(result))],t,domain=sp.ZZ)
    if back!=poly: raise ValueError("ascending FLINT coefficient round-trip failed")
    return result


def build_exact_sturm_chain(poly):
    first=_sturm_input(poly)
    previous,current=first,first.derivative()
    members=[_pack_sturm_poly(previous),_pack_sturm_poly(current)]; quotients=[]; scales=[]
    while True:
        quotient,remainder=divmod(previous,current)
        following=-remainder
        if previous!=quotient*current-following:
            raise ValueError("raw Sturm recurrence failed")
        quotients.append(_pack_sturm_poly(quotient))
        if not following:
            if current.degree()!=0: raise ValueError("Sturm terminal is not a nonzero constant")
            scales.append((1,1))
            break
        if following.degree()>=current.degree(): raise ValueError("Sturm degree did not decrease")
        # Only after checking the raw recurrence: a canonical POSITIVE scale
        # clears denominators/content without changing any Sturm sign.
        scale=_require_flint().fmpq(following.denom(),abs(following.numer().content()))
        if scale<=0: raise ValueError("Sturm scale must be positive")
        scales.append((int(scale.numerator),int(scale.denominator)))
        following=scale*following
        members.append(_pack_sturm_poly(following))
        previous,current=current,following
    chain=ExactSturmChain(tuple(members),tuple(quotients),tuple(scales))
    verify_exact_sturm_chain(poly,chain)
    return chain


def verify_exact_sturm_chain(poly,chain):
    first=_sturm_input(poly)
    if (not isinstance(chain,ExactSturmChain) or len(chain.members)<2
            or len(chain.quotients)!=len(chain.members)-1 or len(chain.scales)!=len(chain.quotients)):
        raise ValueError("invalid Sturm chain shape")
    members=[_unpack_sturm_poly(row) for row in chain.members]
    if members[0]!=first or members[1]!=first.derivative():
        raise ValueError("Sturm input or derivative mismatch")
    if not members[-1] or members[-1].degree()!=0:
        raise ValueError("Sturm terminal must be a nonzero constant")
    for i,quotient in enumerate(chain.quotients,1):
        previous,current=members[i-1:i+1]
        if not current or current.degree()>=previous.degree():
            raise ValueError("Sturm degrees must strictly decrease")
        following=members[i+1] if i+1<len(members) else _require_flint().fmpq_poly([])
        numerator,denominator=chain.scales[i-1]
        if numerator<=0 or denominator<=0: raise ValueError("Sturm scale witness must be positive")
        scale=_require_flint().fmpq(numerator,denominator)
        raw_following=following/scale
        if previous!=_unpack_sturm_poly(quotient)*current-raw_following:
            raise ValueError("Sturm quotient/remainder recurrence mismatch")
        expected=(_require_flint().fmpq(raw_following.denom(),abs(raw_following.numer().content()))
                  if raw_following else _require_flint().fmpq(1))
        if scale!=expected: raise ValueError("Sturm scale witness is not canonical")
    return True


def sturm_variations_at(chain,endpoint):
    endpoint=Fraction(endpoint)
    value=_require_flint().fmpq(endpoint.numerator,endpoint.denominator)
    signs=[]
    for row in chain.members:
        evaluated=_unpack_sturm_poly(row)(value)
        if not evaluated: raise ValueError("zero Sturm member at endpoint boundary")
        signs.append(1 if evaluated>0 else -1)
    return sum(left!=right for left,right in zip(signs,signs[1:]))


def exact_count_in_real_interval(poly,lo,hi):
    """Execute a local count using a reusable, exactly checked Sturm chain."""
    lo,hi=map(Fraction,(lo,hi))
    if lo>=hi: raise ValueError("degenerate or reversed real interval")
    chain=poly if isinstance(poly,ExactSturmChain) else build_exact_sturm_chain(poly)
    count=sturm_variations_at(chain,lo)-sturm_variations_at(chain,hi)
    if count<0 or count>len(chain.members)-1:
        raise ValueError("invalid negative or excessive Sturm variation difference")
    return count


def exact_count_in_box(poly,box):
    """Execute an exact Sturm/Cauchy-index count, rejecting boundary roots."""
    r0,r1=map(sp.Rational,box.real); i0,i1=map(sp.Rational,box.imag)
    if r0>r1 or i0>i1: raise ValueError("reversed root box")
    if i0==i1==0:
        return exact_count_in_real_interval(poly,r0,r1)
    if r0==r1 or i0==i1: raise ValueError("degenerate complex root box")
    corners=(r0+sp.I*i0,r1+sp.I*i0,r1+sp.I*i1,r0+sp.I*i1)
    # On each edge p(start+(end-start)*u), boundary zeros are common
    # real roots of its exact real and imaginary polynomials.
    for start,end in zip(corners,corners[1:]+corners[:1]):
        edge=sp.Poly(poly.as_expr(),t,domain=sp.QQ_I).compose(
            sp.Poly(start+(end-start)*t,t,domain=sp.QQ_I))
        real=sp.Poly.from_list([sp.re(c) for c in edge.all_coeffs()],t,domain=sp.QQ)
        imag=sp.Poly.from_list([sp.im(c) for c in edge.all_coeffs()],t,domain=sp.QQ)
        common=sp.gcd(real,imag)
        if common.is_zero or (common.degree()>0 and common.count_roots(0,1)>0):
            raise ValueError("root on box boundary")
    return int(poly.count_roots(corners[0],corners[2]))


def _count_job(args):
    started=time.perf_counter()
    count=exact_count_in_box(*args)
    kind="real" if args[1].imag==(0,0) else "complex"
    print(f"local_{kind}_count={count} elapsed={time.perf_counter()-started:.3f}s",flush=True)
    return count


def certify_a2_boxes(poly,proposals,dyadic_bits=64,workers=1,sturm_chain=None,*,capability=None):
    """Local exact counts only; FLINT supplies no evidence or multiplicities."""
    _validate_workers(workers)
    if (not isinstance(poly,sp.Poly) or poly.gens!=(t,) or poly.domain!=sp.ZZ
            or poly.degree()<1 or sp.gcd(poly,poly.diff()).degree()!=0):
        raise ValueError("A2 must be a squarefree nonconstant ZZ[t] polynomial")
    if len(proposals)!=poly.degree() or len(set(proposals))!=len(proposals):
        raise ValueError("missing or duplicate proposal")
    real=[]; upper=[]; lower=[]
    for ball in proposals:
        box=outward_dyadic_box(ball,dyadic_bits)
        _check_containment(ball,box)
        if box.imag==(0,0): real.append(box)
        elif box.imag[0]>0: upper.append(box)
        elif box.imag[1]<0: lower.append(box)
        else: raise ValueError("unresolved proposal classification touches real axis")
    if len(lower)!=len(upper): raise ValueError("unbalanced nonreal proposal inventory")
    supplied=real+upper+lower
    if any(not _box_disjoint(a,b) for j,a in enumerate(supplied) for b in supplied[j+1:]):
        raise ValueError("proposal closures overlap")
    proposed=real+upper+[_box_conjugate(box) for box in upper]
    if len(proposed)!=poly.degree(): raise ValueError("proposal inventory does not exhaust degree")
    if any(not _box_disjoint(a,b) for j,a in enumerate(proposed) for b in proposed[j+1:]):
        raise ValueError("proposal closures overlap")
    print(f"local certification: proposed_real={len(real)} proposed_upper={len(upper)} "
          f"dyadic_depth={dyadic_bits} workers={workers}; global_isolation=disabled",flush=True)
    if real:
        sturm_chain=build_exact_sturm_chain(poly) if sturm_chain is None else sturm_chain
        verify_exact_sturm_chain(poly,sturm_chain)
    nonreal=upper+[_box_conjugate(box) for box in upper]
    jobs=[(sturm_chain,box) for box in real]+[(poly,box) for box in nonreal]
    if workers>1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            counts=list(executor.map(_count_job,jobs))
    else:
        counts=[_count_job(job) for job in jobs]
    if any(count!=1 for count in counts): raise ValueError("local exact root count is not one")
    result=[ExactRootBox(box.real,box.imag,count) for box,count in zip(real+nonreal,counts)]
    if sum(counts)!=poly.degree():
        raise ValueError("executed exact counts do not exhaust degree")
    certificate=_assemble_unsigned_box_certificate(EBoxCertificate,poly,result,
                                       ["sturm"]*len(real)+["complex-direct"]*len(nonreal))
    if type(capability) is not bytes or len(capability)!=32:
        raise ValueError("explicit run-local count transport capability required")
    certificate=replace(certificate,seal=hmac.new(capability,_carrier_material(certificate),hashlib.sha256).hexdigest())
    verify_box_certificate(poly,certificate,EBoxCertificate,capability=capability)
    print(f"E certificate: real={len(real)} nonreal={len(nonreal)} executed_sum={sum(counts)}",flush=True)
    return certificate


def certify_odd_transport(a2,e_certificate,workers=1,sturm_chain=None,*,capability=None):
    """Own only O counts; derive every box by the exact parameter involution."""
    _validate_workers(workers)
    verify_box_certificate(a2,e_certificate,EBoxCertificate,capability=capability)
    odd=sp.Poly(a2.as_expr().subs(t,-t),t,domain=sp.ZZ)
    boxes=tuple(sorted((_box_neg(box) for box in e_certificate),key=lambda box:(box.real,box.imag)))
    sturm_chain=build_exact_sturm_chain(a2) if sturm_chain is None else sturm_chain
    verify_exact_sturm_chain(a2,sturm_chain)
    jobs=[(sturm_chain,_box_neg(box)) if box.imag==(0,0) else (odd,box) for box in boxes]
    if workers>1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            counts=list(executor.map(_count_job,jobs))
    else: counts=[_count_job(job) for job in jobs]
    if any(count!=1 for count in counts) or sum(counts)!=a2.degree():
        raise ValueError("O transport exact counts fail")
    certified=[ExactRootBox(box.real,box.imag,count) for box,count in zip(boxes,counts)]
    methods=["sturm-pullback" if box.imag==(0,0) else "complex-direct" for box in boxes]
    result=_assemble_unsigned_box_certificate(OBoxCertificate,odd,certified,methods,e_certificate)
    result=replace(result,seal=hmac.new(capability,_carrier_material(result),hashlib.sha256).hexdigest())
    verify_box_certificate(odd,result,OBoxCertificate,e_certificate,capability=capability)
    real=sum(box.imag==(0,0) for box in boxes)
    print(f"O certificate: real={real} nonreal={len(boxes)-real} executed_sum={sum(counts)}",flush=True)
    return result


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _keys(value, expected):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError(f"expected exactly the fields {sorted(expected)}")


def _parse_poly(rows, degree):
    # Each exporter row is individually trimmed; rectangular padding is invalid.
    if not isinstance(rows, list) or len(rows) != degree + 1:
        raise ValueError("wrong Lambda degree / coefficient row count")
    terms = {}
    for i, row in enumerate(rows):
        if not isinstance(row, list) or not row or len(row) > degree-i+1:
            raise ValueError("malformed trimmed coefficient row")
        for j, value in enumerate(row):
            if not isinstance(value, str) or re.fullmatch(r"0|-?[1-9][0-9]*", value) is None:
                raise ValueError("coefficients must be real decimal integer strings")
            terms[i, j] = int(value)
        if len(row) > 1 and row[-1] == "0":
            raise ValueError("noncanonical trailing zero in trimmed row")
    poly = sp.Poly.from_dict(terms, (Lambda, t), domain=sp.ZZ)
    if poly.degree(Lambda) != degree or rows[-1] != ["1"]:
        raise ValueError("polynomial must have the specified degree and be monic in Lambda")
    return poly


def load_exact_pencils(path):
    """Parse schema 1 strictly; check independently exported E/O transports."""
    data = json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=_unique_object)
    _keys(data, {"schemaVersion", "n", "lambdaConvention", "parameterConvention",
                 "coefficientOrder", "rEven", "rOdd"})
    if (type(data["schemaVersion"]) is not int or data["schemaVersion"] != 1
            or type(data["n"]) is not int or data["n"] != 6
            or data["lambdaConvention"] != "Lambda=2*lambda"
            or data["parameterConvention"] != "t=i*qCSharp"
            or data["coefficientOrder"] != "lambda-lowest-first,t-lowest-first"):
        raise ValueError("unsupported schema, N, conventions or coefficient order")
    result = {}
    for key, parity, odd in (("rEven", "E", False), ("rOdd", "O", True)):
        obj = data[key]
        _keys(obj, {"rOdd", "sectorDimension", "residualInT", "atFactorInT"})
        if obj["rOdd"] is not odd or type(obj["sectorDimension"]) is not int or obj["sectorDimension"] != 45:
            raise ValueError("parity or sector dimension mismatch")
        result[parity] = ParityPencil(parity, _parse_poly(obj["residualInT"],32),
                                    _parse_poly(obj["atFactorInT"],13))
    for field in ("residual", "at"):
        e, o = getattr(result["E"],field), getattr(result["O"],field)
        expected = sp.Poly.from_dict({(i,j): (-1)**j*c for (i,j),c in e.terms()},
                                     (Lambda,t), domain=sp.ZZ)
        if o != expected:
            raise ValueError(f"E/O transport fails for {field}")
    digest=_source_pencil_digest(result["E"])
    if digest!=N6_SOURCE_PENCIL_DIGEST:
        raise ValueError("export differs from the canonical N6 source pencil")
    return {parity:replace(pencil,source_digest=digest) for parity,pencil in result.items()}


def _source_pencil_digest(pencil):
    """Canonical semantic encoding: exact sparse terms plus model conventions."""
    if not isinstance(pencil,ParityPencil) or pencil.parity!="E":
        raise ValueError("canonical N6 source requires the E pencil")
    material={"n":6,"sectorDimension":45,"parameterConvention":"t=i*qCSharp",
              "lambdaConvention":"Lambda=2*lambda","coefficientOrder":"lambda-lowest-first,t-lowest-first",
              "parities":{}}
    for field,degree in (("residual",32),("at",13)):
        poly=getattr(pencil,field)
        if (not isinstance(poly,sp.Poly) or poly.gens!=(Lambda,t) or poly.domain!=sp.ZZ
                or poly.degree(Lambda)!=degree):
            raise ValueError("canonical N6 source polynomial dimensions differ")
    for parity,sign in (("E",1),("O",-1)):
        material["parities"][parity]={field:[[i,j,str(sign**j*c)] for (i,j),c in getattr(pencil,field).terms()]
                                        for field in ("residual","at")}
    return hashlib.sha256(json.dumps(material,sort_keys=True,separators=(",",":")).encode("ascii")).hexdigest()


def _require_canonical_n6_source(pencil):
    if (not isinstance(pencil,ParityPencil) or pencil.source_digest!=N6_SOURCE_PENCIL_DIGEST
            or _source_pencil_digest(pencil)!=N6_SOURCE_PENCIL_DIGEST):
        raise ValueError("writer requires the canonical N6 source pencil digest")


def _matrix(f):
    if not isinstance(f, sp.Poly) or f.gens != (Lambda,t) or f.domain != sp.ZZ:
        raise ValueError("expected Poly(Lambda,t) over ZZ")
    m, d = f.degree(Lambda), max(0,f.degree(t))
    if m < 2 or any(i == m and (j != 0 or c != 1) for (i,j),c in f.terms()):
        raise ValueError("expected a monic polynomial of Lambda degree at least two")
    rows = [[0]*(d+1) for _ in range(m+1)]
    for (i,j),c in f.terms():
        rows[i][j] = int(c)
    return rows


def discriminant_degree_bound(f):
    """Sylvester support bound independent of proposed discriminant degree.

    Each F row shifted by Lambda^s has total degree <= D+s, and each F'
    row <= D-1+s. A determinant term uses all columns 0..2m-2 once.
    Subtracting their Lambda degrees gives (2m-1)D-m^2. The ordinary
    row bound (2m-1)*deg_t(F) is also valid; use the smaller one.
    """
    _matrix(f)
    m = f.degree(Lambda)
    d = max(i+j for (i,j),c in f.terms())
    return min((2*m-1)*max(0,f.degree(t)), (2*m-1)*d-m*m)


def _det_bound(f):
    rows = _matrix(f)
    m = len(rows)-1
    norms = [sum(abs(c) for c in row) for row in rows]
    # Polynomial l1 is submultiplicative: det <= permanent <= row-sum product.
    return sum(norms)**(m-1)*sum(i*v for i,v in enumerate(norms))**m


def coefficient_difference_bound(f, constant, a1, a2):
    """Bound EVERY coefficient of disc(F)-C*t^v*A1*A2^2."""
    return (_det_bound(f)+abs(int(constant))*sum(abs(int(c)) for c in a1.all_coeffs())
            *sum(abs(int(c)) for c in a2.all_coeffs())**2)


def _disc_stream(f, primes, workers):
    _validate_workers(workers)
    m = f.degree(Lambda)
    sign = (-1)**(m*(m-1)//2)
    nodes = list(range(discriminant_degree_bound(f)+1))
    stream = arithmetic.disc_stream(_matrix(f), nodes, primes, workers)
    try:
        for p, resultant in stream:
            yield p, (resultant*sign)%p
    finally:
        stream.close()


def _canonical_layer(poly):
    if (not isinstance(poly,sp.Poly) or poly.gens != (t,) or poly.domain != sp.ZZ
            or poly.is_zero or poly.LC() <= 0 or poly.content() != 1 or poly.eval(0) == 0):
        raise ValueError("layers must be primitive, positive-leading ZZ[t], nonzero at zero")


def certify_layer_identity(f, valuation, a1, a2, constant=None, workers=1):
    """Public fail-closed judge: finite-prime agreement is never a certificate.

    Missing C is proposed by centered CRT of leading ratios past 2*B_D.
    Final equality with modulus > 2*(B_D+B_R) proves the integer identity.
    Every proof prime preserves degree/valuation. Squarefree coprime reduction
    with leading coefficients intact proves shape over Q by Gauss's lemma.
    """
    _matrix(f)
    if type(valuation) is not int or valuation < 0 or type(workers) is not int or workers < 1:
        raise ValueError("invalid valuation or worker count")
    _canonical_layer(a1)
    _canonical_layer(a2)
    if constant is not None and (type(constant) is not int or constant == 0):
        raise ValueError("C must be a nonzero integer")
    degree = valuation+a1.degree()+2*a2.degree()
    if degree > discriminant_degree_bound(f):
        raise ValueError("proposed identity exceeds discriminant support")
    lc = int(a1.LC()*a2.LC()**2)
    bound_d = _det_bound(f)
    pool, pending = [], []
    modulus, rc, mc = 1, 0, 1
    discriminant_residues = [0]*(degree+1)
    shape = False
    started = time.perf_counter()
    stream = _disc_stream(f, arithmetic.prime_stream(), workers)
    try:
        for p, disc in stream:
            if p in pool:
                raise ValueError("duplicate proof prime")
            if len(disc)-1 != degree or arithmetic.valuation_at_zero(disc) != valuation or lc%p == 0:
                raise ValueError(f"proof prime {p} does not preserve degree/valuation")
            pool.append(p)
            pending.append((p,disc))
            discriminant_residues = [arithmetic.crt_pair(old,modulus,int(new),p)
                                     for old,new in zip(discriminant_residues,disc)]
            modulus *= p
            if constant is None:
                ratio = int(disc[0])*pow(lc%p,-1,p)%p
                rc = arithmetic.crt_pair(rc,mc,ratio,p)
                mc *= p
                if mc <= 2*bound_d:
                    if len(pool)%32 == 0:
                        print(f"proof C: {len(pool)} primes, {modulus.bit_length()} bits, {time.perf_counter()-started:.1f}s",flush=True)
                    continue
                constant = rc if rc <= mc//2 else rc-mc
                if constant == 0:
                    raise ValueError("zero reconstructed constant")
            for q, dp in pending:
                x = arithmetic.sympy_poly_desc_modp(a1,q)
                y = arithmetic.sympy_poly_desc_modp(a2,q)
                rhs = arithmetic.polymul(arithmetic.polymul(y,y,q),x,q)
                rhs = np.concatenate([(rhs*(constant%q))%q,np.zeros(valuation,dtype=np.int64)])
                if not np.array_equal(arithmetic.polytrim(rhs),dp):
                    raise ValueError(f"integer layer identity fails modulo {q}")
                if not shape:
                    shape = (len(arithmetic.polygcd(x,arithmetic.polyderiv(x,q),q)) == 1
                             and len(arithmetic.polygcd(y,arithmetic.polyderiv(y,q),q)) == 1
                             and len(arithmetic.polygcd(x,y,q)) == 1)
            pending.clear()
            bound = coefficient_difference_bound(f,constant,a1,a2)
            if modulus > 2*bound:
                if not shape:
                    raise ValueError("no squarefree coprime layer shape in the completed proof pool")
                # These coefficients come from the determinant residues, never RHS.
                # modulus > 2*bound >= 2*B_D proves the centered lift over ZZ.
                lhs = sp.Poly.from_list([c if c <= modulus//2 else c-modulus
                                         for c in discriminant_residues],t,domain=sp.ZZ)
                rhs = constant*sp.Poly(t**valuation,t)*a1*a2**2
                if lhs != rhs:
                    raise ValueError("independently reconstructed discriminant differs from RHS")
                return LayerCertificate(degree,valuation,a1,a2,constant,modulus,bound,
                                        tuple(pool),lhs.as_expr(),rhs.as_expr())
    finally:
        stream.close()
    raise ValueError("prime stream exhausted before the integer proof")


def _lift_layers(f, workers):
    table = {"v":N6_TABLE["valuation"], "A1":N6_TABLE["A1"], "A2":N6_TABLE["A2"]}
    residues, modulus, candidate, stable = None,1,None,0
    started = time.perf_counter()
    stream = _disc_stream(f,arithmetic.prime_stream(),workers)
    try:
        for used,(p,disc) in enumerate(stream,1):
            if len(disc)-1 != N6_TABLE["degD"]:
                raise ValueError("modular discriminant degree differs from N6 proposal")
            try:
                first,second,_ = arithmetic.layers_from_disc(disc,p,table)
            except AssertionError as exc:
                raise ValueError(f"N6 layer proposal fails at prime {p}") from exc
            incoming = [list(map(int,first)),list(map(int,second))]
            if residues is None:
                residues = incoming
            else:
                residues = [[arithmetic.crt_pair(a,modulus,b,p) for a,b in zip(old,new)]
                            for old,new in zip(residues,incoming)]
            modulus *= p
            if used%8:
                continue
            print(f"lift: {used} primes, {modulus.bit_length()} bits, {time.perf_counter()-started:.1f}s",flush=True)
            proposed = []
            for layer in residues:
                rationals = [arithmetic.rational_reconstruct(c,modulus) for c in layer]
                if any(c is None for c in rationals):
                    break
                poly = sp.Poly.from_list([sp.Rational(a,b) for a,b in rationals],t,domain=sp.QQ)
                poly = poly.clear_denoms(convert=True)[1].primitive()[1]
                proposed.append(poly if poly.LC()>0 else -poly)
            if len(proposed) != 2:
                candidate,stable = None,0
            elif proposed == candidate:
                stable += 1
                if stable == 2:
                    return tuple(proposed)
            else:
                candidate,stable = proposed,0
    finally:
        stream.close()


def prove_layer_identity(parity_pencil, workers=1):
    """Lift E's proposed layers, then prove the full integer identity."""
    if parity_pencil.parity != "E":
        raise ValueError("prove E; O is transported only by the strict loader")
    a1,a2 = _lift_layers(parity_pencil.residual,workers)
    if (a1.degree(),a2.degree()) != (124,133):
        raise ValueError("lifted layer degree mismatch")
    return certify_layer_identity(parity_pencil.residual,536,a1,a2,workers=workers)


def _fixed_s1_rows(f):
    """Return the amended fixed Sylvester-1840 first-subresultant row carrier."""
    _matrix(f)
    m = f.degree(Lambda)
    fc = [sp.Poly(c,t,domain=sp.ZZ) for c in sp.Poly(f.as_expr(),Lambda).all_coeffs()]
    derivative = sp.Poly(sp.diff(f.as_expr(),Lambda),Lambda)
    gc = [sp.Poly(c,t,domain=sp.ZZ) for c in derivative.all_coeffs()]
    zero = sp.Poly(0,t,domain=sp.ZZ)
    f_rows = [[zero]*shift+fc+[zero]*(m-2-shift) for shift in range(m-1)]
    g_rows = [[zero]*shift+gc+[zero]*(m-1-shift) for shift in range(m)]
    # SymPy's Sylvester-1840 ordering, with the final row of each family deleted.
    rows = f_rows[:-1]+g_rows[:-1]
    if len(rows) != 2*m-3 or any(len(row) != 2*m-1 for row in rows):
        raise ValueError("invalid fixed first-subresultant matrix shape")
    return rows


def _minor_columns(size):
    return tuple(range(size)), tuple(range(size-1))+(size,)


def _assignment_degree(rows, columns):
    from scipy.optimize import linear_sum_assignment
    impossible = -10**9
    weights = np.array([[impossible if rows[i][j].is_zero else rows[i][j].degree()
                         for j in columns] for i in range(len(rows))],dtype=np.int64)
    rr,cc = linear_sum_assignment(weights,maximize=True)
    if any(weights[i,j] == impossible for i,j in zip(rr,cc)):
        return 0  # the defining minor is structurally the zero polynomial
    return int(weights[rr,cc].sum())


def _row_degree_bound(rows, columns):
    """Independent determinant-degree bound: sum of each row's maximum degree."""
    total = 0
    for row in rows:
        present = [row[j].degree() for j in columns if not row[j].is_zero]
        if not present:
            raise ValueError("zero row in fixed minor")
        total += max(present)
    return total


def _verification_degree_bound(rows,columns):
    """Every determinant product takes one entry per row: sum row maxima.

    No assignment optimizer is used here, even transposed. The looser bound
    deliberately keeps a candidate assignment defect out of the verifier.
    """
    return _row_degree_bound(rows,columns)


def _minor_matches_degree_bound(actual,poly,prime,bound):
    """Compare the complete bound+1 coefficient vector; never truncate a tail."""
    if poly.degree()>bound or len(actual)!=bound+1:
        return False
    coefficients=[int(c)%prime for c in poly.all_coeffs()]
    return list(actual)==[0]*(bound+1-len(coefficients))+coefficients


def _minor_coefficient_bound(rows, columns):
    """Coefficient l1(det) <= permanent(row entry l1) <= product row l1."""
    result = 1
    for row in rows:
        norm = sum(sum(abs(int(c)) for c in row[j].all_coeffs()) for j in columns)
        if norm == 0:
            return 0  # determinant is structurally zero
        result *= norm
    return result


def _fixed_minor_data(f):
    rows = _fixed_s1_rows(f)
    specs = _minor_columns(len(rows))
    degree_bounds = tuple(_assignment_degree(rows,columns) for columns in specs)
    coefficient_bounds = tuple(_minor_coefficient_bound(rows,columns) for columns in specs)
    return rows,specs,degree_bounds,coefficient_bounds


def _poly_matrix(f):
    m = f.degree(Lambda)
    d = max(0,f.degree(t))
    result = [[0]*(d+1) for _ in range(m+1)]
    for (i,j),coefficient in f.terms():
        result[i][j] = int(coefficient)
    return result


def _numeric_fixed_rows(fdesc,p):
    m = len(fdesc)-1
    fdesc=[int(value)%p for value in fdesc]
    gdesc = [fdesc[i]*(m-i)%p for i in range(m)]
    return ([([0]*shift+fdesc+[0]*(m-2-shift)) for shift in range(m-2)]
            +[([0]*shift+gdesc+[0]*(m-1-shift)) for shift in range(m-1)])


def _fixed_minor_values(coefficient_matrix,nodes,p,method):
    values = arithmetic.eval_nodes_modp(coefficient_matrix,nodes,p)
    columns = _minor_columns(2*len(coefficient_matrix)-5)
    output = [[],[]]
    if method == "flint":
        from flint import nmod_mat
    for value in values:
        rows = _numeric_fixed_rows(value[::-1].copy()%p,p)
        for index,chosen in enumerate(columns):
            matrix = [[row[j] for j in chosen] for row in rows]
            if method == "flint":
                determinant = int(nmod_mat(matrix,p).det())
            elif method == "numpy":
                determinant = arithmetic.det_modp(np.array(matrix,dtype=np.int64),p)
            else:
                raise ValueError("unknown fixed-minor implementation")
            output[index].append(determinant)
    return output


def _minor_polys_modp(f,p,degree_bounds,method):
    if p <= max(degree_bounds)+3:
        raise ValueError("prime is too small for distinct interpolation nodes")
    coefficient_matrix = _poly_matrix(f)
    needed = max(degree_bounds)+2  # strictly more nodes than every asserted bound
    candidates = list(range(1,needed+8))
    values = arithmetic.eval_nodes_modp(coefficient_matrix,candidates,p)
    nodes=[]
    raw=[[],[]]
    for node,value in zip(candidates,values):
        rows = _numeric_fixed_rows(value[::-1].copy()%p,p)
        columns = _minor_columns(len(rows))
        probe=[]
        if method == "flint":
            from flint import nmod_mat
            for chosen in columns:
                probe.append(int(nmod_mat([[row[j] for j in chosen] for row in rows],p).det()))
        else:
            for chosen in columns:
                probe.append(arithmetic.det_modp(
                    np.array([[row[j] for j in chosen] for row in rows],dtype=np.int64),p))
        # A simultaneous zero is a specialized higher-gcd index; never use it as a node.
        if probe != [0,0]:
            nodes.append(node)
            raw[0].append(probe[0]); raw[1].append(probe[1])
        if len(nodes) == needed:
            break
    if len(nodes) != needed:
        raise ValueError("too many rejected subresultant specializations")
    result=[]
    for index,ys in enumerate(raw):
        if any(nodes[i]+1!=nodes[i+1] for i in range(len(nodes)-1)):
            polynomial = arithmetic.interp_modp(nodes,ys,p)
        else:
            polynomial = _interp_consecutive_modp(nodes[0],ys,p)
        if len(polynomial)-1 > degree_bounds[index]:
            raise ValueError("interpolation contradicts fixed-minor degree bound")
        result.append([0]*(degree_bounds[index]+1-len(polynomial))+list(map(int,polynomial)))
    return tuple(result)


def _minor_prime_job(arguments):
    f,p,bounds,method=arguments
    return p,_minor_polys_modp(f,p,bounds,method)


def _interp_consecutive_modp(start,ys,p):
    """Fast exact Newton interpolation at start,start+1,...; descending output."""
    work=np.array(ys,dtype=np.int64)%p
    n=len(work)
    newton=np.zeros(n,dtype=np.int64)
    factorial=1
    for order in range(n):
        if order:
            factorial=factorial*order%p
        newton[order]=int(work[0])*pow(factorial,-1,p)%p
        work=(work[1:]-work[:-1])%p
    out=np.zeros(n,dtype=np.int64)
    basis=np.zeros(n,dtype=np.int64); basis[0]=1
    degree=0
    for order,coefficient in enumerate(newton):
        out[:degree+1]=(out[:degree+1]+coefficient*basis[:degree+1])%p
        if order==n-1: break
        shifted=np.zeros(n,dtype=np.int64)
        shifted[1:degree+2]=basis[:degree+1]
        shifted[:degree+1]=(shifted[:degree+1]-(start+order)*basis[:degree+1])%p
        basis=shifted; degree+=1
    return arithmetic.polytrim(out[::-1].copy())


def _exact_small_s1(f,rows,specs,degree_bounds,coefficient_bounds):
    determinants=[]
    for columns in specs:
        matrix=sp.Matrix([[row[j].as_expr() for j in columns] for row in rows])
        determinants.append(sp.expand(matrix.det(method="domain-ge")))
    a,b = [sp.Poly(value,t,domain=sp.ZZ) for value in determinants]
    bound=max(coefficient_bounds)
    modulus=int(sp.nextprime(2*(bound+max(max(map(abs,p.all_coeffs())) for p in (a,b)))+1))
    s1=sp.Poly(a.as_expr()*Lambda+b.as_expr(),Lambda,t,domain=sp.ZZ)
    return SubresultantCertificate(s1,a,b,a,degree_bounds,
        (coefficient_bounds[0],coefficient_bounds[1],coefficient_bounds[0]),
        (modulus,),modulus,(modulus,),modulus,bound+max(max(map(abs,p.all_coeffs())) for p in (a,b)))


def first_subresultant_linear_exact(f,workers=1):
    """Reconstruct raw S1 from fixed minors, then prove it with a disjoint stream.

    Candidate FLINT determinants and verification NumPy determinants share only the
    stated Sylvester convention. No t-content is divided out.
    """
    _validate_workers(workers)
    rows,specs,degree_bounds,coefficient_bounds = _fixed_minor_data(f)
    if f.degree(Lambda) <= 6:
        return _exact_small_s1(f,rows,specs,degree_bounds,coefficient_bounds)
    prime_source=arithmetic.prime_stream()
    executor=ProcessPoolExecutor(max_workers=workers) if workers>1 else None
    def batch(bounds,method):
        primes=[next(prime_source) for _ in range(max(1,workers*2))]
        jobs=((f,p,bounds,method) for p in primes)
        return (list(executor.map(_minor_prime_job,jobs)) if executor is not None
                else [_minor_prime_job(job) for job in jobs])
    candidate_primes=[]
    residues=[None,None]
    modulus=1
    previous=None
    stable=0
    # Candidate reconstruction is proposal-only; the later disjoint stream owns proof.
    try:
        while True:
            done=False
            for p,values in batch(degree_bounds,"flint"):
                candidate_primes.append(p)
                if residues[0] is None:
                    residues=[list(value) for value in values]
                else:
                    residues=[[arithmetic.crt_pair(old,modulus,new,p) for old,new in zip(a,b)]
                              for a,b in zip(residues,values)]
                modulus*=p
                centered=[[x if x<=modulus//2 else x-modulus for x in value] for value in residues]
                candidate=tuple(sp.Poly.from_list(value,t,domain=sp.ZZ) for value in centered)
                if candidate == previous:
                    stable+=1
                else:
                    previous,stable=candidate,0
                if stable >= 2 and modulus > 2*max(max(map(abs,q.all_coeffs())) for q in candidate):
                    done=True; break
            if done: break
        a,b=candidate
        candidate_bound=max(max(map(abs,p.all_coeffs())) for p in (a,b))
        verification_degree_bounds=tuple(_verification_degree_bound(rows,columns) for columns in specs)
        verification_bound=max(coefficient_bounds)+candidate_bound
        verification_primes=[]
        verification_modulus=1
        # These primes are disjoint because candidate primes were consumed from the stream.
        while verification_modulus <= 2*verification_bound:
            for p,measured in batch(verification_degree_bounds,"numpy"):
                for actual,poly,bound in zip(measured,(a,b),verification_degree_bounds):
                    if not _minor_matches_degree_bound(actual,poly,p,bound):
                        raise ValueError(f"fixed-minor verification fails modulo {p}")
                verification_primes.append(p)
                verification_modulus*=p
                if verification_modulus > 2*verification_bound: break
    finally:
        if executor is not None: executor.shutdown()
    if verification_modulus <= 2*verification_bound:
        raise ValueError("verification modulus does not exceed coefficient-difference bound")
    s1=sp.Poly(a.as_expr()*Lambda+b.as_expr(),Lambda,t,domain=sp.ZZ)
    return SubresultantCertificate(s1,a,b,a,degree_bounds,
        (coefficient_bounds[0],coefficient_bounds[1],coefficient_bounds[0]),
        tuple(candidate_primes),modulus,tuple(verification_primes),verification_modulus,
        verification_bound)


def prove_pair_uniqueness(a2,certificate):
    if (sp.gcd(a2,certificate.psc1_raw).degree()!=0
            or sp.gcd(a2,certificate.a_raw).degree()!=0):
        raise ValueError("pair uniqueness fails: triple root or simultaneous double pairs")
    return True


def verify_subresultant_certificate(f,certificate,workers=1):
    """Re-execute fixed-minor evidence, rejecting corrupt bound metadata."""
    _validate_workers(workers)
    rows,specs,degrees,bounds=_fixed_minor_data(f)
    a,b=certificate.a_raw,certificate.b_raw
    if (certificate.degree_bounds!=degrees
            or certificate.coefficient_bounds!=(bounds[0],bounds[1],bounds[0])
            or certificate.psc1_raw!=a
            or certificate.s1!=sp.Poly(a.as_expr()*Lambda+b.as_expr(),Lambda,t,domain=sp.ZZ)):
        raise ValueError("raw fixed-minor coefficient or bounds mismatch")
    difference=max(bounds)+max(max(map(abs,q.all_coeffs())) for q in (a,b))
    if (certificate.verification_bound!=difference
            or certificate.verification_modulus<=2*difference
            or certificate.verification_modulus!=int(sp.prod(certificate.verification_primes))):
        raise ValueError("raw fixed-minor verification modulus/bound mismatch")
    if f.degree(Lambda)<=6:
        exact=_exact_small_s1(f,rows,specs,degrees,bounds)
        if certificate!=exact: raise ValueError("raw fixed-minor specialization mismatch")
        return True
    if set(certificate.candidate_primes)&set(certificate.verification_primes):
        raise ValueError("candidate and verifier prime streams overlap")
    if len(set(certificate.verification_primes))!=len(certificate.verification_primes):
        raise ValueError("duplicate verification prime")
    verification_degrees=tuple(_verification_degree_bound(rows,columns) for columns in specs)
    jobs=((f,p,verification_degrees,"numpy") for p in certificate.verification_primes)
    def check(values):
        for p,measured in values:
            for actual,poly,bound in zip(measured,(a,b),verification_degrees):
                if not _minor_matches_degree_bound(actual,poly,p,bound):
                    raise ValueError("raw fixed-minor specialization mismatch")
    if workers>1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            check(executor.map(_minor_prime_job,jobs))
    else: check(map(_minor_prime_job,jobs))
    return True


def prove_at_nonoverlap(a2,at,certificate):
    a,b=certificate.a_raw,certificate.b_raw
    # Homogeneous Horner evaluation modulo A2 avoids constructing the enormous
    # uncancelled expression a**degree * AT(-b/a) in ZZ(t).
    coefficients=sp.Poly(at.as_expr(),Lambda).all_coeffs()
    am,bm=a.rem(a2),(-b).rem(a2)
    remainder=sp.Poly(coefficients[0],t,domain=sp.ZZ).rem(a2)
    apower=sp.Poly(1,t,domain=sp.ZZ)
    for coefficient in coefficients[1:]:
        apower=(apower*am).rem(a2)
        remainder=(remainder*bm+sp.Poly(coefficient,t,domain=sp.ZZ)*apower).rem(a2)
    if sp.gcd(a2,remainder).degree()!=0:
        raise ValueError("AT overlap at an A2 repeated-root seed")
    return True


def _fraction(value):
    q=sp.Rational(value)
    return Fraction(int(q.p),int(q.q))


def _box_disjoint(a,b):
    return (a.real[1] < b.real[0] or b.real[1] < a.real[0]
            or a.imag[1] < b.imag[0] or b.imag[1] < a.imag[0])


def _box_neg(a):
    return ExactRootBox((-a.real[1],-a.real[0]),(-a.imag[1],-a.imag[0]),a.root_count)


def _box_conjugate(a):
    return ExactRootBox(a.real,(-a.imag[1],-a.imag[0]),a.root_count)


def _box_rotate_t_to_q(a):
    return ExactRootBox(a.imag,(-a.real[1],-a.real[0]))


def _box_half(a):
    return ExactRootBox((a.real[0]/2,a.real[1]/2),(a.imag[0]/2,a.imag[1]/2))


def isolate_a2_roots(poly,eps_den=10**24,*,capability=None):
    """Hybrid local certification in direct t; never global complex isolation."""
    proposals=propose_a2_boxes_flint(poly,128)
    return certify_a2_boxes(poly,proposals,max(16,int(eps_den).bit_length()),capability=capability)




def fujiwara_integer_radius(poly):
    """Ceil of 2*max(|a[n-i]/an|**(1/i), |a0/(2*an)|**(1/n)).

    Raising to i costs 2**i, not 2; the constant term costs 2**(n-1).
    All comparisons and the independent verification use integers only.
    """
    if (not isinstance(poly,sp.Poly) or poly.gens!=(t,) or poly.domain!=sp.ZZ
            or poly.is_zero or poly.degree()<1 or poly.LC()==0):
        raise ValueError("Fujiwara bound requires a nonconstant ZZ[t] polynomial")
    coefficients=[int(c) for c in poly.all_coeffs()]
    leading=abs(coefficients[0])
    radius=1
    degree=poly.degree()
    while any(leading*radius**i < 2**(i-(i==degree))*abs(coefficient)
              for i,coefficient in enumerate(coefficients[1:],1)):
        radius+=1
    verify_fujiwara_radius(poly,radius)
    return radius


def verify_fujiwara_radius(poly,radius):
    if (not isinstance(poly,sp.Poly) or poly.gens!=(t,) or poly.domain!=sp.ZZ
            or poly.is_zero or poly.degree()<1):
        raise ValueError("Fujiwara verification requires nonconstant ZZ[t]")
    if type(radius) is not int or radius<1:
        raise ValueError("invalid Fujiwara radius")
    coefficients=[int(c) for c in poly.all_coeffs()]
    # Separately evaluate rational powers of R/2, including the a0/2 term.
    half=Fraction(radius,2)
    n=len(coefficients)-1
    for i in range(1,n+1):
        term=Fraction(abs(coefficients[i]),abs(coefficients[0]))
        if i==n: term/=2
        if half**i<term:
            raise ValueError("radius does not satisfy the exact Fujiwara inequalities")
    return True


def _bounded_phase_worker(kind,args,queue,capability):
    try:
        if kind=="layer": result=prove_layer_identity(*args)
        elif kind=="layer_candidates": result=_lift_layers(*args)
        elif kind=="s1": result=first_subresultant_linear_exact(*args)
        elif kind=="verify_source": result=verify_inventory_source(*args)
        elif kind=="isolate": result=isolate_a2_roots(*args,capability=capability)
        elif kind=="propose": result=propose_a2_boxes_flint(*args)
        elif kind=="sturm": result=build_exact_sturm_chain(*args)
        elif kind=="certify": result=certify_a2_boxes(*args,capability=capability)
        elif kind=="certify_odd": result=certify_odd_transport(*args,capability=capability)
        elif kind=="quotient":
            a,b,boxes=args
            result=[_lambda_box(a,b,box) for box in boxes]
        elif kind=="localized_quotient":
            localized,boxes=args
            result=[localized_lambda_box(localized,box) for box in boxes]
        elif kind=="quotient_ring": result=build_quotient_ring_certificate(*args)
        elif kind=="quotient_ring_eval":
            certificate,boxes=args
            result=[quotient_ring_lambda_box(certificate,box) for box in boxes]
        elif kind=="pair_at":
            a2,at,s1=args
            result=(prove_pair_uniqueness(a2,s1),prove_at_nonoverlap(a2,at,s1))
        elif kind=="write": result=write_inventory(*args,capability=capability)
        elif kind=="write_generic": result=_write_generic_inventory(*args,capability=capability)
        elif kind=="count":
            poly,boxes,workers=args
            _validate_workers(workers)
            if workers>1:
                with ProcessPoolExecutor(max_workers=workers) as executor:
                    result=list(executor.map(_count_job,((poly,box) for box in boxes)))
            else: result=[_count_job((poly,box)) for box in boxes]
        else: raise ValueError("unknown bounded exact phase")
        queue.put(("ok",result))
    except BaseException as exc:
        queue.put(("error",f"{type(exc).__name__}: {exc}"))


def _bounded_exact_phase(kind,args,timeout_seconds=None,memory_bytes=None,*,capability=None):
    """Run one exact phase in a monitored process tree; fail on its cost boundary."""
    import psutil
    if timeout_seconds is None:
        timeout_seconds=_positive_integer(os.environ.get("ROUTE_B_N6_PHASE_TIMEOUT_SECONDS","600"))
    if memory_bytes is None:
        memory_bytes=_positive_integer(os.environ.get("ROUTE_B_N6_PHASE_MEMORY_GIB","16"))*1024**3
    if (type(timeout_seconds) not in (int,float) or not 0<timeout_seconds<=600
            or type(memory_bytes) is not int or not 0<memory_bytes<=16*1024**3):
        raise ValueError("phase caps must be positive and at most 600s/16GiB")
    queue=mp.Queue()
    needs_capability=kind in {"certify","certify_odd","isolate","write","write_generic"}
    if needs_capability and (type(capability) is not bytes or len(capability)!=32):
        raise ValueError("explicit run-local count transport capability required")
    if not needs_capability and capability is not None:
        raise ValueError("count capability must not enter unrelated phases")
    process=mp.Process(target=_bounded_phase_worker,args=(kind,args,queue,capability),name="n6-"+kind)
    process.start(); started=time.perf_counter(); peak=0
    identity=psutil.Process(process.pid)
    created,command=identity.create_time(),identity.cmdline()
    print(f"phase {kind}: pid={process.pid} command={command} timeout_seconds={timeout_seconds} memory_bytes={memory_bytes}",flush=True)
    try:
        response=None
        while response is None:
            if time.perf_counter()-started>timeout_seconds:
                raise RuntimeError(f"{kind} exact phase exceeded {timeout_seconds}s")
            try:
                tree=[psutil.Process(process.pid)]+psutil.Process(process.pid).children(recursive=True)
                rss=sum(p.memory_info().rss for p in tree if p.is_running())
                peak=max(peak,rss)
            except psutil.Error:
                rss=0
            if rss>memory_bytes:
                raise MemoryError(f"{kind} exact phase exceeded {memory_bytes} bytes RSS")
            try: response=queue.get(timeout=0.05)
            except queue_module.Empty:
                if not process.is_alive():
                    response=queue.get(timeout=3)
        process.join(3)
        status,value=response
        print(f"phase {kind}: {time.perf_counter()-started:.3f}s peak_tree_rss={peak}",flush=True)
        if status!="ok": raise RuntimeError(value)
        return value
    except BaseException:
        print(f"phase {kind}: stopped elapsed={time.perf_counter()-started:.3f}s peak_tree_rss={peak}",flush=True)
        try:
            root=psutil.Process(process.pid)
            if root.create_time()!=created or root.cmdline()!=command:
                raise RuntimeError("child PID identity changed; refusing termination")
            children=root.children(recursive=True)
        except psutil.Error:
            children=[]
        for child in reversed(children):
            try: child.terminate()
            except psutil.Error: pass
        if process.is_alive(): process.terminate()
        process.join()
        raise
    finally:
        queue.close()


def _iv_add(x,y): return (x[0]+y[0],x[1]+y[1])
def _iv_neg(x): return (-x[1],-x[0])
def _iv_mul(x,y):
    values=(x[0]*y[0],x[0]*y[1],x[1]*y[0],x[1]*y[1])
    return min(values),max(values)
def _iv_square(x):
    if x[0]<=0<=x[1]: return Fraction(0),max(x[0]*x[0],x[1]*x[1])
    return min(x[0]*x[0],x[1]*x[1]),max(x[0]*x[0],x[1]*x[1])


def _complex_mul(x,y):
    return (_iv_add(_iv_mul(x.real,y.real),_iv_neg(_iv_mul(x.imag,y.imag))),
            _iv_add(_iv_mul(x.real,y.imag),_iv_mul(x.imag,y.real)))


def _poly_box(poly,box):
    value=ExactRootBox((Fraction(0),Fraction(0)),(Fraction(0),Fraction(0)))
    bits=max(256,max(Fraction(v).denominator.bit_length() for v in box.real+box.imag)+64)
    scale=2**bits
    def outward(bounds):
        return (Fraction((bounds[0]*scale).__floor__(),scale),
                Fraction((bounds[1]*scale).__ceil__(),scale))
    def multiply(left,right):
        real,imag=_complex_mul(left,right)
        return ExactRootBox(outward(real),outward(imag))
    coefficients=list(poly.all_coeffs())
    valuation=0
    while len(coefficients)>1 and coefficients[-1]==0:
        coefficients.pop()
        valuation+=1
    for coefficient in coefficients:
        real,imag=_complex_mul(value,box)
        value=ExactRootBox(outward(_iv_add(real,(Fraction(int(coefficient)),)*2)),outward(imag))
    # Evaluate the original raw t-content by repeated squaring and multiply it
    # back in. No factor is cancelled in the polynomial, quotient or gcd gates.
    power=ExactRootBox((Fraction(1),Fraction(1)),(Fraction(0),Fraction(0)))
    base=box
    while valuation:
        if valuation&1: power=multiply(power,base)
        valuation//=2
        if valuation: base=multiply(base,base)
    value=multiply(value,power)
    return value


def _lambda_box(a,b,box):
    av=_poly_box(a,box)
    norm=_iv_add(_iv_square(av.real),_iv_square(av.imag))
    if norm[0]<=0:
        raise ValueError("quotient denominator interval contains zero")
    bv=_poly_box(b,box)
    inverse=(Fraction(1,norm[1]),Fraction(1,norm[0]))
    reciprocal=ExactRootBox(_iv_mul(av.real,inverse),_iv_mul(_iv_neg(av.imag),inverse))
    real,imag=_complex_mul(bv,reciprocal)
    # The interval operations above are exact. Outward dyadic rounding keeps the
    # exported enclosure compact without discarding any member of that interval.
    scale=2**96
    def outward(bounds):
        return (Fraction((bounds[0]*scale).__floor__(),scale),
                Fraction((bounds[1]*scale).__ceil__(),scale))
    return ExactRootBox(outward(_iv_neg(real)),outward(_iv_neg(imag)))


def _rat_json(value):
    return {"numerator":str(value.numerator),"denominator":str(value.denominator)}


def _box_json(box):
    return {"real":{"lower":_rat_json(box.real[0]),"upper":_rat_json(box.real[1])},
            "imag":{"lower":_rat_json(box.imag[0]),"upper":_rat_json(box.imag[1])}}


def _decimal_in_interval(bounds,target=None,digits=80):
    if bounds[0]==bounds[1]:
        value=bounds[0]
        denominator=value.denominator
        twos=fives=0
        while denominator%2==0: denominator//=2; twos+=1
        while denominator%5==0: denominator//=5; fives+=1
        if denominator!=1:
            raise ValueError("nonterminating point interval cannot carry a decimal seed")
        places=max(twos,fives)
        scaled=value.numerator*2**(places-twos)*5**(places-fives)
        if places==0: return str(scaled),value
        sign="-" if scaled<0 else ""
        raw=str(abs(scaled)).rjust(places+1,"0")
        return sign+raw[:-places]+"."+raw[-places:],value
    target=(bounds[0]+bounds[1])/2 if target is None else target
    scale=10**digits
    base=(target.numerator*scale)//target.denominator
    for k in (base,base+1,base-1):
        value=Fraction(k,scale)
        if bounds[0]<=value<=bounds[1]:
            sign="-" if k<0 else ""
            raw=str(abs(k)).rjust(digits+1,"0")
            text=sign+raw[:-digits]+"."+raw[-digits:]
            return text,value
    raise ValueError("could not place decimal seed inside exact interval")


def _seed_text(real,imag): return {"real":real,"imag":imag}


def _find_exact_partner(boxes,target):
    matches=[i for i,box in enumerate(boxes) if not _box_disjoint(target,box)]
    if len(matches)!=1:
        raise ValueError("exact transport does not identify one partner")
    return matches[0]


def _projective_seed_target(a,b,real,imag):
    """Exact Gaussian-rational scalar target, separate from the H/D enclosure."""
    z=sp.QQ_I.convert(sp.Rational(real)+sp.I*sp.Rational(imag))
    av=a.set_domain(sp.QQ_I).rep.eval(z)
    if not av:
        raise ValueError("projective scalar seed has a zero denominator")
    value=-b.set_domain(sp.QQ_I).rep.eval(z)/av
    return _fraction(value.x),_fraction(value.y)


def _verify_generic_inventory_source(source_pencil,layer,s1,workers=1):
    """Bind both algebraic certificates to the supplied authoritative E pencil.

    Box-carrier seals check transport integrity only. This independently reruns
    the integer discriminant proof and raw fixed-minor proof; it counts no roots.
    """
    if not isinstance(source_pencil,ParityPencil) or source_pencil.parity!="E":
        raise ValueError("writer requires the E source pencil")
    if not isinstance(layer,LayerCertificate):
        raise ValueError("writer requires an exact layer certificate")
    rebuilt=certify_layer_identity(source_pencil.residual,layer.valuation,layer.a1,
                                  layer.a2,layer.constant,workers)
    if rebuilt!=layer:
        raise ValueError("layer certificate metadata differs from source-pencil proof")
    verify_subresultant_certificate(source_pencil.residual,s1,workers)
    prove_pair_uniqueness(layer.a2,s1)
    prove_at_nonoverlap(layer.a2,source_pencil.at,s1)
    return True


def verify_inventory_source(source_pencil,layer,s1,workers=1):
    """N6 source pin plus independently executed algebraic certificate proofs."""
    _require_canonical_n6_source(source_pencil)
    return _verify_generic_inventory_source(source_pencil,layer,s1,workers)


def _make_inventory(layer,s1,a2,roots,odd_certificate,quotient_certificate=None,source_pencil=None,workers=1,*,capability=None,_n6=False):
    verify_box_certificate(a2,roots,EBoxCertificate,capability=capability)
    odd=sp.Poly(a2.as_expr().subs(t,-t),t,domain=sp.ZZ)
    verify_box_certificate(odd,odd_certificate,OBoxCertificate,roots,capability=capability)
    if len(roots)!=a2.degree() or a2!=layer.a2:
        raise ValueError("root inventory is incomplete")
    verifier=verify_inventory_source if _n6 else _verify_generic_inventory_source
    verifier(source_pencil,layer,s1,workers)
    source=(localize_monomial_quotient(a2,s1.a_raw,s1.b_raw,expected_valuation=496)
            if a2.degree()==133 else s1)
    ring=(build_quotient_ring_certificate(a2,source)
          if quotient_certificate is None else quotient_certificate)
    verify_quotient_ring_certificate(a2,source,ring)
    print("writer: completed E/O carriers and H/D identity checked; no count APIs",flush=True)
    seed_a,seed_b=_quotient_source_polynomials(source)
    if any(box.root_count!=1 for box in roots):
        raise ValueError("root inventory lacks an executed exact count of one")
    if any(not _box_disjoint(a,b) for j,a in enumerate(roots) for b in roots[j+1:]):
        raise ValueError("root inventory closures overlap")
    e=roots.boxes; o=odd_certificate.boxes
    prefix="N6" if _n6 else "GENERIC"
    ids={parity:[f"{prefix}-{parity}-A2-T-{i:03d}" for i in range(len(boxes))]
         for parity,boxes in (("E",e),("O",o))}
    lambda_e=[]
    for box in e:
        lambda_e.append(quotient_ring_lambda_box(ring,box))
    print(f"writer: {len(lambda_e)} H/D enclosures complete; building seeds",flush=True)
    loci=[]
    for parity,boxes in (("E",e),("O",o)):
        for index,box in enumerate(boxes):
            e_index=index if parity=="E" else _find_exact_partner(e,_box_neg(box))
            lambda_box=lambda_e[e_index]
            conj_index=_find_exact_partner(boxes,_box_conjugate(box))
            other="O" if parity=="E" else "E"
            other_boxes=o if parity=="E" else e
            parity_index=_find_exact_partner(other_boxes,_box_neg(box))
            tr,t_real=_decimal_in_interval(box.real)
            ti,t_imag=_decimal_in_interval(box.imag)
            # Evaluate the exact rational projective seed at the decimal t seed.
            sign=-1 if parity=="O" else 1
            lr_target,li_target=_projective_seed_target(seed_a,seed_b,sign*t_real,sign*t_imag)
            lr,lr_value=_decimal_in_interval(lambda_box.real,lr_target)
            li,li_value=_decimal_in_interval(lambda_box.imag,li_target)
            # Exact transforms are performed on the decimal rationals, before serialization.
            qr,qi=ti,("-"+tr if not tr.startswith("-") and t_real else tr[1:] if tr.startswith("-") else tr)
            if t_real==0: qi="0"
            pr,pr_value=_decimal_in_interval((lr_value/2,lr_value/2))
            pi,pi_value=_decimal_in_interval((li_value/2,li_value/2))
            loci.append({"id":ids[parity][index],"parity":parity,"algebraicMultiplicity":2,
                "tBox":_box_json(box),"qPhysicalCSharpBox":_box_json(_box_rotate_t_to_q(box)),
                "lambdaClearedBox":_box_json(lambda_box),
                "lambdaPhysicalBox":_box_json(_box_half(lambda_box)),
                "tSeed":_seed_text(tr,ti),"qPhysicalCSharpSeed":_seed_text(qr,qi),
                "lambdaClearedSeed":_seed_text(lr,li),
                "lambdaPhysicalSeed":_seed_text(pr,pi),
                "conjugationPartnerId":ids[parity][conj_index],
                "parityPartnerId":ids[other][parity_index]})
            if len(loci)%32==0:
                print(f"writer: seeds={len(loci)}/{2*len(e)}",flush=True)
    payload={"schemaVersion":3 if _n6 else "generic-polynomial-inventory-v1","n":6 if _n6 else None,
        "model":"open-uniform-XY-delta0-field0-gamma1-SEket-DEbra" if _n6 else "generic-polynomial",
        "conventions":{"parameter":"t=i*qCSharp;qCSharp=-i*t",
            "eigenvalue":"Lambda=2*lambda",
            "coefficientOrder":"lambda-lowest-first,t-lowest-first"},
        "a2Degrees":{"E":a2.degree(),"O":a2.degree()},
        "layerIdentity":{"discriminantDegree":layer.deg_d,"valuation":layer.valuation,
            "a1Degree":layer.a1.degree(),"a2Degree":layer.a2.degree(),
            "constant":str(layer.constant),"proofModulus":str(layer.proof_modulus),
            "proofBound":str(layer.proof_bound)},"loci":loci,"exactRankCertificates":[]}
    if _n6: payload["sourcePencilDigest"]=source_pencil.source_digest
    return payload


def write_inventory(path,layer,s1,a2,roots,odd_certificate=None,payload=None,workers=1,quotient_certificate=None,source_pencil=None,*,capability=None):
    """Public canonical N6 schema-3 writer. Generic pencils are never N6-labelled."""
    _require_canonical_n6_source(source_pencil)
    return _write_generic_inventory(path,layer,s1,a2,roots,odd_certificate,payload,workers,
                                    quotient_certificate,source_pencil,capability=capability,_n6=True)


def _write_generic_inventory(path,layer,s1,a2,roots,odd_certificate=None,payload=None,workers=1,quotient_certificate=None,source_pencil=None,*,capability=None,_n6=False):
    """Internal polynomial writer; generic output has no N6 model claim."""
    canonical=_make_inventory(layer,s1,a2,roots,odd_certificate,quotient_certificate,source_pencil,workers,
                              capability=capability,_n6=_n6)
    if payload is not None and payload!=canonical:
        raise ValueError("supplied inventory differs from exact canonical reconstruction")
    target=Path(path)
    text=json.dumps(canonical,indent=2,ensure_ascii=False)+"\n"
    descriptor,temporary=tempfile.mkstemp(prefix=target.name+".",suffix=".tmp",dir=target.parent)
    os.close(descriptor)
    staging=Path(temporary)
    try:
        staging.write_text(text,encoding="utf-8",newline="\n")
        os.replace(staging,target)
    finally:
        staging.unlink(missing_ok=True)
    return canonical


def verify_inventory_artifact_structure(path):
    """Cheap persisted-schema check, NOT a root-count or source-pencil proof.

    Mathematical source binding belongs to verify_inventory_source during the
    writer phase; no executed-count provenance is recoverable from JSON alone.
    """
    def require(condition):
        if not condition: raise ValueError("invalid persisted N6 inventory structure")
    def integer(text):
        require(type(text) is str)
        value=int(text)
        require(str(value)==text)
        return value
    def rational(value):
        require(type(value) is dict and set(value)=={"numerator","denominator"})
        n,d=integer(value["numerator"]),integer(value["denominator"])
        require(d>0)
        q=Fraction(n,d)
        require((q.numerator,q.denominator)==(n,d))
        return q
    def box(value):
        require(type(value) is dict and set(value)=={"real","imag"})
        axes=[]
        for axis in ("real","imag"):
            require(set(value[axis])=={"lower","upper"})
            lo,hi=(rational(value[axis][edge]) for edge in ("lower","upper"))
            require(lo<=hi); axes.append((lo,hi))
        return ExactRootBox(*axes)
    def seed(value,b):
        require(type(value) is dict and set(value)=={"real","imag"})
        for axis in ("real","imag"):
            require(type(value[axis]) is str and "/" not in value[axis])
        z=Fraction(value["real"]),Fraction(value["imag"])
        require(all(lo<=v<=hi for v,(lo,hi) in zip(z,(b.real,b.imag))))
        return z
    try:
        raw=Path(path).read_bytes()
        data=json.loads(raw.decode("utf-8"))
        require(raw==(json.dumps(data,indent=2,ensure_ascii=False)+"\n").encode("utf-8"))
        require(set(data)=={"schemaVersion","n","model","conventions","a2Degrees","layerIdentity","loci","exactRankCertificates","sourcePencilDigest"})
        require(type(data["schemaVersion"]) is int and data["schemaVersion"]==3)
        require(data["sourcePencilDigest"]==N6_SOURCE_PENCIL_DIGEST)
        fixture=Path(__file__).parent/"tests/fixtures/route_b_a2_n6_residual.json"
        require(load_exact_pencils(fixture)["E"].source_digest==data["sourcePencilDigest"])
        require(type(data["n"]) is int and data["n"]==6)
        require(data["model"]=="open-uniform-XY-delta0-field0-gamma1-SEket-DEbra")
        require(data["conventions"]=={"parameter":"t=i*qCSharp;qCSharp=-i*t","eigenvalue":"Lambda=2*lambda",
                                      "coefficientOrder":"lambda-lowest-first,t-lowest-first"})
        require(data["a2Degrees"]=={"E":133,"O":133} and data["exactRankCertificates"]==[])
        proof=data["layerIdentity"]
        require(set(proof)=={"discriminantDegree","valuation","a1Degree","a2Degree","constant","proofModulus","proofBound"})
        for key,expected in (("discriminantDegree",926),("valuation",536),("a1Degree",124),("a2Degree",133)):
            require(type(proof[key]) is int and proof[key]==expected)
        require(integer(proof["constant"])!=0)
        bound=integer(proof["proofBound"])
        require(bound>0 and integer(proof["proofModulus"])>2*bound)
        rows=data["loci"]
        require(type(rows) is list and len(rows)==266)
        expected_ids=[f"N6-{parity}-A2-T-{i:03d}" for parity in ("E","O") for i in range(133)]
        require([row["id"] for row in rows]==expected_ids)
        parsed={}
        for index,row in enumerate(rows):
            require(set(row)=={"id","parity","algebraicMultiplicity","tBox","qPhysicalCSharpBox","lambdaClearedBox",
                "lambdaPhysicalBox","tSeed","qPhysicalCSharpSeed","lambdaClearedSeed","lambdaPhysicalSeed",
                "conjugationPartnerId","parityPartnerId"})
            require(row["parity"]==("E" if index<133 else "O"))
            require(type(row["algebraicMultiplicity"]) is int and row["algebraicMultiplicity"]==2)
            tb,qb,lb,pb=(box(row[key]) for key in ("tBox","qPhysicalCSharpBox","lambdaClearedBox","lambdaPhysicalBox"))
            require(qb==_box_rotate_t_to_q(tb) and pb==_box_half(lb))
            require(not (tb.real[0]<=0<=tb.real[1] and tb.imag[0]<=0<=tb.imag[1]))
            require(tb.imag==(0,0) or tb.imag[1]<0 or tb.imag[0]>0)
            ts,qs,ls,ps=(seed(row[key],b) for key,b in zip(
                ("tSeed","qPhysicalCSharpSeed","lambdaClearedSeed","lambdaPhysicalSeed"),(tb,qb,lb,pb)))
            require(qs==(ts[1],-ts[0]) and ps==(ls[0]/2,ls[1]/2))
            parsed[row["id"]]=(tb,lb)
        by_id={row["id"]:row for row in rows}
        for parity in ("E","O"):
            boxes=[parsed[row["id"]][0] for row in rows if row["parity"]==parity]
            require(boxes==sorted(boxes,key=lambda b:(b.real,b.imag)))
            require(sum(b.imag==(0,0) for b in boxes)==59 and sum(b.imag[0]>0 for b in boxes)==37)
            require(all(_box_disjoint(a,b) for i,a in enumerate(boxes) for b in boxes[i+1:]))
        for row in rows:
            tb,lb=parsed[row["id"]]
            cp,pp=by_id[row["conjugationPartnerId"]],by_id[row["parityPartnerId"]]
            require(cp["conjugationPartnerId"]==row["id"] and cp["parity"]==row["parity"])
            require(pp["parityPartnerId"]==row["id"] and pp["parity"]!=row["parity"])
            require(parsed[cp["id"]]==(_box_conjugate(tb),_box_conjugate(lb)))
            require(parsed[pp["id"]]==(_box_neg(tb),lb))
        return data
    except (KeyError,TypeError,ZeroDivisionError,UnicodeError) as exc:
        raise ValueError("invalid persisted N6 inventory structure") from exc


class _ProcessTreeSampler:
    """Sample simultaneous process-tree RSS and OS-reported per-process peaks.

    psutil is a benchmark-only dependency. Polling can miss short-lived children;
    the JSON states the sampling interval instead of claiming an exact tree peak.
    """
    def __init__(self, pid):
        import psutil
        self.psutil, self.pid = psutil, pid
        self.peak_tree = 0
        self.peaks = {}
        self.stop_event = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)

    def _sample(self):
        try:
            root = self.psutil.Process(self.pid)
            processes = [root]+root.children(recursive=True)
        except self.psutil.Error:
            return
        current = 0
        for process in processes:
            try:
                info = process.memory_info()
                current += info.rss
                key = f"{process.pid}:{process.create_time()}"
                peak = getattr(info,"peak_wset",info.rss)
                self.peaks[key] = max(self.peaks.get(key,0),peak)
            except self.psutil.Error:
                continue
        self.peak_tree = max(self.peak_tree,current)

    def _run(self):
        while not self.stop_event.is_set():
            self._sample()
            self.stop_event.wait(0.02)

    def __enter__(self):
        self._sample()
        self.thread.start()
        return self

    def __exit__(self, *args):
        self._sample()
        self.stop_event.set()
        self.thread.join()

    def report(self):
        return {"peak_tree_rss_bytes":self.peak_tree,
                "observed_process_peak_wset_bytes":self.peaks,
                "memory_sampling_seconds":0.02}


def recommend_workers(samples, available_memory):
    """Choose measured projected speed among runs using <= half available RAM."""
    eligible = [s for s in samples if s["peak_tree_rss_bytes"] <= available_memory//2]
    if not eligible:
        raise ValueError("no measured worker configuration fits the memory budget")
    return min(eligible,key=lambda s:(s.get("projected_proof_seconds",s["seconds_per_prime"]),
                                      s["workers"]))["workers"]


def _benchmark_probe(kind, poly, queue):
    try:
        if kind == "s1":
            sequence = sp.subresultants(poly.as_expr(),poly.diff(Lambda).as_expr(),Lambda)
            result = {"linear_entries":sum(sp.degree(s,Lambda)==1 for s in sequence)}
        elif kind == "real":
            proposals=propose_a2_boxes_flint(poly,128)
            boxes=[outward_dyadic_box(ball,96) for ball in proposals if _ball_bounds(ball)[1]==(0,0)]
            if not boxes:
                raise ValueError("N6 A2 candidate has no real root to probe")
            count=exact_count_in_real_interval(poly,*boxes[0].real)
            if count!=1: raise ValueError("real local exact count is not one")
            result = {"local_interval":str(boxes[0].real),"exact_root_count":count,
                      "dyadic_depth":96,"global_isolation":False}
        else:
            # Floating roots propose boxes only; exact argument-principle counts
            # decide whether a box isolates a root of this integer polynomial.
            coefficients = np.array([float(c/poly.LC()) for c in poly.all_coeffs()])
            guesses = sorted((z for z in np.roots(coefficients) if z.imag > 0.001),key=abs)
            result = None
            for z in guesses[:8]:
                re = sp.Rational(str(round(float(z.real),10)))
                im = sp.Rational(str(round(float(z.imag),10)))
                radius = sp.Rational(1,10**4)
                lower,upper = re-radius+sp.I*(im-radius),re+radius+sp.I*(im+radius)
                if poly.count_roots(lower,upper) != 1:
                    continue
                tight = sp.Rational(1,10**8)
                lo,hi = re-tight+sp.I*(im-tight),re+tight+sp.I*(im+tight)
                if poly.count_roots(lo,hi) == 1:
                    result = {"isolation_box":[str(lower),str(upper)],
                              "refined_box":[str(lo),str(hi)],"exact_root_count":1}
                    break
            if result is None:
                raise ValueError("numeric proposals did not yield an exact isolated refined root")
        queue.put({"status":"complete",**result})
    except Exception as exc:
        queue.put({"status":"error","error":str(exc)})


def _timed_probe(kind, poly, budget=30):
    queue = mp.Queue()
    process = mp.Process(target=_benchmark_probe,args=(kind,poly,queue))
    started = time.perf_counter()
    process.start()
    with _ProcessTreeSampler(process.pid) as memory:
        process.join(budget)
        if process.is_alive():
            process.terminate()
            process.join()
            report = {"status":"timeout","budget_seconds":budget}
        else:
            try:
                report = queue.get(timeout=2)
            except Exception as exc:
                report = {"status":"error","error":str(exc),"exitcode":process.exitcode}
    queue.close()
    return {**report,"seconds":time.perf_counter()-started,**memory.report()}


def benchmark(path):
    """Measure actual export, N6 arithmetic, child memory and unresolved S1 work.

    The fresh A2 lift is only a timing candidate. This benchmark never promotes
    it to an exact N6 layer; prove_layer_identity is the separate certificate.
    """
    import psutil
    root = Path(__file__).resolve().parents[1]
    output = root/"simulations/results/route_b_a2_n6_benchmark.txt"
    fresh = root/"simulations/results/route_b_a2_n6_benchmark_residual.tmp.json"
    report = {"schemaVersion":1}
    command = ["dotnet","run","--project","compute/RCPsiSquared.Cli","-c","Release",
               "--","route-b-n6-residual","--out",str(fresh)]
    started = time.perf_counter()
    exporter = subprocess.Popen(command,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    with _ProcessTreeSampler(exporter.pid) as memory:
        stdout,stderr = exporter.communicate()
    report["export"] = {"seconds":time.perf_counter()-started,"returncode":exporter.returncode,
                        "stdout":stdout,"stderr":stderr,**memory.report()}
    if exporter.returncode:
        output.write_text(json.dumps(report,indent=2),encoding="utf-8")
        raise RuntimeError(f"fresh C# export failed: {stderr}")
    print(f"benchmark export: {report['export']['seconds']:.3f}s",flush=True)
    started = time.perf_counter()
    pencils = load_exact_pencils(fresh)
    if pencils != load_exact_pencils(path):
        raise ValueError("fresh export differs from the supplied exact pencils")
    e = pencils["E"]
    report["load_seconds"] = time.perf_counter()-started
    started = time.perf_counter()
    with _ProcessTreeSampler(os.getpid()) as memory:
        stream = _disc_stream(e.residual,arithmetic.prime_stream(),1)
        p,disc = next(stream)
        stream.close()
        layers = arithmetic.layers_from_disc(disc,p,{"v":536,"A1":124,"A2":133})
    report["one_prime"] = {"seconds":time.perf_counter()-started,"prime":p,
                           "degree":len(disc)-1,"valuation":layers[2],**memory.report()}
    samples = []
    for workers in (1,2,4,8):
        if workers > (psutil.cpu_count(logical=False) or os.cpu_count() or 1):
            continue
        started = time.perf_counter()
        with _ProcessTreeSampler(os.getpid()) as memory:
            stream = _disc_stream(e.residual,arithmetic.prime_stream(),workers)
            next(stream)
            startup = time.perf_counter()-started
            steady = time.perf_counter()
            for _ in range(32):
                next(stream)
            per_prime = (time.perf_counter()-steady)/32
            stream.close()
        # The row-sum determinant bound estimates the actual proof pool size.
        count = ((_det_bound(e.residual)*2).bit_length()+24)//25
        samples.append({"workers":workers,"startup_seconds":startup,
                        "seconds_per_prime":per_prime,"sample_primes":32,
                        "projected_proof_seconds":startup+count*per_prime,**memory.report()})
        print(f"benchmark workers={workers}: steady {per_prime:.3f}s/prime",flush=True)
    available = psutil.virtual_memory().available
    workers = recommend_workers(samples,available)
    report.update(worker_samples=samples,recommended_workers=workers,available_memory_bytes=available)
    report["s1"] = _timed_probe("s1",e.residual)
    print(f"benchmark exact S1: {report['s1']['status']}",flush=True)
    unresolved = report["s1"]["status"] != "complete"
    report["task4_checkpoint"] = {
        "status":"REQUIRES_REVIEWED_ALTERNATIVE" if unresolved else "REQUIRES_S1_VALIDATION",
        "required_route":"CRT_WITH_RIGOROUS_COEFFICIENT_BOUND" if unresolved else "EXACT_S1_VALIDATION",
        "s1_proved":False,
        "reason":"A timing probe is not the Task4 uniqueness/nonvanishing certificate."}
    started = time.perf_counter()
    with _ProcessTreeSampler(os.getpid()) as memory:
        _,a2 = _lift_layers(e.residual,workers)
    report["roots"] = {"source":"fresh N6 E A2 CRT candidate; timing only",
                       "candidate_degree":a2.degree(),"lift_seconds":time.perf_counter()-started,
                       "lift_memory":memory.report(),"probes":{}}
    for kind in ("real","nonreal"):
        report["roots"]["probes"][kind] = _timed_probe(kind,a2)
        print(f"benchmark N6 A2 {kind}: {report['roots']['probes'][kind]['status']}",flush=True)
    output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2),flush=True)
    return report


def _positive_integer(value):
    if type(value) is not str or re.fullmatch(r"[1-9][0-9]*",value) is None:
        raise argparse.ArgumentTypeError("expected a strictly positive decimal integer")
    return int(value)


def parse_runtime_arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--residual-json",type=Path,required=True)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--benchmark",action="store_true")
    action.add_argument("--prove-layers",action="store_true")
    action.add_argument("--probe-isolation",action="store_true")
    action.add_argument("--out",type=Path)
    for flag in ("workers","phase-timeout-seconds","phase-memory-gib"):
        parser.add_argument("--"+flag,type=_positive_integer)
    args = parser.parse_args(argv)
    for field,default in (("workers","1"),("phase_timeout_seconds","600"),("phase_memory_gib","16")):
        if getattr(args,field) is None:
            try: value=_positive_integer(os.environ.get("ROUTE_B_N6_"+field.upper(),default))
            except argparse.ArgumentTypeError as exc: parser.error(str(exc))
            setattr(args,field,value)
    if args.phase_timeout_seconds>600 or args.phase_memory_gib>16:
        parser.error("phase caps cannot exceed 600 seconds / 16 GiB")
    try: _validate_workers(args.workers)
    except ValueError as exc: parser.error(str(exc))
    return args


def main():
    args=parse_runtime_arguments()
    # Private to this run; never stored on a carrier or in module state.
    capability=os.urandom(32)
    def phase(kind,arguments):
        credential=capability if kind in {"certify","certify_odd","isolate","write"} else None
        return _bounded_exact_phase(kind,arguments,args.phase_timeout_seconds,
                                    args.phase_memory_gib*1024**3,capability=credential)
    if args.out is not None or args.probe_isolation:
        _require_flint()
    if args.benchmark:
        benchmark(args.residual_json)
        return
    if args.probe_isolation:
        started=time.perf_counter()
        _,candidate_a2=phase("layer_candidates",(load_exact_pencils(args.residual_json)["E"].residual,args.workers))
        print("isolation probe source: fresh N6 CRT candidate; timing only, not layer proof",flush=True)
        roots=phase("isolate",(candidate_a2,10**18))
        print(f"isolation probe roots={len(roots)} radius={fujiwara_integer_radius(candidate_a2)} "
              f"total_seconds={time.perf_counter()-started:.3f}",flush=True)
        return
    started = time.perf_counter()
    cert = phase("layer",(load_exact_pencils(args.residual_json)["E"],args.workers))
    if args.out is not None:
        pencils=load_exact_pencils(args.residual_json)
        s1=phase("s1",(pencils["E"].residual,args.workers))
        phase("pair_at",(cert.a2,pencils["E"].at,s1))
        localized=localize_monomial_quotient(cert.a2,s1.a_raw,s1.b_raw,expected_valuation=496)
        ring=phase("quotient_ring",(cert.a2,localized))
        sturm_chain=phase("sturm",(cert.a2,))
        print(f"exact Sturm chain: length={len(sturm_chain.members)} "
              f"degrees={[len(row)-1 for row in sturm_chain.members]}",flush=True)
        print(f"quotient representation: raw_valuations={localized.valuation},{localized.valuation} "
              f"unit_degrees={localized.a_unit.degree()},{localized.b_unit.degree()} "
              f"degree_H={ring.H.degree()} denominator_bits={ring.D.bit_length()}",flush=True)
        roots=None
        for precision,bits in ((128,96),(256,128),(512,192)):
            print(f"hybrid precision={precision} dyadic_depth={bits} workers={args.workers}",flush=True)
            try:
                proposals=phase("propose",(cert.a2,precision))
                candidate=phase("certify",(cert.a2,proposals,bits,args.workers,sturm_chain))
                phase("quotient_ring_eval",(ring,candidate))
                roots=candidate; break
            except RuntimeError as exc:
                # Resource-limit failures stop the run. Only an unresolved
                # proposal/count or denominator exclusion permits refinement.
                if not str(exc).startswith("ValueError:"): raise
                print(f"unresolved: {exc}",flush=True)
        if roots is None:
            raise ValueError("could not certify all A2 root boxes and H/D enclosures")
        odd_certificate=phase("certify_odd",(cert.a2,roots,args.workers,sturm_chain))
        payload=phase("write",(args.out,cert,s1,cert.a2,roots,odd_certificate,None,args.workers,ring,pencils["E"]))
        print(f"EXACT: loci={len(payload['loci'])} S1deg=({s1.a_raw.degree()},{s1.b_raw.degree()}) "
              f"candidate_primes={len(s1.candidate_primes)} verification_primes={len(s1.verification_primes)}")
        print(f"S1_verification_modulus={s1.verification_modulus}")
        print(f"S1_verification_bound={s1.verification_bound}")
        print(f"runtime_seconds={time.perf_counter()-started:.3f}")
        return
    print(f"EXACT: degD={cert.deg_d} valuation={cert.valuation} A1={cert.a1.degree()} A2={cert.a2.degree()}")
    print(f"C={cert.constant}")
    print(f"proof_primes={len(cert.proof_primes)} proof_modulus={cert.proof_modulus} proof_bound={cert.proof_bound}")
    print(f"A1={cert.a1.as_expr()}\nA2={cert.a2.as_expr()}")
    print(f"runtime_seconds={time.perf_counter()-started:.3f}")


if __name__ == "__main__":
    main()
