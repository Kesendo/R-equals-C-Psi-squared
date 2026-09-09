"""Validated algebraic base for F163, conditional on its stored exact certificate.

The exact rank/semisimplicity premise comes from F163; ball residuals are
consistency checks, never proofs that a matrix expression vanishes exactly.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from scipy.linalg import qr
from flint import acb, arb, acb_poly, acb_mat, ctx
import route_b_a2_n6 as exact
from route_b_n6_exact_unfolding import integer_pencil

ROOT=Path(__file__).resolve().parents[1]

def select_near(roots, seed):
    return min(roots, key=lambda z: abs(complex(z)-seed))

def ball_matrix(rows):
    return acb_mat([[acb(x) for x in row] for row in rows])

def identity(n):
    return ball_matrix([[int(i==j) for j in range(n)] for i in range(n)])

def all_contain_zero(matrix):
    return all(matrix[i,j].contains(0) for i in range(matrix.nrows()) for j in range(matrix.ncols()))

def build_base(precision=512):
    ctx.prec=precision
    fixture_path=ROOT/'simulations/tests/fixtures/route_b_a2_n6_residual.json'
    producer=ROOT/'simulations/route_b_n6_exact_unfolding.py'
    report=json.loads((ROOT/'simulations/results/route_b_n6_exact_unfolding.json').read_text())
    assert report['provenance']['fixture_sha256']==hashlib.sha256(fixture_path.read_bytes()).hexdigest()
    assert report['provenance']['script_sha256']==hashlib.sha256(producer.read_bytes()).hexdigest()
    exact.load_exact_pencils(fixture_path) # verifies the canonical semantic digest
    assert report['provenance']['source_pencil_digest']==exact.N6_SOURCE_PENCIL_DIGEST
    assert int(report['layer']['proof_modulus'])>2*int(report['layer']['proof_bound'])
    fixture=json.loads(fixture_path.read_text())
    a=acb_poly([int(c) for c in report['layer']['a2_coefficients_lowest_first']])
    roots=a.roots(tol=arb('1e-100'), maxprec=4096)
    assert len(roots)==133
    t0=select_near(roots, complex(-1.0968139484731172,1.6725941457514577))
    assert abs(complex(t0)-complex(-1.0968139484731172,1.6725941457514577))<1e-10
    atlas=json.loads((ROOT/'simulations/results/route_b_a2_n6.json').read_text())
    locus=next(x for x in atlas['loci'] if x['id']=='N6-E-A2-T-007')
    for coordinate,z in (('real',t0.real),('imag',t0.imag)):
        endpoints=locus['tBox'][coordinate]
        lo=arb(endpoints['lower']['numerator'])/arb(endpoints['lower']['denominator'])
        hi=arb(endpoints['upper']['numerator'])/arb(endpoints['upper']['denominator'])
        assert z>lo and z<hi
    rows=fixture['rEven']['residualInT']
    f=acb_poly([acb_poly([int(c) for c in row])(t0) for row in rows])
    critical=f.derivative().roots(tol=arb('1e-70'), maxprec=4096)
    cleared=select_near(critical,complex(-8.55925566129688,-1.5204627789359494))
    assert abs(complex(cleared)-complex(-8.55925566129688,-1.5204627789359494))<1e-10
    assert not f.derivative().derivative()(cleared).contains(0)
    assert f(cleared).contains(0)
    # F163 supplies one repeated root. Excluding every other critical root
    # identifies the selected one without treating zero containment as equality.
    others=[z for z in critical if z is not cleared]
    assert len(others)==30 and all(not f(z).contains(0) for z in others)
    lam0=cleared/2
    d,bonds,reflection,_=integer_pencil()
    hops=[[sum(b[i][j] for b in bonds) for j in range(90)] for i in range(90)]
    reps=[i for i in range(90) if i<reflection[i]]
    ae=acb_mat([[d[i][j]+d[i][reflection[j]]+t0*(hops[i][j]+hops[i][reflection[j]])-int(ii==jj)*lam0
        for jj,j in enumerate(reps)] for ii,i in enumerate(reps)])
    midpoint=np.array([[complex(ae[i,j]) for j in range(45)] for i in range(45)])
    _,_,pivots=qr(midpoint,pivoting=True)
    chosen=list(map(int,pivots[:43])); free=list(map(int,pivots[43:]))
    minor=acb_mat([[ae[i,j] for j in chosen] for i in chosen])
    rhs=acb_mat([[ae[i,j] for j in free] for i in chosen])
    solved=-minor.solve(rhs)
    ue=acb_mat(45,2)
    for ii,i in enumerate(chosen):
        for j in range(2): ue[i,j]=solved[ii,j]
    for j,i in enumerate(free): ue[i,j]=1
    u=acb_mat(90,2)
    for k,i in enumerate(reps):
        for j in range(2): u[i,j]=u[reflection[i],j]=ue[k,j]
    w=(u.transpose()*u).inv()*u.transpose()
    p=u*w; ident=identity(90); q=ident-p
    l0=acb_mat([[d[i][j]+t0*hops[i][j] for j in range(90)] for i in range(90)])
    s=q*(ident*lam0-l0+p).inv()*q
    residuals={'eigenbasis': all_contain_zero((l0-ident*lam0)*u),
        'dual':all_contain_zero(w*u-identity(2)),
        'projector':all_contain_zero(p*p-p),
        'resolvent':all_contain_zero((ident*lam0-l0)*s-q)}
    assert all(residuals.values())
    q_bonds=[ball_matrix(b)*acb(0,1) for b in bonds]
    return dict(C=ball_matrix(hops)*acb(0,1),left=q_bonds[0],right=q_bonds[-1],
        base_matrix=ident*lam0-l0+p,t0=t0,q0=-acb(0,1)*t0,lambda0=lam0,U=u,W=w,P=p,Q=q,S=s,L0=l0,
        D=ball_matrix(d),T=ball_matrix(hops),bonds=[ball_matrix(b) for b in bonds],
        reflection=reflection,checks=residuals,precision=precision)

if __name__=='__main__':
    start=time.perf_counter(); base=build_base()
    for key in ('t0','q0','lambda0','checks'): print(key,base[key])
    print('seconds',time.perf_counter()-start)
