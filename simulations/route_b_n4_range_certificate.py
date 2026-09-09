"""Real-ball certificates of isolated N4 EP-curve singularities.

Certifies singularity existence/type, not absence of earlier obstructions
along a numerically tracked path and not a maximal Taylor disk.
"""
import hashlib,json
from pathlib import Path
import sympy as s
import mpmath as mp
import numpy as np
from scipy.linalg import qr
from flint import arb,acb,acb_mat,ctx
from route_b_n4_range_polynomial import factor_expression,even_frequency_factor
from route_b_other_n_unfolding import parts

ROOT=Path(__file__).resolve().parents[1]
TARGETS=[
 ('one','triple',('-1.72811602403820664625','.44008951151008331695','-.005138273760692603364')),
 ('one','triple',('-1.66617275439538828710','.426866181753943982808','-.013736634404318642028')),
 ('one','fold',('-9.351510723066957174276','.500823550468681932645','.2647726098223231646961')),
 ('odd','triple',('-0.9441597816697135180265','.24201564006955527180','.5804402904296053012113')),
 ('odd','fold',('-1.781118999662660009548','.446584835339167607503','.0181112843381541134791')),
 ('even','triple',('1.9603436027347338485162','.6947395054024396774154','.540020222277755101278')),
 ('even','fold',('2','2','-.5857864376269049511983')),
]


def profile(name):
    if name=='even':return even_frequency_factor()
    p,(x,y,a,b)=factor_expression();e=s.Symbol('e')
    av,bv=(1+e,1) if name=='one' else (1+e/2,1-e/2)
    return s.Poly(p.subs({a:av,b:bv}),x,y,e).as_expr(),(x,y,e)


def ball_function(expr,variables):
    terms=s.Poly(expr,*variables).terms()
    def evaluate(point):
        powers=[[v**j for j in range(max(p[i] for p,c in terms)+1)] for i,v in enumerate(point)]
        total=arb(0)
        for p,c in terms:
            v=arb(int(s.numer(c)))/arb(int(s.denom(c)))
            for i,k in enumerate(p):v*=powers[i][k]
            total+=v
        return total
    return evaluate


def certify(name,kind,guess):
    p,variables=profile(name);x,y,e=variables
    selected=s.diff(p,x,2) if kind=='triple' else s.diff(p,y)
    equations=s.Matrix([p,s.diff(p,x),selected])
    jacobian=equations.jacobian(variables)
    mp.mp.dps=90
    fn=s.lambdify(variables,equations,'mpmath',cse=True)
    jn=s.lambdify(variables,jacobian,'mpmath',cse=True)
    root=mp.findroot(lambda *v:tuple(fn(*v)),tuple(mp.mpf(v) for v in guess),J=jn,tol=mp.mpf('1e-75'),maxsteps=70)
    center=[arb(str(v)).mid() for v in root]
    funcs=[ball_function(v,variables) for v in equations]
    jfunc=[[ball_function(jacobian[i,j],variables) for j in range(3)] for i in range(3)]
    j0=acb_mat([[f(center) for f in row] for row in jfunc]).inv()
    Y=acb_mat([[acb(j0[i,j].real.mid()) for j in range(3)] for i in range(3)])
    assert not Y.det().contains(0)
    radius=arb('1e-25')
    def check(c):
        box=[v+arb(0,radius) for v in c]
        J=acb_mat([[f(box) for f in row] for row in jfunc])
        defect=acb_mat([[int(i==j) for j in range(3)] for i in range(3)])-Y*J
        F=Y*acb_mat([[f(c)] for f in funcs])
        rows=[sum((abs(defect[i,j]).upper() for j in range(3)),arb(0)) for i in range(3)]
        images=[abs(F[i,0]).upper()+rows[i]*radius for i in range(3)]
        return all(v<1 for v in rows) and all(v<radius for v in images),box,rows,images
    success,box,rows,images=check(center)
    assert success
    displaced=center.copy();displaced[2]+=arb('0.001')
    assert not check(displaced)[0]
    nonzero={}
    for label,expr in [('second',s.diff(p,x,2)),('third',s.diff(p,x,3)),('parameter',s.diff(p,y))]:
        value=ball_function(expr,variables)(box)
        nonzero[label]=str(value)
        if (kind=='triple' and label in ('third','parameter')) or (kind=='fold' and label=='second'):
            assert not value.contains(0)
    ev=lambda expr:ball_function(expr,variables)(box)
    if kind=='triple':
        curvature=-ev(s.diff(p,y))*ev(s.diff(p,x,3))/(ev(s.diff(p,y))*ev(s.diff(p,x,e))-ev(s.diff(p,e))*ev(s.diff(p,x,y)))
    else:
        curvature=-(ev(s.diff(p,y,2))-ev(s.diff(p,x,y))**2/ev(s.diff(p,x,2)))/ev(s.diff(p,e))
    assert not curvature.contains(0)
    q=box[1] if name=='even' else box[1].sqrt()
    omega=box[0] if name=='even' else (-box[0]).sqrt()
    lam=acb(-4,omega)
    d,c,left,right,_=parts(4)
    v={'one':left,'even':(left+right)/2,'odd':(left-right)/2}[name]
    matrix=acb_mat([[acb(int(d[i,j].real))+acb(0,q)*(int(c[i,j].imag)+box[2]*arb(float(v[i,j].imag)))-int(i==j)*lam for j in range(24)] for i in range(24)])
    numeric=np.array([[complex(matrix[i,j]) for j in range(24)] for i in range(24)])
    cols=qr(numeric,pivoting=True)[2][:23]
    rowids=qr(numeric[:,cols].T,pivoting=True)[2][:23]
    minor=acb_mat([[matrix[int(i),int(j)] for j in cols] for i in rowids]).det()
    assert not minor.contains(0)
    # The other exact factor must not share this eigenvalue.
    p4,(xx,yy,a,b)=factor_expression(4)
    av,bv=(1+box[2],arb(1)) if name=='one' else ((1+box[2]/2,1+box[2]/2) if name=='even' else (1+box[2]/2,1-box[2]/2))
    spectator=ball_function(p4,(xx,yy,a,b))([-omega*omega,q*q,av,bv])
    assert not spectator.contains(0)
    if name=='even':
        other=ball_function(p,variables)([-box[0],box[1],box[2]])
        assert not other.contains(0)
    return dict(profile=name,kind=kind,center=[str(v) for v in center],radius=str(radius),q=str(q),omega=str(omega),epsilon=str(box[2]),epsilon_curvature=str(curvature),allowed_real_side='epsilon > endpoint' if curvature>0 else 'epsilon < endpoint',contraction_rows=[str(v) for v in rows],image_bounds=[str(v) for v in images],nonzero_derivatives=nonzero,rank23_minor=str(minor),spectator_factor=str(spectator),success=True,displaced_control_rejected=True)


def main():
    ctx.prec=384
    output=dict(scope=__doc__,precision_bits=ctx.prec,roots=[])
    for name,kind,guess in TARGETS:
        result=certify(name,kind,guess);output['roots'].append(result)
        print(name,kind,result['epsilon'],'certified',flush=True)
    paths=[Path(__file__),ROOT/'simulations/route_b_n4_range_polynomial.py',ROOT/'simulations/results/route_b_n4_range_polynomial.json']
    output['source_hashes']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (ROOT/'simulations/results/route_b_n4_range_certificate.json').write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':main()
