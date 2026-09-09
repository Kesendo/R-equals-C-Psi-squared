"""Isolated real crossing of the two equal-end N4 frequency octics.

384-bit contraction certifies F=F_w=F(-w)=0 locally. This is not a
completeness certificate for the incoming arm, nor an EP3 assertion.
"""
import hashlib,json
from pathlib import Path
import mpmath as mp
import sympy as s
import numpy as np
from scipy.linalg import qr
from flint import arb,acb,acb_mat,ctx
from route_b_n4_range_polynomial import even_frequency_factor,factor_expression,expression
from route_b_n4_range_certificate import ball_function
from route_b_n4_bridge_certificate import pack_ball
from route_b_other_n_unfolding import parts
ROOT=Path(__file__).resolve().parents[1]

def certify(guess):
    ctx.prec=384;mp.mp.dps=100
    f,variables=even_frequency_factor();w,q,e=variables
    other=f.subs(w,-w)
    eq=s.Matrix([f,s.diff(f,w),other]);jac=eq.jacobian(variables)
    fn=s.lambdify(variables,eq,'mpmath',cse=True)
    jn=s.lambdify(variables,jac,'mpmath',cse=True)
    root=mp.findroot(lambda *v:tuple(fn(*v)),tuple(map(mp.mpf,guess)),J=jn,tol=mp.mpf('1e-90'))
    center=[arb(str(v)).mid()for v in root];radius=arb('1e-30')
    funcs=[ball_function(expr,variables)for expr in eq]
    jf=[[ball_function(jac[i,j],variables)for j in range(3)]for i in range(3)]
    inverse=acb_mat([[v(center)for v in row]for row in jf]).inv()
    Y=acb_mat([[inverse[i,j].real.mid()for j in range(3)]for i in range(3)])
    assert not Y.det().contains(0)
    def check(c):
        box=[v+arb(0,radius)for v in c]
        J=acb_mat([[v(box)for v in row]for row in jf])
        defect=acb_mat([[int(i==j)for j in range(3)]for i in range(3)])-Y*J
        residual=Y*acb_mat([[v(c)]for v in funcs])
        rows=[sum((abs(defect[i,j]).upper()for j in range(3)),arb(0))for i in range(3)]
        images=[abs(residual[i,0]).upper()+rows[i]*radius for i in range(3)]
        return all(v<1 for v in rows)and all(v<radius for v in images),box,rows,images
    success,box,rows,images=check(center);assert success
    displaced=center.copy();displaced[2]+=arb('.001');assert not check(displaced)[0]
    p4,(x,y,a,b)=factor_expression(4)
    spectator=p4.subs({x:-w*w,y:q*q,a:1+e/2,b:1+e/2})
    full,vs=expression();xx,yy,aa,bb=vs
    full_frequency=full.subs({xx:-w*w,yy:q*q,aa:1+e/2,bb:1+e/2})
    assert s.Poly(full_frequency-spectator*f*other,*variables).is_zero
    gates={name:ball_function(expr,variables)(box)for name,expr in {
        'F_ww':s.diff(f,w,2),'other_w':s.diff(other,w),'P4':spectator,
        'system_Jacobian_determinant':jac.det(),'F_epsilon':s.diff(f,e)}.items()}
    assert all(not v.contains(0)for v in gates.values())
    d,c,left,right,reflection=parts(4)
    matrix=acb_mat([[acb(int(d[i,j].real)+4*int(i==j))+acb(0,box[1])*(int(c[i,j].imag)+box[2]*arb(int((left[i,j]+right[i,j]).imag))/2)-acb(0,box[0])*int(i==j)for j in range(24)]for i in range(24)])
    numeric=np.array([[complex(matrix[i,j])for j in range(24)]for i in range(24)])
    cols=qr(numeric,pivoting=True)[2][:22];rowids=qr(numeric[:,cols].T,pivoting=True)[2][:22]
    minor=acb_mat([[matrix[int(i),int(j)]for j in cols]for i in rowids]).det();assert not minor.contains(0)
    data=dict(scope=__doc__,precision_bits=ctx.prec,coordinates=['omega','qCSharp','epsilon'],center=[str(v)for v in root],box_dyadic=[pack_ball(v)for v in box],radius_dyadic=pack_ball(radius),contraction_rows_dyadic=[pack_ball(v)for v in rows],image_bounds_dyadic=[pack_ball(v)for v in images],nonzero_gates={k:dict(display=str(v),dyadic=pack_ball(v))for k,v in gates.items()},full_characteristic_factorization_exact=True,full_eigenvalue_algebraic_multiplicity=3,rank_lower_bound=22,rank22_minor=dict(rows=list(map(int,rowids)),columns=list(map(int,cols)),real=pack_ball(minor.real),imag=pack_ball(minor.imag),display=str(minor)),jordan_claim='J2 plus J1; exact parity factors give one eigenvector in each sector, and the nonzero rank22 minor excludes any additional eigenvector.',full_nullity=2,jordan_block_sizes=[2,1],displaced_control_rejected=True,success=True)
    paths=[Path(__file__),ROOT/'simulations/route_b_n4_range_polynomial.py',ROOT/'simulations/results/route_b_n4_range_polynomial.json',ROOT/'simulations/route_b_n4_range_certificate.py',ROOT/'simulations/route_b_n4_bridge_certificate.py',ROOT/'simulations/route_b_other_n_unfolding.py',ROOT/'simulations/framework/weight_coherence_block.py']
    data['source_hashes']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    print('PASS crossing',data['center'],flush=True);print({k:str(v)for k,v in gates.items()},flush=True)
    return data

def main():
    from route_b_n4_parity_factors import derive
    ownership=derive()
    roots=[certify(guess)for guess in [('1.3136606858010466953','.66512116752370169947','-.011642506556512168505'),('1.231','.99','-.375')]]
    path=ROOT/'simulations/route_b_n4_parity_factors.py'
    out=dict(scope=__doc__,roots=roots,parity_ownership=ownership,parity_producer_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),completeness_claim=False)
    (ROOT/'simulations/results/route_b_n4_incoming_crossing.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
