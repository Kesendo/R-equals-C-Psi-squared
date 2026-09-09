"""Uniform real-ball enclosure of the negative-sqrt(5) N4 infinity germ.

Certifies a tail at large q, not its connection to the finite fold.
"""
import hashlib,json
from pathlib import Path
import sympy as s
from flint import arb,acb,acb_mat,ctx
from route_b_n4_fold_infinity import derive
from route_b_n4_range_certificate import ball_function
from route_b_n4_range_polynomial import factor_expression,even_frequency_factor


def main():
    ctx.prec=384;data=derive();r,z,k,u,v,h=s.symbols('r z k u v h')
    g=s.sympify(data['compactified_polynomial'])
    K=s.Poly(s.cancel(g.subs({r:1+z*v,k:3+z*u})/(32*z**3)),u,v,z).as_expr()
    variables=(u,v,z);eq=[K,s.diff(K,v)]
    funcs=[ball_function(a,variables) for a in eq]
    jf=[[ball_function(s.diff(a,b),variables) for b in [u,v]] for a in eq]
    center=[(-5*arb(5).sqrt()/4).mid(),(-arb(5).sqrt()/16).mid()]
    J0=acb_mat([[f(center+[arb(0)]) for f in row] for row in jf]).inv()
    Y=acb_mat([[acb(J0[i,j].real.mid()) for j in range(2)] for i in range(2)])
    assert not Y.det().contains(0)
    R=arb('1e-6');rho=arb('1e-4');zz=R/2+arb(0,R/2)
    def check(c):
        box=[x+arb(0,rho) for x in c]+[zz]
        J=acb_mat([[f(box) for f in row] for row in jf])
        defect=acb_mat([[1,0],[0,1]])-Y*J
        residual=Y*acb_mat([[f(c+[zz])] for f in funcs])
        rows=[sum((abs(defect[i,j]).upper() for j in range(2)),arb(0)) for i in range(2)]
        images=[abs(residual[i,0]).upper()+rho*rows[i] for i in range(2)]
        return all(x<1 for x in rows) and all(x<rho for x in images),box,rows,images
    ok,box,rows,images=check(center);assert ok
    assert not check([center[0]+arb('0.1'),center[1]])[0]
    def even_h(expr):
        poly=s.Poly(s.expand(expr),h)
        assert all(p[0]%2==0 for p,c in poly.terms())
        return s.expand(sum(c*(3+z*u)**(p[0]//2) for p,c in poly.terms()))
    p4,(x,y,a,b)=factor_expression(4)
    spectator=s.cancel(even_h(z**8*p4.subs({x:-(1/z+v)**2,y:1/z**2,a:h/2,b:h/2}))/z**2)
    f,(w,q,e)=even_frequency_factor()
    other=even_h(z**8*f.subs({w:-(1/z+v),q:1/z,e:h-2}))
    bounds={}
    for name,expr in [('K_u',s.diff(K,u)),('K_vv',s.diff(K,v,2)),('scaled_P4',spectator),('scaled_other_octic',other)]:
        value=ball_function(expr,variables)(box);assert not value.contains(0)
        bounds[name]=str(value)
    root=Path(__file__).resolve().parents[1]
    paths=[Path(__file__),root/'simulations/route_b_n4_fold_infinity.py',root/'simulations/route_b_n4_range_polynomial.py',root/'simulations/results/route_b_n4_range_polynomial.json',root/'simulations/route_b_n4_range_certificate.py']
    out=dict(scope=__doc__,precision_bits=ctx.prec,z_interval='[0, 1e-6]',q_range='q >= 1000000',u_center=str(center[0]),v_center=str(center[1]),radius_uv=str(rho),
        row_bounds=[str(a) for a in rows],image_bounds=[str(a) for a in images],nonzero_bounds=bounds,success=True,displaced_center_control_rejected=True,
        source_hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (root/'simulations/results/route_b_n4_infinity_ball.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
