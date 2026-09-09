"""Exact spatial-reflection ownership of N4 equal-end characteristic factors.

Gamma=1, XY, end weights (1+epsilon/2,1,1+epsilon/2), lambda=-4+i*omega.
Build the full coherence matrix and exact integer reflection orbit bases.
This factor identity is not a numerical root or Jordan-character certificate.
"""
import hashlib,json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from route_b_other_n_unfolding import parts
from route_b_n4_range_polynomial import even_frequency_factor,factor_expression

def derive():
    F,(w,q,e)=even_frequency_factor();d,c,left,right,reflection=parts(4)
    assert all(value.imag==0 and value.real==int(value.real) for row in d for value in row)
    assert all(value.real==0 and value.imag==int(value.imag)
               for matrix in [c,left,right] for row in matrix for value in row)
    # B=-i(L+4I) has eigenvalue omega. All entries built exactly from integers.
    B=s.Matrix(24,24,lambda i,j:
        -s.I*int((d[i,j]+4*(i==j)).real)
        +q*(int(c[i,j].imag)+e*s.Rational(1,2)*int((left[i,j]+right[i,j]).imag)))
    R=s.zeros(24)
    for j,i in enumerate(reflection):R[i,j]=1
    assert R*B==B*R
    reps=[j for j in range(24) if j<reflection[j]]
    assert len(reps)==12
    domain=s.QQ_I.poly_ring(q,e);records=[];products=[];quartics=[]
    for parity in [1,-1]:
        U=s.zeros(24,12)
        for col,j in enumerate(reps):U[j,col]=1;U[reflection[j],col]=parity
        assert U.T*U==2*s.eye(12) and R*U==parity*U
        small=s.expand(s.Rational(1,2)*U.T*B*U)
        assert (B*U-U*small).applyfunc(s.expand)==s.zeros(24,12)
        coeff=DomainMatrix.from_Matrix(small).convert_to(domain).charpoly()
        chi=s.Poly.from_list([domain.to_sympy(a) for a in coeff],w).as_expr()
        chosen=F if parity==1 else F.subs(w,-w)
        quotient,remainder=s.div(chi,chosen,w)
        assert s.expand(remainder)==0 and s.degree(quotient,w)==4
        assert s.Poly(chi-chosen*quotient,w,q,e).is_zero
        records.append(dict(parity=parity,dimension=12,frequency_octic='F(omega,q,epsilon)' if parity==1 else 'F(-omega,q,epsilon)',
                            quartic=str(s.factor(quotient)),exact_intertwining=True,exact_division_remainder='0'))
        products.append(chi);quartics.append(quotient)
    P4,(x,y,a,b)=factor_expression(4)
    assert s.Poly(quartics[0]*quartics[1]-P4.subs({x:-w*w,y:q*q,a:1+e/2,b:1+e/2}),w,q,e).is_zero
    # Explicit wrong-parity assignment fails on both sectors.
    assert all(s.Poly(s.rem(products[j],F.subs(w,-w) if j==0 else F,w),w,q,e).as_expr()!=0 for j in [0,1])
    return dict(scope=__doc__,sectors=records,quartic_product_equals_P4=True,wrong_parity_control_rejected=True,
                character_rule='At an F=Fw=0, Fww!=0, Fe!=0 point with F(-w)=0 and its w derivative nonzero, and both quartics nonzero: R+ has EP2 and R- a simple eigenvalue; full character is J2 plus J1, not J3.')

def main():
    out=derive();root=Path(__file__).resolve().parents[1]
    paths=[Path(__file__),root/'simulations/route_b_other_n_unfolding.py',
           root/'simulations/route_b_n4_range_polynomial.py',
           root/'simulations/results/route_b_n4_range_polynomial.json',
           root/'simulations/framework/weight_coherence_block.py']
    out['source_hashes']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (root/'simulations/results/route_b_n4_parity_factors.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
