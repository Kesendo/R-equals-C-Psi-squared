"""Exact full N4 (1,2) characteristic polynomial for end weights (a,1,b).

P(x,y,a,b)=det(lambda I-L), x=(lambda+4)^2, y=qCSharp^2,
gamma=1, XY. Sparse records are [x power,y power,a power,b power,integer].
No branch count or endpoint certificate is asserted by this producer.
"""
import hashlib,json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from route_b_other_n_unfolding import parts,independent_hop

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'simulations/results/route_b_n4_range_polynomial.json'

def derive():
    d,c,left,right,_=parts(4)
    h,j,k=s.symbols('h j k')
    ring=s.ZZ.poly_ring(h,j,k)
    centered=s.Matrix(24,24,lambda r,col:
        int((d[r,col]+4*(r==col)).real)
        +int((left[r,col]/1j).real)*h
        +int(((c-left-right)[r,col]/1j).real)*j
        +int((right[r,col]/1j).real)*k)
    coeff=DomainMatrix.from_Matrix(centered).convert_to(ring).charpoly()
    records=[]
    for index,p in enumerate(coeff):
        assert index%2==0 or not p
        for powers,value in p.to_dict().items():
            assert all(v%2==0 for v in powers)
            ph,pj,pk=(v//2 for v in powers)
            degree=ph+pj+pk
            records.append([(24-index)//2,degree,2*ph,2*pk,int(value)*(-1)**degree])
    return sorted(records,reverse=True),centered

def expression(records=None):
    """Return exact SymPy P and ordered symbols (x,y,a,b)."""
    if records is None:records=json.loads(OUTPUT.read_text())['terms']
    variables=s.symbols('x y a b')
    return sum(v*s.prod(z**p for z,p in zip(variables,powers))
               for *powers,v in records),variables

def evaluator(records=None):
    """NumPy evaluator returning P,P_x,P_xx,P_y; use exact expression for proofs.

    Double precision can suffer cancellation near multiple roots. Evaluation
    is an endpoint locator, not an interval or exact root certificate.
    """
    p,(x,y,a,b)=expression(records)
    return s.lambdify((x,y,a,b),(p,s.diff(p,x),s.diff(p,x,2),s.diff(p,y)),
                      modules='numpy',cse=True)

def factor_expression(degree=8):
    """Return exact generic factor of requested x degree and (x,y,a,b)."""
    data=json.loads(OUTPUT.read_text())
    return expression(next(f['terms'] for f in data['factors'] if f['degree_x']==degree))

def profile_evaluator(name,modules='numpy',degree=8):
    """Return (values,jacobian) at (x,y,e), for P,Px,Pxx,Py of chosen factor.

    Jacobian is the 4x3 derivative matrix; mpmath supports high precision.
    """
    p,(x,y,a,b)=factor_expression(degree);e=s.Symbol('e')
    av,bv={'one':(1+e,1),'even':(1+e/2,1+e/2),'odd':(1+e/2,1-e/2)}[name]
    p=s.Poly(p.subs({a:av,b:bv}),x,y,e).as_expr()
    values=s.Matrix([p,s.diff(p,x),s.diff(p,x,2),s.diff(p,y)])
    return (s.lambdify((x,y,e),values,modules=modules,cse=True),
            s.lambdify((x,y,e),values.jacobian([x,y,e]),modules=modules,cse=True))

def even_frequency_factor():
    """Return equal-end octic carrying lambda=-4+2iq0, with (w,q,e).

    w is frequency: lambda=-4+i*w. The other octic is its w-sign image.
    """
    p,(x,y,a,b)=factor_expression();w,q,e=s.symbols('w q e')
    poly=s.Poly(p.subs({x:-w*w,y:q*q,a:1+e/2,b:1+e/2}),w,q,e).as_expr()
    constant,factors=s.factor_list(poly)
    assert constant==1 and len(factors)==2
    assert all(m==1 and s.degree(f,w)==8 for f,m in factors)
    seed=s.Poly(3*q**4+q**2-1,q)
    selected=[f for f,m in factors
              if s.rem(s.Poly(f.subs({w:2*q,e:0}),q),seed).is_zero]
    assert len(selected)==1
    assert s.Poly(poly-selected[0]*selected[0].subs(w,-w),w,q,e).is_zero
    return selected[0],(w,q,e)

def main():
    records,centered=derive();p,(x,y,a,b)=expression(records)
    z=s.Symbol('z');checks=[];mutant_detected=False
    d,_,_,_,_=parts(4)
    # Independent full Hilbert-space Hamiltonian followed by coherence restriction.
    for qv,av,bv,zv in [(1,1,1,3),(2,2,3,5),(1,0,2,7),(3,-1,2,2)]:
        hop=independent_hop(4,[av,1,bv])
        full=s.Matrix(24,24,lambda r,col:
            s.Integer(int(d[r,col].real))+s.I*qv*int(hop[r,col].imag))
        actual=(zv*s.eye(24)-full).det(method='domain-ge')
        expected=p.subs({x:(zv+4)**2,y:qv*qv,a:av,b:bv})
        assert actual==expected
        checks.append(dict(q=qv,a=av,b=bv,lambda_value=zv,residual='0'))
        mutant=p.subs({x:(zv+4)**2,y:qv*qv,a:0,b:bv})
        mutant_detected|=mutant!=actual
    assert mutant_detected
    assert s.Poly(p,x).degree()==12
    # Exact reflection and independent bond-sign symmetries.
    assert s.Poly(p-p.xreplace({a:b,b:a}),x,y,a,b).is_zero
    constant,factors=s.factor_list(p)
    assert constant==1 and sorted(s.degree(f,x) for f,m in factors)==[4,8]
    assert s.Poly(s.prod(f**m for f,m in factors)-p,x,y,a,b).is_zero
    factor_records=[dict(degree_x=int(s.degree(f,x)),multiplicity=m,
                         terms=[list(powers)+[int(v)] for powers,v in s.Poly(f,x,y,a,b).terms()])
                    for f,m in factors]
    output=dict(scope=__doc__,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        variables=['x=(lambda+4)^2','y=qCSharp^2','a=left end weight','b=right end weight'],
        degree_x=12,term_count=len(records),terms=records,factors=factor_records,
        exact_determinant_checks=checks,missing_left_bond_control_rejected=bool(mutant_detected))
    OUTPUT.write_text(json.dumps(output,indent=2)+'\n')
    print('Exact degree12 characteristic; terms',len(records),'independent determinants',len(checks))

if __name__=='__main__':main()
