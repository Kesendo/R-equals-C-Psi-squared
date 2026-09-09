"""Exact N4 end-response coefficients over Q(t), 3t^4-t^2-1=0.

The embedding is t=i*q0, q0=sqrt((sqrt(13)-1)/6), gamma=1.
Run: python simulations/route_b_n4_self_fold.py
"""
import json
import hashlib
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from route_b_other_n_unfolding import parts, independent_hop

ROOT = Path(__file__).resolve().parents[1]


def main():
    q0 = s.sqrt((s.sqrt(13)-1)/6)
    field = s.QQ.algebraic_field(s.I*q0)
    t = field.unit
    zero = field.zero
    def mat(a):
        return DM.from_Matrix(s.Matrix(a)).convert_to(field)
    def scalar(a, x):
        return a.scalarmul(x)
    def trace(a):
        rows = a.to_list()
        return sum((rows[j][j] for j in range(a.shape[0])), zero)
    def delta(a):
        return 2*trace(a*a)-trace(a)**2
    def polar(a,b):
        return 4*trace(a*b)-2*trace(a)*trace(b)
    def polynomial(x):
        # Field elements are polynomials in this one chosen generator t.
        coeffs = x.to_list()
        return s.Poly.from_list([s.Rational(v.numerator,v.denominator) for v in coeffs],s.Symbol('t')).as_expr()
    def real_form(x):
        p = s.Poly(polynomial(x),s.Symbol('t'))
        assert all(power[0] % 2 == 0 for power,_ in p.terms())
        return s.simplify(p.as_expr().subs(s.Symbol('t')**2,(1-s.sqrt(13))/6))
    d,c,left,right,reflection = parts(4)
    def integer(a):
        assert all(v.imag == 0 and v.real == int(v.real) for row in a for v in row)
        return [[int(v.real) for v in row] for row in a]
    assert (c+left/8-right/4 == independent_hop(4,[1.125,1,.75])).all()
    D = mat(integer(d))
    T = mat(integer(c/1j))
    Vl = mat(integer(left/1j)); Vr = mat(integer(right/1j))
    I = DM.eye(24,field).to_dense()
    L = D+scalar(T,t)
    lam = -field(4)+2*t
    M = L-scalar(I,lam)
    U = M.nullspace().transpose().to_dense()
    assert U.shape == (24,2)
    gram = U.transpose()*U
    W = gram.inv()*U.transpose()
    P = U*W; Q = I-P
    S = Q*(scalar(I,lam)-L+P).inv()*Q
    assert (M*U).is_zero_matrix
    assert (P*P-P).is_zero_matrix
    assert (M*M).rank() == 22
    # Bra complement is an internal permutation precisely at this N.
    from framework.weight_coherence_block import weight_block_configs
    basis = [(a,b) for a in weight_block_configs(4,1) for b in weight_block_configs(4,2)]
    idx = {v:j for j,v in enumerate(basis)}
    B = mat([[int(i == idx[(a,b^15)]) for a,b in basis] for i in range(24)])
    assert (B*D*B+D+scalar(I,field(8))).is_zero_matrix
    assert all((B*v*B-v).is_zero_matrix for v in (T,Vl,Vr))
    # Mutating the dissipator at one coherence must break that same fold equation.
    wrong = D+mat([[int(i==j==0) for j in range(24)] for i in range(24)])
    assert not (B*wrong*B+wrong+scalar(I,field(8))).is_zero_matrix
    R = mat([[int(i==reflection[j]) for j in range(24)] for i in range(24)])
    assert (R*U-U).is_zero_matrix
    G = mat([[int(i==j)*(-1)**((a&10).bit_count()+(b&10).bit_count()) for j,(a,b) in enumerate(basis)] for i in range(24)])
    assert (G*D*G-D).is_zero_matrix
    assert all((G*v*G+v).is_zero_matrix for v in (T,Vl,Vr))
    A = W*scalar(T,t)*U  # Relative q displacement: q=q0*(1+r*eta).
    alpha = delta(A)
    records = {}
    root13 = s.sqrt(13)
    expected = {
        'one': ((663-131*root13)/299, (87-25*root13)/26),
        'even': ((663-131*root13)/299, (87-25*root13)/26),
        'odd': (-(11323+557*root13)/1196, -(12373651+3475859*root13)/9568),
    }
    for name,V in [('one',Vl),('even',scalar(Vl+Vr,field(s.Rational(1,2)))),('odd',scalar(Vl-Vr,field(s.Rational(1,2))))]:
        E = scalar(V,t)
        B1 = W*E*U
        response = W*E*S*E*U if name=='odd' else B1
        if name=='odd': assert B1.is_zero_matrix
        beta = polar(A,response); gamma = delta(response)
        b = real_form(beta/alpha); g = real_form(gamma/alpha)
        assert s.simplify(b-expected[name][0]) == 0
        assert s.simplify(g-expected[name][1]) == 0
        disc = s.simplify(b*b-4*g)
        assert alpha != zero and disc.is_positive is True and g.is_negative is True
        slopes = [s.simplify(q0*(-b-sign*s.sqrt(disc))/2) for sign in (1,-1)]
        records[name] = dict(relative_polynomial=dict(b=str(b),g=str(g),discriminant=str(disc)),q_slopes=[str(v) for v in slopes],q_slopes_decimal=[str(s.N(v,16)) for v in slopes])
    # A self-folded complex-symmetric family can instead exchange EP germs.
    x,e = s.symbols('x e',real=True)
    control = s.I*s.Matrix([[x,e],[e,-x]])
    assert control == -s.conjugate(control)
    control_disc = s.expand(2*s.trace(control**2)-s.trace(control)**2)
    assert control_disc == -4*(x*x+e*e)
    assert s.discriminant(s.Symbol('r')**2+1,s.Symbol('r')) == -4
    dependencies = [ROOT/'simulations/route_b_other_n_unfolding.py',ROOT/'simulations/framework/weight_coherence_block.py']
    out = dict(scope=__doc__,q0=str(q0),lambda0='-4+2*i*q0',dimension=24,nullity=2,squared_nullity=2,alpha_relative=str(polynomial(alpha)),profiles=records,fold_mutation_rejected=True,swapped_branch_control='q=1 +/- i*epsilon, lambda=-4',dependencies={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies},script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    path=ROOT/'simulations/results/route_b_n4_self_fold.json'
    path.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
