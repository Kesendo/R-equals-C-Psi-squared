"""Exact local infinity germs of the equal-end N4 frequency octic.

Source: range_polynomial's full-block factor identity; F89d supplies the
self-fold and F163 the distinction between local existence and path transport.
These are local full-block EP2 germs. Association with a finite seed
or its second fold arm is a separate continuation question.
Run: python simulations/route_b_n4_fold_infinity.py
"""
import hashlib
import json
from pathlib import Path
import sympy as s
from route_b_n4_range_polynomial import even_frequency_factor, factor_expression

ROOT = Path(__file__).resolve().parents[1]


def derive():
    f, (w, q, e) = even_frequency_factor()
    r, z, h, k, u, v = s.symbols('r z h k u v')
    gh = s.Poly(s.expand(z**8*f.subs({w:r/z, q:1/z, e:h-2})), h)
    assert all(p[0] % 2 == 0 for p, _ in gh.terms())
    g = s.expand(sum(c*k**(p[0]//2) for p,c in gh.terms()))
    assert s.Poly(g.subs(k,h*h)-gh.as_expr(), r,z,h).is_zero
    local = s.Poly(s.expand(g.subs({r:1+z*v,k:3+z*u})),z)
    assert min(p[0] for p,c in local.terms()) == 3
    blown = s.cancel(local.as_expr()/(32*z**3))
    H = (u+4*v)*(u-4*v)**2+40*v-8*u
    assert s.expand(blown.subs(z,0)-H) == 0
    equations = s.Matrix([blown,s.diff(blown,v)])
    J = equations.jacobian([u,v])
    res = s.factor(s.resultant(H,s.diff(H,v),v))
    assert s.expand(res-1048576*(u*u-8)*(16*u*u-125)) == 0
    p4,(xx,yy,aa,bb)=factor_expression(4)
    def substitute_even_h(expression):
        poly=s.Poly(s.expand(expression),h)
        assert all(p[0]%2==0 for p,c in poly.terms())
        return s.expand(sum(c*(3+z*u)**(p[0]//2) for p,c in poly.terms()))
    spectator=substitute_even_h(z**8*p4.subs({xx:-(1/z+v)**2,yy:1/z**2,aa:h/2,bb:h/2}))
    other=substitute_even_h(z**8*f.subs({w:-(1/z+v),q:1/z,e:h-2}))
    spectator_series=s.Poly(spectator,z)
    assert min(p[0] for p,c in spectator_series.terms())==2
    spectator_lead=s.factor(spectator_series.coeff_monomial(z**2))
    assert s.expand(spectator_lead-16*((u-4*v)**2+64))==0
    assert other.subs(z,0)==-6912
    rows=[]
    for a,b in [(-2*s.sqrt(2),-s.sqrt(2)/4),
                (2*s.sqrt(2),s.sqrt(2)/4),
                (-5*s.sqrt(5)/4,-s.sqrt(5)/16),
                (5*s.sqrt(5)/4,s.sqrt(5)/16)]:
        pt={u:a,v:b,z:0}
        assert all(s.simplify(expr.subs(pt)) == 0 for expr in equations)
        j0=J.subs(pt).applyfunc(s.simplify)
        determinant=s.simplify(j0.det())
        assert determinant != 0
        spectator_value=s.simplify(spectator_lead.subs(pt))
        assert spectator_value>0
        assert s.simplify(s.diff(H,u).subs(pt)) != 0
        assert s.simplify(s.diff(H,v,2).subs(pt)) != 0
        derivative=(-j0.inv()*equations.diff(z).subs(pt)).applyfunc(s.simplify)
        assert all(s.simplify(expr)==0 for expr in j0*derivative+equations.diff(z).subs(pt))
        e1=s.simplify(a/(2*s.sqrt(3)))
        e2=s.simplify(derivative[0]/(2*s.sqrt(3))-a*a/(24*s.sqrt(3)))
        rows.append(dict(u0=str(a),v0=str(b),jacobian_determinant=str(determinant),
            H_u=str(s.simplify(s.diff(H,u).subs(pt))),
            H_vv=str(s.simplify(s.diff(H,v,2).subs(pt))),
            u1=str(derivative[0]),v1=str(derivative[1]),
            epsilon_1_over_q=str(e1),epsilon_1_over_q_squared=str(e2),
            omega_constant=str(b),omega_1_over_q=str(derivative[1])))
        rows[-1]['P4_leading_coefficient']=str(spectator_value)
    # A dropped leading return term must fail the same leading-system roots.
    mutant=H-40*v
    assert all(s.simplify(mutant.subs({u:s.sympify(row['u0']),v:s.sympify(row['v0'])})) != 0 for row in rows)
    return dict(scope=__doc__,coordinates='z=1/q; r=omega/q; k=(epsilon+2)^2; r=1+z*v; k=3+z*u',
        compactified_polynomial=str(s.collect(g,z)),leading_polynomial=str(H),resultant=str(res),
        limiting_epsilon='sqrt(3)-2',germs=rows,mutation_rejected=True,
        spectator_P4_leading=str(spectator_lead),spectator_other_octic_constant='-6912',
        interpretation='Implicit function theorem: each nonzero Jacobian gives a unique real analytic (u(z),v(z)) germ. omega=q+v0+v1/q+O(q^-2); epsilon=sqrt(3)-2+e1/q+e2/q^2+O(q^-3). Nonzero Hvv gives multiplicity two; the spectator factors stay nonzero. Nonzero Hu gives nonzero epsilon derivative of the full determinant at its root, excluding matrix nullity >=2 by the adjugate derivative identity, hence EP2. No finite-path connection or explicit radius certified.')


if __name__ == '__main__':
    output=derive()
    paths=[Path(__file__),ROOT/'simulations/route_b_n4_range_polynomial.py',ROOT/'simulations/results/route_b_n4_range_polynomial.json']
    output['source_hashes']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (ROOT/'simulations/results/route_b_n4_fold_infinity.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output['germs'],indent=2))
