"""Validated N4 equal-end EP2 bridge: q=2 to the negative-sqrt(5) infinity germ.

The current repository's F163/F89d, self-fold and second-arm notes supply the
local objects; this adds validated connection, not a new spectral mechanism.
The named-store survey also includes proofs, experiments, Confirmations,
GLOSSARY, OpenArcs and CAUGHT_ERRORS: endpoint existence alone is not a bridge.

Every tile proves a contraction uniformly in z. The preceding endpoint's
certified enclosure must be contained in the next uniqueness box. Mere box
overlap is never used as connection evidence. Removing a tile or substituting
the other negative infinity germ must fail the connection audit.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
import sympy as s
from flint import arb, acb, acb_mat, ctx
from route_b_n4_fold_infinity import derive
from route_b_n4_range_certificate import ball_function
from route_b_n4_range_polynomial import factor_expression, even_frequency_factor

ROOT=Path(__file__).resolve().parents[1]


def ab(value):
    if isinstance(value,Fraction):return arb(value.numerator)/arb(value.denominator)
    return arb(str(value))


def pack_ball(value):
    """Lossless dyadic endpoints; decimal display strings are never audited."""
    mid=value.mid().man_exp();rad=value.rad().man_exp()
    record={'mid':[str(mid[0]),int(mid[1])],'rad':[str(rad[0]),int(rad[1])]}
    return record


def rational_bounds(record):
    def dyadic(pair):return Fraction(int(pair[0]))*Fraction(2)**pair[1]
    mid=dyadic(record['mid']);rad=dyadic(record['rad'])
    return mid-rad,mid+rad


def contains_dyadic(outer,inner):
    a,b=rational_bounds(outer);c,d=rational_bounds(inner)
    return a<=c and d<=b


def setup():
    data=derive();r,k,z,u,v,h=s.symbols('r k z u v h')
    g=s.sympify(data['compactified_polynomial'])
    K=s.Poly(s.cancel(g.subs({r:1+z*v,k:3+z*u})/(32*z**3)),u,v,z).as_expr()
    eq=s.Matrix([K,s.diff(K,v)]);jac=eq.jacobian([u,v])
    assert eq.subs({u:-2,v:0,z:s.Rational(1,2)})==s.zeros(2,1)
    assert eq.subs({u:-5*s.sqrt(5)/4,v:-s.sqrt(5)/16,z:0}).applyfunc(s.simplify)==s.zeros(2,1)
    def even_h(expr):
        p=s.Poly(s.expand(expr),h)
        assert all(pw[0]%2==0 for pw,c in p.terms())
        return s.expand(sum(c*(3+z*u)**(pw[0]//2) for pw,c in p.terms()))
    p4,(x,y,a,b)=factor_expression(4)
    spectator=s.cancel(even_h(z**8*p4.subs({x:-(1/z+v)**2,y:1/z**2,a:h/2,b:h/2}))/z**2)
    f,(w,q,e)=even_frequency_factor()
    other=even_h(z**8*f.subs({w:-(1/z+v),q:1/z,e:h-2}))
    expressions={'K_u':s.diff(K,u),'K_vv':s.diff(K,v,2),
                 'scaled_P4':spectator,'scaled_other_octic':other,
                 'k':3+z*u,'r':1+z*v,'u':u}
    variables=(u,v,z)
    return ([ball_function(a,variables) for a in eq],
            [[ball_function(jac[i,j],variables) for j in range(2)]for i in range(2)],
            {name:ball_function(a,variables) for name,a in expressions.items()},
            s.lambdify(variables,eq,'mpmath',cse=True),s.lambdify(variables,jac,'mpmath',cse=True))


def certify_tile(high,low,known,functions):
    funcs,jf,gates,nf,nj=functions
    mid=(high+low)/2
    def solve(zz,guess):
        zm=mp.mpf(zz.numerator)/zz.denominator
        return mp.findroot(lambda a,b:tuple(nf(a,b,zm)),guess,
                           J=lambda a,b:nj(a,b,zm),tol=mp.mpf('1e-55'),maxsteps=35)
    guess=tuple(mp.mpf(str(float(x.mid()))) for x in known)
    values=[]
    for zz in [high,mid,low]:
        value=solve(zz,guess);values.append(value);guess=tuple(value)
    center=[ab(values[1][i]).mid() for i in range(2)]
    spread=max(abs(value[i]-values[1][i]) for value in values for i in range(2))
    radius=ab(spread*mp.mpf('1.6')+mp.mpf('1e-12'))
    zz=ab(mid)+arb(0,ab((high-low)/2))
    J0=acb_mat([[f(center+[ab(mid)]) for f in row] for row in jf]).inv()
    Y=acb_mat([[acb(J0[i,j].real.mid())for j in range(2)]for i in range(2)])
    if Y.det().contains(0):raise ValueError('singular preconditioner')
    identity=acb_mat([[1,0],[0,1]])
    residual=Y*acb_mat([[f(center+[zz])]for f in funcs])
    # Direct interval polynomial evaluation can overestimate z-dependence;
    # enlarge the radius to cover that rigorous residual, not only the
    # approximate endpoint displacement used for initial box selection.
    radius=max(radius,2*max(abs(residual[i,0]).upper()for i in range(2)))
    uv=[c+arb(0,radius) for c in center]
    if not all(uv[i].contains(known[i]) for i in range(2)):
        raise ValueError('previous endpoint enclosure not inside new box')
    box=uv+[zz]
    defect=identity-Y*acb_mat([[f(box) for f in row] for row in jf])
    rows=[sum((abs(defect[i,j]).upper()for j in range(2)),arb(0))for i in range(2)]
    images=[abs(residual[i,0]).upper()+radius*rows[i]for i in range(2)]
    if not(all(x<1 for x in rows)and all(x<radius for x in images)):
        raise ValueError('uniform contraction/inclusion failed')
    bounds={name:fn(box)for name,fn in gates.items()}
    for name in ['K_u','K_vv','scaled_P4','scaled_other_octic']:
        if bounds[name].contains(0):raise ValueError('nonzero gate '+name)
    if not(bounds['k']>0 and bounds['k']<4 and bounds['r']>0 and bounds['u']<0):
        raise ValueError('physical domain gate')
    # Refine an enclosure of the already proved lower-endpoint root. Each
    # mean-value image contains it; intersecting with its previous enclosure
    # preserves that fact. No independent floating root identity is assumed.
    endpoint=uv
    fcenter=Y*acb_mat([[f(center+[ab(low)])]for f in funcs])
    for _ in range(8):
        # The mean-value segment begins at the original fixed center, which
        # need not belong to the shrinking endpoint enclosure. The full uv
        # box contains that segment at every iteration.
        de=identity-Y*acb_mat([[f(uv+[ab(low)])for f in row]for row in jf])
        delta=acb_mat([[endpoint[i]-center[i]]for i in range(2)])
        value=acb_mat([[center[i]]for i in range(2)])-fcenter+de*delta
        endpoint=[endpoint[i].intersection(value[i,0].real)for i in range(2)]
    record=dict(z_high=str(high),z_low=str(low),center=[str(c)for c in center],
                radius=str(radius),box=[str(x)for x in uv],
                incoming_endpoint=[str(x)for x in known],endpoint=[str(x)for x in endpoint],
                box_dyadic=[pack_ball(x)for x in uv],
                endpoint_dyadic=[pack_ball(x)for x in endpoint],
                center_dyadic=[pack_ball(x)for x in center],radius_dyadic=pack_ball(radius),
                preconditioner=[[str(Y[i,j].real)for j in range(2)]for i in range(2)],
                preconditioner_dyadic=[[pack_ball(Y[i,j].real)for j in range(2)]for i in range(2)],
                row_bounds=[str(x)for x in rows],image_bounds=[str(x)for x in images],
                row_bounds_dyadic=[pack_ball(x)for x in rows],
                image_bounds_dyadic=[pack_ball(x)for x in images],
                nonzero_and_domain_bounds={name:value.str(18,more=True)for name,value in bounds.items()},
                nonzero_and_domain_bounds_dyadic={name:pack_ball(value)for name,value in bounds.items()},
                incoming_endpoint_contained=True)
    return record,endpoint


def audit_connection(tiles,wrong_germ=False):
    if not tiles or Fraction(tiles[0]['z_high'])!=Fraction(1,2):return False
    endpoint=[pack_ball(arb(-2)),pack_ball(arb(0))];current=Fraction(1,2)
    for tile in tiles:
        if Fraction(tile['z_high'])!=current:return False
        low=Fraction(tile['z_low'])
        if not 0<=low<current:return False
        box=tile['box_dyadic']
        if not all(contains_dyadic(box[i],endpoint[i])for i in range(2)):return False
        radius_lower=rational_bounds(tile['radius_dyadic'])[0]
        if not all(rational_bounds(x)[1]<1 for x in tile['row_bounds_dyadic']):return False
        if not all(rational_bounds(x)[1]<radius_lower for x in tile['image_bounds_dyadic']):return False
        gates={name:rational_bounds(value)for name,value in tile['nonzero_and_domain_bounds_dyadic'].items()}
        if not all(gates[name][0]*gates[name][1]>0 for name in ['K_u','K_vv','scaled_P4','scaled_other_octic']):return False
        if not(gates['k'][0]>0 and gates['k'][1]<4 and gates['r'][0]>0 and gates['u'][1]<0):return False
        endpoint=tile['endpoint_dyadic'];current=low
    if current!=0:return False
    exact=([-2*arb(2).sqrt(),-arb(2).sqrt()/4]if wrong_germ else
           [-5*arb(5).sqrt()/4,-arb(5).sqrt()/16])
    return all(contains_dyadic(tiles[-1]['box_dyadic'][i],pack_ball(exact[i]))for i in range(2))


def main():
    ctx.prec=384;mp.mp.dps=75
    functions=setup();current=Fraction(1,2);step=Fraction(1,256)
    known=[arb(-2),arb(0)];tiles=[];rejections=0
    while current>0:
        low=max(Fraction(0),current-step)
        try:record,endpoint=certify_tile(current,low,known,functions)
        except (ValueError,ZeroDivisionError):
            step/=2;rejections+=1
            if step<Fraction(1,2**24):raise
            continue
        tiles.append(record);known=endpoint;current=low
        if len(tiles)%25==0:print('tiles',len(tiles),'z',float(current),'rejections',rejections,flush=True)
        step=min(step*2,Fraction(1,32))
    assert audit_connection(tiles)
    assert not audit_connection(tiles[:len(tiles)//2]+tiles[len(tiles)//2+1:])
    assert not audit_connection(tiles,wrong_germ=True)
    summaries={}
    for name in tiles[0]['nonzero_and_domain_bounds_dyadic']:
        bounds=[rational_bounds(t['nonzero_and_domain_bounds_dyadic'][name])for t in tiles]
        lower=min(a for a,b in bounds);upper=max(b for a,b in bounds)
        summaries[name]={'lower_exact':str(lower),'upper_exact':str(upper),
                         'lower_decimal':float(lower),'upper_decimal':float(upper)}
    paths=[Path(__file__),ROOT/'simulations/route_b_n4_fold_infinity.py',ROOT/'simulations/route_b_n4_range_polynomial.py',ROOT/'simulations/results/route_b_n4_range_polynomial.json',ROOT/'simulations/route_b_n4_range_certificate.py']
    out=dict(scope=__doc__,precision_bits=ctx.prec,z_interval='[0,1/2]',q_interval='[2,infinity)',
             tile_count=len(tiles),rejected_attempts=rejections,tiles=tiles,
             exact_start='z=1/2,u=-2,v=0',exact_end='z=0,u=-5sqrt(5)/4,v=-sqrt(5)/16',
             global_gate_bounds=summaries,
             gap_control_rejected=True,wrong_infinity_germ_control_rejected=True,success=True,
             proof='Uniform contraction gives a unique root in each tile box for each z. The preceding endpoint enclosure is wholly inside the next box, so uniqueness identifies the branches. Endpoint enclosures are mean-value images of the proved root, never floating identities. At z=0 the exact negative-sqrt(5) root belongs to the last uniqueness box. Nonzero derivatives and spectators prove full EP2 for finite q; u<0,k>0,r>0 keep positive end bonds and epsilon<sqrt(3)-2<0.',
             source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})
    (ROOT/'simulations/results/route_b_n4_bridge_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS',len(tiles),'tiles, gap and wrong-germ controls rejected',flush=True)


if __name__=='__main__':main()
