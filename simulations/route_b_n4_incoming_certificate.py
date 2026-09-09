"""Validated equal-end N4 incoming branch: exact diabolic seed to q=2 fold.

Stage 0: F163/F89d, PROOF_ROUTE_B_N4_SELF_FOLD, ROUTE_B_N4_RANGE and
ROUTE_B_N4_SECOND_ARM own the local objects. Proofs, experiments,
Confirmations, GLOSSARY, OpenArcs and CAUGHT_ERRORS distinguish them from a
validated connecting path. This producer supplies that incoming connection.
The seed itself is semisimple. The selected frequency factor has a double
root along the open segment; opposite-parity coincidences are retained.

Run: python simulations/route_b_n4_incoming_certificate.py
"""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
import mpmath as mp
import sympy as s
from flint import arb,acb,acb_mat,ctx
from route_b_n4_bridge_certificate import ab,pack_ball,contains_dyadic,rational_bounds
from route_b_n4_range_certificate import ball_function
from route_b_n4_range_polynomial import even_frequency_factor,factor_expression

ROOT=Path(__file__).resolve().parents[1]


def source_hashes():
    names=['route_b_n4_bridge_certificate.py','route_b_n4_range_polynomial.py',
           'results/route_b_n4_range_polynomial.json','route_b_n4_range_certificate.py',
           'route_b_n4_self_fold.py','results/route_b_n4_self_fold.json',
           'route_b_n4_parity_factors.py','results/route_b_n4_parity_factors.json',
           'route_b_n4_incoming_crossing.py','results/route_b_n4_incoming_crossing.json']
    paths=[Path(__file__)]+[ROOT/'simulations'/name for name in names]
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}


def setup():
    f,(w,q,e)=even_frequency_factor()
    a,t,u,v=s.symbols('a t u v');mod=3*a**4+a*a-1
    def reduce(expression):
        p=s.Poly(s.expand(expression),t,u,v)
        return s.expand(sum(s.rem(co,mod,a)*t**pw[0]*u**pw[1]*v**pw[2]for pw,co in p.terms()))
    raw=reduce(f.subs({q:a+(2-a)*t,w:2*a+t*v,e:t*u}))
    poly=s.Poly(raw,t,u,v)
    assert min(pw[0]for pw,co in poly.terms())==2
    K=s.expand(s.Rational(27,256)*s.cancel(raw/t**2))
    H=s.expand(K.subs(t,0))
    vr=s.solve(s.diff(H,v),v)[0]
    A=4639*a*a-1999
    B=-7032*a**3-2510*a*a+2676*a+1172
    C=4524*a**3-3952*a*a-780*a+403
    seed_polynomial=A*u*u+B*u+C
    old_b=(532-786*a*a)/299;old_g=(31-75*a*a)/13;distance=2-a
    inherited=a*a*old_g*u*u+a*old_b*distance*u+distance**2
    assert s.Poly(reduce(seed_polynomial*(a*a*old_g)-inherited*A),t,u,v).is_zero
    numerator,denominator=s.fraction(s.factor(H.subs(v,vr)))
    assert s.expand(s.rem(numerator,mod,a)+12*seed_polynomial)==0
    assert s.cancel(s.diff(H,v).subs(v,vr))==0
    # Exact endpoint identities in Q(a,sqrt(2)), a the positive seed root.
    fold={t:1,u:s.sqrt(2)-2,v:2-2*a}
    assert s.simplify(reduce(K.subs(fold)))==0
    assert s.simplify(reduce(s.diff(K,v).subs(fold)))==0
    aq=((arb(13).sqrt()-1)/6).sqrt()
    am=mp.sqrt((mp.sqrt(13)-1)/6)
    ae=lambda expr:ball_function(expr,(a,))([aq])
    aa,bb,cc=map(ae,[A,B,C]);disc=bb*bb-4*aa*cc
    assert aa>0 and cc<0 and disc>0 and not ae(denominator).contains(0)
    us=(-bb-disc.sqrt())/(2*aa)
    vs=ball_function(vr.as_numer_denom()[0],(a,u))([aq,us])/ball_function(vr.as_numer_denom()[1],(a,u))([aq,us])
    assert us<0
    anchor=[us,vs]
    other_us=(-bb+disc.sqrt())/(2*aa)
    other_vs=ball_function(vr.as_numer_denom()[0],(a,u))([aq,other_us])/ball_function(vr.as_numer_denom()[1],(a,u))([aq,other_us])
    assert other_us>0
    target=[arb(2).sqrt()-2,2-2*aq]
    p4,(x,y,al,be)=factor_expression(4)
    spectator=reduce(p4.subs({x:-(2*a+t*v)**2,y:(a+(2-a)*t)**2,al:1+t*u/2,be:1+t*u/2}))
    other=reduce(f.subs({w:-(2*a+t*v),q:a+(2-a)*t,e:t*u}))
    eq=s.Matrix([K,s.diff(K,v)]);jac=eq.jacobian([u,v])
    variables=(u,v,t,a)
    def compiled(expr):
        # Combine each algebraic coefficient before interval evaluation in
        # u,v,t. Expanding Q as a fourth varying monomial would introduce
        # severe artificial dependency in near-cancelling coefficients.
        terms=[(pw,ball_function(co,(a,))([aq]))for pw,co in s.Poly(expr,u,v,t).terms()]
        def tree(items,index):
            if index<0:return sum((co for pw,co in items),arb(0))
            powers=sorted({pw[index]for pw,co in items},reverse=True)
            return [(power,tree([(pw,co)for pw,co in items if pw[index]==power],index-1))for power in powers]
        structure=tree(terms,2)
        def evaluate(coords):
            def walk(node,index):
                if index<0:return node
                value=arb(0);previous=node[0][0]
                for power,child in node:
                    value=value*coords[index]**(previous-power)+walk(child,index-1)
                    previous=power
                return value*coords[index]**previous
            return walk(structure,2)
        return evaluate
    def wrapped(expr,centered_t=False):
        fn=compiled(expr);derivative=compiled(s.diff(expr,t))if centered_t else None
        def evaluate(point):
            coords=[point[0],point[1],1-point[2]]
            if not centered_t:return fn(coords)
            midpoint=coords[2].mid()
            return fn(coords[:2]+[midpoint])+derivative(coords)*(coords[2]-midpoint)
        return evaluate
    funcs=[wrapped(expr,True)for expr in eq]
    jf=[[wrapped(jac[i,j])for j in range(2)]for i in range(2)]
    # The shared tile engine's k and r slots mean twice the end weight and
    # frequency here, respectively. They are not the outgoing chart's k,r.
    gates={name:wrapped(expr)for name,expr in {
        'K_u':s.diff(K,u),'K_vv':s.diff(K,v,2),'scaled_P4':spectator,
        'other_octic_diagnostic':other,'k':2+t*u,'r':2*a+t*v,'u':u}.items()}
    derivative_parts=[wrapped(expr)for expr in [s.diff(K,t),s.diff(K,u,v),s.diff(K,v,t),s.diff(other,t),s.diff(other,u),s.diff(other,v)]]
    def other_derivative(point):
        kt,kuv,kvt,ht,hu,hv=[fn(point)for fn in derivative_parts]
        up=-kt/gates['K_u'](point)
        vp=-(kvt+kuv*up)/gates['K_vv'](point)
        return ht+hu*up+hv*vp
    gates['other_path_derivative']=other_derivative
    nf=s.lambdify(variables,eq,'mpmath',cse=True)
    nj=s.lambdify(variables,jac,'mpmath',cse=True)
    functions=(funcs,jf,gates,lambda uu,vv,zz:nf(uu,vv,1-zz,am),lambda uu,vv,zz:nj(uu,vv,1-zz,am))
    local_bounds={name:fn(anchor+[arb(1)])for name,fn in gates.items()}
    assert not local_bounds['K_u'].contains(0) and not local_bounds['K_vv'].contains(0)
    data=dict(seed_field='3Q^4+Q^2-1=0; Q=sqrt((sqrt(13)-1)/6)',
              coordinates='t in [0,1]; q=Q+(2-Q)t; omega=2Q+t*v; epsilon=t*u; tile engine z=1-t',
              normalized_polynomial=str(K),normalization='K=27F/(256t^2)',
              seed_u_polynomial=str(seed_polynomial),seed_v_formula=str(vr),
              inherited_slope_polynomial=str(inherited),exact_slope_cross_pin=True,
              seed_anchor_dyadic=[pack_ball(x)for x in anchor],
              other_seed_anchor_dyadic=[pack_ball(other_us),pack_ball(other_vs)],
              fold_anchor_dyadic=[pack_ball(x)for x in target],
              seed_anchor=[x.str(25,more=True)for x in anchor],
              initial_nonzero_bounds={name:pack_ball(value)for name,value in local_bounds.items()},
              domain_slot_meanings={'k':'2+t*u, twice the end bond weight','r':'2Q+t*v, frequency omega'})
    return functions,anchor,target,data


def certify_tile(high,low,known,functions):
    """Incoming version: no exclusion of the other reflection sector."""
    funcs,jf,gates,nf,nj=functions;mid=(high+low)/2
    def solve(zz,guess):
        zm=mp.mpf(zz.numerator)/zz.denominator
        return mp.findroot(lambda a,b:tuple(nf(a,b,zm)),guess,
                           J=lambda a,b:nj(a,b,zm),tol=mp.mpf('1e-55'),maxsteps=35)
    guess=tuple(mp.mpf(str(float(x.mid())))for x in known);values=[]
    for zz in [high,mid,low]:
        value=solve(zz,guess);values.append(value);guess=tuple(value)
    center=[ab(values[1][i]).mid()for i in range(2)]
    spread=max(abs(value[i]-values[1][i])for value in values for i in range(2))
    radius=ab(spread*mp.mpf('1.6')+mp.mpf('1e-12'))
    zz=ab(mid)+arb(0,ab((high-low)/2))
    J0=acb_mat([[f(center+[ab(mid)])for f in row]for row in jf]).inv()
    Y=acb_mat([[acb(J0[i,j].real.mid())for j in range(2)]for i in range(2)])
    if Y.det().contains(0):raise ValueError('singular preconditioner')
    identity=acb_mat([[1,0],[0,1]])
    residual=Y*acb_mat([[f(center+[zz])]for f in funcs])
    radius=max(radius,2*max(abs(residual[i,0]).upper()for i in range(2)))
    uv=[c+arb(0,radius)for c in center]
    if not all(uv[i].contains(known[i])for i in range(2)):raise ValueError('incoming enclosure')
    box=uv+[zz];defect=identity-Y*acb_mat([[f(box)for f in row]for row in jf])
    rows=[sum((abs(defect[i,j]).upper()for j in range(2)),arb(0))for i in range(2)]
    images=[abs(residual[i,0]).upper()+radius*rows[i]for i in range(2)]
    if not(all(x<1 for x in rows)and all(x<radius for x in images)):raise ValueError('uniform contraction')
    bounds={name:fn(box)for name,fn in gates.items()}
    if any(bounds[name].contains(0)for name in ['K_u','K_vv','scaled_P4']):raise ValueError('selected-sector nonzero')
    if bounds['other_octic_diagnostic'].contains(0) and bounds['other_path_derivative'].contains(0):raise ValueError('cross-sector monotonicity')
    if not(bounds['k']>0 and bounds['k']<4 and bounds['r']>0 and bounds['u']<0):raise ValueError('domain')
    endpoint=uv
    fcenter=Y*acb_mat([[f(center+[ab(low)])]for f in funcs])
    de=identity-Y*acb_mat([[f(uv+[ab(low)])for f in row]for row in jf])
    for _ in range(8):
        delta=acb_mat([[endpoint[i]-center[i]]for i in range(2)])
        value=acb_mat([[center[i]]for i in range(2)])-fcenter+de*delta
        endpoint=[endpoint[i].intersection(value[i,0].real)for i in range(2)]
    record=dict(z_high=str(high),z_low=str(low),center=[x.str(20,more=True)for x in center],
                radius=str(radius),box_dyadic=[pack_ball(x)for x in uv],
                endpoint_dyadic=[pack_ball(x)for x in endpoint],center_dyadic=[pack_ball(x)for x in center],
                radius_dyadic=pack_ball(radius),
                preconditioner_dyadic=[[pack_ball(Y[i,j].real)for j in range(2)]for i in range(2)],
                row_bounds_dyadic=[pack_ball(x)for x in rows],image_bounds_dyadic=[pack_ball(x)for x in images],
                nonzero_and_domain_bounds_dyadic={name:pack_ball(value)for name,value in bounds.items()},
                other_start_dyadic=pack_ball(gates['other_octic_diagnostic'](known+[ab(high)])),
                other_end_dyadic=pack_ball(gates['other_octic_diagnostic'](endpoint+[ab(low)])),
                other_sector_zero_possible=bounds['other_octic_diagnostic'].contains(0),incoming_endpoint_contained=True)
    return record,endpoint


def audit(tiles,initial,target):
    if not tiles:return False
    known=initial;current=Fraction(1)
    for tile in tiles:
        if Fraction(tile['z_high'])!=current:return False
        low=Fraction(tile['z_low'])
        if not 0<=low<current:return False
        box=tile['box_dyadic']
        if not all(contains_dyadic(box[i],known[i])for i in range(2)):return False
        if not all(rational_bounds(x)[1]<1 for x in tile['row_bounds_dyadic']):return False
        if not all(rational_bounds(x)[1]<rational_bounds(tile['radius_dyadic'])[0]for x in tile['image_bounds_dyadic']):return False
        bounds={name:rational_bounds(value)for name,value in tile['nonzero_and_domain_bounds_dyadic'].items()}
        if not all(bounds[name][0]*bounds[name][1]>0 for name in ['K_u','K_vv','scaled_P4']):return False
        if not(bounds['k'][0]>0 and bounds['k'][1]<4 and bounds['r'][0]>0 and bounds['u'][1]<0):return False
        if bounds['other_octic_diagnostic'][0]<=0<=bounds['other_octic_diagnostic'][1]:
            if bounds['other_path_derivative'][0]*bounds['other_path_derivative'][1]<=0:return False
        known=tile['endpoint_dyadic'];current=low
    return current==0 and all(contains_dyadic(tiles[-1]['box_dyadic'][i],target[i])for i in range(2))


def close_crossings(data):
    """Certify all exceptional coincidences and identify the external roots."""
    tiles=data['tiles'];possible=[]
    for i,tile in enumerate(tiles):
        lo,hi=rational_bounds(tile['nonzero_and_domain_bounds_dyadic']['other_octic_diagnostic'])
        if lo<=0<=hi:possible.append(i)
    groups=[]
    for index in possible:
        if not groups or index!=groups[-1][-1]+1:groups.append([])
        groups[-1].append(index)
    def sign(record):
        lo,hi=rational_bounds(record)
        assert lo*hi>0
        return 1 if lo>0 else -1
    strips=[]
    for group in groups:
        signs={sign(tiles[i]['nonzero_and_domain_bounds_dyadic']['other_path_derivative'])for i in group}
        assert len(signs)==1
        first,last=tiles[group[0]],tiles[group[-1]]
        start_sign=sign(first['other_start_dyadic']);end_sign=sign(last['other_end_dyadic'])
        assert start_sign==-end_sign and end_sign==next(iter(signs))
        strips.append(dict(tile_indices=group,t_low=str(1-Fraction(first['z_high'])),
                           t_high=str(1-Fraction(last['z_low'])),derivative_sign=end_sign,
                           start_sign=start_sign,end_sign=end_sign,unique_crossing=True))
    crossing_path=ROOT/'simulations/results/route_b_n4_incoming_crossing.json'
    external=json.loads(crossing_path.read_text())
    def enclose(record):
        def dyadic(pair):return arb(int(pair[0]))*arb(2)**int(pair[1])
        # Reconstruction may enlarge the source ball; subsequent strict
        # containment proves inclusion of that larger, hence original, box.
        return dyadic(record['mid'])+arb(0,dyadic(record['rad']))
    aq=((arb(13).sqrt()-1)/6).sqrt();associations=[]
    for index,root in enumerate(external['roots']):
        assert root['success'] and sorted(root['jordan_block_sizes'])==[1,2]
        omega,q,e=[enclose(value)for value in root['box_dyadic']]
        t=(q-aq)/(2-aq);z=1-t;uv=[e/t,(omega-2*aq)/t]
        zlo,zhi=rational_bounds(pack_ball(z));uvd=[pack_ball(x)for x in uv]
        matches=[i for i,tile in enumerate(tiles)
                 if Fraction(tile['z_low'])<=zlo and zhi<=Fraction(tile['z_high'])
                 and all(contains_dyadic(tile['box_dyadic'][j],uvd[j])for j in range(2))]
        assert matches
        groupmatches=[j for j,group in enumerate(groups)if any(i in group for i in matches)]
        assert len(groupmatches)==1
        associations.append(dict(external_root_index=index,tile_indices=matches,strip_index=groupmatches[0],
                                 mapped_t_dyadic=pack_ball(t),mapped_uv_dyadic=uvd,
                                 exact_root_box_contained_in_uniqueness_tile=True))
    def bijective(records):
        return len(records)==len(strips) and sorted(row['strip_index']for row in records)==list(range(len(strips)))
    assert len(strips)==2 and bijective(associations)
    assert not bijective(associations[:-1])
    assert not bijective([associations[0],associations[0]])
    data.update(crossing_strips=strips,crossing_associations=associations,
                opposite_sector_crossing_count=2,crossing_completeness_certified=True,
                crossing_record_omission_control_rejected=True,duplicate_crossing_record_control_rejected=True,
                crossing_certificate_sha256=hashlib.sha256(crossing_path.read_bytes()).hexdigest(),
                crossing_completeness_proof='Every tile outside the listed strips excludes a zero of the other octic. Within each contiguous strip the derivative along the certified path has one strict sign; certified endpoint signs are opposite. Each strip therefore contains exactly one crossing. Both independent three-dimensional root boxes map wholly into path uniqueness tiles, identifying them as those crossings. Their full Jordan type is J2 direct sum J1.')
    return data


def main():
    ctx.prec=384;mp.mp.dps=75
    functions,known,target,data=setup()
    current=Fraction(1);step=Fraction(1,512);tiles=[];rejected=0
    while current>0:
        low=max(Fraction(0),current-step)
        try:record,endpoint=certify_tile(current,low,known,functions)
        except(ValueError,ZeroDivisionError) as error:
            step/=2;rejected+=1
            if step<Fraction(1,2**24):
                print('FAILED at t',float(1-current),'tiles',len(tiles),'cause',str(error),flush=True)
                print('known',known,flush=True)
                raise
            continue
        tiles.append(record);current=low;known=endpoint
        if len(tiles)%50==0:print('tiles',len(tiles),'t',float(1-current),'rejected',rejected,flush=True)
        step=min(2*step,Fraction(1,32))
    initial=data['seed_anchor_dyadic'];end=data['fold_anchor_dyadic']
    assert audit(tiles,initial,end)
    assert not audit(tiles[:len(tiles)//2]+tiles[len(tiles)//2+1:],initial,end)
    wrong=[pack_ball(target[0]+arb('0.1')),pack_ball(target[1])]
    assert not audit(tiles,initial,wrong)
    assert not audit(tiles,data['other_seed_anchor_dyadic'],end)
    summaries={}
    for name in tiles[0]['nonzero_and_domain_bounds_dyadic']:
        values=[rational_bounds(tile['nonzero_and_domain_bounds_dyadic'][name])for tile in tiles]
        lower=min(a for a,b in values);upper=max(b for a,b in values)
        summaries[name]=dict(lower_exact=str(lower),upper_exact=str(upper),lower_decimal=float(lower),upper_decimal=float(upper))
    data.update(scope=__doc__,precision_bits=ctx.prec,tile_count=len(tiles),rejected_attempts=rejected,
                tiles=tiles,global_gate_bounds=summaries,success=True,gap_control_rejected=True,
                wrong_fold_control_rejected=True,
                wrong_seed_branch_control_rejected=True,
                proof='At t=0 the exact quadratic root and eliminated v solve K=Kv=0. Uniform contraction tiles identify the same root via contained endpoint enclosures. The final uniqueness box contains the exact fold. For t>0, F_epsilon=(256/27)t K_u and F_omegaomega=(256/27)K_vv are nonzero. P4 is excluded, but the other octic is tracked without exclusion: cross-parity coincidences can raise full spectral multiplicity. The original t=0 eigenvalue remains the independently certified semisimple seed.',
                source_hashes=source_hashes())
    close_crossings(data)
    (ROOT/'simulations/results/route_b_n4_incoming_certificate.json').write_text(json.dumps(data,indent=2)+'\n')
    print('PASS',len(tiles),'incoming tiles, gap and wrong-fold controls rejected',flush=True)


if __name__=='__main__':
    if '--finalize-crossings' in sys.argv:
        ctx.prec=384;mp.mp.dps=75
        path=ROOT/'simulations/results/route_b_n4_incoming_certificate.json'
        data=json.loads(path.read_text())
        _,_,_,exact_data=setup()
        data.update(exact_data)
        assert audit(data['tiles'],data['seed_anchor_dyadic'],data['fold_anchor_dyadic'])
        assert not audit(data['tiles'],data['other_seed_anchor_dyadic'],data['fold_anchor_dyadic'])
        data['wrong_seed_branch_control_rejected']=True
        close_crossings(data)
        data['source_hashes']=source_hashes()
        path.write_text(json.dumps(data,indent=2)+'\n')
        print('PASS exact crossing completeness and association')
    else:main()
