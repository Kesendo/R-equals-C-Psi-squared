"""High-precision continuation of equal-end N4 EPs in the selected octic.
Numerical paths and endpoint candidates, not certified interval bounds.
"""
import hashlib,json,sys
from pathlib import Path
import mpmath as mp
import sympy as s
from route_b_n4_range_polynomial import even_frequency_factor

def functions():
    f,(w,q,e)=even_frequency_factor()
    v=s.Matrix([f,s.diff(f,w),s.diff(f,w,2),s.diff(f,q)])
    return (s.lambdify((w,q,e),v,'mpmath',cse=True),
            s.lambdify((w,q,e),v.jacobian([w,q,e]),'mpmath',cse=True))

def follow(branch,direction,limit,funcs):
    mp.mp.dps=55;f,jac=funcs
    q0=mp.sqrt((mp.sqrt(13)-1)/6)
    e=direction*mp.mpf('.0001');slopes=[mp.mpf('-.520880800087'),mp.mpf('.100646300073')]
    state=mp.matrix([2*q0,q0+slopes[branch]*e])
    step=mp.mpf('.00001');rows=[];old=None;tol=mp.mpf('1e-40')
    while len(rows)<2500:
        def F(w,q):return tuple(f(w,q,e)[:2])
        def J(w,q):return jac(w,q,e)[:2,:2]
        try:
            new=mp.findroot(F,tuple(state),J=J,tol=tol,maxsteps=40)
            if old is not None and mp.norm(new-state)>mp.mpf('.001'):
                raise ValueError('predictor deviation')
        except (ValueError,ZeroDivisionError):
            if old is None:raise
            step/=2
            if step<mp.mpf('1e-10'):break
            e=old[0]+direction*step;state=old[1]+direction*step*old[2]
            continue
        values=f(*new,e);der=jac(*new,e)
        tangent=mp.lu_solve(der[:2,:2],-der[:2,2])
        rows.append(dict(epsilon=str(e),omega=str(new[0]),q=str(new[1]),
                         residual=str(max(abs(v) for v in values[:2])),Fww=str(values[2]),Fq=str(values[3])))
        old=(e,new,tangent)
        if abs(e)>=limit:break
        step=min(step*mp.mpf('1.2'),mp.mpf('.005'))
        step=min(step,mp.mpf(str(limit))-abs(e))
        e+=direction*step;state=new+direction*step*tangent
    candidates=[]
    if abs(old[0])<limit:
        for index,label in [(2,'Fww'),(3,'Fq')]:
            try:
                def endpoint(w,q,e):
                    v=f(w,q,e);return v[0],v[1],v[index]
                point=mp.findroot(endpoint,(*old[1],old[0]),tol=tol,maxsteps=75)
                candidates.append(dict(kind=label,coordinates=[str(v) for v in point],
                                       near_terminal=bool(mp.norm(point-mp.matrix([*old[1],old[0]]))<mp.mpf('.001')),
                                       residual=str(max(abs(v) for v in endpoint(*point)))))
            except (ValueError,ZeroDivisionError):pass
    return dict(branch=branch,direction=direction,limit=limit,rows=rows,candidates=candidates)

def main():
    funcs=functions();out=dict(scope=__doc__,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),paths={})
    f,(w,q,e)=even_frequency_factor();point={w:2,q:2,e:s.sqrt(2)-2}
    exact={name:str(s.simplify(value.subs(point))) for name,value in
           [('F',f),('Fw',s.diff(f,w)),('Fq',s.diff(f,q)),
            ('Fww',s.diff(f,w,2)),('Fe',s.diff(f,e))]}
    assert [exact[k] for k in ['F','Fw','Fq']]==['0','0','0']
    out['exact_projection_fold']=dict(omega='2',q='2',epsilon='sqrt(2)-2',derivatives=exact,
        boundary='Exact polynomial point; association with continued branch is numerical.')
    target=Path(__file__).parent/'results/route_b_n4_range_even.json'
    for direction in [1,-1]:
        for branch in [0,1]:
            result=follow(branch,direction,1 if direction==1 else 1.99,funcs)
            out['paths'][f'{branch}_{direction}']=result
            print('branch',branch,'direction',direction,'last',result['rows'][-1],'candidates',result['candidates'],flush=True)
            target.write_text(json.dumps(out,indent=2)+'\n')

def enrich_algebra():
    """Recompute exact supplemental identities without repeating stored paths."""
    target=Path(__file__).parent/'results/route_b_n4_range_even.json'
    out=json.loads(target.read_text());f,(w,q,e)=even_frequency_factor()
    point={w:2,q:2,e:s.sqrt(2)-2}
    fwq=s.simplify(s.diff(f,w,q).subs(point));fqq=s.simplify(s.diff(f,q,2).subs(point))
    fww=s.simplify(s.diff(f,w,2).subs(point))
    out['exact_projection_fold']['derivatives'].update(Fwq=str(fwq),Fqq=str(fqq),
        reduced_Fqq=str(s.simplify(fqq-fwq**2/fww)))
    resultant=s.factor(s.resultant(f.subs(w,0),s.diff(f,w).subs(w,0),q))
    vals=s.lambdify((q,e),(f.subs(w,0),s.diff(f,w).subs(w,0)),'mpmath',cse=True)
    mp.mp.dps=55;crossings=[]
    for guess in [('.512','-1.326'),('1.047','-.593')]:
        root=mp.findroot(vals,tuple(map(mp.mpf,guess)),tol=mp.mpf('1e-45'))
        crossings.append(dict(q=str(root[0]),epsilon=str(root[1]),
            on_tracked_upper_negative_branch=guess[0]=='.512',
            residual=str(max(abs(v) for v in vals(*root)))))
    out['zero_frequency']=dict(exact_resultant=str(resultant),numerical_roots=crossings,
        meaning='Positive/negative frequency factors meet; not necessarily a selected-octic branch endpoint.')
    out['decoupled_end_factor']=str(s.factor(f.subs(e,-2)))
    from route_b_other_n_unfolding import parts,independent_hop
    d,*_=parts(4);hop=independent_hop(4,[0,1,0])
    L=s.Matrix(24,24,lambda r,c:int(d[r,c].real)+s.I*s.Rational(1,2)*int(hop[r,c].imag))
    M=L+4*s.eye(24);z=s.Symbol('lambda')
    nullities=[24-(M**power).rank() for power in [1,2,3]]
    assert nullities==[2,4,4]
    out['decoupled_full_block']=dict(q='1/2',epsilon='-2',characteristic=str(s.factor(L.charpoly(z).as_expr())),
        nullities_M_M2_M3=nullities,character='At lambda=-4: two Jordan blocks of size2, not EP4.')
    out['algebra_enrichment_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    target.write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':
    if '--algebra-only' not in sys.argv:main()
    enrich_algebra()
