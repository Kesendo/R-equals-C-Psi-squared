"""High-precision numerical EP paths on the exact N4 characteristic factor.

Root residuals are numerical evidence, not interval certificates.
"""
import json
from pathlib import Path
import mpmath as mp
from route_b_n4_range_polynomial import profile_evaluator


def follow(name,branch,direction,limit,functions=None):
    mp.mp.dps=55
    f,jac=functions or profile_evaluator(name,'mpmath')
    q0=mp.sqrt((mp.sqrt(13)-1)/6)
    slopes={'one':[-.520880800087,.100646300073],'even':[-.520880800087,.100646300073],'odd':[-30.148714776875,37.494111387670]}
    order=2 if name=='odd' else 1
    e=mp.mpf(direction)*mp.mpf('.0001')
    q=q0+mp.mpf(slopes[name][branch])*e**order
    state=mp.matrix([-4*q0*q0,q*q]);step=mp.mpf('.00001')
    rows=[];old=None
    tolerance=mp.eps**mp.mpf('.75')
    while len(rows)<3000:
        def fun(x,y):return tuple(f(x,y,e)[:2])
        def J(x,y):return jac(x,y,e)[:2,:2]
        try:
            new=mp.findroot(fun,tuple(state),J=J,tol=tolerance,maxsteps=30)
            if old is not None and mp.norm(new-state)>mp.mpf('.001'):
                raise ValueError('Predictor deviation')
        except (ValueError,ZeroDivisionError):
            if old is None:raise
            step/=2
            if step<mp.mpf('1e-10'):break
            e=old[0]+direction*step
            state=old[1]+direction*step*old[2]
            continue
        values=f(*new,e);der=jac(*new,e)
        tangent=mp.lu_solve(der[:2,:2],-der[:2,2])
        rows.append(dict(epsilon=str(e),x=str(new[0]),y=str(new[1]),q=str(mp.sqrt(new[1])),omega=str(mp.sqrt(-new[0])),residual=str(max(abs(v) for v in values[:2])),Pxx=str(values[2]),Py=str(values[3])))
        old=(e,new,tangent)
        if abs(e)>=limit:break
        step=min(step*mp.mpf('1.15'),mp.mpf('.002'))
        e+=direction*step;state=new+direction*step*tangent
    endpoints=[]
    if abs(old[0])<limit:
        for index,kind in [(2,'Pxx'),(3,'Py')]:
            try:
                def F(x,y,e):
                    vals=f(x,y,e);return vals[0],vals[1],vals[index]
                candidate=mp.findroot(F,(old[1][0],old[1][1],old[0]),tol=tolerance,maxsteps=60)
                endpoints.append(dict(kind=kind,coordinates=[str(v) for v in candidate],residual=str(max(abs(v) for v in F(*candidate)))))
            except (ValueError,ZeroDivisionError):pass
    return dict(name=name,branch=branch,direction=direction,limit=limit,rows=rows,candidates=endpoints)


if __name__=='__main__':
    out={}
    for name in ['odd']:
        funcs=profile_evaluator(name,'mpmath')
        for branch in [0,1]:
            key=f'{name}_{branch}'
            out[key]=follow(name,branch,1,.7,funcs)
            print(key,'last',out[key]['rows'][-1],'candidates',out[key]['candidates'],flush=True)
            (Path(__file__).parent/'results/route_b_n4_range_precise.json').write_text(json.dumps(out,indent=2)+'\n')
