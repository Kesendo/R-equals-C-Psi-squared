"""Pseudo-arclength continuation of the second real arm at the N4 equal-end turn.

Numerical path association, not a certified global connection.
"""
import json,hashlib
from pathlib import Path
import mpmath as mp
import sympy as s
from route_b_n4_range_even import functions
from route_b_n4_fold_infinity import derive


def cross(a,b):
    return mp.matrix([a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]])


def compact_track(last):
    data=derive();r,z,k,u,v=s.symbols('r z k u v')
    g=s.sympify(data['compactified_polynomial'])
    blown=s.Poly(s.cancel(g.subs({r:1+z*v,k:3+z*u})/(32*z**3)),u,v,z).as_expr()
    eq=s.Matrix([blown,s.diff(blown,v)])
    f=s.lambdify((u,v,z),eq,'mpmath',cse=True)
    J=s.lambdify((u,v,z),eq.jacobian([u,v]),'mpmath',cse=True)
    q=mp.mpf(last['q']);e=mp.mpf(last['epsilon']);omega=mp.mpf(last['omega'])
    zz=1/q;state=mp.matrix([q*((e+2)**2-3),omega-q]);rows=[]
    while zz>mp.mpf('1e-7'):
        state=mp.findroot(lambda a,b:tuple(f(a,b,zz)),tuple(state),J=lambda a,b:J(a,b,zz),tol=mp.mpf('1e-40'))
        rows.append(dict(z=str(zz),u=str(state[0]),v=str(state[1]),epsilon=str(mp.sqrt(3+zz*state[0])-2),q=str(1/zz),residual=str(max(abs(w) for w in f(*state,zz)))))
        zz*=mp.mpf('.8')
    distances=[]
    for germ in data['germs']:
        u0=mp.mpf(str(s.N(s.sympify(germ['u0']),55)));v0=mp.mpf(str(s.N(s.sympify(germ['v0']),55)))
        distances.append(dict(u0=germ['u0'],v0=germ['v0'],distance=str(mp.norm(state-mp.matrix([u0,v0])))))
    selected=min(distances,key=lambda a:mp.mpf(a['distance']))
    return dict(rows=rows,limiting_germ_comparison=distances,nearest_germ=selected,
                scope='Numerical continuation from final finite point to small z; exact germs are supplied separately.')


def main():
    mp.mp.dps=55;f,jac=functions()
    point=mp.matrix([2,2,mp.sqrt(2)-2])
    tangent=mp.matrix([1,1,0])/mp.sqrt(2)
    step=mp.mpf('.005');rows=[];stop='step_limit';crossing=None
    for iteration in range(3000):
        prediction=point+step*tangent
        scale=(1+mp.norm(prediction))**8
        def F(w,q,e):
            vals=f(w,q,e)
            return vals[0]/scale,vals[1]/scale,mp.fdot(tangent,mp.matrix([w,q,e])-prediction)
        def J(w,q,e):
            a=jac(w,q,e)
            return mp.matrix([[a[i,j]/scale for j in range(3)] for i in range(2)]+[list(tangent)])
        try:
            new=mp.findroot(F,tuple(prediction),J=J,tol=mp.mpf('1e-40'),maxsteps=35)
            if mp.norm(new-prediction)>step/3:raise ValueError('predictor correction')
        except (ValueError,ZeroDivisionError):
            step/=2
            if step<mp.mpf('1e-9'):stop='small_step';break
            continue
        d=jac(*new);nt=cross(list(d[0,:]),list(d[1,:]));nt/=mp.norm(nt)
        if mp.fdot(nt,tangent)<0:nt=-nt
        values=f(*new)
        rows.append(dict(omega=str(new[0]),q=str(new[1]),epsilon=str(new[2]),
                         tangent=[str(v) for v in nt],residual=str(max(abs(v) for v in F(*new))),Fww=str(values[2]),Fq=str(values[3])))
        if iteration%50==0:print(iteration,[float(v) for v in new],flush=True)
        if point[2]<0<=new[2]:
            root=mp.findroot(lambda w,q:tuple(f(w,q,0)[:2]),tuple(new[:2]),tol=mp.mpf('1e-40'))
            crossing=dict(omega=str(root[0]),q=str(root[1]),epsilon='0');stop='uniform_crossing';break
        point=new;tangent=nt;step=min(step*mp.mpf('1.12'),mp.mpf('.04'))
        if mp.norm(point)>200:stop='norm_limit';break
    compact=compact_track(rows[-1]) if stop=='step_limit' else None
    out=dict(scope=__doc__,rows=rows,stop=stop,uniform_crossing=crossing,compactification=compact,
             source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).parent/'route_b_n4_range_polynomial.py',Path(__file__).parent/'results/route_b_n4_range_polynomial.json',Path(__file__).parent/'route_b_n4_fold_infinity.py']})
    (Path(__file__).parent/'results/route_b_n4_fold_path.json').write_text(json.dumps(out,indent=2)+'\n')
    print(stop,crossing,'last',rows[-1],flush=True)
    if compact:print('compact limiting germ',compact['nearest_germ'],flush=True)


if __name__=='__main__':main()
