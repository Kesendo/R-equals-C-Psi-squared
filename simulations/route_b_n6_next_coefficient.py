"""F163 next EP-location coefficient from a fixed-plane graph recurrence.
Numerical evaluation of an algebraic formula, not an interval remainder bound.
Run: python simulations/route_b_n6_next_coefficient.py
"""
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy.linalg as la
import route_b_n6_end_profile_probe as ep
ROOT=Path(__file__).resolve().parents[1]
def z(x): return complex(x['real'],x['imag'])
def polar(a,b): return 4*np.trace(a@b)-2*np.trace(a)*np.trace(b)

def graph_coefficients(d,c,v,q0,lam0,order,slope):
    l0=d+q0*c
    isolation=np.sort(abs(la.eigvals(l0)-lam0))[2]
    u,_=ep.base.plane(l0,lam0,isolation/2)
    w=la.solve(u.T@u,u.T)
    p=u@w; Q=np.eye(len(d))-p
    S=Q@la.solve(lam0*np.eye(len(d))-l0+p,Q)
    nmax=2*order
    E=[np.zeros_like(c) for _ in range(nmax+1)]
    E[1]+=q0*v
    E[order]+=slope*c
    if order+1<=nmax: E[order+1]+=slope*v
    X=[u]+[np.zeros_like(u) for _ in range(nmax)]
    K=[np.zeros((2,2),complex) for _ in range(nmax+1)]
    for n in range(1,nmax+1):
        forcing=sum((E[j]@X[n-j] for j in range(1,n+1)),np.zeros_like(u))
        K[n]=w@forcing
        feedback=sum((X[j]@K[n-j] for j in range(1,n)),np.zeros_like(u))
        X[n]=S@(Q@forcing-feedback)
    A=w@c@u
    derivative=polar(K[order],A)
    correction=-polar(K[order],K[2*order])/derivative
    return correction,dict(implicit_derivative=ep.base.encode(derivative),leading_discriminant=ep.base.encode(ep.base.discriminant(K[order])),odd_effective_norms=[float(la.norm(K[i])) for i in range(1,nmax+1,2)],graph_norms=[float(la.norm(x)) for x in X[1:]])

def main():
    original_path=ROOT/'simulations/results/route_b_n6_end_profile_probe.json'
    check_path=ROOT/'simulations/results/route_b_n6_validity_probe.json'
    original=json.loads(original_path.read_text()); measured=json.loads(check_path.read_text())
    q0,lam0=z(original['q0']),z(original['lambda0'])
    d,c,left=ep.base.parts(); right,_=ep.reflected_end(left)
    result=dict(scope=__doc__,seed_id=original['seed_id'],script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),source_sha256=hashlib.sha256(original_path.read_bytes()).hexdigest(),validation_sha256=hashlib.sha256(check_path.read_bytes()).hexdigest(),profiles={})
    for name,(a,b) in ep.DIRECTIONS.items():
        old=original['profiles'][name]; order=old['order']
        slopes=np.array(old['leading_slopes']['real'])+1j*np.array(old['leading_slopes']['imag'])
        entries=[]
        for branch,slope in enumerate(slopes):
            correction,diagnostics=graph_coefficients(d,c,a*left+b*right,q0,lam0,order,slope)
            rows=[]
            for row in measured['profiles'][name]['rows']:
                if row['branch']!=branch: continue
                e=row['epsilon']; actual=z(row['q'])-q0
                remainder=actual-slope*e**order-correction*e**(2*order)
                rows.append(dict(epsilon=e,leading_error=row['relative_displacement_error'],corrected_error=float(abs(remainder)/abs(actual)),scaled_remainder=ep.base.encode(remainder/e**(3*order)),observed_next_coefficient=ep.base.encode((actual-slope*e**order)/e**(2*order))))
            entries.append(dict(branch=branch,slope=ep.base.encode(slope),next_coefficient=ep.base.encode(correction),diagnostics=diagnostics,rows=rows))
            print(name,branch,'next',correction,'errors',[(r['epsilon'],r['corrected_error']) for r in rows if r['epsilon']>=0.01],flush=True)
        result['profiles'][name]=dict(order=order,branches=entries)
    differences=[z(a['next_coefficient'])-z(b['next_coefficient']) for a,b in zip(result['profiles']['one']['branches'],result['profiles']['even']['branches'])]
    odd=[z(x['slope']) for x in result['profiles']['odd']['branches']]
    set_residual=min(max(abs(differences[i]-odd[permutation[i]]) for i in range(2)) for permutation in ((0,1),(1,0)))
    # Exact integer-valued counterexample: reflection-like cancellation by itself
    # cannot replace common bilinear symmetry. The predicted ratios are both zero,
    # but zero is not a root of this antisymmetric B2 discriminant.
    A=np.diag([1,-1]); Bplus=np.array([[0,1],[1,0]]); B2=np.array([[0,1],[-1,0]])
    control=[-polar(c*A+Bplus,B2)/polar(c*A+Bplus,A) for c in (1j,-1j)]
    result['profile_identity']=dict(differences=ep.base.encode(differences),unordered_residual=float(set_residual),nonsymmetric_control_ratios=ep.base.encode(control),nonsymmetric_control_discriminants=ep.base.encode([ep.base.discriminant(x*A+B2) for x in control]))
    (ROOT/'simulations/results/route_b_n6_next_coefficient.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__': main()
