"""F89d palindrome and F131 spatial-reflection action on F163 response objects.
N=4,5,6 selected seeds; XY, gamma=1, q=qCSharp. Real-coefficient end profiles.
Exact dyadic matrix identities plus numerical independent partner-plane readings.
Run: python simulations/route_b_mirror_response.py
"""
import hashlib,json
from pathlib import Path
import numpy as np
import scipy.linalg as la
from framework.weight_coherence_block import weight_block_build,weight_block_configs
import route_b_n6_end_profile_probe as ep
from route_b_n6_next_coefficient import graph_coefficients
ROOT=Path(__file__).resolve().parents[1]
def z(x):return complex(x['real'],x['imag'])
def setup(n,wb):
    basis=[(a,b) for a in weight_block_configs(n,1) for b in weight_block_configs(n,wb)]
    index={x:i for i,x in enumerate(basis)}
    rev=lambda x:int(f'{x:0{n}b}'[::-1],2)
    R=[index[(rev(a),rev(b))] for a,b in basis]
    d=weight_block_build(n,1,wb,0); c=weight_block_build(n,1,wb,1,gamma=0)
    left=np.zeros_like(c)
    for i,(a,b) in enumerate(basis):
        for j,(aa,bb) in enumerate(basis):
            if (a==aa and b^bb==3) or (b==bb and a^aa==3):left[i,j]=c[i,j]
    return basis,d,c,left,left[np.ix_(R,R)],R

def reading(d,c,V,q,lam):
    L=d+q*c; gap=np.sort(abs(la.eigvals(L)-lam))[2]
    U,K=ep.base.plane(L,lam,gap/2);W=la.solve(U.T@U,U.T)
    P=U@W;Q=np.eye(len(d))-P
    S=Q@la.solve(lam*np.eye(len(d))-L+P,Q)
    A=W@c@U; B1=W@(q*V)@U;B2=W@(q*V)@S@(q*V)@U
    return dict(P=P,S=S,U=U,W=W,A=A,B1=B1,B2=B2,gap=float(gap),scalar_residual=float(la.norm(K-lam*np.eye(2))))

def main():
    sources=[ROOT/'simulations/results/route_b_other_n_unfolding.json',ROOT/'simulations/results/route_b_n6_end_profile_probe.json']
    earlier=json.loads(sources[0].read_text());six=json.loads(sources[1].read_text())
    seeds=[next(r for r in earlier['loci'] if r['id']==key) for key in ('N4-real-lambda-plus','N5-O-A2-W-006-Q-plus')]
    seeds.append(dict(n=6,id=six['seed_id'],q=six['q0'],lambda0=six['lambda0']))
    out=dict(scope=__doc__,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'simulations/route_b_n6_end_profile_probe.py',ROOT/'simulations/route_b_n6_next_coefficient.py',ROOT/'simulations/route_b_n6_leakage_probe.py',ROOT/'simulations/framework/weight_coherence_block.py']},sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},readings=[])
    for seed in seeds:
        n=seed['n'];bs,d,c,left,right,R=setup(n,2);bt,dp,cp,lp,rp,Rp=setup(n,n-2)
        index={x:i for i,x in enumerate(bt)};mask=(1<<n)-1
        perm=[index[(a,b^mask)] for a,b in bs]
        assert len(set(perm))==len(bs)==len(bt)
        assert all(Rp[perm[i]]==perm[R[i]] for i in range(len(bs)))
        # Target matrix pulled back to source ordering.
        pull=lambda M:M[np.ix_(perm,perm)]
        rung_pairs=sorted(set(((a^b).bit_count(),(a^(b^mask)).bit_count()) for a,b in bs))
        assert all(k+kp==n for k,kp in rung_pairs)
        exact_checks=[]
        for name,(a,b) in ep.DIRECTIONS.items():
            for q,e in ((.75,.125),(.75+.25j,.125+.0625j)):
                L=d+q*(c+e*(a*left+b*right))
                partner=dp+q.conjugate()*(cp+e.conjugate()*(a*lp+b*rp))
                assert np.array_equal(pull(partner),-L.conjugate()-2*n*np.eye(len(d)))
                reflected=d+q*(c+e*(b*left+a*right))
                assert np.array_equal(L[np.ix_(R,R)],reflected)
                exact_checks.append(dict(profile=name,q=ep.base.encode(q),epsilon=ep.base.encode(e),fold_residual=0,spatial_residual=0))
        qtest=.75+.25j;etest=.125+.0625j
        wrong=dp+qtest.conjugate()*(cp+etest*lp)
        expected=-(d+qtest*(c+etest*left)).conjugate()-2*n*np.eye(len(d))
        control=float(la.norm(pull(wrong)-expected));assert control>0
        q=z(seed['q']);lam=z(seed['lambda0']);qp=q.conjugate();lamp=-lam.conjugate()-2*n
        record=dict(n=n,seed_id=seed['id'],source_block=[1,2],partner_block=[1,n-2],dimension=len(d),disagreement_rung_pairs=rung_pairs,q=seed['q'],lambda0=seed['lambda0'],partner_q=ep.base.encode(qp),partner_lambda=ep.base.encode(lamp),exact_checks=exact_checks,unconjugated_epsilon_control=control,profiles={})
        for name,(a,b) in ep.DIRECTIONS.items():
            V=a*left+b*right;Vp=a*lp+b*rp
            x=reading(d,c,V,q,lam);y=reading(dp,cp,Vp,qp,lamp)
            # These full-space products avoid arbitrary basis alignment.
            residuals={key:float(la.norm(pull(y[key])-sign*x[key].conjugate())) for key,sign in [('P',1),('S',-1)]}
            for tag in ('B1','B2'):
                residuals[tag]=float(la.norm(pull(y['U']@y[tag]@y['W'])+(x['U']@x[tag]@x['W']).conjugate()))
            m=2 if name=='odd' else 1
            cs=ep.leading_slopes(x['A'],x['B2'] if m==2 else x['B1']);ct=ep.leading_slopes(y['A'],y['B2'] if m==2 else y['B1'])
            match=min(((0,1),(1,0)),key=lambda p:max(abs(ct[p[j]]-cs[j].conjugate()) for j in range(2)))
            ds=[graph_coefficients(d,c,V,q,lam,m,s)[0] for s in cs]
            dt=[graph_coefficients(dp,cp,Vp,qp,lamp,m,t)[0] for t in ct]
            record['profiles'][name]=dict(transport_residuals=residuals,source_gap=x['gap'],partner_gap=y['gap'],source_resolvent_norm=float(la.norm(x['S'],2)),partner_resolvent_norm=float(la.norm(y['S'],2)),coefficient_conjugation_error=float(max(abs(ct[match[j]]-cs[j].conjugate()) for j in range(2))),next_coefficient_conjugation_error=float(max(abs(dt[match[j]]-ds[j].conjugate()) for j in range(2))),next_coefficient_relative_error=float(max(abs(dt[match[j]]-ds[j].conjugate())/max(1,abs(ds[j])) for j in range(2))))
        out['readings'].append(record)
        print(n,'partner',lamp,'worst normalized next error',max(x['next_coefficient_relative_error'] for x in record['profiles'].values()),flush=True)
    (ROOT/'simulations/results/route_b_mirror_response.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
