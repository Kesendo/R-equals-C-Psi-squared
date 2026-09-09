"""N4/N5 full (1,2) block tests of the local end-bond unfolding structure.
Gamma=1, XY delta=0, q=qCSharp, hopping2q. Numerical instances, not new
all-locus rank or nonzero certificates. All58 N5 inventory q loci scanned.
Run: python simulations/route_b_other_n_unfolding.py
"""
import hashlib,json
from pathlib import Path
import numpy as np
import scipy.linalg as la
from scipy.optimize import root
from framework.weight_coherence_block import weight_block_build,weight_block_configs
import route_b_n6_end_profile_probe as ep
from route_b_n6_next_coefficient import graph_coefficients,polar
ROOT=Path(__file__).resolve().parents[1]
def z(x):return complex(float(x['real']),float(x['imag']))
def parts(n):
    d=weight_block_build(n,1,2,0); c=weight_block_build(n,1,2,1,gamma=0)
    basis=[(a,b) for a in weight_block_configs(n,1) for b in weight_block_configs(n,2)]
    index={p:i for i,p in enumerate(basis)}
    rev=lambda a:int(f'{a:0{n}b}'[::-1],2)
    reflection=[index[(rev(a),rev(b))] for a,b in basis]
    left=np.zeros_like(c)
    for i,(a,b) in enumerate(basis):
        for j,(aa,bb) in enumerate(basis):
            if (a==aa and b^bb==3) or (b==bb and a^aa==3):left[i,j]=c[i,j]
    right=left[np.ix_(reflection,reflection)]
    return d,c,left,right,reflection

def independent_hop(n,weights):
    h=np.zeros((2**n,2**n))
    for s,w in enumerate(weights):
        for a in range(2**n):
            if ((a>>s)&1)!=((a>>(s+1))&1):h[a^(3<<s),a]=2*w
    a,b=weight_block_configs(n,1),weight_block_configs(n,2)
    return -1j*(np.kron(h[np.ix_(a,a)],np.eye(len(b)))-np.kron(np.eye(len(a)),h[np.ix_(b,b)]))

def coefficients(n,q,lam):
    d,c,left,right,reflection=parts(n); dim=len(d)
    l=d+q*c; distances=np.sort(abs(la.eigvals(l)-lam)); gap=distances[2]
    u,k=ep.base.plane(l,lam,gap/2)
    w=la.solve(u.T@u,u.T); p=u@w; Q=np.eye(dim)-p
    S=Q@la.solve(lam*np.eye(dim)-l+p,Q)
    A=w@c@u
    result=dict(n=n,dimension=dim,q=ep.base.encode(q),lambda0=ep.base.encode(lam),third_distance=float(gap),
        scalar_plane_residual=float(la.norm(k-lam*np.eye(2))),reflection_on_plane=ep.base.encode(w@u[reflection]),profiles={})
    for name,(a,b) in ep.DIRECTIONS.items():
        V=a*left+b*right; B1=w@(q*V)@u; B2=w@(q*V)@S@(q*V)@u
        m=2 if name=='odd' else 1; B=B2 if m==2 else B1
        alpha=ep.base.discriminant(A); beta=polar(A,B); gamma=ep.base.discriminant(B)
        omega=beta**2-4*alpha*gamma
        slopes=ep.leading_slopes(A,B)
        ds=[graph_coefficients(d,c,V,q,lam,m,slope)[0] for slope in slopes]
        independent=independent_hop(n,[1+.125*a]+[1]*(n-3)+[1+.125*b])
        assert np.array_equal(c+.125*V,independent)
        result['profiles'][name]=dict(order=m,alpha=ep.base.encode(alpha),omega=ep.base.encode(omega),first_order_norm=float(la.norm(B1)),slopes=ep.base.encode(slopes),next_coefficients=ep.base.encode(ds))
    one=result['profiles']['one']; even=result['profiles']['even']; odd=result['profiles']['odd']
    unpack=lambda x:np.array(x['real'])+1j*np.array(x['imag'])
    cs1,cs2=unpack(one['slopes']),unpack(even['slopes'])
    perm=min(((0,1),(1,0)),key=lambda pp:max(abs(cs1[j]-cs2[pp[j]]) for j in range(2)))
    diffs=unpack(one['next_coefficients'])-unpack(even['next_coefficients'])[list(perm)]
    os=unpack(odd['slopes'])
    mismatch=min(max(abs(diffs[j]-os[pp[j]]) for j in range(2)) for pp in ((0,1),(1,0)))
    result['identity_absolute_residual']=float(mismatch)
    result['identity_relative_residual']=float(mismatch/max(1,max(abs(os))))
    return result

def follow(n,q0,lam0,record):
    d,c,left,right,_=parts(n); radius=record['third_distance']*.4
    result={}
    for name,(a,b) in ep.DIRECTIONS.items():
        p=record['profiles'][name]; m=p['order']; rows=[]
        for branch in range(2):
            slope=complex(p['slopes']['real'][branch],p['slopes']['imag'][branch])
            previous=slope
            for e in (((.001,.002,.003) if m==2 else (.0001,.0003,.001)) if n==4 else (.001,.003,.01)):
                eta=e**m; hop=c+e*(a*left+b*right)
                def f(x):
                    _,k=ep.base.plane(d+(q0+eta*complex(*x))*hop,lam0,radius)
                    val=ep.base.discriminant(k)/eta**2
                    return np.array([val.real,val.imag])
                def jac(x):
                    h=.001
                    return np.column_stack([(f(x+h*np.eye(2)[j])-f(x-h*np.eye(2)[j]))/(2*h) for j in range(2)])
                sol=root(f,[previous.real,previous.imag],jac=jac,tol=1e-7)
                if not sol.success:raise ValueError((n,name,branch,e,str(sol.message)))
                previous=complex(*sol.x); q=q0+eta*previous
                _,k=ep.base.plane(d+q*hop,lam0,radius)
                lam=np.trace(k)/2
                rows.append(dict(epsilon=e,branch=branch,q=ep.base.encode(q),lambda0=ep.base.encode(lam),scaled_discriminant_residual=float(la.norm(f(sol.x))),traceless_norm=float(la.norm(k-lam*np.eye(2))),third_distance=float(np.sort(abs(la.eigvals(d+q*hop)-lam))[2]),q_displacement_error=float(abs(q-q0-eta*slope)/abs(q-q0)),smallest_singular_values=la.svdvals(d+q*hop-lam*np.eye(len(d)))[-3:].tolist()))
        for e in sorted(set(row['epsilon'] for row in rows)):
            pair=[row for row in rows if row['epsilon']==e]
            separation=abs(z(pair[0]['q'])-z(pair[1]['q']))
            for row in pair:row['branch_separation']=float(separation)
        result[name]=rows
    return result

def residue_control(record):
    n=record['n']; q=z(record['q']); lam=z(record['lambda0'])
    d,c,left,right,reflection=parts(n); l=d+q*c; V=(left-right)/2
    u,_=ep.base.plane(l,lam,record['third_distance']/2)
    w=la.solve(u.T@u,u.T); p=u@w; Q=np.eye(len(d))-p
    S=Q@la.solve(lam*np.eye(len(d))-l+p,Q)
    vals,vecs=la.eig(l); order=np.argsort(abs(vals-lam))
    closest=abs(vals[order[2]]-lam)
    # Numerical cluster selection includes both conjugate members at N5.
    chosen=[j for j in order[2:] if abs(abs(vals[j]-lam)-closest)<1e-10*max(1,closest)]
    pole=np.zeros_like(l); modes=[]
    for j in chosen:
        v=vecs[:,j]; pj=np.outer(v,v)/(v.T@v)
        pole+=pj/(lam-vals[j])
        modes.append(dict(eigenvalue=ep.base.encode(vals[j]),reflection_character=ep.base.encode((v.conj()@v[reflection])/(v.conj()@v)),projector_norm=float(la.norm(pj,2)),eigenprojector_residual=float(la.norm(l@pj-vals[j]*pj)),idempotency_residual=float(la.norm(pj@pj-pj)),seed_overlap=float(la.norm(p@pj))))
    A=w@c@u
    responses={}
    for name,operator in [('full',S),('nearest_only',pole),('nearest_removed',S-pole)]:
        B=w@(q*V)@operator@(q*V)@u
        responses[name]=dict(effective_norm=float(la.norm(B)),odd_slopes=ep.base.encode(ep.leading_slopes(A,B)))
    return dict(resolvent_norm=float(la.norm(S,2)),nearest_modes=modes,responses=responses,scope='Numerical decomposition of reduced resolvent; removing a residue is a diagnostic, not a modified physical Hamiltonian.')

def main():
    source=ROOT/'simulations/results/route_b_a2_n5.json'; data=json.loads(source.read_text())
    rows=[]
    q=np.sqrt((-1+np.sqrt(13))/6)
    for sign in (1,-1):
        r=coefficients(4,q,-4+sign*2j*q);r['id']='N4-real-lambda-'+('plus' if sign==1 else 'minus');rows.append(r)
    for sector in data['sectors']:
        for a in sector['a2Roots']:
            for loc in a['qLoci']:
                r=coefficients(5,z(loc['qPhysicalCSharpSeed']),z(a['lambdaSeed']));r['id']=loc['id'];rows.append(r)
    selected=[rows[0],next(r for r in rows if r['id']=='N5-O-A2-W-006-Q-plus')]
    output=dict(scope=__doc__,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'simulations/route_b_n6_next_coefficient.py',ROOT/'simulations/route_b_n6_end_profile_probe.py',ROOT/'simulations/route_b_n6_leakage_probe.py',ROOT/'simulations/framework/weight_coherence_block.py']},loci=rows,continuations={})
    # Persist coefficient survey even if the finer continuation fails.
    path=ROOT/'simulations/results/route_b_other_n_unfolding.json'
    path.write_text(json.dumps(output,indent=2)+'\n')
    for r in selected:
        print('selected',r['id'],'gap',r['third_distance'],flush=True)
        r['residue_control']=residue_control(r)
        output['continuations'][r['id']]=follow(r['n'],z(r['q']),z(r['lambda0']),r)
    path.write_text(json.dumps(output,indent=2)+'\n')
    print('count',len(rows),'worst relative identity',max(r['identity_relative_residual'] for r in rows),flush=True)
if __name__=='__main__':main()
