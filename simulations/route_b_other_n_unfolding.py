"""N4/N5 full (1,2) block tests of the local end-bond unfolding structure.
Gamma=1, XY delta=0, q=qCSharp, hopping2q. Numerical instances, not new
all-locus rank or nonzero certificates. All58 N5 inventory q loci scanned.
Run: python simulations/route_b_other_n_unfolding.py
"""
import json
import os
import re
from pathlib import Path
for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'

import numpy as np
import scipy.linalg as la
from scipy.optimize import root
from framework.weight_coherence_block import weight_block_build,weight_block_configs
import route_b_n6_end_profile_probe as ep
from route_b_n6_next_coefficient import graph_coefficients,polar
from route_b_artifact_provenance import artifact_provenance, verify_artifact_provenance
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'simulations/results/route_b_a2_n5.json'
PROFILE_ORDERS={'one':1,'even':1,'odd':2}
SELECTED=('N4-real-lambda-plus','N5-O-A2-W-006-Q-plus')
SCALE_FACTORS=(2.**-12,1.,2.**12)
EPS=np.finfo(float).eps
def _finite_number(value):
    """Accept JSON numbers or their canonical decimal spelling, never NaN/Inf."""
    if isinstance(value,str):
        if re.fullmatch(r'-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?',value) is None:
            raise ValueError('Expected canonical finite numeric component')
    elif isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,float,np.integer,np.floating)):
        raise ValueError('Expected finite numeric component')
    number=float(value)
    if not np.isfinite(number):
        raise ValueError('Nonfinite numeric component')
    return number


def z(x):
    return complex(_decode_complex(x,()))
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

def coefficients(n,q,lam,*,generator_scale=1.):
    d,c,left,right,reflection=parts(n); dim=len(d)
    d,c,left,right=(generator_scale*x for x in (d,c,left,right))
    lam=generator_scale*lam
    l=d+q*c; distances=np.sort(abs(la.eigvals(l)-lam)); gap=distances[2]
    u,k=ep.base.plane(l,lam,gap/2)
    w=la.solve(u.T@u,u.T); p=u@w; Q=np.eye(dim)-p
    S=Q@la.solve(lam*np.eye(dim)-l+p,Q)
    A=w@c@u
    error=dict(generator_norm=float(la.norm(l,2)),dual_norm=float(la.norm(w,2)),
               reduced_resolvent_norm=float(la.norm(S,2)))
    unit=EPS*dim*error['generator_norm']*error['reduced_resolvent_norm']*error['dual_norm']**2
    error['identity_relative_bound']=16*unit
    result=dict(n=n,dimension=dim,q=ep.base.encode(q),lambda0=ep.base.encode(lam),third_distance=float(gap),
        scalar_plane_residual=float(la.norm(k-lam*np.eye(2))),reflection_on_plane=ep.base.encode(w@u[reflection]),profiles={},roundoff=error)
    for name,(a,b) in ep.DIRECTIONS.items():
        V=a*left+b*right; B1=w@(q*V)@u; B2=w@(q*V)@S@(q*V)@u
        m=2 if name=='odd' else 1; B=B2 if m==2 else B1
        alpha=ep.base.discriminant(A); beta=polar(A,B); gamma=ep.base.discriminant(B)
        omega=beta**2-4*alpha*gamma
        slopes=ep.leading_slopes(A,B)
        ds=[graph_coefficients(d,c,V,q,lam,m,slope)[0] for slope in slopes]
        independent=generator_scale*independent_hop(n,[1+.125*a]+[1]*(n-3)+[1+.125*b])
        assert np.array_equal(c+.125*V,independent)
        result['profiles'][name]=dict(order=m,alpha=ep.base.encode(alpha),omega=ep.base.encode(omega),first_order_norm=float(la.norm(B1)),slopes=ep.base.encode(slopes),next_coefficients=ep.base.encode(ds),
            a_norm=float(la.norm(A,2)),effective_norm=float(la.norm(B,2)),perturbation_norm=float(la.norm(q*V,2)))
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
                rows.append(dict(epsilon=e,branch=branch,q=ep.base.encode(q),lambda0=ep.base.encode(lam),scaled_discriminant_residual=float(la.norm(f(sol.x))),traceless_norm=float(la.norm(k-lam*np.eye(2))),third_distance=float(np.sort(abs(la.eigvals(d+q*hop)-lam))[2]),q_displacement_error=float(abs(q-q0-eta*slope)/abs(q-q0)),smallest_singular_values=la.svdvals(d+q*hop-lam*np.eye(len(d)))[-3:].tolist(),generator_norm=float(la.norm(d+q*hop,2))))
        for e in sorted(set(row['epsilon'] for row in rows)):
            pair=[row for row in rows if row['epsilon']==e]
            separation=abs(z(pair[0]['q'])-z(pair[1]['q']))
            for row in pair:row['branch_separation']=float(separation)
        result[name]=rows
    return result


def build_coefficient_survey(n5_payload):
    """Compute exactly the two N4 anchors and all 58 N5 inventory loci."""
    rows=[]
    for identifier,(n,q,lam) in _expected_seeds(n5_payload).items():
        row=coefficients(n,q,lam)
        row['id']=identifier
        calibration=[]
        for scale in SCALE_FACTORS:
            measured=row if scale==1 else coefficients(n,q,lam,generator_scale=scale)
            _validate_coefficient_row(measured)
            unit=_roundoff_unit(measured)
            calibration.append(dict(
                generator_scale=scale,**measured['roundoff'],
                identity_relative_residual=measured['identity_relative_residual'],
                identity_ratio=measured['identity_relative_residual']/unit,
                alpha_min_clearance=min(abs(z(p['alpha']))/(16*unit*p['a_norm']**2)
                                        for p in measured['profiles'].values()),
                omega_min_clearance=min(abs(z(p['omega']))/(128*unit*p['a_norm']**2*p['effective_norm']**2)
                                        for p in measured['profiles'].values())))
        row['roundoff']['scale_checks']=calibration
        rows.append(row)
    validate_coefficient_survey(rows,n5_payload)
    return rows


def _expected_seeds(n5_payload):
    q=np.sqrt((-1+np.sqrt(13))/6)
    seeds={'N4-real-lambda-plus':(4,q,-4+2j*q),
           'N4-real-lambda-minus':(4,q,-4-2j*q)}
    try:
        loci=[(loc['id'],z(loc['qPhysicalCSharpSeed']),z(a['lambdaSeed']))
              for sector in n5_payload['sectors'] for a in sector['a2Roots']
              for loc in a['qLoci']]
    except (KeyError,TypeError,ValueError) as error:
        raise ValueError('N5 inventory coverage/schema is invalid') from error
    if len(loci)!=58 or len({item[0] for item in loci})!=58:
        raise ValueError('N5 inventory coverage must contain exactly 58 distinct loci')
    for identifier,q,lam in loci:
        if not identifier.startswith('N5-') or not np.isfinite([q,lam]).all():
            raise ValueError('Invalid N5 locus identity/seed')
        seeds[identifier]=(5,q,lam)
    return seeds


def _finite(value):
    """Reject NaN/Infinity anywhere in the numerical payload, including arrays."""
    if isinstance(value,dict):
        if 'real' in value or 'imag' in value:
            _decode_complex(value)
        for child in value.values(): _finite(child)
    elif isinstance(value,(list,tuple)):
        for child in value: _finite(child)
    elif isinstance(value,(float,int,np.number)) and not np.isfinite(value):
        raise ValueError('Nonfinite survey value')


def _numeric_component(value):
    if isinstance(value,(list,tuple)):
        return [_numeric_component(child) for child in value]
    return _finite_number(value)


def _decode_complex(value,shape=None):
    """An encoded scalar/vector/matrix has exactly real and imag components."""
    if not isinstance(value,dict) or set(value)!={'real','imag'}:
        raise ValueError('Invalid complex encoding schema: expected exactly real/imag keys')
    real=np.asarray(_numeric_component(value['real']))
    imag=np.asarray(_numeric_component(value['imag']))
    if real.shape!=imag.shape or (shape is not None and real.shape!=shape):
        raise ValueError(f'Invalid complex encoding shape: expected {shape}')
    return real+1j*imag


def _unpack(value):
    return _decode_complex(value,(2,))


def _roundoff_unit(row):
    """Dimensionless backward-error model, independent of the identity residual.

    Dense Schur/solve errors start at eps*dim*||L||, then the complementary
    inverse contributes ||S|| and the bilinear dual contributes ||w||^2.
    A 16-operation envelope covers the identity/alpha and 128 covers omega's
    degree-four polar products. These are empirical floating-point guards,
    not certified error bounds. The calibration is rerun over seven decades.
    """
    error=row['roundoff']
    scales=[error[name] for name in ('generator_norm','reduced_resolvent_norm','dual_norm')]
    if min(scales)<=0:
        raise ValueError('Nonpositive roundoff scale')
    unit=EPS*row['dimension']*scales[0]*scales[1]*scales[2]**2
    if error['identity_relative_bound']!=16*unit:
        raise ValueError('Identity bound disagrees with recorded roundoff law')
    return unit


def _validate_coefficient_row(row):
    _finite(row)
    z(row['q']); z(row['lambda0'])
    _decode_complex(row['reflection_on_plane'],(2,2))
    if 'residue_control' in row:
        control=row['residue_control']
        for mode in control['nearest_modes']:
            z(mode['eigenvalue']); z(mode['reflection_character'])
        for response in control['responses'].values():
            _unpack(response['odd_slopes'])
    unit=_roundoff_unit(row)
    if row['third_distance']<=16*EPS*row['dimension']*row['roundoff']['generator_norm']:
        raise ValueError('Pair lacks numerical isolation')
    if not 0<=row['scalar_plane_residual']<=16*unit*row['roundoff']['generator_norm']:
        raise ValueError('Scalar-plane residual exceeds roundoff law')
    if set(row['profiles'])!=set(PROFILE_ORDERS):
        raise ValueError('Expected one/even/odd profiles')
    for name,order in PROFILE_ORDERS.items():
        profile=row['profiles'][name]
        if profile['order']!=order:
            raise ValueError(f'Wrong {name} profile order: expected {order}')
        a,b,v=(profile[key] for key in ('a_norm','effective_norm','perturbation_norm'))
        if min(a,b,v)<=0:
            raise ValueError('Effective coefficient norm is zero')
        if abs(z(profile['alpha']))<=16*unit*a**2:
            raise ValueError('alpha does not clear measured roundoff law')
        if abs(z(profile['omega']))<=128*unit*a**2*b**2:
            raise ValueError('omega does not clear measured roundoff law')
        first_order=profile['first_order_norm']
        if first_order<0 or (first_order>16*unit*v if order==2 else first_order<=16*unit*v):
            raise ValueError('First-order norm disagrees with profile order')
        _unpack(profile['slopes']); _unpack(profile['next_coefficients'])
    absolute=row['identity_absolute_residual']; relative=row['identity_relative_residual']
    if not 0<=relative<=16*unit or absolute<0:
        raise ValueError('Profile identity exceeds independent roundoff law')
    one,even,odd=(row['profiles'][name] for name in ('one','even','odd'))
    cs1,cs2=_unpack(one['slopes']),_unpack(even['slopes'])
    permutation=min(((0,1),(1,0)),key=lambda pp:max(abs(cs1[j]-cs2[pp[j]]) for j in range(2)))
    differences=_unpack(one['next_coefficients'])-_unpack(even['next_coefficients'])[list(permutation)]
    slopes=_unpack(odd['slopes'])
    observed=float(min(max(abs(differences[j]-slopes[pp[j]]) for j in range(2))
                       for pp in ((0,1),(1,0))))
    if absolute!=observed or relative!=float(observed/max(1,max(abs(slopes)))):
        raise ValueError('Stored identity residual disagrees with coefficient data')


def validate_coefficient_survey(rows,n5_payload):
    """Reject missing/duplicate loci, changed seeds and failed numerical laws."""
    expected=_expected_seeds(n5_payload)
    try:
        if len(rows)!=60 or len({r['id'] for r in rows})!=60 or {r['id'] for r in rows}!=set(expected):
            raise ValueError('Coefficient survey coverage must be 2 N4 + 58 N5 = 60 distinct loci')
        for row in rows:
            n,q,lam=expected[row['id']]
            if (row['n']!=n or row['dimension']!={4:24,5:50}[n]
                    or z(row['q'])!=q or z(row['lambda0'])!=lam):
                raise ValueError('Coefficient locus disagrees with its inventory seed/dimension')
            _validate_coefficient_row(row)
            readings=row['roundoff']['scale_checks']
            if [r['generator_scale'] for r in readings]!=list(SCALE_FACTORS):
                raise ValueError('Missing roundoff calibration scales')
            for reading in readings:
                unit=_roundoff_unit(dict(dimension=row['dimension'],roundoff=reading))
                residual=reading['identity_relative_residual']
                if (not 0<=residual<=16*unit or reading['identity_ratio']!=residual/unit
                        or reading['alpha_min_clearance']<=1 or reading['omega_min_clearance']<=1):
                    raise ValueError('Scale calibration fails independent roundoff law')
    except (KeyError,TypeError,IndexError) as error:
        raise ValueError('Invalid coefficient survey schema') from error


def validate_continuations(continuations):
    """Require both distinct branches of every selected epsilon to clear noise.

    Schur/SVD noise scales as eps*dim*||L||. The discriminant of a near-Jordan
    2x2 block has first-order error proportional to its traceless norm; after
    dividing by epsilon^(2m), its bound carries that same denominator. Rebuild
    H from spin flips and D from the disagreement count, independently of the
    end-bond pencil used by follow(). Recompute the lambda-dependent SVD and
    Schur diagnostics. Weyl's absolute singular-value error is bounded here by
    16*eps*dim*(||L||+|lambda|); two such error intervals must fit inside the
    gap between the second and third smallest singular values. The Schur
    cluster projector norm supplies its non-normal conditioning factor. Each
    key fixes (n,q0,lambda0) from the coefficient-survey inventory. Its lambda
    disk is the same fixed isolated disk used by follow(); each q must lie in
    its own leading-predictor disk of radius half the predictor separation.
    These disjoint disks identify the measured local branches, not an
    all-epsilon continuation or an analytic remainder certificate.
    """
    try:
        _finite(continuations)
        if set(continuations)!=set(SELECTED):
            raise ValueError('Wrong selected continuation keys')
        anchors=_expected_seeds(json.loads(SOURCE.read_text(encoding='utf-8')))
        for identifier,profiles in continuations.items():
            n,q_seed,lambda_seed=anchors[identifier]
            seed=coefficients(n,q_seed,lambda_seed)
            seed_unit=_roundoff_unit(seed)
            dim={4:24,5:50}[n]
            basis=[(a,b) for a in weight_block_configs(n,1) for b in weight_block_configs(n,2)]
            decay=np.diag([-2*(a^b).bit_count() for a,b in basis])
            anchor_matrix=decay+q_seed*independent_hop(n,[1]*(n-1))
            anchor_radius=.4*np.sort(abs(la.eigvals(anchor_matrix)-lambda_seed))[2]
            if set(profiles)!=set(PROFILE_ORDERS):
                raise ValueError('Wrong continuation profiles')
            for name,m in PROFILE_ORDERS.items():
                rows=profiles[name]
                epsilons=((.001,.002,.003) if m==2 else (.0001,.0003,.001)) if n==4 else (.001,.003,.01)
                expected={(epsilon,branch) for epsilon in epsilons for branch in (0,1)}
                if len(rows)!=6 or {(r['epsilon'],r['branch']) for r in rows}!=expected:
                    raise ValueError('Continuation needs two branches per selected epsilon')
                for epsilon in epsilons:
                    pair=[r for r in rows if r['epsilon']==epsilon]
                    separation=float(abs(z(pair[0]['q'])-z(pair[1]['q'])))
                    q_floor=16*EPS*max(1,*(abs(z(r['q'])) for r in pair))
                    if separation<=q_floor or any(r['branch_separation']!=separation for r in pair):
                        raise ValueError('Continuation branches coalesced or separation is stale')
                a,b=ep.DIRECTIONS[name]
                slopes=_unpack(seed['profiles'][name]['slopes'])
                for row in rows:
                    q=z(row['q']); lam=z(row['lambda0'])
                    eta=row['epsilon']**m
                    slope=slopes[row['branch']]
                    q_error=(16*EPS*max(1,abs(q_seed),abs(q))
                             +16*seed_unit*eta*max(abs(slopes)))
                    predictor_residual=abs(q-q_seed-eta*slope)
                    if predictor_residual+q_error>=eta*abs(slopes[0]-slopes[1])/2:
                        raise ValueError('Continuation disagrees with selected seed q predictor')
                    if abs(lam-lambda_seed)>=anchor_radius:
                        raise ValueError('Continuation leaves selected seed spectral cluster')
                    if abs(q-q_seed)<=q_error:
                        raise ValueError('Continuation displacement from selected seed is unresolved')
                    relative=float(predictor_residual/abs(q-q_seed))
                    if abs(relative-row['q_displacement_error'])>2*q_error/abs(q-q_seed):
                        raise ValueError('Stored displacement disagrees with selected seed')
                    weights=[1+row['epsilon']*a]+[1]*(n-3)+[1+row['epsilon']*b]
                    matrix=decay+q*independent_hop(n,weights)
                    norm=float(la.norm(matrix,2))
                    if row['generator_norm']!=norm:
                        raise ValueError('Continuation norm disagrees with generator')
                    unit=EPS*dim*norm
                    bound=32*unit*(row['traceless_norm']+unit)/row['epsilon']**(2*m)
                    if not 0<=row['scaled_discriminant_residual']<=bound:
                        raise ValueError('Continuation discriminant exceeds roundoff law')
                    singular=row['smallest_singular_values']
                    if (len(singular)!=3 or singular!=sorted(singular,reverse=True)
                            or not 0<=singular[-1]<=16*unit or singular[-2]<=16*unit):
                        raise ValueError('Continuation singular-value clearance failed')
                    sv_error=16*EPS*dim*(norm+abs(lam))
                    if singular[-3]-singular[-2]<=2*sv_error:
                        raise ValueError('Continuation two-value singular cluster isolation failed')
                    actual_singular=la.svdvals(matrix-lam*np.eye(dim))[-3:]
                    if (actual_singular[-1]>sv_error or actual_singular[-2]<=sv_error
                            or actual_singular[-3]-actual_singular[-2]<=2*sv_error):
                        raise ValueError('Recomputed lambda-dependent singular-value clearance failed')
                    if np.max(abs(actual_singular-np.asarray(singular)))>2*sv_error:
                        raise ValueError('Stored singular values disagree with recorded lambda')
                    values,vl,vr=la.eig(matrix,left=True,right=True)
                    distances=abs(values-lam); ordered=np.argsort(distances)
                    third=distances[ordered[2]]
                    plane,k=ep.base.plane(matrix,lambda_seed,anchor_radius)
                    dual=la.solve(plane.T@plane,plane.T)
                    condition=max(1.,float(la.norm(dual,2)))
                    schur_error=16*unit*condition
                    if abs(np.trace(k)/2-lam)>schur_error:
                        raise ValueError('Recorded lambda disagrees with Schur cluster center')
                    actual_traceless=float(la.norm(k-lam*np.eye(2)))
                    if abs(actual_traceless-row['traceless_norm'])>2*schur_error:
                        raise ValueError('Stored traceless norm disagrees with recorded lambda')
                    actual_discriminant=float(abs(ep.base.discriminant(k))/row['epsilon']**(2*m))
                    disc_error=32*unit*condition*(actual_traceless+unit*condition)/row['epsilon']**(2*m)
                    if (actual_discriminant>disc_error
                            or abs(actual_discriminant-row['scaled_discriminant_residual'])>2*disc_error):
                        raise ValueError('Recomputed discriminant disagrees with recorded continuation')
                    third_index=ordered[2]
                    overlap=abs(vl[:,third_index].conj()@vr[:,third_index])
                    if overlap==0:
                        raise ValueError('Third eigenvalue has unresolved conditioning')
                    third_error=32*unit/overlap
                    if abs(third-row['third_distance'])>third_error:
                        raise ValueError('Stored third distance disagrees with recorded lambda')
                    if row['traceless_norm']<=16*unit or row['third_distance']<=16*unit:
                        raise ValueError('Continuation pair lacks isolation/Jordan clearance')
                    if row['q_displacement_error']<0:
                        raise ValueError('Negative displacement error')
    except (KeyError,TypeError,IndexError) as error:
        raise ValueError('Invalid continuation schema') from error


def validate_provenance(payload):
    """Compare canonical source/script/transitive hashes against current files."""
    verify_artifact_provenance(payload,Path(__file__),SOURCE)

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
    data=json.loads(SOURCE.read_text(encoding='utf-8'))
    rows=build_coefficient_survey(data)
    validate_coefficient_survey(rows,data)
    selected=[next(r for r in rows if r['id']==identifier) for identifier in SELECTED]
    output=dict(scope=__doc__,**artifact_provenance(Path(__file__),SOURCE),loci=rows,continuations={})
    for r in selected:
        print('selected',r['id'],'gap',r['third_distance'],flush=True)
        r['residue_control']=residue_control(r)
        output['continuations'][r['id']]=follow(r['n'],z(r['q']),z(r['lambda0']),r)
    validate_coefficient_survey(rows,data)
    validate_continuations(output['continuations'])
    validate_provenance(output)
    path=ROOT/'simulations/results/route_b_other_n_unfolding.json'
    path.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('count',len(rows),'worst relative identity',max(r['identity_relative_residual'] for r in rows),flush=True)
if __name__=='__main__':main()
