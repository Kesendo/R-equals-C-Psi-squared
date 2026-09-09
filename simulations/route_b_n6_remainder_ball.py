"""Uniform complex-ball contraction for F163 EP branches for three end profiles.
Uses F163 exact premises and acb-enclosed seed. eta=epsilon or epsilon^2.
Run: python simulations/route_b_n6_remainder_ball.py
"""
import hashlib,json,time
import flint
from pathlib import Path
from flint import acb,arb,acb_mat
from route_b_n6_ball_base import build_base,identity
ROOT=Path(__file__).resolve().parents[1]
def midpoint(x): return acb(x.real.mid(),x.imag.mid())
def diskbox(center,r): return center+acb(arb(0,r),arb(0,r))
def polar_det(a,b):
    return a[0,0]*b[1,1]+a[1,1]*b[0,0]-a[0,1]*b[1,0]-a[1,0]*b[0,1]

def equations(base,V,eta,c,l):
    U,W,Q=base['U'],base['W'],base['Q']
    I=identity(90); I2=identity(2)
    key=id(V)
    cache=base.setdefault('_projections',{})
    if key not in cache:
        cache[key]=dict(A=W*base['C']*U,Vp=W*V*U,BC=W*base['C']*Q,BV=W*V*Q,
            DC=Q*base['C']*U,DV=Q*V*U,HC=Q*base['C']*Q,HV=Q*V*Q)
    base.setdefault('_projection_inputs',{})[key]=V
    k=cache[key]; vcoef=base['q0']+eta*c
    B=c*k['BC']+vcoef*k['BV']; D=c*k['DC']+vcoef*k['DV']
    Bc=k['BC']+eta*k['BV']; Dc=k['DC']+eta*k['DV']; H=k['HC']+eta*k['HV']
    R=(base['base_matrix']+eta*(l*I-c*k['HC']-vcoef*k['HV'])).inv()
    RD=R*D; R2D=R*RD; R3D=R*R2D
    M=l*I2-c*k['A']-vcoef*k['Vp']-eta*B*RD
    Ml=I2+eta**2*B*R2D
    Mc=-k['A']-eta*k['Vp']-eta*(Bc*RD+B*R*Dc+eta*B*R*H*RD)
    Mll=-2*eta**3*B*R3D
    Mlc=eta**2*(Bc*R2D+B*R*R*Dc+eta*B*(R*H*R2D+R*R*H*RD))
    F=acb_mat([[M.det()],[polar_det(M,Ml)]])
    J=acb_mat([[polar_det(M,Mc),polar_det(M,Ml)],
        [polar_det(Mc,Ml)+polar_det(M,Mlc),polar_det(Ml,Ml)+polar_det(M,Mll)]])
    return F,J

def check(base,V,c0,l0,R,rho):
    eta=diskbox(acb(0),R)
    _,J0=equations(base,V,acb(0),c0,l0)
    Y0=J0.inv(); Y=acb_mat([[midpoint(Y0[i,j]) for j in range(2)] for i in range(2)])
    assert not Y.det().contains(0)
    F0,_=equations(base,V,eta,c0,l0)
    _,J=equations(base,V,eta,diskbox(c0,rho),diskbox(l0,rho))
    defect=identity(2)-Y*J
    residual=Y*F0
    rowbounds=[sum((abs(defect[i,j]).upper() for j in range(2)),arb(0)) for i in range(2)]
    images=[abs(residual[i,0]).upper()+rowbounds[i]*rho for i in range(2)]
    contraction=all(x<1 for x in rowbounds)
    inclusion=all(x<rho for x in images)
    return dict(success=contraction and inclusion,contraction=contraction,inclusion=inclusion,
        jacobian_row_bounds=[str(x) for x in rowbounds],newton_image_bounds=[str(x) for x in images],
        radius_eta=str(R),radius_variables=str(rho),c_center=str(c0),l_center=str(l0),
        c_abs_upper=str(abs(c0).upper()+rho),c_abs_lower=str(abs(c0).lower()-rho))

def odd_base(base):
    n=90; I=identity(n)
    reflection=acb_mat([[int(j==base['reflection'][i]) for j in range(n)] for i in range(n)])
    even=(I+reflection)/2; odd=(I-reflection)/2
    V=(base['left']-base['right'])/2
    low=odd*V*even; high=even*V*odd
    transformed=dict(base)
    transformed.pop('_projections',None)
    transformed.pop('_projection_inputs',None)
    transformed['L0']=base['L0']+base['q0']*low
    transformed['C']=base['C']+low
    transformed['U']=base['U']+base['q0']*base['S']*V*base['U']
    transformed['W']=base['W']
    transformed['P']=transformed['U']*transformed['W']
    transformed['Q']=I-transformed['P']
    transformed['base_matrix']=base['lambda0']*I-transformed['L0']+transformed['P']
    return transformed,high

def main():
    start=time.perf_counter(); base=build_base()
    source=ROOT/'simulations/results/route_b_n6_end_profile_probe.json'
    old=json.loads(source.read_text())
    output=dict(scope=__doc__,precision=base['precision'],seed_id=old['seed_id'],
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        base_script_sha256=hashlib.sha256((ROOT/'simulations/route_b_n6_ball_base.py').read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),flint_version=flint.__version__,
        premises={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'simulations/results/route_b_n6_exact_unfolding.json',ROOT/'simulations/results/route_b_a2_n6.json',ROOT/'simulations/tests/fixtures/route_b_a2_n6_residual.json']},profiles={})
    transformed,high=odd_base(base)
    profiles=[('one',base,base['left']),('even',base,(base['left']+base['right'])/2),('odd',transformed,high)]
    for name,current,V in profiles:
        entries=[]
        slopes=old['profiles'][name]['leading_slopes']
        for branch in range(2):
            c0=acb(slopes['real'][branch],slopes['imag'][branch])
            leading=current['W']*(c0*current['C']+current['q0']*V)*current['U']
            l0=midpoint((leading[0,0]+leading[1,1])/2)
            attempts=[]
            for exponent in (12,16,18,20,24,28):
                try: result=check(current,V,c0,l0,arb(2)**(-exponent),arb(2)**(-8))
                except (ValueError,ZeroDivisionError) as exc: result=dict(success=False,error=str(exc),radius_exponent=exponent)
                attempts.append(result)
                print(name,branch,exponent,'success',result['success'],'rows',[float(arb(x)) for x in result.get('jacobian_row_bounds',[])],flush=True)
                if result['success']: break
            if attempts[-1]['success']:
                rho=arb(2)**(-8)
                lower=abs(c0).lower()-rho
                # Cauchy applied to c(eta)-c_center, rather than the much larger q.
                result['half_radius_relative_remainder_upper']=str((rho/(2*lower)).upper())
                result['eta_order']=2 if name=='odd' else 1
                result['epsilon_radius_display']=str((arb(result['radius_eta'])**(arb(1)/result['eta_order'])))
                result['cauchy_formula']='relative remainder <= rho/(abs(c_center)-rho) * s^2/(1-s), s=abs(eta)/R < 1; exact Taylor coefficients'
            entries.append(dict(branch=branch,attempts=attempts))
        output['profiles'][name]=entries
    control_profile=profiles[0]
    control=check(control_profile[1],control_profile[2],acb(-4.006341325180786,-1.9824779247218287),acb(-1.81375382667802,2.427436411180607),arb(2)**(-20),arb(2)**(-8))
    assert not control['success'], 'Displaced center control unexpectedly certified'
    output['displaced_center_control']=control
    output['all_branches_certified']=all(any(a['success'] for a in row['attempts']) for entries in output['profiles'].values() for row in entries)
    output['elapsed_seconds']=time.perf_counter()-start
    (ROOT/'simulations/results/route_b_n6_remainder_ball.json').write_text(json.dumps(output,indent=2)+'\n')
    if not output['all_branches_certified']:
        raise RuntimeError('At least one branch has no certified radius; inspect failed attempts')
if __name__=='__main__': main()
