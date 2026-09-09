"""Complex eta-circle exploration of the corrected F163 remainder.
Not a validated bound: sampled maxima, floating-point solves and closure tests
cannot prove holomorphy inside a circle. eta=epsilon for one/even, epsilon^2 for odd.
Run: python simulations/route_b_n6_complex_remainder.py
"""
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy.linalg as la
from scipy.optimize import root
import route_b_n6_end_profile_probe as ep
ROOT=Path(__file__).resolve().parents[1]
def z(x): return complex(x['real'],x['imag'])
def main():
    source=ROOT/'simulations/results/route_b_n6_next_coefficient.json'
    coefficients=json.loads(source.read_text())
    original=json.loads((ROOT/'simulations/results/route_b_n6_end_profile_probe.json').read_text())
    q0,lam0=z(original['q0']),z(original['lambda0'])
    d,c,left=ep.base.parts(); right,_=ep.reflected_end(left)
    radius=original['initial_third_distance']*0.4
    output=dict(scope=__doc__,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),profiles={})
    for name,(a,b) in ep.DIRECTIONS.items():
        v=a*left+b*right
        m=coefficients['profiles'][name]['order']
        circles=[]
        for R in ((0.02,0.05) if m==1 else (0.0025,0.01)):
            for branch in coefficients['profiles'][name]['branches']:
                slope=z(branch['slope']); correction=z(branch['next_coefficient'])
                rows=[]; failure=None
                previous=slope+correction*R
                for k in range(65):
                    angle=2*np.pi*k/64
                    eta=R*np.exp(1j*angle)
                    # Continuous square root, not principal-branch jumps.
                    epsilon=eta if m==1 else np.sqrt(R)*np.exp(0.5j*angle)
                    hop=c+epsilon*v
                    def objective(x):
                        _,block=ep.base.plane(d+(q0+eta*complex(*x))*hop,lam0,radius)
                        value=ep.base.discriminant(block)/eta**2
                        return np.array([value.real,value.imag])
                    def jac(x):
                        h=1e-3
                        return np.column_stack([(objective(x+h*np.eye(2)[j])-objective(x-h*np.eye(2)[j]))/(2*h) for j in range(2)])
                    try:
                        sol=root(objective,[previous.real,previous.imag],jac=jac,tol=1e-8)
                        if not sol.success: raise ValueError(str(sol.message))
                        previous=complex(*sol.x)
                        q=q0+eta*previous
                        u,block=ep.base.plane(d+q*hop,lam0,radius)
                        lam=np.trace(block)/2
                        distances=np.sort(abs(la.eigvals(d+q*hop)-lam))
                        g=(previous-slope-correction*eta)/eta**2
                        rows.append(dict(k=k,eta=ep.base.encode(eta),epsilon=ep.base.encode(epsilon),q=ep.base.encode(q),scaled_remainder=ep.base.encode(g),third_distance=float(distances[2]),scaled_discriminant_residual=float(la.norm(objective(sol.x)))))
                    except ValueError as exc:
                        failure=dict(k=k,message=str(exc)); break
                entry=dict(eta_radius=R,epsilon_radius=R**(1/m),branch=branch['branch'],rows=rows,failure=failure)
                if failure is None:
                    values=np.array([z(row['scaled_remainder']) for row in rows[:-1]])
                    fourier=np.fft.fft(values)/64
                    entry.update(sampled_max_64=float(max(abs(values))),sampled_max_32=float(max(abs(values[::2]))),sampled_min_third_distance=min(row['third_distance'] for row in rows),q_closure=float(abs(z(rows[-1]['q'])-z(rows[0]['q']))),negative_fourier_max=float(max(abs(fourier[33:]))),estimated_next_coefficient=ep.base.encode(fourier[0]))
                    print(name,R,branch['branch'],'max',entry['sampled_max_64'],'closure',entry['q_closure'],'gap',entry['sampled_min_third_distance'],flush=True)
                else: print(name,R,branch['branch'],'FAILED',failure,flush=True)
                circles.append(entry)
        output['profiles'][name]=circles
    (ROOT/'simulations/results/route_b_n6_complex_remainder.json').write_text(json.dumps(output,indent=2)+'\n')
if __name__=='__main__': main()
