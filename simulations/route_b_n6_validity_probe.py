"""Numerical validity of F163 at N6-E-A2-T-007; no certified epsilon radius.
Run with OPENBLAS_NUM_THREADS=1: python simulations/route_b_n6_validity_probe.py
Errors are relative to actual branch displacement, not total q.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy.linalg as la
import route_b_n6_end_profile_probe as ep

ROOT=Path(__file__).resolve().parents[1]
def z(x): return complex(x['real'],x['imag'])
def main():
    source=ROOT/'simulations/results/route_b_n6_end_profile_probe.json'
    old=json.loads(source.read_text())
    q0,lam0=z(old['q0']),z(old['lambda0'])
    d,c,left=ep.base.parts()
    right,_=ep.reflected_end(left)
    l0=d+q0*c
    isolation=old['initial_third_distance']
    v0,_=ep.base.plane(l0,lam0,isolation/2)
    # Singular values are 1-Lipschitz in z. Subtract the maximum distance
    # to the nearest angular node; this covers the entire circle in exact arithmetic.
    # Floating point SVD values here are not interval-enclosed.
    contour=[]
    for fraction in (0.25,0.4,0.5,0.6,0.75):
        r=isolation*fraction
        n=512
        sampled=min(la.svdvals((lam0+r*np.exp(2j*np.pi*k/n))*np.eye(90)-l0)[-1] for k in range(n))
        lower=sampled-2*r*np.sin(np.pi/(2*n))
        contour.append(dict(radius=r,nodes=n,sampled_min_sigma=float(sampled),angular_correction=float(2*r*np.sin(np.pi/(2*n))),continuum_lower_estimate=float(lower)))
    best=max(contour,key=lambda x:x['continuum_lower_estimate'])
    report=dict(seed_id=old['seed_id'],source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope=__doc__,contours=contour,profiles={})
    ep.EPSILONS=(0.01,0.02,0.03,0.05,0.075,0.1)
    for name,(a,b) in ep.DIRECTIONS.items():
        original=old['profiles'][name]
        slopes=np.array(original['leading_slopes']['real'])+1j*np.array(original['leading_slopes']['imag'])
        perturb=a*left+b*right
        result=ep.trace_profile(d,c,perturb,q0,lam0,isolation,v0,original['order'],slopes)
        out=[]
        for row in result['rows']:
            e=row['epsilon']; branch=row['branch']; q=z(row['q'])
            prediction=q0+slopes[branch]*e**original['order']
            err=abs(q-prediction)/abs(q-q0)
            norm=la.norm((q-q0)*c+q*e*perturb,2)
            out.append(dict(epsilon=e,branch=branch,q=row['q'],lambda_value=row['lambda'],local_signature=row.get('local_signature'),contours=row.get('contours'),relative_displacement_error=float(err),third_distance=row['third_distance'],scaled_discriminant_residual=row['scaled_discriminant_residual'],eigenvalue_gap=row['eigenvalue_gap'],traceless_compression_norm=row['traceless_compression_norm'],perturbation_norm=float(norm),fixed_contour_neumann_ratio=float(norm/best['continuum_lower_estimate'])))
        report['profiles'][name]=dict(order=original['order'],rows=out)
        print(name,[(r['epsilon'],r['branch'],round(r['relative_displacement_error'],5),round(r['third_distance'],5)) for r in out if r['epsilon']>=0.01],flush=True)
    dest=ROOT/'simulations/results/route_b_n6_validity_probe.json'
    dest.write_text(json.dumps(report,indent=2)+'\n')
    print('best continuum sigma estimate',best,flush=True)
if __name__=='__main__': main()
