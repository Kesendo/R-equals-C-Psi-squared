"""Numerical one-end N4 EP continuation and candidate endpoint refinement.

Forty-digit polynomial continuation is not a maximal-interval certificate.
Labels 0/1 are the negative/positive leading q slope at epsilon=0.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
import numpy as np
import scipy.linalg as la
import sympy as s
from route_b_n4_range_polynomial import expression
from route_b_other_n_unfolding import parts

ROOT = Path(__file__).resolve().parents[1]


def main():
    mp.mp.dps = 40
    p, (x, y, a, b) = expression()
    factors = s.factor_list(p.subs(b, 1))[1]
    assert sorted((s.degree(f, x), k) for f, k in factors) == [(4, 1), (8, 1)]
    f4 = next(f for f, k in factors if s.degree(f, x) == 4)
    f = next(f for f, k in factors if s.degree(f, x) == 8)
    ev = s.lambdify((x, y, a), (f, s.diff(f, x), s.diff(f, x, 2),
                                  s.diff(f, y), s.diff(f, x, 3), f4), 'mpmath', cse=True)
    text = lambda v: mp.nstr(v, 32)
    encode = lambda v: dict(real=text(mp.re(v)), imag=text(mp.im(v)))
    q0 = mp.sqrt((mp.sqrt(13)-1)/6)
    slopes = [mp.mpf('-.5208808000870880'), mp.mpf('.1006463000731454')]
    tracks = []
    for branch, slope in enumerate(slopes):
        for sign in (-1, 1):
            xx, yy, ee = -4*q0*q0, q0*q0, mp.mpf(0)
            step = mp.mpf('.0001')
            history = []
            for iteration in range(420):
                en = ee+sign*step
                if sign > 0 and en > 1:
                    en = mp.mpf(1)
                if not history:
                    guess = (xx, (q0+slope*en)**2)
                elif len(history) > 1:
                    old = history[-2]
                    ratio = (en-ee)/(ee-old[0])
                    guess = (xx+ratio*(xx-old[1]), yy+ratio*(yy-old[2]))
                else:
                    guess = (xx, yy)
                try:
                    root = mp.findroot(lambda X,Y: ev(X,Y,1+en)[:2], guess,
                                       tol=mp.mpf('1e-30'), maxsteps=35)
                    if max(abs(root[j]-guess[j]) for j in (0, 1)) > mp.mpf('.08'):
                        raise ValueError('predictor jump')
                except (ValueError, ZeroDivisionError):
                    step /= 2
                    if step < mp.mpf('1e-8'):
                        break
                    continue
                xx, yy = root
                ee = en
                history.append((ee, xx, yy))
                step = min(step*mp.mpf('1.2'), mp.mpf('.005'))
                if yy <= 0 or xx >= 0 or abs(ee) >= 1:
                    break
            row = dict(branch=branch,epsilon_sign=sign,points=[
                dict(epsilon=text(e),x=text(X),y=text(Y)) for e,X,Y in history],
                termination='epsilon_1_reached' if ee == 1 else 'adaptive_step_floor_or_budget')
            assert history
            tracks.append(row)
            print('track', branch, sign, text(ee), flush=True)

    # Seeds from the preceding tracked exploration. Refined from equations here;
    # no claim that a nearby solution is automatically the first obstruction.
    proposals = [
        (0,-1,2,'-.00513827','-1.728116','.4400895'),
        (1,-1,2,'-.01373663','-1.666173','.4268662'),
        (1,1,3,'.26477261','-9.35151','.50082355')]
    d,c,left,_,_ = parts(4)
    endpoints = []
    for branch,sign,kind,e0,X0,Y0 in proposals:
        X,Y,A = mp.findroot(lambda X,Y,A: (ev(X,Y,A)[0],ev(X,Y,A)[1],ev(X,Y,A)[kind]),
                            tuple(map(mp.mpf,(X0,Y0,str(mp.mpf(e0)+1)))),
                            tol=mp.mpf('1e-32'),maxsteps=60)
        e = A-1
        vals = ev(X,Y,A)
        q, lam = mp.sqrt(Y), -4+1j*mp.sqrt(-X)
        matrix = d+complex(q)*(c+float(e)*left)-complex(lam)*np.eye(24)
        singular = la.svdvals(matrix)
        beyond_e = e+sign*mp.mpf('1e-6')
        beyond = []
        for imag_sign in (-1,1):
            BX,BY = mp.findroot(lambda X,Y: ev(X,Y,1+beyond_e)[:2],
                               (X+imag_sign*mp.mpc(0,'.003'),Y+imag_sign*mp.mpc(0,'.0001')),
                               tol=mp.mpf('1e-30'),maxsteps=80)
            residual=max(abs(v) for v in ev(BX,BY,1+beyond_e)[:2])
            beyond.append(dict(x=encode(BX),y=encode(BY),q=encode(mp.sqrt(BY)),
                               lambda_value=encode(-4+1j*mp.sqrt(-BX)),residual=text(residual)))
        endpoints.append(dict(branch=branch,epsilon_sign=sign,
            obstruction='Pxx_zero' if kind==2 else 'Py_zero',epsilon=text(e),x=text(X),y=text(Y),
            q=text(q),omega=text(mp.sqrt(-X)),equation_residual=text(max(abs(vals[j]) for j in (0,1,kind))),
            Pxx=text(vals[2]),Py=text(vals[3]),Pxxx=text(vals[4]),spectator_factor=text(vals[5]),
            smallest_full24_singular_values=singular[-4:].tolist(),
            beyond_epsilon=text(beyond_e),beyond=beyond))
        print('endpoint',branch,sign,text(e),flush=True)
    files=[Path(__file__),ROOT/'simulations/route_b_n4_range_polynomial.py',
           ROOT/'simulations/results/route_b_n4_range_polynomial.json',
           ROOT/'simulations/route_b_other_n_unfolding.py',
           ROOT/'simulations/framework/weight_coherence_block.py']
    output=dict(scope=__doc__,precision_digits=mp.mp.dps,tracks=tracks,endpoints=endpoints,
                source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (ROOT/'simulations/results/route_b_n4_range_one.json').write_text(json.dumps(output,indent=2)+'\n')


if __name__ == '__main__':
    main()
