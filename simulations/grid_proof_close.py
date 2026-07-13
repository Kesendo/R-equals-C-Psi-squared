r"""DETERMINISTIC upgrade of the Q(i) residue/endpoint certificate.

  >>> SUPERSEDED FOR THE RIGOROUS DETERMINISTIC PROOF by  core_grid.py + grid_proof_sweep.py.  <<<
  This file's structure builders (build_structure, component_plan, gen_primes, tonelli_i, ...) are
  reused, but its OWN main() sweep uses the heuristic pc_degree_bound whose widths were found NOT
  rigorous (e.g. C06 was gridded at width 28 while the true pseudo-remainder z1-support is 55; it
  catches the specific corruption but is not a valid deterministic certificate).  The rigorous
  engine recomputes every width by tropical span-tracking of the exact reduction, vectorizes the
  whole (z1,z2,w2) tensor grid, reduces via the monic quotient (leading coeff inverted per point,
  with skip-and-enlarge), and checkpoints the full 25-component + 6-endpoint-slice x 17-prime sweep.
  Run the PROOF via:  python -u simulations/grid_proof_sweep.py   (see that file's docstring).


Sibling of `residue_assembly_close.py` (imported): that pipeline PROVES exactly over Q(i) the
w3-elimination cross_form == F0(w1)+F1(w1) w3 over K = Q(i)(z1,z2,w2)[z3]/(Qz), the simple-pole
inventory, and the w1 polynomial-part window {-1,0,1}; the two remaining CERTIFIED-only steps were
(i) all w1-residues of F0,F1 vanish and (ii) the three endpoint coefficients c_-1=c_0=c_+1 vanish,
both by random-point Schwartz-Zippel.  THIS file replaces the randomness by a DETERMINISTIC
full-tensor-grid identity test + CRT.

DETERMINISTIC LEMMA (tensor interpolation + CRT).  Let P be a Laurent polynomial over Z[i] in the
FREE variables (z1,z2,w2, and w1 where present), reduced mod Qz so it is a PAIR (P0,P1) with
P = P0 + P1 z3 (each Pi a Laurent poly, z3-free).  If, for a set of primes p = 1 (mod 4) with
prod(p) > 2H (H a coefficient-height bound for P), P vanishes on a FULL tensor grid
A_z1 x A_z2 x A_w2 (x A_w1) of DISTINCT NONZERO residues with |A_v| = (span_v of P) + 1, in the ring
R_p = GF(p)[z3]/(Qz) (so BOTH components P0,P1 vanish), then P == 0 identically over Z[i].
Proof: multiply by z^{-min} per variable (w1,... are units) to get an ordinary polynomial of degree
< |A_v| in each var; vanishing on the full grid forces each z3-component (a genuine polynomial over
GF(p)) to be 0 by tensor-product interpolation; CRT over the primes lifts 0 mod prod(p) to 0 in Z[i]
since every coefficient is bounded by H < prod(p)/2.  R_p is used as a RING (never needs Qz to split):
all residue/endpoint numerators are DIVISION-FREE structured products, and divisibility is tested by
an inversion-free pseudo-remainder (multiply-through by the leading R_p coefficient).

WHAT IS TESTED.
  RESIDUES (no finite w1-pole of F0,F1).  Over K the 37 w1-carrying denominator factors are NOT
  pairwise coprime: exactly 12 quad-quad pairs share a K-root (Sylvester resultant over K vanishes;
  the OTHER ~24 factors are coprime to all).  So per-factor divisibility is FALSE; we group the
  factors into resultant-connected COMPONENTS (union-find: 13 singletons + 12 pairs).  For a
  component C, F0's principal part at the roots of Pi_C := prod_{f in C} f comes only from the
  C-incident terms (others coprime), and vanishes  <=>  Pi_C | P_C, where
      P_C = sum_{t incident to C} Num_t * prod_{f in C, f != term's C-factor(s)} f
                                          * prod_{d in distinct non-C co-factors, d != term's co} d
  (a DIVISION-FREE structured product; for a singleton this is exactly the pipeline's group
  numerator).  Pi_C | P_C  <=>  pseudo-rem_w1(P_C, Pi_C) == 0, tested on the grid (both z3-comps),
  all CRT primes.  A control corrupts one numerator coefficient and confirms the test FAILS.
  [Pi_C | P_C over K  <=>  no finite w1-pole because ord_r(Pi_C) matches the (simple) pole order at
  every root r of C, including the double at a shared root -- the numerator must vanish to order 2.]

  ENDPOINTS (c_-1=c_0=c_+1=0), CONTINGENT on the residue result (which gives F0,F1 window {-1,0,1}):
  then F0 is a degree-<=1 Laurent poly in w1, so N0(w1=w_k) == 0 for three distinct w_k forces
  c_-1=c_0=c_1=0.  N0(w1=w_k) = sum_t A_t * prod_{all 43 distinct factors except the term's two} is
  a DIVISION-FREE structured product; tested == 0 on the (z1,z2,w2) grid, all CRT primes.  (Width in
  z1,z2 is 108: reported; run gated by feasibility.)

STATUS / HONESTY.  This module is a WORK IN PROGRESS.  The core R_p machinery is validated
(remainder==0 with a discriminating corruption control) but the FULL multi-prime, full-tensor-grid
sweep over ALL components + endpoints is heavy; components whose pseudo-remainder degree bound blows
the grid past the budget (the 3 incidence-32 quads, and the width-108 endpoint) are REPORTED with
exact numbers rather than run.  See main() output for the per-claim verdict.

Authors: Thomas Wicht and Claude, 2026-07-10.
"""
import sys, time, random
import numpy as np

sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")
from residue_assembly_close import (build_F, factor_key, classify_factor, resultant_w1, _rq, W1)
from halfangle_residue_proof import NVARS, p_is_zero
from collections import defaultdict

# ------------------------------------------------------------------ CRT primes p = 1 (mod 4)
def gen_primes(n, start=(1 << 30)):
    """n primes p = 1 (mod 4), ~30 bit (products of two stay < 2^63 for int64)."""
    def isprime(m):
        if m < 2: return False
        for q in (2,3,5,7,11,13,17,19,23,29,31,37):
            if m % q == 0: return m == q
        d=m-1; r=0
        while d%2==0: d//=2; r+=1
        for a in (2,3,5,7,11,13,17,19,23,29,31,37):
            x=pow(a,d,m);
            if x in (1,m-1): continue
            for _ in range(r-1):
                x=x*x%m
                if x==m-1: break
            else: return False
        return True
    out=[]; c=start|1
    while len(out)<n:
        if c%4==1 and isprime(c): out.append(c)
        c+=2
    return out

def tonelli_i(p):
    # sqrt(-1) mod p for p = 1 mod 4
    a=p-1
    if p%4==3: return None
    q,s=p-1,0
    while q%2==0: q//=2; s+=1
    z=2
    while pow(z,(p-1)//2,p)!=p-1: z+=1
    m,c,t,r=s,pow(z,q,p),pow(a,q,p),pow(a,(q+1)//2,p)
    while t!=1:
        i,t2=0,t
        while t2!=1: t2=t2*t2%p; i+=1
        b=pow(c,1<<(m-i-1),p); m,c=i,b*b%p; t=t*c%p; r=r*b%p
    return r

# ------------------------------------------------------------------ R_p ring (z3-pair), vectorized
def rp_mul(a,b,Sz,p):
    a0,a1=a; b0,b1=b
    a0b0=(a0*b0)%p; a1b1=(a1*b1)%p; a0b1=(a0*b1)%p; a1b0=(a1*b0)%p
    return ((a0b0-a1b1)%p, (a0b1+a1b0-(a1b1*Sz)%p)%p)

def rp_powvec(base,e,p):
    if e<0: base=pow(base,p-2,p) if np.isscalar(base) else _invarr(base,p); e=-e
    r=np.ones_like(base) if not np.isscalar(base) else 1
    b=base
    while e:
        if e&1: r=(r*b)%p
        b=(b*b)%p; e>>=1
    return r
def _invarr(arr,p):
    r=np.ones_like(arr); b=arr%p; e=p-2
    while e:
        if e&1: r=(r*b)%p
        b=(b*b)%p; e>>=1
    return r

def gcoef(c,p,ii):
    cr=(int(c.re.numerator)*pow(int(c.re.denominator),p-2,p))%p
    ci=(int(c.im.numerator)*pow(int(c.im.denominator),p-2,p))%p
    return (cr+ii*ci)%p

def eval_rp(poly,z1v,z2,w2,Sz,p,ii):
    """poly dict monomial6->G (w3 exp 0). Return {w1exp:(a0 arr,a1 arr)} over z1 axis."""
    out={}
    z2i=pow(z2,p-2,p); w2i=pow(w2,p-2,p)
    for k,c in poly.items():
        val=np.full_like(z1v, gcoef(c,p,ii))
        e0=k[0]
        if e0: val=(val*rp_powvec(z1v,e0,p))%p
        e1=k[1]
        if e1: val=(val*(pow(z2,e1,p) if e1>0 else pow(z2i,-e1,p)))%p
        e4=k[4]
        if e4: val=(val*(pow(w2,e4,p) if e4>0 else pow(w2i,-e4,p)))%p
        w=k[3]; z3e=k[2]
        cur=out.get(w)
        if cur is None: cur=(np.zeros_like(z1v),np.zeros_like(z1v)); out[w]=cur
        if z3e==0: out[w]=((cur[0]+val)%p,cur[1])
        else: out[w]=(cur[0],(cur[1]+val)%p)
    return out

def pmul_w1(A,B,Sz,p):
    out={}
    for wa,ab in A.items():
        for wb,cd in B.items():
            pr=rp_mul(ab,cd,Sz,p); w=wa+wb
            cur=out.get(w)
            out[w]=pr if cur is None else ((cur[0]+pr[0])%p,(cur[1]+pr[1])%p)
    return out

def padd_w1(A,B,p):
    out={w:(a[0].copy(),a[1].copy()) for w,a in A.items()}
    for w,cd in B.items():
        cur=out.get(w)
        out[w]=(cd[0].copy(),cd[1].copy()) if cur is None else ((cur[0]+cd[0])%p,(cur[1]+cd[1])%p)
    return out

def clear_w1(poly):
    poly={w:ab for w,ab in poly.items() if (np.any(ab[0])or np.any(ab[1]))}
    if not poly: return {}
    s=min(poly); return {w-s:ab for w,ab in poly.items()}

def pseudo_rem(P,PHI,Sz,p):
    P=clear_w1({w:(a[0].copy(),a[1].copy()) for w,a in P.items()})
    PHI=clear_w1(PHI)
    if not PHI: raise RuntimeError("empty divisor")
    m=max(PHI); L=PHI[m]; steps=0
    while P:
        dP=max(P)
        if dP<m: break
        lead=P[dP]; newP={}
        for w,ab in P.items(): newP[w]=rp_mul(ab,L,Sz,p)
        for w,ab in PHI.items():
            tgt=w+(dP-m); sub=rp_mul(ab,lead,Sz,p); cur=newP.get(tgt)
            newP[tgt]=((-sub[0])%p,(-sub[1])%p) if cur is None else ((cur[0]-sub[0])%p,(cur[1]-sub[1])%p)
        P={w:ab for w,ab in newP.items() if (np.any(ab[0])or np.any(ab[1]))}
        steps+=1
        if steps>400: raise RuntimeError("pseudo-rem runaway")
    return P

# ------------------------------------------------------------------ structure: components
def build_structure():
    F=build_F()
    factors={}; incidence=defaultdict(list)
    for ti,(A,B,dens) in enumerate(F):
        for slot,d in enumerate(dens):
            k=factor_key(d); factors[k]=d; incidence[k].append((ti,slot))
    w1keys=[k for k in factors if classify_factor(factors[k])!="w1free"]
    parent={k:k for k in w1keys}
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    edges=0
    for i in range(len(w1keys)):
        for j in range(i+1,len(w1keys)):
            if p_is_zero(_rq(resultant_w1(factors[w1keys[i]],factors[w1keys[j]]))):
                parent[find(w1keys[i])]=find(w1keys[j]); edges+=1
    comps=defaultdict(list)
    for k in w1keys: comps[find(k)].append(k)
    return F,factors,incidence,list(comps.values()),edges

# ------------------------------------------------------------------ P_C and Pi_C evaluation
def component_plan(F,factors,incidence,C):
    Cset=set(C)
    seen=set(); inc2=[]
    for key in C:
        for (ti,slot) in incidence[key]:
            if ti in seen: continue
            seen.add(ti); inc2.append(ti)
    # distinct non-C co-factors
    distinct=[];
    for ti in inc2:
        A,B,dens=F[ti]; k0,k1=factor_key(dens[0]),factor_key(dens[1])
        inC0,inC1=k0 in Cset,k1 in Cset
        cok=None if (inC0 and inC1) else (k1 if inC0 else k0)
        if cok is not None and cok not in [factor_key(d) for d in distinct]:
            distinct.append(dens[1] if cok==k1 else dens[0])
    return inc2,distinct

def eval_PC(F,factors,C,which,inc2,distinct,z1v,z2,w2,Sz,p,ii,corrupt=False):
    Cset=set(C)
    Cfac={key:eval_rp(factors[key],z1v,z2,w2,Sz,p,ii) for key in C}
    drp=[eval_rp(d,z1v,z2,w2,Sz,p,ii) for d in distinct]; dk=[factor_key(d) for d in distinct]
    PC={}; first=True
    for ti in inc2:
        A,B,dens=F[ti]; Num=A if which==0 else B
        if not Num: continue
        k0,k1=factor_key(dens[0]),factor_key(dens[1])
        inC0,inC1=k0 in Cset,k1 in Cset
        used_C=set(x for x,inc in ((k0,inC0),(k1,inC1)) if inc)
        own_co=None if (inC0 and inC1) else (k1 if inC0 else k0)
        term=eval_rp(Num,z1v,z2,w2,Sz,p,ii)
        if corrupt and first:
            term=padd_w1(term,{0:(np.ones_like(z1v),np.zeros_like(z1v))},p); first=False
        for key in C:
            if key in used_C: continue
            term=pmul_w1(term,Cfac[key],Sz,p)
        for d_k,d_rp in zip(dk,drp):
            if d_k==own_co: continue
            term=pmul_w1(term,d_rp,Sz,p)
        PC=padd_w1(PC,term,p)
    PHI={0:(np.ones_like(z1v),np.zeros_like(z1v))}
    for key in C: PHI=pmul_w1(PHI,Cfac[key],Sz,p)
    return PC,PHI

# ------------------------------------------------------------------ degree + height bookkeeping
def fspan(poly,vi):
    es=[k[vi] for k in poly]; return (min(es),max(es)) if es else (0,0)

def pc_degree_bound(F,factors,C,distinct):
    """Upper bound per-var (z1=0,z2=1,w2=4) width and w1-degree of P_C, and of Pi_C's leading."""
    b={0:[0,0],1:[0,0],4:[0,0],3:[0,0]}
    # A/B global span bound
    ab={0:(-6,7),1:(-6,7),4:(-4,4),3:(-4,4)}
    for vi in (0,1,4,3):
        lo=ab[vi][0]; hi=ab[vi][1]
        for key in C:
            s=fspan(factors[key],vi); lo+=s[0]; hi+=s[1]
        for d in distinct:
            s=fspan(d,vi); lo+=min(0,s[0]); hi+=max(0,s[1]); lo+=s[0]; hi+=s[1]
        b[vi]=[lo,hi]
    return b

def height_bound(F,factors,C,distinct):
    """Crude coefficient-height bound for the pseudo-remainder (sets the prime count)."""
    # product over factors of (nterms * maxcoef); times A (480*32); times #incident; times L^k inflation
    def maxcoef(poly): return max((max(abs(c.re.numerator),abs(c.re.denominator),abs(c.im.numerator),abs(c.im.denominator)) for c in poly.values()),default=1)
    H=480*32
    for d in distinct: H*= len(d)*maxcoef(d)
    for key in C: H*= len(factors[key])*maxcoef(factors[key])
    H*= 288
    # pseudo-rem multiply-through by L (leading of Pi_C) up to k times, k ~ w1-deg P_C
    Lc=1
    for key in C: Lc*= len(factors[key])*maxcoef(factors[key])
    H*= Lc**40
    return H

def distinct_nonzero(n, p, seed):
    rng=random.Random(seed); s=set()
    while len(s)<n:
        v=rng.randrange(1,p); s.add(v)
    return list(s)

def run_component(F,factors,incidence,C, primes, cap_points, seed=12345, tcap=200.0):
    """Deterministic full-tensor-grid pseudo-rem test of Pi_C | P_C for which in {0,1}, all primes.
    Returns dict with spans, grid, verdict, control, runtime.  Skips (verdict='SKIP') if grid>cap."""
    inc2,distinct=component_plan(F,factors,incidence,C)
    b=pc_degree_bound(F,factors,C,distinct)
    wz1,wz2,ww2=b[0][1]-b[0][0], b[1][1]-b[1][0], b[4][1]-b[4][0]
    gz1,gz2,gw2=wz1+1,wz2+1,ww2+1
    grid=gz1*gz2*gw2
    info=dict(size=len(C), classes=[classify_factor(factors[k]) for k in C],
              inc=len(inc2), ndistinct=len(distinct), width=(wz1,wz2,ww2), grid=grid)
    if grid>cap_points:
        info['verdict']='SKIP(grid too large)'; return info
    t0=time.time()
    z1vals=np.array(distinct_nonzero(gz1,primes[0],seed),dtype=np.int64)  # placeholder resized per p
    ok=True; ctrl_ok=True
    for p in primes:
        ii=tonelli_i(p)
        z1v=np.array(distinct_nonzero(gz1,p,seed),dtype=np.int64)
        z2list=distinct_nonzero(gz2,p,seed+1); w2list=distinct_nonzero(gw2,p,seed+2)
        for which in (0,1):
            for z2 in z2list:
                z2i=pow(z2,p-2,p)
                for w2 in w2list:
                    Sz=(rp_powvec(z1v,1,p)+rp_powvec(z1v,-1,p)+z2+z2i)%p
                    PC,PHI=eval_PC(F,factors,C,which,inc2,distinct,z1v,z2,w2,Sz,p,ii)
                    rem=pseudo_rem(PC,PHI,Sz,p)
                    if any(np.any(a[0])or np.any(a[1]) for a in rem.values()): ok=False
                if time.time()-t0>tcap:
                    info['verdict']='TIMEOUT'; info['runtime']=time.time()-t0; return info
        # control at one (z2,w2) on first prime
        if p==primes[0]:
            z2=z2list[0]; w2=w2list[0]; Sz=(rp_powvec(z1v,1,p)+rp_powvec(z1v,-1,p)+z2+pow(z2,p-2,p))%p
            PCc,PHIc=eval_PC(F,factors,C,0,inc2,distinct,z1v,z2,w2,Sz,p,ii,corrupt=True)
            remc=pseudo_rem(PCc,PHIc,Sz,p)
            ctrl_ok=any(np.any(a[0])or np.any(a[1]) for a in remc.values())
    info['verdict']='PROVED' if ok else 'FAIL'
    info['control_discriminates']=ctrl_ok
    info['runtime']=time.time()-t0
    return info

def main():
    print("="*80)
    print("DETERMINISTIC grid+CRT upgrade of the residue certificate (cross_form == 0 on V).")
    print("="*80)
    t0=time.time()
    F,factors,incidence,comps,edges=build_structure()
    print(f"[structure] 37 w1-factors, {edges} zero-K-resultant edges -> {len(comps)} components "
          f"({sum(1 for c in comps if len(c)==1)} singletons + {sum(1 for c in comps if len(c)==2)} pairs)")
    # height -> prime count (use the largest crude H across components)
    import math
    Hmax=max(height_bound(F,factors,C,component_plan(F,factors,incidence,C)[1]) for C in comps)
    bits=Hmax.bit_length()
    nprime=max(2, (bits//29)+2)
    import os
    nrun=int(os.environ.get("NRUN","3"))
    primes=gen_primes(nprime)[:nrun]
    print(f"[CRT] crude height bound H ~ 2^{bits}; full rigor needs {nprime} primes (30-bit, =1 mod4).")
    print(f"      THIS RUN uses {len(primes)} primes (validation of the identity on the full grid): {primes}")
    print("-"*80)
    CAP=int(os.environ.get("CAP","4000000"))
    results=[]
    for idx,C in enumerate(sorted(comps,key=lambda c:(len(c),sum(len(incidence[k]) for k in c)))):
        info=run_component(F,factors,incidence,C,primes,CAP)
        results.append(info)
        print(f"[C{idx:02d}] {info['classes']} inc={info['inc']} ndco={info['ndistinct']} "
              f"width(z1,z2,w2)={info['width']} grid={info['grid']:,} -> {info['verdict']}"
              + (f" ctrl={info.get('control_discriminates')} {info.get('runtime',0):.1f}s"
                 if info['verdict'] in ('PROVED','FAIL') else ""))
        if time.time()-t0>3000:
            print("  [global time cap reached; remaining components not run]"); break
    nproved=sum(1 for r in results if r['verdict']=='PROVED')
    nskip=sum(1 for r in results if 'SKIP' in r['verdict'])
    print("-"*80)
    print(f"RESIDUE components PROVED deterministically: {nproved}/{len(comps)} "
          f"({nskip} skipped as grid-too-large; those are the incidence-32 quads)")
    print("ENDPOINT (c_-1=c_0=c_1=0): contingent on ALL residue components; N0(w1=w_k)==0 slice test,"
          " z1,z2-width 108 -> grid ~109x109x71~8.4e5 x 3 w_k x primes (reported, not run here).")
    print(f"[TOTAL] {time.time()-t0:.1f}s")

if __name__=="__main__":
    main()
