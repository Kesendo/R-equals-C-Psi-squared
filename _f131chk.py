import numpy as np
from numpy import kron
I2=np.eye(2); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1]).astype(complex)
n=3; d=8; d2=64
def op(s,site):
    m=s if site==0 else I2
    for k in range(1,n): m=kron(m, s if k==site else I2)
    return m
Zs=[op(Z,l) for l in range(n)]
def superop(js,gs,hs):
    H=np.zeros((d,d),complex)
    for b in range(n-1): H+=js[b]*(op(X,b)@op(X,b+1)+op(Y,b)@op(Y,b+1))
    for l in range(n): H+=hs[l]*Zs[l]
    L=-1j*(kron(H,np.eye(d))-kron(np.eye(d),H.T))
    for l in range(n): L+=gs[l]*(kron(Zs[l],Zs[l].T)-np.eye(d2))
    return L
# R = bit reversal
R=np.zeros((d,d))
for c in range(d):
    rev=0
    for s in range(n):
        if (c>>s)&1: rev|=1<<(n-1-s)
    R[rev,c]=1
RR=kron(R,R)
t=0.37
jB=[1,1]; gB=[0.05,0.08,0.05]; gO=[-0.1,0,0.1]; hB=[0.3,0.3,0.3]
def scan(b,dr,s): return [b[i]+s*dr[i] for i in range(len(b))]
# conjugation identity gamma axis
Lp=superop(jB,scan(gB,gO,t),hB); Lm=superop(jB,scan(gB,gO,-t),hB)
print("conj gamma dev:", np.max(np.abs(RR@Lp@RR-Lm)))
# prep rho_even, readouts
u,v=1,1<<(n-1)
re=np.zeros(d2,complex)
for a,b in[(u,u),(v,v),(u,v),(v,u)]: re[a*d+b]=0.5
Oev=Zs[0]+Zs[n-1]; Oodd=Zs[0]-Zs[n-1]
def moments(L,rho,O,K):
    out=[]; w=rho.copy()
    for k in range(K+1):
        tr=sum(O[i,j]*w[j*d+i] for i in range(d) for j in range(d))
        out.append(tr.real)
        if k<K: w=L@w
    return np.array(out)
mep=moments(Lp,re,Oev,6); mem=moments(Lm,re,Oev,6)
mop=moments(Lp,re,Oodd,6); mom=moments(Lm,re,Oodd,6)
print("even-cell max|m(+t)-m(-t)|:", np.max(np.abs(mep-mem)))
print("odd-cell  max|m(+t)+m(-t)|:", np.max(np.abs(mop+mom)), " max|m|:",np.max(np.abs(mop)))
# zero cell: base generator, R-odd readout
Lb=superop(jB,gB,hB); mz=moments(Lb,re,Oodd,6)
print("zero-cell max|m|:", np.max(np.abs(mz)))
