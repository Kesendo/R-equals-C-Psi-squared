"""Hunt the closed form of Q*(N): identify the EP mode that freezes at the horizon.
At Q* a slow complex pair merges into a real pair. Dump that mode's structure: its
<n_XY>, and the n_diff histogram of its computational-basis coherences |a><b| (popcount
of a^b). The half-lit clue (<n_XY>~1/2) should resolve into which sectors it bridges."""
import numpy as np
I2=np.eye(2); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1]).astype(complex)
def site(op,l,N):
    m=np.array([[1.+0j]])
    for k in range(N): m=np.kron(m, op if k==l else I2)
    return m
def L_super(N,J,g):
    d=2**N; Id=np.eye(d); H=np.zeros((d,d),complex)
    for b in range(N-1): H+=(J/2)*(site(X,b,N)@site(X,b+1,N)+site(Y,b,N)@site(Y,b+1,N))
    L=-1j*(np.kron(H,Id)-np.kron(Id,H.T))
    for l in range(N):
        Zl=site(Z,l,N); L+=g*(np.kron(Zl,Zl.conj())-np.kron(Id,Id))
    return L
def slowest_osc_pair(N,J,g):
    ev,evec=np.linalg.eig(L_super(N,J,g))
    osc=np.where((ev.real<-1e-7)&(np.abs(ev.imag)>1e-6))[0]
    i=osc[np.argmax(ev[osc].real)]; return ev[i],evec[:,i]
def Qstar(N,J=1.0):
    lo,hi=0.3,1.2
    for _ in range(40):
        m=0.5*(lo+hi)
        ev=np.linalg.eigvals(L_super(N,J,m)); nz=ev[ev.real<-1e-7]
        gap=nz.real.max(); band=nz[np.abs(nz.real-gap)<1e-6]
        if np.abs(band.imag).max()>1e-6: lo=m
        else: hi=m
    return J/(0.5*(lo+hi))
def ndiff_hist(N,vec):
    d=2**N; M=vec.reshape(d,d); w=np.abs(M)**2; w=w/(w.sum()+1e-30)
    h={}
    for a in range(d):
        for b in range(d):
            if w[a,b]>1e-6:
                nd=bin(a^b).count("1"); h[nd]=h.get(nd,0)+w[a,b]
    nxy=sum(nd*p for nd,p in h.items())
    return nxy,{k:round(v,3) for k,v in sorted(h.items())}
print(f"{'N':>2} {'Q*':>8} {'g_eff=2/Q*':>10} {'EP lam (just>Q*)':>22} {'<n_XY>':>7}  n_diff-histogram")
for N in [2,3,4,5]:
    Qs=Qstar(N); g=1.0/(Qs*1.02)                    # just above Q* (oscillating side)
    lam,v=slowest_osc_pair(N,1.0,g); nxy,h=ndiff_hist(N,v)
    print(f"{N:>2} {Qs:>8.4f} {2/Qs:>10.4f}  {lam.real:>8.4f}{lam.imag:+.4f}j   {nxy:>7.3f}  {h}")
