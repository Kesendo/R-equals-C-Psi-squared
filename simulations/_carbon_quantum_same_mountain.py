"""Nachsehen: is carbon really different? Compute the coherent<->incoherent threshold Q*(N)
for the XY chain + Z-dephasing (our quantum clock) and compare to the carbon (Frost/Hueckel)
values from FROST_CIRCLE_AS_THE_CLOCK_FACE.md (N=3 -> sqrt2, N=4 -> 1.879, N=5 -> 2.372).
Q* = where the slowest non-zero Liouvillian mode stops oscillating (the band-edge coherence
critically damps). N=2 is our EP base rung (the polyene layer starts at N>=3)."""
import numpy as np

I2=np.eye(2); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1]).astype(complex)

def site(op,l,N):
    m=np.array([[1.+0j]])
    for k in range(N): m=np.kron(m, op if k==l else I2)
    return m

def L_super(N,J,gamma):
    d=2**N; Id=np.eye(d)
    H=np.zeros((d,d),complex)
    for b in range(N-1):
        H += (J/2)*(site(X,b,N)@site(X,b+1,N)+site(Y,b,N)@site(Y,b+1,N))
    L=-1j*(np.kron(H,Id)-np.kron(Id,H.T))          # commutator [H,.]
    for l in range(N):                              # Z-dephasing: gamma*(Z rho Z - rho), gap=2gamma
        Zl=site(Z,l,N); L += gamma*(np.kron(Zl,Zl.conj()) - np.kron(Id,Id))
    return L

def slowest_imag(N,J,gamma):
    ev=np.linalg.eigvals(L_super(N,J,gamma))
    nz=ev[ev.real < -1e-7]
    gap=nz.real.max()                               # closest to 0 = slowest
    band=nz[np.abs(nz.real-gap)<1e-6]               # modes at the gap
    return np.abs(band.imag).max()                  # |Im| of the slowest = oscillating?

def Qstar(N,J=1.0):
    # bisect on gamma: oscillating (small gamma, high Q) vs frozen (large gamma, low Q)
    lo,hi=0.2,2.0                                    # gamma; Q=J/gamma in [0.5,5]
    if slowest_imag(N,J,hi)>1e-6: hi=5.0
    for _ in range(40):
        mid=0.5*(lo+hi)
        if slowest_imag(N,J,mid)>1e-6: lo=mid       # still oscillating -> need more gamma
        else: hi=mid
    g=0.5*(lo+hi); return J/g

carbon={3:np.sqrt(2),4:1.879,5:2.372}
print(f"{'N':>2} {'Q*(quantum)':>12} {'Q*(carbon)':>11} {'match':>7}")
for N in [2,3,4,5]:
    q=Qstar(N)
    c=carbon.get(N)
    cs=f"{c:.3f}" if c else "  (none)"
    mark = "EP base" if N==2 else ("YES" if c and abs(q-c)<0.02 else "diff")
    print(f"{N:>2} {q:>12.4f} {cs:>11} {mark:>7}")
print("\nsqrt2 =", round(np.sqrt(2),4), " (carbon's exact N=3 guess)")


# --- our label: is the slowest mode BELOW Q* dark ({I,Z}), and the rate ~ J^2/gamma? ---
def n_xy_of_mode(N, vec):
    """light content <n_XY> of a vectorized coherence operator (column-stacked rho)."""
    d=2**N; M=vec.reshape(d,d)                       # rho from vec(rho)
    w=np.abs(M)**2; w=w/ (w.sum()+1e-30)
    nxy=0.0
    for a in range(d):
        for b in range(d):
            if w[a,b]>1e-9:
                nxy += w[a,b]*bin(a^b).count("1")     # differing bits = n_XY of |a><b|
    return nxy

def slowest_vec(N,J,gamma):
    L=L_super(N,J,gamma); ev,evec=np.linalg.eig(L)
    mask=ev.real< -1e-7; idx=np.where(mask)[0]
    slow=idx[np.argmax(ev[idx].real)]
    return ev[slow], evec[:,slow]

print("\n--- the mechanism: slowest mode just BELOW Q* (low Q, frozen side) ---")
print(f"{'N':>2} {'Q':>5} {'slowest Re':>11} {'|Im|':>8} {'<n_XY>':>7} {'rate/(J^2/g)':>12}")
for N,Qs in [(3,1.414),(4,1.879),(5,2.372)]:
    Q=0.85*Qs; g=1.0/Q                               # just below the horizon
    lam,v=slowest_vec(N,1.0,g)
    nxy=n_xy_of_mode(N,v); rate=-lam.real
    print(f"{N:>2} {Q:>5.2f} {lam.real:>11.4f} {abs(lam.imag):>8.4f} {nxy:>7.3f} {rate/(1.0/g):>12.3f}")
