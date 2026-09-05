"""WHICH of the two events the three decimals in FROST_CIRCLE_AS_THE_CLOCK_FACE.md belong to
(N=3 -> sqrt2, N=4 -> 1.879, N=5 -> 2.372). Those digits are not an independent carbon measurement:
that note's own script computed them from this same XY + Z-dephasing Liouvillian under the labels
beta <-> J, so this comparison has our numbers on both sides of it. What it can still decide, and
does, is which of OUR two thresholds the note tabulated.

WHICH threshold: the bisection below asks where the slowest non-zero mode of the FULL Liouvillian
stops oscillating, and that is the HANDOVER Q_h(N) -- the point where the {0,2}-coherence pair's
darker real branch reaches the floor Re = -2g and the band-edge survivor takes the gap. The
single-excitation EP Q*(N), where the pair coalesces, sits just ABOVE it: Q_h = 1.878541 / 2.372174
against Q* = 1.87874 / 2.37367 at N=4/5, equal only at N=2,3 where the pair is a clean 2x2. The
note's three decimals are Q_h's: at N=5 it carries 2.372 = Q_h, while the EP rounds to 2.374.
N=2 (Q*=1) is the EP itself and the note's ladder does not reach it."""
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

Qhandover = Qstar                                    # what the bisection measures, named for it
note={3:np.sqrt(2),4:1.879,5:2.372}                  # the three decimals FROST_CIRCLE tabulates
EP={4:1.87874,5:2.37367}                             # the SE coalescence, just above Q_h
print(f"{'N':>2} {'Q_h(quantum)':>13} {'the note':>8} {'carries them':>13} {'EP Q*':>9}")
measured={}
for N in [2,3,4,5]:
    q=Qhandover(N); measured[N]=q
    c=note.get(N)
    cs=f"{c:.3f}" if c else " (none)"
    mark = "EP base" if N==2 else ("YES" if c and round(q,3)==round(c,3) else "NO")
    eps=f"{EP[N]:.5f}" if N in EP else "    =Q_h"
    print(f"{N:>2} {q:>13.6f} {cs:>8} {mark:>13} {eps:>9}")
print("\nsqrt2 =", round(np.sqrt(2),4), " (the note's exact N=3 rung)")

# --- gate: the tabulated digits are the HANDOVER's, and the comparison can tell the two apart ------
# The note reports three decimals, so the comparison's law is its own half-ulp, 5e-4. Do NOT
# gate on round(): the N=4 rounding boundary sits 4e-5 from Q_h and would flip red on a correct
# value if the bisection moved in its sixth decimal.
HALF_ULP = 5e-4
for N in (3,4,5):
    d = abs(measured[N]-note[N])
    assert d < HALF_ULP, f"N={N}: Q_h={measured[N]:.6f} is {d:.2e} from the tabulated digits {note[N]:.3f}"
# Anti-vacuity, reading the MEASUREMENT and not two literals: at N=5 the handover is nearer the
# tabulated digits than the EP is, by nearly a factor 10, so the rung genuinely picks one of the two.
near, far = abs(measured[5]-note[5]), abs(EP[5]-note[5])
assert far > 5*near, f"the N=5 rung no longer separates Q_h ({near:.2e}) from the EP ({far:.2e})"
# At N=4 it cannot separate them: BOTH sit inside the same half-ulp. That rung is a match only, and
# Q_h uses 92% of the note's precision to be one, which is the honest size of the agreement there.
assert abs(EP[4]-note[4]) < HALF_ULP and abs(measured[4]-note[4]) < HALF_ULP,     "N=4 was expected to be the rung that cannot separate the handover from the EP"
print(f"[gate] the three tabulated decimals are the handover's: at N=5 it sits {near:.2e} away and the EP")
print(f"       {far:.2e}, a factor {far/near:.1f}; at N=4 both are inside the note's own 5e-4 and that")
print( "       rung separates nothing.")


# --- our label: is the slowest mode BELOW the handover dark ({I,Z}), and the rate ~ J^2/gamma? ---
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

print("\n--- the mechanism: slowest mode just BELOW the handover (low Q, frozen side) ---")
print(f"{'N':>2} {'Q':>5} {'slowest Re':>11} {'|Im|':>8} {'<n_XY>':>7} {'rate/(J^2/g)':>12}")
for N,Qs in [(3,1.414),(4,1.879),(5,2.372)]:
    Q=0.85*Qs; g=1.0/Q                               # just below the horizon
    lam,v=slowest_vec(N,1.0,g)
    nxy=n_xy_of_mode(N,v); rate=-lam.real
    print(f"{N:>2} {Q:>5.2f} {lam.real:>11.4f} {abs(lam.imag):>8.4f} {nxy:>7.3f} {rate/(1.0/g):>12.3f}")
