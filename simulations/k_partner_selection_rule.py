"""Verify the K-partner selection rule <psi_N|V_b|psi_1>=0 (the carrier never leaks into its
K-partner under any bond defect) and the resulting location-channel rank N-2."""
import numpy as np
def psi(k,i,N): return np.sqrt(2/(N+1))*np.sin(np.pi*k*(i+1)/(N+1))
print(f"{'N':>2} {'<psi_N|V_b|psi_1> max':>21} {'psi_N==(-1)^i psi_1':>19} {'rank(loc chans k=2..N)':>22} {'N-2':>4}")
for N in range(3,9):
    # bond operator on the single-excitation sector: <psi_k|V_b|psi_1> = psi_k(b)psi_1(b+1)+psi_k(b+1)psi_1(b)
    M=np.zeros((N-1,N))                       # rows = bonds 0..N-2, cols = mode k=1..N
    for b in range(N-1):
        for k in range(1,N+1):
            M[b,k-1]=psi(k,b,N)*psi(1,b+1,N)+psi(k,b+1,N)*psi(1,b,N)
    sel = np.abs(M[:,N-1]).max()              # the selection rule: column k=N
    kpart = max(abs(psi(N,i,N)-((-1)**i)*psi(1,i,N)) for i in range(N))   # psi_N = K psi_1?
    rank_loc = np.linalg.matrix_rank(M[:,1:], tol=1e-9)   # location channels k=2..N (k=1 is the strength channel)
    print(f"{N:>2} {sel:>21.2e} {kpart:>19.2e} {rank_loc:>22} {N-2:>4}")
