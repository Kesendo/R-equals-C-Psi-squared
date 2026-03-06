"""
STAR TOPOLOGY v2: RK4 Integration
==================================
Fixes Euler-stepping artifacts (purity > 1.0).
Runs key experiments from v1, validates qualitative findings.

See: experiments/STAR_TOPOLOGY_OBSERVERS.md
"""
import numpy as np

# ─── Pauli / helpers ───
I2 = np.eye(2, dtype=complex)
sx = np.array([[0,1],[1,0]], dtype=complex)
sy = np.array([[0,-1j],[1j,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)

def kron3(A, B, C):
    return np.kron(np.kron(A, B), C)

def star_hamiltonian(J_SA=1.0, J_SB=1.0):
    H = np.zeros((8,8), dtype=complex)
    for sigma in [sx, sy, sz]:
        H += J_SA * kron3(sigma, sigma, I2)
        H += J_SB * kron3(sigma, I2, sigma)
    return H

def lindblad_rhs(rho, H, L_ops):
    drho = -1j * (H @ rho - rho @ H)
    for L in L_ops:
        Ld = L.conj().T
        drho += L @ rho @ Ld - 0.5 * (Ld @ L @ rho + rho @ Ld @ L)
    return drho

def rk4_step(rho, H, L_ops, dt):
    k1 = lindblad_rhs(rho, H, L_ops)
    k2 = lindblad_rhs(rho + 0.5*dt*k1, H, L_ops)
    k3 = lindblad_rhs(rho + 0.5*dt*k2, H, L_ops)
    k4 = lindblad_rhs(rho + dt*k3, H, L_ops)
    rho_new = rho + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    rho_new = 0.5 * (rho_new + rho_new.conj().T)
    rho_new /= np.trace(rho_new).real
    return rho_new

def dephasing_ops(gamma_S, gamma_A, gamma_B):
    ops = []
    if gamma_S > 0: ops.append(np.sqrt(gamma_S) * kron3(sz, I2, I2))
    if gamma_A > 0: ops.append(np.sqrt(gamma_A) * kron3(I2, sz, I2))
    if gamma_B > 0: ops.append(np.sqrt(gamma_B) * kron3(I2, I2, sz))
    return ops

def ptrace(rho8, trace_out):
    """Partial trace: trace_out=0 traces S, 1 traces A, 2 traces B."""
    rho4 = np.zeros((4,4), dtype=complex)
    for i in range(4):
        for j in range(4):
            for k in range(2):
                if trace_out == 0:
                    a,b = divmod(i,2); ap,bp = divmod(j,2)
                    rho4[i,j] += rho8[k*4+a*2+b, k*4+ap*2+bp]
                elif trace_out == 1:
                    s,b = divmod(i,2); sp,bp = divmod(j,2)
                    rho4[i,j] += rho8[s*4+k*2+b, sp*4+k*2+bp]
                elif trace_out == 2:
                    s,a = divmod(i,2); sp,ap = divmod(j,2)
                    rho4[i,j] += rho8[s*4+a*2+k, sp*4+ap*2+k]
    return rho4

def l1_coh(rho):
    d = rho.shape[0]
    return sum(abs(rho[i,j]) for i in range(d) for j in range(d) if i!=j)

def psi_norm(rho):
    return l1_coh(rho) / (rho.shape[0] - 1)

def purity(rho):
    return np.trace(rho @ rho).real

def concurrence(rho4):
    sy2 = np.kron(sy, sy)
    rho_t = sy2 @ rho4.conj() @ sy2
    R = rho4 @ rho_t
    ev = np.sort(np.real(np.sqrt(np.maximum(np.linalg.eigvals(R), 0))))[::-1]
    return max(0, ev[0]-ev[1]-ev[2]-ev[3])

def measure_A_z(rho):
    """Projective Z-measurement on A (qubit 1)."""
    P0 = kron3(I2, np.array([[1,0],[0,0]],dtype=complex), I2)
    P1 = kron3(I2, np.array([[0,0],[0,1]],dtype=complex), I2)
    return P0 @ rho @ P0 + P1 @ rho @ P1

def make_state(name):
    if name == "GHZ":
        psi = np.zeros(8, dtype=complex); psi[0]=psi[7]=1/np.sqrt(2)
    elif name == "W":
        psi = np.zeros(8, dtype=complex); psi[1]=psi[2]=psi[4]=1/np.sqrt(3)
    elif name == "Bell_SA+B":
        bell = np.array([1,0,0,1], dtype=complex)/np.sqrt(2)
        psi = np.kron(bell, np.array([1,1], dtype=complex)/np.sqrt(2))
    elif name == "0++":
        psi = np.kron(np.kron(np.array([1,0],dtype=complex),
              np.array([1,1],dtype=complex)/np.sqrt(2)),
              np.array([1,1],dtype=complex)/np.sqrt(2))
    elif name == "+++":
        p = np.array([1,1],dtype=complex)/np.sqrt(2)
        psi = np.kron(np.kron(p,p),p)
    else:
        raise ValueError(name)
    return np.outer(psi, psi.conj())

def run(state_name, J_SA=1.0, J_SB=1.0, gS=0.05, gA=0.05, gB=0.05,
        dt=0.005, t_max=5.0, meas_at=None):
    """Run star topology simulation. Returns dict of time series."""
    H = star_hamiltonian(J_SA, J_SB)
    L = dephasing_ops(gS, gA, gB)
    rho = make_state(state_name)
    steps = int(t_max/dt)
    measured = False

    rec = {'t':[], 'SA_cpsi':[], 'SB_cpsi':[], 'AB_cpsi':[],
           'SA_R':[], 'SB_R':[], 'AB_R':[], 'R_sum':[], 'purity_full':[],
           'SA_conc':[], 'SB_conc':[], 'AB_conc':[]}

    for step in range(steps+1):
        t = step * dt
        if meas_at and not measured and t >= meas_at:
            rho = measure_A_z(rho); measured = True

        if step % 20 == 0:
            rec['t'].append(round(t,4))
            rec['purity_full'].append(round(purity(rho),6))
            for pair, tr_out, key in [('SA',2,'SA'), ('SB',1,'SB'), ('AB',0,'AB')]:
                rp = ptrace(rho, tr_out)
                p = psi_norm(rp)
                c = concurrence(rp)
                cpsi = c * p
                r = c * p**2
                rec[f'{key}_cpsi'].append(round(cpsi, 6))
                rec[f'{key}_R'].append(round(r, 6))
                rec[f'{key}_conc'].append(round(c, 6))
            rec['R_sum'].append(round(rec['SA_R'][-1]+rec['SB_R'][-1], 6))

        if step < steps:
            rho = rk4_step(rho, H, L, dt)

    return rec


if __name__ == "__main__":
    def crossings(times, vals, thr=0.25):
        cx = []
        for i in range(1, len(vals)):
            if vals[i-1]>thr and vals[i]<=thr:
                tc = times[i-1]+(thr-vals[i-1])/(vals[i]-vals[i-1])*(times[i]-times[i-1])
                cx.append(('down', round(tc,4)))
            elif vals[i-1]<thr and vals[i]>=thr:
                tc = times[i-1]+(thr-vals[i-1])/(vals[i]-vals[i-1])*(times[i]-times[i-1])
                cx.append(('up', round(tc,4)))
        return cx

    def report(name, r, extra=""):
        print(f"\n{'='*60}")
        print(f"  {name} {extra}")
        print(f"{'='*60}")
        print(f"  Purity: t=0 {r['purity_full'][0]:.6f} -> t=end {r['purity_full'][-1]:.6f}")
        for pair in ['SA','SB','AB']:
            cx = crossings(r['t'], r[f'{pair}_cpsi'])
            print(f"  {pair} CPsi crossings: {cx if cx else 'NEVER'}")
        print(f"  R_SA+R_SB: t=0={r['R_sum'][0]:.6f} max={max(r['R_sum']):.6f} end={r['R_sum'][-1]:.6f}")

    print("STAR TOPOLOGY v2 -- RK4, dt=0.005")
    for s in ["GHZ", "W", "Bell_SA+B", "0++", "+++"]:
        r = run(s); report(s, r)
