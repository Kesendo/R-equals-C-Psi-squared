"""USE R=CΨ² (gate from below): is it an invertible decoherence readout?

TAKE 1 (N=2 |++>): REFUTED -- the PERFECT-MIRROR NULL. |++> has stationary populations on the N=2
Heisenberg chain, so R = Psi^2 EXACTLY, C = 1 for every outcome, dev ~ 1e-16: nothing to read, the
inversion fails. The readout only has content where H rotates coherence into populations and dephasing
intercepts it -- see _rcpsi_readout_rabi.py (Rabi qubit) and _rcpsi_readout_f94.py (F94 substrate, the
Tier-1-grounded take), and experiments/RCPSI_DECOHERENCE_READOUT.md. This file is kept as the documented
null; its concluding line is now GATED on the (failing) monotonicity check so it reports honestly.

The grounded meaning (agents, 2026-06-19): R=CΨ² is the generalized Born rule R_i = C_i·Ψ_i².
  Psi^2_i = closed-system Born population  <i|U rho0 U^dag|i>   (unitary; what a closed system gives)
  R_i     = open-system population         <i|rho_open(t)|i>     (Lindblad; what the open system gives)
  C_i     = R_i / Psi^2_i                   (the mirror quality / effective coupling for outcome i)
Born is the C=1 (perfect-mirror) limit. The DECOMPOSITION R=C*Psi^2 is definitional; the CONTENT --
the claim being gated -- is that C_i is a clean, MONOTONE, INVERTIBLE function of the decoherence,
so that measuring R (and knowing the closed-system Psi^2) READS the decoherence from inside.
That is the 'use': the deviation from Born (the '3%') is the signal, not noise.
"""
import sys
import numpy as np
from scipy.linalg import expm

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)


def op(P, l, N):
    m = np.array([[1.0 + 0j]])
    for k in range(N):
        m = np.kron(m, P if k == l else I2)
    return m


N = 2
J = 1.0
bonds = [(0, 1)]
H = sum(J * (op(X, i, N) @ op(X, j, N) + op(Y, i, N) @ op(Y, j, N) + op(Z, i, N) @ op(Z, j, N))
        for (i, j) in bonds)
d = 2 ** N
plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
psi = plus
for _ in range(N - 1):
    psi = np.kron(psi, plus)
rho0 = np.outer(psi, psi.conj())
Id = np.eye(d, dtype=complex)
labels = [format(i, f"0{N}b") for i in range(d)]


def liouvillian(g):
    Lh = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    Ld = np.zeros((d * d, d * d), dtype=complex)
    for l in range(N):
        Zl = op(Z, l, N)
        Ld += g * (np.kron(Zl.T, Zl) - np.kron(Id, Id))
    return Lh + Ld


def R_open(g, t):
    rho = (expm(liouvillian(g) * t) @ rho0.reshape(-1, order="F")).reshape(d, d, order="F")
    return np.real(np.diag(rho))


def Psi2_closed(t):
    U = expm(-1j * H * t)
    rho = U @ rho0 @ U.conj().T
    return np.real(np.diag(rho))


def main():
    t = 0.7
    Psi2 = Psi2_closed(t)
    dom = int(np.argmax(Psi2))                       # the dominant outcome (cleanest C)
    print(f"=== R=CΨ² as a decoherence readout  (N=2 Heisenberg, |++>, t={t}, J={J}) ===\n", flush=True)

    # (1) the decomposition at one gamma: R = C * Psi^2, with C != 1
    g0 = 0.3
    R0 = R_open(g0, t)
    print(f"[1] At γ={g0} (Q=J/γ={J/g0:.2f}, K=γt={g0*t:.3f}):  R_i = C_i · Ψ²_i", flush=True)
    print(f"    {'outcome':>8} {'Ψ² (Born)':>12} {'R (open)':>12} {'C=R/Ψ²':>10}", flush=True)
    for i in range(d):
        if Psi2[i] > 1e-9:
            print(f"    {labels[i]:>8} {Psi2[i]:>12.5f} {R0[i]:>12.5f} {R0[i]/Psi2[i]:>10.5f}", flush=True)

    # (2) the content: is C_dom(γ) a clean, monotone, invertible readout?  And how does the deviation scale?
    print(f"\n[2] dominant outcome |{labels[dom]}>:  C(γ) and the Born-deviation vs decoherence", flush=True)
    print(f"    {'γ':>8} {'Q=J/γ':>8} {'K=γt':>8} {'C=R/Ψ²':>10} {'dev=C−1':>12}", flush=True)
    gs = [0.02, 0.05, 0.1, 0.2, 0.4, 0.8]
    devs = []
    for g in gs:
        C = R_open(g, t)[dom] / Psi2[dom]
        devs.append(C - 1.0)
        print(f"    {g:>8.3g} {J/g:>8.2f} {g*t:>8.3f} {C:>10.5f} {C-1.0:>12.3e}", flush=True)
    monotone = all(devs[k] != devs[k+1] for k in range(len(devs)-1)) and \
        (all(devs[k] < devs[k+1] for k in range(len(devs)-1)) or
         all(devs[k] > devs[k+1] for k in range(len(devs)-1)))
    print(f"    monotone in γ (=> invertible)? {monotone}", flush=True)

    # scaling: fit dev ~ a * γ^p in the small-γ regime (the F94 structure is dev~Q^2 K^3 = J^2 γ t^3 ~ γ^1)
    g_small = np.array([0.02, 0.05, 0.1, 0.2])
    d_small = np.array([abs(R_open(g, t)[dom] / Psi2[dom] - 1.0) for g in g_small])
    p = np.polyfit(np.log(g_small), np.log(d_small), 1)[0]
    print(f"    small-γ power law:  |dev| ∝ γ^{p:.2f}   (F94 form dev∝Q²K³=J²γt³ ⇒ p=1 at fixed t)", flush=True)

    # (3) THE USE -- inversion: a 'measured' open outcome at a HIDDEN γ; recover γ from R and Ψ² alone
    print(f"\n[3] THE USE — invert it: given a measured R (hidden γ) and the closed-system Ψ², recover γ", flush=True)
    cal_g = np.linspace(0.01, 0.9, 400)                       # the C(γ) calibration curve (Ψ² known)
    cal_C = np.array([R_open(g, t)[dom] / Psi2[dom] for g in cal_g])
    for true_g in (0.13, 0.37, 0.61):
        meas_C = R_open(true_g, t)[dom] / Psi2[dom]           # the 'measurement': R_dom / Ψ²_dom
        rec_g = float(np.interp(meas_C, cal_C[::-1], cal_g[::-1]) if cal_C[0] > cal_C[-1]
                      else np.interp(meas_C, cal_C, cal_g))
        print(f"    hidden γ={true_g:.3f}  →  measured C={meas_C:.5f}  →  recovered γ={rec_g:.3f}   "
              f"(err {abs(rec_g-true_g):.1e})", flush=True)

    print("", flush=True)
    if monotone:
        print("=> R=CΨ² used: divide the measured outcome by the closed-system Born amplitude, read C,", flush=True)
        print("   invert the calibration, recover the decoherence. The deviation IS the readout.", flush=True)
    else:
        print("=> REFUTED here (perfect-mirror null): C≡1, the populations are stationary, the deviation is", flush=True)
        print("   ~1e-16 and not monotone, so there is nothing to invert. The readout needs a substrate where", flush=True)
        print("   H rotates coherence↔population (Rabi: _rcpsi_readout_rabi.py; F94: _rcpsi_readout_f94.py).", flush=True)


if __name__ == "__main__":
    main()
