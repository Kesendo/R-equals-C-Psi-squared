"""USE R=CΨ², take 2 (gate from below). Take 1 (N=2 |++>) was REFUTED: stationary populations,
R=Ψ² exactly, C=1, nothing to read (the perfect-mirror limit). The readout only works where the
Born-deviation is nonzero -- where H rotates coherence into populations and dephasing intercepts.
The cleanest such case is the textbook decoherence experiment: a Rabi-driven qubit under Z-dephasing.

  Psi^2_1 = closed Born population of |1>   = <1|U|0><0|U^dag|1> = sin^2(Jt)   (unitary)
  R_1     = open population of |1>          = <1|rho_open(t)|1>                (Lindblad, Z-dephasing)
  C_1     = R_1 / Psi^2_1                   (the mirror quality / coherence factor for outcome |1>)
R=CΨ²: open outcome = closed Born x C. THE USE being gated: is C_1 a clean, monotone, INVERTIBLE
function of the decoherence gamma, so that measuring R_1 (and knowing the closed-system Psi^2_1)
RECOVERS gamma -- i.e. the deviation from Born IS a readout of the environment, from inside.
"""
import sys
import numpy as np
from scipy.linalg import expm

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1.0, -1.0]).astype(complex)
Id = np.eye(2, dtype=complex)
J = 1.0
H = J * X                                  # transverse field: rotates Z-population, makes coherence
rho0 = np.array([[1, 0], [0, 0]], dtype=complex)   # |0>
d = 2


def L(g):
    Lh = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    Ld = g * (np.kron(Z.T, Z) - np.kron(Id, Id))   # Z-dephasing
    return Lh + Ld


def R_open(g, t):
    rho = (expm(L(g) * t) @ rho0.reshape(-1, order="F")).reshape(d, d, order="F")
    return float(np.real(rho[1, 1]))           # P_1 open


def Psi2_closed(t):
    U = expm(-1j * H * t)
    rho = U @ rho0 @ U.conj().T
    return float(np.real(rho[1, 1]))           # P_1 closed = sin^2(Jt)


def main():
    t = 1.0
    P = Psi2_closed(t)
    print(f"=== R=CΨ² as a decoherence readout  (Rabi qubit, H=J·X, |0>, t={t}, J={J}) ===", flush=True)
    print(f"closed Born  Ψ²_|1> = sin²(Jt) = {P:.5f}\n", flush=True)

    print("[1] C_|1>(γ) = R_|1> / Ψ²_|1>  -- is it a clean, monotone, invertible readout of γ?", flush=True)
    print(f"    {'γ':>8} {'Q=J/γ':>8} {'R_|1> (open)':>14} {'C=R/Ψ²':>10} {'dev=1−C':>10}", flush=True)
    gs = [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6]
    Cs = []
    for g in gs:
        R = R_open(g, t) if g > 0 else P
        C = R / P
        Cs.append(C)
        print(f"    {g:>8.3g} {('inf' if g==0 else f'{J/g:.2f}'):>8} {R:>14.5f} {C:>10.5f} {1-C:>10.5f}", flush=True)
    mono = all(Cs[k] > Cs[k + 1] for k in range(len(Cs) - 1))
    print(f"    strictly monotone decreasing in γ (=> invertible)?  {mono}", flush=True)

    print("\n[2] THE USE — invert it: a 'measured' open count R_|1> at a HIDDEN γ; recover γ", flush=True)
    print("    (knowing only the closed-system Born Ψ²_|1> and the C(γ) calibration)", flush=True)
    cal_g = np.linspace(0.0, 2.0, 600)
    cal_C = np.array([(R_open(g, t) if g > 0 else P) / P for g in cal_g])
    ok = True
    for true_g in (0.13, 0.37, 0.61, 1.05):
        meas = R_open(true_g, t)                     # the measurement
        meas_C = meas / P
        rec_g = float(np.interp(meas_C, cal_C[::-1], cal_g[::-1]))   # cal_C decreasing -> reverse
        err = abs(rec_g - true_g)
        ok &= err < 5e-3
        print(f"    hidden γ={true_g:.3f}  →  measured R_|1>={meas:.5f}  (C={meas_C:.5f})  "
              f"→  recovered γ={rec_g:.3f}   (err {err:.1e})", flush=True)

    print(f"\nGATE: monotone-invertible readout? {mono}.   γ recovered from the Born-deviation? {ok}.", flush=True)
    if mono and ok:
        print("=> R=CΨ² USED: divide the measured outcome by the closed-system Born amplitude, read C,", flush=True)
        print("   invert, and the decoherence γ falls out. The deviation from Born IS the readout.", flush=True)
    else:
        print("=> not clean -- diagnose, do not dress it up.", flush=True)


if __name__ == "__main__":
    main()
