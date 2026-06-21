"""CRITICAL gate: is the SINGLE-EXCITATION walk EP (the HARDWARE object,
ep_onset_kingston_may2026) a GENUINE real-axis EP, or also diabolic?

The 5-agent review checked the (n,n+1) COHERENCE block (rate channels) → diabolic.
The hardware measured the SE walk (per-site populations, revival liftoff at Q≈1.5).
These are DIFFERENT blocks. The SE walk's overdamped→underdamped (critical-damping)
transition could be a genuine defective EP. If so, retracting 'the EP' would be an
OVER-CORRECTION — the EP is real (just mis-described by F86a's rate-channel toy).

Test: build the full L (N=3,4), sweep Q=J/γ0, and for the SE-relevant slow modes
check whether two REAL decay rates coalesce into a COMPLEX pair (critical damping),
and whether that coalescence is DEFECTIVE (eigenvector matrix singular / Petermann→∞
CONVERGENTLY) = genuine EP, vs a bounded/diabolic touching.
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
P1 = np.array([[0, 0], [0, 1]], dtype=complex)   # |1><1|, number operator on a site


def op(P, i, N):
    m = np.array([[1]], dtype=complex)
    for k in range(N):
        m = np.kron(m, P if k == i else I2)
    return m


def L_super(J, N, g0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for b in range(N - 1):
        H += 0.5 * J * (op(X, b, N) @ op(X, b + 1, N) + op(Y, b, N) @ op(Y, b + 1, N))
    L = -1j * (np.kron(np.eye(d), H) - np.kron(H.T, np.eye(d)))
    for l in range(N):
        Zl = op(Z, l, N)
        L += g0 * (np.kron(Zl, Zl) - np.eye(d * d))
    return L


def se_coherence_indices(N):
    """Liouville indices |a><b| with popcount(a)=1 and popcount(b)=1 (the SE
    population/coherence block that governs <n_l>(t) of a single-excitation walk)."""
    d = 2 ** N
    se = [i for i in range(d) if bin(i).count("1") == 1]
    idx = []
    for a in se:
        for b in se:
            idx.append(a * d + b)   # row-major vec(rho): index = a*d + b
    return np.array(idx)


def analyze(N, g0=0.05):
    print("=" * 72)
    print(f"N={N}: single-excitation (1,1)-popcount block - overdamped to revival EP?")
    idx = se_coherence_indices(N)
    print(f"   SE (1,1)-block dim = {len(idx)} (= ({N})^2)")
    prevcomplex = None
    for Q in [0.5, 0.8, 1.0, 1.2, 1.5, 1.879, 2.0, 2.5, 3.0, 5.0]:
        J = Q * g0
        L = L_super(J, N, g0)
        Lse = L[np.ix_(idx, idx)]               # the SE (1,1) coherence block
        w, VR = np.linalg.eig(Lse)
        # off-zero modes sorted by |Im| then Re; detect any complex (oscillatory) mode
        n_complex = np.sum(np.abs(w.imag) > 1e-9 * g0 / 0.05)
        # eigenvector-matrix condition number (defectiveness probe)
        condV = np.linalg.cond(VR)
        # max Petermann over modes
        VL = np.linalg.inv(VR).conj().T
        Ks = [(np.linalg.norm(VL[:, i]) * np.linalg.norm(VR[:, i]) /
               max(abs(np.vdot(VL[:, i], VR[:, i])), 1e-300)) ** 2 for i in range(len(w))]
        # slowest two nonzero real parts
        nz = w[np.abs(w) > 1e-9]
        slow = np.sort(nz.real)[-3:] / g0 if len(nz) else []
        print(f"   Q={Q:6}: #oscillatory(Im!=0)={n_complex:2d}  cond(V)={condV:9.2e}  "
              f"maxK={max(Ks):8.2e}  slow Re/g0={np.round(slow,3)}")
    print("   -> overdamped (all real) below the transition, oscillatory (complex pairs) above;")
    print("      a GENUINE EP shows cond(V)/maxK rising to a sharp CONVERGENT peak at the")
    print("      transition Q (defective). A bounded peak = diabolic / no EP.")


for N in (3, 4):
    analyze(N)
