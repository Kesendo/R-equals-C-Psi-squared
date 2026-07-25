# Staircase null-test: feasibility simulation.
# Lineage note: this is the v0 feasibility scout (10/10 checks). The committed verdict
# construction lives in staircase_nulltest_gate.py (v2, matched-pair log-ratio estimator);
# this file additionally holds the collective/correlated-Z cancellation check (V3), which
# the pre-registration cites from HERE, not from the gate.
#
# The test: two coherence cells related by flipping a SPECTATOR site in both patterns
# (the staircase step) decay at rates whose DIFFERENCE contains no Z-dephasing at all,
# uneven or correlated; the difference is exactly the flipped site's amplitude-damping
# rate Gamma_m = 1/T1(m). ZZ crosstalk enters the FREQUENCY difference, not the rate.
#
# Derivation (this session): cell (A,B) decays at 2*sum_{l in A^B} gamma_l (Z part)
#   + 1/2 * sum_l Gamma_l * (n_l^A + n_l^B) (T1 part). Staircase step leaves A^B fixed
#   and raises n^A + n^B by 2 at site m => rate difference = Gamma_m, additive over m.
#
# Protocol shape (idle): conditional Ramsey. Target q0 in |+>, spectators in |s>,
# joint readout (X0/Y0 basis + Z on spectators), post-select spectator pattern s,
# fit slope of ln|C_s(t)|; verdict statistic = slope(s) - slope(00...0).
# Estimator lineage: IBM_CONCENTRATOR_RELOADED matched-pair slope difference;
# grid + shots lineage: PRICE_PAIR_HARDWARE_PREDICTION (8192 shots, 0..80 us grid).
#
# Numbers: ibm_kingston day-of snapshot (concentrator flight 2026-07-11) + measured
# Marrakesh ZZ (price-pair run 4, -3.9 kHz).
#
# Sections:
#   V1  estimator validation: exact-trajectory slope difference vs Sum Gamma_m
#   V2  dephasing blindness: wildly different uneven Tphi profiles -> same difference
#   V3  correlated Z noise: collective-Z Lindblad op -> difference unchanged
#   V4  ZZ: rate difference untouched, frequency difference reads the ZZ
#   V5  free motion on (XY, J on): difference still Sum Gamma_m (T1 perturbation)
#   V6  shot-noise Monte Carlo at hardware budget -> sigma(Delta), z-scores

import sys
import numpy as np

sys.path.insert(0, 'simulations')
import framework as fw
from framework.workflows.propagate_with_hardware_noise import propagate_with_hardware_noise
from framework.lindblad import lindbladian_general
from framework.pauli import _site_op_kron, pauli_matrix

rng = np.random.default_rng(20260725)

checks = []
def check(name, ok, detail=""):
    checks.append((name, bool(ok)))
    print(f"  [{'OK' if ok else 'FAIL'}] {name}  {detail}")

# ---------------------------------------------------------------- configuration
N = 3                      # q0 = target, q1/q2 = spectators
T1_us   = np.array([322.0, 142.0, 237.0])       # kingston day-of chain [109,108,107]
T2s_us  = np.array([136.0,  90.0, 110.0])       # T2* assumed (price-pair-like spread)
GAMMA   = 1.0 / T1_us                            # amplitude damping rates (1/us)
# pure dephasing: 1/Tphi = 1/T2* - 1/(2 T1); cell decays at 2*gamma_Z per site
inv_tphi = 1.0 / T2s_us - 1.0 / (2.0 * T1_us)
GAMMA_Z = inv_tphi / 2.0                         # c-op prefactor convention sqrt(g)*Z
ZZ_RAD_US = -2 * np.pi * 0.0039 / 4.0            # Hamiltonian coefficient zeta_c: the
# measured Marrakesh -3.9 kHz (price-pair run 4) is the CONDITIONAL SHIFT = 4*zeta_c
DETUNE = 2 * np.pi * np.array([0.030, 0.0, 0.0]) # 30 kHz virtual detuning on q0
T_GRID = np.array([0, 1, 2, 4, 7, 11, 16, 24, 36, 54, 80.0])  # price-pair grid, us

ARMS = {"00": (0, 0), "10": (1, 0), "11": (1, 1)}   # spectator patterns

def chain_idle():
    return fw.ChainSystem(N=N, gamma_0=0.0, J=0.0, topology='chain', H_type='xy')

def rho_init(s):
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    ket0 = np.array([1, 0], dtype=complex)
    ket1 = np.array([0, 1], dtype=complex)
    psi = plus
    for b in s:
        psi = np.kron(psi, ket1 if b else ket0)
    return np.outer(psi, psi.conj())

def cell_indices(s):
    # site 0 = leftmost kron factor (most significant bit)
    j = int('1' + ''.join(str(b) for b in s), 2)   # ket pattern with q0 = 1
    i = int('0' + ''.join(str(b) for b in s), 2)   # bra pattern with q0 = 0
    return i, j

def cell_traj(chain, s, t_grid, **noise):
    r0 = rho_init(s)
    i, j = cell_indices(s)
    out = np.empty(len(t_grid), dtype=complex)
    for k, t in enumerate(t_grid):
        rho_t = propagate_with_hardware_noise(chain, r0, t, **noise)
        out[k] = rho_t[j, i]     # <j|rho|i> = coherence |1s><0s| content read as rho[j,i]
    return out

def fit_rate(t, c, drop_below=1e-6):
    m = np.abs(c) > drop_below
    A = np.vstack([np.ones(m.sum()), t[m]]).T
    slope = np.linalg.lstsq(A, np.log(np.abs(c[m])), rcond=None)[0][1]
    return -slope    # decay rate, 1/us

def fit_freq(t, c):
    ph = np.unwrap(np.angle(c))
    A = np.vstack([np.ones(len(t)), t]).T
    return np.linalg.lstsq(A, ph, rcond=None)[0][1]  # rad/us

# sanity: initial cell = 0.5
_r = rho_init((1, 0)); _i, _j = cell_indices((1, 0))
check("V0 cell indexing (initial coherence = 0.5)", abs(_r[_j, _i] - 0.5) < 1e-12)

BASE = dict(T1_l=list(GAMMA), Tphi_l=list(GAMMA_Z), h_z_l=list(DETUNE / 2.0))
# NOTE h_z_l: H term h*Z_l gives cell frequency 2h per differing site -> detune/2

print("\n== V1: exact-trajectory slope differences vs Sum Gamma_m ==")
ch = chain_idle()
rates = {}
for name, s in ARMS.items():
    rates[name] = fit_rate(T_GRID, cell_traj(ch, s, T_GRID, **BASE))
d10 = rates["10"] - rates["00"]
d11 = rates["11"] - rates["00"]
check("V1a flip q1: Delta = Gamma_1", abs(d10 - GAMMA[1]) < 5e-5,
      f"Delta={d10:.6f} vs Gamma_1={GAMMA[1]:.6f} /us")
check("V1b flip q1+q2: Delta = Gamma_1+Gamma_2 (additivity)",
      abs(d11 - (GAMMA[1] + GAMMA[2])) < 1e-4,
      f"Delta={d11:.6f} vs sum={GAMMA[1]+GAMMA[2]:.6f} /us")

print("\n== V2: dephasing blindness (three wildly different uneven Tphi profiles) ==")
deltas = []
for trial in range(3):
    gz = rng.uniform(0.0, 0.05, size=N)   # up to ~10x the physical dephasing
    nz = dict(BASE); nz['Tphi_l'] = list(gz)
    r00 = fit_rate(T_GRID, cell_traj(ch, ARMS["00"], T_GRID, **nz))
    r10 = fit_rate(T_GRID, cell_traj(ch, ARMS["10"], T_GRID, **nz))
    deltas.append(r10 - r00)
spread = max(deltas) - min(deltas)
check("V2 Delta independent of Tphi profile", spread < 1e-6,
      f"spread over profiles = {spread:.2e} /us")

print("\n== V3: correlated (collective) Z noise cancels in the difference ==")
# build by hand: H idle + local Z + T1 + one collective sqrt(g)*(sum c_l Z_l)
Z = pauli_matrix('Z'); SM = np.array([[0, 1], [0, 0]], dtype=complex)
H0 = sum(0.5 * DETUNE[l] / (2 * np.pi) * 0 for l in range(N))  # placeholder 0
H0 = np.zeros((2**N, 2**N), dtype=complex)
for l in range(N):
    H0 += (DETUNE[l] / 2.0) * _site_op_kron(Z, l, N)
c_ops = [np.sqrt(g) * _site_op_kron(Z, l, N) for l, g in enumerate(GAMMA_Z)]
c_ops += [np.sqrt(g) * _site_op_kron(SM, l, N) for l, g in enumerate(GAMMA)]
cvec = rng.uniform(-1, 1, size=N)
c_ops_coll = c_ops + [0.15 * sum(cvec[l] * _site_op_kron(Z, l, N) for l in range(N))]

def cell_traj_manual(s, c_ops_use):
    from scipy.linalg import expm
    L = lindbladian_general(H0, c_ops_use)
    r0 = rho_init(s).flatten('F')
    i, j = cell_indices(s)
    out = np.empty(len(T_GRID), dtype=complex)
    for k, t in enumerate(T_GRID):
        rho_t = (expm(L * t) @ r0).reshape(2**N, 2**N, order='F')
        out[k] = rho_t[j, i]
    return out

d_loc = fit_rate(T_GRID, cell_traj_manual(ARMS["10"], c_ops)) \
      - fit_rate(T_GRID, cell_traj_manual(ARMS["00"], c_ops))
d_col = fit_rate(T_GRID, cell_traj_manual(ARMS["10"], c_ops_coll)) \
      - fit_rate(T_GRID, cell_traj_manual(ARMS["00"], c_ops_coll))
check("V3 collective Z shifts both rates, difference unchanged",
      abs(d_col - d_loc) < 1e-6, f"|d_col - d_loc| = {abs(d_col-d_loc):.2e} /us")

print("\n== V4: ZZ goes into frequency difference, not rate difference ==")
nz = dict(BASE); nz['J_zz'] = ZZ_RAD_US
r00 = cell_traj(ch, ARMS["00"], T_GRID, **nz)
r10 = cell_traj(ch, ARMS["10"], T_GRID, **nz)
d_rate_zz = fit_rate(T_GRID, r10) - fit_rate(T_GRID, r00)
check("V4a rate difference still = Gamma_1 with ZZ on",
      abs(d_rate_zz - GAMMA[1]) < 5e-5, f"Delta={d_rate_zz:.6f} /us")
# frequency fit needs a DENSE grid: the hardware delay grid undersamples the
# 30 kHz fringe (steps up to 26 us vs phase period ~33 us) and unwrap aliases.
T_DENSE = np.linspace(0.0, 12.0, 25)
f00 = fit_freq(T_DENSE, cell_traj(ch, ARMS["00"], T_DENSE, **nz))
f10 = fit_freq(T_DENSE, cell_traj(ch, ARMS["10"], T_DENSE, **nz))
d_freq = f10 - f00
# q0's cell frequency shift when neighbour q1 flips: |Delta_freq| = 4 |J_zz|
check("V4b frequency difference reads the ZZ (dense grid)",
      abs(abs(d_freq) - 4 * abs(ZZ_RAD_US)) < 2e-3,
      f"|d_freq|={abs(d_freq):.5f} vs 4|Jzz|={4*abs(ZZ_RAD_US):.5f} rad/us")

print("\n== V5: free motion on (XY, J>0): difference still Sum Gamma (perturbative) ==")
ch_j = fw.ChainSystem(N=N, gamma_0=0.0, J=0.02, topology='chain', H_type='xy')
nz = dict(BASE)
r00 = fit_rate(T_GRID, cell_traj(ch_j, ARMS["00"], T_GRID, **nz))
r11 = fit_rate(T_GRID, cell_traj(ch_j, ARMS["11"], T_GRID, **nz))
d_j = r11 - r00
check("V5 XY on (J=0.02 rad/us): |Delta - (G1+G2)| small",
      abs(d_j - (GAMMA[1] + GAMMA[2])) < 0.2 * (GAMMA[1] + GAMMA[2]),
      f"Delta={d_j:.6f} vs {GAMMA[1]+GAMMA[2]:.6f} /us")

print("\n== V6: shot-noise Monte Carlo at hardware budget ==")
SHOTS = 8192          # per circuit (per delay, per arm, per basis), price-pair budget
EPS_RO = 0.01         # symmetric readout error per qubit
N_MC = 200

def measurement_probs(rho_t, basis):
    # basis: 'X' or 'Y' on q0, Z on spectators; returns prob over 2^N outcomes
    if basis == 'X':
        U0 = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)      # H
    else:
        U0 = np.array([[1, -1j], [1, 1j]], dtype=complex) / np.sqrt(2)    # H*Sdg
    U = U0
    for _ in range(N - 1):
        U = np.kron(U, np.eye(2))
    p = np.real(np.diag(U @ rho_t @ U.conj().T)).clip(0, None)
    p = p / p.sum()
    # readout confusion, independent per qubit
    M1 = np.array([[1 - EPS_RO, EPS_RO], [EPS_RO, 1 - EPS_RO]])
    M = M1
    for _ in range(N - 1):
        M = np.kron(M, M1)
    return M @ p

# precompute exact rho(t) per arm/delay (with ZZ on, full realistic model)
noise_full = dict(BASE); noise_full['J_zz'] = ZZ_RAD_US
rho_store = {a: [propagate_with_hardware_noise(chain_idle(), rho_init(s), t, **noise_full)
                 for t in T_GRID] for a, s in ARMS.items()}
probs = {(a, b): [measurement_probs(r, b) for r in rho_store[a]]
         for a in ARMS for b in ('X', 'Y')}

def mc_estimate_delta(arm):
    s = ARMS[arm]; s0 = ARMS["00"]
    def arm_rate(a, spat):
        want = int(''.join(str(b) for b in spat), 2)
        cs = np.empty(len(T_GRID), dtype=complex)
        for k in range(len(T_GRID)):
            est = {}
            for b in ('X', 'Y'):
                counts = rng.multinomial(SHOTS, probs[(a, b)][k])
                # post-select spectator outcome == prepared pattern
                sel_p = sel_m = 0
                for out in range(2**N):
                    if (out & (2**(N-1) - 1)) == want:      # low bits = spectators
                        if out >> (N - 1):  sel_m += counts[out]
                        else:               sel_p += counts[out]
                # RAW estimator: divide by TOTAL shots, not by the kept count.
                # Renormalizing by the post-selected count divides out exactly
                # e^{-Gamma_m t} (the spectator's own survival), which is the
                # signal: the normalized conditional coherence is Gamma_1-BLIND
                # (that cancellation is why standard conditional-Ramsey reads
                # only frequency). Post-selection stays in the numerator only.
                est[b] = (sel_p - sel_m) / SHOTS
            cs[k] = est['X'] + 1j * est['Y']
        return fit_rate(T_GRID, cs, drop_below=0.02)
    return arm_rate(arm, s) - arm_rate("00", s0)

d_mc = np.array([mc_estimate_delta("10") for _ in range(N_MC)])
sig = d_mc.std(ddof=1)
bias = d_mc.mean() - GAMMA[1]
z_detect = GAMMA[1] / sig
# what fraction of the spectator's dephasing rate could hide in the difference?
frac_null = sig / (2 * GAMMA_Z[1])
print(f"  shots/circuit={SHOTS}, delays={len(T_GRID)}, MC={N_MC}")
print(f"  Delta_hat = {d_mc.mean():.6f} +- {sig:.6f} /us   (true Gamma_1 = {GAMMA[1]:.6f})")
print(f"  bias = {bias:.2e} /us,  z(detect Gamma_1) = {z_detect:.1f} sigma")
print(f"  null certification: sigma(Delta) = {frac_null:.3f} x spectator pure-dephasing rate")
check("V6a bias small vs sigma", abs(bias) < 0.5 * sig, f"bias/sigma={abs(bias)/sig:.2f}")
check("V6b Gamma_1 detectable at >5 sigma", z_detect > 5, f"z={z_detect:.1f}")

n_ok = sum(1 for _, ok in checks if ok)
print(f"\n{'ALL GREEN' if n_ok == len(checks) else 'FAILURES'}: {n_ok}/{len(checks)} checks")
