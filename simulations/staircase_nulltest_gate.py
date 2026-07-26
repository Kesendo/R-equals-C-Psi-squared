# Staircase null-test: PRE-REGISTRATION GATE v4 (flight 2, Amendment 2).
# v3 -> v4: B-block hardening after flight 1's B-BLOCK-INVALID (mid-batch T1
# telegraphing): per-spectator SINGLE-EXCITED preparations |010> / |001> (was |011>
# both-excited), each prep flown as two half-sets INTERLEAVED through the batch
# (was two end brackets); 142 circuits (was 128). Estimator, arms, grids, verdict
# machinery unchanged.
# v2 -> v3 (round-2 review folds): ZZ convention fixed (the measured -3.9 kHz is the
# CONDITIONAL SHIFT = 4*zeta_c, not the coefficient; v2 predicted 16*zeta_c), SIGMA_DELTA
# sqrt(2) overshoot fixed, B-block CLEAN self-void rate and min-kept-points floor printed.
# Simulates the FLOWN construction end-to-end (circuits -> counts -> committed estimator)
# and freezes the verdict bands. v1 -> v2 folds the three design-review rounds:
#   - estimator: MATCHED-PAIR LOG-RATIO fit ln(|C_s|/|C_00|) (v1's separate per-arm
#     exponential fits carry a -1.5 sigma systematic under Gaussian/quasi-static
#     dephasing envelopes; the ratio is exactly linear in t under the claim)
#   - noise model: adds a QUASI-STATIC (Gaussian-envelope) dephasing component on the
#     target qubit (Gauss-Hermite averaged), which v1's Markovian-only model missed
#   - B-block: delays extended to 600 us (z_inf identifiability to T1 = 500),
#     G-grid widened to [1/700, 1/50], TWO brackets (batch start + end) pooled,
#     with a bracket-consistency statistic (T1 telegraphing guard)
#   - verdicts: CENTERED on the MC means m_i (the residual B-block extrapolation
#     offset is an analysis artifact, removed by construction); bands from
#     N_MC = 2000 (~1.6% MC error), quoted to 2 significant figures
# Run: python simulations/staircase_nulltest_gate.py            (~ minutes)
# Day-of use: edit the DAY-OF INPUT block only, re-run, commit output as addendum.

import sys
import numpy as np
from scipy.linalg import expm

sys.path.insert(0, 'simulations')
from framework.lindblad import lindbladian_general
from framework.pauli import _site_op_kron, pauli_matrix

rng = np.random.default_rng(31415)
N = 3
D = 2**N

# ---- DAY-OF INPUT block (planning values; replace with day-of measurements) ----
T1_us   = np.array([322.0, 142.0, 237.0])
T2s_us  = np.array([136.0,  90.0, 110.0])
P_INF   = np.array([0.015, 0.020, 0.015])
EPS_RO  = np.array([[0.006, 0.012], [0.010, 0.018], [0.008, 0.015]])  # (e01, e10)/qubit
ZZ_SHIFT = -2 * np.pi * 0.0039         # measured CONDITIONAL-RAMSEY SHIFT (price-pair run 4)
ZZ      = ZZ_SHIFT / 4.0               # Hamiltonian coefficient zeta_c (shift = 4*zeta_c);
                                       # v2 wrongly planted the shift in the coefficient slot
                                       # and predicted 16*zeta_c; caught in round-2 review
DETUNE0 = 2 * np.pi * 0.030            # rad/us virtual detuning on q0
SHOTS   = 8192
N_MC    = 2000
QS_FRAC = 0.5                          # fraction of q0 pure dephasing that is quasi-static
T_REF   = 30.0                         # us, envelope-matching reference for sigma_delta

GAMMA_TOT = 1.0 / T1_us
G_DOWN = GAMMA_TOT * (1 - P_INF)
G_UP   = GAMMA_TOT * P_INF
GAMMA_Z_FULL = (1.0 / T2s_us - 1.0 / (2.0 * T1_us)) / 2.0     # c-op prefactor sqrt(g)*Z
GAMMA_Z = GAMMA_Z_FULL.copy()
GAMMA_Z[0] = (1 - QS_FRAC) * GAMMA_Z_FULL[0]                   # Markovian part on q0
G_QS = QS_FRAC * GAMMA_Z_FULL[0]
SIGMA_DELTA = np.sqrt(4.0 * G_QS / T_REF)
# a sqrt(g)*Z c-op decays the one-differing-site cell at 2g, so the removed Markovian half
# decays it at 2*G_QS; match exp(-sd^2 t^2/2) = exp(-2*G_QS*t) at T_REF -> sd^2 = 4 g_qs/T_REF
# (v2 had 8 g_qs/T_REF, a sqrt(2) overshoot; caught in round-2 review)

T_MAIN  = np.array([0, 1, 2, 4, 7, 11, 16, 24, 36, 54, 80.0])
T_DENSE = np.array([0, 0.5, 1, 1.5, 2, 3, 4.0])
T_ALL   = np.unique(np.concatenate([T_MAIN, T_DENSE]))
T_B     = np.array([0, 40, 80, 160, 240, 400, 600.0])

ARMS = {"00": (0, 0), "10": (1, 0), "01": (0, 1), "11": (1, 1)}

Z = pauli_matrix('Z')
SM = np.array([[0, 1], [0, 0]], dtype=complex); SP = SM.T.conj()

def build_L(delta0):
    H = ((DETUNE0 + delta0) / 2.0) * _site_op_kron(Z, 0, N)
    for (a, b) in [(0, 1), (1, 2)]:
        H = H + ZZ * (_site_op_kron(Z, a, N) @ _site_op_kron(Z, b, N))
    c_ops  = [np.sqrt(g) * _site_op_kron(Z,  l, N) for l, g in enumerate(GAMMA_Z) if g > 0]
    c_ops += [np.sqrt(g) * _site_op_kron(SM, l, N) for l, g in enumerate(G_DOWN)]
    c_ops += [np.sqrt(g) * _site_op_kron(SP, l, N) for l, g in enumerate(G_UP)]
    return lindbladian_general(H, c_ops)

def rho0_arm(s):
    psi = np.array([1, 1], dtype=complex) / np.sqrt(2)
    for b in s:
        psi = np.kron(psi, np.array([0, 1], dtype=complex) if b
                      else np.array([1, 0], dtype=complex))
    return np.outer(psi, psi.conj())

def rho0_B(prep):
    # Amendment 2 (flight 2): SINGLE-EXCITED B-block, one preparation per spectator.
    # prep = 1 -> |010> (near spectator excited), prep = 2 -> |001> (far spectator excited).
    # Flight 1's |011> both-excited block showed a structured residual on the short-T1
    # spectator; each reference is now read in exactly its arm's configuration.
    bits = (0, 1, 0) if prep == 1 else (0, 0, 1)
    psi = np.array([1.0 + 0j])
    for b in bits:
        psi = np.kron(psi, np.array([0, 1], dtype=complex) if b else np.array([1, 0], dtype=complex))
    return np.outer(psi, psi.conj())

HAD = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
SDG_H = HAD @ np.array([[1, 0], [0, -1j]], dtype=complex)

def confusion():
    M = np.eye(1)
    for q in range(N):
        Mq = np.array([[1 - EPS_RO[q][0], EPS_RO[q][1]], [EPS_RO[q][0], 1 - EPS_RO[q][1]]])
        M = np.kron(M, Mq)
    return M

M_CONF = confusion()

def probs_from_rho(rho, basis_q0):
    U = (HAD if basis_q0 == 'X' else SDG_H)
    for _ in range(N - 1):
        U = np.kron(U, np.eye(2))
    p = np.real(np.diag(U @ rho @ U.conj().T)).clip(0, None); p /= p.sum()
    return M_CONF @ p

def probs_z(rho):
    p = np.real(np.diag(rho)).clip(0, None); p /= p.sum()
    return M_CONF @ p

# quasi-static average: Gauss-Hermite nodes over delta0 ~ N(0, SIGMA_DELTA^2)
GH_X, GH_W = np.polynomial.hermite_e.hermegauss(15)
NODES = GH_X * SIGMA_DELTA
WEIGHTS = GH_W / GH_W.sum()

print("== gate v4: building probability tables (quasi-static averaged) ==")
P_ARM = {(a, b): {t: np.zeros(D) for t in T_ALL} for a in ARMS for b in ('X', 'Y')}
for delta, w in zip(NODES, WEIGHTS):
    L = build_L(delta)
    for a, s in ARMS.items():
        r0 = rho0_arm(s)
        for t in T_ALL:
            rho_t = (expm(L * t) @ r0.flatten('F')).reshape(D, D, order='F')
            for b in ('X', 'Y'):
                P_ARM[(a, b)][t] += w * probs_from_rho(rho_t, b)
L0 = build_L(0.0)   # B-block populations are detuning-independent
P_B = {}
for prep in (1, 2):
    r0b = rho0_B(prep).flatten('F')
    for t in T_B:
        P_B[(prep, t)] = probs_z((expm(L0 * t) @ r0b).reshape(D, D, order='F'))
e0 = np.zeros(D); e0[0] = 1
e1 = np.zeros(D); e1[-1] = 1
P_CAL0 = probs_z(np.outer(e0, e0).astype(complex))
P_CAL1 = probs_z(np.outer(e1, e1).astype(complex))

def draw(p):
    return rng.multinomial(SHOTS, p).astype(float)

def estimate_eps(cal0, cal1):
    eps = []
    for q in range(N):
        bit = N - 1 - q
        eps.append((sum(cal0[o] for o in range(D) if (o >> bit) & 1) / SHOTS,
                    sum(cal1[o] for o in range(D) if not (o >> bit) & 1) / SHOTS))
    return eps

def minv(eps_hat):
    Mi = np.eye(1)
    for q in range(N):
        Mq = np.array([[1 - eps_hat[q][0], eps_hat[q][1]],
                       [eps_hat[q][0], 1 - eps_hat[q][1]]])
        Mi = np.kron(Mi, np.linalg.inv(Mq))
    return Mi

def ols_slope(t, y):
    A = np.vstack([np.ones(len(t)), t]).T
    return np.linalg.lstsq(A, y, rcond=None)[0][1]

G_GRID = np.linspace(1.0 / 700, 1.0 / 50, 400)

def fit_B(z, t):
    best = None
    for G in G_GRID:
        e = np.exp(-G * t)
        A = np.vstack([1 - e, e]).T
        beta = np.linalg.lstsq(A, z, rcond=None)[0]
        res = ((A @ beta - z)**2).sum()
        if best is None or res < best[0]:
            best = (res, G, beta)
    res, G, (zinf, _z0) = best
    pinf = min(max((1 - zinf) / 2.0, 0.0), 0.5)
    return G, pinf, res, G in (G_GRID[0], G_GRID[-1])

def one_mc():
    eps_hat = estimate_eps(draw(P_CAL0), draw(P_CAL1))
    Mi = minv(eps_hat)
    cells = {}
    for a, s in ARMS.items():
        want = int(f"{s[0]}{s[1]}", 2)
        c = {}
        for t in T_ALL:
            est = {}
            for b in ('X', 'Y'):
                p = Mi @ (draw(P_ARM[(a, b)][t]) / SHOTS)
                v = 0.0
                for out in range(D):
                    if (out & (D // 2 - 1)) == want:
                        v += p[out] * (1 - 2 * (out >> (N - 1)))
                est[b] = v
            c[t] = est['X'] + 1j * est['Y']
        cells[a] = c
    # matched-pair log-ratio rate differences on the main grid
    delta = {}
    min_kept = len(T_MAIN)
    for a in ("10", "01", "11"):
        num = np.array([cells[a][t] for t in T_MAIN])
        den = np.array([cells["00"][t] for t in T_MAIN])
        m = (np.abs(num) > 0.02) & (np.abs(den) > 0.02)
        min_kept = min(min_kept, int(m.sum()))
        delta[a] = -ols_slope(T_MAIN[m], np.log(np.abs(num[m]) / np.abs(den[m])))
    # matched-pair phase-ratio frequency differences on the dense grid
    dfreq = {}
    for a in ("10", "01"):
        num = np.array([cells[a][t] for t in T_DENSE])
        den = np.array([cells["00"][t] for t in T_DENSE])
        dfreq[a] = ols_slope(T_DENSE, np.unwrap(np.angle(num / den)))
    # B-block, Amendment 2: per-spectator SINGLE-EXCITED preps, each run twice
    # (two half-sets, INTERLEAVED through the batch on hardware; the gate has no
    # time axis, so a half-set is statistically a bracket). Spectator q reads its
    # own prep's data only; per-half fits give the bracket-consistency statistic.
    brackets = {1: [], 2: []}
    zs = {1: [], 2: []}
    for _half in (0, 1):
        for q in (1, 2):
            zb = []
            for t in T_B:
                p = Mi @ (draw(P_B[(q, t)]) / SHOTS)
                bit = N - 1 - q
                zb.append(sum(p[o] * (1 - 2 * ((o >> bit) & 1)) for o in range(D)))
            G, pinf, _, _ = fit_B(np.array(zb), T_B)
            brackets[q].append(G * (1 - 2 * pinf))
            zs[q].append(zb)
    pred = {}
    bclean = {}
    for q in (1, 2):
        z_pool = np.array(zs[q]).mean(axis=0)
        G, pinf, res, railed = fit_B(z_pool, T_B)
        pred[q] = G * (1 - 2 * pinf)
        bclean[q] = dict(railed=railed, pinf=pinf, rms=np.sqrt(res / len(T_B)))
    return dict(
        r1=delta["10"] - pred[1],
        r2=delta["11"] - delta["10"] - delta["01"],
        r3=delta["01"] - pred[2],
        d10=delta["10"], d01=delta["01"],
        f10=dfreq["10"], f01=dfreq["01"],
        bc1=brackets[1][0] - brackets[1][1],
        bc2=brackets[2][0] - brackets[2][1],
        min_kept=min_kept,
        rail1=bclean[1]['railed'], pinf1=bclean[1]['pinf'], rms1=bclean[1]['rms'],
        rail2=bclean[2]['railed'], pinf2=bclean[2]['pinf'], rms2=bclean[2]['rms'],
    )

print(f"== gate v4: MC x {N_MC} over the flown construction ==")
runs = [one_mc() for _ in range(N_MC)]
def col(k): return np.array([r[k] for r in runs])

true_d10 = GAMMA_TOT[1] * (1 - 2 * P_INF[1])
true_d01 = GAMMA_TOT[2] * (1 - 2 * P_INF[2])
res = {}
for k, lbl in (('r1', 'P1'), ('r2', 'P2'), ('r3', 'P3')):
    v = col(k)
    res[k] = (v.mean(), v.std(ddof=1))
    print(f"  {lbl} {k}: mean {v.mean():+.2e}  sigma {v.std(ddof=1):.2e} /us"
          f"   (mean/sigma {abs(v.mean())/v.std(ddof=1):.2f})")
mc_rel = 1.0 / np.sqrt(2 * (N_MC - 1))
print(f"  MC relative error on each sigma: {100*mc_rel:.1f}%")
print(f"  Delta(10): {col('d10').mean():.6f} +- {col('d10').std(ddof=1):.6f}"
      f"  (pred {true_d10:.6f})  z_detect {true_d10/col('d10').std(ddof=1):.1f}")
print(f"  Delta(01): {col('d01').mean():.6f} +- {col('d01').std(ddof=1):.6f}"
      f"  (pred {true_d01:.6f})  z_detect {true_d01/col('d01').std(ddof=1):.1f}")
print(f"  dfreq(10): {col('f10').mean():+.5f} +- {col('f10').std(ddof=1):.5f} rad/us (pred {4*ZZ:+.5f})")
print(f"  dfreq(01): {col('f01').mean():+.5f} +- {col('f01').std(ddof=1):.5f} rad/us (pred 0)")
print(f"  bracket consistency sigma: bc1 {col('bc1').std(ddof=1):.2e}  bc2 {col('bc2').std(ddof=1):.2e} /us")
n_circ = len(T_ALL) * len(ARMS) * 2 + 2 * 2 * len(T_B) + 2   # Amendment 2: 2 preps x 2 halves
print(f"  circuits: {n_circ}  shots/circuit: {SHOTS}  total shots: {n_circ*SHOTS:,}")

marg = viol = 0
for r in runs:
    zvals = [abs(r[k] - res[k][0]) / res[k][1] for k in ('r1', 'r2', 'r3')]
    if any(zv > 3 for zv in zvals): viol += 1
    elif any(zv > 2 for zv in zvals): marg += 1
print(f"  H0-true family rates: P(>=1 MARGINAL) = {marg/N_MC:.3f}, P(>=1 VIOLATED) = {viol/N_MC:.3f}")

# B-block CLEAN operating characteristic under H0 (self-void rate) + estimator floor
sb1 = col('bc1').std(ddof=1); sb2 = col('bc2').std(ddof=1)
notclean = sum(1 for r in runs if (
    r['rail1'] or r['rail2'] or r['pinf1'] > 0.10 or r['pinf2'] > 0.10
    or r['rms1'] > 0.02 or r['rms2'] > 0.02
    or abs(r['bc1']) > 3 * sb1 or abs(r['bc2']) > 3 * sb2))
mk = min(r['min_kept'] for r in runs)
print(f"  P(B-block NOT clean | H0) = {notclean/N_MC:.3f}   (self-void rate)")
print(f"  min kept main-grid points over all arms/replicas: {mk}  "
      f"(day-of re-gate HARD abort if < 5)")

print("\nFROZEN CONSTANTS v4 (this gate, seed 31415, N_MC=2000; centers m_i and bands s_i;")
print("quote to 2 significant figures; the verdict statistic is r_i - m_i;")
print("m_i are a small day-of-rate-dependent curvature systematic, removed by centering):")
for k, lbl in (('r1', '1'), ('r2', '2'), ('r3', '3')):
    m, s = res[k]
    print(f"  m{lbl} = {m:+.2e}   s{lbl} = {s:.2e}   /us")
print(f"  sf_near = {col('f10').std(ddof=1):.2e}   sf_far = {col('f01').std(ddof=1):.2e}   rad/us")
print(f"  s_bracket1 = {col('bc1').std(ddof=1):.2e}   s_bracket2 = {col('bc2').std(ddof=1):.2e}   /us")

ok = (abs(res['r1'][0]) < res['r1'][1] and abs(res['r2'][0]) < res['r2'][1]
      and abs(res['r3'][0]) < res['r3'][1]
      and true_d10 / col('d10').std(ddof=1) > 5 and true_d01 / col('d01').std(ddof=1) > 5)
print("GATE:", "PASS (centered verdicts; both splittings > 5 sigma)" if ok else "FAIL")
