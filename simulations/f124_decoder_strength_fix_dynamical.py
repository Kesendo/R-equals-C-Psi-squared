"""GATE-FIRST Stage B: at the DYNAMICAL level, does the carrier's own (strength-channel k=1)
response resolve the DefectDecoder's sign-location ambiguity that Stage 1 showed is dynamical
(living below the exact transition matrix M)?

The carrier pair-state (|vac> + |psi_1>)/sqrt2 lives in the (vacuum + single-excitation) sector,
which XY hopping AND Z-dephasing both preserve (excitation number conserved). So the dynamics is
EXACT in an (N+1)-dim sector: basis {vac, 1_0, .., 1_{N-1}}. rho is (N+1)x(N+1).
  - H_sector = blockdiag(0, h), h[a,a'] the XY hopping matrix (J on bonds; one bond shifted by dJ)
  - Z-dephasing: coherence (mu,nu) decays at rate 2*gamma*popcount(occ(mu) XOR occ(nu))
    (vac<->1_j at 2gamma, 1_i<->1_j at 4gamma), populations undephased.

For each single-bond J-defect we build:
  - LOCATION signature = the per-site purity-deviation profile [P_B(i,t)-P_A(i,t)]_i (the decoder's
    actual SITE-space observable; P_i = (1-p_i)^2 + p_i^2 + 2|v_i|^2, p_i=rho[1+i,1+i], v_i=rho[0,1+i])
  - STRENGTH scalar = the carrier's OWN response Delta<psi_1|rho|psi_1> (the k=1 channel; F124's
    M[b,1]=2 c_a c_{a+1} separates edge 0.289 from interior 0.577)

Mirror pair (N=5): bond 3 weakened (dJ=-0.02) vs bond 1 strengthened (dJ=+0.02). Gate-first:
  G1 - the location-only per-site profiles reproduce the DYNAMICAL ambiguity (a least-squares
       decode on the per-site dictionary does NOT cleanly separate them) -- the thing the exact
       linear M could not reproduce.
  G2 - the strength scalar SEPARATES bond 3 from bond 1 (sign and/or magnitude), tracking M[b,1].
  G3 - augmenting the dictionary with the strength scalar RESOLVES the decode (winner=truth, ratio
       jumps). If it does NOT, the gate fires (F124's strength column is not the dynamical fix).
"""
import sys
import numpy as np
from scipy.linalg import expm

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=4, suppress=True, linewidth=130)


def carrier(N, k=1):
    j = np.arange(N)
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (j + 1) / (N + 1))


def occ(mu, N):
    """site-occupation bitmask of sector basis state mu (0=vac, 1+j = excitation at site j)."""
    return 0 if mu == 0 else (1 << (mu - 1))


def hopping(N, J, defect_bond=None, dJ=0.0):
    h = np.zeros((N, N))
    for b in range(N - 1):
        Jb = J * (1.0 + dJ) if b == defect_bond else J
        h[b, b + 1] = h[b + 1, b] = Jb
    return h


def liouvillian(N, J, gamma, defect_bond=None, dJ=0.0):
    dim = N + 1
    H = np.zeros((dim, dim), complex)
    H[1:, 1:] = hopping(N, J, defect_bond, dJ)
    Id = np.eye(dim)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    # Z-dephasing: diagonal in the vec basis, rate -2 gamma popcount(occ(mu) ^ occ(nu))
    deph = np.zeros(dim * dim)
    for mu in range(dim):
        for nu in range(dim):
            r = bin(occ(mu, N) ^ occ(nu, N)).count("1")
            deph[mu * dim + nu] = -2.0 * gamma * r
    return L + np.diag(deph)


def pair_rho0(N, k=1):
    dim = N + 1
    psi = np.zeros(dim, complex)
    psi[0] = 1.0
    psi[1:] = carrier(N, k)
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def evolve(rho0, L, t):
    dim = int(round(np.sqrt(L.shape[0])))
    return (expm(L * t) @ rho0.reshape(-1)).reshape(dim, dim)


def site_purity_profile(rho, N):
    P = np.zeros(N)
    for i in range(N):
        p_i = rho[1 + i, 1 + i].real
        v_i = rho[0, 1 + i]
        P[i] = (1 - p_i) ** 2 + p_i ** 2 + 2 * abs(v_i) ** 2
    return P


def carrier_survival(rho, N, k=1):
    """<psi_1|rho|psi_1> over the single-excitation block = the carrier's own response."""
    c = carrier(N, k)
    blk = rho[1:, 1:]
    return (c @ blk @ c).real


def main():
    GATE = {"fired": []}

    def gate(name, cond, detail=""):
        flag = "ok " if cond else "GATE-FIRE"
        if not cond:
            GATE["fired"].append(name)
        print(f"   [{flag}] {name}" + (f"   {detail}" if detail else ""))

    N, J, gamma, dJ = 5, 1.0, 0.05, 0.02
    # integrate the signature over a time window (robust vs a single-t sampling artifact)
    ts = np.linspace(1.0, 12.0, 12)

    print("=" * 96)
    print(f"STAGE B (dynamical, N={N}, J={J}, gamma={gamma}, Q={J/gamma:.0f}, dJ={dJ}) - does the carrier")
    print("strength channel resolve the decoder's dynamical sign-location ambiguity?")
    print("=" * 96)

    L_clean = liouvillian(N, J, gamma)
    rho0 = pair_rho0(N)

    def signatures(defect_bond, dJ_signed):
        """time-integrated per-site purity-deviation profile + carrier-survival deviation."""
        Ld = liouvillian(N, J, gamma, defect_bond, dJ_signed)
        prof = np.zeros(N)
        surv = 0.0
        for t in ts:
            rc = evolve(rho0, L_clean, t)
            rd = evolve(rho0, Ld, t)
            prof += site_purity_profile(rd, N) - site_purity_profile(rc, N)
            surv += carrier_survival(rd, N) - carrier_survival(rc, N)
        return prof / len(ts), surv / len(ts)

    # dictionary: per-site profile for each bond at +dJ (the decoder calibrates at one sign)
    print("\nPer-site purity-deviation profiles (the decoder's dynamical dictionary), +dJ per bond:")
    dict_prof = []
    dict_surv = []
    for b in range(N - 1):
        p, s = signatures(b, +dJ)
        dict_prof.append(p); dict_surv.append(s)
        print(f"   bond {b}: profile={np.array2string(p, precision=4)}  strength={s:+.5f}")
    dict_prof = np.array(dict_prof)
    dict_surv = np.array(dict_surv)

    # --- the mirror pair: truth = bond 3 weakened (-dJ) ---
    truth = 3
    obs_prof, obs_surv = signatures(truth, -dJ)
    print(f"\nObserved (truth = bond {truth} WEAKENED, dJ=-{dJ}): profile={np.array2string(obs_prof, precision=4)}  strength={obs_surv:+.5f}")

    # decode location-only: fit obs profile to each bond's profile (allow a free signed amplitude)
    def decode(obs, dictionary):
        res = []
        for b in range(dictionary.shape[0]):
            f = dictionary[b]
            amp = (f @ obs) / (f @ f)
            res.append((np.linalg.norm(obs - amp * f), amp))
        order = sorted(range(len(res)), key=lambda b: res[b][0])
        win, run = order[0], order[1]
        ratio = res[run][0] / res[win][0] if res[win][0] > 1e-12 else np.inf
        return win, run, ratio, res

    win_loc, run_loc, ratio_loc, res_loc = decode(obs_prof, dict_prof)
    amp_win = res_loc[win_loc][1]
    print(f"\nLOCATION-ONLY decode: winner=bond{win_loc} (amp sign {np.sign(amp_win):+.0f}), runner-up=bond{run_loc}, "
          f"residual ratio={ratio_loc:.3f}")
    # G1: dynamical ambiguity present (wrong winner, OR right bond but wrong sign read, OR low ratio)
    sign_ok = (win_loc == truth and amp_win < 0)
    ambiguous = (win_loc != truth) or (not sign_ok) or (ratio_loc < 3.0)
    gate("G1: location-only shows the DYNAMICAL sign-location ambiguity (the exact M could not)",
         ambiguous, f"winner=bond{win_loc}, sign={'-' if amp_win<0 else '+'}, ratio={ratio_loc:.3f}")

    # G2: the strength scalar separates bond 3 (truth) from bond 1 (the confuser)
    print(f"\nStrength scalars: bond 1 (interior, +dJ) = {dict_surv[1]:+.5f}, bond 3 (edge, +dJ) = {dict_surv[3]:+.5f}, "
          f"observed (bond3, -dJ) = {obs_surv:+.5f}")
    sep = abs(dict_surv[1] - dict_surv[3]) > 0.2 * max(abs(dict_surv[1]), abs(dict_surv[3]))
    gate("G2: strength channel separates edge bond 3 from interior bond 1 (tracks M[b,1])", sep,
         f"|s1-s3|={abs(dict_surv[1]-dict_surv[3]):.5f}")

    # G3: augment each bond's signature with the strength scalar, re-decode -> resolves to truth, sign right
    aug_dict = np.column_stack([dict_prof, dict_surv])          # profile + strength
    aug_obs = np.append(obs_prof, obs_surv)
    win_aug, run_aug, ratio_aug, res_aug = decode(aug_obs, aug_dict)
    amp_aug = res_aug[win_aug][1]
    print(f"\nAUGMENTED (profile + strength) decode: winner=bond{win_aug} (amp sign {np.sign(amp_aug):+.0f}), "
          f"runner-up=bond{run_aug}, residual ratio={ratio_aug:.3f}")
    resolved = (win_aug == truth and amp_aug < 0 and ratio_aug >= max(ratio_loc, 3.0))
    gate("G3: adding the strength channel RESOLVES the decode (winner=truth bond 3, sign weakened, ratio up)",
         resolved, f"winner=bond{win_aug}, sign={'-' if amp_aug<0 else '+'}, ratio {ratio_loc:.2f}->{ratio_aug:.2f}")

    print("\n" + "=" * 96)
    if GATE["fired"]:
        print(f"STAGE B: {len(GATE['fired'])} GATE(S) FIRED -> {GATE['fired']}")
        print("Diagnose, do NOT loosen.")
    else:
        print("STAGE B: ALL GATES PASS -> F124's strength channel IS the decoder's de-blinding observable.")
        print("The location-only per-site dictionary inherits the dynamical sign-location ambiguity; adding")
        print("the carrier's own (strength-channel k=1) response - which F124 proves lifts lambda_min from 0")
        print("to E - resolves it. F124 hands the live handshake_decoder its fix: also read the carrier.")
    print("=" * 96)


if __name__ == "__main__":
    main()
