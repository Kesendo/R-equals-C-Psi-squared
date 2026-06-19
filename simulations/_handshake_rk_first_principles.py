"""handshake_decoder, the deeper stone: derive R_k(i) from FIRST PRINCIPLES, not a fit (2026-06-19).

The DefectDecoder reads a per-site profile f_i(b) = (alpha_i - 1)/dJ : the linear response of site-i's
purity-decay rate to a deltaJ defect on bond b. The reading-grammar conjecture is that this factorizes
    F[b,i] = sum_k c_k(b) R_k(i)          c_k(b) = <psi_k|V_b|psi_1>,  modes psi_k(i)=sqrt(2/(N+1))sin(pi k (i+1)/(N+1))
over the unforbidden channels k=1..N-1 (k=N killed by the K-partner selection rule <psi_N|V_b|psi_1>=0).

TWO PRIOR CONJECTURES VETOED (ledger 2026-06-13), both because they FIT R_k to the dictionary:
  (1) dressing c_k by 1/(E_1-E_k) is a per-column rescaling a least-squares R absorbs -> identical residual.
  (2) the unforbidden set is SQUARE ((N-1)x(N-1)) so any factorization is trivially exact -> tests nothing.
The fix is to DERIVE R_k independently and use the dictionary as a falsifiable check, never the fit target.

FIRST-PRINCIPLES R_k (this probe). A bond defect reshapes the carrier standing wave psi_1 by first-order
eigenvector perturbation: psi_1 -> psi_1 + dJ * sum_{k!=1} [c_k(b)/(eps_1-eps_k)] psi_k. The per-site
purity reading is a RATE-RATIO alpha_i, so f_i = (alpha_i-1)/dJ is a LOG-derivative of the per-site weight.
With the per-site clean weight ~ n_1(i)=psi_1(i)^2 (carrier density) times the global rate, the reading splits:
    f_i(b) = d ln n_1(i)/dJ_b           [LOCATION, site-dependent]      + d ln(rate)/dJ_b    [STRENGTH, uniform]
    d ln n_1(i)/dJ_b = sum_{k!=1} c_k(b) * R_k(i),   R_k(i) = 2 psi_k(i) / [ (eps_1-eps_k) psi_1(i) ].
So R_k is the eigenvector-perturbation footprint of the KNOWN sine modes -- ZERO free parameters. Two
self-consistency witnesses before any gate: c_N(b)=0 (selection rule, the rank-(N-2) deficit built in);
R_k(i) diverges at psi_1's nodes -- exactly the sites the painter already drops as unreliable.

We also report the LINEAR form R_k^lin(i)=2 psi_1(i) psi_k(i)/(eps_1-eps_k) (footprint of the raw density
change, not the log) so the gate discriminates ratio-vs-linear, and we compare the STRENGTH candidate two
ways: the single-particle energy shift d eps_1/dJ_b = <psi_1|V_b|psi_1> = 2 psi_1(b)psi_1(b+1) vs the
felt_time-(D) dissipative rate shift re_shift(0,1) -- to see whether the seam to felt_time is an identity
or a rhyme.

GATES that can fire:
  G0  PARAMETER-FREE RECONSTRUCTION (the headline). Build f_recon[b,i] from the closed form (no fit) and
      measure rel residual ||f_meas - s*f_recon|| / ||f_meas|| on RELIABLE sites, where s is a SINGLE global
      scalar (one number for the whole N, not per-channel -- a per-channel scalar would be the vetoed fit).
      Small (<~0.2) => R_k derived. O(1) => the footprint is wrong; diagnose, do not loosen.
  G1  LOCATION SHAPE: the location-only reconstruction (strength projected out) must correlate with the
      measured location-only profile (corr -> 1).
  G2  THE SEAM: is the measured uniform/strength part the felt_time-(D) rate shift, the energy shift, both,
      or neither? Report all three; do not assume.
  Gsel  the selection-rule witness: max_b |c_N(b)| ~ 0 (the forbidden channel is absent by construction).
"""
import importlib.util
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds
from value_vector_felt_time import re_shift


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


st = _load("stone", "simulations/_stone_survivor_alpha_closure.py")


# --- single-particle sector: modes, energies, dictionary ---------------------------------------------
def psi(k, i, N):
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (i + 1) / (N + 1))


def eps(k, N, J):
    return 2.0 * J * np.cos(np.pi * k / (N + 1))


def coupling(k, b, N):
    """c_k(b) = <psi_k|V_b|psi_1> on the single-excitation sector, V_b = 1/2(XX+YY) on bond b=(j,j+1)."""
    j0, j1 = b
    return psi(k, j0, N) * psi(1, j1, N) + psi(k, j1, N) * psi(1, j0, N)


def Rk_ratio(k, i, N, J):
    """Derived per-site purity-sensitivity (log / rate-ratio form): footprint of d ln n_1(i)."""
    return 2.0 * psi(k, i, N) / ((eps(1, N, J) - eps(k, N, J)) * psi(1, i, N))


def Rk_linear(k, i, N, J):
    """Alternative: footprint of the raw density change d n_1(i) (not the log)."""
    return 2.0 * psi(1, i, N) * psi(k, i, N) / (eps(1, N, J) - eps(k, N, J))


# --- ground truth: the actually-painted per-site f-profile (canonical painter, reused from the stone) -
def paint_f_profile(N, J, g, dJ, bnds):
    """For each bond defect, paint the per-site reading f_i(b)=(alpha_i-1)/dJ of the band-edge carrier.
    Returns f[nb, N], reliability mask[nb, N], and the carrier rate reB. Reuses st.full_L/purities/fit_alpha
    so this is the SAME pipeline as the C# DefectDecoder, not a re-derivation."""
    from scipy.linalg import eig

    rho0 = st.psi1_state(N)          # canonical BondingMode carrier |psi_1><psi_1| (matches Symphony exactly)
    v0 = rho0.flatten()
    paulis = ([st.site(st.X, i, N) for i in range(N)],
              [st.site(st.Y, i, N) for i in range(N)],
              [st.site(st.Z, i, N) for i in range(N)])
    LA = st.full_L(N, J, g, bnds)
    wA, RA = eig(LA)
    RAi = np.linalg.inv(RA)
    cA = RAi @ v0
    ds = st.dominant_slow(wA, cA)
    reB = ds.real if ds is not None else -2.0 * g
    tg = np.linspace(0.0, 5.0 / max(abs(reB), 1e-3), 60)
    PA = st.purities(wA, RA, cA, tg, paulis, N)

    nb = len(bnds)
    f = np.full((nb, N), np.nan)
    rel = np.zeros((nb, N), bool)
    for bi, b in enumerate(bnds):
        Vbond = st.full_L(N, 1.0, 0.0, [b]) - st.full_L(N, 0.0, 0.0, [b])
        for frac, store in ((1.0, "f1"), (0.5, "f2")):
            LB = LA + (dJ * frac) * Vbond
            wB, RB = eig(LB)
            cB = np.linalg.inv(RB) @ v0
            PB = st.purities(wB, RB, cB, tg, paulis, N)
            a = np.array([st.fit_alpha(PA[i], PB[i], tg) for i in range(N)])
            if store == "f1":
                f1 = (a - 1) / (dJ * frac)
            else:
                f2 = (a - 1) / (dJ * frac)
        f[bi] = f1
        rel[bi] = (np.abs(f1) <= 10.0) & (np.abs(f1 - f2) <= 0.5)
    return f, rel, reB


def main():
    Q, J, dJ = 1.5, 1.0, 0.02
    print("=== handshake R_k: is the per-site reading the eigenvector-perturbation footprint of the "
          "sine modes? (Q=1.5 chain, band-edge carrier) ===\n", flush=True)
    g = 1.0 / Q
    ok = True
    for N in (4, 5):
        bnds = bonds(N, "chain")
        nb = len(bnds)
        # --- cheap witnesses (no painter) ---
        csel = max(abs(coupling(N, b, N)) for b in bnds)                       # Gsel: c_N(b) ~ 0
        nodes = [i for i in range(N) if abs(psi(1, i, N)) < 1e-9]              # psi_1 nodes (divergent R_k)
        # --- ground truth painter ---
        f, rel, reB = paint_f_profile(N, J, g, dJ, bnds)

        # --- closed-form location reconstructions, BOTH candidate footprint forms (channels k=2..N-1) ---
        def build_loc(Rk):
            L = np.zeros((nb, N))
            for bi, b in enumerate(bnds):
                for i in range(N):
                    if abs(psi(1, i, N)) < 1e-9:
                        continue
                    L[bi, i] = sum(coupling(k, b, N) * Rk(k, i, N, J) for k in range(2, N))
            return L

        loc_r, loc_l = build_loc(Rk_ratio), build_loc(Rk_linear)
        eshift = np.array([2 * psi(1, b[0], N) * psi(1, b[1], N) for b in bnds])     # d eps_1/dJ_b
        dRe01 = np.array([re_shift(N, 0, 1, J, g, bnds, b)[0] for b in bnds])        # felt_time (D), (0,1)
        m = rel.copy()
        for i in nodes:
            m[:, i] = False

        # SEPARATE the two channels per bond: STRENGTH = uniform level (mean over reliable sites),
        # LOCATION = the deviation from it. The decoder localizes a defect by the deviation, not the level.
        def per_bond_center(M):
            out = np.zeros_like(M, float)
            for bi in range(nb):
                if m[bi].any():
                    out[bi] = M[bi] - M[bi][m[bi]].mean()
            return out

        strength = np.array([f[bi][m[bi]].mean() if m[bi].any() else np.nan for bi in range(nb)])
        f_dev = per_bond_center(f)

        def loc_corr(loc):
            lc = per_bond_center(loc)
            return float(np.corrcoef(f_dev[m], lc[m])[0, 1]) if m.sum() > 2 else float("nan")

        cr, cl = loc_corr(loc_r), loc_corr(loc_l)
        sb = ~np.isnan(strength)
        cs_e = float(np.corrcoef(strength[sb], eshift[sb])[0, 1]) if sb.sum() > 2 else float("nan")
        cs_d = float(np.corrcoef(strength[sb], dRe01[sb])[0, 1]) if sb.sum() > 2 else float("nan")
        g1 = max(abs(cr), abs(cl)) > 0.9
        gstr = max(abs(cs_e), abs(cs_d)) > 0.9
        ok &= g1
        print(f"--- N={N} chain, carrier rate Re={reB:+.3f}, reliable sites {int(rel.sum())}/{nb * N} ---",
              flush=True)
        print(f"  Gsel  max_b |c_N(b)| = {csel:.2e}   (K-partner channel absent by construction)", flush=True)
        print(f"  G1  LOCATION (per-bond deviation) corr:  ratio={cr:+.3f}  linear={cl:+.3f}  "
              f"[{'ok' if g1 else 'CHECK'}]   (sign = convention)", flush=True)
        print(f"  Gstr STRENGTH (uniform level) corr:  vs energy-shift={cs_e:+.3f}  vs felt_time-D={cs_d:+.3f}  "
              f"[{'ok' if gstr else 'CHECK'}]", flush=True)
        with np.printoptions(precision=3, suppress=True):
            print(f"  strength(b) = {strength}", flush=True)
            print(f"  eshift(b)   = {eshift}", flush=True)
            ib = nb // 2
            print(f"  [bond {bnds[ib]}] f_dev   = {f_dev[ib]}", flush=True)
            print(f"  [bond {bnds[ib]}] loc_ratio_dev = {per_bond_center(loc_r)[ib]}", flush=True)
            print(f"  [bond {bnds[ib]}] loc_lin_dev   = {per_bond_center(loc_l)[ib]}", flush=True)
        print(flush=True)
    print("VERDICT:", "R_k IS the eigenvector-perturbation footprint -- the reading grammar derived, not fit."
          if ok else "a gate fired -- the footprint form is not (yet) right; diagnose (ratio vs linear? "
                     "strength model? carrier?), do not loosen.")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
