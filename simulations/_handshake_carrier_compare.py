"""handshake R_k, carrier comparison (2026-06-19, Tom: compare both before building further).

The DefectDecoder paints Symphony's BondingMode = PURE |psi_1><psi_1| (single excitation, so <X_i>=<Y_i>=0,
the per-site purity P_i = 1/2(1+<Z_i>^2) is PURE POPULATION dynamics in the N^2-dim (1,1) sector).
But BondingMode.PairState ((|vac>+|psi_1>)/sqrt2, real vacuum-single COHERENCE) is documented as
"the painter's canvas". They are physically different readings. Does the choice change the per-site
reading f_i(b)=(alpha_i-1)/dJ QUALITATIVELY (different shape -> different R_k) or only by a SCALE?

Paint both carriers through the SAME canonical painter (st.full_L/purities/fit_alpha), same interior
defect, N=4,5. Report: the t=0 coherence content per site (confirming pure=0, pair>0), the f-profiles,
their strength (uniform level) and location (deviation), and the cross-carrier correlation of each.
"""
import importlib.util
import sys

import numpy as np
from scipy.linalg import eig

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "simulations")
sys.path.insert(0, "simulations/carbon")
from incompleteness_survivor import bonds


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


st = _load("stone", "simulations/_stone_survivor_alpha_closure.py")


def pure_psi1(N):
    return st.psi1_state(N)                                    # |psi_1><psi_1|, the decoder's actual carrier


def pair_state(N):
    d = 2 ** N
    psi = np.zeros(d, complex)
    psi[0] = 1.0
    for l in range(N):
        psi[1 << (N - 1 - l)] += np.sin(np.pi * (l + 1) / (N + 1))
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())                          # (|vac>+|psi_1>)/sqrt2, the documented canvas


def coherence_content(rho0, paulis, N):
    """t=0 per-site <X_i>^2+<Y_i>^2 (the single-site coherence the purity reads beyond populations)."""
    return np.array([np.trace(paulis[0][i] @ rho0).real ** 2 + np.trace(paulis[1][i] @ rho0).real ** 2
                     for i in range(N)])


def paint_f(N, J, g, dJ, bnds, defect, rho0):
    """Per-site f_i=(alpha_i-1)/dJ for one defect bond, canonical painter. Returns f, reliability, rate."""
    paulis = ([st.site(st.X, i, N) for i in range(N)],
              [st.site(st.Y, i, N) for i in range(N)],
              [st.site(st.Z, i, N) for i in range(N)])
    v0 = rho0.flatten()
    LA = st.full_L(N, J, g, bnds)
    wA, RA = eig(LA)
    cA = np.linalg.inv(RA) @ v0
    ds = st.dominant_slow(wA, cA)
    re = ds.real if ds is not None else -2.0 * g
    tg = np.linspace(0.0, 5.0 / max(abs(re), 1e-3), 60)
    PA = st.purities(wA, RA, cA, tg, paulis, N)
    Vbond = st.full_L(N, 1.0, 0.0, [defect]) - st.full_L(N, 0.0, 0.0, [defect])

    def f_at(frac):
        wB, RB = eig(LA + (dJ * frac) * Vbond)
        cB = np.linalg.inv(RB) @ v0
        PB = st.purities(wB, RB, cB, tg, paulis, N)
        a = np.array([st.fit_alpha(PA[i], PB[i], tg) for i in range(N)])
        return (a - 1) / (dJ * frac)

    f1, f2 = f_at(1.0), f_at(0.5)
    rel = (np.abs(f1) <= 10.0) & (np.abs(f1 - f2) <= 0.5)
    coh = coherence_content(rho0, paulis, N)
    return f1, rel, float(re), coh


def main():
    Q, J, dJ = 1.5, 1.0, 0.02
    g = 1.0 / Q
    print("=== handshake R_k: pure |psi_1> (decoder) vs PairState (documented canvas), same defect ===\n",
          flush=True)
    for N in (4, 5):
        bnds = bonds(N, "chain")
        defect = bnds[len(bnds) // 2]                          # interior bond
        fP, relP, reP, cohP = paint_f(N, J, g, dJ, bnds, defect, pure_psi1(N))
        fS, relS, reS, cohS = paint_f(N, J, g, dJ, bnds, defect, pair_state(N))
        both = relP & relS
        with np.printoptions(precision=3, suppress=True):
            print(f"--- N={N} chain, defect bond {defect} ---", flush=True)
            print(f"  pure |psi_1>:  rate={reP:+.3f}  coh-content(t0)={cohP}  (expect ~0)", flush=True)
            print(f"  PairState:     rate={reS:+.3f}  coh-content(t0)={cohS}  (expect >0)", flush=True)
            print(f"  f pure = {fP}", flush=True)
            print(f"  f pair = {fS}", flush=True)
            sP = fP[both].mean() if both.any() else float("nan")
            sS = fS[both].mean() if both.any() else float("nan")
            print(f"  strength (uniform level):  pure={sP:+.3f}   pair={sS:+.3f}", flush=True)
            devP, devS = fP - sP, fS - sS
            print(f"  location (deviation) pure = {devP}", flush=True)
            print(f"  location (deviation) pair = {devS}", flush=True)
        if both.sum() > 2:
            cf = float(np.corrcoef(fP[both], fS[both])[0, 1])
            cd = float(np.corrcoef((fP - sP)[both], (fS - sS)[both])[0, 1])
            print(f"  corr(f_pure, f_pair) = {cf:+.3f}   corr(location_pure, location_pair) = {cd:+.3f}",
                  flush=True)
            print(f"  => {'SAME shape (carrier choice = scale only)' if abs(cd) > 0.95 else 'DIFFERENT shape (carrier choice changes R_k qualitatively)'}",
                  flush=True)
        print(flush=True)


if __name__ == "__main__":
    main()
